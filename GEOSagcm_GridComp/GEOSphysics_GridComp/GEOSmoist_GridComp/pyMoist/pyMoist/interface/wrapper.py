"""
Wraps pyMoist for GEOS interface use.
"""

import dataclasses
import os
from typing import Callable, Optional, get_type_hints

import numpy as np
from MAPLpyish import CVoidPointer
from mpi4py import MPI

import pyMoist.GFDL_1M as pyGFDL_1M
from ndsl import (
    CompilationConfig,
    CubedSphereCommunicator,
    CubedSpherePartitioner,
    DaceConfig,
    GridIndexing,
    MPIComm,
    NullComm,
    PerformanceCollector,
    QuantityFactory,
    StencilConfig,
    StencilFactory,
    SubtileGridSizer,
    TilePartitioner,
)
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py_utils import backend_is_fortran_aligned, is_gpu_backend
from ndsl.dsl.typing import get_precision
from ndsl.logging import ndsl_log_on_rank_0
from ndsl.optional_imports import cupy as cp
from pyMoist.aer_activation import AerActivation
from pyMoist.interface import (
    InterfaceTransferType,
    MAPLMemoryRepository,
    MemorySpace,
    StencilBackendCompilerOverride,
)
from pyMoist.interface.flags import GFDL1MFlags, MoistFlags
from pyMoist.interface.gfdl_1m import GFDL1MInterface
from pyMoist.UW.compute_uwshcu import ComputeUwshcuInv, UWConfiguration


@dataclasses.dataclass
class MAPLStates:
    import_: CVoidPointer
    export: CVoidPointer
    mapl_comp: CVoidPointer


class GEOSPyMoistWrapper:
    def __init__(
        self,
        mapl_states: MAPLStates,
        flags: MoistFlags,
        backend="numpy",
        fortran_mem_space: MemorySpace = MemorySpace.CPU,
    ) -> None:
        # Look for an override to run on a single node
        single_rank_override = int(os.getenv("GEOS_PYFV3_SINGLE_RANK_OVERRIDE", -1))
        comm = MPIComm()
        if single_rank_override >= 0:
            comm = NullComm(single_rank_override, 6, 42)

        self.backend = backend
        self.flags = flags
        layout = (self.flags.layout_x, self.flags.layout_y)

        # Make a custom performance collector for the GEOS wrapper
        self.perf_collector = PerformanceCollector("GEOS Moist", comm)
        partitioner = CubedSpherePartitioner(TilePartitioner(layout))
        self.communicator = CubedSphereCommunicator(
            comm,
            partitioner,
            timer=self.perf_collector.timestep_timer,
        )
        sizer = SubtileGridSizer.from_tile_params(
            nx_tile=self.flags.npx * self.flags.layout_x,
            ny_tile=self.flags.npy * self.flags.layout_y,
            nz=self.flags.npz,
            n_halo=0,
            layout=layout,
            tile_partitioner=partitioner.tile,
            tile_rank=self.communicator.tile.rank,
        )
        self.quantity_factory = QuantityFactory.from_backend(sizer=sizer, backend=backend)

        self.stencil_config = StencilConfig(
            compilation_config=CompilationConfig(backend=backend, rebuild=False, validate_args=True),
        )

        # Build a DaCeConfig for orchestration.
        # This and all orchestration code are transparent when outside
        # configuration deactivate orchestration
        self.stencil_config.dace_config = DaceConfig(
            communicator=self.communicator,
            backend=self.stencil_config.backend,
            tile_nx=self.flags.npx * self.flags.layout_x,
            tile_nz=self.flags.npz,
            time=True,
        )
        self._is_orchestrated = self.stencil_config.dace_config.is_dace_orchestrated()

        # TODO: Orchestrate all code called from this function

        self._grid_indexing = GridIndexing.from_sizer_and_communicator(sizer=sizer, comm=self.communicator)
        self.stencil_factory = StencilFactory(config=self.stencil_config, grid_indexing=self._grid_indexing)

        # Figure out the interface mode
        tmp_quantity = self.quantity_factory.empty([X_DIM, Y_DIM, Z_DIM], units="")
        default_3D_memory_desc = (tmp_quantity.data.shape, tmp_quantity.data.strides)
        if fortran_mem_space != MemorySpace.CPU:
            raise NotImplementedError("Interface cannot stream Fortran memory resident on GPU")
        if is_gpu_backend(backend):
            self._interface_type = InterfaceTransferType.CPU_TO_GPU_TO_CPU
        else:
            if backend_is_fortran_aligned(backend):
                # This is Fortran layout - we can Map the memory
                self._interface_type = InterfaceTransferType.CPU_MAP
            else:
                # All other layout have to copy the data in/out of Fortran layout
                self._interface_type = InterfaceTransferType.CPU_COPY
        del tmp_quantity

        # Feedback information
        device_ordinal_info = (
            f"  Device PCI bus id: {cp.cuda.Device(0).pci_bus_id}" if is_gpu_backend(backend) else "N/A"
        )
        MPS_pipe_directory = os.getenv("CUDA_MPS_PIPE_DIRECTORY", None)
        MPS_is_on = (
            MPS_pipe_directory is not None
            and is_gpu_backend(backend)
            and os.path.exists(f"{MPS_pipe_directory}/log")
        )
        ndsl_log_on_rank_0.info(
            "pyMoist <> GEOS wrapper initialized (Rank 0):\n"
            f"         Bridge : {self._interface_type.name}\n"
            f"        Backend : {backend}\n"
            f"      Precision : {get_precision()} bit\n"
            f"  Orchestration : {self._is_orchestrated}\n"
            f"          Sizer : {sizer.nx}x{sizer.ny}x{sizer.nz}"
            f"(halo: {sizer.n_halo})\n"
            f" Strides for 3D : {default_3D_memory_desc[1]}\n"
            f"     Device ord : {device_ordinal_info}\n"
            f"     Nvidia MPS : {MPS_is_on}\n"
        )

        # Timer result dict
        self._timings: dict[str, list[float]] = {}

        # JIT system for the component of Moist

        # AER
        self._aer_activation: AerActivation | None = None

        # GFDL 1M
        self._gfdl_1m_interface = GFDL1MInterface(
            self.quantity_factory,
            self.stencil_factory,
            self._interface_type,
        )

        # UW
        self._UW_shallow_convection: Optional[ComputeUwshcuInv] = None
        self.UW_config: UWConfiguration | None = None

        # Initalize MAPL Memory Respositories
        self._mapl_import = MAPLMemoryRepository(
            mapl_states.import_,
            self.quantity_factory,
        )
        self._mapl_export = MAPLMemoryRepository(
            mapl_states.export,
            self.quantity_factory,
        )
        self._mapl_comp = MAPLMemoryRepository(
            mapl_states.mapl_comp,
            self.quantity_factory,
        )

    @property
    def aer_activation(self) -> Callable:
        if not self._aer_activation:
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                self.stencil_config.dace_config,
            ):
                self._aer_activation = AerActivation(
                    stencil_factory=self.stencil_factory,
                    quantity_factory=self.quantity_factory,
                    n_modes=self.flags.n_modes,
                    USE_AERSOL_NN=True,
                )
        return self._aer_activation

    def init_gfdl_1m_configuration(
        self,
        flags: GFDL1MFlags,
        internal_state: CVoidPointer,
    ):
        gfdl_1m_mapl_internal = MAPLMemoryRepository(
            internal_state,
            self.quantity_factory,
        )
        # Get namelist/non-constant parameters passed through the interface
        # Dev NOTE:
        # - data will come as default python POD, e.g. not respecting precision
        # - get_type_hints makes sure we resolve the annotation of the dataclass
        # - we force cast to make sure precision is correct and code generation
        #   respects the intended precision
        upper_case_dict = {}
        for name, type_ in get_type_hints(GFDL1MFlags).items():
            # Don't bring the magic number
            if name.upper() != "MN_123456789":
                upper_case_dict[name] = type_(getattr(flags, name.upper()))

        # Get remaining required parameters from MAPL
        HYDROSTATIC = self._mapl_comp.get_resource("HYDROSTATIC:", bool, default=True)
        PHYS_HYDROSTATIC = self._mapl_comp.get_resource("PHYS_HYDROSTATIC:", bool, default=True)
        MELTFRZ = self._mapl_comp.get_resource("MELTFRZ:", bool, default=True)
        TURNRHCRIT = self._mapl_comp.get_resource("TURNRHCRIT:", np.float32, default=-9999.0)
        PDF_SHAPE = self._mapl_comp.get_resource("PDFSHAPE:", np.int32, default=1)
        # PDF_SHAPE = 1
        ANV_ICEFALL = self._mapl_comp.get_resource("ANV_ICEFALL:", np.float32, default=1.0)
        LS_ICEFALL = self._mapl_comp.get_resource("LS_ICEFALL:", np.float32, default=1.0)
        LIQ_RADII_PARAM = self._mapl_comp.get_resource("LIQ_RADII_PARAM:", np.int32, default=2)
        # LIQ_RADII_PARAM = 2
        ICE_RADII_PARAM = self._mapl_comp.get_resource("ICE_RADII_PARAM:", np.int32, default=1)
        # ICE_RADII_PARAM = 1
        FAC_RI = self._mapl_comp.get_resource("FAC_RI:", np.float32, default=1.0)
        MIN_RI = self._mapl_comp.get_resource("MIN_RI:", np.float32, default=5.0e-6)
        MAX_RI = self._mapl_comp.get_resource("MAX_RI:", np.float32, default=100.0e-6)
        FAC_RL = self._mapl_comp.get_resource("FAC_RL:", np.float32, default=1.0)
        MIN_RL = self._mapl_comp.get_resource("MIN_RL:", np.float32, default=2.5e-6)
        MAX_RL = self._mapl_comp.get_resource("MAX_RL:", np.float32, default=60.0e-6)
        CCW_EVAP_EFF = self._mapl_comp.get_resource("CCW_EVAP_EFF:", np.float32, default=1e-2)
        CCI_EVAP_EFF = self._mapl_comp.get_resource("CCI_EVAP_EFF:", np.float32, default=1e-2)

        config = pyGFDL_1M.GFDL1MConfig(
            LHYDROSTATIC=HYDROSTATIC,  # type: ignore # bool are stupid in numpy
            LPHYS_HYDROSTATIC=PHYS_HYDROSTATIC,  # type: ignore
            LMELTFRZ=MELTFRZ,  # type: ignore
            TURNRHCRIT_PARAM=TURNRHCRIT,
            PDFSHAPE=PDF_SHAPE,
            ANV_ICEFALL=ANV_ICEFALL,
            LS_ICEFALL=LS_ICEFALL,
            LIQ_RADII_PARAM=LIQ_RADII_PARAM,
            ICE_RADII_PARAM=ICE_RADII_PARAM,
            FAC_RI=FAC_RI,
            MIN_RI=MIN_RI,
            MAX_RI=MAX_RI,
            FAC_RL=FAC_RL,
            MIN_RL=MIN_RL,
            MAX_RL=MAX_RL,
            CCW_EVAP_EFF=CCW_EVAP_EFF,
            CCI_EVAP_EFF=CCI_EVAP_EFF,
            **upper_case_dict,
        )

        self._gfdl_1m_interface.init(
            config=config,
            mapl_internal=gfdl_1m_mapl_internal,
            mapl_import=self._mapl_import,
            mapl_export=self._mapl_export,
        )

    def GFDL_1M_Microphysics(self):
        self._gfdl_1m_interface.run(self._timings)

    @property
    def UW_shallow_convection(self) -> Callable:
        if self._UW_shallow_convection is None and self.UW_config is not None:
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                self.stencil_config.dace_config,
            ):
                self._UW_shallow_convection = ComputeUwshcuInv(
                    stencil_factory=self.stencil_factory,
                    quantity_factory=self.quantity_factory,
                    UW_config=self.UW_config,
                )
        return self._UW_shallow_convection

    def init_UW_configuration(
        self,
        NCNST,
        k0,
        windsrcavg,
    ):
        self.UW_config = UWConfiguration(
            NCNST=NCNST,
            k0=k0,
            windsrcavg=windsrcavg,
        )

    def finalize(self):
        import json

        rank = MPI.COMM_WORLD.Get_rank()
        with open(f"pymoist_timings_r{rank}.json", "w") as f:
            json.dump(self._timings, f, indent=4)

        with open(f"internal_orchestration_r{rank}.json", "w") as f:
            json.dump(self.stencil_config.dace_config.performance_collector.times_per_step, f, indent=4)

        with open(f"internal_gt4py_timings_r{rank}.json", "w") as f:
            json.dump(self.stencil_factory.timing_collector.exec_info, f, indent=4)
