"""
Wraps pyMoist for GEOS interface use.
"""

import dataclasses
import enum
import logging
import os
from typing import Callable, Optional

import numpy as np
from gt4py.cartesian.config import build_settings as gt_build_settings
from MAPLpyish import CVoidPointer
from mpi4py import MPI

from ndsl import (
    CompilationConfig,
    CubedSphereCommunicator,
    CubedSpherePartitioner,
    DaceConfig,
    DaCeOrchestration,
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
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.dace.build import set_distributed_caches
from ndsl.dsl.gt4py_utils import is_gpu_backend
from ndsl.dsl.typing import get_precision
from ndsl.logging import ndsl_log
from ndsl.optional_imports import cupy as cp
from pyMoist.aer_activation import AerActivation
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.interface.cuda_profiler import TimedCUDAProfiler
from pyMoist.interface.flags import GFDL1MFlags, MoistFlags
from pyMoist.interface.mapl.memory_factory import MAPLManagedMemory, MAPLMemoryRepository


class MemorySpace(enum.Enum):
    HOST = 0
    DEVICE = 1


class StencilBackendCompilerOverride:
    """Override the NDSL global stencil JIT to allow for 9-rank build
    on any setup.

    This is a workaround that requires to know _exactly_ when build is happening.
    Using this as a context manager, we leverage the DaCe build system to override
    the name and build the 9 codepaths required- while every other rank wait.

    This should be removed when we refactor the GT JIT to distribute building
    much more efficiently
    """

    def __init__(self, comm: MPI.Intracomm, config: DaceConfig):
        self.comm = comm
        self.config = config

        # Orchestration or mono-node is not concerned
        self.no_op = self.config.is_dace_orchestrated() or self.comm.Get_size() == 1

        # We abuse the DaCe build system
        if not self.no_op:
            config._orchestrate = DaCeOrchestration.Build
            set_distributed_caches(config)
            config._orchestrate = DaCeOrchestration.Python

        # We remove warnings from the stencils compiling when in critical and/or
        # error
        if ndsl_log.level > logging.WARNING:
            gt_build_settings["extra_compile_args"]["cxx"].append("-w")
            gt_build_settings["extra_compile_args"]["cuda"].append("-w")

    def __enter__(self):
        if self.no_op:
            return
        if self.config.do_compile:
            ndsl_log.info(f"Stencil backend compiles on {self.comm.Get_rank()}")
        else:
            ndsl_log.info(f"Stencil backend waits on {self.comm.Get_rank()}")
            self.comm.Barrier()
            ndsl_log.info(f"Stencil backend released on {self.comm.Get_rank()}")

    def __exit__(self, type, value, traceback):
        if self.no_op:
            return
        if not self.config.do_compile:
            ndsl_log.info(f"Stencil backend read cache on {self.comm.Get_rank()}")
        else:
            ndsl_log.info(
                f"Stencil backend was compiled on {self.comm.Get_rank()} \
                    now waiting for other ranks"
            )
            self.comm.Barrier()
        ndsl_log.info(f"Rank {self.comm.Get_rank()} ready for execution")


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
        fortran_mem_space: MemorySpace = MemorySpace.HOST,
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
            extra_dim_lengths={},
            layout=layout,
            tile_partitioner=partitioner.tile,
            tile_rank=self.communicator.tile.rank,
        )
        self.quantity_factory = QuantityFactory.from_backend(
            sizer=sizer, backend=backend
        )
        self.nmodes_quantity_factory = AerActivation.make_nmodes_quantity_factory(
            self.quantity_factory
        )

        self.stencil_config = StencilConfig(
            compilation_config=CompilationConfig(
                backend=backend, rebuild=False, validate_args=True
            ),
        )

        # Build a DaCeConfig for orchestration.
        # This and all orchestration code are transparent when outside
        # configuration deactivate orchestration
        self.stencil_config.dace_config = DaceConfig(
            communicator=self.communicator,
            backend=self.stencil_config.backend,
            tile_nx=self.flags.npx * self.flags.layout_x,
            tile_nz=self.flags.npz,
        )
        self._is_orchestrated = self.stencil_config.dace_config.is_dace_orchestrated()

        # TODO: Orchestrate all code called from this function

        self._grid_indexing = GridIndexing.from_sizer_and_communicator(
            sizer=sizer, comm=self.communicator
        )
        self.stencil_factory = StencilFactory(
            config=self.stencil_config, grid_indexing=self._grid_indexing
        )

        self._fortran_mem_space = fortran_mem_space
        self._pace_mem_space = (
            MemorySpace.DEVICE if is_gpu_backend(backend) else MemorySpace.HOST
        )

        # Feedback information
        device_ordinal_info = (
            f"  Device PCI bus id: {cp.cuda.Device(0).pci_bus_id}\n"
            if is_gpu_backend(backend)
            else "N/A"
        )
        MPS_pipe_directory = os.getenv("CUDA_MPS_PIPE_DIRECTORY", None)
        MPS_is_on = (
            MPS_pipe_directory is not None
            and is_gpu_backend(backend)
            and os.path.exists(f"{MPS_pipe_directory}/log")
        )
        ndsl_log.info(
            "pyMoist <> GEOS wrapper initialized: \n"
            f"         bridge : {self._fortran_mem_space}"
            f" > {self._pace_mem_space}\n"
            f"        backend : {backend}\n"
            f"          float : {get_precision()}bit\n"
            f"  orchestration : {self._is_orchestrated}\n"
            f"          sizer : {sizer.nx}x{sizer.ny}x{sizer.nz}"
            f"(halo: {sizer.n_halo})\n"
            f"     Device ord : {device_ordinal_info}\n"
            f"     Nvidia MPS : {MPS_is_on}"
        )

        # Timer result dict
        self._timings = {}

        # JIT system for the component of Moist
        self._aer_activation: Optional[AerActivation] = None
        self._GFDL_1M_driver: Optional[MicrophysicsDriver] = None
        self._GFDL_1M_ready: Optional[bool] = False

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
    def driver(self) -> Callable:
        if not self._GFDL_1M_driver:
            if self.GFDL_1M_config is None:
                raise RuntimeError("GFDL_1M configuration not initalized")
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                self.stencil_config.dace_config,
            ):
                self._GFDL_1M_driver = MicrophysicsDriver(
                    self.stencil_factory,
                    self.quantity_factory,
                    self.GFDL_1M_config,
                )
        return self._GFDL_1M_driver

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

    def make_nmmodes_quantity(self, data):
        qty = self.nmodes_quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            "n/a",
        )
        qty.view[:, :, :, :] = qty.np.asarray(data[:, :, :, :])
        return qty

    def init_gfdl_1m_configuration(
        self,
        flags: GFDL1MFlags,
        internal_state: CVoidPointer,
    ):
        from pyMoist.GFDL_1M.config import GFDL1MConfig
        from pyMoist.GFDL_1M.GFDL_1M import GFDL1M

        self._mapl_internal = MAPLMemoryRepository(
            internal_state,
            self.quantity_factory,
        )
        # Get namelist/non-constant parameters passed through the interface
        upper_case_dict = {}
        for field in dataclasses.fields(GFDL1MFlags):
            # Don't bring the magic number
            if field.name.upper() != "MN_123456789":
                upper_case_dict[field.name] = getattr(flags, field.name.upper())

        # Get remaining required parameters from MAPL
        HYDROSTATIC = self._mapl_comp.get_resource("HYDROSTATIC:", bool, default=True)
        PHYS_HYDROSTATIC = self._mapl_comp.get_resource(
            "PHYS_HYDROSTATIC:", bool, default=True
        )
        MELTFRZ = self._mapl_comp.get_resource("MELTFRZ:", bool, default=True)
        TURNRHCRIT = self._mapl_comp.get_resource(
            "TURNRHCRIT:", np.float32, default=-9999.0
        )
        # PDF_SHAPE = self._mapl_comp.get_resource("PDFSHAPE:", np.int32, default=1)
        PDF_SHAPE = 1
        ANV_ICEFALL = self._mapl_comp.get_resource(
            "ANV_ICEFALL:", np.float32, default=1.0
        )
        LS_ICEFALL = self._mapl_comp.get_resource(
            "LS_ICEFALL:", np.float32, default=1.0
        )
        # LIQ_RADII_PARAM = self._mapl_comp.get_resource("LIQ_RADII_PARAM:", np.int32, default=2)
        LIQ_RADII_PARAM = 2
        # ICE_RADII_PARAM = self._mapl_comp.get_resource("ICE_RADII_PARAM:", np.int32, default=1)
        ICE_RADII_PARAM = 1
        FAC_RI = self._mapl_comp.get_resource("FAC_RI:", np.float32, default=1.0)
        MIN_RI = self._mapl_comp.get_resource("MIN_RI:", np.float32, default=5.0e-6)
        MAX_RI = self._mapl_comp.get_resource("MAX_RI:", np.float32, default=100.0e-6)
        FAC_RL = self._mapl_comp.get_resource("FAC_RL:", np.float32, default=1.0)
        MIN_RL = self._mapl_comp.get_resource("MIN_RL:", np.float32, default=2.5e-6)
        MAX_RL = self._mapl_comp.get_resource("MAX_RL:", np.float32, default=60.0e-6)
        CCW_EVAP_EFF = self._mapl_comp.get_resource(
            "CCW_EVAP_EFF:", np.float32, default=60.0e-6
        )
        CCI_EVAP_EFF = self._mapl_comp.get_resource(
            "CCI_EVAP_EFF:", np.float32, default=60.0e-6
        )

        self.GFDL_1M_config = GFDL1MConfig(
            HYDROSTATIC=HYDROSTATIC,
            PHYS_HYDROSTATIC=PHYS_HYDROSTATIC,
            MELTFRZ=MELTFRZ,
            TURNRHCRIT_PARAM=TURNRHCRIT,
            PDF_SHAPE=PDF_SHAPE,
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

        # Initalize the module
        with StencilBackendCompilerOverride(
            MPI.COMM_WORLD,
            self.stencil_config.dace_config,
        ):
            self.gfdl_1m = GFDL1M(
                self.stencil_factory, self.quantity_factory, self.GFDL_1M_config
            )

        # Link Fortran memory to Python memory #####
        # Fortran memory will only be modified if GFDL1M.__call__
        # is called from within a "with MAPLManagedMemory" statement #####
        # Not all linked fields are modified #####
        self._mapl_internal.register("Q", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QRAIN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QSNOW", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QGRAUPEL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QLCN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QICN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QLLS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("QILS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("CLCN", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("CLLS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("NACTL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_internal.register("NACTI", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("AREA", np.float32, [X_DIM, Y_DIM])
        self._mapl_import.register("PLE", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        self._mapl_import.register("ZLE", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM])
        self._mapl_import.register("T", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("U", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("V", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("FRLAND", np.float32, [X_DIM, Y_DIM])
        self._mapl_import.register("W", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("W2", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("W3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("WSL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("SL2", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("SL3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("WQT", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("QT2", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_import.register("QT3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("CNV_FRC", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("SRF_TYPE", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("SHLW_PRC3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("SHLW_SNO3", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("RHCRIT", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        self._mapl_export.register("RL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("RI", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("EVAPC", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("SUBLC", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("PRCP_RAIN", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("PRCP_SNOW", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("PRCP_ICE", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("PRCP_GRAUPEL", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("FCLD", np.float32, [X_DIM, Y_DIM, Z_DIM])

        self._mapl_export.register("QV", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("QL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("QI", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("QR", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("QS", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("QG", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("LTS", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("EIS", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("ZLCL", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register(
            "DUDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DVDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DTDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQVDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQLDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQIDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQADT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQRDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQSDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQGDT_macro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DUDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DVDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DTDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQVDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQLDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQIDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQADT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQRDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQSDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register(
            "DQGDT_micro", np.float32, [X_DIM, Y_DIM, Z_DIM], True
        )
        self._mapl_export.register("LS_PRCP", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("LS_SNR", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("ICE", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("FRZR", np.float32, [X_DIM, Y_DIM], True)
        self._mapl_export.register("RHX", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        self._mapl_export.register("REV_LS", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        self._mapl_export.register("RSU_LS", np.float32, [X_DIM, Y_DIM, Z_DIM], True)
        self._mapl_export.register(
            "PFL_LS", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True
        )
        self._mapl_export.register(
            "PFI_LS", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True
        )
        self._mapl_export.register(
            "PFL_AN", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True
        )
        self._mapl_export.register(
            "PFI_AN", np.float32, [X_DIM, Y_DIM, Z_INTERFACE_DIM], True
        )
        self._mapl_export.register("DQRL", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("DTDTFRIC", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("DBZ", np.float32, [X_DIM, Y_DIM, Z_DIM])
        self._mapl_export.register("DBZ_MAX", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("DBZ_1KM", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("DBZ_TOP", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("DBZ_M10C", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("CN_PRCP", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("AN_PRCP", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("SC_PRCP", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("CN_SNR", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("AN_SNR", np.float32, [X_DIM, Y_DIM])
        self._mapl_export.register("SC_SNR", np.float32, [X_DIM, Y_DIM])

    def GFDL_1M_Microphysics(self):
        from pyMoist.GFDL_1M.state import (
            CloudFractions,
            LiquidWaterStaticEnergy,
            MixingRatios,
            TotalWater,
            VericalMotion,
        )

        with (
            MAPLManagedMemory(self._mapl_internal) as mapl_internal,
            MAPLManagedMemory(self._mapl_import) as mapl_import,
            MAPLManagedMemory(self._mapl_export) as mapl_export,
        ):
            # Pull the data from the linked Fortran memory
            self.gfdl_1m.mixing_ratios = MixingRatios(
                vapor=mapl_internal.Q,
                rain=mapl_internal.QRAIN,
                snow=mapl_internal.QSNOW,
                graupel=mapl_internal.QGRAUPEL,
                convective_liquid=mapl_internal.QLCN,
                convective_ice=mapl_internal.QICN,
                large_scale_liquid=mapl_internal.QLLS,
                large_scale_ice=mapl_internal.QILS,
            )
            self.gfdl_1m.cloud_fractions = CloudFractions(
                convective=mapl_internal.CLCN,
                large_scale=mapl_internal.CLLS,
            )

            self.gfdl_1m.liquid_concentration = mapl_internal.NACTL
            self.gfdl_1m.ice_concentration = mapl_internal.NACTI
            self.gfdl_1m.area = mapl_import.AREA
            self.gfdl_1m.geopotential_height_interface = mapl_import.ZLE
            self.gfdl_1m.p_interface = mapl_import.PLE
            self.gfdl_1m.t = mapl_import.T
            self.gfdl_1m.u = mapl_import.U
            self.gfdl_1m.v = mapl_import.V
            self.gfdl_1m.land_fraction = mapl_import.FRLAND
            self.gfdl_1m.vertical_motion = VericalMotion(
                velocity=mapl_import.W,
                variance=mapl_import.W2,
                third_moment=mapl_import.W3,
            )
            self.gfdl_1m.liquid_water_static_energy = LiquidWaterStaticEnergy(
                flux=mapl_import.WSL,
                variance=mapl_import.SL2,
                third_moment=mapl_import.SL3,
            )
            self.gfdl_1m.total_water = TotalWater(
                flux=mapl_import.WQT,
                variance=mapl_import.QT2,
                third_moment=mapl_import.QT3,
            )

            self.gfdl_1m.convection_fraction = mapl_export.CNV_FRC
            self.gfdl_1m.surface_type = mapl_export.SRF_TYPE
            if mapl_export.associated("SHLW_PRC3"):
                self.gfdl_1m.shallow_convective_rain = mapl_export.SHLW_PRC3
            else:
                self.gfdl_1m.shallow_convective_rain = None
            if mapl_export.associated("SHLW_SNO3"):
                self.gfdl_1m.shallow_convective_snow = mapl_export.SHLW_SNO3
            else:
                self.gfdl_1m.shallow_convective_rain = None
            self.gfdl_1m.rh_crit = mapl_export.RHCRIT
            # Outputs: model fields originating from within GFDL
            self.gfdl_1m.outputs.liquid_radius = mapl_export.RL
            self.gfdl_1m.outputs.ice_radius = mapl_export.RI
            self.gfdl_1m.outputs.large_scale_nonanvil_precipitation_evaporation = (
                mapl_export.EVAPC
            )
            self.gfdl_1m.outputs.large_scale_nonanvil_precipitation_sublimation = (
                mapl_export.SUBLC
            )
            self.gfdl_1m.outputs.precipitated_rain = mapl_export.PRCP_RAIN
            self.gfdl_1m.outputs.precipitated_snow = mapl_export.PRCP_SNOW
            self.gfdl_1m.outputs.precipitated_ice = mapl_export.PRCP_ICE
            self.gfdl_1m.outputs.precipitated_graupel = mapl_export.PRCP_GRAUPEL

            # Outputs: model fields originating from within GFDL; radiation fields
            self.gfdl_1m.outputs.radiation_cloud_fraction = mapl_export.FCLD
            self.gfdl_1m.outputs.radiation_vapor = mapl_export.QV
            self.gfdl_1m.outputs.radiation_liquid = mapl_export.QL
            self.gfdl_1m.outputs.radiation_ice = mapl_export.QI
            self.gfdl_1m.outputs.radiation_rain = mapl_export.QR
            self.gfdl_1m.outputs.radiation_snow = mapl_export.QS
            self.gfdl_1m.outputs.radiation_graupel = mapl_export.QG
            self.gfdl_1m.outputs.lower_tropospheric_stability = mapl_export.LTS
            self.gfdl_1m.outputs.estimated_inversion_strength = mapl_export.EIS
            if mapl_export.associated("ZLCL"):
                self.gfdl_1m.outputs.z_lcl = mapl_export.ZLCL
            else:
                self.gfdl_1m.outputs.z_lcl = None

            # Outputs: model fields originating from within GFDL; macrophysics/microphysics tendencies
            self.gfdl_1m.outputs.du_dt_macro = mapl_export.DUDT_macro
            self.gfdl_1m.outputs.dv_dt_macro = mapl_export.DVDT_macro
            self.gfdl_1m.outputs.dt_dt_macro = mapl_export.DTDT_macro
            self.gfdl_1m.outputs.dvapor_dt_macro = mapl_export.DQVDT_macro
            self.gfdl_1m.outputs.dliquid_dt_macro = mapl_export.DQLDT_macro
            self.gfdl_1m.outputs.dice_dt_macro = mapl_export.DQIDT_macro
            self.gfdl_1m.outputs.dcloud_fraction_dt_macro = mapl_export.DQADT_macro
            self.gfdl_1m.outputs.drain_dt_macro = mapl_export.DQRDT_macro
            self.gfdl_1m.outputs.dsnow_dt_macro = mapl_export.DQSDT_macro
            self.gfdl_1m.outputs.dgraupel_dt_macro = mapl_export.DQGDT_macro
            self.gfdl_1m.outputs.du_dt_micro = mapl_export.DUDT_micro
            self.gfdl_1m.outputs.dv_dt_micro = mapl_export.DVDT_micro
            self.gfdl_1m.outputs.dt_dt_micro = mapl_export.DTDT_micro
            self.gfdl_1m.outputs.dvapor_dt_micro = mapl_export.DQVDT_micro
            self.gfdl_1m.outputs.dliquid_dt_micro = mapl_export.DQLDT_micro
            self.gfdl_1m.outputs.dice_dt_micro = mapl_export.DQIDT_micro
            self.gfdl_1m.outputs.dcloud_fraction_dt_micro = mapl_export.DQADT_micro
            self.gfdl_1m.outputs.drain_dt_micro = mapl_export.DQRDT_micro
            self.gfdl_1m.outputs.dsnow_dt_micro = mapl_export.DQSDT_micro
            self.gfdl_1m.outputs.dgraupel_dt_micro = mapl_export.DQGDT_micro
            # Outputs: Exports to be filled
            self.gfdl_1m.outputs.large_scale_precip = mapl_export.LS_PRCP
            self.gfdl_1m.outputs.large_scale_snow = mapl_export.LS_SNR
            self.gfdl_1m.outputs.icefall = mapl_export.ICE
            self.gfdl_1m.outputs.freezing_rainfall = mapl_export.FRZR
            self.gfdl_1m.outputs.relative_humidity_after_pdf = mapl_export.RHX
            self.gfdl_1m.outputs.large_scale_nonanvil_precipitation_evaporation = (
                mapl_export.REV_LS
            )
            self.gfdl_1m.outputs.large_scale_nonanvil_precipitation_sublimation = (
                mapl_export.RSU_LS
            )
            self.gfdl_1m.outputs.large_scale_nonanvil_liquid_flux = mapl_export.PFL_LS
            self.gfdl_1m.outputs.large_scale_nonanvil_ice_flux = mapl_export.PFI_LS
            self.gfdl_1m.outputs.anvil_liquid_flux = mapl_export.PFL_AN
            self.gfdl_1m.outputs.anvil_ice_flux = mapl_export.PFI_AN
            if mapl_export.associated("DQRL"):
                self.gfdl_1m.outputs.large_scale_rainwater_source = mapl_export.DQRL
            else:
                self.gfdl_1m.outputs.large_scale_rainwater_source = None
            if mapl_export.associated("DTDTFRIC"):
                self.gfdl_1m.outputs.moist_friction_temperature_tendency = (
                    mapl_export.DTDTFRIC
                )
            else:
                self.gfdl_1m.outputs.moist_friction_temperature_tendency = None
            if mapl_export.associated("DBZ"):
                self.gfdl_1m.outputs.simulated_reflectivity = mapl_export.DBZ
            else:
                self.gfdl_1m.outputs.simulated_reflectivity = None
            if mapl_export.associated("DBZ_MAX"):
                self.gfdl_1m.outputs.maximum_reflectivity = mapl_export.DBZ_MAX
            else:
                self.gfdl_1m.outputs.maximum_reflectivity = None
            if mapl_export.associated("DBZ_1KM"):
                self.gfdl_1m.outputs.one_km_agl_reflectivity = mapl_export.DBZ_1KM
            else:
                self.gfdl_1m.outputs.one_km_agl_reflectivity = None
            if mapl_export.associated("DBZ_TOP"):
                self.gfdl_1m.outputs.echo_top_reflectivity = mapl_export.DBZ_TOP
            else:
                self.gfdl_1m.outputs.echo_top_reflectivity = None
            if mapl_export.associated("DBZ_M10C"):
                self.gfdl_1m.outputs.minus_10c_reflectivity = mapl_export.DBZ_M10C
            else:
                self.gfdl_1m.outputs.minus_10c_reflectivity = None
            # Unused fields, force to zero
            self.gfdl_1m.outputs.deep_convective_precipitation = mapl_export.CN_PRCP
            self.gfdl_1m.outputs.anvil_precipitation = mapl_export.AN_PRCP
            self.gfdl_1m.outputs.shallow_convective_precipitation = mapl_export.SC_PRCP
            self.gfdl_1m.outputs.deep_convective_snow = mapl_export.CN_SNR
            self.gfdl_1m.outputs.anvil_snow = mapl_export.AN_SNR
            self.gfdl_1m.outputs.shallow_convective_snow = mapl_export.SC_SNR

            # Call the module
            with TimedCUDAProfiler("GFDL 1M", self._timings, True):
                self.gfdl_1m()
