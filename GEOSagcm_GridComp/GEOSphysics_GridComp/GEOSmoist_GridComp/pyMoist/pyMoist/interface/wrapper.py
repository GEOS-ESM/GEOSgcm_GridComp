"""
Wraps pyMoist for GEOS interface use.
"""

import dataclasses
import enum
import logging
import numpy as np
import os
from typing import Callable, Optional

from gt4py.cartesian.config import build_settings as gt_build_settings
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
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.dace.build import set_distributed_caches
from ndsl.dsl.gt4py_utils import is_gpu_backend
from ndsl.dsl.typing import get_precision
from ndsl.logging import ndsl_log
from ndsl.optional_imports import cupy as cp
from pyMoist.aer_activation import AerActivation
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.interface.flags import GFDL1MFlags, MoistFlags
from MAPLpyish import CVoidPointer
from pyMoist.interface.mapl.memory_factory import (
    MAPLMemoryRepository,
    MAPLManagedMemory,
)


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
    internal: CVoidPointer
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

        # JIT system for the component of Moist
        self._aer_activation: Optional[AerActivation] = None
        self._GFDL_1M_driver: Optional[MicrophysicsDriver] = None

        # Initalize flags later
        self.gfdl_microphysics_config = None

        self._mapl_internal = MAPLMemoryRepository(
            mapl_states.internal,
            self.quantity_factory,
        )
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
            if self.microphysics_config is None:
                raise RuntimeError("GFDL_1M flags not initalized")
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                self.stencil_config.dace_config,
            ):
                self._GFDL_1M_driver = MicrophysicsDriver(
                    self.stencil_factory,
                    self.quantity_factory,
                    self.microphysics_config,
                )
        return self._GFDL_1M_driver

    def init_gfdl_1m_configuration(
        self,
        flags: GFDL1MFlags,
    ):
        upper_case_dict = {}
        for field in dataclasses.fields(GFDL1MConfig):
            upper_case_dict[field.name] = getattr(flags, field.name.lower())
        self.microphysics_config = GFDL1MConfig(**upper_case_dict)

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

    def GFDL_Single_Moment_Microphysics(self):
        # call MAPL_GetResource( MAPL, MAX_RI , 'MAX_RI:' , DEFAULT=100.e-6, RC=STATUS); VERIFY_(STATUS)
        MAX_RI = self._mapl_comp.get_resource("MAX_RI:", np.float32, default=100.0e-6)
        print(f"MAX_RI From Python = {MAX_RI}, expected 9.99999975E-05")

        # call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
        # __init__
        self._mapl_import.register("T", np.float32, [X_DIM, Y_DIM, Z_DIM])

        # __call__
        with MAPLManagedMemory(self._mapl_import) as mmm:
            print(f"T From Python ({mmm.associated('T')}) = {mmm.T[7, 11, 29]}")
            mmm.T[11, 8, 29] = -300.00
