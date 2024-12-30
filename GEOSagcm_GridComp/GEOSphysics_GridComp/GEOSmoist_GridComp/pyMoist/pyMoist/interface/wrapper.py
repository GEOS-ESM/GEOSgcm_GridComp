"""
Wraps pyMoist for GEOS interface use.
"""

import enum
import logging
import os
from typing import Callable

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
from ndsl.dsl.typing import floating_point_precision
from ndsl.logging import ndsl_log
from ndsl.optional_imports import cupy as cp
from pyMoist.aer_activation import AerActivation
from pyMoist.GFDL_1M.GFDL_1M import GFDL_1M
from pyMoist.interface.flags import MoistFlags


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

    def __exit__(self, type, value, traceback):
        if self.no_op:
            return
        if not self.config.do_compile:
            ndsl_log.info(f"Stencil backend read cache on {self.comm.Get_rank()}")
        else:
            ndsl_log.info(f"Stencil backend compiled on {self.comm.Get_rank()}")
            self.comm.Barrier()


class GEOSPyMoistWrapper:
    def __init__(
        self,
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
            f"          float : {floating_point_precision()}bit\n"
            f"  orchestration : {self._is_orchestrated}\n"
            f"          sizer : {sizer.nx}x{sizer.ny}x{sizer.nz}"
            f"(halo: {sizer.n_halo})\n"
            f"     Device ord : {device_ordinal_info}\n"
            f"     Nvidia MPS : {MPS_is_on}"
        )

        # JIT system for the component of Moist
        self._aer_activation = None
        self._GFDL_1M_evap = None
        self._GFDL_1M_driver = None

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

    @property
    def GFDL_1M_evap(self) -> Callable:
        if not self._GFDL_1M_evap:
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                self.stencil_config.dace_config,
            ):
                self._GFDL_1M_evap = GFDL_1M(
                    stencil_factory=self.stencil_factory,
                    quantity_factory=self.quantity_factory,
                )
        return self._GFDL_1M_evap

    @property
    def GFDL_1M_driver(self) -> Callable:
        if not self._GFDL_1M_driver:
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                self.stencil_config.dace_config,
            ):
                print("_GFDL_1M_driver: BUILD ME PLEASE")
        return self._GFDL_1M_driver

    def make_nmmodes_quantity(self, data):
        qty = self.nmodes_quantity_factory.empty(
            [X_DIM, Y_DIM, Z_DIM, "n_modes"],
            "n/a",
        )
        qty.view[:, :, :, :] = qty.np.asarray(data[:, :, :, :])
        return qty
