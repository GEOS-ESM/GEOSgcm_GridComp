"""
Wraps pyMoist for GEOS interface use.
"""

import dataclasses
import os

from MAPL_PythonBridge import get_MAPLPy
from ndsl import (
    Backend,
    CompilationConfig,
    CubedSphereCommunicator,
    CubedSpherePartitioner,
    DaceConfig,
    GridIndexing,
    LocalComm,
    MPIComm,
    PerformanceCollector,
    QuantityFactory,
    StencilConfig,
    StencilFactory,
    SubtileGridSizer,
    TileCommunicator,
    TilePartitioner,
)
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.typing import get_precision
from ndsl.logging import ndsl_log_on_rank_0
from ndsl.optional_imports import cupy as cp

from pyMoist.fortran.build_helper import InterfaceTransferType, MemorySpace


@dataclasses.dataclass
class NDSLPhysicsConfiguration:
    # Grid layout
    npx: int = 0
    npy: int = 0
    npz: int = 0
    layout_x: int = 1
    layout_y: int = 1
    single_column: bool = False
    backend: str = "st:dace:cpu:KJI"


class NDSLPhysicsStack:
    def __init__(
        self,
        flags: NDSLPhysicsConfiguration,
        fortran_mem_space: MemorySpace = MemorySpace.CPU,
    ) -> None:
        # Look for an override to run on a single node
        single_rank_override = int(os.getenv("GEOS_PYFV3_SINGLE_RANK_OVERRIDE", -1))
        comm = MPIComm()
        if single_rank_override >= 0:
            comm = LocalComm(rank=single_rank_override, total_ranks=6, buffer_dict={})

        self.backend = Backend(flags.backend)
        self.flags = flags
        layout = (self.flags.layout_x, self.flags.layout_y)

        # Make a custom performance collector for the GEOS wrapper
        self.perf_collector = PerformanceCollector("GEOS Moist", comm)

        if flags.single_column:
            partitioner = TilePartitioner(layout)
            self.communicator = TileCommunicator(
                comm,
                partitioner,
                timer=self.perf_collector.timestep_timer,
            )
        else:
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
            backend=self.backend,
        )
        self.quantity_factory = QuantityFactory(sizer=sizer, backend=self.backend)

        self.stencil_config = StencilConfig(
            compilation_config=CompilationConfig(backend=self.backend, rebuild=False, validate_args=True),
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

        self._grid_indexing = GridIndexing.from_sizer_and_communicator(sizer=sizer, comm=self.communicator)
        self.stencil_factory = StencilFactory(config=self.stencil_config, grid_indexing=self._grid_indexing)

        # Figure out the interface mode
        tmp_quantity = self.quantity_factory.empty([I_DIM, J_DIM, K_DIM], units="")
        default_3D_memory_desc = (tmp_quantity.data.shape, tmp_quantity.data.strides)
        if fortran_mem_space != MemorySpace.CPU:
            raise NotImplementedError("Interface cannot stream Fortran memory resident on GPU")
        if self.backend.is_gpu_backend():
            self._interface_type = InterfaceTransferType.CPU_TO_GPU_TO_CPU
        else:
            if self.backend.is_fortran_aligned():
                # This is Fortran layout - we can Map the memory
                self._interface_type = InterfaceTransferType.CPU_MAP
            else:
                # All other layout have to copy the data in/out of Fortran layout
                self._interface_type = InterfaceTransferType.CPU_COPY
        del tmp_quantity

        # Feedback information
        device_ordinal_info = "N/A"
        if cp is not None:
            device_ordinal_info = f"  Device PCI bus id: {cp.cuda.Device(0).pci_bus_id}" if self.backend.is_gpu_backend() else "N/A"
        MPS_pipe_directory = os.getenv("CUDA_MPS_PIPE_DIRECTORY", None)
        MPS_is_on = MPS_pipe_directory is not None and self.backend.is_gpu_backend() and os.path.exists(f"{MPS_pipe_directory}/log")
        ndsl_log_on_rank_0.info(
            "pyMoist <> GEOS wrapper initialized (Rank 0):\n"
            f"         Bridge : {self._interface_type.name}\n"
            f"        Backend : {self.backend}\n"
            f"      Precision : {get_precision()} bit\n"
            f"  Orchestration : {self._is_orchestrated}\n"
            f"          Sizer : {sizer.nx}x{sizer.ny}x{sizer.nz}"
            f"(halo: {sizer.n_halo})\n"
            f" Strides for 3D : {default_3D_memory_desc[1]}\n"
            f"     Device ord : {device_ordinal_info}\n"
            f"     Nvidia MPS : {MPS_is_on}\n"
        )

    @property
    def interface_type(self) -> InterfaceTransferType:
        return self._interface_type


NDSL_PHYSICS: NDSLPhysicsStack | None = None


def _set_NDSL_physics(mapl_state) -> NDSLPhysicsStack:
    """Initialization of the global NDSL stack for Physics parametrization"""
    MAPLPy = get_MAPLPy()
    grid_infos = MAPLPy.get_grid_infos(mapl_state)
    single_column = False
    ny = grid_infos.ny // 6
    if grid_infos.ny < 6:
        single_column = True
        ny = 1

    override_backend = os.getenv("GEOS_DSL_BACKEND", "st:dace:cpu:KJI")

    return NDSLPhysicsStack(
        NDSLPhysicsConfiguration(
            grid_infos.im * grid_infos.nx,
            grid_infos.jm * ny,
            grid_infos.lm,
            grid_infos.nx,
            ny,
            single_column,
            override_backend,
        )
    )


def get_NDSL_physics(mapl_state) -> NDSLPhysicsStack:
    """Retrieve the global MAPLPy accessor"""
    global NDSL_PHYSICS
    if NDSL_PHYSICS is None:
        NDSL_PHYSICS = _set_NDSL_physics(mapl_state)
    return NDSL_PHYSICS
