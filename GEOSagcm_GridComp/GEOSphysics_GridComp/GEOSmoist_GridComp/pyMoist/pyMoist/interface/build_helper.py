import logging

from mpi4py import MPI

from ndsl import DaceConfig, DaCeOrchestration, ndsl_log
from ndsl.dsl.dace.build import set_distributed_caches


from gt4py.cartesian.config import build_settings as gt_build_settings  # isort: skip


class StencilBackendCompilerOverride:
    """Override the NDSL global stencil JIT to allow for 9-rank build
    on any setup.

    This is a workaround that requires to know _exactly_ when build is happening.
    Using this as a context manager, we leverage the DaCe build system to override
    the name and build the 9 codepaths required- while every other rank wait.

    This should be removed when we refactor the GT JIT to distribute building
    much more efficiently
    """

    def __init__(self, comm: MPI.Intracomm, config: DaceConfig | None):
        self.comm = comm

        if config is None:
            raise ValueError("DaCe Config can't be None")

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
