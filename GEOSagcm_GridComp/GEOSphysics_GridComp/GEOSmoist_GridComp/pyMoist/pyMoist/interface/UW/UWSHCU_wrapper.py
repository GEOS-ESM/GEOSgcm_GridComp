import enum
import logging
import os
from datetime import timedelta
from typing import Dict, List, Tuple

import f90nml
import numpy as np
from gt4py.cartesian.config import build_settings as gt_build_settings
from mpi4py import MPI

import pyFV3
from ndsl import (
    CompilationConfig,
    CubedSphereCommunicator,
    CubedSpherePartitioner,
    DaceConfig,
    DaCeOrchestration,
    GridIndexing,
    NullComm,
    PerformanceCollector,
    QuantityFactory,
    StencilConfig,
    StencilFactory,
    SubtileGridSizer,
    TilePartitioner,
    orchestrate,
)
import ndsl.constants
from ndsl.comm.comm_abc import Comm
from ndsl.dsl.dace.build import set_distributed_caches
from ndsl.dsl.gt4py_utils import is_gpu_backend
from ndsl.dsl.typing import floating_point_precision, Float, Int
from ndsl.grid import DampingCoefficients, GridData, MetricTerms
from ndsl.logging import ndsl_log
from ndsl.optional_imports import cupy as cp
from ndsl.utils import safe_assign_array
from ndsl.comm.mpi import MPIComm
from pyFV3.tracers import Tracers

from pyMoist.UW.compute_uwshcu import ComputeUwshcuInv
class UWSHCUwrapper:

    def __init__(
        self,
    ):
        print("UWSHCUwrapper.__init__")

        BACKEND = os.environ.get("GEOS_PYFV3_BACKEND", "gt:gpu")

        # Look for an override to run on a single node
        single_rank_override = int(os.getenv("GEOS_PYFV3_SINGLE_RANK_OVERRIDE", -1))
        if single_rank_override >= 0:
            comm = NullComm(single_rank_override, 6, 42)
        comm = MPIComm()

        # Make a custom performance collector for the UW wrapper
        self.perf_collector = PerformanceCollector("UWSHCU wrapper", comm)

        self.backend = BACKEND
        self.dycore_config = pyFV3._config.DynamicalCoreConfig()
        self.layout = self.dycore_config.layout
        
        partitioner = CubedSpherePartitioner(TilePartitioner(self.layout))
        self.communicator = CubedSphereCommunicator(
            comm,
            partitioner,
            timer=self.perf_collector.timestep_timer,
        )

        sizer = SubtileGridSizer.from_tile_params(
            nx_tile=self.dycore_config.npx - 1,  # NX/NY from config are cell-centers
            ny_tile=self.dycore_config.npy - 1,  # NX/NY from config are cell-centers
            nz=self.dycore_config.npz,
            n_halo=ndsl.constants.N_HALO_DEFAULT,
            extra_dim_lengths={},
            layout=self.dycore_config.layout,
            tile_partitioner=partitioner.tile,
            tile_rank=self.communicator.tile.rank,
        )

        quantity_factory = QuantityFactory.from_backend(sizer=sizer, backend=BACKEND)

        stencil_config = StencilConfig(
            compilation_config=CompilationConfig(
                backend=BACKEND, rebuild=False, validate_args=False
            ),
        )

        # Build a DaCeConfig for orchestration.
        # This and all orchestration code are transparent when outside
        # configuration deactivate orchestration
        stencil_config.dace_config = DaceConfig(
            communicator=self.communicator,
            backend=stencil_config.backend,
            tile_nx=self.dycore_config.npx,
            tile_nz=self.dycore_config.npz,
        )
        self._is_orchestrated = stencil_config.dace_config.is_dace_orchestrated()

        self._grid_indexing = GridIndexing.from_sizer_and_communicator(
            sizer=sizer, comm=self.communicator
        )

        stencil_factory = StencilFactory(
            config=stencil_config, grid_indexing=self._grid_indexing
        )

        self.compute_uwshcu = ComputeUwshcuInv(
            stencil_factory,
            quantity_factory,
        )

    def __call__(
        self,
        # input
        dotransport: Int,
        ncnst: Int,
        k0: Int,
        windsrcavg: Int,
        qtsrchgt: Float,
        qtsrc_fac: Float,
        thlsrc_fac: Float,
        frc_rasn: Float,
        rbuoy: Float,
        epsvarw: Float,
        use_CINcin: Int,
        mumin1: Float,
        rmaxfrac: Float,
        PGFc: Float,
        dt: Float,
        niter_xc: Int,
        criqc: Float,
        rle: Float,
        cridist_opt: Int,
        mixscale: Float,
        rkm: Float,
        detrhgt: Float,
        rdrag: Float,
        use_self_detrain: Int,
        use_cumpenent: Int,
        rpen: Float,
        use_momenflx: Int,
        rdrop: Float,
        pifc0_inv: np.ndarray,
        zifc0_inv: np.ndarray,
        pmid0_inv: np.ndarray,
        zmid0_inv: np.ndarray,
        kpbl_inv: np.ndarray,
        exnmid0_inv: np.ndarray,
        exnifc0_inv: np.ndarray,
        dp0_inv: np.ndarray,
        u0_inv: np.ndarray,
        v0_inv: np.ndarray,
        qv0_inv: np.ndarray,
        ql0_inv: np.ndarray,
        qi0_inv: np.ndarray,
        t0_inv: np.ndarray,
        frland_in: np.ndarray,
        tke_inv: np.ndarray,
        rkfre: np.ndarray,
        cush: np.ndarray,
        shfx: np.ndarray,
        evap: np.ndarray,
        cnvtr: np.ndarray,
        CNV_Tracers: np.ndarray,
        # input-output
        # output
        umf_inv: np.ndarray,
        dcm_inv: np.ndarray,
        qtflx_inv: np.ndarray,
        slflx_inv: np.ndarray,
        uflx_inv: np.ndarray,
        vflx_inv: np.ndarray,
        qvten_inv: np.ndarray,
        qlten_inv: np.ndarray,
        qiten_inv: np.ndarray,
        tten_inv: np.ndarray,
        uten_inv: np.ndarray,
        vten_inv: np.ndarray,
        qrten_inv: np.ndarray,
        qsten_inv: np.ndarray,
        cufrc_inv: np.ndarray,
        fer_inv: np.ndarray,
        fdr_inv: np.ndarray,
        ndrop_inv: np.ndarray,
        nice_inv: np.ndarray,
        qldet_inv: np.ndarray,
        qlsub_inv: np.ndarray,
        qidet_inv: np.ndarray,
        qisub_inv: np.ndarray,
        tpert_out: np.ndarray,
        qpert_out: np.ndarray,
        ):
        print("UWSHCUwrapper.__call__")
        # This is where the actual work would be done
        self.compute_uwshcu(
            # Field inputs
            pifc0_inv=pifc0_inv,
            zifc0_inv=zifc0_inv,
            pmid0_inv=pmid0_inv,
            zmid0_inv=zmid0_inv,
            kpbl_inv=kpbl_inv,
            exnmid0_inv=exnmid0_inv,
            exnifc0_inv=exnifc0_inv,
            dp0_inv=dp0_inv,
            u0_inv=u0_inv,
            v0_inv=v0_inv,
            qv0_inv=qv0_inv,
            ql0_inv=ql0_inv,
            qi0_inv=qi0_inv,
            t0_inv=t0_inv,
            frland=frland_in,
            tke_inv=tke_inv,
            rkfre=rkfre,
            cush=cush,
            shfx=shfx,
            evap=evap,
            cnvtr=cnvtr,
            CNV_Tracers=CNV_Tracers,
            # Float/Int inputs
            dotransport=dotransport,
            ncnst=ncnst,
            k0=k0,
            windsrcavg=windsrcavg,
            qtsrchgt=qtsrchgt,
            qtsrc_fac=qtsrc_fac,
            thlsrc_fac=thlsrc_fac,
            frc_rasn=frc_rasn,
            rbuoy=rbuoy,
            epsvarw=epsvarw,
            use_CINcin=use_CINcin,
            mumin1=mumin1,
            rmaxfrac=rmaxfrac,
            PGFc=PGFc,
            dt=dt,
            niter_xc=niter_xc,
            criqc=criqc,
            rle=rle,
            cridist_opt=cridist_opt,
            mixscale=mixscale,
            rdrag=rdrag,
            rkm=rkm,
            use_self_detrain=use_self_detrain,
            detrhgt=detrhgt,
            use_cumpenent=use_cumpenent,
            rpen=rpen,
            use_momenflx=use_momenflx,
            rdrop=rdrop,
            # Outputs
            umf_inv=umf_inv,
            dcm_inv=dcm_inv,
            qtflx_inv=qtflx_inv,
            slflx_inv=slflx_inv,
            uflx_inv=uflx_inv,
            vflx_inv=vflx_inv,
            qvten_inv=qvten_inv,
            qlten_inv=qlten_inv,
            qiten_inv=qiten_inv,
            tten_inv=tten_inv,
            uten_inv=uten_inv,
            vten_inv=vten_inv,
            qrten_inv=qrten_inv,
            qsten_inv=qsten_inv,
            cufrc_inv=cufrc_inv,
            fer_inv=fer_inv,
            fdr_inv=fdr_inv,
            ndrop_inv=ndrop_inv,
            nice_inv=nice_inv,
            qldet_inv=qldet_inv,
            qlsub_inv=qlsub_inv,
            qidet_inv=qidet_inv,
            qisub_inv=qisub_inv,
            tpert_out=tpert_out,
            qpert_out=qpert_out,
        )