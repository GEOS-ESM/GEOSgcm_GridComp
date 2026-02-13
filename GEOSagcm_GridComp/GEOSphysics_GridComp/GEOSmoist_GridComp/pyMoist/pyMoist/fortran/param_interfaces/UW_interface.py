import numpy as np
from ndsl.dsl.typing import Float, Int
from pyGEOSBridge import GEOSInterfaceCode, get_MAPLPy

from pyMoist.fortran import get_NDSL_physics
from pyMoist.UW import ComputeUwshcuInv, UWConfiguration


class UWGEOSInterface(GEOSInterfaceCode):
    def __init__(self, name: str) -> None:
        # Real initialize happening in `init`
        pass

    def init(self, mapl_state, import_state, export_state) -> None:
        # Make sure we have our NDSL stack setup
        MAPLPy = get_MAPLPy()
        ndsl_stack = get_NDSL_physics(mapl_state)

        if ndsl_stack.quantity_factory.sizer.nz == 72:
            jason_uw = MAPLPy.get_resource("JASON_UW:", mapl_state, default=True)
        else:
            jason_uw = MAPLPy.get_resource("JASON_UW:", mapl_state, default=False)

        # Compile the configuration for UW
        config = UWConfiguration(
            NCNST=23,
            k0=ndsl_stack.quantity_factory.sizer.nz,
            dotransport=MAPLPy.get_resource("USE_TRACER_TRANSP_UW:", mapl_state, default=True),
            dt=MAPLPy.get_resource("DSL__UW_DT:", mapl_state, default=Float(0)),
            windsrcavg=MAPLPy.get_resource("WINDSRCAVG:", mapl_state, default=Int(0)),
            mixscale=MAPLPy.get_resource("MIXSCALE:", mapl_state, default=Float(0.0)),
            criqc=MAPLPy.get_resource(
                "CRIQC:", mapl_state, default=Float(1.0e-3) if jason_uw else Float(0.9e-3)
            ),
            thlsrc_fac=MAPLPy.get_resource(
                "THLSRC_FAC:", mapl_state, default=Float(0.0) if jason_uw else Float(1.0)
            ),
            frc_rasn=MAPLPy.get_resource(
                "THLSRC_FAC:", mapl_state, default=Float(0.0) if jason_uw else Float(1.0)
            ),
            rkm=MAPLPy.get_resource("RKM:", mapl_state, default=Float(12.0) if jason_uw else Float(8.0)),
            rpen=MAPLPy.get_resource("RPEN:", mapl_state, default=Float(3.0)),
            SCLM_SHALLOW=MAPLPy.get_resource("SCLM_SHALLOW:", mapl_state, default=Float(1.0)),
            niter_xc=MAPLPy.get_resource("NITER_XC:", mapl_state, default=Float(0.7)),
            iter_cin=MAPLPy.get_resource("ITER_CIN:", mapl_state, default=Int(2)),
            use_CINcin=MAPLPy.get_resource("USE_CINCIN:", mapl_state, default=True),
            cridist_opt=MAPLPy.get_resource("CRIDIST_OPT:", mapl_state, default=Int(0)),
            use_self_detrain=MAPLPy.get_resource("USE_SELF_DETRAIN:", mapl_state, default=False),
            use_momenflx=MAPLPy.get_resource("USE_MOMENFLX:", mapl_state, default=True),
            use_cumpenent=MAPLPy.get_resource("USE_CUMPENENT:", mapl_state, default=True),
            rle=MAPLPy.get_resource("RLE:", mapl_state, default=Float(0.1)),
            rmaxfrac=MAPLPy.get_resource("RMAXFRAC:", mapl_state, default=Float(0.1)),
            mumin1=MAPLPy.get_resource("MUMIN1:", mapl_state, default=Float(0.906)),
            rbuoy=MAPLPy.get_resource("RBUOY:", mapl_state, default=Float(1.0)),
            rdrag=MAPLPy.get_resource("RDRAG:", mapl_state, default=Float(1.0)),
            epsvarw=MAPLPy.get_resource("EPSVARW:", mapl_state, default=Float(5.0e-4)),
            PGFc=MAPLPy.get_resource("PGFC:", mapl_state, default=Float(0.7)),
            rdrop=MAPLPy.get_resource("SHLW_RDROP:", mapl_state, default=Float(8.0e-6)),
            detrhgt=MAPLPy.get_resource("DETRHGT:", mapl_state, default=Float(1800.0)),
            qtsrc_fac=MAPLPy.get_resource("QTSRC_FAC:", mapl_state, default=Float(0.0)),
            qtsrchgt=MAPLPy.get_resource("QTSRCHGT:", mapl_state, default=Float(40.0)),
            status_PTR2D_1=False,
            status_PTR2D_2=False,
            status_PTR2D_3=False,
            status_PTR3D_1=False,
            status_PTR3D_2=False,
        )

        # Build UW
        self._uw = ComputeUwshcuInv(ndsl_stack.stencil_factory, ndsl_stack.quantity_factory, config)

    def run(self, mapl_state, import_state, export_state) -> None:
        raise RuntimeError("UW requires pyMoist integration requires `run_with_internal`")

    def run_with_internal(self, mapl_state, import_state, export_state, internal_state) -> None:
        MAPLPy = get_MAPLPy()

        # Internals
        Q = MAPLPy.get_pointer("Q", internal_state, dtype=np.float32)  # ??
        QLLS = MAPLPy.get_pointer("QLLS", internal_state, dtype=np.float32)
        QLCN = MAPLPy.get_pointer("QLCN", internal_state, dtype=np.float32)
        CLCN = MAPLPy.get_pointer("CLCN", internal_state, dtype=np.float32)  # ??
        CLLS = MAPLPy.get_pointer("CLLS", internal_state, dtype=np.float32)  # ❌
        QILS = MAPLPy.get_pointer("QILS", internal_state, dtype=np.float32)
        QICN = MAPLPy.get_pointer("QICN", internal_state, dtype=np.float32)
        CUSH = MAPLPy.get_pointer("CUSH", internal_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])

        # Imports
        FRLAND = MAPLPy.get_pointer("FRLAND", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        ZLE = MAPLPy.get_pointer("ZLE", import_state, dtype=np.float32)
        PLE = MAPLPy.get_pointer("PLE", import_state, dtype=np.float32)
        T = MAPLPy.get_pointer("T", import_state, dtype=np.float32)
        U = MAPLPy.get_pointer("U", import_state, dtype=np.float32)
        V = MAPLPy.get_pointer("V", import_state, dtype=np.float32)
        SH = MAPLPy.get_pointer("SH", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        EVAP = MAPLPy.get_pointer("EVAP", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        KPBL_SC = MAPLPy.get_pointer("KPBL_SC", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        TKE = MAPLPy.get_pointer("TKE", import_state, dtype=np.float32)

        # Required Exports (connectivities to moist siblings)
        MFD_SC = MAPLPy.get_pointer("MFD_SC", export_state, dtype=np.float32, alloc=True)
        QLDET_SC = MAPLPy.get_pointer("QLDET_SC", export_state, dtype=np.float32, alloc=True)
        QIDET_SC = MAPLPy.get_pointer("QIDET_SC", export_state, dtype=np.float32, alloc=True)
        CUFRC_SC = MAPLPy.get_pointer("CUFRC_SC", export_state, dtype=np.float32, alloc=True)
        CNPCPRATE = MAPLPy.get_pointer("CNPCPRATE", export_state, dtype=np.float32, alloc=True)
        CNV_FRC = MAPLPy.get_pointer("CNV_FRC", export_state, dtype=np.float32, alloc=True)
        SRF_TYPE = MAPLPy.get_pointer("SRF_TYPE", export_state, dtype=np.float32, alloc=True)
        # Tendency Export
        DUDT_SC = MAPLPy.get_pointer("DUDT_SC", export_state, dtype=np.float32, alloc=True)
        DVDT_SC = MAPLPy.get_pointer("DVDT_SC", export_state, dtype=np.float32, alloc=True)
        DTDT_SC = MAPLPy.get_pointer("DTDT_SC", export_state, dtype=np.float32, alloc=True)
        DQVDT_SC = MAPLPy.get_pointer("DQVDT_SC", export_state, dtype=np.float32, alloc=True)
        DQIDT_SC = MAPLPy.get_pointer("DQIDT_SC", export_state, dtype=np.float32, alloc=True)
        DQLDT_SC = MAPLPy.get_pointer("DQLDT_SC", export_state, dtype=np.float32, alloc=True)
        DQRDT_SC = MAPLPy.get_pointer("DQRDT_SC", export_state, dtype=np.float32, alloc=True)
        DQSDT_SC = MAPLPy.get_pointer("DQSDT_SC", export_state, dtype=np.float32, alloc=True)
        DQADT_SC = MAPLPy.get_pointer("DQADT_SC", export_state, dtype=np.float32, alloc=True)
        # Exports
        UMF_SC = MAPLPy.get_pointer("UMF_SC", export_state, dtype=np.float32, alloc=True)
        DCM_SC = MAPLPy.get_pointer("DCM_SC", export_state, dtype=np.float32, alloc=True)
        ENTR_SC = MAPLPy.get_pointer("ENTR_SC", export_state, dtype=np.float32, alloc=True)
        DETR_SC = MAPLPy.get_pointer("DETR_SC", export_state, dtype=np.float32, alloc=True)
        QLSUB_SC = MAPLPy.get_pointer("QLSUB_SC", export_state, dtype=np.float32, alloc=True)
        QISUB_SC = MAPLPy.get_pointer("QISUB_SC", export_state, dtype=np.float32, alloc=True)
        SC_NDROP = MAPLPy.get_pointer("SC_NDROP", export_state, dtype=np.float32, alloc=True)
        SC_NICE = MAPLPy.get_pointer("SC_NICE", export_state, dtype=np.float32, alloc=True)
        TPERT_SC = MAPLPy.get_pointer("TPERT_SC", export_state, dtype=np.float32, alloc=True)
        QPERT_SC = MAPLPy.get_pointer("QPERT_SC", export_state, dtype=np.float32, alloc=True)
        QTFLX_SC = MAPLPy.get_pointer("QTFLX_SC", export_state, dtype=np.float32, alloc=True)
        SLFLX_SC = MAPLPy.get_pointer("SLFLX_SC", export_state, dtype=np.float32, alloc=True)
        UFLX_SC = MAPLPy.get_pointer("UFLX_SC", export_state, dtype=np.float32, alloc=True)
        VFLX_SC = MAPLPy.get_pointer("VFLX_SC", export_state, dtype=np.float32, alloc=True)
        RKFRE = MAPLPy.get_pointer(
            "RKFRE", export_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2], alloc=True
        )

        # self._uw(...)

    def finalize(self, mapl_state, import_state, export_state) -> None:
        # No finalize call from UW
        pass


CODE = UWGEOSInterface("UW")
