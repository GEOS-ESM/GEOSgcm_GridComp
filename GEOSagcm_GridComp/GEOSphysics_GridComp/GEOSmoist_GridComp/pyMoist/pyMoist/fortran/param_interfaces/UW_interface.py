import numpy as np
from ndsl.dsl.typing import Float, Int
from pyGEOSBridge import GEOSInterfaceCode, get_MAPLPy

from pyMoist.fortran import get_NDSL_physics
from pyMoist.UW import ComputeUwshcuInv, UWConfiguration, UWState


class UWGEOSInterface(GEOSInterfaceCode):
    def __init__(self, name: str) -> None:
        self._state = None

    def init(self, mapl_state, import_state, export_state) -> None:
        # Make sure we have our NDSL stack setup
        MAPLPy = get_MAPLPy()
        ndsl_stack = get_NDSL_physics(mapl_state)

        # Init NDSL state
        self._state = UWState.empty(ndsl_stack.quantity_factory)

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
        )

        # Build UW
        self._uw = ComputeUwshcuInv(ndsl_stack.stencil_factory, ndsl_stack.quantity_factory, config)

    def run(self, mapl_state, import_state, export_state) -> None:
        raise RuntimeError("UW requires pyMoist integration requires `run_with_internal`")

    def run_with_internal(self, mapl_state, import_state, export_state, internal_state) -> None:
        MAPLPy = get_MAPLPy()

        # Internals
        Q = MAPLPy.get_pointer("Q", internal_state, dtype=np.float32)  # ??
        self._state.input_output.qv0_inv.data[:] = Q
        QLLS = MAPLPy.get_pointer("QLLS", internal_state, dtype=np.float32)
        self._state.input.QLLS.data[:] = QLLS
        QLCN = MAPLPy.get_pointer("QLCN", internal_state, dtype=np.float32)
        self._state.input.QLCN.data[:] = QLCN
        CLCN = MAPLPy.get_pointer("CLCN", internal_state, dtype=np.float32)  # ??
        CLLS = MAPLPy.get_pointer("CLLS", internal_state, dtype=np.float32)  # ❌
        QILS = MAPLPy.get_pointer("QILS", internal_state, dtype=np.float32)
        self._state.input.QILS.data[:] = QILS
        QICN = MAPLPy.get_pointer("QICN", internal_state, dtype=np.float32)
        self._state.input.QICN.data[:] = QICN
        CUSH = MAPLPy.get_pointer("CUSH", internal_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        self._state.input_output.cush.data[:] = CUSH

        # Imports
        FRLAND = MAPLPy.get_pointer("FRLAND", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        self._state.input.frland.data[:] = FRLAND
        ZLE = MAPLPy.get_pointer("ZLE", import_state, dtype=np.float32)
        self._state.input.ZLE.data[:] = ZLE
        PLE = MAPLPy.get_pointer("PLE", import_state, dtype=np.float32)
        self._state.input.PLE.data[:] = PLE
        T = MAPLPy.get_pointer("T", import_state, dtype=np.float32)
        self._state.input_output.t0_inv.data[:] = T
        U = MAPLPy.get_pointer("U", import_state, dtype=np.float32)
        self._state.input_output.u0_inv.data[:] = U
        V = MAPLPy.get_pointer("V", import_state, dtype=np.float32)
        self._state.input_output.v0_inv.data[:] = V
        SH = MAPLPy.get_pointer("SH", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        self._state.input.shfx.data[:] = SH
        EVAP = MAPLPy.get_pointer("EVAP", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        self._state.input.evap.data[:] = EVAP
        KPBL_SC = MAPLPy.get_pointer("KPBL_SC", import_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2])
        self._state.input.kpbl_inv.data[:] = KPBL_SC
        TKE = MAPLPy.get_pointer("TKE", import_state, dtype=np.float32)
        self._state.input.tke_inv.data[:] = TKE

        # Required Exports (connectivities to moist siblings)
        MFD_SC = MAPLPy.get_pointer("MFD_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.MFD_SC.data[:] = MFD_SC
        QLDET_SC = MAPLPy.get_pointer("QLDET_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qldet_inv.data[:] = QLDET_SC
        QIDET_SC = MAPLPy.get_pointer("QIDET_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qidet_inv.data[:] = QIDET_SC
        CUFRC_SC = MAPLPy.get_pointer("CUFRC_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.cufrc_inv.data[:] = CUFRC_SC
        CNPCPRATE = MAPLPy.get_pointer("CNPCPRATE", export_state, dtype=np.float32, alloc=True)
        self._state.input_output.cnvtr.data[:] = CNPCPRATE
        CNV_FRC = MAPLPy.get_pointer("CNV_FRC", export_state, dtype=np.float32, alloc=True)
        SRF_TYPE = MAPLPy.get_pointer("SRF_TYPE", export_state, dtype=np.float32, alloc=True)
        # Tendency Export
        DUDT_SC = MAPLPy.get_pointer("DUDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.uten_inv.data[:] = DUDT_SC
        DVDT_SC = MAPLPy.get_pointer("DVDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.vten_inv.data[:] = DVDT_SC
        DTDT_SC = MAPLPy.get_pointer("DTDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.tten_inv.data[:] = DTDT_SC
        DQVDT_SC = MAPLPy.get_pointer("DQVDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qvten_inv.data[:] = DQVDT_SC
        DQIDT_SC = MAPLPy.get_pointer("DQIDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qiten_inv.data[:] = DQIDT_SC
        DQLDT_SC = MAPLPy.get_pointer("DQLDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qlten_inv.data[:] = DQLDT_SC
        DQRDT_SC = MAPLPy.get_pointer("DQRDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qrten_inv.data[:] = DQRDT_SC
        DQSDT_SC = MAPLPy.get_pointer("DQSDT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qsten_inv.data[:] = DQSDT_SC
        DQADT_SC = MAPLPy.get_pointer("DQADT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.DQADT_SC.data[:] = DQADT_SC
        # Exports
        UMF_SC = MAPLPy.get_pointer("UMF_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.umf_inv.data[:] = UMF_SC
        DCM_SC = MAPLPy.get_pointer("DCM_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.dcm_inv.data[:] = DCM_SC
        ENTR_SC = MAPLPy.get_pointer("ENTR_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.fer_inv.data[:] = ENTR_SC
        DETR_SC = MAPLPy.get_pointer("DETR_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.fdr_inv.data[:] = DETR_SC
        QLSUB_SC = MAPLPy.get_pointer("QLSUB_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qlsub_inv.data[:] = QLSUB_SC
        QISUB_SC = MAPLPy.get_pointer("QISUB_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qisub_inv.data[:] = QISUB_SC
        SC_NDROP = MAPLPy.get_pointer("SC_NDROP", export_state, dtype=np.float32, alloc=True)
        self._state.output.ndrop_inv.data[:] = SC_NDROP
        SC_NICE = MAPLPy.get_pointer("SC_NICE", export_state, dtype=np.float32, alloc=True)
        self._state.output.nice_inv.data[:] = SC_NICE
        TPERT_SC = MAPLPy.get_pointer("TPERT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.tpert_out.data[:] = TPERT_SC
        QPERT_SC = MAPLPy.get_pointer("QPERT_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qpert_out.data[:] = QPERT_SC
        QTFLX_SC = MAPLPy.get_pointer("QTFLX_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.qtflx_inv.data[:] = QTFLX_SC
        SLFLX_SC = MAPLPy.get_pointer("SLFLX_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.slflx_inv.data[:] = SLFLX_SC
        UFLX_SC = MAPLPy.get_pointer("UFLX_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.uflx_inv.data[:] = UFLX_SC
        VFLX_SC = MAPLPy.get_pointer("VFLX_SC", export_state, dtype=np.float32, alloc=True)
        self._state.output.vflx_inv.data[:] = VFLX_SC
        RKFRE = MAPLPy.get_pointer(
            "RKFRE", export_state, dtype=np.float32, dims=MAPLPy.grid_dims[0:2], alloc=True
        )
        self._state.output.RKFRE.data[:] = RKFRE

        self._uw(self._state)

    def finalize(self, mapl_state, import_state, export_state) -> None:
        # No finalize call from UW
        pass


CODE = UWGEOSInterface("UW")
