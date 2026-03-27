from MAPL_PythonBridge import UserCode, get_MAPLPy
from MAPL_PythonBridge.types import CVoidPointer
from mpi4py import MPI
from ndsl.dsl.typing import Float, Int
from ndsl.utils import safe_assign_array

from pyMoist.constants import NCNST
from pyMoist.fortran import get_NDSL_physics
from pyMoist.fortran.build_helper import StencilBackendCompilerOverride
from pyMoist.fortran.managed_state import MAPLManagedState
from pyMoist.fortran.memory_factory import MAPLMemoryRepository
from pyMoist.fortran.moist_workarounds import MOIST_WORKAROUNDS
from pyMoist.fortran.profiler import TimedCUDAProfiler
from pyMoist.convection.UW import ComputeUwshcuInv, UWConfiguration, UWState


class UWGEOSInterface(UserCode):
    def __init__(self, name: str) -> None:
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
            NCNST=NCNST,
            k0=ndsl_stack.quantity_factory.sizer.nz,
            dotransport=1 if MAPLPy.get_resource("USE_TRACER_TRANSP_UW:", mapl_state, default=True) else 0,
            dt=MAPLPy.get_resource("DSL__UW_DT:", mapl_state, default=Float(0)),
            windsrcavg=MAPLPy.get_resource("WINDSRCAVG:", mapl_state, default=Int(0) if jason_uw else Int(1)),
            mixscale=MAPLPy.get_resource(
                "MIXSCALE:", mapl_state, default=Float(0.0) if jason_uw else Float(3000.0)
            ),
            criqc=MAPLPy.get_resource(
                "CRIQC:", mapl_state, default=Float(1.0e-3) if jason_uw else Float(0.9e-3)
            ),
            thlsrc_fac=MAPLPy.get_resource(
                "THLSRC_FAC:", mapl_state, default=Float(0.0) if jason_uw else Float(1.0)
            ),
            frc_rasn=MAPLPy.get_resource(
                "FRC_RASN:", mapl_state, default=Float(0.0) if jason_uw else Float(1.0)
            ),
            rkm=MAPLPy.get_resource("RKM:", mapl_state, default=Float(12.0) if jason_uw else Float(8.0)),
            rpen=MAPLPy.get_resource("RPEN:", mapl_state, default=Float(3.0)),
            SCLM_SHALLOW=MAPLPy.get_resource("SCLM_SHALLOW:", mapl_state, default=Float(1.0)),
            niter_xc=MAPLPy.get_resource("NITER_XC:", mapl_state, default=Int(2)),
            iter_cin=MAPLPy.get_resource("ITER_CIN:", mapl_state, default=Int(2)),
            use_CINcin=1 if MAPLPy.get_resource("USE_CINCIN:", mapl_state, default=True) else 0,
            cridist_opt=MAPLPy.get_resource("CRIDIST_OPT:", mapl_state, default=Int(0)),
            use_self_detrain=1 if MAPLPy.get_resource("USE_SELF_DETRAIN:", mapl_state, default=False) else 0,
            use_momenflx=1 if MAPLPy.get_resource("USE_MOMENFLX:", mapl_state, default=True) else 0,
            use_cumpenent=1 if MAPLPy.get_resource("USE_CUMPENENT:", mapl_state, default=True) else 0,
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
        with StencilBackendCompilerOverride(
            MPI.COMM_WORLD,
            ndsl_stack.stencil_factory.config.dace_config,
        ):
            self._uw = ComputeUwshcuInv(ndsl_stack.stencil_factory, ndsl_stack.quantity_factory, config)

        # Init NDSL state
        self._managed_state = MAPLManagedState(
            UWState.empty(
                ndsl_stack.quantity_factory,
                data_dimensions=ndsl_stack.quantity_factory.sizer.data_dimensions,
            ),
            ndsl_stack.interface_type,
        )

    def run(self, mapl_state, import_state, export_state) -> None:
        raise RuntimeError("UW requires pyMoist integration requires `run_with_internal`")

    def run_with_internal(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
        internal_state: CVoidPointer,
    ) -> None:
        ndsl_stack = get_NDSL_physics(mapl_state)
        import_repository = MAPLMemoryRepository(import_state, ndsl_stack.quantity_factory)
        export_repository = MAPLMemoryRepository(export_state, ndsl_stack.quantity_factory)
        internal_repository = MAPLMemoryRepository(internal_state, ndsl_stack.quantity_factory)

        self._managed_state.register_K_interface("input.PLE", "PLE", import_repository)
        self._managed_state.register_K_interface("input.ZLE", "ZLE", import_repository)
        self._managed_state.register("input.QLLS", "QLLS", internal_repository)
        self._managed_state.register("input.QILS", "QILS", internal_repository)
        self._managed_state.register("input.QLCN", "QLCN", internal_repository)
        self._managed_state.register("input.QICN", "QICN", internal_repository)
        self._managed_state.register_2D("input.kpbl_inv", "KPBL_SC", import_repository)
        self._managed_state.register_2D("input.frland", "FRLAND", import_repository)
        self._managed_state.register_K_interface("input.tke_inv", "TKE", import_repository)
        self._managed_state.register_2D("input.shfx", "SH", import_repository)
        self._managed_state.register_2D("input.evap", "EVAP", import_repository)

        self._managed_state.register("input_output.u0_inv", "U", import_repository)
        self._managed_state.register("input_output.v0_inv", "V", import_repository)
        self._managed_state.register("input_output.qv0_inv", "Q", internal_repository)
        self._managed_state.register("input_output.t0_inv", "T", import_repository)
        self._managed_state.register_2D("input_output.cush", "CUSH", internal_repository)
        self._managed_state.register_2D("input_output.cnvtr", "CNPCPRATE", export_repository, alloc=True)

        self._managed_state.register_2D("output.RKFRE", "RKFRE", export_repository, alloc=True)
        self._managed_state.register("output.MFD_SC", "MFD_SC", export_repository, alloc=True)
        self._managed_state.register("output.QLENT_SC", "QLENT_SC", export_repository, alloc=True)
        self._managed_state.register("output.QIENT_SC", "QIENT_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface("output.umf_inv", "UMF_SC", export_repository, alloc=True)
        self._managed_state.register("output.dcm_inv", "DCM_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface(
            "output.qtflx_inv", "QTFLX_SC", export_repository, alloc=True
        )
        self._managed_state.register_K_interface(
            "output.slflx_inv", "SLFLX_SC", export_repository, alloc=True
        )
        self._managed_state.register_K_interface("output.uflx_inv", "UFLX_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface("output.vflx_inv", "VFLX_SC", export_repository, alloc=True)
        self._managed_state.register("output.DQADT_SC", "DQADT_SC", export_repository, alloc=True)
        self._managed_state.register("output.qvten_inv", "DQVDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.qlten_inv", "DQLDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.qiten_inv", "DQIDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.tten_inv", "DTDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.uten_inv", "DUDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.vten_inv", "DVDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.qrten_inv", "DQRDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.qsten_inv", "DQSDT_SC", export_repository, alloc=True)
        self._managed_state.register("output.cufrc_inv", "CUFRC_SC", export_repository, alloc=True)
        self._managed_state.register("output.fer_inv", "ENTR_SC", export_repository, alloc=True)
        self._managed_state.register("output.fdr_inv", "DETR_SC", export_repository, alloc=True)
        self._managed_state.register("output.ndrop_inv", "SC_NDROP", export_repository, alloc=True)
        self._managed_state.register("output.nice_inv", "SC_NICE", export_repository, alloc=True)
        self._managed_state.register("output.qlsub_inv", "QLSUB_SC", export_repository, alloc=True)
        self._managed_state.register("output.qisub_inv", "QISUB_SC", export_repository, alloc=True)
        self._managed_state.register("output.ql0_inv", "QLTOT", export_repository, alloc=True)
        self._managed_state.register("output.qi0_inv", "QITOT", export_repository, alloc=True)
        self._managed_state.register_2D("output.tpert_out", "TPERT_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.qpert_out", "QPERT_SC", export_repository, alloc=True)
        self._managed_state.register("output.qidet_inv", "QIDET_SC", export_repository, alloc=True)
        self._managed_state.register("output.qldet_inv", "QLDET_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface("output.CNV_MFC", "CNV_MFC", export_repository, alloc=True)
        self._managed_state.register("output.CNV_MFD", "CNV_MFD", export_repository, alloc=True)
        self._managed_state.register("output.SHLW_PRC3", "SHLW_PRC3", export_repository, alloc=True)
        self._managed_state.register("output.SHLW_SNO3", "SHLW_SNO3", export_repository, alloc=True)
        self._managed_state.register_2D("output.SC_QT", "SC_QT", export_repository, alloc=True)
        self._managed_state.register_2D("output.SC_MSE", "SC_MSE", export_repository, alloc=True)
        self._managed_state.register_2D("output.CUSH_SC", "CUSH_SC", export_repository, alloc=True)
        self._managed_state.register("input_output.CLCN", "CLCN", internal_repository)

        # Unused from GEOS ?!
        # CLLS = MAPLPy.get_pointer("CLLS", internal_state, dtype=np.float32)
        # CNV_FRC = MAPLPy.get_pointer("CNV_FRC", export_state, dtype=np.float32, alloc=True)
        # SRF_TYPE = MAPLPy.get_pointer("SRF_TYPE", export_state, dtype=np.float32, alloc=True)

        with TimedCUDAProfiler("UW", {}):
            with TimedCUDAProfiler("UW - State copy", {}):
                self._managed_state.fortran_to_ndsl()
                safe_assign_array(
                    self._managed_state.ndsl_state.input_output.CNV_Tracers.data[:],
                    MOIST_WORKAROUNDS.CNV_Tracers().Q[:],
                )

            with TimedCUDAProfiler("UW Numerics", {}):
                self._uw(self._managed_state.ndsl_state)

            with TimedCUDAProfiler("UW - State copy-back", {}):
                safe_assign_array(
                    MOIST_WORKAROUNDS.CNV_Tracers().Q[:],
                    self._managed_state.ndsl_state.input_output.CNV_Tracers.data[:],
                )
                self._managed_state.ndsl_to_fortran()

    def finalize(self, mapl_state, import_state, export_state) -> None:
        self._managed_state.save_recorded()


CODE = UWGEOSInterface("UW")
