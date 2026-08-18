from MAPL_PythonBridge import UserCode, get_MAPLPy
from MAPL_PythonBridge.types import CVoidPointer
from mpi4py import MPI
from ndsl.dsl.typing import Float, Int
from ndsl.utils import safe_assign_array

from pyMoist.convection.UW import ComputeUwshcuInv, UWConfiguration, UWState
from pyMoist.fortran import get_NDSL_physics
from pyMoist.fortran.build_helper import StencilBackendCompilerOverride
from pyMoist.fortran.cuda_profiler import TimedCUDAProfiler
from pyMoist.fortran.managed_state import MAPLManagedState
from pyMoist.fortran.memory_factory import MAPLMemoryRepository
from pyMoist.fortran.moist_workarounds import MOIST_WORKAROUNDS
from pyMoist.convection_tracers import CONVECTION_TRACER_DIM, FloatField_ConvectionTracers, FloatFieldIJ_ConvectionTracers
from ndsl.dsl.gt4py import IJ, IJK
from ndsl.quantity.data_dimensions_field import DataDimensionsField


class UWGEOSInterface(UserCode):
    def __init__(self, name: str) -> None:
        # these must be defined during the first run call because they require the number of convection tracers
        self.config = None
        self._managed_state = None
        self._uw = None

    def init(self, mapl_state, import_state, export_state) -> None:
        # Make sure we have our NDSL stack setup
        MAPLPy = get_MAPLPy()
        ndsl_stack = get_NDSL_physics(mapl_state)

        # Initialize the configuration for UW
        if self.config is None:
            if ndsl_stack.quantity_factory.sizer.nz == 72:
                jason_uw = MAPLPy.get_resource("JASON_UW:", mapl_state, default=True)
                jason_mfd_sc = MAPLPy.get_resource("JASON_MFD_SC:", mapl_state, default=True)
            else:
                jason_uw = MAPLPy.get_resource("JASON_UW:", mapl_state, default=False)
                jason_mfd_sc = MAPLPy.get_resource("JASON_MFD_SC:", mapl_state, default=True)

        self.config = UWConfiguration(
            JASON=True if ndsl_stack.quantity_factory.sizer.nz == 72 else False,
            JASON_MFD_SC=True if ndsl_stack.quantity_factory.sizer.nz == 72 else False,
            REPORT_UW_NEGATIVES=MAPLPy.get_resource("REPORT_UW_NEGATIVES:", mapl_state, default=False),
            UW_USE_EIS=MAPLPy.get_resource("UW_USE_EIS:", mapl_state, default=False if jason_uw else True),
            NCNST=0,  # will be updated during first run call, once tracer packet is built in fortran
            rkfre=MAPLPy.get_resource("RKFRE:", mapl_state, default=Float(1.0) if jason_uw else Float(1.5)),
            k0=ndsl_stack.quantity_factory.sizer.nz,
            dotransport=1 if MAPLPy.get_resource("USE_TRACER_TRANSP_UW:", mapl_state, default=True) else 0,
            dt=MAPLPy.get_resource("DSL__UW_DT:", mapl_state, default=Float(0)),
            windsrcavg=MAPLPy.get_resource("WINDSRCAVG:", mapl_state, default=Int(0) if jason_uw else Int(1)),
            mixscale=MAPLPy.get_resource("MIXSCALE:", mapl_state, default=Float(0.0) if jason_uw else Float(3000.0)),
            criqc=MAPLPy.get_resource("CRIQC:", mapl_state, default=Float(1.0e-3) if jason_uw else Float(3.0e-3)),
            thlsrc_fac=MAPLPy.get_resource("THLSRC_FAC:", mapl_state, default=Float(0.0) if jason_uw else Float(1.0)),
            frc_rasn=MAPLPy.get_resource("FRC_RASN:", mapl_state, default=Float(0.0)),
            rkm=MAPLPy.get_resource("RKM:", mapl_state, default=Float(12.0) if jason_uw else Float(8.0)),
            rpen=MAPLPy.get_resource("RPEN:", mapl_state, default=Float(3.0) if jason_uw else Float(1.5)),
            SCLM_SHALLOW=MAPLPy.get_resource("SCLM_SHALLOW:", mapl_state, default=Float(1.0)),
            niter_xc=MAPLPy.get_resource("NITER_XC:", mapl_state, default=Int(2)),
            iter_cin=MAPLPy.get_resource("ITER_CIN:", mapl_state, default=Int(2)),
            use_CINcin=1 if MAPLPy.get_resource("USE_CINCIN:", mapl_state, default=True) else 0,
            cridist_opt=MAPLPy.get_resource("CRIDIST_OPT:", mapl_state, default=Int(0)),
            use_self_detrain=1 if MAPLPy.get_resource("USE_SELF_DETRAIN:", mapl_state, default=False) else 0,
            use_momenflx=1 if MAPLPy.get_resource("USE_MOMENFLX:", mapl_state, default=True) else 0,
            use_cumpenent=1 if MAPLPy.get_resource("USE_CUMPENENT:", mapl_state, default=True) else 0,
            rle=MAPLPy.get_resource("RLE:", mapl_state, default=Float(0.1)),
            rmaxfrac=MAPLPy.get_resource("RMAXFRAC:", mapl_state, default=Float(0.1) if jason_uw else Float(0.25)),
            mumin1=MAPLPy.get_resource("MUMIN1:", mapl_state, default=Float(0.906)),
            rbuoy=MAPLPy.get_resource("RBUOY:", mapl_state, default=Float(1.0)),
            rdrag=MAPLPy.get_resource("RDRAG:", mapl_state, default=Float(1.0)),
            epsvarw=MAPLPy.get_resource("EPSVARW:", mapl_state, default=Float(5.0e-4)),
            PGFc=MAPLPy.get_resource("PGFC:", mapl_state, default=Float(0.7)),
            rdrop=MAPLPy.get_resource("SHLW_RDROP:", mapl_state, default=Float(8.0e-6)),
            detrhgt=MAPLPy.get_resource("DETRHGT:", mapl_state, default=Float(1800.0)),
            qtsrc_fac=MAPLPy.get_resource("QTSRC_FAC:", mapl_state, default=Float(0.0)),
            qtsrchgt=MAPLPy.get_resource("QTSRCHGT:", mapl_state, default=Float(0.0)),
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

        # set correct size for convection tracer data dimension
        NUMBER_CONVECTION_TRACERS = MOIST_WORKAROUNDS.CNV_Tracers().Q[:].shape[-1]

        # create the data dimensions within the QuantityFactory
        if CONVECTION_TRACER_DIM not in ndsl_stack.quantity_factory.sizer.data_dimensions:
            ndsl_stack.quantity_factory.add_data_dimensions({CONVECTION_TRACER_DIM: MOIST_WORKAROUNDS.CNV_Tracers().Q[:].shape[-1]})
        elif ndsl_stack.quantity_factory.sizer.data_dimensions[CONVECTION_TRACER_DIM] != MOIST_WORKAROUNDS.CNV_Tracers().Q[:].shape[-1]:
            raise ValueError(f"Convection tracer count has been modified since initialization timesteps. If this is intentional, you must re-initialize the NDSL stack.")

        # register the field types (declard in convection_tracers.py) with the correct ddim size
        if not DataDimensionsField.exists("FloatFieldIJ_ConvectionTracers"):
            DataDimensionsField.register(FloatFieldIJ_ConvectionTracers, ndsl_stack.quantity_factory, [CONVECTION_TRACER_DIM], axes=IJ, dtype=Float)
        if not DataDimensionsField.exists("FloatField_ConvectionTracers"):
            DataDimensionsField.register(FloatField_ConvectionTracers, ndsl_stack.quantity_factory, [CONVECTION_TRACER_DIM], axes=IJK, dtype=Float)

        if self._managed_state is None:
            self.config.NCNST = NUMBER_CONVECTION_TRACERS
            # Initialize NDSL state
            self._managed_state = MAPLManagedState(
                UWState.empty(
                    ndsl_stack.quantity_factory,
                    data_dimensions=ndsl_stack.quantity_factory.sizer.data_dimensions,
                ),
                ndsl_stack.interface_type,
            )

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

        self._managed_state.register_2D("output.plcl_out", "PLCL_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.plfc_out", "PLFC_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.pinv_out", "PINV_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.prel_out", "PREL_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.pbup_out", "PBUP_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.cbmf_out", "CBMF_SC", export_repository, alloc=True)
        self._managed_state.register_2D("output.cldhgt_out", "CLDTOP_SC", export_repository, alloc=True)

        self._managed_state.register_2D("output.LTS", "LTS", export_repository, alloc=True)
        self._managed_state.register_2D("output.EIS", "EIS", export_repository, alloc=True)
        self._managed_state.register("output.MFD_SC", "MFD_SC", export_repository, alloc=True)
        self._managed_state.register("output.QLENT_SC", "QLENT_SC", export_repository, alloc=True)
        self._managed_state.register("output.QIENT_SC", "QIENT_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface("output.umf_inv", "UMF_SC", export_repository, alloc=True)
        self._managed_state.register("output.dcm_inv", "DCM_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface("output.qtflx_inv", "QTFLX_SC", export_repository, alloc=True)
        self._managed_state.register_K_interface("output.slflx_inv", "SLFLX_SC", export_repository, alloc=True)
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
        #self._managed_state.register_K_interface("output.CNV_MFC", "CNV_MFC", export_repository, alloc=True)
        #self._managed_state.register("output.CNV_MFD", "CNV_MFD", export_repository, alloc=True)
        self._managed_state.register("output.SHLW_PRC3", "SHLW_PRC3", export_repository, alloc=True)
        self._managed_state.register("output.SHLW_SNO3", "SHLW_SNO3", export_repository, alloc=True)
        self._managed_state.register_2D("output.SC_QT", "SC_QT", export_repository)
        self._managed_state.register_2D("output.SC_MSE", "SC_MSE", export_repository)
        self._managed_state.register_2D("output.CUSH_SC", "CUSH_SC", export_repository)
        self._managed_state.register("input_output.CLCN", "CLCN", internal_repository)

        self._managed_state.register("output.DQVDT_FILL", "DQVDT_FILL_SC", export_repository)
        self._managed_state.register("output.DQLLSDT_FILL", "DQLLSDT_FILL_SC", export_repository)
        self._managed_state.register("output.DQLCNDT_FILL", "DQLCNDT_FILL_SC", export_repository)
        self._managed_state.register("output.DQILSDT_FILL", "DQILSDT_FILL_SC", export_repository)
        self._managed_state.register("output.DQICNDT_FILL", "DQICNDT_FILL_SC", export_repository)

        if self._uw is None:
            # Build UW
            with StencilBackendCompilerOverride(
                MPI.COMM_WORLD,
                ndsl_stack.stencil_factory.config.dace_config,
            ):
                self._uw = ComputeUwshcuInv(ndsl_stack.stencil_factory, ndsl_stack.quantity_factory, self.config)

        with TimedCUDAProfiler("UW", {}):
            with TimedCUDAProfiler("UW - State copy", {}):
                self._managed_state.fortran_to_ndsl()
                safe_assign_array(
                    self._managed_state.ndsl_state.input_output.CNV_Tracers.field[:],
                    MOIST_WORKAROUNDS.CNV_Tracers().Q,
                )

            with TimedCUDAProfiler("UW Numerics", {}):
                self._uw(self._managed_state.ndsl_state)

            with TimedCUDAProfiler("UW - State copy-back", {}):
                safe_assign_array(
                    MOIST_WORKAROUNDS.CNV_Tracers().Q,
                    self._managed_state.ndsl_state.input_output.CNV_Tracers.field[:],
                )
                self._managed_state.ndsl_to_fortran()

    def finalize(self, mapl_state, import_state, export_state) -> None:
        self._managed_state.save_recorded()


CODE = UWGEOSInterface("UW")
