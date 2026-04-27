from typing import Any

from MAPL_PythonBridge import UserCode, get_MAPLPy
from MAPL_PythonBridge.types import CVoidPointer
from mpi4py import MPI
from ndsl.constants import I_DIM, J_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int
from ndsl.utils import safe_assign_array

from pyMoist.constants import NUMBER_OF_TRACERS
from pyMoist.convection.GF_2020 import GF2020, GF2020Config, GF2020CumulusParameterizationConfig, GF2020State
from pyMoist.convection_tracers import ConvectionTracers
from pyMoist.fortran import get_NDSL_physics
from pyMoist.fortran.build_helper import StencilBackendCompilerOverride
from pyMoist.fortran.managed_state import MAPLManagedState
from pyMoist.fortran.memory_factory import MAPLMemoryRepository
from pyMoist.fortran.moist_workarounds import MOIST_WORKAROUNDS
from pyMoist.fortran.profiler import TimedCUDAProfiler
from pyMoist.saturation_tables import SaturationVaporPressureTable


def _default_or_get_from_namelist(default, name_in_namelist: str, namelist: dict[str, Any]) -> Any:
    return default if name_in_namelist not in namelist else namelist[name_in_namelist]


class GF2020Interface(UserCode):
    def __init__(self) -> None:
        pass

    def init(self, mapl_state: CVoidPointer, import_state: CVoidPointer, export_state: CVoidPointer):
        maplpy = get_MAPLPy()
        ndsl_stack = get_NDSL_physics(mapl_state)

        gf_2020_env_setting = maplpy.get_resource("DSL__GF_ENV_SETTING", mapl_state, default=Int(1))

        zero_diff = maplpy.get_resource("ZERO_DIFF:", mapl_state, default=Int(0))
        hydrostatic = maplpy.get_resource("HYDROSTATIC:", mapl_state, default=True)

        sh_md_dp = maplpy.get_resource("SH_MD_DP:", mapl_state, default=True)

        # NOTE for some reason, a bunch of these fields come back as 64 bit, despite MAPLPy detecting the
        # correct 32 bit precision and calling the correct portion of the F to Py bridge
        # to be safe, all config constants have an extra cast
        if zero_diff == 0:
            entrversion = maplpy.get_resource("ENTRVERSION:", mapl_state, default=Int(0))

            # for cumulus parameterizaiton
            gf_min_area = Float(maplpy.get_resource("GF_MIN_AREA:", mapl_state, default=Float(0.0)))
            tau_mid = Int(maplpy.get_resource("TAU_MID:", mapl_state, default=Int(5400)))
            tau_deep = Int(maplpy.get_resource("TAU_DEEP:", mapl_state, default=Int(10800)))
            downdraft_max_height_land_shallow = Float(maplpy.get_resource("HEI_DOWN_LAND_SH:", mapl_state, default=Float(0.0)))
            downdraft_max_height_land_mid = Float(maplpy.get_resource("HEI_DOWN_LAND_MD:", mapl_state, default=Float(0.3)))
            downdraft_max_height_land_deep = Float(maplpy.get_resource("HEI_DOWN_LAND_DP:", mapl_state, default=Float(0.3)))
            downdraft_max_height_ocean_shallow = Float(maplpy.get_resource("HEI_DOWN_OCEAN_SH:", mapl_state, default=Float(0.0)))
            downdraft_max_height_ocean_mid = Float(maplpy.get_resource("HEI_DOWN_OCEAN_MD:", mapl_state, default=Float(0.6)))
            downdraft_max_height_ocean_deep = Float(maplpy.get_resource("HEI_DOWN_OCEAN_DP:", mapl_state, default=Float(0.6)))
            updraft_max_height_land_shallow = Float(maplpy.get_resource("HEI_UPDF_LAND_SH:", mapl_state, default=Float(0.2)))
            updraft_max_height_land_mid = Float(maplpy.get_resource("HEI_UPDF_LAND_MD:", mapl_state, default=Float(0.4)))
            updraft_max_height_land_deep = Float(maplpy.get_resource("HEI_UPDF_LAND_DP:", mapl_state, default=Float(0.4)))
            updraft_max_height_ocean_shallow = Float(maplpy.get_resource("HEI_UPDF_OCEAN_SJ:", mapl_state, default=Float(0.2)))
            updraft_max_height_ocean_mid = Float(maplpy.get_resource("HEI_UPDF_OCEAN_MD:", mapl_state, default=Float(0.4)))
            updraft_max_height_ocean_deep = Float(maplpy.get_resource("HEI_UPDF_OCEAN_DP:", mapl_state, default=Float(0.4)))
            minimum_evap_fraction_land_shallow = Float(maplpy.get_resource("MIN_EDT_LAND_SH:", mapl_state, default=Float(0.0)))
            minimum_evap_fraction_land_mid = Float(maplpy.get_resource("MIN_EDT_LAND_MD:", mapl_state, default=Float(0.1)))
            minimum_evap_fraction_land_deep = Float(maplpy.get_resource("MIN_EDT_LAND_DP:", mapl_state, default=Float(0.1)))
            minimum_evap_fraction_ocean_shallow = Float(maplpy.get_resource("MIN_EDT_OCEAN_SH:", mapl_state, default=Float(0.0)))
            minimum_evap_fraction_ocean_mid = Float(maplpy.get_resource("MIN_EDT_OCEAN_MD:", mapl_state, default=Float(0.1)))
            minimum_evap_fraction_ocean_deep = Float(maplpy.get_resource("MIN_EDT_OCEAN_DP:", mapl_state, default=Float(0.1)))
            maximum_evap_fraction_land_shallow = Float(maplpy.get_resource("MAX_EDT_LAND_SH:", mapl_state, default=Float(0.0)))
            maximum_evap_fraction_land_mid = Float(maplpy.get_resource("MAX_EDT_LAND_MD:", mapl_state, default=Float(0.4)))
            maximum_evap_fraction_land_deep = Float(maplpy.get_resource("MAX_EDT_LAND_DP:", mapl_state, default=Float(0.4)))
            maximum_evap_fraction_ocean_shallow = Float(maplpy.get_resource("MAX_EDT_OCEAN_SH:", mapl_state, default=Float(0.0)))
            maximum_evap_fraction_ocean_mid = Float(maplpy.get_resource("MAX_EDT_OCEAN_MD:", mapl_state, default=Float(0.3)))
            maximum_evap_fraction_ocean_deep = Float(maplpy.get_resource("MAX_EDT_OCEAN_DP:", mapl_state, default=Float(0.3)))
            if hydrostatic:
                sgs_w_timescale = Int(maplpy.get_resource("SGS_W_TIMESCALE:", mapl_state, default=Int(0)))
            else:
                sgs_w_timescale = Int(maplpy.get_resource("SGS_W_TIMESCALE:", mapl_state, default=Int(1)))
            min_entrainment_rate = Float(maplpy.get_resource("MIN_ENTR_RATE:", mapl_state, default=Float(0.3e-4)))
            entrainment_rate_shallow = Float(maplpy.get_resource("ENTR_SH:", mapl_state, default=Float(1.0e-4)))
            entrainment_rate_mid = Float(maplpy.get_resource("ENTR_MD:", mapl_state, default=Float(1.0e-4)))
            entrainment_rate_deep = Float(maplpy.get_resource("ENTR_DP:", mapl_state, default=Float(1.0e-4)))
        else:
            entrversion = Int(maplpy.get_resource("ENTRVERSION:", mapl_state, default=Int(1)))

            # for cumulus parameterization
            gf_min_area = Float(maplpy.get_resource("GF_MIN_AREA:", mapl_state, default=Float(1.0e6)))
            tau_mid = Int(maplpy.get_resource("TAU_MID:", mapl_state, default=Int(3600)))
            tau_deep = Int(maplpy.get_resource("TAU_DEEP:", mapl_state, default=Int(5400)))
            downdraft_max_height_land_shallow = Float(maplpy.get_resource("HEI_DOWN_LAND_SH:", mapl_state, default=Float(0.0)))
            downdraft_max_height_land_mid = Float(maplpy.get_resource("HEI_DOWN_LAND_MD:", mapl_state, default=Float(0.5)))
            downdraft_max_height_land_deep = Float(maplpy.get_resource("HEI_DOWN_LAND_DP:", mapl_state, default=Float(0.5)))
            downdraft_max_height_ocean_shallow = Float(maplpy.get_resource("HEI_DOWN_OCEAN_SH:", mapl_state, default=Float(0.0)))
            downdraft_max_height_ocean_mid = Float(maplpy.get_resource("HEI_DOWN_OCEAN_MD:", mapl_state, default=Float(0.5)))
            downdraft_max_height_ocean_deep = Float(maplpy.get_resource("HEI_DOWN_OCEAN_DP:", mapl_state, default=Float(0.5)))
            updraft_max_height_land_shallow = Float(maplpy.get_resource("HEI_UPDF_LAND_SH:", mapl_state, default=Float(0.2)))
            updraft_max_height_land_mid = Float(maplpy.get_resource("HEI_UPDF_LAND_MD:", mapl_state, default=Float(0.65)))
            updraft_max_height_land_deep = Float(maplpy.get_resource("HEI_UPDF_LAND_DP:", mapl_state, default=Float(0.65)))
            updraft_max_height_ocean_shallow = Float(maplpy.get_resource("HEI_UPDF_OCEAN_SJ:", mapl_state, default=Float(0.2)))
            updraft_max_height_ocean_mid = Float(maplpy.get_resource("HEI_UPDF_OCEAN_MD:", mapl_state, default=Float(0.65)))
            updraft_max_height_ocean_deep = Float(maplpy.get_resource("HEI_UPDF_OCEAN_DP:", mapl_state, default=Float(0.65)))
            minimum_evap_fraction_land_shallow = Float(maplpy.get_resource("MIN_EDT_LAND_SH:", mapl_state, default=Float(0.0)))
            minimum_evap_fraction_land_mid = Float(maplpy.get_resource("MIN_EDT_LAND_MD:", mapl_state, default=Float(0.1)))
            minimum_evap_fraction_land_deep = Float(maplpy.get_resource("MIN_EDT_LAND_DP:", mapl_state, default=Float(0.1)))
            minimum_evap_fraction_ocean_shallow = Float(maplpy.get_resource("MIN_EDT_OCEAN_SH:", mapl_state, default=Float(0.0)))
            minimum_evap_fraction_ocean_mid = Float(maplpy.get_resource("MIN_EDT_OCEAN_MD:", mapl_state, default=Float(0.1)))
            minimum_evap_fraction_ocean_deep = Float(maplpy.get_resource("MIN_EDT_OCEAN_DP:", mapl_state, default=Float(0.1)))
            maximum_evap_fraction_land_shallow = Float(maplpy.get_resource("MAX_EDT_LAND_SH:", mapl_state, default=Float(0.0)))
            maximum_evap_fraction_land_mid = Float(maplpy.get_resource("MAX_EDT_LAND_MD:", mapl_state, default=Float(0.9)))
            maximum_evap_fraction_land_deep = Float(maplpy.get_resource("MAX_EDT_LAND_DP:", mapl_state, default=Float(0.9)))
            maximum_evap_fraction_ocean_shallow = Float(maplpy.get_resource("MAX_EDT_OCEAN_SH:", mapl_state, default=Float(0.0)))
            maximum_evap_fraction_ocean_mid = Float(maplpy.get_resource("MAX_EDT_OCEAN_MD:", mapl_state, default=Float(0.9)))
            maximum_evap_fraction_ocean_deep = Float(maplpy.get_resource("MAX_EDT_OCEAN_DP:", mapl_state, default=Float(0.9)))
            sgs_w_timescale = Int(maplpy.get_resource("SGS_W_TIMESCALE:", mapl_state, default=Int(0)))
            min_entrainment_rate = Float(maplpy.get_resource("MIN_ENTR_RATE:", mapl_state, default=Float(0.1e-4)))
            entrainment_rate_shallow = Float(maplpy.get_resource("ENTR_SH:", mapl_state, default=Float(1.0e-4)))
            entrainment_rate_mid = Float(maplpy.get_resource("ENTR_MD:", mapl_state, default=Float(9.0e-4)))
            entrainment_rate_deep = Float(maplpy.get_resource("ENTR_DP:", mapl_state, default=Float(1.0e-3)))

        config = GF2020Config(
            DT_MOIST=Float(maplpy.get_resource("DSL__GF2020_DT", mapl_state, default=Float(0.0))),
            LHYDROSTATIC=hydrostatic,
            STOCHASTIC_CNV=maplpy.get_resource("STOCHASTIC_CNV:", mapl_state, default=False),
            STOCH_TOP=Float(maplpy.get_resource("STOCH_TOP:", mapl_state, default=Float(2.5))),
            STOCH_BOT=Float(maplpy.get_resource("STOCH_BOT:", mapl_state, default=Float(0.75))),
            GF_MIN_AREA=gf_min_area,
            GF_ENV_SETTING=gf_2020_env_setting,
            ENTRVERSION=entrversion,
            CONVECTION_TRACER=Int(maplpy.get_resource("CONVECTION_TRACER:", mapl_state, default=Int(0))),
            C1=Float(maplpy.get_resource("C1:", mapl_state, default=Float(0.0))),
            ADV_TRIGGER=Int(maplpy.get_resource("ADV_TRIGGER:", mapl_state, default=Int(0))),
            AUTOCONV=Int(maplpy.get_resource("AUTOCONV:", mapl_state, default=Int(1))),
            USE_TRACER_TRANSPORT=Int(maplpy.get_resource("USE_TRACER_TRANSP:", mapl_state, default=Int(1))),
            SCLM_DEEP=Float(maplpy.get_resource("SCLM_DEEP:", mapl_state, default=Float(1.0))),
            FIX_CONVECTIVE_CLOUD=maplpy.get_resource("FIX_CNV_CLOUD:", mapl_state, default=False),
            APPLY_SUBSIDENCE_MICROPHYSICS=Int(maplpy.get_resource("APPLY_SUB_MP:", mapl_state, default=Int(0))),
            NUMBER_OF_TRACERS=NUMBER_OF_TRACERS,
            USE_MOMENTUM_TRANSPORT=Int(maplpy.get_resource("USE_MOMENTUM_TRANSP:", mapl_state, default=Int(1))),
        )

        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(
            # plume dependent
            DOWNDRAFT_MAX_HEIGHT_LAND_SHALLOW=downdraft_max_height_land_shallow,
            DOWNDRAFT_MAX_HEIGHT_LAND_MID=downdraft_max_height_land_mid,
            DOWNDRAFT_MAX_HEIGHT_LAND_DEEP=downdraft_max_height_land_deep,
            DOWNDRAFT_MAX_HEIGHT_OCEAN_SHALLOW=downdraft_max_height_ocean_shallow,
            DOWNDRAFT_MAX_HEIGHT_OCEAN_MID=downdraft_max_height_ocean_mid,
            DOWNDRAFT_MAX_HEIGHT_OCEAN_DEEP=downdraft_max_height_ocean_deep,
            UPDRAFT_MAX_HEIGHT_LAND_SHALLOW=updraft_max_height_land_shallow,
            UPDRAFT_MAX_HEIGHT_LAND_MID=updraft_max_height_land_mid,
            UPDRAFT_MAX_HEIGHT_LAND_DEEP=updraft_max_height_land_deep,
            UPDRAFT_MAX_HEIGHT_OCEAN_SHALLOW=updraft_max_height_ocean_shallow,
            UPDRAFT_MAX_HEIGHT_OCEAN_MID=updraft_max_height_ocean_mid,
            UPDRAFT_MAX_HEIGHT_OCEAN_DEEP=updraft_max_height_ocean_deep,
            MINIMUM_EVAP_FRACTION_LAND_SHALLOW=minimum_evap_fraction_land_shallow,
            MINIMUM_EVAP_FRACTION_LAND_MID=minimum_evap_fraction_land_mid,
            MINIMUM_EVAP_FRACTION_LAND_DEEP=minimum_evap_fraction_land_deep,
            MINIMUM_EVAP_FRACTION_OCEAN_SHALLOW=minimum_evap_fraction_ocean_shallow,
            MINIMUM_EVAP_FRACTION_OCEAN_MID=minimum_evap_fraction_ocean_mid,
            MINIMUM_EVAP_FRACTION_OCEAN_DEEP=minimum_evap_fraction_ocean_deep,
            MAXIMUM_EVAP_FRACTION_LAND_SHALLOW=maximum_evap_fraction_land_shallow,
            MAXIMUM_EVAP_FRACTION_LAND_MID=maximum_evap_fraction_land_mid,
            MAXIMUM_EVAP_FRACTION_LAND_DEEP=maximum_evap_fraction_land_deep,
            MAXIMUM_EVAP_FRACTION_OCEAN_SHALLOW=maximum_evap_fraction_ocean_shallow,
            MAXIMUM_EVAP_FRACTION_OCEAN_MID=maximum_evap_fraction_ocean_mid,
            MAXIMUM_EVAP_FRACTION_OCEAN_DEEP=maximum_evap_fraction_ocean_deep,
            CLOUD_BASE_MASS_FLUX_FACTOR_SHALLOW=Float(maplpy.get_resource("FADJ_MASSFLX_SH:", mapl_state, default=Float(1.0))),
            CLOUD_BASE_MASS_FLUX_FACTOR_MID=Float(maplpy.get_resource("FADJ_MASSFLX_MD:", mapl_state, default=Float(1.0))),
            CLOUD_BASE_MASS_FLUX_FACTOR_DEEP=Float(maplpy.get_resource("FADJ_MASSFLX_DP:", mapl_state, default=Float(1.0))),
            USE_EXCESS_SHALLOW=Int(maplpy.get_resource("USE_EXCESS_SH:", mapl_state, default=Int(3))),
            USE_EXCESS_MID=Int(maplpy.get_resource("USE_EXCESS_MD:", mapl_state, default=Int(2))),
            USE_EXCESS_DEEP=Int(maplpy.get_resource("USE_EXCESS_DP:", mapl_state, default=Int(2))),
            AVERAGE_LAYER_DEPTH_SHALLOW=Float(maplpy.get_resource("AVE_LAYER_SH:", mapl_state, default=Float(30.0))),
            AVERAGE_LAYER_DEPTH_MID=Float(maplpy.get_resource("AVE_LAYER_MD:", mapl_state, default=Float(40.0))),
            AVERAGE_LAYER_DEPTH_DEEP=Float(maplpy.get_resource("AVE_LAYER_DP:", mapl_state, default=Float(40.0))),
            ENABLE_SHALLOW=Int(maplpy.get_resource("SHALLOW:", mapl_state, default=Int(0))),
            ENABLE_MID=Int(maplpy.get_resource("CONGESTUS:", mapl_state, default=Int(1))),
            ENABLE_DEEP=Int(maplpy.get_resource("DEEP:", mapl_state, default=Int(1))),
            ENTRAINMENT_RATE_SHALLOW=entrainment_rate_shallow,
            ENTRAINMENT_RATE_MID=entrainment_rate_mid,
            ENTRAINMENT_RATE_DEEP=entrainment_rate_deep,
            C0_SHAL=Float(maplpy.get_resource("C0_SHAL:", mapl_state, default=Float(0.0))),
            C0_MID=Float(maplpy.get_resource("C0_MID:", mapl_state, default=Float(2.0e-3))),
            C0_DEEP=Float(maplpy.get_resource("C0_DEEP:", mapl_state, default=Float(2.0e-3))),
            TAU_MID=tau_mid,
            TAU_DEEP=tau_deep,
            CLOSURE_CHOICE_SHALLOW=Int(maplpy.get_resource("CLOSURE_SHALLOW:", mapl_state, default=Int(7))),
            CLOSURE_CHOICE_MID=Int(maplpy.get_resource("CLOSURE_CONGESTUS:", mapl_state, default=Int(3))),
            CLOSURE_CHOICE_DEEP=Int(maplpy.get_resource("CLOSURE_DEEP:", mapl_state, default=Int(0))),
            # plume independent
            SHALLOW_MID_DEEP=maplpy.get_resource("SH_MD_DP:", mapl_state, default=True),
            ZERO_DIFF=zero_diff,
            MOIST_TRIGGER=Int(maplpy.get_resource("MOIST_TRIGGER:", mapl_state, default=Int(0))),
            LAMBDA_DEEP=Float(maplpy.get_resource("LAMBDAU_DEEP:", mapl_state, default=Float(0.0))),
            LAMBDA_SHALLOW_DOWN=Float(maplpy.get_resource("LAMBAU_SHDN:", mapl_state, default=Float(2.0))),
            CAP_MAXS=Float(maplpy.get_resource("CAP_MAXS:", mapl_state, default=Float(50.0))),
            OUTPUT_SOUNDING=Int(maplpy.get_resource("OUTPUT_SOUND:", mapl_state, default=Int(0))),
            USE_SCALE_DEP=Int(maplpy.get_resource("USE_SCALE_DEP:", mapl_state, default=Int(1))),
            SATURATION_CALCULATION_CHOICE=Int(maplpy.get_resource("SATUR_CALC:", mapl_state, default=Int(1))),
            CLOUD_LEVEL_GRID=Int(maplpy.get_resource("CLEV_GRID:", mapl_state, default=Int(1))),
            FRAC_MODIS=Int(maplpy.get_resource("FRAC_MODIS:", mapl_state, default=Int(1))),
            BOUNDARY_CONDITION_METHOD=Int(maplpy.get_resource("BC_METH:", mapl_state, default=Int(1))),
            OVERSHOOT=Float(maplpy.get_resource("OVERSHOOT:", mapl_state, default=Float(0.0))),
            USE_MEMORY=Int(maplpy.get_resource("USE_MEMORY:", mapl_state, default=Int(-1))),
            DOWNDRAFT=Int(maplpy.get_resource("DOWNDRAFT:", mapl_state, default=Int(1))),
            USE_WETBULB=Int(maplpy.get_resource("USE_WETBULB:", mapl_state, default=Int(0))),
            DIURNAL_CYCLE=Int(maplpy.get_resource("DICYCL:", mapl_state, default=Int(1))),
            USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES=Int(maplpy.get_resource("USE_LINEAR_SUBCL_MF:", mapl_state, default=Int(0))),
            CRITICAL_MIXING_RATIO_OVER_OCEAN=Float(maplpy.get_resource("QRC_CRIT_OCN:", mapl_state, default=Float(2.0e-4))),
            CRITICAL_MIXING_RATIO_OVER_LAND=Float(maplpy.get_resource("QRC_CRIT_LND:", mapl_state, default=Float(2.0e-4))),
            BETA_SHALLOW=Float(maplpy.get_resource("BETA_SH:", mapl_state, default=Float(2.2))),
            EVAP_FIX=Int(maplpy.get_resource("EVAP_FIX:", mapl_state, default=Int(1))),
            SGS_W_TIMESCALE=sgs_w_timescale,
            VERTICAL_DISCRETIZATION_OPTION=Int(maplpy.get_resource("VERT_DISCR:", mapl_state, default=Int(1))),
            ALP1=Float(maplpy.get_resource("ALP1:", mapl_state, default=Float(1.0))),
            USE_FCT=Int(maplpy.get_resource("USE_FCT:", mapl_state, default=Int(0))),
            MIN_ENTRAINMENT_RATE=min_entrainment_rate,
            USE_SMOOTH_TENDENCIES=Int(maplpy.get_resource("USE_SMOOTH_TEND:", mapl_state, default=Int(1))),
            USE_RAIN_EVAP_BELOW_CLOUD_BASE=Int(maplpy.get_resource("USE_REBCB:", mapl_state, default=Int(1))),
            USE_CLOUD_DISSIPATION=Float(maplpy.get_resource("USE_CLOUD_DISSIPATION:", mapl_state, default=Float(1.0))),
            LIGHTNING_DIAGNOSTICS=Int(maplpy.get_resource("LIGHTNING_DIAG:", mapl_state, default=Int(0))),
            USE_TRACER_SCAVENGE=Int(maplpy.get_resource("USE_TRACER_SCAVEN:", mapl_state, default=Int(1))),
            USE_TRACER_EVAPORATION=Int(maplpy.get_resource("USE_TRACER_EVAP:", mapl_state, default=Int(1))),
            USE_FLUX_FORM=Int(maplpy.get_resource("USE_FLUX_FORM:", mapl_state, default=Int(1))),
            MAX_TEMP_VAPOR_TENDENCY=Float(maplpy.get_resource("MAX_TQ_TEND:", mapl_state, default=Float(100.0))),
        )

        saturation_tables = SaturationVaporPressureTable(ndsl_stack.stencil_factory.backend)

        # Initialize the module
        with StencilBackendCompilerOverride(
            MPI.COMM_WORLD,
            ndsl_stack.stencil_factory.config.dace_config,
        ):
            self._gf_2020 = GF2020(
                stencil_factory=ndsl_stack.stencil_factory,
                quantity_factory=ndsl_stack.quantity_factory,
                config=config,
                cumulus_parameterization_config=cumulus_parameterization_config,
                saturation_tables=saturation_tables,
            )

        # Make the state
        self._managed_state = MAPLManagedState(
            GF2020State.empty(ndsl_stack.quantity_factory),
            ndsl_stack.interface_type,
        )

        self._managed_convection_tracers = MAPLManagedState(
            ConvectionTracers.empty(
                ndsl_stack.quantity_factory,
                data_dimensions={
                    "convection_tracers": config.NUMBER_OF_TRACERS,
                    "size_three_dimension": 3,
                    "size_four_dimension": 4,
                },
            ),
            ndsl_stack.interface_type,
        )

    def run(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
    ):
        pass

    def run_with_internal(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
        internal_state: CVoidPointer,
    ):
        ndsl_stack = get_NDSL_physics(mapl_state)
        import_repository = MAPLMemoryRepository(import_state, ndsl_stack.quantity_factory)
        internal_repository = MAPLMemoryRepository(internal_state, ndsl_stack.quantity_factory)
        export_repository = MAPLMemoryRepository(export_state, ndsl_stack.quantity_factory)

        self._managed_state.register("latitude", "DSL__GF2020_LATS", internal_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("longitude", "DSL__GF2020_LONS", internal_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("p_interface", "PLE", import_repository, dims=[I_DIM, J_DIM, K_INTERFACE_DIM])
        self._managed_state.register("t", "T", import_repository)
        self._managed_state.register("u", "U", import_repository)
        self._managed_state.register("v", "V", import_repository)
        self._managed_state.register("w", "W", import_repository)
        self._managed_state.register("omega", "OMEGA", import_repository)
        self._managed_state.register("t_2m", "T2M", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("specific_humidity_2m", "Q2M", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("t_surface", "TA", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("specific_humidity_surface", "QA", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("vapor", "Q", internal_repository)
        self._managed_state.register("convective_liquid", "QLCN", internal_repository)
        self._managed_state.register("convective_ice", "QICN", internal_repository)
        self._managed_state.register("convective_cloud_fraction", "CLCN", internal_repository)
        self._managed_state.register("large_scale_liquid", "QLLS", internal_repository)
        self._managed_state.register("large_scale_ice", "QILS", internal_repository)
        self._managed_state.register("large_scale_cloud_fraction", "CLLS", internal_repository)
        self._managed_state.register("ice_fraction_in_convective_tower", "CNV_FICE", export_repository)
        self._managed_state.register(
            "p_interface_timestep_start",
            "PLE_DYN_IN",
            import_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
        )
        self._managed_state.register("t_timestep_start", "T_DYN_IN", import_repository)
        self._managed_state.register("u_timestep_start", "U_DYN_IN", import_repository)
        self._managed_state.register("v_timestep_start", "V_DYN_IN", import_repository)
        self._managed_state.register("vapor_timestep_start", "QV_DYN_IN", import_repository)
        self._managed_state.register("geopotential_height_interface", "ZLE", import_repository, dims=[I_DIM, J_DIM, K_INTERFACE_DIM])
        self._managed_state.register("geopotential_height_surface", "PHIS", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("area", "AREA", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("pbl_level", "KPBL", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("convection_fraction", "CNV_FRC", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("surface_type", "SRF_TYPE", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("seed_convection", "STOCH_CNV", export_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("land_fraction", "FRLAND", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("scalar_diffusivity", "KH", import_repository, dims=[I_DIM, J_DIM, K_INTERFACE_DIM])
        self._managed_state.register("buoyancy", "BYNCY", export_repository, alloc=True)
        self._managed_state.register("convective_precipitation_GF", "CNPCPRATE", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("convective_precipitation_RAS", "CNV_PRC3", export_repository, alloc=True)
        self._managed_state.register("convective_rainwater_source", "DQRC", export_repository)
        self._managed_state.register("sensible_heat_flux", "SH", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register(
            "total_water_flux_deep_convection_interface",
            "WQT_DC",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register("sublimation_of_convective_precipitation", "RSU_CN", export_repository, alloc=True)
        self._managed_state.register("evaporation_of_convective_precipitation", "REV_CN", export_repository, alloc=True)
        self._managed_state.register(
            "ice_precip_flux_interface",
            "PFI_CN",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "liquid_precip_flux_interface",
            "PFL_CN",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register("evaporation", "EVAP", import_repository, dims=[I_DIM, J_DIM])
        self._managed_state.register("convective_condensate_source", "CNV_DQCDT", export_repository, alloc=True)
        self._managed_state.register("convective_condensate_grid_mean", "CNV_QC", export_repository, alloc=True)
        self._managed_state.register("entrainment_parameter", "ENTLAM", export_repository, alloc=True)
        self._managed_state.register("lateral_entrainment_rate", "ENTR", export_repository, alloc=True)
        self._managed_state.register("lateral_entrainment_rate_shallow", "ENTR_SH", export_repository, alloc=True)
        self._managed_state.register("lateral_entrainment_rate_mid", "ENTR_MD", export_repository, alloc=True)
        self._managed_state.register("lateral_entrainment_rate_deep", "ENTR_DP", export_repository, alloc=True)
        self._managed_state.register("updraft_areal_fraction", "CNV_UPDF", export_repository, alloc=True)
        self._managed_state.register("updraft_vertical_velocity", "CNV_CVW", export_repository, alloc=True)
        self._managed_state.register("dtdt_shortwave", "RADSW", import_repository)
        self._managed_state.register("dtdt_longwave", "RADLW", import_repository)
        self._managed_state.register("dspecific_humiditydt_pbl", "DQDT_BL", import_repository)
        self._managed_state.register("dtdt_pbl", "DTDT_BL", import_repository)
        self._managed_state.register("dtdt_from_dynamics", "DTDTDYN", import_repository)
        self._managed_state.register("dvapordt_from_dynamics", "DQVDTDYN", import_repository)
        self._managed_state.register("sigma_mid", "SIGMA_MID", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("sigma_deep", "SIGMA_DEEP", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("total_precipitable_water_initial", "TPWI", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register(
            "saturation_total_precipitable_water_initial",
            "TPWI_star",
            export_repository,
            dims=[I_DIM, J_DIM],
            alloc=True,
        )
        self._managed_state.register("dvapordt_deep_convection", "DQVDT_DC", export_repository, alloc=True)
        self._managed_state.register("dtdt_deep_convection", "DTDT_DC", export_repository, alloc=True)
        self._managed_state.register("dudt_deep_convection", "DUDT_DC", export_repository, alloc=True)
        self._managed_state.register("dvdt_deep_convection", "DVDT_DC", export_repository, alloc=True)
        self._managed_state.register("dliquiddt_deep_convection", "DQLDT_DC", export_repository, alloc=True)
        self._managed_state.register("dicedt_deep_convection", "DQIDT_DC", export_repository, alloc=True)
        self._managed_state.register("dcloudfractiondt_deep_convection", "DQADT_DC", export_repository, alloc=True)
        self._managed_state.register(
            "pressure_shallow_convective_cloud_top",
            "CNV_TOPP_SH",
            export_repository,
            dims=[I_DIM, J_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "pressure_mid_convective_cloud_top",
            "CNV_TOPP_MD",
            export_repository,
            dims=[I_DIM, J_DIM],
            alloc=True,
        )
        self._managed_state.register(
            "pressure_deep_convective_cloud_top",
            "CNV_TOPP_DP",
            export_repository,
            dims=[I_DIM, J_DIM],
            alloc=True,
        )
        self._managed_state.register("mass_flux_shallow", "MUPSH", export_repository, alloc=True)
        self._managed_state.register("mass_flux_mid", "MUPMD", export_repository, alloc=True)
        self._managed_state.register("mass_flux_deep_updraft", "MUPDP", export_repository, alloc=True)
        self._managed_state.register(
            "mass_flux_deep_updraft_interface",
            "UMF_DC",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register("mass_flux_deep_updraft_detrained", "MFD_DC", export_repository, alloc=True)
        self._managed_state.register("mass_flux_deep_downdraft", "MDNDP", export_repository, alloc=True)
        self._managed_state.register("mass_flux_cloud_base", "CNV_MF0", export_repository, alloc=True)
        self._managed_state.register("mass_flux_cloud_base_shallow", "MFSH", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("mass_flux_cloud_base_mid", "MFMD", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("mass_flux_cloud_base_deep", "MFDP", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register(
            "total_cumulative_mass_flux_interface",
            "CNV_MFC",
            export_repository,
            dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            alloc=True,
        )
        self._managed_state.register("total_detraining_mass_flux", "CNV_MFD", export_repository, alloc=True)
        self._managed_state.register("convection_code_shallow", "ERRSH", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("convection_code_mid", "ERRMD", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("convection_code_deep", "ERRDP", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cloud_workfunction_0", "AA0", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cloud_workfunction_1", "AA1", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cloud_workfunction_2", "AA2", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cloud_workfunction_3", "AA3", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cloud_workfunction_1_pbl", "AA1_BL", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cloud_workfunction_1_cin", "AA1_CIN", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("pbl_time_scale", "TAU_BL", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("cape_removal_time_scale", "TAU_EC", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("lightning_density", "LFR_GF", export_repository, dims=[I_DIM, J_DIM], alloc=True)
        self._managed_state.register("convection_tracer", "CNV_TR", internal_repository)

        if self._gf_2020 is None:
            raise RuntimeError("GF2020 Runtime called before initialization was done. Abort.")

        with TimedCUDAProfiler("GF 2020 Convection", {}):
            with TimedCUDAProfiler("GF 2020 Convection - State copy", {}):
                self._managed_state.fortran_to_ndsl()
                safe_assign_array(
                    self._managed_convection_tracers.ndsl_state.tracers.data[:],
                    MOIST_WORKAROUNDS.CNV_Tracers().Q[:],
                )
                safe_assign_array(
                    self._managed_convection_tracers.ndsl_state.fscav.data[:],
                    MOIST_WORKAROUNDS.CNV_Tracers().fscav[:],
                )
                safe_assign_array(
                    self._managed_convection_tracers.ndsl_state.vect_hcts.data[:],
                    MOIST_WORKAROUNDS.CNV_Tracers().Vect_Hcts[:],
                )
                safe_assign_array(
                    self._managed_convection_tracers.ndsl_state.use_gcc_washout.data[:],
                    MOIST_WORKAROUNDS.CNV_Tracers().use_gcc_washout[:],
                )

            with TimedCUDAProfiler("GF 2020 Convection Numerics", {}):
                self._managed_state.record("GF2020-In_python")
                # adjust pbl_level from fortran indexing to python indexing
                print("SHIFTING PBL LEVEL F TO P")
                self._managed_state.ndsl_state.pbl_level.field[:] = self._managed_state.ndsl_state.pbl_level.field[:] - 1

                # run GF 2020 Convection
                self._gf_2020(
                    state=self._managed_state.ndsl_state,
                    convection_tracers=self._managed_convection_tracers.ndsl_state,
                )
                # adjust pbl_level from python indexing to fortran indexing
                print("SHIFTING PBL LEVEL P TO F")
                self._managed_state.ndsl_state.pbl_level.field[:] = self._managed_state.ndsl_state.pbl_level.field[:] + 1

                self._managed_state.record("GF2020-Out_python")

            with TimedCUDAProfiler("GF 2020 Convection - State copy-back", {}):
                self._managed_state.ndsl_to_fortran()

    def finalize(
        self,
        mapl_state: CVoidPointer,
        import_state: CVoidPointer,
        export_state: CVoidPointer,
    ):
        self._managed_state.save_recorded()


CODE = GF2020Interface()
