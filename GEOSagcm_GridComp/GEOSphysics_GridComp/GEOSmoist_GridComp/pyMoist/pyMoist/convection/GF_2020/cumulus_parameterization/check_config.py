from ndsl import ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import MELT_GLAC, AEROEVAP, COUPLE_MICROPHYSICS, WRTGRADS


def check_config(
    config: GF2020Config,
    cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
):
    """
    Check config and cumulus_parameterization_config for unexpected options.

    Warn about untested options, block execution of unimplemented options.
    """

    if cumulus_parameterization_config.SATURATION_CALCULATION_CHOICE != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->environmental_conditions construced with"
            "untested SATURATION_CALCULATION_CHOICE option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.CLOUD_LEVEL_GRID != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->environmental_cloud_levels construced with "
            "untested CLOUD_LEVEL_GRID option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.FRAC_MODIS != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->partition_liquid_ice constructed with "
            "untested FRAC_MODIS option. Running untested code... proceed with caution"
        )

    if not MELT_GLAC:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->partition_liquid_ice constructed with "
            "untested MELT_GLAC option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->find_lcl constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if config.ADV_TRIGGER != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->find_lcl constructed with "
            "untested ADV_TRIGGER option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->parcel_moist_static_energy constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->entrainment_rates constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.OVERSHOOT != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->get_convective_cloud_base_level constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->get_convective_cloud_base_level constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.MOIST_TRIGGER != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->get_convective_cloud_base_level constructed with "
            "untested MOIST_TRIGGER option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.USE_MEMORY != -1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->get_convective_cloud_base_level constructed with "
            "untested USE_MEMORY option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->get_convective_cloud_base_level constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.OVERSHOOT != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->updraft_rates_pdf constructed with "
            "untested OVERSHOOT option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->compute_uc_vc constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->updraft_moisture constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.FRAC_MODIS != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->updraft_moisture constructed with "
            "untested FRAC_MODIS option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->updraft_moisture constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->updraft_vertical_velocity constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.USE_WETBULB != 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->downdraft_moisture constructed with "
            "untested USE_WETBULB option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.SGS_W_TIMESCALE != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->DiurnalCycle initialized with "
            "untested SGS_W_TIMESCALE option. Running untested code... proceed with caution"
        )

    if AEROEVAP != 1:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->DowndraftWindshear initialized with "
            "untested AEROEVAP option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.USE_TRACER_EVAPORATION == 0:
        ndsl_log.warning(
            " GF2020-->CumulusParameterization-->AtmosphericComposition-->downdraft_chemistry initialized "
            "with untested USE_TRACER_EVAPORATION option. Running untested code... proceed with caution"
        )

    # generate errors after all warnings
    # TODO find a way to generate all then fail at once, printing all simultaneously
    if cumulus_parameterization_config.OUTPUT_SOUNDING == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with OUTPUT_SOUNDING = 1. This requires"
            "an unimplemented output sounding capabilities (in the Sounding class). Please implement, then"
            "disable this error manually to proceed."
        )

    if config.CONVECTION_TRACER == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with CONVECTION_TRACER = 1. This requires"
            "an unimplemented class ColdPoolParameterization. Please implement, then disable this error"
            "manually to proceed."
        )

    if (
        cumulus_parameterization_config.USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES == 1
        and cumulus_parameterization_config.ENABLE_SHALLOW == 1
    ):
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with"
            "USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES == 1 and shallow plume enabled. This combination requires"
            "a call to the unimplemented function get_delmix in CumulusParameterization and updraft_moisture."
            "Please implement, then disable this error manually to proceed."
        )

    if cumulus_parameterization_config.ZERO_DIFF == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with ZERO_DIFF = 1. This combination"
            "requires an unimplemented porion of UpdraftMassFlux. Please implement, then disable this error"
            "manually to proceed."
        )

    if cumulus_parameterization_config.ENABLE_SHALLOW == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with shallow plume enabled. This requires"
            "an unimplemented porion of UpdraftMassFlux. Please implement, then disable this error"
            "manually to proceed."
        )

    if config.AUTOCONV != 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with AUTOCONV != 1. This requires"
            "an unimplemented option in updraft_moisture. Please implement, then disable this"
            "error manually to proceed."
        )

    if cumulus_parameterization_config.USE_WETBULB == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with USE_WETBULB = 1. This setting requires"
            "the unimplemented function get_wetbulb and an unimplemented option in"
            "downdraft_moist_static_energy_and_moisture_budget. Please implement, then disable this error"
            "manually to proceed."
        )

    if cumulus_parameterization_config.DIURNAL_CYCLE != 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with DIURNAL_CYCLE != 1. This setting"
            "requires an unimplemented option for multiple stencils in DiurnalCycle. Please implement, then"
            "disable this error manually to proceed."
        )

    if config.ADV_TRIGGER == 3:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with ADV_TRIGGER = 3. This setting requires"
            "the unimplemented XieTriggerFunction. Please implement, then disable this error"
            "manually to proceed."
        )

    if cumulus_parameterization_config.VERTICAL_DISCRETIZATION_OPTION != 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with VERTICAL_DISCRETIZATION_OPTION != 1."
            "This setting requires an unimplemented option in VerticalDiscretization. Please implement,"
            "then disable this error manually to proceed."
        )

    if cumulus_parameterization_config.USE_FCT != 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with USE_FCT != 0. This setting requires"
            "an unimplemented option for multiple stencils in VerticalDiscretization. Please implement,"
            "then disable this error manually to proceed."
        )

    if config.APPLY_SUBSIDENCE_MICROPHYSICS != 0:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with APPLY_SUBSIDENCE_MICROPHYSICS != 0."
            "This setting requires an unimeplemented option in EnvironmentalSubsidence. Please implement,"
            "then disable this error manually to proceed."
        )

    if cumulus_parameterization_config.DIURNAL_CYCLE not in (1, 6):
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with DIURNAL_CYCLE != 1 or 6. This"
            "setting requires one or more unimeplemented options in LargeScaleForcing. Please"
            "implement, then disable this error manually to proceed."
        )

    if cumulus_parameterization_config.ENABLE_SHALLOW == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with shallow plume enabled. This requires"
            "an unimeplemented functions in LargeScaleForcing. Please implement, then disable"
            "this error manually to proceed."
        )

    if cumulus_parameterization_config.ENABLE_SHALLOW == 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with shallow plume enabled. This requires"
            "an umplemented portion of ensemble_output_and_feedback. Please impelment, then disable this"
            "error manually to proceed."
        )

    if not COUPLE_MICROPHYSICS:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with COUPLE_MICROPHYSICS != True. This"
            "setting requires an unimeplemented option in cloud_dissipation. Please implement, then"
            "disable this error manually to proceed."
        )

    if not cumulus_parameterization_config.LIGHTNING_DIAGNOSTICS:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with LIGHTNING_DIAGNOSTICS != True. This"
            "setting requires unimeplemented functions. Please implement, then disable this error manually"
            "to proceed."
        )

    if not cumulus_parameterization_config.USE_TRACER_SCAVENGE:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition initialized with"
            "USE_TRACER_SCAVENGE != 1. This setting requires unimeplemented options in updraft_chemistry."
            "Please implement, then disable this error manually to proceed."
        )

    if cumulus_parameterization_config.USE_FLUX_FORM != 1:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition initialized with"
            "USE_FLUX_FORM != 1. This setting requires one or more unimeplemented options in"
            "vertical_transport_part_1. Please implement, then disable this error manually to proceed."
        )

    if cumulus_parameterization_config.ALP1 == 0:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization-->AtmosphericComposition initialized with"
            "ALP1 == 0. This setting requires an unimeplemented option in vertical_transport_part_1."
            "Please implement, then disable this error manually to proceed."
        )

    if not WRTGRADS:
        raise NotImplementedError(
            "[NDSL] GF2020-->CumulusParameterization initialized with WRTGRADS = True. This setting requires"
            "the unimeplemented GATE sounding capabilities (in the GATESounding class). Please implement,"
            "then disable this error manually to proceed."
        )
