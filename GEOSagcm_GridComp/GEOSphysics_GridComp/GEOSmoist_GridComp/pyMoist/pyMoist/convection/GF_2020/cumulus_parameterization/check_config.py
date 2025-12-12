from ndsl import StencilFactory, QuantityFactory, ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)


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
            " GF2020 cumulus parameterization: environmental_conditions construced with"
            "untested SATURATION_CALCULATION_CHOICE option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.CLOUD_LEVEL_GRID != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: environmental_cloud_levels construced with "
            "untested CLOUD_LEVEL_GRID option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.FRAC_MODIS != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: partition_liquid_ice constructed with "
            "untested FRAC_MODIS option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.MELT_GLAC != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: partition_liquid_ice constructed with "
            "untested MELT_GLAC option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: find_lcl constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if config.ADV_TRIGGER != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization find_lcl constructed with "
            "untested ADV_TRIGGER option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: parcel_moist_static_energy constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: entrainment_rates constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.OVERSHOOT != 0:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: get_convective_cloud_base_level constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 0:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: get_convective_cloud_base_level constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.MOIST_TRIGGER != 0:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: get_convective_cloud_base_level constructed with "
            "untested MOIST_TRIGGER option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.USE_MEMORY != -1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: get_convective_cloud_base_level constructed with "
            "untested USE_MEMORY option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: get_convective_cloud_base_level constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.OVERSHOOT != 0:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: updraft_rates_pdf constructed with "
            "untested OVERSHOOT option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.BOUNDARY_CONDITION_METHOD != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: compute_uc_vc constructed with "
            "untested BOUNDARY_CONDITION_METHOD option. Running untested code... proceed with caution"
        )

    if cumulus_parameterization_config.ZERO_DIFF != 1:
        ndsl_log.warning(
            " GF2020 cumulus parameterization: updraft_vertical_velocity constructed with "
            "untested ZERO_DIFF option. Running untested code... proceed with caution"
        )
