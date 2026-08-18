from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants


def check_flags(
    GFDL_1M_config: GFDL1MConfig,
    config_dependent_constants: GFDL1MDriverConfigDependentConstants,
):
    """Checks for any flags that are no meeting the expected value.
    Failing flags are likely not implemented,
    or at the very least not fully implemented

    Args:
        GFDL_1M_config (GFDL1MConfig): dataclass of all constants needed for GFDL Single Moment microphysics,
            gathered either from the namelist or another piece of the model (outside of this module)
        config_dependent_constants (GFDL1MDriverConfigDependentConstants): config dependent constants computed
            within the module

    Raises:
        ValueError: list of non-compliant constants with actual and expected values
    """
    failures = []

    # Helper function to evaluate and record non-compliant flags
    def check_param(param_name, actual, expected):
        if actual != expected:
            failures.append({"parameter": param_name, "actual": actual, "expected": expected})

    # Flag checks (Boolean flags)
    check_param("PHYS_HYDROSTATIC", GFDL_1M_config.LPHYS_HYDROSTATIC, True)
    check_param("HYDROSTATIC", GFDL_1M_config.LHYDROSTATIC, False)
    check_param("CONST_VI", GFDL_1M_config.CONST_VI, False)
    check_param("CONST_VS", GFDL_1M_config.CONST_VS, False)
    check_param("CONST_VG", GFDL_1M_config.CONST_VG, False)
    check_param("CONST_VR", GFDL_1M_config.CONST_VR, False)
    check_param("USE_PPM", GFDL_1M_config.USE_PPM, False)
    check_param("USE_CCN", GFDL_1M_config.USE_CCN, True)
    check_param("DO_QA", GFDL_1M_config.DO_QA, False)
    check_param("FIX_NEGATIVE", GFDL_1M_config.FIX_NEGATIVE, True)
    check_param("FAST_SAT_ADJ", GFDL_1M_config.FAST_SAT_ADJ, False)
    check_param("DO_BIGG", GFDL_1M_config.DO_BIGG, False)
    check_param("DO_EVAP", GFDL_1M_config.DO_EVAP, False)
    check_param("DO_SUBL", GFDL_1M_config.DO_SUBL, False)
    check_param("Z_SLOPE_LIQ", GFDL_1M_config.Z_SLOPE_LIQ, True)
    check_param("Z_SLOPE_ICE", GFDL_1M_config.Z_SLOPE_ICE, True)
    check_param("PROG_CCN", GFDL_1M_config.PROG_CCN, True)
    check_param("PRECIPRAD", GFDL_1M_config.PRECIPRAD, True)
    check_param("MONO_PROF", GFDL_1M_config.MONO_PROF, True)
    check_param("DO_SEDI_HEAT", GFDL_1M_config.DO_SEDI_HEAT, False)
    check_param("SEDI_TRANSPORT", GFDL_1M_config.SEDI_TRANSPORT, True)
    check_param("DO_SEDI_W", GFDL_1M_config.DO_SEDI_W, False)
    check_param("DE_ICE", GFDL_1M_config.DE_ICE, False)
    check_param("MP_PRINT", GFDL_1M_config.MP_PRINT, False)

    # Threshold checks
    if config_dependent_constants.DTS >= 300.0:
        failures.append({"parameter": "DTS", "actual": config_dependent_constants.DTS, "expected": "< 300.0"})

    if failures:
        formatted_failures = "\n".join(f"  - Parameter: {f['parameter']}, Actual: {f['actual']}, Expected: {f['expected']}" for f in failures)
        raise ValueError(f"One or more namelist parameters do not meet expected values:\n{formatted_failures}")
