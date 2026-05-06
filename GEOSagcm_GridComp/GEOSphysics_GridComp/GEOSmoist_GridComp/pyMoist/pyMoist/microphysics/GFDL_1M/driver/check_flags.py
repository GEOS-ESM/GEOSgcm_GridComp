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
        ValueError: list of non-compliant constants
    """
    failed_keywords = []
    if not GFDL_1M_config.LPHYS_HYDROSTATIC:
        failed_keywords.append("PHYS_HYDROSTATIC")
    if GFDL_1M_config.LHYDROSTATIC:
        failed_keywords.append("HYDROSTATIC")
    if GFDL_1M_config.CONST_VI:
        failed_keywords.append("CONST_VI")
    if GFDL_1M_config.CONST_VS:
        failed_keywords.append("CONST_VS")
    if GFDL_1M_config.CONST_VG:
        failed_keywords.append("CONST_VG")
    if GFDL_1M_config.CONST_VR:
        failed_keywords.append("CONST_VR")
    if GFDL_1M_config.USE_PPM:
        failed_keywords.append("USE_PPM")
    if not GFDL_1M_config.USE_CCN:
        failed_keywords.append("USE_CCN")
    if GFDL_1M_config.DO_QA:
        failed_keywords.append("DO_QA")
    if not GFDL_1M_config.FIX_NEGATIVE:
        failed_keywords.append("FIX_NEGATIVE")
    if GFDL_1M_config.FAST_SAT_ADJ:
        failed_keywords.append("FAST_SAT_ADJ")
    if GFDL_1M_config.DO_BIGG:
        failed_keywords.append("DO_BIGG")
    if GFDL_1M_config.DO_EVAP:
        failed_keywords.append("DO_EVAP")
    if GFDL_1M_config.DO_SUBL:
        failed_keywords.append("DO_SUBL")
    if not GFDL_1M_config.Z_SLOPE_LIQ:
        failed_keywords.append("Z_SLOPE_LIQ")
    if not GFDL_1M_config.Z_SLOPE_ICE:
        failed_keywords.append("Z_SLOPE_ICE")
    if not GFDL_1M_config.PROG_CCN:
        failed_keywords.append("PROG_CCN")
    if not GFDL_1M_config.PRECIPRAD:
        failed_keywords.append("PRECIPRAD")
    if not GFDL_1M_config.MONO_PROF:
        failed_keywords.append("MONO_PROF")
    if GFDL_1M_config.DO_SEDI_HEAT:
        failed_keywords.append("DO_SEDI_HEAT")
    if not GFDL_1M_config.SEDI_TRANSPORT:
        failed_keywords.append("SEDI_TRANSPORT")
    if GFDL_1M_config.DO_SEDI_W:
        failed_keywords.append("DO_SEDI_W")
    if GFDL_1M_config.DE_ICE:
        failed_keywords.append("DE_ICE")
    if GFDL_1M_config.MP_PRINT:
        failed_keywords.append("MP_PRINT")
    if config_dependent_constants.DTS >= 300.0:
        failed_keywords.append("DTS")

    if len(failed_keywords) > 0:
        raise ValueError(
            "One or more namelist parameters do not meet \
                expected values. Failing parameters: ",
            failed_keywords,
        )
