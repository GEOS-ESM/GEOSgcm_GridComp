from ndsl.dsl.typing import Float
from pyMoist.GFDL_1M.driver.config import config


def check_flags(
    GFDL_1M_config: config,
    dts: Float,
):
    failed_keywords = []
    if not GFDL_1M_config.PHYS_HYDROSTATIC:
        failed_keywords.append("phys_hydrostatic")
    if GFDL_1M_config.HYDROSTATIC:
        failed_keywords.append("hydrostatic")
    if GFDL_1M_config.CONST_VI:
        failed_keywords.append("const_vi")
    if GFDL_1M_config.CONST_VS:
        failed_keywords.append("const_vs")
    if GFDL_1M_config.CONST_VG:
        failed_keywords.append("const_vg")
    if GFDL_1M_config.CONST_VR:
        failed_keywords.append("const_vr")
    if GFDL_1M_config.USE_PPM:
        failed_keywords.append("use_ppm")
    if not GFDL_1M_config.USE_CCN:
        failed_keywords.append("use_ccn")
    if GFDL_1M_config.DO_QA:
        failed_keywords.append("do_qa")
    if not GFDL_1M_config.FIX_NEGATIVE:
        failed_keywords.append("fix_negative")
    if GFDL_1M_config.FAST_SAT_ADJ:
        failed_keywords.append("fast_sat_adj")
    if GFDL_1M_config.DO_BIGG:
        failed_keywords.append("do_bigg")
    if GFDL_1M_config.DO_EVAP:
        failed_keywords.append("do_evap")
    if GFDL_1M_config.DO_SUBL:
        failed_keywords.append("do_subl")
    if not GFDL_1M_config.Z_SLOPE_LIQ:
        failed_keywords.append("z_slope_liq")
    if not GFDL_1M_config.Z_SLOPE_ICE:
        failed_keywords.append("z_slope_ice")
    if not GFDL_1M_config.PROG_CCN:
        failed_keywords.append("prog_ccn")
    if not GFDL_1M_config.PRECIPRAD:
        failed_keywords.append("preciprad")
    if not GFDL_1M_config.MONO_PROF:
        failed_keywords.append("mono_prof")
    if GFDL_1M_config.DO_SEDI_HEAT:
        failed_keywords.append("do_sedi_heat")
    if not GFDL_1M_config.SEDI_TRANSPORT:
        failed_keywords.append("sedi_transport")
    if GFDL_1M_config.DO_SEDI_W:
        failed_keywords.append("do_sedi_w")
    if GFDL_1M_config.DE_ICE:
        failed_keywords.append("de_ice")
    if GFDL_1M_config.MP_PRINT:
        failed_keywords.append("mp_print")
    if dts >= 300:
        failed_keywords.append("dts")

    if len(failed_keywords) > 0:
        raise ValueError(
            "One or more namelist parameters do not meet \
                expected values. Failing parameters: ",
            failed_keywords,
        )
