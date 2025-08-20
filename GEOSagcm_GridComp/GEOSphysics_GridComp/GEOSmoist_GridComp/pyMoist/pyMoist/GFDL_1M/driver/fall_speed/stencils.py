import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import PARALLEL, computation, exp, interval, log, log10, sqrt

from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.constants import constants


@gtscript.function
def get_top_speed(
    p_dry: Float,
    cnv_frc: Float,
    anv_icefall: Float,
    lsc_icefall: Float,
    den: Float,
    qs: Float,
    qi: Float,
    qg: Float,
    ql: Float,
    t: Float,
):
    """
    Calculate the vertical fall speed of precipitation

    reference Fortran: gfdl_cloud_microphys.F90: subroutine fall_speed
    """
    from __externals__ import const_vg, const_vi, const_vs, vg_fac, vg_max, vi_fac, vi_max, vs_fac, vs_max

    rhof = sqrt(min(10.0, constants.SFCRHO / den))
    if const_vi == True:  # noqa
        vti = vi_fac
    else:
        if qi < constants.THI:
            vti = constants.VF_MIN
        else:
            # -----------------------------------------------------------------------
            # ice:
            # -----------------------------------------------------------------------

            vi1 = 0.01 * vi_fac
            tc = t - constants.TICE  # deg C
            IWC = qi * den * 1.0e3  # Units are g/m3
            # -----------------------------------------------------------------------
            # use deng and mace (2008, grl)
            # https://doi.org/10.1029/2008GL035054
            # -----------------------------------------------------------------------
            viLSC = lsc_icefall * 10.0 ** (
                log10(IWC) * (tc * (constants.AAL * tc + constants.BBL) + constants.CCL)
                + constants.DDL * tc
                + constants.EEL
            )
            viCNV = anv_icefall * 10.0 ** (
                log10(IWC) * (tc * (constants.AAC * tc + constants.BBC) + constants.CCC)
                + constants.DDC * tc
                + constants.EEC
            )
            # Combine
            vti = viLSC * (1.0 - cnv_frc) + viCNV * (cnv_frc)
            # Update units from cm/s to m/s
            vti = vi1 * vti
            # Limits
            vti = min(vi_max, max(constants.VF_MIN, vti))

    # -----------------------------------------------------------------------
    # snow:
    # -----------------------------------------------------------------------

    if const_vs == True:  # noqa
        vts = vs_fac  # 1. ifs_2016
    else:
        if qs < constants.THS:
            vts = constants.VF_MIN
        else:
            vts = vs_fac * constants.VCONS * rhof * exp(0.0625 * log(qs * den / constants.NORMS))
            vts = min(vs_max, max(constants.VF_MIN, vts))

    # -----------------------------------------------------------------------
    # graupel:
    # -----------------------------------------------------------------------

    if const_vg == True:  # noqa
        vtg = vg_fac  # 2.
    else:
        if qg < constants.THG:
            vtg = constants.VF_MIN
        else:
            vtg = vg_fac * constants.VCONG * rhof * sqrt(sqrt(sqrt(qg * den / constants.NORMG)))
            vtg = min(vg_max, max(constants.VF_MIN, vtg))

    return vti, vts, vtg


def fall_speed_core(
    ql1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    t: FloatField,
    t1: FloatField,
    dz: FloatField,
    dz1: FloatField,
    den: FloatField,
    den1: FloatField,
    denfac: FloatField,
    p_dry: FloatField,
    vti: FloatField,
    vts: FloatField,
    vtg: FloatField,
    cnv_frc: FloatFieldIJ,
    anv_icefall: Float,
    lsc_icefall: Float,
):
    """
    Stencil wrapper for fall_speed_core

    calculate the vertical fall speed of precipitation

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    from __externals__ import p_nonhydro

    with computation(PARALLEL), interval(...):
        if p_nonhydro:
            dz1 = dz
            den1 = den  # dry air density remains the same
            denfac = sqrt(constants.SFCRHO / den1)
        else:
            dz1 = dz * t1 / t  # hydrostatic balance
            den1 = den * dz / dz1
            denfac = sqrt(constants.SFCRHO / den1)

        vti, vts, vtg = get_top_speed(p_dry, cnv_frc, anv_icefall, lsc_icefall, den1, qs1, qi1, qg1, ql1, t1)
