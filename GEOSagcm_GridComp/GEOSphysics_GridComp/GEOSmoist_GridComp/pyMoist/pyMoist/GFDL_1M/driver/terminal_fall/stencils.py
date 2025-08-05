from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, computation, exp, interval, function
from ndsl.dsl.typing import BoolField, Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.constants import constants


def check_precip_get_zt(
    q: FloatField,
    vt: FloatField,
    ze: FloatField,
    zt: FloatField,
    precip: FloatFieldIJ,
    precip_fall: FloatFieldIJ,
):
    from __externals__ import dts

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # melting of falling snow into rain
    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    # reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column
    # determine if any precip falls in the column
    # if it falls anywhere in the column, the entire column becomes true
    # initalized to 0 (false), potentially changed to 1 (true)
    with computation(FORWARD), interval(...):
        if q > constants.QPMIN:
            precip_fall = 1
    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine check_column

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 0:
            precip = 0

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            zt = ze - dts * (vt[0, 0, -1] + vt) / 2.0

    with computation(FORWARD), interval(-1, None):
        if precip_fall == 1:
            zt[0, 0, 1] = constants.ZS - dts * vt

    with computation(FORWARD), interval(...):
        if precip_fall == 1:
            if zt[0, 0, 1] >= zt:
                zt[0, 0, 1] = zt - constants.DZ_MIN


def update_dm(
    dm: FloatField,
    dp1: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    precip_fall: FloatFieldIJ,
):
    from __externals__ import do_sedi_w

    with computation(PARALLEL), interval(...):
        if precip_fall == 1:
            if do_sedi_w == True:  # noqa
                dm = dp1 * (1.0 + qv1 + ql1 + qr1 + qi1 + qs1 + qg1)


def update_w1(
    w1: FloatField,
    dm: FloatField,
    m1: FloatField,
    vt: FloatField,
    precip_fall: FloatFieldIJ,
):
    from __externals__ import do_sedi_w

    with computation(FORWARD), interval(0, 1):
        if precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 + m1 * vt) / (dm - m1)

    with computation(FORWARD), interval(1, None):
        if precip_fall == 1:
            if do_sedi_w:
                w1 = (dm * w1 - m1[0, 0, -1] * vt[0, 0, -1] + m1 * vt) / (dm + m1[0, 0, -1] - m1)


def reset(
    m1: FloatField,
    precip_fall: FloatFieldIJ,
):
    # reset masks and temporaries
    with computation(FORWARD), interval(0, 1):
        precip_fall = 0

    with computation(PARALLEL), interval(...):
        # must be reset to zero everywhere in case there is no precipitate
        # in a column and no calculations are performed
        m1 = 0.0


@function
def prefall_melting(
    t: Float,
    qv: Float,
    ql: Float,
    qr: Float,
    qg: Float,
    qs: Float,
    qi: Float,
    m1_sol: Float,
    is_frozen: bool,
    c_air: Float,
    c_vap: Float,
    d0_vap: Float,
    dts: Float,
    lv00: Float,
    ql_mlt: Float,
    tau_imlt: Float,
):
    """
    Melt cloud ice before fall

    reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall
    """

    fac_imlt = 1.0 - exp(-dts / tau_imlt)

    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    m1_sol = 0.0
    lhl = lv00 + d0_vap * t
    lhi = constants.LI00 + constants.DC_ICE * t
    q_liq = ql + qr
    q_sol = qi + qs + qg
    cvm = c_air + qv * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
    lcpk = lhl / cvm
    icpk = lhi / cvm

    # -----------------------------------------------------------------------
    # melting of cloud_ice (before fall) :
    # -----------------------------------------------------------------------

    tc = t - constants.TICE
    if is_frozen == False and (qi > constants.QCMIN and tc > 0.0):  # noqa
        sink = min(qi, fac_imlt * tc / icpk)
        if ql_mlt - ql > 0:
            ans = ql_mlt - ql
        else:
            ans = 0
        tmp = min(sink, ans)
        ql = ql + tmp
        qr = qr + sink - tmp
        qi = qi - sink
        q_liq = q_liq + sink
        q_sol = q_sol - sink
        cvm = c_air + qv * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t - sink * lhi / cvm

    return t, ql, qr, qi, cvm, lhi, icpk, m1_sol


def setup(
    t1: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qg1: FloatField,
    qs1: FloatField,
    qi1: FloatField,
    dz1: FloatField,
    m1_sol: FloatField,
    ze: FloatField,
    zt: FloatField,
    lhi: FloatField,
    icpk: FloatField,
    cvm: FloatField,
    is_frozen: BoolField,
    precip_rain: FloatFieldIJ,
    precip_graupel: FloatFieldIJ,
    precip_snow: FloatFieldIJ,
    precip_ice: FloatFieldIJ,
):
    """
    Calculate terminal fall speed, accounting for
    melting of ice, snow, and graupel during fall

    reference Fortran: gfdl_cloud_microphys.F90: subroutine terminal_fall
    """
    from __externals__ import c_air, c_vap, d0_vap, dts, lv00, ql_mlt, tau_imlt

    # determine frozen levels
    # later operations will only be executed if frozen/melted
    # initalized to is_frozen = False, True = frozen, False = melted
    with computation(PARALLEL), interval(...):
        if t1 <= constants.TICE:
            is_frozen = True

    # we only want the melting layer closest to the surface
    with computation(BACKWARD), interval(0, -1):
        if is_frozen[0, 0, 1] == True and is_frozen[0, 0, 0] == False:  # noqa
            is_frozen = True

    # force surface to "melt" for later calculations
    with computation(PARALLEL), interval(-1, None):
        is_frozen = False

    with computation(PARALLEL), interval(...):
        t1, ql1, qr1, qi1, cvm, lhi, icpk, m1_sol = prefall_melting(
            t1,
            qv1,
            ql1,
            qr1,
            qg1,
            qs1,
            qi1,
            m1_sol,
            is_frozen,
            c_air,
            c_vap,
            d0_vap,
            dts,
            lv00,
            ql_mlt,
            tau_imlt,
        )

    # if timestep is too small turn off melting (only calculate at surface)
    with computation(PARALLEL), interval(0, -1):
        if dts < 300.0:
            is_frozen = True

    with computation(PARALLEL), interval(-1, None):
        if dts < 300.0:
            is_frozen = False

    with computation(BACKWARD), interval(...):
        ze = ze[0, 0, 1] - dz1  # dz < 0

    with computation(FORWARD), interval(0, 1):
        zt = ze

    # -----------------------------------------------------------------------
    # update capacity heat and latend heat coefficient
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        if is_frozen == False:  # noqa
            lhi = constants.LI00 + constants.DC_ICE * t1
            icpk = lhi / cvm

    with computation(FORWARD), interval(-1, None):
        lhi = constants.LI00 + constants.DC_ICE * t1
        icpk = lhi / cvm

    # zero local precipitaton values
    with computation(FORWARD), interval(0, 1):
        precip_rain = 0
        precip_graupel = 0
        precip_snow = 0
        precip_ice = 0


def update_outputs(
    rain: FloatFieldIJ,
    graupel: FloatFieldIJ,
    snow: FloatFieldIJ,
    ice: FloatFieldIJ,
    precip_rain: FloatFieldIJ,
    precip_graupel: FloatFieldIJ,
    precip_snow: FloatFieldIJ,
    precip_ice: FloatFieldIJ,
):
    """
    Update precipitation totals with results of terminal_fall stencil

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    with computation(FORWARD), interval(0, 1):
        rain = rain + precip_rain  # from melted snow & ice that reached the ground
        graupel = graupel + precip_graupel
        snow = snow + precip_snow
        ice = ice + precip_ice
