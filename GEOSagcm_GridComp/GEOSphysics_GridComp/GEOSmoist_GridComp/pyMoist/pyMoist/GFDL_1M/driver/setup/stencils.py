import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import FORWARD, PARALLEL, computation, interval, sqrt

from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.constants import constants
from pyMoist.shared_generic_math import sigma


def init_temporaries(
    t: FloatField,
    dp: FloatField,
    rhcrit3d: FloatField,
    qv: FloatField,
    ql: FloatField,
    qi: FloatField,
    qr: FloatField,
    qs: FloatField,
    qg: FloatField,
    qa: FloatField,
    qn: FloatField,
    qv0: FloatField,
    ql0: FloatField,
    qr0: FloatField,
    qi0: FloatField,
    qs0: FloatField,
    qg0: FloatField,
    qa0: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qa1: FloatField,
    dz: FloatField,
    uin: FloatField,
    vin: FloatField,
    w: FloatField,
    area: FloatFieldIJ,
    t1: FloatField,
    dp1: FloatField,
    omq: FloatField,
    den: FloatField,
    p_dry: FloatField,
    m1: FloatField,
    u1: FloatField,
    v1: FloatField,
    w1: FloatField,
    onemsig: FloatFieldIJ,
    ccn: FloatField,
    c_praut: FloatField,
    rh_limited: FloatField,
    rain: FloatFieldIJ,
    snow: FloatFieldIJ,
    graupel: FloatFieldIJ,
    ice: FloatFieldIJ,
    m2_rain: FloatField,
    m2_sol: FloatField,
    revap: FloatField,
    isubl: FloatField,
):
    """
    Initalize temporary copies of many quantities

    modification to quantities (t, p, q, etc.) made inside of the driver
    are not returned outside of the driver. these copies are necessary
    to ensure that no changes make it to the rest of the model

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    from __externals__ import cpaut

    with computation(PARALLEL), interval(...):
        t1 = t
        dp1 = dp  # moist air mass * grav

        # -----------------------------------------------------------------------
        # import horizontal subgrid variability with pressure dependence
        # total water subgrid deviation in horizontal direction
        # default area dependent form: use dx ~ 100 km as the base
        # -----------------------------------------------------------------------
        rh_limited = min(0.30, 1.0 - rhcrit3d)  # restricted to 70%

        # -----------------------------------------------------------------------
        # convert moist mixing ratios to dry mixing ratios
        # -----------------------------------------------------------------------

        dp1 = dp1 * (1.0 - qv)  # gfs
        omq = dp / dp1

        qv0 = qv * omq
        ql0 = ql * omq
        qr0 = qr * omq
        qi0 = qi * omq
        qs0 = qs * omq
        qg0 = qg * omq

        qv1 = qv0
        ql1 = ql0
        qr1 = qr0
        qi1 = qi0
        qs1 = qs0
        qg1 = qg0

        qa0 = qa
        qa1 = qa

        den = -dp1 / (constants.GRAV * dz)  # density of dry air
        p_dry = den * constants.RDGAS * t  # dry air pressure

        # -----------------------------------------------------------------------
        # for sedi_momentum
        # -----------------------------------------------------------------------

        m1 = 0.0
        u1 = uin
        v1 = vin
        w1 = w

        # ccn needs units #/m^3
        ccn = qn
        c_praut = cpaut * (ccn * constants.RHOR) ** (-1.0 / 3.0)

        # Reset precipitation aggregates to zero
        m2_rain = 0
        m2_sol = 0
        revap = 0
        isubl = 0

    with computation(FORWARD), interval(0, 1):
        # 1 minus sigma used to control minimum cloud
        # fraction needed to autoconvert ql->qr
        onemsig = 1.0 - sigma(sqrt(area))

        # Reset precipitation aggregates to zero
        rain = 0
        snow = 0
        graupel = 0
        ice = 0


@gtscript.function
def fix_negative_core(
    t: Float,
    qv: Float,
    ql: Float,
    qr: Float,
    qi: Float,
    qs: Float,
    qg: Float,
    c_air: Float,
    c_vap: Float,
    lv00: Float,
    d0_vap: Float,
):
    """
    Adjusts/removes negative mixing ratios

    reference Fortran: gfdl_cloud_microphys.F90: subroutine neg_adj
    """
    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    cvm = c_air + qv * c_vap + (qr + ql) * constants.C_LIQ + (qi + qs + qg) * constants.C_ICE
    lcpk = (lv00 + d0_vap * t) / cvm
    icpk = (constants.LI00 + constants.DC_ICE * t) / cvm

    # -----------------------------------------------------------------------
    # ice phase:
    # -----------------------------------------------------------------------

    # if cloud ice < 0, borrow from snow
    if qi < 0.0:
        qs = qs + qi
        qi = 0.0
    # if snow < 0, borrow from graupel
    if qs < 0.0:
        qg = qg + qs
        qs = 0.0
    # if graupel < 0, borrow from rain
    if qg < 0.0:
        qr = qr + qg
        t = t - qg * icpk  # heating
        qg = 0.0

    # -----------------------------------------------------------------------
    # liquid phase:
    # -----------------------------------------------------------------------

    # if rain < 0, borrow from cloud water
    if qr < 0.0:
        ql = ql + qr
        qr = 0.0
    # if cloud water < 0, borrow from water vapor
    if ql < 0.0:
        qv = qv + ql
        t = t - ql * lcpk  # heating
        ql = 0.0

    return t, qv, ql, qr, qi, qs, qg


def fix_negative_values(
    t: FloatField,
    qv: FloatField,
    ql: FloatField,
    qr: FloatField,
    qi: FloatField,
    qs: FloatField,
    qg: FloatField,
    dp: FloatField,
):
    """
    Stencil wrapper for fix_negative_core

    adjusts/removes negative mixing ratios
    updates qv based on new values

    refernce Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    from __externals__ import c_air, c_vap, d0_vap, lv00

    # -----------------------------------------------------------------------
    # fix all negative water species
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(0, -1):
        t, qv, ql, qr, qi, qs, qg = fix_negative_core(t, qv, ql, qr, qi, qs, qg, c_air, c_vap, lv00, d0_vap)
        if qv < 0.0:
            qv[0, 0, 1] = qv[0, 0, 1] + qv * dp / dp[0, 0, 1]
            qv = 0.0

    with computation(FORWARD), interval(-1, None):
        t, qv, ql, qr, qi, qs, qg = fix_negative_core(t, qv, ql, qr, qi, qs, qg, c_air, c_vap, lv00, d0_vap)

        if qv < 0.0 and qv[0, 0, -1] > 0.0:
            dq = min(-qv * dp, qv[0, 0, -1] * dp[0, 0, -1])
            qv[0, 0, -1] = qv[0, 0, -1] - dq / dp[0, 0, -1]
            qv = qv + dq / dp
