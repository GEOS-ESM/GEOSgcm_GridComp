"""Core functions and stencils of the GFDL_1M driver"""

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    BACKWARD,
    computation,
    interval,
    sqrt,
    log,
    log10,
    exp,
    trunc,
    i32,
)
from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, Int, BoolField
import pyMoist.GFDL_1M.GFDL_1M_driver.GFDL_1M_driver_constants as driver_constants
from pyMoist.shared_generic_math import sigma
from pyMoist.shared_incloud_processes import ice_fraction

length = 2621
GlobalTable_driver_qsat = gtscript.GlobalTable[(Float, (length))]


def create_temporaries(
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
    dz: FloatField,
    uin: FloatField,
    vin: FloatField,
    w: FloatField,
    area: FloatFieldIJ,
    t1: FloatField,
    dp1: FloatField,
    omq: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qa1: FloatField,
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
):
    from __externals__ import cpaut

    with computation(PARALLEL), interval(...):
        # -----------------------------------------------------------------------
        # major cloud microphysics
        # mpdrv in Fortran, inlined to avoid another massive function call
        # -----------------------------------------------------------------------

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

        qv1 = qv
        ql1 = ql
        qi1 = qi
        qr1 = qr
        qs1 = qs
        qg1 = qg

        dp1 = dp1 * (1.0 - qv1)  # gfs
        omq = dp / dp1

        qv1 = qv1 * omq
        ql1 = ql1 * omq
        qr1 = qr1 * omq
        qi1 = qi1 * omq
        qs1 = qs1 * omq
        qg1 = qg1 * omq

        qa1 = qa

        den = -dp1 / (driver_constants.grav * dz)  # density of dry air
        p_dry = den * driver_constants.rdgas * t  # dry air pressure

        # -----------------------------------------------------------------------
        # for sedi_momentum
        # -----------------------------------------------------------------------

        m1 = 0.0
        u1 = uin
        v1 = vin
        w1 = w

        # # 1 minus sigma used to control minimum cloud fraction needed to autoconvert ql->qr
        # onemsig = 1.0 - sigma(sqrt(area))

        # ccn needs units #/m^3
        ccn = qn
        c_praut = cpaut * (ccn * driver_constants.rhor) ** (-1.0 / 3.0)

    with computation(FORWARD), interval(0, 1):
        # 1 minus sigma used to control minimum cloud fraction needed to autoconvert ql->qr
        onemsig = 1.0 - sigma(sqrt(area))


@gtscript.function
def fix_negative_values(
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
    Adjusts/removes negative mixing ratios.

    Reference Fortran: gfdl_cloud_microphys.F90: subroutine neg_adj
    """
    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    cvm = (
        c_air
        + qv * c_vap
        + (qr + ql) * driver_constants.c_liq
        + (qi + qs + qg) * driver_constants.c_ice
    )
    lcpk = (lv00 + d0_vap * t) / cvm
    icpk = (driver_constants.li00 + driver_constants.dc_ice * t) / cvm

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


def gfdl_1m_driver_preloop(
    t: FloatField,
    qv: FloatField,
    ql: FloatField,
    qr: FloatField,
    qi: FloatField,
    qs: FloatField,
    qg: FloatField,
    dp: FloatField,
):
    from __externals__ import (
        c_air,
        c_vap,
        d0_vap,
        lv00,
    )

    # -----------------------------------------------------------------------
    # fix all negative water species
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(0, -1):
        t, qv, ql, qr, qi, qs, qg = fix_negative_values(
            t, qv, ql, qr, qi, qs, qg, c_air, c_vap, lv00, d0_vap
        )
        if qv < 0.0:
            qv[0, 0, 1] = qv[0, 0, 1] + qv * dp / dp[0, 0, 1]
            qv = 0.0

    with computation(FORWARD), interval(-1, None):
        t, qv, ql, qr, qi, qs, qg = fix_negative_values(
            t, qv, ql, qr, qi, qs, qg, c_air, c_vap, lv00, d0_vap
        )

        if qv < 0.0 and qv[0, 0, -1] > 0.0:
            dq = min(-qv * dp, qv[0, 0, -1] * dp[0, 0, -1])
            qv[0, 0, -1] = qv[0, 0, -1] - dq / dp[0, 0, -1]
            qv = qv + dq / dp


@gtscript.function
def fall_speed(
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
    Calculated the vertical fall speed of precipitation.

    Reference Fortran: gfdl_cloud_microphys.F90: subroutine fall_speed
    """
    from __externals__ import (
        const_vi,
        const_vs,
        const_vg,
        vi_fac,
        vi_max,
        vs_fac,
        vs_max,
        vg_fac,
        vg_max,
    )

    rhof = sqrt(min(10.0, driver_constants.sfcrho / den))
    if const_vi == True:
        vti = vi_fac
    else:
        if qi < driver_constants.thi:
            vti = driver_constants.vf_min
        else:
            # -----------------------------------------------------------------------
            # ice:
            # -----------------------------------------------------------------------

            vi1 = 0.01 * vi_fac
            tc = t - driver_constants.tice  # deg C
            IWC = qi * den * 1.0e3  # Units are g/m3
            # -----------------------------------------------------------------------
            # use deng and mace (2008, grl)
            # https://doi.org/10.1029/2008GL035054
            # -----------------------------------------------------------------------
            viLSC = lsc_icefall * 10.0 ** (
                log10(IWC)
                * (
                    tc * (driver_constants.aaL * tc + driver_constants.bbL)
                    + driver_constants.ccL
                )
                + driver_constants.ddL * tc
                + driver_constants.eeL
            )
            viCNV = anv_icefall * 10.0 ** (
                log10(IWC)
                * (
                    tc * (driver_constants.aaC * tc + driver_constants.bbC)
                    + driver_constants.ccC
                )
                + driver_constants.ddC * tc
                + driver_constants.eeC
            )
            # Combine
            vti = viLSC * (1.0 - cnv_frc) + viCNV * (cnv_frc)
            # Update units from cm/s to m/s
            vti = vi1 * vti
            # Limits
            vti = min(vi_max, max(driver_constants.vf_min, vti))

    # -----------------------------------------------------------------------
    # snow:
    # -----------------------------------------------------------------------

    if const_vs == True:
        vts = vs_fac  # 1. ifs_2016
    else:
        if qs < driver_constants.ths:
            vts = driver_constants.vf_min
        else:
            vts = (
                vs_fac
                * driver_constants.vcons
                * rhof
                * exp(0.0625 * log(qs * den / driver_constants.norms))
            )
            vts = min(vs_max, max(driver_constants.vf_min, vts))

    # -----------------------------------------------------------------------
    # graupel:
    # -----------------------------------------------------------------------

    if const_vg == True:
        vtg = vg_fac  # 2.
    else:
        if qg < driver_constants.thg:
            vtg = driver_constants.vf_min
        else:
            vtg = (
                vg_fac
                * driver_constants.vcong
                * rhof
                * sqrt(sqrt(sqrt(qg * den / driver_constants.normg)))
            )
            vtg = min(vg_max, max(driver_constants.vf_min, vtg))

    return vti, vts, vtg


def gfdl_1m_driver_loop_1(
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
    from __externals__ import (
        p_nonhydro,
        const_vi,
        const_vs,
        const_vg,
        vi_fac,
        vi_max,
        vs_fac,
        vs_max,
        vg_fac,
        vg_max,
    )

    with computation(PARALLEL), interval(...):
        if p_nonhydro:
            dz1 = dz
            den1 = den  # dry air density remains the same
            denfac = sqrt(driver_constants.sfcrho / den1)
        else:
            dz1 = dz * t1 / t  # hydrostatic balance
            den1 = den * dz / dz1
            denfac = sqrt(driver_constants.sfcrho / den1)

        vti, vts, vtg = fall_speed(
            p_dry, cnv_frc, anv_icefall, lsc_icefall, den1, qs1, qi1, qg1, ql1, t1
        )


def gfdl_1m_driver_loop_2(
    rain: FloatFieldIJ,
    graupel: FloatFieldIJ,
    snow: FloatFieldIJ,
    ice: FloatFieldIJ,
    rain1: FloatFieldIJ,
    graupel1: FloatFieldIJ,
    snow1: FloatFieldIJ,
    ice1: FloatFieldIJ,
):
    with computation(FORWARD), interval(0, 1):
        rain = rain + rain1  # from melted snow & ice that reached the ground
        snow = snow + snow1
        graupel = graupel + graupel1
        ice = ice + ice1


def gfdl_1m_driver_loop_3(
    rain: FloatFieldIJ,
    rain1: FloatFieldIJ,
    evap1: FloatField,
    revap: FloatField,
    m1_rain: FloatField,
    m2_rain: FloatField,
    m1_sol: FloatField,
    m2_sol: FloatField,
    m1: FloatField,
):
    with computation(PARALLEL), interval(...):
        revap = revap + evap1
        m2_rain = m2_rain + m1_rain
        m2_sol = m2_sol + m1_sol
        m1 = m1 + m1_rain + m1_sol

    with computation(FORWARD), interval(0, 1):
        rain = rain + rain1


def gfdl_1m_driver_loop_4(
    isubl: FloatField,
    subl1: FloatField,
):
    with computation(PARALLEL), interval(...):
        isubl = isubl + subl1


def gfdl_1m_driver_postloop(
    qv: FloatField,
    ql: FloatField,
    qr: FloatField,
    qi: FloatField,
    qs: FloatField,
    qg: FloatField,
    qa: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qa1: FloatField,
    qn: FloatField,  # NACTL + NACTI
    qv_dt: FloatField,
    ql_dt: FloatField,
    qr_dt: FloatField,
    qi_dt: FloatField,
    qs_dt: FloatField,
    qg_dt: FloatField,
    qa_dt: FloatField,
    t: FloatField,
    t1: FloatField,
    t_dt: FloatField,
    w: FloatField,
    w1: FloatField,
    uin: FloatField,
    u1: FloatField,
    udt: FloatField,
    vin: FloatField,
    v1: FloatField,
    vdt: FloatField,
    dz: FloatField,
    dp: FloatField,
    dp1: FloatField,
    den: FloatField,
    p_dry: FloatField,
    area: FloatFieldIJ,
    dt_moist: Float,
    fr_land: FloatFieldIJ,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    eis: FloatFieldIJ,
    rh_limited: FloatField,
    m1: FloatField,
    anv_icefall: Float,
    lsc_icefall: Float,
    revap: FloatField,  # strict output
    isubl: FloatField,  # strict output
    rain: FloatFieldIJ,  # strict output
    snow: FloatFieldIJ,  # strict output
    ice: FloatFieldIJ,  # strict output
    graupel: FloatFieldIJ,  # strict output
    m2_rain: FloatField,  # strict output
    m2_sol: FloatField,  # strict output
):
    from __externals__ import (
        c_air,
        c_vap,
        rdt,
        do_sedi_w,
        sedi_transport,
        do_qa,
    )

    # -----------------------------------------------------------------------
    # momentum transportation during sedimentation
    # note: dp1 is dry mass; dp0 is the old moist (total) mass
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(1, None):
        if sedi_transport == True:
            u1 = (dp * u1 + m1[0, 0, -1] * u1[0, 0, -1]) / (dp + m1[0, 0, -1])
            v1 = (dp * v1 + m1[0, 0, -1] * v1[0, 0, -1]) / (dp + m1[0, 0, -1])
            udt = udt + (u1 - uin) * rdt
            vdt = vdt + (v1 - vin) * rdt

    with computation(PARALLEL), interval(...):
        if do_sedi_w:
            w = w1

    # -----------------------------------------------------------------------
    # update moist air mass (actually hydrostatic pressure)
    # convert to dry mixing ratios
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        omq = dp1 / dp
        qv_dt = qv_dt + rdt * (qv1 - qv) * omq
        ql_dt = ql_dt + rdt * (ql1 - ql) * omq
        qr_dt = qr_dt + rdt * (qr1 - qr) * omq
        qi_dt = qi_dt + rdt * (qi1 - qi) * omq
        qs_dt = qs_dt + rdt * (qs1 - qs) * omq
        qg_dt = qg_dt + rdt * (qg1 - qg) * omq
        cvm = (
            c_air
            + qv1 * c_vap
            + (qr1 + ql1) * driver_constants.c_liq
            + (qi1 + qs1 + qg1) * driver_constants.c_ice
        )
        t_dt = t_dt + rdt * (t1 - t) * cvm / driver_constants.cp_air

        # -----------------------------------------------------------------------
        # update cloud fraction tendency
        # -----------------------------------------------------------------------
        if do_qa == False:
            qa_dt = qa_dt + rdt * (
                qa * sqrt((qi1 + ql1) / max(qi + ql, driver_constants.qcmin)) - qa
            )  # New Cloud - Old CloudCloud

    with computation(FORWARD), interval(0, 1):
        # convert to mm / day
        conversion_factor = 86400.0 * rdt * driver_constants.rgrav
        rain = rain * conversion_factor
        snow = snow * conversion_factor
        ice = ice * conversion_factor
        graupel = graupel * conversion_factor
