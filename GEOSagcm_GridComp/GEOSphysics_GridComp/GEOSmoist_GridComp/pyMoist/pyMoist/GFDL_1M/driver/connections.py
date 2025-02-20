"""Core functions and stencils of the GFDL_1M driver"""

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    computation,
    exp,
    interval,
    log,
    log10,
    sqrt,
)

import pyMoist.GFDL_1M.driver.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.shared_generic_math import sigma


GlobalTable_driver_qsat = gtscript.GlobalTable[(Float, (constants.LENGTH))]


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
    initalize temporary copies of many quantities

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
    adjusts/removes negative mixing ratios

    reference Fortran: gfdl_cloud_microphys.F90: subroutine neg_adj
    """
    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    cvm = (
        c_air
        + qv * c_vap
        + (qr + ql) * constants.C_LIQ
        + (qi + qs + qg) * constants.C_ICE
    )
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
    stencil wrapper for fix_negative_core

    adjusts/removes negative mixing ratios
    updates qv based on new values

    refernce Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    from __externals__ import c_air, c_vap, d0_vap, lv00

    # -----------------------------------------------------------------------
    # fix all negative water species
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(0, -1):
        t, qv, ql, qr, qi, qs, qg = fix_negative_core(
            t, qv, ql, qr, qi, qs, qg, c_air, c_vap, lv00, d0_vap
        )
        if qv < 0.0:
            qv[0, 0, 1] = qv[0, 0, 1] + qv * dp / dp[0, 0, 1]
            qv = 0.0

    with computation(FORWARD), interval(-1, None):
        t, qv, ql, qr, qi, qs, qg = fix_negative_core(
            t, qv, ql, qr, qi, qs, qg, c_air, c_vap, lv00, d0_vap
        )

        if qv < 0.0 and qv[0, 0, -1] > 0.0:
            dq = min(-qv * dp, qv[0, 0, -1] * dp[0, 0, -1])
            qv[0, 0, -1] = qv[0, 0, -1] - dq / dp[0, 0, -1]
            qv = qv + dq / dp


@gtscript.function
def fall_speed_core(
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
    calculate the vertical fall speed of precipitation

    reference Fortran: gfdl_cloud_microphys.F90: subroutine fall_speed
    """
    from __externals__ import (
        const_vg,
        const_vi,
        const_vs,
        vg_fac,
        vg_max,
        vi_fac,
        vi_max,
        vs_fac,
        vs_max,
    )

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
            vts = (
                vs_fac
                * constants.VCONS
                * rhof
                * exp(0.0625 * log(qs * den / constants.NORMS))
            )
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
            vtg = (
                vg_fac
                * constants.VCONG
                * rhof
                * sqrt(sqrt(sqrt(qg * den / constants.NORMG)))
            )
            vtg = min(vg_max, max(constants.VF_MIN, vtg))

    return vti, vts, vtg


def fall_speed(
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
    stencil wrapper for fall_speed_core

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

        vti, vts, vtg = fall_speed_core(
            p_dry, cnv_frc, anv_icefall, lsc_icefall, den1, qs1, qi1, qg1, ql1, t1
        )


def warm_rain_update(
    m1_rain: FloatField,
    m1_sol: FloatField,
    rain1: FloatFieldIJ,
    evap1: FloatField,
    revap: FloatField,
    m2_rain: FloatField,
    m2_sol: FloatField,
    m1: FloatField,
    rain: FloatFieldIJ,
):
    """
    update precipitation totals with results of warm_rain stencil

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    with computation(PARALLEL), interval(...):
        revap = revap + evap1
        m2_rain = m2_rain + m1_rain
        m2_sol = m2_sol + m1_sol
        m1 = m1 + m1_rain + m1_sol

        evap1 = 0
        m1_rain = 0
        m1_sol = 0

    with computation(FORWARD), interval(0, 1):
        rain = rain + rain1

        rain1 = 0


def update_tendencies(
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
    """
    compute output tendencies of the microphysics driver

    reference Fortran: gfdl_cloud_microphys.F90:
    subroutine mpdrv, subroutine gfdl_cloud_microphys_driver
    """
    from __externals__ import c_air, c_vap, do_qa, do_sedi_w, rdt, sedi_transport

    # -----------------------------------------------------------------------
    # momentum transportation during sedimentation
    # note: dp1 is dry mass; dp0 is the old moist (total) mass
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(1, None):
        if sedi_transport == True:  # noqa
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
        qv_dt = qv_dt + rdt * (qv1 - qv0) * omq
        ql_dt = ql_dt + rdt * (ql1 - ql0) * omq
        qr_dt = qr_dt + rdt * (qr1 - qr0) * omq
        qi_dt = qi_dt + rdt * (qi1 - qi0) * omq
        qs_dt = qs_dt + rdt * (qs1 - qs0) * omq
        qg_dt = qg_dt + rdt * (qg1 - qg0) * omq
        cvm = (
            c_air
            + qv1 * c_vap
            + (qr1 + ql1) * constants.C_LIQ
            + (qi1 + qs1 + qg1) * constants.C_ICE
        )
        t_dt = t_dt + rdt * (t1 - t) * cvm / constants.CP_AIR

        # -----------------------------------------------------------------------
        # update cloud fraction tendency
        # -----------------------------------------------------------------------
        if do_qa == False:  # noqa
            qa_dt = qa_dt + rdt * (
                qa0 * sqrt((qi1 + ql1) / max(qi0 + ql0, constants.QCMIN)) - qa0
            )  # New Cloud - Old CloudCloud

    with computation(FORWARD), interval(0, 1):
        # convert to mm / day
        conversion_factor = 86400.0 * rdt * constants.RGRAV
        rain = rain * conversion_factor
        snow = snow * conversion_factor
        ice = ice * conversion_factor
        graupel = graupel * conversion_factor
