import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import FORWARD, PARALLEL, computation, exp, i32, interval, log, max, sqrt, trunc

from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.constants import constants
from pyMoist.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat
from pyMoist.GFDL_1M.driver.stencils import wqs2
from pyMoist.shared_incloud_processes import ice_fraction


@gtscript.function
def new_ice_condensate(t, ql, qi, cnv_frc, srf_type):
    """
    Calculate amount of new ice to be frozen at a given point

    reference Fortran: gfdl_cloud_microphys.F90: function new_ice_condensate
    """
    ifrac = ice_fraction(t, cnv_frc, srf_type)
    new_ice_condensate = min(max(0.0, ifrac * (ql + qi) - qi), ql)

    return new_ice_condensate


@gtscript.function
def icloud_melt_freeze(
    t: Float,
    qv: Float,
    ql: Float,
    qr: Float,
    qi: Float,
    qs: Float,
    qg: Float,
    qa: Float,
    den: Float,
    cnv_frc: Float,
    srf_type: Float,
    c_air: Float,
    c_vap: Float,
    fac_frz: Float,
    fac_imlt: Float,
    qi0_crt: Float,
    ql_mlt: Float,
):
    """
    Melting and freezing of cloud ice/water within the icloud subroutine

    reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud
    """

    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    lhi = constants.LI00 + constants.DC_ICE * t
    q_liq = ql + qr
    q_sol = qi + qs + qg
    cvm = c_air + qv * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
    icpk = lhi / cvm

    newice = max(0.0, qi + new_ice_condensate(t, ql, qi, cnv_frc, srf_type))
    newliq = max(0.0, ql + qi - newice)

    melt = fac_imlt * max(0.0, newliq - ql)
    frez = fac_frz * max(0.0, newice - qi)

    if melt > 0.0 and t > constants.TICE and qi > constants.QCMIN:
        # -----------------------------------------------------------------------
        # pimlt: melting of cloud ice
        # -----------------------------------------------------------------------
        if ql_mlt - ql > 0:
            ans = ql_mlt - ql
        else:
            ans = 0
        tmp = min(melt, ans)  # max ql amount

        # new total condensate / old condensate
        qa = max(
            0.0,
            min(
                1.0,
                qa * max(qi + ql - melt + tmp, 0.0) / max(qi + ql, constants.QCMIN),
            ),
        )

        ql = ql + tmp
        qr = qr + melt - tmp
        qi = qi - melt
        q_liq = q_liq + melt
        q_sol = q_sol - melt
        cvm = c_air + qv * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t - melt * lhi / cvm
    elif frez > 0.0 and t <= constants.TICE and ql > constants.QCMIN:
        # -----------------------------------------------------------------------
        # pihom: homogeneous freezing of cloud water into cloud ice
        # this is the 1st occurance of liquid water freezing in the split mp process
        # -----------------------------------------------------------------------
        qi_crt = ice_fraction(t, cnv_frc, srf_type) * qi0_crt / den
        if qi_crt - qi > 0:
            ans = qi_crt
        else:
            ans = 0
        tmp = min(frez, ans)

        # new total condensate / old condensate
        qa = max(
            0.0,
            min(
                1.0,
                qa * max(qi + ql - frez + tmp, 0.0) / max(qi + ql, constants.QCMIN),
            ),
        )

        ql = ql - frez
        qs = qs + frez - tmp
        qi = qi + tmp
        q_liq = q_liq - frez
        q_sol = q_sol + frez
        cvm = c_air + qv * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t + frez * lhi / cvm

    return t, qv, ql, qr, qi, qs, qg, qa, cvm, q_liq, q_sol


@gtscript.function
def acr3d(
    v1: Float,
    v2: Float,
    q1: Float,
    q2: Float,
    c: Float,
    cac_1: Float,
    cac_2: Float,
    cac_3: Float,
    rho: Float,
):
    """
    Compute accretion according to Lin et al. 1983

    reference Fortran: gfdl_cloud_microphys.F90: function acr3d
    """
    t1 = sqrt(q1 * rho)
    s1 = sqrt(q2 * rho)
    s2 = sqrt(s1)  # s1 = s2 ** 2
    acr3d = c * abs(v1 - v2) * q1 * s2 * (cac_1 * t1 + cac_2 * sqrt(t1) * s2 + cac_3 * s1)

    return acr3d


@gtscript.function
def smow_melt(
    tc: Float,
    dqs: Float,
    qsrho: Float,
    psacw: Float,
    psacr: Float,
    c_0: Float,
    c_1: Float,
    c_2: Float,
    c_3: Float,
    c_4: Float,
    rho: Float,
    rhofac: Float,
):
    """
    Melting of snow function (lin et al. 1983)
    note: psacw and psacr must be calc before smlt is called

    reference Fortran: gfdl_cloud_microphys.F90: function smlt
    """
    smow_melt = (c_0 * tc / rho - c_1 * dqs) * (
        c_2 * sqrt(qsrho) + c_3 * qsrho ** 0.65625 * sqrt(rhofac)
    ) + c_4 * tc * (psacw + psacr)

    return smow_melt


@gtscript.function
def graupel_melt(
    tc: Float,
    dqs: Float,
    qgrho: Float,
    pgacw: Float,
    pgacr: Float,
    c_0: Float,
    c_1: Float,
    c_2: Float,
    c_3: Float,
    c_4: Float,
    rho: Float,
):
    """
    Melting of graupel function (lin et al. 1983)
    note: pgacw and pgacr must be calc before gmlt is called

    reference Fortran: gfdl_cloud_microphys.F90: function gmlt
    """
    graupel_melt = (c_0 * tc / rho - c_1 * dqs) * (
        c_2 * sqrt(qgrho) + c_3 * qgrho ** 0.6875 / rho ** 0.25
    ) + c_4 * tc * (pgacw + pgacr)

    return graupel_melt


@gtscript.function
def snow_graupel_coldrain(
    t1: Float,
    qv1: Float,
    ql1: Float,
    qi1: Float,
    qr1: Float,
    qs1: Float,
    qg1: Float,
    qa1: Float,
    p_dry: Float,
    den1: Float,
    denfac: Float,
    vtr: Float,
    vts: Float,
    vtg: Float,
    cvm: Float,
    lhl: Float,
    lhi: Float,
    lcpk: Float,
    icpk: Float,
    tcpk: Float,
    q_liq: Float,
    q_sol: Float,
    di: Float,
    cnv_frc: Float,
    srf_type: Float,
    c_air: Float,
    c_vap: Float,
    cgaci: Float,
    cgacs: Float,
    cgacw: Float,
    cgfr_0: Float,
    cgfr_1: Float,
    cgmlt_0: Float,
    cgmlt_1: Float,
    cgmlt_2: Float,
    cgmlt_3: Float,
    cgmlt_4: Float,
    const_vi: bool,
    csaci: Float,
    csacw: Float,
    csmlt_0: Float,
    csmlt_1: Float,
    csmlt_2: Float,
    csmlt_3: Float,
    csmlt_4: Float,
    dts: Float,
    fac_i2s: Float,
    qi0_crt: Float,
    qs0_crt: Float,
    qs_mlt: Float,
    rdts: Float,
):
    """
    Snow, graupel, cold rain microphysics
    melting, freezing, accretion

    reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud
    """

    t2 = t1
    qv2 = qv1
    ql2 = ql1
    qi2 = qi1
    qr2 = qr1
    qs2 = qs1
    qg2 = qg1

    pgacr = 0.0
    pgacw = 0.0
    tc = t2 - constants.TICE

    if tc >= 0.0:
        # -----------------------------------------------------------------------
        # melting of snow
        # -----------------------------------------------------------------------

        dqs0 = constants.CES0 / p_dry - qv2

        if qs2 > constants.QPMIN:

            # -----------------------------------------------------------------------
            # psacw: accretion of cloud water by snow
            # only rate is used (for snow melt) since tc > 0.
            # -----------------------------------------------------------------------

            if ql2 > constants.QCMIN:
                factor = denfac * csacw * exp(0.8125 * log(qs2 * den1))
                psacw = factor / (1.0 + dts * factor) * ql2  # rate
            else:
                psacw = 0.0

            # -----------------------------------------------------------------------
            # psacr: accretion of rain by melted snow
            # pracs: accretion of snow by rain
            # -----------------------------------------------------------------------

            if qr2 > constants.QPMIN:
                psacr = min(
                    acr3d(
                        vts,
                        vtr,
                        qr2,
                        qs2,
                        constants.CSARC,
                        constants.ACCO_01,
                        constants.ACCO_11,
                        constants.ACCO_21,
                        den1,
                    ),
                    qr2 * rdts,
                )
                pracs = acr3d(
                    vtr,
                    vts,
                    qs2,
                    qr2,
                    constants.CRACS,
                    constants.ACCO_00,
                    constants.ACCO_10,
                    constants.ACCO_20,
                    den1,
                )
            else:
                psacr = 0.0
                pracs = 0.0

            # -----------------------------------------------------------------------
            # total snow sink:
            # psmlt: snow melt (due to rain accretion)
            # -----------------------------------------------------------------------

            psmlt = max(
                0.0,
                smow_melt(
                    tc,
                    dqs0,
                    qs2 * den1,
                    psacw,
                    psacr,
                    csmlt_0,
                    csmlt_1,
                    csmlt_2,
                    csmlt_3,
                    csmlt_4,
                    den1,
                    denfac,
                ),
            )
            sink = min(qs2, min(dts * (psmlt + pracs), tc / icpk))
            qs2 = qs2 - sink
            # sjl, 20170321:
            if qs_mlt - ql2 > 0:
                ans = qs_mlt - ql2
            else:
                ans = 0
            tmp = min(sink, ans)  # max ql due to snow melt

            # new total condensate / old condensate
            qa1 = max(
                0.0,
                min(
                    1.0,
                    qa1 * max(qi2 + ql2 + tmp, 0.0) / max(qi2 + ql2, constants.QCMIN),
                ),
            )

            ql2 = ql2 + tmp
            qr2 = qr2 + sink - tmp
            # sjl, 20170321:
            q_liq = q_liq + sink
            q_sol = q_sol - sink
            cvm = c_air + qv2 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t2 = t2 - sink * lhi / cvm
            tc = t2 - constants.TICE

        # -----------------------------------------------------------------------
        # update capacity heat and latent heat coefficient
        # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * t2
        icpk = lhi / cvm

        # -----------------------------------------------------------------------
        # melting of graupel
        # -----------------------------------------------------------------------

        if qg2 > constants.QPMIN and tc > 0.0:

            # -----------------------------------------------------------------------
            # pgacr: accretion of rain by graupel
            # -----------------------------------------------------------------------

            if qr2 > constants.QPMIN:
                pgacr = min(
                    acr3d(
                        vtg,
                        vtr,
                        qr2,
                        qg2,
                        constants.CGARC,
                        constants.ACCO_02,
                        constants.ACCO_12,
                        constants.ACCO_22,
                        den1,
                    ),
                    rdts * qr2,
                )

            # -----------------------------------------------------------------------
            # pgacw: accretion of cloud water by graupel
            # -----------------------------------------------------------------------

            qden = qg2 * den1
            if ql2 > constants.QCMIN:
                factor = cgacw * qden / sqrt(den1 * sqrt(sqrt(qden)))
                pgacw = factor / (1.0 + dts * factor) * ql2  # rate

            # -----------------------------------------------------------------------
            # pgmlt: graupel melt
            # -----------------------------------------------------------------------

            pgmlt = dts * graupel_melt(
                tc,
                dqs0,
                qden,
                pgacw,
                pgacr,
                cgmlt_0,
                cgmlt_1,
                cgmlt_2,
                cgmlt_3,
                cgmlt_4,
                den1,
            )
            pgmlt = min(max(0.0, pgmlt), min(qg2, tc / icpk))
            qg2 = qg2 - pgmlt
            qr2 = qr2 + pgmlt
            q_liq = q_liq + pgmlt
            q_sol = q_sol - pgmlt
            cvm = c_air + qv2 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t2 = t2 - pgmlt * lhi / cvm

    else:

        # -----------------------------------------------------------------------
        # cloud ice proc:
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # psaci: accretion of cloud ice by snow
        # -----------------------------------------------------------------------

        if qi2 > 3.0e-7:  # cloud ice sink terms
            if qs2 > constants.QPMIN:
                # -----------------------------------------------------------------------
                # sjl added (following lin eq. 23) the temperature dependency
                # to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
                # -----------------------------------------------------------------------
                factor = dts * denfac * csaci * exp(0.05 * tc + 0.8125 * log(qs2 * den1))
                psaci = factor / (1.0 + factor) * qi2
            else:
                psaci = 0.0

            # -----------------------------------------------------------------------
            # psaut: autoconversion: cloud ice -- > snow
            # -----------------------------------------------------------------------

            # -----------------------------------------------------------------------
            # similar to lfo 1983: eq. 21 solved implicitly
            # threshold from wsm6 scheme, hong et al 2004, eq (13) : qi0_crt ~0.8e-4
            # -----------------------------------------------------------------------

            qim = ice_fraction(t2, cnv_frc, srf_type) * qi0_crt / den1

            # -----------------------------------------------------------------------
            # assuming linear subgrid vertical distribution of cloud ice
            # the mismatch computation following lin et al. 1994, mwr
            # -----------------------------------------------------------------------

            if const_vi:
                tmp = fac_i2s
            else:
                tmp = fac_i2s * exp(0.025 * tc)

            di = max(di, constants.QCMIN)
            q_plus = qi2 + di
            if q_plus > (qim + constants.QCMIN):
                if qim > (qi2 - di):
                    dq = (0.25 * (q_plus - qim) ** 2) / di
                else:
                    dq = qi2 - qim
                psaut = tmp * dq
            else:
                psaut = 0.0
            sink = min(qi2, psaci + psaut)

            # new total condensate / old condensate
            qa1 = max(
                0.0,
                min(
                    1.0,
                    qa1 * max(qi2 + ql2 - sink + tmp, 0.0) / max(qi2 + ql2, constants.QCMIN),
                ),
            )

            qi2 = qi2 - sink
            qs2 = qs2 + sink

            # -----------------------------------------------------------------------
            # pgaci: accretion of cloud ice by graupel
            # -----------------------------------------------------------------------

            if qg2 > constants.QPMIN:
                # -----------------------------------------------------------------------
                # factor = dts * cgaci / sqrt (den (k)) *
                # exp (0.05 * tc + 0.875 * log (qg * den (k)))
                # simplified form: remove temp dependency & set the exponent 0.875 -> 1
                # -----------------------------------------------------------------------
                factor = dts * cgaci * sqrt(den1) * qg2
                pgaci = factor / (1.0 + factor) * qi2
                qi2 = qi2 - pgaci
                qg2 = qg2 + pgaci

        # -----------------------------------------------------------------------
        # cold - rain proc:
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # rain to ice, snow, graupel processes:
        # -----------------------------------------------------------------------

        tc = t2 - constants.TICE

        if qr2 > constants.QPMIN and tc < 0.0:

            # -----------------------------------------------------------------------
            # * sink * terms to qr: psacr + pgfr
            # source terms to qs: psacr
            # source terms to qg: pgfr
            # -----------------------------------------------------------------------

            # -----------------------------------------------------------------------
            # psacr accretion of rain by snow
            # -----------------------------------------------------------------------

            if qs2 > constants.QPMIN:  # if snow exists
                psacr = dts * acr3d(
                    vts,
                    vtr,
                    qr2,
                    qs2,
                    constants.CSARC,
                    constants.ACCO_01,
                    constants.ACCO_11,
                    constants.ACCO_21,
                    den1,
                )
            else:
                psacr = 0.0

            # -----------------------------------------------------------------------
            # pgfr: rain freezing -- > graupel
            # -----------------------------------------------------------------------

            pgfr = dts * cgfr_0 / den1 * (exp(-cgfr_1 * tc) - 1.0) * exp(1.75 * log(qr2 * den1))

            # -----------------------------------------------------------------------
            # total sink to qr
            # -----------------------------------------------------------------------

            sink = psacr + pgfr
            factor = min(sink, min(qr2, -tc / icpk)) / max(sink, constants.QPMIN)

            psacr = factor * psacr
            pgfr = factor * pgfr

            sink = psacr + pgfr
            qr2 = qr2 - sink
            qs2 = qs2 + psacr
            qg2 = qg2 + pgfr
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = c_air + qv2 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t2 = t2 + sink * lhi / cvm

        # # -----------------------------------------------------------------------
        # # update capacity heat and latent heat coefficient
        # # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * t2
        icpk = lhi / cvm

        # # -----------------------------------------------------------------------
        # # graupel production terms:
        # # -----------------------------------------------------------------------

        if qs2 > constants.QPMIN:

            # -----------------------------------------------------------------------
            # accretion: snow -- > graupel
            # -----------------------------------------------------------------------

            if qg2 > constants.QPMIN:
                sink = dts * acr3d(
                    vtg,
                    vts,
                    qs2,
                    qs2,
                    cgacs,
                    constants.ACCO_03,
                    constants.ACCO_13,
                    constants.ACCO_23,
                    den1,
                )
            else:
                sink = 0.0

            # -----------------------------------------------------------------------
            # autoconversion snow -- > graupel
            # -----------------------------------------------------------------------

            qsm = qs0_crt / den1
            if qs2 > qsm:
                factor = dts * 1.0e-3 * exp(0.09 * (t2 - constants.TICE))
                sink = sink + factor / (1.0 + factor) * (qs2 - qsm)
            sink = min(qs2, sink)

            # snow existed
            qs2 = qs2 - sink
            qg2 = qg2 + sink

        if qg2 > constants.QPMIN and t2 < (constants.TICE - 0.01):

            # -----------------------------------------------------------------------
            # pgacw: accretion of cloud water by graupel
            # -----------------------------------------------------------------------

            if ql2 > constants.QCMIN:
                qden = qg2 * den1
                factor = dts * cgacw * qden / sqrt(den1 * sqrt(sqrt(qden)))
                pgacw = factor / (1.0 + factor) * ql2
            else:
                pgacw = 0.0

            # -----------------------------------------------------------------------
            # pgacr: accretion of rain by graupel
            # -----------------------------------------------------------------------

            if qr2 > constants.QPMIN:
                pgacr = min(
                    dts
                    * acr3d(
                        vtg,
                        vtr,
                        qr2,
                        qg2,
                        constants.CGARC,
                        constants.ACCO_02,
                        constants.ACCO_12,
                        constants.ACCO_22,
                        den1,
                    ),
                    qr2,
                )
            else:
                pgacr = 0.0

            sink = pgacr + pgacw
            if constants.TICE - t2 > 0:
                ans = constants.TICE - t2
            else:
                ans = 0
            factor = min(sink, ans / icpk) / max(sink, constants.QPMIN)
            pgacr = factor * pgacr
            pgacw = factor * pgacw

            sink = pgacr + pgacw
            qg2 = qg2 + sink
            qr2 = qr2 - pgacr
            ql2 = ql2 - pgacw
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = c_air + qv2 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t2 = t2 + sink * lhi / cvm

    t1 = t2
    qv1 = qv2
    ql1 = ql2
    qi1 = qi2
    qr1 = qr2
    qs1 = qs2
    qg1 = qg2

    return (
        t1,
        qv1,
        ql1,
        qi1,
        qr1,
        qs1,
        qg1,
        qa1,
        p_dry,
        den1,
        denfac,
        vtr,
        vts,
        vtg,
        cvm,
        lhl,
        lhi,
        lcpk,
        icpk,
        tcpk,
    )


@gtscript.function
def iqs1(
    ta: Float,
    den: Float,
    table3: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
):
    """
    Compute saturation specific humidity from table3

    water - ice phase; universal dry / moist formular using air density
    input "den" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function iqs1
    """
    tmin = constants.TABLE_ICE - 160.0
    if ta - tmin > 0:
        ans = ta - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = i32(trunc(ap1))
    es = table3.A[it - 1] + (ap1 - it) * des3.A[it - 1]
    iqs1 = es / (constants.RVGAS * ta * den)

    return iqs1


@gtscript.function
def iqs2(
    ta: Float,
    den: Float,
    table3: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
):
    """
    Compute saturation specific humidity from table3
    with additional calculation of gradient (dq/dt)

    water - ice phase; universal dry / moist formular using air density
    input "den" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function iqs2
    """
    tmin = constants.TABLE_ICE - 160.0
    if ta - tmin > 0:
        ans = ta - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = i32(trunc(ap1))
    es = table3.A[it - 1] + (ap1 - it) * des3.A[it - 1]
    iqs2 = es / (constants.RVGAS * ta * den)
    it = i32(trunc(ap1 - 0.5))
    dqdt = 10.0 * (des3.A[it - 1] + (ap1 - it) * (des3.A[it] - des3.A[it - 1])) / (constants.RVGAS * ta * den)

    return iqs2, dqdt


@gtscript.function
def wqs1(
    ta: Float,
    den: Float,
    table2: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
):
    """
    Compute the saturated specific humidity for table2

    pure water phase; universal dry / moist formular using air density
    input "den" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function wqs1
    """
    tmin = constants.TABLE_ICE - 160.0
    if ta - tmin > 0:
        ans = ta - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = i32(trunc(ap1))
    es = table2.A[it - 1] + (ap1 - it) * des2.A[it - 1]
    wqs1 = es / (constants.RVGAS * ta * den)

    return wqs1


@gtscript.function
def subgrid_z_proc(
    p_dry: Float,
    den1: Float,
    denfac: Float,
    t1: Float,
    qv1: Float,
    ql1: Float,
    qr1: Float,
    qi1: Float,
    qs1: Float,
    qg1: Float,
    qa1: Float,
    rh_limited: Float,
    ccn: Float,
    cnv_frc: Float,
    srf_type: Float,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
    c_air: Float,
    c_vap: Float,
    cssub_0: Float,
    cssub_1: Float,
    cssub_2: Float,
    cssub_3: Float,
    cssub_4: Float,
    d0_vap: Float,
    do_bigg: bool,
    do_evap: bool,
    do_qa: bool,
    dts: Float,
    fac_frz: Float,
    fac_g2v: Float,
    fac_l2v: Float,
    fac_s2v: Float,
    fac_v2g: Float,
    fac_v2s: Float,
    icloud_f: Float,
    lat2: Float,
    lv00: Float,
    preciprad: bool,
    qc_crt: Float,
    qi_lim: Float,
    rh_inc: Float,
    rh_inr: Float,
    t_min: Float,
    t_sub: Float,
):
    """
    Temperature sensitive high vertical resolution processes

    reference Fortran: gfdl_cloud_microphys.F90: subroutine subgrid_z_proc
    """

    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    lhl = lv00 + d0_vap * t1
    lhi = constants.LI00 + constants.DC_ICE * t1
    q_liq = ql1 + qr1
    q_sol = qi1 + qs1 + qg1
    cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
    lcpk = lhl / cvm
    icpk = lhi / cvm
    tcpk = lcpk + icpk
    if constants.TICE - t1 > 0:
        ans = constants.TICE - t1
    else:
        ans = 0
    tcp3 = lcpk + icpk * min(1.0, ans / (constants.TICE - constants.T_WFR))

    rh_adj = 1.0 - rh_limited - rh_inc
    rh_rain = max(0.35, rh_adj - rh_inr)

    subl1 = 0.0

    cycle = False
    if p_dry < constants.P_MIN:
        cycle = True

    # -----------------------------------------------------------------------
    # instant deposit all water vapor to cloud ice when temperature is super low
    # -----------------------------------------------------------------------

    if t1 < t_min and cycle == False:  # noqa
        if qv1 - constants.QVMIN > 0:
            sink = qv1 - constants.QVMIN
        else:
            sink = 0
        qv1 = qv1 - sink
        qi1 = qi1 + sink
        q_sol = q_sol + sink
        cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t1 = t1 + sink * (lhl + lhi) / cvm
        if do_qa == True:  # noqa
            qa1 = 1.0  # air fully saturated; 100 % cloud cover
        cycle = True

    if cycle == False:  # noqa
        # -----------------------------------------------------------------------
        # update heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------
        lhl = lv00 + d0_vap * t1
        lhi = constants.LI00 + constants.DC_ICE * t1
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk
        if constants.TICE - t1 > 0:
            ans = constants.TICE - t1
        else:
            ans = 0
        tcp3 = lcpk + icpk * min(1.0, ans / (constants.TICE - constants.T_WFR))

        # -----------------------------------------------------------------------
        # instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
        # -----------------------------------------------------------------------
        qpz = qv1 + ql1 + qi1
        tin = t1 - (lhl * (ql1 + qi1) + lhi * qi1) / (
            c_air + qpz * c_vap + qr1 * constants.C_LIQ + (qs1 + qg1) * constants.C_ICE
        )
        if tin > t_sub + 6.0:
            rh = qpz / iqs1(tin, den1, table3, des3)
            if rh < rh_adj:  # qpz / rh_adj < qs
                t1 = tin
                qv1 = qpz
                ql1 = 0.0
                qi1 = 0.0
                if do_qa == True:  # noqa
                    qa = 0.0
                cycle = True  # cloud free

    if cycle == False:  # noqa
        # -----------------------------------------------------------------------
        # cloud water < -- > vapor adjustment: LS evaporation
        # -----------------------------------------------------------------------
        if do_evap == True:  # noqa
            qsw, dwsdt = wqs2(t1, den1, table2, des2)
            dq0 = qsw - qv1
            if dq0 > constants.QVMIN:
                factor = min(1.0, fac_l2v * (10.0 * dq0 / qsw))
                evap = min(ql1, factor * ql1 / (1.0 + tcp3 * dwsdt))
            else:
                evap = 0.0
            qv1 = qv1 + evap
            ql1 = ql1 - evap
            q_liq = q_liq - evap
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 - evap * lhl / cvm

        # -----------------------------------------------------------------------
        # update heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * t1
        icpk = lhi / cvm

        # -----------------------------------------------------------------------
        # enforce complete freezing when ice_fraction==1
        # -----------------------------------------------------------------------

        ifrac = ice_fraction(t1, cnv_frc, srf_type)
        if ifrac == 1.0 and ql1 > constants.QCMIN:
            sink = ql1
            ql1 = ql1 - sink
            qi1 = qi1 + sink
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 + sink * lhi / cvm

        # -----------------------------------------------------------------------
        # update heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * t1
        icpk = lhi / cvm

        # -----------------------------------------------------------------------
        # bigg mechanism heterogeneous freezing on existing cloud nuclei
        # -----------------------------------------------------------------------

        tc = constants.TICE - t1
        if do_bigg == True and ql1 > constants.QCMIN and tc > 0.0:  # noqa
            sink = fac_frz * (100.0 / constants.RHOR / ccn) * dts * (exp(0.66 * tc) - 1.0) * den1 * ql1 * ql1
            sink = min(ql1, min(tc / icpk, sink))
            ql1 = ql1 - sink
            qi1 = qi1 + sink
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 + sink * lhi / cvm
            # significant ql existed

        # -----------------------------------------------------------------------
        # update capacity heat and latent heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        lhi = constants.LI00 + constants.DC_ICE * t1
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # sublimation / deposition of LS ice
        # -----------------------------------------------------------------------

        if t1 < constants.TICE:
            qsi, dqsdt = iqs2(t1, den1, table3, des3)
            dq = qv1 - qsi
            sink = min(qi1, dq / (1.0 + tcpk * dqsdt))
            if qi1 > constants.QCMIN:
                # eq 9, hong et al. 2004, mwr
                # for a and b, see dudhia 1989: page 3103 eq (b7) and (b8)
                pidep = (
                    dts
                    * dq
                    * 349138.78
                    * exp(0.875 * log(qi1 * den1))
                    / (qsi * den1 * lat2 / (0.0243 * constants.RVGAS * t1 ** 2) + 4.42478e4)
                )
            else:
                pidep = 0.0
            if dq > 0.0:  # vapor - > ice
                # deposition
                ifrac = ice_fraction(t1, cnv_frc, srf_type)
                tmp = constants.TICE - t1
                qi_crt = 4.92e-11 * exp(1.33 * log(1.0e3 * exp(0.1 * tmp)))
                qi_crt = max(qi_crt, 1.82e-6) * qi_lim * ifrac / den1
                sink = min(sink, min(max(qi_crt - qi1, pidep), tmp / tcpk))
            else:  # ice -- > vapor
                # NOTE sublimation not implemented
                # trigger checked in driver `check_flags` function

                # dev NOTE: unsure how to handle pssub. In Fortran this variable is
                # initalized to nan then used here (at least when do_subl is False,
                # maybe do_subl has other unknown effects)
                # # sublimation
                # if do_subl == True: #noqa
                #     if t1 - t_sub > 0:
                #         ans = t1 - t_sub
                #     else:
                #         ans = 0
                #     pidep = pidep * min(1.0, ans * 0.2)
                #     sink = fac_i2v * max(pidep, sink, -qi1)
                #     subl1 = subl1 + pssub / dts
                # else:
                #     sink = 0.0
                sink = 0
            qv1 = qv1 - sink
            qi1 = qi1 + sink
            q_sol = q_sol + sink
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 + sink * (lhl + lhi) / cvm

        # -----------------------------------------------------------------------
        # update capacity heat and latend heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        lhi = constants.LI00 + constants.DC_ICE * t1
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # sublimation / deposition of snow
        # this process happens for all temp rage
        # -----------------------------------------------------------------------

        if qs1 > constants.QPMIN:
            qsi, dqsdt = iqs2(t1, den1, table3, des3)
            qden = qs1 * den1
            tmp = exp(0.65625 * log(qden))
            tsq = t1 * t1
            dq = (qsi - qv1) / (1.0 + tcpk * dqsdt)
            pssub = (
                cssub_0
                * tsq
                * (cssub_1 * sqrt(qden) + cssub_2 * tmp * sqrt(denfac))
                / (cssub_3 * tsq + cssub_4 * qsi * den1)
            )
            pssub = (qsi - qv1) * dts * pssub
            if pssub > 0.0:  # qs -- > qv, sublimation
                if t1 - t_sub > 0:
                    ans = t1 - t_sub
                else:
                    ans = 0
                pssub = min(fac_s2v * pssub * min(1.0, ans * 0.2), qs1)
                subl1 = subl1 + pssub / dts
            else:
                if t1 > constants.TICE:
                    pssub = 0.0  # no deposition
                else:
                    pssub = max(fac_v2s * pssub, max(dq, (t1 - constants.TICE) / tcpk))
            qs1 = qs1 - pssub
            qv1 = qv1 + pssub
            q_sol = q_sol - pssub
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 - pssub * (lhl + lhi) / cvm

        # -----------------------------------------------------------------------
        # update capacity heat and latend heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        lhi = constants.LI00 + constants.DC_ICE * t1
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # simplified 2 - way grapuel sublimation - deposition mechanism
        # -----------------------------------------------------------------------

        if qg1 > constants.QPMIN:
            qsi, dqsdt = iqs2(t1, den1, table3, des3)
            dq = (qv1 - qsi) / (1.0 + tcpk * dqsdt)
            pgsub = (qv1 / qsi - 1.0) * qg1
            if pgsub > 0.0:  # deposition
                if t1 > constants.TICE:
                    pgsub = 0.0  # no deposition
                else:
                    pgsub = min(
                        fac_v2g * pgsub,
                        min(
                            0.2 * dq,
                            min(ql1 + qr1, (constants.TICE - t1) / tcpk),
                        ),
                    )
            else:  # submilation
                if t1 - t_sub > 0:
                    ans = t1 - t_sub
                else:
                    ans = 0
                pgsub = max(fac_g2v * pgsub, dq) * min(1.0, ans * 0.1)
                subl1 = subl1 + pgsub / dts
            qg1 = qg1 + pgsub
            qv1 = qv1 - pgsub
            q_sol = q_sol + pgsub
            cvm = c_air + qv1 * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t1 = t1 + pgsub * (lhl + lhi) / cvm

            # Fortran ifdef USE_MIN_EVAP goes in here.
            # Not currently executed, so not included

        # -----------------------------------------------------------------------
        # update capacity heat and latend heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t1
        cvm = c_air + (qv1 + q_liq + q_sol) * c_vap
        lcpk = lhl / cvm

        # -----------------------------------------------------------------------
        # compute cloud fraction
        # -----------------------------------------------------------------------
        if do_qa == False:  # noqa
            cycle = True

    if cycle == False:  # noqa
        # -----------------------------------------------------------------------
        # combine water species
        # -----------------------------------------------------------------------
        if preciprad == True:  # noqa
            q_sol = qi1 + qs1 + qg1
            q_liq = ql1 + qr1
        else:
            q_sol = qi1
            q_liq = ql1
        q_cond = q_liq + q_sol

        qpz = qv1 + q_cond  # qpz is conserved

        # -----------------------------------------------------------------------
        # use the "liquid - frozen water temperature" (tin)
        # to compute saturated specific humidity
        # -----------------------------------------------------------------------

        tin = t1 - (lcpk * q_cond + icpk * q_sol)  # minimum temperature

        # -----------------------------------------------------------------------
        # determine saturated specific humidity
        # -----------------------------------------------------------------------

        if tin <= constants.T_WFR:
            # ice phase:
            qstar = iqs1(tin, den1, table3, des3)
        elif tin >= constants.TICE:
            # liquid phase:
            qstar = wqs1(tin, den1, table2, des2)
        else:
            # mixed phase:
            qsi = iqs1(tin, den1, table3, des3)
            qsw = wqs1(tin, den1, table2, des2)
            if q_cond > 3.0e-6:
                rqi = q_sol / q_cond
            else:
                # WMP impose CALIPSO ice polynomial from 0 C to -40 C
                rqi = ice_fraction(tin, cnv_frc, srf_type)
            qstar = rqi * qsi + (1.0 - rqi) * qsw

        # -----------------------------------------------------------------------
        # assuming subgrid linear distribution in horizontal;
        # this is effectively a smoother for the binary cloud scheme
        # -----------------------------------------------------------------------
        if qpz > constants.QCMIN:
            # partial cloudiness by pdf:
            dq = max(constants.QCMIN, rh_limited * qpz)
            q_plus = qpz + dq  # cloud free if qstar > q_plus
            q_minus = qpz - dq
            if icloud_f == 3:
                # triangular
                if q_plus <= qstar:
                    # little/no cloud cover
                    do_nothing = True
                elif qpz <= qstar and qstar < q_plus:  # partial cloud cover
                    qa1 = max(
                        constants.QCMIN,
                        min(
                            1.0,
                            qa1 + (q_plus - qstar) * (q_plus - qstar) / ((q_plus - q_minus) * (q_plus - qpz)),
                        ),
                    )
                elif q_minus <= qstar and qstar < qpz:  # partial cloud cover
                    qa1 = max(
                        constants.QCMIN,
                        min(
                            1.0,
                            qa1
                            + 1.0
                            - (
                                (qstar - q_minus) * (qstar - q_minus) / ((q_plus - q_minus) * (qpz - q_minus))
                            ),
                        ),
                    )
                elif qstar <= q_minus:
                    qa1 = 1.0  # air fully saturated; 100 % cloud cover
            else:
                # top-hat
                if q_plus <= qstar:
                    # little/no cloud cover
                    do_nothing = True
                elif qstar < q_plus and q_cond > qc_crt:
                    qa1 = max(
                        constants.QCMIN,
                        min(1.0, qa1 + (q_plus - qstar) / (dq + dq)),
                    )  # partial cloud cover
                elif qstar <= q_minus:
                    qa1 = 1.0  # air fully saturated; 100 % cloud cover

    return t1, qv1, ql1, qr1, qi1, qs1, qg1, qa1, subl1


def icloud_core(
    t1: FloatField,
    p_dry: FloatField,
    dp1: FloatField,
    qv1: FloatField,
    ql1: FloatField,
    qr1: FloatField,
    qi1: FloatField,
    qs1: FloatField,
    qg1: FloatField,
    qa1: FloatField,
    den1: FloatField,
    denfac: FloatField,
    vts: FloatField,
    vtg: FloatField,
    vtr: FloatField,
    subl1: FloatField,
    rh_limited: FloatField,
    ccn: FloatField,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
):
    """
    Ice cloud processes
    reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud
    Fortran author: Shian-Jiann lin, gfdl
    """

    from __externals__ import (
        c_air,
        c_vap,
        cgaci,
        cgacs,
        cgacw,
        cgfr_0,
        cgfr_1,
        cgmlt_0,
        cgmlt_1,
        cgmlt_2,
        cgmlt_3,
        cgmlt_4,
        const_vi,
        csaci,
        csacw,
        csmlt_0,
        csmlt_1,
        csmlt_2,
        csmlt_3,
        csmlt_4,
        cssub_0,
        cssub_1,
        cssub_2,
        cssub_3,
        cssub_4,
        d0_vap,
        do_bigg,
        do_evap,
        do_qa,
        dts,
        fac_frz,
        fac_g2v,
        fac_i2s,
        fac_imlt,
        fac_l2v,
        fac_s2v,
        fac_v2g,
        fac_v2s,
        icloud_f,
        lat2,
        lv00,
        preciprad,
        qc_crt,
        qi0_crt,
        qi_lim,
        ql_mlt,
        qs0_crt,
        qs_mlt,
        rdts,
        rh_inc,
        rh_inr,
        t_min,
        t_sub,
        z_slope_ice,
    )

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud
    with computation(FORWARD), interval(0, 1):
        # initalize vtr in place of the warm rain calculations
        vtr = 0

    # -----------------------------------------------------------------------
    # sources of cloud ice: pihom, cold rain, and the sat_adj
    # (initiation plus deposition)
    # sources of snow: cold rain, auto conversion + accretion (from cloud ice)
    # sat_adj (deposition; requires pre - existing snow);
    # initial snow comes from auto conversion
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        t1, qv1, ql1, qr1, qi1, qs1, qg1, qa1, cvm, q_liq, q_sol = icloud_melt_freeze(
            t1,
            qv1,
            ql1,
            qr1,
            qi1,
            qs1,
            qg1,
            qa1,
            den1,
            cnv_frc,
            srf_type,
            c_air,
            c_vap,
            fac_frz,
            fac_imlt,
            qi0_crt,
            ql_mlt,
        )

    # begin reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof
    # still within reference Fortran gfdl_cloud_microphys.F90: subroutine icloud
    # set up inputs to "function"
    with computation(PARALLEL), interval(...):
        if z_slope_ice == True:  # noqa
            q_linear_prof = qi1
            h_var_linear_prof = rh_limited
            dm_linear_prof = q_linear_prof  # initalized here to ensure it is created as a 3d field

    with computation(FORWARD), interval(1, None):
        if z_slope_ice == True:  # noqa
            dq_linear_prof = 0.5 * (q_linear_prof - q_linear_prof[0, 0, -1])

    # -----------------------------------------------------------------------
    # use twice the strength of the positive definiteness limiter (lin et al 1994)
    # -----------------------------------------------------------------------
    with computation(FORWARD), interval(0, 1):
        if z_slope_ice == True:  # noqa
            dm_linear_prof = 0

    with computation(FORWARD), interval(1, -1):
        if z_slope_ice == True:  # noqa
            dm_linear_prof = 0.5 * min(abs(dq_linear_prof + dq_linear_prof[0, 0, 1]), 0.5 * q_linear_prof)
            if dq_linear_prof * dq_linear_prof[0, 0, 1] <= 0.0:
                if dq_linear_prof > 0.0:  # local max
                    dm_linear_prof = min(dm_linear_prof, min(dq_linear_prof, -dq_linear_prof[0, 0, 1]))
                else:
                    dm_linear_prof = 0.0

    with computation(FORWARD), interval(-1, None):
        if z_slope_ice == True:  # noqa
            dm_linear_prof = 0

    # -----------------------------------------------------------------------
    # impose a presumed background horizontal variability
    # that is proportional to the value itself
    # -----------------------------------------------------------------------
    with computation(PARALLEL), interval(...):
        if z_slope_ice == True:  # noqa
            dm_linear_prof = max(
                dm_linear_prof,
                max(constants.QVMIN, h_var_linear_prof * q_linear_prof),
            )
        if z_slope_ice == False:  # noqa
            dm_linear_prof = max(constants.QVMIN, h_var_linear_prof * q_linear_prof)

    # handle outputs of "function"
    with computation(PARALLEL), interval(...):
        di = dm_linear_prof

    # end reference Fortran: gfdl_cloud_microphys.F90: subroutine linear_prof

    with computation(PARALLEL), interval(...):
        # -----------------------------------------------------------------------
        # update capacity heat and latent heat coefficient
        # -----------------------------------------------------------------------
        lhl = lv00 + d0_vap * t1
        lhi = constants.LI00 + constants.DC_ICE * t1
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # do nothing above p_min
        # -----------------------------------------------------------------------
        if p_dry >= constants.P_MIN:
            (
                t1,
                qv1,
                ql1,
                qi1,
                qr1,
                qs1,
                qg1,
                qa1,
                p_dry,
                den1,
                denfac,
                vtr,
                vts,
                vtg,
                cvm,
                lhl,
                lhi,
                lcpk,
                icpk,
                tcpk,
            ) = snow_graupel_coldrain(
                t1=t1,
                qv1=qv1,
                ql1=ql1,
                qi1=qi1,
                qr1=qr1,
                qs1=qs1,
                qg1=qg1,
                qa1=qa1,
                p_dry=p_dry,
                den1=den1,
                denfac=denfac,
                vtr=vtr,
                vts=vts,
                vtg=vtg,
                cvm=cvm,
                lhl=lhl,
                lhi=lhi,
                lcpk=lcpk,
                icpk=icpk,
                tcpk=tcpk,
                q_liq=q_liq,
                q_sol=q_sol,
                di=di,
                cnv_frc=cnv_frc,
                srf_type=srf_type,
                c_air=c_air,
                c_vap=c_vap,
                cgaci=cgaci,
                cgacs=cgacs,
                cgacw=cgacw,
                cgfr_0=cgfr_0,
                cgfr_1=cgfr_1,
                cgmlt_0=cgmlt_0,
                cgmlt_1=cgmlt_1,
                cgmlt_2=cgmlt_2,
                cgmlt_3=cgmlt_3,
                cgmlt_4=cgmlt_4,
                const_vi=const_vi,
                csaci=csaci,
                csacw=csacw,
                csmlt_0=csmlt_0,
                csmlt_1=csmlt_1,
                csmlt_2=csmlt_2,
                csmlt_3=csmlt_3,
                csmlt_4=csmlt_4,
                dts=dts,
                fac_i2s=fac_i2s,
                qi0_crt=qi0_crt,
                qs0_crt=qs0_crt,
                qs_mlt=qs_mlt,
                rdts=rdts,
            )

        t1, qv1, ql1, qr1, qi1, qs1, qg1, qa1, subl1 = subgrid_z_proc(
            p_dry=p_dry,
            den1=den1,
            denfac=denfac,
            t1=t1,
            qv1=qv1,
            ql1=ql1,
            qr1=qr1,
            qi1=qi1,
            qs1=qs1,
            qg1=qg1,
            qa1=qa1,
            rh_limited=rh_limited,
            ccn=ccn,
            cnv_frc=cnv_frc,
            srf_type=srf_type,
            table2=table2,
            table3=table3,
            des2=des2,
            des3=des3,
            c_air=c_air,
            c_vap=c_vap,
            cssub_0=cssub_0,
            cssub_1=cssub_1,
            cssub_2=cssub_2,
            cssub_3=cssub_3,
            cssub_4=cssub_4,
            d0_vap=d0_vap,
            do_bigg=do_bigg,
            do_evap=do_evap,
            do_qa=do_qa,
            dts=dts,
            fac_frz=fac_frz,
            fac_g2v=fac_g2v,
            fac_l2v=fac_l2v,
            fac_s2v=fac_s2v,
            fac_v2g=fac_v2g,
            fac_v2s=fac_v2s,
            icloud_f=icloud_f,
            lat2=lat2,
            lv00=lv00,
            preciprad=preciprad,
            qc_crt=qc_crt,
            qi_lim=qi_lim,
            rh_inc=rh_inc,
            rh_inr=rh_inr,
            t_min=t_min,
            t_sub=t_sub,
        )
        # end reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud


def update_precip_total(
    isubl: FloatField,
    subl1: FloatField,
):
    """
    Update precipitation totals with results of icloud stencil

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv
    """
    with computation(PARALLEL), interval(...):
        isubl = isubl + subl1

        subl1 = 0
