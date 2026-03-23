from ndsl import NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import (
    FORWARD,
    PARALLEL,
    computation,
    exp,
    function,
    int32,
    interval,
    log,
    max,
    sqrt,
    trunc,
)
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int

from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.driver.config_constants import GFDL1MDriverConfigDependentConstants
from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.microphysics.GFDL_1M.driver.sat_tables import GFDL_driver_tables, GlobalTable_driver_qsat
from pyMoist.microphysics.GFDL_1M.driver.stencils import wqs2
from pyMoist.shared.incloud_processes import ice_fraction


@function
def new_ice_condensate(
    t: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_ice: Float,
    convection_fraction: Float,
    surface_type: Int,
):
    """Calculate amount of new ice to be frozen at a given point

    reference Fortran: gfdl_cloud_microphys.F90: function new_ice_condensate

    Args:
        t (Float)
        mixing_ratio_liquid (Float)
        mixing_ratio_ice (Float)
        convection_fraction (Float)
        surface_type (Int)
    """
    ifrac = ice_fraction(t, convection_fraction, surface_type)
    new_ice_condensate = min(
        max(0.0, ifrac * (mixing_ratio_liquid + mixing_ratio_ice) - mixing_ratio_ice), mixing_ratio_liquid
    )

    return new_ice_condensate


@function
def icloud_melt_freeze(
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_rain: Float,
    mixing_ratio_ice: Float,
    mixing_ratio_snow: Float,
    mixing_ratio_graupel: Float,
    cloud_fraction: Float,
    density: Float,
    convection_fraction: Float,
    surface_type: Float,
    c_air: Float,
    c_vap: Float,
    fac_frz: Float,
    fac_imlt: Float,
    qi0_crt: Float,
    ql_mlt: Float,
):
    """Melting and freezing of cloud ice/water.

    reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud

    Args:
        t (Float)
        mixing_ratio_vapor (Float)
        mixing_ratio_liquid (Float)
        mixing_ratio_rain (Float)
        mixing_ratio_ice (Float)
        mixing_ratio_snow (Float)
        mixing_ratio_graupel (Float)
        cloud_fraction (Float)
        density (Float)
        convection_fraction (Float)
        surface_type (Float)
        c_air (Float)
        c_vap (Float)
        fac_frz (Float)
        fac_imlt (Float)
        qi0_crt (Float)
        ql_mlt (Float)

    Returns:
        Float: t
        Float: vapor
        Float: liquid
        Float: rain
        Float: ice
        Float: snow
        Float: graupel
        Float: cloud_fraction
        Float: cvm
        Float: q_liq
        Float: q_sol
    """

    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    lhi = constants.LI00 + constants.DC_ICE * t
    q_liq = mixing_ratio_liquid + mixing_ratio_rain
    q_sol = mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel
    cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
    icpk = lhi / cvm

    newice = max(
        0.0,
        mixing_ratio_ice
        + new_ice_condensate(t, mixing_ratio_liquid, mixing_ratio_ice, convection_fraction, surface_type),
    )
    newliq = max(0.0, mixing_ratio_liquid + mixing_ratio_ice - newice)

    melt = fac_imlt * max(0.0, newliq - mixing_ratio_liquid)
    frez = fac_frz * max(0.0, newice - mixing_ratio_ice)

    if melt > 0.0 and t > constants.TICE and mixing_ratio_ice > constants.QCMIN:
        # -----------------------------------------------------------------------
        # pimlt: melting of cloud ice
        # -----------------------------------------------------------------------
        if ql_mlt - mixing_ratio_liquid > 0:
            ans = ql_mlt - mixing_ratio_liquid
        else:
            ans = 0
        tmp = min(melt, ans)  # max mixing_ratio_liquid amount

        # new total condensate / old condensate
        cloud_fraction = max(
            0.0,
            min(
                1.0,
                cloud_fraction
                * max(mixing_ratio_ice + mixing_ratio_liquid - melt + tmp, 0.0)
                / max(mixing_ratio_ice + mixing_ratio_liquid, constants.QCMIN),
            ),
        )

        mixing_ratio_liquid = mixing_ratio_liquid + tmp
        mixing_ratio_rain = mixing_ratio_rain + melt - tmp
        mixing_ratio_ice = mixing_ratio_ice - melt
        q_liq = q_liq + melt
        q_sol = q_sol - melt
        cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t - melt * lhi / cvm
    elif frez > 0.0 and t <= constants.TICE and mixing_ratio_liquid > constants.QCMIN:
        # -----------------------------------------------------------------------
        # pihom: homogeneous freezing of cloud water into cloud ice
        # this is the 1st occurrence of liquid water freezing in the split mp process
        # -----------------------------------------------------------------------
        qi_crt = ice_fraction(t, convection_fraction, surface_type) * qi0_crt / density
        if qi_crt - mixing_ratio_ice > 0:
            ans = qi_crt - mixing_ratio_ice
        else:
            ans = 0
        tmp = min(frez, ans)

        # new total condensate / old condensate
        cloud_fraction = max(
            0.0,
            min(
                1.0,
                cloud_fraction
                * max(mixing_ratio_ice + mixing_ratio_liquid - frez + tmp, 0.0)
                / max(mixing_ratio_ice + mixing_ratio_liquid, constants.QCMIN),
            ),
        )

        mixing_ratio_liquid = mixing_ratio_liquid - frez
        mixing_ratio_snow = mixing_ratio_snow + frez - tmp
        mixing_ratio_ice = mixing_ratio_ice + tmp
        q_liq = q_liq - frez
        q_sol = q_sol + frez
        cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t + frez * lhi / cvm

    return (
        t,
        mixing_ratio_vapor,
        mixing_ratio_liquid,
        mixing_ratio_rain,
        mixing_ratio_ice,
        mixing_ratio_snow,
        mixing_ratio_graupel,
        cloud_fraction,
        cvm,
        q_liq,
        q_sol,
    )


@function
def acr3d(
    terminal_speed_1: Float,
    terminal_speed_2: Float,
    mixing_ratio_1: Float,
    mixing_ratio_2: Float,
    c: Float,
    cac_1: Float,
    cac_2: Float,
    cac_3: Float,
    rho: Float,
):
    """Compute accretion according to Lin et al. 1983

    reference Fortran: gfdl_cloud_microphys.F90: function acr3d

    Args:
        terminal_speed_1 (Float)
        terminal_speed_2 (Float)
        mixing_ratio_1 (Float)
        mixing_ratio_2 (Float)
        c (Float)
        cac_1 (Float)
        cac_2 (Float)
        cac_3 (Float)
        rho (Float)

    Returns:
        Float: acr3d
    """
    t1 = sqrt(mixing_ratio_1 * rho)
    s1 = sqrt(mixing_ratio_2 * rho)
    s2 = sqrt(s1)  # s1 = s2 ** 2
    acr3d = (
        c
        * abs(terminal_speed_1 - terminal_speed_2)
        * mixing_ratio_1
        * s2
        * (cac_1 * t1 + cac_2 * sqrt(t1) * s2 + cac_3 * s1)
    )

    return acr3d


@function
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
    """Melting of snow function (lin et al. 1983)
    note: psacw and psacr must be calc before smlt is called

    reference Fortran: gfdl_cloud_microphys.F90: function smlt

    Args:
        tc (Float)
        dqs (Float)
        qsrho (Float)
        psacw (Float)
        psacr (Float)
        c_0 (Float)
        c_1 (Float)
        c_2 (Float)
        c_3 (Float)
        c_4 (Float)
        rho (Float)
        rhofac (Float)

    Returns:
        Float: smow_melt
    """
    smow_melt = (c_0 * tc / rho - c_1 * dqs) * (
        c_2 * sqrt(qsrho) + c_3 * qsrho**0.65625 * sqrt(rhofac)
    ) + c_4 * tc * (psacw + psacr)

    return smow_melt


@function
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
    """Melting of graupel function (lin et al. 1983)
    note: pgacw and pgacr must be calc before gmlt is called

    reference Fortran: gfdl_cloud_microphys.F90: function gmlt

    Args:
        tc: (Float)
        dqs: (Float)
        qgrho: (Float)
        pgacw: (Float)
        pgacr: (Float)
        c_0: (Float)
        c_1: (Float)
        c_2: (Float)
        c_3: (Float)
        c_4: (Float)
        rho: (Float)

    Returns:
        Float: graupel_melt
    """
    graupel_melt = (c_0 * tc / rho - c_1 * dqs) * (
        c_2 * sqrt(qgrho) + c_3 * qgrho**0.6875 / rho**0.25
    ) + c_4 * tc * (pgacw + pgacr)

    return graupel_melt


@function
def snow_graupel_coldrain(
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_ice: Float,
    mixing_ratio_rain: Float,
    mixing_ratio_snow: Float,
    mixing_ratio_graupel: Float,
    cloud_fraction: Float,
    p_dry: Float,
    density: Float,
    density_factor: Float,
    terminal_fall_rain: Float,
    terminal_fall_snow: Float,
    terminal_fall_graupel: Float,
    cvm: Float,
    lhl: Float,
    lhi: Float,
    lcpk: Float,
    icpk: Float,
    tcpk: Float,
    q_liq: Float,
    q_sol: Float,
    di: Float,
    convection_fraction: Float,
    surface_type: Float,
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
    """Snow, graupel, cold rain microphysics. melting, freezing, accretion

    reference Fortran: gfdl_cloud_microphys.F90: subroutine icloud

    Args:
        t (Float)
        mixing_ratio_vapor (Float)
        mixing_ratio_liquid (Float)
        mixing_ratio_ice (Float)
        mixing_ratio_rain (Float)
        mixing_ratio_snow (Float)
        mixing_ratio_graupel (Float)
        cloud_fraction (Float)
        p_dry (Float)
        density (Float)
        density_factor (Float)
        terminal_fall_rain (Float)
        terminal_fall_snow (Float)
        terminal_fall_graupel (Float)
        cvm (Float)
        lhl (Float)
        lhi (Float)
        lcpk (Float)
        icpk (Float)
        tcpk (Float)
        q_liq (Float)
        q_sol (Float)
        di (Float)
        convection_fraction (Float)
        surface_type (Float)
        c_air (Float)
        c_vap (Float)
        cgaci (Float)
        cgacs (Float)
        cgacw (Float)
        cgfr_0 (Float)
        cgfr_1 (Float)
        cgmlt_0 (Float)
        cgmlt_1 (Float)
        cgmlt_2 (Float)
        cgmlt_3 (Float)
        cgmlt_4 (Float)
        const_vi (bool)
        csaci (Float)
        csacw (Float)
        csmlt_0 (Float)
        csmlt_1 (Float)
        csmlt_2 (Float)
        csmlt_3 (Float)
        csmlt_4 (Float)
        dts (Float)
        fac_i2s (Float)
        qi0_crt (Float)
        qs0_crt (Float)
        qs_mlt (Float)
        rdts (Float)

    Returns:
        Float: t
        Float: vapor
        Float: liquid
        Float: ice
        Float: rain
        Float: snow
        Float: graupel
        Float: cloud_fraction
        Float: p_dry
        Float: density
        Float: density_factor
        Float: terminal_fall_rain
        Float: terminal_fall_snow
        Float: terminal_fall_graupel
        Float: cvm
        Float: lhl
        Float: lhi
        Float: lcpk
        Float: icpk
        Float: tcpk
    """

    internal_t = t
    internal_mixing_ratio_vapor = mixing_ratio_vapor
    internal_mixing_ratio_liquid = mixing_ratio_liquid
    internal_mixing_ratio_ice = mixing_ratio_ice
    internal_mixing_ratio_rain = mixing_ratio_rain
    internal_mixing_ratio_snow = mixing_ratio_snow
    internal_mixing_ratio_graupel = mixing_ratio_graupel

    pgacr = 0.0
    pgacw = 0.0
    tc = internal_t - constants.TICE

    if tc >= 0.0:
        # -----------------------------------------------------------------------
        # melting of snow
        # -----------------------------------------------------------------------

        dqs0 = constants.CES0 / p_dry - internal_mixing_ratio_vapor

        if internal_mixing_ratio_snow > constants.QPMIN:
            # -----------------------------------------------------------------------
            # psacw: accretion of cloud water by snow
            # only rate is used (for snow melt) since tc > 0.
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_liquid > constants.QCMIN:
                factor = density_factor * csacw * exp(0.8125 * log(internal_mixing_ratio_snow * density))
                psacw = factor / (1.0 + dts * factor) * internal_mixing_ratio_liquid  # rate
            else:
                psacw = 0.0

            # -----------------------------------------------------------------------
            # psacr: accretion of rain by melted snow
            # pracs: accretion of snow by rain
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_rain > constants.QPMIN:
                psacr = min(
                    acr3d(
                        terminal_fall_snow,
                        terminal_fall_rain,
                        internal_mixing_ratio_rain,
                        internal_mixing_ratio_snow,
                        constants.CSARC,
                        constants.ACCO_01,
                        constants.ACCO_11,
                        constants.ACCO_21,
                        density,
                    ),
                    internal_mixing_ratio_rain * rdts,
                )
                pracs = acr3d(
                    terminal_fall_rain,
                    terminal_fall_snow,
                    internal_mixing_ratio_snow,
                    internal_mixing_ratio_rain,
                    constants.CRACS,
                    constants.ACCO_00,
                    constants.ACCO_10,
                    constants.ACCO_20,
                    density,
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
                    internal_mixing_ratio_snow * density,
                    psacw,
                    psacr,
                    csmlt_0,
                    csmlt_1,
                    csmlt_2,
                    csmlt_3,
                    csmlt_4,
                    density,
                    density_factor,
                ),
            )
            sink = min(internal_mixing_ratio_snow, min(dts * (psmlt + pracs), tc / icpk))
            internal_mixing_ratio_snow = internal_mixing_ratio_snow - sink
            # sjl, 20170321:
            if qs_mlt - internal_mixing_ratio_liquid > 0:
                ans = qs_mlt - internal_mixing_ratio_liquid
            else:
                ans = 0
            tmp = min(sink, ans)  # max mixing_ratio_liquid due to snow melt

            # new total condensate / old condensate
            cloud_fraction = max(
                0.0,
                min(
                    1.0,
                    cloud_fraction
                    * max(internal_mixing_ratio_ice + internal_mixing_ratio_liquid + tmp, 0.0)
                    / max(internal_mixing_ratio_ice + internal_mixing_ratio_liquid, constants.QCMIN),
                ),
            )

            internal_mixing_ratio_liquid = internal_mixing_ratio_liquid + tmp
            internal_mixing_ratio_rain = internal_mixing_ratio_rain + sink - tmp
            # sjl, 20170321:
            q_liq = q_liq + sink
            q_sol = q_sol - sink
            cvm = (
                c_air
                + internal_mixing_ratio_vapor * c_vap
                + q_liq * constants.C_LIQ
                + q_sol * constants.C_ICE
            )
            internal_t = internal_t - sink * lhi / cvm
            tc = internal_t - constants.TICE

        # -----------------------------------------------------------------------
        # update capacity heat and latent heat coefficient
        # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * internal_t
        icpk = lhi / cvm

        # -----------------------------------------------------------------------
        # melting of graupel
        # -----------------------------------------------------------------------

        if internal_mixing_ratio_graupel > constants.QPMIN and tc > 0.0:
            # -----------------------------------------------------------------------
            # pgacr: accretion of rain by graupel
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_rain > constants.QPMIN:
                pgacr = min(
                    acr3d(
                        terminal_fall_graupel,
                        terminal_fall_rain,
                        internal_mixing_ratio_rain,
                        internal_mixing_ratio_graupel,
                        constants.CGARC,
                        constants.ACCO_02,
                        constants.ACCO_12,
                        constants.ACCO_22,
                        density,
                    ),
                    rdts * internal_mixing_ratio_rain,
                )

            # -----------------------------------------------------------------------
            # pgacw: accretion of cloud water by graupel
            # -----------------------------------------------------------------------

            qden = internal_mixing_ratio_graupel * density
            if internal_mixing_ratio_liquid > constants.QCMIN:
                factor = cgacw * qden / sqrt(density * sqrt(sqrt(qden)))
                pgacw = factor / (1.0 + dts * factor) * internal_mixing_ratio_liquid  # rate

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
                density,
            )
            pgmlt = min(max(0.0, pgmlt), min(internal_mixing_ratio_graupel, tc / icpk))
            internal_mixing_ratio_graupel = internal_mixing_ratio_graupel - pgmlt
            internal_mixing_ratio_rain = internal_mixing_ratio_rain + pgmlt
            q_liq = q_liq + pgmlt
            q_sol = q_sol - pgmlt
            cvm = (
                c_air
                + internal_mixing_ratio_vapor * c_vap
                + q_liq * constants.C_LIQ
                + q_sol * constants.C_ICE
            )
            internal_t = internal_t - pgmlt * lhi / cvm

    else:
        # -----------------------------------------------------------------------
        # cloud ice proc:
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # psaci: accretion of cloud ice by snow
        # -----------------------------------------------------------------------

        if internal_mixing_ratio_ice > 3.0e-7:  # cloud ice sink terms
            if internal_mixing_ratio_snow > constants.QPMIN:
                # -----------------------------------------------------------------------
                # sjl added (following lin eq. 23) the temperature dependency
                # to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
                # -----------------------------------------------------------------------
                factor = (
                    dts
                    * density_factor
                    * csaci
                    * exp(0.05 * tc + 0.8125 * log(internal_mixing_ratio_snow * density))
                )
                psaci = factor / (1.0 + factor) * internal_mixing_ratio_ice
            else:
                psaci = 0.0

            # -----------------------------------------------------------------------
            # psaut: autoconversion: cloud ice -- > snow
            # -----------------------------------------------------------------------

            # -----------------------------------------------------------------------
            # similar to lfo 1983: eq. 21 solved implicitly
            # threshold from wsm6 scheme, hong et al 2004, eq (13) : qi0_crt ~0.8e-4
            # -----------------------------------------------------------------------

            qim = ice_fraction(internal_t, convection_fraction, surface_type) * qi0_crt / density

            # -----------------------------------------------------------------------
            # assuming linear subgrid vertical distribution of cloud ice
            # the mismatch computation following lin et al. 1994, mwr
            # -----------------------------------------------------------------------

            if const_vi:
                tmp = fac_i2s
            else:
                tmp = fac_i2s * exp(0.025 * tc)

            di = max(di, constants.QCMIN)
            q_plus = internal_mixing_ratio_ice + di
            if q_plus > (qim + constants.QCMIN):
                if qim > (internal_mixing_ratio_ice - di):
                    dq = (0.25 * (q_plus - qim) ** 2) / di
                else:
                    dq = internal_mixing_ratio_ice - qim
                psaut = tmp * dq
            else:
                psaut = 0.0
            sink = min(internal_mixing_ratio_ice, psaci + psaut)

            # new total condensate / old condensate
            cloud_fraction = max(
                0.0,
                min(
                    1.0,
                    cloud_fraction
                    * max(internal_mixing_ratio_ice + internal_mixing_ratio_liquid - sink + tmp, 0.0)
                    / max(internal_mixing_ratio_ice + internal_mixing_ratio_liquid, constants.QCMIN),
                ),
            )

            internal_mixing_ratio_ice = internal_mixing_ratio_ice - sink
            internal_mixing_ratio_snow = internal_mixing_ratio_snow + sink

            # -----------------------------------------------------------------------
            # pgaci: accretion of cloud ice by graupel
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_graupel > constants.QPMIN:
                # -----------------------------------------------------------------------
                # factor = dts * cgaci / sqrt (den (k)) *
                # exp (0.05 * tc + 0.875 * log (qg * den (k)))
                # simplified form: remove temp dependency & set the exponent 0.875 -> 1
                # -----------------------------------------------------------------------
                factor = dts * cgaci * sqrt(density) * internal_mixing_ratio_graupel
                pgaci = factor / (1.0 + factor) * internal_mixing_ratio_ice
                internal_mixing_ratio_ice = internal_mixing_ratio_ice - pgaci
                internal_mixing_ratio_graupel = internal_mixing_ratio_graupel + pgaci

        # -----------------------------------------------------------------------
        # cold - rain proc:
        # -----------------------------------------------------------------------

        # -----------------------------------------------------------------------
        # rain to ice, snow, graupel processes:
        # -----------------------------------------------------------------------

        tc = internal_t - constants.TICE

        if internal_mixing_ratio_rain > constants.QPMIN and tc < 0.0:
            # -----------------------------------------------------------------------
            # * sink * terms to qr: psacr + pgfr
            # source terms to mixing_ratio_snow: psacr
            # source terms to mixing_ratio_graupel: pgfr
            # -----------------------------------------------------------------------

            # -----------------------------------------------------------------------
            # psacr accretion of rain by snow
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_snow > constants.QPMIN:  # if snow exists
                psacr = dts * acr3d(
                    terminal_fall_snow,
                    terminal_fall_rain,
                    internal_mixing_ratio_rain,
                    internal_mixing_ratio_snow,
                    constants.CSARC,
                    constants.ACCO_01,
                    constants.ACCO_11,
                    constants.ACCO_21,
                    density,
                )
            else:
                psacr = 0.0

            # -----------------------------------------------------------------------
            # pgfr: rain freezing -- > graupel
            # -----------------------------------------------------------------------

            pgfr = (
                dts
                * cgfr_0
                / density
                * (exp(-cgfr_1 * tc) - 1.0)
                * exp(1.75 * log(internal_mixing_ratio_rain * density))
            )

            # -----------------------------------------------------------------------
            # total sink to qr
            # -----------------------------------------------------------------------

            sink = psacr + pgfr
            factor = min(sink, min(internal_mixing_ratio_rain, -tc / icpk)) / max(sink, constants.QPMIN)

            psacr = factor * psacr
            pgfr = factor * pgfr

            sink = psacr + pgfr
            internal_mixing_ratio_rain = internal_mixing_ratio_rain - sink
            internal_mixing_ratio_snow = internal_mixing_ratio_snow + psacr
            internal_mixing_ratio_graupel = internal_mixing_ratio_graupel + pgfr
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = (
                c_air
                + internal_mixing_ratio_vapor * c_vap
                + q_liq * constants.C_LIQ
                + q_sol * constants.C_ICE
            )
            internal_t = internal_t + sink * lhi / cvm

        # # -----------------------------------------------------------------------
        # # update capacity heat and latent heat coefficient
        # # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * internal_t
        icpk = lhi / cvm

        # # -----------------------------------------------------------------------
        # # graupel production terms:
        # # -----------------------------------------------------------------------

        if internal_mixing_ratio_snow > constants.QPMIN:
            # -----------------------------------------------------------------------
            # accretion: snow -- > graupel
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_graupel > constants.QPMIN:
                sink = dts * acr3d(
                    terminal_fall_graupel,
                    terminal_fall_snow,
                    internal_mixing_ratio_snow,
                    internal_mixing_ratio_graupel,
                    cgacs,
                    constants.ACCO_03,
                    constants.ACCO_13,
                    constants.ACCO_23,
                    density,
                )
            else:
                sink = 0.0

            # -----------------------------------------------------------------------
            # autoconversion snow -- > graupel
            # -----------------------------------------------------------------------

            qsm = qs0_crt / density
            if internal_mixing_ratio_snow > qsm:
                factor = dts * 1.0e-3 * exp(0.09 * (internal_t - constants.TICE))
                sink = sink + factor / (1.0 + factor) * (internal_mixing_ratio_snow - qsm)
            sink = min(internal_mixing_ratio_snow, sink)

            # snow existed
            internal_mixing_ratio_snow = internal_mixing_ratio_snow - sink
            internal_mixing_ratio_graupel = internal_mixing_ratio_graupel + sink

        if internal_mixing_ratio_graupel > constants.QPMIN and internal_t < (constants.TICE - 0.01):
            # -----------------------------------------------------------------------
            # pgacw: accretion of cloud water by graupel
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_liquid > constants.QCMIN:
                qden = internal_mixing_ratio_graupel * density
                factor = dts * cgacw * qden / sqrt(density * sqrt(sqrt(qden)))
                pgacw = factor / (1.0 + factor) * internal_mixing_ratio_liquid
            else:
                pgacw = 0.0

            # -----------------------------------------------------------------------
            # pgacr: accretion of rain by graupel
            # -----------------------------------------------------------------------

            if internal_mixing_ratio_rain > constants.QPMIN:
                pgacr = min(
                    dts
                    * acr3d(
                        terminal_fall_graupel,
                        terminal_fall_rain,
                        internal_mixing_ratio_rain,
                        internal_mixing_ratio_graupel,
                        constants.CGARC,
                        constants.ACCO_02,
                        constants.ACCO_12,
                        constants.ACCO_22,
                        density,
                    ),
                    internal_mixing_ratio_rain,
                )
            else:
                pgacr = 0.0

            sink = pgacr + pgacw
            if constants.TICE - internal_t > 0:
                ans = constants.TICE - internal_t
            else:
                ans = 0
            factor = min(sink, ans / icpk) / max(sink, constants.QPMIN)
            pgacr = factor * pgacr
            pgacw = factor * pgacw

            sink = pgacr + pgacw
            internal_mixing_ratio_graupel = internal_mixing_ratio_graupel + sink
            internal_mixing_ratio_rain = internal_mixing_ratio_rain - pgacr
            internal_mixing_ratio_liquid = internal_mixing_ratio_liquid - pgacw
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = (
                c_air
                + internal_mixing_ratio_vapor * c_vap
                + q_liq * constants.C_LIQ
                + q_sol * constants.C_ICE
            )
            internal_t = internal_t + sink * lhi / cvm

    t = internal_t
    mixing_ratio_vapor = internal_mixing_ratio_vapor
    mixing_ratio_liquid = internal_mixing_ratio_liquid
    mixing_ratio_ice = internal_mixing_ratio_ice
    mixing_ratio_rain = internal_mixing_ratio_rain
    mixing_ratio_snow = internal_mixing_ratio_snow
    mixing_ratio_graupel = internal_mixing_ratio_graupel

    return (
        t,
        mixing_ratio_vapor,
        mixing_ratio_liquid,
        mixing_ratio_ice,
        mixing_ratio_rain,
        mixing_ratio_snow,
        mixing_ratio_graupel,
        cloud_fraction,
        p_dry,
        density,
        density_factor,
        terminal_fall_rain,
        terminal_fall_snow,
        terminal_fall_graupel,
        cvm,
        lhl,
        lhi,
        lcpk,
        icpk,
        tcpk,
    )


@function
def iqs1(
    t: Float,
    density: Float,
    table3: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
):
    """
    Compute saturation specific humidity from table3

    water - ice phase; universal dry / moist formular using air density
    "density" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function iqs1

    Args:
        t (Float)
        density (Float)
        table3 (GlobalTable_driver_qsat)
        des3 (GlobalTable_driver_qsat)
    """
    tmin = constants.TABLE_ICE - 160.0
    if t - tmin > 0:
        ans = t - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = int32(trunc(ap1))
    es = table3.A[it - 1] + (ap1 - it) * des3.A[it - 1]
    iqs1 = es / (constants.RVGAS * t * density)

    return iqs1


@function
def iqs2(
    t: Float,
    density: Float,
    table3: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
):
    """
    Compute saturation specific humidity from table3
    with additional calculation of gradient (dq/dt)

    water - ice phase; universal dry / moist formular using air density
    "density" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function iqs2

    Args:
        t (Float)
        density (Float)
        table3 (GlobalTable_driver_qsat)
        des3 (GlobalTable_driver_qsat)
    """
    tmin = constants.TABLE_ICE - 160.0
    if t - tmin > 0:
        ans = t - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = int32(trunc(ap1))
    es = table3.A[it - 1] + (ap1 - it) * des3.A[it - 1]
    iqs2 = es / (constants.RVGAS * t * density)
    it = int32(trunc(ap1 - 0.5))
    dqdt = (
        10.0 * (des3.A[it - 1] + (ap1 - it) * (des3.A[it] - des3.A[it - 1])) / (constants.RVGAS * t * density)
    )

    return iqs2, dqdt


@function
def wqs1(
    t: Float,
    density: Float,
    table2: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
):
    """Compute the saturated specific humidity for table2

    pure water phase; universal dry / moist formular using air density
    "density" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function wqs1

    Args:
        t (Float)
        density (Float)
        table2 (GlobalTable_driver_qsat)
        des2 (GlobalTable_driver_qsat)
    """
    tmin = constants.TABLE_ICE - 160.0
    if t - tmin > 0:
        ans = t - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = int32(trunc(ap1))
    es = table2.A[it - 1] + (ap1 - it) * des2.A[it - 1]
    wqs1 = es / (constants.RVGAS * t * density)

    return wqs1


@function
def subgrid_z_proc(
    p_dry: Float,
    density: Float,
    density_factor: Float,
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_rain: Float,
    mixing_ratio_ice: Float,
    mixing_ratio_snow: Float,
    mixing_ratio_graupel: Float,
    cloud_fraction: Float,
    rh_limited: Float,
    ccn: Float,
    convection_fraction: Float,
    surface_type: Float,
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
    """Temperature sensitive high vertical resolution processes

    reference Fortran: gfdl_cloud_microphys.F90: subroutine subgrid_z_proc

    Args:
        p_dry (Float)
        density (Float)
        density_factor (Float)
        t (Float)
        mixing_ratio_vapor (Float)
        mixing_ratio_liquid (Float)
        mixing_ratio_rain (Float)
        mixing_ratio_ice (Float)
        mixing_ratio_snow (Float)
        mixing_ratio_graupel (Float)
        cloud_fraction (Float)
        rh_limited (Float)
        ccn (Float)
        convection_fraction (Float)
        surface_type (Float)
        table2 (GlobalTable_driver_qsat)
        table3 (GlobalTable_driver_qsat)
        des2 (GlobalTable_driver_qsat)
        des3 (GlobalTable_driver_qsat)
        c_air (Float)
        c_vap (Float)
        cssub_0 (Float)
        cssub_1 (Float)
        cssub_2 (Float)
        cssub_3 (Float)
        cssub_4 (Float)
        d0_vap (Float)
        do_bigg (bool)
        do_evap (bool)
        do_qa (bool)
        dts (Float)
        fac_frz (Float)
        fac_g2v (Float)
        fac_l2v (Float)
        fac_s2v (Float)
        fac_v2g (Float)
        fac_v2s (Float)
        icloud_f (Float)
        lat2 (Float)
        lv00 (Float)
        preciprad (bool)
        qc_crt (Float)
        qi_lim (Float)
        rh_inc (Float)
        rh_inr (Float)
        t_min (Float)
        t_sub (Float)

    Returns:
        Float: t
        Float: mixing_ratio_vapor
        Float: mixing_ratio_liquid
        Float: mixing_ratio_rain
        Float: mixing_ratio_ice
        Float: mixing_ratio_snow
        Float: mixing_ratio_graupel
        Float: cloud_fraction
        Float: subl1
    """
    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    lhl = lv00 + d0_vap * t
    lhi = constants.LI00 + constants.DC_ICE * t
    q_liq = mixing_ratio_liquid + mixing_ratio_rain
    q_sol = mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel
    cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
    lcpk = lhl / cvm
    icpk = lhi / cvm
    tcpk = lcpk + icpk
    if constants.TICE - t > 0:
        ans = constants.TICE - t
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

    if t < t_min and cycle == False:  # noqa
        if mixing_ratio_vapor - constants.QVMIN > 0:
            sink = mixing_ratio_vapor - constants.QVMIN
        else:
            sink = 0
        mixing_ratio_vapor = mixing_ratio_vapor - sink
        mixing_ratio_ice = mixing_ratio_ice + sink
        q_sol = q_sol + sink
        cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
        t = t + sink * (lhl + lhi) / cvm
        if do_qa == True:  # noqa
            cloud_fraction = 1.0  # air fully saturated; 100 % cloud cover
        cycle = True

    if cycle == False:  # noqa
        # -----------------------------------------------------------------------
        # update heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------
        lhl = lv00 + d0_vap * t
        lhi = constants.LI00 + constants.DC_ICE * t
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk
        if constants.TICE - t > 0:
            ans = constants.TICE - t
        else:
            ans = 0
        tcp3 = lcpk + icpk * min(1.0, ans / (constants.TICE - constants.T_WFR))

        # -----------------------------------------------------------------------
        # instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
        # -----------------------------------------------------------------------
        qpz = mixing_ratio_vapor + mixing_ratio_liquid + mixing_ratio_ice
        tin = t - (lhl * (mixing_ratio_liquid + mixing_ratio_ice) + lhi * mixing_ratio_ice) / (
            c_air
            + qpz * c_vap
            + mixing_ratio_rain * constants.C_LIQ
            + (mixing_ratio_snow + mixing_ratio_graupel) * constants.C_ICE
        )
        if tin > t_sub + 6.0:
            rh = qpz / iqs1(tin, density, table3, des3)
            if rh < rh_adj:  # qpz / rh_adj < mixing_ratio_snow
                t = tin
                mixing_ratio_vapor = qpz
                mixing_ratio_liquid = 0.0
                mixing_ratio_ice = 0.0
                if do_qa == True:  # noqa
                    qa = 0.0
                cycle = True  # cloud free

    if cycle == False:  # noqa
        # -----------------------------------------------------------------------
        # cloud water < -- > vapor adjustment: LS evaporation
        # -----------------------------------------------------------------------
        if do_evap == True:  # noqa
            qsw, dwsdt = wqs2(t, density, table2, des2)
            dq0 = qsw - mixing_ratio_vapor
            if dq0 > constants.QVMIN:
                factor = min(1.0, fac_l2v * (10.0 * dq0 / qsw))
                evap = min(mixing_ratio_liquid, factor * mixing_ratio_liquid / (1.0 + tcp3 * dwsdt))
            else:
                evap = 0.0
            mixing_ratio_vapor = mixing_ratio_vapor + evap
            mixing_ratio_liquid = mixing_ratio_liquid - evap
            q_liq = q_liq - evap
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t - evap * lhl / cvm

        # -----------------------------------------------------------------------
        # update heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * t
        icpk = lhi / cvm

        # -----------------------------------------------------------------------
        # enforce complete freezing when ice_fraction==1
        # -----------------------------------------------------------------------

        ifrac = ice_fraction(t, convection_fraction, surface_type)
        if ifrac == 1.0 and mixing_ratio_liquid > constants.QCMIN:
            sink = mixing_ratio_liquid
            mixing_ratio_liquid = mixing_ratio_liquid - sink
            mixing_ratio_ice = mixing_ratio_ice + sink
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t + sink * lhi / cvm

        # -----------------------------------------------------------------------
        # update heat capacity and latent heat coefficient
        # -----------------------------------------------------------------------

        lhi = constants.LI00 + constants.DC_ICE * t
        icpk = lhi / cvm

        # -----------------------------------------------------------------------
        # bigg mechanism heterogeneous freezing on existing cloud nuclei
        # -----------------------------------------------------------------------

        tc = constants.TICE - t
        if do_bigg == True and mixing_ratio_liquid > constants.QCMIN and tc > 0.0:  # noqa
            sink = (
                fac_frz
                * (100.0 / constants.RHOR / ccn)
                * dts
                * (exp(0.66 * tc) - 1.0)
                * density
                * mixing_ratio_liquid
                * mixing_ratio_liquid
            )
            sink = min(mixing_ratio_liquid, min(tc / icpk, sink))
            mixing_ratio_liquid = mixing_ratio_liquid - sink
            mixing_ratio_ice = mixing_ratio_ice + sink
            q_liq = q_liq - sink
            q_sol = q_sol + sink
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t + sink * lhi / cvm
            # significant mixing_ratio_liquid existed

        # -----------------------------------------------------------------------
        # update capacity heat and latent heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t
        lhi = constants.LI00 + constants.DC_ICE * t
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # sublimation / deposition of LS ice
        # -----------------------------------------------------------------------

        if t < constants.TICE:
            qsi, dqsdt = iqs2(t, density, table3, des3)
            dq = mixing_ratio_vapor - qsi
            sink = min(mixing_ratio_ice, dq / (1.0 + tcpk * dqsdt))
            if mixing_ratio_ice > constants.QCMIN:
                # eq 9, hong et al. 2004, mwr
                # for a and b, see dudhia 1989: page 3103 eq (b7) and (b8)
                pidep = (
                    dts
                    * dq
                    * 349138.78
                    * exp(0.875 * log(mixing_ratio_ice * density))
                    / (qsi * density * lat2 / (0.0243 * constants.RVGAS * t**2) + 4.42478e4)
                )
            else:
                pidep = 0.0
            if dq > 0.0:  # vapor - > ice
                # deposition
                ifrac = ice_fraction(t, convection_fraction, surface_type)
                tmp = constants.TICE - t
                qi_crt = 4.92e-11 * exp(1.33 * log(1.0e3 * exp(0.1 * tmp)))
                qi_crt = max(qi_crt, 1.82e-6) * qi_lim * ifrac / density
                sink = min(sink, min(max(qi_crt - mixing_ratio_ice, pidep), tmp / tcpk))
            else:  # ice -- > vapor
                # NOTE sublimation not implemented
                # trigger checked in driver `check_flags` function

                # dev NOTE: unsure how to handle pssub. In Fortran this variable is
                # initialized to nan then used here (at least when do_subl is False,
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
            mixing_ratio_vapor = mixing_ratio_vapor - sink
            mixing_ratio_ice = mixing_ratio_ice + sink
            q_sol = q_sol + sink
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t + sink * (lhl + lhi) / cvm

        # -----------------------------------------------------------------------
        # update capacity heat and latend heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t
        lhi = constants.LI00 + constants.DC_ICE * t
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # sublimation / deposition of snow
        # this process happens for all temp rage
        # -----------------------------------------------------------------------

        if mixing_ratio_snow > constants.QPMIN:
            qsi, dqsdt = iqs2(t, density, table3, des3)
            qden = mixing_ratio_snow * density
            tmp = exp(0.65625 * log(qden))
            tsq = t * t
            dq = (qsi - mixing_ratio_vapor) / (1.0 + tcpk * dqsdt)
            pssub = (
                cssub_0
                * tsq
                * (cssub_1 * sqrt(qden) + cssub_2 * tmp * sqrt(density_factor))
                / (cssub_3 * tsq + cssub_4 * qsi * density)
            )
            pssub = (qsi - mixing_ratio_vapor) * dts * pssub
            if pssub > 0.0:  # snow -- > vapor, sublimation
                if t - t_sub > 0:
                    ans = t - t_sub
                else:
                    ans = 0
                pssub = min(fac_s2v * pssub * min(1.0, ans * 0.2), mixing_ratio_snow)
                subl1 = subl1 + pssub / dts
            else:
                if t > constants.TICE:
                    pssub = 0.0  # no deposition
                else:
                    pssub = max(fac_v2s * pssub, max(dq, (t - constants.TICE) / tcpk))
            mixing_ratio_snow = mixing_ratio_snow - pssub
            mixing_ratio_vapor = mixing_ratio_vapor + pssub
            q_sol = q_sol - pssub
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t - pssub * (lhl + lhi) / cvm

        # -----------------------------------------------------------------------
        # update capacity heat and latend heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t
        lhi = constants.LI00 + constants.DC_ICE * t
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # -----------------------------------------------------------------------
        # simplified 2 - way grapuel sublimation - deposition mechanism
        # -----------------------------------------------------------------------

        if mixing_ratio_graupel > constants.QPMIN:
            qsi, dqsdt = iqs2(t, density, table3, des3)
            dq = (mixing_ratio_vapor - qsi) / (1.0 + tcpk * dqsdt)
            pgsub = (mixing_ratio_vapor / qsi - 1.0) * mixing_ratio_graupel
            if pgsub > 0.0:  # deposition
                if t > constants.TICE:
                    pgsub = 0.0  # no deposition
                else:
                    pgsub = min(
                        fac_v2g * pgsub,
                        min(
                            0.2 * dq,
                            min(mixing_ratio_liquid + mixing_ratio_rain, (constants.TICE - t) / tcpk),
                        ),
                    )
            else:  # submilation
                if t - t_sub > 0:
                    ans = t - t_sub
                else:
                    ans = 0
                pgsub = max(fac_g2v * pgsub, dq) * min(1.0, ans * 0.1)
                subl1 = subl1 + pgsub / dts
            mixing_ratio_graupel = mixing_ratio_graupel + pgsub
            mixing_ratio_vapor = mixing_ratio_vapor - pgsub
            q_sol = q_sol + pgsub
            cvm = c_air + mixing_ratio_vapor * c_vap + q_liq * constants.C_LIQ + q_sol * constants.C_ICE
            t = t + pgsub * (lhl + lhi) / cvm

            # Fortran ifdef USE_MIN_EVAP goes in here.
            # Not currently executed, so not included

        # -----------------------------------------------------------------------
        # update capacity heat and latend heat coefficient
        # -----------------------------------------------------------------------

        lhl = lv00 + d0_vap * t
        cvm = c_air + (mixing_ratio_vapor + q_liq + q_sol) * c_vap
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
            q_sol = mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel
            q_liq = mixing_ratio_liquid + mixing_ratio_rain
        else:
            q_sol = mixing_ratio_ice
            q_liq = mixing_ratio_liquid
        q_cond = q_liq + q_sol

        qpz = mixing_ratio_vapor + q_cond  # qpz is conserved

        # -----------------------------------------------------------------------
        # use the "liquid - frozen water temperature" (tin)
        # to compute saturated specific humidity
        # -----------------------------------------------------------------------

        tin = t - (lcpk * q_cond + icpk * q_sol)  # minimum temperature

        # -----------------------------------------------------------------------
        # determine saturated specific humidity
        # -----------------------------------------------------------------------

        if tin <= constants.T_WFR:
            # ice phase:
            qstar = iqs1(tin, density, table3, des3)
        elif tin >= constants.TICE:
            # liquid phase:
            qstar = wqs1(tin, density, table2, des2)
        else:
            # mixed phase:
            qsi = iqs1(tin, density, table3, des3)
            qsw = wqs1(tin, density, table2, des2)
            if q_cond > 3.0e-6:
                rqi = q_sol / q_cond
            else:
                # WMP impose CALIPSO ice polynomial from 0 C to -40 C
                rqi = ice_fraction(tin, convection_fraction, surface_type)
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
                    cloud_fraction = max(
                        constants.QCMIN,
                        min(
                            1.0,
                            cloud_fraction
                            + (q_plus - qstar) * (q_plus - qstar) / ((q_plus - q_minus) * (q_plus - qpz)),
                        ),
                    )
                elif q_minus <= qstar and qstar < qpz:  # partial cloud cover
                    cloud_fraction = max(
                        constants.QCMIN,
                        min(
                            1.0,
                            cloud_fraction
                            + 1.0
                            - (
                                (qstar - q_minus) * (qstar - q_minus) / ((q_plus - q_minus) * (qpz - q_minus))
                            ),
                        ),
                    )
                elif qstar <= q_minus:
                    cloud_fraction = 1.0  # air fully saturated; 100 % cloud cover
            else:
                # top-hat
                if q_plus <= qstar:
                    # little/no cloud cover
                    do_nothing = True
                elif qstar < q_plus and q_cond > qc_crt:
                    cloud_fraction = max(
                        constants.QCMIN,
                        min(1.0, cloud_fraction + (q_plus - qstar) / (dq + dq)),
                    )  # partial cloud cover
                elif qstar <= q_minus:
                    cloud_fraction = 1.0  # air fully saturated; 100 % cloud cover

    return (
        t,
        mixing_ratio_vapor,
        mixing_ratio_liquid,
        mixing_ratio_rain,
        mixing_ratio_ice,
        mixing_ratio_snow,
        mixing_ratio_graupel,
        cloud_fraction,
        subl1,
    )


def icloud_core(
    t: FloatField,
    p_dry: FloatField,
    dp: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_liquid: FloatField,
    mixing_ratio_rain: FloatField,
    mixing_ratio_ice: FloatField,
    mixing_ratio_snow: FloatField,
    mixing_ratio_graupel: FloatField,
    cloud_fraction: FloatField,
    density: FloatField,
    density_factor: FloatField,
    terminal_fall_snow: FloatField,
    terminal_fall_graupel: FloatField,
    terminal_fall_rain: FloatField,
    sublimation: FloatField,
    rh_limited: FloatField,
    ccn: FloatField,
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    table2: GlobalTable_driver_qsat,
    table3: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
    des3: GlobalTable_driver_qsat,
):
    """sources of cloud ice: pihom, cold rain, and the sat_adj
    (initiation plus deposition)
    sources of snow: cold rain, auto conversion + accretion (from cloud ice)
    sat_adj (deposition; requires pre - existing snow);
    initial snow comes from auto conversion

    Args:
        t (FloatField)
        p_dry (FloatField)
        dp (FloatField)
        mixing_ratio_vapor (FloatField)
        mixing_ratio_liquid (FloatField)
        mixing_ratio_rain (FloatField)
        mixing_ratio_ice (FloatField)
        mixing_ratio_snow (FloatField)
        mixing_ratio_graupel (FloatField)
        cloud_fraction (FloatField)
        density (FloatField)
        density_factor (FloatField)
        terminal_fall_snow (FloatField)
        terminal_fall_graupel (FloatField)
        terminal_fall_rain (FloatField)
        sublimation (FloatField)
        rh_limited (FloatField)
        ccn (FloatField)
        convection_fraction (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        table2 (GlobalTable_driver_qsat)
        table3 (GlobalTable_driver_qsat)
        des2 (GlobalTable_driver_qsat)
        des3 (GlobalTable_driver_qsat)
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

    with computation(PARALLEL), interval(...):
        (
            t,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            cvm,
            q_liq,
            q_sol,
        ) = icloud_melt_freeze(
            t,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            density,
            convection_fraction,
            surface_type,
            c_air,
            c_vap,
            fac_frz,
            fac_imlt,
            qi0_crt,
            ql_mlt,
        )

    with computation(PARALLEL), interval(...):
        if z_slope_ice == True:  # noqa
            q_linear_prof = mixing_ratio_ice
            h_var_linear_prof = rh_limited
            dm_linear_prof = q_linear_prof  # initialized here to ensure it is created as a 3d field

    with computation(FORWARD), interval(1, None):
        if z_slope_ice == True:  # noqa
            dq_linear_prof = 0.5 * (q_linear_prof - q_linear_prof[0, 0, -1])

    # use twice the strength of the positive definiteness limiter (lin et al 1994)
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

    # impose a presumed background horizontal variability
    # that is proportional to the value itself
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

    with computation(PARALLEL), interval(...):
        # update capacity heat and latent heat coefficient
        lhl = lv00 + d0_vap * t
        lhi = constants.LI00 + constants.DC_ICE * t
        lcpk = lhl / cvm
        icpk = lhi / cvm
        tcpk = lcpk + icpk

        # do nothing above p_min
        if p_dry >= constants.P_MIN:
            (
                t,
                mixing_ratio_vapor,
                mixing_ratio_liquid,
                mixing_ratio_ice,
                mixing_ratio_rain,
                mixing_ratio_snow,
                mixing_ratio_graupel,
                cloud_fraction,
                p_dry,
                density,
                density_factor,
                terminal_fall_rain,
                terminal_fall_snow,
                terminal_fall_graupel,
                cvm,
                lhl,
                lhi,
                lcpk,
                icpk,
                tcpk,
            ) = snow_graupel_coldrain(
                t=t,
                mixing_ratio_vapor=mixing_ratio_vapor,
                mixing_ratio_liquid=mixing_ratio_liquid,
                mixing_ratio_ice=mixing_ratio_ice,
                mixing_ratio_rain=mixing_ratio_rain,
                mixing_ratio_snow=mixing_ratio_snow,
                mixing_ratio_graupel=mixing_ratio_graupel,
                cloud_fraction=cloud_fraction,
                p_dry=p_dry,
                density=density,
                density_factor=density_factor,
                terminal_fall_rain=terminal_fall_rain,
                terminal_fall_snow=terminal_fall_snow,
                terminal_fall_graupel=terminal_fall_graupel,
                cvm=cvm,
                lhl=lhl,
                lhi=lhi,
                lcpk=lcpk,
                icpk=icpk,
                tcpk=tcpk,
                q_liq=q_liq,
                q_sol=q_sol,
                di=di,
                convection_fraction=convection_fraction,
                surface_type=surface_type,
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

        (
            t,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            sublimation,
        ) = subgrid_z_proc(
            p_dry=p_dry,
            density=density,
            density_factor=density_factor,
            t=t,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_liquid=mixing_ratio_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_ice=mixing_ratio_ice,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_graupel=mixing_ratio_graupel,
            cloud_fraction=cloud_fraction,
            rh_limited=rh_limited,
            ccn=ccn,
            convection_fraction=convection_fraction,
            surface_type=surface_type,
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


def update_precip_total(
    sublimation: FloatField,
    driver_sublimation: FloatField,
):
    """
    ensure information is passed back to the rest of the model
    """
    with computation(PARALLEL), interval(...):
        sublimation = sublimation + driver_sublimation

        driver_sublimation = 0


class GFDL1MIceCloud(NDSLRuntime):
    """
    Ice cloud microphysics processes
    bulk cloud micro - physics; processes splitting
    with some un - split sub - grouping
    time implicit (when possible) accretion and autoconversion
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        config_dependent_constants: GFDL1MDriverConfigDependentConstants,
        saturation_tables: GFDL_driver_tables,
    ):
        """Initialize the ice cloud module

        Args:
            stencil_factory (StencilFactory)
            quantity_factory (QuantityFactory)
            config (GFDL1MConfig)
            config_dependent_constants (GFDL1MDriverConfigDependentConstants)
            saturation_tables (GFDL_driver_tables)
        """
        # initialize NDSLRuntime
        super().__init__(stencil_factory)

        # make saturation tables visible at runtime
        self.saturation_tables = saturation_tables

        # initialize locals
        self._sublimation = self.make_local(quantity_factory, [I_DIM, J_DIM, K_DIM], Float)

        # construct stencils
        self._icloud_core = stencil_factory.from_dims_halo(
            func=icloud_core,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "c_air": config_dependent_constants.C_AIR,
                "c_vap": config_dependent_constants.C_VAP,
                "dts": config_dependent_constants.DTS,
                "rdts": config_dependent_constants.RDTS,
                "const_vi": config.CONST_VI,
                "fac_g2v": config_dependent_constants.FAC_G2V,
                "fac_i2s": config_dependent_constants.FAC_I2S,
                "fac_imlt": config_dependent_constants.FAC_IMLT,
                "fac_frz": config_dependent_constants.FAC_FRZ,
                "fac_l2v": config_dependent_constants.FAC_L2V,
                "fac_s2v": config_dependent_constants.FAC_S2V,
                "fac_v2s": config_dependent_constants.FAC_V2S,
                "fac_v2g": config_dependent_constants.FAC_V2G,
                "cgacs": config_dependent_constants.CGACS,
                "csacw": config_dependent_constants.CSACW,
                "csaci": config_dependent_constants.CSACI,
                "cgacw": config_dependent_constants.CGACW,
                "cgaci": config_dependent_constants.CGACI,
                "cgfr_0": config_dependent_constants.CGFR_0,
                "cgfr_1": config_dependent_constants.CGFR_1,
                "csmlt_0": config_dependent_constants.CSMLT_0,
                "csmlt_1": config_dependent_constants.CSMLT_1,
                "csmlt_2": config_dependent_constants.CSMLT_2,
                "csmlt_3": config_dependent_constants.CSMLT_3,
                "csmlt_4": config_dependent_constants.CSMLT_4,
                "cgmlt_0": config_dependent_constants.CGMLT_0,
                "cgmlt_1": config_dependent_constants.CGMLT_1,
                "cgmlt_2": config_dependent_constants.CGMLT_2,
                "cgmlt_3": config_dependent_constants.CGMLT_3,
                "cgmlt_4": config_dependent_constants.CGMLT_4,
                "cssub_0": config_dependent_constants.CSSUB_0,
                "cssub_1": config_dependent_constants.CSSUB_1,
                "cssub_2": config_dependent_constants.CSSUB_2,
                "cssub_3": config_dependent_constants.CSSUB_3,
                "cssub_4": config_dependent_constants.CSSUB_4,
                "qi0_crt": config.QI0_CRT,
                "qs0_crt": config.QS0_CRT,
                "qs_mlt": config.QS_MLT,
                "ql_mlt": config.QL_MLT,
                "z_slope_ice": config.Z_SLOPE_ICE,
                "lv00": config_dependent_constants.LV00,
                "d0_vap": config_dependent_constants.D0_VAP,
                "lat2": config_dependent_constants.LAT2,
                "do_qa": config.DO_QA,
                "do_evap": config.DO_EVAP,
                "do_bigg": config.DO_BIGG,
                "qc_crt": config.QC_CRT,
                "qi_lim": config.QI_LIM,
                "rh_inc": config.RH_INC,
                "rh_inr": config.RH_INR,
                "t_min": config.T_MIN,
                "t_sub": config.T_SUB,
                "preciprad": config.PRECIPRAD,
                "icloud_f": config.ICLOUD_F,
            },
        )

        self._update_output = stencil_factory.from_dims_halo(
            func=update_precip_total,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

    def __call__(
        self,
        t: FloatField,
        p_dry: FloatField,
        dp: FloatField,
        mixing_ratio_vapor: FloatField,
        mixing_ratio_liquid: FloatField,
        mixing_ratio_rain: FloatField,
        mixing_ratio_ice: FloatField,
        mixing_ratio_snow: FloatField,
        mixing_ratio_graupel: FloatField,
        cloud_fraction: FloatField,
        density: FloatField,
        density_factor: FloatField,
        terminal_fall_snow: FloatField,
        terminal_fall_graupel: FloatField,
        terminal_fall_rain: FloatField,
        sublimation: FloatField,
        rh_limited: FloatField,
        ccn: FloatField,
        convection_fraction: FloatFieldIJ,
        surface_type: FloatFieldIJ,
    ):
        """
        Ice cloud microphysics processes
        bulk cloud micro - physics; processes splitting
        with some un - split sub - grouping
        time implicit (when possible) accretion and autoconversion

        Args:
            t (inout): temperature (K)
            p_dry (in): dry air pressure (Pa)
            dp (in): change in pressure between model layers (Pa)
            mixing_ratio_vapor (inout): mixing ratio vapor (kg/kg)
            mixing_ratio_liquid (inout): mixing ratio liquid (kg/kg)
            mixing_ratio_rain (inout): mixing ratio rain (kg/kg)
            mixing_ratio_ice (inout): mixing ratio ice (kg/kg)
            mixing_ratio_snow (inout): mixing ratio snow (kg/kg)
            mixing_ratio_graupel (inout): mixing ratio graupel (kg/kg)
            cloud_fraction (inout): cloud fraction
            density (in): density of the grid cell (kg m^-3)
            density_factor (in): details unknown
            terminal_fall_snow (in): terminal speed of snow (m/s)
            terminal_fall_graupel (in): terminal speed of graupel (m/s)
            terminal_fall_rain (in): terminal speed of rain (m/s)
            sublimation (out): model at-large sublimation (kg kg-1 s-1)
            rh_limited (in): relative humidity with limits imposed
            ccn (in): cloud condensation nuclei
            convection_fraction (in): convection fraction
            surface_type (in): surface type
        """
        self._icloud_core(
            t,
            p_dry,
            dp,
            mixing_ratio_vapor,
            mixing_ratio_liquid,
            mixing_ratio_rain,
            mixing_ratio_ice,
            mixing_ratio_snow,
            mixing_ratio_graupel,
            cloud_fraction,
            density,
            density_factor,
            terminal_fall_snow,
            terminal_fall_graupel,
            terminal_fall_rain,
            self._sublimation,
            rh_limited,
            ccn,
            convection_fraction,
            surface_type,
            self.saturation_tables.table2,
            self.saturation_tables.table3,
            self.saturation_tables.des2,
            self.saturation_tables.des3,
        )

        self._update_output(
            sublimation,
            self._sublimation,
        )
