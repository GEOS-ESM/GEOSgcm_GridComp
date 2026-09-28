from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, function, int32, interval, trunc
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ

from pyMoist.microphysics.GFDL_1M.driver.constants import constants
from pyMoist.microphysics.GFDL_1M.driver.sat_tables import GlobalTable_driver_qsat


@function
def fix_negative_core(
    t: Float,
    mixing_ratio_vapor: Float,
    mixing_ratio_liquid: Float,
    mixing_ratio_rain: Float,
    mixing_ratio_ice: Float,
    mixing_ratio_snow: Float,
    mixing_ratio_graupel: Float,
    c_air: Float,
    c_vap: Float,
    lv00: Float,
    d0_vap: Float,
):
    """Adjusts/removes negative mixing ratios

    reference Fortran: gfdl_cloud_microphys.F90: subroutine neg_adj

    Args:
        t (Float): temperature (Kelvin)
        mixing_ratio_vapor (Float): water vapor mixing ratio (kg/kg)
        mixing_ratio_liquid (Float): liquid water mixing ratio (kg/kg)
        mixing_ratio_rain (Float): rain mixing ratio (kg/kg)
        mixing_ratio_ice (Float): ice mixing ratio (kg/kg)
        mixing_ratio_snow (Float): snow mixing ratio (kg/kg)
        mixing_ratio_graupel (Float): graupel mixing ratio (kg/kg)
        c_air (Float)
        c_vap (Float)
        lv00 (Float)
        d0_vap (Float)

    Returns:
        Float: t
        Float: mixing_ratio_vapor
        Float: mixing_ratio_liquid
        Float: mixing_ratio_rain
        Float: mixing_ratio_ice
        Float: mixing_ratio_snow
        Float: mixing_ratio_graupel
    """
    # -----------------------------------------------------------------------
    # define heat capacity and latent heat coefficient
    # -----------------------------------------------------------------------

    cvm = (
        c_air
        + mixing_ratio_vapor * c_vap
        + (mixing_ratio_rain + mixing_ratio_liquid) * constants.C_LIQ
        + (mixing_ratio_ice + mixing_ratio_snow + mixing_ratio_graupel) * constants.C_ICE
    )
    lcpk = (lv00 + d0_vap * t) / cvm
    icpk = (constants.LI00 + constants.DC_ICE * t) / cvm

    # -----------------------------------------------------------------------
    # ice phase:
    # -----------------------------------------------------------------------

    # if cloud ice < 0, borrow from snow
    if mixing_ratio_ice < 0.0:
        mixing_ratio_snow = mixing_ratio_snow + mixing_ratio_ice
        mixing_ratio_ice = 0.0
    # if snow < 0, borrow from graupel
    if mixing_ratio_snow < 0.0:
        mixing_ratio_graupel = mixing_ratio_graupel + mixing_ratio_snow
        mixing_ratio_snow = 0.0
    # if graupel < 0, borrow from rain
    if mixing_ratio_graupel < 0.0:
        mixing_ratio_rain = mixing_ratio_rain + mixing_ratio_graupel
        t = t - mixing_ratio_graupel * icpk  # heating
        mixing_ratio_graupel = 0.0

    # -----------------------------------------------------------------------
    # liquid phase:
    # -----------------------------------------------------------------------

    # if rain < 0, borrow from cloud water
    if mixing_ratio_rain < 0.0:
        mixing_ratio_liquid = mixing_ratio_liquid + mixing_ratio_rain
        mixing_ratio_rain = 0.0
    # if cloud water < 0, borrow from water vapor
    if mixing_ratio_liquid < 0.0:
        mixing_ratio_vapor = mixing_ratio_vapor + mixing_ratio_liquid
        t = t - mixing_ratio_liquid * lcpk  # heating
        mixing_ratio_liquid = 0.0

    return (
        t,
        mixing_ratio_vapor,
        mixing_ratio_liquid,
        mixing_ratio_rain,
        mixing_ratio_ice,
        mixing_ratio_snow,
        mixing_ratio_graupel,
    )


def fix_negative_values(
    t: FloatField,
    dry_air_mixing_ratio_vapor: FloatField,
    dry_air_mixing_ratio_liquid: FloatField,
    dry_air_mixing_ratio_rain: FloatField,
    dry_air_mixing_ratio_ice: FloatField,
    dry_air_mixing_ratio_snow: FloatField,
    dry_air_mixing_ratio_graupel: FloatField,
    dp: FloatField,
):
    """Adjusts/removes negative mixing ratios and updates vapor and temperature
    where necessary. Core math is done in fix_negative_core

    reference Fortran: gfdl_cloud_microphys.F90: subroutine mpdrv

    Args:
        t (FloatField): temperature (Kelvin)
        dry_air_mixing_ratio_vapor (FloatField): water vapor mixing ratio (kg/kg)
        dry_air_mixing_ratio_liquid (FloatField): liquid water mixing ratio (kg/kg)
        dry_air_mixing_ratio_rain (FloatField): rain mixing ratio (kg/kg)
        dry_air_mixing_ratio_ice (FloatField): ice mixing ratio (kg/kg)
        dry_air_mixing_ratio_snow (FloatField): snow mixing ratio (kg/kg)
        dry_air_mixing_ratio_graupel (FloatField): graupel mixing ratio (kg/kg)
        dp (FloatField): pressure change between model layers (Pa)
    """
    from __externals__ import c_air, c_vap, d0_vap, lv00

    # -----------------------------------------------------------------------
    # fix all negative water species
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(0, -1):
        (
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
        ) = fix_negative_core(
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
            c_air,
            c_vap,
            lv00,
            d0_vap,
        )
        if dry_air_mixing_ratio_vapor < 0.0:
            dry_air_mixing_ratio_vapor[0, 0, 1] = dry_air_mixing_ratio_vapor[0, 0, 1] + dry_air_mixing_ratio_vapor * dp / dp[0, 0, 1]
            dry_air_mixing_ratio_vapor = 0.0

    with computation(FORWARD), interval(-1, None):
        (
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
        ) = fix_negative_core(
            t,
            dry_air_mixing_ratio_vapor,
            dry_air_mixing_ratio_liquid,
            dry_air_mixing_ratio_rain,
            dry_air_mixing_ratio_ice,
            dry_air_mixing_ratio_snow,
            dry_air_mixing_ratio_graupel,
            c_air,
            c_vap,
            lv00,
            d0_vap,
        )

        if dry_air_mixing_ratio_vapor < 0.0 and dry_air_mixing_ratio_vapor[0, 0, -1] > 0.0:
            dq = min(-dry_air_mixing_ratio_vapor * dp, dry_air_mixing_ratio_vapor[0, 0, -1] * dp[0, 0, -1])
            dry_air_mixing_ratio_vapor[0, 0, -1] = dry_air_mixing_ratio_vapor[0, 0, -1] - dq / dp[0, 0, -1]
            dry_air_mixing_ratio_vapor = dry_air_mixing_ratio_vapor + dq / dp


@function
def wqs2(
    ta: Float,
    den: Float,
    table2: GlobalTable_driver_qsat,
    des2: GlobalTable_driver_qsat,
):
    """
    Compute the saturated specific humidity for table2
    with additional calculation of gradient (dq/dt)

    pure water phase; universal dry / moist formula using air density
    input "den" can be either dry or moist air density

    reference Fortran: gfdl_cloud_microphys.F90: function wqs2
    """
    tmin = constants.TABLE_ICE - 160.0

    if ta - tmin > 0:
        ans = ta - tmin
    else:
        ans = 0
    ap1 = 10.0 * ans + 1.0
    ap1 = min(2621.0, ap1)
    it = int32(trunc(ap1))
    es = table2.A[it - 1] + (ap1 - it) * des2.A[it - 1]
    qsat = es / (constants.RVGAS * ta * den)
    it = int32(trunc(ap1 - 0.5))  # check if this rounds or truncates. need truncation here
    # finite diff, del_t = 0.1:
    dqdt = 10.0 * (des2.A[it - 1] + (ap1 - it) * (des2.A[it] - des2.A[it - 1])) / (constants.RVGAS * ta * den)

    return qsat, dqdt


def implicit_fall(
    mixing_ratio: FloatField,
    terminal_speed: FloatField,
    z_interface: FloatField,
    dp: FloatField,
    mass: FloatField,
    precip_flux: FloatField,
    precip: FloatFieldIJ,
    precip_fall: BoolFieldIJ,
):
    """
    Compute the time-implicit monotonic scheme
    """
    from __externals__ import dts

    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            height_diff = z_interface - z_interface[0, 0, 1]
            dd = dts * terminal_speed
            mixing_ratio = mixing_ratio * dp

    # sedimentation: non - vectorizable loop
    with computation(FORWARD), interval(0, 1):
        if precip_fall == True:  # noqa
            qm = mixing_ratio / (height_diff + dd)

    with computation(FORWARD), interval(1, None):
        if precip_fall == True:  # noqa
            qm = (mixing_ratio + dd[0, 0, -1] * qm[0, 0, -1]) / (height_diff + dd)

    # qm is density at this stage
    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            qm = qm * height_diff

    # output mass fluxes: non - vectorizable loop
    with computation(FORWARD), interval(0, 1):
        if precip_fall == True:  # noqa
            mass = mixing_ratio - qm

    with computation(FORWARD), interval(1, None):
        if precip_fall == True:  # noqa
            mass = mass[0, 0, -1] + mixing_ratio - qm

    with computation(FORWARD), interval(-1, None):
        if precip_fall == True:  # noqa
            precip = mass
        else:
            precip = 0

    # update:
    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            mixing_ratio = qm / dp

    with computation(PARALLEL), interval(...):
        if precip_fall == True:  # noqa
            precip_flux = precip_flux + mass
