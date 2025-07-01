"""Stencils and functions called by multiple pyMoist modules.
These functions evaluate various in-cloud microphysical
processes/quantities."""

import gt4py.cartesian.gtscript as gtscript
from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval, log10, sin

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_generic_math import air_density


@gtscript.function
def ice_fraction_modis(
    temp: Float,
):
    # Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
    tc = max(-46.0, min(temp - constants.MAPL_TICE, 46.0))  # convert to celcius and limit range from -46:46 C
    ptc = 7.6725 + 1.0118 * tc + 0.1422 * tc**2 + 0.0106 * tc**3 + 0.000339 * tc**4 + 0.00000395 * tc**5
    ice_frct = 1.0 - (1.0 / (1.0 + exp(-1 * ptc)))
    return ice_frct


@gtscript.function
def ice_fraction(
    temp: Float,
    cnv_frc: Float,
    srf_type: Float,
):
    # Anvil clouds
    # Anvil-Convective sigmoidal function like figure 6(right)
    # Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
    if temp <= constants.JaT_ICE_ALL:
        icefrct_c = 1.000
    elif temp > constants.JaT_ICE_ALL and temp <= constants.JaT_ICE_MAX:
        icefrct_c = sin(
            0.5
            * constants.MAPL_PI
            * (1.00 - (temp - constants.JaT_ICE_ALL) / (constants.JaT_ICE_MAX - constants.JaT_ICE_ALL))
        )
    else:
        icefrct_c = 0.00
    icefrct_c = max(min(icefrct_c, 1.00), 0.00) ** constants.aICEFRPWR
    # Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if srf_type == 2.0:
        if temp <= constants.JiT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.JiT_ICE_ALL and temp <= constants.JiT_ICE_MAX:
            icefrct_m = 1.00 - (temp - constants.JiT_ICE_ALL) / (
                constants.JiT_ICE_MAX - constants.JiT_ICE_ALL
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.iICEFRPWR
    elif srf_type > 1.0:
        if temp <= constants.lT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.lT_ICE_ALL and temp <= constants.lT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * constants.MAPL_PI
                * (1.00 - (temp - constants.lT_ICE_ALL) / (constants.lT_ICE_MAX - constants.lT_ICE_ALL))
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.lICEFRPWR
    else:
        if temp <= constants.oT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.oT_ICE_ALL and temp <= constants.oT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * constants.MAPL_PI
                * (1.00 - (temp - constants.oT_ICE_ALL) / (constants.oT_ICE_MAX - constants.oT_ICE_ALL))
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.oICEFRPWR
    ice_frac = icefrct_m * (1.0 - cnv_frc) + icefrct_c * cnv_frc
    return ice_frac


@gtscript.function
def cloud_effective_radius_liquid(
    pressure: Float,
    temperature: Float,
    liquid_mixing_ratio: Float,
    liquid_concentration: Float,
) -> Float:
    """
    Calculate the effective radius of liquid droplets clouds

    Arguments:
        pressure (in): pressure (millibars)
        temperature (in): temperature (Kelvin)
        liquid_mixing_ratio (in): liquid mixing ratio (kg/kg)
        liquid_concentration (in): liquid cloud droplet concentration (m^-3)

    Returns:
        radius (Float): drop radius
    """
    # Calculate liquid water content
    wc = (
        1.0e3 * air_density(pressure, temperature) * liquid_mixing_ratio
    )  # air density [g/m3] * liquid cloud mixing ratio [kg/kg]
    # Calculate cloud drop number concentration from the aerosol model + ....
    nnx = max(liquid_concentration * 1.0e-6, 10.0)
    # Calculate Radius in meters [m]
    if constants.LIQ_RADII_PARAM == 1:
        # Jason Version
        radius = min(
            60.0e-6,
            max(
                2.5e-6,
                1.0e-6 * constants.BX * (wc / nnx) ** constants.R13BBETA * constants.ABETA * 6.92,
            ),
        )
    else:
        # [liu&daum, 2000 and 2005. liu et al 2008]
        radius = min(
            60.0e-6,
            max(2.5e-6, 1.0e-6 * constants.LBX * (wc / nnx) ** constants.LBE),
        )
    return radius


@gtscript.function
def cloud_effective_radius_ice(
    pressure: Float,
    temperature: Float,
    ice_mixing_ratio: Float,
) -> Float:
    """
    Calculate the effective radius of ice particles in clouds

    Arguments:
        pressure (in): pressure (millibars)
        temperature (in): temperature (Kelvin)
        ice_mixing_ratio (in): liquid mixing ratio (kg/kg)

    Returns:
        radius (Float): ice particle radius
    """
    # Calculate ice water content
    wc = (
        1.0e3 * air_density(pressure, temperature) * ice_mixing_ratio
    )  # air density [g/m3] * ice cloud mixing ratio [kg/kg]
    # Calculate radius in meters [m]
    if constants.ICE_RADII_PARAM == 1:
        # Ice cloud effective radius -- [klaus wyser, 1998]
        if temperature > constants.MAPL_TICE or ice_mixing_ratio <= 0.0:
            bb = -2.0
        else:
            bb = -2.0 + log10(wc / 50.0) * (1.0e-3 * (constants.MAPL_TICE - temperature) ** 1.5)
            # NOTE: there is an issue in this line which causes differences between Fortran and Python
            # the multiplication "-2.0 * log'd result" is performed differently (~60 ULP), despite the log
            # being correct. Needs to be looked into at some point, but not critital for overall performance.
        bb = min(max(bb, -6.0), -2.0)
        radius = 377.4 + 203.3 * bb + 37.91 * bb**2 + 2.3696 * bb**3
        radius = min(150.0e-6, max(5.0e-6, 1.0e-6 * radius))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        tc = temperature - constants.MAPL_TICE
        zfsr = 1.2351 + 0.0105 * tc
        aa = 45.8966 * (wc**0.2214)
        bb = 0.79570 * (wc**0.2535)
        radius = zfsr * (aa + bb * (temperature - 83.15))
        radius = min(150.0e-6, max(5.0e-6, 1.0e-6 * radius * 0.64952))
    return radius


def fix_up_clouds(
    vapor: FloatField,
    t: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    large_scale_cloud_fraction: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    convective_cloud_fraction: FloatField,
) -> None:
    """
    Modify various cloud variables to ensure physical consistency.

    Performed in this order:
        If cloud fraction is too low, move all liquid and frozen water to vapor form
            and remove cloud.
        If liquid water is too low, move all liquid water to vapor form.
        If frozen water is too low, move all frozen water to vapor form.
        If total liquid + frozen water is too low, move all water to vapor form
            and remove cloud.

    Parameters:
    vapor (inout): water vapor mixing ratio
    t (inout): temperature.
    large_scale_liquid (inout): large scale cloud liquid water mixing ratio
    large_scale_ice (inout): large scale cloud frozen water mixing ratio
    large_scale_cloud_fraction (inout): large scale cloud fraction
    convective_liquid (inout): convective cloud liquid water mixing ratio
    convective_ice (inout): convective cloud frozen water mixing ratio
    convective_cloud_fraction (inout): convective cloud fraction
    """
    with computation(PARALLEL), interval(...):
        # fix small convective cloud fraction
        if convective_cloud_fraction < 1.0e-5:
            vapor = vapor + convective_liquid + convective_ice
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * convective_liquid
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * convective_ice
            )
            convective_cloud_fraction = 0.0
            convective_liquid = 0.0
            convective_ice = 0.0
        # fix small large scale cloud fraction
        if large_scale_cloud_fraction < 1.0e-5:
            vapor = vapor + large_scale_liquid + large_scale_ice
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * large_scale_liquid
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * large_scale_ice
            )
            large_scale_cloud_fraction = 0.0
            large_scale_liquid = 0.0
            large_scale_ice = 0.0
        # if large scale liquid water conentration is too low
        if large_scale_liquid < 1.0e-8:
            vapor = vapor + large_scale_liquid
            t = t - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * large_scale_liquid
            large_scale_liquid = 0.0
        # if large scale frozen water conentration is too low
        if large_scale_ice < 1.0e-8:
            vapor = vapor + large_scale_ice
            t = t - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * large_scale_ice
            large_scale_ice = 0.0
        # if convective liquid water conentration is too low
        if convective_liquid < 1.0e-8:
            vapor = vapor + convective_liquid
            t = t - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * convective_liquid
            convective_liquid = 0.0
        # if convective frozen water conentration is too low
        if convective_ice < 1.0e-8:
            vapor = vapor + convective_ice
            t = t - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * convective_ice
            convective_ice = 0.0
        # if total convective water is too low
        if (convective_liquid + convective_ice) < 1.0e-8:
            vapor = vapor + convective_liquid + convective_ice
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * convective_liquid
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * convective_ice
            )
            convective_cloud_fraction = 0.0
            convective_liquid = 0.0
            convective_ice = 0.0
        # if total large scale water is too low
        if (large_scale_liquid + large_scale_ice) < 1.0e-8:
            vapor = vapor + large_scale_liquid + large_scale_ice
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * large_scale_liquid
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * large_scale_ice
            )
            large_scale_cloud_fraction = 0.0
            large_scale_liquid = 0.0
            large_scale_ice = 0.0
