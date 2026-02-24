"""Stencils and functions called by multiple pyMoist modules.
These functions evaluate various in-cloud microphysical
processes/quantities."""

import pyMoist.constants as constants
from ndsl.dsl.gt4py import (
    PARALLEL,
    GlobalTable,
    computation,
    exp,
    float32,
    float64,
    floor,
    function,
    interval,
    log10,
    round_away_from_zero,
    sin,
)
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_generic_math import air_density


@function
def ice_fraction_modis(
    temp: Float,
):
    # Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
    tc = max(-46.0, min(temp - constants.MAPL_TICE, 46.0))  # convert to celcius and limit range from -46:46 C
    ptc = (
        7.6725 + 1.0118 * tc + 0.1422 * tc ** 2 + 0.0106 * tc ** 3 + 0.000339 * tc ** 4 + 0.00000395 * tc ** 5
    )
    ice_frct = 1.0 - (1.0 / (1.0 + exp(-1 * ptc)))
    return ice_frct


@function
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


@function
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


@function
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
        radius = 377.4 + 203.3 * bb + 37.91 * bb ** 2 + 2.3696 * bb ** 3
        radius = min(150.0e-6, max(5.0e-6, 1.0e-6 * radius))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        tc = temperature - constants.MAPL_TICE
        zfsr = 1.2351 + 0.0105 * tc
        aa = 45.8966 * (wc ** 0.2214)
        bb = 0.79570 * (wc ** 0.2535)
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


# able of lookup values of radiative effective radius of ice crystals as a function of temperature from
# -94C to 0C for make_ice_number. Taken from WRF RRTMG radiation code where it is attributed to
# Jon Egill Kristjansson and coauthors. This must be built into a custom shape off-grid quantity,
# passed into the stencil which calls make_ice_number with a custom field type (defined below), then
# passed into the function make_ice_number and accessed with .A[index] to function properly
RADIATIVE_EFFECTIVE_RADIUS = [
    5.92779,
    6.26422,
    6.61973,
    6.99539,
    7.39234,
    7.81177,
    8.25496,
    8.72323,
    9.21800,
    9.74075,
    10.2930,
    10.8765,
    11.4929,
    12.1440,
    12.8317,
    13.5581,
    14.2319,
    15.0351,
    15.8799,
    16.7674,
    17.6986,
    18.6744,
    19.6955,
    20.7623,
    21.8757,
    23.0364,
    24.2452,
    25.5034,
    26.8125,
    27.7895,
    28.6450,
    29.4167,
    30.1088,
    30.7306,
    31.2943,
    31.8151,
    32.3077,
    32.7870,
    33.2657,
    33.7540,
    34.2601,
    34.7892,
    35.3442,
    35.9255,
    36.5316,
    37.1602,
    37.8078,
    38.4720,
    39.1508,
    39.8442,
    40.5552,
    41.2912,
    42.0635,
    42.8876,
    43.7863,
    44.7853,
    45.9170,
    47.2165,
    48.7221,
    50.4710,
    52.4980,
    54.8315,
    57.4898,
    60.4785,
    63.7898,
    65.5604,
    71.2885,
    75.4113,
    79.7368,
    84.2351,
    88.8833,
    93.6658,
    98.5739,
    103.603,
    108.752,
    114.025,
    119.424,
    124.954,
    130.630,
    136.457,
    142.446,
    148.608,
    154.956,
    161.503,
    168.262,
    175.248,
    182.473,
    189.952,
    197.699,
    205.728,
    214.055,
    222.694,
    231.661,
    240.971,
    250.639,
]
RADIATIVE_EFFECTIVE_RADIUS_Table_Type = GlobalTable[(Float, len(RADIATIVE_EFFECTIVE_RADIUS))]


@function
def make_ice_number(
    cloud_ice_mixing_ratio,
    t,
    RADIATIVE_EFFECTIVE_RADIUS,
):
    """
    Get the ice crystal number given cloud ice mixing ratio and temperature.
    Returns the number of droplets per kg per m3.

    Args:
        cloud_ice_mixing_ratio (in): units kg/m3
        t (in): units K
        RADIATIVE_EFFECTIVE_RADIUS: table used for calculations

    Returns:
        crystal_number: units number/(kg*m3)

    Developed by H. Barnes @ NOAA/OAR/ESRL/GSL Earth Prediction Advancement Division
    """

    # DEBUG
    crystal_number = 0.0
    # internal constant
    ice_density = 890.0

    if cloud_ice_mixing_ratio == 0.0:
        crystal_number = 0.0

    else:
        # From the model 3D temperature field, subtract 180K for which
        # index value of RADIATIVE_EFFECTIVE_RADIUS as a start.  Value of corr is for
        # interpolating between neighboring values in the table.

        idx_rei = int(t - 180.0)
        idx_rei = min(max(idx_rei, 0), 93)
        corr = t - floor(t)
        reice = (
            RADIATIVE_EFFECTIVE_RADIUS.A[idx_rei] * (1.0 - corr)
            + RADIATIVE_EFFECTIVE_RADIUS.A[idx_rei + 1] * corr
        )
        deice = 2.0 * reice * 1.0e-6

        internal_lambda = float64(3.0 / deice)

        # value of the dispersion parameter according to Heymsfield et al 2002, Table3.
        t_celcius = t - 273.15

        t_celcius = min(max(t_celcius, -70.0), -15.0)

        if t_celcius > -27.0:
            lambdai = 6.8 * exp(-0.096 * t_celcius)
        else:
            lambdai = 24.8 * exp(-0.049 * t_celcius)

        mui = (0.13 * (lambdai ** 0.64)) - 2.0

        k = (mui + 3) * (mui * 3) / (mui + 2) / (mui + 1)

        crystal_number = (
            k
            * cloud_ice_mixing_ratio
            * internal_lambda
            * internal_lambda
            * internal_lambda
            / (constants.MAPL_PI * ice_density)
        )

    return crystal_number


# table of constants for make_droplet_number, this must be built into a custom shape off-grid quantity,
# passed into the stencil which calls make_droplet_number with a custom field type (defined below), then
# passed into the function make_droplet_number and accessed with .A[index] to function properly
G_RATIO = [24, 60, 120, 210, 336, 504, 720, 990, 1320, 1716, 2184, 2730, 3360, 4080, 4896]
G_RATIO_Table_Type = GlobalTable[(Float, len(G_RATIO))]


@function
def make_droplet_number(
    cloud_water_mixing_ratio,
    num_water_friendly_aerosols,
    G_RATIO,
):
    """
    Get the droplet number given cloud water mixing ratio and number of water-friendly aerosols.
    Returns the number of droplets per kg per m3.

    Args:
        cloud_water_mixing_ratio (in): units kg/m3
        num_water_friendly_aerosols (in): units number/kg
        G_RATIO: table used for calculations

    Returns:
        droplet_number: units number/(kg*m3)

    Developed by H. Barnes @ NOAA/OAR/ESRL/GSL Earth Prediction Advancement Division
    """
    am_r = constants.MAPL_PI * 1000.0 / 6.0

    if cloud_water_mixing_ratio <= 0.0:
        droplet_number = 0.0
    else:
        internal_num_water_friendly_aerosols = max(99.0e6, min(num_water_friendly_aerosols, 5.0e10))
        nu_c = max(2, min(round_away_from_zero(2.5e10 / internal_num_water_friendly_aerosols), 15))

        x1 = max(1.0, min(internal_num_water_friendly_aerosols * 1.0e-9, 10.0)) - 1.0
        xDc = (30.0 - x1 * 20.0 / 9.0) * 1.0e-6

        internal_lambda = (float64(4.0) + nu_c) / xDc

        qnc = (
            cloud_water_mixing_ratio
            / G_RATIO.A[int(nu_c - 1)]
            * internal_lambda
            * internal_lambda
            * internal_lambda
            / am_r
        )
        droplet_number = float32(qnc)

    return droplet_number
