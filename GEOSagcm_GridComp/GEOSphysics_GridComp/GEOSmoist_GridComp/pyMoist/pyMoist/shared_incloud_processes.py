"""Stencils and functions called by multiple pyMoist modules.
These functions evaluate various in-cloud microphysical
processes/quantities."""

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import PARALLEL, computation, exp, interval, log10, sin

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_generic_math import air_density


@gtscript.function
def ice_fraction_modis(
    temp: Float,
):
    # Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
    tc = max(-46.0, min(temp - constants.MAPL_TICE, 46.0))  # convert to celcius and limit range from -46:46 C
    ptc = (
        7.6725
        + 1.0118 * tc
        + 0.1422 * tc**2
        + 0.0106 * tc**3
        + 0.000339 * tc**4
        + 0.00000395 * tc**5
    )
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
    PL: Float,
    TE: Float,
    QC: Float,
    NNL: Float,
    NNI: Float,
) -> Float:
    """
    Calculate the effective radius of liquid clouds [m]
    Implementation of LDRADIUS4 for liquid clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    QC (Float): Liquid cloud mixing ratio.
    NNL (Float): Number concentration of liquid cloud droplets.
    NNI (Float): Number concentration of ice cloud crystals. Not used in function body.

    Returns:
    Float: Effective radius of liquid clouds.
    """
    # Calculate liquid water content
    WC = 1.0e3 * air_density(PL, TE) * QC  # air density [g/m3] * liquid cloud mixing ratio [kg/kg]
    # Calculate cloud drop number concentration from the aerosol model + ....
    NNX = max(NNL * 1.0e-6, 10.0)
    # Calculate Radius in meters [m]
    if constants.LIQ_RADII_PARAM == 1:
        # Jason Version
        RADIUS = min(
            60.0e-6,
            max(
                2.5e-6,
                1.0e-6 * constants.BX * (WC / NNX) ** constants.R13BBETA * constants.ABETA * 6.92,
            ),
        )
    else:
        # [liu&daum, 2000 and 2005. liu et al 2008]
        RADIUS = min(
            60.0e-6,
            max(2.5e-6, 1.0e-6 * constants.LBX * (WC / NNX) ** constants.LBE),
        )
    return RADIUS


@gtscript.function
def cloud_effective_radius_ice(
    PL: Float,
    TE: Float,
    QC: Float,
    NNL: Float,
    NNI: Float,
) -> Float:
    """
    Calculate the effective radius of ice clouds [m]
    Implementation of LDRADIUS4 for Ice clouds

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.
    QC (Float): Ice cloud mixing ratio.
    NNL (Float): Number concentration of liquid cloud droplets.
        Not used in function body, but included in the Fortran code.
    NNI (Float): Number concentration of ice cloud crystals.
        Not used in function body, but included in the Fortran code.

    Returns:
    Float: Effective radius of ice clouds.
    """
    # Calculate ice water content
    WC = 1.0e3 * air_density(PL, TE) * QC  # air density [g/m3] * ice cloud mixing ratio [kg/kg]
    # Calculate radius in meters [m]
    if constants.ICE_RADII_PARAM == 1:
        # Ice cloud effective radius -- [klaus wyser, 1998]
        if TE > constants.MAPL_TICE or QC <= 0.0:
            BB = -2.0
        else:
            BB = -2.0 + log10(WC / 50.0) * (1.0e-3 * (constants.MAPL_TICE - TE) ** 1.5)
        BB = min(max(BB, -6.0), -2.0)
        RADIUS = 377.4 + 203.3 * BB + 37.91 * BB**2 + 2.3696 * BB**3
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        TC = TE - constants.MAPL_TICE
        ZFSR = 1.2351 + 0.0105 * TC
        AA = 45.8966 * (WC**0.2214)
        BB = 0.79570 * (WC**0.2535)
        RADIUS = ZFSR * (AA + BB * (TE - 83.15))
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    return RADIUS


def fix_up_clouds(
    q: FloatField,
    t: FloatField,
    qlls: FloatField,
    qils: FloatField,
    clls: FloatField,
    qlcn: FloatField,
    qicn: FloatField,
    clcn: FloatField,
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
    q (inout): water vapor mixing ratio
    t (inout): temperature.
    qlls (inout): large scale cloud liquid water mixing ratio
    qils (inout): large scale cloud frozen water mixing ratio
    clls (inout): large scale cloud fraction
    qlcn (inout): convective cloud liquid water mixing ratio
    qicn (inout): convective cloud frozen water mixing ratio
    clcn (inout): convective cloud fraction
    """
    with computation(PARALLEL), interval(...):
        # fix small convective cloud fraction
        if clcn < 1.0e-5:
            q = q + qlcn + qicn
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qlcn
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qicn
            )
            clcn = 0.0
            qlcn = 0.0
            qicn = 0.0
        # fix small large scale cloud fraction
        if clls < 1.0e-5:
            q = q + qlls + qils
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qlls
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * qils
            )
            clls = 0.0
            qlls = 0.0
            qils = 0.0
        # if large scale liquid water conentration is too low
        if qlls < 1.0e-8:
            q = q + qlls
            t = t - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qlls
            qlls = 0.0
        # if large scale frozen water conentration is too low
        if qils < 1.0e-8:
            q = q + qils
            t = t - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * qils
            qils = 0.0
        # if convective liquid water conentration is too low
        if qlcn < 1.0e-8:
            q = q + qlcn
            t = t - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qlcn
            qlcn = 0.0
        # if convective frozen water conentration is too low
        if qicn < 1.0e-8:
            q = q + qicn
            t = t - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * qicn
            qicn = 0.0
        # if total convective water is too low
        if (qlcn + qicn) < 1.0e-8:
            q = q + qlcn + qicn
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qlcn
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * qicn
            )
            clcn = 0.0
            qlcn = 0.0
            qicn = 0.0
        # if total large scale water is too low
        if (qlls + qils) < 1.0e-8:
            q = q + qlls + qils
            t = (
                t
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * qlls
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * qils
            )
            clls = 0.0
            qlls = 0.0
            qils = 0.0
