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
        7.6725 + 1.0118 * tc + 0.1422 * tc ** 2 + 0.0106 * tc ** 3 + 0.000339 * tc ** 4 + 0.00000395 * tc ** 5
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
        RADIUS = 377.4 + 203.3 * BB + 37.91 * BB ** 2 + 2.3696 * BB ** 3
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS))
    else:
        # Ice cloud effective radius ----- [Sun, 2001]
        TC = TE - constants.MAPL_TICE
        ZFSR = 1.2351 + 0.0105 * TC
        AA = 45.8966 * (WC ** 0.2214)
        BB = 0.79570 * (WC ** 0.2535)
        RADIUS = ZFSR * (AA + BB * (TE - 83.15))
        RADIUS = min(150.0e-6, max(5.0e-6, 1.0e-6 * RADIUS * 0.64952))
    return RADIUS


def fix_up_clouds(
    QV: FloatField,
    TE: FloatField,
    QLC: FloatField,
    QIC: FloatField,
    CF: FloatField,
    QLA: FloatField,
    QIA: FloatField,
    AF: FloatField,
) -> None:
    """
    Fix cloud variables to ensure physical consistency.

    Parameters:
    QV (3D inout): Water vapor mixing ratio.
    TE (3D inout): Temperature.
    QLC (3D inout): Liquid cloud mixing ratio.
    QIC (3D inout): Ice cloud mixing ratio.
    CF (3D inout): Cloud fraction.
    QLA (3D inout): Anvil liquid cloud mixing ratio.
    QIA (3D inout): Anvil ice cloud mixing ratio.
    AF (3D inout): Anvil cloud fraction.
    """
    with computation(PARALLEL), interval(...):
        # Fix if Anvil cloud fraction too small
        if AF < 1.0e-5:
            QV = QV + QLA + QIA
            TE = (
                TE
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QLA
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QIA
            )
            AF = 0.0
            QLA = 0.0
            QIA = 0.0
        # Fix if LS cloud fraction too small
        if CF < 1.0e-5:
            QV = QV + QLC + QIC
            TE = (
                TE
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QLC
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * QIC
            )
            CF = 0.0
            QLC = 0.0
            QIC = 0.0
        # LS LIQUID too small
        if QLC < 1.0e-8:
            QV = QV + QLC
            TE = TE - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QLC
            QLC = 0.0
        # LS ICE too small
        if QIC < 1.0e-8:
            QV = QV + QIC
            TE = TE - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * QIC
            QIC = 0.0
        # Anvil LIQUID too small
        if QLA < 1.0e-8:
            QV = QV + QLA
            TE = TE - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QLA
            QLA = 0.0
        # Anvil ICE too small
        if QIA < 1.0e-8:
            QV = QV + QIA
            TE = TE - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * QIA
            QIA = 0.0
        # Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
        if (QLA + QIA) < 1.0e-8:
            QV = QV + QLA + QIA
            TE = (
                TE
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QLA
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * QIA
            )
            AF = 0.0
            QLA = 0.0
            QIA = 0.0
        # Fix ALL cloud quants if LS cloud LIQUID+ICE too small
        if (QLC + QIC) < 1.0e-8:
            QV = QV + QLC + QIC
            TE = (
                TE
                - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP) * QLC
                - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CP) * QIC
            )
            CF = 0.0
            QLC = 0.0
            QIC = 0.0
