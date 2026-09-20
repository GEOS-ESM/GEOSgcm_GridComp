"""Stencils and functions called by multiple pyMoist modules.
These functions perform basic math and calculate fundamental
meteorological quantities"""

from ndsl.dsl.gt4py import exp, function
from ndsl.dsl.typing import Float

import pyMoist.constants as constants


@function
def air_density(PL: Float, TE: Float) -> Float:
    """
    Calculate air density [kg/m^3]

    Parameters:
    PL (Float): Pressure level.
    TE (Float): Temperature.

    Returns:
    Float: Calculated air density.
    """
    air_density = (100.0 * PL) / (constants.MAPL_RDRY * TE)
    return air_density


@function
def sigma(dx) -> Float:
    """Arakawa 2011 sigma"""
    sigma = 1.0 - 0.9839 * exp(-0.09835 * (dx / 1000.0))
    return sigma


@function
def compute_estimated_inversion_strength_factor(estimated_inversion_strength) -> Float:
    if estimated_inversion_strength >= 10.0:
        # Very stable regime
        eis_factor = 1.0
    elif estimated_inversion_strength <= 0.0:
        # Very unstable regime
        eis_factor = 0.0
    else:
        # Smooth function from 0 to 1
        eis_factor = (estimated_inversion_strength / 10.0) ** 2

    return eis_factor
