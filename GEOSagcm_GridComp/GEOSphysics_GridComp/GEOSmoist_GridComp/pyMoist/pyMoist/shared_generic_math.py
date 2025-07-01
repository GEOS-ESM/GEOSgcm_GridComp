"""Stencils and functions called by multiple pyMoist modules.
These functions perform basic math and calculate fundamental
meteorological quantities"""

import gt4py.cartesian.gtscript as gtscript

import pyMoist.constants as constants
from ndsl.dsl.gt4py import exp
from ndsl.dsl.typing import Float


@gtscript.function
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


@gtscript.function
def sigma(dx) -> Float:
    sigma = 1.0 - 0.9839 * exp(-0.09835 * (dx / 1000.0))  # Arakawa 2011 sigma
    return sigma
