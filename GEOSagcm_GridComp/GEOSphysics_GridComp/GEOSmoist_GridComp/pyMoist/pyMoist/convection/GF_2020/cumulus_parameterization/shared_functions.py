from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, IntField, Float, Int
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
    int32,
    log,
    exp,
)
from ndsl.dsl.gt4py import function


@function
def saturation_vapor_pressure(t: Float):
    """
    Compute saturation vapor pressure

    Args:
        t (in): temperature in Kelvin
    """
    # CAUTION: This function has not been verified!!!
    t_celcius = t - 273.155
    if t_celcius < -20.0:
        toot = 273.16 / t
        toto = 1 / toot
        eilog = (
            -9.09718 * (toot - 1)
            - 3.56654 * (log(toot) / log(10.0))
            + 0.876793 * (1 - toto)
            + (log(6.1071) / log(10.0))
        )
        saturation_vapor_pressure = 10**eilog
    else:
        tsot = 373.16 / t
        ewlog = -7.90298 * (tsot - 1) + 5.02808 * (log(tsot) / log(10.0))
        ewlog2 = ewlog - 1.3816e-07 * (10 ** (11.344 * (1 - (1 / tsot))) - 1)
        ewlog3 = ewlog2 + 0.0081328 * (10 ** (-3.49149 * (tsot - 1)) - 1)
        ewlog4 = ewlog3 + (log(1013.246) / log(10.0))
        saturation_vapor_pressure = 10**ewlog4

    return saturation_vapor_pressure
