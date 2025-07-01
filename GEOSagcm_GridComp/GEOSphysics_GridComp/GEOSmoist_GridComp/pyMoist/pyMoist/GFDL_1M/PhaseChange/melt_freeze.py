from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval

import pyMoist.constants as constants
from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.shared_incloud_processes import ice_fraction


def melt_freeze(
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    t: FloatField,
    liquid: FloatField,
    ice: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        if t <= constants.MAPL_TICE:
            f_qi = ice_fraction(t, convection_fraction, surface_type)
            d_qil = liquid * (1.0 - exp(-DT_MOIST * f_qi / constants.TAUFRZ))
            d_qil = max(0.0, d_qil)
            ice = ice + d_qil
            liquid = liquid - d_qil
            t = (
                t
                + (constants.MAPL_LATENT_HEAT_SUBLIMATION - constants.MAPL_LATENT_HEAT_VAPORIZATION)
                * d_qil
                / constants.MAPL_CP
            )
        else:
            d_qil = -ice
            d_qil = min(0.0, d_qil)
            ice = ice + d_qil
            liquid = liquid - d_qil
            t = (
                t
                + (constants.MAPL_LATENT_HEAT_SUBLIMATION - constants.MAPL_LATENT_HEAT_VAPORIZATION)
                * d_qil
                / constants.MAPL_CP
            )
