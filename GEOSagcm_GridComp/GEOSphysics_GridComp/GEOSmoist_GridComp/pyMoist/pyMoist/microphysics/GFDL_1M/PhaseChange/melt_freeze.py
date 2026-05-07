from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval
from ndsl.dsl.typing import FloatField, FloatFieldIJ

import pyMoist.constants as constants
from pyMoist.shared.incloud_processes import ice_fraction


def melt_freeze(
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    t: FloatField,
    mixing_ratio_liquid: FloatField,
    mixing_ratio_ice: FloatField,
):
    """Melting/freezing of condensates

    Args:
        convection_fraction (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        t (FloatField)
        mixing_ratio_liquid (FloatField)
        mixing_ratio_ice (FloatField)
    """
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        if t <= constants.MAPL_TICE:
            f_qi = ice_fraction(t, convection_fraction, surface_type)
            d_qil = mixing_ratio_liquid * (1.0 - exp(-DT_MOIST * f_qi / constants.TAUFRZ))
            d_qil = max(0.0, d_qil)
            mixing_ratio_ice = mixing_ratio_ice + d_qil
            mixing_ratio_liquid = mixing_ratio_liquid - d_qil
            t = t + (constants.MAPL_LATENT_HEAT_SUBLIMATION - constants.MAPL_LATENT_HEAT_VAPORIZATION) * d_qil / constants.MAPL_CP
        else:
            d_qil = -mixing_ratio_ice
            d_qil = min(0.0, d_qil)
            mixing_ratio_ice = mixing_ratio_ice + d_qil
            mixing_ratio_liquid = mixing_ratio_liquid - d_qil
            t = t + (constants.MAPL_LATENT_HEAT_SUBLIMATION - constants.MAPL_LATENT_HEAT_VAPORIZATION) * d_qil / constants.MAPL_CP
