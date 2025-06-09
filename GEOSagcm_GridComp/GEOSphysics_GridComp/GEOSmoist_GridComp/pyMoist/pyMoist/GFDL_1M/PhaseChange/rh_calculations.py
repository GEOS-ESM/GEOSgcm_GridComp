from gt4py.cartesian.gtscript import PARALLEL, atan, computation, interval, sqrt, tan

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, IntFieldIJ


def rh_calculations(
    estimated_inversion_strength: FloatFieldIJ,
    minrhcrit: FloatField,
    p_mb: FloatField,
    p_interface_mb: FloatField,
    area: FloatFieldIJ,
    alpha: FloatField,
    k_lcl: IntFieldIJ,
    rh_crit_3d: FloatField,
):
    from __externals__ import DW_LAND, DW_OCEAN, TURNRHCRIT_PARAM, k_end

    with computation(PARALLEL), interval(...):
        # Send the condensates through the pdf after convection
        fac_eis = max(0.0, min(1.0, estimated_inversion_strength / 10.0)) ** 2
        # determine combined minrhcrit in stable/unstable regimes
        minrhcrit = (1.0 - DW_OCEAN) * (1.0 - fac_eis) + (1.0 - DW_LAND) * fac_eis
        # determine the turn pressure using the LCL
        if TURNRHCRIT_PARAM <= 0:
            turnrhcrit = p_mb.at(K=k_lcl - 1) - 250
        else:
            turnrhcrit = TURNRHCRIT_PARAM

    with computation(PARALLEL), interval(0, -1):
        # Use Slingo-Ritter (1985) formulation for critical rel ative humidity
        rh_crit = 1.0
        # lower turn from maxrhcrit=1.0
        if p_mb <= turnrhcrit:
            rh_crit = minrhcrit
        else:
            rh_crit = minrhcrit + (1.0 - minrhcrit) / (19.0) * (
                (
                    atan(
                        (2.0 * (p_mb - turnrhcrit) / (p_interface_mb.at(K=k_end) - turnrhcrit) - 1.0)
                        * tan(20.0 * constants.MAPL_PI / 21.0 - 0.5 * constants.MAPL_PI)
                    )
                    + 0.5 * constants.MAPL_PI
                )
                * 21.0
                / constants.MAPL_PI
                - 1.0
            )
    with computation(PARALLEL), interval(-1, None):
        # lower turn from maxrhcrit=1.0
        if p_mb <= turnrhcrit:
            rh_crit = minrhcrit
        else:
            rh_crit = 1.0
    with computation(PARALLEL), interval(...):
        # include grid cell area scaling and limit RHcrit to > 70%\
        alpha = max(0.0, min(0.30, (1.0 - rh_crit) * sqrt(sqrt(area / 1.0e10))))
        rh_crit_3d = 1 - alpha
