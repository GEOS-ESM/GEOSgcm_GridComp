from ndsl.dsl.gt4py import FORWARD, PARALLEL, atan, computation, interval, sqrt, tan
from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntFieldIJ

import pyMoist.constants as constants


def rh_calculations(
    estimated_inversion_strength: FloatFieldIJ,
    p_mb: FloatField,
    p_interface_mb: FloatField,
    area: FloatFieldIJ,
    lcl_level: IntFieldIJ,
    alpha: FloatField,
):
    """Compute relative humidity for use in the PDF

    Args:
        estimated_inversion_strength (FloatFieldIJ)
        p_mb (FloatField)
        p_interface_mb (FloatField)
        area (FloatFieldIJ)
        lcl_level (IntFieldIJ)
        alpha (FloatField)
    """
    from __externals__ import DW_LAND, DW_OCEAN, TURNRHCRIT_PARAM, k_end

    with computation(FORWARD), interval(0, 1):
        # Send the condensates through the pdf after convection
        fac_eis: FloatFieldIJ = max(0.0, min(1.0, estimated_inversion_strength / 10.0)) ** 2
        # determine combined min_rh_crit in stable/unstable regimes
        min_rh_crit: FloatFieldIJ = (1.0 - DW_OCEAN) * (1.0 - fac_eis) + (1.0 - DW_LAND) * fac_eis

    with computation(PARALLEL), interval(...):
        # determine the turn pressure using the LCL
        if TURNRHCRIT_PARAM <= 0:
            turnrhcrit = p_mb.at(K=lcl_level) - 250
        else:
            turnrhcrit = TURNRHCRIT_PARAM

    with computation(PARALLEL), interval(0, -1):
        # Use Slingo-Ritter (1985) formulation for critical relative humidity
        rh_crit = 1.0
        # lower turn from maxrhcrit=1.0
        if p_mb <= turnrhcrit:
            rh_crit = min_rh_crit
        else:
            rh_crit = min_rh_crit + (1.0 - min_rh_crit) / (19.0) * (
                (
                    atan(
                        (2.0 * (p_mb - turnrhcrit) / (p_interface_mb.at(K=k_end + 1) - turnrhcrit) - 1.0)
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
            rh_crit = min_rh_crit
        else:
            rh_crit = 1.0

    with computation(PARALLEL), interval(...):
        # include grid cell area scaling and limit RHcrit to > 70%\
        # NOTE there is a fundamental mathematical difference between the python and fortran here
        # the double square root is resolved differently, producing errors on the order of 1e-7
        # despite all inputs (area & rh_crit) being identical (0 ULP difference)
        alpha = max(0.0, min(0.30, (1.0 - rh_crit) * sqrt(sqrt(area / 1.0e10))))


def fill_rh_crit_export(alpha: FloatField, rh_crit: FloatField):
    """Export relative humidity to the model - only called if rh_crit is associated in the Fortran

    Args:
        alpha (FloatField)
        rh_crit (FloatField)
    """
    with computation(PARALLEL), interval(...):
        rh_crit = 1 - alpha
