from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField, Int
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
import pyMoist.constants as constants


def hydrostatic_air_density(
    p: FloatField_Plume,
    geopotential_height: FloatField,
    error_code: IntFieldIJ_Plume,
    air_density: FloatField,
    plume: Int,
):
    """Compute air density, assuming hydrostatic balance.

    Args:
        p (FloatField_Plume)
        geopotential_height (FloatField)
        error_code (IntFieldIJ_Plume)
        air_density (FloatField)
        plume (Int)
    """
    with computation(PARALLEL), interval(...):
        # prefil with 0
        air_density = 0.0

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            air_density = (
                100.0
                * (p[0, 0, 0][plume] - p[0, 0, 1][plume])
                / (geopotential_height[0, 0, 1] - geopotential_height)
                / constants.MAPL_GRAV
            )
