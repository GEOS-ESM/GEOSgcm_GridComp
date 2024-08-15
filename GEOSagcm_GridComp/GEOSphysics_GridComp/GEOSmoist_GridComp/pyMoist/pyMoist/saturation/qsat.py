from pyMoist.saturation.formulation import SaturationFormulation
from typing import Any
from ndsl.dsl.typing import Float
import gt4py.cartesian.gtscript as gtscript


class QSat:
    """Traditional satuation specific humidity.
    In Fortran: GEOS_Utilities:QSat

    Uses various formulations of the saturation vapor pressure to compute the saturation specific
    humidity for temperature TL and pressure PL.

    For temperatures <= TMIX (-20C)
    the calculation is done over ice; for temperatures >= ZEROC (0C) the calculation
    is done over liquid water; and in between these values,
    it interpolates linearly between the two.

    The optional RAMP is the width of this
    ice/water ramp (i.e., TMIX = ZEROC-RAMP); its default is 20.

    If PASCALS is true, PL is
    assumed to be in Pa; if false or not present, it is assumed to be in mb.

    Another compile time choice is whether to use the exact formulation
    or a table look-up.
    If UTBL is true, tabled values of the saturation vapor pressures
    are used. These tables are automatically generated at a 0.1K resolution
    for whatever vapor pressure formulation is being used.
    """

    def __init__(
        self,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
        use_table_lookup: bool = True,
    ) -> None:
        if not use_table_lookup:
            raise NotImplementedError(
                "Saturation calculation: exact formulation not available, only table look up"
            )

    @gtscript.function
    def __call__(self, RAMP: Float, pascals: bool, dqsat: Float) -> Float:
        pass
