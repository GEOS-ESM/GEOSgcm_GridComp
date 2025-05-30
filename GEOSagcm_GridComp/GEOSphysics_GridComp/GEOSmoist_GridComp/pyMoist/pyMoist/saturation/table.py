from typing import Dict, Optional

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation.constants import DELTA_T, MAPL_TICE, TABLESIZE, TMINLQU, TMINTBL, TMIX
from pyMoist.saturation.formulation import SaturationFormulation
from pyMoist.saturation.qsat_ice import qsat_ice_scalar_exact
from pyMoist.saturation.qsat_liquid import qsat_liquid_scalar_exact


class SaturationVaporPressureTable:
    """
    Saturation vapor pressure table initialization. Tables are in Pa.

    These tables are automatically generated at a 0.1K resolution
    for whatever vapor pressure formulation is being used.
    """

    def __init__(
        self,
        formulation: SaturationFormulation,
    ) -> None:
        self._estimated_esw = np.empty(TABLESIZE, dtype=Float)
        self._estimated_ese = np.empty(TABLESIZE, dtype=Float)
        self._estimated_esx = np.empty(TABLESIZE, dtype=Float)
        self._TI = np.empty(TABLESIZE, dtype=Float)
        self._LOC = np.empty(TABLESIZE, dtype=Float)

        for i in range(TABLESIZE):
            t = Float(i * DELTA_T) + TMINTBL
            (
                self._estimated_esw[i],
                self._TI[i],
                _,
            ) = qsat_liquid_scalar_exact(t, formulation)

            if t > MAPL_TICE:
                self._estimated_ese[i] = self._estimated_esw[i]
            else:
                (
                    self._estimated_ese[i],
                    self._TI[i],
                    _,
                ) = qsat_ice_scalar_exact(t, formulation)

            t = t - MAPL_TICE

            if t >= TMIX and t < 0.0:
                self._estimated_esx[i] = (t / TMIX) * (
                    self._estimated_ese[i] - self._estimated_esw[i]
                ) + self._estimated_esw[i]
            else:
                self._estimated_esx[i] = self._estimated_ese[i]

        self._estimated_frz, _, _ = qsat_liquid_scalar_exact(MAPL_TICE, formulation)
        self._estimated_lqu, _, _ = qsat_liquid_scalar_exact(TMINLQU, formulation)

    @property
    def ese(self):
        return self._estimated_ese

    @property
    def esw(self):
        return self._estimated_esw

    @property
    def esx(self):
        return self._estimated_esx

    @property
    def frz(self):
        return self._estimated_frz

    @property
    def lqu(self):
        return self._estimated_lqu


# Table needs to be calculated only once
_cached_estimated_saturation: Dict[
    SaturationFormulation, Optional[SaturationVaporPressureTable]
] = {
    SaturationFormulation.MurphyAndKoop: None,
    SaturationFormulation.CAM: None,
    SaturationFormulation.Staars: None,
}


def get_table(
    formulation: SaturationFormulation = SaturationFormulation.Staars,
) -> SaturationVaporPressureTable:
    if _cached_estimated_saturation[formulation] is None:
        _cached_estimated_saturation[formulation] = SaturationVaporPressureTable(
            formulation
        )
    return _cached_estimated_saturation[formulation]
