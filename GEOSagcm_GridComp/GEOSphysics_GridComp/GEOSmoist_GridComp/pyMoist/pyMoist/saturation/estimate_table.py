from ndsl.dsl.typing import Float
import numpy as np
from pyMoist.saturation.constants import (
    MAPL_TICE,
    TMINLQU,
    TMINTBL,
    TABLESIZE,
    DELTA_T,
    TMIX,
)
from pyMoist.saturation.qsat_liquid import qsat_liquid_scalar_exact
from pyMoist.saturation.formulation import SaturationFormulation


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
        self._estimated_blw = np.empty(TABLESIZE, dtype=Float)
        self._estimated_ble = np.empty(TABLESIZE, dtype=Float)
        self._estimated_blx = np.empty(TABLESIZE, dtype=Float)

        for i in range(TABLESIZE):
            t = (i - 1) * DELTA_T + TMINTBL

            self._estimated_blw[i], _ = qsat_liquid_scalar_exact(t, formulation)

            if t > MAPL_TICE:
                self._estimated_ble[i] = self._estimated_blw[i]
            else:
                self._estimated_ble[i] = QSATICE0(t, formulation)

            t = t - MAPL_TICE

            if t >= TMIX and t < 0.0:
                self._estimated_blx[i] = (t / TMIX) * (
                    self._estimated_ble[i] - self._estimated_blw[i]
                ) + self._estimated_blw[i]
            else:
                self._estimated_blx[i] = self._estimated_ble[i]

        self.estimated_frz, _ = qsat_liquid_scalar_exact(MAPL_TICE, formulation)
        self.estimated_lqu, _ = qsat_liquid_scalar_exact(TMINLQU, formulation)

    @property
    def blw(self):
        return self._estimated_blw

    @property
    def ble(self):
        return self._estimated_ble

    @property
    def blx(self):
        return self._estimated_blx

    @property
    def frz(self):
        return self._estimated_frz

    @property
    def lqu(self):
        return self._estimated_lqu


# Table needs to be calculated only once
_cached_estimated_saturation = {
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
