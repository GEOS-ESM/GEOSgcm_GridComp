from typing import Dict, Optional

from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import Z_DIM
from ndsl.dsl.typing import Float
from pyMoist.saturation_tables.constants import DELTA_T, MAPL_TICE, TABLESIZE, TMINLQU, TMINTBL, TMIX
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.ice_exact import ice_exact
from pyMoist.saturation_tables.tables.liquid_exact import liquid_exact


class SaturationVaporPressureTable:
    """
    Saturation vapor pressure table initialization. Tables are in Pa.

    These tables are automatically generated at a 0.1K resolution
    for whatever vapor pressure formulation is being used.
    """

    def __init__(
        self,
        backend,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ) -> None:
        table_compute_domain = (1, 1, TABLESIZE)

        stencil_factory, quantity_factory = get_factories_single_tile(
            table_compute_domain[0],
            table_compute_domain[1],
            table_compute_domain[2],
            0,
            backend,
        )

        if formulation == SaturationFormulation.Staars:
            formulation_int = 1
        elif formulation == SaturationFormulation.CAM:
            formulation_int = 2
        elif formulation == SaturationFormulation.MurphyAndKoop:
            formulation_int = 3

        # self._estimated_esw = np.empty(TABLESIZE, dtype=Float)
        # self._estimated_ese = np.empty(TABLESIZE, dtype=Float)
        # self._estimated_esx = np.empty(TABLESIZE, dtype=Float)

        self._estimated_ese = quantity_factory.zeros([Z_DIM], "n/a")
        self._estimated_esw = quantity_factory.zeros([Z_DIM], "n/a")
        self._estimated_esx = quantity_factory.zeros([Z_DIM], "n/a")

        for i in range(TABLESIZE):
            t = Float(i * DELTA_T) + TMINTBL
            self._estimated_esw.field[i], _ = liquid_exact(t, formulation)

            if t > MAPL_TICE:
                self._estimated_ese.field[i] = self._estimated_esw.field[i]
            else:
                self._estimated_ese.field[i], _ = ice_exact(t, formulation)

            t = t - MAPL_TICE

            if t >= TMIX and t < 0.0:
                self._estimated_esx.field[i] = (t / TMIX) * (
                    self._estimated_ese.field[i] - self._estimated_esw.field[i]
                ) + self._estimated_esw.field[i]
            else:
                self._estimated_esx.field[i] = self._estimated_ese.field[i]

        self._estimated_frz, _ = liquid_exact(MAPL_TICE, formulation)
        self._estimated_lqu, _ = liquid_exact(TMINLQU, formulation)

    @property
    def ese(self):
        return self._estimated_ese.field

    @property
    def esw(self):
        return self._estimated_esw.field

    @property
    def esx(self):
        return self._estimated_esx.field

    @property
    def frz(self):
        return self._estimated_frz

    @property
    def lqu(self):
        return self._estimated_lqu


# Table needs to be calculated only once
_cached_estimated_saturation: Dict[SaturationFormulation, Optional[SaturationVaporPressureTable]] = {
    SaturationFormulation.MurphyAndKoop: None,
    SaturationFormulation.CAM: None,
    SaturationFormulation.Staars: None,
}


def get_table(
    backend,
    formulation: SaturationFormulation = SaturationFormulation.Staars,
) -> SaturationVaporPressureTable:
    if _cached_estimated_saturation[formulation] is None:
        _cached_estimated_saturation[formulation] = SaturationVaporPressureTable(backend, formulation)
    return _cached_estimated_saturation[formulation]
