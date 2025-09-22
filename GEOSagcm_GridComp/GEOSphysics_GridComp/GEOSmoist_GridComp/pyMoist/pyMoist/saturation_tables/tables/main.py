from typing import Dict, Optional

from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import PARALLEL, K, computation, interval
from ndsl.dsl.typing import Float, FloatField, Int
from pyMoist.saturation_tables.constants import DELTA_T, MAPL_TICE, TABLESIZE, TMINLQU, TMINTBL, TMIX
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.constants import IceExactConstants, LiquidExactConstants
from pyMoist.saturation_tables.tables.ice_exact import ice_exact, ice_exact_no_stencil
from pyMoist.saturation_tables.tables.liquid_exact import liquid_exact, liquid_exact_no_stencil


def _compute_tables(ese: FloatField, esw: FloatField, esx: FloatField, formulation: Int):
    with computation(PARALLEL), interval(...):
        t = K * DELTA_T + TMINTBL
        esw, _ = liquid_exact(
            t,
            formulation,
            LiquidExactConstants.B6,
            LiquidExactConstants.B5,
            LiquidExactConstants.B4,
            LiquidExactConstants.B3,
            LiquidExactConstants.B2,
            LiquidExactConstants.B1,
            LiquidExactConstants.B0,
            IceExactConstants.BI6,
            IceExactConstants.BI5,
            IceExactConstants.BI4,
            IceExactConstants.BI3,
            IceExactConstants.BI2,
            IceExactConstants.BI1,
            IceExactConstants.BI0,
            IceExactConstants.S16,
            IceExactConstants.S15,
            IceExactConstants.S14,
            IceExactConstants.S13,
            IceExactConstants.S12,
            IceExactConstants.S11,
            IceExactConstants.S10,
            IceExactConstants.S26,
            IceExactConstants.S25,
            IceExactConstants.S24,
            IceExactConstants.S23,
            IceExactConstants.S22,
            IceExactConstants.S21,
            IceExactConstants.S20,
            LiquidExactConstants.DL_0,
            LiquidExactConstants.DL_1,
            LiquidExactConstants.DL_2,
            LiquidExactConstants.DL_3,
            LiquidExactConstants.DL_4,
            LiquidExactConstants.DL_5,
            LiquidExactConstants.TS,
            LiquidExactConstants.LOGPS,
            LiquidExactConstants.CL_0,
            LiquidExactConstants.CL_1,
            LiquidExactConstants.CL_2,
            LiquidExactConstants.CL_3,
            LiquidExactConstants.CL_4,
            LiquidExactConstants.CL_5,
            LiquidExactConstants.CL_6,
            LiquidExactConstants.CL_7,
            LiquidExactConstants.CL_8,
            LiquidExactConstants.CL_9,
        )

        if t > MAPL_TICE:
            ese = esw
        else:
            ese, _ = ice_exact(
                t,
                formulation,
                IceExactConstants.TMINSTR,
                IceExactConstants.TMINICE,
                IceExactConstants.TSTARR1,
                IceExactConstants.TSTARR2,
                IceExactConstants.TSTARR3,
                IceExactConstants.TSTARR4,
                IceExactConstants.TMAXSTR,
                IceExactConstants.DI_0,
                IceExactConstants.DI_1,
                IceExactConstants.DI_2,
                IceExactConstants.DI_3,
                IceExactConstants.CI_0,
                IceExactConstants.CI_1,
                IceExactConstants.CI_2,
                IceExactConstants.CI_3,
                IceExactConstants.S16,
                IceExactConstants.S15,
                IceExactConstants.S14,
                IceExactConstants.S13,
                IceExactConstants.S12,
                IceExactConstants.S11,
                IceExactConstants.S10,
                IceExactConstants.S26,
                IceExactConstants.S25,
                IceExactConstants.S24,
                IceExactConstants.S23,
                IceExactConstants.S22,
                IceExactConstants.S21,
                IceExactConstants.S20,
                IceExactConstants.BI6,
                IceExactConstants.BI5,
                IceExactConstants.BI4,
                IceExactConstants.BI3,
                IceExactConstants.BI2,
                IceExactConstants.BI1,
                IceExactConstants.BI0,
            )
        t = t - MAPL_TICE

        if t >= TMIX and t < 0.0:
            esx = (t / TMIX) * (ese - esw) + esw
        else:
            esx = ese


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
            formulation_int = Int(1)
        elif formulation == SaturationFormulation.CAM:
            formulation_int = Int(2)
        elif formulation == SaturationFormulation.MurphyAndKoop:
            formulation_int = Int(3)

        self._estimated_ese = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._estimated_esw = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        self._estimated_esx = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        compute_tables = stencil_factory.from_dims_halo(
            func=_compute_tables,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        NDSL_tables = False
        # NOTE Setting NDSL_tables to True will the tables using gt4py stencils and functions.
        # This is turned off because they do not currently verify.
        # Does not verify because:
        #   - Need the ability to pass 32 or 64 bit externals flexibly. Right now it appears they all
        #       get forced to the same precision

        if NDSL_tables == True:
            compute_tables(self._estimated_ese, self._estimated_esw, self._estimated_esx, formulation_int)
        else:
            for i in range(TABLESIZE):
                t = Float(i * DELTA_T) + TMINTBL
                self._estimated_esw.field[0, 0, i], _ = liquid_exact_no_stencil(t, formulation)

                if t > MAPL_TICE:
                    self._estimated_ese.field[0, 0, i] = self._estimated_esw.field[0, 0, i]
                else:
                    self._estimated_ese.field[0, 0, i], _ = ice_exact_no_stencil(t, formulation)

                t = t - MAPL_TICE

                if t >= TMIX and t < 0.0:
                    self._estimated_esx.field[0, 0, i] = (t / TMIX) * (
                        self._estimated_ese.field[0, 0, i] - self._estimated_esw.field[0, 0, i]
                    ) + self._estimated_esw.field[0, 0, i]
                else:
                    self._estimated_esx.field[0, 0, i] = self._estimated_ese.field[0, 0, i]

        self._estimated_frz, _ = liquid_exact_no_stencil(MAPL_TICE, formulation)
        self._estimated_lqu, _ = liquid_exact_no_stencil(TMINLQU, formulation)

    @property
    def ese(self):
        return self._estimated_ese.field[0, 0, :]

    @property
    def esw(self):
        return self._estimated_esw.field[0, 0, :]

    @property
    def esx(self):
        return self._estimated_esx.field[0, 0, :]

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
