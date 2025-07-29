from typing import Dict, Optional

from ndsl.boilerplate import get_factories_single_tile
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, FloatField, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval
from pyMoist.saturation_tables.constants import DELTA_T, MAPL_TICE, TABLESIZE, TMINLQU, TMINTBL, TMIX
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.ice_exact import ice_exact_no_stencil, ice_exact
from pyMoist.saturation_tables.tables.liquid_exact import liquid_exact_no_stencil, liquid_exact
from gt4py.cartesian.gtscript import THIS_K
from pyMoist.saturation_tables.tables.constants import IceExactConstatns
from pyMoist.saturation_tables.tables.constants import LiquidExactConstants


def _compute_tables(ese: FloatField, esw: FloatField, esx: FloatField, formulation: Int):
    from __externals__ import (
        TMINSTR,
        TMINICE,
        TSTARR1,
        TSTARR2,
        TSTARR3,
        TSTARR4,
        TMAXSTR,
        DI_0,
        DI_1,
        DI_2,
        DI_3,
        CI_0,
        CI_1,
        CI_2,
        CI_3,
        S16,
        S15,
        S14,
        S13,
        S12,
        S11,
        S10,
        S26,
        S25,
        S24,
        S23,
        S22,
        S21,
        S20,
        BI6,
        BI5,
        BI4,
        BI3,
        BI2,
        BI1,
        BI0,
        B6,
        B5,
        B4,
        B3,
        B2,
        B1,
        B0,
        BI6,
        BI5,
        BI4,
        BI3,
        BI2,
        BI1,
        BI0,
        S16,
        S15,
        S14,
        S13,
        S12,
        S11,
        S10,
        S26,
        S25,
        S24,
        S23,
        S22,
        S21,
        S20,
        DL_0,
        DL_1,
        DL_2,
        DL_3,
        DL_4,
        DL_5,
        TS,
        LOGPS,
        CL_0,
        CL_1,
        CL_2,
        CL_3,
        CL_4,
        CL_5,
        CL_6,
        CL_7,
        CL_8,
        CL_9,
    )

    with computation(PARALLEL), interval(...):
        t = THIS_K * DELTA_T + TMINTBL
        esw = (
            liquid_exact(
                t,
                formulation,
                B6,
                B5,
                B4,
                B3,
                B2,
                B1,
                B0,
                BI6,
                BI5,
                BI4,
                BI3,
                BI2,
                BI1,
                BI0,
                S16,
                S15,
                S14,
                S13,
                S12,
                S11,
                S10,
                S26,
                S25,
                S24,
                S23,
                S22,
                S21,
                S20,
                DL_0,
                DL_1,
                DL_2,
                DL_3,
                DL_4,
                DL_5,
                TS,
                LOGPS,
                CL_0,
                CL_1,
                CL_2,
                CL_3,
                CL_4,
                CL_5,
                CL_6,
                CL_7,
                CL_8,
                CL_9,
            ),
        )

        if t > MAPL_TICE:
            ese = esw
        else:
            ese = ice_exact(
                t,
                formulation,
                TMINSTR,
                TMINICE,
                TSTARR1,
                TSTARR2,
                TSTARR3,
                TSTARR4,
                TMAXSTR,
                DI_0,
                DI_1,
                DI_2,
                DI_3,
                CI_0,
                CI_1,
                CI_2,
                CI_3,
                S16,
                S15,
                S14,
                S13,
                S12,
                S11,
                S10,
                S26,
                S25,
                S24,
                S23,
                S22,
                S21,
                S20,
                BI6,
                BI5,
                BI4,
                BI3,
                BI2,
                BI1,
                BI0,
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
            externals={
                "TMINSTR": IceExactConstatns.TMINSTR,
                "TMINICE": IceExactConstatns.TMINICE,
                "TSTARR1": IceExactConstatns.TSTARR1,
                "TSTARR2": IceExactConstatns.TSTARR2,
                "TSTARR3": IceExactConstatns.TSTARR3,
                "TSTARR4": IceExactConstatns.TSTARR4,
                "TMAXSTR": IceExactConstatns.TMAXSTR,
                "DI_0": IceExactConstatns.DI[0],
                "DI_1": IceExactConstatns.DI[1],
                "DI_2": IceExactConstatns.DI[2],
                "DI_3": IceExactConstatns.DI[3],
                "CI_0": IceExactConstatns.CI[0],
                "CI_1": IceExactConstatns.CI[1],
                "CI_2": IceExactConstatns.CI[2],
                "CI_3": IceExactConstatns.CI[3],
                "S16": IceExactConstatns.S16,
                "S15": IceExactConstatns.S15,
                "S14": IceExactConstatns.S14,
                "S13": IceExactConstatns.S13,
                "S12": IceExactConstatns.S12,
                "S11": IceExactConstatns.S11,
                "S10": IceExactConstatns.S10,
                "S26": IceExactConstatns.S26,
                "S25": IceExactConstatns.S25,
                "S24": IceExactConstatns.S24,
                "S23": IceExactConstatns.S23,
                "S22": IceExactConstatns.S22,
                "S21": IceExactConstatns.S21,
                "S20": IceExactConstatns.S20,
                "BI6": IceExactConstatns.BI6,
                "BI5": IceExactConstatns.BI5,
                "BI4": IceExactConstatns.BI4,
                "BI3": IceExactConstatns.BI3,
                "BI2": IceExactConstatns.BI2,
                "BI1": IceExactConstatns.BI1,
                "BI0": IceExactConstatns.BI0,
                "B6": LiquidExactConstants.B6,
                "B5": LiquidExactConstants.B5,
                "B4": LiquidExactConstants.B4,
                "B3": LiquidExactConstants.B3,
                "B2": LiquidExactConstants.B2,
                "B1": LiquidExactConstants.B1,
                "B0": LiquidExactConstants.B0,
                "BI6": LiquidExactConstants.BI6,
                "BI5": LiquidExactConstants.BI5,
                "BI4": LiquidExactConstants.BI4,
                "BI3": LiquidExactConstants.BI3,
                "BI2": LiquidExactConstants.BI2,
                "BI1": LiquidExactConstants.BI1,
                "BI0": LiquidExactConstants.BI0,
                "S16": LiquidExactConstants.S16,
                "S15": LiquidExactConstants.S15,
                "S14": LiquidExactConstants.S14,
                "S13": LiquidExactConstants.S13,
                "S12": LiquidExactConstants.S12,
                "S11": LiquidExactConstants.S11,
                "S10": LiquidExactConstants.S10,
                "S26": LiquidExactConstants.S26,
                "S25": LiquidExactConstants.S25,
                "S24": LiquidExactConstants.S24,
                "S23": LiquidExactConstants.S23,
                "S22": LiquidExactConstants.S22,
                "S21": LiquidExactConstants.S21,
                "S20": LiquidExactConstants.S20,
                "DL_0": LiquidExactConstants.DL[0],
                "DL_1": LiquidExactConstants.DL[1],
                "DL_2": LiquidExactConstants.DL[2],
                "DL_3": LiquidExactConstants.DL[3],
                "DL_4": LiquidExactConstants.DL[4],
                "DL_5": LiquidExactConstants.DL[5],
                "TS": LiquidExactConstants.TS,
                "LOGPS": LiquidExactConstants.LOGPS,
                "CL_0": LiquidExactConstants.CL[0],
                "CL_1": LiquidExactConstants.CL[1],
                "CL_2": LiquidExactConstants.CL[2],
                "CL_3": LiquidExactConstants.CL[3],
                "CL_4": LiquidExactConstants.CL[4],
                "CL_5": LiquidExactConstants.CL[5],
                "CL_6": LiquidExactConstants.CL[6],
                "CL_7": LiquidExactConstants.CL[7],
                "CL_8": LiquidExactConstants.CL[8],
                "CL_9": LiquidExactConstants.CL[9],
            },
        )

        # compute_tables(self._estimated_ese, self._estimated_esw, self._estimated_esx, formulation_int)
        # NOTE This function will create the tables using gt4py stencils and functions.
        # This is turned off because they do not currently verify.
        # Two reasons for failures:
        #   - Need the ability to pass 32 OR 64 bit externals flexibly
        #   - Need the ability to return multiple values from a function
        #       (this has worked successfully before, but for some reason is triggering an error here)

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
