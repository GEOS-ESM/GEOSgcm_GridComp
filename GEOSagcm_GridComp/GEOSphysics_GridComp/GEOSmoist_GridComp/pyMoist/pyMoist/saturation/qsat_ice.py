from typing import Optional

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation.constants import (
    DELTA_T,
    ERFAC,
    ESFAC,
    MAPL_TICE,
    MAX_MIXING_RATIO,
)
from pyMoist.saturation.formulation import SaturationFormulation


TMINSTR = Float(-95.0)
TMINICE = MAPL_TICE + TMINSTR


TMINSTR = Float(-95.0)
TSTARR1 = Float(-75.0)
TSTARR2 = Float(-65.0)
TSTARR3 = Float(-50.0)
TSTARR4 = Float(-40.0)
TMAXSTR = Float(+60.0)

DI = [57518.5606e08, 2.01889049, 3.56654, 20.947031]
CI = [9.550426, -5723.265, 3.53068, -0.00728332]

# 64-bit float in Fortran
S16 = 0.516000335e-11 * 100.0
S15 = 0.276961083e-8 * 100.0
S14 = 0.623439266e-6 * 100.0
S13 = 0.754129933e-4 * 100.0
S12 = 0.517609116e-2 * 100.0
S11 = 0.191372282e0 * 100.0
S10 = 0.298152339e1 * 100.0
S26 = 0.314296723e-10 * 100.0
S25 = 0.132243858e-7 * 100.0
S24 = 0.236279781e-5 * 100.0
S23 = 0.230325039e-3 * 100.0
S22 = 0.129690326e-1 * 100.0
S21 = 0.401390832e0 * 100.0
S20 = 0.535098336e1 * 100.0
BI6 = 1.838826904e-10 * 100.0
BI5 = 4.838803174e-8 * 100.0
BI4 = 5.824720280e-6 * 100.0
BI3 = 4.176223716e-4 * 100.0
BI2 = 1.886013408e-2 * 100.0
BI1 = 5.034698970e-1 * 100.0
BI0 = 6.109177956e0 * 100.0


def _saturation_formulation(
    formulation: SaturationFormulation,
    t: Float,
):
    if formulation == SaturationFormulation.Staars:
        TT = t - MAPL_TICE
        if TT < TSTARR1:
            LOC = 1.1
            EX = (
                TT
                * (TT * (TT * (TT * (TT * (TT * S16 + S15) + S14) + S13) + S12) + S11)
                + S10
            )
        elif TT >= TSTARR1 and TT < TSTARR2:
            LOC = 1.2
            W = (TSTARR2 - TT) / (TSTARR2 - TSTARR1)
            EX = W * (
                TT
                * (TT * (TT * (TT * (TT * (TT * S16 + S15) + S14) + S13) + S12) + S11)
                + S10
            ) + (1.0 - W) * (
                TT
                * (TT * (TT * (TT * (TT * (TT * S26 + S25) + S24) + S23) + S22) + S21)
                + S20
            )
        elif TT >= TSTARR2 and TT < TSTARR3:
            LOC = 1.3
            EX = (
                TT
                * (TT * (TT * (TT * (TT * (TT * S26 + S25) + S24) + S23) + S22) + S21)
                + S20
            )
        elif TT >= TSTARR3 and TT < TSTARR4:
            LOC = 1.4
            W = (TSTARR4 - TT) / (TSTARR4 - TSTARR3)
            EX = W * (
                TT
                * (TT * (TT * (TT * (TT * (TT * S26 + S25) + S24) + S23) + S22) + S21)
                + S20
            ) + (1.0 - W) * (
                TT
                * (TT * (TT * (TT * (TT * (TT * BI6 + BI5) + BI4) + BI3) + BI2) + BI1)
                + BI0
            )
        else:
            LOC = 1.5
            EX = (
                TT
                * (TT * (TT * (TT * (TT * (TT * BI6 + BI5) + BI4) + BI3) + BI2) + BI1)
                + BI0
            )
    elif formulation == SaturationFormulation.CAM:
        LOC = 2.0
        TT = MAPL_TICE / t
        EX = DI[0] * np.exp(-(DI[1] / TT + DI[2] * np.log(TT) + DI[3] * TT))
    elif formulation == SaturationFormulation.MurphyAndKoop:
        LOC = 3.0
        EX = np.exp(CI[0] + CI[1] / t + CI[2] * np.log(t) + CI[3] * t)
    return Float(EX), LOC


def qsat_ice_scalar_exact(
    temperature: Float,
    formulation: SaturationFormulation = SaturationFormulation.Staars,
    PL: Optional[Float] = None,
    DQ: Optional[Float] = None,
):
    """Reference Fortran: QSATICE0 w/ UTBL=False"""
    if temperature < TMINICE:
        TI = TMINICE
    elif temperature > MAPL_TICE:
        TI = MAPL_TICE
    else:
        TI = temperature

    DX = 0.0  # only calulcated when DQ is not none
    EX, LOC = _saturation_formulation(formulation, TI)

    if DQ is not None:
        if temperature < TMINICE:
            DDQ = 0.0
        elif temperature > MAPL_TICE:
            DDQ = 0.0
        else:
            if PL > EX:
                DD = EX
                TI = temperature + DELTA_T
                EX, _ = _saturation_formulation(formulation, TI)
                DDQ = EX - DD
                EX = DD
    if PL is not None:
        if PL > EX:
            DD = ESFAC / (PL - (1.0 - ESFAC) * EX)
            EX = EX * DD
            if DQ is not None:
                DX = DDQ * ERFAC * PL * DD * DD
        else:
            EX = MAX_MIXING_RATIO
            if DQ is not None:
                DX = 0.0
    else:
        if DQ is not None:
            DX = DDQ * (1.0 / DELTA_T)

    return EX, TI, DX, LOC
