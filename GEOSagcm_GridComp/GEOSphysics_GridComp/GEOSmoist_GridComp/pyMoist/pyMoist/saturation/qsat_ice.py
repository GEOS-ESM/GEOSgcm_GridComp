from typing import Optional

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation.constants import DELTA_T, ERFAC, ESFAC, MAPL_TICE, MAX_MIXING_RATIO
from pyMoist.saturation.formulation import SaturationFormulation


f64 = np.float64

TMINSTR = Float(-95.0)
TMINICE = MAPL_TICE + TMINSTR

TMINSTR = Float(-95.0)
TSTARR1 = Float(-75.0)
TSTARR2 = Float(-65.0)
TSTARR3 = Float(-50.0)
TSTARR4 = Float(-40.0)
TMAXSTR = Float(+60.0)

DI = [
    f64(Float(57518.5606e08)),
    f64(Float(2.01889049)),
    f64(Float(3.56654)),
    f64(Float(20.947031)),
]
CI = [
    f64(Float(9.550426)),
    f64(Float(-5723.265)),
    f64(Float(3.53068)),
    f64(Float(-0.00728332)),
]

# 64-bit float in Fortran
S16 = f64(Float(0.516000335e-11) * Float(100.0))
S15 = f64(Float(0.276961083e-8) * Float(100.0))
S14 = f64(Float(0.623439266e-6) * Float(100.0))
S13 = f64(Float(0.754129933e-4) * Float(100.0))
S12 = f64(Float(0.517609116e-2) * Float(100.0))
S11 = f64(Float(0.191372282e0) * Float(100.0))
S10 = f64(Float(0.298152339e1) * Float(100.0))
S26 = f64(Float(0.314296723e-10) * Float(100.0))
S25 = f64(Float(0.132243858e-7) * Float(100.0))
S24 = f64(Float(0.236279781e-5) * Float(100.0))
S23 = f64(Float(0.230325039e-3) * Float(100.0))
S22 = f64(Float(0.129690326e-1) * Float(100.0))
S21 = f64(Float(0.401390832e0) * Float(100.0))
S20 = f64(Float(0.535098336e1) * Float(100.0))
BI6 = f64(Float(1.838826904e-10) * Float(100.0))
BI5 = f64(Float(4.838803174e-8) * Float(100.0))
BI4 = f64(Float(5.824720280e-6) * Float(100.0))
BI3 = f64(Float(4.176223716e-4) * Float(100.0))
BI2 = f64(Float(1.886013408e-2) * Float(100.0))
BI1 = f64(Float(5.034698970e-1) * Float(100.0))
BI0 = f64(Float(6.109177956e0) * Float(100.0))


def _saturation_formulation(
    formulation: SaturationFormulation,
    t: Float,
):
    if formulation == SaturationFormulation.Staars:
        TT = t - MAPL_TICE
        if TT < TSTARR1:
            EX = (
                TT
                * (TT * (TT * (TT * (TT * (TT * S16 + S15) + S14) + S13) + S12) + S11)
                + S10
            )
        elif TT >= TSTARR1 and TT < TSTARR2:
            W = (TSTARR2 - TT) / (TSTARR2 - TSTARR1)
            EX = W * (
                TT
                * (TT * (TT * (TT * (TT * (TT * S16 + S15) + S14) + S13) + S12) + S11)
                + S10
            ) + (Float(1.0) - W) * (
                TT
                * (TT * (TT * (TT * (TT * (TT * S26 + S25) + S24) + S23) + S22) + S21)
                + S20
            )
        elif TT >= TSTARR2 and TT < TSTARR3:
            EX = (
                TT
                * (TT * (TT * (TT * (TT * (TT * S26 + S25) + S24) + S23) + S22) + S21)
                + S20
            )
        elif TT >= TSTARR3 and TT < TSTARR4:
            W = (TSTARR4 - TT) / (TSTARR4 - TSTARR3)
            EX = W * (
                TT
                * (TT * (TT * (TT * (TT * (TT * S26 + S25) + S24) + S23) + S22) + S21)
                + S20
            ) + (Float(1.0) - W) * (
                TT
                * (TT * (TT * (TT * (TT * (TT * BI6 + BI5) + BI4) + BI3) + BI2) + BI1)
                + BI0
            )
        else:
            EX = (
                TT
                * (TT * (TT * (TT * (TT * (TT * BI6 + BI5) + BI4) + BI3) + BI2) + BI1)
                + BI0
            )
    elif formulation == SaturationFormulation.CAM:
        TT = MAPL_TICE / t
        EX = DI[0] * np.exp(-(DI[1] / TT + DI[2] * np.log(TT) + DI[3] * TT))
    elif formulation == SaturationFormulation.MurphyAndKoop:
        EX = np.exp(CI[0] + CI[1] / t + CI[2] * np.log(t) + CI[3] * t)
    return Float(EX)


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

    DX = Float(0.0)  # only calulcated when DQ is not none
    EX = _saturation_formulation(formulation, TI)

    if DQ is not None:
        if temperature < TMINICE:
            DDQ = Float(0.0)
        elif temperature > MAPL_TICE:
            DDQ = Float(0.0)
        else:
            if PL > EX:
                DD = EX
                TI = temperature + DELTA_T
                EX, _ = _saturation_formulation(formulation, TI)
                DDQ = EX - DD
                EX = DD
    if PL is not None:
        if PL > EX:
            DD = ESFAC / (PL - (Float(1.0) - ESFAC) * EX)
            EX = EX * DD
            if DQ is not None:
                DX = DDQ * ERFAC * PL * DD * DD
        else:
            EX = MAX_MIXING_RATIO
            if DQ is not None:
                DX = Float(0.0)
    else:
        if DQ is not None:
            DX = DDQ * (Float(1.0) / DELTA_T)

    return EX, TI, DX
