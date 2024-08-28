from ndsl.dsl.typing import Float
from typing import Optional
from pyMoist.saturation.constants import (
    TMINLQU,
    DELTA_T,
    MAPL_TICE,
    ESFAC,
    MAX_MIXING_RATIO,
    TMAXTBL,
    DEGSUBS,
)
from pyMoist.saturation.formulation import SaturationFormulation
import numpy as np


ERFAC = DEGSUBS / ESFAC


# Below are actual 64-bit float in Fortran
B6 = 6.136820929e-11 * 100.0
B5 = 2.034080948e-8 * 100.0
B4 = 3.031240396e-6 * 100.0
B3 = 2.650648471e-4 * 100.0
B2 = 1.428945805e-2 * 100.0
B1 = 4.436518521e-1 * 100.0
B0 = 6.107799961e0 * 100.0
BI6 = 1.838826904e-10 * 100.0
BI5 = 4.838803174e-8 * 100.0
BI4 = 5.824720280e-6 * 100.0
BI3 = 4.176223716e-4 * 100.0
BI2 = 1.886013408e-2 * 100.0
BI1 = 5.034698970e-1 * 100.0
BI0 = 6.109177956e0 * 100.0
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

DL = [-7.902980, 5.02808, -1.3816, 11.344, 8.1328, -3.49149]
TS = 373.16
LOGPS = 3.005714898  # log10(1013.246)
CL = [
    54.842763,
    -6763.22,
    -4.21000,
    0.000367,
    0.0415,
    218.8,
    53.878000,
    -1331.22,
    -9.44523,
    0.014025,
]


def _saturation_formulation(formulation: SaturationFormulation, t: Float):
    if formulation == SaturationFormulation.Staars:
        TT = t - MAPL_TICE
        EX = TT * (TT * (TT * (TT * (TT * (TT * B6 + B5) + B4) + B3) + B2) + B1) + B0
    elif formulation == SaturationFormulation.CAM:
        TT = TS / t
        EX = 10.0 ** (
            DL[0] * (TT - 1.0)
            + DL[1] * np.log10(TT)
            + DL[2] * (10.0 ** (DL[3] * (1.0 - (1.0 / TT))) - 1.0) / 10000000.0
            + DL[4] * (10.0 ** (DL[5] * (TT - 1.0)) - 1.0) / 1000.0
            + LOGPS
            + 2.0
        )
    elif formulation == SaturationFormulation.MurphyAndKoop:
        EX = np.exp(
            (CL[0] + CL[1] / t + CL[2] * np.log(t) + CL[3] * t)
            + np.tanh(CL[4] * (t - CL[5]))
            * (CL[6] + CL[7] / t + CL[8] * np.log(t) + CL[9] * t)
        )
    return EX


def qsat_liquid_scalar_exact(
    temperature: Float,
    formulation: SaturationFormulation = SaturationFormulation.Staars,
    PL: Optional[Float] = None,
    DQ: Optional[Float] = None,
):
    """Reference Fortran: QSATLQU0 w/ UTBL=False"""

    if temperature < TMINLQU:
        TI = TMINLQU
    elif temperature > TMAXTBL:
        TI = TMAXTBL
    else:
        TI = temperature

    DX = 0  # DX only calculated when DQ is present
    EX = _saturation_formulation(formulation, TI)

    if DQ is not None:
        if temperature < TMINLQU:
            DDQ = 0.0
        elif temperature > TMAXTBL:
            DDQ = 0.0
        else:
            if PL > EX:
                DD = EX
                TI = temperature + DELTA_T
                EX = _saturation_formulation(formulation, TI)
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
    elif DQ is not None:
        DX = DDQ * (1.0 / DELTA_T)

    return EX, TI, DX


def qsat_liquid_scalar_table():
    """Reference Fortran: QSATLQU0 w/ UTBL=True"""
    raise NotImplementedError("Nope.")
