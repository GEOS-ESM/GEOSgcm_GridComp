from typing import Optional

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation.constants import DELTA_T, ERFAC, ESFAC, MAPL_TICE, MAX_MIXING_RATIO, TMAXTBL, TMINLQU
from pyMoist.saturation.formulation import SaturationFormulation


f64 = np.float64

# Below are actual 64-bit float in Fortran
B6 = f64(Float(6.136820929e-11) * Float(100.0))
B5 = f64(Float(2.034080948e-8) * Float(100.0))
B4 = f64(Float(3.031240396e-6) * Float(100.0))
B3 = f64(Float(2.650648471e-4) * Float(100.0))
B2 = f64(Float(1.428945805e-2) * Float(100.0))
B1 = f64(Float(4.436518521e-1) * Float(100.0))
B0 = f64(Float(6.107799961e0) * Float(100.0))
BI6 = f64(Float(1.838826904e-10) * Float(100.0))
BI5 = f64(Float(4.838803174e-8) * Float(100.0))
BI4 = f64(Float(5.824720280e-6) * Float(100.0))
BI3 = f64(Float(4.176223716e-4) * Float(100.0))
BI2 = f64(Float(1.886013408e-2) * Float(100.0))
BI1 = f64(Float(5.034698970e-1) * Float(100.0))
BI0 = f64(Float(6.109177956e0) * Float(100.0))
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

DL = [
    f64(Float(-7.902980)),
    f64(Float(5.02808)),
    f64(Float(-1.3816)),
    f64(Float(11.344)),
    f64(Float(8.1328)),
    f64(Float(-3.49149)),
]
TS = f64(Float(373.16))
LOGPS = f64(Float(3.005714898))  # log10(1013.246)
CL = [
    f64(Float(54.842763)),
    f64(Float(-6763.22)),
    f64(Float(-4.21000)),
    f64(Float(0.000367)),
    f64(Float(0.0415)),
    f64(Float(218.8)),
    f64(Float(53.878000)),
    f64(Float(-1331.22)),
    f64(Float(-9.44523)),
    f64(Float(0.014025)),
]


def _saturation_formulation(formulation: SaturationFormulation, t: Float):
    if formulation == SaturationFormulation.Staars:
        TT = t - MAPL_TICE
        EX = TT * (TT * (TT * (TT * (TT * (TT * B6 + B5) + B4) + B3) + B2) + B1) + B0
    elif formulation == SaturationFormulation.CAM:
        TT = TS / t
        EX = Float(10.0) ** (
            DL[0] * (TT - Float(1.0))
            + DL[1] * np.log10(TT)
            + DL[2]
            * (Float(10.0) ** (DL[3] * (Float(1.0) - (Float(1.0) / TT))) - Float(1.0))
            / Float(10000000.0)
            + DL[4]
            * (Float(10.0) ** (DL[5] * (TT - Float(1.0))) - Float(1.0))
            / Float(1000.0)
            + LOGPS
            + Float(2.0)
        )
    elif formulation == SaturationFormulation.MurphyAndKoop:
        EX = np.exp(
            (CL[0] + CL[1] / t + CL[2] * np.log(t) + CL[3] * t)
            + np.tanh(CL[4] * (t - CL[5]))
            * (CL[6] + CL[7] / t + CL[8] * np.log(t) + CL[9] * t)
        )
    return Float(EX)


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

    DX = Float(0.0)  # DX only calculated when DQ is present
    EX = _saturation_formulation(formulation, TI)

    if DQ is not None:
        if temperature < TMINLQU:
            DDQ = Float(0.0)
        elif temperature > TMAXTBL:
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
    elif DQ is not None:
        DX = DDQ * (Float(1.0) / DELTA_T)

    return EX, TI, DX
