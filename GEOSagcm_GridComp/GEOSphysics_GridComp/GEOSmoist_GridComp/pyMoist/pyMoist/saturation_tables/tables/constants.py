from dataclasses import dataclass

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation_tables.constants import MAPL_TICE


@dataclass
class IceExactConstatns:
    TMINSTR = Float(-95.0)
    TMINICE = MAPL_TICE + TMINSTR

    TSTARR1 = Float(-75.0)
    TSTARR2 = Float(-65.0)
    TSTARR3 = Float(-50.0)
    TSTARR4 = Float(-40.0)
    TMAXSTR = Float(+60.0)

    DI = [
        np.float64(Float(57518.5606e08)),
        np.float64(Float(2.01889049)),
        np.float64(Float(3.56654)),
        np.float64(Float(20.947031)),
    ]
    CI = [
        np.float64(Float(9.550426)),
        np.float64(Float(-5723.265)),
        np.float64(Float(3.53068)),
        np.float64(Float(-0.00728332)),
    ]

    # 64-bit float in Fortran
    S16 = np.float64(Float(0.516000335e-11) * Float(100.0))
    S15 = np.float64(Float(0.276961083e-8) * Float(100.0))
    S14 = np.float64(Float(0.623439266e-6) * Float(100.0))
    S13 = np.float64(Float(0.754129933e-4) * Float(100.0))
    S12 = np.float64(Float(0.517609116e-2) * Float(100.0))
    S11 = np.float64(Float(0.191372282e0) * Float(100.0))
    S10 = np.float64(Float(0.298152339e1) * Float(100.0))
    S26 = np.float64(Float(0.314296723e-10) * Float(100.0))
    S25 = np.float64(Float(0.132243858e-7) * Float(100.0))
    S24 = np.float64(Float(0.236279781e-5) * Float(100.0))
    S23 = np.float64(Float(0.230325039e-3) * Float(100.0))
    S22 = np.float64(Float(0.129690326e-1) * Float(100.0))
    S21 = np.float64(Float(0.401390832e0) * Float(100.0))
    S20 = np.float64(Float(0.535098336e1) * Float(100.0))
    BI6 = np.float64(Float(1.838826904e-10) * Float(100.0))
    BI5 = np.float64(Float(4.838803174e-8) * Float(100.0))
    BI4 = np.float64(Float(5.824720280e-6) * Float(100.0))
    BI3 = np.float64(Float(4.176223716e-4) * Float(100.0))
    BI2 = np.float64(Float(1.886013408e-2) * Float(100.0))
    BI1 = np.float64(Float(5.034698970e-1) * Float(100.0))
    BI0 = np.float64(Float(6.109177956e0) * Float(100.0))


@dataclass
class LiquidExactConstants:
    # Below are actual 64-bit float in Fortran
    B6 = np.float64(Float(6.136820929e-11) * Float(100.0))
    B5 = np.float64(Float(2.034080948e-8) * Float(100.0))
    B4 = np.float64(Float(3.031240396e-6) * Float(100.0))
    B3 = np.float64(Float(2.650648471e-4) * Float(100.0))
    B2 = np.float64(Float(1.428945805e-2) * Float(100.0))
    B1 = np.float64(Float(4.436518521e-1) * Float(100.0))
    B0 = np.float64(Float(6.107799961e0) * Float(100.0))
    BI6 = np.float64(Float(1.838826904e-10) * Float(100.0))
    BI5 = np.float64(Float(4.838803174e-8) * Float(100.0))
    BI4 = np.float64(Float(5.824720280e-6) * Float(100.0))
    BI3 = np.float64(Float(4.176223716e-4) * Float(100.0))
    BI2 = np.float64(Float(1.886013408e-2) * Float(100.0))
    BI1 = np.float64(Float(5.034698970e-1) * Float(100.0))
    BI0 = np.float64(Float(6.109177956e0) * Float(100.0))
    S16 = np.float64(Float(0.516000335e-11) * Float(100.0))
    S15 = np.float64(Float(0.276961083e-8) * Float(100.0))
    S14 = np.float64(Float(0.623439266e-6) * Float(100.0))
    S13 = np.float64(Float(0.754129933e-4) * Float(100.0))
    S12 = np.float64(Float(0.517609116e-2) * Float(100.0))
    S11 = np.float64(Float(0.191372282e0) * Float(100.0))
    S10 = np.float64(Float(0.298152339e1) * Float(100.0))
    S26 = np.float64(Float(0.314296723e-10) * Float(100.0))
    S25 = np.float64(Float(0.132243858e-7) * Float(100.0))
    S24 = np.float64(Float(0.236279781e-5) * Float(100.0))
    S23 = np.float64(Float(0.230325039e-3) * Float(100.0))
    S22 = np.float64(Float(0.129690326e-1) * Float(100.0))
    S21 = np.float64(Float(0.401390832e0) * Float(100.0))
    S20 = np.float64(Float(0.535098336e1) * Float(100.0))

    DL = [
        np.float64(Float(-7.902980)),
        np.float64(Float(5.02808)),
        np.float64(Float(-1.3816)),
        np.float64(Float(11.344)),
        np.float64(Float(8.1328)),
        np.float64(Float(-3.49149)),
    ]
    TS = np.float64(Float(373.16))
    LOGPS = np.float64(Float(3.005714898))  # log10(1013.246)
    CL = [
        np.float64(Float(54.842763)),
        np.float64(Float(-6763.22)),
        np.float64(Float(-4.21000)),
        np.float64(Float(0.000367)),
        np.float64(Float(0.0415)),
        np.float64(Float(218.8)),
        np.float64(Float(53.878000)),
        np.float64(Float(-1331.22)),
        np.float64(Float(-9.44523)),
        np.float64(Float(0.014025)),
    ]
