from typing import Optional

import numpy as np

from ndsl.dsl.typing import Float, Int
from pyMoist.saturation_tables.constants import (
    DELTA_T,
    ERFAC,
    ESFAC,
    MAPL_TICE,
    MAX_MIXING_RATIO,
    TMAXTBL,
    TMINLQU,
)
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.constants import LiquidExactConstants
from ndsl.dsl.gt4py import function, exp, tanh, log10, log


@function
def _saturation_formulation(
    t: Float,
    formulation: Int,
    B6: Float,
    B5: Float,
    B4: Float,
    B3: Float,
    B2: Float,
    B1: Float,
    B0: Float,
    BI6: Float,
    BI5: Float,
    BI4: Float,
    BI3: Float,
    BI2: Float,
    BI1: Float,
    BI0: Float,
    S16: Float,
    S15: Float,
    S14: Float,
    S13: Float,
    S12: Float,
    S11: Float,
    S10: Float,
    S26: Float,
    S25: Float,
    S24: Float,
    S23: Float,
    S22: Float,
    S21: Float,
    S20: Float,
    DL_0: Float,
    DL_1: Float,
    DL_2: Float,
    DL_3: Float,
    DL_4: Float,
    DL_5: Float,
    TS: Float,
    LOGPS: Float,
    CL_0: Float,
    CL_1: Float,
    CL_2: Float,
    CL_3: Float,
    CL_4: Float,
    CL_5: Float,
    CL_6: Float,
    CL_7: Float,
    CL_8: Float,
    CL_9: Float,
):
    if formulation == 1:
        tt = t - MAPL_TICE
        ex = tt * (tt * (tt * (tt * (tt * (tt * B6 + B5) + B4) + B3) + B2) + B1) + B0
    elif formulation == 2:
        tt = TS / t
        ex = 10.0 ** (
            DL_0 * (tt - 1.0)
            + DL_1 * log10(tt)
            + DL_2 * (10.0 ** (DL_3 * (1.0 - (1.0 / tt))) - 1.0) / 10000000.0
            + DL_4 * (10.0 ** (DL_5 * (tt - 1.0)) - 1.0) / 1000.0
            + LOGPS
            + 2.0
        )
    elif formulation == 3:
        ex = exp(
            (CL_0 + CL_1 / t + CL_2 * log(t) + CL_3 * t)
            + tanh(CL_4 * (t - CL_5)) * (CL_6 + CL_7 / t + CL_8 * log(t) + CL_9 * t)
        )
    return ex


def _saturation_formulation_no_stencil(
    t: Float, formulation: SaturationFormulation = SaturationFormulation.Staars
):
    if formulation == SaturationFormulation.Staars:
        tt = t - MAPL_TICE
        ex = (
            tt
            * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt * (tt * LiquidExactConstants.B6 + LiquidExactConstants.B5)
                            + LiquidExactConstants.B4
                        )
                        + LiquidExactConstants.B3
                    )
                    + LiquidExactConstants.B2
                )
                + LiquidExactConstants.B1
            )
            + LiquidExactConstants.B0
        )
    elif formulation == SaturationFormulation.CAM:
        tt = LiquidExactConstants.TS / t
        ex = Float(10.0) ** (
            LiquidExactConstants.DL[0] * (tt - Float(1.0))
            + LiquidExactConstants.DL[1] * np.log10(tt)
            + LiquidExactConstants.DL[2]
            * (Float(10.0) ** (LiquidExactConstants.DL[3] * (Float(1.0) - (Float(1.0) / tt))) - Float(1.0))
            / Float(10000000.0)
            + LiquidExactConstants.DL[4]
            * (Float(10.0) ** (LiquidExactConstants.DL[5] * (tt - Float(1.0))) - Float(1.0))
            / Float(1000.0)
            + LiquidExactConstants.LOGPS
            + Float(2.0)
        )
    elif formulation == SaturationFormulation.MurphyAndKoop:
        ex = np.exp(
            (
                LiquidExactConstants.CL[0]
                + LiquidExactConstants.CL[1] / t
                + LiquidExactConstants.CL[2] * np.log(t)
                + LiquidExactConstants.CL[3] * t
            )
            + np.tanh(LiquidExactConstants.CL[4] * (t - LiquidExactConstants.CL[5]))
            * (
                LiquidExactConstants.CL[6]
                + LiquidExactConstants.CL[7] / t
                + LiquidExactConstants.CL[8] * np.log(t)
                + LiquidExactConstants.CL[9] * t
            )
        )
    return Float(ex)


@function
def liquid_exact(
    t_in: Float,
    formulation: Int,
    B6: Float,
    B5: Float,
    B4: Float,
    B3: Float,
    B2: Float,
    B1: Float,
    B0: Float,
    BI6: Float,
    BI5: Float,
    BI4: Float,
    BI3: Float,
    BI2: Float,
    BI1: Float,
    BI0: Float,
    S16: Float,
    S15: Float,
    S14: Float,
    S13: Float,
    S12: Float,
    S11: Float,
    S10: Float,
    S26: Float,
    S25: Float,
    S24: Float,
    S23: Float,
    S22: Float,
    S21: Float,
    S20: Float,
    DL_0: Float,
    DL_1: Float,
    DL_2: Float,
    DL_3: Float,
    DL_4: Float,
    DL_5: Float,
    TS: Float,
    LOGPS: Float,
    CL_0: Float,
    CL_1: Float,
    CL_2: Float,
    CL_3: Float,
    CL_4: Float,
    CL_5: Float,
    CL_6: Float,
    CL_7: Float,
    CL_8: Float,
    CL_9: Float,
    p: Float = 1e15,
    dq: Float = 1e15,
):
    """Reference Fortran: QSATLQU0 w/ UTBL=False"""

    if t_in < TMINLQU:
        t = TMINLQU
    elif t_in > TMAXTBL:
        t = TMAXTBL
    else:
        t = t_in

    dx = 0.0  # DX only calculated when DQ is present
    ex = _saturation_formulation(
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
    )

    if dq != 1e15:
        if t_in < TMINLQU:
            ddq = 0.0
        elif t_in > TMAXTBL:
            ddq = 0.0
        else:
            if p > ex:
                dd = ex
                t = t_in + DELTA_T
                ex = _saturation_formulation(
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
                )
                ddq = ex - dd
                ex = dd

    if p != 1e15:
        if p > ex:
            dd = ESFAC / (p - (1.0 - ESFAC) * ex)
            ex = ex * dd
            if dq != 1e15:
                dx = ddq * ERFAC * p * dd * dd
        else:
            ex = MAX_MIXING_RATIO
            if dq != 1e15:
                dx = 0.0
    elif dq != 1e15:
        dx = ddq * (1.0 / DELTA_T)

    return ex, dx


def liquid_exact_no_stencil(
    t_in: Float,
    formulation: SaturationFormulation = SaturationFormulation.Staars,
    p: Optional[Float] = None,
    dq: Optional[Float] = None,
):
    """Reference Fortran: QSATLQU0 w/ UTBL=False"""

    if t_in < TMINLQU:
        t = TMINLQU
    elif t_in > TMAXTBL:
        t = TMAXTBL
    else:
        t = t_in

    dx = Float(0.0)  # DX only calculated when DQ is present
    ex = _saturation_formulation_no_stencil(t=t, formulation=formulation)

    if dq is not None:
        if t_in < TMINLQU:
            ddq = Float(0.0)
        elif t_in > TMAXTBL:
            ddq = Float(0.0)
        else:
            if p > ex:
                dd = ex
                t = t_in + DELTA_T
                ex, _ = _saturation_formulation_no_stencil(t=t, formulation=formulation)
                ddq = ex - dd
                ex = dd

    if p is not None:
        if p > ex:
            dd = ESFAC / (p - (Float(1.0) - ESFAC) * ex)
            ex = ex * dd
            if dq is not None:
                dx = ddq * ERFAC * p * dd * dd
        else:
            ex = MAX_MIXING_RATIO
            if dq is not None:
                dx = Float(0.0)
    elif dq is not None:
        dx = ddq * (Float(1.0) / DELTA_T)

    return ex, dx
