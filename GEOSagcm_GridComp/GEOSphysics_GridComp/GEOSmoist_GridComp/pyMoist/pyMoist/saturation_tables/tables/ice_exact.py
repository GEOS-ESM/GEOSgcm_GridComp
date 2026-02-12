from typing import Optional

import numpy as np
from ndsl.dsl.gt4py import exp, function, log
from ndsl.dsl.typing import Float, Int

from pyMoist.saturation_tables.constants import DELTA_T, ERFAC, ESFAC, MAPL_TICE, MAX_MIXING_RATIO
from pyMoist.saturation_tables.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.constants import IceExactConstants


@function
def _saturation_formulation(
    t: Float,
    formulation: Int,
    TMINSTR: Float,
    TMINICE: Float,
    TSTARR1: Float,
    TSTARR2: Float,
    TSTARR3: Float,
    TSTARR4: Float,
    TMAXSTR: Float,
    DI_0: Float,
    DI_1: Float,
    DI_2: Float,
    DI_3: Float,
    CI_0: Float,
    CI_1: Float,
    CI_2: Float,
    CI_3: Float,
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
    BI6: Float,
    BI5: Float,
    BI4: Float,
    BI3: Float,
    BI2: Float,
    BI1: Float,
    BI0: Float,
):
    if formulation == 1:
        tt = t - MAPL_TICE
        if tt < TSTARR1:
            ex = tt * (tt * (tt * (tt * (tt * (tt * S16 + S15) + S14) + S13) + S12) + S11) + S10
        elif tt >= TSTARR1 and tt < TSTARR2:
            W = (TSTARR2 - tt) / (TSTARR2 - TSTARR1)
            ex = W * (tt * (tt * (tt * (tt * (tt * (tt * S16 + S15) + S14) + S13) + S12) + S11) + S10) + (
                1.0 - W
            ) * (tt * (tt * (tt * (tt * (tt * (tt * S26 + S25) + S24) + S23) + S22) + S21) + S20)
        elif tt >= TSTARR2 and tt < TSTARR3:
            ex = tt * (tt * (tt * (tt * (tt * (tt * S26 + S25) + S24) + S23) + S22) + S21) + S20
        elif tt >= TSTARR3 and tt < TSTARR4:
            W = (TSTARR4 - tt) / (TSTARR4 - TSTARR3)
            ex = W * (tt * (tt * (tt * (tt * (tt * (tt * S26 + S25) + S24) + S23) + S22) + S21) + S20) + (
                1.0 - W
            ) * (tt * (tt * (tt * (tt * (tt * (tt * BI6 + BI5) + BI4) + BI3) + BI2) + BI1) + BI0)
        else:
            ex = tt * (tt * (tt * (tt * (tt * (tt * BI6 + BI5) + BI4) + BI3) + BI2) + BI1) + BI0
    elif formulation == 2:
        tt = MAPL_TICE / t
        ex = DI_0 * exp(-(DI_1 / tt + DI_2 * log(tt) + DI_3 * tt))
    elif formulation == 3:
        ex = exp(CI_0 + CI_1 / t + CI_2 * log(t) + CI_3 * t)
    return ex


def _saturation_formulation_no_stencil(
    t: Float,
    formulation: SaturationFormulation,
):
    if formulation == SaturationFormulation.Staars:
        tt = t - MAPL_TICE
        if tt < IceExactConstants.TSTARR1:
            ex = (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.S16 + IceExactConstants.S15)
                                + IceExactConstants.S14
                            )
                            + IceExactConstants.S13
                        )
                        + IceExactConstants.S12
                    )
                    + IceExactConstants.S11
                )
                + IceExactConstants.S10
            )
        elif tt >= IceExactConstants.TSTARR1 and tt < IceExactConstants.TSTARR2:
            W = (IceExactConstants.TSTARR2 - tt) / (IceExactConstants.TSTARR2 - IceExactConstants.TSTARR1)
            ex = W * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.S16 + IceExactConstants.S15)
                                + IceExactConstants.S14
                            )
                            + IceExactConstants.S13
                        )
                        + IceExactConstants.S12
                    )
                    + IceExactConstants.S11
                )
                + IceExactConstants.S10
            ) + (Float(1.0) - W) * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.S26 + IceExactConstants.S25)
                                + IceExactConstants.S24
                            )
                            + IceExactConstants.S23
                        )
                        + IceExactConstants.S22
                    )
                    + IceExactConstants.S21
                )
                + IceExactConstants.S20
            )
        elif tt >= IceExactConstants.TSTARR2 and tt < IceExactConstants.TSTARR3:
            ex = (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.S26 + IceExactConstants.S25)
                                + IceExactConstants.S24
                            )
                            + IceExactConstants.S23
                        )
                        + IceExactConstants.S22
                    )
                    + IceExactConstants.S21
                )
                + IceExactConstants.S20
            )
        elif tt >= IceExactConstants.TSTARR3 and tt < IceExactConstants.TSTARR4:
            W = (IceExactConstants.TSTARR4 - tt) / (IceExactConstants.TSTARR4 - IceExactConstants.TSTARR3)
            ex = W * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.S26 + IceExactConstants.S25)
                                + IceExactConstants.S24
                            )
                            + IceExactConstants.S23
                        )
                        + IceExactConstants.S22
                    )
                    + IceExactConstants.S21
                )
                + IceExactConstants.S20
            ) + (Float(1.0) - W) * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.BI6 + IceExactConstants.BI5)
                                + IceExactConstants.BI4
                            )
                            + IceExactConstants.BI3
                        )
                        + IceExactConstants.BI2
                    )
                    + IceExactConstants.BI1
                )
                + IceExactConstants.BI0
            )
        else:
            ex = (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt * (tt * IceExactConstants.BI6 + IceExactConstants.BI5)
                                + IceExactConstants.BI4
                            )
                            + IceExactConstants.BI3
                        )
                        + IceExactConstants.BI2
                    )
                    + IceExactConstants.BI1
                )
                + IceExactConstants.BI0
            )
    elif formulation == SaturationFormulation.CAM:
        tt = MAPL_TICE / t
        ex = IceExactConstants.DI[0] * np.exp(
            -(
                IceExactConstants.DI[1] / tt
                + IceExactConstants.DI[2] * np.log(tt)
                + IceExactConstants.DI[3] * tt
            )
        )
    elif formulation == SaturationFormulation.MurphyAndKoop:
        ex = np.exp(
            IceExactConstants.CI[0]
            + IceExactConstants.CI[1] / t
            + IceExactConstants.CI[2] * np.log(t)
            + IceExactConstants.CI[3] * t
        )
    return Float(ex)


@function
def ice_exact(
    t_in: Float,
    formulation: Int,
    TMINSTR: Float,
    TMINICE: Float,
    TSTARR1: Float,
    TSTARR2: Float,
    TSTARR3: Float,
    TSTARR4: Float,
    TMAXSTR: Float,
    DI_0: Float,
    DI_1: Float,
    DI_2: Float,
    DI_3: Float,
    CI_0: Float,
    CI_1: Float,
    CI_2: Float,
    CI_3: Float,
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
    BI6: Float,
    BI5: Float,
    BI4: Float,
    BI3: Float,
    BI2: Float,
    BI1: Float,
    BI0: Float,
    p: Float = 1e15,
    dq: Float = 1e15,
):
    """Reference Fortran: QSATICE0 w/ UTBL=False"""
    if t_in < TMINICE:
        t = TMINICE
    elif t_in > MAPL_TICE:
        t = MAPL_TICE
    else:
        t = t_in

    dx = 0.0  # only calculated when DQ is not none
    ex = _saturation_formulation(
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

    if dq != 1e15:
        if t_in < TMINICE:
            ddq = 0.0
        elif t_in > MAPL_TICE:
            ddq = 0.0
        else:
            if p > ex:
                dd = ex
                t = t_in + DELTA_T
                ex = _saturation_formulation(
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
    else:
        if dq != 1e15:
            dx = ddq * (1.0 / DELTA_T)

    return ex, dx


def ice_exact_no_stencil(
    t_in: Float,
    formulation: SaturationFormulation = SaturationFormulation.Staars,
    p: Optional[Float] = None,
    dq: Optional[Float] = None,
):
    """Reference Fortran: QSATICE0 w/ UTBL=False"""
    if t_in < IceExactConstants.TMINICE:
        t = IceExactConstants.TMINICE
    elif t_in > MAPL_TICE:
        t = MAPL_TICE
    else:
        t = t_in

    dx = Float(0.0)  # only calculated when DQ is not none
    ex = _saturation_formulation_no_stencil(t=t, formulation=formulation)

    if dq is not None:
        if t_in < IceExactConstants.TMINICE:
            ddq = Float(0.0)
        elif t_in > MAPL_TICE:
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
    else:
        if dq is not None:
            dx = ddq * (Float(1.0) / DELTA_T)

    return ex, dx
