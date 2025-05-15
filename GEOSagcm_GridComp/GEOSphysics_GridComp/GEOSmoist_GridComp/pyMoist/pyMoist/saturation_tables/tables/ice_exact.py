from typing import Optional

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation_old.constants import (
    DELTA_T,
    ERFAC,
    ESFAC,
    MAPL_TICE,
    MAX_MIXING_RATIO,
)
from pyMoist.saturation_old.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.constants import IceExactConstatns


def _saturation_formulation(
    t: Float,
    formulation: SaturationFormulation,
):
    if formulation == SaturationFormulation.Staars:
        tt = t - MAPL_TICE
        if tt < IceExactConstatns.TSTARR1:
            ex = (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt
                                * (tt * IceExactConstatns.S16 + IceExactConstatns.S15)
                                + IceExactConstatns.S14
                            )
                            + IceExactConstatns.S13
                        )
                        + IceExactConstatns.S12
                    )
                    + IceExactConstatns.S11
                )
                + IceExactConstatns.S10
            )
        elif tt >= IceExactConstatns.TSTARR1 and tt < IceExactConstatns.TSTARR2:
            W = (IceExactConstatns.TSTARR2 - tt) / (
                IceExactConstatns.TSTARR2 - IceExactConstatns.TSTARR1
            )
            ex = W * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt
                                * (tt * IceExactConstatns.S16 + IceExactConstatns.S15)
                                + IceExactConstatns.S14
                            )
                            + IceExactConstatns.S13
                        )
                        + IceExactConstatns.S12
                    )
                    + IceExactConstatns.S11
                )
                + IceExactConstatns.S10
            ) + (Float(1.0) - W) * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt
                                * (tt * IceExactConstatns.S26 + IceExactConstatns.S25)
                                + IceExactConstatns.S24
                            )
                            + IceExactConstatns.S23
                        )
                        + IceExactConstatns.S22
                    )
                    + IceExactConstatns.S21
                )
                + IceExactConstatns.S20
            )
        elif tt >= IceExactConstatns.TSTARR2 and tt < IceExactConstatns.TSTARR3:
            ex = (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt
                                * (tt * IceExactConstatns.S26 + IceExactConstatns.S25)
                                + IceExactConstatns.S24
                            )
                            + IceExactConstatns.S23
                        )
                        + IceExactConstatns.S22
                    )
                    + IceExactConstatns.S21
                )
                + IceExactConstatns.S20
            )
        elif tt >= IceExactConstatns.TSTARR3 and tt < IceExactConstatns.TSTARR4:
            W = (IceExactConstatns.TSTARR4 - tt) / (
                IceExactConstatns.TSTARR4 - IceExactConstatns.TSTARR3
            )
            ex = W * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt
                                * (tt * IceExactConstatns.S26 + IceExactConstatns.S25)
                                + IceExactConstatns.S24
                            )
                            + IceExactConstatns.S23
                        )
                        + IceExactConstatns.S22
                    )
                    + IceExactConstatns.S21
                )
                + IceExactConstatns.S20
            ) + (Float(1.0) - W) * (
                tt
                * (
                    tt
                    * (
                        tt
                        * (
                            tt
                            * (
                                tt
                                * (tt * IceExactConstatns.BI6 + IceExactConstatns.BI5)
                                + IceExactConstatns.BI4
                            )
                            + IceExactConstatns.BI3
                        )
                        + IceExactConstatns.BI2
                    )
                    + IceExactConstatns.BI1
                )
                + IceExactConstatns.BI0
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
                                tt
                                * (tt * IceExactConstatns.BI6 + IceExactConstatns.BI5)
                                + IceExactConstatns.BI4
                            )
                            + IceExactConstatns.BI3
                        )
                        + IceExactConstatns.BI2
                    )
                    + IceExactConstatns.BI1
                )
                + IceExactConstatns.BI0
            )
    elif formulation == SaturationFormulation.CAM:
        tt = MAPL_TICE / t
        ex = IceExactConstatns.DI[0] * np.exp(
            -(
                IceExactConstatns.DI[1] / tt
                + IceExactConstatns.DI[2] * np.log(tt)
                + IceExactConstatns.DI[3] * tt
            )
        )
    elif formulation == SaturationFormulation.MurphyAndKoop:
        ex = np.exp(
            IceExactConstatns.CI[0]
            + IceExactConstatns.CI[1] / t
            + IceExactConstatns.CI[2] * np.log(t)
            + IceExactConstatns.CI[3] * t
        )
    return Float(ex)


def ice_exact(
    t_in: Float,
    formulation: SaturationFormulation = SaturationFormulation.Staars,
    p: Optional[Float] = None,
    dq: Optional[Float] = None,
):
    """Reference Fortran: QSATICE0 w/ UTBL=False"""
    if t_in < IceExactConstatns.TMINICE:
        t = IceExactConstatns.TMINICE
    elif t_in > MAPL_TICE:
        t = MAPL_TICE
    else:
        t = t_in

    dx = Float(0.0)  # only calulcated when DQ is not none
    ex = _saturation_formulation(t=t, formulation=formulation)

    if dq is not None:
        if t_in < IceExactConstatns.TMINICE:
            ddq = Float(0.0)
        elif t_in > MAPL_TICE:
            ddq = Float(0.0)
        else:
            if p > ex:
                dd = ex
                t = t_in + DELTA_T
                ex, _ = _saturation_formulation(t=t, formulation=formulation)
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
