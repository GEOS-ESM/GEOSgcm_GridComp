from typing import Optional

import numpy as np

from ndsl.dsl.typing import Float
from pyMoist.saturation_old.constants import (
    DELTA_T,
    ERFAC,
    ESFAC,
    MAPL_TICE,
    MAX_MIXING_RATIO,
    TMAXTBL,
    TMINLQU,
)
from pyMoist.saturation_old.formulation import SaturationFormulation
from pyMoist.saturation_tables.tables.constants import LiquidExactConstants


def _saturation_formulation(
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
                            tt
                            * (tt * LiquidExactConstants.B6 + LiquidExactConstants.B5)
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
            * (
                Float(10.0)
                ** (LiquidExactConstants.DL[3] * (Float(1.0) - (Float(1.0) / tt)))
                - Float(1.0)
            )
            / Float(10000000.0)
            + LiquidExactConstants.DL[4]
            * (
                Float(10.0) ** (LiquidExactConstants.DL[5] * (tt - Float(1.0)))
                - Float(1.0)
            )
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


def liquid_exact(
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
    ex = _saturation_formulation(t=t, formulation=formulation)

    if dq is not None:
        if t_in < TMINLQU:
            ddq = Float(0.0)
        elif t_in > TMAXTBL:
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
    elif dq is not None:
        dx = ddq * (Float(1.0) / DELTA_T)

    return ex, dx
