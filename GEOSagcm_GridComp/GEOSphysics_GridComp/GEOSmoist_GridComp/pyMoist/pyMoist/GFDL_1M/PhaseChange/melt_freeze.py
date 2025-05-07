from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    exp,
    interval,
)

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.shared_incloud_processes import ice_fraction


def melt_freeze(
    dt: Float,
    cnv_frc: FloatFieldIJ,
    srf_type: FloatFieldIJ,
    T: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
):
    with computation(PARALLEL), interval(...):
        if T <= constants.MAPL_TICE:
            fQi = ice_fraction(T, cnv_frc, srf_type)
            dQil = QLCN * (1.0 - exp(-dt * fQi / constants.TAUFRZ))
            dQil = max(0.0, dQil)
            QICN = QICN + dQil
            QLCN = QLCN - dQil
            T = (
                T
                + (
                    constants.MAPL_LATENT_HEAT_SUBLIMATION
                    - constants.MAPL_LATENT_HEAT_VAPORIZATION
                )
                * dQil
                / constants.MAPL_CP
            )
        else:
            dQil = -QICN
            dQil = min(0.0, dQil)
            QICN = QICN + dQil
            QLCN = QLCN - dQil
            T = (
                T
                + (
                    constants.MAPL_LATENT_HEAT_SUBLIMATION
                    - constants.MAPL_LATENT_HEAT_VAPORIZATION
                )
                * dQil
                / constants.MAPL_CP
            )
