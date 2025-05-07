import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_incloud_processes import cloud_effective_radius_ice


def sublimate(
    PLmb: FloatField,
    T: FloatField,
    Q: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    QST: FloatField,
    SUBLC: FloatField,
):
    from __externals__ import DT_MOIST, CCI_EVAP_EFF

    with computation(PARALLEL), interval(...):
        SUBLC = Q
        RHCRIT = 1
        # Sublimation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        ES = (
            100.0 * PLmb * QST / (constants.EPSILON + (1.0 - constants.EPSILON) * QST)
        )  # (100s <-^ convert from mbar to Pa)
        RHx = min(Q / QST, 1.00)
        K1 = (
            (constants.MAPL_LATENT_HEAT_VAPORIZATION**2)
            * constants.RHO_I
            / (constants.K_COND * constants.MAPL_RVAP * (T**2))
        )
        K2 = (
            constants.MAPL_RVAP
            * T
            * constants.RHO_I
            / (constants.DIFFU * (1000.0 / PLmb) * ES)
        )
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if CLCN > 0.0 and QICN > 0.0:
            QCm = QICN / CLCN
        else:
            QCm = 0.0
        radius = cloud_effective_radius_ice(PLmb, T, QCm, NACTL, NACTI)
        if RHx < RHCRIT and radius > 0.0:
            SUBL = (
                CCI_EVAP_EFF
                * QICN
                * DT_MOIST
                * (RHCRIT - RHx)
                / ((K1 + K2) * radius**2)
            )
            SUBL = min(SUBL, QICN)
        else:
            SUBL = 0.0
        QC = QLCN + QICN
        if QC > 0.0:
            CLCN = CLCN * (QC - SUBL) / QC
        Q = Q + SUBL
        QICN = QICN - SUBL
        T = T - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CPDRY) * SUBL
        SUBLC = (Q - SUBLC) / DT_MOIST
