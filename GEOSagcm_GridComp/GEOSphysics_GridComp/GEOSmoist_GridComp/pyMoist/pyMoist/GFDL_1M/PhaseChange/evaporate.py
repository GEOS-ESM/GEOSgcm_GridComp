from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_incloud_processes import cloud_effective_radius_liquid


def evaporate(
    PLmb: FloatField,
    T: FloatField,
    Q: FloatField,
    QLCN: FloatField,
    QICN: FloatField,
    CLCN: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    QST: FloatField,
    EVAPC: FloatField,
):
    from __externals__ import DT_MOIST, CCW_EVAP_EFF

    with computation(PARALLEL), interval(...):
        EVAPC = Q
        RHCRIT = 1
        # Evaporation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        ES = (
            100.0 * PLmb * QST / (constants.EPSILON + (1.0 - constants.EPSILON) * QST)
        )  # (100's <-^ convert from mbar to Pa)
        RHx = min(Q / QST, 1.00)
        K1 = (
            (constants.MAPL_LATENT_HEAT_VAPORIZATION**2)
            * constants.RHO_W
            / (constants.K_COND * constants.MAPL_RVAP * (T**2))
        )
        K2 = (
            constants.MAPL_RVAP
            * T
            * constants.RHO_W
            / (constants.DIFFU * (1000.0 / PLmb) * ES)
        )
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if CLCN > 0.0 and QLCN > 0.0:
            QCm = QLCN / CLCN
        else:
            QCm = 0.0
        RADIUS = cloud_effective_radius_liquid(PLmb, T, QCm, NACTL, NACTI)
        if RHx < RHCRIT and RADIUS > 0.0:
            EVAP = (
                CCW_EVAP_EFF
                * QLCN
                * DT_MOIST
                * (RHCRIT - RHx)
                / ((K1 + K2) * RADIUS**2)
            )
            EVAP = min(EVAP, QLCN)
        else:
            EVAP = 0.0
        QC = QLCN + QICN
        if QC > 0.0:
            CLCN = CLCN * (QC - EVAP) / QC
        Q = Q + EVAP
        QLCN = QLCN - EVAP
        T = T - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CPDRY) * EVAP
        EVAPC = (Q - EVAPC) / DT_MOIST
