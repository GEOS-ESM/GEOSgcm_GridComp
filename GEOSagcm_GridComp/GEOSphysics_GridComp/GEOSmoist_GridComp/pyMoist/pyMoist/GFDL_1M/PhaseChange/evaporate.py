from gt4py.cartesian.gtscript import (
    PARALLEL,
    computation,
    interval,
)

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_incloud_processes import cloud_effective_radius_liquid


def evaporate(
    p_mb: FloatField,
    t: FloatField,
    q: FloatField,
    qlcn: FloatField,
    qicn: FloatField,
    clcn: FloatField,
    nactl: FloatField,
    nacti: FloatField,
    qst: FloatField,
    evapc: FloatField,
):
    from __externals__ import DT_MOIST, CCW_EVAP_EFF

    with computation(PARALLEL), interval(...):
        evapc = q
        rh_crit = 1
        # Evaporation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        es = (
            100.0 * p_mb * qst / (constants.EPSILON + (1.0 - constants.EPSILON) * qst)
        )  # (100's <-^ convert from mbar to Pa)
        rhx = min(q / qst, 1.00)
        k1 = (
            (constants.MAPL_LATENT_HEAT_VAPORIZATION**2)
            * constants.RHO_W
            / (constants.K_COND * constants.MAPL_RVAP * (t**2))
        )
        k2 = (
            constants.MAPL_RVAP
            * t
            * constants.RHO_W
            / (constants.DIFFU * (1000.0 / p_mb) * es)
        )
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if clcn > 0.0 and qlcn > 0.0:
            qcm = qlcn / clcn
        else:
            qcm = 0.0
        radius = cloud_effective_radius_liquid(p_mb, t, qcm, nactl, nacti)
        if rhx < rh_crit and radius > 0.0:
            evap = (
                CCW_EVAP_EFF
                * qlcn
                * DT_MOIST
                * (rh_crit - rhx)
                / ((k1 + k2) * radius**2)
            )
            evap = min(evap, qlcn)
        else:
            evap = 0.0
        qc = qlcn + qicn
        if qc > 0.0:
            clcn = clcn * (qc - evap) / qc
        q = q + evap
        qlcn = qlcn - evap
        t = t - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CPDRY) * evap
        evapc = (q - evapc) / DT_MOIST
