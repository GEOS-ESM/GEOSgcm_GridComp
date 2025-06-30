from gt4py.cartesian.gtscript import PARALLEL, computation, interval

import pyMoist.constants as constants
from ndsl.dsl.typing import Float, FloatField
from pyMoist.shared_incloud_processes import cloud_effective_radius_liquid


def evaporate(
    p_mb: FloatField,
    t: FloatField,
    vapor: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    convective_cloud_fraction: FloatField,
    nactl: FloatField,
    nacti: FloatField,
    qsat: FloatField,
    evapc: FloatField,
):
    from __externals__ import CCW_EVAP_EFF, DT_MOIST

    with computation(PARALLEL), interval(...):
        evapc = vapor
        rh_crit = 1
        # Evaporation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        es = (
            100.0 * p_mb * qsat / (constants.EPSILON + (1.0 - constants.EPSILON) * qsat)
        )  # (100's <-^ convert from mbar to Pa)
        rhx = min(vapor / qsat, 1.00)
        k1 = (
            (constants.MAPL_LATENT_HEAT_VAPORIZATION ** 2)
            * constants.RHO_W
            / (constants.K_COND * constants.MAPL_RVAP * (t ** 2))
        )
        k2 = constants.MAPL_RVAP * t * constants.RHO_W / (constants.DIFFU * (1000.0 / p_mb) * es)
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if convective_cloud_fraction > 0.0 and convective_liquid > 0.0:
            qcm = convective_liquid / convective_cloud_fraction
        else:
            qcm = 0.0
        radius = cloud_effective_radius_liquid(p_mb, t, qcm, nactl)
        if rhx < rh_crit and radius > 0.0:
            evap = CCW_EVAP_EFF * convective_liquid * DT_MOIST * (rh_crit - rhx) / ((k1 + k2) * radius ** 2)
            evap = min(evap, convective_liquid)
        else:
            evap = 0.0
        qc = convective_liquid + convective_ice
        if qc > 0.0:
            convective_cloud_fraction = convective_cloud_fraction * (qc - evap) / qc
        vapor = vapor + evap
        convective_liquid = convective_liquid - evap
        t = t - (constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CPDRY) * evap
        evapc = (vapor - evapc) / DT_MOIST
