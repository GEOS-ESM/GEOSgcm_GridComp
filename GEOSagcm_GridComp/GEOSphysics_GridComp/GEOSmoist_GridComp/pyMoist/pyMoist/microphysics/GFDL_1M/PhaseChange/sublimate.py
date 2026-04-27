from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField

import pyMoist.constants as constants
from pyMoist.shared.incloud_processes import cloud_effective_radius_ice


def sublimate(
    p_mb: FloatField,
    t: FloatField,
    mixing_ratio_vapor: FloatField,
    mixing_ratio_convective_liquid: FloatField,
    mixing_ratio_convective_ice: FloatField,
    convective_cloud_fraction: FloatField,
    liquid_concentration: FloatField,
    ice_concentration: FloatField,
    saturation_specific_humidity: FloatField,
    sublimation: FloatField,
):
    """Quantify evaporation of excess ice prior to microphysics driver run

    Args:
        p_mb (FloatField): pressure (mb)
        t (FloatField): temperature (Kelvin)
        mixing_ratio_vapor (FloatField): water vapor mixing ratio (kg/kg)
        mixing_ratio_convective_liquid (FloatField): unitless
        mixing_ratio_convective_ice (FloatField): unitless
        convective_cloud_fraction (FloatField): unitless
        liquid_concentration (FloatField): liquid particle concentration (m^-3)
        ice_concentration (FloatField): ice particle concentration (m^-3)
        saturation_specific_humidity (FloatField)
        sublimation (FloatField): sublimation of cloud ice (kg kg-1 s-1)
    """
    from __externals__ import CCI_EVAP_EFF, DT_MOIST

    with computation(PARALLEL), interval(...):
        sublimation = mixing_ratio_vapor
        rh_crit = 1
        # Sublimation of cloud water. DelGenio et al formulation
        # (Eq.s 15-17, 1996, J. Clim., 9, 270-303)
        es = (
            100.0 * p_mb * saturation_specific_humidity / (constants.EPSILON + (1.0 - constants.EPSILON) * saturation_specific_humidity)
        )  # (100s <-^ convert from mbar to Pa)
        rhx = min(mixing_ratio_vapor / saturation_specific_humidity, 1.00)
        k1 = (constants.MAPL_LATENT_HEAT_VAPORIZATION**2) * constants.RHO_I / (constants.K_COND * constants.MAPL_RVAP * (t**2))
        k2 = constants.MAPL_RVAP * t * constants.RHO_I / (constants.DIFFU * (1000.0 / p_mb) * es)
        # Here, DIFFU is given for 1000 mb so 1000./PLmb accounts
        # for increased diffusivity at lower pressure
        if convective_cloud_fraction > 0.0 and mixing_ratio_convective_ice > 0.0:
            qcm = mixing_ratio_convective_ice / convective_cloud_fraction
        else:
            qcm = 0.0
        radius = cloud_effective_radius_ice(p_mb, t, qcm)
        if rhx < rh_crit and radius > 0.0:
            subl = CCI_EVAP_EFF * mixing_ratio_convective_ice * DT_MOIST * (rh_crit - rhx) / ((k1 + k2) * radius**2)
            subl = min(subl, mixing_ratio_convective_ice)
        else:
            subl = 0.0
        qc = mixing_ratio_convective_liquid + mixing_ratio_convective_ice
        if qc > 0.0:
            convective_cloud_fraction = convective_cloud_fraction * (qc - subl) / qc
        mixing_ratio_vapor = mixing_ratio_vapor + subl
        mixing_ratio_convective_ice = mixing_ratio_convective_ice - subl
        t = t - (constants.MAPL_LATENT_HEAT_SUBLIMATION / constants.MAPL_CPDRY) * subl
        sublimation = (mixing_ratio_vapor - sublimation) / DT_MOIST
