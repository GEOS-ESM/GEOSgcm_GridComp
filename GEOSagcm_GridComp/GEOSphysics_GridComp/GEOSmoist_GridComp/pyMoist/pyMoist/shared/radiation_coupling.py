from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField

from pyMoist.shared.incloud_processes import cloud_effective_radius_ice, cloud_effective_radius_liquid


def radiation_coupling(
    temperature: FloatField,
    pressure: FloatField,
    large_scale_cloud_fraction: FloatField,
    convective_cloud_fraction: FloatField,
    vapor: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    liquid_concentration: FloatField,
    ice_concentration: FloatField,
    radiation_vapor: FloatField,
    radiation_liquid: FloatField,
    radiation_ice: FloatField,
    radiation_rain: FloatField,
    radiation_snow: FloatField,
    radiation_graupel: FloatField,
    radiation_cloud_fraction: FloatField,
    liquid_radius: FloatField,
    ice_radius: FloatField,
) -> None:
    """Couple radiation with cloud variables to ensure physical consistency.

    Args:
        temperature (FloatField)
        pressure (FloatField)
        large_scale_cloud_fraction (FloatField)
        convective_cloud_fraction (FloatField)
        vapor (FloatField)
        large_scale_liquid (FloatField)
        large_scale_ice (FloatField)
        convective_liquid (FloatField)
        convective_ice (FloatField)
        rain (FloatField)
        snow (FloatField)
        graupel (FloatField)
        liquid_concentration (FloatField)
        ice_concentration (FloatField)
        radiation_vapor (FloatField)
        radiation_liquid (FloatField)
        radiation_ice (FloatField)
        radiation_rain (FloatField)
        radiation_snow (FloatField)
        radiation_graupel (FloatField)
        radiation_cloud_fraction (FloatField)
        liquid_radius (FloatField)
        ice_radius (FloatField)
    """

    from __externals__ import FAC_RI, FAC_RL, MAX_RI, MAX_RL, MIN_RI, MIN_RL

    with computation(PARALLEL), interval(...):
        # water vapor
        radiation_vapor = vapor

        # total cloud fraction
        radiation_cloud_fraction = max(min(large_scale_cloud_fraction + convective_cloud_fraction, 1.0), 0.0)
        if radiation_cloud_fraction >= 1.0e-5:
            radiation_liquid = (large_scale_liquid + convective_liquid) / radiation_cloud_fraction if (large_scale_liquid + convective_liquid) >= 1.0e-8 else 0.0
            radiation_ice = (large_scale_ice + convective_ice) / radiation_cloud_fraction if (large_scale_ice + convective_ice) >= 1.0e-8 else 0.0
            radiation_rain = rain / radiation_cloud_fraction if rain >= 1.0e-8 else 0.0
            radiation_snow = snow / radiation_cloud_fraction if snow >= 1.0e-8 else 0.0
            radiation_graupel = graupel / radiation_cloud_fraction if graupel >= 1.0e-8 else 0.0
        else:
            radiation_cloud_fraction = 0.0
            radiation_liquid = 0.0
            radiation_ice = 0.0
            radiation_rain = 0.0
            radiation_snow = 0.0
            radiation_graupel = 0.0

        # Cap the high end of condensates
        radiation_liquid = min(radiation_liquid, 0.01)
        radiation_ice = min(radiation_ice, 0.01)
        radiation_rain = min(radiation_rain, 0.01)
        radiation_snow = min(radiation_snow, 0.01)
        radiation_graupel = min(radiation_graupel, 0.01)

        # Liquid radii - Brams formulation with limits
        liquid_radius = max(
            MIN_RL,
            min(
                cloud_effective_radius_liquid(pressure, temperature, radiation_liquid, liquid_concentration) * FAC_RL,
                MAX_RL,
            ),
        )
        # Ice radii - Brams formulation with limits
        ice_radius = max(
            MIN_RI,
            min(
                cloud_effective_radius_ice(pressure, temperature, radiation_ice) * FAC_RI,
                MAX_RI,
            ),
        )
