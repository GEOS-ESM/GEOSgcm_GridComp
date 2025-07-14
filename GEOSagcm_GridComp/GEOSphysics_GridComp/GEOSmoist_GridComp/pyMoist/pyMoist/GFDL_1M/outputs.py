from dataclasses import dataclass

import numpy.typing as npt

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM


@dataclass(frozen=True)
class Outputs:
    """
    Collection of all fields computed locally and returned to the rest of the model
    """

    lower_tropospheric_stability: Quantity | npt.NDArray
    estimated_inversion_strength: Quantity | npt.NDArray
    du_dt_macro: Quantity | npt.NDArray
    dv_dt_macro: Quantity | npt.NDArray
    dt_dt_macro: Quantity | npt.NDArray
    dvapor_dt_macro: Quantity | npt.NDArray
    dliquid_dt_macro: Quantity | npt.NDArray
    dice_dt_macro: Quantity | npt.NDArray
    dcloud_fraction_dt_macro: Quantity | npt.NDArray
    drain_dt_macro: Quantity | npt.NDArray
    dsnow_dt_macro: Quantity | npt.NDArray
    dgraupel_dt_macro: Quantity | npt.NDArray
    du_dt_micro: Quantity | npt.NDArray
    dv_dt_micro: Quantity | npt.NDArray
    dt_dt_micro: Quantity | npt.NDArray
    dvapor_dt_micro: Quantity | npt.NDArray
    dliquid_dt_micro: Quantity | npt.NDArray
    dice_dt_micro: Quantity | npt.NDArray
    dcloud_fraction_dt_micro: Quantity | npt.NDArray
    drain_dt_micro: Quantity | npt.NDArray
    dsnow_dt_micro: Quantity | npt.NDArray
    dgraupel_dt_micro: Quantity | npt.NDArray
    radiation_cloud_fraction: Quantity | npt.NDArray
    radiation_ice: Quantity | npt.NDArray
    radiation_liquid: Quantity | npt.NDArray
    radiation_vapor: Quantity | npt.NDArray
    radiation_rain: Quantity | npt.NDArray
    radiation_snow: Quantity | npt.NDArray
    radiation_graupel: Quantity | npt.NDArray
    ice_radius: Quantity | npt.NDArray
    liquid_radius: Quantity | npt.NDArray
    large_scale_nonanvil_precipitation_evaporation: Quantity | npt.NDArray
    large_scale_nonanvil_precipitation_sublimation: Quantity | npt.NDArray
    relative_humidity_after_pdf: Quantity | npt.NDArray
    precipitated_rain: Quantity | npt.NDArray
    precipitated_snow: Quantity | npt.NDArray
    precipitated_ice: Quantity | npt.NDArray
    precipitated_graupel: Quantity | npt.NDArray
    large_scale_precip: Quantity | npt.NDArray
    large_scale_snow: Quantity | npt.NDArray
    icefall: Quantity | npt.NDArray
    freezing_rainfall: Quantity | npt.NDArray
    large_scale_nonanvil_ice_flux: Quantity | npt.NDArray
    large_scale_nonanvil_liquid_flux: Quantity | npt.NDArray
    anvil_liquid_flux: Quantity | npt.NDArray
    anvil_ice_flux: Quantity | npt.NDArray
    deep_convective_precipitation: Quantity | npt.NDArray
    anvil_precipitation: Quantity | npt.NDArray
    shallow_convective_precipitation: Quantity | npt.NDArray
    deep_convective_snow: Quantity | npt.NDArray
    anvil_snow: Quantity | npt.NDArray
    shallow_convective_snow: Quantity | npt.NDArray
    # Optional outputs
    z_lcl: Quantity | npt.NDArray | None
    large_scale_rainwater_source: Quantity | npt.NDArray | None
    moist_friction_temperature_tendency: Quantity | npt.NDArray | None
    simulated_reflectivity: Quantity | npt.NDArray | None
    maximum_reflectivity: Quantity | npt.NDArray | None
    one_km_agl_reflectivity: Quantity | npt.NDArray | None
    echo_top_reflectivity: Quantity | npt.NDArray | None
    minus_10c_reflectivity: Quantity | npt.NDArray | None

    @classmethod
    def zeros(cls, quantity_factory: QuantityFactory):
        attrs = {}
        _3D_fields = [
            "du_dt_macro",
            "dv_dt_macro",
            "dt_dt_macro",
            "dvapor_dt_macro",
            "dliquid_dt_macro",
            "dice_dt_macro",
            "dcloud_fraction_dt_macro",
            "drain_dt_macro",
            "dsnow_dt_macro",
            "dgraupel_dt_macro",
            "du_dt_micro",
            "dv_dt_micro",
            "dt_dt_micro",
            "dvapor_dt_micro",
            "dliquid_dt_micro",
            "dice_dt_micro",
            "dcloud_fraction_dt_micro",
            "drain_dt_micro",
            "dsnow_dt_micro",
            "dgraupel_dt_micro",
            "radiation_cloud_fraction",
            "radiation_ice",
            "radiation_liquid",
            "radiation_vapor",
            "radiation_rain",
            "radiation_snow",
            "radiation_graupel",
            "ice_radius",
            "liquid_radius",
            "large_scale_nonanvil_precipitation_evaporation",
            "large_scale_nonanvil_precipitation_sublimation",
            "relative_humidity_after_pdf",
            "large_scale_rainwater_source",
            "moist_friction_temperature_tendency",
            "simulated_reflectivity",
        ]
        for field in _3D_fields:
            attrs[field] = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        _2D_fields = [
            "lower_tropospheric_stability",
            "estimated_inversion_strength",
            "z_lcl",
            "precipitated_rain",
            "precipitated_snow",
            "precipitated_ice",
            "precipitated_graupel",
            "large_scale_precip",
            "large_scale_snow",
            "icefall",
            "freezing_rainfall",
            "maximum_reflectivity",
            "one_km_agl_reflectivity",
            "echo_top_reflectivity",
            "minus_10c_reflectivity",
            "deep_convective_precipitation",
            "anvil_precipitation",
            "shallow_convective_precipitation",
            "deep_convective_snow",
            "anvil_snow",
            "shallow_convective_snow",
        ]
        for field in _2D_fields:
            attrs[field] = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        _3D_fields_K_interface = [
            "large_scale_nonanvil_ice_flux",
            "large_scale_nonanvil_liquid_flux",
            "anvil_liquid_flux",
            "anvil_ice_flux",
        ]
        for field in _3D_fields_K_interface:
            attrs[field] = quantity_factory.zeros(
                [X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a"
            )

        return cls(**attrs)
