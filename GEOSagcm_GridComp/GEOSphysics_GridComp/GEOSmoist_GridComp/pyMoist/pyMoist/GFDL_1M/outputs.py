from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM


@dataclass
class Outputs:
    lower_tropospheric_stability: Quantity
    estimated_inversion_stregnth: Quantity
    z_lcl: Quantity
    du_dt_macro: Quantity
    dv_dt_macro: Quantity
    dt_dt_macro: Quantity
    dvapor_dt_macro: Quantity
    dliquid_dt_macro: Quantity
    dice_dt_macro: Quantity
    dcloud_fraction_dt_macro: Quantity
    drain_dt_macro: Quantity
    dsnow_dt_macro: Quantity
    dgraupel_dt_macro: Quantity
    du_dt_micro: Quantity
    dv_dt_micro: Quantity
    dt_dt_micro: Quantity
    dvapor_dt_micro: Quantity
    dliquid_dt_micro: Quantity
    dice_dt_micro: Quantity
    dcloud_fraction_dt_micro: Quantity
    drain_dt_micro: Quantity
    dsnow_dt_micro: Quantity
    dgraupel_dt_micro: Quantity
    radiation_cloud_fraction: Quantity
    radiation_ice: Quantity
    radiation_liquid: Quantity
    radiation_vapor: Quantity
    radiation_rain: Quantity
    radiation_snow: Quantity
    radiation_graupel: Quantity
    ice_radius: Quantity
    liquid_radius: Quantity
    large_scale_nonanvil_precipitation_evaporation: Quantity
    large_scale_nonanvil_precipitation_sublimation: Quantity
    relative_humidity_after_pdf: Quantity
    precipitated_rain: Quantity
    precipitated_snow: Quantity
    precipitated_ice: Quantity
    precipitated_graupel: Quantity
    large_scale_precip: Quantity
    large_scale_snow: Quantity
    icefall: Quantity
    freezing_rainfall: Quantity
    large_scale_nonanvil_ice_flux: Quantity
    large_scale_nonanvil_liquid_flux: Quantity
    anvil_liquid_flux: Quantity
    anvil_ice_flux: Quantity
    large_scale_rainwater_source: Quantity
    moist_friction_temperature_tendency: Quantity
    simulated_reflectivity: Quantity
    maximum_reflectivity: Quantity
    one_km_agl_reflectivity: Quantity
    echo_top_reflectivity: Quantity
    minus_10c_reflectivity: Quantity
    deep_convective_precipitation: Quantity
    anvil_precipitation: Quantity
    shallow_convective_precipitation: Quantity
    deep_convective_snow: Quantity
    anvil_snow: Quantity
    shallow_convective_snow: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):
        lower_tropospheric_stability = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        estimated_inversion_stregnth = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        z_lcl = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        du_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dv_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dt_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapor_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dliquid_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dice_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcloud_fraction_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        drain_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dsnow_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dgraupel_dt_macro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        du_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dv_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dt_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapor_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dliquid_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dice_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcloud_fraction_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        drain_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dsnow_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dgraupel_dt_micro = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_cloud_fraction = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_ice = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_liquid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_vapor = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_rain = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_snow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        radiation_graupel = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ice_radius = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        liquid_radius = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        large_scale_nonanvil_precipitation_evaporation = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        large_scale_nonanvil_precipitation_sublimation = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        relative_humidity_after_pdf = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        precipitated_rain = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        precipitated_snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        precipitated_ice = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        precipitated_graupel = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        large_scale_precip = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        large_scale_snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        icefall = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        freezing_rainfall = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        large_scale_nonanvil_ice_flux = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        large_scale_nonanvil_liquid_flux = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        anvil_liquid_flux = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        anvil_ice_flux = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        large_scale_rainwater_source = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        moist_friction_temperature_tendency = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        simulated_reflectivity = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        maximum_reflectivity = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        one_km_agl_reflectivity = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        echo_top_reflectivity = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        minus_10c_reflectivity = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        deep_convective_precipitation = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        anvil_precipitation = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        shallow_convective_precipitation = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        deep_convective_snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        anvil_snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        shallow_convective_snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        return cls(
            lower_tropospheric_stability,
            estimated_inversion_stregnth,
            z_lcl,
            du_dt_macro,
            dv_dt_macro,
            dt_dt_macro,
            dvapor_dt_macro,
            dliquid_dt_macro,
            dice_dt_macro,
            dcloud_fraction_dt_macro,
            drain_dt_macro,
            dsnow_dt_macro,
            dgraupel_dt_macro,
            du_dt_micro,
            dv_dt_micro,
            dt_dt_micro,
            dvapor_dt_micro,
            dliquid_dt_micro,
            dice_dt_micro,
            dcloud_fraction_dt_micro,
            drain_dt_micro,
            dsnow_dt_micro,
            dgraupel_dt_micro,
            radiation_cloud_fraction,
            radiation_ice,
            radiation_liquid,
            radiation_vapor,
            radiation_rain,
            radiation_snow,
            radiation_graupel,
            ice_radius,
            liquid_radius,
            large_scale_nonanvil_precipitation_evaporation,
            large_scale_nonanvil_precipitation_sublimation,
            relative_humidity_after_pdf,
            precipitated_rain,
            precipitated_snow,
            precipitated_ice,
            precipitated_graupel,
            large_scale_precip,
            large_scale_snow,
            icefall,
            freezing_rainfall,
            large_scale_nonanvil_ice_flux,
            large_scale_nonanvil_liquid_flux,
            anvil_liquid_flux,
            anvil_ice_flux,
            large_scale_rainwater_source,
            moist_friction_temperature_tendency,
            simulated_reflectivity,
            maximum_reflectivity,
            one_km_agl_reflectivity,
            echo_top_reflectivity,
            minus_10c_reflectivity,
            deep_convective_precipitation,
            anvil_precipitation,
            shallow_convective_precipitation,
            deep_convective_snow,
            anvil_snow,
            shallow_convective_snow,
        )
