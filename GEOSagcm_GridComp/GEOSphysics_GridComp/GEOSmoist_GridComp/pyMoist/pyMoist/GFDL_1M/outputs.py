from dataclasses import dataclass

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM


@dataclass
class Outputs:
    lower_tropospheric_stability: Quantity
    estimated_inversion_strength: Quantity
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
    evaporation: Quantity
    sublimation: Quantity
    precipitated_rain: Quantity
    precipitated_snow: Quantity
    precipitated_ice: Quantity
    precipitated_graupel: Quantity

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
        evaporation = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        sublimation = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        precipitated_rain = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        precipitated_snow = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        precipitated_ice = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        precipitated_graupel = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

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
            evaporation,
            sublimation,
            precipitated_rain,
            precipitated_snow,
            precipitated_ice,
            precipitated_graupel,
        )
