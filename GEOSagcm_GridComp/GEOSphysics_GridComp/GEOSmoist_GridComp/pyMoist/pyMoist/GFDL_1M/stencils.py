from ndsl.dsl.gt4py import PARALLEL, computation, interval
from ndsl.dsl.typing import FloatField


def prepare_tendencies(
    u: FloatField,
    v: FloatField,
    t: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    du_dt: FloatField,
    dv_dt: FloatField,
    dt_dt: FloatField,
    dvapor_dt: FloatField,
    dliquid_dt: FloatField,
    dice_dt: FloatField,
    dcloud_fraction_dt: FloatField,
    drain_dt: FloatField,
    dsnow_dt: FloatField,
    dgraupel_dt: FloatField,
):
    with computation(PARALLEL), interval(...):
        du_dt = u
        dv_dt = v
        dt_dt = t
        dvapor_dt = vapor
        dliquid_dt = convective_liquid + large_scale_liquid
        dice_dt = convective_ice + large_scale_ice
        dcloud_fraction_dt = convective_cloud_fraction + large_scale_cloud_fraction
        drain_dt = rain
        dsnow_dt = snow
        dgraupel_dt = graupel


def update_radiation(
    t: FloatField,
    u: FloatField,
    v: FloatField,
    radiation_cloud_fraction: FloatField,
    radiation_ice: FloatField,
    radiation_liquid: FloatField,
    radiation_vapor: FloatField,
    radiation_rain: FloatField,
    radiation_snow: FloatField,
    radiation_graupel: FloatField,
    dcloud_fraction_dt: FloatField,
    dtdt: FloatField,
    dudt: FloatField,
    dvdt: FloatField,
    dicedt: FloatField,
    dliquiddt: FloatField,
    dvapordt: FloatField,
    draindt: FloatField,
    dsnowdt: FloatField,
    dgraupeldt: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        t = t + dtdt * DT_MOIST
        u = u + dudt * DT_MOIST
        v = v + dvdt * DT_MOIST
        radiation_cloud_fraction = min(
            1.0, max(0.0, radiation_cloud_fraction + dcloud_fraction_dt * DT_MOIST)
        )
        radiation_ice = radiation_ice + dicedt * DT_MOIST
        radiation_liquid = radiation_liquid + dliquiddt * DT_MOIST
        radiation_vapor = radiation_vapor + dvapordt * DT_MOIST
        radiation_rain = radiation_rain + draindt * DT_MOIST
        radiation_snow = radiation_snow + dsnowdt * DT_MOIST
        radiation_graupel = radiation_graupel + dgraupeldt * DT_MOIST


def update_tendencies(
    u: FloatField,
    v: FloatField,
    t: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    du_dt: FloatField,
    dv_dt: FloatField,
    dt_dt: FloatField,
    dvapor_dt: FloatField,
    dliquid_dt: FloatField,
    dice_dt: FloatField,
    dcloud_fraction_dt: FloatField,
    drain_dt: FloatField,
    dsnow_dt: FloatField,
    dgraupel_dt: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        du_dt = (u - du_dt) / DT_MOIST
        dv_dt = (v - dv_dt) / DT_MOIST
        dt_dt = (t - dt_dt) / DT_MOIST
        dvapor_dt = (vapor - dvapor_dt) / DT_MOIST
        dliquid_dt = ((convective_liquid + large_scale_liquid) - dliquid_dt) / DT_MOIST
        dice_dt = ((convective_ice + large_scale_ice) - dice_dt) / DT_MOIST
        dcloud_fraction_dt = (
            (convective_cloud_fraction + large_scale_cloud_fraction) - dcloud_fraction_dt
        ) / DT_MOIST
        drain_dt = (rain - drain_dt) / DT_MOIST
        dsnow_dt = (snow - dsnow_dt) / DT_MOIST
        dgraupel_dt = (graupel - dgraupel_dt) / DT_MOIST


def get_total_concentration(
    ice_concentration: FloatField,
    liquid_concentration: FloatField,
    total_concentration: FloatField,
):
    with computation(PARALLEL), interval(...):
        total_concentration = ice_concentration + liquid_concentration


def prepare_radiation(
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    radiation_cloud_fraction: FloatField,
    convective_liquid: FloatField,
    large_scale_liquid: FloatField,
    radiation_liquid: FloatField,
    convective_ice: FloatField,
    large_scale_ice: FloatField,
    radiation_ice: FloatField,
    vapor: FloatField,
    radiation_vapor: FloatField,
    rain: FloatField,
    radiation_rain: FloatField,
    snow: FloatField,
    radiation_snow: FloatField,
    graupel: FloatField,
    radiation_graupel: FloatField,
):
    with computation(PARALLEL), interval(...):
        # cloud fraction
        radiation_cloud_fraction = min(convective_cloud_fraction + large_scale_cloud_fraction, 1.0)
        # liquid
        radiation_liquid = convective_liquid + large_scale_liquid
        # ice
        radiation_ice = convective_ice + large_scale_ice
        # vapor
        radiation_vapor = vapor
        # RAIN
        radiation_rain = rain
        # snow
        radiation_snow = snow
        # graupel
        radiation_graupel = graupel


def reset_micro_tendencies(
    dvapordt: FloatField,
    dliquiddt: FloatField,
    draindt: FloatField,
    dicedt: FloatField,
    dsnowdt: FloatField,
    dgraupeldt: FloatField,
    dcloudfractiondt: FloatField,
    dtdt: FloatField,
    dudt: FloatField,
    dvdt: FloatField,
):
    with computation(PARALLEL), interval(...):
        dvapordt = 0
        dliquiddt = 0
        draindt = 0
        dicedt = 0
        dsnowdt = 0
        dgraupeldt = 0
        dcloudfractiondt = 0
        dtdt = 0
        dudt = 0
        dvdt = 0
