from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, interval, sqrt
from ndsl.dsl.typing import FloatField, FloatFieldIJ

from pyMoist.microphysics.GFDL_1M.driver.constants import constants


def update_tendencies(
    mixing_ratio_vapor_unmodified: FloatField,
    mixing_ratio_liquid_unmodified: FloatField,
    mixing_ratio_rain_unmodified: FloatField,
    mixing_ratio_ice_unmodified: FloatField,
    mixing_ratio_snow_unmodified: FloatField,
    mixing_ratio_graupel_unmodified: FloatField,
    cloud_fraction_unmodified: FloatField,
    mixing_ratio_driver_vapor: FloatField,
    mixing_ratio_driver_liquid: FloatField,
    mixing_ratio_driver_rain: FloatField,
    mixing_ratio_driver_ice: FloatField,
    mixing_ratio_driver_snow: FloatField,
    mixing_ratio_driver_graupel: FloatField,
    dvapordt: FloatField,
    dliquiddt: FloatField,
    draindt: FloatField,
    dicedt: FloatField,
    dsnowdt: FloatField,
    dgraupeldt: FloatField,
    dcloudfractiondt: FloatField,
    t_unmodified: FloatField,
    driver_t: FloatField,
    dtdt: FloatField,
    w_unmodified: FloatField,
    driver_w: FloatField,
    u_unmodified: FloatField,
    driver_u: FloatField,
    dudt: FloatField,
    v_unmodified: FloatField,
    driver_v: FloatField,
    dvdt: FloatField,
    dp_unmodified: FloatField,
    driver_dp: FloatField,
    driver_mass: FloatField,
    rain: FloatFieldIJ,
    snow: FloatFieldIJ,
    ice: FloatFieldIJ,
    graupel: FloatFieldIJ,
):
    """
    Compute output tendencies of the microphysics driver

    Arguments:
        mixing_ratio_vapor_unmodified (in): water vapor mixing ratio, unmodified within driver (kg/kg)
        mixing_ratio_liquid_unmodified (in): in cloud liquid water, unmodified within driver (kg/kg)
        mixing_ratio_rain_unmodified (in): falling rain, unmodified within driver (kg/kg)
        mixing_ratio_ice_unmodified (in): in cloud frozen water, unmodified within driver (kg/kg)
        mixing_ratio_snow_unmodified (in): in cloud snow, unmodified within driver (kg/kg)
        mixing_ratio_graupel_unmodified (in): in cloud graupel, unmodified within driver (kg/kg)
        cloud_fraction_unmodified (in): cloud fraction (convective + large scale), unmodified within driver
        mixing_ratio_driver_vapor (in): water vapor mixing ratio, driver modified
        mixing_ratio_driver_liquid (in): in cloud liquid water, driver modified
        mixing_ratio_driver_rain (in): falling rain mixing ratio, driver modified
        mixing_ratio_driver_ice (in): in cloud frozen water, driver modified
        mixing_ratio_driver_snow (in): in cloud snow, driver modified
        mixing_ratio_driver_graupel (in): in cloud graupel, driver modified
        dvapordt (out): water vapor tendency
        dliquiddt (out): in cloud liquid water tendency
        draindt (out): falling rain tendency
        dicedt (out): in cloud frozen water tendency
        dsnowdt (out): in cloud snow tendency
        dgraupeldt (out): in cloud graupel tendency
        dcloudfractiondt (out): cloud fraction (convective + large scale) tendency
        t_unmodified (in): atmospheric temperature, unmodified within driver (K)
        driver_t (in): atmospheric temperature, driver modified (K)
        dtdt (out): atmospheric temperature tendency
        w_unmodified (in): vertical velocity, unmodified within driver (m/s)
        driver_w (in): vertical velocity, driver modified (m/s)
        u_unmodified (in): eastward winds, unmodified within driver (m/s)
        driver_u: eastward winds, driver modified (m/s)
        dudt: eastward wind tendency
        v_unmodified: northward winds, unmodified within driver (m/s)
        driver_v: northward winds, driver modified (m/s)
        dvdt: northward wind tendency
        dp_unmodified: change in pressure between model levels, unmodified by driver (mb)
        driver_dp: change in pressure between model levels, driver modified (mb)
        driver_mass: details unknown
        rain: precipitated rain at surface (kg/m^2/s)
        snow: precipitated snow at surface (kg/m^2/s)
        ice: precipitated ice at surface (kg/m^2/s)
        graupel: precipitated graupel at surface (kg/m^2/s)

    reference Fortran: gfdl_cloud_microphys.F90:
    subroutines mpdrv, gfdl_cloud_microphys_driver
    """
    from __externals__ import c_air, c_vap, do_qa, do_sedi_w, rdt, sedi_transport

    # -----------------------------------------------------------------------
    # momentum transportation during sedimentation
    # note: dp1 is dry mass; dp0 is the old moist (total) mass
    # -----------------------------------------------------------------------

    with computation(FORWARD), interval(1, None):
        if sedi_transport == True:  # noqa
            driver_u = (dp_unmodified * driver_u + driver_mass[0, 0, -1] * driver_u[0, 0, -1]) / (dp_unmodified + driver_mass[0, 0, -1])
            driver_v = (dp_unmodified * driver_v + driver_mass[0, 0, -1] * driver_v[0, 0, -1]) / (dp_unmodified + driver_mass[0, 0, -1])
            dudt = dudt + (driver_u - u_unmodified) * rdt
            dvdt = dvdt + (driver_v - v_unmodified) * rdt

    with computation(PARALLEL), interval(...):
        if do_sedi_w:
            w_unmodified = driver_w

    # -----------------------------------------------------------------------
    # update moist air mass (actually hydrostatic pressure)
    # convert to dry mixing ratios
    # -----------------------------------------------------------------------

    with computation(PARALLEL), interval(...):
        omq = driver_dp / dp_unmodified
        dvapordt = dvapordt + rdt * (mixing_ratio_driver_vapor - mixing_ratio_vapor_unmodified) * omq
        dliquiddt = dliquiddt + rdt * (mixing_ratio_driver_liquid - mixing_ratio_liquid_unmodified) * omq
        draindt = draindt + rdt * (mixing_ratio_driver_rain - mixing_ratio_rain_unmodified) * omq
        dicedt = dicedt + rdt * (mixing_ratio_driver_ice - mixing_ratio_ice_unmodified) * omq
        dsnowdt = dsnowdt + rdt * (mixing_ratio_driver_snow - mixing_ratio_snow_unmodified) * omq
        dgraupeldt = dgraupeldt + rdt * (mixing_ratio_driver_graupel - mixing_ratio_graupel_unmodified) * omq
        cvm = (
            c_air
            + mixing_ratio_driver_vapor * c_vap
            + (mixing_ratio_driver_rain + mixing_ratio_driver_liquid) * constants.C_LIQ
            + (mixing_ratio_driver_ice + mixing_ratio_driver_snow + mixing_ratio_driver_graupel) * constants.C_ICE
        )
        dtdt = dtdt + rdt * (driver_t - t_unmodified) * cvm / constants.CP_AIR

        # -----------------------------------------------------------------------
        # update cloud fraction tendency
        # -----------------------------------------------------------------------
        if do_qa == False:  # noqa
            dcloudfractiondt = dcloudfractiondt + rdt * (
                cloud_fraction_unmodified
                * sqrt((mixing_ratio_driver_ice + mixing_ratio_driver_liquid) / max(mixing_ratio_ice_unmodified + mixing_ratio_liquid_unmodified, constants.QCMIN))
                - cloud_fraction_unmodified
            )  # New Cloud - Old CloudCloud

    with computation(FORWARD), interval(0, 1):
        # convert to mm / day
        conversion_factor = 86400.0 * rdt * constants.RGRAV
        rain = rain * conversion_factor
        snow = snow * conversion_factor
        ice = ice * conversion_factor
        graupel = graupel * conversion_factor
