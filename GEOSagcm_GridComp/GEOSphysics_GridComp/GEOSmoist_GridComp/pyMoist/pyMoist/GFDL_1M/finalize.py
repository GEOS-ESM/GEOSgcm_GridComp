from ndsl import StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, function, interval, sqrt
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.constants import MAPL_CP, MAPL_GRAV
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.GFDL_1M.outputs import Outputs
from pyMoist.GFDL_1M.state import CloudFractions, MixingRatios
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.radiation_coupling import GFDL1MRadiationCoupling
from pyMoist.redistribute_clouds import RedistributeClouds
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


@function
def fix_negative_precip(
    precip: Float,
):
    if precip < 1.0e-8:
        precip = 0.0

    return precip


def finalize_precip(
    rain_from_driver: FloatFieldIJ,
    snow_from_driver: FloatFieldIJ,
    ice_from_driver: FloatFieldIJ,
    graupel_from_driver: FloatFieldIJ,
    precipitated_rain: FloatFieldIJ,
    precipitated_snow: FloatFieldIJ,
    precipitated_ice: FloatFieldIJ,
    precipitated_graupel: FloatFieldIJ,
    large_scale_precip: FloatFieldIJ,
    large_scale_snow: FloatFieldIJ,
    icefall: FloatFieldIJ,
    freezing_rainfall: FloatFieldIJ,
    large_scale_nonanvil_ice_flux: FloatField,
    large_scale_nonanvil_ice_flux_from_driver: FloatField,
    large_scale_nonanvil_liquid_flux: FloatField,
    large_scale_nonanvil_liquid_flux_from_driver: FloatField,
    convective_liquid: FloatField,
    liquid_for_radiation: FloatField,
    anvil_liquid_flux: FloatField,
    convective_ice: FloatField,
    ice_for_radiation: FloatField,
    anvil_ice_flux: FloatField,
    evaporation: FloatField,
    evaporation_from_driver: FloatField,
    sublimation: FloatField,
    sublimation_from_driver: FloatField,
    radiation_vapor: FloatField,
    radiation_rain: FloatField,
    radiation_snow: FloatField,
    radiation_graupel: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
):
    """
    Must be constructed using Z_INTERFACE_DIM
    """
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(0, -1):
        # send driver evaporation and sublimation outputs back to the rest of the model
        evaporation = evaporation_from_driver
        sublimation = sublimation_from_driver

    with computation(FORWARD), interval(0, 1):
        # send precip diagnostics back to the rest of the model
        # and convert from mm/day to kg m-2 s-1
        precipitated_rain = max(rain_from_driver / 86400.0, 0.0)
        precipitated_snow = max(snow_from_driver / 86400.0, 0.0)
        precipitated_ice = max(ice_from_driver / 86400.0, 0.0)
        precipitated_graupel = max(graupel_from_driver / 86400.0, 0.0)
        # Fill GEOS precip diagnostics
        large_scale_precip = precipitated_rain
        large_scale_snow = precipitated_snow
        icefall = precipitated_ice + precipitated_graupel
        freezing_rainfall = 0.0

    with computation(PARALLEL), interval(1, None):
        # Bring in precipitation fluxes from driver
        large_scale_nonanvil_ice_flux = large_scale_nonanvil_ice_flux_from_driver
        large_scale_nonanvil_liquid_flux = large_scale_nonanvil_liquid_flux_from_driver

    with computation(PARALLEL), interval(...):
        # Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
        large_scale_nonanvil_ice_flux = large_scale_nonanvil_ice_flux / (MAPL_GRAV * DT_MOIST)
        large_scale_nonanvil_liquid_flux = large_scale_nonanvil_liquid_flux / (MAPL_GRAV * DT_MOIST)

    with computation(PARALLEL), interval(1, None):
        # Redistribute precipitation fluxes for chemistry
        anvil_ice_flux = large_scale_nonanvil_ice_flux * min(
            1.0, max(convective_ice[0, 0, -1] / max(ice_for_radiation[0, 0, -1], 1.0e-8), 0.0)
        )
        large_scale_nonanvil_ice_flux = large_scale_nonanvil_ice_flux - anvil_ice_flux

        anvil_liquid_flux = large_scale_nonanvil_liquid_flux * min(
            1.0, max(convective_liquid[0, 0, -1] / max(liquid_for_radiation[0, 0, -1], 1.0e-8), 0.0)
        )
        large_scale_nonanvil_liquid_flux = large_scale_nonanvil_liquid_flux - anvil_liquid_flux

    with computation(PARALLEL), interval(0, -1):
        # cleanup suspended precipitation condensates
        radiation_rain = fix_negative_precip(radiation_rain)
        radiation_snow = fix_negative_precip(radiation_snow)
        radiation_graupel = fix_negative_precip(radiation_graupel)

    with computation(PARALLEL), interval(0, -1):
        vapor = radiation_vapor
        rain = radiation_rain
        snow = radiation_snow
        graupel = radiation_graupel


def fix_humidity(
    relative_humidity: FloatField,
    vapor: FloatField,
    t: FloatField,
    p_mb: FloatField,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
):
    with computation(PARALLEL), interval(...):
        qsat, _ = saturation_specific_humidity(t, p_mb * 100, ese, esx)
        relative_humidity = vapor / qsat


def fix_mixing_ratio(
    mixing_ratio: FloatField,
    mass: FloatField,
    k_sum_1: FloatFieldIJ,
    k_sum_2: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        k_sum_1 = k_sum_1 + (mixing_ratio * mass)

    with computation(PARALLEL), interval(...):
        if mixing_ratio < 0.0:
            mixing_ratio = 0.0

    with computation(FORWARD), interval(...):
        k_sum_2 = k_sum_2 + (mixing_ratio * mass)

    with computation(PARALLEL), interval(...):
        if k_sum_2 > 0.0:
            factor = (k_sum_2 - k_sum_1) / k_sum_2

        # reduce Q proportionally to the increase in TPW
        mixing_ratio = mixing_ratio * (1.0 - factor)

    # reset temporaries to zero for later uses
    with computation(FORWARD), interval(0, 1):
        k_sum_1 = 0
        k_sum_2 = 0


def minimum_mixing_ratio(
    mixing_ratio: FloatField,
    minimum: Float,
):
    with computation(PARALLEL), interval(...):
        mixing_ratio = min(mixing_ratio, minimum)


def fix_radii(
    convective_ice: FloatField,
    convective_liquid: FloatField,
    large_scale_ice: FloatField,
    large_scale_liquid: FloatField,
    ice_radius: FloatField,
    liquid_radius: FloatField,
):
    with computation(PARALLEL), interval(...):
        if large_scale_ice + convective_ice <= 0.0:
            ice_radius = 36.0e-6
        if large_scale_liquid + convective_liquid <= 0.0:
            liquid_radius = 14.0e-6


def update_rainwater_source(
    large_scale_rainwater_source: FloatField,
    drain_dt_macro: FloatField,
    drain_dt_micro: FloatField,
):
    with computation(PARALLEL), interval(...):
        large_scale_rainwater_source = drain_dt_macro + drain_dt_micro


def dissipative_ke_heating(
    mass: FloatField,
    u0: FloatField,
    v0: FloatField,
    du_dt_macro: FloatField,
    du_dt_micro: FloatField,
    dv_dt_macro: FloatField,
    dv_dt_micro: FloatField,
    t_tendency: FloatField,
    dts: FloatFieldIJ,
    fpi: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        # total KE dissipation estimate
        dts = dts - ((du_dt_macro + du_dt_micro) * u0 + (dv_dt_macro + dv_dt_micro) * v0) * mass
        # [sic] fpi needed for calcualtion of conversion to pot. energyintegrated
        ke = sqrt(
            (du_dt_macro + du_dt_micro) * (du_dt_macro + du_dt_micro)
            + (dv_dt_macro + dv_dt_micro) * (dv_dt_macro + dv_dt_micro)
        )
        fpi = fpi + ke * mass

    with computation(PARALLEL), interval(...):
        if fpi > 0.0:
            t_tendency = (ke / fpi) * dts * (1.0 / MAPL_CP)

    with computation(FORWARD), interval(0, 1):
        # reset temporaries to zero for future uses
        dts = 0
        fpi = 0


class Finalize:
    """
    Computes tendencies and output diagnostic fields using the following functions:

    redistribute_clouds: ensure in-cloud fields have physically reasonable values
    finalize_precip: ensure precipitaiton values are physically readonable and have the correct units
    radiation_coupling: send data to relevant fields for subsequent radiation calculations
    fix_humidity (conditional): recompute humidity post phase_change and driver components
    fix_mixing_ratio: remove negative mixing ratio values, maintaining mass conservation
    minimum_mixing_ratio: enforce minimum mixing ratios
    fix_radii: fix ice/liquid radii values where ice/liquid is absent
    update_tendencies: update microphysics tendencies
    large_scale_rainwater_source (conditional): update rainwater if output is enabled
    dissipative_ke_heating (condtional): compute temperature tendency due to friction

    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        GFDL_1M_config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
        update_tendencies,
    ):
        self.update_tendencies = update_tendencies
        self.GFDL_1M_config = GFDL_1M_config
        self.saturation_tables = saturation_tables

        # construct stencils
        self.redistribute_clouds = RedistributeClouds(
            stencil_factory,
        )

        self.finalize_precip = stencil_factory.from_dims_halo(
            func=finalize_precip,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            externals={
                "DT_MOIST": GFDL_1M_config.DT_MOIST,
            },
        )

        self.radiation_coupling = GFDL1MRadiationCoupling(
            stencil_factory,
            GFDL_1M_config.DO_QA,
            FAC_RL=GFDL_1M_config.FAC_RL,
            MIN_RL=GFDL_1M_config.MIN_RL,
            MAX_RL=GFDL_1M_config.MAX_RL,
            FAC_RI=GFDL_1M_config.FAC_RI,
            MIN_RI=GFDL_1M_config.MIN_RI,
            MAX_RI=GFDL_1M_config.MAX_RI,
        )

        self.fix_humidity = stencil_factory.from_dims_halo(
            func=fix_humidity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.fix_mixing_ratio = stencil_factory.from_dims_halo(
            func=fix_mixing_ratio,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.minimum_mixing_ratio = stencil_factory.from_dims_halo(
            func=minimum_mixing_ratio,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.fix_radii = stencil_factory.from_dims_halo(
            func=fix_radii,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.update_rainwater_source = stencil_factory.from_dims_halo(
            func=update_rainwater_source,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self.dissipative_ke_heating = stencil_factory.from_dims_halo(
            func=dissipative_ke_heating,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        t: FloatField,
        u: FloatField,
        v: FloatField,
        ice_concentration: FloatField,
        liquid_concentration: FloatField,
        mixing_ratios: MixingRatios,
        cloud_fractions: CloudFractions,
        masks: Masks,
        outputs: Outputs,
        temporaries: Temporaries,
        driver: MicrophysicsDriver,
    ):
        self.redistribute_clouds(
            cloud_fraction=outputs.radiation_cloud_fraction,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            liquid=outputs.radiation_liquid,
            convective_liquid=mixing_ratios.convective_liquid,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            ice=outputs.radiation_ice,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_ice=mixing_ratios.large_scale_ice,
            vapor=outputs.radiation_vapor,
            temperature=t,
        )

        self.finalize_precip(
            rain_from_driver=driver.outputs.rain,
            snow_from_driver=driver.outputs.snow,
            ice_from_driver=driver.outputs.ice,
            graupel_from_driver=driver.outputs.graupel,
            precipitated_rain=outputs.precipitated_rain,
            precipitated_snow=outputs.precipitated_snow,
            precipitated_ice=outputs.precipitated_ice,
            precipitated_graupel=outputs.precipitated_graupel,
            large_scale_precip=outputs.large_scale_precip,
            large_scale_snow=outputs.large_scale_snow,
            icefall=outputs.icefall,
            freezing_rainfall=outputs.freezing_rainfall,
            large_scale_nonanvil_ice_flux=outputs.large_scale_nonanvil_ice_flux,
            large_scale_nonanvil_ice_flux_from_driver=driver.outputs.m2_sol,
            large_scale_nonanvil_liquid_flux=outputs.large_scale_nonanvil_liquid_flux,
            large_scale_nonanvil_liquid_flux_from_driver=driver.outputs.m2_rain,
            convective_liquid=mixing_ratios.convective_liquid,
            liquid_for_radiation=outputs.radiation_liquid,
            anvil_liquid_flux=outputs.anvil_liquid_flux,
            convective_ice=mixing_ratios.convective_ice,
            ice_for_radiation=outputs.radiation_ice,
            anvil_ice_flux=outputs.anvil_ice_flux,
            evaporation=outputs.large_scale_nonanvil_precipitation_evaporation,
            evaporation_from_driver=driver.outputs.revap,
            sublimation=outputs.large_scale_nonanvil_precipitation_sublimation,
            sublimation_from_driver=driver.outputs.isubl,
            radiation_vapor=outputs.radiation_vapor,
            radiation_rain=outputs.radiation_rain,
            radiation_snow=outputs.radiation_snow,
            radiation_graupel=outputs.radiation_graupel,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
        )

        self.radiation_coupling(
            mixing_ratios.vapor,
            t,
            mixing_ratios.large_scale_liquid,
            mixing_ratios.large_scale_ice,
            cloud_fractions.large_scale,
            mixing_ratios.convective_liquid,
            mixing_ratios.convective_ice,
            cloud_fractions.convective,
            temporaries.p_mb,
            mixing_ratios.rain,
            mixing_ratios.snow,
            mixing_ratios.graupel,
            liquid_concentration,
            ice_concentration,
            outputs.radiation_vapor,
            outputs.radiation_liquid,
            outputs.radiation_ice,
            outputs.radiation_rain,
            outputs.radiation_snow,
            outputs.radiation_graupel,
            outputs.radiation_cloud_fraction,
            outputs.liquid_radius,
            outputs.ice_radius,
            outputs.relative_humidity_after_pdf,
            self.saturation_tables.ese,
            self.saturation_tables.esx,
        )

        if self.GFDL_1M_config.DO_QA is True:
            self.fix_humidity(
                relative_humidity=outputs.relative_humidity_after_pdf,
                vapor=mixing_ratios.vapor,
                t=t,
                p_mb=temporaries.p_mb,
                ese=self.saturation_tables.ese,
                esx=self.saturation_tables.esx,
            )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_vapor,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_liquid,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_ice,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_rain,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_snow,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_graupel,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.fix_mixing_ratio(
            mixing_ratio=outputs.radiation_cloud_fraction,
            mass=temporaries.mass,
            k_sum_1=temporaries.temporary_2d_1,
            k_sum_2=temporaries.temporary_2d_2,
        )

        self.minimum_mixing_ratio(
            mixing_ratio=outputs.radiation_liquid,
            minimum=Float(0.001),
        )

        self.minimum_mixing_ratio(
            mixing_ratio=outputs.radiation_ice,
            minimum=Float(0.001),
        )

        self.minimum_mixing_ratio(
            mixing_ratio=outputs.radiation_rain,
            minimum=Float(0.01),
        )

        self.minimum_mixing_ratio(
            mixing_ratio=outputs.radiation_snow,
            minimum=Float(0.01),
        )

        self.minimum_mixing_ratio(
            mixing_ratio=outputs.radiation_graupel,
            minimum=Float(0.01),
        )

        self.fix_radii(
            convective_ice=mixing_ratios.convective_ice,
            convective_liquid=mixing_ratios.convective_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            ice_radius=outputs.ice_radius,
            liquid_radius=outputs.liquid_radius,
        )

        self.update_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratios.vapor,
            rain=mixing_ratios.rain,
            snow=mixing_ratios.snow,
            graupel=mixing_ratios.graupel,
            convective_liquid=mixing_ratios.convective_liquid,
            convective_ice=mixing_ratios.convective_ice,
            large_scale_liquid=mixing_ratios.large_scale_liquid,
            large_scale_ice=mixing_ratios.large_scale_ice,
            convective_cloud_fraction=cloud_fractions.convective,
            large_scale_cloud_fraction=cloud_fractions.large_scale,
            du_dt=outputs.du_dt_micro,
            dv_dt=outputs.dv_dt_micro,
            dt_dt=outputs.dt_dt_micro,
            dvapor_dt=outputs.dvapor_dt_micro,
            dliquid_dt=outputs.dliquid_dt_micro,
            dice_dt=outputs.dice_dt_micro,
            dcloud_fraction_dt=outputs.dcloud_fraction_dt_micro,
            drain_dt=outputs.drain_dt_micro,
            dsnow_dt=outputs.dsnow_dt_micro,
            dgraupel_dt=outputs.dgraupel_dt_micro,
        )

        if outputs.large_scale_rainwater_source is not None:
            self.update_rainwater_source(
                outputs.large_scale_rainwater_source, outputs.drain_dt_macro, outputs.drain_dt_micro
            )

        if outputs.moist_friction_temperature_tendency is not None:
            self.dissipative_ke_heating(
                mass=temporaries.mass,
                u0=temporaries.u_unmodified,
                v0=temporaries.v_unmodified,
                du_dt_macro=outputs.du_dt_macro,
                du_dt_micro=outputs.du_dt_micro,
                dv_dt_macro=outputs.dv_dt_macro,
                dv_dt_micro=outputs.dv_dt_micro,
                t_tendency=outputs.moist_friction_temperature_tendency,
                dts=temporaries.temporary_2d_1,
                fpi=temporaries.temporary_2d_2,
            )

        if (
            outputs.simulated_reflectivity is not None
            or outputs.maximum_reflectivity is not None
            or outputs.one_km_agl_reflectivity is not None
            or outputs.echo_top_reflectivity is not None
            or outputs.minus_10c_reflectivity is not None
        ):
            raise NotImplementedError("Diagnostic radar output not implemented yet.")
