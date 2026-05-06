from ndsl import NDSLRuntime, QuantityFactory, StencilFactory, ndsl_log
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, function, interval, sqrt
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from ndsl.stencils.basic_operations import copy

from pyMoist.constants import MAPL_CP, MAPL_GRAV
from pyMoist.microphysics.GFDL_1M.config import GFDL1MConfig
from pyMoist.microphysics.GFDL_1M.radiation_coupling import GFDL1MRadiationCoupling
from pyMoist.saturation_tables import GlobalTable_saturation_tables, SaturationVaporPressureTable, saturation_specific_humidity
from pyMoist.shared.redistribute_clouds import redistribute_clouds


@function
def fix_negative_precip(
    precip: Float,
):
    if precip < 1.0e-8:
        precip = 0.0

    return precip


def finalize_precip(
    precipitated_rain: FloatFieldIJ,
    precipitated_snow: FloatFieldIJ,
    precipitated_ice: FloatFieldIJ,
    precipitated_graupel: FloatFieldIJ,
    large_scale_precip: FloatFieldIJ,
    large_scale_snow: FloatFieldIJ,
    icefall: FloatFieldIJ,
    freezing_rainfall: FloatFieldIJ,
    large_scale_nonanvil_ice_flux: FloatField,
    large_scale_nonanvil_liquid_flux: FloatField,
    convective_liquid: FloatField,
    liquid_for_radiation: FloatField,
    anvil_liquid_flux: FloatField,
    convective_ice: FloatField,
    ice_for_radiation: FloatField,
    anvil_ice_flux: FloatField,
    radiation_vapor: FloatField,
    radiation_rain: FloatField,
    radiation_snow: FloatField,
    radiation_graupel: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(FORWARD), interval(0, 1):
        # send precip diagnostics back to the rest of the model
        # and convert from mm/day to kg m-2 s-1
        precipitated_rain = max(precipitated_rain / 86400.0, 0.0)
        precipitated_snow = max(precipitated_snow / 86400.0, 0.0)
        precipitated_ice = max(precipitated_ice / 86400.0, 0.0)
        precipitated_graupel = max(precipitated_graupel / 86400.0, 0.0)
        # Fill GEOS precip diagnostics
        large_scale_precip = precipitated_rain
        large_scale_snow = precipitated_snow
        icefall = precipitated_ice + precipitated_graupel
        freezing_rainfall = 0.0

    with computation(FORWARD), interval(...):
        # Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
        large_scale_nonanvil_ice_flux[0, 0, 1] = large_scale_nonanvil_ice_flux[0, 0, 1] / (MAPL_GRAV * DT_MOIST)
        large_scale_nonanvil_liquid_flux[0, 0, 1] = large_scale_nonanvil_liquid_flux[0, 0, 1] / (MAPL_GRAV * DT_MOIST)

    with computation(FORWARD), interval(...):
        # Redistribute precipitation fluxes for chemistry
        anvil_ice_flux[0, 0, 1] = large_scale_nonanvil_ice_flux[0, 0, 1] * min(
            1.0,
            max(convective_ice / max(ice_for_radiation, 1.0e-8), 0.0),
        )
        large_scale_nonanvil_ice_flux[0, 0, 1] = large_scale_nonanvil_ice_flux[0, 0, 1] - anvil_ice_flux[0, 0, 1]

        anvil_liquid_flux[0, 0, 1] = large_scale_nonanvil_liquid_flux[0, 0, 1] * min(
            1.0,
            max(
                convective_liquid / max(liquid_for_radiation, 1.0e-8),
                0.0,
            ),
        )
        large_scale_nonanvil_liquid_flux[0, 0, 1] = large_scale_nonanvil_liquid_flux[0, 0, 1] - anvil_liquid_flux[0, 0, 1]

    with computation(PARALLEL), interval(...):
        # cleanup suspended precipitation condensates
        radiation_rain = fix_negative_precip(radiation_rain)
        radiation_snow = fix_negative_precip(radiation_snow)
        radiation_graupel = fix_negative_precip(radiation_graupel)

    with computation(PARALLEL), interval(...):
        vapor = radiation_vapor
        rain = radiation_rain
        snow = radiation_snow
        graupel = radiation_graupel


def fix_humidity(
    relative_humidity: FloatField,
    vapor: FloatField,
    t: FloatField,
    p_mb: FloatField,
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
):
    with computation(PARALLEL), interval(...):
        qsat, _ = saturation_specific_humidity(t, p_mb * 100.0, ese, esx)
        relative_humidity = vapor / qsat


def fix_mixing_ratio(
    mixing_ratio: FloatField,
    mass: FloatField,
):
    # predefine two FloatFieldIJ internal fields
    with computation(FORWARD), interval(0, 1):
        k_sum_1: FloatFieldIJ = 0.0
        k_sum_2: FloatFieldIJ = 0.0

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
):
    # predefine two FloatFieldIJ internal fields
    with computation(FORWARD), interval(0, 1):
        dts: FloatFieldIJ = 0.0
        fpi: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        # total KE dissipation estimate
        dts = dts - ((du_dt_macro + du_dt_micro) * u0 + (dv_dt_macro + dv_dt_micro) * v0) * mass
        # [sic] fpi needed for calculation of conversion to pot. energy integrated
        ke = sqrt((du_dt_macro + du_dt_micro) * (du_dt_macro + du_dt_micro) + (dv_dt_macro + dv_dt_micro) * (dv_dt_macro + dv_dt_micro))
        fpi = fpi + ke * mass

    with computation(PARALLEL), interval(...):
        if fpi > 0.0:
            t_tendency = (ke / fpi) * dts * (1.0 / MAPL_CP)

    with computation(FORWARD), interval(0, 1):
        # reset temporaries to zero for future uses
        dts = 0
        fpi = 0


class GFDL1MFinalize(NDSLRuntime):
    """
    Computes tendencies and output diagnostic fields using the following functions:

    redistribute_clouds: ensure in-cloud fields have physically reasonable values
    finalize_precip: ensure precipitation values are physically reasonable and have the correct units
    radiation_coupling: send data to relevant fields for subsequent radiation calculations
    fix_humidity (conditional): recompute humidity post phase_change and driver components
    fix_mixing_ratio: remove negative mixing ratio values, maintaining mass conservation
    minimum_mixing_ratio: enforce minimum mixing ratios
    fix_radii: fix ice/liquid radii values where ice/liquid is absent
    update_tendencies: update microphysics tendencies
    large_scale_rainwater_source (conditional): update rainwater if output is enabled
    dissipative_ke_heating (conditional): compute temperature tendency due to friction

    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
        update_tendencies,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make the config, pre-build stencil, and saturation tables visible at runtime
        self.config = config
        self.update_tendencies = update_tendencies
        self.saturation_tables = saturation_tables

        # construct stencils
        self._redistribute_clouds = stencil_factory.from_dims_halo(
            func=redistribute_clouds,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._finalize_precip = stencil_factory.from_dims_halo(
            func=finalize_precip,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )

        self._radiation_coupling = GFDL1MRadiationCoupling(
            stencil_factory=stencil_factory,
            config=config,
            saturation_tables=saturation_tables,
        )

        self._fix_humidity = stencil_factory.from_dims_halo(
            func=fix_humidity,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._fix_mixing_ratio = stencil_factory.from_dims_halo(
            func=fix_mixing_ratio,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._minimum_mixing_ratio = stencil_factory.from_dims_halo(
            func=minimum_mixing_ratio,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._fix_radii = stencil_factory.from_dims_halo(
            func=fix_radii,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._update_rainwater_source = stencil_factory.from_dims_halo(
            func=update_rainwater_source,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._dissipative_ke_heating = stencil_factory.from_dims_halo(
            func=dissipative_ke_heating,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._copy = stencil_factory.from_dims_halo(
            func=copy,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        # Dev NOTE: this is an orchestration workaround. Direct call to
        #           `self.saturation_tables.X` fails closure capture for
        #           argument reconstruction at call time
        self._ese = self.saturation_tables.ese
        self._esx = self.saturation_tables.esx

    def __call__(
        self,
        t,
        u,
        v,
        mixing_ratio_vapor,
        mixing_ratio_convective_liquid,
        mixing_ratio_large_scale_liquid,
        mixing_ratio_convective_ice,
        mixing_ratio_large_scale_ice,
        mixing_ratio_rain,
        mixing_ratio_snow,
        mixing_ratio_graupel,
        cloud_fraction_convective,
        cloud_fraction_large_scale,
        non_anvil_large_scale_precip,
        non_anvil_large_scale_snow,
        non_anvil_large_scale_ice_precip_flux,
        non_anvil_large_scale_liquid_precip_flux,
        anvil_liquid_precip_flux,
        anvil_ice_precip_flux,
        surface_rain,
        surface_snow,
        surface_ice,
        surface_graupel,
        icefall,
        freezing_rainfall,
        concentration_liquid,
        concentration_ice,
        cloud_particle_effective_radius_liquid,
        cloud_particle_effective_radius_ice,
        relative_humidity_after_pdf,
        large_scale_rainwater_source,
        radiation_vapor,
        radiation_liquid,
        radiation_rain,
        radiation_snow,
        radiation_graupel,
        radiation_ice,
        radiation_cloud_fraction,
        dudt_micro,
        dvdt_micro,
        dtdt_micro,
        dvapordt_micro,
        dliquiddt_micro,
        dicedt_micro,
        dcloud_fractiondt_micro,
        draindt_micro,
        dsnowdt_micro,
        dgraupeldt_micro,
        dudt_macro,
        dvdt_macro,
        draindt_macro,
        dtdt_friction_pressure_weighted,
        local_p_mb,
        local_mass,
        local_u_unmodified,
        local_v_unmodified,
        simulated_reflectivity,
        maximum_composite_reflectivity,
        base_1km_agl_reflectivity,
        echo_top_reflectivity,
        minus_10c_reflectivity,
        mass_fraction_suspended_rain,
        mass_fraction_suspended_snow,
        mass_fraction_suspended_graupel,
    ):
        self._redistribute_clouds(
            cloud_fraction=radiation_cloud_fraction,
            convective_cloud_fraction=cloud_fraction_convective,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            liquid=radiation_liquid,
            convective_liquid=mixing_ratio_convective_liquid,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            ice=radiation_ice,
            convective_ice=mixing_ratio_convective_ice,
            large_scale_ice=mixing_ratio_large_scale_ice,
            vapor=radiation_vapor,
            temperature=t,
        )

        self._finalize_precip(
            precipitated_rain=surface_rain,
            precipitated_snow=surface_snow,
            precipitated_ice=surface_ice,
            precipitated_graupel=surface_graupel,
            large_scale_precip=non_anvil_large_scale_precip,
            large_scale_snow=non_anvil_large_scale_snow,
            icefall=icefall,
            freezing_rainfall=freezing_rainfall,
            large_scale_nonanvil_ice_flux=non_anvil_large_scale_ice_precip_flux,
            large_scale_nonanvil_liquid_flux=non_anvil_large_scale_liquid_precip_flux,
            convective_liquid=mixing_ratio_convective_liquid,
            liquid_for_radiation=radiation_liquid,
            anvil_liquid_flux=anvil_liquid_precip_flux,
            convective_ice=mixing_ratio_convective_ice,
            ice_for_radiation=radiation_ice,
            anvil_ice_flux=anvil_ice_precip_flux,
            radiation_vapor=radiation_vapor,
            radiation_rain=radiation_rain,
            radiation_snow=radiation_snow,
            radiation_graupel=radiation_graupel,
            vapor=mixing_ratio_vapor,
            rain=mixing_ratio_rain,
            snow=mixing_ratio_snow,
            graupel=mixing_ratio_graupel,
        )

        self._radiation_coupling(
            t=t,
            mixing_ratio_vapor=mixing_ratio_vapor,
            mixing_ratio_large_scale_liquid=mixing_ratio_large_scale_liquid,
            mixing_ratio_large_scale_ice=mixing_ratio_large_scale_ice,
            mixing_ratio_convective_liquid=mixing_ratio_convective_liquid,
            mixing_ratio_rain=mixing_ratio_rain,
            mixing_ratio_snow=mixing_ratio_snow,
            mixing_ratio_graupel=mixing_ratio_graupel,
            mixing_ratio_convective_ice=mixing_ratio_convective_ice,
            cloud_fraction_large_scale=cloud_fraction_large_scale,
            cloud_fraction_convective=cloud_fraction_convective,
            concentration_liquid=concentration_liquid,
            concentration_ice=concentration_ice,
            liquid_radius=cloud_particle_effective_radius_liquid,
            ice_radius=cloud_particle_effective_radius_ice,
            relative_humidity_after_pdf=relative_humidity_after_pdf,
            radiation_vapor=radiation_vapor,
            radiation_liquid=radiation_liquid,
            radiation_ice=radiation_ice,
            radiation_rain=radiation_rain,
            radiation_snow=radiation_snow,
            radiation_graupel=radiation_graupel,
            radiation_cloud_fraction=radiation_cloud_fraction,
            local_p_mb=local_p_mb,
        )

        if self.config.DO_QA is True:
            self._fix_humidity(
                relative_humidity=relative_humidity_after_pdf,
                vapor=mixing_ratio_vapor,
                t=t,
                p_mb=local_p_mb,
                ese=self._ese,
                esx=self._esx,
            )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_vapor,
            mass=local_mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_liquid,
            mass=local_mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_ice,
            mass=local_mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_rain,
            mass=local_mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_snow,
            mass=local_mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_graupel,
            mass=local_mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=radiation_cloud_fraction,
            mass=local_mass,
        )

        self._minimum_mixing_ratio(
            mixing_ratio=radiation_liquid,
            minimum=Float(0.001),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=radiation_ice,
            minimum=Float(0.001),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=radiation_rain,
            minimum=Float(0.01),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=radiation_snow,
            minimum=Float(0.01),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=radiation_graupel,
            minimum=Float(0.01),
        )

        self._fix_radii(
            convective_ice=mixing_ratio_convective_ice,
            convective_liquid=mixing_ratio_convective_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            ice_radius=cloud_particle_effective_radius_ice,
            liquid_radius=cloud_particle_effective_radius_liquid,
        )

        self.update_tendencies(
            u=u,
            v=v,
            t=t,
            vapor=mixing_ratio_vapor,
            rain=mixing_ratio_rain,
            snow=mixing_ratio_snow,
            graupel=mixing_ratio_graupel,
            convective_liquid=mixing_ratio_convective_liquid,
            convective_ice=mixing_ratio_convective_ice,
            large_scale_liquid=mixing_ratio_large_scale_liquid,
            large_scale_ice=mixing_ratio_large_scale_ice,
            convective_cloud_fraction=cloud_fraction_convective,
            large_scale_cloud_fraction=cloud_fraction_large_scale,
            du_dt=dudt_micro,
            dv_dt=dvdt_micro,
            dt_dt=dtdt_micro,
            dvapor_dt=dvapordt_micro,
            dliquid_dt=dliquiddt_micro,
            dice_dt=dicedt_micro,
            dcloud_fraction_dt=dcloud_fractiondt_micro,
            drain_dt=draindt_micro,
            dsnow_dt=dsnowdt_micro,
            dgraupel_dt=dgraupeldt_micro,
        )

        if large_scale_rainwater_source is not None:
            self._update_rainwater_source(
                large_scale_rainwater_source,
                draindt_macro,
                draindt_micro,
            )

        if dtdt_friction_pressure_weighted is not None:
            self._dissipative_ke_heating(
                mass=local_mass,
                u0=local_u_unmodified,
                v0=local_v_unmodified,
                du_dt_macro=dudt_macro,
                du_dt_micro=dudt_micro,
                dv_dt_macro=dvdt_macro,
                dv_dt_micro=dvdt_micro,
                t_tendency=dtdt_friction_pressure_weighted,
            )

        if (
            simulated_reflectivity is not None
            or maximum_composite_reflectivity is not None
            or base_1km_agl_reflectivity is not None
            or echo_top_reflectivity is not None
            or minus_10c_reflectivity is not None
        ):
            ndsl_log.warning("Diagnostic radar output not implemented yet.")

        # new code from v11.8.1, is not tested (translate tests are based on data from v11.5.2)
        if mass_fraction_suspended_rain is not None:
            self._copy(mixing_ratio_rain, mass_fraction_suspended_rain)

        if mass_fraction_suspended_snow is not None:
            self._copy(mixing_ratio_snow, mass_fraction_suspended_snow)

        if mass_fraction_suspended_graupel is not None:
            self._copy(mixing_ratio_graupel, mass_fraction_suspended_graupel)
