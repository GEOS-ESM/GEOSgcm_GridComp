from ndsl import StencilFactory, QuantityFactory, NDSLRuntime, ndsl_log, orchestrate
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, computation, function, interval, sqrt
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ
from pyMoist.constants import MAPL_CP, MAPL_GRAV
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import MicrophysicsDriver
from pyMoist.GFDL_1M.masks import Masks
from pyMoist.GFDL_1M.state_old import MicrophysicState, Outputs
from pyMoist.GFDL_1M.temporaries import Temporaries
from pyMoist.radiation_coupling import GFDL1MRadiationCoupling
from pyMoist.redistribute_clouds import redistribute_clouds
from pyMoist.saturation_tables import (
    GlobalTable_saturation_tables,
    SaturationVaporPressureTable,
    saturation_specific_humidity,
)
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.locals import GFDL1MLocals


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
    """

    Dev note on `PFL_LS`, `PFI_LS`, `PFL_AN`, `PFI_AN `:
        Fortran use large_scale_nonanvil_ice_flux(PFI_LS) & anvil_ice_flux(PFI_AN) etc.,
        as Z_INTERFACE_DIMS fields, but level == 0 is never touched.
        It's reset every time to 0 and calculation are done on [1:] everywhere (before and in this code).
        Therefore, it's not a Z_INTERFACE_DIM and we treat it as a Z_DIM.
    """
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

    with computation(PARALLEL), interval(...):
        # Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
        large_scale_nonanvil_ice_flux = large_scale_nonanvil_ice_flux / (MAPL_GRAV * DT_MOIST)
        large_scale_nonanvil_liquid_flux = large_scale_nonanvil_liquid_flux / (MAPL_GRAV * DT_MOIST)

    with computation(PARALLEL), interval(1, None):
        # Redistribute precipitation fluxes for chemistry
        anvil_ice_flux = large_scale_nonanvil_ice_flux * min(
            1.0,
            max(convective_ice / max(ice_for_radiation, 1.0e-8), 0.0),
        )
        large_scale_nonanvil_ice_flux = large_scale_nonanvil_ice_flux - anvil_ice_flux

        anvil_liquid_flux = large_scale_nonanvil_liquid_flux * min(
            1.0,
            max(
                convective_liquid / max(liquid_for_radiation, 1.0e-8),
                0.0,
            ),
        )
        large_scale_nonanvil_liquid_flux = large_scale_nonanvil_liquid_flux - anvil_liquid_flux

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
        qsat, _ = saturation_specific_humidity(t, p_mb * 100, ese, esx)
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
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
        saturation_tables: SaturationVaporPressureTable,
        update_tendencies,
    ):
        # NOTE disabled because state has changed
        # orchestrate(
        #     obj=self,
        #     config=stencil_factory.config.dace_config,
        #     dace_compiletime_args=[
        #         "mixing_ratios",
        #         "cloud_fractions",
        #         "masks",
        #         "outputs",
        #         "temporaries",
        #         "driver",
        #     ],
        # )

        # make the config, pre-build stencil, and saturation tables visible at runtime
        self.config = config
        self.update_tendencies = update_tendencies
        self.saturation_tables = saturation_tables

        # construct stencils
        self._redistribute_clouds = stencil_factory.from_dims_halo(
            func=redistribute_clouds,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._finalize_precip = stencil_factory.from_dims_halo(
            func=finalize_precip,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"DT_MOIST": config.DT_MOIST},
        )

        self._radiation_coupling = GFDL1MRadiationCoupling(
            stencil_factory=stencil_factory,
            config=config,
            saturation_tables=saturation_tables,
        )

        self._fix_humidity = stencil_factory.from_dims_halo(
            func=fix_humidity,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._fix_mixing_ratio = stencil_factory.from_dims_halo(
            func=fix_mixing_ratio,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._minimum_mixing_ratio = stencil_factory.from_dims_halo(
            func=minimum_mixing_ratio,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._fix_radii = stencil_factory.from_dims_halo(
            func=fix_radii,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_rainwater_source = stencil_factory.from_dims_halo(
            func=update_rainwater_source,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._dissipative_ke_heating = stencil_factory.from_dims_halo(
            func=dissipative_ke_heating,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        # Dev NOTE: this is an orchestration workaround. Direct call to
        #           `self.saturation_tables.X` fails closure capture for
        #           argument reconstruction at call time
        self._ese = self.saturation_tables.ese
        self._esx = self.saturation_tables.esx

    def __call__(
        self,
        state: GFDL1MState,
        locals: GFDL1MLocals,
    ):
        self._redistribute_clouds(
            cloud_fraction=state.radiation_field.cloud_fraction,
            convective_cloud_fraction=state.cloud_fraction.convective,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            liquid=state.radiation_field.liquid,
            convective_liquid=state.mixing_ratio.convective_liquid,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            ice=state.radiation_field.ice,
            convective_ice=state.mixing_ratio.convective_ice,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            vapor=state.radiation_field.vapor,
            temperature=state.t,
        )

        self._finalize_precip(
            precipitated_rain=state.precipitation_at_surface.rain,
            precipitated_snow=state.precipitation_at_surface.snow,
            precipitated_ice=state.precipitation_at_surface.ice,
            precipitated_graupel=state.precipitation_at_surface.graupel,
            large_scale_precip=state.non_anvil_large_scale.precip,
            large_scale_snow=state.non_anvil_large_scale.snow,
            icefall=state.icefall,
            freezing_rainfall=state.freezing_rainfall,
            large_scale_nonanvil_ice_flux=state.non_anvil_large_scale.ice_precip_flux,
            large_scale_nonanvil_liquid_flux=state.non_anvil_large_scale.liquid_precip_flux,
            convective_liquid=state.mixing_ratio.convective_liquid,
            liquid_for_radiation=state.radiation_field.liquid,
            anvil_liquid_flux=state.anvil.liquid_precip_flux,
            convective_ice=state.mixing_ratio.convective_ice,
            ice_for_radiation=state.radiation_field.ice,
            anvil_ice_flux=state.anvil.ice_precip_flux,
            radiation_vapor=state.radiation_field.vapor,
            radiation_rain=state.radiation_field.rain,
            radiation_snow=state.radiation_field.snow,
            radiation_graupel=state.radiation_field.graupel,
            vapor=state.mixing_ratio.vapor,
            rain=state.mixing_ratio.rain,
            snow=state.mixing_ratio.snow,
            graupel=state.mixing_ratio.graupel,
        )

        self._radiation_coupling(state=state, locals=locals)

        if self.config.DO_QA is True:
            self._fix_humidity(
                relative_humidity=state.relative_humidity_after_pdf,
                vapor=state.mixing_ratio.vapor,
                t=state.t,
                p_mb=locals.p_mb,
                ese=self._ese,
                esx=self._esx,
            )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.vapor,
            mass=locals.mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.liquid,
            mass=locals.mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.ice,
            mass=locals.mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.rain,
            mass=locals.mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.snow,
            mass=locals.mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.graupel,
            mass=locals.mass,
        )

        self._fix_mixing_ratio(
            mixing_ratio=state.radiation_field.cloud_fraction,
            mass=locals.mass,
        )

        self._minimum_mixing_ratio(
            mixing_ratio=state.radiation_field.liquid,
            minimum=Float(0.001),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=state.radiation_field.ice,
            minimum=Float(0.001),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=state.radiation_field.rain,
            minimum=Float(0.01),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=state.radiation_field.snow,
            minimum=Float(0.01),
        )

        self._minimum_mixing_ratio(
            mixing_ratio=state.radiation_field.graupel,
            minimum=Float(0.01),
        )

        self._fix_radii(
            convective_ice=state.mixing_ratio.convective_ice,
            convective_liquid=state.mixing_ratio.convective_liquid,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            ice_radius=state.cloud_particle_effective_radius.ice,
            liquid_radius=state.cloud_particle_effective_radius.liquid,
        )

        self.update_tendencies(
            u=state.u,
            v=state.v,
            t=state.t,
            vapor=state.mixing_ratio.vapor,
            rain=state.mixing_ratio.rain,
            snow=state.mixing_ratio.snow,
            graupel=state.mixing_ratio.graupel,
            convective_liquid=state.mixing_ratio.convective_liquid,
            convective_ice=state.mixing_ratio.convective_ice,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            convective_cloud_fraction=state.cloud_fraction.convective,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            du_dt=state.tendencies.dudt_micro,
            dv_dt=state.tendencies.dvdt_micro,
            dt_dt=state.tendencies.dtdt_micro,
            dvapor_dt=state.tendencies.dvapordt_micro,
            dliquid_dt=state.tendencies.dliquiddt_micro,
            dice_dt=state.tendencies.dicedt_micro,
            dcloud_fraction_dt=state.tendencies.dcloud_fractiondt_micro,
            drain_dt=state.tendencies.draindt_micro,
            dsnow_dt=state.tendencies.dsnowdt_micro,
            dgraupel_dt=state.tendencies.dgraupeldt_micro,
        )

        # NOTE need a better way to figure out if this is associated, this will never be none with the new state
        if state.large_scale_rainwater_source is not None:
            self._update_rainwater_source(
                state.large_scale_rainwater_source,
                state.tendencies.draindt_macro,
                state.tendencies.draindt_micro,
            )

        # NOTE need a better way to figure out if this is associated, this will never be none with the new state
        if state.tendencies.dtdt_friction_pressure_weighted is not None:
            self._dissipative_ke_heating(
                mass=locals.mass,
                u0=locals.u_unmodified,
                v0=locals.v_unmodified,
                du_dt_macro=state.tendencies.dudt_macro,
                du_dt_micro=state.tendencies.dudt_micro,
                dv_dt_macro=state.tendencies.dvdt_macro,
                dv_dt_micro=state.tendencies.dvdt_micro,
                t_tendency=state.tendencies.dtdt_friction_pressure_weighted,
            )

        # NOTE DISABLED DURING DEVELOPMENT.
        # NOTE need a better way to figure out if this is associated, this will never be none with the new state
        # if (
        #     outputs.simulated_reflectivity is not None
        #     or outputs.maximum_reflectivity is not None
        #     or outputs.one_km_agl_reflectivity is not None
        #     or outputs.echo_top_reflectivity is not None
        #     or outputs.minus_10c_reflectivity is not None
        # ):
        #     ndsl_log.warning("Diagnostic radar output not implemented yet.")
