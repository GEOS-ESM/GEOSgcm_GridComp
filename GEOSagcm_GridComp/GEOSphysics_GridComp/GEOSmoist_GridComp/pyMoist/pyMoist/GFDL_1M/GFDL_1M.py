from ndsl import NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.GFDL_1M.config import GFDL1MConfig
from pyMoist.GFDL_1M.driver.driver import GFDL1MDriver
from pyMoist.GFDL_1M.finalize import GFDL1MFinalize
from pyMoist.GFDL_1M.locals import GFDL1MLocals
from pyMoist.GFDL_1M.PhaseChange.phase_change import PhaseChange
from pyMoist.GFDL_1M.setup import GFDL1MSetup
from pyMoist.GFDL_1M.state import GFDL1MState
from pyMoist.GFDL_1M.stencils import (
    get_total_concentration,
    prepare_radiation,
    prepare_tendencies,
    reset_micro_tendencies,
    update_radiation,
    update_tendencies,
)
from pyMoist.saturation_tables import get_saturation_vapor_pressure_table


class GFDL1M(NDSLRuntime):
    """
    GFDL Single Moment microphysics

    The primary purpose of this code is to compute macro/microphysical tendencies to be applied to state
    variables (p, t, wind, etc.). This code requires all fields to be preloaded with Fortran memory or
    otherwise supplied between the __init__ and __call__ steps.

    Performs the following functions to achieve this goal:
    __init__
        - initalize saturaiton vapor pressure tables, initalize temporary/output fields, construct stencils
        Arguments: StencilFactory, QuantityFactory, GFDL1MConfig

    __call__
        - setup: compute additional required fields, create pristine copies of input variables
        - phase_change: create new condensates, perform phase change operations
        - driver: precipitate condensates
        - finalize: compute tendencies, prepare fields to be returned to the larger model
        Arguments: none (data needs to be pre-loaded)
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GFDL1MConfig,
    ):
        super().__init__(stencil_factory)

        if stencil_factory.grid_indexing.n_halo != 0:
            raise ValueError("halo needs to be zero for GFDL Single Moment microphysics")

        # Initalize saturation tables
        saturation_tables = get_saturation_vapor_pressure_table(stencil_factory.backend)

        # Locals
        self._locals = GFDL1MLocals.make_locals(quantity_factory)

        # Build components
        self.prepare_tendencies = stencil_factory.from_dims_halo(
            func=prepare_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._setup = GFDL1MSetup(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            saturation_tables=saturation_tables,
        )

        self._phase_change = PhaseChange(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            saturation_tables=saturation_tables,
        )

        self._update_tendencies = stencil_factory.from_dims_halo(
            func=update_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
            },
        )

        self._reset_micro_tendencies = stencil_factory.from_dims_halo(
            func=reset_micro_tendencies,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._prepare_radiation = stencil_factory.from_dims_halo(
            func=prepare_radiation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._get_total_concentration = stencil_factory.from_dims_halo(
            func=get_total_concentration,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._driver = GFDL1MDriver(
            stencil_factory,
            quantity_factory,
            config,
        )

        self._update_radiation = stencil_factory.from_dims_halo(
            func=update_radiation,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "DT_MOIST": config.DT_MOIST,
            },
        )

        self._finalize = GFDL1MFinalize(
            stencil_factory=stencil_factory,
            quantity_factory=quantity_factory,
            config=config,
            saturation_tables=saturation_tables,
            update_tendencies=self._update_tendencies,
        )

    def __call__(
        self,
        state: GFDL1MState,
    ):
        # NOTE test GFDL_1M_setup: ❌ fails with tiny error (operational differences?)
        self._setup(
            p_interface=state.p_interface,
            z_interface=state.z_interface,
            u=state.u,
            v=state.v,
            t=state.t,
            lcl_height=state.lcl_height,
            lower_tropospheric_stability=state.lower_tropospheric_stability,
            estimated_inversion_strength=state.estimated_inversion_strength,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_rain=state.mixing_ratio.rain,
            mixing_ratio_snow=state.mixing_ratio.snow,
            mixing_ratio_graupel=state.mixing_ratio.graupel,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            cloud_fraction_convective=state.cloud_fraction.convective,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            shallow_convection_rain=state.shallow_convection_rain,
            shallow_convection_snow=state.shallow_convection_snow,
            dudt_macro=state.tendencies.dudt_macro,
            dvdt_macro=state.tendencies.dvdt_macro,
            dtdt_macro=state.tendencies.dtdt_macro,
            dvapordt_macro=state.tendencies.dvapordt_macro,
            dliquiddt_macro=state.tendencies.dliquiddt_macro,
            dicedt_macro=state.tendencies.dicedt_macro,
            dcloud_fractiondt_macro=state.tendencies.dcloud_fractiondt_macro,
            draindt_macro=state.tendencies.draindt_macro,
            dsnowdt_macro=state.tendencies.dsnowdt_macro,
            dgraupeldt_macro=state.tendencies.dgraupeldt_macro,
            local_p_mb=self._locals.p_mb,
            local_p_interface_mb=self._locals.p_interface_mb,
            local_edge_height_above_surface=self._locals.edge_height_above_surface,
            local_layer_height_above_surface=self._locals.layer_height_above_surface,
            local_layer_thickness=self._locals.layer_thickness,
            local_layer_thickness_negative=self._locals.layer_thickness_negative,
            local_dp=self._locals.dp,
            local_mass=self._locals.mass,
            local_mass_inverse=self._locals.mass_inverse,
            local_saturation_specific_humidity=self._locals.saturation_specific_humidity,
            local_dsaturation_specific_humidity=self._locals.dsaturation_specific_humidity,
            local_u_unmodified=self._locals.u_unmodified,
            local_v_unmodified=self._locals.v_unmodified,
            local_lcl_level=self._locals.lcl_level,
        )

        self._phase_change(
            t=state.t,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            cloud_fraction_convective=state.cloud_fraction.convective,
            concentration_ice=state.concentration.ice,
            concentration_liquid=state.concentration.liquid,
            relative_humidity_after_pdf=state.relative_humidity_after_pdf,
            estimated_inversion_strength=state.estimated_inversion_strength,
            area=state.area,
            critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
            pdf_iters=state.hydrostatic_pdf_iterations,
            cloud_liquid_evaporation=state.cloud_liquid_evaporation,
            cloud_ice_sublimation=state.cloud_ice_sublimation,
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            local_lcl_level=self._locals.lcl_level,
            local_p_mb=self._locals.p_mb,
            local_p_interface_mb=self._locals.p_interface_mb,
            local_saturation_specific_humidity=self._locals.saturation_specific_humidity,
        )

        self._update_tendencies(
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
            du_dt=state.tendencies.dudt_macro,
            dv_dt=state.tendencies.dvdt_macro,
            dt_dt=state.tendencies.dtdt_macro,
            dvapor_dt=state.tendencies.dvapordt_macro,
            dliquid_dt=state.tendencies.dliquiddt_macro,
            dice_dt=state.tendencies.dicedt_macro,
            dcloud_fraction_dt=state.tendencies.dcloud_fractiondt_macro,
            drain_dt=state.tendencies.draindt_macro,
            dsnow_dt=state.tendencies.dsnowdt_macro,
            dgraupel_dt=state.tendencies.dgraupeldt_macro,
        )

        # prepare microphysics tendencies
        self.prepare_tendencies(
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

        self._reset_micro_tendencies(
            dvapordt=self._locals.driver_tencencies.dvapordt,
            dliquiddt=self._locals.driver_tencencies.dliquiddt,
            draindt=self._locals.driver_tencencies.draindt,
            dicedt=self._locals.driver_tencencies.dicedt,
            dsnowdt=self._locals.driver_tencencies.dsnowdt,
            dgraupeldt=self._locals.driver_tencencies.dgraupeldt,
            dcloudfractiondt=self._locals.driver_tencencies.dcloudfractiondt,
            dtdt=self._locals.driver_tencencies.dtdt,
            dudt=self._locals.driver_tencencies.dudt,
            dvdt=self._locals.driver_tencencies.dvdt,
        )

        self._prepare_radiation(
            convective_cloud_fraction=state.cloud_fraction.convective,
            large_scale_cloud_fraction=state.cloud_fraction.large_scale,
            radiation_cloud_fraction=state.radiation_field.cloud_fraction,
            convective_liquid=state.mixing_ratio.convective_liquid,
            large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            radiation_liquid=state.radiation_field.liquid,
            convective_ice=state.mixing_ratio.convective_ice,
            large_scale_ice=state.mixing_ratio.large_scale_ice,
            radiation_ice=state.radiation_field.ice,
            vapor=state.mixing_ratio.vapor,
            radiation_vapor=state.radiation_field.vapor,
            rain=state.mixing_ratio.rain,
            radiation_rain=state.radiation_field.rain,
            snow=state.mixing_ratio.snow,
            radiation_snow=state.radiation_field.snow,
            graupel=state.mixing_ratio.graupel,
            radiation_graupel=state.radiation_field.graupel,
        )

        self._get_total_concentration(
            ice_concentration=state.concentration.ice,
            liquid_concentration=state.concentration.liquid,
            total_concentration=self._locals.total_concentration,
        )

        self._driver(
            t=state.t,
            u=state.u,
            v=state.v,
            w=state.vertical_motion.velocity,
            dz=self._locals.layer_thickness_negative,
            dp=self._locals.dp,
            area=state.area,
            land_fraction=state.land_fraction,
            convection_fraction=state.convection_fraction,
            surface_type=state.surface_type,
            estimated_inversion_strength=state.estimated_inversion_strength,
            critical_relative_humidity_for_pdf=state.critical_relative_humidity_for_pdf,
            vapor=state.radiation_field.vapor,
            liquid=state.radiation_field.liquid,
            rain=state.radiation_field.rain,
            ice=state.radiation_field.ice,
            snow=state.radiation_field.snow,
            graupel=state.radiation_field.graupel,
            cloud_fraction=state.radiation_field.cloud_fraction,
            total_concentration=self._locals.total_concentration,
            dvapordt=self._locals.driver_tencencies.dvapordt,
            dliquiddt=self._locals.driver_tencencies.dliquiddt,
            draindt=self._locals.driver_tencencies.draindt,
            dicedt=self._locals.driver_tencencies.dicedt,
            dsnowdt=self._locals.driver_tencencies.dsnowdt,
            dgraupeldt=self._locals.driver_tencencies.dgraupeldt,
            dcloudfractiondt=self._locals.driver_tencencies.dcloudfractiondt,
            dtdt=self._locals.driver_tencencies.dtdt,
            dudt=self._locals.driver_tencencies.dudt,
            dvdt=self._locals.driver_tencencies.dvdt,
            liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
            ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
            evaporation=state.non_anvil_large_scale.evaporation,
            sublimation=state.non_anvil_large_scale.sublimation,
            surface_precip_rain=state.precipitation_at_surface.rain,
            surface_precip_snow=state.precipitation_at_surface.snow,
            surface_precip_ice=state.precipitation_at_surface.ice,
            surface_precip_graupel=state.precipitation_at_surface.graupel,
        )

        self._update_radiation(
            t=state.t,
            u=state.u,
            v=state.v,
            radiation_cloud_fraction=state.radiation_field.cloud_fraction,
            radiation_ice=state.radiation_field.ice,
            radiation_liquid=state.radiation_field.liquid,
            radiation_vapor=state.radiation_field.vapor,
            radiation_rain=state.radiation_field.rain,
            radiation_snow=state.radiation_field.snow,
            radiation_graupel=state.radiation_field.graupel,
            dcloud_fraction_dt=self._locals.driver_tencencies.dcloudfractiondt,
            dtdt=self._locals.driver_tencencies.dtdt,
            dudt=self._locals.driver_tencencies.dudt,
            dvdt=self._locals.driver_tencencies.dvdt,
            dicedt=self._locals.driver_tencencies.dicedt,
            dliquiddt=self._locals.driver_tencencies.dliquiddt,
            dvapordt=self._locals.driver_tencencies.dvapordt,
            draindt=self._locals.driver_tencencies.draindt,
            dsnowdt=self._locals.driver_tencencies.dsnowdt,
            dgraupeldt=self._locals.driver_tencencies.dgraupeldt,
        )

        self._finalize(
            t=state.t,
            u=state.u,
            v=state.v,
            mixing_ratio_vapor=state.mixing_ratio.vapor,
            mixing_ratio_convective_liquid=state.mixing_ratio.convective_liquid,
            mixing_ratio_large_scale_liquid=state.mixing_ratio.large_scale_liquid,
            mixing_ratio_convective_ice=state.mixing_ratio.convective_ice,
            mixing_ratio_large_scale_ice=state.mixing_ratio.large_scale_ice,
            mixing_ratio_rain=state.mixing_ratio.rain,
            mixing_ratio_snow=state.mixing_ratio.snow,
            mixing_ratio_graupel=state.mixing_ratio.graupel,
            cloud_fraction_convective=state.cloud_fraction.convective,
            cloud_fraction_large_scale=state.cloud_fraction.large_scale,
            non_anvil_large_scale_precip=state.non_anvil_large_scale.precip,
            non_anvil_large_scale_snow=state.non_anvil_large_scale.snow,
            non_anvil_large_scale_ice_precip_flux=state.non_anvil_large_scale.ice_precip_flux,
            non_anvil_large_scale_liquid_precip_flux=state.non_anvil_large_scale.liquid_precip_flux,
            anvil_liquid_precip_flux=state.anvil.liquid_precip_flux,
            anvil_ice_precip_flux=state.anvil.ice_precip_flux,
            surface_rain=state.precipitation_at_surface.rain,
            surface_snow=state.precipitation_at_surface.snow,
            surface_ice=state.precipitation_at_surface.ice,
            surface_graupel=state.precipitation_at_surface.graupel,
            icefall=state.icefall,
            freezing_rainfall=state.freezing_rainfall,
            concentration_liquid=state.concentration.liquid,
            concentration_ice=state.concentration.ice,
            cloud_particle_effective_radius_liquid=state.cloud_particle_effective_radius.liquid,
            cloud_particle_effective_radius_ice=state.cloud_particle_effective_radius.ice,
            relative_humidity_after_pdf=state.relative_humidity_after_pdf,
            large_scale_rainwater_source=state.large_scale_rainwater_source,
            radiation_vapor=state.radiation_field.vapor,
            radiation_liquid=state.radiation_field.liquid,
            radiation_rain=state.radiation_field.rain,
            radiation_snow=state.radiation_field.snow,
            radiation_graupel=state.radiation_field.graupel,
            radiation_ice=state.radiation_field.ice,
            radiation_cloud_fraction=state.radiation_field.cloud_fraction,
            dudt_micro=state.tendencies.dudt_micro,
            dvdt_micro=state.tendencies.dvdt_micro,
            dtdt_micro=state.tendencies.dtdt_micro,
            dvapordt_micro=state.tendencies.dvapordt_micro,
            dliquiddt_micro=state.tendencies.dliquiddt_micro,
            dicedt_micro=state.tendencies.dicedt_micro,
            dcloud_fractiondt_micro=state.tendencies.dcloud_fractiondt_micro,
            draindt_micro=state.tendencies.draindt_micro,
            dsnowdt_micro=state.tendencies.dsnowdt_micro,
            dgraupeldt_micro=state.tendencies.dgraupeldt_micro,
            dudt_macro=state.tendencies.dudt_macro,
            dvdt_macro=state.tendencies.dvdt_macro,
            draindt_macro=state.tendencies.draindt_macro,
            dtdt_friction_pressure_weighted=state.tendencies.dtdt_friction_pressure_weighted,
            local_p_mb=self._locals.p_mb,
            local_mass=self._locals.mass,
            local_u_unmodified=self._locals.u_unmodified,
            local_v_unmodified=self._locals.v_unmodified,
        )
