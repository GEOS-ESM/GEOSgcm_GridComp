from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist import GFDL1M, GFDL1MConfig, GFDL1MState
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGFDL_1M(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        # NOTE not all fields are associated in v11.5.2. fields which are not associated are not serialized,
        # and not output by the parameterization, and have therefore been manually disabled within this test
        self.in_vars["data_vars"] = {
            "precipitation_at_surface_deep_convective_precipitation": {},
            "precipitation_at_surface_anvil_precipitation": {},
            "precipitation_at_surface_shallow_convective_precipitation": {},
            "precipitation_at_surface_deep_convective_snow": {},
            "precipitation_at_surface_anvil_snow": {},
            "precipitation_at_surface_shallow_convective_snow": {},
            # "lcl_height": {}, not associated in v11.5.2
            "shallow_convection_rain": {},
            "shallow_convection_snow": {},
            # "large_scale_rainwater_source": {}, not associated in v11.5.2
            "tendencies_dtdt_friction_pressure_weighted": {},
            "mixing_ratio_vapor": {},
            "mixing_ratio_rain": {},
            "mixing_ratio_snow": {},
            "mixing_ratio_graupel": {},
            "mixing_ratio_large_scale_liquid": {},
            "mixing_ratio_large_scale_ice": {},
            "mixing_ratio_convective_liquid": {},
            "mixing_ratio_convective_ice": {},
            "cloud_fraction_large_scale": {},
            "cloud_fraction_convective": {},
            "concentration_liquid": {},
            "concentration_ice": {},
            "area": {},
            "p_interface": {},
            "z_interface": {},
            "t": {},
            "u": {},
            "v": {},
            "land_fraction": {},
            "covariance_liquid_water_static_energy_and_total_water_specific_humidity": {},
            "omega": {},
            "pdf_first_plume_fractional_area": {},
            "vertical_motion_velocity": {},
            "vertical_motion_variance": {},
            "vertical_motion_third_moment": {},
            "liquid_water_static_energy_flux": {},
            "liquid_water_static_energy_variance": {},
            "liquid_water_static_energy_third_moment": {},
            "total_water_flux": {},
            "total_water_variance": {},
            "total_water_third_moment": {},
            "lower_tropospheric_stability": {},
            "estimated_inversion_strength": {},
            "tendencies_dcloud_fractiondt_macro": {},
            "tendencies_dvapordt_macro": {},
            "tendencies_dicedt_macro": {},
            "tendencies_dliquiddt_macro": {},
            "tendencies_draindt_macro": {},
            "tendencies_dgraupeldt_macro": {},
            "tendencies_dsnowdt_macro": {},
            "tendencies_dudt_macro": {},
            "tendencies_dvdt_macro": {},
            "tendencies_dtdt_macro": {},
            "convection_fraction": {},
            "surface_type": {},
            "cloud_liquid_evaporation": {},
            "cloud_ice_sublimation": {},
            "icefall": {},
            "freezing_rainfall": {},
            "relative_humidity_after_pdf": {},
            "buoyancy_flux": {},
            "liquid_water_flux": {},
            "hydrostatic_pdf_iterations": {},
            "radiation_field_cloud_fraction": {},
            "radiation_field_vapor": {},
            "radiation_field_liquid": {},
            "radiation_field_ice": {},
            "radiation_field_rain": {},
            "radiation_field_snow": {},
            "radiation_field_graupel": {},
            "cloud_particle_effective_radius_liquid": {},
            "cloud_particle_effective_radius_ice": {},
            "precipitation_at_surface_rain": {},
            "precipitation_at_surface_snow": {},
            "precipitation_at_surface_ice": {},
            "precipitation_at_surface_graupel": {},
            "non_anvil_large_scale_precip": {},
            "non_anvil_large_scale_snow": {},
            "non_anvil_large_scale_evaporation": {},
            "non_anvil_large_scale_sublimation": {},
            "non_anvil_large_scale_liquid_precip_flux": {},
            "non_anvil_large_scale_ice_precip_flux": {},
            "anvil_liquid_precip_flux": {},
            "anvil_ice_precip_flux": {},
            "critical_relative_humidity_for_pdf": {},
            "tendencies_dcloud_fractiondt_micro": {},
            "tendencies_dvapordt_micro": {},
            "tendencies_dicedt_micro": {},
            "tendencies_dliquiddt_micro": {},
            "tendencies_draindt_micro": {},
            "tendencies_dgraupeldt_micro": {},
            "tendencies_dsnowdt_micro": {},
            "tendencies_dudt_micro": {},
            "tendencies_dvdt_micro": {},
            "tendencies_dtdt_micro": {},
            # "radar_simulated_reflectivity": {}, not associated in v11.5.2
            # "radar_maximum_composite_reflectivity": {}, not associated in v11.5.2
            # "radar_base_1km_agl_reflectivity": {}, not associated in v11.5.2
            # "radar_echo_top_reflectivity": {}, not associated in v11.5.2
            # "radar_minus_10c_reflectivity": {}, not associated in v11.5.2
        }

        self.out_vars = self.in_vars["data_vars"].copy()

        # Initialize saturation tables
        self.saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GFDL_1M-constants")

    def compute(self, inputs):
        # initialize constants
        config = GFDL1MConfig(**self.constants)

        # initialize state
        state = GFDL1MState.zeros(self.quantity_factory)

        # fill state with input data
        state.precipitation_at_surface.deep_convective_precipitation.field[:] = inputs["precipitation_at_surface_deep_convective_precipitation"]
        state.precipitation_at_surface.anvil_precipitation.field[:] = inputs["precipitation_at_surface_anvil_precipitation"]
        state.precipitation_at_surface.shallow_convective_precipitation.field[:] = inputs["precipitation_at_surface_shallow_convective_precipitation"]
        state.precipitation_at_surface.deep_convective_snow.field[:] = inputs["precipitation_at_surface_deep_convective_snow"]
        state.precipitation_at_surface.anvil_snow.field[:] = inputs["precipitation_at_surface_anvil_snow"]
        state.precipitation_at_surface.shallow_convective_snow.field[:] = inputs["precipitation_at_surface_shallow_convective_snow"]
        # state.lcl_height.field[:] = inputs["lcl_height"]
        state.shallow_convection_rain.field[:] = inputs["shallow_convection_rain"]
        state.shallow_convection_snow.field[:] = inputs["shallow_convection_snow"]
        # state.large_scale_rainwater_source.field[:] = inputs["large_scale_rainwater_source"]
        state.tendencies.dtdt_friction_pressure_weighted.field[:] = inputs["tendencies_dtdt_friction_pressure_weighted"]
        state.mixing_ratio.vapor.field[:] = inputs["mixing_ratio_vapor"]
        state.mixing_ratio.rain.field[:] = inputs["mixing_ratio_rain"]
        state.mixing_ratio.snow.field[:] = inputs["mixing_ratio_snow"]
        state.mixing_ratio.graupel.field[:] = inputs["mixing_ratio_graupel"]
        state.mixing_ratio.large_scale_liquid.field[:] = inputs["mixing_ratio_large_scale_liquid"]
        state.mixing_ratio.large_scale_ice.field[:] = inputs["mixing_ratio_large_scale_ice"]
        state.mixing_ratio.convective_liquid.field[:] = inputs["mixing_ratio_convective_liquid"]
        state.mixing_ratio.convective_ice.field[:] = inputs["mixing_ratio_convective_ice"]
        state.cloud_fraction.large_scale.field[:] = inputs["cloud_fraction_large_scale"]
        state.cloud_fraction.convective.field[:] = inputs["cloud_fraction_convective"]
        state.concentration.liquid.field[:] = inputs["concentration_liquid"]
        state.concentration.ice.field[:] = inputs["concentration_ice"]
        state.area.field[:] = inputs["area"]
        state.p_interface.field[:] = inputs["p_interface"]
        state.z_interface.field[:] = inputs["z_interface"]
        state.t.field[:] = inputs["t"]
        state.u.field[:] = inputs["u"]
        state.v.field[:] = inputs["v"]
        state.land_fraction.field[:] = inputs["land_fraction"]
        state.covariance_liquid_water_static_energy_and_total_water_specific_humidity.field[:] = inputs[
            "covariance_liquid_water_static_energy_and_total_water_specific_humidity"
        ]
        state.omega.field[:] = inputs["omega"]
        state.pdf_first_plume_fractional_area.field[:] = inputs["pdf_first_plume_fractional_area"]
        state.vertical_motion.velocity.field[:] = inputs["vertical_motion_velocity"]
        state.vertical_motion.variance.field[:] = inputs["vertical_motion_variance"]
        state.vertical_motion.third_moment.field[:] = inputs["vertical_motion_third_moment"]
        state.liquid_water_static_energy.flux.field[:] = inputs["liquid_water_static_energy_flux"]
        state.liquid_water_static_energy.variance.field[:] = inputs["liquid_water_static_energy_variance"]
        state.liquid_water_static_energy.third_moment.field[:] = inputs["liquid_water_static_energy_third_moment"]
        state.total_water.flux.field[:] = inputs["total_water_flux"]
        state.total_water.variance.field[:] = inputs["total_water_variance"]
        state.total_water.third_moment.field[:] = inputs["total_water_third_moment"]
        state.lower_tropospheric_stability.field[:] = inputs["lower_tropospheric_stability"]
        state.estimated_inversion_strength.field[:] = inputs["estimated_inversion_strength"]
        state.tendencies.dcloud_fractiondt_macro.field[:] = inputs["tendencies_dcloud_fractiondt_macro"]
        state.tendencies.dvapordt_macro.field[:] = inputs["tendencies_dvapordt_macro"]
        state.tendencies.dicedt_macro.field[:] = inputs["tendencies_dicedt_macro"]
        state.tendencies.dliquiddt_macro.field[:] = inputs["tendencies_dliquiddt_macro"]
        state.tendencies.draindt_macro.field[:] = inputs["tendencies_draindt_macro"]
        state.tendencies.dgraupeldt_macro.field[:] = inputs["tendencies_dgraupeldt_macro"]
        state.tendencies.dsnowdt_macro.field[:] = inputs["tendencies_dsnowdt_macro"]
        state.tendencies.dudt_macro.field[:] = inputs["tendencies_dudt_macro"]
        state.tendencies.dvdt_macro.field[:] = inputs["tendencies_dvdt_macro"]
        state.tendencies.dtdt_macro.field[:] = inputs["tendencies_dtdt_macro"]
        state.convection_fraction.field[:] = inputs["convection_fraction"]
        state.surface_type.field[:] = inputs["surface_type"]
        state.cloud_liquid_evaporation.field[:] = inputs["cloud_liquid_evaporation"]
        state.cloud_ice_sublimation.field[:] = inputs["cloud_ice_sublimation"]
        state.icefall.field[:] = inputs["icefall"]
        state.freezing_rainfall.field[:] = inputs["freezing_rainfall"]
        state.relative_humidity_after_pdf.field[:] = inputs["relative_humidity_after_pdf"]
        state.buoyancy_flux.field[:] = inputs["buoyancy_flux"]
        state.liquid_water_flux.field[:] = inputs["liquid_water_flux"]
        state.hydrostatic_pdf_iterations.field[:] = inputs["hydrostatic_pdf_iterations"]
        state.radiation_field.cloud_fraction.field[:] = inputs["radiation_field_cloud_fraction"]
        state.radiation_field.vapor.field[:] = inputs["radiation_field_vapor"]
        state.radiation_field.liquid.field[:] = inputs["radiation_field_liquid"]
        state.radiation_field.ice.field[:] = inputs["radiation_field_ice"]
        state.radiation_field.rain.field[:] = inputs["radiation_field_rain"]
        state.radiation_field.snow.field[:] = inputs["radiation_field_snow"]
        state.radiation_field.graupel.field[:] = inputs["radiation_field_graupel"]
        state.cloud_particle_effective_radius.liquid.field[:] = inputs["cloud_particle_effective_radius_liquid"]
        state.cloud_particle_effective_radius.ice.field[:] = inputs["cloud_particle_effective_radius_ice"]
        state.precipitation_at_surface.rain.field[:] = inputs["precipitation_at_surface_rain"]
        state.precipitation_at_surface.snow.field[:] = inputs["precipitation_at_surface_snow"]
        state.precipitation_at_surface.ice.field[:] = inputs["precipitation_at_surface_ice"]
        state.precipitation_at_surface.graupel.field[:] = inputs["precipitation_at_surface_graupel"]
        state.non_anvil_large_scale.precip.field[:] = inputs["non_anvil_large_scale_precip"]
        state.non_anvil_large_scale.snow.field[:] = inputs["non_anvil_large_scale_snow"]
        state.non_anvil_large_scale.evaporation.field[:] = inputs["non_anvil_large_scale_evaporation"]
        state.non_anvil_large_scale.sublimation.field[:] = inputs["non_anvil_large_scale_sublimation"]
        state.non_anvil_large_scale.liquid_precip_flux.field[:] = inputs["non_anvil_large_scale_liquid_precip_flux"]
        state.non_anvil_large_scale.ice_precip_flux.field[:] = inputs["non_anvil_large_scale_ice_precip_flux"]
        state.anvil.liquid_precip_flux.field[:] = inputs["anvil_liquid_precip_flux"]
        state.anvil.ice_precip_flux.field[:] = inputs["anvil_ice_precip_flux"]
        state.critical_relative_humidity_for_pdf.field[:] = inputs["critical_relative_humidity_for_pdf"]
        state.tendencies.dcloud_fractiondt_micro.field[:] = inputs["tendencies_dcloud_fractiondt_micro"]
        state.tendencies.dvapordt_micro.field[:] = inputs["tendencies_dvapordt_micro"]
        state.tendencies.dicedt_micro.field[:] = inputs["tendencies_dicedt_micro"]
        state.tendencies.dliquiddt_micro.field[:] = inputs["tendencies_dliquiddt_micro"]
        state.tendencies.draindt_micro.field[:] = inputs["tendencies_draindt_micro"]
        state.tendencies.dgraupeldt_micro.field[:] = inputs["tendencies_dgraupeldt_micro"]
        state.tendencies.dsnowdt_micro.field[:] = inputs["tendencies_dsnowdt_micro"]
        state.tendencies.dudt_micro.field[:] = inputs["tendencies_dudt_micro"]
        state.tendencies.dvdt_micro.field[:] = inputs["tendencies_dvdt_micro"]
        state.tendencies.dtdt_micro.field[:] = inputs["tendencies_dtdt_micro"]
        # state.radar.simulated_reflectivity.field[:] = inputs["radar_simulated_reflectivity"]
        # state.radar.maximum_composite_reflectivity.field[:] = inputs["radar_maximum_composite_reflectivity"]
        # state.radar.base_1km_agl_reflectivity.field[:] = inputs["radar_base_1km_agl_reflectivity"]
        # state.radar.echo_top_reflectivity.field[:] = inputs["radar_echo_top_reflectivity"]
        # state.radar.minus_10c_reflectivity.field[:] = inputs["radar_minus_10c_reflectivity"]

        # initialize test class
        code = GFDL1M(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
        )

        # execute test code
        code(state=state)

        return {
            "precipitation_at_surface_deep_convective_precipitation": state.precipitation_at_surface.deep_convective_precipitation.field[:],  # noqa
            "precipitation_at_surface_anvil_precipitation": state.precipitation_at_surface.anvil_precipitation.field[:],  # noqa
            "precipitation_at_surface_shallow_convective_precipitation": state.precipitation_at_surface.shallow_convective_precipitation.field[:],  # noqa
            "precipitation_at_surface_deep_convective_snow": state.precipitation_at_surface.deep_convective_snow.field[:],  # noqa
            "precipitation_at_surface_anvil_snow": state.precipitation_at_surface.anvil_snow.field[:],
            "precipitation_at_surface_shallow_convective_snow": state.precipitation_at_surface.shallow_convective_snow.field[:],  # noqa
            # "lcl_height": state.lcl_height.field[:],
            "shallow_convection_rain": state.shallow_convection_rain.field[:],
            "shallow_convection_snow": state.shallow_convection_snow.field[:],
            # "large_scale_rainwater_source": state.large_scale_rainwater_source.field[:],
            "tendencies_dtdt_friction_pressure_weighted": state.tendencies.dtdt_friction_pressure_weighted.field[:],  # noqa
            "mixing_ratio_vapor": state.mixing_ratio.vapor.field[:],
            "mixing_ratio_rain": state.mixing_ratio.rain.field[:],
            "mixing_ratio_snow": state.mixing_ratio.snow.field[:],
            "mixing_ratio_graupel": state.mixing_ratio.graupel.field[:],
            "mixing_ratio_large_scale_liquid": state.mixing_ratio.large_scale_liquid.field[:],
            "mixing_ratio_large_scale_ice": state.mixing_ratio.large_scale_ice.field[:],
            "mixing_ratio_convective_liquid": state.mixing_ratio.convective_liquid.field[:],
            "mixing_ratio_convective_ice": state.mixing_ratio.convective_ice.field[:],
            "cloud_fraction_large_scale": state.cloud_fraction.large_scale.field[:],
            "cloud_fraction_convective": state.cloud_fraction.convective.field[:],
            "concentration_liquid": state.concentration.liquid.field[:],
            "concentration_ice": state.concentration.ice.field[:],
            "area": state.area.field[:],
            "p_interface": state.p_interface.field[:],
            "z_interface": state.z_interface.field[:],
            "t": state.t.field[:],
            "u": state.u.field[:],
            "v": state.v.field[:],
            "land_fraction": state.land_fraction.field[:],
            "covariance_liquid_water_static_energy_and_total_water_specific_humidity": state.covariance_liquid_water_static_energy_and_total_water_specific_humidity.field[  # noqa
                :
            ],
            "omega": state.omega.field[:],
            "pdf_first_plume_fractional_area": state.pdf_first_plume_fractional_area.field[:],
            "vertical_motion_velocity": state.vertical_motion.velocity.field[:],
            "vertical_motion_variance": state.vertical_motion.variance.field[:],
            "vertical_motion_third_moment": state.vertical_motion.third_moment.field[:],
            "liquid_water_static_energy_flux": state.liquid_water_static_energy.flux.field[:],
            "liquid_water_static_energy_variance": state.liquid_water_static_energy.variance.field[:],
            "liquid_water_static_energy_third_moment": state.liquid_water_static_energy.third_moment.field[:],
            "total_water_flux": state.total_water.flux.field[:],
            "total_water_variance": state.total_water.variance.field[:],
            "total_water_third_moment": state.total_water.third_moment.field[:],
            "lower_tropospheric_stability": state.lower_tropospheric_stability.field[:],
            "estimated_inversion_strength": state.estimated_inversion_strength.field[:],
            "tendencies_dcloud_fractiondt_macro": state.tendencies.dcloud_fractiondt_macro.field[:],
            "tendencies_dvapordt_macro": state.tendencies.dvapordt_macro.field[:],
            "tendencies_dicedt_macro": state.tendencies.dicedt_macro.field[:],
            "tendencies_dliquiddt_macro": state.tendencies.dliquiddt_macro.field[:],
            "tendencies_draindt_macro": state.tendencies.draindt_macro.field[:],
            "tendencies_dgraupeldt_macro": state.tendencies.dgraupeldt_macro.field[:],
            "tendencies_dsnowdt_macro": state.tendencies.dsnowdt_macro.field[:],
            "tendencies_dudt_macro": state.tendencies.dudt_macro.field[:],
            "tendencies_dvdt_macro": state.tendencies.dvdt_macro.field[:],
            "tendencies_dtdt_macro": state.tendencies.dtdt_macro.field[:],
            "convection_fraction": state.convection_fraction.field[:],
            "surface_type": state.surface_type.field[:],
            "cloud_liquid_evaporation": state.cloud_liquid_evaporation.field[:],
            "cloud_ice_sublimation": state.cloud_ice_sublimation.field[:],
            "icefall": state.icefall.field[:],
            "freezing_rainfall": state.freezing_rainfall.field[:],
            "relative_humidity_after_pdf": state.relative_humidity_after_pdf.field[:],
            "buoyancy_flux": state.buoyancy_flux.field[:],
            "liquid_water_flux": state.liquid_water_flux.field[:],
            "hydrostatic_pdf_iterations": state.hydrostatic_pdf_iterations.field[:],
            "radiation_field_cloud_fraction": state.radiation_field.cloud_fraction.field[:],
            "radiation_field_vapor": state.radiation_field.vapor.field[:],
            "radiation_field_liquid": state.radiation_field.liquid.field[:],
            "radiation_field_ice": state.radiation_field.ice.field[:],
            "radiation_field_rain": state.radiation_field.rain.field[:],
            "radiation_field_snow": state.radiation_field.snow.field[:],
            "radiation_field_graupel": state.radiation_field.graupel.field[:],
            "cloud_particle_effective_radius_liquid": state.cloud_particle_effective_radius.liquid.field[:],
            "cloud_particle_effective_radius_ice": state.cloud_particle_effective_radius.ice.field[:],
            "precipitation_at_surface_rain": state.precipitation_at_surface.rain.field[:],
            "precipitation_at_surface_snow": state.precipitation_at_surface.snow.field[:],
            "precipitation_at_surface_ice": state.precipitation_at_surface.ice.field[:],
            "precipitation_at_surface_graupel": state.precipitation_at_surface.graupel.field[:],
            "non_anvil_large_scale_precip": state.non_anvil_large_scale.precip.field[:],
            "non_anvil_large_scale_snow": state.non_anvil_large_scale.snow.field[:],
            "non_anvil_large_scale_evaporation": state.non_anvil_large_scale.evaporation.field[:],
            "non_anvil_large_scale_sublimation": state.non_anvil_large_scale.sublimation.field[:],
            "non_anvil_large_scale_liquid_precip_flux": state.non_anvil_large_scale.liquid_precip_flux.field[:],
            "non_anvil_large_scale_ice_precip_flux": state.non_anvil_large_scale.ice_precip_flux.field[:],
            "anvil_liquid_precip_flux": state.anvil.liquid_precip_flux.field[:],
            "anvil_ice_precip_flux": state.anvil.ice_precip_flux.field[:],
            "critical_relative_humidity_for_pdf": state.critical_relative_humidity_for_pdf.field[:],
            "tendencies_dcloud_fractiondt_micro": state.tendencies.dcloud_fractiondt_micro.field[:],
            "tendencies_dvapordt_micro": state.tendencies.dvapordt_micro.field[:],
            "tendencies_dicedt_micro": state.tendencies.dicedt_micro.field[:],
            "tendencies_dliquiddt_micro": state.tendencies.dliquiddt_micro.field[:],
            "tendencies_draindt_micro": state.tendencies.draindt_micro.field[:],
            "tendencies_dgraupeldt_micro": state.tendencies.dgraupeldt_micro.field[:],
            "tendencies_dsnowdt_micro": state.tendencies.dsnowdt_micro.field[:],
            "tendencies_dudt_micro": state.tendencies.dudt_micro.field[:],
            "tendencies_dvdt_micro": state.tendencies.dvdt_micro.field[:],
            "tendencies_dtdt_micro": state.tendencies.dtdt_micro.field[:],
            # inputs["radar_simulated_reflectivity"]: state.radar.simulated_reflectivity.field[:],
            # inputs["radar_maximum_composite_reflectivity"]: state.radar.maximum_composite_reflectivity.field[:], # noqa
            # inputs["radar_base_1km_agl_reflectivity"]: state.radar.base_1km_agl_reflectivity.field[:],
            # inputs["radar_echo_top_reflectivity"]: state.radar.echo_top_reflectivity.field[:],
            # inputs["radar_minus_10c_reflectivity"]: state.radar.minus_10c_reflectivity.field[:],
        }
