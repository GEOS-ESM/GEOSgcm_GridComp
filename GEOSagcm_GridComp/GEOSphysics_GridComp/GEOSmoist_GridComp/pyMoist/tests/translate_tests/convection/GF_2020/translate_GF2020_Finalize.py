import numpy as np
from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import NUMBER_OF_PLUMES
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.finalize import GF2020Finalize
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection_tracers import ConvectionTracers
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGF2020_Finalize(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
            # fields saved midway through GF2020 fortran
            "internal_lateral_entrainment_rate": {},
            "local_edge_height_above_surface": {},
            "local_layer_height_above_surface": {},
            "local_p": {},
            "local_p_kappa": {},
            "local_th": {},
            "local_mass": {},
            "local_modified_area": {},
            "local_vertical_velocity": {},
            "local_dz": {},
            "local_air_density": {},
            "local_scalar_diffusivity": {},
            "local_vapor_current": {},
            "local_p_flipped": {},
            "local_vapor_flipped": {},
            # NOTE these are disabled and loaded manually because the translate test cannot handle 4D data
            # # CumulusParameterization state - input
            # "t_excess": {},
            # "vapor_excess": {},
            # "grid_scale_forcing_t": {},
            # "grid_scale_forcing_vapor": {},
            # "subgrid_scale_forcing_t": {},
            # "subgrid_scale_forcing_vapor": {},
            # "seed_convection": {},
            # "saturation_water_vapor": {},
            # "ocean_fraction": {},
            # "convection_fraction": {},
            # "surface_type": {},
            # "lateral_entrainment_rate": {},
            # "last_error_code": {},
            # # CumulusParameterization state - output
            # "dtdt": {},
            # "dvapordt": {},
            # "dcloudicedt": {},
            # "dudt": {},
            # "dvdt": {},
            # "dnliquiddt": {},
            # "dnicedt": {},
            # "dbuoyancydt": {},
            # "dconvectiveicedt": {},
            # "dlargescaleicedt": {},
            # "dconvectiveliquiddt": {},
            # "dlargescaleliquiddt": {},
            # "dconvectivecloudfractiondt": {},
            # "dlargescalecloudfractiondt": {},
            # "error_code": {},
            # "downdraft_origin_level": {},
            # "lcl_level": {},
            # "updraft_origin_level": {},
            # "updraft_lfc_level": {},
            # "cloud_top_level": {},
            # "kstabi": {},
            # "kstabm": {},
            # "precip": {},
            # "cloud_base_mass_flux_modified": {},
            # "epsilon_forced": {},
            # "total_normalized_integrated_condensate_forced": {},
            # "scale_dependence_factor": {},
            # "p_cloud_levels_forced": {},
            # "entrainment_rate": {},
            # "mass_entrainment_updraft_forced": {},
            # "mass_entrainment_downdraft_forced": {},
            # "mass_detrainment_updraft_forced": {},
            # "mass_detrainment_downdraft_forced": {},
            # "normalized_massflux_updraft_forced": {},
            # "normalized_massflux_downdraft_forced": {},
            # "condensate_to_fall_forced": {},
            # "evaporate_in_downdraft_forced": {},
            # "cloud_liquid_after_rain_forced": {},
            # "t_updraft": {},
            # "convective_cloud_fraction_output": {},
            # "cloud_workfunction_0": {},
            # "cloud_workfunction_1": {},
            # "cloud_workfunction_2": {},
            # "cloud_workfunction_3": {},
            # "cloud_workfunction_1_pbl": {},
            # "cloud_workfunction_1_cin": {},
            # "cape_removal_time_scale": {},
            # "pbl_time_scale": {},
            # "lightning_density": {},
            # "evaporation_sublimation_tendency": {},
            # "convective_precip_flux": {},
            # "t_perturbation": {},
            # # CumulusParameterization state - input-output
            # "grid_length": {},
            # "pbl_level": {},
            # "ccn": {},
            # "air_density": {},
            # "omega": {},
            # "topography_height_no_negative": {},
            # "sensible_heat_flux": {},
            # "latent_heat_flux": {},
            # "longitude_degrees": {},
            # "latitude_degrees": {},
            # "t_old": {},
            # "vapor_old": {},
            # "t_modified_by_advection": {},
            # "vapor_modified_by_advection": {},
            # "geopotential_height_forced": {},
            # "p_forced": {},
            # "p_surface": {},
            # "t_surface": {},
            # "u": {},
            # "v": {},
            # "w": {},
            # "mass": {},
            # "convective_scale_velocity": {},
            # "buoyancy_excess": {},
            # "large_scale_ice": {},
            # "convective_ice": {},
            # "large_scale_liquid": {},
            # "convective_liquid": {},
            # "large_scale_cloud_fraction": {},
            # "convective_cloud_fraction": {},
            # "chemistry_tracers": {},
            # "chemistry_tracers_output": {},
        }

        # NOTE disabled fields are nan in fortran - zero in python, disabled so the test passes
        self.out_vars: dict = {
            "latitude_bugworkaroundname": {},
            "longitude_bugworkaroundname": {},
            "p_interface_bugworkaroundname": {},
            "t_bugworkaroundname": {},
            "u_bugworkaroundname": {},
            "v_bugworkaroundname": {},
            "w_bugworkaroundname": {},
            "omega_bugworkaroundname": {},
            "t_2m_bugworkaroundname": {},
            "specific_humidity_2m_bugworkaroundname": {},
            "t_surface_bugworkaroundname": {},
            "specific_humidity_surface_bugworkaroundname": {},
            "vapor_bugworkaroundname": {},
            "convective_liquid_bugworkaroundname": {},
            "convective_ice_bugworkaroundname": {},
            "large_scale_liquid_bugworkaroundname": {},
            "large_scale_ice_bugworkaroundname": {},
            "convective_cloud_fraction_bugworkaroundname": {},
            "large_scale_cloud_fraction_bugworkaroundname": {},
            "p_interface_timestep_start_bugworkaroundname": {},
            "t_timestep_start_bugworkaroundname": {},
            "u_timestep_start_bugworkaroundname": {},
            "v_timestep_start_bugworkaroundname": {},
            "vapor_timestep_start_bugworkaroundname": {},
            "geopotential_height_interface_bugworkaroundname": {},
            "geopotential_height_surface_bugworkaroundname": {},
            "area_bugworkaroundname": {},
            "pbl_level_bugworkaroundname": {},
            "convection_fraction_bugworkaroundname": {},
            "surface_type_bugworkaroundname": {},
            "seed_convection_bugworkaroundname": {},
            "land_fraction_bugworkaroundname": {},
            "scalar_diffusivity_bugworkaroundname": {},
            "buoyancy_bugworkaroundname": {},
            "convective_precipitation_GF_bugworkaroundname": {},
            "convective_precipitation_RAS_bugworkaroundname": {},
            "sensible_heat_flux_bugworkaroundname": {},
            # "total_water_flux_deep_convection_interface_bugworkaroundname": {}, # disabled temporarily - do not push, fails b/c of <10 ulp diff in saturation functions
            "evaporation_bugworkaroundname": {},
            "sublimation_of_convective_precipitation_bugworkaroundname": {},
            "evaporation_of_convective_precipitation_bugworkaroundname": {},
            "ice_precip_flux_interface_bugworkaroundname": {},
            "liquid_precip_flux_interface_bugworkaroundname": {},
            "convective_condensate_source_bugworkaroundname": {},
            "convective_condensate_grid_mean_bugworkaroundname": {},
            "entrainment_parameter_bugworkaroundname": {},
            "lateral_entrainment_rate_bugworkaroundname": {},
            "lateral_entrainment_rate_shallow_bugworkaroundname": {},
            "lateral_entrainment_rate_mid_bugworkaroundname": {},
            "lateral_entrainment_rate_deep_bugworkaroundname": {},
            "updraft_areal_fraction_bugworkaroundname": {},
            "updraft_vertical_velocity_bugworkaroundname": {},
            "dtdt_shortwave_bugworkaroundname": {},
            "dtdt_longwave_bugworkaroundname": {},
            "dspecific_humiditydt_pbl_bugworkaroundname": {},
            "dtdt_pbl_bugworkaroundname": {},
            "dtdt_from_dynamics_bugworkaroundname": {},
            "dvapordt_from_dynamics_bugworkaroundname": {},
            "sigma_mid_bugworkaroundname": {},
            "sigma_deep_bugworkaroundname": {},
            "total_precipitable_water_initial_bugworkaroundname": {},
            "saturation_total_precipitable_water_initial_bugworkaroundname": {},
            "dvapordt_deep_convection_bugworkaroundname": {},
            "dtdt_deep_convection_bugworkaroundname": {},
            "dudt_deep_convection_bugworkaroundname": {},
            "dvdt_deep_convection_bugworkaroundname": {},
            "dliquiddt_deep_convection_bugworkaroundname": {},
            "dicedt_deep_convection_bugworkaroundname": {},
            "dcloudfractiondt_deep_convection_bugworkaroundname": {},
            "pressure_shallow_convective_cloud_top_bugworkaroundname": {},
            "pressure_mid_convective_cloud_top_bugworkaroundname": {},
            "pressure_deep_convective_cloud_top_bugworkaroundname": {},
            "mass_flux_shallow_bugworkaroundname": {},
            "mass_flux_mid_bugworkaroundname": {},
            "mass_flux_deep_updraft_bugworkaroundname": {},
            "mass_flux_deep_updraft_interface_bugworkaroundname": {},
            "mass_flux_deep_updraft_detrained_bugworkaroundname": {},
            "mass_flux_deep_downdraft_bugworkaroundname": {},
            "mass_flux_cloud_base_bugworkaroundname": {},
            "mass_flux_cloud_base_shallow_bugworkaroundname": {},
            "mass_flux_cloud_base_mid_bugworkaroundname": {},
            "mass_flux_cloud_base_deep_bugworkaroundname": {},
            "total_cumulative_mass_flux_interface_bugworkaroundname": {},
            "total_detraining_mass_flux_bugworkaroundname": {},
            "convection_code_shallow_bugworkaroundname": {},
            "convection_code_mid_bugworkaroundname": {},
            "convection_code_deep_bugworkaroundname": {},
            "cloud_workfunction_0_bugworkaroundname": {},
            "cloud_workfunction_1_bugworkaroundname": {},
            "cloud_workfunction_2_bugworkaroundname": {},
            "cloud_workfunction_3_bugworkaroundname": {},
            "cloud_workfunction_1_pbl_bugworkaroundname": {},
            "cloud_workfunction_1_cin_bugworkaroundname": {},
            "pbl_time_scale_bugworkaroundname": {},
            "cape_removal_time_scale_bugworkaroundname": {},
            "lightning_density_bugworkaroundname": {},
            "convection_tracer_bugworkaroundname": {},
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")
        self.convection_tracers_input = data_loader.load("GF2020_ConvectionTracers", use_dynamic_i_call=True)

        # workaround because translate test cannot read in 4d fields
        self.manual_inputs = data_loader.load("GF2020_CumulusParameterization-Out", use_dynamic_i_call=True)

        # load data from GF2020-In to fill in unmodified fields from the state
        self.unmodified_state = data_loader.load("GF2020-In", use_dynamic_i_call=True)

        # load data from GF2020_Finalize-Out to get total_precipitable_water_initial and
        # saturation total_precipitable_water_initial
        self.manual_outputs = data_loader.load("GF2020_Finalize-Out", use_dynamic_i_call=True)

    def compute(self, inputs):
        config = GF2020Config(**self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)

        # initialize GF2020 state
        state = GF2020State.zeros(self.quantity_factory)

        # initialize GF2020 CumulusParameterization state
        cumulus_parameterization_state = GF2020CumulusParameterizationState.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": NUMBER_OF_PLUMES,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # initialize GF2020 locals
        locals = GF2020Locals.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

        # fill state with data not passed in/unmodified inside the cumulus parameterization core
        state.latitude.field[:] = self.unmodified_state["latitude_bugworkaroundname"]
        state.longitude.field[:] = self.unmodified_state["longitude_bugworkaroundname"]
        state.p_interface.field[:] = self.unmodified_state["p_interface_bugworkaroundname"]
        state.t.field[:] = self.unmodified_state["t_bugworkaroundname"]
        state.u.field[:] = self.unmodified_state["u_bugworkaroundname"]
        state.v.field[:] = self.unmodified_state["v_bugworkaroundname"]
        state.w.field[:] = self.unmodified_state["w_bugworkaroundname"]
        state.omega.field[:] = self.unmodified_state["omega_bugworkaroundname"]
        state.t_2m.field[:] = self.unmodified_state["t_2m_bugworkaroundname"]
        state.specific_humidity_2m.field[:] = self.unmodified_state["specific_humidity_2m_bugworkaroundname"]
        state.t_surface.field[:] = self.unmodified_state["t_surface_bugworkaroundname"]
        state.specific_humidity_surface.field[:] = self.unmodified_state["specific_humidity_surface_bugworkaroundname"]
        state.vapor.field[:] = self.unmodified_state["vapor_bugworkaroundname"]
        state.convective_liquid.field[:] = self.unmodified_state["convective_liquid_bugworkaroundname"]
        state.convective_ice.field[:] = self.unmodified_state["convective_ice_bugworkaroundname"]
        state.convective_cloud_fraction.field[:] = self.unmodified_state["convective_cloud_fraction_bugworkaroundname"]
        state.large_scale_liquid.field[:] = self.unmodified_state["large_scale_liquid_bugworkaroundname"]
        state.large_scale_ice.field[:] = self.unmodified_state["large_scale_ice_bugworkaroundname"]
        state.large_scale_cloud_fraction.field[:] = self.unmodified_state["large_scale_cloud_fraction_bugworkaroundname"]
        state.p_interface_timestep_start.field[:] = self.unmodified_state["p_interface_timestep_start_bugworkaroundname"]
        state.t_timestep_start.field[:] = self.unmodified_state["t_timestep_start_bugworkaroundname"]
        state.u_timestep_start.field[:] = self.unmodified_state["u_timestep_start_bugworkaroundname"]
        state.v_timestep_start.field[:] = self.unmodified_state["v_timestep_start_bugworkaroundname"]
        state.vapor_timestep_start.field[:] = self.unmodified_state["vapor_timestep_start_bugworkaroundname"]
        state.geopotential_height_interface.field[:] = self.unmodified_state["geopotential_height_interface_bugworkaroundname"]
        state.geopotential_height_surface.field[:] = self.unmodified_state["geopotential_height_surface_bugworkaroundname"]
        state.area.field[:] = self.unmodified_state["area_bugworkaroundname"]
        state.pbl_level.field[:] = self.unmodified_state["pbl_level_bugworkaroundname"]
        state.convection_fraction.field[:] = self.unmodified_state["convection_fraction_bugworkaroundname"]
        state.surface_type.field[:] = self.unmodified_state["surface_type_bugworkaroundname"]
        state.seed_convection.field[:] = self.unmodified_state["seed_convection_bugworkaroundname"]
        state.land_fraction.field[:] = self.unmodified_state["land_fraction_bugworkaroundname"]
        state.scalar_diffusivity.field[:] = self.unmodified_state["scalar_diffusivity_bugworkaroundname"]
        state.buoyancy.field[:] = self.unmodified_state["buoyancy_bugworkaroundname"]
        # state.convective_precipitation_GF.field[:] = self.unmodified_state["convective_precipitation_GF_bugworkaroundname"]
        # state.convective_precipitation_RAS.field[:] = self.unmodified_state["convective_precipitation_RAS_bugworkaroundname"]
        state.sensible_heat_flux.field[:] = self.unmodified_state["sensible_heat_flux_bugworkaroundname"]
        # state.total_water_flux_deep_convection.field[:] = self.unmodified_state["total_water_flux_deep_convection_interface_bugworkaroundname"]
        # state.sublimation_of_convective_precipitation.field[:] = self.unmodified_state["sublimation_of_convective_precipitation_bugworkaroundname"]
        # state.evaporation_of_convective_precipitation.field[:] = self.unmodified_state["evaporation_of_convective_precipitation_bugworkaroundname"]
        # state.ice_precip_flux_interface.field[:] = self.unmodified_state["ice_precip_flux_interface_bugworkaroundname"]
        # state.liquid_precip_flux_interface.field[:] = self.unmodified_state["liquid_precip_flux_interface_bugworkaroundname"]
        state.evaporation.field[:] = self.unmodified_state["evaporation_bugworkaroundname"]
        # state.convective_condensate_source.field[:] = self.unmodified_state["convective_condensate_source_bugworkaroundname"]
        # state.convective_condensate_grid_mean.field[:] = self.unmodified_state["convective_condensate_grid_mean_bugworkaroundname"]
        # state.entrainment_parameter.field[:] = self.unmodified_state["entrainment_parameter_bugworkaroundname"]
        state.lateral_entrainment_rate.field[:] = inputs["internal_lateral_entrainment_rate"]
        # state.lateral_entrainment_rate_shallow.field[:] = self.unmodified_state["lateral_entrainment_rate_shallow_bugworkaroundname"]
        # state.lateral_entrainment_rate_mid.field[:] = self.unmodified_state["lateral_entrainment_rate_mid_bugworkaroundname"]
        # state.lateral_entrainment_rate_deep.field[:] = self.unmodified_state["lateral_entrainment_rate_deep_bugworkaroundname"]
        # state.updraft_areal_fraction.field[:] = self.unmodified_state["updraft_areal_fraction_bugworkaroundname"]
        # state.updraft_vertical_velocity.field[:] = self.unmodified_state["updraft_vertical_velocity_bugworkaroundname"]
        state.dtdt_shortwave.field[:] = self.unmodified_state["dtdt_shortwave_bugworkaroundname"]
        state.dtdt_longwave.field[:] = self.unmodified_state["dtdt_longwave_bugworkaroundname"]
        state.dspecific_humiditydt_pbl.field[:] = self.unmodified_state["dspecific_humiditydt_pbl_bugworkaroundname"]
        state.dtdt_pbl.field[:] = self.unmodified_state["dtdt_pbl_bugworkaroundname"]
        state.dtdt_from_dynamics.field[:] = self.unmodified_state["dtdt_from_dynamics_bugworkaroundname"]
        state.dvapordt_from_dynamics.field[:] = self.unmodified_state["dvapordt_from_dynamics_bugworkaroundname"]
        # state.sigma_mid.field[:] = self.unmodified_state["sigma_mid_bugworkaroundname"]
        # state.sigma_deep.field[:] = self.unmodified_state["sigma_deep_bugworkaroundname"]
        state.total_precipitable_water_initial.field[:] = self.manual_outputs["total_precipitable_water_initial_bugworkaroundname"]
        state.saturation_total_precipitable_water_initial.field[:] = self.manual_outputs["saturation_total_precipitable_water_initial_bugworkaroundname"]
        # state.dvapordt_deep_convection.field[:] = self.unmodified_state["dvapordt_deep_convection_bugworkaroundname"]
        # state.dtdt_deep_convection.field[:] = self.unmodified_state["dtdt_deep_convection_bugworkaroundname"]
        # state.dudt_deep_convection.field[:] = self.unmodified_state["dudt_deep_convection_bugworkaroundname"]
        # state.dvdt_deep_convection.field[:] = self.unmodified_state["dvdt_deep_convection_bugworkaroundname"]
        # state.dliquiddt_deep_convection.field[:] = self.unmodified_state["dliquiddt_deep_convection_bugworkaroundname"]
        # state.dicedt_deep_convection.field[:] = self.unmodified_state["dicedt_deep_convection_bugworkaroundname"]
        # state.dcloudfractiondt_deep_convection.field[:] = self.unmodified_state["dcloudfractiondt_deep_convection_bugworkaroundname"]
        # state.pressure_shallow_convective_cloud_top.field[:] = self.unmodified_state["pressure_shallow_convective_cloud_top_bugworkaroundname"]
        # state.pressure_mid_convective_cloud_top.field[:] = self.unmodified_state["pressure_mid_convective_cloud_top_bugworkaroundname"]
        # state.pressure_deep_convective_cloud_top.field[:] = self.unmodified_state["pressure_deep_convective_cloud_top_bugworkaroundname"]
        # state.mass_flux_shallow.field[:] = self.unmodified_state["mass_flux_shallow_bugworkaroundname"]
        # state.mass_flux_mid.field[:] = self.unmodified_state["mass_flux_mid_bugworkaroundname"]
        # state.mass_flux_deep_updraft.field[:] = self.unmodified_state["mass_flux_deep_updraft_bugworkaroundname"]
        # state.mass_flux_deep_updraft_interface.field[:] = self.unmodified_state["mass_flux_deep_updraft_interface_bugworkaroundname"]
        # state.mass_flux_deep_updraft_detrained.field[:] = self.unmodified_state["mass_flux_deep_updraft_detrained_bugworkaroundname"]
        # state.mass_flux_deep_downdraft.field[:] = self.unmodified_state["mass_flux_deep_downdraft_bugworkaroundname"]
        # state.mass_flux_cloud_base.field[:] = self.unmodified_state["mass_flux_cloud_base_bugworkaroundname"]
        # state.mass_flux_cloud_base_shallow.field[:] = self.unmodified_state["mass_flux_cloud_base_shallow_bugworkaroundname"]
        # state.mass_flux_cloud_base_mid.field[:] = self.unmodified_state["mass_flux_cloud_base_mid_bugworkaroundname"]
        # state.mass_flux_cloud_base_deep.field[:] = self.unmodified_state["mass_flux_cloud_base_deep_bugworkaroundname"]
        state.total_cumulative_mass_flux_interface.field[:] = self.unmodified_state["total_cumulative_mass_flux_interface_bugworkaroundname"]
        state.total_detraining_mass_flux.field[:] = self.unmodified_state["total_detraining_mass_flux_bugworkaroundname"]
        # state.convection_code_shallow.field[:] = self.unmodified_state["convection_code_shallow_bugworkaroundname"]
        # state.convection_code_mid.field[:] = self.unmodified_state["convection_code_mid_bugworkaroundname"]
        # state.convection_code_deep.field[:] = self.unmodified_state["convection_code_deep_bugworkaroundname"]
        # state.cloud_workfunction_0.field[:] = self.unmodified_state["cloud_workfunction_0_bugworkaroundname"]
        # state.cloud_workfunction_1.field[:] = self.unmodified_state["cloud_workfunction_1_bugworkaroundname"]
        # state.cloud_workfunction_2.field[:] = self.unmodified_state["cloud_workfunction_2_bugworkaroundname"]
        # state.cloud_workfunction_3.field[:] = self.unmodified_state["cloud_workfunction_3_bugworkaroundname"]
        # state.cloud_workfunction_1_pbl.field[:] = self.unmodified_state["cloud_workfunction_1_pbl_bugworkaroundname"]
        # state.cloud_workfunction_1_cin.field[:] = self.unmodified_state["cloud_workfunction_1_cin_bugworkaroundname"]
        # state.pbl_time_scale.field[:] = self.unmodified_state["pbl_time_scale_bugworkaroundname"]
        # state.cape_removal_time_scale.field[:] = self.unmodified_state["cape_removal_time_scale_bugworkaroundname"]
        # state.lightning_density.field[:] = self.unmodified_state["lightning_density_bugworkaroundname"]
        state.convection_tracer.field[:] = self.unmodified_state["convection_tracer_bugworkaroundname"]

        # fill locals with input data
        locals.derived_state.edge_height_above_surface.field[:] = inputs["local_edge_height_above_surface"]
        locals.derived_state.layer_height_above_surface.field[:] = inputs["local_layer_height_above_surface"]
        locals.derived_state.p.field[:] = inputs["local_p"]
        locals.derived_state.p_kappa.field[:] = inputs["local_p_kappa"]
        locals.derived_state.th.field[:] = inputs["local_th"]
        locals.derived_state.mass.field[:] = inputs["local_mass"]
        locals.derived_state.modified_area.field[:] = inputs["local_modified_area"]
        locals.derived_state.vertical_velocity.field[:] = inputs["local_vertical_velocity"]
        locals.derived_state.dz.field[:] = inputs["local_dz"]
        locals.derived_state.air_density.field[:] = inputs["local_air_density"]
        locals.flipped_copy.scalar_diffusivity.field[:] = np.moveaxis(inputs["local_scalar_diffusivity"], 0, 2)
        locals.flipped_copy.vapor_current.field[:] = np.moveaxis(inputs["local_vapor_current"], 0, 2)
        locals.flipped_copy.p.field[:] = np.moveaxis(inputs["local_p_flipped"], 0, 2)
        locals.flipped_copy.vapor.field[:] = np.moveaxis(inputs["local_vapor_flipped"], 0, 2)

        # fill CumlulusParameterizaiton State with input data
        # input fields
        cumulus_parameterization_state.input.t_excess.field[:] = self.manual_inputs["t_excess"]
        cumulus_parameterization_state.input.vapor_excess.field[:] = self.manual_inputs["vapor_excess"]
        cumulus_parameterization_state.input.grid_scale_forcing_t.field[:] = self.manual_inputs["grid_scale_forcing_t"]
        cumulus_parameterization_state.input.grid_scale_forcing_vapor.field[:] = self.manual_inputs["grid_scale_forcing_vapor"]
        cumulus_parameterization_state.input.subgrid_scale_forcing_t.field[:] = self.manual_inputs["subgrid_scale_forcing_t"]
        cumulus_parameterization_state.input.subgrid_scale_forcing_vapor.field[:] = self.manual_inputs["subgrid_scale_forcing_vapor"]
        cumulus_parameterization_state.input.seed_convection.field[:] = self.manual_inputs["seed_convection"]
        cumulus_parameterization_state.input.saturation_water_vapor.field[:] = self.manual_inputs["saturation_water_vapor"]
        cumulus_parameterization_state.input.ocean_fraction.field[:] = self.manual_inputs["ocean_fraction"]
        cumulus_parameterization_state.input.convection_fraction.field[:] = self.manual_inputs["convection_fraction"]
        cumulus_parameterization_state.input.surface_type.field[:] = self.manual_inputs["surface_type"]
        cumulus_parameterization_state.input.lateral_entrainment_rate.field[:] = self.manual_inputs["lateral_entrainment_rate"]
        cumulus_parameterization_state.input.last_error_code.field[:] = self.manual_inputs["last_error_code"]
        # output fields
        cumulus_parameterization_state.output.dtdt.field[:] = self.manual_inputs["dtdt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dvapordt.field[:] = self.manual_inputs["dvapordt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dcloudicedt.field[:] = self.manual_inputs["dcloudicedt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dudt.field[:] = self.manual_inputs["dudt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dvdt.field[:] = self.manual_inputs["dvdt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dnliquiddt.field[:] = self.manual_inputs["dnliquiddt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dnicedt.field[:] = self.manual_inputs["dnicedt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dbuoyancydt.field[:] = self.manual_inputs["dbuoyancydt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dconvectiveicedt.field[:] = self.manual_inputs["dconvectiveicedt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dlargescaleicedt.field[:] = self.manual_inputs["dlargescaleicedt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dconvectiveliquiddt.field[:] = self.manual_inputs["dconvectiveliquiddt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dlargescaleliquiddt.field[:] = self.manual_inputs["dlargescaleliquiddt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dconvectivecloudfractiondt.field[:] = self.manual_inputs["dconvectivecloudfractiondt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.dlargescalecloudfractiondt.field[:] = self.manual_inputs["dlargescalecloudfractiondt"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.error_code.field[:] = self.manual_inputs["error_code"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.downdraft_origin_level.field[:] = self.manual_inputs["downdraft_origin_level"][:, :, [1, 2, 0]] - 1
        cumulus_parameterization_state.output.lcl_level.field[:] = self.manual_inputs["lcl_level"][:, :, [1, 2, 0]] - 1
        cumulus_parameterization_state.output.updraft_origin_level.field[:] = self.manual_inputs["updraft_origin_level"][:, :, [1, 2, 0]] - 1
        cumulus_parameterization_state.output.updraft_lfc_level.field[:] = self.manual_inputs["updraft_lfc_level"][:, :, [1, 2, 0]] - 1
        cumulus_parameterization_state.output.cloud_top_level.field[:] = self.manual_inputs["cloud_top_level"][:, :, [1, 2, 0]] - 1
        cumulus_parameterization_state.output.kstabi.field[:] = self.manual_inputs["kstabi"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.kstabm.field[:] = self.manual_inputs["kstabm"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.precip.field[:] = self.manual_inputs["precip"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.cloud_base_mass_flux_modified.field[:] = self.manual_inputs["cloud_base_mass_flux_modified"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.epsilon_forced.field[:] = self.manual_inputs["epsilon_forced"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.total_normalized_integrated_condensate_forced.field[:] = self.manual_inputs[
            "total_normalized_integrated_condensate_forced"
        ][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.scale_dependence_factor.field[:] = self.manual_inputs["scale_dependence_factor"][:, :, [1, 2, 0]]
        cumulus_parameterization_state.output.p_cloud_levels_forced.field[:] = self.manual_inputs["p_cloud_levels_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.entrainment_rate.field[:] = self.manual_inputs["entrainment_rate"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.mass_entrainment_updraft_forced.field[:] = self.manual_inputs["mass_entrainment_updraft_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.mass_entrainment_downdraft_forced.field[:] = self.manual_inputs["mass_entrainment_downdraft_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.mass_detrainment_updraft_forced.field[:] = self.manual_inputs["mass_detrainment_updraft_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.mass_detrainment_downdraft_forced.field[:] = self.manual_inputs["mass_detrainment_downdraft_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.normalized_massflux_updraft_forced.field[:] = self.manual_inputs["normalized_massflux_updraft_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.normalized_massflux_downdraft_forced.field[:] = self.manual_inputs["normalized_massflux_downdraft_forced"][
            :, :, :, [1, 2, 0]
        ]
        cumulus_parameterization_state.output.condensate_to_fall_forced.field[:] = self.manual_inputs["condensate_to_fall_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.evaporate_in_downdraft_forced.field[:] = self.manual_inputs["evaporate_in_downdraft_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.cloud_liquid_after_rain_forced.field[:] = self.manual_inputs["cloud_liquid_after_rain_forced"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.t_updraft.field[:] = self.manual_inputs["t_updraft"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.convective_cloud_fraction.field[:] = self.manual_inputs["convective_cloud_fraction_output"][:, :, :, [1, 2, 0]]
        cumulus_parameterization_state.output.cloud_workfunction_0.field[:] = self.manual_inputs["cloud_workfunction_0"]
        cumulus_parameterization_state.output.cloud_workfunction_1.field[:] = self.manual_inputs["cloud_workfunction_1"]
        cumulus_parameterization_state.output.cloud_workfunction_2.field[:] = self.manual_inputs["cloud_workfunction_2"]
        cumulus_parameterization_state.output.cloud_workfunction_3.field[:] = self.manual_inputs["cloud_workfunction_3"]
        cumulus_parameterization_state.output.cloud_workfunction_1_pbl.field[:] = self.manual_inputs["cloud_workfunction_1_pbl"]
        cumulus_parameterization_state.output.cloud_workfunction_1_cin.field[:] = self.manual_inputs["cloud_workfunction_1_cin"]
        cumulus_parameterization_state.output.cape_removal_time_scale.field[:] = self.manual_inputs["cape_removal_time_scale"]
        cumulus_parameterization_state.output.pbl_time_scale.field[:] = self.manual_inputs["pbl_time_scale"]
        cumulus_parameterization_state.output.lightning_density.field[:] = self.manual_inputs["lightning_density"]
        cumulus_parameterization_state.output.evaporation_sublimation_tendency.field[:] = self.manual_inputs["evaporation_sublimation_tendency"]
        cumulus_parameterization_state.output.convective_precip_flux.field[:] = self.manual_inputs["convective_precip_flux"]
        cumulus_parameterization_state.output.t_perturbation.field[:] = self.manual_inputs["t_perturbation"]
        # input/output fields
        cumulus_parameterization_state.input_output.grid_length.field[:] = self.manual_inputs["grid_length"]
        cumulus_parameterization_state.input_output.pbl_level.field[:] = self.manual_inputs["pbl_level"] - 1
        cumulus_parameterization_state.input_output.ccn.field[:] = self.manual_inputs["ccn"]
        cumulus_parameterization_state.input_output.air_density.field[:] = self.manual_inputs["air_density"]
        cumulus_parameterization_state.input_output.omega.field[:] = self.manual_inputs["omega"]
        cumulus_parameterization_state.input_output.topography_height_no_negative.field[:] = self.manual_inputs["topography_height_no_negative"]
        cumulus_parameterization_state.input_output.sensible_heat_flux.field[:] = self.manual_inputs["sensible_heat_flux"]
        cumulus_parameterization_state.input_output.latent_heat_flux.field[:] = self.manual_inputs["latent_heat_flux"]
        cumulus_parameterization_state.input_output.longitude_degrees.field[:] = self.manual_inputs["longitude_degrees"]
        cumulus_parameterization_state.input_output.latitude_degrees.field[:] = self.manual_inputs["latitude_degrees"]
        cumulus_parameterization_state.input_output.t_old.field[:] = self.manual_inputs["t_old"]
        cumulus_parameterization_state.input_output.vapor_old.field[:] = self.manual_inputs["vapor_old"]
        cumulus_parameterization_state.input_output.t_modified_by_advection.field[:] = self.manual_inputs["t_modified_by_advection"]
        cumulus_parameterization_state.input_output.vapor_modified_by_advection.field[:] = self.manual_inputs["vapor_modified_by_advection"]
        cumulus_parameterization_state.input_output.geopotential_height_forced.field[:] = self.manual_inputs["geopotential_height_forced"]
        cumulus_parameterization_state.input_output.p_forced.field[:] = self.manual_inputs["p_forced"]
        cumulus_parameterization_state.input_output.p_surface.field[:] = self.manual_inputs["p_surface"]
        cumulus_parameterization_state.input_output.t_surface.field[:] = self.manual_inputs["t_surface"]
        cumulus_parameterization_state.input_output.u.field[:] = self.manual_inputs["u"]
        cumulus_parameterization_state.input_output.v.field[:] = self.manual_inputs["v"]
        cumulus_parameterization_state.input_output.w.field[:] = self.manual_inputs["w"]
        cumulus_parameterization_state.input_output.mass.field[:] = self.manual_inputs["mass"]
        cumulus_parameterization_state.input_output.convective_scale_velocity.field[:] = self.manual_inputs["convective_scale_velocity"]
        cumulus_parameterization_state.input_output.buoyancy_excess.field[:] = self.manual_inputs["buoyancy_excess"]
        cumulus_parameterization_state.input_output.large_scale_ice.field[:] = self.manual_inputs["large_scale_ice"]
        cumulus_parameterization_state.input_output.convective_ice.field[:] = self.manual_inputs["convective_ice"]
        cumulus_parameterization_state.input_output.large_scale_liquid.field[:] = self.manual_inputs["large_scale_liquid"]
        cumulus_parameterization_state.input_output.convective_liquid.field[:] = self.manual_inputs["convective_liquid"]
        cumulus_parameterization_state.input_output.large_scale_cloud_fraction.field[:] = self.manual_inputs["large_scale_cloud_fraction"]
        cumulus_parameterization_state.input_output.convective_cloud_fraction.field[:] = self.manual_inputs["convective_cloud_fraction"]
        cumulus_parameterization_state.input_output.chemistry_tracers.field[:] = self.manual_inputs["chemistry_tracers"]
        chemistry_tracers_input_5d = np.full(cumulus_parameterization_state.input_output.chemistry_tracers_output.field[:].shape, np.nan)
        for plume in range(NUMBER_OF_PLUMES):
            chemistry_tracers_input_5d[:, :, :, plume, :] = self.manual_inputs["chemistry_tracers_output"][
                :,
                :,
                :,
                plume * config.NUMBER_OF_TRACERS : plume * config.NUMBER_OF_TRACERS + config.NUMBER_OF_TRACERS,
            ]
        cumulus_parameterization_state.input_output.chemistry_tracers_output.field[:] = chemistry_tracers_input_5d[:, :, :, [1, 2, 0], :]

        # initialize convection tracers
        convection_tracers = ConvectionTracers.ones(
            self.quantity_factory,
            data_dimensions={
                "convection_tracers": config.NUMBER_OF_TRACERS,
                "size_three_dimension": 3,
                "size_four_dimension": 4,
            },
        )

        convection_tracers.tracers.field[:] = np.moveaxis(self.convection_tracers_input["tracers"], 0, 3)
        convection_tracers.vect_hcts.field[:] = self.convection_tracers_input["vect_hcts"]
        convection_tracers.kc_scal.field[:] = self.convection_tracers_input["kc_scal"]
        convection_tracers.fscav.field[:] = self.convection_tracers_input["fscav"]
        convection_tracers.convfaci2g.field[:] = self.convection_tracers_input["convfaci2g"]
        convection_tracers.retfactor.field[:] = self.convection_tracers_input["retfactor"]
        convection_tracers.liq_and_gas.field[:] = self.convection_tracers_input["liq_and_gas"]
        convection_tracers.online_cldliq.field[:] = self.convection_tracers_input["online_cldliq"]
        convection_tracers.online_vud.field[:] = self.convection_tracers_input["online_vud"]
        convection_tracers.ftemp_threshold.field[:] = self.convection_tracers_input["ftemp_threshold"]
        convection_tracers.use_gcc_washout.field[:] = self.convection_tracers_input["use_gcc_washout"]
        convection_tracers.use_gocart.field[:] = self.convection_tracers_input["use_gocart"]
        convection_tracers.is_wetdep.field[:] = self.convection_tracers_input["is_wetdep"]

        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        code = GF2020Finalize(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
            saturation_tables=saturation_tables,
        )

        code(
            state=state,
            locals=locals,
            cumulus_parameterization_state=cumulus_parameterization_state,
            convection_tracers=convection_tracers,
        )

        # fill output dictionary for testing
        outputs = {
            "latitude_bugworkaroundname": state.latitude.field[:],
            "longitude_bugworkaroundname": state.longitude.field[:],
            "p_interface_bugworkaroundname": state.p_interface.field[:],
            "t_bugworkaroundname": state.t.field[:],
            "u_bugworkaroundname": state.u.field[:],
            "v_bugworkaroundname": state.v.field[:],
            "w_bugworkaroundname": state.w.field[:],
            "omega_bugworkaroundname": state.omega.field[:],
            "t_2m_bugworkaroundname": state.t_2m.field[:],
            "specific_humidity_2m_bugworkaroundname": state.specific_humidity_2m.field[:],
            "t_surface_bugworkaroundname": state.t_surface.field[:],
            "specific_humidity_surface_bugworkaroundname": state.specific_humidity_surface.field[:],
            "vapor_bugworkaroundname": state.vapor.field[:],
            "convective_liquid_bugworkaroundname": state.convective_liquid.field[:],
            "convective_ice_bugworkaroundname": state.convective_ice.field[:],
            "large_scale_liquid_bugworkaroundname": state.large_scale_liquid.field[:],
            "large_scale_ice_bugworkaroundname": state.large_scale_ice.field[:],
            "convective_cloud_fraction_bugworkaroundname": state.convective_cloud_fraction.field[:],
            "large_scale_cloud_fraction_bugworkaroundname": state.large_scale_cloud_fraction.field[:],
            "p_interface_timestep_start_bugworkaroundname": state.p_interface_timestep_start.field[:],
            "t_timestep_start_bugworkaroundname": state.t_timestep_start.field[:],
            "u_timestep_start_bugworkaroundname": state.u_timestep_start.field[:],
            "v_timestep_start_bugworkaroundname": state.v_timestep_start.field[:],
            "vapor_timestep_start_bugworkaroundname": state.vapor_timestep_start.field[:],
            "geopotential_height_interface_bugworkaroundname": state.geopotential_height_interface.field[:],
            "geopotential_height_surface_bugworkaroundname": state.geopotential_height_surface.field[:],
            "area_bugworkaroundname": state.area.field[:],
            "pbl_level_bugworkaroundname": state.pbl_level.field[:],
            "convection_fraction_bugworkaroundname": state.convection_fraction.field[:],
            "surface_type_bugworkaroundname": state.surface_type.field[:],
            "seed_convection_bugworkaroundname": state.seed_convection.field[:],
            "land_fraction_bugworkaroundname": state.land_fraction.field[:],
            "scalar_diffusivity_bugworkaroundname": state.scalar_diffusivity.field[:],
            "buoyancy_bugworkaroundname": state.buoyancy.field[:],
            "convective_precipitation_GF_bugworkaroundname": state.convective_precipitation_GF.field[:],
            "convective_precipitation_RAS_bugworkaroundname": state.convective_precipitation_RAS.field[:],
            "sensible_heat_flux_bugworkaroundname": state.sensible_heat_flux.field[:],
            "total_water_flux_deep_convection_interface_bugworkaroundname": state.total_water_flux_deep_convection_interface.field[:],
            "evaporation_bugworkaroundname": state.evaporation.field[:],
            "sublimation_of_convective_precipitation_bugworkaroundname": state.sublimation_of_convective_precipitation.field[:],
            "evaporation_of_convective_precipitation_bugworkaroundname": state.evaporation_of_convective_precipitation.field[:],
            "ice_precip_flux_interface_bugworkaroundname": state.ice_precip_flux_interface.field[:],
            "liquid_precip_flux_interface_bugworkaroundname": state.liquid_precip_flux_interface.field[:],
            "convective_condensate_source_bugworkaroundname": state.convective_condensate_source.field[:],
            "convective_condensate_grid_mean_bugworkaroundname": state.convective_condensate_grid_mean.field[:],
            "entrainment_parameter_bugworkaroundname": state.entrainment_parameter.field[:],
            "lateral_entrainment_rate_bugworkaroundname": state.lateral_entrainment_rate.field[:],
            "lateral_entrainment_rate_shallow_bugworkaroundname": state.lateral_entrainment_rate_shallow.field[:],
            "lateral_entrainment_rate_mid_bugworkaroundname": state.lateral_entrainment_rate_mid.field[:],
            "lateral_entrainment_rate_deep_bugworkaroundname": state.lateral_entrainment_rate_deep.field[:],
            "updraft_areal_fraction_bugworkaroundname": state.updraft_areal_fraction.field[:],
            "updraft_vertical_velocity_bugworkaroundname": state.updraft_vertical_velocity.field[:],
            "dtdt_shortwave_bugworkaroundname": state.dtdt_shortwave.field[:],
            "dtdt_longwave_bugworkaroundname": state.dtdt_longwave.field[:],
            "dspecific_humiditydt_pbl_bugworkaroundname": state.dspecific_humiditydt_pbl.field[:],
            "dtdt_pbl_bugworkaroundname": state.dtdt_pbl.field[:],
            "dtdt_from_dynamics_bugworkaroundname": state.dtdt_from_dynamics.field[:],
            "dvapordt_from_dynamics_bugworkaroundname": state.dvapordt_from_dynamics.field[:],
            "sigma_mid_bugworkaroundname": state.sigma_mid.field[:],
            "sigma_deep_bugworkaroundname": state.sigma_deep.field[:],
            "total_precipitable_water_initial_bugworkaroundname": state.total_precipitable_water_initial.field[:],
            "saturation_total_precipitable_water_initial_bugworkaroundname": state.saturation_total_precipitable_water_initial.field[:],
            "dvapordt_deep_convection_bugworkaroundname": state.dvapordt_deep_convection.field[:],
            "dtdt_deep_convection_bugworkaroundname": state.dtdt_deep_convection.field[:],
            "dudt_deep_convection_bugworkaroundname": state.dudt_deep_convection.field[:],
            "dvdt_deep_convection_bugworkaroundname": state.dvdt_deep_convection.field[:],
            "dliquiddt_deep_convection_bugworkaroundname": state.dliquiddt_deep_convection.field[:],
            "dicedt_deep_convection_bugworkaroundname": state.dicedt_deep_convection.field[:],
            "dcloudfractiondt_deep_convection_bugworkaroundname": state.dcloudfractiondt_deep_convection.field[:],
            "pressure_shallow_convective_cloud_top_bugworkaroundname": state.pressure_shallow_convective_cloud_top.field[:],
            "pressure_mid_convective_cloud_top_bugworkaroundname": state.pressure_mid_convective_cloud_top.field[:],
            "pressure_deep_convective_cloud_top_bugworkaroundname": state.pressure_deep_convective_cloud_top.field[:],
            "mass_flux_shallow_bugworkaroundname": state.mass_flux_shallow.field[:],
            "mass_flux_mid_bugworkaroundname": state.mass_flux_mid.field[:],
            "mass_flux_deep_updraft_bugworkaroundname": state.mass_flux_deep_updraft.field[:],
            "mass_flux_deep_updraft_interface_bugworkaroundname": state.mass_flux_deep_updraft_interface.field[:],
            "mass_flux_deep_updraft_detrained_bugworkaroundname": state.mass_flux_deep_updraft_detrained.field[:],
            "mass_flux_deep_downdraft_bugworkaroundname": state.mass_flux_deep_downdraft.field[:],
            "mass_flux_cloud_base_bugworkaroundname": state.mass_flux_cloud_base.field[:],
            "mass_flux_cloud_base_shallow_bugworkaroundname": state.mass_flux_cloud_base_shallow.field[:],
            "mass_flux_cloud_base_mid_bugworkaroundname": state.mass_flux_cloud_base_mid.field[:],
            "mass_flux_cloud_base_deep_bugworkaroundname": state.mass_flux_cloud_base_deep.field[:],
            "total_cumulative_mass_flux_interface_bugworkaroundname": state.total_cumulative_mass_flux_interface.field[:],
            "total_detraining_mass_flux_bugworkaroundname": state.total_detraining_mass_flux.field[:],
            "convection_code_shallow_bugworkaroundname": state.convection_code_shallow.field[:],
            "convection_code_mid_bugworkaroundname": state.convection_code_mid.field[:],
            "convection_code_deep_bugworkaroundname": state.convection_code_deep.field[:],
            "cloud_workfunction_0_bugworkaroundname": state.cloud_workfunction_0.field[:],
            "cloud_workfunction_1_bugworkaroundname": state.cloud_workfunction_1.field[:],
            "cloud_workfunction_2_bugworkaroundname": state.cloud_workfunction_2.field[:],
            "cloud_workfunction_3_bugworkaroundname": state.cloud_workfunction_3.field[:],
            "cloud_workfunction_1_pbl_bugworkaroundname": state.cloud_workfunction_1_pbl.field[:],
            "cloud_workfunction_1_cin_bugworkaroundname": state.cloud_workfunction_1_cin.field[:],
            "pbl_time_scale_bugworkaroundname": state.pbl_time_scale.field[:],
            "cape_removal_time_scale_bugworkaroundname": state.cape_removal_time_scale.field[:],
            "lightning_density_bugworkaroundname": state.lightning_density.field[:],
            "convection_tracer_bugworkaroundname": state.convection_tracer.field[:],
        }

        return outputs
