import numpy as np
from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py

from pyMoist.convection.GF_2020 import GF2020, GF2020Config, GF2020CumulusParameterizationConfig, GF2020State
from pyMoist.convection_tracers import ConvectionTracers
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGF2020(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        _namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, stencil_factory)
        self.quantity_factory = grid.quantity_factory

        self.in_vars["data_vars"] = {
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
            "total_water_flux_deep_convection_interface_bugworkaroundname": {},
            "evaporation_bugworkaroundname": {},
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
            "dvapordt_deep_convection_bugworkaroundname": {},
            "dtdt_deep_convection_bugworkaroundname": {},
            "dudt_deep_convection_bugworkaroundname": {},
            "dvdt_deep_convection_bugworkaroundname": {},
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

        self.out_vars = self.in_vars["data_vars"].copy()
        del self.out_vars[
            "total_water_flux_deep_convection_interface_bugworkaroundname"
        ]  # disabled temporarily - do not push, fails b/c of <10 ulp diff in saturation functions

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.cu_param_constants = data_loader.load("GF2020_CumulusParameterization-constants")
        self.convection_tracers_input = data_loader.load("GF2020_ConvectionTracers")

        # workaround because translate test cannot read in 4d fields
        self.manual_inputs = data_loader.load("GF2020_CumulusParameterization-Out", use_dynamic_i_call=True)

        # load data from GF2020-In to fill in unmodified fields from the state
        self.unmodified_state = data_loader.load("GF2020-In", use_dynamic_i_call=True)

    def compute(self, inputs):
        config = GF2020Config(**self.constants)
        cumulus_parameterization_config = GF2020CumulusParameterizationConfig(**self.cu_param_constants)

        # initialize saturation tables
        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

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

        # initialize GF2020 state
        state = GF2020State.zeros(self.quantity_factory)

        state.latitude.field[:] = inputs["latitude_bugworkaroundname"]
        state.longitude.field[:] = inputs["longitude_bugworkaroundname"]
        state.p_interface.field[:] = inputs["p_interface_bugworkaroundname"]
        state.t.field[:] = inputs["t_bugworkaroundname"]
        state.u.field[:] = inputs["u_bugworkaroundname"]
        state.v.field[:] = inputs["v_bugworkaroundname"]
        state.w.field[:] = inputs["w_bugworkaroundname"]
        state.omega.field[:] = inputs["omega_bugworkaroundname"]
        state.t_2m.field[:] = inputs["t_2m_bugworkaroundname"]
        state.specific_humidity_2m.field[:] = inputs["specific_humidity_2m_bugworkaroundname"]
        state.t_surface.field[:] = inputs["t_surface_bugworkaroundname"]
        state.specific_humidity_surface.field[:] = inputs["specific_humidity_surface_bugworkaroundname"]
        state.vapor.field[:] = inputs["vapor_bugworkaroundname"]
        state.convective_liquid.field[:] = inputs["convective_liquid_bugworkaroundname"]
        state.convective_ice.field[:] = inputs["convective_ice_bugworkaroundname"]
        state.large_scale_liquid.field[:] = inputs["large_scale_liquid_bugworkaroundname"]
        state.large_scale_ice.field[:] = inputs["large_scale_ice_bugworkaroundname"]
        state.convective_cloud_fraction.field[:] = inputs["convective_cloud_fraction_bugworkaroundname"]
        state.large_scale_cloud_fraction.field[:] = inputs["large_scale_cloud_fraction_bugworkaroundname"]
        state.p_interface_timestep_start.field[:] = inputs["p_interface_timestep_start_bugworkaroundname"]
        state.t_timestep_start.field[:] = inputs["t_timestep_start_bugworkaroundname"]
        state.u_timestep_start.field[:] = inputs["u_timestep_start_bugworkaroundname"]
        state.v_timestep_start.field[:] = inputs["v_timestep_start_bugworkaroundname"]
        state.vapor_timestep_start.field[:] = inputs["vapor_timestep_start_bugworkaroundname"]
        state.geopotential_height_interface.field[:] = inputs["geopotential_height_interface_bugworkaroundname"]
        state.geopotential_height_surface.field[:] = inputs["geopotential_height_surface_bugworkaroundname"]
        state.area.field[:] = inputs["area_bugworkaroundname"]
        state.pbl_level.field[:] = inputs["pbl_level_bugworkaroundname"] - 1
        state.convection_fraction.field[:] = inputs["convection_fraction_bugworkaroundname"]
        state.surface_type.field[:] = inputs["surface_type_bugworkaroundname"]
        state.land_fraction.field[:] = inputs["land_fraction_bugworkaroundname"]
        state.scalar_diffusivity.field[:] = inputs["scalar_diffusivity_bugworkaroundname"]
        state.buoyancy.field[:] = inputs["buoyancy_bugworkaroundname"]
        state.convective_precipitation_GF.field[:] = inputs["convective_precipitation_GF_bugworkaroundname"]
        state.convective_precipitation_RAS.field[:] = inputs["convective_precipitation_RAS_bugworkaroundname"]
        state.sensible_heat_flux.field[:] = inputs["sensible_heat_flux_bugworkaroundname"]
        state.total_water_flux_deep_convection_interface.field[:] = inputs["total_water_flux_deep_convection_interface_bugworkaroundname"]
        state.evaporation.field[:] = inputs["evaporation_bugworkaroundname"]
        state.convective_condensate_source.field[:] = inputs["convective_condensate_source_bugworkaroundname"]
        state.convective_condensate_grid_mean.field[:] = inputs["convective_condensate_grid_mean_bugworkaroundname"]
        state.entrainment_parameter.field[:] = inputs["entrainment_parameter_bugworkaroundname"]
        state.lateral_entrainment_rate.field[:] = inputs["lateral_entrainment_rate_bugworkaroundname"]
        state.lateral_entrainment_rate_shallow.field[:] = inputs["lateral_entrainment_rate_shallow_bugworkaroundname"]
        state.lateral_entrainment_rate_mid.field[:] = inputs["lateral_entrainment_rate_mid_bugworkaroundname"]
        state.lateral_entrainment_rate_deep.field[:] = inputs["lateral_entrainment_rate_deep_bugworkaroundname"]
        state.updraft_areal_fraction.field[:] = inputs["updraft_areal_fraction_bugworkaroundname"]
        state.updraft_vertical_velocity.field[:] = inputs["updraft_vertical_velocity_bugworkaroundname"]
        state.dtdt_shortwave.field[:] = inputs["dtdt_shortwave_bugworkaroundname"]
        state.dtdt_longwave.field[:] = inputs["dtdt_longwave_bugworkaroundname"]
        state.dspecific_humiditydt_pbl.field[:] = inputs["dspecific_humiditydt_pbl_bugworkaroundname"]
        state.dtdt_pbl.field[:] = inputs["dtdt_pbl_bugworkaroundname"]
        state.dtdt_from_dynamics.field[:] = inputs["dtdt_from_dynamics_bugworkaroundname"]
        state.dvapordt_from_dynamics.field[:] = inputs["dvapordt_from_dynamics_bugworkaroundname"]
        state.sigma_mid.field[:] = inputs["sigma_mid_bugworkaroundname"]
        state.sigma_deep.field[:] = inputs["sigma_deep_bugworkaroundname"]
        state.dvapordt_deep_convection.field[:] = inputs["dvapordt_deep_convection_bugworkaroundname"]
        state.dtdt_deep_convection.field[:] = inputs["dtdt_deep_convection_bugworkaroundname"]
        state.dudt_deep_convection.field[:] = inputs["dudt_deep_convection_bugworkaroundname"]
        state.dvdt_deep_convection.field[:] = inputs["dvdt_deep_convection_bugworkaroundname"]
        state.pressure_shallow_convective_cloud_top.field[:] = inputs["pressure_shallow_convective_cloud_top_bugworkaroundname"]
        state.pressure_mid_convective_cloud_top.field[:] = inputs["pressure_mid_convective_cloud_top_bugworkaroundname"]
        state.pressure_deep_convective_cloud_top.field[:] = inputs["pressure_deep_convective_cloud_top_bugworkaroundname"]
        state.mass_flux_shallow.field[:] = inputs["mass_flux_shallow_bugworkaroundname"]
        state.mass_flux_mid.field[:] = inputs["mass_flux_mid_bugworkaroundname"]
        state.mass_flux_deep_updraft.field[:] = inputs["mass_flux_deep_updraft_bugworkaroundname"]
        state.mass_flux_deep_updraft_interface.field[:] = inputs["mass_flux_deep_updraft_interface_bugworkaroundname"]
        state.mass_flux_deep_updraft_detrained.field[:] = inputs["mass_flux_deep_updraft_detrained_bugworkaroundname"]
        state.mass_flux_deep_downdraft.field[:] = inputs["mass_flux_deep_downdraft_bugworkaroundname"]
        state.mass_flux_cloud_base.field[:] = inputs["mass_flux_cloud_base_bugworkaroundname"]
        state.mass_flux_cloud_base_shallow.field[:] = inputs["mass_flux_cloud_base_shallow_bugworkaroundname"]
        state.mass_flux_cloud_base_mid.field[:] = inputs["mass_flux_cloud_base_mid_bugworkaroundname"]
        state.mass_flux_cloud_base_deep.field[:] = inputs["mass_flux_cloud_base_deep_bugworkaroundname"]
        state.convection_code_shallow.field[:] = inputs["convection_code_shallow_bugworkaroundname"]
        state.convection_code_mid.field[:] = inputs["convection_code_mid_bugworkaroundname"]
        state.convection_code_deep.field[:] = inputs["convection_code_deep_bugworkaroundname"]
        state.cloud_workfunction_0.field[:] = inputs["cloud_workfunction_0_bugworkaroundname"]
        state.cloud_workfunction_1.field[:] = inputs["cloud_workfunction_1_bugworkaroundname"]
        state.cloud_workfunction_2.field[:] = inputs["cloud_workfunction_2_bugworkaroundname"]
        state.cloud_workfunction_3.field[:] = inputs["cloud_workfunction_3_bugworkaroundname"]
        state.cloud_workfunction_1_pbl.field[:] = inputs["cloud_workfunction_1_pbl_bugworkaroundname"]
        state.cloud_workfunction_1_cin.field[:] = inputs["cloud_workfunction_1_cin_bugworkaroundname"]
        state.pbl_time_scale.field[:] = inputs["pbl_time_scale_bugworkaroundname"]
        state.cape_removal_time_scale.field[:] = inputs["cape_removal_time_scale_bugworkaroundname"]
        state.lightning_density.field[:] = inputs["lightning_density_bugworkaroundname"]
        state.convection_tracer.field[:] = inputs["convection_tracer_bugworkaroundname"]

        # initialize test code
        code = GF2020(
            stencil_factory=self.stencil_factory,
            quantity_factory=self.quantity_factory,
            config=config,
            cumulus_parameterization_config=cumulus_parameterization_config,
            saturation_tables=saturation_tables,
        )

        # run test code
        code(
            state=state,
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
            "pbl_level_bugworkaroundname": state.pbl_level.field[:] + 1,
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
