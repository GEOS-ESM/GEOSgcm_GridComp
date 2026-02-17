from f90nml import Namelist
from ndsl import StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.setup import GF2020Setup
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.constants import (
    NUMBER_OF_PLUMES,
)
from pyMoist.constants import NUMBER_OF_TRACERS
import numpy as np
from pyMoist.convection_tracers import ConvectionTracers
from ndsl.constants import X_DIM, Y_DIM, Z_DIM


class TranslateGF2020_Setup(TranslateFortranData2Py):
    def __init__(
        self,
        grid: Grid,
        namelist: Namelist,
        stencil_factory: StencilFactory,
    ):
        super().__init__(grid, namelist, stencil_factory)
        self.stencil_factory = stencil_factory
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
            "total_water_flux_deep_convection_bugworkaroundname": {},
            "evaporation_bugworkaroundname": {},
            "convective_condensate_source_bugworkaroundname": {},
            "convective_condensate_grid_mean_bugworkaroundname": {},
            "entrainment_parameter_bugworkaroundname": {},
            "lateral_entrainment_rate_bugworkaroundname": {},
            "lateral_entrainment_rate_shallow_bugworkaroundname": {},
            "lateral_entrainment_rate_mid_bugworkaroundname": {},
            "lateral_entrainment_rate_deep_bugworkaroundname": {},
            "updraft_area_fraction_bugworkaroundname": {},
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
            "pressure_shallow_convective_cloud_top_bugworkaroundname": {},
            "pressure_mid_convective_cloud_top_bugworkaroundname": {},
            "pressure_deep_convective_cloud_top_bugworkaroundname": {},
            "mass_flux_shalow_bugworkaroundname": {},
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
            "cloud_work_function_0_bugworkaroundname": {},
            "cloud_work_function_1_bugworkaroundname": {},
            "cloud_work_function_2_bugworkaroundname": {},
            "cloud_work_function_3_bugworkaroundname": {},
            "cloud_work_function_1_pbl_bugworkaroundname": {},
            "cloud_work_function_1_cin_bugworkaroundname": {},
            "pbl_time_scale_bugworkaroundname": {},
            "cape_removal_time_scale_bugworkaroundname": {},
            "lighting_density_bugworkaroundname": {},
            "convection_tracer_bugworkaroundname": {},
        }

        # NOTE disabled fields are nan in fortran - zero in python, disabled so the test passes
        self.out_vars = {
            # internal locals
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
            # CumulusParameterization state - input
            "t_excess": {},
            "vapor_excess": {},
            "grid_scale_forcing_t": {},
            "grid_scale_forcing_vapor": {},
            "subgrid_scale_forcing_t": {},
            "subgrid_scale_forcing_vapor": {},
            "seed_convection": {},
            # "saturation_water_vapor": {},
            "ocean_fraction": {},
            "convection_fraction": {},
            "surface_type": {},
            "lateral_entrainment_rate": {},
            "last_error_code": {},
            # CumulusParameterization state - output
            "dtdt": {},
            "dvapordt": {},
            "dcloudicedt": {},
            "dudt": {},
            "dvdt": {},
            "dnliquiddt": {},
            "dnicedt": {},
            "dbuoyancydt": {},
            # "dconvectiveicedt": {},
            # "dlargescaleicedt": {},
            # "dconvectiveliquiddt": {},
            # "dlargescaleliquiddt": {},
            # "dconvectivecloudfractiondt": {},
            # "dlargescalecloudfractiondt": {},
            "error_code": {},
            "downdraft_origin_level": {},
            "lcl_level": {},
            "updraft_origin_level": {},
            "updraft_lfc_level": {},
            "cloud_top_level": {},
            "kstabi": {},
            "kstabm": {},
            "precip": {},
            "cloud_base_mass_flux_modified": {},
            "epsilon_forced": {},
            "total_normalized_integrated_condensate_forced": {},
            "scale_dependence_factor": {},
            "p_cloud_levels_forced": {},
            "entrainment_rate": {},
            "mass_entrainment_updraft_forced": {},
            "mass_entrainment_downdraft_forced": {},
            "mass_detrainment_updraft_forced": {},
            "mass_detrainment_downdraft_forced": {},
            "normalized_massflux_updraft_forced": {},
            "normalized_massflux_downdraft_forced": {},
            "condensate_to_fall_forced": {},
            "evaporate_in_downdraft_forced": {},
            "cloud_liquid_after_rain_forced": {},
            "t_updraft": {},
            "convective_cloud_fraction_output": {},
            "cloud_workfunction_0": {},
            "cloud_workfunction_1": {},
            "cloud_workfunction_2": {},
            "cloud_workfunction_3": {},
            "cloud_workfunction_1_pbl": {},
            "cloud_workfunction_1_cin": {},
            "cape_removal_time_scale": {},
            "pbl_time_scale": {},
            "lightning_density": {},
            "evaporation_sublimation_tendency": {},
            "convective_precip_flux": {},
            "t_perturbation": {},
            # CumulusParameterization state - input-output
            "grid_length": {},
            "pbl_level": {},
            "ccn": {},
            "air_density": {},
            "omega": {},
            "topography_height_no_negative": {},
            "sensible_heat_flux": {},
            "latent_heat_flux": {},
            "longitude_degrees": {},
            "latitude_degrees": {},
            "t_old": {},
            "vapor_old": {},
            "t_modified_by_advection": {},
            "vapor_modified_by_advection": {},
            "geopotential_height_forced": {},
            "p_forced": {},
            "p_surface": {},
            "t_surface": {},
            "u": {},
            "v": {},
            "w": {},
            "mass": {},
            "convective_scale_velocity": {},
            "buoyancy_excess": {},
            # "large_scale_ice": {},
            # "convective_ice": {},
            # "large_scale_liquid": {},
            # "convective_liquid": {},
            # "large_scale_cloud_fraction": {},
            # "convective_cloud_fraction": {},
            "chemistry_tracers": {},
            "chemistry_tracers_output": {},
        }

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")
        self.convection_tracers_input = data_loader.load("GF2020_ConvectionTracers")

    def compute(self, inputs):
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)

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

        # fill GF2020 state with input data
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
        state.geopotential_height_interface.field[:] = inputs[
            "geopotential_height_interface_bugworkaroundname"
        ]
        state.geopotential_height_surface.field[:] = inputs["geopotential_height_surface_bugworkaroundname"]
        state.area.field[:] = inputs["area_bugworkaroundname"]
        state.pbl_level.field[:] = inputs["pbl_level_bugworkaroundname"]
        state.convection_fraction.field[:] = inputs["convection_fraction_bugworkaroundname"]
        state.surface_type.field[:] = inputs["surface_type_bugworkaroundname"]
        state.seed_convection.field[:] = inputs["seed_convection_bugworkaroundname"]
        state.land_fraction.field[:] = inputs["land_fraction_bugworkaroundname"]
        state.scalar_diffusivity.field[:] = inputs["scalar_diffusivity_bugworkaroundname"]
        state.buoyancy.field[:] = inputs["buoyancy_bugworkaroundname"]
        state.convective_precipitation_GF.field[:] = inputs["convective_precipitation_GF_bugworkaroundname"]
        state.convective_precipitation_RAS.field[:] = inputs["convective_precipitation_RAS_bugworkaroundname"]
        state.sensible_heat_flux.field[:] = inputs["sensible_heat_flux_bugworkaroundname"]
        state.total_water_flux_deep_convection.field[:] = inputs[
            "total_water_flux_deep_convection_bugworkaroundname"
        ]
        state.evaporation.field[:] = inputs["evaporation_bugworkaroundname"]
        state.convective_condensate_source.field[:] = inputs["convective_condensate_source_bugworkaroundname"]
        state.convective_condensate_grid_mean.field[:] = inputs[
            "convective_condensate_grid_mean_bugworkaroundname"
        ]
        state.entrainment_parameter.field[:] = inputs["entrainment_parameter_bugworkaroundname"]
        state.lateral_entrainment_rate.field[:] = inputs["lateral_entrainment_rate_bugworkaroundname"]
        state.lateral_entrainment_rate_shallow.field[:] = inputs[
            "lateral_entrainment_rate_shallow_bugworkaroundname"
        ]
        state.lateral_entrainment_rate_mid.field[:] = inputs["lateral_entrainment_rate_mid_bugworkaroundname"]
        state.lateral_entrainment_rate_deep.field[:] = inputs[
            "lateral_entrainment_rate_deep_bugworkaroundname"
        ]
        state.updraft_area_fraction.field[:] = inputs["updraft_area_fraction_bugworkaroundname"]
        state.updraft_vertical_velocity.field[:] = inputs["updraft_vertical_velocity_bugworkaroundname"]
        state.dtdt_shortwave.field[:] = inputs["dtdt_shortwave_bugworkaroundname"]
        state.dtdt_longwave.field[:] = inputs["dtdt_longwave_bugworkaroundname"]
        state.dspecific_humiditydt_pbl.field[:] = inputs["dspecific_humiditydt_pbl_bugworkaroundname"]
        state.dtdt_pbl.field[:] = inputs["dtdt_pbl_bugworkaroundname"]
        state.dtdt_from_dynamics.field[:] = inputs["dtdt_from_dynamics_bugworkaroundname"]
        state.dvapordt_from_dynamics.field[:] = inputs["dvapordt_from_dynamics_bugworkaroundname"]
        state.sigma_mid.field[:] = inputs["sigma_mid_bugworkaroundname"]
        state.sigma_deep.field[:] = inputs["sigma_deep_bugworkaroundname"]
        state.total_precipitable_water_initial.field[:] = inputs[
            "total_precipitable_water_initial_bugworkaroundname"
        ]
        state.saturation_total_precipitable_water_initial.field[:] = inputs[
            "saturation_total_precipitable_water_initial_bugworkaroundname"
        ]
        state.dvapordt_deep_convection.field[:] = inputs["dvapordt_deep_convection_bugworkaroundname"]
        state.dtdt_deep_convection.field[:] = inputs["dtdt_deep_convection_bugworkaroundname"]
        state.dudt_deep_convection.field[:] = inputs["dudt_deep_convection_bugworkaroundname"]
        state.dvdt_deep_convection.field[:] = inputs["dvdt_deep_convection_bugworkaroundname"]
        state.pressure_shallow_convective_cloud_top.field[:] = inputs[
            "pressure_shallow_convective_cloud_top_bugworkaroundname"
        ]
        state.pressure_mid_convective_cloud_top.field[:] = inputs[
            "pressure_mid_convective_cloud_top_bugworkaroundname"
        ]
        state.pressure_deep_convective_cloud_top.field[:] = inputs[
            "pressure_deep_convective_cloud_top_bugworkaroundname"
        ]
        state.mass_flux_shalow.field[:] = inputs["mass_flux_shalow_bugworkaroundname"]
        state.mass_flux_mid.field[:] = inputs["mass_flux_mid_bugworkaroundname"]
        state.mass_flux_deep_updraft.field[:] = inputs["mass_flux_deep_updraft_bugworkaroundname"]
        state.mass_flux_deep_updraft_interface.field[:] = inputs[
            "mass_flux_deep_updraft_interface_bugworkaroundname"
        ]
        state.mass_flux_deep_updraft_detrained.field[:] = inputs[
            "mass_flux_deep_updraft_detrained_bugworkaroundname"
        ]
        state.mass_flux_deep_downdraft.field[:] = inputs["mass_flux_deep_downdraft_bugworkaroundname"]
        state.mass_flux_cloud_base.field[:] = inputs["mass_flux_cloud_base_bugworkaroundname"]
        state.mass_flux_cloud_base_shallow.field[:] = inputs["mass_flux_cloud_base_shallow_bugworkaroundname"]
        state.mass_flux_cloud_base_mid.field[:] = inputs["mass_flux_cloud_base_mid_bugworkaroundname"]
        state.mass_flux_cloud_base_deep.field[:] = inputs["mass_flux_cloud_base_deep_bugworkaroundname"]
        state.convection_code_shallow.field[:] = inputs["convection_code_shallow_bugworkaroundname"]
        state.convection_code_mid.field[:] = inputs["convection_code_mid_bugworkaroundname"]
        state.convection_code_deep.field[:] = inputs["convection_code_deep_bugworkaroundname"]
        state.cloud_workfunction_0.field[:] = inputs["cloud_work_function_0_bugworkaroundname"]
        state.cloud_workfunction_1.field[:] = inputs["cloud_work_function_1_bugworkaroundname"]
        state.cloud_workfunction_2.field[:] = inputs["cloud_work_function_2_bugworkaroundname"]
        state.cloud_workfunction_3.field[:] = inputs["cloud_work_function_3_bugworkaroundname"]
        state.cloud_workfunction_1_pbl.field[:] = inputs["cloud_work_function_1_pbl_bugworkaroundname"]
        state.cloud_workfunction_1_cin.field[:] = inputs["cloud_work_function_1_cin_bugworkaroundname"]
        state.pbl_time_scale.field[:] = inputs["pbl_time_scale_bugworkaroundname"]
        state.cape_removal_time_scale.field[:] = inputs["cape_removal_time_scale_bugworkaroundname"]
        state.lightning_density.field[:] = inputs["lighting_density_bugworkaroundname"]
        state.convection_tracer.field[:] = inputs["convection_tracer_bugworkaroundname"]

        # initialize GF2020 locals
        locals = GF2020Locals.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
                "convection_tracers": config.NUMBER_OF_TRACERS,
            },
        )

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

        setup = GF2020Setup(
            stencil_factory=self.stencil_factory, quantity_factory=self.quantity_factory, config=config
        )

        setup(
            state=state,
            locals=locals,
            cumulus_parameterization_state=cumulus_parameterization_state,
            saturation_tables=saturation_tables,
            convection_tracers=convection_tracers,
            scm_stop=False,
        )

        # fill output dictionary for testing

        # collapse plume dim for chemistry_tracers_output
        # NOTE ideally this has no numpy dependency
        chemistry_tracers_output_5d_reordered = (
            cumulus_parameterization_state.input_output.chemistry_tracers_output.field[:, :, :, [2, 0, 1], :]
        )
        grid_size = self.stencil_factory.grid_indexing.get_shape([X_DIM, Y_DIM, Z_DIM])
        chemistry_tracers_output_4d = np.full(
            [grid_size[0], grid_size[1], grid_size[2], NUMBER_OF_PLUMES * NUMBER_OF_TRACERS], np.nan
        )
        for plume in range(NUMBER_OF_PLUMES):
            chemistry_tracers_output_4d[
                :,
                :,
                :,
                plume * config.NUMBER_OF_TRACERS : plume * config.NUMBER_OF_TRACERS
                + config.NUMBER_OF_TRACERS,
            ] = chemistry_tracers_output_5d_reordered[:, :, :, plume, :]

        # fill top level of a few fields with nans to make test pass
        # these levels are NEVER read, so they don't need to be tested anyway
        cumulus_parameterization_state.input_output.t_old.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.vapor_old.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.air_density.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.t_modified_by_advection.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.vapor_modified_by_advection.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.geopotential_height_forced.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.p_forced.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.u.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.v.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.w.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.mass.field[:, :, -1] = np.nan
        cumulus_parameterization_state.input_output.buoyancy_excess.field[:, :, -1] = np.nan

        outputs = {
            # GF2020 locals
            "local_edge_height_above_surface": locals.derived_state.edge_height_above_surface.field[:],
            "local_layer_height_above_surface": locals.derived_state.layer_height_above_surface.field[:],
            "local_p": locals.derived_state.p.field[:],
            "local_p_kappa": locals.derived_state.p_kappa.field[:],
            "local_th": locals.derived_state.th.field[:],
            "local_mass": locals.derived_state.mass.field[:],
            "local_modified_area": locals.derived_state.modified_area.field[:],
            "local_vertical_velocity": locals.derived_state.vertical_velocity.field[:],
            "local_dz": locals.derived_state.dz.field[:],
            "local_air_density": locals.derived_state.air_density.field[:],
            "local_scalar_diffusivity": np.moveaxis(locals.derived_state.scalar_diffusivity.field[:], 2, 0),
            # GF2020 CumulusParameterization fields
            # input fields
            "t_excess": cumulus_parameterization_state.input.t_excess.field[:],
            "vapor_excess": cumulus_parameterization_state.input.vapor_excess.field[:],
            "grid_scale_forcing_t": cumulus_parameterization_state.input.grid_scale_forcing_t.field[:],
            "grid_scale_forcing_vapor": cumulus_parameterization_state.input.grid_scale_forcing_vapor.field[
                :
            ],
            "subgrid_scale_forcing_t": cumulus_parameterization_state.input.subgrid_scale_forcing_t.field[:],
            "subgrid_scale_forcing_vapor": cumulus_parameterization_state.input.subgrid_scale_forcing_vapor.field[
                :
            ],
            "seed_convection": cumulus_parameterization_state.input.seed_convection.field[:],
            "saturation_water_vapor": cumulus_parameterization_state.input.saturation_water_vapor.field[:],
            "ocean_fraction": cumulus_parameterization_state.input.ocean_fraction.field[:],
            "convection_fraction": cumulus_parameterization_state.input.convection_fraction.field[:],
            "surface_type": cumulus_parameterization_state.input.surface_type.field[:],
            "lateral_entrainment_rate": cumulus_parameterization_state.input.lateral_entrainment_rate.field[
                :
            ],
            "last_error_code": cumulus_parameterization_state.input.last_error_code.field[:],
            # output fields
            "dtdt": cumulus_parameterization_state.output.dtdt.field[:, :, :, [2, 0, 1]],
            "dvapordt": cumulus_parameterization_state.output.dvapordt.field[:, :, :, [2, 0, 1]],
            "dcloudicedt": cumulus_parameterization_state.output.dcloudicedt.field[:, :, :, [2, 0, 1]],
            "dudt": cumulus_parameterization_state.output.dudt.field[:, :, :, [2, 0, 1]],
            "dvdt": cumulus_parameterization_state.output.dvdt.field[:, :, :, [2, 0, 1]],
            "dnliquiddt": cumulus_parameterization_state.output.dnliquiddt.field[:, :, :, [2, 0, 1]],
            "dnicedt": cumulus_parameterization_state.output.dnicedt.field[:, :, :, [2, 0, 1]],
            "dbuoyancydt": cumulus_parameterization_state.output.dbuoyancydt.field[:, :, :, [2, 0, 1]],
            "dconvectiveicedt": cumulus_parameterization_state.output.dconvectiveicedt.field[
                :, :, :, [2, 0, 1]
            ],
            "dlargescaleicedt": cumulus_parameterization_state.output.dlargescaleicedt.field[
                :, :, :, [2, 0, 1]
            ],
            "dconvectiveliquiddt": cumulus_parameterization_state.output.dconvectiveliquiddt.field[
                :, :, :, [2, 0, 1]
            ],
            "dlargescaleliquiddt": cumulus_parameterization_state.output.dlargescaleliquiddt.field[
                :, :, :, [2, 0, 1]
            ],
            "dconvectivecloudfractiondt": cumulus_parameterization_state.output.dconvectivecloudfractiondt.field[
                :, :, :, [2, 0, 1]
            ],
            "dlargescalecloudfractiondt": cumulus_parameterization_state.output.dlargescalecloudfractiondt.field[
                :, :, :, [2, 0, 1]
            ],
            "error_code": cumulus_parameterization_state.output.error_code.field[:, :, [2, 0, 1]],
            "downdraft_origin_level": cumulus_parameterization_state.output.downdraft_origin_level.field[
                :, :, [2, 0, 1]
            ],
            "lcl_level": cumulus_parameterization_state.output.lcl_level.field[:, :, [2, 0, 1]],
            "updraft_origin_level": cumulus_parameterization_state.output.updraft_origin_level.field[
                :, :, [2, 0, 1]
            ],
            "updraft_lfc_level": cumulus_parameterization_state.output.updraft_lfc_level.field[
                :, :, [2, 0, 1]
            ],
            "cloud_top_level": cumulus_parameterization_state.output.cloud_top_level.field[:, :, [2, 0, 1]],
            "kstabi": cumulus_parameterization_state.output.kstabi.field[:, :, [2, 0, 1]],
            "kstabm": cumulus_parameterization_state.output.kstabm.field[:, :, [2, 0, 1]],
            "precip": cumulus_parameterization_state.output.precip.field[:, :, [2, 0, 1]],
            "cloud_base_mass_flux_modified": cumulus_parameterization_state.output.cloud_base_mass_flux_modified.field[
                :, :, [2, 0, 1]
            ],
            "epsilon_forced": cumulus_parameterization_state.output.epsilon_forced.field[:, :, [2, 0, 1]],
            "total_normalized_integrated_condensate_forced": cumulus_parameterization_state.output.total_normalized_integrated_condensate_forced.field[
                :, :, [2, 0, 1]
            ],
            "scale_dependence_factor": cumulus_parameterization_state.output.scale_dependence_factor.field[
                :, :, [2, 0, 1]
            ],
            "p_cloud_levels_forced": cumulus_parameterization_state.output.p_cloud_levels_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "entrainment_rate": cumulus_parameterization_state.output.entrainment_rate.field[
                :, :, :, [2, 0, 1]
            ],
            "mass_entrainment_updraft_forced": cumulus_parameterization_state.output.mass_entrainment_updraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "mass_entrainment_downdraft_forced": cumulus_parameterization_state.output.mass_entrainment_downdraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "mass_detrainment_updraft_forced": cumulus_parameterization_state.output.mass_detrainment_updraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "mass_detrainment_downdraft_forced": cumulus_parameterization_state.output.mass_detrainment_downdraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "normalized_massflux_updraft_forced": cumulus_parameterization_state.output.normalized_massflux_updraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "normalized_massflux_downdraft_forced": cumulus_parameterization_state.output.normalized_massflux_downdraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "condensate_to_fall_forced": cumulus_parameterization_state.output.condensate_to_fall_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "evaporate_in_downdraft_forced": cumulus_parameterization_state.output.evaporate_in_downdraft_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "cloud_liquid_after_rain_forced": cumulus_parameterization_state.output.cloud_liquid_after_rain_forced.field[
                :, :, :, [2, 0, 1]
            ],
            "t_updraft": cumulus_parameterization_state.output.t_updraft.field[:, :, :, [2, 0, 1]],
            "convective_cloud_fraction_output": cumulus_parameterization_state.output.convective_cloud_fraction.field[
                :, :, :, [2, 0, 1]
            ],
            "cloud_workfunction_0": cumulus_parameterization_state.output.cloud_workfunction_0.field[:],
            "cloud_workfunction_1": cumulus_parameterization_state.output.cloud_workfunction_1.field[:],
            "cloud_workfunction_2": cumulus_parameterization_state.output.cloud_workfunction_2.field[:],
            "cloud_workfunction_3": cumulus_parameterization_state.output.cloud_workfunction_3.field[:],
            "cloud_workfunction_1_pbl": cumulus_parameterization_state.output.cloud_workfunction_1_pbl.field[
                :
            ],
            "cloud_workfunction_1_cin": cumulus_parameterization_state.output.cloud_workfunction_1_cin.field[
                :
            ],
            "cape_removal_time_scale": cumulus_parameterization_state.output.cape_removal_time_scale.field[:],
            "pbl_time_scale": cumulus_parameterization_state.output.pbl_time_scale.field[:],
            "lightning_density": cumulus_parameterization_state.output.lightning_density.field[:],
            "evaporation_sublimation_tendency": cumulus_parameterization_state.output.evaporation_sublimation_tendency.field[
                :
            ],
            "convective_precip_flux": cumulus_parameterization_state.output.convective_precip_flux.field[:],
            "t_perturbation": cumulus_parameterization_state.output.t_perturbation.field[:],
            # input/output fields
            "grid_length": cumulus_parameterization_state.input_output.grid_length.field[:],
            "pbl_level": cumulus_parameterization_state.input_output.pbl_level.field[:] + 1,
            "ccn": cumulus_parameterization_state.input_output.ccn.field[:],
            "air_density": cumulus_parameterization_state.input_output.air_density.field[:],
            "omega": cumulus_parameterization_state.input_output.omega.field[:],
            "topography_height_no_negative": cumulus_parameterization_state.input_output.topography_height_no_negative.field[
                :
            ],
            "sensible_heat_flux": cumulus_parameterization_state.input_output.sensible_heat_flux.field[:],
            "latent_heat_flux": cumulus_parameterization_state.input_output.latent_heat_flux.field[:],
            "longitude_degrees": cumulus_parameterization_state.input_output.longitude_degrees.field[:],
            "latitude_degrees": cumulus_parameterization_state.input_output.latitude_degrees.field[:],
            "t_old": cumulus_parameterization_state.input_output.t_old.field[:],
            "vapor_old": cumulus_parameterization_state.input_output.vapor_old.field[:],
            "t_modified_by_advection": cumulus_parameterization_state.input_output.t_modified_by_advection.field[
                :
            ],
            "vapor_modified_by_advection": cumulus_parameterization_state.input_output.vapor_modified_by_advection.field[
                :
            ],
            "geopotential_height_forced": cumulus_parameterization_state.input_output.geopotential_height_forced.field[
                :
            ],
            "p_forced": cumulus_parameterization_state.input_output.p_forced.field[:],
            "p_surface": cumulus_parameterization_state.input_output.p_surface.field[:],
            "t_surface": cumulus_parameterization_state.input_output.t_surface.field[:],
            "u": cumulus_parameterization_state.input_output.u.field[:],
            "v": cumulus_parameterization_state.input_output.v.field[:],
            "w": cumulus_parameterization_state.input_output.w.field[:],
            "mass": cumulus_parameterization_state.input_output.mass.field[:],
            "convective_scale_velocity": cumulus_parameterization_state.input_output.convective_scale_velocity.field[
                :
            ],
            "buoyancy_excess": cumulus_parameterization_state.input_output.buoyancy_excess.field[:],
            "large_scale_ice": cumulus_parameterization_state.input_output.large_scale_ice.field[:],
            "convective_ice": cumulus_parameterization_state.input_output.convective_ice.field[:],
            "large_scale_liquid": cumulus_parameterization_state.input_output.large_scale_liquid.field[:],
            "convective_liquid": cumulus_parameterization_state.input_output.convective_liquid.field[:],
            "large_scale_cloud_fraction": cumulus_parameterization_state.input_output.large_scale_cloud_fraction.field[
                :
            ],
            "convective_cloud_fraction": cumulus_parameterization_state.input_output.convective_cloud_fraction.field[
                :
            ],
            "chemistry_tracers": cumulus_parameterization_state.input_output.chemistry_tracers.field[:],
            "chemistry_tracers_output": chemistry_tracers_output_4d,
        }

        return outputs
