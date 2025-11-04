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
            "latitude": {},
            "longitude": {},
            "p_interface": {},
            "t": {},
            "u": {},
            "v": {},
            "w": {},
            "omega": {},
            "t_2m": {},
            "specific_humidity_2m": {},
            "t_surface": {},
            "specific_humidity_surface": {},
            "vapor": {},
            "convective_liquid": {},
            "convective_ice": {},
            "large_scale_liquid": {},
            "large_scale_ice": {},
            "convective_cloud_fraction": {},
            "large_scale_cloud_fraction": {},
            "p_interface_timestep_start": {},
            "t_timestep_start": {},
            "u_timestep_start": {},
            "v_timestep_start": {},
            "vapor_timestep_start": {},
            "geopotential_height_interface": {},
            "geopotential_height_surface": {},
            "area": {},
            "pbl_level": {},
            "convection_fraction": {},
            "surface_type": {},
            "land_fraction": {},
            "scalar_diffusivity": {},
            "buoyancy": {},
            "convective_precipitation_GF": {},
            "convective_precipitation_RAS": {},
            "sensible_heat_flux": {},
            "total_water_flux_deep_convection": {},
            "evaporation": {},
            "convective_condensate_source": {},
            "convective_condensate_grid_mean": {},
            "entrainment_parameter": {},
            "lateral_entrainment_rate": {},
            "lateral_entrainment_rate_shallow": {},
            "lateral_entrainment_rate_mid": {},
            "lateral_entrainment_rate_deep": {},
            "updraft_area_fraction": {},
            "updraft_vertical_velocity": {},
            "dtdt_shortwave": {},
            "dtdt_longwave": {},
            "dspecific_humiditydt_pbl": {},
            "dtdt_pbl": {},
            "dtdt_from_dynamics": {},
            "dvapordt_from_dynamics": {},
            "sigma_mid": {},
            "sigma_deep": {},
            "total_precipitable_water_initial": {},
            "saturation_total_precipitable_water_initial": {},
            "dvapordt_deep_convection": {},
            "dtdt_deep_convection": {},
            "dudt_deep_convection": {},
            "dvdt_deep_convection": {},
            "pressure_shallow_convective_cloud_top": {},
            "pressure_mid_convective_cloud_top": {},
            "pressure_deep_convective_cloud_top": {},
            "mass_flux_shalow": {},
            "mass_flux_mid": {},
            "mass_flux_deep_updraft": {},
            "mass_flux_deep_updraft_interface": {},
            "mass_flux_deep_updraft_detrained": {},
            "mass_flux_deep_downdraft": {},
            "mass_flux_cloud_base": {},
            "mass_flux_cloud_base_shallow": {},
            "mass_flux_cloud_base_mid": {},
            "mass_flux_cloud_base_deep": {},
            "convection_code_shallow": {},
            "convection_code_mid": {},
            "convection_code_deep": {},
            "cloud_work_function_0": {},
            "cloud_work_function_1": {},
            "cloud_work_function_2": {},
            "cloud_work_function_3": {},
            "cloud_work_function_1_pbl": {},
            "cloud_work_function_1_cin": {},
            "pbl_time_scale": {},
            "cape_removal_time_scale": {},
            "lighting_density": {},
            "convection_tracer": {},
        }

        self.out_vars = {}
        self.out_vars.update(
            {
                "edge_height_above_surface": {},
                "layer_height_above_surface": {},
                "p": {},
                "p_kappa": {},
                "th": {},
                "mass": {},
                "vertical_velocity": {},
                "seed_convection": {},
                "modified_area": {},
                "t_2m_local": {},
                "evaporation_local": {},
                "sensible_heat_flux_local": {},
                "topography_height": {},
                "ocean_fraction": {},
                "grid_length": {},
                "pbl_level_cu_param_input": {},
                "t_local": {},
                "p_local": {},
                "vapor_local": {},
                "vapor_current_local": {},
                "u_local": {},
                "v_local": {},
                "vertical_velocity_local": {},
                "layer_height_above_surface_local": {},
                "edge_height_above_surface_local": {},
                "mass_local": {},
                "scalar_diffusivity_local": {},
                "lateral_entrainment_rate_local": {},
                "convective_liquid_local": {},
                "convective_ice_local": {},
                "convective_cloud_fraction_local": {},
                "large_scale_liquid_local": {},
                "large_scale_ice_local": {},
                "large_scale_cloud_fraction_local": {},
                "t_surface_cu_param_input": {},
                "p_surface_cu_param_input": {},
                "grid_scale_forcing_t_cu_param_input": {},
                "grid_scale_forcing_vapor_cu_param_input": {},
                "subgrid_scale_forcing_t_cu_param_input": {},
                "subgrid_scale_forcing_vapor_cu_param_input": {},
                "advective_forcing_t_cu_param_input": {},
                "buoyancy_excess": {},
                "t_excess_cu_param_input": {},
                "vapor_excess_cu_param_input": {},
                "last_ierr": {},
                "fix_out_vapor": {},
                "conprr": {},
                "evap_subl_tendency": {},
                "convective_precip_flux": {},
                "t_perturbation": {},
                "omega_cu_param_input": {},
                "ccn": {},
                "topography_height_no_negative": {},
                "latitude_degrees": {},
                "longitude_degrees": {},
                "geopotential_height_cu_param_input": {},
                "p_cu_param_input": {},
                "t_cu_param_input": {},
                "vapor_timestep_start_cu_param_input": {},
                "vapor_current_cu_param_input": {},
                "air_density_cu_param_input": {},
                "u_cu_param_input": {},
                "v_cu_param_input": {},
                "w_cu_param_input": {},
                "mass_cu_param_input": {},
                "t_modified_by_advection": {},
                "vapor_modified_by_advection": {},
                # "convective_liquid_cu_param_input": {}, # all nan first timestep
                # "convective_ice_cu_param_input": {}, # all nan first timestep
                # "convective_cloud_fraction_cu_param_input": {}, # all nan first timestep
                # "large_scale_liquid_cu_param_input": {}, # all nan first timestep
                # "large_scale_ice_cu_param_input": {}, # all nan first timestep
                # "large_scale_cloud_fraction_cu_param_input": {}, # all nan first timestep
                "pbl_height_cu_param_input": {},
                "sensible_heat_flux_cu_param_inputs": {},
                "latent_heat_flux_cu_param_inputs": {},
                "convective_scale_velosity_cu_param_input": {},
                "ocean_fraction_cu_param_input": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF2020-constants")

    def compute(self, inputs):
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)

        state = GF2020State(**inputs)

        locals = GF2020Locals.zeros(
            self.quantity_factory,
            data_dimensions={
                "plumes": 3,
            },
        )

        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        setup = GF2020Setup(self.stencil_factory, self.quantity_factory, config)

        setup(state, locals, saturation_tables)

        import numpy as np

        # top rows are not computed in Fortran, retains initalized value (nan)
        # Python initalizes with zero, fill top row with nan for test passage
        locals.cumulus_parameterization_input.p_mb.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.t.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.vapor_timestep_start.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.vapor_current.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.u.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.v.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.w.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.buoyancy_excess.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.geopotential_height.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.air_density.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.mass.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.t_modified_by_advection.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.vapor_modified_by_advection.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.convective_liquid.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.convective_ice.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.convective_cloud_fraction.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.large_scale_liquid.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.large_scale_ice.field[:, :, -1] = np.nan
        locals.cumulus_parameterization_input.large_scale_cloud_fraction.field[:, :, -1] = np.nan
        print(f"README: {locals.cumulus_parameterization_input.pbl_level.field[:].shape}")
        inputs.update(
            {
                "edge_height_above_surface": locals.derived_state.edge_height_above_surface.field[:],
                "layer_height_above_surface": locals.derived_state.layer_height_above_surface.field[:],
                "p": locals.derived_state.p.field[:],
                "p_kappa": locals.derived_state.p_kappa.field[:],
                "th": locals.derived_state.th.field[:],
                "mass": locals.derived_state.mass.field[:],
                "vertical_velocity": locals.derived_state.vertical_velocity.field[:],
                "seed_convection": locals.derived_state.seed_convection.field[:],
                "modified_area": locals.derived_state.modified_area.field[:],
                "t_2m_local": locals.local_copy.t_2m.field[:],
                "evaporation_local": locals.local_copy.evaporation.field[:],
                "sensible_heat_flux_local": locals.local_copy.sensible_heat_flux.field[:],
                "topography_height": locals.derived_state.topography_height.field[:],
                "ocean_fraction": locals.cumulus_parameterization_input.ocean_fraction.field[:],
                "grid_length": locals.cumulus_parameterization_input.grid_length.field[:],
                "pbl_level_cu_param_input": locals.cumulus_parameterization_input.pbl_level.field[:] + 1,
                "t_local": np.moveaxis(locals.local_copy.t.field[:], 2, 0),
                "p_local": np.moveaxis(locals.local_copy.p.field[:], 2, 0),
                "vapor_local": np.moveaxis(locals.local_copy.vapor.field[:], 2, 0),
                "vapor_current_local": np.moveaxis(locals.local_copy.vapor_current.field[:], 2, 0),
                "u_local": np.moveaxis(locals.local_copy.u.field[:], 2, 0),
                "v_local": np.moveaxis(locals.local_copy.v.field[:], 2, 0),
                "vertical_velocity_local": np.moveaxis(locals.local_copy.vertical_velocity.field[:], 2, 0),
                "layer_height_above_surface_local": np.moveaxis(
                    locals.local_copy.layer_height_above_surface.field[:], 2, 0
                ),
                "edge_height_above_surface_local": np.moveaxis(
                    locals.local_copy.edge_height_above_surface.field[:], 2, 0
                ),
                "mass_local": np.moveaxis(locals.local_copy.mass.field[:], 2, 0),
                "scalar_diffusivity_local": np.moveaxis(locals.local_copy.scalar_diffusivity.field[:], 2, 0),
                "lateral_entrainment_rate_local": np.moveaxis(
                    locals.local_copy.lateral_entrainment_rate.field[:], 2, 0
                ),
                "convective_liquid_local": np.moveaxis(locals.local_copy.convective_liquid.field[:], 2, 0),
                "convective_ice_local": np.moveaxis(locals.local_copy.convective_ice.field[:], 2, 0),
                "convective_cloud_fraction_local": np.moveaxis(
                    locals.local_copy.convective_cloud_fraction.field[:], 2, 0
                ),
                "large_scale_liquid_local": np.moveaxis(locals.local_copy.large_scale_liquid.field[:], 2, 0),
                "large_scale_ice_local": np.moveaxis(locals.local_copy.large_scale_ice.field[:], 2, 0),
                "large_scale_cloud_fraction_local": np.moveaxis(
                    locals.local_copy.large_scale_cloud_fraction.field[:], 2, 0
                ),
                "t_surface_cu_param_input": locals.cumulus_parameterization_input.t_surface.field[:],
                "p_surface_cu_param_input": locals.cumulus_parameterization_input.p_surface.field[:],
                "grid_scale_forcing_t_cu_param_input": np.moveaxis(
                    locals.cumulus_parameterization_input.grid_scale_forcing_t.field[:], 2, 0
                ),
                "grid_scale_forcing_vapor_cu_param_input": np.moveaxis(
                    locals.cumulus_parameterization_input.grid_scale_forcing_vapor.field[:], 2, 0
                ),
                "subgrid_scale_forcing_t_cu_param_input": np.moveaxis(
                    locals.cumulus_parameterization_input.subgrid_scale_forcing_t.field[:], 2, 0
                ),
                "subgrid_scale_forcing_vapor_cu_param_input": np.moveaxis(
                    locals.cumulus_parameterization_input.subgrid_scale_forcing_vapor.field[:], 2, 0
                ),
                "advective_forcing_t_cu_param_input": np.moveaxis(
                    locals.cumulus_parameterization_input.advective_forcing_t.field[:], 2, 0
                ),
                "buoyancy_excess": locals.cumulus_parameterization_input.buoyancy_excess.field[:],
                "t_excess_cu_param_input": locals.cumulus_parameterization_input.t_excess.field[:],
                "vapor_excess_cu_param_input": locals.cumulus_parameterization_input.vapor_excess.field[:],
                "last_ierr": locals.miscelaneous_diagnostic.last_ierr.field[:],
                "fix_out_vapor": locals.miscelaneous_diagnostic.fix_out_vapor.field[:],
                "conprr": locals.miscelaneous_diagnostic.conprr.field[:],
                "evap_subl_tendency": locals.cumulus_parameterization_output.evap_subl_tendency.field[:],
                "convective_precip_flux": locals.cumulus_parameterization_output.convective_precip_flux.field[
                    :
                ],
                "t_perturbation": locals.cumulus_parameterization_input.t_perturbation.field[:],
                "omega_cu_param_input": locals.cumulus_parameterization_input.omega.field[:],
                "ccn": locals.cumulus_parameterization_input.ccn.field[:],
                "topography_height_no_negative": locals.cumulus_parameterization_input.topography_height.field[
                    :
                ],
                "latitude_degrees": locals.cumulus_parameterization_input.latitude_degrees.field[:],
                "longitude_degrees": locals.cumulus_parameterization_input.longitude_degrees.field[:],
                "geopotential_height_cu_param_input": locals.cumulus_parameterization_input.geopotential_height.field[
                    :
                ],
                "p_cu_param_input": locals.cumulus_parameterization_input.p_mb.field[:],
                "t_cu_param_input": locals.cumulus_parameterization_input.t.field[:],
                "vapor_timestep_start_cu_param_input": locals.cumulus_parameterization_input.vapor_timestep_start.field[
                    :
                ],
                "vapor_current_cu_param_input": locals.cumulus_parameterization_input.vapor_current.field[:],
                "air_density_cu_param_input": locals.cumulus_parameterization_input.air_density.field[:],
                "u_cu_param_input": locals.cumulus_parameterization_input.u.field[:],
                "v_cu_param_input": locals.cumulus_parameterization_input.v.field[:],
                "w_cu_param_input": locals.cumulus_parameterization_input.w.field[:],
                "mass_cu_param_input": locals.cumulus_parameterization_input.mass.field[:],
                "t_modified_by_advection": locals.cumulus_parameterization_input.t_modified_by_advection.field[
                    :
                ],
                "vapor_modified_by_advection": locals.cumulus_parameterization_input.vapor_modified_by_advection.field[
                    :
                ],
                "convective_liquid_cu_param_input": locals.cumulus_parameterization_input.convective_liquid.field[
                    :
                ],
                "convective_ice_cu_param_input": locals.cumulus_parameterization_input.convective_ice.field[
                    :
                ],
                "convective_cloud_fraction_cu_param_input": locals.cumulus_parameterization_input.convective_cloud_fraction.field[
                    :
                ],
                "large_scale_liquid_cu_param_input": locals.cumulus_parameterization_input.large_scale_liquid.field[
                    :
                ],
                "large_scale_ice_cu_param_input": locals.cumulus_parameterization_input.large_scale_ice.field[
                    :
                ],
                "large_scale_cloud_fraction_cu_param_input": locals.cumulus_parameterization_input.large_scale_cloud_fraction.field[
                    :
                ],
                "pbl_height_cu_param_input": locals.cumulus_parameterization_input.pbl_height.field[:],
                "sensible_heat_flux_cu_param_inputs": locals.cumulus_parameterization_input.sensible_heat_flux.field[
                    :
                ],
                "latent_heat_flux_cu_param_inputs": locals.cumulus_parameterization_input.latent_heat_flux.field[
                    :
                ],
                "convective_scale_velosity_cu_param_input": locals.cumulus_parameterization_input.convective_scale_velosity.field[
                    :
                ],
                "ocean_fraction_cu_param_input": locals.cumulus_parameterization_input.ocean_fraction.field[
                    :
                ],
            }
        )

        return inputs
