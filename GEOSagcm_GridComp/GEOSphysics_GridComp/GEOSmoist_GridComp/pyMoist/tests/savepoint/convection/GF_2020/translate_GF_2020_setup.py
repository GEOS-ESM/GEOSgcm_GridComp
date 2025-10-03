from ndsl import Namelist, StencilFactory
from ndsl.stencils.testing.grid import Grid
from ndsl.stencils.testing.savepoint import DataLoader
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.setup import GF2020Setup
from pyMoist.convection.GF_2020.temporaries import GF2020Temporaries
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable


class TranslateGF_2020_setup(TranslateFortranData2Py):
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
            "pressure_deep_convection_cloud_top": {},
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

        self.out_vars = self.in_vars["data_vars"].copy()
        del (
            self.out_vars["latitude"],
            self.out_vars["longitude"],
            self.out_vars["w"],
            self.out_vars["omega"],
            self.out_vars["geopotential_height_interface"],
            self.out_vars["area"],
        )
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
                "pbl_level_local": {},
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
                "p_surface": {},
                "grid_scale_forcing_t": {},
                "grid_scale_forcing_vapor": {},
                "subgrid_scale_forcing_t": {},
                "subgrid_scale_forcing_vapor": {},
                "advective_forcing_t": {},
                "buoyancy_excess": {},
                "t_excess": {},
                "vapor_excess": {},
                "last_ierr": {},
                "fix_out_vapor": {},
                "conprr": {},
                "evap_subl_tendency_cu_param": {},
                "convective_precip_flux_cu_param": {},
                "t_perturbation_cu_param": {},
                "omega_cu_param": {},
                "ccn": {},
                "dtdt_cu_param_shallow": {},
                "dtdt_cu_param_mid": {},
                "dtdt_cu_param_deep": {},
                "dudt_cu_param_shallow": {},
                "dudt_cu_param_mid": {},
                "dudt_cu_param_deep": {},
                "dvdt_cu_param_shallow": {},
                "dvdt_cu_param_mid": {},
                "dvdt_cu_param_deep": {},
                "dvapordt_cu_param_shallow": {},
                "dvapordt_cu_param_mid": {},
                "dvapordt_cu_param_deep": {},
                "dvapordt_cu_param_combined": {},
                "dcloudicedt_cu_param_shallow": {},
                "dcloudicedt_cu_param_mid": {},
                "dcloudicedt_cu_param_deep": {},
                "dnicedt_cu_param_shallow": {},
                "dnicedt_cu_param_mid": {},
                "dnicedt_cu_param_deep": {},
                "dnliquiddt_cu_param_shallow": {},
                "dnliquiddt_cu_param_mid": {},
                "dnliquiddt_cu_param_deep": {},
                "dbuoyancydt_cu_param_shallow": {},
                "dbuoyancydt_cu_param_mid": {},
                "dbuoyancydt_cu_param_deep": {},
                # "dconvectiveicedt_cu_param_shallow": {},
                # "dconvectiveicedt_cu_param_mid": {},
                # "dconvectiveicedt_cu_param_deep": {},
                # "dlargescaleicedt_cu_param_shallow": {},
                # "dlargescaleicedt_cu_param_mid": {},
                # "dlargescaleicedt_cu_param_deep": {},
                # "dconvectiveliquiddt_cu_param_shallow": {},
                # "dconvectiveliquiddt_cu_param_mid": {},
                # "dconvectiveliquiddt_cu_param_deep": {},
                # "dlargescaleliquiddt_cu_param_shallow": {},
                # "dlargescaleliquiddt_cu_param_mid": {},
                # "dlargescaleliquiddt_cu_param_deep": {},
                # "dconvectivecloudfractiondt_cu_param_shallow": {},
                # "dconvectivecloudfractiondt_cu_param_mid": {},
                # "dconvectivecloudfractiondt_cu_param_deep": {},
                # "dlargescalecloudfractiondt_cu_param_shallow": {},
                # "dlargescalecloudfractiondt_cu_param_mid": {},
                # "dlargescalecloudfractiondt_cu_param_deep": {},
                "topography_height_no_negative": {},
                "latitude_degrees": {},
                "longitude_degrees": {},
            }
        )

    def extra_data_load(self, data_loader: DataLoader):
        self.constants = data_loader.load("GF_2020-constants")

    def compute(self, inputs):
        config = GF2020Config(SINGLE_COLUMN_MODE=False, **self.constants)

        state = GF2020State(**inputs)

        temporaries = GF2020Temporaries.make(self.quantity_factory)

        saturation_tables = SaturationVaporPressureTable(self.stencil_factory.backend)

        setup = GF2020Setup(self.stencil_factory, self.quantity_factory, config)

        setup(state, saturation_tables, temporaries)

        from numpy import moveaxis

        inputs.update(
            {
                "edge_height_above_surface": temporaries.edge_height_above_surface.field[:],
                "layer_height_above_surface": temporaries.layer_height_above_surface.field[:],
                "p": temporaries.p.field[:],
                "p_kappa": temporaries.p_kappa.field[:],
                "th": temporaries.th.field[:],
                "mass": temporaries.mass.field[:],
                "vertical_velocity": temporaries.vertical_velocity.field[:],
                "seed_convection": temporaries.seed_convection.field[:],
                "modified_area": temporaries.modified_area.field[:],
                "t_2m_local": temporaries.t_2m_local.field[:],
                "evaporation_local": temporaries.evaporation_local.field[:],
                "sensible_heat_flux_local": temporaries.sensible_heat_flux_local.field[:],
                "topography_height": temporaries.topography_height.field[:],
                "ocean_fraction": temporaries.ocean_fraction.field[:],
                "grid_length": temporaries.grid_length.field[:],
                "pbl_level_local": temporaries.pbl_level_local.field[:],
                "t_local": moveaxis(temporaries.t_local.field[:], 2, 0),
                "p_local": moveaxis(temporaries.p_local.field[:], 2, 0),
                "vapor_local": moveaxis(temporaries.vapor_local.field[:], 2, 0),
                "vapor_current_local": moveaxis(temporaries.vapor_current_local.field[:], 2, 0),
                "u_local": moveaxis(temporaries.u_local.field[:], 2, 0),
                "v_local": moveaxis(temporaries.v_local.field[:], 2, 0),
                "vertical_velocity_local": moveaxis(temporaries.vertical_velocity_local.field[:], 2, 0),
                "layer_height_above_surface_local": moveaxis(
                    temporaries.layer_height_above_surface_local.field[:], 2, 0
                ),
                "edge_height_above_surface_local": moveaxis(
                    temporaries.edge_height_above_surface_local.field[:], 2, 0
                ),
                "mass_local": moveaxis(temporaries.mass_local.field[:], 2, 0),
                "scalar_diffusivity_local": moveaxis(temporaries.scalar_diffusivity_local.field[:], 2, 0),
                "lateral_entrainment_rate_local": moveaxis(
                    temporaries.lateral_entrainment_rate_local.field[:], 2, 0
                ),
                "convective_liquid_local": moveaxis(temporaries.convective_liquid_local.field[:], 2, 0),
                "convective_ice_local": moveaxis(temporaries.convective_ice_local.field[:], 2, 0),
                "convective_cloud_fraction_local": moveaxis(
                    temporaries.convective_cloud_fraction_local.field[:], 2, 0
                ),
                "large_scale_liquid_local": moveaxis(temporaries.large_scale_liquid_local.field[:], 2, 0),
                "large_scale_ice_local": moveaxis(temporaries.large_scale_ice_local.field[:], 2, 0),
                "large_scale_cloud_fraction_local": moveaxis(
                    temporaries.large_scale_cloud_fraction_local.field[:], 2, 0
                ),
                "p_surface": temporaries.p_surface.field[:],
                "grid_scale_forcing_t": moveaxis(temporaries.grid_scale_forcing_t.field[:], 2, 0),
                "grid_scale_forcing_vapor": moveaxis(temporaries.grid_scale_forcing_vapor.field[:], 2, 0),
                "subgrid_scale_forcing_t": moveaxis(temporaries.subgrid_scale_forcing_t.field[:], 2, 0),
                "subgrid_scale_forcing_vapor": moveaxis(
                    temporaries.subgrid_scale_forcing_vapor.field[:], 2, 0
                ),
                "advective_forcing_t": moveaxis(temporaries.advective_forcing_t.field[:], 2, 0),
                "buoyancy_excess": moveaxis(temporaries.buoyancy_excess.field[:], 2, 0),
                "t_excess": temporaries.t_excess.field[:],
                "vapor_excess": temporaries.vapor_excess.field[:],
                "last_ierr": temporaries.last_ierr.field[:],
                "fix_out_vapor": temporaries.fix_out_vapor.field[:],
                "conprr": temporaries.conprr.field[:],
                "evap_subl_tendency_cu_param": temporaries.evap_subl_tendency_cu_param.field[:],
                "convective_precip_flux_cu_param": temporaries.convective_precip_flux_cu_param.field[:],
                "t_perturbation_cu_param": temporaries.t_perturbation_cu_param.field[:],
                "omega_cu_param": temporaries.omega_cu_param.field[:],
                "ccn": temporaries.ccn.field[:],
                "dtdt_cu_param_shallow": temporaries.dtdt_cu_param_shallow.field[:],
                "dtdt_cu_param_mid": temporaries.dtdt_cu_param_mid.field[:],
                "dtdt_cu_param_deep": temporaries.dtdt_cu_param_deep.field[:],
                "dudt_cu_param_shallow": temporaries.dudt_cu_param_shallow.field[:],
                "dudt_cu_param_mid": temporaries.dudt_cu_param_mid.field[:],
                "dudt_cu_param_deep": temporaries.dudt_cu_param_deep.field[:],
                "dvdt_cu_param_shallow": temporaries.dvdt_cu_param_shallow.field[:],
                "dvdt_cu_param_mid": temporaries.dvdt_cu_param_mid.field[:],
                "dvdt_cu_param_deep": temporaries.dvdt_cu_param_deep.field[:],
                "dvapordt_cu_param_shallow": temporaries.dvapordt_cu_param_shallow.field[:],
                "dvapordt_cu_param_mid": temporaries.dvapordt_cu_param_mid.field[:],
                "dvapordt_cu_param_deep": temporaries.dvapordt_cu_param_deep.field[:],
                "dvapordt_cu_param_combined": temporaries.dvapordt_cu_param_combined.field[:],
                "dcloudicedt_cu_param_shallow": temporaries.dcloudicedt_cu_param_shallow.field[:],
                "dcloudicedt_cu_param_mid": temporaries.dcloudicedt_cu_param_mid.field[:],
                "dcloudicedt_cu_param_deep": temporaries.dcloudicedt_cu_param_deep.field[:],
                "dnicedt_cu_param_shallow": temporaries.dnicedt_cu_param_shallow.field[:],
                "dnicedt_cu_param_mid": temporaries.dnicedt_cu_param_mid.field[:],
                "dnicedt_cu_param_deep": temporaries.dnicedt_cu_param_deep.field[:],
                "dnliquiddt_cu_param_shallow": temporaries.dnliquiddt_cu_param_shallow.field[:],
                "dnliquiddt_cu_param_mid": temporaries.dnliquiddt_cu_param_mid.field[:],
                "dnliquiddt_cu_param_deep": temporaries.dnliquiddt_cu_param_deep.field[:],
                "dbuoyancydt_cu_param_shallow": temporaries.dbuoyancydt_cu_param_shallow.field[:],
                "dbuoyancydt_cu_param_mid": temporaries.dbuoyancydt_cu_param_mid.field[:],
                "dbuoyancydt_cu_param_deep": temporaries.dbuoyancydt_cu_param_deep.field[:],
                "dconvectiveicedt_cu_param_shallow": temporaries.dconvectiveicedt_cu_param_shallow.field[:],
                "dconvectiveicedt_cu_param_mid": temporaries.dconvectiveicedt_cu_param_mid.field[:],
                "dconvectiveicedt_cu_param_deep": temporaries.dconvectiveicedt_cu_param_deep.field[:],
                "dlargescaleicedt_cu_param_shallow": temporaries.dlargescaleicedt_cu_param_shallow.field[:],
                "dlargescaleicedt_cu_param_mid": temporaries.dlargescaleicedt_cu_param_mid.field[:],
                "dlargescaleicedt_cu_param_deep": temporaries.dlargescaleicedt_cu_param_deep.field[:],
                "dconvectiveliquiddt_cu_param_shallow": temporaries.dconvectiveliquiddt_cu_param_shallow.field[
                    :
                ],
                "dconvectiveliquiddt_cu_param_mid": temporaries.dconvectiveliquiddt_cu_param_mid.field[:],
                "dconvectiveliquiddt_cu_param_deep": temporaries.dconvectiveliquiddt_cu_param_deep.field[:],
                "dlargescaleliquiddt_cu_param_shallow": temporaries.dlargescaleliquiddt_cu_param_shallow.field[
                    :
                ],
                "dlargescaleliquiddt_cu_param_mid": temporaries.dlargescaleliquiddt_cu_param_mid.field[:],
                "dlargescaleliquiddt_cu_param_deep": temporaries.dlargescaleliquiddt_cu_param_deep.field[:],
                "dconvectivecloudfractiondt_cu_param_shallow": temporaries.dconvectivecloudfractiondt_cu_param_shallow.field[
                    :
                ],
                "dconvectivecloudfractiondt_cu_param_mid": temporaries.dconvectivecloudfractiondt_cu_param_mid.field[
                    :
                ],
                "dconvectivecloudfractiondt_cu_param_deep": temporaries.dconvectivecloudfractiondt_cu_param_deep.field[
                    :
                ],
                "dlargescalecloudfractiondt_cu_param_shallow": temporaries.dlargescalecloudfractiondt_cu_param_shallow.field[
                    :
                ],
                "dlargescalecloudfractiondt_cu_param_mid": temporaries.dlargescalecloudfractiondt_cu_param_mid.field[
                    :
                ],
                "dlargescalecloudfractiondt_cu_param_deep": temporaries.dlargescalecloudfractiondt_cu_param_deep.field[
                    :
                ],
                "topography_height_no_negative": temporaries.topography_height_no_negative.field[:],
                "latitude_degrees": temporaries.latitude_degrees.field[:],
                "longitude_degrees": temporaries.longitude_degrees.field[:],
            }
        )

        return inputs
