import dataclasses

from ndsl import Quantity, State
from ndsl.constants import I_XIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Float


@dataclasses.dataclass
class GF2020State(State):
    latitude: Quantity = dataclasses.field(
        metadata={
            "name": "latitude",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    longitude: Quantity = dataclasses.field(
        metadata={
            "name": "longitude",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_interface: Quantity = dataclasses.field(
        metadata={
            "name": "p_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t: Quantity = dataclasses.field(
        metadata={
            "name": "t",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    u: Quantity = dataclasses.field(
        metadata={
            "name": "u",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    v: Quantity = dataclasses.field(
        metadata={
            "name": "v",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    w: Quantity = dataclasses.field(
        metadata={
            "name": "w",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    omega: Quantity = dataclasses.field(
        metadata={
            "name": "omega",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_2m: Quantity = dataclasses.field(
        metadata={
            "name": "t_2m",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    specific_humidity_2m: Quantity = dataclasses.field(
        metadata={
            "name": "specific_humidity_2m",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_surface: Quantity = dataclasses.field(
        metadata={
            "name": "t_surface",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    specific_humidity_surface: Quantity = dataclasses.field(
        metadata={
            "name": "specific_humidity_surface",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor: Quantity = dataclasses.field(
        metadata={
            "name": "vapor",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_liquid: Quantity = dataclasses.field(
        metadata={
            "name": "convective_liquid",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_ice: Quantity = dataclasses.field(
        metadata={
            "name": "convective_ice",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_cloud_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "convective_cloud_fraction",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_liquid: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_liquid",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_ice: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_ice",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    large_scale_cloud_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "large_scale_cloud_fraction",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice_fraction_in_convective_tower: Quantity = dataclasses.field(
        metadata={
            "name": "ice_fraction_in_convective_tower",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_interface_timestep_start: Quantity = dataclasses.field(
        metadata={
            "name": "p_interface_timestep_start",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_timestep_start: Quantity = dataclasses.field(
        metadata={
            "name": "t_timestep_start",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    u_timestep_start: Quantity = dataclasses.field(
        metadata={
            "name": "u_timestep_start",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    v_timestep_start: Quantity = dataclasses.field(
        metadata={
            "name": "v_timestep_start",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_timestep_start: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_timestep_start",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height_interface: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height_surface: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_surface",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    area: Quantity = dataclasses.field(
        metadata={
            "name": "area",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pbl_level: Quantity = dataclasses.field(
        metadata={
            "name": "pbl_level",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convection_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "convection_fraction",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    surface_type: Quantity = dataclasses.field(
        metadata={
            "name": "surface_type",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    seed_convection: Quantity = dataclasses.field(
        metadata={
            "name": "seed_convection",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    land_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "land_fraction",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    scalar_diffusivity: Quantity = dataclasses.field(
        metadata={
            "name": "scalar_diffusivity",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    buoyancy: Quantity = dataclasses.field(
        metadata={
            "name": "buoyancy",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_precipitation_GF: Quantity = dataclasses.field(
        metadata={
            "name": "convective_precipitation_GF",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_precipitation_RAS: Quantity = dataclasses.field(
        metadata={
            "name": "convective_precipitation_RAS",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_rainwater_source: Quantity = dataclasses.field(
        metadata={
            "name": "convective_rainwater_source",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    sensible_heat_flux: Quantity = dataclasses.field(
        metadata={
            "name": "sensible_heat_flux",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    total_water_flux_deep_convection_interface: Quantity = dataclasses.field(
        metadata={
            "name": "total_water_flux_deep_convection_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    sublimation_of_convective_precipitation: Quantity = dataclasses.field(
        metadata={
            "name": "sublimation_of_convective_precipitation",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation_of_convective_precipitation: Quantity = dataclasses.field(
        metadata={
            "name": "evaporation_of_convective_precipitation",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice_precip_flux_interface: Quantity = dataclasses.field(
        metadata={
            "name": "ice_precip_flux_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liquid_precip_flux_interface: Quantity = dataclasses.field(
        metadata={
            "name": "liquid_precip_flux_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation: Quantity = dataclasses.field(
        metadata={
            "name": "evaporation",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_condensate_source: Quantity = dataclasses.field(
        metadata={
            "name": "convective_condensate_source",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_condensate_grid_mean: Quantity = dataclasses.field(
        metadata={
            "name": "convective_condensate_grid_mean",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    entrainment_parameter: Quantity = dataclasses.field(
        metadata={
            "name": "entrainment_parameter",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lateral_entrainment_rate: Quantity = dataclasses.field(
        metadata={
            "name": "lateral_entrainment_rate",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lateral_entrainment_rate_shallow: Quantity = dataclasses.field(
        metadata={
            "name": "lateral_entrainment_rate_shallow",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lateral_entrainment_rate_mid: Quantity = dataclasses.field(
        metadata={
            "name": "lateral_entrainment_rate_mid",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lateral_entrainment_rate_deep: Quantity = dataclasses.field(
        metadata={
            "name": "lateral_entrainment_rate_deep",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    updraft_areal_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "updraft_areal_fraction",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    updraft_vertical_velocity: Quantity = dataclasses.field(
        metadata={
            "name": "updraft_vertical_velocity",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt_shortwave: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt_shortwave",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt_longwave: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt_longwave",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dspecific_humiditydt_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "dspecific_humiditydt_pbl",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt_pbl",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt_from_dynamics: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt_from_dynamics",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvapordt_from_dynamics: Quantity = dataclasses.field(
        metadata={
            "name": "dvapordt_from_dynamics",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    sigma_mid: Quantity = dataclasses.field(
        metadata={
            "name": "sigma_mid",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    sigma_deep: Quantity = dataclasses.field(
        metadata={
            "name": "sigma_deep",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    total_precipitable_water_initial: Quantity = dataclasses.field(
        metadata={
            "name": "total_precipitable_water_initial",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    saturation_total_precipitable_water_initial: Quantity = dataclasses.field(
        metadata={
            "name": "saturation_total_precipitable_water_initial",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvapordt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dvapordt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dudt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dudt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvdt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dvdt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dliquiddt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dliquiddt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dicedt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dicedt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dcloudfractiondt_deep_convection: Quantity = dataclasses.field(
        metadata={
            "name": "dcloudfractiondt_deep_convection",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pressure_shallow_convective_cloud_top: Quantity = dataclasses.field(
        metadata={
            "name": "pressure_shallow_convective_cloud_top",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pressure_mid_convective_cloud_top: Quantity = dataclasses.field(
        metadata={
            "name": "pressure_mid_convective_cloud_top",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pressure_deep_convective_cloud_top: Quantity = dataclasses.field(
        metadata={
            "name": "pressure_deep_convective_cloud_top",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_shallow: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_shallow",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_mid: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_mid",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_deep_updraft: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_deep_updraft",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_deep_updraft_interface: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_deep_updraft_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_deep_updraft_detrained: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_deep_updraft_detrained",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_deep_downdraft: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_deep_downdraft",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_cloud_base: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_cloud_base",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_cloud_base_shallow: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_cloud_base_shallow",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_cloud_base_mid: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_cloud_base_mid",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_cloud_base_deep: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_cloud_base_deep",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    total_cumulative_mass_flux_interface: Quantity = dataclasses.field(
        metadata={
            "name": "total_cumulative_mass_flux_interface",
            "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    total_detraining_mass_flux: Quantity = dataclasses.field(
        metadata={
            "name": "total_detraining_mass_flux",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convection_code_shallow: Quantity = dataclasses.field(
        metadata={
            "name": "convection_code_shallow",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convection_code_mid: Quantity = dataclasses.field(
        metadata={
            "name": "convection_code_mid",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convection_code_deep: Quantity = dataclasses.field(
        metadata={
            "name": "convection_code_deep",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_workfunction_0: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_0",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_workfunction_1: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_workfunction_2: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_2",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_workfunction_3: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_3",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_workfunction_1_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1_pbl",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_workfunction_1_cin: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1_cin",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pbl_time_scale: Quantity = dataclasses.field(
        metadata={
            "name": "pbl_time_scale",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cape_removal_time_scale: Quantity = dataclasses.field(
        metadata={
            "name": "cape_removal_time_scale",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lightning_density: Quantity = dataclasses.field(
        metadata={
            "name": "lighting_density",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convection_tracer: Quantity = dataclasses.field(
        metadata={
            "name": "convection_tracer",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
