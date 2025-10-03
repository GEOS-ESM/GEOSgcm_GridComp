from dataclasses import dataclass

import numpy.typing as npt

from ndsl import Quantity


@dataclass
class GF2020State:
    """
    Units:
        latitude
        longitude
        p_interface
        t
        u
        v
        w
        omega
        t_2m
        specific_humidity_2m
        t_surface
        specific_humidity_surface
        vapor
        convective_liquid
        convective_ice
        convective_cloud_fraction
        large_scale_liquid
        large_scale_ice
        large_scale_cloud_fraction
        p_interface_timestep_start
        t_timestep_start
        u_timestep_start
        v_timestep_start
        vapor_timestep_start
        geopotential_height_interface
        geopotential_height_surface
        area
        pbl_level
        convection_fraction
        surface_type
        land_fraction
        scalar_diffusivity
        buoyancy
        convective_precipitation_GF
        convective_precipitation_RAS
        sensible_heat_flux
        total_water_flux_deep_convection
        evaporation
        convective_condensate_source
        convective_condensate_grid_mean
        entrainment_parameter
        lateral_entrainment_rate
        lateral_entrainment_rate_shallow
        lateral_entrainment_rate_mid
        lateral_entrainment_rate_deep
        updraft_area_fraction
        updraft_vertical_velocity
        dtdt_shortwave
        dtdt_longwave
        dspecific_humiditydt_pbl
        dtdt_pbl
        dtdt_from_dynamics
        dvapordt_from_dynamics
        sigma_mid
        sigma_deep
        total_precipitable_water_initial
        saturation_total_precipitable_water_initial
        dvapordt_deep_convection
        dtdt_deep_convection
        dudt_deep_convection
        dvdt_deep_convection
        pressure_shallow_convective_cloud_top
        pressure_mid_convective_cloud_top
        pressure_deep_convection_cloud_top
        mass_flux_shalow
        mass_flux_mid
        mass_flux_deep_updraft
        mass_flux_deep_updraft_interface
        mass_flux_deep_updraft_detrained
        mass_flux_deep_downdraft
        mass_flux_cloud_base
        mass_flux_cloud_base_shallow
        mass_flux_cloud_base_mid
        mass_flux_cloud_base_deep
        convection_code_shallow
        convection_code_mid
        convection_code_deep
        cloud_work_function_0
        cloud_work_function_1
        cloud_work_function_2
        cloud_work_function_3
        cloud_work_function_1_pbl
        cloud_work_function_1_cin
        pbl_time_scale
        cape_removal_time_scale
        lighting_density
        convection_tracer
    """

    latitude: Quantity | npt.NDArray
    longitude: Quantity | npt.NDArray
    p_interface: Quantity | npt.NDArray
    t: Quantity | npt.NDArray
    u: Quantity | npt.NDArray
    v: Quantity | npt.NDArray
    w: Quantity | npt.NDArray
    omega: Quantity | npt.NDArray
    t_2m: Quantity | npt.NDArray
    specific_humidity_2m: Quantity | npt.NDArray
    t_surface: Quantity | npt.NDArray
    specific_humidity_surface: Quantity | npt.NDArray
    vapor: Quantity | npt.NDArray
    convective_liquid: Quantity | npt.NDArray
    convective_ice: Quantity | npt.NDArray
    convective_cloud_fraction: Quantity | npt.NDArray
    large_scale_liquid: Quantity | npt.NDArray
    large_scale_ice: Quantity | npt.NDArray
    large_scale_cloud_fraction: Quantity | npt.NDArray
    p_interface_timestep_start: Quantity | npt.NDArray
    t_timestep_start: Quantity | npt.NDArray
    u_timestep_start: Quantity | npt.NDArray
    v_timestep_start: Quantity | npt.NDArray
    vapor_timestep_start: Quantity | npt.NDArray
    geopotential_height_interface: Quantity | npt.NDArray
    geopotential_height_surface: Quantity | npt.NDArray
    area: Quantity | npt.NDArray
    pbl_level: Quantity | npt.NDArray
    convection_fraction: Quantity | npt.NDArray
    surface_type: Quantity | npt.NDArray
    land_fraction: Quantity | npt.NDArray
    scalar_diffusivity: Quantity | npt.NDArray
    buoyancy: Quantity | npt.NDArray
    convective_precipitation_GF: Quantity | npt.NDArray
    convective_precipitation_RAS: Quantity | npt.NDArray
    sensible_heat_flux: Quantity | npt.NDArray
    total_water_flux_deep_convection: Quantity | npt.NDArray
    evaporation: Quantity | npt.NDArray
    convective_condensate_source: Quantity | npt.NDArray
    convective_condensate_grid_mean: Quantity | npt.NDArray
    entrainment_parameter: Quantity | npt.NDArray
    lateral_entrainment_rate: Quantity | npt.NDArray
    lateral_entrainment_rate_shallow: Quantity | npt.NDArray
    lateral_entrainment_rate_mid: Quantity | npt.NDArray
    lateral_entrainment_rate_deep: Quantity | npt.NDArray
    updraft_area_fraction: Quantity | npt.NDArray
    updraft_vertical_velocity: Quantity | npt.NDArray
    dtdt_shortwave: Quantity | npt.NDArray
    dtdt_longwave: Quantity | npt.NDArray
    dspecific_humiditydt_pbl: Quantity | npt.NDArray
    dtdt_pbl: Quantity | npt.NDArray
    dtdt_from_dynamics: Quantity | npt.NDArray
    dvapordt_from_dynamics: Quantity | npt.NDArray
    sigma_mid: Quantity | npt.NDArray
    sigma_deep: Quantity | npt.NDArray
    total_precipitable_water_initial: Quantity | npt.NDArray
    saturation_total_precipitable_water_initial: Quantity | npt.NDArray
    dvapordt_deep_convection: Quantity | npt.NDArray
    dtdt_deep_convection: Quantity | npt.NDArray
    dudt_deep_convection: Quantity | npt.NDArray
    dvdt_deep_convection: Quantity | npt.NDArray
    pressure_shallow_convective_cloud_top: Quantity | npt.NDArray
    pressure_mid_convective_cloud_top: Quantity | npt.NDArray
    pressure_deep_convection_cloud_top: Quantity | npt.NDArray
    mass_flux_shalow: Quantity | npt.NDArray
    mass_flux_mid: Quantity | npt.NDArray
    mass_flux_deep_updraft: Quantity | npt.NDArray
    mass_flux_deep_updraft_interface: Quantity | npt.NDArray
    mass_flux_deep_updraft_detrained: Quantity | npt.NDArray
    mass_flux_deep_downdraft: Quantity | npt.NDArray
    mass_flux_cloud_base: Quantity | npt.NDArray
    mass_flux_cloud_base_shallow: Quantity | npt.NDArray
    mass_flux_cloud_base_mid: Quantity | npt.NDArray
    mass_flux_cloud_base_deep: Quantity | npt.NDArray
    convection_code_shallow: Quantity | npt.NDArray
    convection_code_mid: Quantity | npt.NDArray
    convection_code_deep: Quantity | npt.NDArray
    cloud_work_function_0: Quantity | npt.NDArray
    cloud_work_function_1: Quantity | npt.NDArray
    cloud_work_function_2: Quantity | npt.NDArray
    cloud_work_function_3: Quantity | npt.NDArray
    cloud_work_function_1_pbl: Quantity | npt.NDArray
    cloud_work_function_1_cin: Quantity | npt.NDArray
    pbl_time_scale: Quantity | npt.NDArray
    cape_removal_time_scale: Quantity | npt.NDArray
    lighting_density: Quantity | npt.NDArray
    convection_tracer: Quantity | npt.NDArray
