import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020CumulusParameterizationLocals(State):
    t_new: Quantity = dataclasses.field(
        metadata={
            "name": "t_new",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "t_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "t_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_excess: Quantity = dataclasses.field(
        metadata={
            "name": "t_excess",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_new_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "t_new_pbl",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_forced: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_excess: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_excess",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_forced_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_forced_pbl",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    moist_static_energy: Quantity = dataclasses.field(
        metadata={
            "name": "moist_static_energy",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    maximum_updraft_origin_level: Quantity = dataclasses.field(
        metadata={
            "name": "maximum_updraft_origin_level",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    kstabm: Quantity = dataclasses.field(
        metadata={
            "name": "kstabm",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    ocean_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "ocean_fraction",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cap_max: Quantity = dataclasses.field(
        metadata={
            "name": "cap_max",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    error_code_2: Quantity = dataclasses.field(
        metadata={
            "name": "error_code_2",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    error_code_3: Quantity = dataclasses.field(
        metadata={
            "name": "error_code_3",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    cap_max_increment: Quantity = dataclasses.field(
        metadata={
            "name": "cap_max_increment",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height_modified: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_modified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    geopotential_height_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "geopotential_height_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_0: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_0",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_1: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_2: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_2",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_3: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_3",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_0_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_0_pbl",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_1_pbl: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1_pbl",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_work_function_1_fa: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_work_function_1_fa",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cin1: Quantity = dataclasses.field(
        metadata={
            "name": "cin1",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    k_x_modified: Quantity = dataclasses.field(
        metadata={
            "name": "k_x_modified",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    epsilon: Quantity = dataclasses.field(
        metadata={
            "name": "epsilon",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    pbl_time_scale: Quantity = dataclasses.field(
        metadata={
            "name": "pbl_time_scale",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_wetbulb: Quantity = dataclasses.field(
        metadata={
            "name": "t_wetbulb",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vapor_wetbulb: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_wetbulb",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    tau_ecmwf: Quantity = dataclasses.field(
        metadata={
            "name": "tau_ecmwf",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    f_dicycle_modified: Quantity = dataclasses.field(
        metadata={
            "name": "f_dicycle_modified",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    add_buoyancy: Quantity = dataclasses.field(
        metadata={
            "name": "add_buoyancy",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hcdo: Quantity = dataclasses.field(
        metadata={
            "name": "hcdo",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cupclw: Quantity = dataclasses.field(
        metadata={
            "name": "cupclw",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    qrcdo: Quantity = dataclasses.field(
        metadata={
            "name": "qrcdo",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    c1d: Quantity = dataclasses.field(
        metadata={
            "name": "c1d",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evap_bcb: Quantity = dataclasses.field(
        metadata={
            "name": "evap_bcb",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_flux_ensemble: Quantity = dataclasses.field(
        metadata={
            "name": "mass_flux_ensemble",
            "dims": [X_DIM, Y_DIM, "ensemble_members"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precipitation_ensemble: Quantity = dataclasses.field(
        metadata={
            "name": "precipitation_ensemble",
            "dims": [X_DIM, Y_DIM, "ensemble_members"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    scale_dependence_factor_downdraft: Quantity = dataclasses.field(
        metadata={
            "name": "scale_dependence_factor_downdraft",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    random_number: Quantity = dataclasses.field(
        metadata={
            "name": "random_number",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    detrainment_function_updraft: Quantity = dataclasses.field(
        metadata={
            "name": "detrainment_function_updraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    detrainment_function_downdraft: Quantity = dataclasses.field(
        metadata={
            "name": "detrainment_function_downdraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    entrainment_rate_downdraft: Quantity = dataclasses.field(
        metadata={
            "name": "entrainment_rate_downdraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    epsilon_min: Quantity = dataclasses.field(
        metadata={
            "name": "epsilon_min",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    epsilon_max: Quantity = dataclasses.field(
        metadata={
            "name": "epsilon_max",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    arbitrary_numerical_parameter: Quantity = dataclasses.field(
        metadata={
            "name": "arbitrary_numerical_parameter",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_moist_static_energy: Quantity = dataclasses.field(
        metadata={
            "name": "environment_moist_static_energy",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_moist_static_energy_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "environment_moist_static_energy_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_moist_static_energy_forced: Quantity = dataclasses.field(
        metadata={
            "name": "environment_moist_static_energy_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_moist_static_energy_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "environment_moist_static_energy_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_moist_static_energy: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_moist_static_energy",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_moist_static_energy_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_moist_static_energy_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_moist_static_energy_forced: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_moist_static_energy_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_moist_static_energy_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_moist_static_energy_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_mixing_ratio: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_mixing_ratio",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_mixing_ratio_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_mixing_ratio_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_mixing_ratio_forced: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_mixing_ratio_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    environment_saturation_mixing_ratio_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "environment_saturation_mixing_ratio_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "p_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    u_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "u_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    u_c: Quantity = dataclasses.field(
        metadata={
            "name": "u_c",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    v_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "v_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    v_c: Quantity = dataclasses.field(
        metadata={
            "name": "v_c",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    gamma_cloud_levels: Quantity = dataclasses.field(
        metadata={
            "name": "gamma_cloud_levels",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    gamma_cloud_levels_forced: Quantity = dataclasses.field(
        metadata={
            "name": "gamma_cloud_levels_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hydrostatic_air_density: Quantity = dataclasses.field(
        metadata={
            "name": "hydrostatic_air_density",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    partition_liquid_ice: Quantity = dataclasses.field(
        metadata={
            "name": "partition_liquid_ice",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    melting_layer: Quantity = dataclasses.field(
        metadata={
            "name": "melting_layer",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    detrainment_start_level: Quantity = dataclasses.field(
        metadata={
            "name": "detrainment_start_level",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    moist_static_energy_origin_level: Quantity = dataclasses.field(
        metadata={
            "name": "moist_static_energy_origin_level",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    moist_static_energy_origin_level_forced: Quantity = dataclasses.field(
        metadata={
            "name": "moist_static_energy_origin_level_forced",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    start_level: Quantity = dataclasses.field(
        metadata={
            "name": "start_level",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    cloud_moist_static_energy: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_moist_static_energy",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_moist_static_energy_forced: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_moist_static_energy_forced",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_moist_static_energy_forced_transported: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_moist_static_energy_forced_transported",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    negative_buoyancy_depth: Quantity = dataclasses.field(
        metadata={
            "name": "negative_buoyancy_depth",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    frh_lfc: Quantity = dataclasses.field(
        metadata={
            "name": "frh_lfc",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_entrainment_updraft: Quantity = dataclasses.field(
        metadata={
            "name": "mass_entrainment_updraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_detrainment_updraft: Quantity = dataclasses.field(
        metadata={
            "name": "mass_detrainment_updraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_entrainment_u_updraft: Quantity = dataclasses.field(
        metadata={
            "name": "mass_entrainment_u_updraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_detrainment_u_updraft: Quantity = dataclasses.field(
        metadata={
            "name": "mass_detrainment_u_updraft",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
