import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020CumulusParameterizationState(State):
    @dataclasses.dataclass
    class Input:
        t_excess: Quantity = dataclasses.field(
            metadata={
                "name": "t_excess",
                "dims": [X_DIM, Y_DIM],
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

        grid_scale_forcing_t: Quantity = dataclasses.field(
            metadata={
                "name": "grid_scale_forcing_t",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        grid_scale_forcing_vapor: Quantity = dataclasses.field(
            metadata={
                "name": "grid_scale_forcing_vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        subgrid_scale_forcing_t: Quantity = dataclasses.field(
            metadata={
                "name": "subgrid_scale_forcing_t",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        subgrid_scale_forcing_vapor: Quantity = dataclasses.field(
            metadata={
                "name": "subgrid_scale_forcing_vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        seed_convection: Quantity = dataclasses.field(
            metadata={
                "name": "seed_convection",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        saturation_water_vapor: Quantity = dataclasses.field(
            metadata={
                "name": "saturation_water_vapor",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
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
        convection_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "convection_fraction",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        surface_type: Quantity = dataclasses.field(
            metadata={
                "name": "surface_type",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        lateral_entrainment_rate: Quantity = dataclasses.field(
            metadata={
                "name": "lateral_entrainment_rate",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        last_error_code: Quantity = dataclasses.field(
            metadata={
                "name": "last_error_code",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Output:
        t: Quantity = dataclasses.field(
            metadata={
                "name": "t",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloudicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dcloudicedt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt: Quantity = dataclasses.field(
            metadata={
                "name": "dudt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnliquiddt: Quantity = dataclasses.field(
            metadata={
                "name": "dnliquiddt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dnicedt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dbuoyancydt: Quantity = dataclasses.field(
            metadata={
                "name": "dbuoyancydt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveicedt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleicedt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveliquiddt: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveliquiddt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleliquiddt: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleliquiddt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectivecloudfractiondt: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectivecloudfractiondt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescalecloudfractiondt: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescalecloudfractiondt",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        error_code: Quantity = dataclasses.field(
            metadata={
                "name": "error_code",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Int,
            }
        )
        downdraft_origin_level: Quantity = dataclasses.field(
            metadata={
                "name": "downdraft_origin_level",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        lcl_level: Quantity = dataclasses.field(
            metadata={
                "name": "lcl_level",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        updraft_origin_level: Quantity = dataclasses.field(
            metadata={
                "name": "updraft_origin_level",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        updraft_lfc_level: Quantity = dataclasses.field(
            metadata={
                "name": "updraft_lfc_level",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        cloud_top: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_top",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        kstabi: Quantity = dataclasses.field(
            metadata={
                "name": "kstabi",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        kstabm: Quantity = dataclasses.field(
            metadata={
                "name": "kstabm",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        precip: Quantity = dataclasses.field(
            metadata={
                "name": "precip",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        cloud_base_mass_flux: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_base_mass_flux",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        epsilon: Quantity = dataclasses.field(
            metadata={
                "name": "epsilon",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        pwav: Quantity = dataclasses.field(
            metadata={
                "name": "pwav",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        scale_dependence_factor: Quantity = dataclasses.field(
            metadata={
                "name": "scale_dependence_factor",
                "dims": [X_DIM, Y_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_cloud_levels_forced: Quantity = dataclasses.field(
            metadata={
                "name": "p_cloud_levels_forced",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        entrainment_rate: Quantity = dataclasses.field(
            metadata={
                "name": "entrainment_rate",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        updraft_mass_entrainment: Quantity = dataclasses.field(
            metadata={
                "name": "updraft_mass_entrainment",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        updraft_mass_detrainment: Quantity = dataclasses.field(
            metadata={
                "name": "updraft_mass_detrainment",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        downdraft_mass_entrainment: Quantity = dataclasses.field(
            metadata={
                "name": "downdraft_mass_entrainment",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        downdraft_mass_detrainment: Quantity = dataclasses.field(
            metadata={
                "name": "downdraft_mass_detrainment",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        z_updraft: Quantity = dataclasses.field(
            metadata={
                "name": "z_updraft",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        z_downdraft: Quantity = dataclasses.field(
            metadata={
                "name": "z_downdraft",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_updraft: Quantity = dataclasses.field(
            metadata={
                "name": "p_updraft",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_downdraft: Quantity = dataclasses.field(
            metadata={
                "name": "p_downdraft",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        cloud_liquid_after_rain: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_liquid_after_rain",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_updraft: Quantity = dataclasses.field(
            metadata={
                "name": "t_updraft",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_cloud_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "convective_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
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
        cloud_work_function_1_pbl: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_work_function_1_pbl",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        cloud_work_function_1_cin: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_work_function_1_cin",
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
        cape_removal_time_scale: Quantity = dataclasses.field(
            metadata={
                "name": "cape_removal_time_scale",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        lightning_density: Quantity = dataclasses.field(
            metadata={
                "name": "lightning_density",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        evap_subl_tendency: Quantity = dataclasses.field(
            metadata={
                "name": "evap_subl_tendency",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "convective_precip_flux",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_perturbation: Quantity = dataclasses.field(
            metadata={
                "name": "t_perturbation",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class InputOutput:
        grid_length: Quantity = dataclasses.field(
            metadata={
                "name": "grid_length",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        pbl_level: Quantity = dataclasses.field(
            metadata={
                "name": "pbl_level",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ccn: Quantity = dataclasses.field(
            metadata={
                "name": "ccn",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        air_density: Quantity = dataclasses.field(
            metadata={
                "name": "air_density",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        omega: Quantity = dataclasses.field(
            metadata={
                "name": "omega",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        topography_height_no_negative: Quantity = dataclasses.field(
            metadata={
                "name": "topography_height_no_negative",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        sensible_heat_flux: Quantity = dataclasses.field(
            metadata={
                "name": "sensible_heat_flux",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        latent_heat_flux: Quantity = dataclasses.field(
            metadata={
                "name": "latent_heat_flux",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        longitude_degrees: Quantity = dataclasses.field(
            metadata={
                "name": "longitude_degrees",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        latitude_degrees: Quantity = dataclasses.field(
            metadata={
                "name": "latitude_degrees",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_old: Quantity = dataclasses.field(
            metadata={
                "name": "t_old",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor_old: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_old",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_modified_by_advection: Quantity = dataclasses.field(
            metadata={
                "name": "t_modified_by_advection",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor_modified_by_advection: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_modified_by_advection",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        geopotential_height_forced: Quantity = dataclasses.field(
            metadata={
                "name": "geopotential_height_forced",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_forced: Quantity = dataclasses.field(
            metadata={
                "name": "p_forced",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_surface: Quantity = dataclasses.field(
            metadata={
                "name": "p_surface",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_surface: Quantity = dataclasses.field(
            metadata={
                "name": "t_surface",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        u: Quantity = dataclasses.field(
            metadata={
                "name": "u",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        v: Quantity = dataclasses.field(
            metadata={
                "name": "v",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        w: Quantity = dataclasses.field(
            metadata={
                "name": "w",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        mass: Quantity = dataclasses.field(
            metadata={
                "name": "mass",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_scale_velosity: Quantity = dataclasses.field(
            metadata={
                "name": "convective_scale_velosity",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        buoyancy_excess: Quantity = dataclasses.field(
            metadata={
                "name": "buoyancy_excess",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_ice_input: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_ice_input",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_ice_input: Quantity = dataclasses.field(
            metadata={
                "name": "convective_ice_input",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_liquid_input: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_liquid_input",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_liquid_input: Quantity = dataclasses.field(
            metadata={
                "name": "convective_liquid_input",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_cloud_fraction_input: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_cloud_fraction_input",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_cloud_fraction_input: Quantity = dataclasses.field(
            metadata={
                "name": "convective_cloud_fraction_input",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    input: Input
    output: Output
    input_output: InputOutput
