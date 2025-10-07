import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020Locals(State):
    @dataclasses.dataclass
    class DerivedState:
        edge_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "edge_height_above_surface",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        layer_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "layer_height_above_surface",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p: Quantity = dataclasses.field(
            metadata={
                "name": "p",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_kappa: Quantity = dataclasses.field(
            metadata={
                "name": "p_kappa",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        th: Quantity = dataclasses.field(
            metadata={
                "name": "th",
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
        vertical_velocity: Quantity = dataclasses.field(
            metadata={
                "name": "vertical_velocity",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice_fraciton: Quantity = dataclasses.field(
            metadata={
                "name": "ice_fraciton",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        saturation_specific_humidity: Quantity = dataclasses.field(
            metadata={
                "name": "saturation_specific_humidity",
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
        modified_area: Quantity = dataclasses.field(
            metadata={
                "name": "modified_area",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        topography_height: Quantity = dataclasses.field(
            metadata={
                "name": "topography_height",
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
        advective_forcing_t: Quantity = dataclasses.field(
            metadata={
                "name": "advective_forcing_t",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_perturbation_horizontal: Quantity = dataclasses.field(
            metadata={
                "name": "t_perturbation_horizontal",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_perturbation_vertical: Quantity = dataclasses.field(
            metadata={
                "name": "t_perturbation_vertical",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Output:
        precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "precip_flux",
                "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        evap_subl_tendency: Quantity = dataclasses.field(
            metadata={
                "name": "evap_subl_tendency",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class LocalCopy:
        t_2m: Quantity = dataclasses.field(
            metadata={
                "name": "t_2m",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        evaporation: Quantity = dataclasses.field(
            metadata={
                "name": "evaporation",
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
        t: Quantity = dataclasses.field(
            metadata={
                "name": "t",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p: Quantity = dataclasses.field(
            metadata={
                "name": "p",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_mb: Quantity = dataclasses.field(
            metadata={
                "name": "p_mb",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor: Quantity = dataclasses.field(
            metadata={
                "name": "vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor_current: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_current",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        vertical_velocity: Quantity = dataclasses.field(
            metadata={
                "name": "vertical_velocity",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        layer_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "layer_height_above_surface",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        edge_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "edge_height_above_surface",
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
        scalar_diffusivity: Quantity = dataclasses.field(
            metadata={
                "name": "scalar_diffusivity",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        convective_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "convective_liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_ice: Quantity = dataclasses.field(
            metadata={
                "name": "convective_ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_cloud_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "convective_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_ice: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_cloud_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class CumulusParameterizationInput:
        geopotential_height: Quantity = dataclasses.field(
            metadata={
                "name": "geopotential_height",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_mb: Quantity = dataclasses.field(
            metadata={
                "name": "p_mb",
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
        t: Quantity = dataclasses.field(
            metadata={
                "name": "t",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor_timestep_start: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_timestep_start",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor_current: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_current",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        omega: Quantity = dataclasses.field(
            metadata={
                "name": "omega",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        mass: Quantity = dataclasses.field(
            metadata={
                "name": "mass",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        topography_height: Quantity = dataclasses.field(
            metadata={
                "name": "topography_height",
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
        grid_length: Quantity = dataclasses.field(
            metadata={
                "name": "grid_length",
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
        t_perturbation: Quantity = dataclasses.field(
            metadata={
                "name": "t_perturbation",
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
        vapor_excess: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_excess",
                "dims": [X_DIM, Y_DIM],
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
        ccn: Quantity = dataclasses.field(
            metadata={
                "name": "ccn",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "convective_liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_ice: Quantity = dataclasses.field(
            metadata={
                "name": "convective_ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_cloud_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "convective_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_liquid: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_ice: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        large_scale_cloud_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "large_scale_cloud_fraction",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        longitude_degrees: Quantity = dataclasses.field(
            metadata={
                "name": "longitude_degrees",
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
                "dtype": Int,
            }
        )
        pbl_height: Quantity = dataclasses.field(
            metadata={
                "name": "pbl_height",
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
        convective_scale_velosity: Quantity = dataclasses.field(
            metadata={
                "name": "convective_scale_velosity",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class CumulusParameterizationOutput:
        dtdt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dudt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dudt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dudt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt_combined: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt_combined",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloudicedt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dcloudicedt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloudicedt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dcloudicedt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloudicedt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dcloudicedt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnicedt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dnicedt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnicedt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dnicedt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnicedt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dnicedt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnliquiddt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dnliquiddt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnliquiddt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dnliquiddt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnliquiddt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dnliquiddt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dbuoyancydt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dbuoyancydt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dbuoyancydt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dbuoyancydt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dbuoyancydt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dbuoyancydt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveicedt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveicedt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveicedt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveicedt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveicedt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveicedt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleicedt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleicedt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleicedt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleicedt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleicedt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleicedt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveliquiddt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveliquiddt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveliquiddt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveliquiddt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveliquiddt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveliquiddt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleliquiddt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleliquiddt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleliquiddt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleliquiddt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleliquiddt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleliquiddt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectivecloudfractiondt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectivecloudfractiondt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectivecloudfractiondt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectivecloudfractiondt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectivecloudfractiondt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectivecloudfractiondt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescalecloudfractiondt_shallow: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescalecloudfractiondt_shallow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescalecloudfractiondt_mid: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescalecloudfractiondt_mid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescalecloudfractiondt_deep: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescalecloudfractiondt_deep",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        evap_subl_tendency: Quantity = dataclasses.field(
            metadata={
                "name": "evap_subl_tendency",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        convective_precip_flux: Quantity = dataclasses.field(
            metadata={
                "name": "convective_precip_flux",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class MiscelaneousDiagnostic:
        fix_out_vapor: Quantity = dataclasses.field(
            metadata={
                "name": "fix_out_vapor",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        last_ierr: Quantity = dataclasses.field(
            metadata={
                "name": "last_ierr",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        conprr: Quantity = dataclasses.field(
            metadata={
                "name": "conprr",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    derived_state: DerivedState
    output: Output
    local_copy: LocalCopy
    cumulus_parameterization_input: CumulusParameterizationInput
    cumulus_parameterization_output: CumulusParameterizationOutput
    miscelaneous_diagnostic: MiscelaneousDiagnostic
