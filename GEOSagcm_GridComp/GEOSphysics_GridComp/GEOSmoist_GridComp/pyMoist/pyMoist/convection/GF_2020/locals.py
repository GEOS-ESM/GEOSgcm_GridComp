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
        ocean_fraction: Quantity = dataclasses.field(
            metadata={
                "name": "ocean_fraction",
                "dims": [X_DIM, Y_DIM],
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
        t_surface: Quantity = dataclasses.field(
            metadata={
                "name": "t_surface",
                "dims": [X_DIM, Y_DIM],
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
        p_surface: Quantity = dataclasses.field(
            metadata={
                "name": "p_surface",
                "dims": [X_DIM, Y_DIM],
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
        pbl_level: Quantity = dataclasses.field(
            metadata={
                "name": "pbl_level",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Int,
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
        t_surface: Quantity = dataclasses.field(
            metadata={
                "name": "t_surface",
                "dims": [X_DIM, Y_DIM],
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

    @dataclasses.dataclass
    class CumulusParameterizationOutput:
        dtdt: Quantity = dataclasses.field(
            metadata={
                "name": "dtdt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt: Quantity = dataclasses.field(
            metadata={
                "name": "dudt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt: Quantity = dataclasses.field(
            metadata={
                "name": "dvdt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvapordt: Quantity = dataclasses.field(
            metadata={
                "name": "dvapordt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
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
        dcloudicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dcloudicedt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dnicedt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dnliquiddt: Quantity = dataclasses.field(
            metadata={
                "name": "dnliquiddt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dbuoyancydt: Quantity = dataclasses.field(
            metadata={
                "name": "dbuoyancydt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveicedt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleicedt: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleicedt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectiveliquiddt: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectiveliquiddt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescaleliquiddt: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescaleliquiddt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dconvectivecloudfractiondt: Quantity = dataclasses.field(
            metadata={
                "name": "dconvectivecloudfractiondt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dlargescalecloudfractiondt: Quantity = dataclasses.field(
            metadata={
                "name": "dlargescalecloudfractiondt",
                "dims": [X_DIM, Y_DIM, Z_DIM, "plumes"],
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
    class CumulusParameterizationInternal:
        t: Quantity = dataclasses.field(
            metadata={
                "name": "t",
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
        t_pbl: Quantity = dataclasses.field(
            metadata={
                "name": "t_pbl",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vapor_pbl: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_pbl",
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
        moist_static_energy: Quantity = dataclasses.field(
            metadata={
                "name": "moist_static_energy",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        kbmax: Quantity = dataclasses.field(
            metadata={
                "name": "kbmax",
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
        ierr2: Quantity = dataclasses.field(
            metadata={
                "name": "ierr2",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Int,
            }
        )
        ierr3: Quantity = dataclasses.field(
            metadata={
                "name": "ierr3",
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
        ierrc: Quantity = dataclasses.field(
            metadata={
                "name": "ierrc",
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
        cloud_work_function_1_fa: Quantity = dataclasses.field(
            metadata={
                "name": "cloud_work_function_1_fa",
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
        cin: Quantity = dataclasses.field(
            metadata={
                "name": "cin",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        xk_x: Quantity = dataclasses.field(
            metadata={
                "name": "xk_x",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        edt: Quantity = dataclasses.field(
            metadata={
                "name": "edt",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        edto: Quantity = dataclasses.field(
            metadata={
                "name": "edto",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_wetblub: Quantity = dataclasses.field(
            metadata={
                "name": "t_wetblub",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        q_wetbulb: Quantity = dataclasses.field(
            metadata={
                "name": "q_wetbulb",
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
        xf_dicycle: Quantity = dataclasses.field(
            metadata={
                "name": "xf_dicycle",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        x_add_buoy: Quantity = dataclasses.field(
            metadata={
                "name": "x_add_buoy",
                "dims": [X_DIM, Y_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        z: Quantity = dataclasses.field(
            metadata={
                "name": "z",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        xz: Quantity = dataclasses.field(
            metadata={
                "name": "xz",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        hcot: Quantity = dataclasses.field(
            metadata={
                "name": "hcot",
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
        xf_ens: Quantity = dataclasses.field(
            metadata={
                "name": "xf_ens",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        pr_ens: Quantity = dataclasses.field(
            metadata={
                "name": "pr_ens",
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
        cup_max_incremenet: Quantity = dataclasses.field(
            metadata={
                "name": "cup_max_incremenet",
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
    cumulus_parameterization_internal: CumulusParameterizationInternal
