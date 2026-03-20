import dataclasses

from ndsl import Quantity, State
from ndsl.constants import I_XIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GF2020Locals(State):
    @dataclasses.dataclass
    class DerivedState:
        edge_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "edge_height_above_surface",
                "dims": [I_XIM, J_DIM, K_INTERFACE_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        layer_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "layer_height_above_surface",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p: Quantity = dataclasses.field(
            metadata={
                "name": "p",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_kappa: Quantity = dataclasses.field(
            metadata={
                "name": "p_kappa",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        th: Quantity = dataclasses.field(
            metadata={
                "name": "th",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        mass: Quantity = dataclasses.field(
            metadata={
                "name": "mass",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        vertical_velocity: Quantity = dataclasses.field(
            metadata={
                "name": "vertical_velocity",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice_fraciton: Quantity = dataclasses.field(
            metadata={
                "name": "ice_fraciton",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        saturation_specific_humidity: Quantity = dataclasses.field(
            metadata={
                "name": "saturation_specific_humidity",
                "dims": [I_XIM, J_DIM, K_DIM],
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
        modified_area: Quantity = dataclasses.field(
            metadata={
                "name": "modified_area",
                "dims": [I_XIM, J_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_perturbation_horizontal: Quantity = dataclasses.field(
            metadata={
                "name": "t_perturbation_horizontal",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        t_perturbation_vertical: Quantity = dataclasses.field(
            metadata={
                "name": "t_perturbation_vertical",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dz: Quantity = dataclasses.field(
            metadata={
                "name": "dz",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        air_density: Quantity = dataclasses.field(
            metadata={
                "name": "air_density",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        scalar_diffusivity: Quantity = dataclasses.field(
            metadata={
                "name": "scalar_diffusivity",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class FlippedCopy:
        t_2m: Quantity = dataclasses.field(
            metadata={
                "name": "t_2m",
                "dims": [I_XIM, J_DIM],
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
        t_surface: Quantity = dataclasses.field(
            metadata={
                "name": "t_surface",
                "dims": [I_XIM, J_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p: Quantity = dataclasses.field(
            metadata={
                "name": "p",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        p_surface: Quantity = dataclasses.field(
            metadata={
                "name": "p_surface",
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
        vapor_current: Quantity = dataclasses.field(
            metadata={
                "name": "vapor_current",
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
                "name": "vertical_velocity",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        layer_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "layer_height_above_surface",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        edge_height_above_surface: Quantity = dataclasses.field(
            metadata={
                "name": "edge_height_above_surface",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        mass: Quantity = dataclasses.field(
            metadata={
                "name": "mass",
                "dims": [I_XIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        scalar_diffusivity: Quantity = dataclasses.field(
            metadata={
                "name": "scalar_diffusivity",
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
        pbl_level: Quantity = dataclasses.field(
            metadata={
                "name": "pbl_level",
                "dims": [I_XIM, J_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Int,
            }
        )

    do_this_column: Quantity = dataclasses.field(
        metadata={
            "name": "do_this_column",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    fix_out_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "fix_out_vapor",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_tendency_from_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "t_tendency_from_vapor",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rtgt: Quantity = dataclasses.field(
        metadata={
            "name": "rtgt",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    aot500: Quantity = dataclasses.field(
        metadata={
            "name": "aot500",
            "dims": [I_XIM, J_DIM],
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
    sensible_heat_flux: Quantity = dataclasses.field(
        metadata={
            "name": "sensible_heat_flux",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    topography_height: Quantity = dataclasses.field(
        metadata={
            "name": "topography_height",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ocean_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "ocean_fraction",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    grid_length: Quantity = dataclasses.field(
        metadata={
            "name": "grid_length",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    buoyancy_excess: Quantity = dataclasses.field(
        metadata={
            "name": "buoyancy_excess",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    grid_scale_forcing_t: Quantity = dataclasses.field(
        metadata={
            "name": "grid_scale_forcing_t",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    grid_scale_forcing_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "grid_scale_forcing_vapor",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    subgrid_scale_forcing_t: Quantity = dataclasses.field(
        metadata={
            "name": "subgrid_scale_forcing_t",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    subgrid_scale_forcing_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "subgrid_scale_forcing_vapor",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    advective_forcing_t: Quantity = dataclasses.field(
        metadata={
            "name": "advective_forcing_t",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    saturation_water_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "saturation_water_vapor",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvapordt: Quantity = dataclasses.field(
        metadata={
            "name": "dvapordt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dcloudicedt: Quantity = dataclasses.field(
        metadata={
            "name": "dcloudicedt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dudt: Quantity = dataclasses.field(
        metadata={
            "name": "dudt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvdt: Quantity = dataclasses.field(
        metadata={
            "name": "dvdt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dlarge_scale_icedt: Quantity = dataclasses.field(
        metadata={
            "name": "dlarge_scale_icedt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvective_icedt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvective_icedt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dlarge_scale_liquiddt: Quantity = dataclasses.field(
        metadata={
            "name": "dlarge_scale_liquiddt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvective_liquiddt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvective_liquiddt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dlarge_scale_cloud_fractiondt: Quantity = dataclasses.field(
        metadata={
            "name": "dlarge_scale_cloud_fractiondt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvective_cloud_fractiondt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvective_cloud_fractiondt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dbuoyancydt: Quantity = dataclasses.field(
        metadata={
            "name": "dbuoyancydt",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvection_tracersdt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvection_tracersdt",
            "dims": [I_XIM, J_DIM, K_DIM, "convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation_sublimation_tendency: Quantity = dataclasses.field(
        metadata={
            "name": "evaporation_sublimation_tendency",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convective_precip_flux: Quantity = dataclasses.field(
        metadata={
            "name": "convective_precip_flux",
            "dims": [I_XIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    precip: Quantity = dataclasses.field(
        metadata={
            "name": "precip",
            "dims": [I_XIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )

    derived_state: DerivedState
    flipped_copy: FlippedCopy
