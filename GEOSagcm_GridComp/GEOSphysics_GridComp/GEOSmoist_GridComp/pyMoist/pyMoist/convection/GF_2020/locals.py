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
        dz: Quantity = dataclasses.field(
            metadata={
                "name": "dz",
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
        scalar_diffusivity: Quantity = dataclasses.field(
            metadata={
                "name": "scalar_diffusivity",
                "dims": [X_DIM, Y_DIM, Z_DIM],
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
        w: Quantity = dataclasses.field(
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

    do_this_column: Quantity = dataclasses.field(
        metadata={
            "name": "do_this_column",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    total_dbuoyancydt: Quantity = dataclasses.field(
        metadata={
            "name": "total_dbuoyancydt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    fix_out_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "fix_out_vapor",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t_tendency_from_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "t_tendency_from_vapor",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rtgt: Quantity = dataclasses.field(
        metadata={
            "name": "rtgt",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    aot500: Quantity = dataclasses.field(
        metadata={
            "name": "aot500",
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
    vapor_current: Quantity = dataclasses.field(
        metadata={
            "name": "vapor_current",
            "dims": [X_DIM, Y_DIM, Z_DIM],
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
    saturation_water_vapor: Quantity = dataclasses.field(
        metadata={
            "name": "saturation_water_vapor",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dtdt: Quantity = dataclasses.field(
        metadata={
            "name": "dtdt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvapordt: Quantity = dataclasses.field(
        metadata={
            "name": "dvapordt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dcloudicedt: Quantity = dataclasses.field(
        metadata={
            "name": "dcloudicedt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dudt: Quantity = dataclasses.field(
        metadata={
            "name": "dudt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dvdt: Quantity = dataclasses.field(
        metadata={
            "name": "dvdt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dlarge_scale_icedt: Quantity = dataclasses.field(
        metadata={
            "name": "dlarge_scale_icedt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvective_icedt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvective_icedt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dlarge_scale_liquiddt: Quantity = dataclasses.field(
        metadata={
            "name": "dlarge_scale_liquiddt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvective_liquiddt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvective_liquiddt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dlarge_scale_cloud_fractiondt: Quantity = dataclasses.field(
        metadata={
            "name": "dlarge_scale_cloud_fractiondt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvective_cloud_fractiondt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvective_cloud_fractiondt",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dconvection_tracersdt: Quantity = dataclasses.field(
        metadata={
            "name": "dconvection_tracersdt",
            "dims": [X_DIM, Y_DIM, Z_DIM, "convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation_sublimation_tendency: Quantity = dataclasses.field(
        metadata={
            "name": "evaporation_sublimation_tendency",
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

    derived_state: DerivedState
    flipped_copy: FlippedCopy
