import dataclasses

from ndsl import Local, LocalState
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.typing import Float


@dataclasses.dataclass
class GFDL1MDriverLocals(LocalState):
    p_dry: Local = dataclasses.field(
        metadata={
            "name": "p_dry",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t: Local = dataclasses.field(
        metadata={
            "name": "t",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dz: Local = dataclasses.field(
        metadata={
            "name": "dz",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dp: Local = dataclasses.field(
        metadata={
            "name": "dp",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_fraction: Local = dataclasses.field(
        metadata={
            "name": "cloud_fraction",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density_unmodified: Local = dataclasses.field(
        metadata={
            "name": "density_unmodified",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density: Local = dataclasses.field(
        metadata={
            "name": "density",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density_factor: Local = dataclasses.field(
        metadata={
            "name": "density_factor",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass: Local = dataclasses.field(
        metadata={
            "name": "mass",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    u: Local = dataclasses.field(
        metadata={
            "name": "u",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    v: Local = dataclasses.field(
        metadata={
            "name": "v",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    w: Local = dataclasses.field(
        metadata={
            "name": "w",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    one_minus_sigma: Local = dataclasses.field(
        metadata={
            "name": "one_minus_sigma",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ccn: Local = dataclasses.field(
        metadata={
            "name": "ccn",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    c_praut: Local = dataclasses.field(
        metadata={
            "name": "c_praut",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rh_limited: Local = dataclasses.field(
        metadata={
            "name": "rh_limited",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lhi: Local = dataclasses.field(
        metadata={
            "name": "lhi",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    icpk: Local = dataclasses.field(
        metadata={
            "name": "icpk",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hold_data: Local = dataclasses.field(
        metadata={
            "name": "hold_data",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liquid_precip_flux: Local = dataclasses.field(
        metadata={
            "name": "liquid_precip_flux",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice_precip_flux: Local = dataclasses.field(
        metadata={
            "name": "ice_precip_flux",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rain: Local = dataclasses.field(
        metadata={
            "name": "rain",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice: Local = dataclasses.field(
        metadata={
            "name": "ice",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    snow: Local = dataclasses.field(
        metadata={
            "name": "snow",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    graupel: Local = dataclasses.field(
        metadata={
            "name": "graupel",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation: Local = dataclasses.field(
        metadata={
            "name": "evaporation",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )

    @dataclasses.dataclass
    class TerminalSpeed:
        rain: Local = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Local = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Local = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Local = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class DryAirMixingRatio:
        vapor: Local = dataclasses.field(
            metadata={
                "name": "vapor",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        liquid: Local = dataclasses.field(
            metadata={
                "name": "liquid",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        rain: Local = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Local = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Local = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Local = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Unmodified:
        @dataclasses.dataclass
        class MixingRatio:
            vapor: Local = dataclasses.field(
                metadata={
                    "name": "vapor",
                    "dims": [I_DIM, J_DIM, K_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            liquid: Local = dataclasses.field(
                metadata={
                    "name": "liquid",
                    "dims": [I_DIM, J_DIM, K_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            rain: Local = dataclasses.field(
                metadata={
                    "name": "rain",
                    "dims": [I_DIM, J_DIM, K_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            ice: Local = dataclasses.field(
                metadata={
                    "name": "ice",
                    "dims": [I_DIM, J_DIM, K_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            snow: Local = dataclasses.field(
                metadata={
                    "name": "snow",
                    "dims": [I_DIM, J_DIM, K_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            graupel: Local = dataclasses.field(
                metadata={
                    "name": "graupel",
                    "dims": [I_DIM, J_DIM, K_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )

        mixing_ratio: MixingRatio

    terminal_speed: TerminalSpeed
    dry_air_mixing_ratio: DryAirMixingRatio
    unmodified: Unmodified
