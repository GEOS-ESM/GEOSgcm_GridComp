import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float


@dataclasses.dataclass
class GFDL1MDriverLocals(State):
    p_dry: Quantity = dataclasses.field(
        metadata={
            "name": "p_dry",
            "dims": [X_DIM, Y_DIM, Z_DIM],
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
    dz: Quantity = dataclasses.field(
        metadata={
            "name": "dz",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dp: Quantity = dataclasses.field(
        metadata={
            "name": "dp",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_fraction: Quantity = dataclasses.field(
        metadata={
            "name": "cloud_fraction",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density_unmodified: Quantity = dataclasses.field(
        metadata={
            "name": "density_unmodified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density: Quantity = dataclasses.field(
        metadata={
            "name": "density",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density_factor: Quantity = dataclasses.field(
        metadata={
            "name": "density_factor",
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
    one_minus_sigma: Quantity = dataclasses.field(
        metadata={
            "name": "one_minus_sigma",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ccn: Quantity = dataclasses.field(
        metadata={
            "name": "ccn",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    c_praut: Quantity = dataclasses.field(
        metadata={
            "name": "c_praut",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rh_limited: Quantity = dataclasses.field(
        metadata={
            "name": "rh_limited",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lhi: Quantity = dataclasses.field(
        metadata={
            "name": "lhi",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    icpk: Quantity = dataclasses.field(
        metadata={
            "name": "icpk",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hold_data: Quantity = dataclasses.field(
        metadata={
            "name": "hold_data",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liquid_precip_flux: Quantity = dataclasses.field(
        metadata={
            "name": "liquid_precip_flux",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice_precip_flux: Quantity = dataclasses.field(
        metadata={
            "name": "ice_precip_flux",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rain: Quantity = dataclasses.field(
        metadata={
            "name": "rain",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice: Quantity = dataclasses.field(
        metadata={
            "name": "ice",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    snow: Quantity = dataclasses.field(
        metadata={
            "name": "snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    graupel: Quantity = dataclasses.field(
        metadata={
            "name": "graupel",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation: Quantity = dataclasses.field(
        metadata={
            "name": "evaporation",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )

    @dataclasses.dataclass
    class TerminalSpeed:
        rain: Quantity = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Quantity = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Quantity = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Quantity = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class DryAirMixingRatio:
        vapor: Quantity = dataclasses.field(
            metadata={
                "name": "vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        liquid: Quantity = dataclasses.field(
            metadata={
                "name": "liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        rain: Quantity = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Quantity = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Quantity = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Quantity = dataclasses.field(
            metadata={
                "name": "graupel",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    @dataclasses.dataclass
    class Unmodified:
        @dataclasses.dataclass
        class MixingRatio:
            vapor: Quantity = dataclasses.field(
                metadata={
                    "name": "vapor",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            liquid: Quantity = dataclasses.field(
                metadata={
                    "name": "liquid",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            rain: Quantity = dataclasses.field(
                metadata={
                    "name": "rain",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            ice: Quantity = dataclasses.field(
                metadata={
                    "name": "ice",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            snow: Quantity = dataclasses.field(
                metadata={
                    "name": "snow",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            graupel: Quantity = dataclasses.field(
                metadata={
                    "name": "graupel",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )

        mixing_ratio: MixingRatio

    terminal_speed: TerminalSpeed
    dry_air_mixing_ratio: DryAirMixingRatio
    unmodified: Unmodified
