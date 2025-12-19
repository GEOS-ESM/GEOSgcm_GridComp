import dataclasses

from ndsl import Local, LocalState
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float


@dataclasses.dataclass
class GFDL1MDriverLocals(LocalState):
    p_dry: Local = dataclasses.field(
        metadata={
            "name": "p_dry",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    t: Local = dataclasses.field(
        metadata={
            "name": "t",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dz: Local = dataclasses.field(
        metadata={
            "name": "dz",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dp: Local = dataclasses.field(
        metadata={
            "name": "dp",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    cloud_fraction: Local = dataclasses.field(
        metadata={
            "name": "cloud_fraction",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density_unmodified: Local = dataclasses.field(
        metadata={
            "name": "density_unmodified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density: Local = dataclasses.field(
        metadata={
            "name": "density",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    density_factor: Local = dataclasses.field(
        metadata={
            "name": "density_factor",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass: Local = dataclasses.field(
        metadata={
            "name": "mass",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    u: Local = dataclasses.field(
        metadata={
            "name": "u",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    v: Local = dataclasses.field(
        metadata={
            "name": "v",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    w: Local = dataclasses.field(
        metadata={
            "name": "w",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    one_minus_sigma: Local = dataclasses.field(
        metadata={
            "name": "one_minus_sigma",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ccn: Local = dataclasses.field(
        metadata={
            "name": "ccn",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    c_praut: Local = dataclasses.field(
        metadata={
            "name": "c_praut",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rh_limited: Local = dataclasses.field(
        metadata={
            "name": "rh_limited",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lhi: Local = dataclasses.field(
        metadata={
            "name": "lhi",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    icpk: Local = dataclasses.field(
        metadata={
            "name": "icpk",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    hold_data: Local = dataclasses.field(
        metadata={
            "name": "hold_data",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liquid_precip_flux: Local = dataclasses.field(
        metadata={
            "name": "liquid_precip_flux",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice_precip_flux: Local = dataclasses.field(
        metadata={
            "name": "ice_precip_flux",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    rain: Local = dataclasses.field(
        metadata={
            "name": "rain",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ice: Local = dataclasses.field(
        metadata={
            "name": "ice",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    snow: Local = dataclasses.field(
        metadata={
            "name": "snow",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    graupel: Local = dataclasses.field(
        metadata={
            "name": "graupel",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    evaporation: Local = dataclasses.field(
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
        rain: Local = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Local = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Local = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Local = dataclasses.field(
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
        vapor: Local = dataclasses.field(
            metadata={
                "name": "vapor",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        liquid: Local = dataclasses.field(
            metadata={
                "name": "liquid",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        rain: Local = dataclasses.field(
            metadata={
                "name": "rain",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        ice: Local = dataclasses.field(
            metadata={
                "name": "ice",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        snow: Local = dataclasses.field(
            metadata={
                "name": "snow",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        graupel: Local = dataclasses.field(
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
            vapor: Local = dataclasses.field(
                metadata={
                    "name": "vapor",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            liquid: Local = dataclasses.field(
                metadata={
                    "name": "liquid",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            rain: Local = dataclasses.field(
                metadata={
                    "name": "rain",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            ice: Local = dataclasses.field(
                metadata={
                    "name": "ice",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            snow: Local = dataclasses.field(
                metadata={
                    "name": "snow",
                    "dims": [X_DIM, Y_DIM, Z_DIM],
                    "units": "?",
                    "intent": "?",
                    "dtype": Float,
                }
            )
            graupel: Local = dataclasses.field(
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
