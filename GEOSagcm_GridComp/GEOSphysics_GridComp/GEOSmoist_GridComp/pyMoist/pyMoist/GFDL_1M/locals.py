import dataclasses

from ndsl import Local, LocalState
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GFDL1MLocals(LocalState):
    p_interface_mb: Local = dataclasses.field(
        metadata={
            "name": "p_interface_mb",
            "dims": [I_DIM, J_DIM, K_INTERFACE_DIM],
            "units": "millibars",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_mb: Local = dataclasses.field(
        metadata={
            "name": "p_mb",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "millibars",
            "intent": "?",
            "dtype": Float,
        }
    )
    edge_height_above_surface: Local = dataclasses.field(
        metadata={
            "name": "edge_height_above_surface",
            "dims": [I_DIM, J_DIM, K_INTERFACE_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_height_above_surface: Local = dataclasses.field(
        metadata={
            "name": "layer_height_above_surface",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_thickness: Local = dataclasses.field(
        metadata={
            "name": "layer_thickness",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_thickness_negative: Local = dataclasses.field(
        metadata={
            "name": "layer_thickness_negative",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    dp: Local = dataclasses.field(
        metadata={
            "name": "dp",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "Pa",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass: Local = dataclasses.field(
        metadata={
            "name": "mass",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "kg m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_inverse: Local = dataclasses.field(
        metadata={
            "name": "mass_inverse",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "kg m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    u_unmodified: Local = dataclasses.field(
        metadata={
            "name": "u_unmodified",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    v_unmodified: Local = dataclasses.field(
        metadata={
            "name": "v_unmodified",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    saturation_specific_humidity: Local = dataclasses.field(
        metadata={
            "name": "saturation_specific_humidity",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dsaturation_specific_humidity: Local = dataclasses.field(
        metadata={
            "name": "dsaturation_specific_humidity",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lcl_level: Local = dataclasses.field(
        metadata={
            "name": "lcl_level",
            "dims": [I_DIM, J_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    total_concentration: Local = dataclasses.field(
        metadata={
            "name": "total_concentration",
            "dims": [I_DIM, J_DIM, K_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )

    @dataclasses.dataclass
    class DriverTendencies(LocalState):
        dvapordt: Local = dataclasses.field(
            metadata={
                "name": "dvapordt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dliquiddt: Local = dataclasses.field(
            metadata={
                "name": "dliquiddt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        draindt: Local = dataclasses.field(
            metadata={
                "name": "draindt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dicedt: Local = dataclasses.field(
            metadata={
                "name": "dicedt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dsnowdt: Local = dataclasses.field(
            metadata={
                "name": "dsnowdt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dgraupeldt: Local = dataclasses.field(
            metadata={
                "name": "dgraupeldt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloudfractiondt: Local = dataclasses.field(
            metadata={
                "name": "dcloudfractiondt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt: Local = dataclasses.field(
            metadata={
                "name": "dtdt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt: Local = dataclasses.field(
            metadata={
                "name": "dudt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt: Local = dataclasses.field(
            metadata={
                "name": "dvdt",
                "dims": [I_DIM, J_DIM, K_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    driver_tendencies: DriverTendencies
