import dataclasses

from ndsl import Local, LocalState
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GFDL1MLocals(LocalState):
    p_interface_mb: Local = dataclasses.field(
        metadata={
            "name": "p_interface_mb",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "millibars",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_mb: Local = dataclasses.field(
        metadata={
            "name": "p_mb",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "millibars",
            "intent": "?",
            "dtype": Float,
        }
    )
    edge_height_above_surface: Local = dataclasses.field(
        metadata={
            "name": "edge_height_above_surface",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_height_above_surface: Local = dataclasses.field(
        metadata={
            "name": "layer_height_above_surface",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_thickness: Local = dataclasses.field(
        metadata={
            "name": "layer_thickness",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_thickness_negative: Local = dataclasses.field(
        metadata={
            "name": "layer_thickness_negative",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    dz: Local = dataclasses.field(
        metadata={
            "name": "dz",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "Pa",
            "intent": "?",
            "dtype": Float,
        }
    )
    dp: Local = dataclasses.field(
        metadata={
            "name": "dp",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "Pa",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass: Local = dataclasses.field(
        metadata={
            "name": "mass",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_inverse: Local = dataclasses.field(
        metadata={
            "name": "mass_inverse",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    u_unmodified: Local = dataclasses.field(
        metadata={
            "name": "u_unmodified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    v_unmodified: Local = dataclasses.field(
        metadata={
            "name": "v_unmodified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    saturation_specific_humidity: Local = dataclasses.field(
        metadata={
            "name": "saturation_specific_humidity",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    dsaturation_specific_humidity: Local = dataclasses.field(
        metadata={
            "name": "dsaturation_specific_humidity",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lcl_level: Local = dataclasses.field(
        metadata={
            "name": "lcl_level",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
    total_concentration: Local = dataclasses.field(
        metadata={
            "name": "total_concentration",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )

    @dataclasses.dataclass
    class DriverTendencies:
        dvapordt: Local = dataclasses.field(
            metadata={
                "name": "dvapordt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dliquiddt: Local = dataclasses.field(
            metadata={
                "name": "dliquiddt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        draindt: Local = dataclasses.field(
            metadata={
                "name": "draindt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dicedt: Local = dataclasses.field(
            metadata={
                "name": "dicedt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dsnowdt: Local = dataclasses.field(
            metadata={
                "name": "dsnowdt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dgraupeldt: Local = dataclasses.field(
            metadata={
                "name": "dgraupeldt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dcloudfractiondt: Local = dataclasses.field(
            metadata={
                "name": "dcloudfractiondt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dtdt: Local = dataclasses.field(
            metadata={
                "name": "dtdt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dudt: Local = dataclasses.field(
            metadata={
                "name": "dudt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )
        dvdt: Local = dataclasses.field(
            metadata={
                "name": "dvdt",
                "dims": [X_DIM, Y_DIM, Z_DIM],
                "units": "?",
                "intent": "?",
                "dtype": Float,
            }
        )

    driver_tencencies: DriverTendencies
