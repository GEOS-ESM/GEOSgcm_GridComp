import dataclasses

from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import Float, Int


@dataclasses.dataclass
class GFDL1MLocals(State):
    p_interface_mb: Quantity = dataclasses.field(
        metadata={
            "name": "p_interface_mb",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "millibars",
            "intent": "?",
            "dtype": Float,
        }
    )
    p_mb: Quantity = dataclasses.field(
        metadata={
            "name": "p_mb",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "millibars",
            "intent": "?",
            "dtype": Float,
        }
    )
    edge_height_above_surface: Quantity = dataclasses.field(
        metadata={
            "name": "edge_height_above_surface",
            "dims": [X_DIM, Y_DIM, Z_INTERFACE_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_height_above_surface: Quantity = dataclasses.field(
        metadata={
            "name": "layer_height_above_surface",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_thickness: Quantity = dataclasses.field(
        metadata={
            "name": "layer_thickness",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    layer_thickness_negative: Quantity = dataclasses.field(
        metadata={
            "name": "layer_thickness_negative",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m",
            "intent": "?",
            "dtype": Float,
        }
    )
    dp: Quantity = dataclasses.field(
        metadata={
            "name": "dp",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "Pa",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass: Quantity = dataclasses.field(
        metadata={
            "name": "mass",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    mass_inverse: Quantity = dataclasses.field(
        metadata={
            "name": "mass_inverse",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "kg m-2",
            "intent": "?",
            "dtype": Float,
        }
    )
    u_unmodified: Quantity = dataclasses.field(
        metadata={
            "name": "u_unmodified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m s-1",
            "intent": "?",
            "dtype": Float,
        }
    )
    v_unmodified: Quantity = dataclasses.field(
        metadata={
            "name": "v_unmodified",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "m s-1",
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
    dsaturation_specific_humidity: Quantity = dataclasses.field(
        metadata={
            "name": "dsaturation_specific_humidity",
            "dims": [X_DIM, Y_DIM, Z_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    lcl_level: Quantity = dataclasses.field(
        metadata={
            "name": "lcl_level",
            "dims": [X_DIM, Y_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Int,
        }
    )
