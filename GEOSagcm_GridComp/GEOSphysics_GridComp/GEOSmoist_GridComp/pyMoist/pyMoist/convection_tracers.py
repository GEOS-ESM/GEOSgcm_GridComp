import dataclasses
from ndsl import Quantity, State
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import Float, Bool

import dataclasses
from typing import Dict


@dataclasses.dataclass
class ConvectionTracers(State):
    """
    Dataclass of Convection Tracers, contains both the numerical data of the tracers
    (stored in the "tracer" field) and metadata, each stored in its off-grid field

    Must be initalized with the following extra dimensions:
        "convection_tracers": number of convective tracers, must be defined prior to initalization
        "size_three_dimension": fixed dimension of size three for metadata
        "size_four_dimension": fixed dimension of size four for metadata
    """

    tracers: Quantity = dataclasses.field(
        metadata={
            "name": "tracers",
            "dims": [X_DIM, Y_DIM, Z_DIM, "convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    fscav: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vect_hcts: Quantity = dataclasses.field(
        metadata={
            "name": "vect_hcts",
            "dims": ["convection_tracers", "size_four_dimension"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    kc_scal: Quantity = dataclasses.field(
        metadata={
            "name": "kc_scal",
            "dims": ["convection_tracers", "size_three_dimension"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convfaci2g: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    retfactor: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liq_and_gas: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    online_cldliq: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    online_vud: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ftemp_threshold: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    use_gcc_washout: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )
    use_gocart: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )
    is_wetdep: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": ["convection_tracers"],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )
