import dataclasses

from ndsl import Quantity, State
from ndsl.constants import I_DIM, J_DIM, K_DIM
from ndsl.dsl.typing import Bool, Float
from ndsl.quantity.data_dimensions_field import DataDimensionsField

# dimension names for convection tracers
CONVECTION_TRACER_DIM = "convection_tracers"
SIZE_THREE_DIM = "size_three_dimension"
SIZE_FOUR_DIM = "size_four_dimension"


# field types with convection tracers ddim
FloatFieldIJ_ConvectionTracers = DataDimensionsField.declare()
FloatField_ConvectionTracers = DataDimensionsField.declare()
FloatField_ConvectionTracers_Plume = DataDimensionsField.declare()
ConvectionTracerMetaDataTable_Float = DataDimensionsField.declare()
ConvectionTracerMetaDataTable_Bool = DataDimensionsField.declare()
ConvectionTracerMetaDataTable_x3 = DataDimensionsField.declare()
ConvectionTracerMetaDataTable_x4 = DataDimensionsField.declare()


@dataclasses.dataclass
class ConvectionTracers(State):
    """
    Dataclass of Convection Tracers, contains both the numerical data of the tracers
    (stored in the "tracer" field) and metadata, each stored in its off-grid field

    Must be initialized with the following extra dimensions:
        "convection_tracers": number of convective tracers, must be defined prior to initalization
        "size_three_dimension": fixed dimension of size three for metadata
        "size_four_dimension": fixed dimension of size four for metadata
    """

    tracers: Quantity = dataclasses.field(
        metadata={
            "name": "tracers",
            "dims": [I_DIM, J_DIM, K_DIM, CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    fscav: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    vect_hcts: Quantity = dataclasses.field(
        metadata={
            "name": "vect_hcts",
            "dims": [CONVECTION_TRACER_DIM, "size_four_dimension"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    kc_scal: Quantity = dataclasses.field(
        metadata={
            "name": "kc_scal",
            "dims": [CONVECTION_TRACER_DIM, "size_three_dimension"],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    convfaci2g: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    retfactor: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    liq_and_gas: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    online_cldliq: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    online_vud: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    ftemp_threshold: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Float,
        }
    )
    use_gcc_washout: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )
    use_gocart: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )
    is_wetdep: Quantity = dataclasses.field(
        metadata={
            "name": "fscav",
            "dims": [CONVECTION_TRACER_DIM],
            "units": "?",
            "intent": "?",
            "dtype": Bool,
        }
    )
