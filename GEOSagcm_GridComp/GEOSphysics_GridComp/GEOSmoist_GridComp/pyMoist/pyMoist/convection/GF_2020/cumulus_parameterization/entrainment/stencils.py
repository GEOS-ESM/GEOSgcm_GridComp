from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ, Int
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    BACKWARD,
    computation,
    interval,
    int32,
    log,
    exp,
)
from ndsl.dsl.gt4py import function
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import saturation_vapor_pressure
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume

def entrainment_rates(
    vapor: FloatField,
    environment_saturation_mixing_ratio: FloatField,
    lcl_level: IntField,
    error_code: IntFieldIJ_Plume,
    entrainment_rate: FloatField,
    updraft_detrainment_function: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(...):
        if error_code[0,0][plume] == 0:
            frh = vapor