from ndsl.dsl.typing import FloatField, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K, BACKWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)


def tracer_output(
    error_code: FloatFieldIJ_Plume,
    plume: Int,
    tup: FloatField,
    tempco: FloatField,
    t_cup: FloatField,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume]:
            tup = tempco
    with computation(PARALLEL), interval(-1, None):
        if error_code[0, 0][plume]:
            tup = t_cup
