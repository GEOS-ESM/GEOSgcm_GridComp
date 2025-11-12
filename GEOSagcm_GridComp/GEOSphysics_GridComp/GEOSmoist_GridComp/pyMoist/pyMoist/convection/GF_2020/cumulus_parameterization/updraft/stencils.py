from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, Bool
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
    FloatFieldIJ,
)


def in_cloud_updraft_air_temperature(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    FIRST_GUESS_W: Bool,
    tempco: FloatField_Plume,
    hco: FloatField_Plume,
    zo_cup: FloatField_Plume,
    qco: FloatField_Plume,
    tn_cup: FloatField_Plume,
):
    with computation(PARALLEL), interval(0, -1):
        if FIRST_GUESS_W == False:
            if error_code[0, 0][plume] == 0:
                tempco = (1.0 / cumulus_parameterization_constants.CP) * (
                    hco
                    - constants.MAPL_GRAV * zo_cup
                    - cumulus_parameterization_constants.XLV * qco
                )

    with computation(PARALLEL), interval(-1, None):
        if FIRST_GUESS_W == False:
            if error_code[0, 0][plume] == 0:
                tempco = tn_cup
    with computation(PARALLEL), interval(...):
        if FIRST_GUESS_W == False:
            if error_code[0, 0][plume] != 0:
                tempco = tn_cup
