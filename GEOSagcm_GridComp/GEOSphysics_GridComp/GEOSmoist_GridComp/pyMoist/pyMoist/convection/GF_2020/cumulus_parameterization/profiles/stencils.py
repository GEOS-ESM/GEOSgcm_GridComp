from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
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
    tempcdo: FloatField_Plume,
    hcdo: FloatField_Plume,
    zo_cup: FloatField_Plume,
    qcdo: FloatField_Plume,
    tn_cup: FloatField_Plume,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            tempcdo = (1.0 / cumulus_parameterization_constants.CP) * (
                hcdo
                - constants.MAPL_GRAV * zo_cup
                - cumulus_parameterization_constants.XLV * qcdo
            )
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            tempcdo = tn_cup
