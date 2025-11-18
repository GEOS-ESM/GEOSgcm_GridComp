from ndsl.dsl.typing import FloatField, Int, IntField, Float,
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K, BACKWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)


@function
def get_buoy(
    hc: Float,
    he_cup: Float,
    hes_cup: Float,
    ierr: Int,
    kbcon: Int,
    klcl: Int,
    ktop: Int,
):

    dby = 0.0
    if ierr == 0:
        if K <= klcl:
            dby = hc - he_cup

        if K >= klcl + 1 and K <= ktop + 1:
            dby = hc - hes_cup

    return dby


def get_buoyancy(
    # In
    hc: FloatField,
    he_cup: FloatField,
    hes_cup: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    kbcon: IntField,
    klcl: IntField,
    ktop: IntField,
    # Out
    dby: FloatField,
):

    with computation(PARALLEL), interval(...):
        # Subtract 1 from k-level arrays to account for Python starting at zero
        dby = get_buoy(hc, he_cup, hes_cup, error_code[0,0][plume], kbcon - 1, klcl - 1, ktop - 1)
