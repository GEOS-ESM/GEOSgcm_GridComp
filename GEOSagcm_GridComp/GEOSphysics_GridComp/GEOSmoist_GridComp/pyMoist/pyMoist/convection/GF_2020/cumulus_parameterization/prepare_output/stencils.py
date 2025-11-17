from ndsl.dsl.gt4py import (
    PARALLEL,
    computation,
    interval,
    FORWARD,
    function,
    BACKWARD,
    K,
)
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, Int, IntFieldIJ
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.shared_incloud_processes import ice_fraction


def prepare_output(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    xmb: FloatFieldIJ,
    pwavo: FloatField,
    pwevo: FloatField,
    zuo: FloatField,
    zdo: FloatField,
    pwo: FloatField,
    pwdo: FloatField,
    up_massentro: FloatField,
    up_massdetro: FloatField,
    dd_massentro: FloatField,
    dd_massdetro: FloatField,
    zenv: FloatField,
    subten_Q: FloatField,
    subten_H: FloatField,
    subten_T: FloatField,
):
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            pwavo = xmb * pwavo
            pwevo = xmb * pwevo
            zuo = xmb * zuo
            zdo = xmb * zdo
            pwo = xmb * pwo
            pwdo = xmb * pwdo
            up_massentro = xmb * up_massentro
            up_massdetro = xmb * up_massdetro
            dd_massentro = xmb * dd_massentro
            dd_massdetro = xmb * dd_massdetro
            zenv = xmb * zenv

    with computation(PARALLEL), interval(...):
        subten_Q = xmb * subten_Q
        subten_H = xmb * subten_H
        subten_T = xmb * subten_T
