from ndsl.dsl.typing import (
    FloatField,
    FloatFieldIJ,
    Float,
    IntFieldIJ,
    Int,
    Bool,
    IntField,
)
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
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


def cup_up_aa0(
    dby: FloatField,
    gamma_cup: FloatField,
    ktop: IntField,
    kbcon: IntField,
    k22: IntField,
    z_cup: FloatField,
    t_cup: FloatField,
    zu: FloatField,
    integ: IntField,
    integ_interval: IntField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    aa0: FloatField,
):
    from __externals__ import k_start

    with computation(PARALLEL), interval(...):
        if integ == 1:
            if integ_interval == cumulus_parameterization_constants.BL:
                kbeg = k_start
                kend = kbcon - 2
            elif integ_interval == cumulus_parameterization_constants.CIN:
                kbeg = k22 - 1
                kend = kbcon - 2

        else:
            kbeg = kbcon - 1
            kend = ktop - 1

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= kbeg and K <= kend:
                dz = z_cup[0, 0, 1] - z_cup
                aa_1 = (
                    zu
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * t_cup)
                    )
                    * dby
                    / (1.0 + gamma_cup)
                )
                aa_2 = (
                    zu[0, 0, 1]
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * t_cup[0, 0, 1])
                    )
                    * dby[0, 0, 1]
                    / (1.0 + gamma_cup[0, 0, 1])
                )
                da = 0.5 * (aa_1 + aa_2) * dz
                aa0_below = aa0.at(K=K - 1)
                aa0 = aa0_below + da

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            aa0_temp = aa0.at(K=kend)
            aa0 = aa0_temp


def cloud_work_function_zero(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    aa1: FloatField,
):
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if aa1 == 0.0:
                error_code[0, 0][plume] = 17
                # ierrc[0,0][plume]="cloud work function zero"
