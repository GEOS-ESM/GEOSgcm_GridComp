from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, IntField
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


def get_melting_profile(
    MELT_GLAC: IntField,
    cumulus: IntField,
    edto: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    melting_layer: FloatField,
    p_liq_ice: FloatField,
    po_cup: FloatField,
    pwdo: FloatField,
    pwo: FloatField,
    qrco: FloatField,
    tn_cup: FloatField,
    total_pwo_solid_phase: FloatFieldIJ,
    melting: FloatField,
):
    from __externals__ import k_end

    with computation(FORWARD), interval(...):
        ktf = k_end - 1
        if MELT_GLAC == 1 and cumulus == 1:
            pwo_solid_phase = 0.0
            pwo_eff = 0.0
            melting = 0.0

            if error_code[0, 0][plume] > 0:
                melting = 0.0

            total_pwo_solid_phase = 0.0

    with computation(FORWARD), interval(...):
        if MELT_GLAC == 1 and cumulus == 1:
            if K <= ktf - 1:
                if error_code[0, 0][plume] == 0:
                    dp = 100.0 * (po_cup - po_cup[0, 0, 1])

                    pwo_eff = 0.5 * (pwo + pwo[0, 0, 1])

                    pwo_solid_phase = (1.0 - p_liq_ice) * pwo_eff

                    total_pwo_solid_phase = (
                        total_pwo_solid_phase
                        + pwo_solid_phase * dp / constants.MAPL_GRAV
                    )

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == 1 and cumulus == 1:
            if K <= ktf:
                if error_code[0, 0][plume] == 0:
                    melting = melting_layer * (
                        total_pwo_solid_phase
                        / (
                            100
                            * (po_cup.at(K=0) - po_cup.at(K=ktf))
                            / constants.MAPL_GRAV
                        )
                    )
        else:
            melting = 0.0
