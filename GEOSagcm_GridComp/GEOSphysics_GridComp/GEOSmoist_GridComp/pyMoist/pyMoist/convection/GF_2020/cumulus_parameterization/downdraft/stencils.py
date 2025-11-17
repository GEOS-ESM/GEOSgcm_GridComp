from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K, BACKWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)


def moist_static_energy_and_moisture_budget(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    heso_cup: FloatField_Plume,
    u_cup: FloatField_Plume,
    v_cup: FloatField_Plume,
    cumulus: Int,
    use_wetbulb: Int,
    jmin: IntFieldIJ_Plume,
    t_wetbulb: FloatFieldIJ_Plume,
    q_wetbulb: FloatFieldIJ_Plume,
    zo_cup: FloatField_Plume,
    hc: FloatField_Plume,
    zdo: FloatField_Plume,
    dd_massdetro: FloatField_Plume,
    dd_massentro: FloatField_Plume,
    dd_massdetru: FloatField_Plume,
    dd_massentru: FloatField_Plume,
    us: FloatField_Plume,
    vs: FloatField_Plume,
    pgcon: FloatFieldIJ_Plume,
    heo: FloatField_Plume,
    hcdo: FloatField_Plume,
):
    with computation(PARALLEL), interval(...):
        hcdo = heso_cup
        ucd = u_cup
        vcd = v_cup
        dbydo = 0.0
        bud = 0.0

    with computation(PARALLEL), interval(...):
        if (
            error_code[0, 0][plume] == 0
            or cumulus != cumulus_parameterization_constants.shallow
        ):
            i_wb = 0
            if use_wetbulb == 1:
                if K == jmin:
                    hcdo = 0.5 * (
                        cumulus_parameterization_constants.CP * t_wetbulb
                        + cumulus_parameterization_constants.XLV * q_wetbulb
                        + zo_cup * constants.MAPL_GRAV
                        + hc
                    )
                i_wb = 1
            if K == jmin:
                dbydo = hcdo - heso_cup
                bud: FloatFieldIJ = dbydo * (zo_cup.at(K=jmin + 1) - zo_cup)

    with computation(BACKWARD), interval(...):
        if (
            error_code[0, 0][plume] == 0
            or cumulus != cumulus_parameterization_constants.shallow
        ):
            if K <= jmin - i_wb:
                denom: FloatFieldIJ = zdo[0, 0, 1] - 0.5 * dd_massdetro + dd_massentro
                denomU: FloatFieldIJ = zdo[0, 0, 1] - 0.5 * dd_massdetru + dd_massentru
                if denom > 0.0 and denomU > 0.0:
                    dzo: FloatFieldIJ = zo_cup[0, 0, 1] - zo_cup

                    ucd = (
                        ucd[0, 0, 1] * zdo[0, 0, 1]
                        - 0.5 * dd_massdetru * ucd[0, 0, 1]
                        + dd_massentru * us
                        - pgcon * zdo[0, 0, 1] * (us[0, 0, 1] - us)
                    ) / denomU
                    vcd = (
                        vcd[0, 0, 1] * zdo[0, 0, 1]
                        - 0.5 * dd_massdetru * vcd[0, 0, 1]
                        + dd_massentru * vs
                        - pgcon * zdo[0, 0, 1] * (vs[0, 0, 1] - vs)
                    ) / denomU

                    hcdo = (
                        hcdo[0, 0, 1] * zdo[0, 0, 1]
                        - 0.5 * dd_massdetro * hcdo[0, 0, 1]
                        + dd_massentro * heo
                    ) / denom

                    dbydo = hcdo - heso_cup

                    bud = bud + dbydo * dzo
                else:
                    ucd = ucd[0, 0, 1]
                    vcd = vcd[0, 0, 1]
                    hcdo = hcdo[0, 0, 1]

    with computation(PARALLEL), interval(...):
        if (
            error_code[0, 0][plume] == 0
            or cumulus != cumulus_parameterization_constants.shallow
        ):
            if bud > 0:
                ierr = 7
                # ierrc(i)='downdraft is not negatively buoyant '
