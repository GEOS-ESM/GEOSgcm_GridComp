from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, IntField
from ndsl.dsl.gt4py import (
    computation,
    PARALLEL,
    interval,
    FORWARD,
    K,
    BACKWARD,
    exp,
    function,
    int32,
    float32,
)
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)

# Parameters needed for get_wetbulb
RD = 287.06
RV = 461.52
RCPD = 1004.71
RTT = 273.16
HOH2O = 1000.0
RLVTT = 2.5008e6
RLSTT = 2.8345e6
RETV = RV / RD - 1.0
RLMLT = RLSTT - RLVTT
RCPV = 4.0 * RV
R2ES = 611.21 * RD / RV
R3LES = 17.502
R3IES = 22.587
R4LES = 32.19
R4IES = -0.7
R5LES = R3LES * (RTT - R4LES)
R5IES = R3IES * (RTT - R4IES)
R5ALVCP = R5LES * RLVTT / RCPD
R5ALSCP = R5IES * RLSTT / RCPD
RALVDCP = RLVTT / RCPD
RALSDCP = RLSTT / RCPD
RALFDCP = RLMLT / RCPD
RTWAT = RTT
RTBER = RTT - 5.0
RTBERCU = RTT - 5.0
RTICE = RTT - 23.0
RTICECU = RTT - 23.0
RTWAT_RTICE_R = 1.0 / (RTWAT - RTICE)
RTWAT_RTICECU_R = 1.0 / (RTWAT - RTICECU)
RVTMP2 = RCPV / RCPD - 1.0
ZQMAX = 0.5


@function
def get_wetbulb(
    local_t_cloud_levels,
    local_vapor_cloud_levels_forced,
    p_cloud_levels_forced,
):
    PT = local_t_cloud_levels
    PQ = local_vapor_cloud_levels_forced
    PSP = p_cloud_levels_forced * 100.0

    if PT > RTT:
        Z3ES = R3LES
        Z4ES = R4LES
        Z5ALCP = R5ALVCP
        ZALDCP = RALVDCP
    else:
        Z3ES = R3IES
        Z4ES = R4IES
        Z5ALCP = R5ALSCP
        ZALDCP = RALSDCP

    PTARE = PT
    ZQP = 1.0 / PSP

    FOEALFCU = min(
        1.0, ((max(RTICECU, min(RTWAT, PTARE)) - RTICECU) * RTWAT_RTICECU_R) ** 2
    )
    FOEEWMCU = R2ES * (
        FOEALFCU * exp(R3LES * (PTARE - RTT) / (PTARE - R4LES))
        + (1.0 - FOEALFCU) * exp(R3IES * (PTARE - RTT) / (PTARE - R4IES))
    )
    ZQSAT = FOEEWMCU * ZQP

    ZQSAT = min(cumulus_parameterization_constants.MAX_QSAT, ZQSAT)
    ZCOR = 1.0 / (1.0 - RETV * ZQSAT)
    ZQSAT = ZQSAT * ZCOR

    FOEDEMCU = FOEALFCU * R5ALVCP * (1.0 / (PTARE - R4LES) ** 2) + (
        1.0 - FOEALFCU
    ) * R5ALSCP * (1.0 / (PTARE - R4IES) ** 2)

    ZCOND = (PQ - ZQSAT) / (1.0 + ZQSAT * ZCOR * FOEDEMCU)

    ZCOND = min(ZCOND, 0.0)

    FOELDCPMCU = FOEALFCU * RALVDCP + (1.0 - FOEALFCU) * RALSDCP
    PT = PT + FOELDCPMCU * ZCOND

    PQ = PQ - ZCOND

    PTARE = PT

    FOEALFCU = min(
        1.0, ((max(RTICECU, min(RTWAT, PTARE)) - RTICECU) * RTWAT_RTICECU_R) ** 2
    )
    FOEEWMCU = R2ES * (
        FOEALFCU * exp(R3LES * (PTARE - RTT) / (PTARE - R4LES))
        + (1.0 - FOEALFCU) * exp(R3IES * (PTARE - RTT) / (PTARE - R4IES))
    )
    ZQSAT = FOEEWMCU * ZQP

    ZQSAT = min(0.5, ZQSAT)
    ZCOR = 1.0 / (1.0 - RETV * ZQSAT)
    ZQSAT = ZQSAT * ZCOR

    FOEDEMCU = FOEALFCU * R5ALVCP * (1.0 / (PTARE - R4LES) ** 2) + (
        1.0 - FOEALFCU
    ) * R5ALSCP * (1.0 / (PTARE - R4IES) ** 2)
    ZCOND1 = (PQ - ZQSAT) / (1.0 + ZQSAT * ZCOR * FOEDEMCU)

    if ZCOND == 0.0:
        ZCOND1 = min(ZCOND1, 0.0)

    FOELDCPMCU = FOEALFCU * RALVDCP + (1.0 - FOEALFCU) * RALSDCP
    PT = PT + FOELDCPMCU * ZCOND1
    PQ = PQ - ZCOND1

    local_vapor_wetbulb = PQ
    local_t_wetbulb = PT
    evap = -ZCOND1

    return (local_vapor_wetbulb, local_t_wetbulb)


def downdraft_wet_bulb(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    jmin: IntFieldIJ,
    local_vapor_cloud_levels_forced: FloatField,
    local_t_cloud_levels: FloatField,
    p_cloud_levels_forced: FloatField,
    local_vapor_wetbulb: FloatFieldIJ,
    local_t_wetbulb: FloatFieldIJ,
):
    from __externals__ import USE_WETBULB

    with computation(PARALLEL), interval(...):
        if USE_WETBULB == 1 and plume != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                lev: IntFieldIJ = jmin
                local_vapor_wetbulb, local_t_wetbulb = get_wetbulb(
                    local_vapor_cloud_levels_forced.at(K=lev),
                    local_t_cloud_levels.at(K=lev),
                    p_cloud_levels_forced.at(K=lev),
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


beta3 = Float(-1.13)
alpha3 = Float(1.9)


def cup_dd_edt(
    ccn: FloatField,
    cumulus: IntField,
    edtmax: FloatField,
    edtmin: FloatField,
    kbcon: IntField,
    ktop: IntField,
    maxens2: IntField,
    p: FloatField,
    psum2: FloatField,
    psumh: FloatField,
    pwav: FloatField,
    pwev: FloatField,
    us: FloatField,
    vs: FloatField,
    z: FloatField,
    aeroevap: IntField,
    sdp: FloatFieldIJ,
    vshear: FloatFieldIJ,
    vws: FloatFieldIJ,
    pef: FloatFieldIJ,
    dp: FloatFieldIJ,
    edt: FloatField,
    edtc: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):

    with computation(FORWARD), interval(...):
        edt = 0.0
        vshear = 0.0
        edtc = 0.0
    with computation(FORWARD), interval(...):
        sdp = 0.0
        vws = 0.0
        if cumulus != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                idx = kbcon - 1
                while idx >= kbcon - 1 and idx <= ktop - 1:
                    dp = p.at(K=idx) - p.at(K=idx + 1)
                    vws = vws + (
                        (
                            abs(
                                (us.at(K=idx + 1) - us.at(K=idx))
                                / (z.at(K=idx + 1) - z.at(K=idx))
                            )
                            + abs(
                                (vs.at(K=idx + 1) - vs.at(K=idx))
                                / (z.at(K=idx + 1) - z.at(K=idx))
                            )
                        )
                        * dp
                    )
                    sdp = sdp + dp
                    idx += 1
                vshear = 1.0e3 * vws / sdp

    with computation(FORWARD), interval(...):
        if cumulus != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                pef = (
                    float32(1.591)
                    - float32(0.639) * vshear
                    + float32(0.0953) * (vshear**2)
                    - float32(0.00496) * (vshear**3)
                )
                pef = min(pef, 0.9)
                pef = max(pef, 0.1)
                zkbc = z.at(K=kbcon - 1) * 3.281e-3
                prezk = 0.02

                if zkbc > 3.0:
                    prezk = 0.96729352 + zkbc * (
                        -0.70034167
                        + zkbc
                        * (
                            0.162179896
                            + zkbc
                            * (-1.2569798e-2 + zkbc * (4.2772e-4 - zkbc * 5.44e-6))
                        )
                    )

                if zkbc > 25.0:
                    prezk = 2.4

                pefb = 1.0 / (1.0 + prezk)
                pefb = min(pefb, 0.9)
                pefb = max(pefb, 0.1)

                edt = 1.0 - 0.5 * (pefb + pef)

                if aeroevap > 1:
                    aeroadd = (cumulus_parameterization_constants.CCNCLEAN**beta3) * (
                        (psumh) ** (alpha3 - 1)
                    )

                    prop_c = 0.5 * (pefb + pef) / aeroadd
                    aeroadd = (ccn**beta3) * ((psum2) ** (alpha3 - 1))

                    aeroadd = prop_c * aeroadd
                    pefc = aeroadd
                    if pefc > 0.9:
                        pefc = 0.9
                    if pefc < 0.1:
                        pefc = 0.1
                    edt = 1.0 - pefc
                    if aeroevap == 2:
                        edt = 1.0 - 0.25 * (pefb + pef + 2.0 * pefc)

                einc = 0.2 * edt
                if K <= maxens2 - 1:
                    edtc = edt + float(int32(K) - int32(1)) * einc

    with computation(PARALLEL), interval(...):
        if cumulus != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                if K <= maxens2 - 1:
                    edtc = -edtc * pwav / pwev
                    if edtc > edtmax:
                        edtc = edtmax
                    if edtc < edtmin:
                        edtc = edtmin

                temp = edtc.at(K=0)
                edtc = temp


def update_edto(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    maxens2: Int,
    sigd: FloatFieldIJ,
    edto: FloatFieldIJ,
    # edtc: #FloatFieldIJ_Maxens2??
):
    with computation(PARALLEL), interval(...):
        iedt = 1
        while iedt <= maxens2:
            if error_code[0, 0][plume] == 0:
                edto = sigd * edtc[0, 0][iedt]
                edt = edto
