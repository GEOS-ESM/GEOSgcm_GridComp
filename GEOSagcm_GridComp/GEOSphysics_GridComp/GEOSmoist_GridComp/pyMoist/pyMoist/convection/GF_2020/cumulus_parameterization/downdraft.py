from ndsl import StencilFactory, QuantityFactory, Local, NDSLRuntime, Quantity
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import (
    GF2020CumulusParameterizationConfig,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import (
    GF2020CumulusParameterizationState,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import (
    GF2020CumulusParameterizationLocals,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as constants

from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ, Int, BoolFieldIJ
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
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import unknown_find_level
from ndsl.stencils.column_operations import column_max


def get_critical_level(
    error_code: IntFieldIJ_Plume,
    critical_level: IntFieldIJ,
    cloud_top_level: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    MAX_DOWNDRAFT_ORIGIN_HEIGHt: Float,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            cloud_top_height: FloatFieldIJ = (
                geopotential_height_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume])
                - topography_height_no_negative
            ) * 0.6
            cloud_top_height = min(
                cloud_top_height + topography_height_no_negative,
                MAX_DOWNDRAFT_ORIGIN_HEIGHt + topography_height_no_negative,
            )

        # setup mask to stop the next block
        stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and stop_computation == False:
            if geopotential_height_cloud_levels_forced >= cloud_top_height:
                critical_level = K
                stop_computation = True


def fill_plume_data_dimension(
    data: IntFieldIJ,
    data_dimension_field: IntFieldIJ_Plume,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        data_dimension_field[0, 0][plume] = data


def fill_from_plume_data_dimension(
    data: IntFieldIJ,
    data_dimension_field: IntFieldIJ_Plume,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        data = data_dimension_field[0, 0][plume]


def get_downdraft_origin_level(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    detrainment_start_level: IntFieldIJ,
    downdraft_origin_level: IntFieldIJ,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    geopotential_height_cloud_levels_forced: FloatField,
    melting_layer: FloatField,
    MINIMUM_DEPTH: Float,
    plume: Int,
):
    """
    Get the downdraft origin level

    For shallow plume, return 0 (downdraft is disabled). For mid and deep plume, perform full calculation.

    Args:
        error_code
        cloud_top_level
        updraft_lfc_level
        detrainment_start_level
        downdraft_origin_level
        environment_saturation_moist_static_energy_cloud_levels_forced
        geopotential_height_cloud_levels_forced
        melting_layer
        MINIMUM_DEPTH
        plume

    """
    from __externals__ import MELT_GLAC, k_end

    with computation(FORWARD), interval(0, 1):
        if plume == 0:
            downdraft_origin_level = 0
        elif plume == 1:
            # setup internal constants
            beta: FloatFieldIJ = 0.05
        elif plume == 2:
            # setup internal constants
            beta: FloatFieldIJ = 0.02

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            # predefine field for next block
            moist_static_energy_internal = 0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if plume == 2 and MELT_GLAC == True:  # noqa
                _, max_index = column_max(melting_layer, 0, k_end)
                downdraft_origin_level = max(downdraft_origin_level, max_index)

                downdraft_origin_level_internal = downdraft_origin_level
                keep_going = True
                while keep_going == True:
                    keep_going = False
                    if downdraft_origin_level_internal - 1 < detrainment_start_level:
                        detrainment_start_level = downdraft_origin_level_internal - 1
                    if downdraft_origin_level_internal >= cloud_top_level[0, 0][plume] - 1:
                        downdraft_origin_level_internal = cloud_top_level[0, 0][plume] - 2
                        level_initial = downdraft_origin_level_internal
                        moist_static_energy_internal[0, 0, level_initial] = (
                            environment_saturation_moist_static_energy_cloud_levels_forced.at(K=level_initial)
                        )
                        # dz = geopotential_height_cloud_levels_forced.at(K=level_initial+1) - geopotential_height_cloud_levels_forced.at(K=level_initial)
                        dh = 0
                        level = level_initial
                        stop_while_level = False
                        while level >= 0 and stop_while_level == False:
                            moist_static_energy_internal[0, 0, level] = (
                                environment_saturation_moist_static_energy_cloud_levels_forced.at(
                                    K=level_initial
                                )
                            )
                            dz = geopotential_height_cloud_levels_forced.at(
                                K=level + 1
                            ) - geopotential_height_cloud_levels_forced.at(K=level)
                            dh = dh + dz * (
                                moist_static_energy_internal.at(K=level)
                                - environment_saturation_moist_static_energy_cloud_levels_forced.at(K=level)
                            )
                            if dh >= 0:
                                downdraft_origin_level_internal = downdraft_origin_level_internal - 1
                                if downdraft_origin_level_internal >= 5:
                                    keep_going = True
                                else:
                                    error_code[0, 0][plume] = 9
                                    stop_while_level = True
                            level -= 1

                downdraft_origin_level = downdraft_origin_level_internal
                if downdraft_origin_level_internal <= 5:
                    error_code[0, 0][plume] = 4

    with computation(FORWARD), interval(0, 1):
        # must have at least depth_min m between cloud convective base and cloud top.
        if error_code[0, 0][plume] == 0:
            if downdraft_origin_level - 1 <= detrainment_start_level:
                detrainment_start_level = downdraft_origin_level - 1
            if (
                -geopotential_height_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume])
                + geopotential_height_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume])
                < MINIMUM_DEPTH
            ):
                error_code[0, 0][plume] = 6


class DowndraftOriginLevel(NDSLRuntime):
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make a modifyable copy of the QuantityFactory and add a data dimension
        self.quantity_factory = quantity_factory
        self.quantity_factory.add_data_dimensions(
            {"plume": cumulus_parameterization_constants.NUMBER_OF_PLUMES}
        )

        # make config and cumulus_parameterization_config visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # initalize locals
        self._critical_level: Local = self.make_local(quantity_factory, [X_DIM, Y_DIM], Int)
        self._downdraft_origin_level_ddim: Local = self.make_local(
            self.quantity_factory, [X_DIM, Y_DIM, "plume"], Int
        )

        # construct stencils
        self._get_critical_level = stencil_factory.from_dims_halo(
            func=get_critical_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._fill_plume_data_dimension = stencil_factory.from_dims_halo(
            func=fill_plume_data_dimension,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._unknown_find_level = stencil_factory.from_dims_halo(
            func=unknown_find_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._fill_from_plume_data_dimension = stencil_factory.from_dims_halo(
            func=fill_from_plume_data_dimension,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._get_downdraft_origin_level = stencil_factory.from_dims_halo(
            func=get_downdraft_origin_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"MELT_GLAC": cumulus_parameterization_config.MELT_GLAC},
        )

    def __call__(
        self,
        error_code: Quantity,
        cloud_top_level: Quantity,
        geopotential_height_cloud_levels_forced: Quantity,
        topography_height_no_negative: Quantity,
        environment_saturation_moist_static_energy_cloud_levels_forced: Quantity,
        updraft_origin_level: Quantity,
        downdraft_origin_level: Quantity,
        updraft_lfc_level: Quantity,
        detrainment_start_level: Quantity,
        melting_layer: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._get_critical_level(
            error_code=error_code,
            critical_level=self._critical_level,
            cloud_top_level=cloud_top_level,
            geopotential_height_cloud_levels_forced=geopotential_height_cloud_levels_forced,
            topography_height_no_negative=topography_height_no_negative,
            MAX_DOWNDRAFT_ORIGIN_HEIGHt=plume_dependent_constants.MAX_DOWNDRAFT_ORIGIN_HEIGHt,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._fill_plume_data_dimension(
            data=downdraft_origin_level,
            data_dimension_field=self._downdraft_origin_level_ddim,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._unknown_find_level(
            array=environment_saturation_moist_static_energy_cloud_levels_forced,
            start_index=updraft_origin_level,
            end_index=self._critical_level,
            out_index=self._downdraft_origin_level_ddim,
            error_code=error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._get_downdraft_origin_level(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            updraft_lfc_level=updraft_lfc_level,
            detrainment_start_level=detrainment_start_level,
            downdraft_origin_level=downdraft_origin_level,
            environment_saturation_moist_static_energy_cloud_levels_forced=environment_saturation_moist_static_energy_cloud_levels_forced,
            geopotential_height_cloud_levels_forced=geopotential_height_cloud_levels_forced,
            melting_layer=melting_layer,
            MINIMUM_DEPTH=plume_dependent_constants.MINIMUM_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


######## NOTE TODO NOTE README NOTE TODO TODO NOTE EVERYTHING BELOW HERE NEEDS TO BE REWORKED

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

    FOEALFCU = min(1.0, ((max(RTICECU, min(RTWAT, PTARE)) - RTICECU) * RTWAT_RTICECU_R) ** 2)
    FOEEWMCU = R2ES * (
        FOEALFCU * exp(R3LES * (PTARE - RTT) / (PTARE - R4LES))
        + (1.0 - FOEALFCU) * exp(R3IES * (PTARE - RTT) / (PTARE - R4IES))
    )
    ZQSAT = FOEEWMCU * ZQP

    ZQSAT = min(cumulus_parameterization_constants.MAX_QSAT, ZQSAT)
    ZCOR = 1.0 / (1.0 - RETV * ZQSAT)
    ZQSAT = ZQSAT * ZCOR

    FOEDEMCU = FOEALFCU * R5ALVCP * (1.0 / (PTARE - R4LES) ** 2) + (1.0 - FOEALFCU) * R5ALSCP * (
        1.0 / (PTARE - R4IES) ** 2
    )

    ZCOND = (PQ - ZQSAT) / (1.0 + ZQSAT * ZCOR * FOEDEMCU)

    ZCOND = min(ZCOND, 0.0)

    FOELDCPMCU = FOEALFCU * RALVDCP + (1.0 - FOEALFCU) * RALSDCP
    PT = PT + FOELDCPMCU * ZCOND

    PQ = PQ - ZCOND

    PTARE = PT

    FOEALFCU = min(1.0, ((max(RTICECU, min(RTWAT, PTARE)) - RTICECU) * RTWAT_RTICECU_R) ** 2)
    FOEEWMCU = R2ES * (
        FOEALFCU * exp(R3LES * (PTARE - RTT) / (PTARE - R4LES))
        + (1.0 - FOEALFCU) * exp(R3IES * (PTARE - RTT) / (PTARE - R4IES))
    )
    ZQSAT = FOEEWMCU * ZQP

    ZQSAT = min(0.5, ZQSAT)
    ZCOR = 1.0 / (1.0 - RETV * ZQSAT)
    ZQSAT = ZQSAT * ZCOR

    FOEDEMCU = FOEALFCU * R5ALVCP * (1.0 / (PTARE - R4LES) ** 2) + (1.0 - FOEALFCU) * R5ALSCP * (
        1.0 / (PTARE - R4IES) ** 2
    )
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
        if error_code[0, 0][plume] == 0 or cumulus != cumulus_parameterization_constants.shallow:
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
        if error_code[0, 0][plume] == 0 or cumulus != cumulus_parameterization_constants.shallow:
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
                        hcdo[0, 0, 1] * zdo[0, 0, 1] - 0.5 * dd_massdetro * hcdo[0, 0, 1] + dd_massentro * heo
                    ) / denom

                    dbydo = hcdo - heso_cup

                    bud = bud + dbydo * dzo
                else:
                    ucd = ucd[0, 0, 1]
                    vcd = vcd[0, 0, 1]
                    hcdo = hcdo[0, 0, 1]

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 or cumulus != cumulus_parameterization_constants.shallow:
            if bud > 0:
                ierr = 7
                # ierrc(i)='downdraft is not negatively buoyant '


beta3 = Float(-1.13)
alpha3 = Float(1.9)


def cup_dd_edt(
    ccn: FloatFieldIJ,
    local_epsilon_max: FloatFieldIJ,
    local_epsilon_min: FloatFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    p_forced: FloatField,
    local_psum: FloatFieldIJ,
    local_psumh: FloatFieldIJ,
    local_pwavo: FloatFieldIJ,
    local_pwevo: FloatFieldIJ,
    u: FloatField,
    v: FloatField,
    geopotential_height_forced: FloatField,
    local_epsilon: FloatFieldIJ,
    local_edtc: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):

    with computation(FORWARD), interval(...):
        local_epsilon = 0.0
        vshear: FloatFieldIJ = 0.0
        local_edtc = 0.0
    with computation(FORWARD), interval(...):
        sdp: FloatFieldIJ = 0.0
        vws: FloatFieldIJ = 0.0
        if plume != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                idx = updraft_lfc_level[0, 0][plume] - 1
                while idx >= updraft_lfc_level[0, 0][plume] - 1 and idx <= cloud_top_level[0, 0][plume] - 1:
                    dp: FloatFieldIJ = p_forced.at(K=idx) - p_forced.at(K=idx + 1)
                    vws = vws + (
                        (
                            abs(
                                (u.at(K=idx + 1) - u.at(K=idx))
                                / (
                                    geopotential_height_forced.at(K=idx + 1)
                                    - geopotential_height_forced.at(K=idx)
                                )
                            )
                            + abs(
                                (v.at(K=idx + 1) - v.at(K=idx))
                                / (
                                    geopotential_height_forced.at(K=idx + 1)
                                    - geopotential_height_forced.at(K=idx)
                                )
                            )
                        )
                        * dp
                    )
                    sdp = sdp + dp
                    idx += 1
                vshear = 1.0e3 * vws / sdp

    with computation(FORWARD), interval(...):
        if plume != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                pef: FloatFieldIJ = (
                    float32(1.591)
                    - float32(0.639) * vshear
                    + float32(0.0953) * (vshear**2)
                    - float32(0.00496) * (vshear**3)
                )
                pef = min(pef, 0.9)
                pef = max(pef, 0.1)
                zkbc = geopotential_height_forced.at(K=updraft_lfc_level[0, 0][plume] - 1) * 3.281e-3
                prezk = 0.02

                if zkbc > 3.0:
                    prezk = 0.96729352 + zkbc * (
                        -0.70034167
                        + zkbc * (0.162179896 + zkbc * (-1.2569798e-2 + zkbc * (4.2772e-4 - zkbc * 5.44e-6)))
                    )

                if zkbc > 25.0:
                    prezk = 2.4

                pefb = 1.0 / (1.0 + prezk)
                pefb = min(pefb, 0.9)
                pefb = max(pefb, 0.1)

                local_epsilon = 1.0 - 0.5 * (pefb + pef)

                if cumulus_parameterization_constants.AEROEVAP > 1:
                    aeroadd = (cumulus_parameterization_constants.CCNCLEAN**beta3) * (
                        (local_psumh) ** (alpha3 - 1)
                    )

                    prop_c = 0.5 * (pefb + pef) / aeroadd
                    aeroadd = (ccn**beta3) * ((local_psum) ** (alpha3 - 1))

                    aeroadd = prop_c * aeroadd
                    pefc = aeroadd
                    if pefc > 0.9:
                        pefc = 0.9
                    if pefc < 0.1:
                        pefc = 0.1
                    local_epsilon = 1.0 - pefc
                    if cumulus_parameterization_constants.AEROEVAP == 2:
                        local_epsilon = 1.0 - 0.25 * (pefb + pef + 2.0 * pefc)

                einc = 0.2 * local_epsilon
                if K <= cumulus_parameterization_constants.MAXENS2 - 1:
                    local_edtc = local_epsilon + float(int32(K) - int32(1)) * einc

    with computation(FORWARD), interval(...):
        if plume != cumulus_parameterization_constants.shallow:
            if error_code[0, 0][plume] == 0:
                if K <= cumulus_parameterization_constants.MAXENS2 - 1:
                    local_edtc = -local_edtc * local_pwavo / local_pwevo
                    if local_edtc > local_epsilon_max:
                        local_edtc = local_epsilon_max
                    if local_edtc < local_epsilon_min:
                        local_edtc = local_epsilon_min


def update_edto(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_scale_dependence_factor_downdraft: FloatFieldIJ,
    epsilon: FloatFieldIJ_Plume,
    local_edtc: FloatFieldIJ,
    local_epsilon: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        iedt = 0
        while iedt < cumulus_parameterization_constants.MAXENS2:
            if error_code[0, 0][plume] == 0:
                epsilon[0, 0][plume] = local_scale_dependence_factor_downdraft * local_edtc
                local_epsilon = epsilon[0, 0][plume]
            iedt += 1


class DowndraftNormalizedMassFlux:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftLateralMassFlux:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftWetBlub:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._downdraft_wet_bulb = stencil_factory.from_dims_halo(
            func=downdraft_wet_bulb,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"USE_WETBULB": cumulus_parameterization_config.USE_WETBULB},
        )

        if self.cumulus_parameterization_config.USE_WETBULB == 1:
            raise NotImplementedError(
                f"Coding limitation: USE_WETBULB = 0 is expected, getting USE_WETBULB = 1"
            )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        # self._downdraft_wet_bulb(
        #     error_code= state.output.error_code,
        #     plume=plume_dependent_constants.PLUME_INDEX,
        #     jmin=,
        #     local_vapor_cloud_levels_forced=locals.vapor_cloud_levels_forced,
        #     local_t_cloud_levels=locals.t_cloud_levels,
        #     p_cloud_levels_forced=state.output.p_cloud_levels_forced,
        #     local_vapor_wetbulb=locals.vapor_wetbulb,
        #     local_t_wetbulb=locals.t_wetbulb,
        # )


class DowndraftMoistStaticEnergyAndMoistureBudget:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._moist_static_energy_and_moisture_budget = stencil_factory.from_dims_halo(
            func=moist_static_energy_and_moisture_budget,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._moist_static_energy_and_moisture_budget(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            # heso_cup=,
            # u_cup=,
            # v_cup=,
            # cumulus=,
            # use_wetbulb=,
            # jmin=,
            # t_wetbulb=,
            # q_wetbulb=,
            # zo_cup=,
            # hc=,
            # zdo=,
            # dd_massdetro=,
            # dd_massentro=,
            # dd_massdetru=,
            # dd_massentru=,
            # us=,
            # vs=,
            # pgcon=,
            # heo=,
            # hcdo=,
        )


class DowndraftMoistureProperties:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftWindshear:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._cup_dd_edt = stencil_factory.from_dims_halo(
            func=cup_dd_edt,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._update_edto = stencil_factory.from_dims_halo(
            func=update_edto,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        # NOTE This code will break if MAXENS2 != 1

        self._cup_dd_edt(
            ccn=state.input_output.ccn,
            local_epsilon_max=locals.epsilon_max,
            local_epsilon_min=locals.epsilon_min,
            updraft_lfc_level=state.output.updraft_lfc_level,
            cloud_top_level=state.output.cloud_top_level,
            p_forced=state.input_output.p_forced,
            local_psum=locals.psum,
            local_psumh=locals.psumh,
            local_pwavo=locals.pwavo,
            local_pwevo=locals.pwevo,
            u=state.input_output.u,
            v=state.input_output.v,
            geopotential_height_forced=state.input_output.geopotential_height_forced,
            local_epsilon=locals.epsilon,
            local_edtc=locals.edtc,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._update_edto(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_scale_dependence_factor_downdraft=locals.scale_dependence_factor_downdraft,
            epsilon=state.output.epsilon,
            local_edtc=locals.edtc,
            local_epsilon=locals.epsilon,
        )
