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
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    get_updraft_origin_conditions,
)


def parcel_moist_static_energy(
    error_code: IntFieldIJ_Plume,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    updraft_origin_level: IntFieldIJ,
    p: FloatField,
    environmenet_moist_static_energy: FloatField,
    environmenet_moist_static_energy_forced: FloatField,
    t_perturbation: FloatField,
    moist_static_energy_origin_level: FloatFieldIJ,
    moist_static_energy_origin_level_forced: FloatFieldIJ,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    from __externals__ import k_end, BOUNDARY_CONDITION_METHOD

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            modification = (
                cumulus_parameterization_constants.XLV * vapor_excess
                + cumulus_parameterization_constants.CP * t_excess
            ) + add_buoyancy

            moist_static_energy_origin_level = get_updraft_origin_conditions(
                field=environmenet_moist_static_energy,
                scalar_perturbation=modification,
                p=p,
                updraft_origin_level=updraft_origin_level,
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=t_perturbation,
            )
            moist_static_energy_origin_level_forced = get_updraft_origin_conditions(
                field=environmenet_moist_static_energy_forced,
                scalar_perturbation=modification,
                p=p,
                updraft_origin_level=updraft_origin_level,
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=t_perturbation,
            )


def first_guess_mse(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    start_level: IntFieldIJ_Plume,
    ktop: IntFieldIJ_Plume,
    up_massdetro: FloatField_Plume,
    up_massentro: FloatField_Plume,
    zuo: FloatField_Plume,
    zu: FloatField_Plume,
    heo: FloatField_Plume,
    zqexec: FloatFieldIJ_Plume,
    hco: FloatField_Plume,
    heso_cup: FloatField_Plume,
    ztexec: FloatFieldIJ_Plume,
    x_add_buoy: FloatFieldIJ_Plume,
):

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= ktop + 1:
                denom: FloatFieldIJ = (
                    zu.at(K=K - 1)
                    - 0.5 * up_massdetro.at(K=K - 1)
                    + up_massentro.at(K=K - 1)
                )
                if denom > 0.0:
                    hco = (
                        hco.at(K=K - 1) * zuo.at(K=K - 1)
                        - 0.5 * up_massdetro.at(K=K - 1) * hco.at(K=K - 1)
                        + up_massentro.at(K=K - 1) * heo.at(K=K - 1)
                    ) / denom
                    if K == start_level + 1:
                        x_add: FloatFieldIJ = (
                            cumulus_parameterization_constants.XLV * zqexec
                            + cumulus_parameterization_constants.CP * ztexec
                        ) + x_add_buoy
                        hco = hco + x_add * up_massentro.at(K=K - 1) / denom

            else:
                hco = hco.at(K=K - 1)

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= ktop + 2:
                hco = heso_cup


def moist_static_energy_inside_cloud(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    start_level: IntFieldIJ,
    xhc: FloatField,
    xhkb: FloatFieldIJ,
    ktop: IntFieldIJ,
    up_massdetro: FloatField,
    up_massentro: FloatField,
    xzu: FloatField,
    xhe: FloatField,
    p_liq_ice: FloatField,
    zqexec: FloatFieldIJ,
    ztexec: FloatFieldIJ,
    x_add_buoy: FloatFieldIJ,
    qrco: FloatField,
    xhes_cup: FloatField,
):
    with computation(PARALLEL), interval(...):
        xhc = 0.0

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= start_level:
                xhc = xhkb

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= ktop + 1:
                denom: FloatFieldIJ = (
                    xzu.at(K=K - 1)
                    - 0.5 * up_massdetro.at(K=K - 1)
                    + up_massentro.at(K=K - 1)
                )
                if denom == 0.0:
                    xhc = xhc.at(K=K - 1)
                else:
                    xhc = (
                        xhc.at(K=K - 1) * xzu.at(K=K - 1)
                        - 0.5 * up_massdetro.at(K=K - 1) * xhc.at(K=K - 1)
                        + up_massentro.at(K=K - 1) * xhe.at(K=K - 1)
                    ) / denom
                    if K == start_level + 1:
                        x_add: FloatFieldIJ = (
                            cumulus_parameterization_constants.XLV * zqexec
                            + cumulus_parameterization_constants.CP * ztexec
                        ) + x_add_buoy
                        xhc = xhc + x_add * up_massentro.at(K=K - 1) / denom

                xhc = (
                    xhc
                    + cumulus_parameterization_constants.XLF * (1.0 - p_liq_ice) * qrco
                )

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= ktop + 2:
                xhc = xhes_cup
                xzu = 0.0
