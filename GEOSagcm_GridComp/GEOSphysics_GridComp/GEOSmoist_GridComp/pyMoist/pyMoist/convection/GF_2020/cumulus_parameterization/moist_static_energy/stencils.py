from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import get_cloud_boundary_conditions


def parcel_moist_static_energy(
    error_code: IntFieldIJ_Plume,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    updraft_origin_level: FloatFieldIJ,
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

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            modification = (
                cumulus_parameterization_constants.XLV * vapor_excess
                + cumulus_parameterization_constants.CP
                + t_excess
            ) + add_buoyancy

        moist_static_energy_origin_level = get_cloud_boundary_conditions(
            field=environmenet_moist_static_energy,
            scalar_perturbation=modification,
            p=p,
            updraft_origin_level=updraft_origin_level,
            ocean_fraction=ocean_fraction,
            BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
            AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
            k_end=k_end,
            compute_perturbation=True,
            perturbation_field=t_perturbation,
        )
        moist_static_energy_origin_level_forced = get_cloud_boundary_conditions(
            field=environmenet_moist_static_energy,
            scalar_perturbation=modification,
            p=p,
            updraft_origin_level=updraft_origin_level,
            ocean_fraction=ocean_fraction,
            BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
            AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
            k_end=k_end,
            compute_perturbation=True,
            perturbation_field=t_perturbation,
        )
