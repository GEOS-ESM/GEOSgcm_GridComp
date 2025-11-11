from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntField, IntFieldIJ, Int
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from gt4py.cartesian.gtscript import (
    PARALLEL,
    FORWARD,
    BACKWARD,
    computation,
    interval,
    int32,
    log,
    exp,
    K,
)
from ndsl.dsl.gt4py import function
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import saturation_vapor_pressure
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume


def entrainment_rates(
    vapor: FloatField,
    environment_saturation_mixing_ratio: FloatField,
    lcl_level: IntFieldIJ_Plume,
    error_code: IntFieldIJ_Plume,
    entrainment_rate: FloatField_Plume,
    detrainment_function_updraft: FloatField,
    plume: Int,
):
    from __externals__ import ZERO_DIFF

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            frh = min(
                vapor
                / max(environment_saturation_mixing_ratio, cumulus_parameterization_constants.smaller_qv),
                1.0,
            )
            if K > lcl_level[0, 0][plume]:
                entrainment_rate[0, 0, 0][plume] = (
                    entrainment_rate[0, 0, 0][plume]
                    * (1.3 - frh)
                    * (
                        environment_saturation_mixing_ratio
                        / environment_saturation_mixing_ratio.at(K=lcl_level[0, 0][plume])
                    )
                    ** 3
                )
            else:
                entrainment_rate[0, 0, 0][plume] = entrainment_rate[0, 0, 0][plume] * (1.3 - frh)
            if ZERO_DIFF == 1:
                detrainment_function_updraft = 0.75e-4 * (1.6 - frh)
            else:
                if plume == 0:
                    detrainment_function_updraft = 0.75 * entrainment_rate[0, 0, 0][plume]
                if plume == 1:
                    detrainment_function_updraft = 0.5 * entrainment_rate[0, 0, 0][plume]
                if plume == 2:
                    detrainment_function_updraft = 0.1 * entrainment_rate[0, 0, 0][plume]


def downdraft_entrainment_profiles(
    lateral_entrainment_rate: FloatField,
    entrainment_rate_downdraft: FloatField,
    detrainment_function_downdraft: FloatField,
    scale_dependence_factor_downdraft: FloatFieldIJ,
    plume_entrainment_rate: Float,
):
    """
    Get the entrainment and detrainment profiles for the downdraft

    Args:
        lateral_entrainment_rate (in)
        entrainment_rate_downdraft (out)
        detrainment_function_downdraft (out)
        scale_dependence_factor_downdraft (out)
        plume_entrainment_rate (in)
    """
    from __externals__ import DOWNDRAFT

    with computation(PARALLEL), interval(0, -1):
        entrainment_rate_downdraft = lateral_entrainment_rate * plume_entrainment_rate * 0.3
        detrainment_function_downdraft = entrainment_rate_downdraft

    with computation(FORWARD), interval(0, 1):
        if DOWNDRAFT == 0:
            scale_dependence_factor_downdraft = 1.0
        else:
            scale_dependence_factor_downdraft = 0.0
