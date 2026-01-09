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
    FloatFieldIJ_ensemble_2,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_stencils import unknown_find_level
from ndsl.stencils.column_operations import column_max, column_min, column_max_ddim


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


def downdraft_mass_flux(
    error_code: IntFieldIJ_Plume,
    detrainment_start_level: IntFieldIJ,
    downdraft_origin_level: IntFieldIJ_Plume,
    pbl_level: IntFieldIJ,
    updraft_origin_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    lcl_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    p_surface: FloatFieldIJ,
    normalized_massflux_downdraft: FloatField,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    ocean_fraction: FloatFieldIJ,
    random_number: FloatFieldIJ,
    DOWNDRAFT_MAX_HEIGHT_LAND: Float,
    DOWNDRAFT_MAX_HEIGHT_OCEAN: Float,
    plume: Int,
):
    """
    Handle mass fluxes in the downdraft. For plumes massflux is forced to zero.

    Args:
        error_code
        detrainment_start_level
        downdraft_origin_level
        pbl_level
        updraft_origin_level
        updraft_lfc_level
        lcl_level
        p_cloud_levels_forced
        p_surface
        normalized_massflux_downdraft
        normalized_massflux_downdraft_forced
        ocean_fraction
        random_number
        DOWNDRAFT_MAX_HEIGHT_LAND
        DOWNDRAFT_MAX_HEIGHT_OCEAN
        plume
    """
    from __externals__ import ZERO_DIFF, k_end

    with computation(FORWARD), interval(...):
        normalized_massflux_downdraft = 0.0

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            # set internal constants
            beta: FloatFieldIJ = 2.5

            height_down = (
                1.0 - ocean_fraction
            ) * DOWNDRAFT_MAX_HEIGHT_LAND + ocean_fraction * DOWNDRAFT_MAX_HEIGHT_OCEAN

            # non-zero-diff-APR-08-2020
            if ZERO_DIFF == 1:
                height_down = 0.5
            # non-zero-diff-APR-08-2020

            p_max: FloatFieldIJ = (
                height_down * p_cloud_levels_forced.at(K=downdraft_origin_level[0, 0][plume], ddim=[plume])
                + (1.0 - height_down) * p_surface
            )

    with computation(PARALLEL), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0:
            p_internal = abs(p_cloud_levels_forced[0, 0, 0][plume] - p_max)

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            _, _min_loc = column_min(p_internal, 0, downdraft_origin_level[0, 0][plume])
            min_loc: IntFieldIJ = _min_loc

            # this alpha constrains the location of the maximun ZU to be at "kb_adj" vertical level
            alpha: FloatFieldIJ = 1.0 + (beta - 1.0) * (
                (min_loc + 1) / (downdraft_origin_level[0, 0][plume] + 2)
            ) / (1.0 - ((min_loc + 1) / (downdraft_origin_level[0, 0][plume] + 2)))

    with computation(PARALLEL), interval(1, None):
        if (
            plume != 0
            and error_code[0, 0][plume] == 0
            and K <= min(downdraft_origin_level[0, 0][plume] + 1, k_end - 1)
        ):
            ratio = float(K + 1) / (downdraft_origin_level[0, 0][plume] + 2)
            normalized_massflux_downdraft_forced[0, 0, 0][plume] = ratio ** (alpha - 1.0) * (1.0 - ratio) ** (
                beta - 1.0
            )

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            normalized_massflux_downdraft_forced[0, 0, 0][plume] = 0.0

            # get max value for next block
            _max_val, _ = column_max_ddim(
                normalized_massflux_downdraft_forced,
                plume,
                0,
                min(k_end, downdraft_origin_level[0, 0][plume] + 1),
            )
            max_val: FloatFieldIJ = _max_val

    with computation(FORWARD), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if max_val <= 0:
                normalized_massflux_downdraft_forced[0, 0, 0][plume] = 0.0
                error_code[0, 0][plume] = 51
            else:
                if K <= min(k_end, downdraft_origin_level[0, 0][plume] + 1):
                    normalized_massflux_downdraft_forced[0, 0, 0][plume] = (
                        normalized_massflux_downdraft_forced[0, 0, 0][plume] / (1.0e-9 + max_val)
                    )


def downdraft_lateral_massflux(
    error_code: IntFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    normalized_massflux_downdraft: FloatField,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_modified: FloatField,
    detrainment_function_downdraft: FloatField,
    entrainment_rate_downdraft: FloatField,
    mass_entrainment_downdraft: FloatField,
    mass_detrainment_downdraft: FloatField,
    mass_entrainment_downdraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    mass_entrainment_u_downdraft: FloatField,
    mass_detrainment_u_downdraft: FloatField,
    LAMBDA_DOWN: Float,
    plume: Int,
):
    """
    Get the lateral massfluxes for the downdraft.

    For mid and deep plumes massfluxes are computed, for shallow plumes massfluxes are forced to zero.

    Args:
        error_code
        downdraft_origin_level
        geopotential_height_cloud_levels_forced
        normalized_massflux_downdraft
        normalized_massflux_downdraft_forced
        normalized_massflux_downdraft_modified
        detrainment_function_downdraft
        entrainment_rate_downdraft
        mass_entrainment_downdraft
        mass_detrainment_downdraft
        mass_entrainment_downdraft_forced
        mass_detrainment_downdraft_forced
        mass_entrainment_u_downdraft
        mass_detrainment_u_downdraft
        plume
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        # set entrainment/detrainment to zero
        detrainment_function_downdraft = 0.0
        mass_entrainment_downdraft = 0.0
        mass_detrainment_downdraft = 0.0
        mass_entrainment_downdraft_forced[0, 0, 0][plume] = 0.0
        mass_detrainment_downdraft_forced[0, 0, 0][plume] = 0.0
        mass_entrainment_u_downdraft = 0.0
        mass_detrainment_u_downdraft = 0.0

    # rest of this stencil does not execute for shallow plumes (plume = 0)
    with computation(PARALLEL), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0 and K < downdraft_origin_level[0, 0][plume]:
            detrainment_function_downdraft = entrainment_rate_downdraft

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            entrainment_rate_downdraft = 0.0

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            # get maximum index for next block
            _, _max_loc = column_max_ddim(normalized_massflux_downdraft_forced, plume, 0, k_end)
            max_loc: IntFieldIJ = _max_loc

    with computation(BACKWARD), interval(0, -1):
        if (
            plume != 0
            and error_code[0, 0][plume] == 0
            and K >= max_loc
            and K <= downdraft_origin_level[0, 0][plume]
        ):

            # from downdraft_origin_level to maximum value of normalized_massflux_downdraft, change entrainment
            dzo = geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                detrainment_function_downdraft * dzo * normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                normalized_massflux_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 1][plume]
                + mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = max(
                0.0, mass_entrainment_downdraft_forced[0, 0, 0][plume]
            )
            # check dd_massdetro in case of dd_massentro has been changed above
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                mass_entrainment_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 0][plume]
                + normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

    with computation(BACKWARD), interval(0, -1):
        if plume != 0 and error_code[0, 0][plume] == 0 and K < max_loc:
            # from maximum value zd to surface -> change detrainment
            dzo = geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                entrainment_rate_downdraft * dzo * normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                normalized_massflux_downdraft_forced[0, 0, 1][plume]
                + mass_entrainment_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 0][plume]
            )
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = max(
                0.0, mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )
            # check dd_massentro in case of dd_massdetro has been changed above
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                mass_detrainment_downdraft_forced[0, 0, 0][plume]
                + normalized_massflux_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

    with computation(BACKWARD), interval(0, -1):
        if plume != 0 and error_code[0, 0][plume] == 0 and K <= downdraft_origin_level[0, 0][plume]:
            normalized_massflux_downdraft_modified = normalized_massflux_downdraft_forced[0, 0, 0][plume]
            normalized_massflux_downdraft = normalized_massflux_downdraft_forced[0, 0, 0][plume]
            mass_entrainment_downdraft = mass_entrainment_downdraft_forced[0, 0, 0][plume]
            mass_detrainment_downdraft = mass_detrainment_downdraft_forced[0, 0, 0][plume]
            mass_entrainment_u_downdraft = (
                mass_entrainment_downdraft_forced[0, 0, 0][plume]
                + LAMBDA_DOWN * mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )
            mass_detrainment_u_downdraft = (
                mass_detrainment_downdraft_forced[0, 0, 0][plume]
                + LAMBDA_DOWN * mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )


def downdraft_moist_static_energy_and_buoyancy(
    error_code: IntFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ_Plume,
    u: FloatField,
    u_cloud_levels: FloatField,
    u_c_downdraft: FloatField,
    v: FloatField,
    v_cloud_levels: FloatField,
    v_c_downdraft: FloatField,
    environment_moist_static_energy_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    cloud_moist_static_energy: FloatField,
    cloud_moist_static_energy_downdraft_forced: FloatField,
    buoyancy_downdraft_forced: FloatField,
    t_wetbulb: FloatFieldIJ,
    vapor_wetbulb: FloatFieldIJ,
    geopotential_height_cloud_levels_forced: FloatField,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    mass_entrainment_downdraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    mass_entrainment_u_downdraft: FloatField,
    mass_detrainment_u_downdraft: FloatField,
    plume: Int,
):
    """
    Compute moist static energy and buoyancy for the downdraft.

    For shallow plumes: majority of the code is skipped. buoyancy_downdraft_forced is forced to zero,
    u_c_downdraft and v_c_downdraft are forced to u_cloud_levels and v_cloud_levels, respectively, and
    cloud_moist_static_energy_downdraft_forced is forced to environment_saturation_moist_static_energy.

    Args:
        error_code
        downdraft_origin_level
        u
        u_cloud_levels
        u_c_downdraft
        v
        v_cloud_levels
        v_c_downdraft
        environment_moist_static_energy_forced
        environment_saturation_moist_static_energy_cloud_levels_forced
        cloud_moist_static_energy
        cloud_moist_static_energy_downdraft_forced
        buoyancy_downdraft_forced
        t_wetbulb
        vapor_wetbulb
        geopotential_height_cloud_levels_forced
        normalized_massflux_downdraft_forced
        mass_entrainment_downdraft_forced
        mass_detrainment_downdraft_forced
        mass_entrainment_u_downdraft
        mass_detrainment_u_downdraft
        plume
    """
    from __externals__ import USE_WETBULB, PRESSURE_GRADIENT_CONSTANT

    with computation(PARALLEL), interval(...):
        cloud_moist_static_energy_downdraft_forced = (
            environment_saturation_moist_static_energy_cloud_levels_forced
        )
        u_c_downdraft = u_cloud_levels
        v_c_downdraft = v_cloud_levels
        buoyancy_downdraft_forced = 0.0
        buoyancy_downdraft = 0.0

    with computation(FORWARD), interval(0, 1):
        # set bounds for next block
        lower_bound: IntFieldIJ = downdraft_origin_level[0, 0][plume]
        upper_bound: IntFieldIJ = downdraft_origin_level[0, 0][plume] + 1

    with computation(FORWARD), interval(lower_bound, upper_bound):
        buoyancy_downdraft: FloatFieldIJ = 0.0
        if error_code[0, 0][plume] == 0 and plume != 0:
            wetbulb_adjustment: IntFieldIJ = 0
            # for future test)
            if USE_WETBULB == 1:
                cloud_moist_static_energy_downdraft_forced = 0.5 * (
                    cumulus_parameterization_constants.CP * t_wetbulb
                    + cumulus_parameterization_constants.XLV * vapor_wetbulb
                    + geopotential_height_cloud_levels_forced * constants.MAPL_GRAV
                    + cloud_moist_static_energy
                )
                wetbulb_adjustment = 1

            buoyancy_downdraft_forced = (
                cloud_moist_static_energy_downdraft_forced
                - environment_saturation_moist_static_energy_cloud_levels_forced
            )
            buoyancy_downdraft = buoyancy_downdraft_forced * (
                geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
            )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0:
            upper_bound = downdraft_origin_level[0, 0][plume]

    with computation(BACKWARD), interval(0, upper_bound):
        if error_code[0, 0][plume] == 0 and plume != 0:
            denom = (
                normalized_massflux_downdraft_forced[0, 0, 1][plume]
                - 0.5 * mass_detrainment_downdraft_forced[0, 0, 0][plume]
                + mass_entrainment_downdraft_forced[0, 0, 0][plume]
            )
            denom_u = (
                normalized_massflux_downdraft_forced[0, 0, 1][plume]
                - 0.5 * mass_detrainment_u_downdraft
                + mass_entrainment_u_downdraft
            )

            # tmp fix for denominator being zero
            if denom > 0.0 and denom_u > 0.0:
                dz = (
                    geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
                )

                u_c_downdraft = (
                    u_c_downdraft[0, 0, 1] * normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    - 0.5 * mass_detrainment_u_downdraft * u_c_downdraft[0, 0, 1]
                    + mass_entrainment_u_downdraft * u
                    - PRESSURE_GRADIENT_CONSTANT
                    * normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    * (u[0, 0, 1] - u)
                ) / denom_u
                v_c_downdraft = (
                    v_c_downdraft[0, 0, 1] * normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    - 0.5 * mass_detrainment_u_downdraft * v_c_downdraft[0, 0, 1]
                    + mass_entrainment_u_downdraft * v
                    - PRESSURE_GRADIENT_CONSTANT
                    * normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    * (v[0, 0, 1] - v)
                ) / denom_u

                cloud_moist_static_energy_downdraft_forced = (
                    cloud_moist_static_energy_downdraft_forced[0, 0, 1]
                    * normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    - 0.5
                    * mass_detrainment_downdraft_forced[0, 0, 0][plume]
                    * cloud_moist_static_energy_downdraft_forced[0, 0, 1]
                    + mass_entrainment_downdraft_forced[0, 0, 0][plume]
                    * environment_moist_static_energy_forced
                ) / denom

                buoyancy_downdraft_forced = (
                    cloud_moist_static_energy_downdraft_forced
                    - environment_saturation_moist_static_energy_cloud_levels_forced
                )
                buoyancy_downdraft = buoyancy_downdraft + buoyancy_downdraft_forced * dz
            else:
                u_c_downdraft = u_c_downdraft[0, 0, 1]
                v_c_downdraft = v_c_downdraft[0, 0, 1]
                cloud_moist_static_energy_downdraft_forced = cloud_moist_static_energy_downdraft_forced[
                    0, 0, 1
                ]

    with computation(FORWARD), interval(0, 1):
        if buoyancy_downdraft > 0:
            error_code[0, 0][plume] = 7


def downdraft_moisture(
    error_code: IntFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ_Plume,
    t_cloud_levels_forced: FloatField,
    t_wetbulb: FloatFieldIJ,
    vapor_forced: FloatField,
    vapor_cloud_levels_forced: FloatField,
    environment_saturation_mixing_ratio_cloud_levels_forced: FloatField,
    cloud_total_water_after_entrainment_forced: FloatField,
    cloud_total_water_after_entrainment_downdraft_forced: FloatField,
    downdraft_saturation_vapor_forced: FloatField,
    vapor_wetbulb: FloatFieldIJ,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    environment_moist_static_energy_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    cloud_moist_static_energy_downdraft_forced: FloatField,
    evaporate_in_downdraft_forced: FloatField_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    mass_entrainment_downdraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    gamma_cloud_levels_forced: FloatField,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    total_normalized_integrated_evaporate_forced: FloatFieldIJ_Plume,
    buoyancy: FloatFieldIJ,
    plume: Int,
):
    """
    Compute the moisture profile for the downdraft.

    For shallow plumes outputs are forced to zero and calculation is terminated.

    Args:
        error_code
        downdraft_origin_level
        t_cloud_levels_forced
        t_wetbulb
        vapor_forced
        vapor_cloud_levels_forced
        environment_saturation_mixing_ratio_cloud_levels_forced
        cloud_total_water_after_entrainment_forced
        cloud_total_water_after_entrainment_downdraft_forced
        downdraft_saturation_vapor_forced
        vapor_wetbulb
        normalized_massflux_downdraft_forced
        environment_moist_static_energy_forced
        environment_saturation_moist_static_energy_cloud_levels_forced
        cloud_moist_static_energy_downdraft_forced
        evaporate_in_downdraft_forced
        geopotential_height_cloud_levels_forced
        mass_entrainment_downdraft_forced
        mass_detrainment_downdraft_forced
        gamma_cloud_levels_forced
        total_normalized_integrated_condensate_forced
        total_normalized_integrated_evaporate_forced
        buoyancy
        plume
    """
    from __externals__ import USE_WETBULB, ZERO_DIFF, EVAP_FIX

    with computation(FORWARD), interval(0, 1):
        internal_loop_constant: IntFieldIJ = 1

    with computation(FORWARD), interval(0, 1):
        # prefill outputs with zero
        buoyancy = 0.0
        total_normalized_integrated_evaporate_forced[0, 0][plume] = 0.0

    with computation(PARALLEL), interval(...):
        # prefill outputs with zero
        cloud_total_water_after_entrainment_downdraft_forced = 0.0
        downdraft_saturation_vapor_forced = 0.0
        evaporate_in_downdraft_forced[0, 0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0:
            # set bounds for next block
            lower_bound: IntFieldIJ = downdraft_origin_level[0, 0][plume]
            upper_bound: IntFieldIJ = downdraft_origin_level[0, 0][plume] + 1

    with computation(FORWARD), interval(lower_bound, upper_bound):
        if error_code[0, 0][plume] == 0 and plume != 0:
            # boundary condition at downdraft_origin_level ('level of free sinking')
            dz = geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced

            cloud_total_water_after_entrainment_downdraft_forced = vapor_cloud_levels_forced

            if USE_WETBULB == 1:
                # mixture 50% env air + updraft
                cloud_total_water_after_entrainment_downdraft_forced = 0.5 * (
                    vapor_wetbulb + cloud_total_water_after_entrainment_forced
                )

            d_moist_static_energy = (
                cloud_moist_static_energy_downdraft_forced
                - environment_saturation_moist_static_energy_cloud_levels_forced
            )

            if d_moist_static_energy < 0:
                downdraft_saturation_vapor_forced = (
                    environment_saturation_mixing_ratio_cloud_levels_forced
                    + (1.0 / cumulus_parameterization_constants.XLV)
                    * (gamma_cloud_levels_forced / (1.0 + gamma_cloud_levels_forced))
                    * d_moist_static_energy
                )
            else:
                downdraft_saturation_vapor_forced = environment_saturation_mixing_ratio_cloud_levels_forced

            evaporate_in_downdraft_forced[0, 0, 0][plume] = normalized_massflux_downdraft_forced[0, 0, 0][
                plume
            ] * min(
                0.0, cloud_total_water_after_entrainment_downdraft_forced - downdraft_saturation_vapor_forced
            )
            cloud_total_water_after_entrainment_downdraft_forced = downdraft_saturation_vapor_forced
            total_normalized_integrated_evaporate_forced[0, 0][plume] = (
                total_normalized_integrated_evaporate_forced[0, 0][plume]
                + evaporate_in_downdraft_forced[0, 0, 0][plume]
            )
            buoyancy = dz * d_moist_static_energy

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0:
            # set bounds for next block
            upper_bound = downdraft_origin_level[0, 0][plume]

    with computation(BACKWARD), interval(0, upper_bound):
        if error_code[0, 0][plume] == 0 and plume != 0:
            dz = geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced

            # downward transport + mixing
            denom = (
                normalized_massflux_downdraft_forced[0, 0, 1][plume]
                - 0.5 * mass_detrainment_downdraft_forced[0, 0, 0][plume]
                + mass_entrainment_downdraft_forced[0, 0, 0][plume]
            )
            if denom == 0.0:
                cloud_total_water_after_entrainment_downdraft_forced = (
                    cloud_total_water_after_entrainment_downdraft_forced[0, 0, 1]
                )
            else:
                cloud_total_water_after_entrainment_downdraft_forced = (
                    cloud_total_water_after_entrainment_downdraft_forced[0, 0, 1]
                    * normalized_massflux_downdraft_forced[0, 0, 1][plume]
                    - 0.5
                    * mass_detrainment_downdraft_forced[0, 0, 0][plume]
                    * cloud_total_water_after_entrainment_downdraft_forced[0, 0, 1]
                    + mass_entrainment_downdraft_forced[0, 0, 0][plume] * vapor_forced
                ) / denom

            # to be negatively buoyant, hcd should be smaller than hes!
            # ideally, dh should be negative till dd hits ground, but that is not always the case
            d_moist_static_energy = (
                cloud_moist_static_energy_downdraft_forced
                - environment_saturation_moist_static_energy_cloud_levels_forced
            )
            buoyancy = buoyancy + dz * d_moist_static_energy
            downdraft_saturation_vapor_forced = (
                environment_saturation_mixing_ratio_cloud_levels_forced
                + (1.0 / cumulus_parameterization_constants.XLV)
                * (gamma_cloud_levels_forced / (1.0 + gamma_cloud_levels_forced))
                * d_moist_static_energy
            )

            # rain water evaporation amount at layer k
            dq_eva = cloud_total_water_after_entrainment_downdraft_forced - downdraft_saturation_vapor_forced

            if dq_eva > 0.0:
                dq_eva = 0.0
                downdraft_saturation_vapor_forced = cloud_total_water_after_entrainment_downdraft_forced
            # amount of the evaporated rain water
            evaporate_in_downdraft_forced[0, 0, 0][plume] = (
                normalized_massflux_downdraft_forced[0, 0, 0][plume] * dq_eva
            )  # kg[water vapor]/kg[air]

            # source term for in-downdraft water vapor mixing ratio
            cloud_total_water_after_entrainment_downdraft_forced = downdraft_saturation_vapor_forced  # => equiv to qcd = qcd - dq_eva !( -dq_eva >0 => source term for qcd)

            # total evaporated rain water
            total_normalized_integrated_evaporate_forced[0, 0][plume] = (
                total_normalized_integrated_evaporate_forced[0, 0][plume]
                + evaporate_in_downdraft_forced[0, 0, 0][plume]
            )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0:
            if total_normalized_integrated_evaporate_forced[0, 0][plume] >= 0 and internal_loop_constant == 1:
                error_code[0, 0][plume] = 70

            if buoyancy >= 0 and internal_loop_constant == 1:
                error_code[0, 0][plume] = 73

            if ZERO_DIFF == 0 and EVAP_FIX == 1:
                if (
                    abs(total_normalized_integrated_evaporate_forced[0, 0][plume])
                    > total_normalized_integrated_condensate_forced[0, 0][plume]
                    and error_code[0, 0][plume] == 0
                ):
                    fix_evap = total_normalized_integrated_condensate_forced[0, 0][plume] / (
                        1.0e-16 + abs(total_normalized_integrated_evaporate_forced[0, 0][plume])
                    )
                    total_normalized_integrated_evaporate_forced[0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0:
            # set bounds for next block
            upper_bound = downdraft_origin_level[0, 0][plume] + 1

    with computation(BACKWARD), interval(0, upper_bound):
        if error_code[0, 0][plume] == 0 and plume != 0:
            if ZERO_DIFF == 0 and EVAP_FIX == 1:
                evaporate_in_downdraft_forced[0, 0, 0][plume] = (
                    evaporate_in_downdraft_forced[0, 0, 0][plume] * fix_evap
                )
                total_normalized_integrated_evaporate_forced[0, 0][plume] = (
                    total_normalized_integrated_evaporate_forced[0, 0][plume]
                    + evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                dq_eva = evaporate_in_downdraft_forced[0, 0, 0][plume] / (
                    1.0e-16 + normalized_massflux_downdraft_forced[0, 0, 0][plume]
                )
                cloud_total_water_after_entrainment_downdraft_forced = (
                    downdraft_saturation_vapor_forced + dq_eva
                )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0:
            if ZERO_DIFF == 0 and EVAP_FIX == 1:
                if total_normalized_integrated_evaporate_forced[0, 0][plume] >= 0.0:
                    error_code[0, 0][plume] = 70


def downdraft_temperature(
    error_code: IntFieldIJ_Plume,
    downdraft_column_temperature_forced: FloatField,
    cloud_moist_static_energy_downdraft_forced: FloatField,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_total_water_after_entrainment_downdraft_forced: FloatField,
    t_cloud_levels_forced: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            downdraft_column_temperature_forced = (1.0 / cumulus_parameterization_constants.CP) * (
                cloud_moist_static_energy_downdraft_forced
                - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
                - cumulus_parameterization_constants.XLV
                * cloud_total_water_after_entrainment_downdraft_forced
            )
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            downdraft_column_temperature_forced = t_cloud_levels_forced


def downdraft_windshear(
    error_code: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    geopotential_height_forced: FloatField,
    p_forced: FloatField,
    u: FloatField,
    v: FloatField,
    ccn: FloatFieldIJ,
    psum: FloatFieldIJ,
    psumh: FloatFieldIJ,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    total_normalized_integrated_evaporate_forced: FloatFieldIJ_Plume,
    epsilon: FloatFieldIJ,
    epsilon_min: FloatFieldIJ,
    epsilon_max: FloatFieldIJ,
    epsilon_computed: FloatFieldIJ_ensemble_2,
    plume: Int,
):
    from __externals__ import AEROEVAP

    with computation(FORWARD), interval(0, 1):
        # initalize internal constants
        alpha3: FloatFieldIJ = 1.9
        beta3: FloatFieldIJ = -1.13

        # zero input fields
        epsilon = 0.0

        # zero every part of the ensemble dimension
        count = 0
        while count < cumulus_parameterization_constants.MAXENS2:
            epsilon_computed[0, 0][count] = 0.0
            count += 1

        # initalize internal 2d temporaries
        vshear: FloatFieldIJ = 0.0
        sdp: FloatFieldIJ = 0.0
        vws: FloatFieldIJ = 0.0

        lower_bound: FloatFieldIJ = updraft_lfc_level[0, 0][plume]
        upper_bound: FloatFieldIJ = cloud_top_level[0, 0][plume] + 1

    with computation(FORWARD), interval(lower_bound, upper_bound):
        if plume != 0 and error_code[0, 0][plume] == 0:
            dp = p_forced - p_forced[0, 0, 1]
            vws = (
                vws
                + (
                    abs((u[0, 0, 1] - u) / (geopotential_height_forced[0, 0, 1] - geopotential_height_forced))
                    + abs(
                        (v[0, 0, 1] - v) / (geopotential_height_forced[0, 0, 1] - geopotential_height_forced)
                    )
                )
                * dp
            )
            sdp = sdp + dp

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            vshear = 1.0e3 * vws / sdp

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            precip_efficiency = 1.591 - 0.639 * vshear + 0.0953 * (vshear**2) - 0.00496 * (vshear**3)
            precip_efficiency = min(precip_efficiency, 0.9)
            precip_efficiency = max(precip_efficiency, 0.1)

            # cloud base precip efficiency
            zkbc = geopotential_height_forced.at(K=updraft_lfc_level[0, 0][plume]) * 3.281e-3
            prezk = 0.02
            if zkbc > 3.0:
                prezk = 0.96729352 + zkbc * (
                    -0.70034167
                    + zkbc * (0.162179896 + zkbc * (-1.2569798e-2 + zkbc * (4.2772e-4 - zkbc * 5.44e-6)))
                )
            if zkbc > 25.0:
                prezk = 2.4

            precip_efficiency_b = 1.0 / (1.0 + prezk)
            precip_efficiency_b = min(precip_efficiency_b, 0.9)
            precip_efficiency_b = max(precip_efficiency_b, 0.1)

            epsilon = 1.0 - 0.5 * (precip_efficiency_b + precip_efficiency)

            if AEROEVAP > 1:
                aeroadd = (cumulus_parameterization_constants.CCNCLEAN**beta3) * (
                    (psumh) ** (alpha3 - 1)
                )  # *1.e6
                prop_c = 0.5 * (precip_efficiency_b + precip_efficiency) / aeroadd
                aeroadd = (ccn**beta3) * ((psum) ** (alpha3 - 1))  # *1.e6
                aeroadd = prop_c * aeroadd
                precip_efficiency_c = aeroadd
                if precip_efficiency_c > 0.9:
                    precip_efficiency_c = 0.9
                if precip_efficiency_c < 0.1:
                    precip_efficiency_c = 0.1
                epsilon = 1.0 - precip_efficiency_c
                if AEROEVAP == 2:
                    epsilon = 1.0 - 0.25 * (
                        precip_efficiency_b + precip_efficiency + 2.0 * precip_efficiency_c
                    )

            # epsilon here is 1-precip_efficiency!
            einc = 0.2 * epsilon
            count = 0
            while count < cumulus_parameterization_constants.MAXENS2:
                epsilon_computed[0, 0][count] = epsilon + (count - 1) * einc
                count += 1

            count = 0
            while count < cumulus_parameterization_constants.MAXENS2:
                epsilon_computed[0, 0][count] = (
                    -epsilon_computed[0, 0][count]
                    * total_normalized_integrated_condensate_forced[0, 0][plume]
                    / total_normalized_integrated_evaporate_forced[0, 0][plume]
                )
                if epsilon_computed[0, 0][count] > epsilon_max:
                    epsilon_computed[0, 0][count] = epsilon_max
                if epsilon_computed[0, 0][count] < epsilon_min:
                    epsilon_computed[0, 0][count] = epsilon_min
                count += 1


def update_epsilon_forced(
    error_code: IntFieldIJ_Plume,
    scale_dependence_factor_downdraft: FloatFieldIJ,
    epsilon_computed: FloatFieldIJ_ensemble_2,
    epsilon: FloatFieldIJ,
    epsilon_forced: FloatFieldIJ_Plume,
    plume: Int,
):
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            count = 0
            while count < cumulus_parameterization_constants.MAXENS2:
                epsilon_forced[0, 0][plume] = (
                    scale_dependence_factor_downdraft * epsilon_computed[0, 0][count]
                )
                epsilon = epsilon_forced[0, 0][plume]
                count += 1


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


class DowndraftWetBlub:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftWindShear:
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
        self._downdraft_windshear = stencil_factory.from_dims_halo(
            func=downdraft_windshear,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"AEROEVAP": cumulus_parameterization_config.AEROEVAP},
        )

        self._update_epsilon_forced = stencil_factory.from_dims_halo(
            func=update_epsilon_forced,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        error_code: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        geopotential_height_forced: Quantity,
        p_forced: Quantity,
        u: Quantity,
        v: Quantity,
        ccn: Quantity,
        psum: Quantity,
        psumh: Quantity,
        total_normalized_integrated_condensate_forced: Quantity,
        total_normalized_integrated_evaporate_forced: Quantity,
        scale_dependence_factor_downdraft: Quantity,
        epsilon: Quantity,
        epsilon_min: Quantity,
        epsilon_max: Quantity,
        epsilon_computed: Quantity,
        epsilon_forced: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._downdraft_windshear(
            error_code=error_code,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            geopotential_height_forced=geopotential_height_forced,
            p_forced=p_forced,
            u=u,
            v=v,
            ccn=ccn,
            psum=psum,
            psumh=psumh,
            total_normalized_integrated_condensate_forced=total_normalized_integrated_condensate_forced,
            total_normalized_integrated_evaporate_forced=total_normalized_integrated_evaporate_forced,
            epsilon=epsilon,
            epsilon_min=epsilon_min,
            epsilon_max=epsilon_max,
            epsilon_computed=epsilon_computed,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._update_epsilon_forced(
            error_code=error_code,
            scale_dependence_factor_downdraft=scale_dependence_factor_downdraft,
            epsilon_computed=epsilon_computed,
            epsilon=epsilon,
            epsilon_forced=epsilon_forced,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


######## NOTE TODO NOTE README NOTE TODO TODO NOTE EVERYTHING BELOW HERE NEEDS TO BE REWORKED


class DowndraftLateralMassFlux:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


class DowndraftMoistureProperties:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
