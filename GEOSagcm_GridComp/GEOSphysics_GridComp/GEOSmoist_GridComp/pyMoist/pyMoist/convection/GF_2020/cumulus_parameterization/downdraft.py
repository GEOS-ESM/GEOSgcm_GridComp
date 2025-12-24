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
from ndsl.stencils.column_operations import column_max, column_max_ddim


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


def downdraft_lateral_massflux(
    error_code: IntFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ,
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
        # zero out entrainment/detrainment
        detrainment_function_downdraft = 0.0
        mass_entrainment_downdraft = 0.0
        mass_detrainment_downdraft = 0.0
        mass_entrainment_downdraft_forced[0, 0, 0][plume] = 0.0
        mass_detrainment_downdraft_forced[0, 0, 0][plume] = 0.0
        mass_entrainment_u_downdraft = 0.0
        mass_detrainment_u_downdraft = 0.0

    with computation(PARALLEL), interval(0, downdraft_origin_level - 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            detrainment_function_downdraft = entrainment_rate_downdraft

    with computation(FORWARD), interval(0, 1):
        entrainment_rate_downdraft = 0.0

        # get max_loc for next block
        _, _max_loc = column_max_ddim(normalized_massflux_downdraft_forced, plume, 0, k_end)
        max_loc: IntFieldIJ = _max_loc

    with computation(BACKWARD), interval(max_loc, downdraft_origin_level):
        if error_code[0, 0][plume] == 0:
            # from downdraft_origin_level to maximum value of
            # normalized_massflux_downdraft_forced -> change entrainment
            normalized_massflux_downdraft_forced[0, 0, 0][plume] = (
                geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
            )
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                detrainment_function_downdraft
                * normalized_massflux_downdraft_forced[0, 0, 0][plume]
                * normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                normalized_massflux_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 1][plume]
                + mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = max(
                0.0, mass_entrainment_downdraft_forced[0, 0, 0][plume]
            )

            # check mass_detrainment_downdraft_forced in case
            # mass_entrainment_downdraft_forced has been changed above
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                mass_entrainment_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 0][plume]
                + normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # get max_loc for next block
            _, _max_loc = column_max_ddim(normalized_massflux_downdraft_forced, plume, 0, k_end)
            max_loc: IntFieldIJ = _max_loc

    with computation(BACKWARD), interval(0, max_loc):
        if error_code[0, 0][plume] == 0:
            # from maximum value normalized_massflux_downdraft_forced to surface -> change detrainment
            dz = geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                entrainment_rate_downdraft * dz * normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

            mass_detrainment_downdraft_forced[0, 0, 0][plume] = (
                normalized_massflux_downdraft_forced[0, 0, 1][plume]
                + mass_entrainment_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 0][plume]
            )
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = max(
                0.0, mass_detrainment_downdraft_forced[0, 0, 0][plume]
            )
            # check mass_entrainment_downdraft_forced in case
            # mass_detrainment_downdraft_forced has been changed above
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = (
                mass_detrainment_downdraft_forced[0, 0, 0][plume]
                + normalized_massflux_downdraft_forced[0, 0, 0][plume]
                - normalized_massflux_downdraft_forced[0, 0, 1][plume]
            )

    with computation(BACKWARD), interval(0, downdraft_origin_level):
        if error_code[0, 0][plume] == 0:
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
    downdraft_origin_level: IntFieldIJ,
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

    with computation(FORWARD), interval(downdraft_origin_level, downdraft_origin_level + 1):
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

    with computation(BACKWARD), interval(0, downdraft_origin_level):
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


######## NOTE TODO NOTE README NOTE TODO TODO NOTE EVERYTHING BELOW HERE NEEDS TO BE REWORKED


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
