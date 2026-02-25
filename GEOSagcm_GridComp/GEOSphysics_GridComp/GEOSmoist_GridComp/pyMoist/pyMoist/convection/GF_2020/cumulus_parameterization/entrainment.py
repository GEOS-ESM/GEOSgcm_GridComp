from gt4py.cartesian.gtscript import FORWARD, PARALLEL, K, computation, interval

import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, Int, IntFieldIJ
from ndsl.stencils.column_operations import column_max_ddim
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import get_cloud_boundary_conditions


def entrainment_rates(
    vapor_cloud_levels_forced: FloatField,
    environment_saturation_mixing_ratio_cloud_levels_forced: FloatField,
    lcl_level: IntFieldIJ_Plume,
    error_code: IntFieldIJ_Plume,
    entrainment_rate: FloatField_Plume,
    detrainment_function_updraft: FloatField,
    plume: Int,
):
    """Determine the entrainment/detrainment rates, updating the initial estimate
    using data computed inside the core.

    Args:
        vapor_cloud_levels_forced (FloatField)
        environment_saturation_mixing_ratio_cloud_levels_forced (FloatField)
        lcl_level (IntFieldIJ_Plume)
        error_code (IntFieldIJ_Plume)
        entrainment_rate (FloatField_Plume)
        detrainment_function_updraft (FloatField)
        plume (Int)
    """
    from __externals__ import MIN_ENTRAINMENT_RATE, ZERO_DIFF

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            frh = min(
                vapor_cloud_levels_forced
                / max(
                    environment_saturation_mixing_ratio_cloud_levels_forced,
                    cumulus_parameterization_constants.smaller_qv,
                ),
                1.0,
            )
            if K >= lcl_level[0, 0][plume]:
                entrainment_rate[0, 0, 0][plume] = (
                    entrainment_rate[0, 0, 0][plume]
                    * (1.3 - frh)
                    * (
                        environment_saturation_mixing_ratio_cloud_levels_forced
                        / environment_saturation_mixing_ratio_cloud_levels_forced.at(K=lcl_level[0, 0][plume])
                    )
                    ** 3
                )
            else:
                entrainment_rate[0, 0, 0][plume] = entrainment_rate[0, 0, 0][plume] * (1.3 - frh)
            if ZERO_DIFF == 1:
                detrainment_function_updraft = 0.75e-4 * (1.6 - frh)
            else:
                entrainment_rate[0, 0, 0][plume] = max(entrainment_rate[0, 0, 0][plume], MIN_ENTRAINMENT_RATE)
                if plume == cumulus_parameterization_constants.SHALLOW:
                    detrainment_function_updraft = 0.75 * entrainment_rate[0, 0, 0][plume]
                if plume == cumulus_parameterization_constants.MID:
                    detrainment_function_updraft = 0.5 * entrainment_rate[0, 0, 0][plume]
                if plume == cumulus_parameterization_constants.DEEP:
                    detrainment_function_updraft = 0.1 * entrainment_rate[0, 0, 0][plume]


def downdraft_entrainment_profiles(
    lateral_entrainment_rate: FloatField,
    entrainment_rate_downdraft: FloatField,
    detrainment_function_downdraft: FloatField,
    scale_dependence_factor_downdraft: FloatFieldIJ,
    plume_entrainment_rate: Float,
):
    """Get the entrainment and detrainment profiles for the downdraft

    Args:
        lateral_entrainment_rate (FloatField)
        entrainment_rate_downdraft (FloatField)
        detrainment_function_downdraft (FloatField)
        scale_dependence_factor_downdraft (FloatFieldIJ)
        plume_entrainment_rate (Float)
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


def compute_lateral_massflux(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    geopotential_height: FloatField,
    normalized_massflux_updraft: FloatField_Plume,
    detrainment_function_updraft: FloatField,
    entrainment_rate: FloatField_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_updraft: FloatField,
    mass_detrainment_updraft: FloatField,
    updraft_lfc_level: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    pbl_level: IntFieldIJ,
    mass_entrainment_u_updraft: FloatField,
    mass_detrainment_u_updraft: FloatField,
    LAMBDA_DEEP: Float,
    plume: Int,
):
    """Compute massfluxes associated with entrainment and detrainment.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        geopotential_height (FloatField)
        normalized_massflux_updraft (FloatField_Plume)
        detrainment_function_updraft (FloatField)
        entrainment_rate (FloatField_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        mass_entrainment_updraft_forced (FloatField_Plume)
        mass_detrainment_updraft_forced (FloatField_Plume)
        mass_entrainment_updraft (FloatField)
        mass_detrainment_updraft (FloatField)
        updraft_lfc_level (IntFieldIJ_Plume)
        updraft_origin_level (IntFieldIJ_Plume)
        pbl_level (IntFieldIJ)
        mass_entrainment_u_updraft (FloatField)
        mass_detrainment_u_updraft (FloatField)
        LAMBDA_DEEP (Float)
        plume (Int)
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        mass_entrainment_updraft_forced[0, 0, 0][plume] = 0.0
        mass_detrainment_updraft_forced[0, 0, 0][plume] = 0.0
        mass_entrainment_u_updraft = 0.0
        mass_detrainment_u_updraft = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            _, index = column_max_ddim(normalized_massflux_updraft, plume, 0, k_end)
            max_index: IntFieldIJ = index

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            # will not allow detrainment below cloud base or in the PBL
            if plume == cumulus_parameterization_constants.SHALLOW:
                if (
                    K
                    <= max(
                        updraft_lfc_level[0, 0][plume],
                        updraft_origin_level[0, 0][plume],
                    )
                    + 1
                ):
                    detrainment_function_updraft = 0
            else:
                if K <= max_index + 1:
                    detrainment_function_updraft = 0

    # mass entrainment and detrainment are defined on model levels
    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            # below location of maximum value normalized_massflux_updraft -> change entrainment
            if K <= max_index - 1:
                height_massflux_avg = (
                    geopotential_height[0, 0, 1] - geopotential_height
                ) * normalized_massflux_updraft[0, 0, 0][plume]
                mass_detrainment_updraft_forced[0, 0, 0][plume] = (
                    detrainment_function_updraft * height_massflux_avg
                )
                mass_entrainment_updraft_forced[0, 0, 0][plume] = (
                    normalized_massflux_updraft[0, 0, 1][plume]
                    + -normalized_massflux_updraft[0, 0, 0][plume]
                    + mass_detrainment_updraft_forced[0, 0, 0][plume]
                )
                mass_entrainment_updraft_forced[0, 0, 0][plume] = max(
                    mass_entrainment_updraft_forced[0, 0, 0][plume], 0.0
                )

                # check mass_detrainment_updraft_forced in case it has been changed above
                mass_detrainment_updraft_forced[0, 0, 0][plume] = (
                    -normalized_massflux_updraft[0, 0, 1][plume]
                    + normalized_massflux_updraft[0, 0, 0][plume]
                    + mass_entrainment_updraft_forced[0, 0, 0][plume]
                )

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= max_index and K <= cloud_top_level[0, 0][plume]:
                # above location of maximum value normalized_massflux_updraft -> change entrainment
                height_massflux_avg = (
                    geopotential_height[0, 0, 1] - geopotential_height
                ) * normalized_massflux_updraft[0, 0, 0][plume]
                mass_entrainment_updraft_forced[0, 0, 0][plume] = (
                    entrainment_rate[0, 0, 0][plume] * height_massflux_avg
                )
                mass_detrainment_updraft_forced[0, 0, 0][plume] = (
                    normalized_massflux_updraft[0, 0, 0][plume]
                    + mass_entrainment_updraft_forced[0, 0, 0][plume]
                    - normalized_massflux_updraft[0, 0, 1][plume]
                )
                mass_detrainment_updraft_forced[0, 0, 0][plume] = max(
                    mass_detrainment_updraft_forced[0, 0, 0][plume], 0.0
                )

                # check mass_entrainment_updraft_forced in case it has been changed above
                mass_entrainment_updraft_forced[0, 0, 0][plume] = (
                    -normalized_massflux_updraft[0, 0, 0][plume]
                    + mass_detrainment_updraft_forced[0, 0, 0][plume]
                    + normalized_massflux_updraft[0, 0, 1][plume]
                )

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            mass_entrainment_updraft = mass_entrainment_updraft_forced[0, 0, 0][plume]
            mass_detrainment_updraft = mass_detrainment_updraft_forced[0, 0, 0][plume]

    with computation(FORWARD), interval(1, None):
        if error_code[0, 0][plume] == 0:
            # for weaker mixing
            mass_entrainment_u_updraft[0, 0, -1] = (
                mass_entrainment_updraft_forced[0, 0, -1][plume]
                + LAMBDA_DEEP * mass_detrainment_updraft_forced[0, 0, -1][plume]
            )
            mass_detrainment_u_updraft[0, 0, -1] = (
                mass_detrainment_updraft_forced[0, 0, -1][plume]
                + LAMBDA_DEEP * mass_detrainment_updraft_forced[0, 0, -1][plume]
            )


def compute_uc_vc(
    u_c: FloatField,
    v_c: FloatField,
    cloud_moist_static_energy: FloatField,
    cloud_moist_static_energy_forced: FloatField,
    error_code: IntFieldIJ_Plume,
    start_level: IntFieldIJ,
    moist_static_energy_origin_level: FloatFieldIJ,
    moist_static_energy_origin_level_forced: FloatFieldIJ,
    u_cloud_levels: FloatField,
    v_cloud_levels: FloatField,
    p: FloatField,
    updraft_origin_level: IntFieldIJ_Plume,
    ocean_fraction: FloatFieldIJ,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    """Compute u and v for the c-grid, update cloud moist static energy
    below start_level (nominally LCL level).

    Args:
        u_c (FloatField)
        v_c (FloatField)
        cloud_moist_static_energy (FloatField)
        cloud_moist_static_energy_forced (FloatField)
        error_code (IntFieldIJ_Plume)
        start_level (IntFieldIJ)
        moist_static_energy_origin_level (FloatFieldIJ)
        moist_static_energy_origin_level_forced (FloatFieldIJ)
        u_cloud_levels (FloatField)
        v_cloud_levels (FloatField)
        p (FloatField)
        updraft_origin_level (IntFieldIJ_Plume)
        ocean_fraction (FloatFieldIJ)
        AVERAGE_LAYER_DEPTH (Float)
        plume (Int)
    """
    from __externals__ import BOUNDARY_CONDITION_METHOD, k_end

    with computation(PARALLEL), interval(...):
        u_c = 0.0
        v_c = 0.0
        cloud_moist_static_energy = 0.0
        cloud_moist_static_energy_forced = 0.0

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break
        # this is never touched so long as compute_perturbation=False
        dummy_field_no_read = 0.0

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= start_level:
            cloud_moist_static_energy = moist_static_energy_origin_level
            cloud_moist_static_energy_forced = moist_static_energy_origin_level_forced

            # get uc and vc as average between layers below k22
            u_c = get_cloud_boundary_conditions(
                field=u_cloud_levels,
                scalar_perturbation=0,
                p=p,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )

            v_c = get_cloud_boundary_conditions(
                field=v_cloud_levels,
                scalar_perturbation=0,
                p=p,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )
