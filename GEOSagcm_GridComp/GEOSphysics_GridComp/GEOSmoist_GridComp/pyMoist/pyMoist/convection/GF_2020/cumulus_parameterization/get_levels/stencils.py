from ndsl.dsl.gt4py import computation, interval, FORWARD, K, function, PARALLEL
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, BoolFieldIJ
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    get_updraft_origin_conditions,
    compute_dewpoint,
    column_max,
)
import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants


def find_maximum_updraft_origin_level(
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    maximum_updraft_origin_level: IntFieldIJ,
    MAX_UPDRAFT_ORIGIN_HEIGHT: Float,
    plume: Int,
):
    """
    Set the highest level from which an updraft may originate.

    Args:
        geopotential_height (in)
        topography_height_no_negative (in): topography height, without features extending below sea level
        error_code (in)
        maximum_updraft_origin_level (out)
        MAX_UPDRAFT_ORIGIN_HEIGHT (in): max height(m) above ground where updraft air can originate
        plume (in): specifies the current plume
    """
    with computation(FORWARD), interval(0, 1):
        stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0 and stop_computation == 0:
            if geopotential_height > MAX_UPDRAFT_ORIGIN_HEIGHT + topography_height_no_negative:
                maximum_updraft_origin_level = K
                stop_computation: BoolFieldIJ = True


def find_detrainmet_start_level(
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    detrainment_start_level: IntFieldIJ,
    DETRAINMENT_CRITICAL_DEPTH: Float,
    plume: Int,
):
    """
    Set the highest level from which an updraft may originate.

    Args:
        geopotential_height (in)
        topography_height_no_negative (in): topography height, without features extending below sea level
        error_code (in)
        maximum_updraft_origin_level (out)
        DETRAINMENT_CRITICAL_DEPTH (in): depth(m) over which downdraft detrains all its mass
        plume (in): specifies the current plume
    """
    with computation(FORWARD), interval(0, 1):
        stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0 and stop_computation == 0:
            if geopotential_height > DETRAINMENT_CRITICAL_DEPTH + topography_height_no_negative:
                detrainment_start_level = K
                stop_computation: BoolFieldIJ = True


def find_highest_moist_static_energy_level(
    moist_static_energy: FloatField,
    error_code: IntFieldIJ_Plume,
    maximum_updraft_origin_level: IntFieldIJ,
    updraft_origin_level: IntFieldIJ,
    plume: Int,
):
    """
    Find the level with the highest moist static energy - which will
    be used as the starting point for any subsequent updrafts

    Args:
        moist_static_energy (in)
        error_code (in)
        maximum_updraft_origin_level (in)
        updraft_origin_level (out)
        plume (in): specifies the current plume
    """
    with computation(FORWARD), interval(0, 1):
        # prefil output
        updraft_origin_level = 0

        if plume == 0:
            # start at surface for shallow plume
            start_level: IntFieldIJ = 0
        else:
            # start above surface
            start_level: IntFieldIJ = 1

        if error_code[0, 0][plume] == 0:
            level = start_level
            while level <= maximum_updraft_origin_level + 1:
                if moist_static_energy.at(K=level) > moist_static_energy.at(
                    K=max(start_level, updraft_origin_level)
                ):
                    updraft_origin_level = level
                level += 1

            updraft_origin_level = updraft_origin_level + start_level - 1
            updraft_origin_level = max(updraft_origin_level, start_level)

            if plume == 0:
                # for shallow plumes
                updraft_origin_level = min(2, updraft_origin_level)
                if updraft_origin_level > maximum_updraft_origin_level:
                    error_code[0, 0][plume] = 2
            else:
                # other plumes
                if updraft_origin_level > maximum_updraft_origin_level:
                    updraft_origin_level = start_level


@function
def get_lcl(p_source, t_source, vapor_source):
    # internal constants
    cpg = 102.45
    rgas = 287.0
    cp = 1004.0
    p00 = 1.0e5
    g = 9.80
    rocp = rgas / cp
    p00i = 1.0 / p00
    cpor = cp / rgas
    cpi = 1.0 / cp
    p00k = 26.870941  #  = p00 ** rocp
    p00ki = 1.0 / p00k

    # convert to pascals
    p_source = 100 * p_source

    # simpler, cheaper method
    dewpoint = compute_dewpoint(p_source, vapor_source)
    t_lcl = dewpoint - (0.001296 * dewpoint + 0.1963) * (t_source - dewpoint)
    p_lcl = p_source * (t_lcl / t_source) ** cpor
    dz_lcl = 127 * (t_source - dewpoint)
    if dz_lcl <= 0.0:
        dz_lcl = -999.0

    return t_lcl, p_lcl, dz_lcl


def find_lcl(
    p: FloatField,
    p_cloud_levels: FloatField,
    t_excess: FloatFieldIJ,
    t_cloud_levels_forced: FloatField,
    t_perturbation: FloatField,
    vapor_excess: FloatFieldIJ,
    vapor_cloud_levels_forced: FloatField,
    omega: FloatField,
    air_density: FloatField,
    geopotential_height_cloud_levels: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    updraft_origin_level: IntFieldIJ,
    grid_length: FloatFieldIJ,
    lcl_level: IntFieldIJ_Plume,
    error_code: IntFieldIJ_Plume,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    from __externals__ import k_end, BOUNDARY_CONDITION_METHOD, ADV_TRIGGER

    with computation(PARALLEL), interval(...):
        dummy_field_no_read = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # default value
            lcl_level[0, 0][plume] = updraft_origin_level

            # get conditions for source parcel
            vapor_source = get_updraft_origin_conditions(
                field=vapor_cloud_levels_forced,
                scalar_perturbation=vapor_excess,
                p=p,
                updraft_origin_level=updraft_origin_level,
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )
            t_source = get_updraft_origin_conditions(
                field=t_cloud_levels_forced,
                scalar_perturbation=t_excess,
                p=p,
                updraft_origin_level=updraft_origin_level,
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )
            p_source = get_updraft_origin_conditions(
                field=p_cloud_levels,
                scalar_perturbation=0,
                p=p,
                updraft_origin_level=updraft_origin_level,
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )

        # initalize 2d temporaries
        p_lcl: FloatFieldIJ = 0.0
        t_lcl: FloatFieldIJ = 0.0
        dz_lcl: FloatFieldIJ = 0.0
        z_source: FloatFieldIJ = 0.0
        stop_computation: BoolFieldIJ = False

        p_lcl, t_lcl, dz_lcl = get_lcl(p_source=p_source, t_source=t_source, vapor_source=vapor_source)

        if dz_lcl >= 0.0:
            z_source = get_updraft_origin_conditions(
                field=geopotential_height_cloud_levels,
                scalar_perturbation=0,
                p=p,
                updraft_origin_level=updraft_origin_level,
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if dz_lcl >= 0.0:
                if geopotential_height_cloud_levels >= z_source + dz_lcl and stop_computation == False:
                    lcl_level[0, 0][plume] = max(K, updraft_origin_level)
                    stop_computation = True

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            lcl_level[0, 0][plume] = min(lcl_level[0, 0][plume], k_end - 5)
            if cumulus_parameterization_constants.USE_LCL == True and plume == 1:
                error_code[0, 0][plume] = 21

            if ADV_TRIGGER == 1 and plume == 0:
                wkf: FloatFieldIJ = 0.02  # m/s

                level = lcl_level[0, 0][plume]
                dz_lcl = geopotential_height_cloud_levels - topography_height_no_negative

                ckf = wkf
                if dz_lcl <= 2.0e3:
                    ckf = wkf * dz_lcl / 2000

                    wkflcl = (
                        -(omega.at(K=max(0, level - 1))) / air_density.at(K=max(0, level - 1))
                        + omega / air_density
                        + omega[0, 0, 1] / air_density[0, 0, 1] / (3.0 * constants.MAPL_GRAV)
                    )

                    # check to see if cloud is buoyant using fritsch-chappell trigger
                    # function described in kain and fritsch (1992)...w0avg is an
                    # aproximate value for the running-mean grid-scale vertical
                    # velocity, which gives smoother fields of convective initiation
                    # than the instantaneous value...formula relating temperature
                    # perturbation to vertical velocity has been used with the most
                    # success at grid lengths near 25 km.  for different grid-lengths,
                    # adjust vertical velocity to equivalent value for 25 km grid
                    # length, assuming linear dependence of w on grid length...
                    if grid_length >= 25.0e3:
                        wkflcl = wkflcl * grid_length / 25.0e3 - ckf
                    else:
                        wkflcl = wkflcl - ckf

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and ADV_TRIGGER == 1 and plume == 0:
            # Kain (2004) Eq. 1
            t_perturbation = 4.64 * wkflcl ** (1.0 / 3.0)

            if t_perturbation > 2.0:
                t_perturbation = 2.0
            elif t_perturbation < -2.0:
                t_perturbation = -2.0


def set_start_level(
    updraft_origin_level: IntFieldIJ,
    start_level: IntFieldIJ,
):
    with computation(FORWARD), interval(...):
        updraft_origin_level = start_level


def convective_cloud_base_level(
    error_code: IntFieldIJ_Plume,
    cloud_moist_static_energy_forced_t: FloatField,
    cap_max: FloatFieldIJ,
    updraft_origin_level: IntFieldIJ,
    start_level: IntFieldIJ,
    moist_static_energy_origin_level_forced: FloatFieldIJ,
    updraft_lfc_level: IntFieldIJ,
    maximum_updraft_origin_level: IntFieldIJ,
    negative_buoyancy_depth: FloatFieldIJ,
    frh_lfc: FloatFieldIJ,
    geopotential_height_cloud_levels_forced: FloatField,
    entrainment_rate: FloatField,
    environment_moist_static_energy_forced: FloatField,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    plume: Int,
):
    from __externals__ import OVERSHOOT, ZERO_DIFF

    with computation(PARALLEL), interval(...):
        cloud_moist_static_energy_forced_t = 0.0
        dby = 0.0

    with computation(FORWARD), interval(0, 1):
        start_level = 0
        cap_max_internal = cap_max
        updraft_lfc_level = maximum_updraft_origin_level + 3
        negative_buoyancy_depth = 0.0
        frh_lfc = 0.0
        if error_code[0, 0][plume] == 0:
            if ZERO_DIFF == 1:
                start_level_internal = updraft_origin_level
            else:
                start_level_internal = start_level

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= start_level_internal:
                cloud_moist_static_energy_forced_t = moist_static_energy_origin_level_forced

    # determine the level of convective cloud base
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            while error_code[0, 0][plume] == 0:
                updraft_lfc_level = (
                    start_level_internal  # NOTE 2D FIELD, PROBLEM - WILL CAUSE NO ALL AXIS PRESENT ERROR
                )
                dz = (
                    geopotential_height_cloud_levels_forced
                    - geopotential_height_cloud_levels_forced[0, 0, -1]
                )
                cloud_moist_static_energy_forced_t = (
                    (1.0 - 0.5 * entrainment_rate[0, 0, -1] * dz)
                    * cloud_moist_static_energy_forced_t[0, 0, -1]
                    + entrainment_rate[0, 0, -1] * dz * environment_moist_static_energy_forced[0, 0, -1]
                ) / (1.0 + 0.5 * entrainment_rate[0, 0, -1] * dz)
                if K == start_level_internal + 1:
                    modification = (
                        cumulus_parameterization_constants.XLV * vapor_excess
                        + cumulus_parameterization_constants.CP * t_excess
                    ) + add_buoyancy
                    cloud_moist_static_energy_forced_t = cloud_moist_static_energy_forced_t + modification
