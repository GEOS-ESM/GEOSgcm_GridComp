from ndsl.dsl.gt4py import computation, interval, FORWARD, K, function, PARALLEL, exp
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, BoolFieldIJ
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
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
    updraft_origin_level: IntFieldIJ_Plume,
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
        updraft_origin_level[0, 0][plume] = 0

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
                    K=max(start_level, updraft_origin_level[0, 0][plume])
                ):
                    updraft_origin_level[0, 0][plume] = level
                level += 1

            updraft_origin_level[0, 0][plume] = updraft_origin_level[0, 0][plume] + start_level - 1
            updraft_origin_level[0, 0][plume] = max(updraft_origin_level[0, 0][plume], start_level)

            if plume == 0:
                # for shallow plumes
                updraft_origin_level[0, 0][plume] = min(2, updraft_origin_level[0, 0][plume])
                if updraft_origin_level[0, 0][plume] > maximum_updraft_origin_level:
                    error_code[0, 0][plume] = 2
            else:
                # other plumes
                if updraft_origin_level[0, 0][plume] > maximum_updraft_origin_level:
                    updraft_origin_level[0, 0][plume] = start_level


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
    updraft_origin_level: IntFieldIJ_Plume,
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
            lcl_level[0, 0][plume] = updraft_origin_level[0,0][plume]

            # get conditions for source parcel
            vapor_source = get_updraft_origin_conditions(
                field=vapor_cloud_levels_forced,
                scalar_perturbation=vapor_excess,
                p=p,
                updraft_origin_level=updraft_origin_level[0,0][plume],
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
                updraft_origin_level=updraft_origin_level[0,0][plume],
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
                updraft_origin_level=updraft_origin_level[0,0][plume],
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
                updraft_origin_level=updraft_origin_level[0,0][plume],
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
                    lcl_level[0, 0][plume] = max(K, updraft_origin_level[0,0][plume])
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


def set_start_level(updraft_origin_level: IntFieldIJ_Plume, start_level: IntFieldIJ, plume: Int):
    with computation(FORWARD), interval(...):
        updraft_origin_level[0, 0][plume] = start_level


def convective_cloud_base_level(
    error_code: IntFieldIJ_Plume,
    cloud_moist_static_energy_forced_transported: FloatField,
    cap_max: FloatFieldIJ,
    updraft_origin_level: IntFieldIJ_Plume,
    start_level: IntFieldIJ,
    moist_static_energy_origin_level_forced: FloatFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    maximum_updraft_origin_level: IntFieldIJ,
    negative_buoyancy_depth: FloatFieldIJ,
    frh_lfc: FloatFieldIJ,
    geopotential_height_cloud_levels_forced: FloatField,
    entrainment_rate: FloatField_Plume,
    environment_moist_static_energy_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    p_cloud_levels_forced: FloatField_Plume,
    vapor_forced: FloatField,
    environment_saturation_mixing_ratio_forced: FloatField,
    ocean_fraction: FloatFieldIJ,
    cap_max_increment: FloatFieldIJ,
    t_perturbation: FloatField,
    p_forced: FloatField,
    cloud_top: IntFieldIJ_Plume,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    """
    Determine level of free convection (LFC) for the source parcel.
    This stencil contains an open-ended vertical solver with multiple nested K-intervals.
    To implement this properly, the entire stencil has been constructed on an interval(0, 1),
    and all K read/writes have been done with absolute indexes or relative offsets.

    May also modify cloud_top, depending on LFC location.

    Args:
        error_code
        cloud_moist_static_energy_forced_transported
        cap_max
        updraft_origin_level
        start_level
        moist_static_energy_origin_level_forced
        updraft_lfc_level
        maximum_updraft_origin_level
        negative_buoyancy_depth
        frh_lfc
        geopotential_height_cloud_levels_forced
        entrainment_rate
        environment_moist_static_energy_forced
        environment_saturation_moist_static_energy_cloud_levels_forced
        t_excess
        vapor_excess
        add_buoyancy
        p_cloud_levels_forced
        vapor_forced
        environment_saturation_mixing_ratio_forced
        ocean_fraction
        cap_max_increment
        t_perturbation
        p_forced
        cloud_top
        AVERAGE_LAYER_DEPTH
        plume
    """
    from __externals__ import (
        OVERSHOOT,
        ZERO_DIFF,
        k_end,
        MOIST_TRIGGER,
        USE_MEMORY,
        BOUNDARY_CONDITION_METHOD,
    )

    with computation(PARALLEL), interval(...):
        # prefill some fields
        cloud_moist_static_energy_forced_transported = 0.0
        dby = 0.0

    with computation(FORWARD), interval(0, 1):
        # internal constants
        frh_crit_O = 0.7
        frh_crit_L = 0.7

        # mask for solver
        stop_solver = False

        # prefill some fields
        start_level = 0
        cap_max_internal = cap_max
        # default value
        updraft_lfc_level[0, 0][plume] = maximum_updraft_origin_level + 3
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
                cloud_moist_static_energy_forced_transported = moist_static_energy_origin_level_forced

    # determine the level of convective cloud base
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            while error_code[0, 0][plume] == 0 and stop_solver == False:
                updraft_lfc_level[0, 0][plume] = start_level_internal
                level = 0
                while level <= k_end:
                    dz = geopotential_height_cloud_levels_forced.at(
                        K=level
                    ) - geopotential_height_cloud_levels_forced.at(K=level - 1)
                    temporary = (
                        (1.0 - 0.5 * entrainment_rate.at(K=level - 1, ddim=[plume]) * dz)
                        * cloud_moist_static_energy_forced_transported.at(K=level - 1)
                        + entrainment_rate.at(K=level - 1, ddim=[plume])
                        * dz
                        * environment_moist_static_energy_forced.at(K=level - 1)
                    ) / (1.0 + 0.5 * entrainment_rate.at(K=level - 1, ddim=[plume]) * dz)
                    cloud_moist_static_energy_forced_transported[0, 0, level] = temporary
                    if K == start_level_internal + 1:
                        modification = (
                            cumulus_parameterization_constants.XLV * vapor_excess
                            + cumulus_parameterization_constants.CP * t_excess
                        ) + add_buoyancy
                        temporary = cloud_moist_static_energy_forced_transported.at(K=level) + modification
                        cloud_moist_static_energy_forced_transported[0, 0, level] = temporary
                    level += 1

                while (
                    cloud_moist_static_energy_forced_transported.at(K=updraft_lfc_level[0, 0][plume])
                    < environment_saturation_moist_static_energy_cloud_levels_forced.at(
                        K=updraft_lfc_level[0, 0][plume]
                    )
                    and stop_solver == False
                ):
                    updraft_lfc_level[0, 0][plume] += 1
                    if updraft_lfc_level[0, 0][plume] > maximum_updraft_origin_level + 2:
                        error_code[0, 0][plume] = 3
                        stop_solver = True

                if stop_solver == False:
                    # cloud base pressure and max moist static energy pressure
                    # i.e., the depth (in mb) of the layer of negative buoyancy
                    negative_buoyancy_depth = -(
                        p_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume], ddim=[plume])
                        - p_cloud_levels_forced.at(K=start_level_internal, ddim=[plume])
                    )

                    if MOIST_TRIGGER == 1:
                        frh_lfc = 0.0
                        dzh = 0.0
                        level = updraft_origin_level
                        while level <= updraft_lfc_level[0, 0][plume]:
                            dz = geopotential_height_cloud_levels_forced.at(
                                K=level
                            ) - geopotential_height_cloud_levels_forced.at(K=max(level - 1, 0))
                            frh_lfc = frh_lfc + dz * (
                                vapor_forced.at(K=level)
                                / environment_saturation_mixing_ratio_forced.at(K=level)
                            )
                            dzh = dzh + dz
                            level += 1

                        frh_lfc = frh_lfc / (dzh + 1.0e-16)
                        frh_crit = frh_crit_O * ocean_fraction + frh_crit_L * (1.0 - ocean_fraction)
                        fx = (2.0 / 0.78) * exp(-((frh_lfc - frh_crit) ** 2)) * (frh_lfc - frh_crit)
                        fx = max(-1.0, min(1.0, fx))

                        del_cap_max = fx * cap_max_increment
                        cap_max_internal = min(max(cap_max + del_cap_max, 10.0), 150.0)

                    # test if the air parcel has enough energy to reach the positive buoyant region
                    if cap_max_internal > negative_buoyancy_depth:
                        stop_solver = True

                if stop_solver == False:
                    # if am here -> kbcon not found for air parcels from k22 level
                    updraft_origin_level += 1
                    if USE_MEMORY == 20:
                        # increase capmax
                        cap_max_internal = cap_max_internal + cap_max_increment

                    # get new moist_static_energy_origin_level_forced
                    modification = (
                        cumulus_parameterization_constants.XLV * vapor_excess
                        + cumulus_parameterization_constants.CP * t_excess
                    ) + add_buoyancy
                    moist_static_energy_origin_level_forced = get_updraft_origin_conditions(
                        field=environment_moist_static_energy_forced,
                        scalar_perturbation=modification,
                        p=p_forced,
                        updraft_origin_level=updraft_origin_level,
                        ocean_fraction=ocean_fraction,
                        BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                        AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                        k_end=k_end,
                        compute_perturbation=True,
                        perturbation_field=t_perturbation,
                    )
                    start_level_internal += 1
                    cloud_moist_static_energy_forced_transported[0, 0, start_level_internal] = (
                        moist_static_energy_origin_level_forced
                    )

            if updraft_lfc_level[0, 0][plume] == 0:
                error_code[0, 0][plume] = 33

    with computation(FORWARD), interval(0, 1):
        cloud_top[0, 0][plume] = k_end - 1
        if error_code[0, 0][plume] == 0:
            start_level_internal = updraft_lfc_level[0, 0][plume]

    with computation(FORWARD), interval(1, -1):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1:
                dz = (
                    geopotential_height_cloud_levels_forced
                    - geopotential_height_cloud_levels_forced[0, 0, -1]
                )
                denom = 1.0 + 0.5 * entrainment_rate[0, 0, -1][plume] * dz
                if denom == 0.0:
                    cloud_moist_static_energy_forced_transported = (
                        cloud_moist_static_energy_forced_transported[0, 0, -1]
                    )
                else:
                    cloud_moist_static_energy_forced_transported = (
                        (1.0 - 0.5 * entrainment_rate[0, 0, -1][plume] * dz)
                        * cloud_moist_static_energy_forced_transported[0, 0, -1]
                        + entrainment_rate[0, 0, -1][plume]
                        * dz
                        * environment_moist_static_energy_forced[0, 0, -1]
                    ) / denom

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # mask for next computation block
            top_found = False

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K >= start_level and K <= k_end - 1 and top_found == False:
                if (
                    cloud_moist_static_energy_forced_transported
                    < environment_saturation_moist_static_energy_cloud_levels_forced
                ):
                    cloud_top[0, 0][plume] = K - 1
                    top_found = True

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if cloud_top[0, 0][plume] <= updraft_lfc_level[0, 0][plume] + 1:
                error_code[0, 0][plume] = 41

            if OVERSHOOT > 1.0e-6 and error_code[0, 0][plume] == 0:
                z_overshoot = (1.0 + OVERSHOOT) * geopotential_height_cloud_levels_forced.at(
                    K=cloud_top[0, 0][plume]
                )

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and OVERSHOOT > 1.0e-6:
            if (
                K >= cloud_top[0, 0][plume]
                and K <= k_end - 2
                and z_overshoot < geopotential_height_cloud_levels_forced
            ):
                cloud_top[0, 0][plume] = min(K - 1, k_end - 2)


def updraft_rates_pdf(
    entrainment_rate: FloatField_Plume,
    moist_static_energy: FloatField,
    saturation_moist_static_energy: FloatField,
    moist_static_energy_origin_level: FloatFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    geopotential_height: FloatField,
    cloud_moist_static_energy: FloatField,
    error_code: IntFieldIJ_Plume,
    cloud_top: IntFieldIJ_Plume,
    plume: Int,
):
    from __externals__ import k_end, OVERSHOOT

    with computation(PARALLEL), interval(...):
        # default value
        cloud_moist_static_energy = 0.0

    with computation(FORWARD), interval(0, 1):
        if plume != 0:
            # default value
            cloud_top[0, 0][plume] = k_end - 3

        if plume != 0 and error_code[0, 0][plume] == 0:
            start_level: IntFieldIJ = updraft_lfc_level[0, 0][plume]

    with computation(PARALLEL), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if K <= start_level:
                cloud_moist_static_energy = moist_static_energy_origin_level

    with computation(FORWARD), interval(1, None):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if K > start_level and K < k_end - 2:
                dz = geopotential_height - geopotential_height[0, 0, -1]

                cloud_moist_static_energy = (
                    (1.0 - 0.5 * entrainment_rate[0, 0, -1][plume] * dz) * cloud_moist_static_energy[0, 0, -1]
                    + entrainment_rate[0, 0, -1][plume] * dz * moist_static_energy[0, 0, -1]
                ) / (1.0 + 0.5 * entrainment_rate[0, 0, -1][plume] * dz)

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            # set up mask for next computation
            stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if K > start_level and K < k_end - 2 and stop_computation == False:
                # find the height where the parcel is no longer saturated
                if cloud_moist_static_energy < saturation_moist_static_energy:
                    cloud_top[0, 0][plume] = K - 1
                    stop_computation = True

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if cloud_top[0, 0][plume] <= updraft_lfc_level[0, 0][plume] + 1:
                error_code[0, 0][plume] = 41

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0 and OVERSHOOT > 0:
            z_overshoot = (1.0 + OVERSHOOT) * geopotential_height.at(K=cloud_top[0, 0][plume])

    with computation(FORWARD), interval(0, 1):
        # set up mask for next computation
        stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and plume != 0 and OVERSHOOT > 0:
            if K >= cloud_top[0, 0][plume] and K < k_end - 1 and stop_computation == False:
                if z_overshoot < geopotential_height:
                    cloud_top[0, 0][plume] = min(K - 1, k_end - 3)
                    stop_computation = True


def cloud_top_checks(
    cloud_top: IntFieldIJ_Plume,
    p: FloatField_Plume,
    geopotential_height: FloatField,
    error_code: IntFieldIJ_Plume,
    last_error_code: IntFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    MINIMUM_DEPTH: Float,
    plume: Int,
):
    # check if cloud_top is too low for deep convection
    with computation(FORWARD), interval(0, 1):
        if plume == 2 and error_code[0, 0][plume] == 0:
            if p.at(K=cloud_top[0, 0][plume], ddim=[plume]) > 400:
                error_code[0, 0][plume] = 22

    # check if cloud_top is too high for mid convection
    with computation(FORWARD), interval(0, 1):
        if plume == 1 and error_code[0, 0][plume] == 0:
            if p.at(K=cloud_top[0, 0][plume], ddim=[plume]) < 400:
                error_code[0, 0][plume] = 22

    # check if cloud_top is too high for shallow convection
    with computation(FORWARD), interval(0, 1):
        if plume == 0 and error_code[0, 0][plume] == 0:
            if p.at(K=cloud_top[0, 0][plume], ddim=[plume]) < 400:
                error_code[0, 0][plume] = 23

    # avoid double-counting plumes
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and last_error_code == 0:
            if p.at(K=cloud_top[0, 0][plume], ddim=[plume]) < 700:
                error_code[0, 0][plume] = 27

    # last checks for cloud_top
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if (
                geopotential_height.at(K=cloud_top[0, 0][plume])
                - geopotential_height.at(K=updraft_lfc_level[0, 0][plume])
                < MINIMUM_DEPTH
            ):
                error_code[0, 0][plume] = 5
        if cloud_top[0, 0][plume] <= updraft_lfc_level[0, 0][plume]:
            error_code[0, 0][plume] = 5
