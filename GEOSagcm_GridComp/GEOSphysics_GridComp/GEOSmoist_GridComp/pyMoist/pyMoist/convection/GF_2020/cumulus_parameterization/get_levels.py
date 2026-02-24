import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl import QuantityFactory, StencilFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.gt4py import FORWARD, PARALLEL, K, computation, exp, function, interval
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ, Int, IntFieldIJ
from ndsl.logging import ndsl_log
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
from pyMoist.convection.GF_2020.cumulus_parameterization.locals import GF2020CumulusParameterizationLocals
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    compute_dewpoint,
    get_cloud_boundary_conditions,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState


def find_maximum_updraft_origin_level(
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    maximum_updraft_origin_level: IntFieldIJ,
    MAX_UPDRAFT_ORIGIN_HEIGHT: Float,
    plume: Int,
):
    """Determine the highest level from which an updraft may originate.

    Args:
        geopotential_height (FloatField)
        topography_height_no_negative (FloatFieldIJ)
        error_code (IntFieldIJ_Plume)
        maximum_updraft_origin_level (IntFieldIJ)
        MAX_UPDRAFT_ORIGIN_HEIGHT (Float)
        plume (Int)
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
    """Determine the lowest level at which detrainment may occur.

    Args:
        geopotential_height (FloatField)
        topography_height_no_negative (FloatFieldIJ)
        error_code (IntFieldIJ_Plume)
        detrainment_start_level (IntFieldIJ)
        DETRAINMENT_CRITICAL_DEPTH (Float)
        plume (Int)
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
    """Detemine the level with the highest moist static energy. This level will
        be used as the starting point for any subsequent updrafts

    Args:
        moist_static_energy (FloatField)
        error_code (IntFieldIJ_Plume)
        maximum_updraft_origin_level (IntFieldIJ)
        updraft_origin_level (IntFieldIJ_Plume)
        plume (Int)
    """
    with computation(FORWARD), interval(0, 1):
        # prefil output
        updraft_origin_level[0, 0][plume] = 0

        if plume == cumulus_parameterization_constants.SHALLOW:
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

            if plume == cumulus_parameterization_constants.SHALLOW:
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
    p00k = 26.870941
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
    """Determine the LCL level. For shallow plumes, if LCL is buoyant, calculate a temperature perturbation.

    Args:
        p (FloatField)
        p_cloud_levels (FloatField)
        t_excess (FloatFieldIJ)
        t_cloud_levels_forced (FloatField)
        t_perturbation (FloatField)
        vapor_excess (FloatFieldIJ)
        vapor_cloud_levels_forced (FloatField)
        omega (FloatField)
        air_density (FloatField)
        geopotential_height_cloud_levels (FloatField)
        topography_height_no_negative (FloatFieldIJ)
        ocean_fraction (FloatFieldIJ)
        updraft_origin_level (IntFieldIJ_Plume)
        grid_length (FloatFieldIJ)
        lcl_level (IntFieldIJ_Plume)
        error_code (IntFieldIJ_Plume)
        AVERAGE_LAYER_DEPTH (Float)
        plume (Int)
    """
    from __externals__ import ADV_TRIGGER, BOUNDARY_CONDITION_METHOD, k_end

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break
        # this is never touched so long as compute_perturbation=False
        dummy_field_no_read = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # default value
            lcl_level[0, 0][plume] = updraft_origin_level[0, 0][plume]

            # get conditions for source parcel
            vapor_source = get_cloud_boundary_conditions(
                field=vapor_cloud_levels_forced,
                scalar_perturbation=vapor_excess,
                p=p,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )
            t_source = get_cloud_boundary_conditions(
                field=t_cloud_levels_forced,
                scalar_perturbation=t_excess,
                p=p,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )
            p_source = get_cloud_boundary_conditions(
                field=p_cloud_levels,
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

        # initialize 2d temporaries
        p_lcl: FloatFieldIJ = 0.0
        t_lcl: FloatFieldIJ = 0.0
        dz_lcl: FloatFieldIJ = 0.0
        z_source: FloatFieldIJ = 0.0
        stop_computation: BoolFieldIJ = False

        p_lcl, t_lcl, dz_lcl = get_lcl(p_source=p_source, t_source=t_source, vapor_source=vapor_source)

        if dz_lcl >= 0.0:
            z_source = get_cloud_boundary_conditions(
                field=geopotential_height_cloud_levels,
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

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if dz_lcl >= 0.0:
                if geopotential_height_cloud_levels >= z_source + dz_lcl and stop_computation == False:
                    lcl_level[0, 0][plume] = max(K, updraft_origin_level[0, 0][plume])
                    stop_computation = True

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            lcl_level[0, 0][plume] = min(lcl_level[0, 0][plume], k_end - 5)
            if (
                cumulus_parameterization_constants.USE_LCL == True
                and plume == cumulus_parameterization_constants.MID
            ):
                error_code[0, 0][plume] = 21

            if ADV_TRIGGER == 1 and plume == cumulus_parameterization_constants.SHALLOW:
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
        if (
            error_code[0, 0][plume] == 0
            and ADV_TRIGGER == 1
            and plume == cumulus_parameterization_constants.SHALLOW
        ):
            # Kain (2004) Eq. 1
            t_perturbation = 4.64 * wkflcl ** (1.0 / 3.0)

            if t_perturbation > 2.0:
                t_perturbation = 2.0
            elif t_perturbation < -2.0:
                t_perturbation = -2.0


def set_start_level(lcl_level: IntFieldIJ_Plume, start_level: IntFieldIJ, plume: Int):
    """Copy LCL level data into another field called "start_level", which will modified later when searching
    for an equilibrium level.

    Args:
        lcl_level (IntFieldIJ_Plume)
        start_level (IntFieldIJ)
        plume (Int)
    """
    with computation(FORWARD), interval(...):
        start_level = lcl_level[0, 0][plume]


def get_convective_cloud_base_level(
    error_code: IntFieldIJ_Plume,
    lcl_level: IntFieldIJ_Plume,
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
    environment_moist_static_energy_cloud_levels_forced: FloatField,
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
    cloud_top_level: IntFieldIJ_Plume,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    """Determine level of free convection (LFC) for the source parcel.

    May also modify cloud_top_level, depending on LFC location.

    This stencil contains an open-ended vertical solver with multiple nested K-intervals.
    To implement this properly, the entire stencil has been constructed on an interval(0, 1),
    and all K read/writes have been done with absolute indexes or relative offsets. The alternative is
    to break this into a series of stencils (at least seven, based on current understanding of the structure
    of this code) and pass data between them using a much larger number of locals.

    The current method (single stencil) is horrible for optimization and speed, but was chosen nonetheless
    because it allows for more readable code and should allow for a more seamless implementation of a proper
    solution once the necessary tool/feature has been developed.

    Args:
        error_code (IntFieldIJ_Plume): _description_
        lcl_level (IntFieldIJ_Plume): _description_
        cloud_moist_static_energy_forced_transported (FloatField): _description_
        cap_max (FloatFieldIJ): _description_
        updraft_origin_level (IntFieldIJ_Plume): _description_
        start_level (IntFieldIJ): _description_
        moist_static_energy_origin_level_forced (FloatFieldIJ): _description_
        updraft_lfc_level (IntFieldIJ_Plume): _description_
        maximum_updraft_origin_level (IntFieldIJ): _description_
        negative_buoyancy_depth (FloatFieldIJ): _description_
        frh_lfc (FloatFieldIJ): _description_
        geopotential_height_cloud_levels_forced (FloatField): _description_
        entrainment_rate (FloatField_Plume): _description_
        environment_moist_static_energy_forced (FloatField): _description_
        environment_moist_static_energy_cloud_levels_forced (FloatField): _description_
        environment_saturation_moist_static_energy_cloud_levels_forced (FloatField): _description_
        t_excess (FloatFieldIJ): _description_
        vapor_excess (FloatFieldIJ): _description_
        add_buoyancy (FloatFieldIJ): _description_
        p_cloud_levels_forced (FloatField_Plume): _description_
        vapor_forced (FloatField): _description_
        environment_saturation_mixing_ratio_forced (FloatField): _description_
        ocean_fraction (FloatFieldIJ): _description_
        cap_max_increment (FloatFieldIJ): _description_
        t_perturbation (FloatField): _description_
        p_forced (FloatField): _description_
        cloud_top_level (IntFieldIJ_Plume): _description_
        AVERAGE_LAYER_DEPTH (Float): _description_
        plume (Int): _description_
    """
    from __externals__ import (
        BOUNDARY_CONDITION_METHOD,
        MOIST_TRIGGER,
        OVERSHOOT,
        USE_MEMORY,
        ZERO_DIFF,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        # prefill some fields
        cloud_moist_static_energy_forced_transported = 0.0
        dby = 0.0

        # make garbage field so the get_cloud_boundary_conditions call does not break
        # this is never touched so long as compute_perturbation=False
        dummy_field_no_read = 0.0

    with computation(FORWARD), interval(0, 1):
        # internal constants
        frh_crit_O = 0.7
        frh_crit_L = 0.7

        # prefill some fields
        start_level_internal: IntFieldIJ = 0
        cap_max_internal = cap_max

        # initialize some 2d temporaries
        dzh: FloatFieldIJ = 0.0
        found_level: BoolFieldIJ = False

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and ZERO_DIFF == 1:
            start_level_internal = lcl_level[0, 0][plume]
        else:
            start_level_internal = start_level

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= start_level_internal:
            cloud_moist_static_energy_forced_transported = (
                moist_static_energy_origin_level_forced  # assumed no entraiment between these layers
            )

    # determine the level of convective cloud base (updraft_lfc_level)

    with computation(FORWARD), interval(0, 1):
        skip_last_check = False
        updraft_lfc_level[0, 0][plume] = maximum_updraft_origin_level + 3
        negative_buoyancy_depth = 0.0
        frh_lfc = 0.0
        if error_code[0, 0][plume] == 0:
            continue_outer_while_loop = True
            while error_code[0, 0][plume] == 0 and continue_outer_while_loop == True:
                updraft_lfc_level[0, 0][plume] = start_level_internal
                level = start_level_internal + 1
                while level <= maximum_updraft_origin_level + 3:
                    dz = (
                        geopotential_height_cloud_levels_forced[0, 0, level]
                        - geopotential_height_cloud_levels_forced[0, 0, level - 1]
                    )
                    cloud_moist_static_energy_forced_transported[0, 0, level] = (
                        (1.0 - 0.5 * entrainment_rate[0, 0, level - 1][plume] * dz)
                        * cloud_moist_static_energy_forced_transported[0, 0, level - 1]
                        + entrainment_rate[0, 0, level - 1][plume]
                        * dz
                        * environment_moist_static_energy_forced[0, 0, level - 1]
                    ) / (1.0 + 0.5 * entrainment_rate[0, 0, level - 1][plume] * dz)
                    if level == start_level_internal + 1:
                        modification = (
                            cumulus_parameterization_constants.XLV * vapor_excess
                            + cumulus_parameterization_constants.CP * t_excess
                        ) + add_buoyancy
                        cloud_moist_static_energy_forced_transported[0, 0, level] = (
                            cloud_moist_static_energy_forced_transported[0, 0, level] + modification
                        )
                    level += 1

                continue_inner_while_loop = True
                while (
                    cloud_moist_static_energy_forced_transported.at(K=updraft_lfc_level[0, 0][plume])
                    < environment_saturation_moist_static_energy_cloud_levels_forced.at(
                        K=updraft_lfc_level[0, 0][plume]
                    )
                ) and continue_inner_while_loop == True:
                    updraft_lfc_level[0, 0][plume] = updraft_lfc_level[0, 0][plume] + 1
                    if updraft_lfc_level[0, 0][plume] > maximum_updraft_origin_level + 2:
                        error_code[0, 0][plume] = 3
                        continue_inner_while_loop = False

                if error_code[0, 0][plume] != 0:
                    continue_outer_while_loop = False
                    skip_last_check = True

                if continue_outer_while_loop == True:
                    # cloud base pressure and max moist static energy pressure
                    # i.e., the depth (in mb) of the layer of negative buoyancy
                    negative_buoyancy_depth = -(
                        p_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume], ddim=[plume])
                        - p_cloud_levels_forced.at(K=start_level_internal, ddim=[plume])
                    )

                    if MOIST_TRIGGER == 1:
                        frh_lfc = 0.0
                        dzh = 0
                        level = updraft_origin_level[0, 0][plume]
                        while level <= updraft_lfc_level[0, 0][plume]:
                            dz = (
                                geopotential_height_cloud_levels_forced[0, 0, level]
                                - geopotential_height_cloud_levels_forced[0, 0, max(level - 1, 0)]
                            )
                            frh_lfc = frh_lfc + dz * (
                                vapor_forced[0, 0, level]
                                / environment_moist_static_energy_forced[0, 0, level]
                            )
                            dzh = dzh + dz

                        frh_lfc = frh_lfc / (dzh + 1.0e-16)
                        frh_crit = frh_crit_O * ocean_fraction + frh_crit_L * (1.0 - ocean_fraction)

                        fx = (
                            (2.0 / 0.78) * exp(-((frh_lfc - frh_crit) ** 2)) * (frh_lfc - frh_crit)
                        )  # exponential
                        fx = max(-1.0, min(1.0, fx))

                        del_cap_max = fx * cap_max_increment

                        cap_max_internal = min(max(cap_max + del_cap_max, 10.0), 150.0)

                    # test if the air parcel has enough energy to reach the positive buoyant region
                    if cap_max_internal > negative_buoyancy_depth:
                        continue_outer_while_loop = False
                        skip_last_check = True

                if continue_outer_while_loop == True:
                    # if am here -> updraft_lfc_level not found for air parcels from updraft_origin_level
                    updraft_origin_level[0, 0][plume] = updraft_origin_level[0, 0][plume] + 1

                    if USE_MEMORY == 20:
                        cap_max_internal = cap_max_internal + cap_max_increment

                    # get new moist_static_energy_origin_level_forced
                    modification = (
                        cumulus_parameterization_constants.XLV * vapor_excess
                        + cumulus_parameterization_constants.CP * t_excess
                    ) + add_buoyancy
                    moist_static_energy_origin_level_forced = get_cloud_boundary_conditions(
                        field=environment_moist_static_energy_cloud_levels_forced,
                        scalar_perturbation=modification,
                        p=p_forced,
                        updraft_origin_level=updraft_origin_level[0, 0][plume],
                        ocean_fraction=ocean_fraction,
                        BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                        AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                        k_end=k_end,
                        compute_perturbation=True,
                        perturbation_field=t_perturbation,
                    )

                    start_level_internal = start_level_internal + 1
                    cloud_moist_static_energy_forced_transported[
                        0, 0, start_level_internal
                    ] = moist_static_energy_origin_level_forced

            if skip_last_check == True:
                # last check for updraft_lfc_level
                if updraft_lfc_level[0, 0][plume] == 0:
                    error_code[0, 0][plume] = 33

    # determine the level of neutral buoyancy - ktop

    with computation(FORWARD), interval(0, 1):
        cloud_top_level[0, 0][plume] = k_end - 2

        if error_code[0, 0][plume] == 0:
            start_level_internal = updraft_lfc_level[0, 0][plume]

    with computation(FORWARD), interval(1, -2):
        if error_code[0, 0][plume] == 0 and K > start_level_internal:
            dz = geopotential_height_cloud_levels_forced - geopotential_height_cloud_levels_forced[0, 0, -1]
            denom = 1.0 + 0.5 * entrainment_rate[0, 0, -1][plume] * dz
            if denom == 0.0:
                cloud_moist_static_energy_forced_transported = cloud_moist_static_energy_forced_transported[
                    0, 0, -1
                ]
            else:
                cloud_moist_static_energy_forced_transported = (
                    (1.0 - 0.5 * entrainment_rate[0, 0, -1][plume] * dz)
                    * cloud_moist_static_energy_forced_transported[0, 0, -1]
                    + entrainment_rate[0, 0, -1][plume]
                    * dz
                    * environment_moist_static_energy_forced[0, 0, -1]
                ) / denom

    with computation(FORWARD), interval(0, -2):
        if error_code[0, 0][plume] == 0 and K > start_level_internal and found_level == False:
            if (
                cloud_moist_static_energy_forced_transported
                < environment_saturation_moist_static_energy_cloud_levels_forced
            ):
                cloud_top_level[0, 0][plume] = K - 1
                found_level = True

    with computation(FORWARD), interval(0, 1):
        if (
            error_code[0, 0][plume] == 0
            and cloud_top_level[0, 0][plume] <= updraft_lfc_level[0, 0][plume] + 1
        ):
            error_code[0, 0][plume] = 41

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and OVERSHOOT >= 1.0e-6:
            z_overshoot: FloatFieldIJ = (1.0 + OVERSHOOT) * geopotential_height_cloud_levels_forced.at(
                K=cloud_top_level[0, 0][plume]
            )

    with computation(FORWARD), interval(0, 1):
        # reset mask
        found_level = False

    with computation(FORWARD), interval(0, -3):
        if (
            error_code[0, 0][plume] == 0
            and OVERSHOOT >= 1.0e-6
            and K >= cloud_top_level[0, 0][plume]
            and found_level == False
        ):
            if z_overshoot < geopotential_height_cloud_levels_forced:
                cloud_top_level[0, 0][plume] = min(K - 1, k_end - 2)
                found_level = True


def get_cloud_top(
    entrainment_rate: FloatField_Plume,
    environment_moist_static_energy_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    moist_static_energy_origin_level_forced: FloatFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_moist_static_energy_forced_transported: FloatField,
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    plume: Int,
):
    """Determine the initial estimate for the cloud top level. This can include an overshooting top above the
    equilibrium level, if the configuration option OVERSHOOT == True.

    Args:
        entrainment_rate (FloatField_Plume)
        environment_moist_static_energy_forced (FloatField)
        environment_saturation_moist_static_energy_cloud_levels_forced (FloatField)
        moist_static_energy_origin_level_forced (FloatFieldIJ)
        updraft_lfc_level (IntFieldIJ_Plume)
        geopotential_height_cloud_levels_forced (FloatField)
        cloud_moist_static_energy_forced_transported (FloatField)
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        plume (Int)
    """
    from __externals__ import OVERSHOOT, k_end

    with computation(PARALLEL), interval(...):
        # default value
        cloud_moist_static_energy_forced_transported = 0.0

    with computation(FORWARD), interval(0, 1):
        if plume != 0:
            # default value
            cloud_top_level[0, 0][plume] = k_end - 3

        if plume != 0 and error_code[0, 0][plume] == 0:
            start_level: IntFieldIJ = updraft_lfc_level[0, 0][plume]

    with computation(PARALLEL), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if K <= start_level:
                cloud_moist_static_energy_forced_transported = moist_static_energy_origin_level_forced

    with computation(FORWARD), interval(1, None):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if K > start_level and K < k_end - 2:
                dz = (
                    geopotential_height_cloud_levels_forced
                    - geopotential_height_cloud_levels_forced[0, 0, -1]
                )

                cloud_moist_static_energy_forced_transported = (
                    (1.0 - 0.5 * entrainment_rate[0, 0, -1][plume] * dz)
                    * cloud_moist_static_energy_forced_transported[0, 0, -1]
                    + entrainment_rate[0, 0, -1][plume]
                    * dz
                    * environment_moist_static_energy_forced[0, 0, -1]
                ) / (1.0 + 0.5 * entrainment_rate[0, 0, -1][plume] * dz)

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            # set up mask for next computation
            stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(...):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if K > start_level and K < k_end - 2 and stop_computation == False:
                # find the height where the parcel is no longer saturated
                if (
                    cloud_moist_static_energy_forced_transported
                    < environment_saturation_moist_static_energy_cloud_levels_forced
                ):
                    cloud_top_level[0, 0][plume] = K - 1
                    stop_computation = True

    with computation(FORWARD), interval(0, 1):
        if plume != 0 and error_code[0, 0][plume] == 0:
            if cloud_top_level[0, 0][plume] <= updraft_lfc_level[0, 0][plume] + 1:
                error_code[0, 0][plume] = 41

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and plume != 0 and OVERSHOOT > 0:
            z_overshoot = (1.0 + OVERSHOOT) * geopotential_height_cloud_levels_forced.at(
                K=cloud_top_level[0, 0][plume]
            )

    with computation(FORWARD), interval(0, 1):
        # set up mask for next computation
        stop_computation: BoolFieldIJ = False

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and plume != 0 and OVERSHOOT > 0:
            if K >= cloud_top_level[0, 0][plume] and K < k_end - 1 and stop_computation == False:
                if z_overshoot < geopotential_height_cloud_levels_forced:
                    cloud_top_level[0, 0][plume] = min(K - 1, k_end - 3)
                    stop_computation = True


def cloud_top_checks(
    cloud_top_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    geopotential_height_cloud_levels: FloatField,
    error_code: IntFieldIJ_Plume,
    last_error_code: IntFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    MINIMUM_DEPTH: Float,
    plume: Int,
):
    """Perform a series of checks on the initial estimate of cloud top level.

    Args:
        cloud_top_level (IntFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        geopotential_height_cloud_levels (FloatField)
        error_code (IntFieldIJ_Plume)
        last_error_code (IntFieldIJ)
        updraft_lfc_level (IntFieldIJ_Plume)
        MINIMUM_DEPTH (Float)
        plume (Int)
    """
    # check if cloud_top_level is too low for deep convection
    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.DEEP and error_code[0, 0][plume] == 0:
            if p_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume], ddim=[plume]) > 400:
                error_code[0, 0][plume] = 22

    # check if cloud_top_level is too high for mid convection
    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.MID and error_code[0, 0][plume] == 0:
            if p_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume], ddim=[plume]) < 400:
                error_code[0, 0][plume] = 22

    # check if cloud_top_level is too high for shallow convection
    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.SHALLOW and error_code[0, 0][plume] == 0:
            if p_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume], ddim=[plume]) < 400:
                error_code[0, 0][plume] = 23

    # avoid double-counting plumes
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 and last_error_code == 0:
            if p_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume], ddim=[plume]) < 700:
                error_code[0, 0][plume] = 27

    # last checks for cloud_top_level
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if (
                geopotential_height_cloud_levels.at(K=cloud_top_level[0, 0][plume])
                - geopotential_height_cloud_levels.at(K=updraft_lfc_level[0, 0][plume])
                < MINIMUM_DEPTH
            ):
                error_code[0, 0][plume] = 5
        if cloud_top_level[0, 0][plume] <= updraft_lfc_level[0, 0][plume]:
            error_code[0, 0][plume] = 5


class MaximumUpdraftOriginLevel:
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
        self._find_maximum_updraft_origin_level = stencil_factory.from_dims_halo(
            func=find_maximum_updraft_origin_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_maximum_updraft_origin_level(
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            error_code=state.output.error_code,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            MAX_UPDRAFT_ORIGIN_HEIGHT=plume_dependent_constants.MAX_UPDRAFT_ORIGIN_HEIGHT,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class DowndraftDetrainmentLevel:
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
        self._find_detrainmet_start_level = stencil_factory.from_dims_halo(
            func=find_detrainmet_start_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_detrainmet_start_level(
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            topography_height_no_negative=state.input_output.topography_height_no_negative,
            error_code=state.output.error_code,
            detrainment_start_level=locals.detrainment_start_level,
            DETRAINMENT_CRITICAL_DEPTH=plume_dependent_constants.DETRAINMENT_CRITICAL_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class HighestMoistStaticEnergyLevel:
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
        self._find_highest_moist_static_energy_level = stencil_factory.from_dims_halo(
            func=find_highest_moist_static_energy_level,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._find_highest_moist_static_energy_level(
            moist_static_energy=locals.environment_moist_static_energy_cloud_levels_forced,
            error_code=state.output.error_code,
            maximum_updraft_origin_level=locals.maximum_updraft_origin_level,
            updraft_origin_level=state.output.updraft_origin_level,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class CloudTop:
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
        self._updraft_rates_pdf = stencil_factory.from_dims_halo(
            func=get_cloud_top,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        )

        self._cloud_top_checks = stencil_factory.from_dims_halo(
            func=cloud_top_checks,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):

        if self.cumulus_parameterization_config.OVERSHOOT != 0:
            ndsl_log.warning(
                " GF2020 cumulus parameterization called CloudTop with "
                "untested OVERSHOOT option. Running untested code... proceed with caution"
            )

        self._updraft_rates_pdf(
            entrainment_rate=state.output.entrainment_rate,
            moist_static_energy=locals.environment_moist_static_energy_forced,
            saturation_moist_static_energy=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
            moist_static_energy_origin_level=locals.moist_static_energy_origin_level_forced,
            updraft_lfc_level=state.output.updraft_lfc_level,
            geopotential_height=locals.geopotential_height_cloud_levels_forced,
            cloud_moist_static_energy=locals.cloud_moist_static_energy_forced_transported,
            error_code=state.output.error_code,
            cloud_top_level=state.output.cloud_top_level,
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._cloud_top_checks(
            cloud_top_level=state.output.cloud_top_level,
            p=state.output.p_cloud_levels_forced,
            geopotential_height=locals.geopotential_height_cloud_levels,
            error_code=state.output.error_code,
            last_error_code=state.input.last_error_code,
            updraft_lfc_level=state.output.updraft_lfc_level,
            MINIMUM_DEPTH=plume_dependent_constants.MINIMUM_DEPTH,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
