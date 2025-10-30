from ndsl.dsl.gt4py import computation, interval, FORWARD, K
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, BoolFieldIJ
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import IntFieldIJ_Plume


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
