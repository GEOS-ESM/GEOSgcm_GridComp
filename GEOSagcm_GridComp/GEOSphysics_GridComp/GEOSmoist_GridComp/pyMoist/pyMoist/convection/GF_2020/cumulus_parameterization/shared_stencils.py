from ndsl.dsl.typing import FloatField, FloatFieldIJ, IntFieldIJ, Int, BoolFieldIJ
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import FloatField_Plume, IntFieldIJ_Plume
from gt4py.cartesian.gtscript import (
    FORWARD,
    PARALLEL,
    BACKWARD,
    computation,
    interval,
    K,
    sqrt,
)
import pyMoist.constants as constants
from ndsl.stencils.column_operations import column_max


def unknown_find_level(
    array: FloatField,
    start_index: IntFieldIJ_Plume,
    end_index: IntFieldIJ,
    out_index: IntFieldIJ_Plume,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    Details/purpose unknown

    Args:
        array (in): array to be analyzed
        start_index (in): start index for analysis
        end_index (in): end index for analysis
        out_index (out): output of analysis
        error_code (in): field for stopping flow through the scheme and tracking errors
        plume (in): specifies the current plume
    """
    with computation(FORWARD), interval(0, 1):
        out_index[0, 0][plume] = start_index[0, 0][plume]

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            x = array.at(K=start_index[0, 0][plume])
            stop_index = max(start_index[0, 0][plume] + 1, end_index)

            level = start_index[0, 0][plume] + 1
            while level >= start_index[0, 0][plume] + 1 and level <= stop_index:
                if array.at(K=level) < x:
                    x = array.at(K=level)
                    out_index[0, 0][plume] = level
                level += 1


def updraft_vertical_velocity(
    vertical_velocity_3d: FloatField,
    vertical_velocity_2d: FloatFieldIJ,
    convective_scale_velocity: FloatFieldIJ,
    entrainment_rate: FloatField_Plume,
    detrainment_function_updraft: FloatField,
    geopotential_height_cloud_levels_forced: FloatField,
    t_cloud_levels_forced: FloatField,
    updraft_column_temperature_forced: FloatField,
    cloud_total_water_after_entrainment_forced: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    vapor_forced: FloatField,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    from __externals__ import ZERO_DIFF, k_end

    # internal constants
    with computation(FORWARD), interval(0, 1):
        ctea: FloatFieldIJ = 1.0 / 3.0
        cteb: FloatFieldIJ = 2.0
        visc: FloatFieldIJ = 2000.0
        eps: FloatFieldIJ = 0.622
        f: FloatFieldIJ = 2.0
        C_d: FloatFieldIJ = 0.506
        gam: FloatFieldIJ = 0.5
        beta: FloatFieldIJ = 1.875
        smooth: BoolFieldIJ = True
        n_smooth: IntFieldIJ = 1

        if ZERO_DIFF == 1:
            ftun1: FloatFieldIJ = 1.0
            ftun2: FloatFieldIJ = 0.5
        else:
            ftun1: FloatFieldIJ = 0.25
            ftun2: FloatFieldIJ = 1.0

    # initialize arrays to zero
    with computation(PARALLEL), interval(...):
        vertical_velocity_3d = 0.0

    with computation(FORWARD), interval(0, 1):
        vertical_velocity_2d = 0.0

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= updraft_lfc_level[0, 0][plume]:
                vertical_velocity_3d = max(1.0, convective_scale_velocity**2)

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K >= updraft_lfc_level[0, 0][plume] and K <= cloud_top_level[0, 0][plume]:
                dz = (
                    geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
                )

                t_ve = 0.5 * (
                    t_cloud_levels_forced * (1.0 + (vapor_forced / eps) / (1.0 + vapor_forced))
                    + t_cloud_levels_forced[0, 0, 1]
                    * (1.0 + (vapor_forced[0, 0, 1] / eps) / (1.0 + vapor_forced[0, 0, 1]))
                )

                t_v = 0.5 * (
                    updraft_column_temperature_forced
                    * (
                        1.0
                        + (cloud_total_water_after_entrainment_forced / eps)
                        / (1.0 + cloud_total_water_after_entrainment_forced)
                    )
                    + updraft_column_temperature_forced[0, 0, 1]
                    * (
                        1.0
                        + (cloud_total_water_after_entrainment_forced[0, 0, 1] / eps)
                        / (1.0 + cloud_total_water_after_entrainment_forced[0, 0, 1])
                    )
                )

                bu = constants.MAPL_GRAV * (
                    (t_v - t_ve) / t_ve
                    - ftun2
                    * 0.50
                    * (
                        cloud_liquid_after_rain_forced[0, 0, 1][plume]
                        + cloud_liquid_after_rain_forced[0, 0, 0][plume]
                    )
                )

                dw1 = 2.0 / (f * (1.0 + gam)) * bu * dz
                if ZERO_DIFF == 1:
                    kx = max(entrainment_rate[0, 0, 0][plume], detrainment_function_updraft) * dz
                else:
                    kx = (
                        (1.0 + beta * C_d)
                        * max(entrainment_rate[0, 0, 0][plume], detrainment_function_updraft)
                        * dz
                        * ftun1
                    )

                dw2 = (vertical_velocity_3d) - 2.0 * kx * (vertical_velocity_3d)

                vertical_velocity_3d[0, 0, 1] = (dw1 + dw2) / (1.0 + kx)

                if vertical_velocity_3d[0, 0, 1] < 0.0:
                    vertical_velocity_3d[0, 0, 1] = 0.5 * vertical_velocity_3d

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if smooth == True:  # noqa
                if ZERO_DIFF == 1:
                    if K <= cloud_top_level[0, 0][plume] - 2:
                        nvs: IntFieldIJ = 0
                        vs: FloatFieldIJ = 0.0
                        level_inner_loop = max(K - n_smooth, 0)
                        while level_inner_loop <= min(K + n_smooth, k_end - 1):
                            nvs = nvs + 1
                            vs = vs + vertical_velocity_3d.at(K=level_inner_loop)
                            level_inner_loop += 1
                        vertical_velocity_3d = vs / (1.0e-16 + nvs)
                else:
                    if K <= cloud_top_level[0, 0][plume] + 1:
                        vs: FloatFieldIJ = 0.0
                        dz1m: FloatFieldIJ = 0.0
                        level_inner_loop = max(K - n_smooth, 0)
                        while level_inner_loop <= min(K + n_smooth, k_end - 1):
                            dz = geopotential_height_cloud_levels_forced.at(
                                K=level_inner_loop + 1
                            ) - geopotential_height_cloud_levels_forced.at(K=level_inner_loop)
                            vs = vs + dz * vertical_velocity_3d.at(K=level_inner_loop)
                            dz1m = dz1m + dz
                            level_inner_loop += 1
                        vertical_velocity_3d = vs / (1.0e-16 + dz1m)

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            # convert to vertical velocity
            vertical_velocity_3d = sqrt(max(0.1, vertical_velocity_3d))

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            max_val, _ = column_max(vertical_velocity_3d, 0, k_end)
            if max_val < 1.0:
                error_code[0, 0][plume] = 54

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 or error_code[0, 0][plume] == 54:
            # sanity check
            if vertical_velocity_3d < 1.0:
                vertical_velocity_3d = 1.0
            if vertical_velocity_3d > 20.0:
                vertical_velocity_3d = 20.0

            if K >= cloud_top_level[0, 0][plume] + 1 and ZERO_DIFF == 0:
                vertical_velocity_3d = 0.1

    # get the column average vertical velocity
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 or error_code[0, 0][plume] == 54:
            if K >= updraft_lfc_level[0, 0][plume] and K <= cloud_top_level[0, 0][plume]:
                dz = (
                    geopotential_height_cloud_levels_forced[0, 0, 1] - geopotential_height_cloud_levels_forced
                )
                vertical_velocity_2d = vertical_velocity_2d + vertical_velocity_3d * dz

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0 or error_code[0, 0][plume] == 54:
            vertical_velocity_2d = vertical_velocity_2d / (
                geopotential_height_cloud_levels_forced.at(K=cloud_top_level[0, 0][plume] + 1)
                - geopotential_height_cloud_levels_forced.at(K=updraft_lfc_level[0, 0][plume])
                + 1.0e-16
            )
            vertical_velocity_2d = max(1.0, vertical_velocity_2d)


def tridiag(
    m: IntFieldIJ_Plume,
    a: FloatField,
    b: FloatField,
    c: FloatField,
    f: FloatField,
    error_code: IntFieldIJ_Plume,
    plume: Int,
):
    """
    this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
    an updated "f" at time t+1 is the output

    Args:
        m
        a
        b
        c
        f
    """
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and K == m[0, 0][plume]:
            c = 0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            q = -c / b
            f = f / b

    with computation(FORWARD), interval(1, None):
        if error_code[0, 0][plume] == 0 and K <= m[0, 0][plume]:
            p = 1.0 / (b + a * q[0, 0, -1])
            q = -c * p
            f = p * (f - a * f[0, 0, -1])

    with computation(BACKWARD), interval(...):
        if error_code[0, 0][plume] == 0 and K < m[0, 0][plume]:
            f = f + q * f[0, 0, 1]
