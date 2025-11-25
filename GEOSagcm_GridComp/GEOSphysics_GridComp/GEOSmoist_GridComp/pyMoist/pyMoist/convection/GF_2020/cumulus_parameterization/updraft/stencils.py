from ndsl.dsl.typing import (
    FloatField,
    FloatFieldIJ,
    Float,
    IntFieldIJ,
    Int,
    Bool,
    IntField,
)
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)


def in_cloud_updraft_air_temperature(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_incloud_air_temp: FloatField,
    local_cloud_moist_static_energy_forced: FloatField,
    local_geopotential_height_cloud_levels_forced: FloatField,
    local_incloud_water_vapor_mixing_ratio: FloatField,
    local_t_cloud_levels_forced: FloatField,
):
    from __externals__ import FIRST_GUESS_W

    with computation(PARALLEL), interval(0, -1):
        if FIRST_GUESS_W == 0:
            if error_code[0, 0][plume] == 0:
                local_incloud_air_temp = (
                    1.0 / cumulus_parameterization_constants.CP
                ) * (
                    local_cloud_moist_static_energy_forced
                    - constants.MAPL_GRAV
                    * local_geopotential_height_cloud_levels_forced
                    - cumulus_parameterization_constants.XLV
                    * local_incloud_water_vapor_mixing_ratio
                )

    with computation(PARALLEL), interval(-1, None):
        if FIRST_GUESS_W == 0:
            if error_code[0, 0][plume] == 0:
                local_incloud_air_temp = local_t_cloud_levels_forced

    with computation(PARALLEL), interval(...):
        if FIRST_GUESS_W == 0:
            if error_code[0, 0][plume] != 0:
                local_incloud_air_temp = local_t_cloud_levels_forced


def cup_up_aa0(
    local_buoyancy: FloatField,
    local_gamma_cloud_levels: FloatField,
    cloud_top: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    local_geopotential_height_cloud_levels: FloatField,
    local_t_cloud_levels: FloatField,
    local_normalized_massflux_updraft: FloatField,
    local_integ: IntFieldIJ,
    local_integ_interval: IntFieldIJ,
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_cloud_work_function: FloatFieldIJ,
):
    from __externals__ import k_start

    with computation(FORWARD), interval(...):
        local_cloud_work_function = 0.0

        if local_integ == 1:
            if local_integ_interval == cumulus_parameterization_constants.BL:
                kbeg = k_start
                kend = updraft_lfc_level[0, 0][plume] - 2
            elif local_integ_interval == cumulus_parameterization_constants.CIN:
                kbeg = updraft_origin_level[0, 0][plume] - 1
                kend = updraft_lfc_level[0, 0][plume] - 2

        else:
            kbeg = updraft_lfc_level[0, 0][plume] - 1
            kend = cloud_top[0, 0][plume] - 1

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= kbeg and K <= kend:
                dz = (
                    local_geopotential_height_cloud_levels[0, 0, 1]
                    - local_geopotential_height_cloud_levels
                )
                aa_1 = (
                    local_normalized_massflux_updraft
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * local_t_cloud_levels)
                    )
                    * local_buoyancy
                    / (1.0 + local_gamma_cloud_levels)
                )
                aa_2 = (
                    local_normalized_massflux_updraft[0, 0, 1]
                    * (
                        constants.MAPL_GRAV
                        / (
                            cumulus_parameterization_constants.CP
                            * local_t_cloud_levels[0, 0, 1]
                        )
                    )
                    * local_buoyancy[0, 0, 1]
                    / (1.0 + local_gamma_cloud_levels[0, 0, 1])
                )
                da = 0.5 * (aa_1 + aa_2) * dz

                local_cloud_work_function = local_cloud_work_function + da


def cloud_work_function_zero(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_cloud_work_function: FloatFieldIJ,
):
    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if local_cloud_work_function == 0.0:
                error_code[0, 0][plume] = 17
                # ierrc[0,0][plume]="cloud work function zero"
