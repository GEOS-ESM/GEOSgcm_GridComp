from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, IntFieldIJ, Int, IntField
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    IntFieldIJ_Plume,
    FloatFieldIJ_Plume,
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    get_cloud_boundary_conditions,
    liquid_fraction,
)


def updraft_moisture_light(
    start_level: IntFieldIJ,
    error_code: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_vapor_mixing_ratio_forced: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    precipitable_water_updraft_forced: FloatField_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    cloud_moist_static_energy_forced: FloatField,
    unspecifid_temperature: FloatField,
    ocean_fraction: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    p_forced: FloatField,
    cloud_top_level: IntFieldIJ_Plume,
    buoyancy_forced: FloatField,
    cloud_liquid_before_rain_forced: FloatField,
    t_cloud_levels: FloatField,
    vapor_forced: FloatField,
    gamma_cloud_levels_forced: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
    environment_saturation_mixing_ratio_cloud_levels_forced: FloatField,
    updraft_origin_level: FloatFieldIJ_Plume,
    vapor_cloud_levels_forced: FloatField,
    vapor_excess: FloatFieldIJ,
    mass_entrainment_updraft: FloatField,
    mass_detrainment_updraft: FloatField,
    psum: FloatFieldIJ,
    psumh: FloatFieldIJ,
    c1d: FloatField,
    add_buoyancy: FloatFieldIJ,
    AVERAGE_LAYER_DEPTH: Float,
    C0: Float,
    plume: Int,
):
    from __externals__ import (
        BOUNDARY_CONDITION_METHOD,
        k_end,
        MELT_GLAC,
        FRAC_MODIS,
        QRC_CRIT_OCN,
        QRC_CRIT_LND,
    )

    with computation(FORWARD), interval(0, 1):
        total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
        psum = 0.0
        psumh = 0.0

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break (this is never touched)
        dummy_field_no_read = 0.0 + BOUNDARY_CONDITION_METHOD

    with computation(PARALLEL), interval(...):
        precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
        cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0
        cloud_liquid_before_rain_forced = 0.0
        unspecifid_temperature = t_cloud_levels
        cloud_vapor_mixing_ratio_forced = vapor_cloud_levels_forced

    # with computation(FORWARD), interval(0, 1):
    #     if error_code[0, 0][plume] == 0:
    #         # get boundary condition for cloud_vapor_mixing_ratio_forced
    #         vapor_boundary_condition = get_cloud_boundary_conditions(
    #             field=vapor_cloud_levels_forced,
    #             scalar_perturbation=0,
    #             p=p_forced,
    #             updraft_origin_level=updraft_origin_level[0, 0][plume],
    #             ocean_fraction=ocean_fraction,
    #             BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
    #             AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
    #             k_end=k_end,
    #             compute_perturbation=False,
    #             perturbation_field=dummy_field_no_read,
    #         )
    #         cloud_vapor_mixing_ratio_forced = (
    #             vapor_boundary_condition
    #             + vapor_excess
    #             + 0.5 * add_buoyancy / cumulus_parameterization_constants.XLV
    #         )

    # with computation(FORWARD), interval(...):
    #     if error_code[0, 0][plume] == 0:
    #         if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:
    #             dz = (
    #                 geopotential_height_cloud_levels_forced
    #                 - geopotential_height_cloud_levels_forced[0, 0, -1]
    #             )

    #             # saturation in cloud, this is what is allowed to be in it
    #             saturation_cloud_liquid_rain_forced = (
    #                 environment_saturation_mixing_ratio_cloud_levels_forced
    #                 + (1 / cumulus_parameterization_constants.XLV)
    #                 * (gamma_cloud_levels_forced / (1.0 + gamma_cloud_levels_forced))
    #                 * buoyancy_forced
    #             )

    #             # 1. steady state plume equation, for what could be in cloud without condensation
    #             denom = (
    #                 normalized_massflux_updraft_forced[0, 0, -1][plume]
    #                 - 0.5 * mass_detrainment_updraft[0, 0, -1]
    #                 + mass_entrainment_updraft[0, 0, -1]
    #             )
    #             if denom > 0.0:
    #                 cloud_vapor_mixing_ratio_forced = (
    #                     cloud_vapor_mixing_ratio_forced[0, 0, -1]
    #                     * normalized_massflux_updraft_forced[0, 0, -1][plume]
    #                     - 0.5 * mass_detrainment_updraft[0, 0, -1] * cloud_vapor_mixing_ratio_forced[0, 0, -1]
    #                     + mass_entrainment_updraft[0, 0, -1] * vapor_forced[0, 0, -1]
    #                 ) / denom
    #                 if K == start_level + 1:
    #                     cloud_vapor_mixing_ratio_forced = (
    #                         cloud_vapor_mixing_ratio_forced
    #                         + vapor_excess * mass_entrainment_updraft[0, 0, -1] / denom
    #                     )
    #             else:
    #                 cloud_vapor_mixing_ratio_forced = cloud_vapor_mixing_ratio_forced[0, 0, -1]

    #             # total condensed water before rainout
    #             cloud_liquid_before_rain_forced = max(
    #                 0.0, cloud_vapor_mixing_ratio_forced - saturation_cloud_liquid_rain_forced
    #             )
    #             # updraft temp
    #             unspecifid_temperature = (1.0 / cumulus_parameterization_constants.CP) * (
    #                 cloud_moist_static_energy_forced
    #                 - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
    #                 - cumulus_parameterization_constants.XLV * saturation_cloud_liquid_rain_forced
    #             )

    #             # add glaciation effect on the MSE
    #             if MELT_GLAC == 1:
    #                 delta_cloud_mse_glac = (
    #                     cloud_liquid_before_rain_forced
    #                     * (
    #                         1.0
    #                         - liquid_fraction(
    #                             unspecifid_temperature, convection_fraction, surface_type, FRAC_MODIS
    #                         )
    #                     )
    #                     * cumulus_parameterization_constants.XLF
    #                 )

    #                 unspecifid_temperature = (
    #                     unspecifid_temperature
    #                     + (1.0 / cumulus_parameterization_constants.CP) * delta_cloud_mse_glac
    #                 )

    #             if C0 < 1.0e-6:
    #                 cx0 = 0.0
    #             else:
    #                 cx0 = (c1d + C0) * dz

    #             cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced / (1.0 + cx0)
    #             min_liq = ocean_fraction * QRC_CRIT_OCN + (1.0 - ocean_fraction) * QRC_CRIT_LND
    #             precipitable_water_updraft_forced[0, 0, 0][plume] = cx0 * max(
    #                 0.0, cloud_liquid_after_rain_forced[0, 0, 0][plume] - min_liq
    #             )  # units kg[rain]/kg[air]

    #             # convert pw to normalized pw
    #             precipitable_water_updraft_forced[0, 0, 0][plume] = (
    #                 precipitable_water_updraft_forced[0, 0, 0][plume]
    #                 * normalized_massflux_updraft_forced[0, 0, 0][plume]
    #             )

    #             # total water (vapor + condensed) in updraft after the rainout
    #             cloud_vapor_mixing_ratio_forced = cloud_liquid_after_rain_forced[0, 0, 0][plume] + min(
    #                 cloud_vapor_mixing_ratio_forced, saturation_cloud_liquid_rain_forced
    #             )

    # with computation(PARALLEL), interval(...):
    #     if error_code[0, 0][plume] == 0:
    #         # get back water vapor qc
    #         if K <= cloud_top_level[0, 0][plume]:
    #             cloud_vapor_mixing_ratio_forced = (
    #                 cloud_vapor_mixing_ratio_forced - cloud_liquid_after_rain_forced[0, 0, 0][plume]
    #             )


def in_cloud_updraft_air_temperature(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_incloud_air_temp_forced: FloatField,
    local_hcdo: FloatField,
    local_geopotential_height_cloud_levels_forced: FloatField,
    local_incloud_water_vapor_mixing_ratio_forced: FloatField,
    local_t_cloud_levels_forced: FloatField,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            local_incloud_air_temp_forced = (1.0 / cumulus_parameterization_constants.CP) * (
                local_hcdo
                - constants.MAPL_GRAV * local_geopotential_height_cloud_levels_forced
                - cumulus_parameterization_constants.XLV * local_incloud_water_vapor_mixing_ratio_forced
            )
    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            local_incloud_air_temp_forced = local_t_cloud_levels_forced


def get_melting_profile(
    error_code: IntFieldIJ_Plume,
    plume: Int,
    local_melting_layer: FloatField,
    local_partition_liquid_ice: FloatField,
    p_cloud_levels_forced: FloatField_Plume,
    precipitable_water_updraft_forced: FloatField_Plume,
    local_melting: FloatField,
):
    from __externals__ import k_end, MELT_GLAC

    with computation(FORWARD), interval(...):
        ktf = k_end - 1
        if MELT_GLAC == 1 and plume == 2:
            pwo_solid_phase = 0.0
            pwo_eff = 0.0
            local_melting = 0.0

            if error_code[0, 0][plume] > 0:
                local_melting = 0.0

            total_pwo_solid_phase: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(...):
        if MELT_GLAC == 1 and plume == 2:
            if K <= ktf - 1:
                if error_code[0, 0][plume] == 0:
                    dp = 100.0 * (
                        p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
                    )

                    pwo_eff = 0.5 * (
                        precipitable_water_updraft_forced[0, 0, 0][plume]
                        + precipitable_water_updraft_forced[0, 0, 1][plume]
                    )

                    pwo_solid_phase = (1.0 - local_partition_liquid_ice) * pwo_eff

                    total_pwo_solid_phase = total_pwo_solid_phase + pwo_solid_phase * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == 1 and plume == 2:
            if K <= ktf:
                if error_code[0, 0][plume] == 0:
                    local_melting = local_melting_layer * (
                        total_pwo_solid_phase
                        / (
                            100
                            * (
                                p_cloud_levels_forced.at(K=0, ddim=[plume])
                                - p_cloud_levels_forced.at(K=ktf, ddim=[plume])
                            )
                            / constants.MAPL_GRAV
                        )
                    )
        else:
            local_melting = 0.0
