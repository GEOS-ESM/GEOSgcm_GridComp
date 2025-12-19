from ndsl import StencilFactory, QuantityFactory
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
from ndsl.dsl.typing import (
    FloatField,
    FloatFieldIJ,
    Float,
    IntFieldIJ,
    Int,
)
from ndsl.dsl.gt4py import computation, PARALLEL, interval, FORWARD, K, sqrt, exp
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
import pyMoist.constants as constants
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import (
    get_cloud_boundary_conditions,
    liquid_fraction,
)


def updraft_moisture(
    start_level: IntFieldIJ,
    error_code: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_vapor_mixing_ratio_forced: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    precipitable_water_updraft_forced: FloatField_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    cloud_moist_static_energy_forced: FloatField,
    miscellaneous_temperature: FloatField,
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
    normalized_massflux_updraft_forced: FloatField,
    environment_saturation_mixing_ratio_cloud_levels_forced: FloatField,
    updraft_origin_level: IntFieldIJ_Plume,
    vapor_cloud_levels_forced: FloatField,
    vapor_excess: FloatFieldIJ,
    ccn: FloatFieldIJ,
    mass_entrainment_updraft: FloatField,
    mass_detrainment_updraft: FloatField,
    psum: FloatFieldIJ,
    psumh: FloatFieldIJ,
    c1d: FloatField,
    add_buoyancy: FloatFieldIJ,
    vertical_velocity_3d: FloatField,
    C0: Float,
    AVERAGE_LAYER_DEPTH: Float,
    plume: Int,
):
    from __externals__ import (
        k_end,
        BOUNDARY_CONDITION_METHOD,
        USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
        AUTOCONV,
        CRITICAL_MIXING_RATIO_OVER_OCEAN,
        CRITICAL_MIXING_RATIO_OVER_LAND,
        FRAC_MODIS,
        ZERO_DIFF,
    )

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break (this is never touched)
        dummy_field_no_read = 0.0

    with computation(FORWARD), interval(0, 1):
        # internal constants
        BDISPM: FloatFieldIJ = 0.366  # berry--size dispersion (maritime)
        BDISPC: FloatFieldIJ = 0.146  # berry--size dispersion (continental)
        T_BF: FloatFieldIJ = 268.16
        T_ICE_BF: FloatFieldIJ = 235.16
        RK: FloatFieldIJ = 3.0
        XEXP: FloatFieldIJ = 2.0

    with computation(FORWARD), interval(0, 1):
        total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
        psum = 0.0
        psumh = 0.0

    with computation(PARALLEL), interval(...):
        precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
        miscellaneous_temperature = t_cloud_levels
        cloud_liquid_before_rain_forced = 0.0
        cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0  # liq/ice water
        cloud_vapor_mixing_ratio_forced = 0.0  # total water: liq/ice = vapor water

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            vapor_source: FloatFieldIJ = get_cloud_boundary_conditions(
                field=vapor_cloud_levels_forced,
                scalar_perturbation=0,
                p=p_forced,
                updraft_origin_level=updraft_origin_level[0, 0][plume],
                ocean_fraction=ocean_fraction,
                BOUNDARY_CONDITION_METHOD=BOUNDARY_CONDITION_METHOD,
                AVERAGE_LAYER_DEPTH=AVERAGE_LAYER_DEPTH,
                k_end=k_end,
                compute_perturbation=False,
                perturbation_field=dummy_field_no_read,
            )

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= start_level:
            cloud_vapor_mixing_ratio_forced = (
                vapor_source + vapor_excess + 0.5 * add_buoyancy / cumulus_parameterization_constants.XLV
            )
            cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        if (
            error_code[0, 0][plume] == 0 and USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES == 1 and plume == 0
        ):  # only for shallow plume
            get_delmix_implementation_here = True

    with computation(PARALLEL), interval(...):
        # initalize mask to stop computation in the next block
        stop_current_index = False

    with computation(FORWARD), interval(1, None):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:
                dz = (
                    geopotential_height_cloud_levels_forced
                    - geopotential_height_cloud_levels_forced[0, 0, -1]
                )
                # saturation  in cloud, this is what is allowed to be in it
                qrch = (
                    environment_saturation_mixing_ratio_cloud_levels_forced
                    + (1.0 / cumulus_parameterization_constants.XLV)
                    * (gamma_cloud_levels_forced / (1.0 + gamma_cloud_levels_forced))
                    * buoyancy_forced
                )

                #    1. steady state plume equation, for what could
                #       be in cloud without condensation
                denom = (
                    normalized_massflux_updraft_forced[0, 0, -1]
                    - 0.5 * mass_detrainment_updraft[0, 0, -1]
                    + mass_entrainment_updraft[0, 0, -1]
                )

                if denom > 0.0:
                    cloud_vapor_mixing_ratio_forced = (
                        cloud_vapor_mixing_ratio_forced[0, 0, -1]
                        * normalized_massflux_updraft_forced[0, 0, -1]
                        - 0.5 * mass_detrainment_updraft[0, 0, -1] * cloud_vapor_mixing_ratio_forced[0, 0, -1]
                        + mass_entrainment_updraft[0, 0, -1] * vapor_forced[0, 0, -1]
                    ) / denom

                    if K == start_level + 1:
                        cloud_vapor_mixing_ratio_forced = (
                            cloud_vapor_mixing_ratio_forced
                            + vapor_excess * mass_entrainment_updraft[0, 0, -1] / denom
                        )
                    # assuming no liq/ice water in the environment
                    cloud_liquid_after_rain_forced[0, 0, 0][plume] = (
                        cloud_liquid_after_rain_forced[0, 0, -1][plume]
                        * normalized_massflux_updraft_forced[0, 0, -1]
                        - 0.5
                        * mass_detrainment_updraft[0, 0, -1]
                        * cloud_liquid_after_rain_forced[0, 0, -1][plume]
                    ) / denom

                else:
                    cloud_vapor_mixing_ratio_forced = cloud_vapor_mixing_ratio_forced[0, 0, -1]
                    cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_after_rain_forced[0, 0, -1][
                        plume
                    ]

                # updraft temp
                miscellaneous_temperature = (1.0 / cumulus_parameterization_constants.CP) * (
                    cloud_moist_static_energy_forced
                    - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
                    - cumulus_parameterization_constants.XLV * qrch
                )

                # total condensed water before rainout
                cloud_liquid_before_rain_forced = max(0.0, cloud_vapor_mixing_ratio_forced - qrch)

                cloud_liquid_after_rain_forced[0, 0, 0][plume] = min(
                    cloud_liquid_before_rain_forced, cloud_liquid_after_rain_forced[0, 0, 0][plume]
                )

                # production term => condensation/diffusional growth
                cup = (
                    max(
                        0.0,
                        cloud_vapor_mixing_ratio_forced
                        - qrch
                        - cloud_liquid_after_rain_forced[0, 0, 0][plume],
                    )
                    / dz
                )

                if C0 < 1.0e-6:
                    cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                    cloud_vapor_mixing_ratio_forced = cloud_liquid_after_rain_forced[0, 0, 0][plume] + min(
                        cloud_vapor_mixing_ratio_forced, qrch
                    )
                    total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
                    psum = psum + cloud_liquid_before_rain_forced * normalized_massflux_updraft_forced * dz

                    stop_current_index = True

                if stop_current_index == False:
                    if AUTOCONV == 1:
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        cx0 = (c1d + C0) * dz
                        cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced / (
                            1.0 + cx0
                        )
                        precipitable_water_updraft_forced[0, 0, 0][plume] = cx0 * max(
                            0.0, cloud_liquid_after_rain_forced[0, 0, 0][plume] - min_liq
                        )  # units kg[rain]/kg[air]
                        # convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                        precipitable_water_updraft_forced[0, 0, 0][plume] = (
                            precipitable_water_updraft_forced[0, 0, 0][plume]
                            * normalized_massflux_updraft_forced
                        )

                    elif AUTOCONV == 2:
                        # this is similar to AUTOCONV == 1 with temperature dependence
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        cx0 = (
                            (c1d + C0)
                            * dz
                            * liquid_fraction(
                                miscellaneous_temperature, convection_fraction, surface_type, FRAC_MODIS
                            )
                        )
                        cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced / (
                            1.0 + cx0
                        )
                        precipitable_water_updraft_forced[0, 0, 0][plume] = cx0 * max(
                            0.0, cloud_liquid_after_rain_forced[0, 0, 0][plume] - min_liq
                        )  # units kg[rain]/kg[air]
                        # --- convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                        precipitable_water_updraft_forced[0, 0, 0][plume] = (
                            precipitable_water_updraft_forced[0, 0, 0][plume]
                            * normalized_massflux_updraft_forced
                        )

                    elif AUTOCONV == 3:
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        if cloud_liquid_before_rain_forced <= min_liq:
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
                        else:
                            cx0 = C0 * liquid_fraction(
                                miscellaneous_temperature, convection_fraction, surface_type, FRAC_MODIS
                            )
                            cx0 = max(cx0, 0.50 * C0)
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_after_rain_forced[
                                0, 0, 0
                            ][plume] * exp(-cx0 * dz) + (cup / cx0) * (1.0 - exp(-cx0 * dz))
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = min(
                                cloud_liquid_before_rain_forced,
                                cloud_liquid_after_rain_forced[0, 0, 0][plume],
                            )
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                cloud_liquid_before_rain_forced
                                - cloud_liquid_after_rain_forced[0, 0, 0][plume]
                            )
                            # convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                precipitable_water_updraft_forced[0, 0, 0][plume]
                                * normalized_massflux_updraft_forced
                            )

                    elif AUTOCONV == 4:
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        if cloud_liquid_before_rain_forced <= min_liq:
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
                        else:
                            tem1 = liquid_fraction(
                                miscellaneous_temperature, convection_fraction, surface_type, FRAC_MODIS
                            )
                            cbf = 1.0
                            if miscellaneous_temperature < T_BF:
                                cbf = 1.0 + 0.5 * sqrt(
                                    min(max(T_BF - miscellaneous_temperature, 0.0), T_BF - T_ICE_BF)
                                )
                            qrc_crit_BF = ccn / cbf
                            cx0 = (
                                C0
                                * cbf
                                * (tem1 * 1.3 + (1.0 - tem1))
                                / (0.75 * min(15.0, max(vertical_velocity_3d, 1.0)))
                            )
                            # analytical solution
                            cx0 = cx0 * (
                                1.0
                                - exp(-((cloud_liquid_after_rain_forced[0, 0, 0][plume] / qrc_crit_BF) ** 2))
                            )
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_after_rain_forced[
                                0, 0, 0
                            ][plume] * exp(-cx0 * dz) + (cup / cx0) * (1.0 - exp(-cx0 * dz))
                            precipitable_water_updraft_forced[0, 0, 0][plume] = max(
                                cloud_liquid_before_rain_forced
                                - cloud_liquid_after_rain_forced[0, 0, 0][plume],
                                0.0,
                            )
                            # convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                precipitable_water_updraft_forced[0, 0, 0][plume]
                                * normalized_massflux_updraft_forced
                            )

                    elif AUTOCONV == 5:
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        if cloud_liquid_before_rain_forced <= min_liq:
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
                        else:
                            cx0 = (c1d + C0) * (
                                1.0
                                + 0.33
                                * liquid_fraction(
                                    miscellaneous_temperature, convection_fraction, surface_type, FRAC_MODIS
                                )
                            )
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_after_rain_forced[
                                0, 0, 0
                            ][plume] * exp(-cx0 * dz) + (cup / cx0) * (1.0 - exp(-cx0 * dz))
                            precipitable_water_updraft_forced[0, 0, 0][plume] = max(
                                0.0,
                                cloud_liquid_before_rain_forced
                                - cloud_liquid_after_rain_forced[0, 0, 0][plume],
                            )  # units kg[rain]/kg[air]
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = (
                                cloud_liquid_before_rain_forced
                                - precipitable_water_updraft_forced[0, 0, 0][plume]
                            )
                            # convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                precipitable_water_updraft_forced[0, 0, 0][plume]
                                * normalized_massflux_updraft_forced
                            )

                    elif AUTOCONV == 6:
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        if cloud_liquid_before_rain_forced <= min_liq:
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
                        else:
                            cx0 = (c1d + C0) * dz
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = (
                                cloud_liquid_before_rain_forced
                            ) * exp(-cx0)
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                cloud_liquid_before_rain_forced
                                - cloud_liquid_after_rain_forced[0, 0, 0][plume]
                            )
                            # convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                precipitable_water_updraft_forced[0, 0, 0][plume]
                                * normalized_massflux_updraft_forced
                            )

                    elif AUTOCONV == 7:
                        min_liq = (
                            ocean_fraction * CRITICAL_MIXING_RATIO_OVER_OCEAN
                            + (1.0 - ocean_fraction) * CRITICAL_MIXING_RATIO_OVER_LAND
                        )
                        if cloud_liquid_before_rain_forced <= min_liq:
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = 0.0
                        else:
                            cx0 = c1d + C0
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_after_rain_forced[
                                0, 0, 0
                            ][plume] * exp(-cx0 * dz) + (cup / cx0) * (1.0 - exp(-cx0 * dz))
                            precipitable_water_updraft_forced[0, 0, 0][plume] = max(
                                cloud_liquid_before_rain_forced
                                - cloud_liquid_after_rain_forced[0, 0, 0][plume],
                                0.0,
                            )
                            cloud_liquid_after_rain_forced[0, 0, 0][plume] = (
                                cloud_liquid_before_rain_forced
                                - precipitable_water_updraft_forced[0, 0, 0][plume]
                            )
                            # convert precipitable_water_updraft_forced to normalized precipitable_water_updraft_forced
                            precipitable_water_updraft_forced[0, 0, 0][plume] = (
                                precipitable_water_updraft_forced[0, 0, 0][plume]
                                * normalized_massflux_updraft_forced
                            )

            if ZERO_DIFF == 0 and total_normalized_integrated_condensate_forced[0, 0][plume] < 0.0:
                error_code[0, 0][plume] = 66

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume] + 1:
            # get back cloud_vapor_mixing_ratio_forced
            cloud_vapor_mixing_ratio_forced = (
                cloud_vapor_mixing_ratio_forced - cloud_liquid_after_rain_forced[0, 0, 0][plume]
            )


def updraft_moist_static_energy_and_momentum_budget(
    error_code: IntFieldIJ_Plume,
    start_level: IntFieldIJ,
    cloud_top_level: IntFieldIJ_Plume,
    p_forced: FloatField,
    environment_moist_static_energy: FloatField,
    environment_moist_static_energy_forced: FloatField,
    environment_moist_static_energy_cloud_levels: FloatField,
    environment_moist_static_energy_cloud_levels_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    cloud_moist_static_energy: FloatField,
    cloud_moist_static_energy_forced: FloatField,
    normalized_massflux_updraft: FloatField,
    normalized_massflux_updraft_forced: FloatField,
    mass_entrainment_updraft: FloatField,
    mass_detrainment_updraft: FloatField,
    mass_entrainment_u_updraft: FloatField,
    mass_detrainment_u_updraft: FloatField,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    u: FloatField,
    v: FloatField,
    u_c: FloatField,
    v_c: FloatField,
    u_cloud_levels: FloatField,
    v_cloud_levels: FloatField,
    partition_liquid_ice: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    vapor_excess: FloatFieldIJ,
    t_excess: FloatFieldIJ,
    add_buoyancy: FloatFieldIJ,
    plume: Int,
):
    from __externals__ import USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES, PRESSURE_GRADIENT_CONSTANT

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and plume == 0 and USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES == 1:
            # only for shallow plume
            get_delmix_implementation_here = True

    with computation(FORWARD), interval(1, None):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:
                denom = (
                    normalized_massflux_updraft[0, 0, -1]
                    - 0.5 * mass_detrainment_updraft[0, 0, -1]
                    + mass_entrainment_updraft[0, 0, -1]
                )
                denom_u = (
                    normalized_massflux_updraft[0, 0, -1]
                    - 0.5 * mass_detrainment_u_updraft[0, 0, -1]
                    + mass_entrainment_u_updraft[0, 0, -1]
                )

                if denom > 0.0 and denom_u > 0.0:
                    cloud_moist_static_energy = (
                        cloud_moist_static_energy[0, 0, -1] * normalized_massflux_updraft[0, 0, -1]
                        - (0.5 * mass_detrainment_updraft[0, 0, -1]) * cloud_moist_static_energy[0, 0, -1]
                        + mass_entrainment_updraft[0, 0, -1] * environment_moist_static_energy[0, 0, -1]
                    ) / denom

                    cloud_moist_static_energy_forced = (
                        cloud_moist_static_energy_forced[0, 0, -1]
                        * normalized_massflux_updraft_forced[0, 0, -1]
                        - 0.5
                        * mass_detrainment_updraft_forced[0, 0, -1][plume]
                        * cloud_moist_static_energy_forced[0, 0, -1]
                        + mass_entrainment_updraft_forced[0, 0, -1][plume]
                        * environment_moist_static_energy_forced[0, 0, -1]
                    ) / denom

                    if K == start_level + 1:
                        modification = (
                            cumulus_parameterization_constants.XLV * vapor_excess
                            + cumulus_parameterization_constants.CP * t_excess
                        ) + add_buoyancy
                        cloud_moist_static_energy_forced = (
                            cloud_moist_static_energy_forced
                            + modification * mass_entrainment_updraft_forced[0, 0, -1][plume] / denom
                        )
                        cloud_moist_static_energy = (
                            cloud_moist_static_energy
                            + modification * mass_entrainment_updraft[0, 0, -1] / denom
                        )

                    u_c = (
                        u_c[0, 0, -1] * normalized_massflux_updraft[0, 0, -1]
                        - 0.5 * mass_detrainment_u_updraft[0, 0, -1] * u_c[0, 0, -1]
                        + mass_entrainment_u_updraft[0, 0, -1] * u[0, 0, -1]
                        - PRESSURE_GRADIENT_CONSTANT
                        * 0.5
                        * (normalized_massflux_updraft + normalized_massflux_updraft[0, 0, -1])
                        * (u_cloud_levels - u_cloud_levels[0, 0, -1])
                    ) / denom_u

                    v_c = (
                        v_c[0, 0, -1] * normalized_massflux_updraft[0, 0, -1]
                        - 0.5 * mass_detrainment_u_updraft[0, 0, -1] * v_c[0, 0, -1]
                        + mass_entrainment_u_updraft[0, 0, -1] * v[0, 0, -1]
                        - PRESSURE_GRADIENT_CONSTANT
                        * 0.5
                        * (normalized_massflux_updraft + normalized_massflux_updraft[0, 0, -1])
                        * (v_cloud_levels - v_cloud_levels[0, 0, -1])
                    ) / denom_u

                else:
                    cloud_moist_static_energy = cloud_moist_static_energy[0, 0, -1]
                    cloud_moist_static_energy_forced = cloud_moist_static_energy_forced[0, 0, -1]
                    u_c = u_c[0, 0, -1]
                    v_c = v_c[0, 0, -1]

                cloud_moist_static_energy = (
                    cloud_moist_static_energy
                    + (1.0 - partition_liquid_ice)
                    * cloud_liquid_after_rain_forced[0, 0, 0][plume]
                    * cumulus_parameterization_constants.XLF
                )
                cloud_moist_static_energy_forced = (
                    cloud_moist_static_energy_forced
                    + (1.0 - partition_liquid_ice)
                    * cloud_liquid_after_rain_forced[0, 0, 0][plume]
                    * cumulus_parameterization_constants.XLF
                )

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0 and K >= cloud_top_level[0, 0][plume] + 2:
            cloud_moist_static_energy = environment_saturation_moist_static_energy_cloud_levels
            u_c = u_cloud_levels
            v_c = v_cloud_levels
            cloud_moist_static_energy_forced = environment_saturation_moist_static_energy_cloud_levels_forced
            normalized_massflux_updraft = 0.0
            normalized_massflux_updraft_forced = 0.0


def updraft_temperature(
    error_code: IntFieldIJ_Plume,
    updraft_t: FloatField,
    cloud_moist_static_energy_forced: FloatField,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_vapor_mixing_ratio_forced: FloatField,
    t_cloud_levels_forced: FloatField,
    plume: Int,
):
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            updraft_t = (1.0 / cumulus_parameterization_constants.CP) * (
                cloud_moist_static_energy_forced
                - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
                - cumulus_parameterization_constants.XLV * cloud_vapor_mixing_ratio_forced
            )

    with computation(PARALLEL), interval(-1, None):
        if error_code[0, 0][plume] == 0:
            updraft_t = t_cloud_levels_forced

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            updraft_t = t_cloud_levels_forced


def cup_up_aa0(
    local_buoyancy: FloatField,
    local_gamma_cloud_levels: FloatField,
    cloud_top_level: IntFieldIJ_Plume,
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
            kend = cloud_top_level[0, 0][plume] - 1

    with computation(FORWARD), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            if K >= kbeg and K <= kend:
                dz = local_geopotential_height_cloud_levels[0, 0, 1] - local_geopotential_height_cloud_levels
                aa_1 = (
                    local_normalized_massflux_updraft
                    * (constants.MAPL_GRAV / (cumulus_parameterization_constants.CP * local_t_cloud_levels))
                    * local_buoyancy
                    / (1.0 + local_gamma_cloud_levels)
                )
                aa_2 = (
                    local_normalized_massflux_updraft[0, 0, 1]
                    * (
                        constants.MAPL_GRAV
                        / (cumulus_parameterization_constants.CP * local_t_cloud_levels[0, 0, 1])
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


class UpdraftMassFluxProfile:
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
        # self._updraft_rates_pdf = stencil_factory.from_dims_halo(
        #     func=updraft_rates_pdf,
        #     compute_dims=[X_DIM, Y_DIM, Z_DIM],
        #     externals={"OVERSHOOT": cumulus_parameterization_config.OVERSHOOT},
        # )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        pass
        # self._updraft_rates_pdf(
        #     entrainment_rate=state.output.entrainment_rate,
        #     moist_static_energy=locals.environment_moist_static_energy_forced,
        #     saturation_moist_static_energy=locals.environment_saturation_mixing_ratio_cloud_levels_forced,
        #     moist_static_energy_origin_level=locals.moist_static_energy_origin_level_forced,
        #     updraft_lfc_level=state.output.updraft_lfc_level,
        #     geopotential_height=locals.geopotential_height_cloud_levels_forced,
        #     cloud_moist_static_energy=locals.cloud_moist_static_energy_forced_transported,
        #     error_code=state.output.error_code,
        #     cloud_top_level=state.output.cloud_top_level,
        #     plume=plume_dependent_constants.PLUME_INDEX,
        # )


class UpdraftInitialWorkfunctions:
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
        self._cup_up_aa0 = stencil_factory.from_dims_halo(
            func=cup_up_aa0,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._cloud_work_function_zero = stencil_factory.from_dims_halo(
            func=cloud_work_function_zero,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy,
            local_gamma_cloud_levels=locals.gamma_cloud_levels,
            cloud_top_level=state.output.cloud_top_level,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            local_geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
            local_t_cloud_levels=locals.t_cloud_levels,
            local_normalized_massflux_updraft=locals.normalized_massflux_updraft,
            local_integ=locals.integ,
            local_integ_interval=locals.integ_interval,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_0,
        )

        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy_forced,
            local_gamma_cloud_levels=locals.gamma_cloud_levels_forced,
            cloud_top_level=state.output.cloud_top_level,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            local_geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels_forced,
            local_t_cloud_levels=locals.t_cloud_levels_forced,
            local_normalized_massflux_updraft=locals.normalized_massflux_updraft_forced,
            local_integ=locals.integ,
            local_integ_interval=locals.integ_interval,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_1,
        )

        self._cloud_work_function_zero(
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_1,
        )


class UpdraftCIN:
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
        self._cup_up_aa0 = stencil_factory.from_dims_halo(
            func=cup_up_aa0,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020CumulusParameterizationState,
        locals: GF2020CumulusParameterizationLocals,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy,
            local_gamma_cloud_levels=locals.gamma_cloud_levels,
            cloud_top_level=state.output.cloud_top_level,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            local_geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels,
            local_t_cloud_levels=locals.t_cloud_levels,
            local_normalized_massflux_updraft=locals.normalized_massflux_updraft,
            local_integ=locals.integ,
            local_integ_interval=locals.integ_interval,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_0,
        )

        self._cup_up_aa0(
            local_buoyancy=locals.buoyancy_forced,
            local_gamma_cloud_levels=locals.gamma_cloud_levels_forced,
            cloud_top_level=state.output.cloud_top_level,
            updraft_lfc_level=state.output.updraft_lfc_level,
            updraft_origin_level=state.output.updraft_origin_level,
            local_geopotential_height_cloud_levels=locals.geopotential_height_cloud_levels_forced,
            local_t_cloud_levels=locals.t_cloud_levels_forced,
            local_normalized_massflux_updraft=locals.normalized_massflux_updraft_forced,
            local_integ=locals.integ,
            local_integ_interval=locals.integ_interval,
            error_code=state.output.error_code,
            plume=plume_dependent_constants.PLUME_INDEX,
            local_cloud_work_function=locals.cloud_work_function_1,
        )


class UpdraftUpdateWorkfunctions:
    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass
