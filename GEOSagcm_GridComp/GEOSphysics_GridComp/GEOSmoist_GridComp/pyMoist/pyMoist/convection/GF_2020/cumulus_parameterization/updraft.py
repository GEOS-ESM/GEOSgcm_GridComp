import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl import Local, NDSLRuntime, Quantity, QuantityFactory, StencilFactory
from ndsl.constants import I_XIM, J_DIM, K_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, Field, GlobalTable, K, computation, interval
from ndsl.dsl.typing import BoolFieldIJ, Float, FloatField, FloatFieldIJ, Int, IntFieldIJ
from ndsl.stencils.column_operations import column_max_ddim, column_min
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.config import GF2020CumulusParameterizationConfig
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Ensemble,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.plume_dependent_constants import (
    GF2020PlumeDependentConstants,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.shared_functions import get_cloud_boundary_conditions


# initialize constants and field type for UpdraftMassFlux stencil
_X_ALPHA = [
    3.699999,
    3.699999,
    3.699999,
    3.699999,
    3.024999,
    2.559999,
    2.249999,
    2.028571,
    1.862500,
    1.733333,
    1.630000,
    1.545454,
    1.475000,
    1.415385,
    1.364286,
    1.320000,
    1.281250,
    1.247059,
    1.216667,
    1.189474,
    1.165000,
    1.142857,
    1.122727,
    1.104348,
    1.087500,
    1.075000,
    1.075000,
    1.075000,
    1.075000,
    1.075000,
]

_G_ALPHA = [
    4.1706450,
    4.1706450,
    4.1706450,
    4.1706450,
    2.0469250,
    1.3878370,
    1.1330030,
    1.012418,
    0.9494680,
    0.9153771,
    0.8972442,
    0.8885444,
    0.8856795,
    0.8865333,
    0.8897996,
    0.8946404,
    0.9005030,
    0.9070138,
    0.9139161,
    0.9210315,
    0.9282347,
    0.9354376,
    0.9425780,
    0.9496124,
    0.9565111,
    0.9619183,
    0.9619183,
    0.9619183,
    0.9619183,
    0.9619183,
]
_CONSTANTS_TABLES_TYPE: Field = GlobalTable[(Float, len(_X_ALPHA))]


def updraft_mass_flux(
    error_code: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    pbl_level: IntFieldIJ,
    updraft_lfc_level: IntFieldIJ_Plume,
    lcl_level: IntFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    p_surface: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    normalized_massflux_updraft: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_updraft_modified: FloatField,
    random_number: FloatFieldIJ,
    UPDRAFT_MAX_HEIGHT_LAND: Float,
    UPDRAFT_MAX_HEIGHT_OCEAN: Float,
    plume: Int,
    X_ALPHA: _CONSTANTS_TABLES_TYPE,
    G_ALPHA: _CONSTANTS_TABLES_TYPE,
):
    """Handle mass fluxes in the updraft. This code has a number of potential paths depending on configuration
    settings - some of these paths are not yet implemented.

    This code may execute for all plumes - shallow plume execution is currently not implemented.

    Args:
        error_code (IntFieldIJ_Plume): _description_
        updraft_origin_level (IntFieldIJ_Plume): _description_
        cloud_top_level (IntFieldIJ_Plume): _description_
        pbl_level (IntFieldIJ): _description_
        updraft_lfc_level (IntFieldIJ_Plume): _description_
        lcl_level (IntFieldIJ_Plume): _description_
        p_cloud_levels_forced (FloatField_Plume): _description_
        p_surface (FloatFieldIJ): _description_
        ocean_fraction (FloatFieldIJ): _description_
        normalized_massflux_updraft (FloatField): _description_
        normalized_massflux_updraft_forced (FloatField_Plume): _description_
        normalized_massflux_updraft_modified (FloatField): _description_
        random_number (FloatFieldIJ): _description_
        UPDRAFT_MAX_HEIGHT_LAND (Float): _description_
        UPDRAFT_MAX_HEIGHT_OCEAN (Float): _description_
        plume (Int): _description_
        X_ALPHA (_CONSTANTS_TABLES_TYPE): _description_
        G_ALPHA (_CONSTANTS_TABLES_TYPE): _description_
    """
    from __externals__ import ZERO_DIFF, k_end

    with computation(FORWARD), interval(0, 1):
        # initialize constants
        PX: FloatFieldIJ = 45 / 120  # px sets the pressure level of max zu
        BETA_DEEP: FloatFieldIJ = 1.25
        G_BETA_DEEP: FloatFieldIJ = 0.8974707
        execution_choice: IntFieldIJ = -999

        # set up masks
        stop_computation: BoolFieldIJ = False
        stop_do_loop: BoolFieldIJ = False

        if ZERO_DIFF == 1:
            if plume == cumulus_parameterization_constants.DEEP and ocean_fraction > 0.90:
                # deep plume over ocean
                execution_choice = 11
            if plume == cumulus_parameterization_constants.DEEP and ocean_fraction <= 0.90:
                # deep plume over land
                execution_choice = 12
            if plume == cumulus_parameterization_constants.MID:
                # mid plume
                execution_choice = 5
        else:
            if plume == cumulus_parameterization_constants.DEEP:
                # deep plume
                execution_choice = 20
            if plume == cumulus_parameterization_constants.MID:
                # mid plume
                execution_choice = 20

    with computation(PARALLEL), interval(...):
        # fill momentum with 0
        normalized_massflux_updraft_forced[0, 0, 0][plume] = 0.0
        normalized_massflux_updraft_h = 0.0
        normalized_massflux_updraft_l = 0.0

    # EXECUTION CHOICE 20
    with computation(FORWARD), interval(0, 1):
        # land/ocean
        if error_code[0, 0][plume] == 0:
            if execution_choice == 20:

                height_updraft: FloatFieldIJ = (
                    1.0 - ocean_fraction
                ) * UPDRAFT_MAX_HEIGHT_LAND + ocean_fraction * UPDRAFT_MAX_HEIGHT_OCEAN
                # add a randomic perturbation
                height_updraft = height_updraft + random_number

                # height_updraft parameter goes from 0 to 1 = rainfall decreases with height_updraft
                p_max_normalized_massflux_updraft: FloatFieldIJ = (p_surface - 100.0) * (
                    1.0 - 0.5 * height_updraft
                ) + 0.6 * p_cloud_levels_forced.at(
                    K=cloud_top_level[0, 0][plume], ddim=[plume]
                ) * 0.5 * height_updraft

                # beta parameter: must be larger than 1
                # higher makes the profile sharper around the maximum zu
                beta: FloatFieldIJ = max(1.1, 2.1 - 0.5 * height_updraft)

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if execution_choice == 20:
                p_internal = abs(p_cloud_levels_forced[0, 0, 0][plume] - p_max_normalized_massflux_updraft)

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if execution_choice == 20:
                _, _min_loc = column_min(p_internal, 0, cloud_top_level[0, 0][plume])
                updraft_origin_level_adj: IntFieldIJ = _min_loc + 1
                updraft_origin_level_adj = max(
                    updraft_origin_level[0, 0][plume] + 1, updraft_origin_level_adj
                )
                updraft_origin_level_adj = min(updraft_origin_level_adj, cloud_top_level[0, 0][plume] + 1)

                # this alpha constrains the location of the maximun normalized_massflux_updraft_forced
                # to be at "updraft_origin_level_adj" vertical level
                alpha: FloatFieldIJ = 1.0 + (
                    (beta - 1.0)
                    * ((updraft_origin_level_adj / (cloud_top_level[0, 0][plume] + 2)))
                    / (1.0 - ((updraft_origin_level_adj) / (cloud_top_level[0, 0][plume] + 2)))
                )

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0:
            if (
                execution_choice == 20
                and K >= lcl_level[0, 0][plume] - 1
                and K <= min(k_end, cloud_top_level[0, 0][plume])
            ):
                ratio = float(K + 1) / (cloud_top_level[0, 0][plume] + 2)
                normalized_massflux_updraft_forced[0, 0, 0][plume] = (ratio ** (alpha - 1.0)) * (
                    (1.0 - ratio) ** (beta - 1.0)
                )

    with computation(BACKWARD), interval(1, -1):
        if error_code[0, 0][plume] == 0:
            if execution_choice == 20 and K <= lcl_level[0, 0][plume]:
                # special treatment below updraft_origin_level/updraft_lcl_level
                normalized_massflux_updraft_forced[0, 0, 0][plume] = (
                    normalized_massflux_updraft_forced[0, 0, 1][plume] * 0.5
                )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if execution_choice == 20:
                normalized_massflux_updraft_forced[0, 0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            _max_val, _ = column_max_ddim(
                normalized_massflux_updraft_forced,
                plume,
                0,
                min(k_end, cloud_top_level[0, 0][plume] + 1),
            )
            max_val: FloatFieldIJ = _max_val

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if max_val <= 0:
                normalized_massflux_updraft_forced[0, 0, 0][plume] = 0.0
                error_code[0, 0][plume] = 51
            else:
                if K <= min(k_end, cloud_top_level[0, 0][plume] + 1):
                    normalized_massflux_updraft_forced[0, 0, 0][plume] = normalized_massflux_updraft_forced[
                        0, 0, 0
                    ][plume] / (1.0e-9 + max_val)

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            normalized_massflux_updraft_modified = normalized_massflux_updraft_forced[0, 0, 0][plume]
            normalized_massflux_updraft = normalized_massflux_updraft_forced[0, 0, 0][plume]


def updraft_moisture(
    start_level: IntFieldIJ,
    error_code: IntFieldIJ_Plume,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_total_water_after_entrainment_forced: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    cloud_moist_static_energy_forced: FloatField,
    updraft_column_temperature_forced: FloatField,
    ocean_fraction: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    p_forced: FloatField,
    cloud_top_level: IntFieldIJ_Plume,
    d_buoyancy_forced: FloatField,
    cloud_liquid_before_rain_forced: FloatField,
    t_cloud_levels: FloatField,
    vapor_forced: FloatField,
    gamma_cloud_levels_forced: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
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
    """Compute moisture properties of the updraft and quantify and segregate precipitable water.

    Args:
        start_level (IntFieldIJ)
        error_code (IntFieldIJ_Plume)
        geopotential_height_cloud_levels_forced (FloatField)
        cloud_total_water_after_entrainment_forced (FloatField)
        cloud_liquid_after_rain_forced (FloatField_Plume)
        condensate_to_fall_forced (FloatField_Plume)
        total_normalized_integrated_condensate_forced (FloatFieldIJ_Plume)
        cloud_moist_static_energy_forced (FloatField)
        updraft_column_temperature_forced (FloatField)
        ocean_fraction (FloatFieldIJ)
        convection_fraction (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        p_forced (FloatField)
        cloud_top_level (IntFieldIJ_Plume)
        d_buoyancy_forced (FloatField)
        cloud_liquid_before_rain_forced (FloatField)
        t_cloud_levels (FloatField)
        vapor_forced (FloatField)
        gamma_cloud_levels_forced (FloatField)
        normalized_massflux_updraft_forced (FloatField_Plume)
        environment_saturation_mixing_ratio_cloud_levels_forced (FloatField)
        updraft_origin_level (IntFieldIJ_Plume)
        vapor_cloud_levels_forced (FloatField)
        vapor_excess (FloatFieldIJ)
        ccn (FloatFieldIJ)
        mass_entrainment_updraft (FloatField)
        mass_detrainment_updraft (FloatField)
        psum (FloatFieldIJ)
        psumh (FloatFieldIJ)
        c1d (FloatField)
        add_buoyancy (FloatFieldIJ)
        vertical_velocity_3d (FloatField)
        C0 (Float)
        AVERAGE_LAYER_DEPTH (Float)
        plume (Int)
    """
    from __externals__ import (
        AUTOCONV,
        BOUNDARY_CONDITION_METHOD,
        CRITICAL_MIXING_RATIO_OVER_LAND,
        CRITICAL_MIXING_RATIO_OVER_OCEAN,
        USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
        ZERO_DIFF,
        k_end,
    )

    with computation(PARALLEL), interval(...):
        # make garbage field so the get_cloud_boundary_conditions call does not break
        # this is never touched so long as compute_perturbation=False
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
        # no precip for small clouds
        total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
        psum = 0.0
        psumh = 0.0

    with computation(PARALLEL), interval(0, -1):
        condensate_to_fall_forced[0, 0, 0][plume] = 0.0
        updraft_column_temperature_forced = t_cloud_levels
        cloud_liquid_before_rain_forced = 0.0
        cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0  # liq/ice water
        cloud_total_water_after_entrainment_forced = (
            vapor_cloud_levels_forced  # total water: liq/ice = vapor water
        )

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            # get boundary condition for cloud_total_water_after_entrainment_forced
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

    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0 and K <= start_level:
            # get boundary condition for cloud_total_water_after_entrainment_forced
            cloud_total_water_after_entrainment_forced = (
                vapor_source + vapor_excess + 0.5 * add_buoyancy / cumulus_parameterization_constants.XLV
            )
            cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0

    with computation(FORWARD), interval(0, 1):
        if (
            error_code[0, 0][plume] == 0
            and USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES == 1
            and plume == cumulus_parameterization_constants.SHALLOW
        ):  # only for shallow plume
            # option to produce linear fluxes in the sub-cloud layer
            get_delmix_implementation_here = True

    with computation(PARALLEL), interval(...):
        # initialize mask to stop computation in the next block
        stop_current_index = False

    with computation(FORWARD), interval(1, None):
        if error_code[0, 0][plume] == 0:
            if K >= start_level + 1 and K <= cloud_top_level[0, 0][plume] + 1:
                dz = (
                    geopotential_height_cloud_levels_forced
                    - geopotential_height_cloud_levels_forced[0, 0, -1]
                )
                # saturation  in cloud, this is what is allowed to be in it
                saturation_cloud_liquid = (
                    environment_saturation_mixing_ratio_cloud_levels_forced
                    + (1.0 / cumulus_parameterization_constants.XLV)
                    * (gamma_cloud_levels_forced / (1.0 + gamma_cloud_levels_forced))
                    * d_buoyancy_forced
                )

                # steady state plume equation, for what could be in cloud without condensation
                denom = (
                    normalized_massflux_updraft_forced[0, 0, -1][plume]
                    - 0.5 * mass_detrainment_updraft[0, 0, -1]
                    + mass_entrainment_updraft[0, 0, -1]
                )

                if denom > 0.0:
                    cloud_total_water_after_entrainment_forced = (
                        cloud_total_water_after_entrainment_forced[0, 0, -1]
                        * normalized_massflux_updraft_forced[0, 0, -1][plume]
                        - 0.5
                        * mass_detrainment_updraft[0, 0, -1]
                        * cloud_total_water_after_entrainment_forced[0, 0, -1]
                        + mass_entrainment_updraft[0, 0, -1] * vapor_forced[0, 0, -1]
                    ) / denom

                    if K == start_level + 1:
                        cloud_total_water_after_entrainment_forced = (
                            cloud_total_water_after_entrainment_forced
                            + vapor_excess * mass_entrainment_updraft[0, 0, -1] / denom
                        )
                    # assuming no liq/ice water in the environment
                    cloud_liquid_after_rain_forced[0, 0, 0][plume] = (
                        cloud_liquid_after_rain_forced[0, 0, -1][plume]
                        * normalized_massflux_updraft_forced[0, 0, -1][plume]
                        - 0.5
                        * mass_detrainment_updraft[0, 0, -1]
                        * cloud_liquid_after_rain_forced[0, 0, -1][plume]
                    ) / denom

                else:
                    cloud_total_water_after_entrainment_forced = cloud_total_water_after_entrainment_forced[
                        0, 0, -1
                    ]
                    cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_after_rain_forced[0, 0, -1][
                        plume
                    ]

                # updraft temp
                updraft_column_temperature_forced = (1.0 / cumulus_parameterization_constants.CP) * (
                    cloud_moist_static_energy_forced
                    - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
                    - cumulus_parameterization_constants.XLV * saturation_cloud_liquid
                )

                # total condensed water before rainout
                cloud_liquid_before_rain_forced = max(
                    0.0, cloud_total_water_after_entrainment_forced - saturation_cloud_liquid
                )

                cloud_liquid_after_rain_forced[0, 0, 0][plume] = min(
                    cloud_liquid_before_rain_forced, cloud_liquid_after_rain_forced[0, 0, 0][plume]
                )

                # production term => condensation/diffusional growth
                cup = (
                    max(
                        0.0,
                        cloud_total_water_after_entrainment_forced
                        - saturation_cloud_liquid
                        - cloud_liquid_after_rain_forced[0, 0, 0][plume],
                    )
                    / dz
                )

                if C0 < 1.0e-6:
                    cloud_liquid_after_rain_forced[0, 0, 0][plume] = cloud_liquid_before_rain_forced
                    cloud_total_water_after_entrainment_forced = cloud_liquid_after_rain_forced[0, 0, 0][
                        plume
                    ] + min(cloud_total_water_after_entrainment_forced, saturation_cloud_liquid)
                    total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
                    psum = (
                        psum
                        + cloud_liquid_before_rain_forced
                        * normalized_massflux_updraft_forced[0, 0, 0][plume]
                        * dz
                    )

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
                        condensate_to_fall_forced[0, 0, 0][plume] = cx0 * max(
                            0.0, cloud_liquid_after_rain_forced[0, 0, 0][plume] - min_liq
                        )  # units kg[rain]/kg[air]
                        # convert precipitable_water_updraft_forced to
                        # normalized precipitable_water_updraft_forced
                        condensate_to_fall_forced[0, 0, 0][plume] = (
                            condensate_to_fall_forced[0, 0, 0][plume]
                            * normalized_massflux_updraft_forced[0, 0, 0][plume]
                        )

                    # total water (vapor + condensed) in updraft after the rainout
                    cloud_total_water_after_entrainment_forced = cloud_liquid_after_rain_forced[0, 0, 0][
                        plume
                    ] + min(cloud_total_water_after_entrainment_forced, saturation_cloud_liquid)

                    # integrated normalized condensates
                    total_normalized_integrated_condensate_forced[0, 0][plume] = (
                        total_normalized_integrated_condensate_forced[0, 0][plume]
                        + condensate_to_fall_forced[0, 0, 0][plume]
                    )
                    psum = (
                        psum
                        + cloud_liquid_before_rain_forced
                        * normalized_massflux_updraft_forced[0, 0, 0][plume]
                        * dz
                    )

            if ZERO_DIFF == 0 and total_normalized_integrated_condensate_forced[0, 0][plume] < 0.0:
                error_code[0, 0][plume] = 66

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume] + 1:
            # get back cloud_vapor_mixing_ratio_forced
            cloud_total_water_after_entrainment_forced = (
                cloud_total_water_after_entrainment_forced - cloud_liquid_after_rain_forced[0, 0, 0][plume]
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
    normalized_massflux_updraft_forced: FloatField_Plume,
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
    """Compute the impate of precipitation on the moist static energy, horizontal winds,
    and massflux within the updraft.

    Args:
        error_code (IntFieldIJ_Plume)
        start_level (IntFieldIJ)
        cloud_top_level (IntFieldIJ_Plume)
        p_forced (FloatField)
        environment_moist_static_energy (FloatField)
        environment_moist_static_energy_forced (FloatField)
        environment_moist_static_energy_cloud_levels (FloatField)
        environment_moist_static_energy_cloud_levels_forced (FloatField)
        environment_saturation_moist_static_energy_cloud_levels (FloatField)
        environment_saturation_moist_static_energy_cloud_levels_forced (FloatField)
        cloud_moist_static_energy (FloatField)
        cloud_moist_static_energy_forced (FloatField)
        normalized_massflux_updraft (FloatField)
        normalized_massflux_updraft_forced (FloatField_Plume)
        mass_entrainment_updraft (FloatField)
        mass_detrainment_updraft (FloatField)
        mass_entrainment_u_updraft (FloatField)
        mass_detrainment_u_updraft (FloatField)
        mass_detrainment_updraft_forced (FloatField_Plume)
        mass_entrainment_updraft_forced (FloatField_Plume)
        u (FloatField)
        v (FloatField)
        u_c (FloatField)
        v_c (FloatField)
        u_cloud_levels (FloatField)
        v_cloud_levels (FloatField)
        partition_liquid_ice (FloatField)
        cloud_liquid_after_rain_forced (FloatField_Plume)
        vapor_excess (FloatFieldIJ)
        t_excess (FloatFieldIJ)
        add_buoyancy (FloatFieldIJ)
        plume (Int)
    """
    from __externals__ import PRESSURE_GRADIENT_CONSTANT, USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES

    with computation(PARALLEL), interval(...):
        if (
            error_code[0, 0][plume] == 0
            and plume == cumulus_parameterization_constants.SHALLOW
            and USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES == 1
        ):
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
                        * normalized_massflux_updraft_forced[0, 0, -1][plume]
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
            normalized_massflux_updraft_forced[0, 0, 0][plume] = 0.0


def updraft_temperature(
    error_code: IntFieldIJ_Plume,
    updraft_column_temperature_forced: FloatField,
    cloud_moist_static_energy_forced: FloatField,
    geopotential_height_cloud_levels_forced: FloatField,
    cloud_total_water_after_entrainment_forced: FloatField,
    t_cloud_levels_forced: FloatField,
    plume: Int,
):
    """Compute the temperature within the updraft.

    Args:
        error_code (IntFieldIJ_Plume)
        updraft_column_temperature_forced (FloatField)
        cloud_moist_static_energy_forced (FloatField)
        geopotential_height_cloud_levels_forced (FloatField)
        cloud_total_water_after_entrainment_forced (FloatField)
        t_cloud_levels_forced (FloatField)
        plume (Int)
    """
    with computation(PARALLEL), interval(0, -1):
        if error_code[0, 0][plume] == 0:
            updraft_column_temperature_forced = (1.0 / cumulus_parameterization_constants.CP) * (
                cloud_moist_static_energy_forced
                - constants.MAPL_GRAV * geopotential_height_cloud_levels_forced
                - cumulus_parameterization_constants.XLV * cloud_total_water_after_entrainment_forced
            )

    with computation(PARALLEL), interval(-1, None):
        if error_code[0, 0][plume] == 0:
            updraft_column_temperature_forced = t_cloud_levels_forced

    with computation(PARALLEL), interval(...):
        if error_code[0, 0][plume] != 0:
            updraft_column_temperature_forced = t_cloud_levels_forced


def cloud_workfunction_aa0(
    error_code: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    geopotential_height: FloatField,
    normalized_massflux_updraft: FloatField,
    d_buoyancy: FloatField,
    gamma_cloud_levels: FloatField,
    t_cloud_levels: FloatField,
    workfunction: FloatFieldIJ,
    mode: Int,
    plume: Int,
):
    """Compute cloud workfunction aa0

    Stencil calculation is controlled my "mode" argument:
        mode = 0: default. uses dynamic bounds (updraft_lfc_level --> cloud_top_level)
        mode = 1: boundary layer. uses fixed lower bound (surface --> updraft_lfc_level - 1)
        mode = 2: convective inhibition. uses alternate dynamic bounds
            (updraft_origin_level --> updraft_lfc_level - 1)

    Args:
        error_code (IntFieldIJ_Plume)
        updraft_origin_level (IntFieldIJ_Plume)
        updraft_lfc_level (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        geopotential_height (FloatField)
        normalized_massflux_updraft (FloatField)
        d_buoyancy (FloatField)
        gamma_cloud_levels (FloatField)
        t_cloud_levels (FloatField)
        workfunction (FloatFieldIJ)
        mode (Int)
        plume (Int)
    """
    with computation(FORWARD), interval(0, 1):
        # initialize workfunction to zero
        workfunction = 0.0

        # set up bounds for next block
        if mode == 0:
            lower_bound: IntFieldIJ = updraft_lfc_level[0, 0][plume]
            upper_bound: IntFieldIJ = cloud_top_level[0, 0][plume]
        elif mode == 1:
            lower_bound: IntFieldIJ = 0
            upper_bound: IntFieldIJ = updraft_lfc_level[0, 0][plume] - 1
        elif mode == 2:
            lower_bound: IntFieldIJ = updraft_origin_level[0, 0][plume]
            upper_bound: IntFieldIJ = updraft_lfc_level[0, 0][plume] - 1

        upper_bound = upper_bound + 1

    with computation(FORWARD), interval(lower_bound, upper_bound):
        if error_code[0, 0][plume] == 0:
            dz = geopotential_height[0, 0, 1] - geopotential_height
            workfunction_current_level = (
                normalized_massflux_updraft
                * (constants.MAPL_GRAV / (cumulus_parameterization_constants.CP * t_cloud_levels))
                * d_buoyancy
                / (1.0 + gamma_cloud_levels)
            )
            workfunction_level_above = (
                normalized_massflux_updraft[0, 0, 1]
                * (constants.MAPL_GRAV / (cumulus_parameterization_constants.CP * t_cloud_levels[0, 0, 1]))
                * d_buoyancy[0, 0, 1]
                / (1.0 + gamma_cloud_levels[0, 0, 1])
            )
            d_workfunction = 0.5 * (workfunction_current_level + workfunction_level_above) * dz

            workfunction = workfunction + d_workfunction


def check_cloud_workfunction_1(
    error_code: IntFieldIJ_Plume,
    cloud_workfunction_1: FloatFieldIJ,
    plume: Int,
):
    """Check the validity of cloud_workfunction_1.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_workfunction_1 (FloatFieldIJ)
        plume (Int)
    """
    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if cloud_workfunction_1 == 0.0:
                error_code[0, 0][plume] = 17


def compute_precipitation_ensemble(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    precipitation_ensemble: FloatFieldIJ_Ensemble,
    plume: Int,
):
    from __externals__ import C0_MID, ENSEMBLE_MEMBERS

    with computation(FORWARD), interval(...):
        if error_code[0, 0][plume] == 0 and K <= cloud_top_level[0, 0][plume]:
            member = 0
            while member < ENSEMBLE_MEMBERS:
                precipitation_ensemble[0, 0][member] = (
                    precipitation_ensemble[0, 0][member]
                    + condensate_to_fall_forced[0, 0, 0][plume]
                    + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                member += 1

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            if precipitation_ensemble[0, 0][6] < 1.0e-6 and C0_MID > 0.0 and plume != 0:
                error_code[0, 0][plume] = 18
                member = 0
                while member < ENSEMBLE_MEMBERS:
                    precipitation_ensemble[0, 0][member] = 0.0
                member += 1

            member = 0
            while member < ENSEMBLE_MEMBERS:
                if precipitation_ensemble[0, 0][plume] < 1.0e-5:
                    precipitation_ensemble[0, 0][plume] = 0.0
                member += 1


class UpdraftMassFlux(NDSLRuntime):
    """Handle mass fluxes in the updraft. This code has a number of potential paths depending on configuration
    settings - some of these paths are not yet implemented.

    This code may execute for all plumes - shallow plume execution is currently not implemented.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # add dimension to quantityfactory and create classes for constants
        quantity_factory.add_data_dimensions({"UpdraftMassFlux_constants": len(_X_ALPHA)})

        self._X_ALPHA: Local = quantity_factory.zeros(["UpdraftMassFlux_constants"], "n/a")
        self._G_ALPHA: Local = quantity_factory.zeros(["UpdraftMassFlux_constants"], "n/a")

        self._X_ALPHA.field[:] = _X_ALPHA
        self._G_ALPHA.field[:] = _G_ALPHA

        # construct stencil
        self._updraft_mass_flux = stencil_factory.from_dims_halo(
            func=updraft_mass_flux,
            compute_dims=[I_XIM, J_DIM, K_DIM],
            externals={
                "ZERO_DIFF": cumulus_parameterization_config.ZERO_DIFF,
                "BETA_SHALLOW": cumulus_parameterization_config.BETA_SHALLOW,
                "USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES": cumulus_parameterization_config.USE_LINEAR_SUBCLOUD_MOISTURE_FLUXES,
            },
        )

    def __call__(
        self,
        error_code: Quantity,
        updraft_origin_level: Quantity,
        cloud_top_level: Quantity,
        pbl_level: Quantity,
        updraft_lfc_level: Quantity,
        lcl_level: Quantity,
        p_cloud_levels_forced: Quantity,
        p_surface: Quantity,
        ocean_fraction: Quantity,
        normalized_massflux_updraft: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        normalized_massflux_updraft_modified: Quantity,
        random_number: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._updraft_mass_flux(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            cloud_top_level=cloud_top_level,
            pbl_level=pbl_level,
            updraft_lfc_level=updraft_lfc_level,
            lcl_level=lcl_level,
            p_cloud_levels_forced=p_cloud_levels_forced,
            p_surface=p_surface,
            ocean_fraction=ocean_fraction,
            normalized_massflux_updraft=normalized_massflux_updraft,
            normalized_massflux_updraft_forced=normalized_massflux_updraft_forced,
            normalized_massflux_updraft_modified=normalized_massflux_updraft_modified,
            random_number=random_number,
            UPDRAFT_MAX_HEIGHT_LAND=plume_dependent_constants.UPDRAFT_MAX_HEIGHT_LAND,
            UPDRAFT_MAX_HEIGHT_OCEAN=plume_dependent_constants.UPDRAFT_MAX_HEIGHT_OCEAN,
            plume=plume_dependent_constants.PLUME_INDEX,
            X_ALPHA=self._X_ALPHA,
            G_ALPHA=self._G_ALPHA,
        )


class UpdraftInitialWorkfunctions(NDSLRuntime):
    """Compute initial estimates for cloud_workfunction_0 and cloud_workfunction_1."""

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._cloud_workfunction_aa0 = stencil_factory.from_dims_halo(
            func=cloud_workfunction_aa0,
            compute_dims=[I_XIM, J_DIM, K_DIM],
        )

        self._check_cloud_workfunction_1 = stencil_factory.from_dims_halo(
            func=check_cloud_workfunction_1,
            compute_dims=[I_XIM, J_DIM, K_DIM],
        )

    def __call__(
        self,
        error_code: Quantity,
        updraft_origin_level: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        geopotential_height_cloud_levels: Quantity,
        geopotential_height_cloud_levels_forced: Quantity,
        normalized_massflux_updraft: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        d_buoyancy: Quantity,
        d_buoyancy_forced: Quantity,
        gamma_cloud_levels: Quantity,
        gamma_cloud_levels_forced: Quantity,
        t_cloud_levels: Quantity,
        t_cloud_levels_forced: Quantity,
        cloud_workfunction_0: Quantity,
        cloud_workfunction_1: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._cloud_workfunction_aa0(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            geopotential_height=geopotential_height_cloud_levels,
            normalized_massflux_updraft=normalized_massflux_updraft,
            d_buoyancy=d_buoyancy,
            gamma_cloud_levels=gamma_cloud_levels,
            t_cloud_levels=t_cloud_levels,
            workfunction=cloud_workfunction_0,
            mode=Int(0),
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._cloud_workfunction_aa0(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            geopotential_height=geopotential_height_cloud_levels_forced,
            normalized_massflux_updraft=normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            d_buoyancy=d_buoyancy_forced,
            gamma_cloud_levels=gamma_cloud_levels_forced,
            t_cloud_levels=t_cloud_levels_forced,
            workfunction=cloud_workfunction_1,
            mode=Int(0),
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._check_cloud_workfunction_1(
            error_code=error_code,
            cloud_workfunction_1=cloud_workfunction_1,
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class UpdraftCIN(NDSLRuntime):
    """Compute initial estimate of updraft CIN"""

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._cloud_workfunction_aa0 = stencil_factory.from_dims_halo(
            func=cloud_workfunction_aa0,
            compute_dims=[I_XIM, J_DIM, K_DIM],
        )

    def __call__(
        self,
        error_code: Quantity,
        updraft_origin_level: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        geopotential_height_cloud_levels: Quantity,
        geopotential_height_cloud_levels_forced: Quantity,
        normalized_massflux_updraft: Quantity,
        normalized_massflux_updraft_forced: Quantity,
        d_buoyancy: Quantity,
        d_buoyancy_forced: Quantity,
        gamma_cloud_levels: Quantity,
        gamma_cloud_levels_forced: Quantity,
        t_cloud_levels: Quantity,
        t_cloud_levels_forced: Quantity,
        cin_0: Quantity,
        cin_1: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._cloud_workfunction_aa0(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            geopotential_height=geopotential_height_cloud_levels,
            normalized_massflux_updraft=normalized_massflux_updraft,
            d_buoyancy=d_buoyancy,
            gamma_cloud_levels=gamma_cloud_levels,
            t_cloud_levels=t_cloud_levels,
            workfunction=cin_0,
            mode=Int(2),
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._cloud_workfunction_aa0(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            geopotential_height=geopotential_height_cloud_levels_forced,
            normalized_massflux_updraft=normalized_massflux_updraft_forced.field[
                :, :, :, plume_dependent_constants.PLUME_INDEX
            ],
            d_buoyancy=d_buoyancy_forced,
            gamma_cloud_levels=gamma_cloud_levels_forced,
            t_cloud_levels=t_cloud_levels_forced,
            workfunction=cin_1,
            mode=Int(2),
            plume=plume_dependent_constants.PLUME_INDEX,
        )


class UpdateWorkfunctionAndPrecipitationEnsemble(NDSLRuntime):
    """Update cloud workfunctions and precipitation ensemble with the result of convection."""

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        cumulus_parameterization_config: GF2020CumulusParameterizationConfig,
    ):
        # init NDSLRuntime
        super().__init__(stencil_factory)

        # make configuration visible at runtime
        self.config = config
        self.cumulus_parameterization_config = cumulus_parameterization_config

        # construct stencils and functions
        self._cloud_workfunction_aa0 = stencil_factory.from_dims_halo(
            func=cloud_workfunction_aa0,
            compute_dims=[I_XIM, J_DIM, K_DIM],
        )

        self._compute_precipitation_ensemble = stencil_factory.from_dims_halo(
            func=compute_precipitation_ensemble,
            compute_dims=[I_XIM, J_DIM, K_DIM],
            externals={
                "ENSEMBLE_MEMBERS": cumulus_parameterization_constants.MAXENS1
                * cumulus_parameterization_constants.MAXENS2
                * cumulus_parameterization_constants.MAXENS3,
                "C0_MID": cumulus_parameterization_config.C0_MID,
            },
        )

    def __call__(
        self,
        error_code: Quantity,
        updraft_origin_level: Quantity,
        updraft_lfc_level: Quantity,
        cloud_top_level: Quantity,
        geopotential_height_cloud_levels_modified: Quantity,
        normalized_massflux_updraft_modified: Quantity,
        d_buoyancy_modified: Quantity,
        gamma_cloud_levels: Quantity,
        t_cloud_levels_modified: Quantity,
        cloud_workfunction_0_modified: Quantity,
        condensate_to_fall_forced: Quantity,
        evaporate_in_downdraft_forced: Quantity,
        epsilon_forced: Quantity,
        precipitation_ensemble: Quantity,
        plume_dependent_constants: GF2020PlumeDependentConstants,
    ):
        self._cloud_workfunction_aa0(
            error_code=error_code,
            updraft_origin_level=updraft_origin_level,
            updraft_lfc_level=updraft_lfc_level,
            cloud_top_level=cloud_top_level,
            geopotential_height=geopotential_height_cloud_levels_modified,
            normalized_massflux_updraft=normalized_massflux_updraft_modified,
            d_buoyancy=d_buoyancy_modified,
            gamma_cloud_levels=gamma_cloud_levels,
            t_cloud_levels=t_cloud_levels_modified,
            workfunction=cloud_workfunction_0_modified,
            mode=Int(0),
            plume=plume_dependent_constants.PLUME_INDEX,
        )

        self._compute_precipitation_ensemble(
            error_code=error_code,
            cloud_top_level=cloud_top_level,
            condensate_to_fall_forced=condensate_to_fall_forced,
            evaporate_in_downdraft_forced=evaporate_in_downdraft_forced,
            epsilon_forced=epsilon_forced,
            precipitation_ensemble=precipitation_ensemble,
            plume=plume_dependent_constants.PLUME_INDEX,
        )
