from ndsl import NDSLRuntime, QuantityFactory, StencilFactory
from ndsl.constants import I_DIM, J_DIM, K_DIM, K_INTERFACE_DIM
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, K, abs, computation, floor, interval, max, min, sqrt
from ndsl.dsl.typing import Float, FloatField, FloatFieldIJ, IntFieldIJ

import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_ConvectionTracers,
    FloatField_ConvectionTracers_Plume,
    FloatField_Plume,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.convection.GF_2020.cumulus_parameterization.state import GF2020CumulusParameterizationState
from pyMoist.convection.GF_2020.locals import GF2020Locals
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection_tracers import ConvectionTracers
from pyMoist.saturation_tables.saturation_specific_humidity_functions import saturation_specific_humidity
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.saturation_tables.types import GlobalTable_saturation_tables


def compute_extra_inputs_from_state(
    p_interface: FloatField,
    p: FloatField,
    p_kappa: FloatField,
    edge_height_above_surface: FloatField,
    layer_height_above_surface: FloatField,
    geopotential_height_interface: FloatField,
    t: FloatField,
    th: FloatField,
    vapor: FloatField,
    mass: FloatField,
    w: FloatField,
    omega: FloatField,
    vertical_motion: FloatField,
    tpwi: FloatFieldIJ,
    tpwi_star: FloatFieldIJ,
    seed_convection: FloatFieldIJ,
    area: FloatFieldIJ,
    modified_area: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    ese: GlobalTable_saturation_tables,
    esx: GlobalTable_saturation_tables,
):
    """
    Performs initial setup for the GF 2020 convection scheme:
     - Compute derived states
     - initialize stochastic variability for convection
     - Modify area (m^2) here so GF scale dependence has a convection_fraction dependence

    This stencil MUST be built using K_INTERFACE_DIM to function properly.

    Args:
        p_interface (FloatField)
        p (FloatField)
        p_kappa (FloatField)
        edge_height_above_surface (FloatField)
        layer_height_above_surface (FloatField)
        geopotential_height_interface (FloatField)
        t (FloatField)
        th (FloatField)
        vapor (FloatField)
        mass (FloatField)
        w (FloatField)
        omega (FloatField)
        vertical_motion (FloatField)
        tpwi (FloatFieldIJ)
        tpwi_star (FloatFieldIJ)
        seed_convection (FloatFieldIJ)
        area (FloatFieldIJ)
        modified_area (FloatFieldIJ)
        convection_fraction (FloatFieldIJ)
        ese (GlobalTable_saturation_tables)
        esx (GlobalTable_saturation_tables)
    """
    from __externals__ import GF_MIN_AREA, LHYDROSTATIC, STOCH_BOT, STOCH_TOP, STOCHASTIC_CONVECTION, k_end

    # compute derived states
    with computation(PARALLEL), interval(...):
        edge_height_above_surface = geopotential_height_interface - geopotential_height_interface.at(K=k_end)

    with computation(PARALLEL), interval(0, -1):
        p = 0.5 * (p_interface + p_interface[0, 0, 1])
        p_kappa = (p / constants.MAPL_P00) ** (constants.MAPL_KAPPA)
        layer_height_above_surface = 0.5 * (edge_height_above_surface + edge_height_above_surface[0, 0, 1])
        th = t / p_kappa
        mass = (p_interface[0, 0, 1] - p_interface) / constants.MAPL_GRAV

        if LHYDROSTATIC:
            vertical_motion = -1 * omega / (constants.MAPL_GRAV * p / (constants.MAPL_RDRY * t * (1.0 + constants.MAPL_VIREPS * vapor)))
        else:
            vertical_motion = w

    with computation(FORWARD), interval(0, 1):
        tpwi = vapor * mass
        qsat, _ = saturation_specific_humidity(t, p, ese, esx)
        tpwi_star = qsat * mass

    with computation(FORWARD), interval(1, -1):
        tpwi = tpwi + vapor * mass
        qsat, _ = saturation_specific_humidity(t, p, ese, esx)
        tpwi_star = tpwi_star + qsat * mass

    with computation(FORWARD), interval(0, 1):
        # initialize stochastic variability for convection
        if STOCHASTIC_CONVECTION:
            # Create bit-processor-reproducible random white noise for convection [0:1]
            seedini = 1000000 * (100 * t.at(K=k_end) - floor(100 * t.at(K=k_end)))
            seed_convection = sqrt(max(min(seedini / 1000000.0, 1.0), 0.0))
            # Create stochastic variability to GF sigma
            seed_convection = sqrt(1.0 - (1.0 - seed_convection)) * (STOCH_TOP - STOCH_BOT) + STOCH_BOT
        else:
            seed_convection = 1.0

        # Modify area (m^2) here so GF scale dependence has a convection_fraction dependence
        if GF_MIN_AREA > 0:
            if area > GF_MIN_AREA:
                modified_area = GF_MIN_AREA * convection_fraction + area * (1.0 - convection_fraction)
            else:
                modified_area = area
        elif GF_MIN_AREA < 0:
            if area > abs(GF_MIN_AREA):
                modified_area = area * convection_fraction + abs(GF_MIN_AREA) * (1.0 - convection_fraction)
            else:
                modified_area = area
        else:
            modified_area = area


def pass_back_to_model_state(local_seed_convection: FloatFieldIJ, model_state_seed_convection: FloatFieldIJ):
    with computation(FORWARD), interval(0, 1):
        model_state_seed_convection = local_seed_convection


def zero_state(
    dvapordt_deep_convection: FloatField,
    dtdt_deep_convection: FloatField,
    dudt_deep_convection: FloatField,
    dvdt_deep_convection: FloatField,
    sigma_deep: FloatFieldIJ,
    sigma_mid: FloatFieldIJ,
    mass_flux_shalow: FloatField,
    mass_flux_mid: FloatField,
    mass_flux_deep_updraft: FloatField,
    mass_flux_deep_updraft_interface: FloatField,
    mass_flux_deep_updraft_detrained: FloatField,
    mass_flux_deep_downdraft: FloatField,
    mass_flux_cloud_base: FloatField,
    mass_flux_cloud_base_shallow: FloatFieldIJ,
    mass_flux_cloud_base_mid: FloatFieldIJ,
    mass_flux_cloud_base_deep: FloatFieldIJ,
    convection_code_shallow: FloatFieldIJ,
    convection_code_mid: FloatFieldIJ,
    convection_code_deep: FloatFieldIJ,
    cloud_workfunction_0: FloatFieldIJ,
    cloud_workfunction_1: FloatFieldIJ,
    cloud_workfunction_2: FloatFieldIJ,
    cloud_workfunction_3: FloatFieldIJ,
    cloud_workfunction_1_pbl: FloatFieldIJ,
    cloud_workfunction_1_cin: FloatFieldIJ,
    convective_precipitation_RAS: FloatField,
    convective_precipitation_GF: FloatFieldIJ,
    convective_condensate_source: FloatField,
    convective_condensate_grid_mean: FloatField,
    total_water_flux_deep_convection_interface: FloatField,
    updraft_area_fraction: FloatField,
    updraft_vertical_velocity: FloatField,
    entrainment_parameter: FloatField,
    lightning_density: FloatFieldIJ,
    pbl_time_scale: FloatFieldIJ,
    cape_removal_time_scale: FloatFieldIJ,
):
    """
    Zero fields from the GEOS model state.

    All of these fields are outputs of GF2020 which may be set in the CumulusParameterization core routine,
    so we are really just clearing any data from the previous timestep to ensure nothing leaks through.

    Must be built with K_INTERFACE_DIM.

    Args:
        dvapordt_deep_convection (FloatField)
        dtdt_deep_convection (FloatField)
        dudt_deep_convection (FloatField)
        dvdt_deep_convection (FloatField)
        sigma_deep (FloatFieldIJ)
        sigma_mid (FloatFieldIJ)
        mass_flux_shalow (FloatField)
        mass_flux_mid (FloatField)
        mass_flux_deep_updraft (FloatField)
        mass_flux_deep_updraft_interface (FloatField)
        mass_flux_deep_updraft_detrained (FloatField)
        mass_flux_deep_downdraft (FloatField)
        mass_flux_cloud_base (FloatField)
        mass_flux_cloud_base_shallow (FloatFieldIJ)
        mass_flux_cloud_base_mid (FloatFieldIJ)
        mass_flux_cloud_base_deep (FloatFieldIJ)
        convection_code_shallow (FloatFieldIJ)
        convection_code_mid (FloatFieldIJ)
        convection_code_deep (FloatFieldIJ)
        cloud_work_function_0 (FloatFieldIJ)
        cloud_work_function_1 (FloatFieldIJ)
        cloud_work_function_2 (FloatFieldIJ)
        cloud_work_function_3 (FloatFieldIJ)
        cloud_work_function_1_pbl (FloatFieldIJ)
        cloud_work_function_1_cin (FloatFieldIJ)
        convective_precipitation_RAS (FloatField)
        convective_precipitation_GF (FloatFieldIJ)
        convective_condensate_source (FloatField)
        convective_condensate_grid_mean (FloatField)
        total_water_flux_deep_convection_interface (FloatField)
        updraft_area_fraction (FloatField)
        updraft_vertical_velocity (FloatField)
        entrainment_parameter (FloatField)
        lightning_density (FloatFieldIJ)
        pbl_time_scale (FloatFieldIJ)
        cape_removal_time_scale (FloatFieldIJ)
    """
    with computation(PARALLEL), interval(0, -1):
        dvapordt_deep_convection = 0.0
        dtdt_deep_convection = 0.0
        dudt_deep_convection = 0.0
        dvdt_deep_convection = 0.0

    with computation(FORWARD), interval(0, 1):
        sigma_deep = 0.0
        sigma_mid = 0.0

    with computation(PARALLEL), interval(0, -1):
        mass_flux_shalow = 0.0
        mass_flux_mid = 0.0
        mass_flux_deep_updraft = 0.0

    with computation(PARALLEL), interval(...):
        mass_flux_deep_updraft_interface = 0.0

    with computation(PARALLEL), interval(0, -1):
        mass_flux_deep_updraft_detrained = 0.0
        mass_flux_deep_downdraft = 0.0
        mass_flux_cloud_base = 0.0

    with computation(FORWARD), interval(0, 1):
        mass_flux_cloud_base_shallow = 0.0
        mass_flux_cloud_base_mid = 0.0
        mass_flux_cloud_base_deep = 0.0

        convection_code_shallow = 0.0
        convection_code_mid = 0.0
        convection_code_deep = 0.0

        cloud_workfunction_0 = 0.0
        cloud_workfunction_1 = 0.0
        cloud_workfunction_2 = 0.0
        cloud_workfunction_3 = 0.0
        cloud_workfunction_1_pbl = 0.0
        cloud_workfunction_1_cin = 0.0

    with computation(PARALLEL), interval(0, -1):
        convective_precipitation_RAS = 0.0

    with computation(FORWARD), interval(0, 1):
        convective_precipitation_GF = 0.0

    with computation(PARALLEL), interval(0, -1):
        convective_condensate_source = 0.0
        convective_condensate_grid_mean = 0.0

    with computation(PARALLEL), interval(...):
        total_water_flux_deep_convection_interface = 0.0

    with computation(PARALLEL), interval(0, -1):
        updraft_area_fraction = 0.0
        updraft_vertical_velocity = 0.0
        entrainment_parameter = 0.0

    with computation(FORWARD), interval(0, 1):
        lightning_density = 0.0
        pbl_time_scale = 0.0
        cape_removal_time_scale = 0.0


def prefill_cumulus_parameterization_state(
    error_code: IntFieldIJ_Plume,
    downdraft_origin_level: IntFieldIJ_Plume,
    lcl_level: IntFieldIJ_Plume,
    updraft_origin_level: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    kstabi: IntFieldIJ_Plume,
    kstabm: IntFieldIJ_Plume,
    precip: FloatFieldIJ_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    total_normalized_integrated_condensate_forced: FloatFieldIJ_Plume,
    scale_dependence_factor: FloatFieldIJ_Plume,
    p_cloud_levels_forced: FloatField_Plume,
    entrainment_rate: FloatField_Plume,
    mass_entrainment_updraft_forced: FloatField_Plume,
    mass_entrainment_downdraft_forced: FloatField_Plume,
    mass_detrainment_updraft_forced: FloatField_Plume,
    mass_detrainment_downdraft_forced: FloatField_Plume,
    normalized_massflux_updraft_forced: FloatField_Plume,
    normalized_massflux_downdraft_forced: FloatField_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    t_updraft: FloatField_Plume,
    convective_cloud_fraction_output: FloatField_Plume,
    dtdt: FloatField_Plume,
    dudt: FloatField_Plume,
    dvdt: FloatField_Plume,
    dvapordt: FloatField_Plume,
    dcloudicedt: FloatField_Plume,
    dnicedt: FloatField_Plume,
    dnliquiddt: FloatField_Plume,
    dbuoyancydt: FloatField_Plume,
    chemistry_tracers_output: FloatField_ConvectionTracers_Plume,
    evaporation_sublimation_tendency: FloatField,
    convective_precip_flux: FloatField,
    t_perturbation: FloatField,
    omega: FloatField,
    large_scale_ice: FloatField,
    convective_ice: FloatField,
    large_scale_liquid: FloatField,
    convective_liquid: FloatField,
    large_scale_cloud_fraction: FloatField,
    convective_cloud_fraction: FloatField,
    lightning_density: FloatFieldIJ,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    last_error_code: IntFieldIJ,
):
    """
    Zero fields from the CumulusParameterization state.

    All of these fields may be set in the CumulusParameterization core routine,
    so they need to be reset to ensure no lingering data from the previous
    timestep makes it through.

    Args:
        error_code (IntFieldIJ_Plume)
        downdraft_origin_level (IntFieldIJ_Plume)
        lcl_level (IntFieldIJ_Plume)
        updraft_origin_level (IntFieldIJ_Plume)
        updraft_lfc_level (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        kstabi (IntFieldIJ_Plume)
        kstabm (IntFieldIJ_Plume)
        precip (FloatFieldIJ_Plume)
        cloud_base_mass_flux_modified (FloatFieldIJ_Plume)
        epsilon_forced (FloatFieldIJ_Plume)
        total_normalized_integrated_condensate_forced (FloatFieldIJ_Plume)
        scale_dependence_factor (FloatFieldIJ_Plume)
        p_cloud_levels_forced (FloatField_Plume)
        entrainment_rate (FloatField_Plume)
        mass_entrainment_updraft_forced (FloatField_Plume)
        mass_entrainment_downdraft_forced (FloatField_Plume)
        mass_detrainment_updraft_forced (FloatField_Plume)
        mass_detrainment_downdraft_forced (FloatField_Plume)
        normalized_massflux_updraft_forced (FloatField_Plume)
        normalized_massflux_downdraft_forced (FloatField_Plume)
        condensate_to_fall_forced (FloatField_Plume)
        evaporate_in_downdraft_forced (FloatField_Plume)
        cloud_liquid_after_rain_forced (FloatField_Plume)
        t_updraft (FloatField_Plume)
        convective_cloud_fraction_output (FloatField_Plume)
        dtdt (FloatField_Plume)
        dudt (FloatField_Plume)
        dvdt (FloatField_Plume)
        dvapordt (FloatField_Plume)
        dcloudicedt (FloatField_Plume)
        dnicedt (FloatField_Plume)
        dnliquiddt (FloatField_Plume)
        dbuoyancydt (FloatField_Plume)
        chemistry_tracers_output (FloatField_ConvectionTracers_Plume)
        evaporation_sublimation_tendency (FloatField)
        convective_precip_flux (FloatField)
        t_perturbation (FloatField)
        omega (FloatField)
        large_scale_ice (FloatField)
        convective_ice (FloatField)
        large_scale_liquid (FloatField)
        convective_liquid (FloatField)
        large_scale_cloud_fraction (FloatField)
        convective_cloud_fraction (FloatField)
        lightning_density (FloatFieldIJ)
        t_excess (FloatFieldIJ)
        vapor_excess (FloatFieldIJ)
        last_error_code (IntFieldIJ)
    """
    from __externals__ import APPLY_SUBSIDENCE_MICROPHYSICS, NUMBER_OF_PLUMES

    with computation(FORWARD), interval(0, 1):
        plume = 0
        while plume < NUMBER_OF_PLUMES:
            error_code[0, 0][plume] = 0
            downdraft_origin_level[0, 0][plume] = 0
            lcl_level[0, 0][plume] = 0
            updraft_origin_level[0, 0][plume] = 0
            updraft_lfc_level[0, 0][plume] = 0
            cloud_top_level[0, 0][plume] = 0
            kstabi[0, 0][plume] = 0
            kstabm[0, 0][plume] = 0
            precip[0, 0][plume] = 0.0
            cloud_base_mass_flux_modified[0, 0][plume] = 0.0
            epsilon_forced[0, 0][plume] = 0.0
            total_normalized_integrated_condensate_forced[0, 0][plume] = 0.0
            scale_dependence_factor[0, 0][plume] = 0.0
            plume += 1

    with computation(PARALLEL), interval(...):
        plume = 0
        while plume < NUMBER_OF_PLUMES:
            p_cloud_levels_forced[0, 0, 0][plume] = 0.0
            entrainment_rate[0, 0, 0][plume] = 0.0
            mass_entrainment_updraft_forced[0, 0, 0][plume] = 0.0
            mass_entrainment_downdraft_forced[0, 0, 0][plume] = 0.0
            mass_detrainment_updraft_forced[0, 0, 0][plume] = 0.0
            mass_detrainment_downdraft_forced[0, 0, 0][plume] = 0.0
            normalized_massflux_updraft_forced[0, 0, 0][plume] = 0.0
            normalized_massflux_downdraft_forced[0, 0, 0][plume] = 0.0
            condensate_to_fall_forced[0, 0, 0][plume] = 0.0
            evaporate_in_downdraft_forced[0, 0, 0][plume] = 0.0
            cloud_liquid_after_rain_forced[0, 0, 0][plume] = 0.0
            t_updraft[0, 0, 0][plume] = 0.0
            convective_cloud_fraction_output[0, 0, 0][plume] = 0.0

            dtdt[0, 0, 0][plume] = 0.0
            dudt[0, 0, 0][plume] = 0.0
            dvdt[0, 0, 0][plume] = 0.0
            dvapordt[0, 0, 0][plume] = 0.0
            dcloudicedt[0, 0, 0][plume] = 0.0
            dnicedt[0, 0, 0][plume] = 0.0
            dnliquiddt[0, 0, 0][plume] = 0.0
            dbuoyancydt[0, 0, 0][plume] = 0.0

            tracer = 0
            while tracer < constants.NUMBER_OF_TRACERS:
                chemistry_tracers_output[0, 0, 0][plume, tracer] = 0.0
                tracer += 1

            plume += 1

        evaporation_sublimation_tendency = 0.0
        convective_precip_flux = 0.0
        t_perturbation = 0.0
        omega = 0.0

        if APPLY_SUBSIDENCE_MICROPHYSICS == 1:
            large_scale_ice = 0.0
            convective_ice = 0.0
            large_scale_liquid = 0.0
            convective_liquid = 0.0
            large_scale_cloud_fraction = 0.0
            convective_cloud_fraction = 0.0

    with computation(FORWARD), interval(0, 1):
        lightning_density = 0.0
        t_excess = 0.0
        vapor_excess = 0.0
        last_error_code = -999


def prefill_locals(
    rtgt: FloatFieldIJ,
    precip: FloatFieldIJ,
    t_tendency_from_vapor: FloatField,
    dtdt: FloatField,
    dvapordt: FloatField,
    dcloudicedt: FloatField,
    dudt: FloatField,
    dvdt: FloatField,
    dlarge_scale_icedt: FloatField,
    dconvective_icedt: FloatField,
    dlarge_scale_liquiddt: FloatField,
    dconvective_liquiddt: FloatField,
    dlarge_scale_cloud_fractiondt: FloatField,
    dconvective_cloud_fractiondt: FloatField,
    dbuoyancydt: FloatField,
    evaporation_sublimation_tendency: FloatField,
    convective_precip_flux: FloatField,
):
    """
    Zero local fields which are conditionally written to ensure no data remains from the previous timestep.

    Args:
        rtgt (FloatFieldIJ)
        precip (FloatFieldIJ)
        t_tendency_from_vapor (FloatField)
        dtdt (FloatField)
        dvapordt (FloatField)
        dcloudicedt (FloatField)
        dudt (FloatField)
        dvdt (FloatField)
        dlarge_scale_icedt (FloatField)
        dconvective_icedt (FloatField)
        dlarge_scale_liquiddt (FloatField)
        dconvective_liquiddt (FloatField)
        dlarge_scale_cloud_fractiondt (FloatField)
        dconvective_cloud_fractiondt (FloatField)
        dbuoyancydt (FloatField)
        evaporation_sublimation_tendency (FloatField)
        convective_precip_flux (FloatField)
    """
    with computation(FORWARD), interval(0, 1):
        rtgt = 1.0
        precip = 0.0

    with computation(PARALLEL), interval(...):
        t_tendency_from_vapor = 0.0
        dtdt = 0.0
        dvapordt = 0.0
        dcloudicedt = 0.0
        dudt = 0.0
        dvdt = 0.0
        dlarge_scale_icedt = 0.0
        dconvective_icedt = 0.0
        dlarge_scale_liquiddt = 0.0
        dconvective_liquiddt = 0.0
        dlarge_scale_cloud_fractiondt = 0.0
        dconvective_cloud_fractiondt = 0.0
        dbuoyancydt = 0.0
        evaporation_sublimation_tendency = 0.0
        convective_precip_flux = 0.0


def set_2d_fields(
    aot500: FloatFieldIJ,
    t: FloatField,
    t_2m_max: Float,
    t_2m: FloatFieldIJ,
    t_2m_local: FloatFieldIJ,
    evaporation: FloatFieldIJ,
    evaporation_local: FloatFieldIJ,
    sensible_heat_flux: FloatFieldIJ,
    sensible_heat_flux_local: FloatFieldIJ,
    p_interface: FloatField,
    vapor: FloatField,
    geopotential_height_surface: FloatFieldIJ,
    topography_height: FloatFieldIJ,
    land_fraction: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    area: FloatFieldIJ,
    grid_length: FloatFieldIJ,
    pbl_level: FloatFieldIJ,
    pbl_level_flipped: IntFieldIJ,
):
    """
    Compute 2-dimensional inputs for the CumulusParameterization core routine.

    Args:
        aot500 (FloatFieldIJ)
        t (FloatField)
        t_2m_max (Float)
        t_2m (FloatFieldIJ)
        t_2m_local (FloatFieldIJ)
        evaporation (FloatFieldIJ)
        evaporation_local (FloatFieldIJ)
        sensible_heat_flux (FloatFieldIJ)
        sensible_heat_flux_local (FloatFieldIJ)
        p_interface (FloatField)
        vapor (FloatField)
        geopotential_height_surface (FloatFieldIJ)
        topography_height (FloatFieldIJ)
        land_fraction (FloatFieldIJ)
        ocean_fraction (FloatFieldIJ)
        area (FloatFieldIJ)
        grid_length (FloatFieldIJ)
        pbl_level (FloatFieldIJ)
        pbl_level_flipped (IntFieldIJ)
    """
    from __externals__ import SIZE_I_DIM, SIZE_J_DIM, k_end

    with computation(FORWARD), interval(0, 1):
        aot500 = 0.1

        if t_2m_max < 1.0e-6:
            # in case convection is called before surface data is computed
            t_2m_local = t.at(K=k_end)  # kelvin
        else:
            t_2m_local = t_2m  # kelvin

        # moisture flux from sfc
        evaporation_local = evaporation  # kg m-2 s-1

        # sensible–heat_flux comes in W m-2, below it is converted to K m s-1
        sensible_heat_flux_local = sensible_heat_flux / (1004.0 * p_interface.at(K=k_end + 1) / (287.04 * t.at(K=k_end) * (1.0 + 0.608 * vapor.at(K=k_end))))  # K m s-1

        # topography height  (m)
        topography_height = geopotential_height_surface / constants.MAPL_GRAV

        # land/ocean fraction: land if < 1 ,ocean if = 1
        ocean_fraction = 1.0 - land_fraction

        # grid length for the scale awareness (m)
        grid_length = sqrt(area)
        # special setting for SCM runs
        if SIZE_I_DIM == 1 and SIZE_J_DIM == 1:
            grid_length = 100000.0

        # flip pbl_level
        if pbl_level != -1.0:
            pbl_level_flipped = k_end - int(round(pbl_level))
        else:
            pbl_level_flipped = 0


def choose_environment_and_flip_k_axis(
    dz: FloatField,
    air_density: FloatField,
    geopotential_height_interface: FloatField,
    t: FloatField,
    t_flipped: FloatField,
    t_timestep_start: FloatField,
    p: FloatField,
    p_interface: FloatField,
    p_interface_timestep_start: FloatField,
    p_flipped: FloatField,
    p_surface_flipped: FloatFieldIJ,
    vapor: FloatField,
    vapor_timestep_start: FloatField,
    vapor_flipped: FloatField,
    vapor_current_flipped: FloatField,
    u: FloatField,
    u_timestep_start: FloatField,
    u_flipped: FloatField,
    v: FloatField,
    v_timestep_start: FloatField,
    v_flipped: FloatField,
    w: FloatField,
    w_flipped: FloatField,
    layer_height_above_surface: FloatField,
    layer_height_above_surface_flipped: FloatField,
    edge_height_above_surface: FloatField,
    edge_height_above_surface_flipped: FloatField,
    mass: FloatField,
    mass_flipped: FloatField,
    scalar_diffusivity: FloatField,
    scalar_diffusivity_flipped: FloatField,
    lateral_entrainment_rate: FloatField,
    lateral_entrainment_rate_flipped: FloatField,
    buoyancy: FloatField,
    buoyancy_excess: FloatField,
    dtdt_shortwave: FloatField,
    dtdt_longwave: FloatField,
    dtdt_from_dynamics: FloatField,
    dtdt_pbl: FloatField,
    dvapordt_from_dynamics: FloatField,
    dspecific_humiditydt_pbl: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    advective_forcing_t: FloatField,
    convective_liquid: FloatField,
    convective_liquid_flipped: FloatField,
    convective_ice: FloatField,
    convective_ice_flipped: FloatField,
    convective_cloud_fraction: FloatField,
    convective_cloud_fraction_flipped: FloatField,
    large_scale_liquid: FloatField,
    large_scale_liquid_flipped: FloatField,
    large_scale_ice: FloatField,
    large_scale_ice_flipped: FloatField,
    large_scale_cloud_fraction: FloatField,
    large_scale_cloud_fraction_flipped: FloatField,
    convection_tracer: FloatField,
    total_precipitable_water_initial: FloatFieldIJ,
    saturation_total_precipitable_water_initial: FloatFieldIJ,
    saturation_water_vapor: FloatFieldIJ,
):
    """
    Get the desired state for the convection scheme, controlled by external GF_ENV_SETTING.

    GF_ENV_SETTING = 0:
        use state updated by during current timestep (dynamics and physics)
    GF_ENV_SETTING = 1:
        use state from start of timestep (pre-dynamics)

    Args:
        dz (FloatField)
        air_density (FloatField)
        geopotential_height_interface (FloatField)
        t (FloatField)
        t_flipped (FloatField)
        t_timestep_start (FloatField)
        p (FloatField)
        p_interface (FloatField)
        p_interface_timestep_start (FloatField)
        p_flipped (FloatField)
        p_surface_flipped (FloatFieldIJ)
        vapor (FloatField)
        vapor_timestep_start (FloatField)
        vapor_flipped (FloatField)
        vapor_current_flipped (FloatField)
        u (FloatField)
        u_timestep_start (FloatField)
        u_flipped (FloatField)
        v (FloatField)
        v_timestep_start (FloatField)
        v_flipped (FloatField)
        w (FloatField)
        w_flipped (FloatField)
        layer_height_above_surface (FloatField)
        layer_height_above_surface_flipped (FloatField)
        edge_height_above_surface (FloatField)
        edge_height_above_surface_flipped (FloatField)
        mass (FloatField)
        mass_flipped (FloatField)
        scalar_diffusivity (FloatField)
        scalar_diffusivity_flipped (FloatField)
        lateral_entrainment_rate (FloatField)
        lateral_entrainment_rate_flipped (FloatField)
        buoyancy (FloatField)
        buoyancy_excess (FloatField)
        dtdt_shortwave (FloatField)
        dtdt_longwave (FloatField)
        dtdt_from_dynamics (FloatField)
        dtdt_pbl (FloatField)
        dvapordt_from_dynamics (FloatField)
        dspecific_humiditydt_pbl (FloatField)
        grid_scale_forcing_t (FloatField)
        grid_scale_forcing_vapor (FloatField)
        subgrid_scale_forcing_t (FloatField)
        subgrid_scale_forcing_vapor (FloatField)
        advective_forcing_t (FloatField)
        convective_liquid (FloatField)
        convective_liquid_flipped (FloatField)
        convective_ice (FloatField)
        convective_ice_flipped (FloatField)
        convective_cloud_fraction (FloatField)
        convective_cloud_fraction_flipped (FloatField)
        large_scale_liquid (FloatField)
        large_scale_liquid_flipped (FloatField)
        large_scale_ice (FloatField)
        large_scale_ice_flipped (FloatField)
        large_scale_cloud_fraction (FloatField)
        large_scale_cloud_fraction_flipped (FloatField)
        convection_tracer (FloatField)
        total_precipitable_water_initial (FloatFieldIJ)
        saturation_total_precipitable_water_initial (FloatFieldIJ)
        saturation_water_vapor (FloatFieldIJ)
    """
    from __externals__ import CONVECTION_TRACER, ENTRVERSION, GF_ENV_SETTING, k_end

    # 1st setting: enviromental state is the one already modified by dyn + physics
    with computation(PARALLEL), interval(...):
        if GF_ENV_SETTING == 0:
            dz = -(geopotential_height_interface[0, 0, 1] - geopotential_height_interface)
            air_density = p / (287.04 * t * (1.0 + 0.608 * vapor))
            t_flipped = t.at(K=k_end - K)
            p_flipped = p.at(K=k_end - K)
            vapor_flipped = vapor.at(K=k_end - K)
            vapor_current_flipped = vapor.at(K=k_end - K)
            u_local = u.at(K=k_end - K)
            v_local = v.at(K=k_end - K)
            w_flipped = w.at(K=k_end - K)
            layer_height_above_surface_flipped = layer_height_above_surface.at(K=k_end - K)
            edge_height_above_surface_flipped = edge_height_above_surface.at(K=k_end - K)
            mass_flipped = mass.at(K=k_end - K)
            scalar_diffusivity_flipped = scalar_diffusivity.at(K=k_end - K)

            # Grid and sub-grid scale forcings for convection
            grid_scale_forcing_t = 0.0
            grid_scale_forcing_vapor = 0.0
            subgrid_scale_forcing_t = 0.0
            subgrid_scale_forcing_vapor = 0.0
            advective_forcing_t = 0.0

    with computation(FORWARD), interval(0, 1):
        if GF_ENV_SETTING == 0:
            p_surface_flipped = p_interface.at(K=k_end + 1)

    # 2nd setting: environmental state is that one before any tendency
    # is applied (i.e, at begin of each time step).
    # Get back the model state, heights and others variables at time N
    # (or at the beggining of current time step)
    # In physics, the state vars (t,u,v,PLE) are untouched and represent the
    # model state after dynamics phase 1. But, "Q" is modified by physics, so
    # depending on what was called before this subroutine, "Q" may be already
    # changed from what it was just after dynamics phase 1. To solve this issue,
    # "Q" just after dynamics is saved in the var named "QV_DYN_IN" in "GEOS_AgcmGridComp.F90".
    with computation(PARALLEL), interval(...):
        if GF_ENV_SETTING == 1:
            mass_n = (p_interface_timestep_start[0, 0, 1] - p_interface_timestep_start) * (1.0 / constants.MAPL_GRAV)
            p_n = 0.5 * (p_interface_timestep_start + p_interface_timestep_start[0, 0, 1])
            p_kappa_interface_n = (p_interface_timestep_start / constants.MAPL_P00) ** (constants.MAPL_RGAS / constants.MAPL_CP)
            if K == k_end:
                p_kappa_surface_n = (p_interface_timestep_start[0, 0, 1] / constants.MAPL_P00) ** (constants.MAPL_RGAS / constants.MAPL_CP)
            p_kappa_n = (p_n / constants.MAPL_P00) ** (constants.MAPL_RGAS / constants.MAPL_CP)
            edge_height_above_surface_n = (t_timestep_start / p_kappa_n) * (1.0 + constants.MAPL_VIREPS * vapor_timestep_start)

    with computation(BACKWARD), interval(...):
        if GF_ENV_SETTING == 1:
            if K == k_end:
                layer_height_above_surface_n = 0 + (constants.MAPL_CP / constants.MAPL_GRAV) * (p_kappa_surface_n - p_kappa_n) * edge_height_above_surface_n
            else:
                layer_height_above_surface_n = (
                    edge_height_above_surface_n[0, 0, 1]
                    + (constants.MAPL_CP / constants.MAPL_GRAV) * (p_kappa_interface_n[0, 0, 1] - p_kappa_n) * edge_height_above_surface_n
                )
            edge_height_above_surface_n = (
                layer_height_above_surface_n + (constants.MAPL_CP / constants.MAPL_GRAV) * (p_kappa_n - p_kappa_interface_n) * edge_height_above_surface_n
            )

    with computation(PARALLEL), interval(...):
        if GF_ENV_SETTING == 1:
            if K == k_end:
                dz = -(0 - edge_height_above_surface_n)
            else:
                dz = -(edge_height_above_surface_n[0, 0, 1] - edge_height_above_surface_n)
            air_density = p_n / (287.04 * t_timestep_start * (1.0 + 0.608 * vapor_timestep_start))
            t_flipped = t_timestep_start.at(K=k_end - K)
            p_flipped = p_n.at(K=k_end - K)
            vapor_flipped = vapor_timestep_start.at(K=k_end - K)
            vapor_current_flipped = vapor.at(K=k_end - K)
            u_flipped = u_timestep_start.at(K=k_end - K)
            v_flipped = v_timestep_start.at(K=k_end - K)
            w_flipped = w.at(K=k_end - K)
            layer_height_above_surface_flipped = layer_height_above_surface_n.at(K=k_end - K)
            if K == 0:
                edge_height_above_surface_flipped = 0
            else:
                edge_height_above_surface_flipped = edge_height_above_surface_n.at(K=k_end - K + 1)
            mass_flipped = mass_n.at(K=k_end - K)
            scalar_diffusivity_flipped = scalar_diffusivity.at(K=k_end - K)

            # Grid and sub-grid scale forcings for convection
            grid_scale_forcing_t = dtdt_from_dynamics.at(K=k_end - K) + dtdt_shortwave.at(K=k_end - K) + dtdt_longwave.at(K=k_end - K)
            grid_scale_forcing_vapor = dvapordt_from_dynamics.at(K=k_end - K)
            subgrid_scale_forcing_t = dtdt_pbl.at(K=k_end - K)
            subgrid_scale_forcing_vapor = dspecific_humiditydt_pbl.at(K=k_end - K)
            advective_forcing_t = dtdt_from_dynamics.at(K=k_end - K)

    with computation(FORWARD), interval(0, 1):
        if GF_ENV_SETTING == 1:
            p_surface_flipped = p_interface_timestep_start.at(K=k_end + 1)

    # remainder is the same for both settings
    with computation(PARALLEL), interval(...):
        convective_liquid_flipped = convective_liquid.at(K=k_end - K)
        convective_ice_flipped = convective_ice.at(K=k_end - K)
        convective_cloud_fraction_flipped = convective_cloud_fraction.at(K=k_end - K)
        large_scale_liquid_flipped = large_scale_liquid.at(K=k_end - K)
        large_scale_ice_flipped = large_scale_ice.at(K=k_end - K)
        large_scale_cloud_fraction_flipped = large_scale_cloud_fraction.at(K=k_end - K)

        if ENTRVERSION == 0:
            # eq 6 of https://doi.org/10.1029/2021JD034881
            lateral_entrainment_rate_flipped = 0.71 * max(0.5, w.at(K=k_end - K)) ** (-1.17) * max(0.1, buoyancy.at(K=k_end - K)) ** (-0.36)
        else:
            lateral_entrainment_rate_flipped = 1.0

    with computation(PARALLEL), interval(...):
        # must be in separate computation to ensure lateral_entrainment_rate_local is written before read
        lateral_entrainment_rate = lateral_entrainment_rate_flipped.at(K=k_end - K)

        if CONVECTION_TRACER == 1:
            buoyancy_excess = convection_tracer.at(K=k_end - K)
        else:
            buoyancy_excess = 0.0

    with computation(FORWARD), interval(0, 1):
        # saturation column_water_vapor
        if CONVECTION_TRACER == 1:
            saturation_water_vapor = total_precipitable_water_initial / (1.0e-6 + saturation_total_precipitable_water_initial)
            saturation_water_vapor = min(1.0, max(0.0, saturation_water_vapor))


def copy_into_cumulus_parameterization_state(
    grid_length_local: FloatFieldIJ,
    grid_length: FloatFieldIJ,
    saturation_water_vapor_local: FloatFieldIJ,
    saturation_water_vapor: FloatFieldIJ,
    seed_convection_model_state: FloatFieldIJ,
    seed_convection: FloatFieldIJ,
    convection_fraction_model_state: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    surface_type_model_state: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    grid_scale_forcing_t_local: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor_local: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t_local: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor_local: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    lateral_entrainment_rate_flipped: FloatField,
    lateral_entrainment_rate: FloatField,
):
    """
    One-to-one copy fields into the CumulusParameterization state.

    This division has been inforced for readibility. Some of the earlier functions could
    write directly to the CumulusParameterization state, but local copies have been
    introduced for consistency, so that the data is being written to only one state at a time.

    Args:
        grid_length_local (FloatFieldIJ)
        grid_length (FloatFieldIJ)
        saturation_water_vapor_local (FloatFieldIJ)
        saturation_water_vapor (FloatFieldIJ)
        seed_convection_model_state (FloatFieldIJ)
        seed_convection (FloatFieldIJ)
        convection_fraction_model_state (FloatFieldIJ)
        convection_fraction (FloatFieldIJ)
        surface_type_model_state (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        grid_scale_forcing_t_local (FloatField)
        grid_scale_forcing_t (FloatField)
        grid_scale_forcing_vapor_local (FloatField)
        grid_scale_forcing_vapor (FloatField)
        subgrid_scale_forcing_t_local (FloatField)
        subgrid_scale_forcing_t (FloatField)
        subgrid_scale_forcing_vapor_local (FloatField)
        subgrid_scale_forcing_vapor (FloatField)
        lateral_entrainment_rate_flipped (FloatField)
        lateral_entrainment_rate (FloatField)
    """
    with computation(FORWARD), interval(0, 1):
        grid_length = grid_length_local
        saturation_water_vapor = saturation_water_vapor_local
        seed_convection = seed_convection_model_state
        convection_fraction = convection_fraction_model_state
        surface_type = surface_type_model_state

    with computation(PARALLEL), interval(...):
        grid_scale_forcing_t = grid_scale_forcing_t_local
        grid_scale_forcing_vapor = grid_scale_forcing_vapor_local
        subgrid_scale_forcing_t = subgrid_scale_forcing_t_local
        subgrid_scale_forcing_vapor = subgrid_scale_forcing_vapor_local
        lateral_entrainment_rate = lateral_entrainment_rate_flipped


def prepare_cumulus_paramaterization_state(
    aot500: FloatFieldIJ,
    ccn: FloatFieldIJ,
    ocean_fraction_local: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    p_surface_flipped: FloatFieldIJ,
    p_surface: FloatFieldIJ,
    t_2m_flipped: FloatFieldIJ,
    t_surface: FloatFieldIJ,
    topography_height: FloatFieldIJ,
    topography_height_no_negative: FloatFieldIJ,
    pbl_level_flipped: IntFieldIJ,
    pbl_level: IntFieldIJ,
    latitude_model_state: FloatFieldIJ,
    latitude: FloatFieldIJ,
    longitude_model_state: FloatFieldIJ,
    longitude: FloatFieldIJ,
    rtgt: FloatFieldIJ,
    geopotential_height_forced: FloatField,
    layer_height_above_surface_flipped: FloatField,
    edge_height_above_surface_flipped: FloatField,
    p_flipped: FloatField,
    p_forced: FloatField,
    t_flipped: FloatField,
    t_old: FloatField,
    vapor_flipped: FloatField,
    vapor_old: FloatField,
    air_density: FloatField,
    u_flipped: FloatField,
    u: FloatField,
    v_flipped: FloatField,
    v: FloatField,
    w_flipped: FloatField,
    w: FloatField,
    mass_flipped: FloatField,
    mass: FloatField,
    omega: FloatField,
    buoyancy_excess_local: FloatField,
    buoyancy_excess: FloatField,
    advective_forcing_t: FloatField,
    t_modified_by_advection: FloatField,
    grid_scale_forcing_vapor: FloatField,
    vapor_modified_by_advection: FloatField,
    convective_liquid_flipped: FloatField,
    convective_liquid: FloatField,
    convective_ice_flipped: FloatField,
    convective_ice: FloatField,
    convective_cloud_fraction_flipped: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_liquid_flipped: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice_flipped: FloatField,
    large_scale_ice: FloatField,
    large_scale_cloud_fraction_flipped: FloatField,
    large_scale_cloud_fraction: FloatField,
    convection_tracers: FloatField_ConvectionTracers,
    chemistry_tracers: FloatField_ConvectionTracers,
    sensible_heat_flux_local: FloatFieldIJ,
    sensible_heat_flux: FloatFieldIJ,
    evaporation_local: FloatFieldIJ,
    latent_heat_flux: FloatFieldIJ,
    convective_scale_velocity: FloatFieldIJ,
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
):
    """
    Compute inputs for the cumulus parameterization and fill the rest of the state with non-one-to-one-copies.

    Args:
        aot500 (FloatFieldIJ)
        ccn (FloatFieldIJ)
        ocean_fraction_local (FloatFieldIJ)
        ocean_fraction (FloatFieldIJ)
        p_surface_flipped (FloatFieldIJ)
        p_surface (FloatFieldIJ)
        t_2m_flipped (FloatFieldIJ)
        t_surface (FloatFieldIJ)
        topography_height (FloatFieldIJ)
        topography_height_no_negative (FloatFieldIJ)
        pbl_level_flipped (IntFieldIJ)
        pbl_level (IntFieldIJ)
        latitude_model_state (FloatFieldIJ)
        latitude (FloatFieldIJ)
        longitude_model_state (FloatFieldIJ)
        longitude (FloatFieldIJ)
        rtgt (FloatFieldIJ)
        geopotential_height_forced (FloatField)
        layer_height_above_surface_flipped (FloatField)
        edge_height_above_surface_flipped (FloatField)
        p_flipped (FloatField)
        p_forced (FloatField)
        t_flipped (FloatField)
        t_old (FloatField)
        vapor_flipped (FloatField)
        vapor_old (FloatField)
        air_density (FloatField)
        u_flipped (FloatField)
        u (FloatField)
        v_flipped (FloatField)
        v (FloatField)
        w_flipped (FloatField)
        w (FloatField)
        mass_flipped (FloatField)
        mass (FloatField)
        omega (FloatField)
        buoyancy_excess_local (FloatField)
        buoyancy_excess (FloatField)
        advective_forcing_t (FloatField)
        t_modified_by_advection (FloatField)
        grid_scale_forcing_vapor (FloatField)
        vapor_modified_by_advection (FloatField)
        convective_liquid_flipped (FloatField)
        convective_liquid (FloatField)
        convective_ice_flipped (FloatField)
        convective_ice (FloatField)
        convective_cloud_fraction_flipped (FloatField)
        convective_cloud_fraction (FloatField)
        large_scale_liquid_flipped (FloatField)
        large_scale_liquid (FloatField)
        large_scale_ice_flipped (FloatField)
        large_scale_ice (FloatField)
        large_scale_cloud_fraction_flipped (FloatField)
        large_scale_cloud_fraction (FloatField)
        convection_tracers (FloatField_ConvectionTracers)
        chemistry_tracers (FloatField_ConvectionTracers)
        sensible_heat_flux_local (FloatFieldIJ)
        sensible_heat_flux (FloatFieldIJ)
        evaporation_local (FloatFieldIJ)
        latent_heat_flux (FloatFieldIJ)
        convective_scale_velocity (FloatFieldIJ)
        t_excess (FloatFieldIJ)
        vapor_excess (FloatFieldIJ)
    """
    from __externals__ import APPLY_SUBSIDENCE_MICROPHYSICS, AUTOCONV, DT_MOIST, USE_TRACER_TRANSPORT, k_end

    with computation(FORWARD), interval(0, 1):
        if AUTOCONV == 2:
            ccn = max(100.0, (370.37 * (0.01 + max(0.0, aot500))) ** 1.555)
        else:
            ccn = 100.0

        ocean_fraction = ocean_fraction_local
        p_surface = p_surface_flipped * 1.0e-2  # mbar
        t_surface = t_2m_flipped
        topography_height_no_negative = max(0.0, topography_height)
        pbl_level = pbl_level_flipped
        latitude = latitude_model_state * 180.0 / 3.14159
        longitude = longitude_model_state * 180.0 / 3.14159

    with computation(PARALLEL), interval(0, -1):
        # heights, current pressure, temperature and water vapor mixing ratio
        geopotential_height_forced = layer_height_above_surface_flipped * rtgt + topography_height
        p_forced = p_flipped * 1.0e-2  # mbar
        t_old = t_flipped
        vapor_old = vapor_flipped  # at beginning of the timestep

        # air density, TKE and cloud liquid water mixing ratio
        air_density = 1.0e2 * p_forced / (287.04 * t_old * (1.0 + 0.608 * vapor_old))

        # wind velocities
        u = u_flipped
        v = v_flipped
        w = w_flipped
        mass = mass_flipped
        omega = -constants.MAPL_GRAV * air_density * w

        # buoyancy excess
        buoyancy_excess = buoyancy_excess_local

        # temperature/water vapor modified only by advection
        t_modified_by_advection = t_old + advective_forcing_t * DT_MOIST
        vapor_modified_by_advection = vapor_old + grid_scale_forcing_vapor * DT_MOIST

        if APPLY_SUBSIDENCE_MICROPHYSICS == 1:
            # microphysics ice and liquid mixing ratio, and cloud fraction of the host model
            # (only subsidence is applied)
            convective_liquid = convective_liquid_flipped
            convective_ice = convective_ice_flipped
            convective_cloud_fraction = convective_cloud_fraction_flipped
            large_scale_liquid = large_scale_liquid_flipped
            large_scale_ice = large_scale_ice_flipped
            large_scale_cloud_fraction = large_scale_cloud_fraction_flipped

    with computation(PARALLEL), interval(...):
        if USE_TRACER_TRANSPORT == 1:
            tracer = 0
            while tracer < constants.NUMBER_OF_TRACERS:
                chemistry_tracers[0, 0, 0][tracer] = max(convection_tracers.at(K=k_end - K, ddim=[tracer]), constants.FLOAT_TINY)
                tracer += 1

    with computation(FORWARD), interval(0, 1):
        pbl_internal: FloatFieldIJ = geopotential_height_forced.at(K=pbl_level) - topography_height

    # NOTE variables in this section have unintelligible names
    # could not desipher during Fortran --> Python port, so Fortran names remain
    with computation(FORWARD), interval(0, 1):
        # get execess temperature and water vapor for source air parcels
        pten = t_old.at(K=0)
        pqen = vapor_old.at(K=0)
        paph = 100.0 * p_surface
        zrho = paph / (287.04 * (t_old.at(K=0) * (1.0 + 0.608 * vapor_old.at(K=0))))

        # sensible and latent sfc fluxes for the heat-engine closure
        sensible_heat_flux = zrho * cumulus_parameterization_constants.CP * sensible_heat_flux_local  # W/m^2
        latent_heat_flux = zrho * cumulus_parameterization_constants.XLV * evaporation_local  # W/m^2

        # local le and h fluxes for W*
        pahfs = -sensible_heat_flux_local * zrho * 1004.64  # W/m^2
        pqhfl = -evaporation_local  # kg/m^2/s

        # buoyancy flux (h+le)
        zkhvfl = (pahfs / 1004.64 + 0.608 * pten * pqhfl) / zrho  # K m s-1

        # depth of 1st model layer
        # geopotential_height_forced.at(K=0) - topography_height is ~ 1/2 of the depth
        # of 1st model layer, so multiply by 2
        pgeoh = 2.0 * (geopotential_height_forced.at(K=0) - topography_height) * constants.MAPL_GRAV  # m+2 s-2

        # convective-scale velocity w*
        # in the future, change 0.001 by ustar^3
        convective_scale_velocity = max(0.0, 0.001 - 1.5 * 0.41 * zkhvfl * pgeoh / pten)  # m+3 s-3

        if convective_scale_velocity > constants.FLOAT_TINY:
            # convective-scale velocity w*
            convective_scale_velocity = 1.2 * convective_scale_velocity**0.3333

            # temperature excess
            t_excess = max(0.0, -1.5 * pahfs / (zrho * convective_scale_velocity * 1004.64))  # K

            # moisture excess
            vapor_excess = max(0.0, -1.5 * pqhfl / (zrho * convective_scale_velocity))  # kg kg-1

        # convective_scale_velocity for shallow convection closure (Grant 2001)

        # depth of the pbl
        pgeoh = pbl_internal * constants.MAPL_GRAV

        # convective_scale_velocity W* (m/s)
        convective_scale_velocity = max(0.0, 0.001 - 1.5 * 0.41 * zkhvfl * pgeoh / pten)
        convective_scale_velocity = 1.2 * convective_scale_velocity**0.3333


class GF2020Setup(NDSLRuntime):
    """
    This class performs the entire setup sequence for the GF2020 convection parameterization scheme

    In the source Fortran codee, this code is split across three subroutines nested as follows:

    - GF_Run
        - GF2020_INTERFACE
            - GF2020_DRV up to "------ CALL CUMULUS PARAMETERIZATION"

    This python implementation simplifies this structure by bringing all setup calculations to the same level.
    An effort has been made to reduce duplicate/unnecessary locals where possible, but some have been retained
    for the sake of readibility.
    """

    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        config: GF2020Config,
        saturation_tables: SaturationVaporPressureTable,
    ):
        super().__init__(stencil_factory)

        # make inputs visible at runtime
        self.stencil_factory = stencil_factory
        self.config = config
        self.saturation_tables = saturation_tables

        # check config for unimplemented paths
        if config.ADV_TRIGGER == 2:
            raise NotImplementedError(
                "[NDSL] GF2020-->Setup initialized with shallow plume enabled. This requires"
                "an umplemented portion of ensemble_output_and_feedback. Please impelment, then disable this"
                "error manually to proceed."
            )

        # Construct stencils
        self._compute_extra_inputs_from_state = stencil_factory.from_dims_halo(
            func=compute_extra_inputs_from_state,
            compute_dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
            externals={
                "STOCHASTIC_CONVECTION": config.STOCHASTIC_CNV,
                "STOCH_TOP": config.STOCH_TOP,
                "STOCH_BOT": config.STOCH_BOT,
                "GF_MIN_AREA": config.GF_MIN_AREA,
                "LHYDROSTATIC": config.LHYDROSTATIC,
            },
        )

        self._pass_back_to_model_state = stencil_factory.from_dims_halo(
            func=pass_back_to_model_state,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._zero_state = stencil_factory.from_dims_halo(
            func=zero_state,
            compute_dims=[I_DIM, J_DIM, K_INTERFACE_DIM],
        )

        self._prefill_cumulus_parameterization_state = stencil_factory.from_dims_halo(
            func=prefill_cumulus_parameterization_state,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "NUMBER_OF_PLUMES": cumulus_parameterization_constants.NUMBER_OF_PLUMES,
                "APPLY_SUBSIDENCE_MICROPHYSICS": config.APPLY_SUBSIDENCE_MICROPHYSICS,
            },
        )

        self._prefill_locals = stencil_factory.from_dims_halo(
            func=prefill_locals,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._set_2d_fields = stencil_factory.from_dims_halo(
            func=set_2d_fields,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "SIZE_I_DIM": stencil_factory.grid_indexing.domain[0],
                "SIZE_J_DIM": stencil_factory.grid_indexing.domain[1],
            },
        )

        self._choose_environment_and_flip_k_axis = stencil_factory.from_dims_halo(
            func=choose_environment_and_flip_k_axis,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "GF_ENV_SETTING": config.GF_ENV_SETTING,
                "ENTRVERSION": config.ENTRVERSION,
                "CONVECTION_TRACER": config.CONVECTION_TRACER,
            },
        )

        self._copy_into_cumulus_parameterization_state = stencil_factory.from_dims_halo(
            func=copy_into_cumulus_parameterization_state,
            compute_dims=[I_DIM, J_DIM, K_DIM],
        )

        self._prepare_cumulus_paramaterization_state = stencil_factory.from_dims_halo(
            func=prepare_cumulus_paramaterization_state,
            compute_dims=[I_DIM, J_DIM, K_DIM],
            externals={
                "AUTOCONV": config.AUTOCONV,
                "DT_MOIST": config.DT_MOIST,
                "APPLY_SUBSIDENCE_MICROPHYSICS": config.APPLY_SUBSIDENCE_MICROPHYSICS,
                "USE_TRACER_TRANSPORT": config.USE_TRACER_TRANSPORT,
            },
        )

    def __call__(
        self,
        state: GF2020State,
        locals: GF2020Locals,
        cumulus_parameterization_state: GF2020CumulusParameterizationState,
        convection_tracers: ConvectionTracers,
    ):
        """
        Perform setup calculations

        Args:
            state (GF2020State): NDSL State containing all model fields required for GF2020.
            locals (GF2020Locals): NDSL LocalState containing all locals for GF2020.
            cumulus_parameterization_state (GF2020CumulusParameterizationState): NDSL State containing all
                fields required for the CumulusParameterization.
            convection_tracers (ConvectionTracers): Collection of tracers from the rest of the model which
                will be updated within convection. These may come from a variety of sources, and need to be
                collected into the expected ConvectionTracers data type before being passed down.

        """
        self._compute_extra_inputs_from_state(
            p_interface=state.p_interface,
            p=locals.derived_state.p,
            p_kappa=locals.derived_state.p_kappa,
            edge_height_above_surface=locals.derived_state.edge_height_above_surface,
            layer_height_above_surface=locals.derived_state.layer_height_above_surface,
            geopotential_height_interface=state.geopotential_height_interface,
            t=state.t,
            th=locals.derived_state.th,
            vapor=state.vapor,
            mass=locals.derived_state.mass,
            w=state.w,
            omega=state.omega,
            vertical_motion=locals.derived_state.vertical_velocity,
            tpwi=state.total_precipitable_water_initial,
            tpwi_star=state.saturation_total_precipitable_water_initial,
            seed_convection=locals.derived_state.seed_convection,
            area=state.area,
            modified_area=locals.derived_state.modified_area,
            convection_fraction=state.convection_fraction,
            ese=self.saturation_tables.ese,
            esx=self.saturation_tables.esx,
        )

        if state.seed_convection is not None:
            self._pass_back_to_model_state(
                local_seed_convection=locals.derived_state.seed_convection,
                model_state_seed_convection=state.seed_convection,
            )

        self._zero_state(
            dvapordt_deep_convection=state.dvapordt_deep_convection,
            dtdt_deep_convection=state.dtdt_deep_convection,
            dudt_deep_convection=state.dudt_deep_convection,
            dvdt_deep_convection=state.dvdt_deep_convection,
            sigma_deep=state.sigma_deep,
            sigma_mid=state.sigma_mid,
            mass_flux_shalow=state.mass_flux_shallow,
            mass_flux_mid=state.mass_flux_mid,
            mass_flux_deep_updraft=state.mass_flux_deep_updraft,
            mass_flux_deep_updraft_interface=state.mass_flux_deep_updraft_interface,
            mass_flux_deep_updraft_detrained=state.mass_flux_deep_updraft_detrained,
            mass_flux_deep_downdraft=state.mass_flux_deep_downdraft,
            mass_flux_cloud_base=state.mass_flux_cloud_base,
            mass_flux_cloud_base_shallow=state.mass_flux_cloud_base_shallow,
            mass_flux_cloud_base_mid=state.mass_flux_cloud_base_mid,
            mass_flux_cloud_base_deep=state.mass_flux_cloud_base_deep,
            convection_code_shallow=state.convection_code_shallow,
            convection_code_mid=state.convection_code_mid,
            convection_code_deep=state.convection_code_deep,
            cloud_workfunction_0=state.cloud_workfunction_0,
            cloud_workfunction_1=state.cloud_workfunction_1,
            cloud_workfunction_2=state.cloud_workfunction_2,
            cloud_workfunction_3=state.cloud_workfunction_3,
            cloud_workfunction_1_pbl=state.cloud_workfunction_1_pbl,
            cloud_workfunction_1_cin=state.cloud_workfunction_1_cin,
            convective_precipitation_RAS=state.convective_precipitation_RAS,
            convective_precipitation_GF=state.convective_precipitation_GF,
            convective_condensate_source=state.convective_condensate_source,
            convective_condensate_grid_mean=state.convective_condensate_grid_mean,
            total_water_flux_deep_convection_interface=state.total_water_flux_deep_convection_interface,
            updraft_area_fraction=state.updraft_areal_fraction,
            updraft_vertical_velocity=state.updraft_vertical_velocity,
            entrainment_parameter=state.entrainment_parameter,
            lightning_density=state.lightning_density,
            pbl_time_scale=state.pbl_time_scale,
            cape_removal_time_scale=state.cape_removal_time_scale,
        )

        self._prefill_cumulus_parameterization_state(
            error_code=cumulus_parameterization_state.output.error_code,
            downdraft_origin_level=cumulus_parameterization_state.output.downdraft_origin_level,
            lcl_level=cumulus_parameterization_state.output.lcl_level,
            updraft_origin_level=cumulus_parameterization_state.output.updraft_origin_level,
            updraft_lfc_level=cumulus_parameterization_state.output.updraft_lfc_level,
            cloud_top_level=cumulus_parameterization_state.output.cloud_top_level,
            kstabi=cumulus_parameterization_state.output.kstabi,
            kstabm=cumulus_parameterization_state.output.kstabm,
            precip=cumulus_parameterization_state.output.precip,
            cloud_base_mass_flux_modified=cumulus_parameterization_state.output.cloud_base_mass_flux_modified,
            epsilon_forced=cumulus_parameterization_state.output.epsilon_forced,
            total_normalized_integrated_condensate_forced=cumulus_parameterization_state.output.total_normalized_integrated_condensate_forced,
            scale_dependence_factor=cumulus_parameterization_state.output.scale_dependence_factor,
            p_cloud_levels_forced=cumulus_parameterization_state.output.p_cloud_levels_forced,
            entrainment_rate=cumulus_parameterization_state.output.entrainment_rate,
            mass_entrainment_updraft_forced=cumulus_parameterization_state.output.mass_entrainment_updraft_forced,
            mass_entrainment_downdraft_forced=cumulus_parameterization_state.output.mass_entrainment_downdraft_forced,
            mass_detrainment_updraft_forced=cumulus_parameterization_state.output.mass_detrainment_updraft_forced,
            mass_detrainment_downdraft_forced=cumulus_parameterization_state.output.mass_detrainment_downdraft_forced,
            normalized_massflux_updraft_forced=cumulus_parameterization_state.output.normalized_massflux_updraft_forced,
            normalized_massflux_downdraft_forced=cumulus_parameterization_state.output.normalized_massflux_downdraft_forced,
            condensate_to_fall_forced=cumulus_parameterization_state.output.condensate_to_fall_forced,
            evaporate_in_downdraft_forced=cumulus_parameterization_state.output.evaporate_in_downdraft_forced,
            cloud_liquid_after_rain_forced=cumulus_parameterization_state.output.cloud_liquid_after_rain_forced,
            t_updraft=cumulus_parameterization_state.output.t_updraft,
            convective_cloud_fraction_output=cumulus_parameterization_state.output.convective_cloud_fraction,
            dtdt=cumulus_parameterization_state.output.dtdt,
            dudt=cumulus_parameterization_state.output.dudt,
            dvdt=cumulus_parameterization_state.output.dvdt,
            dvapordt=cumulus_parameterization_state.output.dvapordt,
            dcloudicedt=cumulus_parameterization_state.output.dcloudicedt,
            dnicedt=cumulus_parameterization_state.output.dnicedt,
            dnliquiddt=cumulus_parameterization_state.output.dnliquiddt,
            dbuoyancydt=cumulus_parameterization_state.output.dbuoyancydt,
            chemistry_tracers_output=cumulus_parameterization_state.input_output.chemistry_tracers_output,
            evaporation_sublimation_tendency=cumulus_parameterization_state.output.evaporation_sublimation_tendency,
            convective_precip_flux=cumulus_parameterization_state.output.convective_precip_flux,
            t_perturbation=cumulus_parameterization_state.output.t_perturbation,
            omega=cumulus_parameterization_state.input_output.omega,
            large_scale_ice=cumulus_parameterization_state.input_output.large_scale_ice,
            convective_ice=cumulus_parameterization_state.input_output.convective_ice,
            large_scale_liquid=cumulus_parameterization_state.input_output.large_scale_liquid,
            convective_liquid=cumulus_parameterization_state.input_output.convective_liquid,
            large_scale_cloud_fraction=cumulus_parameterization_state.input_output.large_scale_cloud_fraction,
            convective_cloud_fraction=cumulus_parameterization_state.input_output.convective_cloud_fraction,
            lightning_density=cumulus_parameterization_state.output.lightning_density,
            t_excess=cumulus_parameterization_state.input.t_excess,
            vapor_excess=cumulus_parameterization_state.input.vapor_excess,
            last_error_code=cumulus_parameterization_state.input.last_error_code,
        )

        self._prefill_locals(
            rtgt=locals.rtgt,
            precip=locals.precip,
            t_tendency_from_vapor=locals.t_tendency_from_vapor,
            dtdt=locals.dtdt,
            dvapordt=locals.dvapordt,
            dcloudicedt=locals.dcloudicedt,
            dudt=locals.dudt,
            dvdt=locals.dvdt,
            dlarge_scale_icedt=locals.dlarge_scale_icedt,
            dconvective_icedt=locals.dconvective_icedt,
            dlarge_scale_liquiddt=locals.dlarge_scale_liquiddt,
            dconvective_liquiddt=locals.dconvective_liquiddt,
            dlarge_scale_cloud_fractiondt=locals.dlarge_scale_cloud_fractiondt,
            dconvective_cloud_fractiondt=locals.dconvective_cloud_fractiondt,
            dbuoyancydt=locals.dbuoyancydt,
            evaporation_sublimation_tendency=locals.evaporation_sublimation_tendency,
            convective_precip_flux=locals.convective_precip_flux,
        )

        # workaround because max of full field cannot be determined inside a stencil
        t_2m_max = Float(state.t_2m.field.max().item())

        # # if surface temperature is not yet set in single column mode, stop the entire convection scheme
        if self.stencil_factory.grid_indexing.get_shape([I_DIM, J_DIM]) == (1, 1) and t_2m_max < 1.0e-6:
            # NOTE this value goes into scm_stop - needs to be made a part of the LocalState, but currently
            # LocalStates cannot support scalars
            return True

        self._set_2d_fields(
            aot500=locals.aot500,
            t=state.t,
            t_2m_max=t_2m_max,
            t_2m=state.t_2m,
            t_2m_local=locals.flipped_copy.t_2m,
            evaporation=state.evaporation,
            evaporation_local=locals.evaporation,
            sensible_heat_flux=state.sensible_heat_flux,
            sensible_heat_flux_local=locals.sensible_heat_flux,
            p_interface=state.p_interface,
            vapor=state.vapor,
            geopotential_height_surface=state.geopotential_height_surface,
            topography_height=locals.topography_height,
            land_fraction=state.land_fraction,
            ocean_fraction=locals.ocean_fraction,
            area=state.area,
            grid_length=locals.grid_length,
            pbl_level=state.pbl_level,
            pbl_level_flipped=locals.flipped_copy.pbl_level,
        )

        self._choose_environment_and_flip_k_axis(
            dz=locals.derived_state.dz,
            air_density=locals.derived_state.air_density,
            geopotential_height_interface=state.geopotential_height_interface,
            t=state.t,
            t_flipped=locals.flipped_copy.t,
            t_timestep_start=state.t_timestep_start,
            p=locals.derived_state.p,
            p_interface=state.p_interface,
            p_interface_timestep_start=state.p_interface_timestep_start,
            p_flipped=locals.flipped_copy.p,
            p_surface_flipped=locals.flipped_copy.p_surface,
            vapor=state.vapor,
            vapor_timestep_start=state.vapor_timestep_start,
            vapor_flipped=locals.flipped_copy.vapor,
            vapor_current_flipped=locals.flipped_copy.vapor_current,
            u=state.u,
            u_timestep_start=state.u_timestep_start,
            u_flipped=locals.flipped_copy.u,
            v=state.v,
            v_timestep_start=state.v_timestep_start,
            v_flipped=locals.flipped_copy.v,
            w=state.w,
            w_flipped=locals.flipped_copy.w,
            layer_height_above_surface=locals.derived_state.layer_height_above_surface,
            layer_height_above_surface_flipped=locals.flipped_copy.layer_height_above_surface,
            edge_height_above_surface=locals.derived_state.edge_height_above_surface,
            edge_height_above_surface_flipped=locals.flipped_copy.edge_height_above_surface,
            mass=locals.derived_state.mass,
            mass_flipped=locals.flipped_copy.mass,
            scalar_diffusivity=state.scalar_diffusivity,
            scalar_diffusivity_flipped=locals.flipped_copy.scalar_diffusivity,
            lateral_entrainment_rate=state.lateral_entrainment_rate,
            lateral_entrainment_rate_flipped=locals.flipped_copy.lateral_entrainment_rate,
            buoyancy=state.buoyancy,
            buoyancy_excess=locals.buoyancy_excess,
            dtdt_shortwave=state.dtdt_shortwave,
            dtdt_longwave=state.dtdt_longwave,
            dtdt_from_dynamics=state.dtdt_from_dynamics,
            dtdt_pbl=state.dtdt_pbl,
            dvapordt_from_dynamics=state.dvapordt_from_dynamics,
            dspecific_humiditydt_pbl=state.dspecific_humiditydt_pbl,
            grid_scale_forcing_t=locals.grid_scale_forcing_t,
            grid_scale_forcing_vapor=locals.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t=locals.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor=locals.subgrid_scale_forcing_vapor,
            advective_forcing_t=locals.advective_forcing_t,
            convective_liquid=state.convective_liquid,
            convective_liquid_flipped=locals.flipped_copy.convective_liquid,
            convective_ice=state.convective_ice,
            convective_ice_flipped=locals.flipped_copy.convective_ice,
            convective_cloud_fraction=state.convective_cloud_fraction,
            convective_cloud_fraction_flipped=locals.flipped_copy.convective_cloud_fraction,
            large_scale_liquid=state.large_scale_liquid,
            large_scale_liquid_flipped=locals.flipped_copy.large_scale_liquid,
            large_scale_ice=state.large_scale_ice,
            large_scale_ice_flipped=locals.flipped_copy.large_scale_ice,
            large_scale_cloud_fraction=state.large_scale_cloud_fraction,
            large_scale_cloud_fraction_flipped=locals.flipped_copy.large_scale_cloud_fraction,
            convection_tracer=state.convection_tracer,
            total_precipitable_water_initial=state.total_precipitable_water_initial,
            saturation_total_precipitable_water_initial=state.saturation_total_precipitable_water_initial,
            saturation_water_vapor=locals.saturation_water_vapor,
        )

        if self.config.ADV_TRIGGER == 2:
            raise NotImplementedError("option not implemented, should have been caught at initialization")

        self._copy_into_cumulus_parameterization_state(
            grid_length_local=locals.grid_length,
            grid_length=cumulus_parameterization_state.input_output.grid_length,
            saturation_water_vapor_local=locals.saturation_water_vapor,
            saturation_water_vapor=cumulus_parameterization_state.input.saturation_water_vapor,
            seed_convection_model_state=locals.derived_state.seed_convection,
            seed_convection=cumulus_parameterization_state.input.seed_convection,
            convection_fraction_model_state=state.convection_fraction,
            convection_fraction=cumulus_parameterization_state.input.convection_fraction,
            surface_type_model_state=state.surface_type,
            surface_type=cumulus_parameterization_state.input.surface_type,
            grid_scale_forcing_t_local=locals.grid_scale_forcing_t,
            grid_scale_forcing_t=cumulus_parameterization_state.input.grid_scale_forcing_t,
            grid_scale_forcing_vapor_local=locals.grid_scale_forcing_vapor,
            grid_scale_forcing_vapor=cumulus_parameterization_state.input.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t_local=locals.subgrid_scale_forcing_t,
            subgrid_scale_forcing_t=cumulus_parameterization_state.input.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor_local=locals.subgrid_scale_forcing_vapor,
            subgrid_scale_forcing_vapor=cumulus_parameterization_state.input.subgrid_scale_forcing_vapor,
            lateral_entrainment_rate_flipped=locals.flipped_copy.lateral_entrainment_rate,
            lateral_entrainment_rate=cumulus_parameterization_state.input.lateral_entrainment_rate,
        )

        self._prepare_cumulus_paramaterization_state(
            aot500=locals.aot500,
            ccn=cumulus_parameterization_state.input_output.ccn,
            ocean_fraction_local=locals.ocean_fraction,
            ocean_fraction=cumulus_parameterization_state.input.ocean_fraction,
            p_surface_flipped=locals.flipped_copy.p_surface,
            p_surface=cumulus_parameterization_state.input_output.p_surface,
            t_2m_flipped=locals.flipped_copy.t_2m,
            t_surface=cumulus_parameterization_state.input_output.t_surface,
            topography_height=locals.topography_height,
            topography_height_no_negative=cumulus_parameterization_state.input_output.topography_height_no_negative,
            pbl_level_flipped=locals.flipped_copy.pbl_level,
            pbl_level=cumulus_parameterization_state.input_output.pbl_level,
            latitude_model_state=state.latitude,
            latitude=cumulus_parameterization_state.input_output.latitude_degrees,
            longitude_model_state=state.longitude,
            longitude=cumulus_parameterization_state.input_output.longitude_degrees,
            rtgt=locals.rtgt,
            geopotential_height_forced=cumulus_parameterization_state.input_output.geopotential_height_forced,
            layer_height_above_surface_flipped=locals.flipped_copy.layer_height_above_surface,
            edge_height_above_surface_flipped=locals.flipped_copy.edge_height_above_surface,
            p_flipped=locals.flipped_copy.p,
            p_forced=cumulus_parameterization_state.input_output.p_forced,
            t_flipped=locals.flipped_copy.t,
            t_old=cumulus_parameterization_state.input_output.t_old,
            vapor_flipped=locals.flipped_copy.vapor,
            vapor_old=cumulus_parameterization_state.input_output.vapor_old,
            air_density=cumulus_parameterization_state.input_output.air_density,
            u_flipped=locals.flipped_copy.u,
            u=cumulus_parameterization_state.input_output.u,
            v_flipped=locals.flipped_copy.v,
            v=cumulus_parameterization_state.input_output.v,
            w_flipped=locals.flipped_copy.w,
            w=cumulus_parameterization_state.input_output.w,
            mass_flipped=locals.flipped_copy.mass,
            mass=cumulus_parameterization_state.input_output.mass,
            omega=cumulus_parameterization_state.input_output.omega,
            buoyancy_excess_local=locals.buoyancy_excess,
            buoyancy_excess=cumulus_parameterization_state.input_output.buoyancy_excess,
            advective_forcing_t=locals.advective_forcing_t,
            t_modified_by_advection=cumulus_parameterization_state.input_output.t_modified_by_advection,
            grid_scale_forcing_vapor=cumulus_parameterization_state.input.grid_scale_forcing_vapor,
            vapor_modified_by_advection=cumulus_parameterization_state.input_output.vapor_modified_by_advection,
            convective_liquid_flipped=locals.flipped_copy.convective_liquid,
            convective_liquid=cumulus_parameterization_state.input_output.convective_liquid,
            convective_ice_flipped=locals.flipped_copy.convective_ice,
            convective_ice=cumulus_parameterization_state.input_output.convective_ice,
            convective_cloud_fraction_flipped=locals.flipped_copy.convective_cloud_fraction,
            convective_cloud_fraction=cumulus_parameterization_state.input_output.convective_cloud_fraction,
            large_scale_liquid_flipped=locals.flipped_copy.large_scale_liquid,
            large_scale_liquid=cumulus_parameterization_state.input_output.large_scale_liquid,
            large_scale_ice_flipped=locals.flipped_copy.large_scale_ice,
            large_scale_ice=cumulus_parameterization_state.input_output.large_scale_ice,
            large_scale_cloud_fraction_flipped=locals.flipped_copy.large_scale_cloud_fraction,
            large_scale_cloud_fraction=cumulus_parameterization_state.input_output.large_scale_cloud_fraction,
            convection_tracers=convection_tracers.tracers,
            chemistry_tracers=cumulus_parameterization_state.input_output.chemistry_tracers,
            sensible_heat_flux_local=locals.sensible_heat_flux,
            sensible_heat_flux=cumulus_parameterization_state.input_output.sensible_heat_flux,
            evaporation_local=locals.evaporation,
            latent_heat_flux=cumulus_parameterization_state.input_output.latent_heat_flux,
            convective_scale_velocity=cumulus_parameterization_state.input_output.convective_scale_velocity,
            t_excess=cumulus_parameterization_state.input.t_excess,
            vapor_excess=cumulus_parameterization_state.input.vapor_excess,
        )

        # NOTE this value goes into scm_stop - needs to be made a part of the LocalState, but currently
        # LocalStates cannot support scalars
        return False
