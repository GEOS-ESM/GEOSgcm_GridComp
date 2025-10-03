import copy
from ndsl import StencilFactory, QuantityFactory, ndsl_log
from ndsl.dsl.gt4py import PARALLEL, interval, computation, FORWARD, sqrt, max, min, abs, floor, BACKWARD
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, K
import pyMoist.constants as constants
from pyMoist.saturation_tables.types import GlobalTable_saturation_tables
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.temporaries import GF2020Temporaries


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
     - Initalize stochastic variability for convection
     - Modify area (m^2) here so GF scale dependence has a convection_fraction dependence

    This stencil MUST be built using Z_INTERFACE_DIM to function properly.
    """
    from __externals__ import k_end, STOCHASTIC_CONVECTION, STOCH_TOP, STOCH_BOT, GF_MIN_AREA, LHYDROSTATIC

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
            vertical_motion = (
                -1
                * omega
                / (
                    constants.MAPL_GRAV
                    * p
                    / (constants.MAPL_RDRY * t * (1.0 + constants.MAPL_VIREPS * vapor))
                )
            )
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

        # Initalize stochastic variability for convection
        if STOCHASTIC_CONVECTION == True:  # noqa
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


def zero_state_fields(
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
    cloud_work_function_0: FloatFieldIJ,
    cloud_work_function_1: FloatFieldIJ,
    cloud_work_function_2: FloatFieldIJ,
    cloud_work_function_3: FloatFieldIJ,
    cloud_work_function_1_pbl: FloatFieldIJ,
    cloud_work_function_1_cin: FloatFieldIJ,
    convective_precipitation_RAS: FloatField,
    convective_precipitation_GF: FloatFieldIJ,
    convective_condensate_source: FloatField,
    convective_condensate_grid_mean: FloatField,
    total_water_flux_deep_convection: FloatField,
    updraft_area_fraction: FloatField,
    updraft_vertical_velocity: FloatField,
    lighting_density: FloatFieldIJ,
    entrainment_parameter: FloatField,
    pbl_time_scale: FloatFieldIJ,
    cape_removal_time_scale: FloatFieldIJ,
):
    """
    Must be built with Z_INTERFACE_DIM.
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

        cloud_work_function_0 = 0.0
        cloud_work_function_1 = 0.0
        cloud_work_function_2 = 0.0
        cloud_work_function_3 = 0.0
        cloud_work_function_1_pbl = 0.0
        cloud_work_function_1_cin = 0.0

    with computation(PARALLEL), interval(0, -1):
        convective_precipitation_RAS = 0.0

    with computation(FORWARD), interval(0, 1):
        convective_precipitation_GF = 0.0

    with computation(PARALLEL), interval(0, -1):
        convective_condensate_source = 0.0

    with computation(PARALLEL), interval(...):
        convective_condensate_grid_mean = 0.0
        total_water_flux_deep_convection = 0.0

    with computation(PARALLEL), interval(0, -1):
        updraft_area_fraction = 0.0
        updraft_vertical_velocity = 0.0

    with computation(FORWARD), interval(0, 1):
        lighting_density = 0.0

    with computation(PARALLEL), interval(0, -1):
        entrainment_parameter = 0.0

    with computation(FORWARD), interval(0, 1):
        pbl_time_scale = 0.0
        cape_removal_time_scale = 0.0


def compute_2d_inputs(
    t: FloatField,
    t_2m_max: Float,
    t_2m: FloatFieldIJ,
    t_2m_local: FloatFieldIJ,
    evaporation: FloatFieldIJ,
    evaporation_local: FloatFieldIJ,
    sensible_heat_flux: FloatFieldIJ,
    p_interface: FloatField,
    vapor: FloatField,
    geopotential_height_surface: FloatFieldIJ,
    topography_height: FloatFieldIJ,
    land_fraction: FloatFieldIJ,
    ocean_fraction: FloatFieldIJ,
    area: FloatFieldIJ,
    grid_length: FloatFieldIJ,
    pbl_level: FloatFieldIJ,
    pbl_level_local: FloatFieldIJ,
):
    from __externals__ import k_end, SINGLE_COLUMN_MODE

    with computation(FORWARD), interval(0, 1):
        if t_2m_max < 1.0e-6:
            # in case convection is called before surface data is computed
            t_2m_local = t.at(K=k_end)
        else:
            t_2m_local = t_2m
        # moisture flux from sfc
        evaporation_local = evaporation  # kg m-2 s-1
        # sensible–heat_flux comes in W m-2, below it is converted to K m s-1
        sensible_heat_flux_local = sensible_heat_flux / (
            1004.0 * p_interface.at(K=k_end) / (287.04 * t.at(K=k_end) * (1.0 + 0.608 * vapor.at(K=k_end)))
        )  # K m s-1
        # topography height  (m)
        topography_height = geopotential_height_surface / constants.MAPL_GRAV
        # land/ocean fraction: land if < 1 ,ocean if = 1
        ocean_fraction = 1.0 - land_fraction
        # grid length for the scale awareness (m)
        grid_length = sqrt(area)
        # special setting for SCM runs
        if SINGLE_COLUMN_MODE == True:
            dx2d = 100000.0

        # TODO not confident this is correct
        if pbl_level != 0.0:
            pbl_level_local = k_end - int(round(pbl_level))
        else:
            pbl_level_local = 1


def set_local_state(
    geopotential_height_interface: FloatField,
    t: FloatField,
    t_local: FloatField,
    p: FloatField,
    p_local: FloatField,
    vapor: FloatField,
    vapor_local: FloatField,
    vapor_current_local: FloatField,
    u: FloatField,
    u_local: FloatField,
    v: FloatField,
    v_local: FloatField,
    vertical_velocity: FloatField,
    vertical_velocity_local: FloatField,
    layer_height_above_surface: FloatField,
    layer_height_above_surface_local: FloatField,
    edge_height_above_surface: FloatField,
    edge_height_above_surface_local: FloatField,
    mass: FloatField,
    mass_local: FloatField,
    scalar_diffusivity: FloatField,
    scalar_diffusivity_local: FloatField,
    lateral_entrainment_rate_local: FloatField,
    lateral_entrainment_rate: FloatField,
    p_surface: FloatFieldIJ,
    buoyancy: FloatField,
    p_interface: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    large_scale_cloud_fraction: FloatField,
    p_interface_timestep_start: FloatField,
    t_timestep_start: FloatField,
    u_timestep_start: FloatField,
    v_timestep_start: FloatField,
    vapor_timestep_start: FloatField,
    dtdt_shortwave: FloatField,
    dtdt_longwave: FloatField,
    dtdt_from_dynamics: FloatField,
    dvapordt_from_dynamics: FloatField,
    dtdt_pbl: FloatField,
    dspecific_humiditydt_pbl: FloatField,
    grid_scale_forcing_t: FloatField,
    grid_scale_forcing_vapor: FloatField,
    subgrid_scale_forcing_t: FloatField,
    subgrid_scale_forcing_vapor: FloatField,
    advective_forcing_t: FloatField,
    convective_liquid_local: FloatField,
    convective_ice_local: FloatField,
    convective_cloud_fraction_local: FloatField,
    large_scale_liquid_local: FloatField,
    large_scale_ice_local: FloatField,
    large_scale_cloud_fraction_local: FloatField,
    convection_tracer: FloatField,
    buoyancy_excess: FloatField,
):
    """
    Get the desired state for the convection scheme, controlled by external GF_ENV_SETTING.

    GF_ENV_SETTING = 0:
        use state updated by during current timestep (dynamics and physics)
    GF_ENV_SETTING = 1:
        use state from start of timestep (pre-dynamics)
    """
    from __externals__ import k_end, GF_ENV_SETTING, ENTRVERSION, CONVECTION_TRACER

    with computation(PARALLEL), interval(...):
        if GF_ENV_SETTING == 0:
            dz = -(geopotential_height_interface[0, 0, 1] - geopotential_height_interface)
            air_density = p / (287.04 * t * (1.0 + 0.608 * vapor))
            t_local = t.at(K=k_end - K)
            p_local = p.at(K=k_end - K)
            vapor_local = vapor.at(K=k_end - K)
            vapor_current_local = vapor.at(K=k_end - K)
            u_local = u.at(K=k_end - K)
            v_local = v.at(K=k_end - K)
            vertical_velocity_local = vertical_velocity.at(K=k_end - K)
            layer_height_above_surface_local = layer_height_above_surface.at(K=k_end - K)
            edge_height_above_surface_local = edge_height_above_surface.at(K=k_end - K)
            mass_local = mass.at(K=k_end - K)
            scalar_diffusivity_local = scalar_diffusivity.at(K=k_end - K)

            # Grid and sub-grid scale forcings for convection
            grid_scale_forcing_t = 0.0
            grid_scale_forcing_vapor = 0.0
            subgrid_scale_forcing_t = 0.0
            subgrid_scale_forcing_vapor = 0.0
            advective_forcing_t = 0.0

    with computation(FORWARD), interval(0, 1):
        if GF_ENV_SETTING == 0:
            p_surface = p_interface.at(K=k_end + 1)

    with computation(PARALLEL), interval(...):
        if GF_ENV_SETTING == 1:
            mass_n = (p_interface_timestep_start[0, 0, 1] - p_interface_timestep_start) * (
                1.0 / constants.MAPL_GRAV
            )
            p_n = 0.5 * (p_interface_timestep_start + p_interface_timestep_start[0, 0, 1])
            p_kappa_interface_n = (p_interface_timestep_start / constants.MAPL_P00) ** (
                constants.MAPL_RGAS / constants.MAPL_CP
            )
            if K == k_end:
                p_kappa_surface_n = (p_interface_timestep_start[0, 0, 1] / constants.MAPL_P00) ** (
                    constants.MAPL_RGAS / constants.MAPL_CP
                )
            p_kappa_n = (p_n / constants.MAPL_P00) ** (constants.MAPL_RGAS / constants.MAPL_CP)
            edge_height_above_surface_n = (t_timestep_start / p_kappa_n) * (
                1.0 + constants.MAPL_VIREPS * vapor_timestep_start
            )

    with computation(BACKWARD), interval(...):
        if GF_ENV_SETTING == 1:
            if K == k_end:
                layer_height_above_surface_n = (
                    0
                    + (constants.MAPL_CP / constants.MAPL_GRAV)
                    * (p_kappa_surface_n - p_kappa_n)
                    * edge_height_above_surface_n
                )
            else:
                layer_height_above_surface_n = (
                    edge_height_above_surface_n[0, 0, 1]
                    + (constants.MAPL_CP / constants.MAPL_GRAV)
                    * (p_kappa_interface_n[0, 0, 1] - p_kappa_n)
                    * edge_height_above_surface_n
                )
            edge_height_above_surface_n = (
                layer_height_above_surface_n
                + (constants.MAPL_CP / constants.MAPL_GRAV)
                * (p_kappa_n - p_kappa_interface_n)
                * edge_height_above_surface_n
            )

    with computation(PARALLEL), interval(...):
        if GF_ENV_SETTING == 1:
            if K == k_end:
                dz = -(0 - edge_height_above_surface_n)
            else:
                dz = -(edge_height_above_surface_n[0, 0, 1] - edge_height_above_surface_n)
            air_density = p_n / (287.04 * t_timestep_start * (1.0 + 0.608 * vapor_timestep_start))
            t_local = t_timestep_start.at(K=k_end - K)
            p_local = p_n.at(K=k_end - K)
            vapor_local = vapor_timestep_start.at(K=k_end - K)
            vapor_current_local = vapor.at(K=k_end - K)
            u_local = u_timestep_start.at(K=k_end - K)
            v_local = v_timestep_start.at(K=k_end - K)
            vertical_velocity_local = vertical_velocity.at(K=k_end - K)
            layer_height_above_surface_local = layer_height_above_surface_n.at(K=k_end - K)
            if K == 0:
                edge_height_above_surface_local = 0
            else:
                edge_height_above_surface_local = edge_height_above_surface_n.at(K=k_end - K + 1)
            mass_local = mass_n.at(K=k_end - K)
            scalar_diffusivity_local = scalar_diffusivity.at(K=k_end - K)

            # Grid and sub-grid scale forcings for convection
            grid_scale_forcing_t = (
                dtdt_from_dynamics.at(K=k_end - K)
                + dtdt_shortwave.at(K=k_end - K)
                + dtdt_longwave.at(K=k_end - K)
            )
            grid_scale_forcing_vapor = dvapordt_from_dynamics.at(K=k_end - K)
            subgrid_scale_forcing_t = dtdt_pbl.at(K=k_end - K)
            subgrid_scale_forcing_vapor = dspecific_humiditydt_pbl.at(K=k_end - K)
            advective_forcing_t = dtdt_from_dynamics.at(K=k_end - K)

    with computation(FORWARD), interval(0, 1):
        if GF_ENV_SETTING == 1:
            p_surface = p_interface_timestep_start.at(K=k_end + 1)

    with computation(PARALLEL), interval(...):
        # remainder is the same for both settings
        convective_liquid_local = convective_liquid.at(K=k_end - K)
        convective_ice_local = convective_ice.at(K=k_end - K)
        convective_cloud_fraction_local = convective_cloud_fraction.at(K=k_end - K)
        large_scale_liquid_local = large_scale_liquid.at(K=k_end - K)
        large_scale_ice_local = large_scale_ice.at(K=k_end - K)
        large_scale_cloud_fraction_local = large_scale_cloud_fraction.at(K=k_end - K)

        if ENTRVERSION == 0:
            # eq 6 of https://doi.org/10.1029/2021JD034881
            lateral_entrainment_rate_local = (
                0.71
                * max(0.5, vertical_velocity.at(K=k_end - K)) ** (-1.17)
                * max(0.1, buoyancy.at(K=k_end - K)) ** (-0.36)
            )
        else:
            lateral_entrainment_rate_local = 1.0

    with computation(PARALLEL), interval(...):
        # must be in separate computation to ensure lateral_entrainment_rate_local is written before read
        lateral_entrainment_rate = lateral_entrainment_rate_local.at(K=k_end - K)
        if CONVECTION_TRACER == 1:
            buoyancy_excess = convection_tracer.at(K=k_end - K)
        else:
            buoyancy_excess = 0.0


def prepare_cumulus_parameterization(
    t_excess: FloatFieldIJ,
    vapor_excess: FloatFieldIJ,
    last_ierr: FloatFieldIJ,
    fix_out_vapor: FloatFieldIJ,
    conprr: FloatFieldIJ,
    evap_subl_tendency_cu_param: FloatField,
    convective_precip_flux_cu_param: FloatField,
    t_perturbation_cu_param_horizontal: FloatField,
    t_perturbation_cu_param_vertical: FloatField,
    t_perturbation_cu_param: FloatField,
    omega_cu_param: FloatField,
    ccn: FloatFieldIJ,
    dtdt_cu_param_shallow: FloatField,
    dtdt_cu_param_mid: FloatField,
    dtdt_cu_param_deep: FloatField,
    dudt_cu_param_shallow: FloatField,
    dudt_cu_param_mid: FloatField,
    dudt_cu_param_deep: FloatField,
    dvdt_cu_param_shallow: FloatField,
    dvdt_cu_param_mid: FloatField,
    dvdt_cu_param_deep: FloatField,
    dvapordt_cu_param_shallow: FloatField,
    dvapordt_cu_param_mid: FloatField,
    dvapordt_cu_param_deep: FloatField,
    dvapordt_cu_param_combined: FloatField,
    dcloudicedt_cu_param_shallow: FloatField,
    dcloudicedt_cu_param_mid: FloatField,
    dcloudicedt_cu_param_deep: FloatField,
    dnicedt_cu_param_shallow: FloatField,
    dnicedt_cu_param_mid: FloatField,
    dnicedt_cu_param_deep: FloatField,
    dnliquiddt_cu_param_shallow: FloatField,
    dnliquiddt_cu_param_mid: FloatField,
    dnliquiddt_cu_param_deep: FloatField,
    dbuoyancydt_cu_param_shallow: FloatField,
    dbuoyancydt_cu_param_mid: FloatField,
    dbuoyancydt_cu_param_deep: FloatField,
    dconvectiveicedt_cu_param_shallow: FloatField,
    dconvectiveicedt_cu_param_mid: FloatField,
    dconvectiveicedt_cu_param_deep: FloatField,
    dlargescaleicedt_cu_param_shallow: FloatField,
    dlargescaleicedt_cu_param_mid: FloatField,
    dlargescaleicedt_cu_param_deep: FloatField,
    dconvectiveliquiddt_cu_param_shallow: FloatField,
    dconvectiveliquiddt_cu_param_mid: FloatField,
    dconvectiveliquiddt_cu_param_deep: FloatField,
    dlargescaleliquiddt_cu_param_shallow: FloatField,
    dlargescaleliquiddt_cu_param_mid: FloatField,
    dlargescaleliquiddt_cu_param_deep: FloatField,
    dconvectivecloudfractiondt_cu_param_shallow: FloatField,
    dconvectivecloudfractiondt_cu_param_mid: FloatField,
    dconvectivecloudfractiondt_cu_param_deep: FloatField,
    dlargescalecloudfractiondt_cu_param_shallow: FloatField,
    dlargescalecloudfractiondt_cu_param_mid: FloatField,
    dlargescalecloudfractiondt_cu_param_deep: FloatField,
    p_surface: FloatFieldIJ,
    topography_height: FloatFieldIJ,
    topography_height_no_negative: FloatFieldIJ,
    latitude: FloatFieldIJ,
    longitude: FloatFieldIJ,
    latitude_degrees: FloatFieldIJ,
    longitude_degrees: FloatFieldIJ,
):
    from __externals__ import APPLY_SUB_MP, ADV_TRIGGER, USE_TRACER_TRANSP, AUTOCONV

    with computation(FORWARD), interval(0, 1):
        # pre-fill some fields
        t_excess = 0.0
        vapor_excess = 0.0
        last_ierr = -999
        fix_out_vapor = 1.0
        conprr = 0.0

    with computation(PARALLEL), interval(...):
        # pre-fill some fields
        evap_subl_tendency_cu_param = 0.0
        convective_precip_flux_cu_param = 0.0
        t_perturbation_cu_param = 0.0

        dtdt_cu_param_shallow = 0.0
        dtdt_cu_param_mid = 0.0
        dtdt_cu_param_deep = 0.0
        dudt_cu_param_shallow = 0.0
        dudt_cu_param_mid = 0.0
        dudt_cu_param_deep = 0.0
        dvdt_cu_param_shallow = 0.0
        dvdt_cu_param_mid = 0.0
        dvdt_cu_param_deep = 0.0
        dvapordt_cu_param_shallow = 0.0
        dvapordt_cu_param_mid = 0.0
        dvapordt_cu_param_deep = 0.0
        dvapordt_cu_param_combined = 0.0
        dcloudicedt_cu_param_shallow = 0.0
        dcloudicedt_cu_param_mid = 0.0
        dcloudicedt_cu_param_deep = 0.0
        dnicedt_cu_param_shallow = 0.0
        dnicedt_cu_param_mid = 0.0
        dnicedt_cu_param_deep = 0.0
        dnliquiddt_cu_param_shallow = 0.0
        dnliquiddt_cu_param_mid = 0.0
        dnliquiddt_cu_param_deep = 0.0
        dbuoyancydt_cu_param_shallow = 0.0
        dbuoyancydt_cu_param_mid = 0.0
        dbuoyancydt_cu_param_deep = 0.0

        if APPLY_SUB_MP == 1:
            dconvectiveicedt_cu_param_shallow = 0.0
            dconvectiveicedt_cu_param_mid = 0.0
            dconvectiveicedt_cu_param_deep = 0.0
            dlargescaleicedt_cu_param_shallow = 0.0
            dlargescaleicedt_cu_param_mid = 0.0
            dlargescaleicedt_cu_param_deep = 0.0
            dconvectiveliquiddt_cu_param_shallow = 0.0
            dconvectiveliquiddt_cu_param_mid = 0.0
            dconvectiveliquiddt_cu_param_deep = 0.0
            dlargescaleliquiddt_cu_param_shallow = 0.0
            dlargescaleliquiddt_cu_param_mid = 0.0
            dlargescaleliquiddt_cu_param_deep = 0.0
            dconvectivecloudfractiondt_cu_param_shallow = 0.0
            dconvectivecloudfractiondt_cu_param_mid = 0.0
            dconvectivecloudfractiondt_cu_param_deep = 0.0
            dlargescalecloudfractiondt_cu_param_shallow = 0.0
            dlargescalecloudfractiondt_cu_param_mid = 0.0
            dlargescalecloudfractiondt_cu_param_deep = 0.0

        omega_cu_param = 0.0

        if ADV_TRIGGER == 2:
            t_perturbation_cu_param = t_perturbation_cu_param_horizontal + t_perturbation_cu_param_vertical

        if USE_TRACER_TRANSP == 1:
            tracer_stuff_not_implemented_yet = 0

    with computation(FORWARD), interval(0, 1):
        if AUTOCONV == 2:
            ccn = max(100.0, (370.37 * (0.01 + max(0.0, 0.1))) ** 1.555)
        else:
            ccn = 100.0

        p_surface = p_surface * 1.0e-2

        topography_height_no_negative = max(0, topography_height)
        latitude_degrees = latitude * 180.0 / 3.14159
        longitude_degrees = longitude * 180.0 / 3.14159


class GF2020Setup:
    def __init__(
        self, stencil_factory: StencilFactory, quantity_factory: QuantityFactory, GF_2020_config: GF2020Config
    ):
        self.GF_2020_config = GF_2020_config

        self.plume_order = "SH_MD_DP"
        ndsl_log.warning("plume order currently set manually, need to integrate this into config")

        # Construct stencils
        self._compute_extra_inputs_from_state = stencil_factory.from_dims_halo(
            func=compute_extra_inputs_from_state,
            compute_dims=[X_DIM, Y_DIM, Z_INTERFACE_DIM],
            externals={
                "STOCHASTIC_CONVECTION": GF_2020_config.STOCHASTIC_CNV,
                "STOCH_TOP": GF_2020_config.STOCH_TOP,
                "STOCH_BOT": GF_2020_config.STOCH_BOT,
                "GF_MIN_AREA": GF_2020_config.GF_MIN_AREA,
                "LHYDROSTATIC": GF_2020_config.LHYDROSTATIC,
            },
        )

        self._zero_state_fields = stencil_factory.from_dims_halo(
            func=zero_state_fields,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._compute_2d_inputs = stencil_factory.from_dims_halo(
            func=compute_2d_inputs,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "SINGLE_COLUMN_MODE": GF_2020_config.SINGLE_COLUMN_MODE,
            },
        )

        self._set_local_state = stencil_factory.from_dims_halo(
            func=set_local_state,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "GF_ENV_SETTING": GF_2020_config.GF_ENV_SETTING,
                "ENTRVERSION": GF_2020_config.ENTRVERSION,
                "CONVECTION_TRACER": GF_2020_config.CONVECTION_TRACER,
            },
        )

        self._prepare_cumulus_parameterization = stencil_factory.from_dims_halo(
            func=prepare_cumulus_parameterization,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
            externals={
                "APPLY_SUB_MP": GF_2020_config.APPLY_SUB_MP,
                "ADV_TRIGGER": GF_2020_config.ADV_TRIGGER,
                "USE_TRACER_TRANSP": GF_2020_config.USE_TRACER_TRANSP,
                "AUTOCONV": GF_2020_config.AUTOCONV,
            },
        )

    def __call__(
        self,
        state: GF2020State,
        saturation_tables: SaturationVaporPressureTable,
        temporaries: GF2020Temporaries,
    ):
        self._compute_extra_inputs_from_state(
            p_interface=state.p_interface,
            p=temporaries.p,
            p_kappa=temporaries.p_kappa,
            edge_height_above_surface=temporaries.edge_height_above_surface,
            layer_height_above_surface=temporaries.layer_height_above_surface,
            geopotential_height_interface=state.geopotential_height_interface,
            t=state.t,
            th=temporaries.th,
            vapor=state.vapor,
            mass=temporaries.mass,
            w=state.w,
            omega=state.omega,
            vertical_motion=temporaries.vertical_velocity,
            tpwi=state.total_precipitable_water_initial,
            tpwi_star=state.saturation_total_precipitable_water_initial,
            seed_convection=temporaries.seed_convection,
            area=state.area,
            modified_area=temporaries.modified_area,
            convection_fraction=state.convection_fraction,
            ese=saturation_tables.ese,
            esx=saturation_tables.esx,
        )

        self._zero_state_fields(
            dvapordt_deep_convection=state.dvapordt_deep_convection,
            dtdt_deep_convection=state.dtdt_deep_convection,
            dudt_deep_convection=state.dudt_deep_convection,
            dvdt_deep_convection=state.dvdt_deep_convection,
            sigma_deep=state.sigma_deep,
            sigma_mid=state.sigma_mid,
            mass_flux_shalow=state.mass_flux_shalow,
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
            cloud_work_function_0=state.cloud_work_function_0,
            cloud_work_function_1=state.cloud_work_function_1,
            cloud_work_function_2=state.cloud_work_function_2,
            cloud_work_function_3=state.cloud_work_function_3,
            cloud_work_function_1_pbl=state.cloud_work_function_1_pbl,
            cloud_work_function_1_cin=state.cloud_work_function_1_cin,
            convective_precipitation_RAS=state.convective_precipitation_RAS,
            convective_precipitation_GF=state.convective_precipitation_GF,
            convective_condensate_source=state.convective_condensate_source,
            convective_condensate_grid_mean=state.convective_condensate_grid_mean,
            total_water_flux_deep_convection=state.total_water_flux_deep_convection,
            updraft_area_fraction=state.updraft_area_fraction,
            updraft_vertical_velocity=state.updraft_vertical_velocity,
            lighting_density=state.lighting_density,
            entrainment_parameter=state.entrainment_parameter,
            pbl_time_scale=state.pbl_time_scale,
            cape_removal_time_scale=state.cape_removal_time_scale,
        )

        # workaround because max of full field cannot be determined inside a stencil
        t_2m_max = state.t_2m.max()

        # if surface temperature is not yet set in single column mode, do not run
        if self.GF_2020_config.SINGLE_COLUMN_MODE == True and t_2m_max < 1.0e-6:
            return

        self._compute_2d_inputs(
            t=state.t,
            t_2m_max=t_2m_max,
            t_2m=state.t_2m,
            t_2m_local=temporaries.t_2m_local,
            evaporation=state.evaporation,
            evaporation_local=temporaries.evaporation_local,
            sensible_heat_flux=state.sensible_heat_flux,
            p_interface=state.p_interface,
            vapor=state.vapor,
            geopotential_height_surface=state.geopotential_height_surface,
            topography_height=temporaries.topography_height,
            land_fraction=state.land_fraction,
            ocean_fraction=temporaries.ocean_fraction,
            area=state.area,
            grid_length=temporaries.grid_length,
            pbl_level=state.pbl_level,
            pbl_level_local=temporaries.pbl_level_local,
        )

        self._set_local_state(
            geopotential_height_interface=state.geopotential_height_interface,
            t=state.t,
            t_local=temporaries.t_local,
            p=temporaries.p,
            p_local=temporaries.p_local,
            vapor=state.vapor,
            vapor_local=temporaries.vapor_local,
            vapor_current_local=temporaries.vapor_current_local,
            u=state.u,
            u_local=temporaries.u_local,
            v=state.v,
            v_local=temporaries.v_local,
            vertical_velocity=temporaries.vertical_velocity,
            vertical_velocity_local=temporaries.vertical_velocity_local,
            layer_height_above_surface=temporaries.layer_height_above_surface,
            layer_height_above_surface_local=temporaries.layer_height_above_surface_local,
            edge_height_above_surface=temporaries.edge_height_above_surface,
            edge_height_above_surface_local=temporaries.edge_height_above_surface_local,
            mass=temporaries.mass,
            mass_local=temporaries.mass_local,
            scalar_diffusivity=state.scalar_diffusivity,
            scalar_diffusivity_local=temporaries.scalar_diffusivity_local,
            lateral_entrainment_rate_local=temporaries.lateral_entrainment_rate_local,
            lateral_entrainment_rate=state.lateral_entrainment_rate,
            p_surface=temporaries.p_surface,
            buoyancy=state.buoyancy,
            p_interface=state.p_interface,
            convective_liquid=state.convective_liquid,
            convective_ice=state.convective_ice,
            convective_cloud_fraction=state.convective_cloud_fraction,
            large_scale_liquid=state.large_scale_liquid,
            large_scale_ice=state.large_scale_ice,
            large_scale_cloud_fraction=state.large_scale_cloud_fraction,
            p_interface_timestep_start=state.p_interface_timestep_start,
            t_timestep_start=state.t_timestep_start,
            u_timestep_start=state.u_timestep_start,
            v_timestep_start=state.v_timestep_start,
            vapor_timestep_start=state.vapor_timestep_start,
            dtdt_shortwave=state.dtdt_shortwave,
            dtdt_longwave=state.dtdt_longwave,
            dtdt_from_dynamics=state.dtdt_from_dynamics,
            dvapordt_from_dynamics=state.dvapordt_from_dynamics,
            dtdt_pbl=state.dtdt_pbl,
            dspecific_humiditydt_pbl=state.dspecific_humiditydt_pbl,
            grid_scale_forcing_t=temporaries.grid_scale_forcing_t,
            grid_scale_forcing_vapor=temporaries.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t=temporaries.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor=temporaries.subgrid_scale_forcing_vapor,
            advective_forcing_t=temporaries.advective_forcing_t,
            convective_liquid_local=temporaries.convective_liquid_local,
            convective_ice_local=temporaries.convective_ice_local,
            convective_cloud_fraction_local=temporaries.convective_cloud_fraction_local,
            large_scale_liquid_local=temporaries.large_scale_liquid_local,
            large_scale_ice_local=temporaries.large_scale_ice_local,
            large_scale_cloud_fraction_local=temporaries.large_scale_cloud_fraction_local,
            convection_tracer=state.convection_tracer,
            buoyancy_excess=temporaries.buoyancy_excess,
        )

        if self.GF_2020_config.ADV_TRIGGER == 2:
            raise NotImplementedError("ADV_TRIGGER == 2 not implemented yet")

        self._prepare_cumulus_parameterization(
            t_excess=temporaries.t_excess,
            vapor_excess=temporaries.vapor_excess,
            last_ierr=temporaries.last_ierr,
            fix_out_vapor=temporaries.fix_out_vapor,
            conprr=temporaries.conprr,
            evap_subl_tendency_cu_param=temporaries.evap_subl_tendency_cu_param,
            convective_precip_flux_cu_param=temporaries.convective_precip_flux_cu_param,
            t_perturbation_cu_param_horizontal=temporaries.t_perturbation_cu_param_horizontal,
            t_perturbation_cu_param_vertical=temporaries.t_perturbation_cu_param_vertical,
            t_perturbation_cu_param=temporaries.t_perturbation_cu_param,
            omega_cu_param=temporaries.omega_cu_param,
            ccn=temporaries.ccn,
            dtdt_cu_param_shallow=temporaries.dtdt_cu_param_shallow,
            dtdt_cu_param_mid=temporaries.dtdt_cu_param_mid,
            dtdt_cu_param_deep=temporaries.dtdt_cu_param_deep,
            dudt_cu_param_shallow=temporaries.dudt_cu_param_shallow,
            dudt_cu_param_mid=temporaries.dudt_cu_param_mid,
            dudt_cu_param_deep=temporaries.dudt_cu_param_deep,
            dvdt_cu_param_shallow=temporaries.dvdt_cu_param_shallow,
            dvdt_cu_param_mid=temporaries.dvdt_cu_param_mid,
            dvdt_cu_param_deep=temporaries.dvdt_cu_param_deep,
            dvapordt_cu_param_shallow=temporaries.dvapordt_cu_param_shallow,
            dvapordt_cu_param_mid=temporaries.dvapordt_cu_param_mid,
            dvapordt_cu_param_deep=temporaries.dvapordt_cu_param_deep,
            dvapordt_cu_param_combined=temporaries.dvapordt_cu_param_combined,
            dcloudicedt_cu_param_shallow=temporaries.dcloudicedt_cu_param_shallow,
            dcloudicedt_cu_param_mid=temporaries.dcloudicedt_cu_param_mid,
            dcloudicedt_cu_param_deep=temporaries.dcloudicedt_cu_param_deep,
            dnicedt_cu_param_shallow=temporaries.dnicedt_cu_param_shallow,
            dnicedt_cu_param_mid=temporaries.dnicedt_cu_param_mid,
            dnicedt_cu_param_deep=temporaries.dnicedt_cu_param_deep,
            dnliquiddt_cu_param_shallow=temporaries.dnliquiddt_cu_param_shallow,
            dnliquiddt_cu_param_mid=temporaries.dnliquiddt_cu_param_mid,
            dnliquiddt_cu_param_deep=temporaries.dnliquiddt_cu_param_deep,
            dbuoyancydt_cu_param_shallow=temporaries.dbuoyancydt_cu_param_shallow,
            dbuoyancydt_cu_param_mid=temporaries.dbuoyancydt_cu_param_mid,
            dbuoyancydt_cu_param_deep=temporaries.dbuoyancydt_cu_param_deep,
            dconvectiveicedt_cu_param_shallow=temporaries.dconvectiveicedt_cu_param_shallow,
            dconvectiveicedt_cu_param_mid=temporaries.dconvectiveicedt_cu_param_mid,
            dconvectiveicedt_cu_param_deep=temporaries.dconvectiveicedt_cu_param_deep,
            dlargescaleicedt_cu_param_shallow=temporaries.dlargescaleicedt_cu_param_shallow,
            dlargescaleicedt_cu_param_mid=temporaries.dlargescaleicedt_cu_param_mid,
            dlargescaleicedt_cu_param_deep=temporaries.dlargescaleicedt_cu_param_deep,
            dconvectiveliquiddt_cu_param_shallow=temporaries.dconvectiveliquiddt_cu_param_shallow,
            dconvectiveliquiddt_cu_param_mid=temporaries.dconvectiveliquiddt_cu_param_mid,
            dconvectiveliquiddt_cu_param_deep=temporaries.dconvectiveliquiddt_cu_param_deep,
            dlargescaleliquiddt_cu_param_shallow=temporaries.dlargescaleliquiddt_cu_param_shallow,
            dlargescaleliquiddt_cu_param_mid=temporaries.dlargescaleliquiddt_cu_param_mid,
            dlargescaleliquiddt_cu_param_deep=temporaries.dlargescaleliquiddt_cu_param_deep,
            dconvectivecloudfractiondt_cu_param_shallow=temporaries.dconvectivecloudfractiondt_cu_param_shallow,
            dconvectivecloudfractiondt_cu_param_mid=temporaries.dconvectivecloudfractiondt_cu_param_mid,
            dconvectivecloudfractiondt_cu_param_deep=temporaries.dconvectivecloudfractiondt_cu_param_deep,
            dlargescalecloudfractiondt_cu_param_shallow=temporaries.dlargescalecloudfractiondt_cu_param_shallow,
            dlargescalecloudfractiondt_cu_param_mid=temporaries.dlargescalecloudfractiondt_cu_param_mid,
            dlargescalecloudfractiondt_cu_param_deep=temporaries.dlargescalecloudfractiondt_cu_param_deep,
            p_surface=temporaries.p_surface,
            topography_height=temporaries.topography_height,
            topography_height_no_negative=temporaries.topography_height_no_negative,
            latitude=state.latitude,
            longitude=state.longitude,
            latitude_degrees=temporaries.latitude_degrees,
            longitude_degrees=temporaries.longitude_degrees,
        )
