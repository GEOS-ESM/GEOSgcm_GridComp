import copy
from ndsl import StencilFactory, QuantityFactory, ndsl_log
from ndsl.dsl.gt4py import PARALLEL, interval, computation, FORWARD, sqrt, max, min, abs, floor, BACKWARD
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Float, K, IntFieldIJ
import pyMoist.constants as constants
from pyMoist.saturation_tables.types import GlobalTable_saturation_tables
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.saturation_tables.tables.main import SaturationVaporPressureTable
from pyMoist.convection.GF_2020.config import GF2020Config
from pyMoist.convection.GF_2020.state import GF2020State
from pyMoist.convection.GF_2020.locals import GF2020Locals


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


def set_2d_fields(
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
    pbl_level_cu_param_input: IntFieldIJ,
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

        if pbl_level != 0.0:
            pbl_level_cu_param_input = k_end - int(round(pbl_level))
        else:
            pbl_level_cu_param_input = 0


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
    lateral_entrainment_rate: FloatField,
    lateral_entrainment_rate_local: FloatField,
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
    t_perturbation_horizontal: FloatField,
    t_perturbation_vertical: FloatField,
    t_perturbation: FloatField,
    ccn: FloatFieldIJ,
    p_surface: FloatFieldIJ,
    topography_height: FloatFieldIJ,
    topography_height_no_negative: FloatFieldIJ,
    latitude: FloatFieldIJ,
    longitude: FloatFieldIJ,
    latitude_degrees: FloatFieldIJ,
    longitude_degrees: FloatFieldIJ,
    layer_height_above_surface: FloatField,
    geopotential_height_cu_param_input: FloatField,
    p: FloatField,
    p_mb: FloatField,
    t_local: FloatField,
    t_cu_param_input: FloatField,
    vapor_local: FloatField,
    vapor_current_local: FloatField,
    vapor_timestep_start_cu_param_input: FloatField,
    vapor_current_cu_param_input: FloatField,
    air_density_cu_param_input: FloatField,
    u_local: FloatField,
    v_local: FloatField,
    w_local: FloatField,
    u_cu_param_input: FloatField,
    v_cu_param_input: FloatField,
    w_cu_param_input: FloatField,
    omega_cu_param_input: FloatField,
    mass_local: FloatField,
    mass_cu_param_input: FloatField,
    advective_forcing_t: FloatField,
    t_modified_by_advection: FloatField,
    grid_scale_forcing_vapor: FloatField,
    vapor_modified_by_advection: FloatField,
    pbl_level_cu_param_input: IntFieldIJ,
    pbl_height_cu_param_input: FloatFieldIJ,
    sensible_heat_flux_local: FloatFieldIJ,
    sensible_heat_flux_cu_param_input: FloatFieldIJ,
    evaporation_local: FloatFieldIJ,
    latent_heat_flux_cu_param_input: FloatFieldIJ,
    convective_scale_velosity_cu_param_input: FloatFieldIJ,
    t_excess_cu_param_input: FloatFieldIJ,
    vapor_excess_cu_param_input: FloatFieldIJ,
):
    from __externals__ import ADV_TRIGGER, USE_TRACER_TRANSP, AUTOCONV, DT_MOIST

    with computation(PARALLEL), interval(...):
        if ADV_TRIGGER == 2:
            t_perturbation = t_perturbation_horizontal + t_perturbation_vertical

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

    with computation(PARALLEL), interval(0, -1):
        # NOTE lots of redundancies/unnecessary calculations from here down which can be removed
        # These were maintained during porting to retain as much consistency with the Fortran as possible.

        geopotential_height_cu_param_input = layer_height_above_surface + topography_height
        p_mb = p * 1.0e-2
        t_cu_param_input = t_local
        vapor_timestep_start_cu_param_input = vapor_local
        vapor_current_cu_param_input = vapor_current_local
        air_density_cu_param_input = (
            1.0e2 * p_mb / (287.04 * t_cu_param_input * (1.0 + 0.608 * vapor_timestep_start_cu_param_input))
        )

        u_cu_param_input = u_local
        v_cu_param_input = v_local
        w_cu_param_input = w_local
        omega_cu_param_input = -1 * constants.MAPL_GRAV * air_density_cu_param_input * w_local

        mass_cu_param_input = mass_local

        t_modified_by_advection = t_cu_param_input + advective_forcing_t * DT_MOIST
        vapor_modified_by_advection = (
            vapor_timestep_start_cu_param_input + grid_scale_forcing_vapor * DT_MOIST
        )

    with computation(FORWARD), interval(0, 1):
        pbl_height_cu_param_input = (
            geopotential_height_cu_param_input.at(K=pbl_level_cu_param_input) - topography_height
        )

    with computation(FORWARD), interval(0, 1):
        # get execess T and Q for source air parcels
        density = (
            100.0 * p_surface / (287.04 * (t_cu_param_input * (1.0 + 0.608 * vapor_current_cu_param_input)))
        )
        # sensible and latent sfc fluxes for the heat-engine closure
        sensible_heat_flux_cu_param_input = density * constants.MAPL_CP * sensible_heat_flux_local  # W/m^2
        latent_heat_flux_cu_param_input = density * constants.MAPL_ALHL * evaporation_local  # W/m^2

        buoyancy_flux = (
            -sensible_heat_flux_cu_param_input * density * 1004.64 / 1004.64
            + 0.608 * t_cu_param_input * latent_heat_flux_cu_param_input
        ) / density  # K m s-1

        depth_of_first_model_layer = (
            2.0 * (geopotential_height_cu_param_input.at(K=0) - topography_height) * constants.MAPL_GRAV
        )  # m+2 s-2

        convective_scale_velosity_cu_param_input = max(
            0.0, 0.001 - 1.5 * 0.41 * buoyancy_flux * depth_of_first_model_layer / t_cu_param_input
        )  # m+3 s-3

        if convective_scale_velosity_cu_param_input > constants.FLOAT_TINY:
            # convective-scale velocity w*
            convective_scale_velosity_cu_param_input = 1.2 * convective_scale_velosity_cu_param_input**0.3333
            # temperature excess
            t_excess_cu_param_input = max(
                0.0,
                -1.5
                * -sensible_heat_flux_cu_param_input
                * density
                * 1004.64
                / (density * convective_scale_velosity_cu_param_input * 1004.64),
            )  # K
            # moisture  excess
            vapor_excess_cu_param_input = max(
                0.0,
                -1.5 * latent_heat_flux_cu_param_input / (density * convective_scale_velosity_cu_param_input),
            )  # kg kg-1

        # convective scale velosity for shallow convection closure (Grant 2001)
        # depth of the pbl
        # convective-scale velocity W* (m/s)
        convective_scale_velosity_cu_param_input = max(
            0.0,
            0.001
            - 1.5 * 0.41 * buoyancy_flux * pbl_height_cu_param_input * constants.MAPL_GRAV / t_cu_param_input,
        )
        convective_scale_velosity_cu_param_input = 1.2 * convective_scale_velosity_cu_param_input**0.3333


def prepare_cumulus_paramaterization_microphysics(
    convective_liquid_local: FloatField,
    convective_ice_local: FloatField,
    convective_cloud_fraction_local: FloatField,
    large_scale_liquid_local: FloatField,
    large_scale_ice_local: FloatField,
    large_scale_cloud_fraction_local: FloatField,
    convective_liquid_cu_param_input: FloatField,
    convective_ice_cu_param_input: FloatField,
    convective_cloud_fraction_cu_param_input: FloatField,
    large_scale_liquid_cu_param_input: FloatField,
    large_scale_ice_cu_param_input: FloatField,
    large_scale_cloud_fraction_cu_param_input: FloatField,
):
    with computation(PARALLEL), interval(0, -1):
        convective_liquid_cu_param_input = convective_liquid_local
        convective_ice_cu_param_input = convective_ice_local
        convective_cloud_fraction_cu_param_input = convective_cloud_fraction_local
        large_scale_liquid_cu_param_input = large_scale_liquid_local
        large_scale_ice_cu_param_input = large_scale_ice_local
        large_scale_cloud_fraction_cu_param_input = large_scale_cloud_fraction_local


class GF2020Setup:
    """
    This class performs the entire setup sequence for the GF2020 convection parameterization scheme, based
    on the Fortran code available in GEOS v11.4.2.

    In GEOS v11.4.2, this code is split across three subroutines nested as follows:

    - GF_Run
        - GF2020_INTERFACE
            - GF2020_DRV

    This python implementation simplifies this structure by bringing all setup calculations to the same level,
    reducing some (but likely not all) redundent/duplicate field definitions in the process.

    The call function for this class reads a takes a NDSL state ("state") filled with model state as input and
    returns a NDSL state ("locals") filled with computed quantities prepared to be passed to the cumulus
    parameterization core.
    """

    def __init__(
        self, stencil_factory: StencilFactory, quantity_factory: QuantityFactory, GF_2020_config: GF2020Config
    ):
        """
        Build stencils
        """
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

        self._set_2d_fields = stencil_factory.from_dims_halo(
            func=set_2d_fields,
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
                "ADV_TRIGGER": GF_2020_config.ADV_TRIGGER,
                "USE_TRACER_TRANSP": GF_2020_config.USE_TRACER_TRANSP,
                "AUTOCONV": GF_2020_config.AUTOCONV,
                "DT_MOIST": GF_2020_config.DT_MOIST,
            },
        )

        self._prepare_cumulus_paramaterization_microphysics = stencil_factory.from_dims_halo(
            func=prepare_cumulus_paramaterization_microphysics,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        state: GF2020State,
        saturation_tables: SaturationVaporPressureTable,
        locals: GF2020Locals,
    ):
        """
        Perform setup calculations

        Args:
            state: NDSL data class containing all fields from overarching model used at some point in GF2020
            saturation_tables: saturation vapor pressure tables, for liquid, ice, and dynamic surfaces
            locals: all internal GF2020 fields

        """
        # TODO reset all temporaries to zero
        # TODO reset last_ierr = -999, fix_out_vapor = 1
        locals.miscelaneous_diagnostic.last_ierr.field[:] = -999
        locals.miscelaneous_diagnostic.fix_out_vapor.field[:] = 1

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

        self._set_2d_fields(
            t=state.t,
            t_2m_max=t_2m_max,
            t_2m=state.t_2m,
            t_2m_local=locals.local_copy.t_2m,
            evaporation=state.evaporation,
            evaporation_local=locals.local_copy.evaporation,
            sensible_heat_flux=state.sensible_heat_flux,
            p_interface=state.p_interface,
            vapor=state.vapor,
            geopotential_height_surface=state.geopotential_height_surface,
            topography_height=locals.derived_state.topography_height,
            land_fraction=state.land_fraction,
            ocean_fraction=locals.cumulus_parameterization_input.ocean_fraction,
            area=state.area,
            grid_length=locals.cumulus_parameterization_input.grid_length,
            pbl_level=state.pbl_level,
            pbl_level_cu_param_input=locals.cumulus_parameterization_input.pbl_level,
        )

        self._set_local_state(
            geopotential_height_interface=state.geopotential_height_interface,
            t=state.t,
            t_local=locals.local_copy.t,
            p=locals.derived_state.p,
            p_local=locals.local_copy.p,
            vapor=state.vapor,
            vapor_local=locals.local_copy.vapor,
            vapor_current_local=locals.local_copy.vapor_current,
            u=state.u,
            u_local=locals.local_copy.u,
            v=state.v,
            v_local=locals.local_copy.v,
            vertical_velocity=locals.derived_state.vertical_velocity,
            vertical_velocity_local=locals.local_copy.vertical_velocity,
            layer_height_above_surface=locals.derived_state.layer_height_above_surface,
            layer_height_above_surface_local=locals.local_copy.layer_height_above_surface,
            edge_height_above_surface=locals.derived_state.edge_height_above_surface,
            edge_height_above_surface_local=locals.local_copy.edge_height_above_surface,
            mass=locals.derived_state.mass,
            mass_local=locals.local_copy.mass,
            scalar_diffusivity=state.scalar_diffusivity,
            scalar_diffusivity_local=locals.local_copy.scalar_diffusivity,
            lateral_entrainment_rate=state.lateral_entrainment_rate,
            lateral_entrainment_rate_local=locals.local_copy.lateral_entrainment_rate,
            p_surface=locals.cumulus_parameterization_input.p_surface,
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
            grid_scale_forcing_t=locals.derived_state.grid_scale_forcing_t,
            grid_scale_forcing_vapor=locals.derived_state.grid_scale_forcing_vapor,
            subgrid_scale_forcing_t=locals.derived_state.subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor=locals.derived_state.subgrid_scale_forcing_vapor,
            advective_forcing_t=locals.derived_state.advective_forcing_t,
            convective_liquid_local=locals.local_copy.convective_liquid,
            convective_ice_local=locals.local_copy.convective_ice,
            convective_cloud_fraction_local=locals.local_copy.convective_cloud_fraction,
            large_scale_liquid_local=locals.local_copy.large_scale_liquid,
            large_scale_ice_local=locals.local_copy.large_scale_ice,
            large_scale_cloud_fraction_local=locals.local_copy.large_scale_cloud_fraction,
            convection_tracer=state.convection_tracer,
            buoyancy_excess=locals.cumulus_parameterization_input.buoyancy_excess,
        )

        if self.GF_2020_config.ADV_TRIGGER == 2:
            raise NotImplementedError("ADV_TRIGGER == 2 not implemented yet")

        self._prepare_cumulus_parameterization(
            t_perturbation_horizontal=locals.derived_state.t_perturbation_horizontal,
            t_perturbation_vertical=locals.derived_state.t_perturbation_vertical,
            t_perturbation=locals.cumulus_parameterization_input.t_perturbation,
            ccn=locals.cumulus_parameterization_input.ccn,
            p_surface=locals.cumulus_parameterization_input.p_surface,
            topography_height=locals.derived_state.topography_height,
            topography_height_no_negative=locals.cumulus_parameterization_input.topography_height,
            latitude=state.latitude,
            longitude=state.longitude,
            latitude_degrees=locals.cumulus_parameterization_input.latitude_degrees,
            longitude_degrees=locals.cumulus_parameterization_input.longitude_degrees,
            layer_height_above_surface=locals.local_copy.layer_height_above_surface,
            geopotential_height_cu_param_input=locals.cumulus_parameterization_input.geopotential_height,
            p=locals.local_copy.p,
            p_mb=locals.cumulus_parameterization_input.p_mb,
            t_local=locals.local_copy.t,
            t_cu_param_input=locals.cumulus_parameterization_input.t,
            vapor_local=locals.local_copy.vapor,
            vapor_current_local=locals.local_copy.vapor_current,
            vapor_timestep_start_cu_param_input=locals.cumulus_parameterization_input.vapor_timestep_start,
            vapor_current_cu_param_input=locals.cumulus_parameterization_input.vapor_current,
            air_density_cu_param_input=locals.cumulus_parameterization_input.air_density,
            u_local=locals.local_copy.u,
            v_local=locals.local_copy.v,
            w_local=locals.local_copy.vertical_velocity,
            u_cu_param_input=locals.cumulus_parameterization_input.u,
            v_cu_param_input=locals.cumulus_parameterization_input.v,
            w_cu_param_input=locals.cumulus_parameterization_input.w,
            omega_cu_param_input=locals.cumulus_parameterization_input.omega,
            mass_local=locals.local_copy.mass,
            mass_cu_param_input=locals.cumulus_parameterization_input.mass,
            advective_forcing_t=locals.derived_state.advective_forcing_t,
            t_modified_by_advection=locals.cumulus_parameterization_input.t_modified_by_advection,
            grid_scale_forcing_vapor=locals.derived_state.grid_scale_forcing_vapor,
            vapor_modified_by_advection=locals.cumulus_parameterization_input.vapor_modified_by_advection,
            pbl_level_cu_param_input=locals.cumulus_parameterization_input.pbl_level,
            pbl_height_cu_param_input=locals.cumulus_parameterization_input.pbl_height,
            sensible_heat_flux_local=locals.local_copy.sensible_heat_flux,
            sensible_heat_flux_cu_param_input=locals.cumulus_parameterization_input.sensible_heat_flux,
            evaporation_local=locals.local_copy.evaporation,
            latent_heat_flux_cu_param_input=locals.cumulus_parameterization_input.latent_heat_flux,
            convective_scale_velosity_cu_param_input=locals.cumulus_parameterization_input.convective_scale_velosity,
            t_excess_cu_param_input=locals.cumulus_parameterization_input.t_excess,
            vapor_excess_cu_param_input=locals.cumulus_parameterization_input.vapor_excess,
        )

        if self.GF_2020_config.APPLY_SUB_MP == 1:
            self._prepare_cumulus_paramaterization_microphysics(
                convective_liquid_local=locals.local_copy.convective_liquid,
                convective_ice_local=locals.local_copy.convective_ice,
                convective_cloud_fraction_local=locals.local_copy.convective_cloud_fraction,
                large_scale_liquid_local=locals.local_copy.large_scale_liquid,
                large_scale_ice_local=locals.local_copy.large_scale_ice,
                large_scale_cloud_fraction_local=locals.local_copy.large_scale_cloud_fraction,
                convective_liquid_cu_param_input=locals.cumulus_parameterization_input.convective_liquid,
                convective_ice_cu_param_input=locals.cumulus_parameterization_input.convective_ice,
                convective_cloud_fraction_cu_param_input=locals.cumulus_parameterization_input.convective_cloud_fraction,
                large_scale_liquid_cu_param_input=locals.cumulus_parameterization_input.large_scale_liquid,
                large_scale_ice_cu_param_input=locals.cumulus_parameterization_input.large_scale_ice,
                large_scale_cloud_fraction_cu_param_input=locals.cumulus_parameterization_input.large_scale_cloud_fraction,
            )

        if self.GF_2020_config.USE_TRACER_TRANSP == 1:
            ndsl_log.warning("tracer stuff not yet implemented")
