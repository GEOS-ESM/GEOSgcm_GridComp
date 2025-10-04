import dataclasses
import copy

from ndsl import Quantity, QuantityFactory
from ndsl.constants import X_DIM, Y_DIM, Z_DIM, Z_INTERFACE_DIM
from ndsl.dsl.dace.orchestration import dace_inhibitor
from ndsl.dsl.typing import Int


@dataclasses.dataclass
class GF2020Temporaries:
    edge_height_above_surface: Quantity
    precip_flux: Quantity
    layer_height_above_surface: Quantity
    p: Quantity
    p_kappa: Quantity
    th: Quantity
    mass: Quantity
    evap_subl_tendency: Quantity
    vertical_velocity: Quantity
    ice_fraciton: Quantity
    saturation_specific_humidity: Quantity
    seed_convection: Quantity
    modified_area: Quantity
    t_2m_local: Quantity
    evaporation_local: Quantity
    sensible_heat_flux_local: Quantity
    pbl_level_local: Quantity
    t_local: Quantity
    p_local: Quantity
    p_mb_local: Quantity
    vapor_local: Quantity
    vapor_current_local: Quantity
    u_local: Quantity
    v_local: Quantity
    vertical_velocity_local: Quantity
    layer_height_above_surface_local: Quantity
    edge_height_above_surface_local: Quantity
    mass_local: Quantity
    scalar_diffusivity_local: Quantity
    lateral_entrainment_rate_local: Quantity
    convective_liquid_local: Quantity
    convective_ice_local: Quantity
    convective_cloud_fraction_local: Quantity
    large_scale_liquid_local: Quantity
    large_scale_ice_local: Quantity
    large_scale_cloud_fraction_local: Quantity
    topography_height: Quantity
    topography_height_no_negative: Quantity
    ocean_fraction: Quantity
    grid_length: Quantity
    p_surface: Quantity
    grid_scale_forcing_t: Quantity
    grid_scale_forcing_vapor: Quantity
    subgrid_scale_forcing_t: Quantity
    subgrid_scale_forcing_vapor: Quantity
    advective_forcing_t: Quantity
    buoyancy_excess: Quantity
    t_excess: Quantity
    vapor_excess: Quantity
    last_ierr: Quantity
    fix_out_vapor: Quantity
    conprr: Quantity
    evap_subl_tendency_cu_param: Quantity
    convective_precip_flux_cu_param: Quantity
    t_perturbation_cu_param_horizontal: Quantity
    t_perturbation_cu_param_vertical: Quantity
    t_perturbation_cu_param: Quantity
    omega_cu_param: Quantity
    ccn: Quantity
    dtdt_cu_param_shallow: Quantity
    dtdt_cu_param_mid: Quantity
    dtdt_cu_param_deep: Quantity
    dudt_cu_param_shallow: Quantity
    dudt_cu_param_mid: Quantity
    dudt_cu_param_deep: Quantity
    dvdt_cu_param_shallow: Quantity
    dvdt_cu_param_mid: Quantity
    dvdt_cu_param_deep: Quantity
    dvapordt_cu_param_shallow: Quantity
    dvapordt_cu_param_mid: Quantity
    dvapordt_cu_param_deep: Quantity
    dvapordt_cu_param_combined: Quantity
    dcloudicedt_cu_param_shallow: Quantity
    dcloudicedt_cu_param_mid: Quantity
    dcloudicedt_cu_param_deep: Quantity
    dnicedt_cu_param_shallow: Quantity
    dnicedt_cu_param_mid: Quantity
    dnicedt_cu_param_deep: Quantity
    dnliquiddt_cu_param_shallow: Quantity
    dnliquiddt_cu_param_mid: Quantity
    dnliquiddt_cu_param_deep: Quantity
    dbuoyancydt_cu_param_shallow: Quantity
    dbuoyancydt_cu_param_mid: Quantity
    dbuoyancydt_cu_param_deep: Quantity
    dconvectiveicedt_cu_param_shallow: Quantity
    dconvectiveicedt_cu_param_mid: Quantity
    dconvectiveicedt_cu_param_deep: Quantity
    dlargescaleicedt_cu_param_shallow: Quantity
    dlargescaleicedt_cu_param_mid: Quantity
    dlargescaleicedt_cu_param_deep: Quantity
    dconvectiveliquiddt_cu_param_shallow: Quantity
    dconvectiveliquiddt_cu_param_mid: Quantity
    dconvectiveliquiddt_cu_param_deep: Quantity
    dlargescaleliquiddt_cu_param_shallow: Quantity
    dlargescaleliquiddt_cu_param_mid: Quantity
    dlargescaleliquiddt_cu_param_deep: Quantity
    dconvectivecloudfractiondt_cu_param_shallow: Quantity
    dconvectivecloudfractiondt_cu_param_mid: Quantity
    dconvectivecloudfractiondt_cu_param_deep: Quantity
    dlargescalecloudfractiondt_cu_param_shallow: Quantity
    dlargescalecloudfractiondt_cu_param_mid: Quantity
    dlargescalecloudfractiondt_cu_param_deep: Quantity
    latitude_degrees: Quantity
    longitude_degrees: Quantity

    @classmethod
    def make(cls, quantity_factory: QuantityFactory):

        ## Used to compute/hold extra state-based inputs
        # 3D Z_INTERFACE_DIM Fields
        edge_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        precip_flux = quantity_factory.zeros([X_DIM, Y_DIM, Z_INTERFACE_DIM], "n/a")
        # 3D Z_DIM Fields
        layer_height_above_surface = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p_kappa = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        th = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mass = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        evap_subl_tendency = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vertical_velocity = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ice_fraciton = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        saturation_specific_humidity = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        # 2D Fields
        seed_convection = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        modified_area = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        # Used to hold local (modifiable) copy of the state, passed to cumulus parameterization
        t_2m_local = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        evaporation_local = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        sensible_heat_flux_local = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        pbl_level_local = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        t_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        p_mb_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vapor_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vapor_current_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        u_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        v_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        vertical_velocity_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        layer_height_above_surface_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        edge_height_above_surface_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        mass_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        scalar_diffusivity_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        lateral_entrainment_rate_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        convective_liquid_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        convective_ice_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        convective_cloud_fraction_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        large_scale_liquid_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        large_scale_ice_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        large_scale_cloud_fraction_local = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        # Extra inputs for cumulus parameterization
        topography_height = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        topography_height_no_negative = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        ocean_fraction = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        grid_length = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        p_surface = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        grid_scale_forcing_t = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        grid_scale_forcing_vapor = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        subgrid_scale_forcing_t = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        subgrid_scale_forcing_vapor = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        advective_forcing_t = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        buoyancy_excess = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t_excess = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        vapor_excess = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        last_ierr = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        fix_out_vapor = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        conprr = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        evap_subl_tendency_cu_param = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        convective_precip_flux_cu_param = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t_perturbation_cu_param_horizontal = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t_perturbation_cu_param_vertical = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        t_perturbation_cu_param = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        omega_cu_param = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        ccn = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        latitude_degrees = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")
        longitude_degrees = quantity_factory.zeros([X_DIM, Y_DIM], "n/a")

        # Outputs from cumulus parameterization
        dtdt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dtdt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dtdt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dudt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dudt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dudt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvdt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvdt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvdt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapordt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapordt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapordt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dvapordt_cu_param_combined = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcloudicedt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcloudicedt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dcloudicedt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dnicedt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dnicedt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dnicedt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dnliquiddt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dnliquiddt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dnliquiddt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dbuoyancydt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dbuoyancydt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dbuoyancydt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectiveicedt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectiveicedt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectiveicedt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescaleicedt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescaleicedt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescaleicedt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectiveliquiddt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectiveliquiddt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectiveliquiddt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescaleliquiddt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescaleliquiddt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescaleliquiddt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectivecloudfractiondt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectivecloudfractiondt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dconvectivecloudfractiondt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescalecloudfractiondt_cu_param_shallow = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescalecloudfractiondt_cu_param_mid = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        dlargescalecloudfractiondt_cu_param_deep = quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")

        return cls(
            edge_height_above_surface,
            precip_flux,
            layer_height_above_surface,
            p,
            p_kappa,
            th,
            mass,
            evap_subl_tendency,
            vertical_velocity,
            ice_fraciton,
            saturation_specific_humidity,
            seed_convection,
            modified_area,
            t_2m_local,
            evaporation_local,
            sensible_heat_flux_local,
            pbl_level_local,
            t_local,
            p_local,
            p_mb_local,
            vapor_local,
            vapor_current_local,
            u_local,
            v_local,
            vertical_velocity_local,
            layer_height_above_surface_local,
            edge_height_above_surface_local,
            mass_local,
            scalar_diffusivity_local,
            lateral_entrainment_rate_local,
            convective_liquid_local,
            convective_ice_local,
            convective_cloud_fraction_local,
            large_scale_liquid_local,
            large_scale_ice_local,
            large_scale_cloud_fraction_local,
            topography_height,
            topography_height_no_negative,
            ocean_fraction,
            grid_length,
            p_surface,
            grid_scale_forcing_t,
            grid_scale_forcing_vapor,
            subgrid_scale_forcing_t,
            subgrid_scale_forcing_vapor,
            advective_forcing_t,
            buoyancy_excess,
            t_excess,
            vapor_excess,
            last_ierr,
            fix_out_vapor,
            conprr,
            evap_subl_tendency_cu_param,
            convective_precip_flux_cu_param,
            t_perturbation_cu_param_horizontal,
            t_perturbation_cu_param_vertical,
            t_perturbation_cu_param,
            omega_cu_param,
            ccn,
            dtdt_cu_param_shallow,
            dtdt_cu_param_mid,
            dtdt_cu_param_deep,
            dudt_cu_param_shallow,
            dudt_cu_param_mid,
            dudt_cu_param_deep,
            dvdt_cu_param_shallow,
            dvdt_cu_param_mid,
            dvdt_cu_param_deep,
            dvapordt_cu_param_shallow,
            dvapordt_cu_param_mid,
            dvapordt_cu_param_deep,
            dvapordt_cu_param_combined,
            dcloudicedt_cu_param_shallow,
            dcloudicedt_cu_param_mid,
            dcloudicedt_cu_param_deep,
            dnicedt_cu_param_shallow,
            dnicedt_cu_param_mid,
            dnicedt_cu_param_deep,
            dnliquiddt_cu_param_shallow,
            dnliquiddt_cu_param_mid,
            dnliquiddt_cu_param_deep,
            dbuoyancydt_cu_param_shallow,
            dbuoyancydt_cu_param_mid,
            dbuoyancydt_cu_param_deep,
            dconvectiveicedt_cu_param_shallow,
            dconvectiveicedt_cu_param_mid,
            dconvectiveicedt_cu_param_deep,
            dlargescaleicedt_cu_param_shallow,
            dlargescaleicedt_cu_param_mid,
            dlargescaleicedt_cu_param_deep,
            dconvectiveliquiddt_cu_param_shallow,
            dconvectiveliquiddt_cu_param_mid,
            dconvectiveliquiddt_cu_param_deep,
            dlargescaleliquiddt_cu_param_shallow,
            dlargescaleliquiddt_cu_param_mid,
            dlargescaleliquiddt_cu_param_deep,
            dconvectivecloudfractiondt_cu_param_shallow,
            dconvectivecloudfractiondt_cu_param_mid,
            dconvectivecloudfractiondt_cu_param_deep,
            dlargescalecloudfractiondt_cu_param_shallow,
            dlargescalecloudfractiondt_cu_param_mid,
            dlargescalecloudfractiondt_cu_param_deep,
            latitude_degrees,
            longitude_degrees,
        )
