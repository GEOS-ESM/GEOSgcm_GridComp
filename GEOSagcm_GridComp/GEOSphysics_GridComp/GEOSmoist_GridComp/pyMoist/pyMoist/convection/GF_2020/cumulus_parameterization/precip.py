import pyMoist.constants as constants
import pyMoist.convection.GF_2020.cumulus_parameterization.constants as cumulus_parameterization_constants
from ndsl.dsl.gt4py import BACKWARD, FORWARD, PARALLEL, K, computation, function, interval, sqrt
from ndsl.dsl.typing import FloatField, FloatFieldIJ, Int, IntFieldIJ
from pyMoist.convection.GF_2020.cumulus_parameterization.field_types import (
    FloatField_Plume,
    FloatFieldIJ_Plume,
    IntFieldIJ_Plume,
)
from pyMoist.shared.incloud_processes import ice_fraction


@function
def liquid_fraction(
    t,
    convection_fraction,
    surface_type,
    FRAC_MODIS,
):
    """
    Get the fraction of liquid condensates

    Args:
        t (in): temperature
        convection_fraction (in)
        surface_type (in)
        FRAC_MODIS (in): use fraction liq/ice content derived from MODIS/CALIPO sensors
    """
    if FRAC_MODIS == 1:
        liquid_fraction = 1.0 - ice_fraction(t, convection_fraction, surface_type)
    else:
        liquid_fraction = min(
            1.0,
            (
                max(0.0, (t - cumulus_parameterization_constants.T_ICE))
                / (cumulus_parameterization_constants.T_0 - cumulus_parameterization_constants.T_ICE)
            )
            ** 2,
        )

    return liquid_fraction


def partition_liquid_ice(
    t: FloatField,
    p: FloatField_Plume,
    geopotential_height: FloatField,
    topography_height_no_negative: FloatFieldIJ,
    surface_type: FloatFieldIJ,
    convection_fraction: FloatFieldIJ,
    error_code: IntFieldIJ_Plume,
    melting_layer: FloatField,
    part_liquid_ice: FloatField,
    plume: Int,
):
    """Partition total condensate into liquid and ice phases.

    Args:
        t (FloatField)
        p (FloatField_Plume)
        geopotential_height (FloatField)
        topography_height_no_negative (FloatFieldIJ)
        surface_type (FloatFieldIJ)
        convection_fraction (FloatFieldIJ)
        error_code (IntFieldIJ_Plume)
        melting_layer (FloatField)
        part_liquid_ice (FloatField)
        plume (Int)
    """
    from __externals__ import FRAC_MODIS, MELT_GLAC, k_end

    with computation(PARALLEL), interval(...):
        # constants, set internally because they may differ from global constants
        # and need to only exist inside this stencil
        t1 = 276.16
        z_meltlayer1 = 4000.0
        z_meltlayer2 = 6000.0
        del_t = 3.0

        # prefill some fields
        part_liquid_ice = 1.0
        melting_layer = 0.0

    with computation(PARALLEL), interval(0, -1):
        if MELT_GLAC == True and plume == cumulus_parameterization_constants.DEEP:
            if error_code[0, 0][plume] == 0:
                # get function of T for partition of total condensate into liq and ice phases
                part_liquid_ice = liquid_fraction(t, convection_fraction, surface_type, FRAC_MODIS)

    with computation(PARALLEL), interval(0, -1):
        if MELT_GLAC == True and plume == cumulus_parameterization_constants.DEEP:
            if error_code[0, 0][plume] == 0:
                # define the melting layer (the layer will be between T_0+1 < TEMP < T_1
                if t <= (cumulus_parameterization_constants.T_0 - del_t):
                    melting_layer = 0.0

                elif t < (cumulus_parameterization_constants.T_0 + del_t) and t > (
                    cumulus_parameterization_constants.T_0 - del_t
                ):
                    melting_layer = (
                        (t - (cumulus_parameterization_constants.T_0 - del_t)) / (2.0 * del_t)
                    ) ** 2

                else:
                    melting_layer = 1.0

                melting_layer = melting_layer * (1.0 - melting_layer)

    with computation(FORWARD), interval(0, 1):
        if MELT_GLAC == True and plume == cumulus_parameterization_constants.DEEP:
            # normalize vertical integral of melting_layer to 1
            norm: FloatFieldIJ = 0.0

    with computation(FORWARD), interval(0, -2):
        if MELT_GLAC == True and plume == cumulus_parameterization_constants.DEEP:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                dp = 100.0 * (p[0, 0, 0][plume] - p[0, 0, 1][plume])
                norm = norm + melting_layer * dp / constants.MAPL_GRAV

    with computation(PARALLEL), interval(...):
        if MELT_GLAC == True and plume == cumulus_parameterization_constants.DEEP:
            if error_code[0, 0][plume] == 0:
                # normalize vertical integral of melting_layer to 1
                melting_layer = (
                    melting_layer
                    / (norm + 1.0e-6)
                    * (
                        100
                        * (p.at(K=0, ddim=[plume]) - p.at(K=k_end - 1, ddim=[plume]))
                        / constants.MAPL_GRAV
                    )
                )


class PrecipFactor:
    """
    Get the pickup of ensemble average precipitation, following Neelin et al 2009.

    Fortran --> Python translation note: this runs in the fortran, but it does not modify inputs
    and the output is never used. For this reason, it is not implemented in the Python version.
    """

    def __init__(self):
        pass

    def __call__(self, *args, **kwds):
        pass


def get_precip_fluxes(
    error_code: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    epsilon_forced: FloatFieldIJ_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    precipitation_flux: FloatField,
    evaporation_flux: FloatField,
    plume: Int,
):
    """Compute the evaporation and precipitation flux throughout the entire column.

    Args:
        error_code (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        cloud_base_mass_flux_modified (FloatFieldIJ_Plume)
        epsilon_forced (FloatFieldIJ_Plume)
        condensate_to_fall_forced (FloatField_Plume)
        evaporate_in_downdraft_forced (FloatField_Plume)
        precipitation_flux (FloatField)
        evaporation_flux (FloatField)
        plume (Int)
    """
    with computation(PARALLEL), interval(...):
        precipitation_flux = 0.0
        evaporation_flux = 0.0

    with computation(BACKWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                precipitation_flux = precipitation_flux[0, 0, 1] + cloud_base_mass_flux_modified[0, 0][
                    plume
                ] * (
                    condensate_to_fall_forced[0, 0, 0][plume]
                    + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                precipitation_flux = max(0.0, precipitation_flux)

                evaporation_flux = (
                    evaporation_flux[0, 0, 1]
                    - cloud_base_mass_flux_modified[0, 0][plume]
                    * epsilon_forced[0, 0][plume]
                    * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                evaporation_flux = max(0.0, evaporation_flux)


def rain_evaporation_below_cloud_base(
    error_code: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    ocean_fraction: FloatFieldIJ,
    p_cloud_levels_forced: FloatField_Plume,
    p_surface: FloatFieldIJ,
    t_cloud_levels: FloatField,
    vapor_cloud_levels_forced: FloatField,
    environment_saturation_mixing_ratio_cloud_levels: FloatField,
    epsilon_forced: FloatFieldIJ_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    condensate_to_fall_forced: FloatField_Plume,
    evaporate_in_downdraft_forced: FloatField_Plume,
    precip: FloatFieldIJ_Plume,
    precipitation_flux: FloatField,
    evaporation_flux: FloatField,
    evaporation_below_cloud_base: FloatField,
    dtdt: FloatField_Plume,
    dvapordt: FloatField_Plume,
    dbuoyancydt: FloatField_Plume,
    plume: Int,
):
    """Allow rain the evaporate below the cloud base.

    Args:
        error_code (IntFieldIJ_Plume)
        updraft_lfc_level (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        ocean_fraction (FloatFieldIJ)
        p_cloud_levels_forced (FloatField_Plume)
        p_surface (FloatFieldIJ)
        t_cloud_levels (FloatField)
        vapor_cloud_levels_forced (FloatField)
        environment_saturation_mixing_ratio_cloud_levels (FloatField)
        epsilon_forced (FloatFieldIJ_Plume)
        cloud_base_mass_flux_modified (FloatFieldIJ_Plume)
        condensate_to_fall_forced (FloatField_Plume)
        evaporate_in_downdraft_forced (FloatField_Plume)
        precip (FloatFieldIJ_Plume)
        precipitation_flux (FloatField)
        evaporation_flux (FloatField)
        evaporation_below_cloud_base (FloatField)
        dtdt (FloatField_Plume)
        dvapordt (FloatField_Plume)
        dbuoyancydt (FloatField_Plume)
        plume (Int)
    """
    with computation(FORWARD), interval(0, 1):
        # setup internal constants
        alpha1: FloatFieldIJ = 5.44e-4
        alpha2: FloatFieldIJ = 5.09e-3
        alpha3: FloatFieldIJ = 0.5777
        c_conv: FloatFieldIJ = 0.05

    with computation(FORWARD), interval(0, 1):
        if plume == cumulus_parameterization_constants.SHALLOW:
            critical_rh_ocean: FloatFieldIJ = 1.0
            critical_rh_land: FloatFieldIJ = 1.0
            eff_c_conv: FloatFieldIJ = min(0.2, max(cloud_base_mass_flux_modified[0, 0][plume], c_conv))
        else:
            critical_rh_ocean: FloatFieldIJ = 0.95
            critical_rh_land: FloatFieldIJ = 0.85
            eff_c_conv: FloatFieldIJ = c_conv

        total_evaporation_below_cloud_base: FloatFieldIJ = 0.0

    with computation(PARALLEL), interval(...):
        precipitation_flux = 0.0
        evaporation_flux = 0.0

    with computation(FORWARD), interval(0, 1):
        if error_code[0, 0][plume] == 0:
            critical_rh: FloatFieldIJ = critical_rh_ocean * ocean_fraction + critical_rh_land * (
                1.0 - ocean_fraction
            )

    with computation(BACKWARD), interval(...):
        if error_code[0, 0][plume] == 0:
            if K <= cloud_top_level[0, 0][plume]:
                dp: FloatFieldIJ = 100.0 * (
                    p_cloud_levels_forced[0, 0, 0][plume] - p_cloud_levels_forced[0, 0, 1][plume]
                )

                if K <= updraft_lfc_level[0, 0][plume]:
                    vapor_deficit: FloatFieldIJ = max(
                        0.0,
                        (
                            critical_rh * environment_saturation_mixing_ratio_cloud_levels
                            - vapor_cloud_levels_forced
                        ),
                    )

                    evaporation_below_cloud_base = (
                        eff_c_conv
                        * alpha1
                        * vapor_deficit
                        * (
                            sqrt(p_cloud_levels_forced[0, 0, 0][plume] / p_surface)
                            / alpha2
                            * precipitation_flux[0, 0, 1]
                            / eff_c_conv
                        )
                        ** alpha3
                    )

                    evaporation_below_cloud_base = evaporation_below_cloud_base * dp / constants.MAPL_GRAV

                else:

                    evaporation_below_cloud_base = 0.0

                precipitation_flux = (
                    precipitation_flux[0, 0, 1]
                    - evaporation_below_cloud_base
                    + cloud_base_mass_flux_modified[0, 0][plume]
                    * (
                        condensate_to_fall_forced[0, 0, 0][plume]
                        + epsilon_forced[0, 0][plume] * evaporate_in_downdraft_forced[0, 0, 0][plume]
                    )
                )
                precipitation_flux = max(0.0, precipitation_flux)

                evaporation_flux = (
                    evaporation_flux[0, 0, 1]
                    + evaporation_below_cloud_base
                    - cloud_base_mass_flux_modified[0, 0][plume]
                    * epsilon_forced[0, 0][plume]
                    * evaporate_in_downdraft_forced[0, 0, 0][plume]
                )
                evaporation_flux = max(0.0, evaporation_flux)

                total_evaporation_below_cloud_base = (
                    total_evaporation_below_cloud_base + evaporation_below_cloud_base
                )

                del_vapor = evaporation_below_cloud_base * constants.MAPL_GRAV / dp
                del_t = (
                    -evaporation_below_cloud_base
                    * constants.MAPL_GRAV
                    / dp
                    * (cumulus_parameterization_constants.XLV / cumulus_parameterization_constants.CP)
                )

                dvapordt[0, 0, 0][plume] = dvapordt[0, 0, 0][plume] + del_vapor
                dtdt[0, 0, 0][plume] = dtdt[0, 0, 0][plume] + del_t
                dbuoyancydt[0, 0, 0][plume] = (
                    dbuoyancydt[0, 0, 0][plume]
                    + cumulus_parameterization_constants.CP * del_t
                    + cumulus_parameterization_constants.XLV * del_vapor
                )

                precip[0, 0][plume] = precip[0, 0][plume] - evaporation_below_cloud_base


def cloud_dissapation(
    error_code: IntFieldIJ_Plume,
    updraft_lfc_level: IntFieldIJ_Plume,
    cloud_top_level: IntFieldIJ_Plume,
    hydrostatic_air_density: FloatField,
    geopotential_height_forced: FloatField,
    t_cloud_levels_forced: FloatField,
    normalized_massflux_updraft_forced: FloatField_Plume,
    cloud_base_mass_flux_modified: FloatFieldIJ_Plume,
    vapor_cloud_levels_forced: FloatField,
    cloud_liquid_after_rain_forced: FloatField_Plume,
    environment_saturation_mixing_ratio_cloud_levels_forced: FloatField,
    environment_saturation_moist_static_energy_cloud_levels_forced: FloatField,
    vertical_velocity_3d: FloatField,
    scale_dependence_factor: FloatFieldIJ_Plume,
    dtdt: FloatField_Plume,
    dvapordt: FloatField_Plume,
    dcloudicedt: FloatField_Plume,
    plume: Int,
):
    """After excess water has precipitated, dissipate clouds which are no longer saturated.

    Args:
        error_code (IntFieldIJ_Plume)
        updraft_lfc_level (IntFieldIJ_Plume)
        cloud_top_level (IntFieldIJ_Plume)
        hydrostatic_air_density (FloatField)
        geopotential_height_forced (FloatField)
        t_cloud_levels_forced (FloatField)
        normalized_massflux_updraft_forced (FloatField_Plume)
        cloud_base_mass_flux_modified (FloatFieldIJ_Plume)
        vapor_cloud_levels_forced (FloatField)
        cloud_liquid_after_rain_forced (FloatField_Plume)
        environment_saturation_mixing_ratio_cloud_levels_forced (FloatField)
        environment_saturation_moist_static_energy_cloud_levels_forced (FloatField)
        vertical_velocity_3d (FloatField)
        scale_dependence_factor (FloatFieldIJ_Plume)
        dtdt (FloatField_Plume)
        dvapordt (FloatField_Plume)
        dcloudicedt (FloatField_Plume)
        plume (Int)
    """
    from __externals__ import COUPLE_MICROPHYSICS, DTIME, USE_CLOUD_DISSIPATION

    with computation(FORWARD), interval(0, 1):
        # setup internal constants
        cloud_lifetime: FloatFieldIJ = 1800.0
        version_x: IntFieldIJ = 2

    with computation(BACKWARD), interval(...):
        if (
            error_code[0, 0][plume] == 0
            and USE_CLOUD_DISSIPATION >= 0.0
            and K <= cloud_top_level[0, 0][plume]
            and K >= updraft_lfc_level[0, 0][plume]
        ):
            # cloud liq/ice remained in the convection plume
            precip_dissipation = max(
                0.0, cloud_liquid_after_rain_forced[0, 0, 0][plume] - dcloudicedt[0, 0, 0][plume] * DTIME
            )

            # get relative humidity
            f_rh = 0.0

            # estimation of the fractional area
            fractional_area = (
                (cloud_base_mass_flux_modified[0, 0][plume] / scale_dependence_factor[0, 0][plume])
                * normalized_massflux_updraft_forced[0, 0, 0][plume]
                / (hydrostatic_air_density * vertical_velocity_3d)
            )

            # source of enviroment moistening/cooling due to the 'remained' cloud dissipation into it.
            out_precip_dissipation = (precip_dissipation * (1.0 - f_rh)) / cloud_lifetime

            # NOTE other option (if this is true) is not implemented
            if not (version_x == 1 or COUPLE_MICROPHYSICS == False):
                dcloudicedt[0, 0, 0][plume] = (
                    dcloudicedt[0, 0, 0][plume]
                    + out_precip_dissipation * fractional_area * USE_CLOUD_DISSIPATION
                )

            cloud_liquid_after_rain_forced[0, 0, 0][plume] = max(
                0.0,
                cloud_liquid_after_rain_forced[0, 0, 0][plume]
                - out_precip_dissipation * USE_CLOUD_DISSIPATION * fractional_area * DTIME,
            )
