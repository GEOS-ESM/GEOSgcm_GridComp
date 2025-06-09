from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, BoolFieldIJ, Int, IntFieldIJ
from ndsl.dsl.gt4py import (
    computation,
    interval,
    PARALLEL,
    FORWARD,
    function,
    log,
    BACKWARD,
)
from gt4py.cartesian.gtscript import THIS_K
from pyMoist.saturation_tables.qsat_functions import saturation_specific_humidity
from pyMoist.field_types import GlobalTable_saturaion_tables
from pyMoist.constants import (
    MAPL_GRAV,
    MAPL_RGAS,
    MAPL_RVAP,
    MAPL_CPDRY,
    MAPL_CPVAP,
    MAPL_P00,
    MAPL_KAPPA,
    MAPL_ALHL,
    MAPL_CP,
)


def calculate_derived_states(
    p_interface: FloatField,
    p_interface_mb: FloatField,
    p_mb: FloatField,
    geopotential_height_interface: FloatField,
    edge_height_above_surface: FloatField,
    layer_height_above_surface: FloatField,
    layer_thickness: FloatField,
    layer_thinkness_negative: FloatField,
    dp: FloatField,
    mass: FloatField,
    t: FloatField,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    qsat: FloatField,
    dqsat: FloatField,
    u: FloatField,
    u_unmodified: FloatField,
    v: FloatField,
    v_unmodified: FloatField,
    temporary_3d: FloatField,
    th: FloatField,
):
    """
    Computes derived state fields required for the rest of the GFDL single moment
    microphysics module.

    This stencil MUST be built using Z_INTERFACE_DIM to funciton properly.
    """
    from __externals__ import k_end

    with computation(PARALLEL), interval(...):
        p_interface_mb = p_interface * 0.01
        edge_height_above_surface = geopotential_height_interface - geopotential_height_interface.at(K=k_end)
    with computation(FORWARD), interval(0, -1):
        p_mb = 0.5 * (p_interface_mb + p_interface_mb[0, 0, 1])
        layer_height_above_surface = 0.5 * (edge_height_above_surface + edge_height_above_surface[0, 0, 1])
        layer_thickness = edge_height_above_surface - edge_height_above_surface[0, 0, 1]
        layer_thinkness_negative = -1.0 * layer_thickness
        dp = p_interface[0, 0, 1] - p_interface
        mass = dp / MAPL_GRAV
        qsat, dqsat = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        u_unmodified = u
        v_unmodified = v
        temporary_3d = (100.0 * p_mb / MAPL_P00) ** (MAPL_KAPPA)
        th = t / temporary_3d


@function
def find_tlcl(
    t: Float,
    rh: Float,
):
    """
    Computes the LCL temperature

    Arguments:
        t: temperature at surface (K)
        rh: relative humidity at surface

    Returns:
        tlcl: LCL temperature
    """
    term1 = 1.0 / (t - 55.0)
    term2 = log(max(0.1, rh) / 100.0) / 2840.0
    denom = term1 - term2
    tlcl = (1.0 / denom) + 55.0
    return tlcl


def find_klcl(
    t: FloatField,
    p_mb: FloatField,
    vapor: FloatField,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    found_level: BoolFieldIJ,
    k_lcl: IntFieldIJ,
):
    """
    Find the level of the lifted condensation level (LCL).

    Arguments:
        t (in): Atmospheric temperature (K)
        p_mb (in): pressure (mb)
        vapor (in): water vapor mixing radio (kg/kg)
        ese (in): saturation vapor pressure table, details unknown
        esx (in): saturation vapor pressure table, details unknown
        found_level (in): boolean mask used to stop calculation after LCL is identified
        k_lcl (out): LCL level
    """
    from __externals__ import k_end

    # get LCL pressure
    with computation(FORWARD), interval(-1, None):
        qsat, _ = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        rhsfc = 100 * vapor / qsat
        tlcl = find_tlcl(t=t, rh=rhsfc)
        rm = (1 - vapor) * MAPL_RGAS + vapor * MAPL_RVAP
        cpm = (1.0 - vapor) * MAPL_CPDRY + vapor * MAPL_CPVAP
        plcl = p_mb * ((tlcl / t) ** (cpm / rm))

    # find nearest level <= LCL pressure
    with computation(BACKWARD), interval(...):
        if found_level == False:
            k_lcl = THIS_K
        if p_mb <= plcl.at(K=k_end):
            found_level = True

    # Reset mask for future use
    with computation(FORWARD), interval(0, 1):
        found_level = False


def vertical_interpolation(
    field: FloatField,
    interpolated_field: FloatFieldIJ,
    p_interface_mb: FloatField,
    target_pressure: Float,
    pb: FloatFieldIJ,
    pt: FloatFieldIJ,
    boolean_2d_mask: BoolFieldIJ,
    interface: bool = False,
):
    """
    Interpolate to a specific vertical level.

    Only works for non-interface fields. Must be constructed using Z_INTERFACE_DIM.

    Arguments:
        field (in): three dimensional field to be interpolated to a specific pressure
        interpolated_field (out): output two dimension field of interpolated values
        p_interface_mb (in): interface pressure in mb
        target_pressure (in): target pressure for interpolation in Pascals
        pb (in): placeholder 2d quantity, can be removed onces 2d temporaries are available
        pt (in): placeholder 2d quantity, can be removed onces 2d temporaries are available
        boolean_2d_mask (in): boolean mask to track when each cell is modified
        interface (in): specifies if input 'field' is an interface (True) or non-interface (False) field
    """
    # from __externals__ import k_end

    # mask tracks which points have been touched. check later on ensures that every point has been touched
    with computation(FORWARD), interval(0, 1):
        boolean_2d_mask = False

    with computation(PARALLEL), interval(...):
        p = log(p_interface_mb * 100)

    with computation(FORWARD), interval(-1, None):
        if interface == True:
            pb = p
    with computation(FORWARD), interval(-1, None):
        if interface == False:
            pb = 0.5 * (p[0, 0, -1] + p)

    with computation(BACKWARD), interval(0, -1):
        if interface == True:
            pt = p.at(K=THIS_K)
            if log(target_pressure) > pt and log(target_pressure) <= pb:
                al = (pb - log(target_pressure)) / (pb - pt)
                interpolated_field = field.at(K=THIS_K) * al + field.at(K=THIS_K + 1) * (1.0 - al)
            pb = pt

    with computation(BACKWARD), interval(1, -1):
        if interface == False:
            pt = 0.5 * (p.at(K=THIS_K - 1) + p.at(K=THIS_K))
            if log(target_pressure) > pt and log(target_pressure) <= pb:
                al = (pb - log(target_pressure)) / (pb - pt)
                interpolated_field = field.at(K=THIS_K - 1) * al + field.at(K=THIS_K) * (1.0 - al)
                boolean_2d_mask = True
            pb = pt

    with computation(FORWARD), interval(-2, -1):
        if interface == False:
            pt = 0.5 * (p + p[0, 0, -1])
            pb = 0.5 * (p + p[0, 0, 1])
            if log(target_pressure) > pb and log(target_pressure) <= p[0, 0, 1]:
                interpolated_field = field[0, 0, -1]

            # ensure every point was actually touched
            if boolean_2d_mask == False:
                interpolated_field = p

    # reset masks and temporaries for later use
    with computation(FORWARD), interval(0, 1):
        boolean_2d_mask = False
        pb = 0
        pt = 0


def find_eis(
    t: FloatField,
    th: FloatField,
    layer_height_above_surface: FloatField,
    t700: FloatFieldIJ,
    th700: FloatFieldIJ,
    z700: FloatFieldIJ,
    k_lcl: IntFieldIJ,
    ese: GlobalTable_saturaion_tables,
    esx: GlobalTable_saturaion_tables,
    lower_tropospheric_stability: FloatFieldIJ,
    estimated_inversion_strength: FloatFieldIJ,
):
    """
    Find estimated inversion strength. Returns Estimated Inversion
    Strength (K) according to Wood and Betherton, J.Climate, 2006.
    Based on Fortran code written by Donifan Barahona.
    """
    with computation(FORWARD), interval(-1, None):
        lower_tropospheric_stability = th700 - th
        zlcl = layer_height_above_surface.at(K=k_lcl - 1)

        # Simplified single adiabat eq4 of https://doi.org/10.1175/JCLI3988.1
        t850 = 0.5 * (t + t700)
        qs850, _ = saturation_specific_humidity(t=t850, p=100 * 850, ese=ese, esx=esx)
        gamma850 = (1.0 + (MAPL_ALHL * qs850 / (MAPL_RGAS * t850))) / (
            1.0 + (MAPL_ALHL * MAPL_ALHL * qs850 / (MAPL_CP * MAPL_RVAP * t850 * t850))
        )
        gamma850 = MAPL_GRAV / MAPL_CP * (1.0 - gamma850)
        estimated_inversion_strength = lower_tropospheric_stability - gamma850 * (z700 - zlcl)


def prepare_tendencies(
    u: FloatField,
    v: FloatField,
    t: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    du_dt: FloatField,
    dv_dt: FloatField,
    dt_dt: FloatField,
    dvapor_dt: FloatField,
    dliquid_dt: FloatField,
    dice_dt: FloatField,
    dcloud_fraction_dt: FloatField,
    drain_dt: FloatField,
    dsnow_dt: FloatField,
    dgraupel_dt: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        du_dt = u
        dv_dt = v
        dt_dt = t
        dvapor_dt = vapor
        dliquid_dt = convective_liquid + large_scale_liquid
        dice_dt = convective_ice + large_scale_ice
        dcloud_fraction_dt = convective_cloud_fraction + large_scale_cloud_fraction
        drain_dt = rain
        dsnow_dt = snow
        dgraupel_dt = graupel


def update_tendencies(
    u: FloatField,
    v: FloatField,
    t: FloatField,
    vapor: FloatField,
    rain: FloatField,
    snow: FloatField,
    graupel: FloatField,
    convective_liquid: FloatField,
    convective_ice: FloatField,
    large_scale_liquid: FloatField,
    large_scale_ice: FloatField,
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    du_dt: FloatField,
    dv_dt: FloatField,
    dt_dt: FloatField,
    dvapor_dt: FloatField,
    dliquid_dt: FloatField,
    dice_dt: FloatField,
    dcloud_fraction_dt: FloatField,
    drain_dt: FloatField,
    dsnow_dt: FloatField,
    dgraupel_dt: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        du_dt = (u - du_dt) / DT_MOIST
        dv_dt = (v - dv_dt) / DT_MOIST
        dt_dt = (t - dt_dt) / DT_MOIST
        dvapor_dt = (vapor - dvapor_dt) / DT_MOIST
        dliquid_dt = ((convective_liquid + large_scale_liquid) - dliquid_dt) / DT_MOIST
        dice_dt = ((convective_ice + large_scale_ice) - dice_dt) / DT_MOIST
        dcloud_fraction_dt = (
            (convective_cloud_fraction + large_scale_cloud_fraction) - dcloud_fraction_dt
        ) / DT_MOIST
        drain_dt = (rain - drain_dt) / DT_MOIST
        dsnow_dt = (snow - dsnow_dt) / DT_MOIST
        dgraupel_dt = (graupel - dgraupel_dt) / DT_MOIST


def prepare_radiation_quantities(
    convective_cloud_fraction: FloatField,
    large_scale_cloud_fraction: FloatField,
    radiation_cloud_fraction: FloatField,
    convective_liquid: FloatField,
    large_scale_liquid: FloatField,
    radiation_liquid: FloatField,
    convective_ice: FloatField,
    large_scale_ice: FloatField,
    radiation_ice: FloatField,
    vapor: FloatField,
    radiation_vapor: FloatField,
    rain: FloatField,
    radiation_rain: FloatField,
    snow: FloatField,
    radiation_snow: FloatField,
    graupel: FloatField,
    radiation_graupel: FloatField,
):
    # cloud fraction
    radiation_cloud_fraction = min(convective_cloud_fraction + large_scale_cloud_fraction, 1.0)
    # liquid
    radiation_liquid = convective_liquid + large_scale_liquid
    # ice
    radiation_ice = convective_ice + large_scale_ice
    # vapor
    radiation_vapor = vapor
    # RAIN
    radiation_rain = rain
    # snow
    radiation_snow = snow
    # graupel
    radiation_graupel = graupel


def apply_driver_tendencies(
    t: FloatField,
    u: FloatField,
    v: FloatField,
    radiation_cloud_fraction: FloatField,
    radiation_ice: FloatField,
    radiation_liquid: FloatField,
    radiation_vapor: FloatField,
    radiation_rain: FloatField,
    radiation_snow: FloatField,
    radiation_graupel: FloatField,
    dcloud_fraction_dt: FloatField,
    dt_dt: FloatField,
    du_dt: FloatField,
    dv_dt: FloatField,
    dice_dt: FloatField,
    dliquid_dt: FloatField,
    dvapor_dt: FloatField,
    drain_dt: FloatField,
    dsnow_dt: FloatField,
    dgraupel_dt: FloatField,
):
    from __externals__ import DT_MOIST

    with computation(PARALLEL), interval(...):
        t = t + dt_dt * DT_MOIST
        u = u + du_dt * DT_MOIST
        v = v + dv_dt * DT_MOIST
        radiation_cloud_fraction = min(
            1.0, max(0.0, radiation_cloud_fraction + dcloud_fraction_dt * DT_MOIST)
        )
        radiation_ice = radiation_ice + dice_dt * DT_MOIST
        radiation_liquid = radiation_liquid + dliquid_dt * DT_MOIST
        radiation_vapor = radiation_vapor + dvapor_dt * DT_MOIST
        radiation_rain = radiation_rain + drain_dt * DT_MOIST
        radiation_snow = radiation_snow + dsnow_dt * DT_MOIST
        radiation_graupel = radiation_graupel + dgraupel_dt * DT_MOIST
