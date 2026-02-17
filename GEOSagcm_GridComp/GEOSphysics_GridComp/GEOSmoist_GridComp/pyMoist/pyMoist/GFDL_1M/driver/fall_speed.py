from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval, log, log10, sqrt
from ndsl.dsl.typing import FloatField, FloatFieldIJ

from pyMoist.GFDL_1M.driver.constants import constants


def fall_speed(
    liquid: FloatField,
    ice: FloatField,
    snow: FloatField,
    graupel: FloatField,
    t_unmodified: FloatField,
    t: FloatField,
    dz_unmodified: FloatField,
    dz: FloatField,
    density_unmodified: FloatField,
    density: FloatField,
    density_factor: FloatField,
    ice_terminal_velocity: FloatField,
    snow_terminal_velocity: FloatField,
    graupel_terminal_velocity: FloatField,
    convection_fraction: FloatFieldIJ,
):
    """
    calculate the vertical fall speed of precipitation

    Arguments:
        liquid (inout): in cloud liquid mixing radio (kg/kg)
        ice (inout): in cloud ice mixing radio (kg/kg)
        snow (inout): in cloud snow mixing radio (kg/kg)
        graupel (inout): in cloud graupel mixing radio (kg/kg)
        t_unmodified (in): atmospheric temperature, unmodified throughout the driver (K)
        t (in): atmospheric temperature, modified throughout the driver (K)
        dz_unmodified (in): layer thickness (m)
        dz (out): layer thickness
        density_unmodified (in): density
        density (out): density, potentially corrected for hydrostatic balance
        density_factor (out): details unknown
        ice_terminal_velocity (out): terminal fall speed for ice
        snow_terminal_velocity (out): terminal fall speed for snow
        graupel_terminal_velocity (out): terminal fall speed for graupel
        convection_fraction (in): convection fraction

    Externals:
        const_vg (bool): controls constant vs computed terminal velocity for graupel
        const_vi (bool): controls constant vs computed terminal velocity for ice
        const_vs (bool): controls constant vs computed terminal velocity for snow
        p_nonhydro: performs hydrostatic adjustment on air density
        vg_fac (Float): details unknown
        vg_max (Float): details unknown
        vi_fac (Float): details unknown
        vi_max (Float): details unknown
        vs_fac (Float): details unknown
        vs_max (Float): details unknown
        anv_icefall (in): details unknown
        ls_icefall (in): details unknown

    """
    from __externals__ import (
        anv_icefall,
        const_vg,
        const_vi,
        const_vs,
        ls_icefall,
        p_nonhydro,
        vg_fac,
        vg_max,
        vi_fac,
        vi_max,
        vs_fac,
        vs_max,
    )

    with computation(PARALLEL), interval(...):
        if p_nonhydro:
            dz = dz_unmodified
            density = density_unmodified  # dry air density remains the same
            density_factor = sqrt(constants.SFCRHO / density)
        else:
            dz = dz_unmodified * t / t_unmodified  # hydrostatic balance
            density = density_unmodified * dz_unmodified / dz
            density_factor = sqrt(constants.SFCRHO / density)

        rhof = sqrt(min(10.0, constants.SFCRHO / density))
        if const_vi == True:  # noqa
            ice_terminal_velocity = vi_fac
        else:
            if ice < constants.THI:
                ice_terminal_velocity = constants.VF_MIN
            else:
                # -----------------------------------------------------------------------
                # ice:
                # -----------------------------------------------------------------------

                vi1 = 0.01 * vi_fac
                tc = t - constants.TICE  # deg C
                IWC = ice * density * 1.0e3  # Units are g/m3
                # -----------------------------------------------------------------------
                # use deng and mace (2008, grl)
                # https://doi.org/10.1029/2008GL035054
                # -----------------------------------------------------------------------
                viLSC = ls_icefall * 10.0 ** (
                    log10(IWC) * (tc * (constants.AAL * tc + constants.BBL) + constants.CCL)
                    + constants.DDL * tc
                    + constants.EEL
                )
                viCNV = anv_icefall * 10.0 ** (
                    log10(IWC) * (tc * (constants.AAC * tc + constants.BBC) + constants.CCC)
                    + constants.DDC * tc
                    + constants.EEC
                )
                # Combine
                ice_terminal_velocity = viLSC * (1.0 - convection_fraction) + viCNV * (convection_fraction)
                # Update units from cm/s to m/s
                ice_terminal_velocity = vi1 * ice_terminal_velocity
                # Limits
                ice_terminal_velocity = min(vi_max, max(constants.VF_MIN, ice_terminal_velocity))

        # -----------------------------------------------------------------------
        # snow:
        # -----------------------------------------------------------------------

        if const_vs == True:  # noqa
            snow_terminal_velocity = vs_fac  # 1. ifs_2016
        else:
            if snow < constants.THS:
                snow_terminal_velocity = constants.VF_MIN
            else:
                snow_terminal_velocity = (
                    vs_fac * constants.VCONS * rhof * exp(0.0625 * log(snow * density / constants.NORMS))
                )
                snow_terminal_velocity = min(vs_max, max(constants.VF_MIN, snow_terminal_velocity))

        # -----------------------------------------------------------------------
        # graupel:
        # -----------------------------------------------------------------------

        if const_vg == True:  # noqa
            graupel_terminal_velocity = vg_fac  # 2.
        else:
            if graupel < constants.THG:
                graupel_terminal_velocity = constants.VF_MIN
            else:
                graupel_terminal_velocity = (
                    vg_fac * constants.VCONG * rhof * sqrt(sqrt(sqrt(graupel * density / constants.NORMG)))
                )
                graupel_terminal_velocity = min(vg_max, max(constants.VF_MIN, graupel_terminal_velocity))
