from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval, log, log10, sqrt
from ndsl.dsl.typing import FloatField, FloatFieldIJ

from pyMoist.microphysics.GFDL_1M.driver.constants import constants


def fall_speed(
    mixing_ratio_liquid: FloatField,
    mixing_ratio_ice: FloatField,
    mixing_ratio_snow: FloatField,
    mixing_ratio_graupel: FloatField,
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
    """Calculate the vertical fall speed of precipitation

    Arguments:
        mixing_ratio_liquid (FloatField): (inout) in cloud liquid mixing radio (kg/kg)
        mixing_ratio_ice (FloatField): (inout) in cloud ice mixing radio (kg/kg)
        mixing_ratio_snow (FloatField): (inout) in cloud snow mixing radio (kg/kg)
        mixing_ratio_graupel (FloatField): (inout) in cloud graupel mixing radio (kg/kg)
        t_unmodified (FloatField): (in) atmospheric temperature, unmodified throughout the driver (K)
        t (FloatField): (in) atmospheric temperature, modified throughout the driver (K)
        dz_unmodified (FloatField): (in) layer thickness (m)
        dz (FloatField): (out) layer thickness
        density_unmodified (FloatField): (in) density
        density (FloatField): (out) density, potentially corrected for hydrostatic balance
        density_factor (FloatField): (out) details unknown
        ice_terminal_velocity (FloatField): (out) terminal fall speed for ice
        snow_terminal_velocity (FloatField): (out) terminal fall speed for snow
        graupel_terminal_velocity (FloatField): (out) terminal fall speed for graupel
        convection_fraction (FloatFieldIJ): (in) convection fraction
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
            if mixing_ratio_ice < constants.THI:
                ice_terminal_velocity = constants.VF_MIN
            else:
                # -----------------------------------------------------------------------
                # ice:
                # -----------------------------------------------------------------------

                vi1 = 0.01 * vi_fac
                tc = t - constants.TICE  # deg C
                IWC = mixing_ratio_ice * density * 1.0e3  # Units are g/m3
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
            if mixing_ratio_snow < constants.THS:
                snow_terminal_velocity = constants.VF_MIN
            else:
                snow_terminal_velocity = (
                    vs_fac
                    * constants.VCONS
                    * rhof
                    * exp(0.0625 * log(mixing_ratio_snow * density / constants.NORMS))
                )
                snow_terminal_velocity = min(vs_max, max(constants.VF_MIN, snow_terminal_velocity))

        # -----------------------------------------------------------------------
        # graupel:
        # -----------------------------------------------------------------------

        if const_vg == True:  # noqa
            graupel_terminal_velocity = vg_fac  # 2.
        else:
            if mixing_ratio_graupel < constants.THG:
                graupel_terminal_velocity = constants.VF_MIN
            else:
                graupel_terminal_velocity = (
                    vg_fac
                    * constants.VCONG
                    * rhof
                    * sqrt(sqrt(sqrt(mixing_ratio_graupel * density / constants.NORMG)))
                )
                graupel_terminal_velocity = min(vg_max, max(constants.VF_MIN, graupel_terminal_velocity))
