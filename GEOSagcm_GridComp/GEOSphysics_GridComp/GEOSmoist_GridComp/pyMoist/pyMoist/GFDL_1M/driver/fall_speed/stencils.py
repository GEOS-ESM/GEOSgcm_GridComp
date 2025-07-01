from ndsl.dsl.gt4py import PARALLEL, computation, exp, interval, log, log10, sqrt

from ndsl.dsl.typing import FloatField, FloatFieldIJ
from pyMoist.GFDL_1M.driver.constants import constants


def fall_speed(
    ql: FloatField,
    qi: FloatField,
    qs: FloatField,
    qg: FloatField,
    t: FloatField,
    t1: FloatField,
    dz: FloatField,
    dz1: FloatField,
    den: FloatField,
    den1: FloatField,
    denfac: FloatField,
    vti: FloatField,
    vts: FloatField,
    vtg: FloatField,
    cnv_frc: FloatFieldIJ,
):
    """
    calculate the vertical fall speed of precipitation

    Arguments:
        ql (inout): in cloud liquid mixing radio (kg/kg)
        qi (inout): in cloud ice mixing radio (kg/kg)
        qs (inout): in cloud snow mixing radio (kg/kg)
        qg (inout): in cloud graupel mixing radio (kg/kg)
        t (in): atmopsheric temperature, unmodified throughout the driver (K)
        t1 (in): atmospheric temperature, modified thorughout the driver (K)
        dz (in): layer thickness (m)
        dz1 (out): layer thickness
        den (in): density
        den1 (out): density, potentially corrected for hydrostatic balance
        denfac (out): details unknown
        vti (out): terminal fall speed for ice
        vts (out): terminal fall speed for snow
        vtg (out): terminal fall speed for graupel
        cnv_frc (in): convection fraction

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

    reference Fortran: gfdl_cloud_microphys.F90: subroutines mpdrv and fall_speed
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
            dz1 = dz
            den1 = den  # dry air density remains the same
            denfac = sqrt(constants.SFCRHO / den1)
        else:
            dz1 = dz * t1 / t  # hydrostatic balance
            den1 = den * dz / dz1
            denfac = sqrt(constants.SFCRHO / den1)

        rhof = sqrt(min(10.0, constants.SFCRHO / den1))
        if const_vi == True:  # noqa
            vti = vi_fac
        else:
            if qi < constants.THI:
                vti = constants.VF_MIN
            else:
                # -----------------------------------------------------------------------
                # ice:
                # -----------------------------------------------------------------------

                vi1 = 0.01 * vi_fac
                tc = t1 - constants.TICE  # deg C
                IWC = qi * den1 * 1.0e3  # Units are g/m3
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
                vti = viLSC * (1.0 - cnv_frc) + viCNV * (cnv_frc)
                # Update units from cm/s to m/s
                vti = vi1 * vti
                # Limits
                vti = min(vi_max, max(constants.VF_MIN, vti))

        # -----------------------------------------------------------------------
        # snow:
        # -----------------------------------------------------------------------

        if const_vs == True:  # noqa
            vts = vs_fac  # 1. ifs_2016
        else:
            if qs < constants.THS:
                vts = constants.VF_MIN
            else:
                vts = vs_fac * constants.VCONS * rhof * exp(0.0625 * log(qs * den1 / constants.NORMS))
                vts = min(vs_max, max(constants.VF_MIN, vts))

        # -----------------------------------------------------------------------
        # graupel:
        # -----------------------------------------------------------------------

        if const_vg == True:  # noqa
            vtg = vg_fac  # 2.
        else:
            if qg < constants.THG:
                vtg = constants.VF_MIN
            else:
                vtg = vg_fac * constants.VCONG * rhof * sqrt(sqrt(sqrt(qg * den1 / constants.NORMG)))
                vtg = min(vg_max, max(constants.VF_MIN, vtg))
