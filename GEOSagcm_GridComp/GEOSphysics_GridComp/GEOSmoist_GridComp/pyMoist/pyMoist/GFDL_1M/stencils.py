from ndsl.dsl.typing import Float, FloatFieldIJ, FloatField, BoolFieldIJ
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
)


def dumb_stencil(in_data: FloatField, out_data: FloatField):
    with computation(PARALLEL), interval(...):
        out_data = in_data


def calculate_derived_states(
    p_interface: FloatField,
    p_interface_mb: FloatField,
    p_mb: FloatField,
    geopotential_height_interface: FloatField,
    edge_height_above_surface: FloatField,
    layer_height_above_surface: FloatField,
    layer_thickness: FloatField,
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
        dp = p_interface[0, 0, 1] - p_interface
        mass = dp / MAPL_GRAV
        qsat, dqsat = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        u_unmodified = u
        v_unmodified = v


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
    k_lcl: FloatFieldIJ,
    test_field: FloatFieldIJ,
):
    """
    Find the level of the lifted condensation level (LCL).

    Arguments:
        t (in): Atmospheric temperature (K)
        p_mb (in): pressure (mb)
        vapor: water vapor mixing radio (kg/kg)
        ese (in): saturation vapor pressure table, details unknown
        esx (in): saturation vapor pressure table, details unknown
        found_level (in): boolean mask used to stop calculation after LCL is identified
        k_lcl (out): LCL level
    """
    from __externals__ import k_end

    # get LCL pressure
    with computation(FORWARD), interval(0, 1):
        qsat, _ = saturation_specific_humidity(t=t, p=p_mb * 100, ese=ese, esx=esx)
        rhsfc = 100 * vapor.at(K=k_end) / qsat
        test_field = qsat
        tlcl = find_tlcl(t=t, rh=rhsfc)
        rm = (1 - vapor.at(K=k_end)) * MAPL_RGAS + vapor.at(K=k_end) * MAPL_RVAP
        cpm = (1.0 - vapor.at(K=k_end)) * MAPL_CPDRY + vapor.at(K=k_end) * MAPL_CPVAP
        plcl = p_mb.at(K=k_end) * ((tlcl / t.at(K=k_end)) ** (cpm / rm))

    # find nearest level <= LCL pressure
    with computation(BACKWARD), interval(...):
        if found_level == False:
            k_lcl = THIS_K
        if p_mb <= plcl.at(K=0):
            found_level = True

    # Reset mask for future use
    with computation(FORWARD), interval(0, 1):
        found_level = False


# def find_eis(
#     t: FloatField,
#     p_mb: FloatField,
#     p_interface_mb: FloatField,
#     qsat: FloatField,
#     layer_height_above_surface: FloatField,
#     klcl: FloatFieldIJ,
#     lower_tropospheric_stability: FloatFieldIJ,
#     estimated_inversion_strength: FloatFieldIJ,
# ):
#     """
#     Find estimated inversion strength. Returns Estimated Inversion
#     Strength (K) according to Wood and Betherton, J.Climate, 2006.
#     Based on Fortran code written by Donifan Barahona.
#     """
#     with computation(PARALLEL):
#         with interval(...):
#             temporary = t/(100.0*p_mb/MAPL_P00)**(MAPL_KAPPA)
#         with interval(0,1):
#         lower_tropospheric_stability =
#           LTS(I,J) = TH700(I,J)-TH(I,J,LM)
#           ZLCL = ZL0(I,J,KLCL(I,J)-1)

#           ! Simplified single adiabat eq4 of https://doi.org/10.1175/JCLI3988.1
#            T850 = 0.5*(T(I,J,LM)+T700(I,J))
#           QS850 = GEOS_QSAT(T850, 850.0)

#           GAMMA850 =  (1.0+(          MAPL_ALHL*QS850/(        MAPL_RGAS*T850     )))/ &
#                       (1.0+(MAPL_ALHL*MAPL_ALHL*QS850/(MAPL_CP*MAPL_RVAP*T850*T850)))
#           GAMMA850 =  gravbcp*(1.0-GAMMA850)

#           EIS(I,J) =  LTS(I,J) - GAMMA850*(Z700(I,J) - ZLCL)
