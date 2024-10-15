import copy

import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, FORWARD, sin
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatField,
    Float,
    Int,
    IntFieldIJ,
    IntField,
    FloatFieldIJ,
    BoolField,
)
from ndsl import StencilFactory, QuantityFactory
from pyMoist.field_types import FloatField_NTracers
from pyMoist.saturation.qsat import QSat, QSat_Float, FloatField_Extra_Dim
from pyMoist.saturation.formulation import SaturationFormulation
import pyMoist.constants as const


P00 = Float(1e5)  # Reference pressure
ZVIR = Float(0.609)  # r_H2O/r_air-1
ROVCP = const.MAPL_RGAS / const.MAPL_CP  # Gas constant over specific heat


@gtscript.function
def slope(
    kmask: Float,
    field: Float,
    field_above: Float,
    field_below: Float,
    p0: Float,
    p0_above: Float,
    p0_below: Float,
):
    """
    Function that calculates slope of a given field at each grid cell.

    Inputs:
    kmask (Float): K-level mask (e.g., 1 for k=0, 2 for k=1,71, 3 for k=72)
    field (Float): Field of interest [N/A]
    field_above (Float): The k-level above field (e.g., field[0,0,1]) [N/A]
    field_below (Float): The k-level below field (e.g., field[0,0,-1]) [N/A]
    p0 (Float): Pressure [Pa]
    p0_above (Float): The k-level above p0 (e.g., p0[0,0,1]) [Pa]
    p0_below (Float): The k-level below p0 (e.g., p0[0,0,-1]) [Pa]

    Returns:
    slope (Float): Slope of the field of interest [N/A]
    """
    if kmask == 1:
        value = (field_above - field) / (p0_above - p0)
        if value > 0.0:
            slope = max(0.0, value)
        else:
            slope = min(0.0, value)
    elif kmask == 2:
        above_value = (field_above - field) / (p0_above - p0)
        below_value = (field - field_below) / (p0 - p0_below)
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))
    else:
        above_value = (field_above - field) / (p0_above - p0)
        below_value = (field - field_below) / (p0 - p0_below)
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))

    return slope


@gtscript.function
def exnerfn(
    p: Float,
) -> Float:
    """
    Function that calculates the Exner function for a given pressure.

    Inputs:
    p (Float): Atmospheric pressure [Pa]

    Returns:
    (p / 100000.0) ** (const.MAPL_RGAS / const.MAPL_CP) (Float): Exner function
    """

    return (p / 100000.0) ** (const.MAPL_RGAS / const.MAPL_CP)


@gtscript.function
def ice_fraction(
    temp: Float,
    cnv_frc: Float,
    srf_type: Float,
):
    # Anvil clouds
    # Anvil-Convective sigmoidal function like figure 6(right)
    # Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
    if temp <= const.JaT_ICE_ALL:
        icefrct_c = 1.000
    elif temp > const.JaT_ICE_ALL and temp <= const.JaT_ICE_MAX:
        icefrct_c = sin(
            0.5
            * const.MAPL_PI
            * (
                1.00
                - (temp - const.JaT_ICE_ALL) / (const.JaT_ICE_MAX - const.JaT_ICE_ALL)
            )
        )
    else:
        icefrct_c = 0.00
    icefrct_c = max(min(icefrct_c, 1.00), 0.00) ** const.aICEFRPWR
    # Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if srf_type == 2.0:
        if temp <= const.JiT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > const.JiT_ICE_ALL and temp <= const.JiT_ICE_MAX:
            icefrct_m = 1.00 - (temp - const.JiT_ICE_ALL) / (
                const.JiT_ICE_MAX - const.JiT_ICE_ALL
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** const.iICEFRPWR
    elif srf_type > 1.0:
        if temp <= const.lT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > const.lT_ICE_ALL and temp <= const.lT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * const.MAPL_PI
                * (
                    1.00
                    - (temp - const.lT_ICE_ALL) / (const.lT_ICE_MAX - const.lT_ICE_ALL)
                )
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** const.lICEFRPWR
    else:
        if temp <= const.oT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > const.oT_ICE_ALL and temp <= const.oT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * const.MAPL_PI
                * (
                    1.00
                    - (temp - const.oT_ICE_ALL) / (const.oT_ICE_MAX - const.oT_ICE_ALL)
                )
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** const.oICEFRPWR
    ice_frac = icefrct_m * (1.0 - cnv_frc) + icefrct_c * cnv_frc
    return ice_frac


@gtscript.function
def conden(
    p: Float,
    thl: Float,
    qt: Float,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
):
    """
    Function that determines if condensation process has occurred.

    Inputs:
    p (Float): Pressure [Pa]
    thl (Float): Temperature [K]
    qt (Float): Mixing ratio [kg/kg]
    ese (FloatField_Extra_Dim): Used in QSat_Float [N/A]
    esx (FloatField_Extra_Dim): Used in for QSat_Float [N/A]

    Returns:
    th (Float): Temperature [K]
    qv (Float): Water vapor mixing ratio [kg/kg]
    ql (Float): Liquid water mixing ratio [kg/kg]
    qi (Float): Ice water mixing ratio [kg/kg]
    rvls (Float): Saturation specific humidity [kg/kg]
    id_check (Int): Flag that indicates if condensation occurs (0 for no condensation, 1 for condensation).
    """

    tc = thl * exnerfn(p)  # tc is a real*8 (64-bit) in the Fortran

    nu = ice_fraction(tc, 0.0, 0.0)  # nu is a real*8 (64-bit) in the Fortran
    leff = (1.0 - nu) * const.MAPL_ALHL + (
        nu * const.MAPL_ALHS
    )  # Effective latent heat, nu is a real*8 (64-bit) in the Fortran
    temps = tc
    ps = p
    qs, _ = QSat_Float(
        ese, esx, temps, ps / 100.0
    )  # qs is a real*8 (64-bit) in the Fortran
    rvls = qs

    if qs >= qt:  # no condensation
        id_check = 0
        qv = qt
        qc = 0.0
        ql = 0.0
        qi = 0.0
        th = thl
    else:  # condensation
        iteration = 0
        while iteration < 10:
            temps = temps + ((tc - temps) * const.MAPL_CP / leff + qt - rvls) / (
                const.MAPL_CP / leff
                + const.EPSILON * leff * rvls / (const.MAPL_RGAS * temps * temps)
            )
            qs, _ = QSat_Float(
                ese, esx, temps, ps / 100.0
            )  # qs is a real*8 (64-bit) in the Fortran
            rvls = qs
            iteration += 1
        qc = max(qt - qs, 0.0)  # qc is a real*8 (64-bit) in the Fortran
        qv = qt - qc
        ql = qc * (1.0 - nu)
        qi = nu * qc
        th = temps / exnerfn(p)
        if abs((temps - (leff / const.MAPL_CP) * qc) - tc) >= 1.0:
            id_check = 1
        else:
            id_check = 0

    return th, qv, ql, qi, rvls, id_check


def compute_uwshcu(
    dotransport: Float,
    ncnst: Int,
    k0: Int,
    kpbl_in: IntFieldIJ,
    pifc0_in: FloatField,
    pmid0_in: FloatField,
    exnmid0_in: FloatField,
    u0_in: FloatField,
    v0_in: FloatField,
    qv0_in: FloatField,
    ql0_in: FloatField,
    qi0_in: FloatField,
    th0_in: FloatField,
    tr0_inout: FloatField_NTracers,
    umf_out: FloatField,
    dcm_out: FloatField,
    qvten_out: FloatField,
    qlten_out: FloatField,
    qiten_out: FloatField,
    sten_out: FloatField,
    uten_out: FloatField,
    vten_out: FloatField,
    qrten_out: FloatField,
    qsten_out: FloatField,
    cufrc_out: FloatField,
    fer_out: FloatField,
    fdr_out: FloatField,
    qldet_out: FloatField,
    qidet_out: FloatField,
    qlsub_out: FloatField,
    qisub_out: FloatField,
    ndrop_out: FloatField,
    nice_out: FloatField,
    tpert_out: FloatFieldIJ,
    qpert_out: FloatFieldIJ,
    qtflx_out: FloatField,
    slflx_out: FloatField,
    uflx_out: FloatField,
    vflx_out: FloatField,
    tr0: FloatField_NTracers,
    ssthl0: FloatField,
    ssqt0: FloatField,
    ssu0: FloatField,
    ssv0: FloatField,
    sstr0: FloatField_NTracers,
    thj: FloatField,
    qvj: FloatField,
    qlj: FloatField,
    qij: FloatField,
    qse: FloatField,
    id_check: IntField,
    kmask: FloatField,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
    id_exit: BoolField,
    id_check_flag: BoolField,
):
    """
    University of Washington Shallow Convection Scheme

    Described in Park and Bretherton. 2008. J. Climate :

    'The University of Washington shallow convection and
    moist turbulent schemes and their impact on climate
    simulations with the Community Atmosphere Model'

    Coded in CESM by Sungsu Park. Oct.2005. May.2008.

    Coded in GEOS by Nathan Arnold. July 2016.

    For general questions, email sungsup@ucar.edu or
    sungsu@atmos.washington.edu

    For GEOS-specific questions, email nathan.arnold@nasa.gov

    In the original Fortran some constants where renamed.
    We use the original names, follows is a correspondence,
    between the names we use and the original fortran:
    MAPL_ALHL             = xlv
    MAPL_ALHF             = xlf
    MAPL_ALHS             = xls
    MAPL_CP               = cp
    MAPL_RGAS             = r
    MAPL_GRAV             = g
    EPSILON               = ep2


    Parameters:
    dotransport (Float in): Transport tracers [1 true]
    ncnst (Int in): Number of tracers
    k0 (Int in): Number of vertical levels
    kpbl_in (2D in): Boundary layer top layer index
    pifc0_in (3D in): Environmental pressure at interfaces [Pa]
    pmid0_in (3D in): Environmental pressure at midpoints [Pa]
    exnmid0_in (3D in): Exner function at midpoints
    u0_in (3D in): Environmental zonal wind [m/s]
    v0_in (3D in): Environmental meridional wind [m/s]
    qv0_in (3D in): Environmental specific humidity
    ql0_in (3D in): Environmental liquid water specific humidity
    qi0_in (3D in): Environmental ice specific humidity
    th0_in (3D in): Environmental potential temperature [K]
    tr0_inout (4D in): Environmental tracers [ #, kg/kg ]
    umf_out (3D out): Updraft mass flux at the interfaces [ kg/m2/s ]
    dcm_out (3D out): Detrained cloudy air mass
    qvten_out (3D out): Tendency of water vapor specific humidity [ kg/kg/s ]
    qlten_out (3D out): Tendency of liquid water specific humidity [ kg/kg/s ]
    qiten_out (3D out): Tendency of ice specific humidity [ kg/kg/s ]
    sten_out (3D out): Tendency of dry static energy [ J/kg/s ]
    uten_out (3D out): Tendency of zonal wind [ m/s2 ]
    vten_out (3D out): Tendency of meridional wind [ m/s2 ]
    qrten_out (3D out): Tendency of rain water specific humidity [ kg/kg/s ]
    qsten_out (3D out): Tendency of snow specific humidity [ kg/kg/s ]
    cufrc_out (3D out): Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    fer_out (3D out): Fractional lateral entrainment rate [ 1/Pa ]
    fdr_out (3D out): Fractional lateral detrainment rate [ 1/Pa ]
    qldet_out (3D out): [N/A]
    qidet_out (3D out): [N/A]
    qlsub_out (3D out): [N/A]
    qisub_out (3D out): [N/A]
    ndrop_out (3D out): [N/A]
    nice_out (3D out): [N/A]
    tpert_out (2D out): Temperature perturbation
    qpert_out (2D out): Humidity perturbation
    qtflx_out (3D out): [N/A]
    slflx_out (3D out): [N/A]
    uflx_out (3D out): [N/A]
    vflx_out (3D out): [N/A]
    tr0 (4D in): Environmental tracers [ #, kg/kg ]
    ssthl0 (3D inout): [N/A]
    ssqt0 (3D inout): [N/A]
    ssu0 (3D inout): [N/A]
    ssv0 (3D inout): [N/A]
    sstr0 (4D inout): [N/A]
    thj (3D inout): [N/A]
    qvj (3D inout): [N/A]
    qlj (3D inout): [N/A]
    qij (3D inout): [N/A]
    qse (3D inout): [N/A]
    id_check (3D in): 1 for condensation, 0 for no condensation
    kmask (3D in): K-level mask (e.g., 1 for k=0, 2 for k=1,71, 3 for k=72)
    ese: FloatField_Extra_Dim,
    esx:: FloatField_Extra_Dim,
    id_exit (3D inout): Flag that will raise an error if id_exit == False
    id_check_flag (3D inout): Flag that will raise an error if id_check == 1

    """

    with computation(FORWARD), interval(...):
        # Initialize output variables
        umf_out = 0.0
        dcm_out = 0.0
        cufrc_out = 0.0
        fer_out = const.MAPL_UNDEF
        fdr_out = const.MAPL_UNDEF
        qldet_out = 0.0
        qidet_out = 0.0
        qlsub_out = 0.0
        qisub_out = 0.0
        ndrop_out = 0.0
        nice_out = 0.0
        qtflx_out[0, 0, 1] = 0.0
        slflx_out[0, 0, 1] = 0.0
        uflx_out[0, 0, 1] = 0.0
        vflx_out[0, 0, 1] = 0.0
        tpert_out = 0.0
        qpert_out = 0.0

    # Start Main Calculation
    with computation(PARALLEL), interval(0, 1):
        pmid0 = pmid0_in
        pmid0_above = pmid0_in[0, 0, 1]
        u0 = u0_in
        u0_above = u0_in[0, 0, 1]
        v0 = v0_in
        v0_above = v0_in[0, 0, 1]
        qv0 = qv0_in
        qv0_above = qv0_in[0, 0, 1]
        ql0 = ql0_in
        ql0_above = ql0_in[0, 0, 1]
        qi0 = qi0_in
        qi0_above = qi0_in[0, 0, 1]
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        exnmid0 = exnmid0_in
        exnmid0_above = exnmid0_in[0, 0, 1]
        t0 = th0_in * exnmid0
        t0_above = th0_in[0, 0, 1] * exnmid0_above
        thl0 = (
            t0
            - ((const.MAPL_ALHL * ql0) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0) / const.MAPL_CP)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((const.MAPL_ALHL * ql0_above) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0_above) / const.MAPL_CP)
        ) / exnmid0_above

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        ssthl0 = slope(
            kmask, thl0, thl0_above, thl0_above, pmid0, pmid0_above, pmid0_above
        )
        ssqt0 = slope(kmask, qt0, qt0_above, qt0_above, pmid0, pmid0_above, pmid0_above)
        ssu0 = slope(kmask, u0, u0_above, u0_above, pmid0, pmid0_above, pmid0_above)
        ssv0 = slope(kmask, v0, v0_above, v0_above, pmid0, pmid0_above, pmid0_above)

        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                sstr0[0, 0, 0][n] = slope(
                    kmask,
                    tr0[0, 0, 0][n],
                    tr0[0, 0, 1][n],
                    tr0[0, 0, 1][n],
                    pmid0,
                    pmid0_above,
                    pmid0_above,
                )
                n += 1

    with computation(PARALLEL), interval(1, -1):
        pmid0 = pmid0_in
        pmid0_above = pmid0_in[0, 0, 1]
        pmid0_below = pmid0_in[0, 0, -1]
        u0 = u0_in
        u0_above = u0_in[0, 0, 1]
        u0_below = u0_in[0, 0, -1]
        v0 = v0_in
        v0_above = v0_in[0, 0, 1]
        v0_below = v0_in[0, 0, -1]
        qv0 = qv0_in
        qv0_above = qv0_in[0, 0, 1]
        qv0_below = qv0_in[0, 0, -1]
        ql0 = ql0_in
        ql0_above = ql0_in[0, 0, 1]
        ql0_below = ql0_in[0, 0, -1]
        qi0 = qi0_in
        qi0_above = qi0_in[0, 0, 1]
        qi0_below = qi0_in[0, 0, -1]

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        exnmid0 = exnmid0_in
        exnmid0_above = exnmid0_in[0, 0, 1]
        exnmid0_below = exnmid0_in[0, 0, -1]
        t0 = th0_in * exnmid0
        t0_above = th0_in[0, 0, 1] * exnmid0_above
        t0_below = th0_in[0, 0, -1] * exnmid0_below
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        qt0_below = qv0_below + ql0_below + qi0_below
        thl0 = (
            t0
            - ((const.MAPL_ALHL * ql0) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0) / const.MAPL_CP)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((const.MAPL_ALHL * ql0_above) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0_above) / const.MAPL_CP)
        ) / exnmid0_above
        thl0_below = (
            t0_below
            - ((const.MAPL_ALHL * ql0_below) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0_below) / const.MAPL_CP)
        ) / exnmid0_below

        ssthl0 = slope(
            kmask, thl0, thl0_above, thl0_below, pmid0, pmid0_above, pmid0_below
        )
        ssqt0 = slope(kmask, qt0, qt0_above, qt0_below, pmid0, pmid0_above, pmid0_below)
        ssu0 = slope(kmask, u0, u0_above, u0_below, pmid0, pmid0_above, pmid0_below)
        ssv0 = slope(kmask, v0, v0_above, v0_below, pmid0, pmid0_above, pmid0_below)

        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                sstr0[0, 0, 0][n] = slope(
                    kmask,
                    tr0[0, 0, 0][n],
                    tr0[0, 0, 1][n],
                    tr0[0, 0, -1][n],
                    pmid0,
                    pmid0_above,
                    pmid0_below,
                )
                n += 1

    with computation(PARALLEL), interval(-1, None):
        pmid0 = pmid0_in[0, 0, -1]
        pmid0_above = pmid0_in
        pmid0_below = pmid0_in[0, 0, -2]
        u0 = u0_in[0, 0, -1]
        u0_above = u0_in
        u0_below = u0_in[0, 0, -2]
        v0 = v0_in[0, 0, -1]
        v0_above = v0_in
        v0_below = v0_in[0, 0, -2]
        qv0 = qv0_in[0, 0, -1]
        qv0_above = qv0_in
        qv0_below = qv0_in[0, 0, -2]
        ql0 = ql0_in[0, 0, -1]
        ql0_above = ql0_in
        ql0_below = ql0_in[0, 0, -2]
        qi0 = qi0_in[0, 0, -1]
        qi0_above = qi0_in
        qi0_below = qi0_in[0, 0, -2]

        exnmid0 = exnmid0_in[0, 0, -1]
        exnmid0_above = exnmid0_in
        exnmid0_below = exnmid0_in[0, 0, -2]
        t0 = th0_in[0, 0, -1] * exnmid0
        t0_above = th0_in * exnmid0_above
        t0_below = th0_in[0, 0, -2] * exnmid0_below
        qt0 = qv0 + ql0 + qi0
        qt0_above = qv0_above + ql0_above + qi0_above
        qt0_below = qv0_below + ql0_below + qi0_below
        thl0 = (
            t0
            - ((const.MAPL_ALHL * ql0) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0) / const.MAPL_CP)
        ) / exnmid0
        thl0_above = (
            t0_above
            - ((const.MAPL_ALHL * ql0_above) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0_above) / const.MAPL_CP)
        ) / exnmid0_above
        thl0_below = (
            t0_below
            - ((const.MAPL_ALHL * ql0_below) / const.MAPL_CP)
            - ((const.MAPL_ALHS * qi0_below) / const.MAPL_CP)
        ) / exnmid0_below

        if dotransport == 1.0:
            n = 0
            # Loop over tracers
            while n < ncnst:
                tr0[0, 0, 0][n] = tr0_inout[0, 0, 0][n]
                n += 1

        ssthl0 = slope(
            kmask, thl0, thl0_above, thl0_below, pmid0, pmid0_above, pmid0_below
        )
        ssqt0 = slope(kmask, qt0, qt0_above, qt0_below, pmid0, pmid0_above, pmid0_below)
        ssu0 = slope(kmask, u0, u0_above, u0_below, pmid0, pmid0_above, pmid0_below)
        ssv0 = slope(kmask, v0, v0_above, v0_below, pmid0, pmid0_above, pmid0_below)

        # Calculate slope for each tracer by hand
        if dotransport == 1.0:
            n = 0
            while n < ncnst:
                sstr0[0, 0, 0][n] = slope(
                    kmask,
                    tr0[0, 0, -1][n],
                    tr0[0, 0, 0][n],
                    tr0[0, 0, -2][n],
                    pmid0,
                    pmid0_above,
                    pmid0_below,
                )
                n += 1

    with computation(PARALLEL), interval(...):
        pmid0 = pmid0_in
        exnmid0 = exnmid0_in
        qv0 = qv0_in
        ql0 = ql0_in
        qi0 = qi0_in
        qt0 = qv0 + ql0 + qi0
        t0 = th0_in * exnmid0
        thl0 = (
            t0
            - const.MAPL_ALHL * ql0 / const.MAPL_CP
            - const.MAPL_ALHS * qi0 / const.MAPL_CP
        ) / exnmid0
        zvir = 0.609  # r_H2O/r_air-1
        pifc0 = pifc0_in
        thl0bot = thl0 + ssthl0 * (pifc0 - pmid0)
        qt0bot = qt0 + ssqt0 * (pifc0 - pmid0)

        thj, qvj, qlj, qij, qse, id_check = conden(pifc0, thl0bot, qt0bot, ese, esx)

        if id_check == 1:
            id_check_flag = True

        thv0bot = thj * (1.0 + zvir * qvj - qlj - qij)
        thvl0bot = thl0bot * (1.0 + zvir * qt0bot)

        thl0top = thl0 + ssthl0 * (pifc0_in[0, 0, 1] - pmid0)
        qt0top = qt0 + ssqt0 * (pifc0_in[0, 0, 1] - pmid0)

    with computation(PARALLEL), interval(0, -1):
        thj, qvj, qlj, qij, qse, id_check = conden(
            pifc0_in[0, 0, 1], thl0top, qt0top, ese, esx
        )

        if id_check == 1:
            id_check_flag = True

        thv0top = thj * (1.0 + zvir * qvj - qlj - qij)
        thvl0top = thl0top * (1.0 + zvir * qt0top)

    with computation(PARALLEL), interval(-1, None):
        thv0top = thv0bot
        thvl0top = thvl0bot

    with computation(FORWARD), interval(...):
        """
        Below 'iter' loop is for implicit CIN closure

        It is important to note that this iterative cin loop is located at the outer  
        shell of the code. Thus, source air properties can also be changed during the 
        iterative cin calculation, because cumulus convection induces non-zero fluxes 
        even at interfaces below PBL top height through 'fluxbelowinv' subroutine. 
        """

        iteration = 0
        iter_max = 2
        while iteration < iter_max:
            """
            Cumulus scale height                                                    
            In contrast to the premitive code, cumulus scale height is iteratively 
            calculated at each time step, and at each iterative cin step.          
            It is not clear whether I should locate below two lines within or  out 
            of the iterative cin loop. 
            """

            """
            Find PBL top height interface index, 'kinv-1' where 'kinv' is the layer 
            index with PBLH in it. When PBLH is exactly at interface, 'kinv' is the 
            layer index having PBLH as a lower interface.                          
            In the previous code, I set the lower limit of 'kinv' by 2  in order to
            be consistent with the other parts of the code. However in the modified 
            code, I allowed 'kinv' to be 1 & if 'kinv = 1', I just exit the program 
            without performing cumulus convection. This new approach seems to be    
            more reasonable: if PBL height is within 'kinv=1' layer, surface is STL 
            interface (bflxs <= 0) and interface just above the surface should be   
            either non-turbulent (Ri>0.19) or stably turbulent (0<=Ri<0.19 but this 
            interface is identified as a base external interface of upperlying CL.  
            Thus, when 'kinv=1', PBL scheme guarantees 'bflxs <= 0'.  For this case 
            it is reasonable to assume that cumulus convection does not happen.     
            When these is SBCL, PBL height from the PBL scheme is likely to be very 
            close at 'kinv-1' interface, but not exactly, since 'zi' information is 
            changed between two model time steps. In order to ensure correct identi 
            fication of 'kinv' for general case including SBCL, I imposed an offset 
            of 5 [m] in the below 'kinv' finding block.
            """

            # Invert kpbl index
            if kpbl_in >= k0:
                kinv = k0 - kpbl_in + 1
            else:
                kinv = 5

            if kinv <= 1:
                id_exit = True

            # Cumulus convection goes here
            # This section has NOT BEEN PORTED
            # If id_exit == False the code will be triggered and an error will be raised

            iteration += 1

        # Initialize output variables when cumulus convection was not performed.
        umf_out[0, 0, 1] = 0.0
        dcm_out = 0.0
        qvten_out = 0.0
        qlten_out = 0.0
        qiten_out = 0.0
        sten_out = 0.0
        uten_out = 0.0
        vten_out = 0.0
        qrten_out = 0.0
        qsten_out = 0.0
        cufrc_out = 0.0
        cush_inout = -1.0
        qldet_out = 0.0
        qidet_out = 0.0
        fer_out = const.MAPL_UNDEF
        fdr_out = const.MAPL_UNDEF
        qtflx_out[0, 0, 1] = 0.0
        slflx_out[0, 0, 1] = 0.0
        uflx_out[0, 0, 1] = 0.0
        vflx_out[0, 0, 1] = 0.0

    with computation(PARALLEL), interval(0, 1):
        umf_out = 0.0
        qtflx_out = 0.0
        slflx_out = 0.0
        uflx_out = 0.0
        vflx_out = 0.0


class ComputeUwshcu:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ) -> None:
        """
        Initialize the ComputeUwshcu class.

        Parameters:
        stencil_factory (StencilFactory): Factory for creating stencil computations.
        quantity_factory (QuantityFactory): Factory for creating quantities.
        formulation: Saturation Formulation used for QSat
        """

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._compute_uwshcu = self.stencil_factory.from_dims_halo(
            func=compute_uwshcu,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

        self._k_mask = self.quantity_factory.zeros([X_DIM, Y_DIM, Z_DIM], "n/a")
        for i in range(0, self._k_mask.view[:].shape[0]):
            for j in range(0, self._k_mask.view[:].shape[1]):
                for k in range(0, self._k_mask.view[:].shape[2]):
                    if k == 0:
                        self._k_mask.view[i, j, k] = 1
                    elif k == 71:
                        self._k_mask.view[i, j, k] = 3
                    else:
                        self._k_mask.view[i, j, k] = 2

        self.id_exit = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )

        self.id_check_flag = self.quantity_factory.zeros(
            [X_DIM, Y_DIM, Z_DIM], "n/a", dtype=bool
        )

    @staticmethod
    def make_ntracers_quantity_factory(ijk_quantity_factory: QuantityFactory):
        ntracers_quantity_factory = copy.deepcopy(ijk_quantity_factory)
        ntracers_quantity_factory.set_extra_dim_lengths(
            **{
                "ntracers": const.NCNST,
            }
        )
        return ntracers_quantity_factory

    def __call__(
        self,
        dotransport: Float,
        ncnst: Int,
        k0: Int,
        kpbl_in: IntFieldIJ,
        pifc0_in: FloatField,
        pmid0_in: FloatField,
        exnmid0_in: FloatField,
        u0_in: FloatField,
        v0_in: FloatField,
        qv0_in: FloatField,
        ql0_in: FloatField,
        qi0_in: FloatField,
        th0_in: FloatField,
        tr0_inout: FloatField_NTracers,
        umf_out: FloatField,
        dcm_out: FloatField,
        qvten_out: FloatField,
        qlten_out: FloatField,
        qiten_out: FloatField,
        sten_out: FloatField,
        uten_out: FloatField,
        vten_out: FloatField,
        qrten_out: FloatField,
        qsten_out: FloatField,
        cufrc_out: FloatField,
        fer_out: FloatField,
        fdr_out: FloatField,
        qldet_out: FloatField,
        qidet_out: FloatField,
        qlsub_out: FloatField,
        qisub_out: FloatField,
        ndrop_out: FloatField,
        nice_out: FloatField,
        tpert_out: FloatFieldIJ,
        qpert_out: FloatFieldIJ,
        qtflx_out: FloatField,
        slflx_out: FloatField,
        uflx_out: FloatField,
        vflx_out: FloatField,
        tr0_test: FloatField_NTracers,
        ssthl0_test: FloatField,
        ssqt0_test: FloatField,
        ssu0_test: FloatField,
        ssv0_test: FloatField,
        sstr0_test: FloatField_NTracers,
        thj_test: FloatField,
        qvj_test: FloatField,
        qlj_test: FloatField,
        qij_test: FloatField,
        qse_test: FloatField,
        id_check_test: IntField,
        formulation: SaturationFormulation = SaturationFormulation.Staars,
    ):
        self.qsat = QSat(
            self.stencil_factory,
            self.quantity_factory,
            formulation=formulation,
        )

        self.id_exit.view[:, :, :] = False
        self.id_check_flag.view[:, :, :] = False

        self._compute_uwshcu(
            dotransport=dotransport,
            ncnst=ncnst,
            k0=k0,
            kpbl_in=kpbl_in,
            pifc0_in=pifc0_in,
            pmid0_in=pmid0_in,
            exnmid0_in=exnmid0_in,
            u0_in=u0_in,
            v0_in=v0_in,
            qv0_in=qv0_in,
            ql0_in=ql0_in,
            qi0_in=qi0_in,
            th0_in=th0_in,
            tr0_inout=tr0_inout,
            umf_out=umf_out,
            dcm_out=dcm_out,
            qvten_out=qvten_out,
            qlten_out=qlten_out,
            qiten_out=qiten_out,
            sten_out=sten_out,
            uten_out=uten_out,
            vten_out=vten_out,
            qrten_out=qrten_out,
            qsten_out=qsten_out,
            cufrc_out=cufrc_out,
            fer_out=fer_out,
            fdr_out=fdr_out,
            qldet_out=qldet_out,
            qidet_out=qidet_out,
            qlsub_out=qlsub_out,
            qisub_out=qisub_out,
            ndrop_out=ndrop_out,
            nice_out=nice_out,
            tpert_out=tpert_out,
            qpert_out=qpert_out,
            qtflx_out=qtflx_out,
            slflx_out=slflx_out,
            uflx_out=uflx_out,
            vflx_out=vflx_out,
            tr0=tr0_test,
            ssthl0=ssthl0_test,
            ssqt0=ssqt0_test,
            ssu0=ssu0_test,
            ssv0=ssv0_test,
            sstr0=sstr0_test,
            thj=thj_test,
            qvj=qvj_test,
            qlj=qlj_test,
            qij=qij_test,
            qse=qse_test,
            id_check=id_check_test,
            kmask=self._k_mask,
            ese=self.qsat.ese,
            esx=self.qsat.esx,
            id_exit=self.id_exit,
            id_check_flag=self.id_check_flag,
        )

        if (self.id_exit.view[:, :, :] == False).all():
            raise NotImplementedError(
                "Expected id_exit == True, got id_exit == False! "
                "Cumulus convection was unexpectedly triggered! "
                "This code has not been ported."
            )

        if (self.id_check_flag.view[:, :, :] == True).all():
            raise NotImplementedError(
                "Expected id_check == 0, got id_check == 1! "
                "Code has been triggered that has not been ported."
            )
