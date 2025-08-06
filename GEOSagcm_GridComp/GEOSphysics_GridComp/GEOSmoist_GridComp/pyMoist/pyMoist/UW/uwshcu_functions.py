import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import (
    K,
    erfc,
    exp,
    float32,
    float64,
    log,
    sin,
    sqrt,
    int32,
)

import pyMoist.constants as constants
import pyMoist.pyMoist_constants as py_constants
from ndsl.dsl.typing import Float, FloatField, Int
from pyMoist.saturation.qsat import FloatField_Extra_Dim, QSat_Float


P00 = Float(1e5)  # Reference pressure
zvir = Float(0.609)  # r_H2O/r_air-1
ROVCP = constants.MAPL_RGAS / constants.MAPL_CP  # Gas constant over specific heat


@gtscript.function
def exnerfn(
    p: Float,
) -> Float:
    """
    Function that calculates the Exner function for a given pressure.

    Inputs:
    p [Float]: Atmospheric pressure [Pa]

    Returns:
    (p / 100000.0) ** (constants.MAPL_RGAS / constants.MAPL_CP) [Float]:
    Exner function [unitless]
    """

    return (p / 100000.0) ** (constants.MAPL_RGAS / constants.MAPL_CP)


@gtscript.function
def slope_bot(
    field: FloatField,
    p0: FloatField,
):
    """
    Function that calculates slope at bottom layer of a field.

    Inputs:
    field (FloatField): Field of interest [N/A]
    p0 (FloatField): Pressure [Pa]

    Returns:
    slope (Float): Slope of the field of interest [N/A]
    """
    if K == 0:
        value = (field[0, 0, 1] - field) / (p0[0, 0, 1] - p0)
        if value > 0.0:
            slope = max(0.0, value)
        else:
            slope = min(0.0, value)

    return slope


@gtscript.function
def slope_mid(
    max_k: Int,
    field: FloatField,
    p0: FloatField,
):
    """
    Function that calculates slope at mid layers of a field.

    Inputs:
    max_k [Int]: Max k level (e.g., 71)
    field [FloatField]: Field of interest [n/a]
    p0 [FloatField]: Pressure [Pa]

    Returns:
    slope [Float]: Slope of the field of interest [n/a]
    """

    if K > 0 and K < max_k:
        above_value = (field[0, 0, 1] - field) / (p0[0, 0, 1] - p0)
        below_value = (field - field[0, 0, -1]) / (p0 - p0[0, 0, -1])
        if above_value > 0.0:
            slope = max(0.0, min(above_value, below_value))
        else:
            slope = min(0.0, max(above_value, below_value))

    return slope


# @gtscript.function
# def slope_top(
#     max_k: Int,
#     field: FloatField,
#     p0: FloatField,
# ):
#     """
#     Function that calculates slope at mid layers of a field.

#     Inputs:
#     max_k [Int]: Max k level (e.g., 71)
#     field [FloatField]: Field of interest [n/a]
#     p0 [FloatField]: Pressure [Pa]

#     Returns:
#     slope [Float]: Slope of the field of interest [n/a]
#     """

#     if K == max_k:
#         above_value = (field[0, 0, -1] - field) / (p0[0, 0, -1] - p0)
#         below_value = (field - field[0, 0, -2]) / (p0 - p0[0, 0, -2])
#         if above_value > 0.0:
#             slope = max(0.0, min(above_value, below_value))
#         else:
#             slope = min(0.0, max(above_value, below_value))

#     return slope


@gtscript.function
def ice_fraction(
    temp: Float,
    cnv_frc: Float,
    srf_type: Float,
):
    """
    Function description

    Inputs:
    temp [Float]: Temperature [K]
    cnv_frc [Float]: Convective fraction [?]
    srf_type [Float]: Surface type [?]

    Returns:
    ice_frac [Float]: Ice fraction [unitless]
    """
    # Anvil clouds
    # Anvil-Convective sigmoidal function like figure 6(right)
    # Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
    if temp <= constants.JaT_ICE_ALL:
        icefrct_c = 1.000
    elif temp > constants.JaT_ICE_ALL and temp <= constants.JaT_ICE_MAX:
        icefrct_c = sin(
            0.5
            * constants.MAPL_PI
            * (
                1.00
                - (temp - constants.JaT_ICE_ALL)
                / (constants.JaT_ICE_MAX - constants.JaT_ICE_ALL)
            )
        )
    else:
        icefrct_c = 0.00
    icefrct_c = max(min(icefrct_c, 1.00), 0.00) ** constants.aICEFRPWR
    # Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if srf_type == 2.0:
        if temp <= constants.JiT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.JiT_ICE_ALL and temp <= constants.JiT_ICE_MAX:
            icefrct_m = 1.00 - (temp - constants.JiT_ICE_ALL) / (
                constants.JiT_ICE_MAX - constants.JiT_ICE_ALL
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.iICEFRPWR
    elif srf_type > 1.0:
        if temp <= constants.lT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.lT_ICE_ALL and temp <= constants.lT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * constants.MAPL_PI
                * (
                    1.00
                    - (temp - constants.lT_ICE_ALL)
                    / (constants.lT_ICE_MAX - constants.lT_ICE_ALL)
                )
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.lICEFRPWR
    else:
        if temp <= constants.oT_ICE_ALL:
            icefrct_m = 1.000
        elif temp > constants.oT_ICE_ALL and temp <= constants.oT_ICE_MAX:
            icefrct_m = sin(
                0.5
                * constants.MAPL_PI
                * (
                    1.00
                    - (temp - constants.oT_ICE_ALL)
                    / (constants.oT_ICE_MAX - constants.oT_ICE_ALL)
                )
            )
        else:
            icefrct_m = 0.00
        icefrct_m = max(min(icefrct_m, 1.00), 0.00) ** constants.oICEFRPWR
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
    p [Float]: Pressure [Pa]
    thl [Float]: Liquid potential temperature [K]
    qt [Float]: Mixing ratio [kg/kg]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]

    Returns:
    th [Float]: Temperature [K]
    qv [Float]: Water vapor mixing ratio [kg/kg]
    ql [Float]: Liquid water mixing ratio [kg/kg]
    qi [Float]: Ice water mixing ratio [kg/kg]
    rvls [Float]: Saturation specific humidity [kg/kg]
    id_check (Int): Flag that indicates if condensation occurs
    (0 for no condensation, 1 for condensation).
    """

    tc: float64 = float32(thl) * exnerfn(p)

    nu: float64 = ice_fraction(float32(tc), 0.0, 0.0)
    leff: float64 = (float64(1.0) - nu) * constants.MAPL_LATENT_HEAT_VAPORIZATION + (
        nu * constants.MAPL_LATENT_HEAT_SUBLIMATION
    )
    temps: float32 = tc
    ps: float32 = p
    qs, _ = QSat_Float(ese, esx, temps, ps / 100.0)
    rvls: float32 = float64(qs)

    if qs >= qt:  # no condensation
        id_check = 0
        qv: float32 = qt
        qc: float64 = 0.0
        ql: float32 = 0.0
        qi: float32 = 0.0
        th: float32 = thl
    else:  # condensation
        iteration = 0
        while iteration < 10:
            temps = temps + ((tc - temps) * constants.MAPL_CP / leff + qt - rvls) / (
                constants.MAPL_CP / leff
                + constants.EPSILON
                * leff
                * rvls
                / (constants.MAPL_RGAS * temps * temps)
            )
            qs, _ = QSat_Float(ese, esx, temps, ps / 100.0)
            rvls = qs
            iteration += 1
        qc: float64 = max(qt - qs, float64(0.0))
        qv = qt - qc
        ql = qc * (float64(1.0) - nu)
        qi = nu * qc
        th = temps / exnerfn(p)
        if abs((temps - (leff / constants.MAPL_CP) * qc) - tc) >= float64(1.0):
            id_check = 1
        else:
            id_check = 0

    return (
        float32(th),
        float32(qv),
        float32(ql),
        float32(qi),
        float32(rvls),
        int32(id_check),
    )


@gtscript.function
def compute_alpha(
    del_CIN: Float,
    ke: Float,
):
    """
    Subroutine to compute proportionality factor for
    implicit CIN calculation.

    Inputs:
    del_CIN [Float]: Difference between initial and final CIN calculations [J/kg]
    ke [Float]: Evaporative efficiency [?]

    Returns:
    compute_alpha [Float]: Proportionality factor for CIN calculation [unitless]
    """

    x0: float64 = float64(0.0)
    del_CIN8_f64: float64 = float64(del_CIN)
    ke8_f64: float64 = ke
    iteration = 0
    while iteration < 10:
        x1 = x0 - (exp(-x0 * ke8_f64 * del_CIN8_f64) - x0) / (
            -ke8_f64 * del_CIN8_f64 * exp(-x0 * ke8_f64 * del_CIN8_f64) - 1.0
        )
        x0 = x1
        iteration += 1

    compute_alpha = float32(x0)

    return compute_alpha


@gtscript.function
def compute_mumin2(
    mulcl: Float,
    rmaxfrax: Float,
    mulow: Float,
):
    """
    Subroutine to compute critical 'mu' (normalized CIN) such
    that updraft fraction at the LCL is equal to 'rmaxfrac'.

    Inputs:
    mulcl [Float]: Some var at the LCL [?]
    rmaxfrac [Float]: Maximum core updraft fraction [unitless]
    mulow [Float]: Some var at the bottom interface [?]

    Returns:
    compute_mumin2 [Float]: Critical mu (normalized CIN) [unitless]

    """

    x0: float64 = mulow
    iteration = 0
    while iteration < 10:
        ex: float64 = exp(-(x0**2))
        ef: float64 = erfc(x0)  # Complimentary error fraction function
        exf: float64 = ex / ef
        f: float64 = (
            float64(0.5) * exf**2
            - float64(0.5) * (ex / float64(2.0) / rmaxfrax) ** 2
            - (mulcl * float64(2.5066) / float64(2.0)) ** 2
        )
        fs: float64 = (float64(2.0) * exf**2) * (
            exf / sqrt(constants.MAPL_PI) - x0
        ) + (float64(0.5) * x0 * ex**2) / (rmaxfrax**2)
        x1: float64 = x0 - f / fs
        x0 = x1
        iteration += 1

    compute_mumin2 = float32(x0)

    return compute_mumin2


@gtscript.function
def compute_ppen(
    wtwb: Float,
    drag: Float,
    bogbot: Float,
    bogtop: Float,
    rho0j: Float,
    dpen: Float,
):
    """
    Function to compute critical 'ppen[Pa]<0' ( pressure dis.
    from 'pifc0(kpen-1)' to the cumulus top where cumulus updraft
    vertical velocity is exactly zero ) by considering exact
    non-zero fer(kpen).

    Inputs:
    wtwb [Float]: Updraft vertical velocity at lower interface [m/s]
    drag [Float]: Drag coefficient [unitless]
    bogbot [Float]: Cloud buoyancy at base interface [n/a]
    bogtop [Float]: Cloud buoyancy at top interface [n/a]
    rho0j [Float]: Density of water [kg/m^3] [?]
    dpen [Float]: Environmental layer pressure thickness [Pa] > 0

    Returns:
    compute_ppen [Float]: Critical ppen [Pa]
    """

    # Buoyancy slope
    SB: float64 = (bogtop - bogbot) / dpen

    # Sign of slope, 'f' at x = 0
    # If 's00>0', 'w' increases with height.
    s00: float64 = bogbot / rho0j - drag * wtwb

    if drag * dpen < float64(1.0e-4):
        if s00 >= float64(0.0):
            x0: float64 = dpen
        else:
            x0 = max(float64(0.0), min(dpen, float64(-0.5) * wtwb / s00))
    else:
        if s00 >= float64(0.0):
            x0 = dpen
        else:
            x0 = float64(0.0)

        iteration = 0
        while iteration < 5:
            aux: float64 = min(max(float64(-2.0) * drag * x0, -20.0), 20.0)

            f: float64 = exp(aux) * (
                wtwb - (bogbot - SB / (2.0 * drag)) / (drag * rho0j)
            ) + (SB * x0 + bogbot - SB / (2.0 * drag)) / (drag * rho0j)
            fs: float64 = -2.0 * drag * exp(aux) * (
                wtwb - (bogbot - SB / (2.0 * drag)) / (drag * rho0j)
            ) + (SB) / (drag * rho0j)

            x1: float64 = x0 - f / fs
            x0 = x1
            iteration += 1

    compute_ppen = -max(float64(0.0), min(dpen, x0))

    return compute_ppen


@gtscript.function
def getbuoy(
    pbot: Float,
    thv0bot: Float,
    ptop: Float,
    thv0top: Float,
    thvubot: Float,
    thvutop: Float,
    cin_in: Float,
    plfc_in: Float,
):
    """
    Function to calculate integrated CIN [ J/kg = m2/s2 ] and
    'cinlcl, plfc' if any. Assume 'thv' is linear in each layer
    both for cumulus and environment. Note that this subroutine
    only includes positive CIN in calculation - if there is any
    negative CIN, it is assumed to be zero. This is slightly
    different from 'single_cin' below, where both positive  and
    negative CIN are included.

    Inputs:
    pbot [Float]: Pressure at bottom layer [Pa]
    thv0bot [Float]: Some sort of temperature at bottom [?]
    ptop [Float]: Pressure at top layer [Pa]
    thv0top [Float]: Some sort of temperature at top [?]
    thvubot [Float]: Some sort of temperature at bot [?]
    thvutop [Float]: Some sort of temperature at top [?]
    cin_in [Float]: Convective inhibition [J/kg]
    plfc_in [Float]: Pressure at the level of free convection [Pa]

    Returns:
    plfc [Float]: Pressure at level of free convection [Pa]
    cin [Float]: Integreated CIN [J/kg]
    """
    plfc = plfc_in
    cin = cin_in

    if thvubot > thv0bot and thvutop > thv0top:
        plfc = pbot
    elif thvubot <= thv0bot and thvutop <= thv0top:
        cin = cin_in - ((thvubot / thv0bot - 1.0) + (thvutop / thv0top - 1.0)) * (
            pbot - ptop
        ) / (
            pbot / (constants.MAPL_RGAS * thv0bot * exnerfn(pbot))
            + ptop / (constants.MAPL_RGAS * thv0top * exnerfn(ptop))
        )
    elif thvubot > thv0bot and thvutop <= thv0top:
        frc = (thvutop / thv0top - 1.0) / (
            (thvutop / thv0top - 1.0) - (thvubot / thv0bot - 1.0)
        )
        cin = cin_in - (thvutop / thv0top - 1.0) * (
            (ptop + frc * (pbot - ptop)) - ptop
        ) / (
            pbot / (constants.MAPL_RGAS * thv0bot * exnerfn(pbot))
            + ptop / (constants.MAPL_RGAS * thv0top * exnerfn(ptop))
        )
    else:
        frc = (thvubot / thv0bot - 1.0) / (
            (thvubot / thv0bot - 1.0) - (thvutop / thv0top - 1.0)
        )
        plfc = pbot - frc * (pbot - ptop)
        cin = cin_in - ((thvubot / thv0bot - 1.0) * (pbot - plfc)) / (
            (
                pbot / (constants.MAPL_RGAS * thv0bot * exnerfn(pbot))
                + ptop / (constants.MAPL_RGAS * thv0top * exnerfn(ptop))
            )
        )

    return plfc, cin  # Note: plfc and cin are returned, but not always used


@gtscript.function
def qsinvert(
    qt: Float,
    thl: Float,
    ps_in: Float,
    ese: FloatField_Extra_Dim,
    esx: FloatField_Extra_Dim,
):
    """
    Function calculating saturation pressure ps (or pLCL) from qt and
    thl ( liquid potential temperature,  NOT liquid virtual potential
    temperature) by inverting Bolton formula. I should check later if
    current use of 'leff' instead of 'xlv' here is reasonable or not.

    Inputs:
    qt [Float]: Mixing ratio [kg/kg]
    thl [Float]: Liquid potential temperature [K]
    ps_in [Float]: Pressure [Pa]
    ese [FloatField_Extra_Dim]: Used in QSat_Float [n/a]
    esx [FloatField_Extra_Dim]: Used in QSat_Float [n/a]

    Returns:
    qsinvert [Float]: Saturation pressure [Pa]
    """

    psmin: float64 = float64(
        10000.0
    )  # Default saturation pressure [Pa] if iteration does not converge
    dpsmax: float64 = float64(1.0)  # Tolerance [Pa] for convergence of iteration
    p00 = 1e5
    rovcp = constants.MAPL_RDRY / constants.MAPL_CP

    # Calculate best initial guess of pLCL
    Ti: float64 = thl * (ps_in / p00) ** rovcp
    Tgeos: float32 = Ti
    Pgeos: float32 = float32(ps_in)
    qs, dqsdT = QSat_Float(ese, esx, Tgeos, Pgeos / 100.0)
    es: float64 = (
        ps_in
        * qs
        / (py_constants.ep2 + (float64(1.0) - py_constants.ep2) * float64(qs))
    )
    rhi: float64 = qt / float64(qs)

    if rhi <= float64(0.01):
        qsinvert: float32 = psmin

    else:

        TLCL: float64 = float64(55.0) + float64(1.0) / (
            float64(1.0) / (Ti - float64(55.0)) - log(rhi) / float64(2840.0)
        )  # Bolton's formula. MWR.1980.Eq.(22)
        PiLCL: float64 = TLCL / thl
        ps: float64 = p00 * (PiLCL) ** (float64(1.0) / rovcp)

        iteration = 0
        while iteration < 10:
            Pis: float64 = (ps / p00) ** rovcp  # Exner function
            Ts: float64 = thl * Pis
            Tgeos = Ts
            Pgeos = ps
            qs, dqsdT = QSat_Float(ese, esx, Tgeos, Pgeos / 100.0, DQSAT_trigger=True)
            gam: float64 = (
                constants.MAPL_LATENT_HEAT_VAPORIZATION / constants.MAPL_CP
            ) * float64(dqsdT)
            err: float64 = qt - qs
            nu: float64 = ice_fraction(float32(Ts), 0.0, 0.0)
            leff: float64 = (
                (float64(1.0) - nu) * constants.MAPL_LATENT_HEAT_VAPORIZATION
                + nu * constants.MAPL_LATENT_HEAT_SUBLIMATION
            )
            dlnqsdT: float64 = gam * (constants.MAPL_CP / leff) / qs
            dTdPis: float64 = thl
            dPisdps: float64 = rovcp * Pis / ps
            dlnqsdps: float64 = float64(-1.0) / (ps - (1.0 - py_constants.ep2) * es)
            derrdps: float64 = -qs * (dlnqsdT * dTdPis * dPisdps + dlnqsdps)
            dps: float64 = -err / derrdps
            ps = ps + dps

            if ps < float64(0.0):
                qsinvert = psmin
                iteration = 10

            elif abs(dps) <= dpsmax:
                qsinvert = ps
                iteration = 10

            else:
                qsinvert = psmin

            iteration += 1

    return float32(qsinvert)


@gtscript.function
def sign(
    a: Float,
    b: Float,
):
    """
    Function that returns the magnitude of one argument and the sign of another.

    Inputs:
    a [Float]: Argument of which the magnitude is needed [unitless]
    b [Float]: Argument of which the sign is needed [unitless]

    Returns:
    result [Float]: The magnitude of a and sign of b [unitless]
    """

    if b >= 0.0:
        result = abs(a)
    else:
        result = -abs(a)

    return result


@gtscript.function
def roots(
    a: Float,
    b: Float,
    c: Float,
):
    """
    Function to solve a second order polynomial equation of the
    form [ax^2 + bx + c].

    Inputs:
    a [Float]: Coefficient of the x^2 term [unitless]
    b [Float]: Coefficient of x [unitless]
    c [Float]: Constant term [unitless]

    Returns:
    r1 [Float]: The first root of the polynomial [unitless]
    r2  [Float]: The second root of the polynomial [unitless]
    status [Int]: 0 if roots are found. 1, 2, or 3 if
    there are no roots [unitless]
    """

    status = 0

    if a == 0:  # Form b*x + c = 0
        if b == 0:  # Failure: c = 0
            status = 1
        else:  # b*x + c = 0
            r1 = -c / b
        r2 = r1
    else:
        if b == 0:  # Form a*x**2 + c = 0
            if a * c > 0:  # Failure: x**2 = -c/a < 0
                status = 2
            else:  # x**2 = -c/a
                r1 = sqrt(-c / a)
            r2 = -r1
        else:  # Form a*x**2 + b*x + c = 0
            if (b**2 - 4.0 * a * c) < 0.0:  # Failure, no real roots
                status = 3
            else:
                q = -0.5 * (b + sign(1.0, b) * sqrt(b**2 - 4.0 * a * c))
                r1 = q / a
                r2 = c / q

    return r1, r2, status


@gtscript.function
def single_cin(
    pbot: Float,
    thv0bot: Float,
    ptop: Float,
    thv0top: Float,
    thvubot: Float,
    thvutop: Float,
):
    """
    Function to calculate a single layer CIN by summing all
    positive and negative CIN.

    Inputs:
    pbot [Float]: Pressure at bottom layer [Pa]
    thv0bot [Float]: Some sort of temperature at bottom layer [?]
    ptop [Float]: Pressure at top of layer [Pa]
    thv0top [Float]: Some sort of temperature at top layer [?]
    thvubot [Float]: Some sort of temperature at bottom layer [?]
    thvutop [Float]: Some sort of temperature at top layer [?]

    Returns:
    single_cin [Float]: Convective Inhibition (CIN) of a single layer [J/kg]
    """

    single_cin = (
        ((1.0 - thvubot / thv0bot) + (1.0 - thvutop / thv0top))
        * (pbot - ptop)
        / (
            pbot / (constants.MAPL_RGAS * thv0bot * exnerfn(pbot))
            + ptop / (constants.MAPL_RGAS * thv0top * exnerfn(ptop))
        )
    )

    return single_cin
