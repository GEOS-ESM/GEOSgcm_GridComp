from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float
from ndsl import Quantity, QuantityFactory, StencilFactory
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log, exp, sqrt
import pyMoist.aer_activation_constants as constants

import numpy as np

# Global space
FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (constants.n_modes))]

'''
def Erf(
    x: Float
):
    """
    Compute the error function for a given value x.
    
    Parameters:
    x (Float): Input value for which the error function is to be computed.
    
    Returns:
    Float: The error function value for the input x.
    """
    erf = 0.0
    if x < 0.0e+00:
        return -GammP(0.5, x**2)
    else:
        return GammP(0.5, x**2)

def GammLn(
    xx: Float
    )-> Float:
    """
    Compute the natural logarithm of the gamma function for a given value xx.
    
    Parameters:
    xx (Float): Input value for which the natural logarithm of the gamma function is to be computed.
    
    Returns:
    Float: The natural logarithm of the gamma function value for the input xx.
    """
    cof = [
        76.18009172947146, 
        -86.50532032941677, 
        24.01409824083091, 
        -1.231739572450155, 
        0.001208650973866179, 
        -0.000005395239384953
    ]
    stp = 2.5066282746310005
    
    x = xx
    y = x
    tmp = x + 5.5
    tmp = (x + 0.5) * log(tmp) - tmp
    ser = 1.000000000190015
    
    for j in range(6):
        y += 1.0
        ser += cof[j] / y
    
    gammln = tmp + log(stp * ser / x)
    return gammln
    
def GammP(
    a: Float,
    x: Float 
)-> Float:
    """
    Compute the incomplete gamma function for given values a and x.
    
    Parameters:
    a (Float): Parameter a for the incomplete gamma function.
    x (Float): Parameter x for the incomplete gamma function.
    
    Returns:
    Float: The incomplete gamma function value for the input parameters a and x.
    """
    #Not sure what these variables are
    if x < 0.0 or a <= 0.0:
        raise ValueError("aero_actv: function gammp: bad arguments")
    gamser = 0.0
    gammcf = 0.0
    gln = GammLn(a)
    if x < a + 1.0:
        _gser_stencil(gamser, a, x, gln)
        gammp = gamser
    else:
        _gcf_matrix_stencil(gammcf, a, x, gln)
        gammp = 1.0 - gammcf
    return gammp

def _get_act_frac_stencil(
    nmodes: Int, 
    xnap: FloatField, 
    rg: FloatField, 
    sigmag: FloatField, 
    tkelvin: Float, 
    ptot: Float, 
    wupdraft: Float, 
    nact: FloatField, 
    bibar: FloatField,
):
    """
    Compute the activation fraction stencil for aerosol properties.
    
    Parameters:
    nmodes (Int): Number of modes.
    xnap (FloatField): Number concentration for each mode.
    rg (FloatField): Geometric mean dry radius for each mode.
    sigmag (FloatField): Geometric standard deviation for each mode.
    tkelvin (Float): Absolute temperature in Kelvin.
    ptot (Float): Total ambient pressure in Pa.
    wupdraft (Float): Updraft velocity in m/s.
    nact (FloatField): Activation number concentration for each mode.
    bibar (FloatField): Intermediate computation field.
    
    Returns:
    None
    """
    with computation(PARALLEL), interval(...):
        # Call to _act_frac_mat_stencil
        _act_frac_mat_stencil(nmodes, xnap, rg, sigmag, bibar, tkelvin, ptot, wupdraft, nact)

@gtscript.function
def act_frac_mat_loop(sm, eta, ac, rg, fracactn, nact, sigmag, bibar, dum, gamma, xnap, zeta):
        #These variables must be computed for each mode
        n = 0
        xlogsigm = log(sigmag[0,0,0][n])
        smax = 0.0
        while n < constants.n_modes:
            sm[0,0,0][n] = (2.0 / sqrt(bibar[0,0,0][n])) * (a / (3.0 * rg[0,0,0][n])) ** 1.5
            eta[0,0,0][n] = dum ** 3 / (constants.TWOPI * constants.DENH2O * gamma * xnap[0,0,0][n])
            f1 = 0.5 * exp(2.50 * xlogsigm[0,0,0][n] ** 2)
            f2 = 1.0 + 0.25 * xlogsigm[0,0,0][n]
            smax = smax + (f1*(zeta/eta[0,0,0][n])**1.5+ f2*(sm[0,0,0][n]**2/(eta[0,0,0][n]+3.0*zeta))**0.75)/sm[0,0,0][n]**2
            n +=1
        
        smax = 1.0e+00 / sqrt(smax)  
        n=0
        while n < constants.n_modes:
            #lines 534
            ac[0,0,0][n] = rg[0,0,0][n] * (sm[0,0,0][n] / smax) ** 0.66666666666666667
            u = log(ac[0,0,0][n] / rg[0,0,0][n]) / (constants.SQRT2 * xlogsigm[0,0,0][n])
            fracactn[0,0,0][n] = 0.5 * (1.0 - Erf(u))
            nact[0,0,0][n] = fracactn[0,0,0][n] * xnap[0,0,0][n]
            n+=1

def _act_frac_mat_stencil(
    n_modes: Int,
    xnap: Int,
    rg: Float,
    sigmag: Float,
    bibar: Float,
    tkelvin: Float,
    ptot: Float,
    wupdraft: FloatField, 
    nact: FloatField,
):
    
    """
    Perform matrix computations for the activation fraction stencil.
    
    Parameters:
    n_modes (Int): Number of modes.
    xnap (Int): Number concentration for each mode.
    rg (Float): Geometric mean dry radius for each mode.
    sigmag (Float): Geometric standard deviation for each mode.
    bibar (Float): Intermediate computation field.
    tkelvin (Float): Absolute temperature in Kelvin.
    ptot (Float): Total ambient pressure in Pa.
    wupdraft (Float): Updraft velocity in m/s.
    nact (FloatField): Activation number concentration for each mode.
    
    Returns:
    None
    """
    with computation(PARALLEL), interval(...):
        
        #rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
        #values close to those given in a-z et al. 1998 figure 5. 
        
        rdrp = 0.105e-06 #[m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.
        sm = np.zeros(constants.n_modes) 
        eta = np.zeros(constants.n_modes)
        ac = np.zeros(constants.n_modes)
        fracactn = np.zeros(constants.n_modes)
        #These variables are common to all modes and need only be computed once. 
        dv = constants.DIJH2O0 * (constants.P0DIJ / ptot) * (tkelvin / constants.T0DIJ) ** 1.94e+00 #[m^2/s] (p&k,2nd ed., p.503)
        surten = 76.10e-3 - 0.155e-3 * (tkelvin - 273.15e+00) #[j/m^2]
        wpe = exp(77.34491296 - 7235.424651 / tkelvin - 8.2 * log(tkelvin) + tkelvin * 5.7113e-3) #[pa]
        dumw = sqrt(constants.TWOPI * constants.WMOLMASS / constants.RGASJMOL / tkelvin) #[s/m]
        dvprime = dv / ((rdrp / (rdrp + constants.DELTAV)) + (dv * dumw / (rdrp * constants.ALPHAC)))
        xka = (5.69 + 0.017 * (tkelvin - 273.15)) * 418.4e-5
        duma = sqrt(constants.TWOPI * constants.AMOLMASS / constants.RGASJMOL / tkelvin)
        xkaprime = xka / ((rdrp / (rdrp + constants.DELTAT)) + (xka * duma / (rdrp * constants.ALPHAT * constants.DENH2O * constants.CPAIR)))
        g = 1.0 / ((constants.DENH2O * constants.RGASJMOL * tkelvin) / (wpe * dvprime * constants.WMOLMASS) + 
                    ((constants.HEATVAP * constants.DENH2O) / (xkaprime * tkelvin)) * 
                    ((constants.HEATVAP * constants.WMOLMASS) / (constants.RGASJMOL * tkelvin) - 1.0))
        a = (2.0 * surten * constants.WMOLMASS) / (constants.DENH2O * constants.RGASJMOL * tkelvin)
        alpha = (constants.GRAVITY / (constants.RGASJMOL * tkelvin)) * ((constants.WMOLMASS * constants.HEATVAP) / (constants.CPAIR * tkelvin) - constants.AMOLMASS)
        gamma = (constants.RGASJMOL * tkelvin) / (wpe * constants.WMOLMASS) + (constants.WMOLMASS * constants.HEATVAP * constants.HEATVAP) / (constants.CPAIR * ptot * constants.AMOLMASS * tkelvin)
        dum = sqrt(alpha * wupdraft / g)
        zeta = 2.0 * a * dum / 3.0

        #These variables must be computed for each mode
        act_frac_mat_loop(sm, eta, ac, rg, fracactn, nact, sigmag, bibar, dum, gamma, xnap, zeta)

def _gser_stencil(
    gamser: FloatField, 
    a: Float, 
    x: Float, 
    gln: FloatField,
):
    """
    Compute the series representation of the incomplete gamma function.
    
    Parameters:
    gamser (Float): Output value of the series representation of the incomplete gamma function.
    a (Float): Parameter a for the incomplete gamma function.
    x (Float): Parameter x for the incomplete gamma function.
    gln (Float): Natural logarithm of the gamma function.
    
    Returns:
    gamser (Float): Output value of the series representation of the incomplete gamma function.er 
    """
    with computation(PARALLEL), interval(...):
        eps = 3.0e-9
        itmax = 10000
        if x <= 0:
            if x < 0:
                raise ValueError("x < 0 in gser")
            gamser = 0.0
        ap = a
        sum_ = 1.0 / a
        del_ = sum_
        n = 1
        while n <= itmax:
            ap += 1
            del_ *= x / ap
            sum_ += del_
            if abs(del_) < abs(sum_) * eps:
                gamser = sum_ * exp(-x + a * log(x) - gln)
                return gamser
            n += 1
        gamser = sum_ * exp(-x + a * log(x) - gln)
        return gamser


@gtscript.function
def gcf_matrix_computation(itmax: Int,
                            eps: Float, 
                            fpmin: Float, 
                            b: Float, 
                            c: Float,  
                            d: Float, 
                            h: Float,
):
    i = 1
    while i <= itmax:
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        del_ = d * c
        h *= del_
        if abs(del_ - 1.0) < eps:
            return b, c, d, h
        i += 1
    return b, c, d, h

def _gcf_matrix_stencil(
    gammcf: FloatField, 
    a: Float, 
    x: Float, 
    gln: FloatField,
):
    """
    Compute the continued fraction representation of the incomplete gamma function.
    
    Parameters:
    gammcf (Float): Output value of the continued fraction representation of the incomplete gamma function.
    a (Float): Parameter a for the incomplete gamma function.
    x (Float): Parameter x for the incomplete gamma function.
    gln (Float): Natural logarithm of the gamma function.
    
    Returns:
    None
    """
    with computation(PARALLEL), interval(...):
        itmax = 10000
        eps = 3.0e-7
        fpmin = 1.0e-30
        b = x + 1.0 - a
        c = 1.0 / fpmin
        d = 1.0 / b
        h = d
        gcf_matrix_computation(itmax, eps, fpmin, b, c, d, h)
        gammcf = exp(-x + a * log(x) - gln) * h

'''


def aer_activation_stencil(
    q: FloatField,
    t: FloatField,
    plo: FloatField,
    ple: FloatField,
    zlo: FloatField,
    zle: FloatField,
    qlcn: FloatField,
    qicn: FloatField,
    qlls: FloatField,
    qils: FloatField,
    sh: FloatField,
    evap: FloatField,
    kpbl: FloatField,
    tke: FloatField,
    vvel: FloatField,
    FRLAND: FloatField,
    NACTL: FloatField,
    NACTI: FloatField,
    NWFA: FloatField,
    NN_LAND: Float,
    NN_OCEAN: Float,
):
    """
    Perform aerosol activation computations based on input atmospheric and aerosol properties.

    Parameters:
    q (Floatfield): Specific humidity field.
    t (Floatfield): Temperature field.
    plo (Floatfield): Low-level pressure field.
    ple (Floatfield): Low-level pressure field at the end of the time step.
    zlo (Floatfield): Low-level height field.
    zle (Floatfield): Low-level height field at the end of the time step.
    qlcn (Floatfield): Cloud liquid water mixing ratio field.
    qicn ((Floatfield): Cloud ice mixing ratio field.
    qlls (Floatfield): Large-scale cloud liquid water mixing ratio field.
    qils (Floatfield): Large-scale cloud ice mixing ratio field.
    sh (Floatfield): Specific humidity field.
    evap (Floatfield): Evaporation rate field.
    kpbl (Floatfield): Planetary boundary layer height field.
    tke (Floatfield): Turbulent kinetic energy field.
    vvel (Floatfield): Vertical velocity field.
    FRLAND (Floatfield): Fraction of land field.
    NACTL (Floatfield): Activated cloud droplet number concentration field.
    NACTI (Floatfield): Activated ice crystal number concentration field.
    NWFA (Floatfield): Newly formed aerosol number concentration field.
    NN_LAND (Float): Number concentration over land field.
    NN_OCEAN (Float): Number concentration over ocean field.

    Returns:
    None
    """
    with computation(PARALLEL), interval(...):
        NACTL = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
        NACTI = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)

        # might need to change variable names
        tk = t
        press = plo
        air_den = press * 28.8e-3 / 8.31 / tk
        qi = (qicn + qils) * 1.0e3
        ql = (qlcn + qlls) * 1.0e3
        wupdraft = vvel + sqrt(tke)

        # Liquid Clouds
        """
        if tk >= constants.MAPL_TICE - 40.0 and plo > 10000.0 and 0.1 < wupdraft < 100.0:
            ni[1:n_modes] = max(AeroProps[i, j, k, :, 0] * air_den, constants.ZERO_PAR)
            rg[1:n_modes] = max(AeroProps[i, j, k, :, 1] * 0.5 * 1.e6, constants.ZERO_PAR)
            sig0[1:n_modes] = AeroProps[i, j, k, :, 2]
            bibar[1:n_modes] = max(AeroProps[i, j, k, :, 4], constants.ZERO_PAR)
            _get_act_frac_stencil(n_modes, ni[1:n_modes], rg[1:n_modes], sig0[1:n_modes], tk, press, wupdraft, nact[1:n_modes], bibar[1:n_modes])
            numbinit = 0.0
            NACTL = 0.0ma
        """
        # Ice Clouds
        if tk <= constants.MAPL_TICE and (
            qi > np.finfo(float).eps or ql > np.finfo(float).eps
        ):
            numbinit = 0.0
            n = 0
            while n <= constants.n_modes:
                if AERO_DGN[0, 0, 0][n] >= 0.5e-6:
                    numbinit += AERO_NUM[0, 0, 0][n]
                    n += 1
            numbinit *= air_den
            NACTI = (
                constants.AI
                * ((constants.MAPL_TICE - tk) ** constants.BI)
                * (
                    numbinit
                    ** (constants.CI * (constants.MAPL_TICE - tk) + constants.DI)
                )
            )

            # apply limits for NACTL/NACTI
        if NACTL < constants.NN_MIN:
            NACTL = constants.NN_MIN
        if NACTL > constants.NN_MAX:
            NACTL = constants.NN_MAX
        if NACTI < constants.NN_MIN:
            NACTI = constants.NN_MIN
        if NACTI > constants.NN_MAX:
            NACTI = constants.NN_MAX


def ddim_test(this_is_4d: FloatField_NModes):
    with computation(PARALLEL), interval(...):
        lev = 0
        while lev < 10:
            this_is_4d[0, 0, 0][lev] = 42
            lev += 1


class AerActivation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        n_modes: Int,
        USE_AERSOL_NN: bool,
    ):

        if constants.n_modes != n_modes:
            raise NotImplementedError(
                f"Coding limitation: 14 modes are expected, getting {n_modes}"
            )

        if not USE_AERSOL_NN:
            raise NotImplementedError("Non NN Aerosol not implemented")

        """
        self._get_act_frac = stencil_factory.from_dims_halo(
            func = _get_act_frac_stencil,
            compute_dims = [X_DIM, Y_DIM, Z_DIM],
        )
        self._act_frac_mat = stencil_factory.from_dims_halo(
            func = _act_frac_mat_stencil,
            compute_dims = [X_DIM, Y_DIM, Z_DIM],
        )
        self._gser = stencil_factory.from_dims_halo(
            func = _gser_stencil,
            compute_dims = [X_DIM, Y_DIM, Z_DIM],
        )
        self._gcf_matrix = stencil_factory.from_dims_halo(
            func = _gcf_matrix_stencil,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )
        """

        self.ddim_test_stencil = stencil_factory.from_origin_domain(
            func=ddim_test,
            origin=(0, 0, 0),
            domain=stencil_factory.grid_indexing.domain,
        )

    # GEOS_moistGridComp for aero props line  5400ish

    def __call__(
        self,
        q: FloatField,
        t: FloatField,
        plo: FloatField,
        ple: FloatField,
        zlo: FloatField,
        zle: FloatField,
        qlcn: FloatField,
        qicn: FloatField,
        qlls: FloatField,
        qils: FloatField,
        sh: Float,
        evap: Float,
        kpbl: Float,
        tke: FloatField,
        vvel: FloatField,
        FRLAND: Float,
        NACTL: FloatField,
        NACTI: FloatField,
        NWFA: FloatField,
        NN_LAND: Int,
        NN_OCEAN: Int,
    ):
        """
        Perform aerosol activation calculations.

        Parameters:
        IM (Int): Size of the first dimension.
        JM (Int): Size of the second dimension.
        LM (Int): Size of the third dimension.
        q (FloatField): Specific humidity field.
        t (FloatField): Temperature field.
        plo (FloatField): Low-level pressure field.
        ple (FloatField): Low-level pressure field at the end of the time step.
        zlo (FloatField): Low-level height field.
        zle (FloatField): Low-level height field at the end of the time step.
        qlcn (FloatField): Cloud liquid water mixing ratio field.
        qicn (FloatField): Cloud ice mixing ratio field.
        qlls (FloatField): Large-scale cloud liquid water mixing ratio field.
        qils (FloatField): Large-scale cloud ice mixing ratio field.
        sh (Float): Specific humidity value.
        evap (Float): Evaporation rate value.
        kpbl (Float): Planetary boundary layer height value.
        tke (FloatField): Turbulent kinetic energy field.
        vvel (FloatField): Vertical velocity field.
        FRLAND (Float): Fraction of land value.
        NACTL (FloatField): Activated cloud droplet number concentration field.
        NACTI (FloatField): Activated ice crystal number concentration field.
        NWFA (FloatField): Newly formed aerosol number concentration field.
        NN_LAND (Int): Number concentration over land field.
        NN_OCEAN (Int): Number concentration over ocean field.

        Returns:
        None

        self._get_act_frac(
                            nmodes,
                            xnap,
                            rg,
                            sigmag,
                            tkelvin,
                            ptot,
                            wupdraft,
                            nact,
                            bibar
                            )
        self._act_frac_mat(
                            nmodes,
                            xnap,
                            rg,
                            sigmag,
                            bibar,
                            tkelvin,
                            ptot,
                            wupdraft,
                            nact,
                            )
        self._gser(
                    gamser,
                    a,
                    x,
                    gln,
                    )
        self._gcf_matrix(
                        gamser,
                        a,
                        x,
                        gln,
                        )
        """

    def ddim_debug(
        self,
        aero_f_dust: FloatField_NModes,
    ):
        self.ddim_test_stencil(aero_f_dust)
