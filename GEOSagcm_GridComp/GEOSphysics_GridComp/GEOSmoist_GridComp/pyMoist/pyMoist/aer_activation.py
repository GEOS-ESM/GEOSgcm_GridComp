from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float 
from ndsl import Quantity, QuantityFactory, StencilFactory
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log, exp, sqrt, BACKWARD
import pyMoist.aer_activation_constants as constants

import sys
sys.path.append('serialbox/python')
import serialbox as ser

@gtscript.function
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
        erf =- GammP(0.5, x**2)
    else:
        erf = GammP(0.5, x**2)
    return erf

@gtscript.function
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
    
@gtscript.function
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
    if x < a + 1.0:
        _gser_stencil(gamser, a, x, gln)
        gammp = gamser
    else:
        _gcf_matrix_stencil(gammcf, a, x, gln)
        gammp = 1.0e00 - gammcf
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
        xlogsigm = log(sigmag[n])
        smax = 0.0
        while n < constants.n_modes:
            sm[n] = (2.0 / sqrt(bibar[n])) * (a / (3.0 * rg[n])) ** 1.5
            eta[n] = dum ** 3 / (constants.TWOPI * constants.DENH2O * gamma * xnap[n])
            f1 = 0.5 * exp(2.50 * xlogsigm[n] ** 2)
            f2 = 1.0 + 0.25 * xlogsigm[n]
            smax = smax + (f1*(zeta/eta[n])**1.5+ f2*(sm[n]**2/(eta[n]+3.0*zeta))**0.75)/sm[n]**2
            n +=1
        
        smax = 1.0e+00 / sqrt(smax)  
        n=0
        while n < constants.n_modes:
            #lines 534
            ac[n] = rg[n] * (sm[n] / smax) ** 0.66666666666666667
            u = log(ac[n] / rg[n]) / (constants.SQRT2 * xlogsigm[n])
            fracactn[n] = 0.5 * (1.0 - Erf(u))
            nact[n] = fracactn[n] * xnap[n]
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
        '''
        rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
        values close to those given in a-z et al. 1998 figure 5. 
        '''
        rdrp = 0.105e-06 #[m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.
        sm = [constants.n_modes]
        eta = [constants.n_modes]
        ac = [constants.n_modes]
        fracactn = [constants.n_modes]
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
            gamser = 0
        else:
            ap = a
            sum_ = 1.0 / a
            del_ = sum_
            for n in range(1, itmax + 1):
                ap += 1
                del_ *= x / ap
                sum_ += del_
                if abs(del_) < abs(sum_) * eps:
                    raise ValueError("aero_actv: function _gser_stencil: a too large, itmax too small")
            gamser = sum_ * exp(-x + a * log(x) - gln)

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
            break
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
        
@gtscript.function
def aer_activation_lquid_clouds(ni,rg,sig0,tk,press,wupdraft,bibar,nact,air_den):
    _get_act_frac_stencil(ni, rg, sig0, tk, press, wupdraft, bibar, nact)
    numbinit = 0.0
    NACTL = 0.0
    for n in range(constants.n_modes):
        numbinit += AeroProps[i, j, k, n, 0] * air_den
        NACTL += nact[n]
    NACTL = min(NACTL, 0.99 * numbinit)

@gtscript.function
def aer_activation_ice_clouds():
    
#do not inlcude aero_aci and AeroProps in the inputs into the stencil
def aer_activation_stencil( 
        IM: Int, 
        JM: Int, 
        LM: Int,
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
        USE_AERO_BUFFER: bool, #get rid of. set to True
        AeroProps: FloatField, #get rid of
        aero_aci: Int, #get rid of
        NACTL: FloatField, 
        NACTI: FloatField,
        NWFA: FloatField, 
        NN_LAND: Float, 
        NN_OCEAN: Float,
):
    """
    Perform aerosol activation computations based on input atmospheric and aerosol properties.
    
    Parameters:
    IM (int): Size of the first dimension.
    JM (int): Size of the second dimension.
    LM (int): Size of the third dimension.
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
    USE_AERO_BUFFER (bool): Flag to use aerosol buffer.
    AeroProps (Floatfield): Aerosol properties field.
    aero_aci (Int): Aerosol-cloud interaction field.
    NACTL (Floatfield): Activated cloud droplet number concentration field.
    NACTI (Floatfield): Activated ice crystal number concentration field.
    NWFA (Floatfield): Newly formed aerosol number concentration field.
    NN_LAND (Float): Number concentration over land field.
    NN_OCEAN (Float): Number concentration over ocean field.
    
    Returns:
    None
    """

    with computation(PARALLEL), interval(...):

        #n_modes for loop between 1-8 
        n = 1
        aci_num
        aci_dgn
        aci_sigma
        aci_hygroscopicity
        aci_density
        aci_f_dust
        aci_f_soot
        aci_f_organic

        

        while n <= 8:
            if n == 1:
                num = aci_num
            elif n == 2:
                dpg = aci_dgn
            elif n == 3:
                sig = aci_sigma
            elif n == 4:
                kap = aci_hygroscopicity
            elif n == 5:
                den = aci_density
            elif n == 6:
                fdust = aci_f_dust 
            elif n == 7:
                fsoot = aci_f_soot
            elif n == 8:
                forg = aci_f_organic
            n += 1
        
        USE_AERO_BUFFER = True
        USE_AEROSOL_NN = True
        AeroProps_num = 0.0

        #loop for all buffer dimensions
        AeroPropsBuffer_Loop(aci_num, aci_dgn, aci_sigma, aci_hygroscopicity, aci_density, aci_f_dust, aci_f_soot, aci_f_organic)

        num[n] = aci_num
        dpg[n] = aci_dgn
        sig[n] = aci_sigma
        kap[n] = aci_hygroscopicity
        den[n] = aci_density
        fdust[n] = aci_f_dust
        fsoot[n] = aci_f_soot
        forg[n] =  aci_f_organic
        nmods = n_modes 

        kpbli = max(min(round(kpbl), LM-1), 1).astype(Int) #line 96 of fortran

        #Activated aerosol # concentration for liq/ice phases (units: m^-3)
        numbinit = 0.0
        WC = 0.0
        BB = 0.0
        RAUX = 0.0

        #determining aerosol number concentration at cloud base
        k = kpbli
        tk = t
        press = plo
        air_den = press * 28.8e-3 / 8.31 / tk

        with computation(BACKWARD), interval(...):
            NACTL = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
            NACTI = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)

            with computation(PARALLEL), interval(...): 
                #might need to change variable names
                tk = t
                press = plo
                air_den = press * 28.8e-3 / 8.31 / tk
                qi = (qicn + qils) * 1.e+3 
                ql = (qlcn + qlls) * 1.e+3
                wupdraft = vvel + sqrt(tke)
                
                #Liquid Clouds
                if tk >= constants.MAPL_TICE - 40.0 and plo > 10000.0 and 0.1 < wupdraft < 100.0:
                    ni[1:n_modes] = max(AeroProps[i, j, k, :, 0] * air_den, constants.ZERO_PAR)
                    rg[1:n_modes] = max(AeroProps[i, j, k, :, 1] * 0.5 * 1.e6, constants.ZERO_PAR)
                    sig0[1:n_modes] = AeroProps[i, j, k, :, 2]
                    bibar[1:n_modes] = max(AeroProps[i, j, k, :, 4], constants.ZERO_PAR)
                    _get_act_frac_stencil(n_modes, ni[1:n_modes], rg[1:n_modes], sig0[1:n_modes], tk, press, wupdraft, nact[1:n_modes], bibar[1:n_modes])
                    numbinit = 0.0
                    NACTL = 0.0
                    aer_activation_lquid_clouds(ni,rg,sig0,tk,press,wupdraft,bibar,nact,air_den)

                #Ice Clouds
                if tk <= constants.MAPL_TICE and (qi > finfo(float).eps or ql > finfo(float).eps):
                    aer_activation_ice_clouds()
                    numbinit = 0.0
                    for n in range(constants.n_modes):
                        if AeroProps[i, j, k, n, 1] >= 0.5e-6:
                            numbinit += AeroProps[i, j, k, n, 0]
                    numbinit *= air_den
                    NACTI = constants.AI * ((constants.MAPL_TICE - tk) ** constants.BI) * (numbinit ** (constants.CI * (constants.MAPL_TICE - tk) + constants.DI))

                     #apply limits for NACTL/NACTI
            if NACTL < NN_MIN:
                NACTL = NN_MIN
            if NACTL > NN_MAX:
                NACTL = NN_MAX
            if NACTI < NN_MIN:
                NACTI = NN_MIN
            if NACTI > NN_MAX:
                NACTI = NN_MAX
            else:
                NACTL = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
                NACTI = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)

class AerActivation:
    def __init__(
            self,
            stencil_factory: StencilFactory,
            quantity_factory: QuantityFactory,
            do_qa: bool,
    ):
        
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
        
#need to make literals for all the esmf and mapl calls
#GEOS_moistGridComp for aero props line  5400ish

#I got rid of aero_aci call
    def __call__(
        self,
        IM: Int,
        JM: Int,
        LM: Int,
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
        USE_AERO_BUFFER: bool,
        AeroProps: FloatField, 
        NACTL: FloatField,
        NACTI: FloatField,
        NWFA: FloatField,
        NN_LAND: Int,
        NN_OCEAN: Int
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
        USE_AERO_BUFFER (bool): Flag to use aerosol buffer.
        AeroProps (FloatField): Aerosol properties field.
        NACTL (FloatField): Activated cloud droplet number concentration field.
        NACTI (FloatField): Activated ice crystal number concentration field.
        NWFA (FloatField): Newly formed aerosol number concentration field.
        NN_LAND (Int): Number concentration over land field.
        NN_OCEAN (Int): Number concentration over ocean field.

        Returns:
        None
        """
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

    #start port at line 213