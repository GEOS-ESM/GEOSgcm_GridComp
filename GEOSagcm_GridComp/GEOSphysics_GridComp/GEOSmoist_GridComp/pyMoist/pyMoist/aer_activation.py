from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float, FloatFieldIJ
from ndsl import Quantity, QuantityFactory, StencilFactory
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log, exp, sqrt
import pyMoist.aer_activation_constants as constants
import numpy as np

# Global space
FloatField_NModes = gtscript.Field[gtscript.IJK, (Float, (constants.n_modes))]

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
    stp = 2.5066282746310005
    
    x = xx
    y = x
    tmp = x + 5.5
    tmp = (x + 0.5) * log(tmp) - tmp
    ser = 1.000000000190015
    
    ser += 76.18009172947146 / (y+1)
    ser += -86.50532032941677 / (y+1)
    ser += 24.01409824083091 / (y+1)
    ser += -1.231739572450155 / (y+1)
    ser += 0.001208650973866179 / (y+1)
    ser += -0.000005395239384953 / (y+1)

    gammln = tmp + log(stp * ser / x)
    return gammln

@gtscript.function
def _gser_stencil(
    a: Float, 
    x: Float, 
    gln: Float,
)-> Float:
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
    eps = 3.0e-9
    itmax = 10000
    gln = GammLn(a)
    if x <= 0:
        #if x < 0:
            #raise ValueError("x < 0 in gser")
        gamser = 0.0
    else:
        ap = a
        sum_ = 1.0 / a
        del_ = sum_
        n = 0
        while n < itmax:
            ap += 1.0
            del_ *= x / ap
            sum_ += del_
            if abs(del_) < abs(sum_) * eps: #this might be wrong, might need ()
                gamser = sum_ * exp(-x + a * log(x) - gln)
                n = itmax
            n += 1
        gamser = sum_ * exp(-x + a * log(x) - gln)
    return gamser

@gtscript.function
def _gcf_matrix_stencil( 
    a: Float, 
    x: Float, 
    gln: Float,
)->Float:
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
    itmax = 10000
    eps = 3.0e-7
    fpmin = 1.0e-30
    gln = GammLn(a)
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b
    h = d

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
            i = itmax+1 #breaks the loop
        i += 1
    return exp(-x + a * log(x) - gln) * h

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
    #Fortran messages here potential bad arguments
    #TODO; Allow print in GT4Py
    #if (x < 0.0) or (a <= 0.0):
    #    raise ValueError("aero_actv: function gammp: bad arguments")
    gln = GammLn(a)
    if x < a + 1.0:
        gammp = _gser_stencil(a, x, gln)
    else:
        gammp = 1.0 - _gcf_matrix_stencil(a, x, gln)
    return gammp

@gtscript.function
def Erf(
    x: Float
)->Float:
    """
    Compute the error function for a given value x.
    
    Parameters:
    x (Float): Input value for which the error function is to be computed.
    
    Returns:
    Float: The error function value for the input x.
    """
    erf = 0.0
    if x < 0.0e+00:
        erf = -1.0*GammP(0.5, x**2)
    else:
        erf = GammP(0.5, x**2)
    return erf

def aer_activation_stencil(
        aero_dgn: FloatField_NModes,
        aero_num: FloatField_NModes,
        nacti: FloatField,
        t: FloatField,
        plo: FloatField,
        qicn: FloatField,
        qils: FloatField,
        qlcn: FloatField,
        qlls: FloatField,
        nn_land: Float,
        frland: FloatFieldIJ,
        nn_ocean: Float,
        aero_hygroscopicity: FloatField_NModes,
        nwfa: FloatField,
        nactl:FloatField,
        vvel: FloatField,
        tke: FloatField,
        aero_sigma: FloatField_NModes,
        nact: FloatField_NModes,
        ni: FloatField_NModes,
        rg: FloatField_NModes,
        sig0: FloatField_NModes,
        bibar: FloatField_NModes
): 
    with computation(PARALLEL), interval(...):
        
        nfaux = 0.0
        n = 0
        while n < constants.n_modes:
            if aero_hygroscopicity[0,0,0][n] > 0.4: #aero_kap is aero_hygroscopicity
                nfaux += aero_num[0,0,0][n]
            n += 1
        nwfa =  nfaux

        nactl = nn_land * frland + nn_ocean * (1.0 - frland)
        nacti = nn_land * frland + nn_ocean * (1.0 - frland)

        #determing aerosol number concentration at cloud base
        tk = t
        press = plo
        air_den = press * 28.8e-3 / 8.31 / tk
        qi = (qicn + qils) * 1.0e3
        ql = (qlcn + qlls) * 1.0e3
        wupdraft = vvel + sqrt(tke)

        # Liquid Clouds
        if (tk >= (constants.MAPL_TICE - 40.0)) and (plo > 10000.0) and (0.1 < wupdraft) and (wupdraft < 100.0):
            n = 0
            while n < constants.n_modes:
                ni_temp = aero_num[0,0,0][n]* air_den
                # rg_temp = aero_dgn[0,0,0][n] * 0.5 * 1.e6
                # sig0_temp = aero_sigma[0,0,0][n]
                # bibar_temp = aero_hygroscopicity[0,0,0][n]
                ni = max(ni_temp, constants.ZERO_PAR)
                # rg = max(rg_temp, constants.ZERO_PAR)
                # sig0 = sig0_temp
                # bibar = max(bibar_temp, constants.ZERO_PAR)
                n += 1 

            #rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
            #values close to those given in a-z et al. 1998 figure 5. 
            
            rdrp = 0.105e-06 #[m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.
            #These variables are common to all modes and need only be computed once. 
            dv = constants.DIJH2O0 * (constants.P0DIJ / plo) * (tk / constants.T0DIJ) ** 1.94e+00 #[m^2/s] (p&k,2nd ed., p.503)
            surten = 76.10e-3 - 0.155e-3 * (tk - 273.15e+00) #[j/m^2]
            wpe = exp(77.34491296 - 7235.424651 / tk - 8.2 * log(tk) + tk * 5.7113e-3) #[pa]
            dumw = sqrt(constants.TWOPI * constants.WMOLMASS / constants.RGASJMOL / tk) #[s/m]
            dvprime = dv / ((rdrp / (rdrp + constants.DELTAV)) + (dv * dumw / (rdrp * constants.ALPHAC)))
            xka = (5.69 + 0.017 * (tk - 273.15)) * 418.4e-5
            duma = sqrt(constants.TWOPI * constants.AMOLMASS / constants.RGASJMOL / tk)
            xkaprime = xka / ((rdrp / (rdrp + constants.DELTAT)) + (xka * duma / (rdrp * constants.ALPHAT * constants.DENH2O * constants.CPAIR)))
            g = 1.0 / ((constants.DENH2O * constants.RGASJMOL * tk) / (wpe * dvprime * constants.WMOLMASS) + 
                        ((constants.HEATVAP * constants.DENH2O) / (xkaprime * tk)) * 
                        ((constants.HEATVAP * constants.WMOLMASS) / (constants.RGASJMOL * tk) - 1.0))
            a = (2.0 * surten * constants.WMOLMASS) / (constants.DENH2O * constants.RGASJMOL * tk)
            alpha = (constants.GRAVITY / (constants.RGASJMOL * tk)) * ((constants.WMOLMASS * constants.HEATVAP) / (constants.CPAIR * tk) - constants.AMOLMASS)
            gamma = (constants.RGASJMOL * tk) / (wpe * constants.WMOLMASS) + (constants.WMOLMASS * constants.HEATVAP * constants.HEATVAP) / (constants.CPAIR * plo * constants.AMOLMASS * tk)
            dum = sqrt(alpha * wupdraft / g)
            zeta = 2.0 * a * dum / 3.0

            #These variables must be computed for each mode
            n = 0
            while n < constants.n_modes:
                xlogsigm = log(sig0[0,0,0][n])
                smax = 0.0 #double precision
                sm = (2.0 / sqrt(bibar[0,0,0][n])) * (a / (3.0 * rg[0,0,0][n])) ** 1.5
                eta = dum ** 3 / (constants.TWOPI * constants.DENH2O * gamma * ni[0,0,0][n])
                f1 = 0.5 * exp(2.50 * xlogsigm ** 2)
                f2 = 1.0 + 0.25 * xlogsigm
                smax = smax + (f1*(zeta/eta)**1.5+ f2*(sm**2/(eta+3.0*zeta))**0.75)/sm**2
                n +=1
            
            smax = 1.0e+00 / sqrt(smax)  
            n = 0
            u = 0.0
            while n < constants.n_modes:
                sm = (2.0 / sqrt(bibar[0,0,0][n])) * (a / (3.0 * rg[0,0,0][n])) ** 1.5
                xlogsigm = log(sig0[0,0,0][n])
                ac = rg[0,0,0][n] * (sm / smax) ** 0.66666666666666667
                u = log(ac / rg[0,0,0][n]) / (constants.SQRT2 * xlogsigm)
                fracactn = 0.5 * (1.0 - Erf(u))
                nact[0,0,0][n] = fracactn * ni[0,0,0][n]
                n+=1

            numbinit = 0.0
            nactl = 0.0
            n = 0
            while n < constants.n_modes:
                numbinit += aero_num[0,0,0][n] * air_den
                nactl += nact[0,0,0][n]
                n += 1
            nactl = min(nactl, 0.99*numbinit)
            
        # Ice Clouds
        if (tk <= constants.MAPL_TICE) and (qi > constants.FLOAT_TINY or ql > constants.FLOAT_TINY):
            numbinit = 0.0
            n = 0
            while n < constants.n_modes:
                if aero_dgn[0,0,0][n] >= 0.5e-6:
                    numbinit += aero_num[0,0,0][n]
                n += 1
            numbinit *= air_den
            nacti = (
                constants.AI
                * ((constants.MAPL_TICE - tk) ** constants.BI)
                * (
                    numbinit
                    ** (constants.CI * (constants.MAPL_TICE - tk) + constants.DI)
                )
            )

        # apply limits for NACTL/NACTI
        if nactl < constants.NN_MIN:
            nactl = constants.NN_MIN
        if nactl > constants.NN_MAX:
            nactl = constants.NN_MAX
        if nacti < constants.NN_MIN:
            nacti = constants.NN_MIN
        if nacti > constants.NN_MAX:
            nacti = constants.NN_MAX

class AerActivation:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
        n_modes: Int,
        USE_AERSOL_NN: bool,
    ):
        #Dead lines in aer_actv_single_moment.F90:238-244 for literal kpbli 
        #Variables never used: numbinit, WC, BB, RAUX

        self._nact = quantity_factory._numpy.zeros(
            (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        )
        self._ni = quantity_factory._numpy.zeros(
            (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        )
        self._rg = quantity_factory._numpy.zeros(
            (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        )
        self._sig0 = quantity_factory._numpy.zeros(
            (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        )
        self._bibar = quantity_factory._numpy.zeros(
            (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        )
        #self._ni_temp = quantity_factory._numpy.zeros(
        #    (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        #)
        #self._rg_temp = quantity_factory._numpy.zeros(
        #    (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        #)
        #self._sig0_temp = quantity_factory._numpy.zeros(
        #    (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        #)
        #self._aero_hygroscopicity_temp = quantity_factory._numpy.zeros(
        #    (stencil_factory.grid_indexing.domain[0], stencil_factory.grid_indexing.domain[1], stencil_factory.grid_indexing.domain[2], constants.n_modes), dtype=Float
        #)


        
        if constants.n_modes != n_modes:
            raise NotImplementedError(
                f"Coding limitation: 14 modes are expected, getting {n_modes}"
            )

        if not USE_AERSOL_NN:
            raise NotImplementedError("Non NN Aerosol not implemented")
        
        self.higher_dimensional_storages = stencil_factory.from_origin_domain(
            func=aer_activation_stencil,
            origin=(0, 0, 0),
            domain=stencil_factory.grid_indexing.domain,
        )

    def __call__(
        self,
        aero_dgn: FloatField_NModes,
        aero_num: FloatField_NModes,
        nacti: FloatField,
        t: FloatField,
        plo: FloatField,
        qicn: FloatField,
        qils: FloatField,
        qlcn: FloatField,
        qlls: FloatField,
        nn_land: Float,
        frland: FloatFieldIJ,
        nn_ocean: Float,
        aero_hygroscopicity: FloatField_NModes,
        nwfa: FloatField,
        nactl:FloatField,
        vvel: FloatField,
        tke: FloatField,
        aero_sigma: FloatField_NModes,
    ):
        self.higher_dimensional_storages(
            aero_dgn,
            aero_num,
            nacti,
            t,
            plo,
            qicn,
            qils,
            qlcn,
            qlls,
            nn_land,
            frland,
            nn_ocean,
            aero_hygroscopicity,
            nwfa,
            nactl,
            vvel,
            tke,
            aero_sigma,
            self._nact,
            self._ni,
            self._rg,
            self._sig0,
            self._bibar,
            
        )        

        # self._ni_temp,
        # self._rg_temp,
        # self._sig0_temp,
        # self._aero_hygroscopicity_temp