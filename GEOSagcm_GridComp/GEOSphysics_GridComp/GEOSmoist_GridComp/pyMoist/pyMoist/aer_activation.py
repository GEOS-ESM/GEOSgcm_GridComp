from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float 
from ndsl import Quantity, QuantityFactory, StencilFactory
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log, exp, sqrt
import pyMoist.aer_activation_constants as constants

@gtscript.function
def Erf(
    x:Float
):
    erf = 0.0
    if x < 0.0e+00:
        erf =- GammP(0.5, x**2)
    else:
        erf = GammP(0.5, x**2)
    return erf

@gtscript.function
def GammLn(xx):
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
    #Not sure what these variables are
    gamser = 0.0
    gammcf = 0.0
    gln = 0.0
    if x < 0.0 or a <= 0.0:
        raise ValueError("aero_actv: function gammp: bad arguments")
    if x < a+1.0:
        _gser_stencil(gamser,a,x,gln)
        gammp=gamser
    else:
        _gcf_matrix_stencil(gammcf,a,x,gln)
        gammp=1.0e00-gammcf
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
    bibar: FloatField
):
    with computation(PARALLEL), interval(...):
        # Call to _act_frac_mat_stencil
        _act_frac_mat_stencil(nmodes, xnap, rg, sigmag, bibar, tkelvin, ptot, wupdraft, nact)

def _act_frac_mat_stencil(
    nmodes: Int,
    xnap: Int,
    rg: Float,
    sigmag: Float,
    bibar: Float,
    tkelvin: Float,
    ptot: Float,
    wupdraft: FloatField, 
    nact: FloatField,
):
    with computation(PARALLEL), interval(...):
        '''
        rdrp is the radius value used in eqs.(17) & (18) and was adjusted to yield eta and zeta 
        values close to those given in a-z et al. 1998 figure 5. 
        '''
        rdrp = 0.105e-06 #[m] tuned to approximate the results in figures 1-5 in a-z et al. 1998.
        n_modes = 8 #This might change
        sm = [n_modes]
        eta = [n_modes]
        ac = [n_modes]
        fracactn = [n_modes]
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
        n = 0
        xlogsigm = log(sigmag[n])
        smax = 0.0
        while n < n_modes:
            sm[n] = (2.0 / sqrt(bibar[n])) * (a / (3.0 * rg[n])) ** 1.5
            eta[n] = dum ** 3 / (constants.TWOPI * constants.DENH2O * gamma * xnap[n])
            f1 = 0.5 * exp(2.50 * xlogsigm[n] ** 2)
            f2 = 1.0 + 0.25 * xlogsigm[n]
            smax = smax + (f1*(zeta/eta[n])**1.5+ f2*(sm[n]**2/(eta[n]+3.0*zeta))**0.75)/sm[n]**2
            n +=1
        
        smax = 1.0e+00 / sqrt(smax)  
        n=0
        while n < n_modes:
            #lines 534
            ac[n] = rg[n] * (sm[n] / smax) ** 0.66666666666666667
            u = log(ac[n] / rg[n]) / (constants.SQRT2 * xlogsigm[n])
            fracactn[n] = 0.5 * (1.0 - Erf(u))
            nact[n] = fracactn[n] * xnap[n]
            n+=1

def _gser_stencil(
    gamser: FloatField, 
    a: Float, 
    x: Float, 
    gln: FloatField
):
    eps = 3.0e-9
    itmax = 10000
    with computation(PARALLEL), interval(...):
        if x <= 0:
            if x < 0:
                raise ValueError("x < 0 in gser")
            gamser[...] = 0
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
            gamser[...] = sum_ * exp(-x + a * log(x) - gln)


def _gcf_matrix_stencil(
    gammcf: FloatField, 
    a: Float, 
    x: Float, 
    gln: FloatField
):
    with computation(PARALLEL), interval(...):
        itmax = 10000
        eps = 3.0e-7
        fpmin = 1.0e-30
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
                break
            i += 1
        else:
            raise ValueError("AERO_ACTV: Function GCF: A too large, Itmax too small.")
        gammcf[...] = exp(-x + a * log(x) - gln) * h
        return gammcf


#def aer_activation_lquid_clouds():

#def aer_activation_ice_clourds():
    
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
        NN_OCEAN: Float
):

    with computation(PARALLEL), interval(...):
        USE_AERO_BUFFER = True
        USE_AEROSOL_NN = True
        AeroProps_num = 0.0
        n_modes = 8
        kpbli = max(min(round(kpbl), LM-1), 1).astype(Int) #line 96 of fortran

        k = kpbli
        tk = t
        press = plo
        air_den = press * 28.8e-3 / 8.31 / tk

        for k in range(LM-1, -1, -1):
            NACTL = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
            NACTI = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
            for j in range(JM):
                for i in range(IM):
                    tk = t
                    press = plo
                    air_den = press * 28.8e-3 / 8.31 / tk
                    qi = (qicn+ qils) * 1.e+3
                    ql = (qlcn + qlls) * 1.e+3
                    wupdraft = vvel + sqrt(tke)

                    #Liquid Clouds

                    if tk >= constants.MAPL_TICE - 40.0 and plo[i, j, k] > 10000.0 and 0.1 < wupdraft < 100.0:
                        ni[:] = max(AeroProps[i, j, k, :, 0] * air_den, constants.ZERO_PAR)
                        rg[:] = max(AeroProps[i, j, k, :, 1] * 0.5 * 1.e6, constants.ZERO_PAR)
                        sig0[:] = AeroProps[i, j, k, :, 2]
                        bibar[:] = max(AeroProps[i, j, k, :, 4], constants.ZERO_PAR)
                        self._get_act_frac(ni, rg, sig0, tk, press, wupdraft, bibar, nact)
                        numbinit = 0.0
                        NACTL[i, j, k] = 0.0
                        for n in range(n_modes):
                            numbinit += AeroProps[i, j, k, n, 0] * air_den
                            NACTL[i, j, k] += nact[n]
                        NACTL[i, j, k] = min(NACTL[i, j, k], 0.99 * numbinit)
                    
                    #Ice Clouds
                    if tk <= constants.MAPL_TICE and (qi > np.finfo(float).eps or ql > np.finfo(float).eps):
                        numbinit = 0.0
                        for n in range(n_modes):
                            if AeroProps[i, j, k, n, 1] >= 0.5e-6:
                                numbinit += AeroProps[i, j, k, n, 0]
                        numbinit *= air_den
                        NACTI = ai * ((constants.MAPL_TICE - tk) ** bi) * (numbinit ** (ci * (constants.MAPL_TICE - tk) + di))

                    NACTL = clip(NACTL[i, j, k], constants.NN_MIN, constants.NN_MAX)
                    NACTI = clip(NACTI[i, j, k], constants.NN_MIN, constants.NN_MAX)

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
    def __call__(
        self,
        IM,
        JM,
        LM,
        q,
        t,
        plo,
        ple,
        zlo,
        zle,
        qlcn,
        qicn,
        qlls,
        qils,
        sh,
        evap,
        kpbl,
        tke,
        vvel,
        FRLAND,
        USE_AERO_BUFFER: bool,
        AeroProps: FloatField, 
        aero_aci,
        NACTL,
        NACTI,
        NWFA,
        NN_LAND,
        NN_OCEAN

    ):
        aci_num
        aci_dgn
        aci_sigma
        aci_hygroscopicity
        aci_density
        aci_f_dust
        aci_f_soot
        aci_f_organic

        #n_modes for loop between 1-8 
        n_modes = 1
        while n_modes < 8:
            if n_modes == 1:
                num = aci_num
            elif n_modes == 2:
                dpg = aci_dgn
            elif n_modes == 3:
                sig = aci_sigma
            elif n_modes == 4:
                kap = aci_hygroscopicity
            elif n_modes == 5:
                den = aci_density
            elif n_modes == 6:
                fdust = aci_f_dust 
            elif n_modes == 7:
                fsoot = aci_f_soot
            elif n_modes == 8:
                forg = aci_f_organic
            n_modes += 1

''''
    self._get_act_frac_stencil(

    )
    self._act_frac_mat_stencil(

    )
    self._gser_stencil(

    )
    self._gcf_matrix_stencil(

    )
'''

 #start port at line 213