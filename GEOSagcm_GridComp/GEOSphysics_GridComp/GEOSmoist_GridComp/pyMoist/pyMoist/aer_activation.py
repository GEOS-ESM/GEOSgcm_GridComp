from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import FloatField, Int, Float 
from ndsl import Quantity, QuantityFactory, StencilFactory
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log, exp, sqrt  # type: ignore
import pyMoist.aer_activation_constants as constants

def erf(
        x:Float
)-> None:
    
def get_act_frac_stencil(
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
        # Call to act_frac_mat_stencil
        _act_frac_mat_stencil(nmodes, xnap, rg, sigmag, bibar, tkelvin, ptot, wupdraft, nact)

def _act_frac_mat_stencil(
    nmodes: Int,
    xnap: 
    rg:
    sigmag:
    bibar:
    tkelvin: Float,
    ptot:
    wupdraft: FloatField, 
    nact: FloatField
):
    with computation(PARALLEL), interval(...):
        rdrp = 0.105e-06 
        for n in range(ni.shape[0]):
            dv = constants.DIJH2O0 * (constants.P0DIJ / ptot) * (tkelvin / constants.T0DIJ) ** 1.94
            surten = 76.10e-3 - 0.155e-3 * (tkelvin - 273.15)
            wpe = exp(77.34491296 - 7235.424651 / tkelvin - 8.2 * log(tkelvin) + tkelvin * 5.7113e-3)
            dumw = sqrt(constants.TWOPI * constants.WMOLMASS / constants.RGASJMOL / tkelvin)
            dvprime = dv / ((constants.DELTAV / (constants.DELTAV + constants.DELTAV)) + (dv * dumw / (constants.DELTAV * constants.ALPHAC)))
            xka = (5.69 + 0.017 * (tkelvin - 273.15)) * 418.4e-5
            duma = sqrt(constants.TWOPI * constants.AMOLMASS / constants.RGASJMOL / tkelvin)
            xkaprime = xka / ((constants.DELTAV / (constants.DELTAV + constants.DELTAT)) + (xka * duma / (constants.DELTAV * constants.ALPHAT * constants.DENH2O * constants.CPAIR)))
            g = 1.0 / ((constants.DENH2O * constants.RGASJMOL * tkelvin) / (wpe * dvprime * constants.WMOLMASS) + 
                       ((constants.HEATVAP * constants.DENH2O) / (xkaprime * tkelvin)) * 
                       ((constants.HEATVAP * constants.WMOLMASS) / (constants.RGASJMOL * tkelvin) - 1.0))
            a = (2.0 * surten * constants.WMOLMASS) / (constants.DENH2O * constants.RGASJMOL * tkelvin)
            alpha = (constants.GRAVITY / (constants.RGASJMOL * tkelvin)) * ((constants.WMOLMASS * constants.HEATVAP) / (constants.CPAIR * tkelvin) - constants.AMOLMASS)
            gamma = (constants.RGASJMOL * tkelvin) / (wpe * constants.WMOLMASS) + (constants.WMOLMASS * constants.HEATVAP * constants.HEATVAP) / (constants.CPAIR * ptot * constants.AMOLMASS * tkelvin)
            dum = sqrt(alpha * wupdraft / g)
            zeta = 2.0 * a * dum / 3.0
            #These variables must be computed for each mode
            xlogsigm = log(sig0[n])
            sm[n] = (2.0 / sqrt(bibar[n])) * (a / (3.0 * rg[n])) ** 1.5
            eta = dum ** 3 / (constants.TWOPI * constants.DENH2O * gamma * ni[n])
            f1 = 0.5 * exp(2.50 * xlogsigm ** 2)
            f2 = 1.0 + 0.25 * xlogsigm
            smax = 1.0 / sqrt(f1 * (zeta / eta) ** 1.5 + f2 * (sm[n] ** 2 / (eta + 3.0 * zeta)) ** 0.75)
            ac = rg[n] * (sm[n] / smax) ** 0.66666666666666667
            u = log(ac / rg[n]) / (constants.SQRT2 * xlogsigm)
            fracactn = 0.5 * (1.0 - erf(u))
            nact[n] = fracactn * ni[n]

def _gser_stencil(
    a: Float, 
    x: Float, 
    gamser: FloatField, 
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
                    break
            gamser[...] = sum_ * exp(-x + a * log(x) - gln)

def _gcf_matrix_stencil(
    a: Float, 
    x: Float, 
    gammcf: FloatField, 
    gln: FloatField
):
    itmax = 10000
    eps = 3.0e-7
    fpmin = 1.0e-30
    with computation(PARALLEL), interval(...):
        b = x + 1.0 - a
        c = 1.0 / fpmin
        d = 1.0 / b
        h = d
        for i in range(1, itmax + 1):
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
        gammcf[...] = exp(-x + a * log(x) - gln) * h

def aer_activation_stencil(
        IM: int, 
        JM: int, 
        LM: int,
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
        USE_AERO_BUFFER,
        AeroProps,
        areo_aci,
        NACTL,
        NACTI,
        NWFA,
        NN_LAND: Float, 
        NN_OCEAN: Float):

    with computation(PARALLEL), interval(...):
        AeroProps_num = 0.0

class AerActivation:
    def __init__(
            self,
            stencil_factory: StencilFactory,
            quantity_factory: QuantityFactory,
            do_qa: bool,
    ):
        '''
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
        '''
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
        kpbli = max(min(round(kpbl), LM-1), 1).astype(Int) #line 96 of fortran

        for j in range(JM):
            for i in range(IM):
                k = kpbli[i, j]
                tk = t[i, j, k]
                press = plo[i, j, k]
                air_den = press * 28.8e-3 / 8.31 / tk

        for k in range(LM-1, -1, -1):
            NACTL[:, :, k] = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
            NACTI[:, :, k] = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
            for j in range(JM):
                for i in range(IM):
                    tk = t[i, j, k]
                    press = plo[i, j, k]
                    air_den = press * 28.8e-3 / 8.31 / tk
                    qi = (qicn[i, j, k] + qils[i, j, k]) * 1.e+3
                    ql = (qlcn[i, j, k] + qlls[i, j, k]) * 1.e+3
                    wupdraft = vvel[i, j, k] + sqrt(tke[i, j, k])

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
                        NACTI[i, j, k] = ai * ((constants.MAPL_TICE - tk) ** bi) * (numbinit ** (ci * (constants.MAPL_TICE - tk) + di))

                    NACTL[i, j, k] = clip(NACTL[i, j, k], constants.NN_MIN, constants.NN_MAX)
                    NACTI[i, j, k] = clip(NACTI[i, j, k], constants.NN_MIN, constants.NN_MAX)

        else:
            NACTL[:, :, :] = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)
            NACTI[:, :, :] = NN_LAND * FRLAND + NN_OCEAN * (1.0 - FRLAND)

    '''
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