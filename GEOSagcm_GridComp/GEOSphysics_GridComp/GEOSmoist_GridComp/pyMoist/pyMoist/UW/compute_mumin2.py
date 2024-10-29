import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, sqrt, exp
import pyMoist.constants as constants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    Int,
)
from ndsl import StencilFactory, QuantityFactory

@gtscript.function
def erfc(x: Float):
    '''
    Complimentary error function approximation borowed from 
    Abramowits and Stegun's Handbook of Mathematical Functions
    '''

    # Constants for the approximation
    p = 0.3275911
    a1, a2, a3, a4, a5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429

    # Approximation for the error function complement
    t = 1 / (1 + p * x)
    erf_approx = 1 - (a1 * t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5) * exp(-x**2)

    return erf_approx

@gtscript.function
def compute_mumin2(
    mulcl: Float,
    rmaxfrax: Float,
    mulow: Float,
    ):
    
    #Subroutine to compute critical 'mu' (normalized CIN) such 
    #that updraft fraction at the LCL is equal to 'rmaxfrac'.  
    
    x0 = mulow
    iteration = 0
    while iteration < 10:
       ex = exp(-x0**2)
       ef = erfc(x0) # Fortran error fraction function
       exf = ex/ef
       f  = 0.5*exf**2 - 0.5*(ex/2.0/rmaxfrax)**2 - (mulcl*2.5066/2.0)**2
       fs = (2.*exf**2)*(exf/sqrt(constants.MAPL_PI)-x0) + (0.5*x0*ex**2)/(rmaxfrax**2)
       x1 = x0 - f/fs     
       x0 = x1
       iteration+=1

    compute_mumin2 = x0

    return compute_mumin2
    

