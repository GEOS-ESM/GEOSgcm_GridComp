import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, sqrt, exp, erfc
import pyMoist.constants as constants
from ndsl.dsl.typing import (
    Float,
    Int,
)


@gtscript.function
def compute_mumin2(
    mulcl: Float,
    rmaxfrax: Float,
    mulow: Float,
    ):
    
    #Subroutine to compute critical 'mu' (normalized CIN) such 
    #that updraft fraction at the LCL is equal to 'rmaxfrac'.  
    
    x0:f64 = mulow
    iteration = 0
    while iteration < 10:
       ex:f64 = exp(-x0**2)
       ef:f64 = erfc(x0) # Complimentary error fraction function
       exf:f64 = ex/ef
       f:f64  = f64(0.5)*exf**2 - f64(0.5)*(ex/f64(2.0)/rmaxfrax)**2 - (mulcl*f64(2.5066)/f64(2.0))**2
       fs:f64 = (f64(2.0)*exf**2)*(exf/sqrt(constants.MAPL_PI)-x0) + (f64(0.5)*x0*ex**2)/(rmaxfrax**2)
       x1:f64 = x0 - f/fs     
       x0 = x1
       iteration+=1

    compute_mumin2 = f32(x0)

    return compute_mumin2
    

