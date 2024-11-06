import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import exp
import pyMoist.constants as constants
from ndsl.dsl.typing import (
    Float,
)


@gtscript.function
def compute_alpha(
    del_CIN: Float,
    ke: Float,
    ):
    
    #Subroutine to compute proportionality factor for 
    #implicit CIN calculation. 
    
    x0:f64 = f64(0.0)
    del_CIN8_f64:f64 = f64(del_CIN)
    ke8_f64:f64 = ke
    iteration = 0
    while iteration < 10:
        x1 = x0 - (exp(-x0*ke8_f64*del_CIN8_f64) - x0)/(-ke8_f64*del_CIN8_f64*exp(-x0*ke8_f64*del_CIN8_f64) - 1.0)
        x0 = x1
        iteration+=1
     
    compute_alpha = f32(x0)

    return compute_alpha