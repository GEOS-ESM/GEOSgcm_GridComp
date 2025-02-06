import gt4py.cartesian.gtscript as gtscript
import pyMoist.constants as constants
from ndsl.dsl.typing import (
    Float,
)
from pyMoist.UW.compute_uwshcu import exnerfn

@gtscript.function
def sign(
    a: Float,
    b: Float,
    ):
    
    #Function that returns the magnitude of one argument and the sign of another.
    
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
    
    # Subroutine to solve the second order polynomial equation. 
    # I should check this subroutine later.   
    
    status = 0

    if a == 0:                           # Form b*x + c = 0
        if b == 0:                       # Failure: c = 0
            status = 1
        else:                              # b*x + c = 0
            r1 = -c/b
        r2 = r1
    else:
        if b == 0:                       # Form a*x**2 + c = 0
            if a*c > 0:                  # Failure: x**2 = -c/a < 0
                status = 2  
            else:                          # x**2 = -c/a 
                r1 = sqrt(-c/a)
            r2 = -r1
        else:                              # Form a*x**2 + b*x + c = 0
            if (b**2 - 4.*a*c) < 0.0:      # Failure, no real roots
                status = 3
            else:
                q  = -0.5*(b + sign(1.0,b)*sqrt(b**2 - 4.*a*c))
                r1 =  q/a
                r2 =  c/q

    return r1, r2, status