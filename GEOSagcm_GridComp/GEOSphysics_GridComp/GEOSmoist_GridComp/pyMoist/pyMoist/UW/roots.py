import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, sqrt
import pyMoist.pyMoist_constants as constants
import pyMoist.radiation_coupling_constants as radconstants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    Int,
)
from ndsl import StencilFactory, QuantityFactory
import numpy as np

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
    
    #Subroutine to solve the second order polynomial equation. 
    #I should check this subroutine later.   
    
    status = 0

    if a == 0.0:                           #Form b*x + c = 0
        if b == 0.0:                       # Failure: c = 0
            status = 1
        else:                              # b*x + c = 0
            r1 = -c/b
        r2 = r1
    else:
        if b == 0.0:                       # Form a*x**2 + c = 0
            if a*c > 0.0:                  # Failure: x**2 = -c/a < 0
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



def Roots(
    aquad: Float,
    bquad: Float,
    cquad: Float,
    xc1: Float,
    xc2: Float,
    status: Int,
    ):

    with computation(PARALLEL), interval(...):
        xc1, xc2, status = roots(aquad,bquad,cquad)


class Get_Roots:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._Get_Roots = self.stencil_factory.from_dims_halo(
            func=Roots,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        aquad: Float,
        bquad: Float,
        cquad: Float,
        xc1: Float,
        xc2: Float,
        status: Int,
    ):

        self._Get_Roots(
            aquad = aquad,
            bquad = bquad,
            cquad = cquad,
            xc1 = xc1,
            xc2 = xc2,
            status = status,
        )



