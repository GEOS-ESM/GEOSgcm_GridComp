import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, sqrt, exp
import pyMoist.pyMoist_constants as constants
import pyMoist.radiation_coupling_constants as radconstants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
    Int,
)
from ndsl import StencilFactory, QuantityFactory
import math
import numpy as np

def compute_mumin2(mulcl,rmaxfrac,mulow):
    '''
    Subroutine to compute critical 'mu' (normalized CIN) such 
    that updraft fraction at the LCL is equal to 'rmaxfrac'.  
    '''
    x0 = mulow
    iteration = 0
    while iteration < 10:
       ex = np.exp(-x0**2)
       ef = math.erfc(x0) # Fortran error fraction function
       exf = ex/ef
       f  = 0.5*exf**2 - 0.5*(ex/2.0/rmaxfrac)**2 - (mulcl*2.5066/2.0)**2
       fs = (2.*exf**2)*(exf/np.sqrt(radconstants.MAPL_PI)-x0) + (0.5*x0*ex**2)/(rmaxfrac**2)
       x1 = x0 - f/fs     
       x0 = x1
       iteration+=1

    compute_mumin2 = x0

    return compute_mumin2

mulcl = 0.958077192
rmaxfrac = 0.100000001
mu = 0.983174622

mumin2 = compute_mumin2(mulcl,rmaxfrac,mu)
print(mumin2)

'''
@gtscript.function
def compute_mumin2(
    mulcl: Float,
    rmaxfrac: Float,
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
       f  = 0.5*exf**2 - 0.5*(ex/2._r8/rmaxfrac)**2 - (mulcl*2.5066/2._r8)**2
       fs = (2.*exf**2)*(exf/sqrt(radconstants.MAPL_PI)-x0) + (0.5*x0*ex**2)/(rmaxfrac**2)
       x1 = x0 - f/fs     
       x0 = x1
       iteration+=1

    compute_mumin2 = x0

    return compute_mumin2


def ComputeMumin2(
    mulcl: Float,
    rmaxfrac: Float,
    mulow: Float,
    mumin2: Float,
    ):

    with computation(PARALLEL), interval(...):
        mumin2 = compute_mumin2(mulcl, rmaxfrac, mulow)

class Compute_Mumin2:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._Compute_Mumin2 = self.stencil_factory.from_dims_halo(
            func=ComputeMumin2,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        mulcl: Float,
        rmaxfrac: Float,
        mulow: Float,
        mumin2: Float,
    ):

        self._Get_Roots(
            mulcl = mulcl,
            rmaxfrac = rmaxfrac,
            mulow = mulow,
            mumin2 = mumin2,
        )
'''
    

