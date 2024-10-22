import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, exp, FORWARD
import pyMoist.pyMoist_constants as constants
import pyMoist.radiation_coupling_constants as radconstants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    FloatFieldIJ,
    Float,
)
from ndsl import StencilFactory, QuantityFactory


@gtscript.function
def compute_alpha(
    del_CIN: Float,
    ke: Float,
    ):
    
    #Subroutine to compute proportionality factor for 
    #implicit CIN calculation. 
    
    x0 = 0.0
    del_CIN8_f64 = del_CIN
    ke8_f64 = ke
    iteration = 0
    while iteration < 10:
        x1 = x0 - (exp(-x0*ke8_f64*del_CIN8_f64) - x0)/(-ke8_f64*del_CIN8_f64*exp(-x0*ke8_f64*del_CIN8_f64) - 1.0)
        x0 = x1
        iteration+=1
     
    compute_alpha = x0

    return compute_alpha


def ComputeAlpha(
    del_CIN: FloatFieldIJ,
    ke: FloatFieldIJ,
    alpha: FloatFieldIJ,
    ):

    with computation(FORWARD), interval(...):
        alpha = compute_alpha(del_CIN, ke)

        
class Compute_Alpha:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._Compute_Alpha = self.stencil_factory.from_dims_halo(
            func=ComputeAlpha,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        del_CIN: FloatFieldIJ,
        ke: FloatFieldIJ,
        alpha: FloatFieldIJ,
    ):

        self._Compute_Alpha(
            del_CIN=del_CIN,
            ke=ke,
            alpha=alpha,
        )


  
    

