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
from pyMoist.UW.compute_uwshcu import exnerfn


@gtscript.function
def single_cin(
   pbot,
   thv0bot,
   ptop,
   thv0top,
   thvubot,
   thvutop
    ):
    
    #Function to calculate a single layer CIN by summing all 
    #positive and negative CIN.  
    
    single_cin = ( (1.- thvubot/thv0bot) + (1.- thvutop/thv0top)) * ( pbot - ptop ) / ( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top*exnerfn(ptop)) )

    return single_cin

def SingleCin(
    pbot: FloatFieldIJ,
    thv0bot: FloatFieldIJ,
    ptop: FloatFieldIJ,
    thv0top: FloatFieldIJ,
    thvubot: FloatFieldIJ,
    thvutop: FloatFieldIJ,
    cin1: FloatFieldIJ,
    cin2: FloatFieldIJ,
    ):

    with computation(FORWARD), interval(...):
        cin2 = cin1 + single_cin(pbot, thv0bot,ptop,thv0top,thvubot,thvutop)

class Single_Cin:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._SingleCin = self.stencil_factory.from_dims_halo(
            func=SingleCin,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        pbot: FloatFieldIJ,
        thv0bot: FloatFieldIJ,
        ptop: FloatFieldIJ,
        thv0top: FloatFieldIJ,
        thvubot: FloatFieldIJ,
        thvutop: FloatFieldIJ,
        cin1: FloatFieldIJ,
        cin2: FloatFieldIJ,
    ):

        self._Single_Cin(
            pbot=pbot,
            thv0bot=thv0bot,
            ptop=ptop,
            thv0top=thv0top,
            thvubot=thvubot,
            thvutop=thvutop,
            cin1=cin1,
            cin2=cin2,
        )


  
    


