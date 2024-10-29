import gt4py.cartesian.gtscript as gtscript
import pyMoist.constants as constants
from ndsl.dsl.typing import (
    Float,
)
from pyMoist.UW.compute_uwshcu import exnerfn

@gtscript.function
def exnerfn(
    p: Float,
    ):
    """
    Function that calculates the Exner function for a given pressure.

    Inputs:
    p (Float): Atmospheric pressure [Pa]

    Returns:
    (p / 100000.0) ** (const.MAPL_RGAS / const.MAPL_CP) (Float): Exner function
    """

    return (p / 100000.0) ** (constants.MAPL_RGAS / constants.MAPL_CP)

@gtscript.function
def single_cin(
   pbot: Float,
   thv0bot: Float,
   ptop: Float,
   thv0top: Float,
   thvubot: Float,
   thvutop: Float,
    ):
    
    #Function to calculate a single layer CIN by summing all 
    #positive and negative CIN.  
    
    single_cin = ( (1.- thvubot/thv0bot) + (1.- thvutop/thv0top)) * ( pbot - ptop ) / ( pbot/(constants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(constants.MAPL_RGAS*thv0top*exnerfn(ptop)) )

    return single_cin