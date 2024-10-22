import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, interval, PARALLEL, log
import pyMoist.pyMoist_constants as constants
import pyMoist.radiation_coupling_constants as radconstants
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
)
from ndsl import StencilFactory, QuantityFactory
#from compute_uwshcu import exnerfn

def exnerfn(p):
    """
    Function that calculates the Exner function for a given pressure.

    Inputs:
    p (Float): Atmospheric pressure [Pa]

    Returns:
    (p / 100000.0) ** (radconstants.MAPL_RGAS / constants.cp) (Float): Exner function
    """

    return (p / 100000.0) ** (radconstants.MAPL_RGAS / constants.cp)

def getbuoy(pbot,thv0bot,ptop,thv0top,thvubot,thvutop):
    '''
    Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and 
    'cinlcl, plfc' if any. Assume 'thv' is linear in each layer 
    both for cumulus and environment. Note that this subroutine 
    only includes positive CIN in calculation - if there is any 
    negative CIN, it is assumed to be zero.    This is slightly 
    different from 'single_cin' below, where both positive  and 
    negative CIN are included. 
    '''
    if thvubot > thv0bot and thvutop > thv0top:
        plfc = pbot
        return plfc
    elif thvubot < thv0bot and thvutop < thv0top:
        cin  = cin - ( (thvubot/thv0bot - 1.) + (thvutop/thv0top - 1.)) * (pbot - ptop) / ( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top*exnerfn(ptop)) )
        return cin
    elif thvubot > thv0bot and thvutop < thv0top:
        frc  = ( thvutop/thv0top - 1.) / ( (thvutop/thv0top - 1.) - (thvubot/thv0bot - 1.) )
        cin  = cin - ( thvutop/thv0top - 1.) * ( (ptop + frc*(pbot - ptop)) - ptop ) / ( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top*exnerfn(ptop)) )
        return cin
    else:            
        frc  = ( thvubot/thv0bot - 1.) / ( (thvubot/thv0bot - 1.) - (thvutop/thv0top - 1.) )
        plfc = pbot - frc * ( pbot - ptop )
        cin  = cin - ( thvubot/thv0bot - 1.)*(pbot - plfc)/( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top * exnerfn(ptop)))
        return plfc, cin
    
pbot = 93720.7188
thv0bot = 293.494476
ptop = 92191.5234
thv0top = 293.764496
thvubot = 293.707458
thvutop = 294.277283

plfc = getbuoy(pbot, thv0bot,ptop,thv0top,thvubot,thvutop)

print(plfc)

'''
@gtscript.function
def getbuoy(
    pbot: Float,
    thv0bot: Float,
    ptop: Float,
    thv0top: Float,
    thvubot: Float,
    thvutop: Float,
    plfc: Float,
    cin: Float,
    ):
    
    #Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and 
    #'cinlcl, plfc' if any. Assume 'thv' is linear in each layer 
    #both for cumulus and environment. Note that this subroutine 
    #only includes positive CIN in calculation - if there is any 
    #negative CIN, it is assumed to be zero.    This is slightly 
    #different from 'single_cin' below, where both positive  and 
    #negative CIN are included. 
    
    if thvubot > thv0bot and thvutop > thv0top:
        plfc = pbot
        #return
    elif thvubot < thv0bot and thvutop < thv0top:
        cin  = cin - ( (thvubot/thv0bot - 1.) + (thvutop/thv0top - 1.)) * (pbot - ptop) / ( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top*exnerfn(ptop)) )
    elif thvubot > thv0bot and thvutop < thv0top:
        frc  = ( thvutop/thv0top - 1.) / ( (thvutop/thv0top - 1.) - (thvubot/thv0bot - 1.) )
        cin  = cin - ( thvutop/thv0top - 1.) * ( (ptop + frc*(pbot - ptop)) - ptop ) / ( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top*exnerfn(ptop)) )
    else:            
        frc  = ( thvubot/thv0bot - 1.) / ( (thvubot/thv0bot - 1.) - (thvutop/thv0top - 1.) )
        plfc = pbot - frc * ( pbot - ptop )
        cin  = cin - ( thvubot/thv0bot - 1.)*(pbot - plfc)/( pbot/(radconstants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(radconstants.MAPL_RGAS*thv0top * exnerfn(ptop)))

    return plfc, cin # These returns are optional

def GetBuoy(
    plcl: FloatField_IJ,
    thv0lcl: FloatField_IJ,
    pifc0: FloatField_IJ,
    thv0top: FloatField_IJ,
    thvubot: FloatField_IJ,
    thvutop: FloatField_IJ,
    plfc: FloatField_IJ,
    cin: FloatField_IJ,
    ):

    with computation(PARALLEL), interval(...):
        plfc, cin = getbuoy(plcl,thv0lcl,pifc0,thv0top,thvubot,thvutop)


class Get_Buoy:
    def __init__(
        self,
        stencil_factory: StencilFactory,
        quantity_factory: QuantityFactory,
    ) -> None:

        self.stencil_factory = stencil_factory
        self.quantity_factory = quantity_factory
        grid_indexing = stencil_factory.grid_indexing
        self._Get_Buoy = self.stencil_factory.from_dims_halo(
            func=GetBuoy,
            compute_dims=[X_DIM, Y_DIM, Z_DIM],
        )

    def __call__(
        self,
        plcl: FloatField_IJ,
        thv0bot: FloatField_IJ,
        pifc0: FloatField_IJ,
        thv0top: FloatField_IJ,
        thvubot: FloatField_IJ,
        thvutop: FloatField_IJ,
        plfc: FloatField_IJ,
        cin: FloatField_IJ,
    ):

        self._Get_Buoy(
            plcl=plcl,
            thv0bot=thv0bot,
            pifc0=pifc0,
            thv0top=thv0top,
            thvubot= thvubot,
            thvutop=thvutop,
            plfc=plfc,
            cin=cin,
        )
'''

  
    

