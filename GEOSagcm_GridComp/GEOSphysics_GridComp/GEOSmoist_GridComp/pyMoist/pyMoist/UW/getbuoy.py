import gt4py.cartesian.gtscript as gtscript
import pyMoist.constants as constants
from ndsl.dsl.typing import (
    Float,
)
from pyMoist.UW.compute_uwshcu import exnerfn


@gtscript.function
def getbuoy(
    pbot: Float,
    thv0bot: Float,
    ptop: Float,
    thv0top: Float,
    thvubot: Float,
    thvutop: Float,
    cin_in: Float,
    plfc_in: Float,
    ):
    
    #Subroutine to calculate integrated CIN [ J/kg = m2/s2 ] and 
    #'cinlcl, plfc' if any. Assume 'thv' is linear in each layer 
    #both for cumulus and environment. Note that this subroutine 
    #only includes positive CIN in calculation - if there is any 
    #negative CIN, it is assumed to be zero.    This is slightly 
    #different from 'single_cin' below, where both positive  and 
    #negative CIN are included. 
    plfc = plfc_in
    cin = cin_in

    if thvubot > thv0bot and thvutop > thv0top:
        plfc = pbot
    elif thvubot < thv0bot and thvutop < thv0top:
        cin  = cin_in - ( (thvubot/thv0bot - 1.) + (thvutop/thv0top - 1.)) * (pbot - ptop) / ( pbot/(constants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(constants.MAPL_RGAS*thv0top*exnerfn(ptop)) )
    elif thvubot > thv0bot and thvutop < thv0top:
        frc  = ( thvutop/thv0top - 1.) / ( (thvutop/thv0top - 1.) - (thvubot/thv0bot - 1.) )
        cin  = cin_in - ( thvutop/thv0top - 1.) * ( (ptop + frc*(pbot - ptop)) - ptop ) / ( pbot/(constants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(constants.MAPL_RGAS*thv0top*exnerfn(ptop)) )
    else:          
        frc  = ( thvubot/thv0bot - 1.) / ( (thvubot/thv0bot - 1.) - (thvutop/thv0top - 1.) )
        plfc = pbot - frc * ( pbot - ptop )
        cin  = cin_in - ( thvubot/thv0bot - 1.)*(pbot - plfc)/( pbot/(constants.MAPL_RGAS*thv0bot*exnerfn(pbot)) + ptop/(constants.MAPL_RGAS*thv0top * exnerfn(ptop)))

    return plfc, cin # Note: plfc and cin are returned, but not always used

  
    

