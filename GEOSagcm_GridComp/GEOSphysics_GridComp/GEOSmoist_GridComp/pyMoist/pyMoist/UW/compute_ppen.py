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
import numpy as np

def compute_ppen(
    wtwb,D,bogbot,bogtop,rho0j,dpen
    ):
    '''
    Subroutine to compute critical 'ppen[Pa]<0' ( pressure dis. 
    from 'pifc0(kpen-1)' to the cumulus top where cumulus updraft 
    vertical velocity is exactly zero ) by considering exact    
    non-zero fer(kpen).  
    '''
    # Buoyancy slope
    SB = ( bogtop - bogbot ) / dpen

    # Sign of slope, 'f' at x = 0
    # If 's00>0', 'w' increases with height.
    s00 = bogbot / rho0j - D * wtwb

    if D*dpen < 1.0e-4:
        if s00 >= 0.0:
            x0 = dpen       
        else:
            x0 = max(0.0,min(dpen,-0.5*wtwb/s00))
    else:
        if s00 >= 0.0:
            x0 = dpen
        else:
            x0 = 0.0

        iteration=0
        while iteration < 5:
            aux  = min(max(-2.0*D*x0, -20.0), 20.0)
           
            f  = np.exp(aux)*(wtwb-(bogbot-SB/(2.0*D))/(D*rho0j)) + (SB*x0+bogbot-SB/(2.0*D))/(D*rho0j)
            fs = -2.0*D*np.exp(aux)*(wtwb-(bogbot-SB/(2.0*D))/(D*rho0j)) + (SB)/(D*rho0j)
          
            x1 = x0 - f/fs     
            x0 = x1
            iteration+=1
   
    compute_ppen = -max(0.0,min(dpen,x0))

    return compute_ppen

wtwb = 2.48292089E-02
drage = 9.80464276E-04
bogbot = -1.46985054E-04
bogtop = -2.57849693E-04
rhomid0j = 1.05870783
dp0 = 1522.86719

ppen = compute_ppen(wtwb,drage,bogbot,bogtop,rhomid0j,dp0)
print(ppen)

'''
@gtscript.function
def compute_ppen(
    wtwb: Float,
    D: Float,
    bogbot: Float,
    bogtop: Float,
    rho0j: Float,
    dpen: Float,
    ):
    
    #Subroutine to compute critical 'ppen[Pa]<0' ( pressure dis. 
    #from 'pifc0(kpen-1)' to the cumulus top where cumulus updraft 
    #vertical velocity is exactly zero ) by considering exact    
    #non-zero fer(kpen).  
    
    # Buoyancy slope
    SB = ( bogtop - bogbot ) / dpen

    # Sign of slope, 'f' at x = 0
    # If 's00>0', 'w' increases with height.
    s00 = bogbot / rho0j - D * wtwb

    if D*dpen < 1.0e-4:
        if s00 >= 0.0:
            x0 = dpen       
        else:
            x0 = max(0.0,min(dpen,-0.5*wtwb/s00))
    else:
        if s00 >= 0.0:
            x0 = dpen
        else:
            x0 = 0.0

        iteration=0
        while iteration < 5:
            aux  = min(max(-2.0*D*x0, -20.0), 20.0)
           
            f  = exp(aux)*(wtwb-(bogbot-SB/(2.0*D))/(D*rho0j)) + (SB*x0+bogbot-SB/(2.0*D))/(D*rho0j)
            fs = -2.0*D*exp(aux)*(wtwb-(bogbot-SB/(2.0*D))/(D*rho0j)) + (SB)/(D*rho0j)
          
            x1 = x0 - f/fs     
            x0 = x1
   
    compute_ppen = -max(0.0,min(dpen,x0))

    return compute_ppen
'''