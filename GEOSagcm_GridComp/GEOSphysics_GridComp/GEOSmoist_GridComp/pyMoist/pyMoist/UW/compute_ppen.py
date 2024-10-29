from ndsl import Namelist, Quantity, StencilFactory
from ndsl.dsl.typing import FloatField, Float
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.stencils.testing.translate import TranslateFortranData2Py
from pyMoist.saturation import QSat
import gt4py.cartesian.gtscript as gtscript
from gt4py.cartesian.gtscript import computation, PARALLEL, interval, exp, FORWARD
import gt4py.cartesian.gtscript as gtscript
from ndsl.constants import X_DIM, Y_DIM, Z_DIM
from ndsl.dsl.typing import (
    Float,
)

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
            iteration+=1
   
    compute_ppen = -max(0.0,min(dpen,x0))

    return compute_ppen