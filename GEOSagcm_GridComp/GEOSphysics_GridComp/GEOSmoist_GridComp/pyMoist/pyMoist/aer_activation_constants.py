'''Aer Activation Constants are taken from GEOS/src/Shared/@MAPL/shared/Constants'''

from ndsl.dsl.typing import Float
import numpy as np

MAPL_TICE = 273.16  # K

#32 bit
R_AIR = 3.47e-3 #m3 Pa kg-1K-1

#64 bit
ZERO_PAR = 1.e-6 #small non-zero value
AI = 0.0000594
BI = 3.33
CI =  0.0264
DI = 0.0033
BETAAI = -2.262e+3
GAMAI = 5.113e+6
DELTAAI = 2.809e+3
DENSIC = 917.0 #Ice crystal density in kgm-3
#32 bit
NN_MIN = 100.0e6
NN_MAX = 1000.0e6

#ACTFRAC_Mat constants - all 64 bit
#parameters
PI = 3.141592653589793e+00
TWOPI = 2.0e+00 * PI
SQRT2 = 1.414213562e+00
THREESQRT2BY2 = 1.5e+00 * SQRT2
AVGNUM = 6.0221367e+23
RGASJMOL = 8.31451e+00 
WMOLMASS = 18.01528e-03
AMOLMASS = 28.966e-03
ASMOLMSS = 132.1406e-03
DENH2O   = 1.00e+03 
DENAMSUL = 1.77e+03 
XNUAMSUL = 3.00e+00
PHIAMSUL = 1.000e+00
GRAVITY  = 9.81e+00 
HEATVAP  = 40.66e+03/WMOLMASS
CPAIR    = 1006.0e+00
T0DIJ    = 273.15e+00
P0DIJ    = 101325.0e+00
DIJH2O0  = 0.211e-04
DELTAV = 1.096e-07
DELTAT = 2.160e-07
ALPHAC = 1.000e+00
ALPHAT = 0.960e+00
 
n_modes = 14

FLOAT_EPSILON = np.finfo(Float).eps 