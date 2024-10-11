import numpy as np

from ndsl.dsl.typing import Float


"""
All global constants used for aer_actv_single_moment.F90
"""
# MAPL Constants
MAPL_TICE = 273.16  # K

# 64 bit
ZERO_PAR = 1.0e-6  # small non-zero value
AI = 0.0000594
BI = 3.33
CI = 0.0264
DI = 0.0033

BETAAI = -2.262e3
GAMAI = 5.113e6
DELTAAI = 2.809e3
DENSIC = 917.0  # Ice crystal density in kgm-3

# 32 bit
NN_MIN = 100.0e6
NN_MAX = 1000.0e6

# ACTFRAC_Mat constants - all 64 bit
PI = 3.141592653589793e00
TWOPI = 2.0e00 * PI
SQRT2 = 1.414213562e00
THREESQRT2BY2 = 1.5e00 * SQRT2

AVGNUM = 6.0221367e23  # [1/mol]
RGASJMOL = 8.31451e00  # [j/mol/k]
WMOLMASS = 18.01528e-03  # molar mass of h2o [kg/mol]
AMOLMASS = 28.966e-03  # molar mass of air     [kg/mol]
ASMOLMSS = 132.1406e-03  # molar mass of nh42so4 [kg/mol]
DENH2O = 1.00e03  # density of water [kg/m^3]
DENAMSUL = 1.77e03  # density of pure ammonium sulfate [kg/m^3]
XNUAMSUL = 3.00e00  # number of ions formed when the salt is dissolved in water [1]
PHIAMSUL = 1.000e00  # osmotic coefficient value in a-r 1998. [1]
GRAVITY = 9.81e00  # grav. accel. at the earth's surface [m/s/s]
HEATVAP = 40.66e03 / WMOLMASS  # latent heat of vap. for water and tnbp [j/kg]
CPAIR = 1006.0e00  # heat capacity of air [j/kg/k]
T0DIJ = 273.15e00  # reference temp. for dv [k]
P0DIJ = 101325.0e00  # reference pressure for dv [pa]
DIJH2O0 = 0.211e-04  # reference value of dv [m^2/s] (p&k,2nd ed., p.503)
DELTAV = 1.096e-07  # vapor jump length [m]
DELTAT = 2.160e-07  # thermal jump length [m]
ALPHAC = 1.000e00  # condensation mass accommodation coefficient [1]
ALPHAT = 0.960e00  # thermal accommodation coefficient [1]

# Define how many modes in an aerosol
n_modes = 14

# Python euqivalent of Fortran's tiny(X)
FLOAT_TINY = np.finfo(Float).tiny
