#-----------Constants are taken from GEOS/src/Shared/@MAPL/shared/Constants--------
#Math constants
MAPL_PI_R8 = 3.14159265358979323846e0 
MAPL_PI = MAPL_PI_R8 
MAPL_DEGREES_TO_RADIANS_R8 = MAPL_PI_R8 / 180.e0
MAPL_DEGREES_TO_RADIANS = MAPL_PI / 180.0
MAPL_RADIANS_TO_DEGREES = 180.e0 / MAPL_PI_R8

#Define whether or not to use CODATA 2018 constants
CODATA_2018_CONSTANTS = True

# Universal Constants
if CODATA_2018_CONSTANTS:
    MAPL_STFBOL = 5.670374419E-8  # W/(m^2 K^4)
    MAPL_AVOGAD = 6.02214076E26   # 1/kmol
    MAPL_RUNIV = 8314.462618      # J/(Kmole K)
else: 
    MAPL_STFBOL = 5.6734E-8  # W/(m^2 K^4)
    MAPL_AVOGAD = 6.023E26   # 1/kmol
    MAPL_RUNIV = 8314.47    # J/(Kmole K)

# Earth Constants
MAPL_PSDRY = 98305.0  # Pa
MAPL_SECONDS_PER_SIDEREAL_DAY = 86164.0  # s
MAPL_GRAV = 9.80665  # m^2/s
MAPL_RADIUS = 6371.0E3  # m
MAPL_OMEGA_R8 = 2.0 * MAPL_PI_R8 / MAPL_SECONDS_PER_SIDEREAL_DAY  # 1/s
MAPL_OMEGA = 2.0 * MAPL_PI / MAPL_SECONDS_PER_SIDEREAL_DAY  # 1/s
MAPL_EARTH_ECCENTRICITY = 8.181919084262200e-2  # --
MAPL_EARTH_SEMIMAJOR_AXIS = 6378137.0  # m
MAPL_KM_PER_DEG = (1.0 / (MAPL_RADIUS / 1000.0)) * MAPL_RADIANS_TO_DEGREES
MAPL_DEG_PER_KM = (MAPL_RADIUS / 1000.0) * MAPL_DEGREES_TO_RADIANS_R8

# Physical properties
MAPL_H2OMW = 18.015  # kg/Kmole
MAPL_O3MW = 47.9982  # kg/Kmole
MAPL_LATENT_HEAT_VAPORIZATION = 2.4665E6  # J/kg @15C @1atm
MAPL_ALHL = MAPL_LATENT_HEAT_VAPORIZATION  # J/kg
MAPL_LATENT_HEAT_FUSION = 3.3370E5  # J/kg @1atm
MAPL_ALHF = MAPL_LATENT_HEAT_FUSION  # J/kg
MAPL_LATENT_HEAT_SUBLIMATION = MAPL_ALHL + MAPL_ALHF  # J/kg
MAPL_ALHS = MAPL_LATENT_HEAT_SUBLIMATION  # J/kg

# Earth Specific Chemistry and Thermodynamic Constants
MAPL_AIRMW = 28.965  # kg/Kmole
MAPL_RDRY = MAPL_RUNIV / MAPL_AIRMW  # J/(kg K)
MAPL_CPDRY = 3.5 * MAPL_RDRY  # J/(kg K)
MAPL_CVDRY = MAPL_CPDRY - MAPL_RDRY  # J/(kg K)
MAPL_RVAP = MAPL_RUNIV / MAPL_H2OMW  # J/(kg K)
MAPL_CPVAP = 4.0 * MAPL_RVAP  # J/(kg K)
MAPL_CVVAP = MAPL_CPVAP - MAPL_RVAP  # J/(kg K)
MAPL_KAPPA = MAPL_RDRY / MAPL_CPDRY  # (2.0/7.0)
MAPL_EPSILON = MAPL_H2OMW / MAPL_AIRMW  # --
MAPL_DELTAP = MAPL_CPVAP / MAPL_CPDRY  # --
MAPL_DELTAV = MAPL_CVVAP / MAPL_CVDRY  # --
MAPL_GAMMAD = MAPL_CPDRY / MAPL_CVDRY  # --
MAPL_RGAS = MAPL_RDRY  # J/(kg K) (DEPRECATED)
MAPL_CP = MAPL_RGAS / MAPL_KAPPA  # J/(kg K) (DEPRECATED)
MAPL_VIREPS = 1.0 / MAPL_EPSILON - 1.0  # (DEPRECATED)
MAPL_P00 = 100000.0  # Pa
MAPL_CAPICE = 2000.0  # J/(K kg)
MAPL_CAPWTR = 4218.0  # J/(K kg)
MAPL_RHOWTR = 1000.0  # kg/m^3
MAPL_NUAIR = 1.533E-5  # m^2/S (@ 18C)
MAPL_TICE = 273.16  # K
MAPL_SRFPRS = 98470  # Pa
MAPL_KARMAN = 0.40  # --
MAPL_USMIN = 1.00  # m/s
MAPL_RHO_SEAWATER = 1026.0  # sea water density [kg/m^3]
MAPL_RHO_SEAICE = 917.0  # sea ice density [kg/m^3]
MAPL_RHO_SNOW = 330.0  # snow density [kg/m^3]
MAPL_CELSIUS_TO_KELVIN = 273.15  # K

#Constants used for radiation coupling 
alhlbcp = MAPL_ALHL / MAPL_CP
alhsbcp = MAPL_ALHS / MAPL_CP

#Constants used for LDRADIUS4_ICE and LDRADIUS4_LIQUID
MAPL_RGAS = 287.0
MAPL_TICE = 273.15
LIQ_RADII_PARAM = 1
ICE_RADII_PARAM = 1

bx = 100.0 * (3.0 / (4.0 * MAPL_PI)) ** (1.0 / 3.0)
r13bbeta = 1.0 / 3.0 - 0.14
abeta = 0.07

RHO_W = 1000.0
Ldiss = 0.07
Lk = 0.75
Lbx = Ldiss * 1.0e3 * (3.0 / (4.0 * MAPL_PI * Lk * RHO_W * 1.0e-3)) ** (1.0 / 3.0)
Lbe = 1.0 / 3.0 - 0.14

#Taken from GEOS_GFDL_1M_InterfaceMod.F90:269-274 - uses MAPL_GetResource to get variables. These do not change in the radiation coupling loop
MIN_RL = 2.5e-6
MAX_RL = 60.0e-6
FAC_RL = 1.0
MIN_RI = 5.0e-6
MAX_RI = 100.0e-6
FAC_RI = 1.0

#Aer Activation Constants 
#32 bit
R_AIR = 3.47e-3 #m3 Pa kg-1K-1
#64 bit
zero_par = 1.e-6 #small non-zero value
ai = 0.0000594
bi = 3.33
ci =  0.0264
di = 0.0033
betaai = -2.262e+3
gamai = 5.113e+6
deltaai = 2.809e+3
densic = 917.0 #Ice crystal density in kgm-3
#32 bit
NN_MIN = 100.0e6
NN_MAX = 1000.0e6

#ACTFRAC_Mat constants - all 64 bit
pi = 3.141592653589793e+00
twopi = 2.0e+00 * pi
sqrt2 = 1.414213562e+00
threesqrt2by2 = 1.5e+00 * sqrt2
avgnum = 6.0221367e+23
rgasjmol = 8.31451e+00 
wmolmass = 18.01528e-03
amolmass = 28.966e-03
asmolmss = 132.1406e-03
denh2o   = 1.00e+03 
denamsul = 1.77e+03 
xnuamsul = 3.00e+00
phiamsul = 1.000e+00
gravity  = 9.81e+00 
heatvap  = 40.66e+03/wmolmass
cpair    = 1006.0e+00
t0dij    = 273.15e+00
p0dij    = 101325.0e+00
dijh2o0  = 0.211e-04