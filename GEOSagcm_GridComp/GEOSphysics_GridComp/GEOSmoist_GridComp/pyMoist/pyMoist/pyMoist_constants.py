### File containing all constants used for pyMoist ###

from ndsl.dsl.typing import Float


# Generic Global Constants
taufrz = 450.0


# Following taken from MathConstants.F90
PI = Float(3.14159265358979323846)
deg_2_rad = PI / 180
rad_to_deg = 180 / PI

# Following taken from PhysicalConstants.F90
# Universal Constants
CODATA_2018_CONSTANTS = True  # set for now, needs to be dynamic
if CODATA_2018_CONSTANTS is True:
    stfbol = Float(5.670374419e-8)  # W/(m^2 K^4)
    avogad = Float(6.02214076e26)  # 1/kmol
    runiv = Float(8314.462618)  # J/(Kmole K)
else:
    stfbol = Float(5.6734e-8)  # W/(m^2 K^4)
    avogad = Float(6.023e26)  # 1/kmol
    runiv = Float(8314.47)  # J/(Kmole K)

# Earth Constants
psdry = Float(98305)  # Pa
seconds_per_sidereal_day = Float(86164.0)  # s
gravity = Float(9.80665)  # m^2/s
radius = Float(6371.0e3)  # m
omega = 2 * PI / seconds_per_sidereal_day  # 1/s
eath_eccentricity = Float(8.181919084262200e-2)  # --
earth_semimajor_axis = Float(6378137)  # m
km_per_deg = (1 / (radius / 1000)) * rad_to_deg
deg_per_km = (radius / 1000) * deg_2_rad

# Physical Properties
h2omw = Float(18.015)  # kg/Kmole
o3mw = Float(47.9982)  # kg/Kmole
latent_heat_vaporization = Float(2.4665e6)  # J/kg @15C @1atm
latent_heat_fusion = Float(3.3370e5)  # J/kg @1atm
latent_heat_sublimation = latent_heat_vaporization + latent_heat_fusion  # J/kg

# Earth Specific Chemistry and Thermodynamic Constants
airmw = Float(28.965)  # kg/Kmole
rdry = runiv / airmw  # J/(kg K)
cpdry = 3.5 * rdry  # J/(kg K)
cvdry = cpdry - rdry  # J/(kg K)
rvap = runiv / h2omw  # J/(kg K)
cpvap = 4 * rvap  # J/(kg K)
cvvap = cpvap - rvap  # J/(kg K)
kappa = rdry / cpdry  # (2.0/7.0) # translator note no clue what 2/7 references
cp = rdry / kappa
epsilon = h2omw / airmw  # --
deltap = cpvap / cpdry  # --
deltav = cvvap / cvdry  # --
gammad = cpdry / cvdry  # --
reference_pressure = Float(100000)  # Pa
capice = Float(2000)  # J/(K kg)
capwtr = Float(4218)  # J/(K kg)
rhowtr = Float(1000)  # kg/m^3
nuair = Float(1.533e-5)  # m^2/s (@ 18C)
t_ice = Float(273.16)  # K
surface_pressure = Float(98470)  # Pa
karman = Float(0.40)  # --
usmin = Float(1.00)  # m/s
rho_seawater = Float(1026.0)  # sea water density [kg/m^3]
rho_seaice = Float(917.0)  # sea ice density [kg/m^3]
rho_snow = Float(330.0)  # snow density [kg/m^3]
celsius_to_kelvin = Float(273.15)  # K


# ice_fraction constants
# In anvil/convective clouds
aT_ICE_ALL = Float(252.16)
aT_ICE_MAX = Float(268.16)
aICEFRPWR = Float(2.0)
# Over snow/ice SRF_TYPE = 2
iT_ICE_ALL = Float(236.16)
iT_ICE_MAX = Float(261.16)
iICEFRPWR = Float(6.0)
# Over Land     SRF_TYPE = 1
lT_ICE_ALL = Float(239.16)
lT_ICE_MAX = Float(261.16)
lICEFRPWR = Float(2.0)
# Over Oceans   SRF_TYPE = 0
oT_ICE_ALL = Float(238.16)
oT_ICE_MAX = Float(263.16)
oICEFRPWR = Float(4.0)
# Jason method constants (translator note: I have no clue who jason is)
# In anvil/convective clouds
JaT_ICE_ALL = Float(245.16)
JaT_ICE_MAX = Float(261.16)
JaICEFRPWR = Float(2.0)
# Over snow/ice
JiT_ICE_ALL = Float(t_ice - 40.0)
JiT_ICE_MAX = Float(t_ice)
JiICEFRPWR = Float(4.0)

# Other miscellaneous parameters
k_cond = Float(2.4e-2)  # 1/(J*s*K)
diffu = Float(2.2e-5)  # m^2/s
taufrz = Float(450)
dqcmax = Float(1e-4)
K_COND = Float(2.4e-2)  # J m**-1 s**-1 K**-1
DIFFU = Float(2.2e-5)  # m**2 s**-1
dQCmax = Float(1.0e-4)

# cloud radius Constants based on DOI 10.1088/1748-9326/3/4/045021
LIQ_RADII_PARAM = 2
ICE_RADII_PARAM = 1
RHO_W = Float(1000)
RHO_I = Float(916.8)  # Density of ice crystal in kg/m^3
abeta = Float(0.07)
r13bbeta = Float(1.0 / 3.0 - 0.14)
bx = Float(100.0 * (3.0 / (4.0 * PI)) ** (1.0 / 3.0))
Ldiss = Float(0.07)  # tunable dispersion effect
Lk = Float(0.75)  # tunable shape effect (0.5:1)
Lbe = Float(1.0 / 3.0 - 0.14)
Lbx = Float(Ldiss * 1.0e3 * (3.0 / (4.0 * PI * Lk * RHO_W * 1.0e-3)) ** (1.0 / 3.0))
# cloud radius eqs are in cgs units

# GEOS_Utilities Constatns
TIMX = Float(-20.0)

ESFAC = Float(h2omw / airmw)
MAX_MIXING_RATIO = Float(1.0)
ZEROC = Float(t_ice)

TMINTBL = Float(150.0)
TMAXTBL = Float(333.0)
DEGSUBS = Float(100)
ERFAC = Float(DEGSUBS / ESFAC)
DELTA_T = Float(1.0 / DEGSUBS)
TABLESIZE = Float((TMAXTBL - TMINTBL) * DEGSUBS + 1)
TMIX = Float(-20.0)

UTBL = True
TYPE = Float(1)

FIRST = True

TMINSTR = Float(-95.0)
TSTARR1 = Float(-75.0)
TSTARR2 = Float(-65.0)
TSTARR3 = Float(-50.0)
TSTARR4 = Float(-40.0)
TMAXSTR = Float(60.0)

B6 = Float(6.136820929e-11 * 100.0)
B5 = Float(2.034080948e-8 * 100.0)
B4 = Float(3.031240396e-6 * 100.0)
B3 = Float(2.650648471e-4 * 100.0)
B2 = Float(1.428945805e-2 * 100.0)
B1 = Float(4.436518521e-1 * 100.0)
B0 = Float(6.107799961e0 * 100.0)
BI6 = Float(1.838826904e-10 * 100.0)
BI5 = Float(4.838803174e-8 * 100.0)
BI4 = Float(5.824720280e-6 * 100.0)
BI3 = Float(4.176223716e-4 * 100.0)
BI2 = Float(1.886013408e-2 * 100.0)
BI1 = Float(5.034698970e-1 * 100.0)
BI0 = Float(6.109177956e0 * 100.0)
S16 = Float(0.516000335e-11 * 100.0)
S15 = Float(0.276961083e-8 * 100.0)
S14 = Float(0.623439266e-6 * 100.0)
S13 = Float(0.754129933e-4 * 100.0)
S12 = Float(0.517609116e-2 * 100.0)
S11 = Float(0.191372282e0 * 100.0)
S10 = Float(0.298152339e1 * 100.0)
S26 = Float(0.314296723e-10 * 100.0)
S25 = Float(0.132243858e-7 * 100.0)
S24 = Float(0.236279781e-5 * 100.0)
S23 = Float(0.230325039e-3 * 100.0)
S22 = Float(0.129690326e-1 * 100.0)
S21 = Float(0.401390832e0 * 100.0)
S20 = Float(0.535098336e1 * 100.0)

TMINLQU = Float(ZEROC - 40.0)
TMINICE = Float(ZEROC + TMINSTR)


# aer activation constants
R_AIR = Float(3.47e-3)  # m3 Pa kg-1K-1 # also used in GFDL_1M, but defined in aer

debug = True
if debug is True:
    latent_heat_vaporization = 2466500.00
    latent_heat_fusion = 333700.000
    latent_heat_sublimation = 2800200.00
    cp = 1004.68225
