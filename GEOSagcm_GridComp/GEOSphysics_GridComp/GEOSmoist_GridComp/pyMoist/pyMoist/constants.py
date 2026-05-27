"""File containing all constants used for pyMoist"""

import numpy as np

from ndsl.dsl.typing import Float


_f32 = np.float32
_f64 = np.float64
_i32 = np.int32

# Define number of tracers
NCNST = _i32(23)

# MAPL_UNDEF is set to 1E15 in the Fortran
# We keep it as is for now to match 11.5.2 GEOS
MAPL_UNDEF = Float(1e15)

# Math Constants
MAPL_PI_R8 = _f64(3.14159265358979323846e0)
MAPL_PI = _f32(MAPL_PI_R8)
MAPL_DEGREES_TO_RADIANS_R8 = MAPL_PI_R8 / _f64(180.0e0)
MAPL_DEGREES_TO_RADIANS = _f32(MAPL_PI / Float(180.0))
MAPL_RADIANS_TO_DEGREES = _f64(_f64(180.0e0) / MAPL_PI_R8)

# Following taken from PhysicalConstants.F90
# Universal Constants
CODATA_2018_CONSTANTS = False  # set for now, needs to be dynamic
if CODATA_2018_CONSTANTS:
    MAPL_STFBOL = Float(5.670374419e-8)  # W/(m^2 K^4)
    MAPL_AVOGAD = Float(6.02214076e26)  # 1/kmol
    MAPL_RUNIV = Float(8314.462618)  # J/(Kmole K)
else:
    MAPL_STFBOL = Float(5.6734e-8)  # W/(m^2 K^4)
    MAPL_AVOGAD = Float(6.023e26)  # 1/kmol
    MAPL_RUNIV = Float(8314.47)  # J/(Kmole K)

# Earth Constants
MAPL_PSDRY = _f64(98305)  # Pa
MAPL_SECONDS_PER_SIDEREAL_DAY = Float(86164.0)  # s
MAPL_GRAV = Float(9.80665)  # m^2/s
MAPL_RADIUS = Float(6371.0e3)  # m
MAPL_OMEGA_R8 = _f64(_f64(2.0) * MAPL_PI_R8 / MAPL_SECONDS_PER_SIDEREAL_DAY)  # 1/s
MAPL_OMEGA = _f32(Float(2.0) * MAPL_PI / MAPL_SECONDS_PER_SIDEREAL_DAY)  # 1/s
MAPL_EARTH_ECCENTRICITY = _f64(8.181919084262200e-2)  # --
MAPL_EARTH_SEMIMAJOR_AXIS = _f64(6378137)  # m
MAPL_KM_PER_DEG = _f64(_f64(1.0) / (MAPL_RADIUS / _f64(1000))) * MAPL_RADIANS_TO_DEGREES
MAPL_DEG_PER_KM = _f64((MAPL_RADIUS / _f64(1000)) * MAPL_DEGREES_TO_RADIANS_R8)

# Physical Properties
MAPL_LATENT_HEAT_VAPORIZATION = Float(2.4665e6)  # J/kg @15C @1atm
MAPL_ALHL = MAPL_LATENT_HEAT_VAPORIZATION  # J/kg
MAPL_LATENT_HEAT_FUSION = Float(3.3370e5)  # J/kg @1atm
MAPL_ALHF = MAPL_LATENT_HEAT_FUSION  # J/kg
MAPL_LATENT_HEAT_SUBLIMATION = MAPL_ALHL + MAPL_ALHF  # J/kg
MAPL_ALHS = MAPL_LATENT_HEAT_SUBLIMATION  # J/kg
MAPL_H2OMW = Float(18.015)  # kg/Kmole

# Earth Specific Chemistry and Thermodynamic Constants
MAPL_TICE = Float(273.16)  # K
MAPL_AIRMW = Float(28.965)  # kg/Kmole
MAPL_RDRY = MAPL_RUNIV / MAPL_AIRMW  # J/(kg K)
MAPL_CPDRY = Float(3.5) * MAPL_RDRY  # J/(kg K)
MAPL_KAPPA = MAPL_RDRY / MAPL_CPDRY  # (2.0/7.0)
MAPL_CVDRY = MAPL_CPDRY - MAPL_RDRY  # J/(kg K)
MAPL_RVAP = MAPL_RUNIV / MAPL_H2OMW  # J/(kg K)
MAPL_CPVAP = Float(4) * MAPL_RVAP  # J/(kg K)
MAPL_CVVAP = MAPL_CPVAP - MAPL_RVAP  # J/(kg K)
MAPL_RGAS = MAPL_RDRY  # MAPL_RDRY  # J/(kg K) (DEPRECATED)
MAPL_CP = MAPL_RGAS / MAPL_KAPPA  # J/(kg K) (DEPRECATED)

EPSILON = MAPL_H2OMW / MAPL_AIRMW  # --
MAPL_CELSIUS_TO_KELVIN = Float(273.15)  # K


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
JiT_ICE_ALL = Float(MAPL_TICE - 40.0)
JiT_ICE_MAX = Float(MAPL_TICE)
JiICEFRPWR = Float(4.0)

# Other miscellaneous parameters
TAUFRZ = Float(450)
K_COND = Float(2.4e-2)  # J m**-1 s**-1 K**-1
DIFFU = Float(2.2e-5)  # m**2 s**-1
dQCMAX = Float(1.0e-4)

# cloud radius Constants based on DOI 10.1088/1748-9326/3/4/045021
RHO_I = Float(916.8)  # Density of ice crystal in kg/m^3
# cloud radius eqs are in cgs units

# Constants for _fix_up_clouds_stencil
ALHLBCP = MAPL_ALHL / MAPL_CP
ALHSBCP = MAPL_ALHS / MAPL_CP

# Constants for cloud_effective_radius_liquid and cloud_effective_radius_ice
LIQ_RADII_PARAM = 2
ICE_RADII_PARAM = 1
BX = Float(100.0) * (Float(3.0) / (Float(4.0) * MAPL_PI)) ** (Float(1.0) / Float(3.0))
R13BBETA = Float(1.0) / Float(3.0) - Float(0.14)
ABETA = Float(0.07)
RHO_W = Float(1000.0)
LDISS = Float(0.07)
LK = Float(0.75)
LBX = (
    LDISS
    * Float(1.0e3)
    * (Float(3.0) / (Float(4.0) * MAPL_PI * LK * RHO_W * Float(1.0e-3)))
    ** (Float(1.0) / Float(3.0))
)
LBE = Float(1.0) / Float(3.0) - Float(0.14)

# Aer Activation constants
R_AIR = Float(3.47e-3)  # m3 Pa kg-1K-1, also used in GFDL_1M, but defined in aer


# Python euqivalent of Fortran's tiny(X)
FLOAT_TINY = np.finfo(Float).tiny

# Define how many modes in an Aerosol Activation
N_MODES = 14
