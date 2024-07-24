"""All constants used for Radiation Coupling Port. Global constants taken from GEOS/src/Shared/@MAPL/shared/Constants"""

# Define whether or not to use CODATA 2018 constants
CODATA_2018_CONSTANTS = True

# Math Constants
MAPL_PI_R8 = 3.14159265358979323846e0
MAPL_PI = MAPL_PI_R8
MAPL_DEGREES_TO_RADIANS_R8 = MAPL_PI_R8 / 180.0e0
MAPL_DEGREES_TO_RADIANS = MAPL_PI / 180.0
MAPL_RADIANS_TO_DEGREES = 180.0e0 / MAPL_PI_R8

# Universal Constants
if CODATA_2018_CONSTANTS:
    MAPL_RUNIV = 8314.462618  # J/(Kmole K)
else:
    MAPL_RUNIV = 8314.47  # J/(Kmole K)

MAPL_LATENT_HEAT_VAPORIZATION = 2.4665e6  # J/kg @15C @1atm
MAPL_ALHL = MAPL_LATENT_HEAT_VAPORIZATION  # J/kg
MAPL_LATENT_HEAT_FUSION = 3.3370e5  # J/kg @1atm
MAPL_ALHF = MAPL_LATENT_HEAT_FUSION  # J/kg
MAPL_LATENT_HEAT_SUBLIMATION = MAPL_ALHL + MAPL_ALHF  # J/kg
MAPL_ALHS = MAPL_LATENT_HEAT_SUBLIMATION  # J/kg

# Earth Specific Chemistry and Thermodynamic Constants
MAPL_AIRMW = 28.965  # kg/Kmole
MAPL_RDRY = MAPL_RUNIV / MAPL_AIRMW  # J/(kg K)
MAPL_CPDRY = 3.5 * MAPL_RDRY  # J/(kg K)
MAPL_KAPPA = MAPL_RDRY / MAPL_CPDRY  # (2.0/7.0)
MAPL_RGAS = MAPL_RDRY  # J/(kg K) (DEPRECATED)
MAPL_CP = MAPL_RGAS / MAPL_KAPPA  # J/(kg K) (DEPRECATED)

# Constants for _fix_up_clouds_stencil
ALHLBCP = MAPL_ALHL / MAPL_CP
ALHSBCP = MAPL_ALHS / MAPL_CP

# Constants for cloud_effective_radius_liquid and cloud_effective_radius_ice
MAPL_RGAS = 287.0
MAPL_TICE = 273.16  # K
LIQ_RADII_PARAM = 0
ICE_RADII_PARAM = 1
BX = 100.0 * (3.0 / (4.0 * MAPL_PI)) ** (1.0 / 3.0)
R13BBETA = 1.0 / 3.0 - 0.14
ABETA = 0.07
RHO_W = 1000.0
LDISS = 0.07
LK = 0.75
LBX = LDISS * 1.0e3 * (3.0 / (4.0 * MAPL_PI * LK * RHO_W * 1.0e-3)) ** (1.0 / 3.0)
LBE = 1.0 / 3.0 - 0.14

# Taken from GEOS_GFDL_1M_InterfaceMod.F90:269-274 - uses MAPL_GetResource to get variables. These do not change in the radiation coupling loop
MIN_RL = 2.5e-6
MAX_RL = 60.0e-6
FAC_RL = 1.0
MIN_RI = 5.0e-6
MAX_RI = 100.0e-6
FAC_RI = 1.0
