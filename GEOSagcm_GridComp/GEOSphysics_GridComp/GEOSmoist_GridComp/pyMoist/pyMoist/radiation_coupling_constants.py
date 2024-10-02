"""All constants used for Radiation Coupling Port
Global constants taken from GEOS/src/Shared/@MAPL/shared/Constants"""

from ndsl.dsl.typing import Float


# Define whether or not to use CODATA 2018 constants
CODATA_2018_CONSTANTS = True

# MAPL_UNDEF is set to 1E15 in the Fortran
# We keep it as is for now to match 11.5.2 GEOS
MAPL_UNDEF = 1e15

# Math Constants
MAPL_PI_R8 = Float(3.14159265358979323846e0)
MAPL_PI = MAPL_PI_R8
MAPL_DEGREES_TO_RADIANS_R8 = MAPL_PI_R8 / Float(180.0e0)
MAPL_DEGREES_TO_RADIANS = MAPL_PI / Float(180.0)
MAPL_RADIANS_TO_DEGREES = Float(180.0e0) / MAPL_PI_R8

# Universal Constants
if CODATA_2018_CONSTANTS:
    MAPL_RUNIV = Float(8314.462618)  # J/(Kmole K)
else:
    MAPL_RUNIV = Float(8314.47)  # J/(Kmole K)

MAPL_LATENT_HEAT_VAPORIZATION = 2.4665e6  # J/kg @15C @1atm
MAPL_ALHL = MAPL_LATENT_HEAT_VAPORIZATION  # J/kg
MAPL_LATENT_HEAT_FUSION = 3.3370e5  # J/kg @1atm
MAPL_ALHF = MAPL_LATENT_HEAT_FUSION  # J/kg
MAPL_LATENT_HEAT_SUBLIMATION = MAPL_ALHL + MAPL_ALHF  # J/kg
MAPL_ALHS = MAPL_LATENT_HEAT_SUBLIMATION  # J/kg

# Earth Specific Chemistry and Thermodynamic Constants
MAPL_AIRMW = Float(28.965)  # kg/Kmole
MAPL_RDRY = MAPL_RUNIV / MAPL_AIRMW  # J/(kg K)
MAPL_CPDRY = Float(3.5) * MAPL_RDRY  # J/(kg K)
MAPL_KAPPA = MAPL_RDRY / MAPL_CPDRY  # (2.0/7.0)

####
# MAPL_RGAS
#
# The fortran gives a definition that doesn't link to the actual
# value we print out. There's a misunderstanding on where this value
# comes from on our end. We overwrite here the value with the one
# we read from Fortran
MAPL_RGAS = MAPL_RDRY  # J/(kg K) (DEPRECATED)
MAPL_RGAS = 287.052307  # Fortran f32 expected value
####

MAPL_CP = MAPL_RGAS / MAPL_KAPPA  # J/(kg K) (DEPRECATED)

# Constants for _fix_up_clouds_stencil
ALHLBCP = MAPL_ALHL / MAPL_CP
ALHSBCP = MAPL_ALHS / MAPL_CP

# Constants for cloud_effective_radius_liquid and cloud_effective_radius_ice
MAPL_TICE = Float(273.16)  # K
LIQ_RADII_PARAM = 0
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
# Fortran f32 exepcted value
# LBX = 47.7948799
LBE = Float(1.0) / Float(3.0) - Float(0.14)
# Fortran f32 exepcted value
# LBE = 0.193333343

# Taken from GEOS_GFDL_1M_InterfaceMod.F90:269-274
# - uses MAPL_GetResource to get variables.
# These do not change in the radiation coupling loop
MIN_RL = Float(2.5e-6)
MAX_RL = Float(60.0e-6)
FAC_RL = Float(1.0)
MIN_RI = Float(5.0e-6)
MAX_RI = Float(100.0e-6)
FAC_RI = Float(1.0)
