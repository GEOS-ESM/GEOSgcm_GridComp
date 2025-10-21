"""Shared GEOS5/GF2020 Constants"""

from gt4py.cartesian.gtscript import int32, float32

# plume spectral size
MAXIENS = int32(3)
DEEP = int32(1)
SHAL = int32(2)
MID = int32(3)

# physical constants

MAPL_GRAV = float32(9.80665)  # m^2/s
RGAS = float32(287.0)  # J K-1 kg-1
CP = float32(1004.0)  # J K-1 kg-1
RV = float32(461.0)  # J K-1 kg-1
P00 = float32(1.0e5)  # hPa
TCRIT = float32(258.0)  # K
CPOR = float32(CP / RGAS)
XLV = float32(2.5e6)  # J kg-1
AKMIN = float32(1.0)  # #
TKMIN = float32(1.0e-5)  # m+2 s-2
CCNCLEAN = float32(250.0)  # # cm-3
T_0 = float32(273.16)  # K
T_ICE = float32(235.16)  # K
XLF = float32(0.333e6)  # latent heat of freezing (J K-1 kg-1)
MAX_QSAT = 0.5  # kg/kg
MX_BUOY = (
    CP * 5.0 + XLV * 2.0e-3
)  # temp exc=5 K, q deficit=2 g/kg (=> mx_buoy ~ 10 kJ/kg)

# Default autoconversion parameter for GEOS-Chem species [s-1]
KC_DEFAULT_GCC = float32(5.0e-3)

CNV_2MOM = False
