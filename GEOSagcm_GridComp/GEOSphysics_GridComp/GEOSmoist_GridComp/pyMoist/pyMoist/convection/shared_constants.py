"""Shared GEOS5/GF2020 Constants"""

from gt4py.cartesian.gtscript import i32, f32

# plume spectral size
MAXIENS = i32(3)
DEEP = i32(1)
SHAL = i32(2)
MID = i32(3)

# physical constants

RGAS = f32(287.0)  # J K-1 kg-1
CP = f32(1004.0)  # J K-1 kg-1
RV = f32(461.0)  # J K-1 kg-1
P00 = f32(1.0e5)  # hPa
TCRIT = f32(258.0)  # K
CPOR = f32(CP / RGAS)
XLV = f32(2.5e6)  # J kg-1
AKMIN = f32(1.0)  # #
TKMIN = f32(1.0e-5)  # m+2 s-2
CCNCLEAN = f32(250.0)  # # cm-3
T_0 = f32(273.16)  # K
T_ICE = f32(235.16)  # K
XLF = f32(0.333e6)  # latent heat of freezing (J K-1 kg-1)
MAX_QSAT = 0.5  # kg/kg
MX_BUOY = (
    CP * 5.0 + XLV * 2.0e-3
)  # temp exc=5 K, q deficit=2 g/kg (=> mx_buoy ~ 10 kJ/kg)

# Default autoconversion parameter for GEOS-Chem species [s-1]
KC_DEFAULT_GCC = f32(5.0e-3)

CNV_2MOM = False
