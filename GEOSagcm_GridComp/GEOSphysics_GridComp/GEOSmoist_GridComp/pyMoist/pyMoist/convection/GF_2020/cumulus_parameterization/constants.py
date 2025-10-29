from ndsl.dsl.typing import Float, Int


# physical constants used ONLY for GF_2020, else there may be unintended conflicts with global constants
RGAS = Float(287.0)  # J K-1 kg-1
CP = Float(1004.0)  # J K-1 kg-1
RV = Float(461.0)  # J K-1 kg-1
P00 = Float(1.0e5)  # hPa
TCRIT = Float(258.0)  # K
CPOR = CP / RGAS
XLV = Float(2.5e6)  # J kg-1
AKMIN = Float(1.0)  # #
TKMIN = Float(1.0e-5)  # m+2 s-2
CCNCLEAN = Float(250.0)  # # cm-3
T_0 = Float(273.16)  # K
T_ICE = Float(235.16)  # K
XLF = Float(0.333e6)  # latent heat of freezing (J K-1 kg-1)
MAX_QSAT = Float(0.5)  # kg/kg
MX_BUOY = CP * Float(5.0) + XLV * Float(2.0e-3)  # temp exc=5 K, q deficit=2 g/kg (=> mx_buoy ~ 10 kJ/kg)

smaller_qv = Float(1.0e-16)  # kg/kg

MAXENS1 = Int(1)  # ensemble one on cap_max
MAXENS2 = Int(1)  # ensemble two on precip efficiency
MAXENS3 = Int(16)  # ensemble three done in cup_forcing_ens16 for G3d

deep = Int(1)
shallow = Int(2)
mid = Int(3)
