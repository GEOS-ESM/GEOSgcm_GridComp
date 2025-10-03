from ndsl.dsl.typing import Float
from pyMoist.constants import MAPL_GRAV

from dataclasses import dataclass


@dataclass
class DriverConstants:
    rgas = Float(287.0)  # J K-1 kg-1
    cp = Float(1004.0)  # J K-1 kg-1
    rv = Float(461.0)  # J K-1 kg-1
    p00 = Float(1.0e5)  # hPa
    tcrit = Float(258.0)  # K
    g = MAPL_GRAV  # m s-2
    cpor = cp / rgas
    xlv = Float(2.5e6)  # J kg-1
    akmin = Float(1.0)  #
    tkmin = Float(1.0e-5)  # m+2 s-2
    ccnclean = Float(250.0)  # cm-3
    T_0 = Float(273.16)  # K
    T_ice = Float(235.16)  # K
    xlf = Float(0.333e6)  # latent heat of freezing (J K-1 kg-1)
    max_qsat = Float(0.5)  # kg/kg
    mx_buoy = cp * 5.0 + xlv * 2.0e-3  # temp exc=5 K, q deficit=2 g/kg (=> mx_buoy ~ 10 kJ/kg)


driver_constants = DriverConstants()
