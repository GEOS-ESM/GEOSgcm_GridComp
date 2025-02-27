"""
Constants for GFDL_1M Driver. Some of these are redefined
(with different values) from the generic GFDL constants,
others are duplicated. Both of these need a proper solution
Commented comments are defined in namelist and set at stencil init.
"""

import numpy as np

from ndsl.dsl.typing import Float, Int


# constants from driver module

GRAV = Float(9.80665)  # gfs: acceleration due to gravity
RDGAS = Float(287.05)  # gfs: gas constant for dry air
RVGAS = Float(461.50)  # gfs: gas constant for water vapor
CP_AIR = Float(1004.6)  # gfs: heat capacity of dry air at constant pressure
HLV = Float(2.5e6)  # gfs: latent heat of evaporation
HLF = Float(3.3358e5)  # gfs: latent heat of fusion
PI = Float(3.1415926535897931)  # gfs: ratio of circle circumference to diameter
CP_VAP = (
    Float(4.0) * RVGAS
)  # 1846.0, heat capacity of water vapore at constnat pressure
CV_AIR = CP_AIR - RDGAS  # 717.55, heat capacity of dry air at constant volume
CV_VAP = Float(3.0) * RVGAS  # 1384.5, heat capacity of water vapor at constant volume

C_ICE = Float(1972.0)  # gfdl: heat capacity of ice at - 15 deg c
C_LIQ = Float(4185.5)  # gfdl: heat capacity of water at 15 deg c

EPS = RDGAS / RVGAS  # 0.6219934995
ZVIR = RVGAS / RDGAS - Float(1.0)  # 0.6077338443

T_ICE = Float(273.16)  # freezing temperature
TABLE_ICE = Float(273.16)  # freezing point for qs table

E_00 = Float(611.21)  # ifs: saturation vapor pressure at 0 deg c

DC_VAP = CP_VAP - C_LIQ  # - 2339.5, isobaric heating / cooling
DC_ICE = C_LIQ - C_ICE  # 2213.5, isobaric heating / colling

HLV0 = HLV  # gfs: evaporation latent heat coefficient at 0 deg c
HLF0 = HLF  # gfs: fussion latent heat coefficient at 0 deg c

LV0 = (
    HLV0 - DC_VAP * T_ICE
)  # 3.13905782e6, evaporation latent heat coefficient at 0 deg k
LI00 = HLF0 - DC_ICE * T_ICE  # - 2.7105966e5, fusion latent heat coefficient at 0 deg k

D2ICE = DC_VAP + DC_ICE  # - 126, isobaric heating / cooling
LI2 = LV0 + LI00  # 2.86799816e6, sublimation latent heat coefficient at 0 deg k

QPMIN = Float(1.0e-8)  # min value for suspended rain/snow/liquid/ice precip
QVMIN = Float(1.0e-20)  # min value for water vapor (treated as zero)
QCMIN = Float(1.0e-12)  # min value for cloud condensates

VR_MIN = Float(1.0e-3)  # min fall speed for rain
VF_MIN = Float(1.0e-5)  # min fall speed for cloud ice, snow, graupel

DZ_MIN = Float(1.0e-2)  # use for correcting flipped height

SFCRHO = Float(1.2)  # surface air density
RHOR = Float(1.0e3)  # density of rain water, lin83

RC = (Float(4.0) / Float(3.0)) * PI * RHOR

# cloud microphysics switches

DO_SETUP = True  # setup constants and parameters
P_NONHYDRO = False  # perform hydrosatic adjustment on air density

DT_FR = Float(8.0)  # epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
# minimum temperature water can exist (moore & molinero nov. 2011, nature)
# dt_fr can be considered as the error bar

P_MIN = Float(100.0)  # minimum pressure (pascal) for mp to operate

# -----------------------------------------------------------------------
# namelist parameters
# -----------------------------------------------------------------------

TICE = Float(273.16)  # set tice = 165. to trun off ice - phase phys (kessler emulator)

LOG_10 = np.log(10.0, dtype=Float)
TICE0 = Float(273.16) - Float(0.01)
T_WFR = Float(273.16) - Float(
    40.0
)  # supercooled water can exist down to - 40 c, which is the "absolute"

# Constants moved from setup functions
RGRAV = Float(1.0) / GRAV

# fall velocity constants:

THI = Float(1.0e-8)  # cloud ice threshold for terminal fall
THG = Float(1.0e-8)
THS = Float(1.0e-8)
AAC = Float(-4.18334e-5)
BBC = Float(-0.00525867)
CCC = Float(-0.0486519)
DDC = Float(0.00251197)
EEC = Float(1.91523)
AAL = Float(-1.70704e-5)
BBL = Float(-0.00319109)
CCL = Float(-0.0169876)
DDL = Float(0.00410839)
EEL = Float(1.93644)

# marshall - palmer constants

VCONS = Float(6.6280504)
VCONG = Float(87.2382675)
NORMS = Float(942477796.076938)
NORMG = Float(5026548245.74367)


# From setupm

GAM263 = Float(1.456943)
GAM275 = Float(1.608355)
GAM209 = Float(1.827363)
GAM325 = Float(2.54925)
GAM350 = Float(3.323363)
GAM380 = Float(4.694155)
GAM425 = Float(8.285063)
GAM450 = Float(11.631769)
GAM480 = Float(17.837789)
GAM625 = Float(184.860962)
GAM680 = Float(496.604067)

# intercept parameters

RNZR = Float(8.0e6)  # lin83
RNZS = Float(3.0e6)  # lin83
RNZG = Float(4.0e6)  # rh84

# density parameters

RHOS = Float(0.1e3)  # lin83 (snow density; 1 / 10 of water)
RHOG = Float(0.4e3)  # rh84 (graupel density)
ACC = [Float(5.0), Float(2.0), Float(0.5)]

# computed constants

PIE = Float(4.0) * np.arctan(1.0, dtype=Float)

VDIFU = Float(2.11e-5)
TCOND = Float(2.36e-2)

VISK = Float(1.259e-5)
HLTS = Float(2.8336e6)
HLTC = Float(2.5e6)
HLTF = Float(3.336e5)

CH2O = Float(4.1855e3)
RI50 = Float(1.0e-4)

PISQ = PIE * PIE
SCM3 = (VISK / VDIFU) ** (Float(1.0) / Float(3.0))

CRACS = PISQ * RNZR * RNZS * RHOS
CSARC = PISQ * RNZR * RNZS * RHOR
CGARC = PISQ * RNZR * RNZG * RHOR

ACT = np.zeros(8, dtype=Float)

ACT[0] = PIE * RNZS * RHOS
ACT[1] = PIE * RNZR * RHOR
ACT[5] = PIE * RNZG * RHOG
ACT[2] = ACT[1]
ACT[3] = ACT[0]
ACT[4] = ACT[1]
ACT[6] = ACT[0]
ACT[7] = ACT[5]

ACT_0 = ACT[0]
ACT_1 = ACT[1]
ACT_2 = ACT[2]
ACT_3 = ACT[3]
ACT_4 = ACT[4]
ACT_5 = ACT[5]
ACT_6 = ACT[6]
ACT_7 = ACT[7]

ACCO = np.zeros([3, 4], dtype=Float)
for i in range(1, 4):
    for k in range(1, 5):
        ACCO[i - 1, k - 1] = ACC[i - 1] / (
            ACT[2 * k - 2] ** ((7 - i) * Float(0.25))
            * ACT[2 * k - 1] ** (i * Float(0.25))
        )

ACCO_00 = ACCO[0, 0]
ACCO_01 = ACCO[0, 1]
ACCO_02 = ACCO[0, 2]
ACCO_03 = ACCO[0, 3]
ACCO_10 = ACCO[1, 0]
ACCO_11 = ACCO[1, 1]
ACCO_12 = ACCO[1, 2]
ACCO_13 = ACCO[1, 3]
ACCO_20 = ACCO[2, 0]
ACCO_21 = ACCO[2, 1]
ACCO_22 = ACCO[2, 2]
ACCO_23 = ACCO[2, 3]


GCON = Float(40.74) * np.sqrt(SFCRHO, dtype=Float)  # 44.628

# subl and revp: five constants for three separate processes


ES0 = Float(6.107799961e2)  # ~6.1 mb
CES0 = EPS * ES0


# terinal_fall / warm_rain constants

ZS = Float(0)


# warm_rain constants:

VCONR = Float(2503.23638966667)
NORMR = Float(25132741228.7183)
THR = Float(1.0e-8)

SO3 = Float(7.0) / Float(3.0)


# q table constants

LENGTH = Int(2621)
DELT = Float(0.1)
TMIN = TABLE_ICE - Float(160.0)
ESBASW = Float(1013246.0)
TBASW = TABLE_ICE + Float(100.0)
ESBASI = Float(6107.1)
