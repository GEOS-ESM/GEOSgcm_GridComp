"""
Constants for GFDL_1M Driver. Some of these are redefined
(with different values) from the generic GFDL constants,
others are duplicated. Both of these need a proper solution
Commented comments are defined in namelist and set at stencil init.
"""

import numpy as np
from gt4py.cartesian.gtscript import i32

from ndsl.dsl.typing import Float


# constants from gfdl_1m_driver module

grav = Float(9.80665)  # gfs: acceleration due to gravity
rdgas = Float(287.05)  # gfs: gas constant for dry air
rvgas = Float(461.50)  # gfs: gas constant for water vapor
cp_air = Float(1004.6)  # gfs: heat capacity of dry air at constant pressure
hlv = Float(2.5e6)  # gfs: latent heat of evaporation
hlf = Float(3.3358e5)  # gfs: latent heat of fusion
pi = Float(3.1415926535897931)  # gfs: ratio of circle circumference to diameter
cp_vap = (
    Float(4.0) * rvgas
)  # 1846.0, heat capacity of water vapore at constnat pressure
cv_air = cp_air - rdgas  # 717.55, heat capacity of dry air at constant volume
cv_vap = Float(3.0) * rvgas  # 1384.5, heat capacity of water vapor at constant volume

c_ice = Float(1972.0)  # gfdl: heat capacity of ice at - 15 deg c
c_liq = Float(4185.5)  # gfdl: heat capacity of water at 15 deg c

eps = rdgas / rvgas  # 0.6219934995
zvir = rvgas / rdgas - Float(1.0)  # 0.6077338443

t_ice = Float(273.16)  # freezing temperature
table_ice = Float(273.16)  # freezing point for qs table

e00 = Float(611.21)  # ifs: saturation vapor pressure at 0 deg c

dc_vap = cp_vap - c_liq  # - 2339.5, isobaric heating / cooling
dc_ice = c_liq - c_ice  # 2213.5, isobaric heating / colling

hlv0 = hlv  # gfs: evaporation latent heat coefficient at 0 deg c
hlf0 = hlf  # gfs: fussion latent heat coefficient at 0 deg c

lv0 = (
    hlv0 - dc_vap * t_ice
)  # 3.13905782e6, evaporation latent heat coefficient at 0 deg k
li00 = hlf0 - dc_ice * t_ice  # - 2.7105966e5, fusion latent heat coefficient at 0 deg k

d2ice = dc_vap + dc_ice  # - 126, isobaric heating / cooling
li2 = lv0 + li00  # 2.86799816e6, sublimation latent heat coefficient at 0 deg k

qpmin = Float(1.0e-8)  # min value for suspended rain/snow/liquid/ice precip
qvmin = Float(1.0e-20)  # min value for water vapor (treated as zero)
qcmin = Float(1.0e-12)  # min value for cloud condensates

vr_min = Float(1.0e-3)  # min fall speed for rain
vf_min = Float(1.0e-5)  # min fall speed for cloud ice, snow, graupel

dz_min = Float(1.0e-2)  # use for correcting flipped height

sfcrho = Float(1.2)  # surface air density
rhor = Float(1.0e3)  # density of rain water, lin83

rc = (Float(4.0) / Float(3.0)) * pi * rhor

# cloud microphysics switchers

do_setup = True  # setup constants and parameters
p_nonhydro = False  # perform hydrosatic adjustment on air density

ables_are_initialized = False

dt_fr = Float(8.0)  # epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
# minimum temperature water can exist (moore & molinero nov. 2011, nature)
# dt_fr can be considered as the error bar

p_min = Float(100.0)  # minimum pressure (pascal) for mp to operate

# -----------------------------------------------------------------------
# namelist parameters
# -----------------------------------------------------------------------

tice = Float(273.16)  # set tice = 165. to trun off ice - phase phys (kessler emulator)

log_10 = np.log(10.0, dtype=Float)
tice0 = Float(273.16) - Float(0.01)
t_wfr = Float(273.16) - Float(
    40.0
)  # supercooled water can exist down to - 40 c, which is the "absolute"

# Constants moved from setup functions
rgrav = Float(1.0) / grav

# fall velocity constants:

thi = Float(1.0e-8)  # cloud ice threshold for terminal fall
thg = Float(1.0e-8)
ths = Float(1.0e-8)
aaC = Float(-4.18334e-5)
bbC = Float(-0.00525867)
ccC = Float(-0.0486519)
ddC = Float(0.00251197)
eeC = Float(1.91523)
aaL = Float(-1.70704e-5)
bbL = Float(-0.00319109)
ccL = Float(-0.0169876)
ddL = Float(0.00410839)
eeL = Float(1.93644)

# marshall - palmer constants

vcons = Float(6.6280504)
vcong = Float(87.2382675)
norms = Float(942477796.076938)
normg = Float(5026548245.74367)


# From setupm

gam263 = Float(1.456943)
gam275 = Float(1.608355)
gam290 = Float(1.827363)
gam325 = Float(2.54925)
gam350 = Float(3.323363)
gam380 = Float(4.694155)
gam425 = Float(8.285063)
gam450 = Float(11.631769)
gam480 = Float(17.837789)
gam625 = Float(184.860962)
gam680 = Float(496.604067)

# intercept parameters

rnzr = Float(8.0e6)  # lin83
rnzs = Float(3.0e6)  # lin83
rnzg = Float(4.0e6)  # rh84

# density parameters

rhos = Float(0.1e3)  # lin83 (snow density; 1 / 10 of water)
rhog = Float(0.4e3)  # rh84 (graupel density)
acc = [Float(5.0), Float(2.0), Float(0.5)]

# computed constants

pie = Float(4.0) * np.arctan(1.0, dtype=Float)

vdifu = Float(2.11e-5)
tcond = Float(2.36e-2)

visk = Float(1.259e-5)
hlts = Float(2.8336e6)
hltc = Float(2.5e6)
hltf = Float(3.336e5)

ch2o = Float(4.1855e3)
ri50 = Float(1.0e-4)

pisq = pie * pie
scm3 = (visk / vdifu) ** (Float(1.0) / Float(3.0))

cracs = pisq * rnzr * rnzs * rhos
csacr = pisq * rnzr * rnzs * rhor
cgacr = pisq * rnzr * rnzg * rhor

act = np.zeros(8)

act[0] = pie * rnzs * rhos
act[1] = pie * rnzr * rhor
act[5] = pie * rnzg * rhog
act[2] = act[1]
act[3] = act[0]
act[4] = act[1]
act[6] = act[0]
act[7] = act[5]

act_0 = act[0]
act_1 = act[1]
act_2 = act[2]
act_3 = act[3]
act_4 = act[4]
act_5 = act[5]
act_6 = act[6]
act_7 = act[7]

acco = np.zeros([3, 4])
for i in range(1, 4):
    for k in range(1, 5):
        acco[i - 1, k - 1] = acc[i - 1] / (
            act[2 * k - 2] ** ((7 - i) * Float(0.25))
            * act[2 * k - 1] ** (i * Float(0.25))
        )

acco_00 = acco[0, 0]
acco_01 = acco[0, 1]
acco_02 = acco[0, 2]
acco_03 = acco[0, 3]
acco_10 = acco[1, 0]
acco_11 = acco[1, 1]
acco_12 = acco[1, 2]
acco_13 = acco[1, 3]
acco_20 = acco[2, 0]
acco_21 = acco[2, 1]
acco_22 = acco[2, 2]
acco_23 = acco[2, 3]


gcon = Float(40.74) * np.sqrt(sfcrho)  # 44.628

# subl and revp: five constants for three separate processes


es0 = Float(6.107799961e2)  # ~6.1 mb
ces0 = eps * es0


# terinal_fall / warm_rain constants

zs = Float(0)


# warm_rain constants:

vconr = Float(2503.23638966667)
normr = Float(25132741228.7183)
thr = Float(1.0e-8)

so3 = Float(7.0) / Float(3.0)


# q table constants

length = i32(2621)
