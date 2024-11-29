"""
Constants for GFDL_1M Driver. Some of these are redefined (with different values)
from the generic GFDL constants, others are duplicated. Both of these need a proper solution
Commented comments are defined in namelist and set at stencil init.
"""

import numpy as np

### constants from gfdl_1m_driver module

grav = 9.80665  # gfs: acceleration due to gravity
rdgas = 287.05  # gfs: gas constant for dry air
rvgas = 461.50  # gfs: gas constant for water vapor
cp_air = 1004.6  # gfs: heat capacity of dry air at constant pressure
hlv = 2.5e6  # gfs: latent heat of evaporation
hlf = 3.3358e5  # gfs: latent heat of fusion
pi = 3.1415926535897931  # gfs: ratio of circle circumference to diameter
cp_vap = 4.0 * rvgas  # 1846.0, heat capacity of water vapore at constnat pressure
cv_air = cp_air - rdgas  # 717.55, heat capacity of dry air at constant volume
cv_vap = 3.0 * rvgas  # 1384.5, heat capacity of water vapor at constant volume

c_ice = 1972.0  # gfdl: heat capacity of ice at - 15 deg c
c_liq = 4185.5  # gfdl: heat capacity of water at 15 deg c

eps = rdgas / rvgas  # 0.6219934995
zvir = rvgas / rdgas - 1.0  # 0.6077338443

t_ice = 273.16  # freezing temperature
table_ice = 273.16  # freezing point for qs table

e00 = 611.21  # ifs: saturation vapor pressure at 0 deg c

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

qpmin = 1.0e-8  # min value for suspended rain/snow/liquid/ice precip
qvmin = 1.0e-20  # min value for water vapor (treated as zero)
qcmin = 1.0e-12  # min value for cloud condensates

vr_min = 1.0e-3  # min fall speed for rain
vf_min = 1.0e-5  # min fall speed for cloud ice, snow, graupel

dz_min = 1.0e-2  # use for correcting flipped height

sfcrho = 1.2  # surface air density
rhor = 1.0e3  # density of rain water, lin83

rc = (4.0 / 3.0) * pi * rhor

# cloud microphysics switchers

# icloud_f = 0  # cloud scheme
# irain_f = 0  # cloud water to rain auto conversion scheme

# de_ice = False  # to prevent excessive build - up of cloud ice from external sources
# sedi_transport = False  # transport of momentum in sedimentation
# do_sedi_w = False  # transport of vertical motion in sedimentation
# do_sedi_heat = False  # transport of heat in sedimentation
# prog_ccn = False  # do prognostic ccn (yi ming's method)
# do_bigg = False  # do bigg mechanism freezing of supercooled liquid on aerosol nuclei
# do_evap = False  # do evaporation
# do_subl = False  # do sublimation
# do_qa = False  # do inline cloud fraction (WMP: in FV3 dynamics)
# preciprad = True  # consider precipitates in cloud fraciton calculation
# fix_negative = False  # fix negative water species
do_setup = True  # setup constants and parameters
p_nonhydro = False  # perform hydrosatic adjustment on air density

ables_are_initialized = False

dt_fr = 8.0  # epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
# minimum temperature water can exist (moore & molinero nov. 2011, nature)
# dt_fr can be considered as the error bar

p_min = 100.0  # minimum pressure (pascal) for mp to operate

# -----------------------------------------------------------------------
# namelist parameters
# -----------------------------------------------------------------------

# cld_min = 0.05  # minimum cloud fraction
tice = 273.16  # set tice = 165. to trun off ice - phase phys (kessler emulator)

log_10 = np.log(10.0)
tice0 = 273.16 - 0.01
t_wfr = (
    273.16 - 40.0
)  # supercooled water can exist down to - 40 c, which is the "absolute"

# t_min = 178.0  # min temp to freeze - dry all water vapor
# t_sub = 184.0  # min temp for sublimation of cloud ice
# mp_time = 150.0  # maximum micro - physics time step (sec)

# relative humidity increment

# rh_inc = 0.25  # rh increment for complete evaporation of cloud water and cloud ice
# rh_inr = 0.25  # rh increment for minimum evaporation of rain
# rh_ins = 0.25  # rh increment for sublimation of snow

# conversion time scale

# tau_r2g = 900.0  # rain freezing during fast_sat
# tau_smlt = 900.0  # snow melting
# tau_g2r = 600.0  # graupel melting to rain
# tau_imlt = 600.0  # cloud ice melting
# tau_i2s = 1000.0  # cloud ice to snow auto - conversion
# tau_l2r = 900.0  # cloud water to rain auto - conversion
# tau_v2l = 150.0  # water vapor to cloud water (condensation)
# tau_l2v = 300.0  # cloud water to water vapor (evaporation)
# tau_i2v = 300.0  # cloud ice to water vapor (sublimation)
# tau_s2v = 600.0  # snow sublimation
# tau_v2s = 21600.0  # snow deposition -- make it a slow process
# tau_g2v = 900.0  # graupel sublimation
# tau_v2g = 21600.0  # graupel deposition -- make it a slow process
# tau_revp = 600.0  # rain re-evaporation
# tau_frz = 450.0  # timescale for liquid-ice freezing

# horizontal subgrid variability

# dw_land = 0.20  # base value for subgrid deviation / variability over land
# dw_ocean = 0.10  # base value for ocean

# prescribed ccn

# ccn_o = 90.0  # ccn over ocean (cm^ - 3)
# ccn_l = 270.0  # ccn over land (cm^ - 3)

# rthreshu = 7.0e-6  # critical cloud drop radius (micro m)
# rthreshs = 10.0e-6  # critical cloud drop radius (micro m)

# sat_adj0 = 0.90  # adjustment factor (0: no, 1: full) during fast_sat_adj

# qc_crt = 5.0e-8  # mini condensate mixing ratio to allow partial cloudiness

# qi_lim = 1.0  # cloud ice limiter to prevent large ice build up

# ql_mlt = 2.0e-3  # max value of cloud water allowed from melted cloud ice
# qs_mlt = 1.0e-6  # max cloud water due to snow melt

# ql_gen = 1.0e-3  # max cloud water generation
# qi_gen = 9.82679e-5  # max cloud ice generation at -40 C

# cloud condensate upper bounds: "safety valves" for ql & qi

# ql0_max = 2.0e-3  # max cloud water value (auto converted to rain)
# qi0_max = 1.0e-4  # max cloud ice value (by other sources)

# qi0_crt = 1.0e-4  # cloud ice to snow autoconversion threshold (was 1.e-4)
# qi0_crt is highly dependent on horizontal resolution
# qr0_crt = 1.0e-4  # rain to snow or graupel / hail threshold
# lfo used * mixing ratio * = 1.e-4 (hail in lfo)
# qs0_crt = 1.0e-3  # snow to graupel density threshold (0.6e-3 in purdue lin scheme)

# c_paut = 0.55  # autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
# c_psaci = 0.02  # accretion: cloud ice to snow (was 0.1 in zetac)
# c_piacr = 5.0  # accretion: rain to ice:
# c_cracw = 0.9  # rain accretion efficiency
# c_pgacs = 2.0e-3  # snow to graupel "accretion" eff. (was 0.1 in zetac)
# c_pgaci = 0.05  #  ice to graupel "accretion" eff.

# decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

# alin = 842.0  # "a" in lin1983
# clin = 4.8  # "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)

# fall velocity tuning constants:

# const_vi = False  # if .t. the constants are specified by v * _fac
# const_vs = False  # if .t. the constants are specified by v * _fac
# const_vg = False  # if .t. the constants are specified by v * _fac
# const_vr = False  # if .t. the constants are specified by v * _fac

# good values:

# vi_fac = 1.0  # if const_vi: 1 / 3
# vs_fac = 1.0  # if const_vs: 1.
# vg_fac = 1.0  # if const_vg: 2.
# vr_fac = 1.0  # if const_vr: 4.

# upper bounds of fall speed (with variable speed option)

# vi_max = 1.0  # max fall speed for ice
# vs_max = 2.0  # max fall speed for snow
# vg_max = 12.0  # max fall speed for graupel
# vr_max = 12.0  # max fall speed for rain

# cloud microphysics switchers

# fast_sat_adj = False  # has fast saturation adjustments
# z_slope_liq = True  # use linear mono slope for autocconversions
# z_slope_ice = False  # use linear mono slope for autocconversions
# use_ccn = False  # use input ccn when .T. else use ccn_o/ccn_l
# use_ppm = False  # use ppm fall scheme
# mono_prof = True  # perform terminal fall with mono ppm scheme
# mp_print = False  # cloud microphysics debugging printout

# Constants moved from setup functions
rgrav = 1.0 / grav

### fall velocity constants:

thi = 1.0e-8  # cloud ice threshold for terminal fall
thg = 1.0e-8
ths = 1.0e-8
aaC = -4.18334e-5
bbC = -0.00525867
ccC = -0.0486519
ddC = 0.00251197
eeC = 1.91523
aaL = -1.70704e-5
bbL = -0.00319109
ccL = -0.0169876
ddL = 0.00410839
eeL = 1.93644

# marshall - palmer constants

vcons = 6.6280504
vcong = 87.2382675
norms = 942477796.076938
normg = 5026548245.74367


### From setupm

gam263 = 1.456943
gam275 = 1.608355
gam290 = 1.827363
gam325 = 2.54925
gam350 = 3.323363
gam380 = 4.694155
gam425 = 8.285063
gam450 = 11.631769
gam480 = 17.837789
gam625 = 184.860962
gam680 = 496.604067

# intercept parameters

rnzr = 8.0e6  # lin83
rnzs = 3.0e6  # lin83
rnzg = 4.0e6  # rh84

# density parameters

rhos = 0.1e3  # lin83 (snow density; 1 / 10 of water)
rhog = 0.4e3  # rh84 (graupel density)
acc = [5.0, 2.0, 0.5]

# computed constants

pie = 4.0 * np.arctan(1.0)

vdifu = 2.11e-5
tcond = 2.36e-2

visk = 1.259e-5
hlts = 2.8336e6
hltc = 2.5e6
hltf = 3.336e5

ch2o = 4.1855e3
ri50 = 1.0e-4

pisq = pie * pie
scm3 = (visk / vdifu) ** (1.0 / 3.0)

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
            act[2 * k - 2] ** ((7 - i) * 0.25) * act[2 * k - 1] ** (i * 0.25)
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


gcon = 40.74 * np.sqrt(sfcrho)  # 44.628

# subl and revp: five constants for three separate processes


es0 = 6.107799961e2  # ~6.1 mb
ces0 = eps * es0


# terinal_fall / warm_rain constants

zs = 0


# warm_rain constants:

vconr = 2503.23638966667
normr = 25132741228.7183
thr = 1.0e-8

so3 = 7.0 / 3.0
