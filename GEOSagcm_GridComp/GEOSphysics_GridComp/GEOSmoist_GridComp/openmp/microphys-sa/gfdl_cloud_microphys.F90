!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Cloud Microphysics.
!*
!* The GFDL Cloud Microphysics is free software: you can
!* redistribute it and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The GFDL Cloud Microphysics is distributed in the hope it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the GFDL Cloud Microphysics.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'gfdl_cloud_microphys' contains the full GFDL cloud
!! microphysics \cite chen2013seasonal.
!>@details The module is paired with 'fv_cmp', which performs the "fast"
!! processes
!>author Shian-Jiann Lin, Linjiong Zhou

! =======================================================================
! cloud micro - physics package for gfdl global cloud resolving model
! the algorithms are originally derived from lin et al 1983. most of the
! key elements have been simplified / improved. this code at this stage
! bears little to no similarity to the original lin mp in zetac.
! therefore, it is best to be called gfdl micro - physics (gfdl mp) .
! developer: shian-jiann lin, linjiong zhou
! =======================================================================

module gfdl2_cloud_microphys_mod

  use omp_lib

  implicit none

  private

  public gfdl_cloud_microphys_driver, gfdl_cloud_microphys_init

  real :: missing_value = - 1.e10

  logical :: module_is_initialized = .false.
  logical :: qsmith_tables_initialized = .false.

  character (len = 17) :: mod_name = 'gfdl_cloud_microphys'

  real, parameter :: MAPL_PI = 4.0 * atan(1.0)
  real, parameter :: grav = 9.80665 !< gfs: acceleration due to gravity
  real, parameter :: rdgas = 287.05 !< gfs: gas constant for dry air
  real, parameter :: rvgas = 461.50 !< gfs: gas constant for water vapor
  real, parameter :: cp_air = 1004.6 !< gfs: heat capacity of dry air at constant pressure
  real, parameter :: hlv = 2.5e6 !< gfs: latent heat of evaporation
  real, parameter :: hlf = 3.3358e5 !< gfs: latent heat of fusion
  real, parameter :: pi = 3.1415926535897931 !< gfs: ratio of circle circumference to diameter
  ! real, parameter :: cp_air = rdgas * 7. / 2. ! 1004.675, heat capacity of dry air at constant pressure
  real, parameter :: cp_vap = 4.0 * rvgas !< 1846.0, heat capacity of water vapore at constnat pressure
  ! real, parameter :: cv_air = 717.56 ! satoh value
  real, parameter :: cv_air = cp_air - rdgas !< 717.55, heat capacity of dry air at constant volume
  ! real, parameter :: cv_vap = 1410.0 ! emanuel value
  real, parameter :: cv_vap = 3.0 * rvgas !< 1384.5, heat capacity of water vapor at constant volume

  ! the following two are from emanuel's book "atmospheric convection"
  ! real, parameter :: c_ice = 2106.0 ! heat capacity of ice at 0 deg c: c = c_ice + 7.3 * (t - tice)
  ! real, parameter :: c_liq = 4190.0 ! heat capacity of water at 0 deg c

  real, parameter :: c_ice = 1972.0 !< gfdl: heat capacity of ice at - 15 deg c
  real, parameter :: c_liq = 4185.5 !< gfdl: heat capacity of water at 15 deg c
  ! real, parameter :: c_liq = 4218.0 ! ifs: heat capacity of liquid at 0 deg c

  real, parameter :: eps = rdgas / rvgas ! 0.6219934995
  real, parameter :: zvir = rvgas / rdgas - 1. !< 0.6077338443

  real, parameter :: t_ice = 273.16 !< freezing temperature
  real, parameter :: table_ice = 273.16 !< freezing point for qs table

  ! real, parameter :: e00 = 610.71 ! gfdl: saturation vapor pressure at 0 deg c
  real, parameter :: e00 = 611.21 !< ifs: saturation vapor pressure at 0 deg c

  real, parameter :: dc_vap = cp_vap - c_liq !< - 2339.5, isobaric heating / cooling
  real, parameter :: dc_ice = c_liq - c_ice !< 2213.5, isobaric heating / colling

  real, parameter :: hlv0 = hlv !< gfs: evaporation latent heat coefficient at 0 deg c
  ! real, parameter :: hlv0 = 2.501e6 ! emanuel appendix - 2
  real, parameter :: hlf0 = hlf !< gfs: fussion latent heat coefficient at 0 deg c
  ! real, parameter :: hlf0 = 3.337e5 ! emanuel

  real, parameter :: lv0 = hlv0 - dc_vap * t_ice!< 3.13905782e6, evaporation latent heat coefficient at 0 deg k
  real, parameter :: li00 = hlf0 - dc_ice * t_ice!< - 2.7105966e5, fusion latent heat coefficient at 0 deg k

  real, parameter :: d2ice = dc_vap + dc_ice !< - 126, isobaric heating / cooling
  real, parameter :: li2 = lv0 + li00 !< 2.86799816e6, sublimation latent heat coefficient at 0 deg k

  real, parameter :: qpmin = 1.e-8  !< min value for suspended rain/snow/liquid/ice precip
  real, parameter :: qvmin = 1.e-20 !< min value for water vapor (treated as zero)
  real, parameter :: qcmin = 1.e-12 !< min value for cloud condensates

  real, parameter :: vr_min = 1.e-3 !< min fall speed for rain
  real, parameter :: vf_min = 1.e-5 !< min fall speed for cloud ice, snow, graupel

  real, parameter :: dz_min = 1.e-2 ! use for correcting flipped height

  real, parameter :: sfcrho = 1.2 !< surface air density
  real, parameter :: rhor = 1.e3 !< density of rain water, lin83

  real, parameter :: rc = (4. / 3.) * pi * rhor

  integer, parameter :: TABLE_LENGTH = 2621

  real :: cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw !< constants for accretions
  real :: acco (3, 4) !< constants for accretions
  real :: cssub (5), cgsub (5), crevp (5), cgfr (2), csmlt (5), cgmlt (5)

  real :: es0, ces0
  real :: pie, rgrav
  real :: c_air, c_vap

  real :: lati, latv, lats, lat2, lcp, icp, tcp !< used in bigg mechanism and wet bulk

  real :: d0_vap !< the same as dc_vap, except that cp_vap can be cp_vap or cv_vap
  real :: lv00 !< the same as lv0, except that cp_vap can be cp_vap or cv_vap

  ! cloud microphysics switchers

  integer :: icloud_f = 0 !< cloud scheme
  integer :: irain_f = 0 !< cloud water to rain auto conversion scheme

  logical :: de_ice = .false. !< to prevent excessive build - up of cloud ice from external sources
  logical :: sedi_transport = .false. !< transport of momentum in sedimentation
  logical :: do_sedi_w = .false. !< transport of vertical motion in sedimentation
  logical :: do_sedi_heat = .false. !< transport of heat in sedimentation
  logical :: prog_ccn = .false. !< do prognostic ccn (yi ming's method)
  logical :: do_bigg = .false. !< do bigg mechanism freezing of supercooled liquid on aerosol nuclei
  logical :: do_evap = .false. !< do evaporation
  logical :: do_subl = .false. !< do sublimation
  logical :: do_qa = .false. !< do inline cloud fraction (WMP: in FV3 dynamics)
  logical :: preciprad = .true. !< consider precipitates in cloud fraciton calculation
  logical :: fix_negative = .false. !< fix negative water species
  logical :: do_setup = .true. !< setup constants and parameters
  logical :: p_nonhydro = .false. !< perform hydrosatic adjustment on air density

  real, dimension(TABLE_LENGTH) :: table, table2, table3, tablew
  real, dimension(TABLE_LENGTH) :: des, des2, des3, desw

  logical :: tables_are_initialized = .false.

  ! logical :: root_proc
  ! integer :: id_rh, id_vtr, id_vts, id_vtg, id_vti, id_rain, id_snow, id_graupel, &
  ! id_ice, id_prec, id_cond, id_var, id_droplets
  ! integer :: gfdl_mp_clock ! clock for timing of driver routine

  real :: dt_fr = 8. !< epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
  ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
  ! dt_fr can be considered as the error bar

  real, parameter :: p_min = 100. !< minimum pressure (pascal) for mp to operate

  ! slj, the following parameters are for cloud - resolving resolution: 1 - 5 km

  ! qi0_crt = 0.8e-4
  ! qs0_crt = 0.6e-3
  ! c_psaci = 0.1
  ! c_pgacs = 0.1
  ! c_pgaci = 0.05

  ! -----------------------------------------------------------------------
  !> namelist parameters
  ! -----------------------------------------------------------------------

  real :: cld_min = 0.05 !< minimum cloud fraction
  real :: tice = 273.16 !< set tice = 165. to trun off ice - phase phys (kessler emulator)

  real :: log_10 = log (10.)
  real :: tice0 = 273.16 - 0.01
  real, parameter :: t_wfr = 273.16 - 40.0 ! supercooled water can exist down to - 40 c, which is the "absolute"

  real :: t_min = 178. !< min temp to freeze - dry all water vapor
  real :: t_sub = 184. !< min temp for sublimation of cloud ice
  real :: mp_time = 150. !< maximum micro - physics time step (sec)

  ! relative humidity increment

  real :: rh_inc = 0.25 !< rh increment for complete evaporation of cloud water and cloud ice
  real :: rh_inr = 0.25 !< rh increment for minimum evaporation of rain
  real :: rh_ins = 0.25 !< rh increment for sublimation of snow

  ! conversion time scale

  real :: tau_r2g = 900. !< rain freezing during fast_sat
  real :: tau_smlt = 900. !< snow melting
  real :: tau_g2r = 600. !< graupel melting to rain
  real :: tau_imlt = 600. !< cloud ice melting
  real :: tau_i2s = 1000. !< cloud ice to snow auto - conversion
  real :: tau_l2r = 900. !< cloud water to rain auto - conversion
  real :: tau_v2l = 150. !< water vapor to cloud water (condensation)
  real :: tau_l2v = 300. !< cloud water to water vapor (evaporation)
  real :: tau_i2v = 300. !< cloud ice to water vapor (sublimation)
  real :: tau_s2v = 600. !< snow sublimation
  real :: tau_v2s = 21600. !< snow deposition -- make it a slow process
  real :: tau_g2v = 900. !< graupel sublimation
  real :: tau_v2g = 21600. !< graupel deposition -- make it a slow process
  real :: tau_revp = 600. !< rain re-evaporation
  real :: tau_frz = 450. !, timescale for liquid-ice freezing
  ! horizontal subgrid variability

  real :: dw_land = 0.20 !< base value for subgrid deviation / variability over land
  real :: dw_ocean = 0.10 !< base value for ocean

  ! prescribed ccn

  real :: ccn_o = 90. !< ccn over ocean (cm^ - 3)
  real :: ccn_l = 270. !< ccn over land (cm^ - 3)

  real :: rthreshu =  7.0e-6 !< critical cloud drop radius (micro m)
  real :: rthreshs = 10.0e-6 !< critical cloud drop radius (micro m)

  real :: sat_adj0 = 0.90 !< adjustment factor (0: no, 1: full) during fast_sat_adj

  real :: qc_crt = 5.0e-8 !< mini condensate mixing ratio to allow partial cloudiness

  real :: qi_lim = 1. !< cloud ice limiter to prevent large ice build up

  real :: ql_mlt = 2.0e-3 !< max value of cloud water allowed from melted cloud ice
  real :: qs_mlt = 1.0e-6 !< max cloud water due to snow melt

  real :: ql_gen = 1.0e-3 !< max cloud water generation
  real :: qi_gen = 9.82679e-5 !< max cloud ice generation at -40 C

  ! cloud condensate upper bounds: "safety valves" for ql & qi

  real :: ql0_max = 2.0e-3 !< max cloud water value (auto converted to rain)
  real :: qi0_max = 1.0e-4 !< max cloud ice value (by other sources)

  real :: qi0_crt = 1.0e-4 !< cloud ice to snow autoconversion threshold (was 1.e-4)
  !! qi0_crt is highly dependent on horizontal resolution
  real :: qr0_crt = 1.0e-4 !< rain to snow or graupel / hail threshold
  !! lfo used * mixing ratio * = 1.e-4 (hail in lfo)
  real :: qs0_crt = 1.0e-3 !< snow to graupel density threshold (0.6e-3 in purdue lin scheme)

  real :: c_paut = 0.55 !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion)
  real :: c_psaci = 0.02 !< accretion: cloud ice to snow (was 0.1 in zetac)
  real :: c_piacr = 5.0 !< accretion: rain to ice:
  real :: c_cracw = 0.9 !< rain accretion efficiency
  real :: c_pgacs = 2.0e-3 !< snow to graupel "accretion" eff. (was 0.1 in zetac)
  real :: c_pgaci = 0.05   !<  ice to graupel "accretion" eff.

  ! decreasing clin to reduce csacw (so as to reduce cloud water --- > snow)

  real :: alin = 842.0 !< "a" in lin1983
  real :: clin = 4.8 !< "c" in lin 1983, 4.8 -- > 6. (to ehance ql -- > qs)

  ! fall velocity tuning constants:

  logical :: const_vi = .false. !< if .t. the constants are specified by v * _fac
  logical :: const_vs = .false. !< if .t. the constants are specified by v * _fac
  logical :: const_vg = .false. !< if .t. the constants are specified by v * _fac
  logical :: const_vr = .false. !< if .t. the constants are specified by v * _fac

  ! good values:

  real :: vi_fac = 1. !< if const_vi: 1 / 3
  real :: vs_fac = 1. !< if const_vs: 1.
  real :: vg_fac = 1. !< if const_vg: 2.
  real :: vr_fac = 1. !< if const_vr: 4.

  ! upper bounds of fall speed (with variable speed option)

  real :: vi_max = 1.0 !< max fall speed for ice
  real :: vs_max = 2.0 !< max fall speed for snow
  real :: vg_max = 12. !< max fall speed for graupel
  real :: vr_max = 12. !< max fall speed for rain

  ! cloud microphysics switchers

  logical :: fast_sat_adj = .false. !< has fast saturation adjustments
  logical :: z_slope_liq = .true. !< use linear mono slope for autocconversions
  logical :: z_slope_ice = .false. !< use linear mono slope for autocconversions
  logical :: use_ccn = .false. !< use input ccn when .T. else use ccn_o/ccn_l
  logical :: use_ppm = .false. !< use ppm fall scheme
  logical :: mono_prof = .true. !< perform terminal fall with mono ppm scheme
  logical :: mp_print = .false. !< cloud microphysics debugging printout

  ! real :: global_area = - 1.

  ! -----------------------------------------------------------------------
  ! namelist
  ! -----------------------------------------------------------------------

  namelist / gfdl_cloud_microphysics_nml /                                  &
       mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
       vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max,  &
       vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max,    &
       qi0_crt, qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
       const_vs, const_vg, const_vr, use_ccn, rthreshu, rthreshs, ccn_l, ccn_o, qc_crt, &
       tau_g2v, tau_v2g, tau_s2v, tau_v2s, &
       tau_revp, tau_frz, do_bigg, do_evap, do_subl, &
       sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, &
       tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, c_pgaci,  &
       z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin,              &
       preciprad, cld_min, use_ppm, mono_prof,         &
       do_sedi_heat, sedi_transport, do_sedi_w, dt_fr, de_ice, icloud_f, irain_f, mp_print

  public                                                                    &
       mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, dw_land, dw_ocean, &
       vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max,  &
       vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max,    &
       qi0_crt, qr0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
       const_vs, const_vg, const_vr, use_ccn, rthreshu, rthreshs, ccn_l, ccn_o, qc_crt, &
       tau_g2v, tau_v2g, tau_s2v, tau_v2s, &
       tau_revp, tau_frz, do_bigg, do_evap, do_subl, &
       sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, &
       tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, c_pgaci,  &
       z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin,              &
       preciprad, cld_min, use_ppm, mono_prof,         &
       do_sedi_heat, sedi_transport, do_sedi_w, dt_fr, de_ice, icloud_f, irain_f, mp_print

  !$omp declare target( &
  !$omp   des2, desw, table2, tablew, &

  !$omp   d0_vap, lv00, c_vap, c_air, tau_revp, &
  !$omp   tau_v2l, tau_l2v, tau_i2v, tau_s2v, tau_v2s, tau_g2v, &
  !$omp   tau_v2g, tau_frz, tau_imlt, tau_smlt, tau_i2s, tau_g2r, &
  !$omp   tice, tice0, rh_inc, rh_inr, t_min, do_qa, t_sub, do_evap, &
  !$omp   do_bigg, qi_lim, do_subl, preciprad, icloud_f, qc_crt, lat2, z_slope_ice, &
  !$omp   c_paut, prog_ccn, fix_negative, p_nonhydro, sedi_transport, ql_mlt, qs_mlt, qi0_crt, qs0_crt, &
  !$omp   const_vi, vi_fac, vi_max, const_vs, vs_fac, vs_max, const_vg, vg_fac, vg_max, const_vr, vr_fac, vr_max, &
  !$omp   do_sedi_w, use_ppm, mono_prof, rthreshs, rthreshu, irain_f, z_slope_liq, do_sedi_heat, &
  !$omp   ql0_max, dt_fr, sat_adj0, dw_land, dw_ocean, c_psaci, c_pgacs, &
  !$omp   ccn_l, ccn_o, c_cracw, use_ccn, de_ice, mp_time, &

  !$omp   ces0, cracs, cracw, &
  !$omp   csaci, csacr, csacw, cgaci, cgacr, cgacs, cgacw, &
  !$omp   cssub, crevp, csmlt, cgmlt, cgfr, acco)

contains

  ! -----------------------------------------------------------------------
  ! the driver of the gfdl cloud microphysics
  ! -----------------------------------------------------------------------

  !>@brief The subroutine 'gfdl_cloud_microphys_driver' executes the full GFDL
  !! cloud microphysics.
  subroutine gfdl_cloud_microphys_driver ( &
       qv, ql, qr, qi, qs, qg, qa, qn, &
       qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, &
       pt_dt, pt, w, &
       uin, vin, udt, vdt, &
       dz, delp, area, dt_in, &
       land, cnv_fraction, srf_type, eis,                                &
       rhcrit, anv_icefall, lsc_icefall,                                 &
       revap, isubl,                                                     &
       rain, snow, ice, graupel, &
       m2_rain, m2_sol, &
       hydrostatic, phys_hydrostatic, &
       iis, iie, jjs, jje, kks, kke, ktop, kbot)

    implicit none

    logical, intent (in) :: hydrostatic, phys_hydrostatic
    integer, intent (in) :: iis, iie, jjs, jje !< physics window
    integer, intent (in) :: kks, kke !< vertical dimension
    integer, intent (in) :: ktop, kbot !< vertical compute domain

    real, intent (in) :: dt_in !< physics time step

    real, intent (in), dimension (:, :) :: area !< cell area
    real, intent (in), dimension (:, :) :: land !< land fraction
    real, intent (in), dimension (:, :) :: cnv_fraction !< diagnosed convective fraction
    real, intent (in), dimension (:, :) :: srf_type
    real, intent (in), dimension (:, :) :: eis  !< estimated inversion strength
    real, intent (in), dimension (:, :, :) :: rhcrit

    real, intent (in) :: anv_icefall, lsc_icefall

    real, intent (in), dimension (:, :, :) :: delp, dz, uin, vin
    real, intent (in), dimension (:, :, :) :: pt, qv, ql, qr, qg, qa, qn

    real, intent (inout), dimension (:, :, :) :: qi, qs
    real, intent (inout), dimension (:, :, :) :: pt_dt, qa_dt, udt, vdt, w
    real, intent (inout), dimension (:, :, :) :: qv_dt, ql_dt, qr_dt
    real, intent (inout), dimension (:, :, :) :: qi_dt, qs_dt, qg_dt

    real, intent (out), dimension (:, :) :: rain, snow, ice, graupel
    real, intent (out), dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
    real, intent (out), dimension (:, :, :) :: revap ! Rain evaporation
    real, intent (out), dimension (:, :, :) :: isubl ! Ice sublimation

    ! logical :: used

    real :: mpdt, rdt, dts, convt, tot_prec

    integer :: i, j, k
    integer :: is, ie, js, je !< physics window
    integer :: ks, ke !< vertical dimension
    integer :: days, ntimes

    real, dimension (iie - iis + 1, jje - jjs + 1) :: prec_mp, prec1, cond, w_var, rh0

    real, dimension (iie - iis + 1, jje - jjs + 1, kke - kks + 1) :: vt_r, vt_s, vt_g, vt_i, qn2

    real :: allmax

    is = 1
    js = 1
    ks = 1
    ie = iie - iis + 1
    je = jje - jjs + 1
    ke = kke - kks + 1

    ! call mpp_clock_begin (gfdl_mp_clock)

    ! -----------------------------------------------------------------------
    ! define heat capacity of dry air and water vapor based on hydrostatical property
    ! -----------------------------------------------------------------------

    if (phys_hydrostatic .or. hydrostatic) then
       c_air = cp_air
       c_vap = cp_vap
       p_nonhydro = .false.
    else
       c_air = cv_air
       c_vap = cv_vap
       p_nonhydro = .true.
    endif

    d0_vap = c_vap - c_liq
    lv00 = hlv0 - d0_vap * t_ice

    if (hydrostatic) do_sedi_w = .false.

    ! -----------------------------------------------------------------------
    ! define latent heat coefficient used in wet bulb and bigg mechanism
    ! -----------------------------------------------------------------------

    latv = hlv
    lati = hlf
    lats = latv + lati
    lat2 = lats * lats

    !$omp target update to(c_air, c_vap, p_nonhydro, d0_vap, lv00, do_sedi_w, lat2)

    lcp = latv / cp_air
    icp = lati / cp_air
    tcp = (latv + lati) / cp_air

    ! tendency zero out for am moist processes should be done outside the driver

    ! -----------------------------------------------------------------------
    ! define cloud microphysics sub time step
    ! -----------------------------------------------------------------------

    mpdt   = min (dt_in, mp_time)
    rdt    = 1. / dt_in
    ntimes = nint (dt_in / mpdt)

    ! small time step:
    dts = dt_in / real (ntimes)

    ! call get_time (time, seconds, days)

    ! -----------------------------------------------------------------------
    ! initialize precipitation
    ! -----------------------------------------------------------------------

    do j = js, je
       do i = is, ie
          graupel (i, j) = 0.
          rain (i, j) = 0.
          snow (i, j) = 0.
          ice (i, j) = 0.
          cond (i, j) = 0.
       enddo
    enddo

    ! -----------------------------------------------------------------------
    ! major cloud microphysics
    ! -----------------------------------------------------------------------

    ! print *, 'gfdl_cloud_microphys_driver - calling mpdrv'
    call mpdrv ( &
         hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, qg, &
         qa, qn, dz, is, ie, js, je, ks, ke, ktop, kbot, dt_in, ntimes, &
         rain(:, js:je), snow(:, js:je), graupel(:, js:je), ice(:, js:je), m2_rain, &
         m2_sol, cond(:, js:je), area(:, js:je), &
         land(:, js:je), cnv_fraction(:, js:je), srf_type(:, js:je), eis(:, js:je), &
         rhcrit, anv_icefall, lsc_icefall, &
         revap, isubl, &
         udt, vdt, pt_dt, &
         qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, w_var, vt_r, &
         vt_s, vt_g, vt_i, qn2)
    ! print *, 'gfdl_cloud_microphys_driver - completed mpdrv'

    ! -----------------------------------------------------------------------
    ! no clouds allowed above ktop
    ! -----------------------------------------------------------------------

    if (ks < ktop) then
       do k = ks, ktop
          do j = js, je
             do i = is, ie
                qa_dt (i, j, k) = 0.
             enddo
          enddo
       enddo
    endif

    ! -----------------------------------------------------------------------
    ! diagnostic output
    ! -----------------------------------------------------------------------

    ! if (id_vtr > 0) then
    ! used = send_data (id_vtr, vt_r, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_vts > 0) then
    ! used = send_data (id_vts, vt_s, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_vtg > 0) then
    ! used = send_data (id_vtg, vt_g, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_vti > 0) then
    ! used = send_data (id_vti, vt_i, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_droplets > 0) then
    ! used = send_data (id_droplets, qn2, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_var > 0) then
    ! used = send_data (id_var, w_var, time, is_in = iis, js_in = jjs)
    ! endif

    ! convert to mm / day

    convt = 86400. * rdt * rgrav
    do j = js, je
       do i = is, ie
          rain (i, j) = rain (i, j) * convt
          snow (i, j) = snow (i, j) * convt
          ice (i, j) = ice (i, j) * convt
          graupel (i, j) = graupel (i, j) * convt
          prec_mp (i, j) = rain (i, j) + snow (i, j) + ice (i, j) + graupel (i, j)
       enddo
    enddo

    ! if (id_cond > 0) then
    ! do j = js, je
    ! do i = is, ie
    ! cond (i, j) = cond (i, j) * rgrav
    ! enddo
    ! enddo
    ! used = send_data (id_cond, cond, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (id_snow > 0) then
    ! used = send_data (id_snow, snow, time, iis, jjs)
    ! used = send_data (id_snow, snow, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (snow, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean snow = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_graupel > 0) then
    ! used = send_data (id_graupel, graupel, time, iis, jjs)
    ! used = send_data (id_graupel, graupel, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (graupel, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean graupel = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_ice > 0) then
    ! used = send_data (id_ice, ice, time, iis, jjs)
    ! used = send_data (id_ice, ice, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (ice, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean ice_mp = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_rain > 0) then
    ! used = send_data (id_rain, rain, time, iis, jjs)
    ! used = send_data (id_rain, rain, time, is_in = iis, js_in = jjs)
    ! if (mp_print .and. seconds == 0) then
    ! tot_prec = g_sum (rain, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'mean rain = ', tot_prec
    ! endif
    ! endif
    !
    ! if (id_rh > 0) then !not used?
    ! used = send_data (id_rh, rh0, time, iis, jjs)
    ! used = send_data (id_rh, rh0, time, is_in = iis, js_in = jjs)
    ! endif
    !
    !
    ! if (id_prec > 0) then
    ! used = send_data (id_prec, prec_mp, time, iis, jjs)
    ! used = send_data (id_prec, prec_mp, time, is_in = iis, js_in = jjs)
    ! endif

    ! if (mp_print) then
    ! prec1 (:, :) = prec1 (:, :) + prec_mp (:, :)
    ! if (seconds == 0) then
    ! prec1 (:, :) = prec1 (:, :) * dt_in / 86400.
    ! tot_prec = g_sum (prec1, is, ie, js, je, area, 1)
    ! if (root_proc) write (*, *) 'daily prec_mp = ', tot_prec
    ! prec1 (:, :) = 0.
    ! endif
    ! endif

    ! call mpp_clock_end (gfdl_mp_clock)

  end subroutine gfdl_cloud_microphys_driver

  ! -----------------------------------------------------------------------
  !>@brief gfdl cloud microphysics, major program
  !>@details lin et al., 1983, jam, 1065 - 1092, and
  !! rutledge and hobbs, 1984, jas, 2949 - 2972
  !! terminal fall is handled lagrangianly by conservative fv algorithm
  !>@param pt: temperature (k)
  !>@param 6 water species:
  !>@param 1) qv: water vapor (kg / kg)
  !>@param 2) ql: cloud water (kg / kg)
  !>@param 3) qr: rain (kg / kg)
  !>@param 4) qi: cloud ice (kg / kg)
  !>@param 5) qs: snow (kg / kg)
  !>@param 6) qg: graupel (kg / kg)
  ! -----------------------------------------------------------------------
  subroutine mpdrv (hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, &
       qg, qa, qn, dz, is, ie, js, je, ks, ke, ktop, kbot, dt_in, ntimes, &
       rain, snow, graupel, ice, m2_rain, m2_sol, cond, area1, land, &
       cnv_fraction, srf_type, eis, rhcrit, anv_icefall, lsc_icefall, revap, isubl, &
       u_dt, v_dt, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, &
       w_var, vt_r, vt_s, vt_g, vt_i, qn2)

    implicit none

    logical, intent (in) :: hydrostatic

    integer, intent (in) :: is, ie, js, je, ks, ke
    integer, intent (in) :: ntimes, ktop, kbot

    real, intent (in) :: dt_in

    real, intent (in), dimension (is:, js:) :: area1, land
    real, intent (in), dimension (is:, js:) :: cnv_fraction
    real, intent (in), dimension (is:, js:) :: srf_type
    real, intent (in), dimension (is:, js:) :: eis

    real, intent (in), dimension (is:, js:, ks:) :: rhcrit

    real, intent (in) :: anv_icefall, lsc_icefall

    real, intent (in), dimension (is:, js:, ks:) :: uin, vin, delp, pt, dz
    real, intent (in), dimension (is:, js:, ks:) :: qv, qi, ql, qr, qs, qg, qa, qn

    real, intent (inout), dimension (is:, js:, ks:) :: u_dt, v_dt, w, pt_dt, qa_dt
    real, intent (inout), dimension (is:, js:, ks:) :: qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
    real, intent (  out), dimension (is:, js:, ks:) :: revap, isubl

    real, intent (inout), dimension (is:, js:) :: rain, snow, ice, graupel, cond

    real, intent (out), dimension (is:, js:) :: w_var

    real, intent (out), dimension (is:, js:, ks:) :: vt_r, vt_s, vt_g, vt_i, qn2

    real, intent (out), dimension (is:, js:, ks:) :: m2_rain, m2_sol

    ! real, dimension (ktop:kbot) :: h_var1d
    ! real, dimension (ktop:kbot) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    ! real, dimension (ktop:kbot) :: vtiz, vtsz, vtgz, vtrz
    ! real, dimension (ktop:kbot) :: dp1, dz1
    ! real, dimension (ktop:kbot) :: qv0, ql0, qr0, qi0, qs0, qg0
    ! real, dimension (ktop:kbot) :: den, den0, tz, p1, denfac
    ! real, dimension (ktop:kbot) :: ccn, c_praut, m1_rain, m1_sol, m1, evap1, subl1
    ! real, dimension (ktop:kbot) :: w1, tmp1, tmp2

    real, dimension (is:ie, js:je, ktop:kbot) :: h_var1d
    real, dimension (is:ie, js:je, ktop:kbot) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (is:ie, js:je, ktop:kbot) :: vtiz, vtsz, vtgz, vtrz
    real, dimension (is:ie, js:je, ktop:kbot) :: dp1, dz1
    real, dimension (is:ie, js:je, ktop:kbot) :: qv0, ql0, qr0, qi0, qs0, qg0
    real, dimension (is:ie, js:je, ktop:kbot) :: den, tz, p1, denfac
    real, dimension (is:ie, js:je, ktop:kbot) :: ccn, c_praut, m1_rain, m1_sol, m1, evap1, subl1
    real, dimension (is:ie, js:je, ktop:kbot) :: w1

    real :: cpaut, t0, dts, rdt, den0
    real, dimension(is:ie, js:je) :: r1, s1, i1, g1
    real :: cvm, omq
    real :: u1_k, u1_km1, v1_k, v1_km1

    integer :: i, j, k, n

    integer :: num_devices, nteams, nthreads
    logical :: initial_device

    ! num_devices = omp_get_num_devices()
    ! print *, "Number of available devices", num_devices

    ! !$omp target map(initial_device, nteams, nthreads)
    ! initial_device = omp_is_initial_device()
    ! nteams = omp_get_num_teams()
    ! nthreads = omp_get_num_threads()
    ! !$omp end target
    ! if (initial_device) then
    !    print *, "Running on host"
    ! else
    !    print *, "Running on device with ", nteams, " teams in total and ", nthreads, " threads per team"
    ! end if

    ! -----------------------------------------------------------------------
    ! calculate cloud condensation nuclei (ccn)
    ! the following is based on klein eq. 15
    ! -----------------------------------------------------------------------

    ! c_air = 0
    ! print *, 'c_air (host): ', c_air
    ! !$omp target update from(c_air)
    ! print *, 'c_air (device): ', c_air

    ! print *, 'pt: ', minval(pt), maxval(pt), sum(pt)
    ! print *, 'delp: ', minval(delp), maxval(delp), sum(delp)
    ! print *, 'qv: ', minval(qv), maxval(qv), sum(qv)
    ! print *, 'ql: ', minval(ql), maxval(ql), sum(ql)
    ! print *, 'qi: ', minval(qi), maxval(qi), sum(qi)
    ! print *, 'qr: ', minval(qr), maxval(qr), sum(qr)
    ! print *, 'qs: ', minval(qs), maxval(qs), sum(qs)
    ! print *, 'qg: ', minval(qg), maxval(qg), sum(qg)
    ! print *, 'qn: ', minval(qn), maxval(qn), sum(qn)
    ! print *, 'w min/max/sum: ', minval(w), maxval(w), sum(w)

    dts = dt_in / real (ntimes)
    rdt = 1. / dt_in

    ! -----------------------------------------------------------------------
    ! calculate cloud condensation nuclei (ccn)
    ! the following is based on klein eq. 15
    ! -----------------------------------------------------------------------

    !$omp target data &
    !$omp   map(to: &
    !$omp     is, ie, js, je, ktop, kbot, ntimes, dt_in, dts, rdt, &
    !$omp     area1, land, cnv_fraction, srf_type, eis, &
    !$omp     rhcrit, anv_icefall, lsc_icefall, &
    !$omp     uin, vin, delp, pt, dz, &
    !$omp     qv, qi, ql, qr, qs, qg, qa, qn) &

    !$omp   map(tofrom: &
    !$omp     u_dt, v_dt, w, pt_dt, qa_dt, &
    !$omp     qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, &
    !$omp     rain, snow, ice, graupel, cond) &

    !$omp   map(from: &
    !$omp     revap, isubl, w_var, &
    !$omp     vt_r, vt_s, vt_g, vt_i, qn2, m2_rain, m2_sol) &

    !$omp   map(alloc: &
    !$omp     tz, dp1, h_var1d, &
    !$omp     qvz, qlz, qiz, qrz, qsz, qgz, &
    !$omp     qv0, ql0, qi0, qr0, qs0, qg0, &
    !$omp     m1, w1, p1, ccn, c_praut, dz1, den, denfac, &
    !$omp     vtrz, vtsz, vtiz, vtgz, m1_sol)

    !$omp target teams distribute parallel do simd collapse(3) private(t0)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             ! Initialize
             m2_rain (i, j, k) = 0.
             m2_sol (i, j, k) = 0.
             revap (i, j, k) = 0.
             isubl (i, j, k) = 0.

             t0 = pt (i, j, k)
             tz (i, j, k) = t0
             dp1 (i, j, k) = delp (i, j, k)
             ! dp0 (k) = dp1 (k) ! moist air mass * grav

             ! import horizontal subgrid variability with pressure dependence
             ! total water subgrid deviation in horizontal direction
             ! default area dependent form: use dx ~ 100 km as the base
             h_var1d (i, j, k) = min(0.30,1.0 - rhcrit(i,j,k)) ! restricted to 70%
#ifdef __GFORTRAN__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! Convert moist mixing ratios to dry mixing ratios

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
#endif
             qvz (i, j, k) = qv (i, j, k)
             qlz (i, j, k) = ql (i, j, k)
             qiz (i, j, k) = qi (i, j, k)
             qrz (i, j, k) = qr (i, j, k)
             qsz (i, j, k) = qs (i, j, k)
             qgz (i, j, k) = qg (i, j, k)
#ifdef __GFORTRAN__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
#endif
             qaz (i, j, k) = qa (i, j, k)
#ifdef __GFORTRAN__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(omq, den0)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
#endif
             ! dp1: dry air_mass
             ! dp1 (k) = dp1 (k) * (1. - (qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)))
             dp1 (i, j, k) = dp1 (i, j, k) * (1. - qvz (i, j, k)) ! gfs
             omq = delp (i, j, k) / dp1 (i, j, k)
             qvz (i, j, k) = qvz (i, j, k) * omq
             qlz (i, j, k) = qlz (i, j, k) * omq
             qrz (i, j, k) = qrz (i, j, k) * omq
             qiz (i, j, k) = qiz (i, j, k) * omq
             qsz (i, j, k) = qsz (i, j, k) * omq
             qgz (i, j, k) = qgz (i, j, k) * omq

             den0 = - dp1 (i, j, k) / (grav * dz (i, j, k)) ! density of dry air
             p1 (i, j, k) = den0 * rdgas * t0 ! dry air pressure

             ! -----------------------------------------------------------------------
             ! save a copy of old value for computing tendencies
             ! -----------------------------------------------------------------------

             qv0 (i, j, k) = qvz (i, j, k)
             ql0 (i, j, k) = qlz (i, j, k)
             qr0 (i, j, k) = qrz (i, j, k)
             qi0 (i, j, k) = qiz (i, j, k)
             qs0 (i, j, k) = qsz (i, j, k)
             qg0 (i, j, k) = qgz (i, j, k)
#ifdef __GFORTRAN__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! for sedi_momentum

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(cpaut)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
#endif
             m1 (i, j, k) = 0.
             if (do_sedi_w) w1 (i, j, k) = w (i, j, k)
#ifdef __GFORTRAN__
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(3) default(shared) private(cpaut)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
#endif
             cpaut = c_paut * 0.104 * grav / 1.717e-5
             ! ccn needs units #/m^3
             if (prog_ccn) then
                ! qn has units # / m^3
                ccn (i, j, k) = qn (i, j, k)
                c_praut (i, j, k) = cpaut * (ccn (i, j, k) * rhor) ** (- 1. / 3.)
             else
                ! qn has units # / m^3
                ccn (i, j, k) = qn (i, j, k)
                !!! use GEOS ccn: ccn (i, j, k) = (ccn_l * land (i) + ccn_o * (1. - land (i))) * 1.e6
                c_praut (i, j, k) = cpaut * (ccn (i, j, k) * rhor) ** (- 1. / 3.)
             endif
             m2_rain (i, j, k) = c_praut (i,j, k)
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! fix all negative water species
    ! -----------------------------------------------------------------------

    if (fix_negative) then
       call neg_adj (ie, ie, js, je, ktop, kbot, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)
    endif

    do n = 1, ntimes ! TODO: do n = 1, ntimes

       ! dry air density

       !$omp target teams distribute parallel do simd collapse(3) private(t0)
       do k = ktop, kbot
          do j = js, je
             do i = is, ie
                if (p_nonhydro) then
                   dz1 (i, j, k) = dz (i, j, k)
                   den(i,j,k) = - dp1 (i,j,k) / (grav * dz (i, j, k)) ! density of dry air
                   denfac (i, j, k) = sqrt (sfcrho / den (i, j, k))
                else
                   t0 = pt (i, j, k)
                   dz1 (i, j, k) = dz (i, j, k) * tz (i, j, k) / t0 ! hydrostatic balance
                   den(i,j,k) = (- dp1 (i,j,k) / (grav * dz (i, j, k))) * dz(i,j,k) / dz1(i,j,k) ! density of dry air
                   denfac (i, j, k) = sqrt (sfcrho / den (i, j, k))
                end if
#ifdef __GFORTRAN__
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

       !$omp target teams distribute parallel do simd collapse(3) private(t0)
       do k = ktop, kbot
          do j = js, je
             do i = is, ie
#endif
                ! sedimentation of cloud ice, snow, and graupel
                call fall_speed ( &
                     ktop, kbot, p1(i, j, k), cnv_fraction(i, j), anv_icefall, lsc_icefall, &
                     den(i, j, k), qsz(i, j, k), qiz(i, j, k), qgz(i, j, k), qlz(i, j, k), tz(i, j, k), &
                     vtsz(i, j, k), vtiz(i, j, k), vtgz(i, j, k)) ! output
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

       call terminal_fall_3d ( &
            dts, is, ie, js, je, ktop, kbot, &
            tz, qvz, qlz, qrz, qgz, qsz, qiz, &
            dz1, dp1, den, vtgz, vtsz, vtiz, &
            r1, g1, s1, i1, & ! output
            m1_sol, w1)

       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             rain (i, j) = rain (i, j) + r1 (i, j) ! from melted snow & ice that reached the ground
             snow (i, j) = snow (i, j) + s1 (i, j)
             graupel (i, j) = graupel (i, j) + g1 (i, j)
             ice (i, j) = ice (i, j) + i1 (i, j)
          end do
       end do

       ! heat transportation during sedimentation

       if (do_sedi_heat) then
          call sedi_heat_3d (is, ie, js, je, ktop, kbot, dp1, m1_sol, dz1, tz, qvz, qlz, qrz, qiz, qsz, qgz, c_ice)
       endif

       ! warm rain processes

       call warm_rain_3d ( &
            dts, is, ie, js, je, ktop, kbot, &
            dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
            eis, den, denfac, ccn, c_praut, vtrz, &
            r1, evap1, m1_rain, w1, h_var1d)

       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             rain (i, j) = rain (i, j) + r1 (i, j)
          end do
       end do
       !$omp end target teams distribute parallel do simd

       !$omp target teams distribute parallel do simd collapse(3)
       do k = ktop, kbot
          do i = is, ie
             do j = js, je
                revap (i,j,k) = revap (i,j,k) + evap1(i, j, k)
                m2_rain (i, j, k) = m2_rain (i, j, k) + m1_rain (i, j, k)
                m2_sol (i, j, k) = m2_sol (i, j, k) + m1_sol (i, j, k)
                m1 (i, j, k) = m1 (i, j, k) + m1_rain (i, j, k) + m1_sol (i, j, k)
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

       ! -----------------------------------------------------------------------
       ! ice - phase microphysics
       ! -----------------------------------------------------------------------

       call icloud_3d ( &
            is, ie, js, je, ktop, kbot, tz, p1, &
            qvz, qlz, qrz, qiz, qsz, qgz, dp1, den, denfac, &
            vtsz, vtgz, vtrz, qaz, dts, subl1, h_var1d, ccn, cnv_fraction, srf_type)

       !$omp target teams distribute parallel do simd collapse(3)
       do k = ktop, kbot
          do i = is, ie
             do j = js, je
                isubl (i,j,k) = isubl (i,j,k) + subl1(i, j, k)
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

    end do ! ntimes

    !$omp end target data

    ! print *, 'm2_rain: ', minval(m2_rain), maxval(m2_rain), sum(m2_rain)
    ! print *, 'm2_sol: ', minval(m2_sol), maxval(m2_sol), sum(m2_sol)
    ! print *, 'rain: ', minval(rain), maxval(rain), sum(rain)
    ! print *, 'snow: ', minval(snow), maxval(snow), sum(snow)
    ! print *, 'graupel: ', minval(graupel), maxval(graupel), sum(graupel)
    ! print *, 'ice: ', minval(ice), maxval(ice), sum(ice)

    ! do j = js, je

    !    do i = is, ie

    !       ! -----------------------------------------------------------------------
    !       ! momentum transportation during sedimentation
    !       ! note: dp1 is dry mass; delp(i, j, :) is the old moist (total) mass
    !       ! -----------------------------------------------------------------------

    !       if (sedi_transport) then
    !          v1_km1 = vin (i, j, ktop)
    !          u1_km1 = uin (i, j, ktop)
    !          !$acc loop seq
    !          do k = ktop + 1, kbot
    !             u1_k = uin (i, j, k)
    !             v1_k = vin (i, j, k)
    !             u1_k = (delp (i, j, k) * u1_k + m1 (k - 1) * u1_km1) / (delp (i, j, k) + m1 (k - 1))
    !             v1_k = (delp (i, j, k) * v1_k + m1 (k - 1) * v1_km1) / (delp (i, j, k) + m1 (k - 1))
    !             u_dt (i, j, k) = u_dt (i, j, k) + (u1_k - uin (i, j, k)) * rdt
    !             v_dt (i, j, k) = v_dt (i, j, k) + (v1_k - vin (i, j, k)) * rdt
    !             u1_km1 = u1_k ! store for next iteration
    !             v1_km1 = v1_k
    !          enddo
    !       endif

    !       if (do_sedi_w) then
    !          !$acc loop vector
    !          do k = ktop, kbot
    !             w (i, j, k) = w1 (k)
    !          enddo
    !       endif

    !       ! -----------------------------------------------------------------------
    !       ! update moist air mass (actually hydrostatic pressure)
    !       ! convert to dry mixing ratios
    !       ! -----------------------------------------------------------------------

    !       !$acc loop vector private(omq, cvm, t0)
    !       do k = ktop, kbot
    !          t0 = pt (i, j, k)
    !          omq = dp1 (k) / delp (i, j, k)
    !          qv_dt (i, j, k) = qv_dt (i, j, k) + rdt * (qvz (k) - qv0 (k)) * omq
    !          ql_dt (i, j, k) = ql_dt (i, j, k) + rdt * (qlz (k) - ql0 (k)) * omq
    !          qr_dt (i, j, k) = qr_dt (i, j, k) + rdt * (qrz (k) - qr0 (k)) * omq
    !          qi_dt (i, j, k) = qi_dt (i, j, k) + rdt * (qiz (k) - qi0 (k)) * omq
    !          qs_dt (i, j, k) = qs_dt (i, j, k) + rdt * (qsz (k) - qs0 (k)) * omq
    !          qg_dt (i, j, k) = qg_dt (i, j, k) + rdt * (qgz (k) - qg0 (k)) * omq
    !          cvm = c_air + qvz (k) * c_vap + (qrz (k) + qlz (k)) * c_liq + (qiz (k) + qsz (k) + qgz (k)) * c_ice
    !          pt_dt (i, j, k) = pt_dt (i, j, k) + rdt * (tz (k) - t0) * cvm / cp_air
    !       enddo

    !       ! -----------------------------------------------------------------------
    !       ! update cloud fraction tendency
    !       ! -----------------------------------------------------------------------
    !       if (.not. do_qa) then
    !          !$acc loop vector
    !          do k = ktop, kbot
    !             qa_dt (i, j, k) = qa_dt (i, j, k) + rdt * (                           &
    !                  qa(i, j, k)*SQRT( (qiz(k)+qlz(k)) / max(qi0(k)+ql0(k),qcmin) ) - & ! New Cloud -
    !                  qa(i, j, k) )                                                      ! Old Cloud
    !          enddo
    !       endif

    !       ! -----------------------------------------------------------------------
    !       ! fms diagnostics:
    !       ! -----------------------------------------------------------------------

    !       ! if (id_cond > 0) then
    !       ! do k = ktop, kbot ! total condensate
    !       ! cond (i) = cond (i) + dp1 (k) * (qlz (k) + qrz (k) + qsz (k) + qiz (k) + qgz (k))
    !       ! enddo
    !       ! endif
    !       !
    !       ! if (id_vtr > 0) then
    !       ! do k = ktop, kbot
    !       ! vt_r (i, j, k) = vtrz (k)
    !       ! enddo
    !       ! endif
    !       !
    !       ! if (id_vts > 0) then
    !       ! do k = ktop, kbot
    !       ! vt_s (i, j, k) = vtsz (k)
    !       ! enddo
    !       ! endif
    !       !
    !       ! if (id_vtg > 0) then
    !       ! do k = ktop, kbot
    !       ! vt_g (i, j, k) = vtgz (k)
    !       ! enddo
    !       ! endif
    !       !
    !       ! if (id_vts > 0) then
    !       ! do k = ktop, kbot
    !       ! vt_i (i, j, k) = vtiz (k)
    !       ! enddo
    !       ! endif
    !       !
    !       ! if (id_droplets > 0) then
    !       ! do k = ktop, kbot
    !       ! qn2 (i, j, k) = ccn (k)
    !       ! enddo
    !       ! endif

    !    enddo

    ! enddo
    ! ! !$omp end target teams distribute

  end subroutine mpdrv

  ! -----------------------------------------------------------------------
  !> sedimentation of heat
  ! -----------------------------------------------------------------------

  subroutine sedi_heat_3d (is, ie, js, je, ktop, kbot, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)

    implicit none

    ! TODO: This routine has not been optimized since we run with do_sedi_heat = .false.
    !$omp declare target

    ! input q fields are dry mixing ratios, and dm is dry air mass

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: dm, m1, dz, qv, ql, qr, qi, qs, qg

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: tz

    real, intent (in) :: cw ! heat capacity

    real :: dgz, cvn

    real :: tmp

    integer :: i, j, k

    ! -----------------------------------------------------------------------
    ! sjl, july 2014
    ! assumption: the ke in the falling condensates is negligible compared to the potential energy
    ! that was unaccounted for. local thermal equilibrium is assumed, and the loss in pe is transformed
    ! into internal energy (to heat the whole grid box)
    ! backward time - implicit upwind transport scheme:
    ! dm here is dry air mass
    ! -----------------------------------------------------------------------

    k = ktop
    !$omp target teams distribute parallel do collapse(2) private(cvn, dgz)
    do i = is, ie
       do j = js, je
          cvn = &
               dm (i, j, k) * (cv_air + qv (i, j, k) * cv_vap + &
               (qr (i, j, k) + ql (i, j, k)) * c_liq + &
               (qi (i, j, k) + qs (i, j, k) + qg (i, j, k)) * c_ice)
          dgz = - 0.5 * grav * dz (i, j, k) ! > 0
          tz (i, j, k) = ( &
               (cvn + m1 (i, j, k) * cw) * tz (i, j, k) + &
               m1 (i, j, k) * dgz &
               ) / (cvn + m1 (i, j, k) * cw)
       end do
    end do
    !$omp end target teams distribute parallel do

    ! -----------------------------------------------------------------------
    ! implicit algorithm: can't be vectorized
    ! needs an inner i - loop for vectorization
    ! -----------------------------------------------------------------------

    do k = ktop + 1, kbot
       !$omp target teams distribute parallel do collapse(2) private(cvn, dgz)
       do j = js, je
          do i = is, ie
             dgz = - 0.5 * grav * dz (i, j, k) ! > 0
             cvn = dm (i, j, k) * ( &
                  cv_air + &
                  qv (i, j, k) * cv_vap + &
                  (qr (i, j, k) + ql (i, j, k)) * c_liq + &
                  (qi (i, j, k) + qs (i, j, k) + qg (i, j, k)) * c_ice &
                  )
             tz (i, j, k) = ( &
                  (cvn + cw * (m1 (i, j, k) - m1 (i, j, k - 1))) * tz (i, j, k) + &
                  m1 (i, j, k - 1) * cw * tz (i, j, k - 1) + &
                  dgz * (m1 (i, j, k - 1) + m1 (i, j, k)) &
                  ) / (cvn + cw * m1 (i, j, k))
          end do
       end do
       !$omp end target teams distribute parallel do
    end do

  end subroutine sedi_heat_3d

  ! -----------------------------------------------------------------------
  !> warm rain cloud microphysics
  ! -----------------------------------------------------------------------

  subroutine warm_rain_3d ( &
       dt, is, ie, js, je, ktop, kbot, &
       dp, dz, tz, qv, ql, qr, qi, qs, qg, qa, &
       eis, den, denfac, ccn, c_praut, vtr, &
       r1, evap1, m1_rain, w1, h_var)

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in) :: dt !< time step (s)

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: h_var

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: dp, dz, den
    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: denfac, ccn, c_praut

    real, intent (in), dimension (is:ie, js:je) :: eis !< estimated inversion strength

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: tz, vtr
    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: evap1, m1_rain, w1

    real, intent (out), dimension (is:ie, js:je) :: r1

    real, parameter :: so3 = 7. / 3.

    real, dimension (is:ie, js:je, ktop:kbot) :: dl, dm, revap, isubl, qadum
    real, dimension (is:ie, js:je, ktop:kbot + 1) :: ze, zt

    real :: sink, dq, qc0, qc
    real :: fac_rc, qden
    real :: zs
    real :: dt5

    integer :: i, j, k

    ! fall velocity constants:

    real, parameter :: vconr = 2503.23638966667
    real, parameter :: normr = 25132741228.7183
    real, parameter :: thr = 1.e-8

    logical, dimension(is:ie, js:je) :: no_fall

    zs = 0.
    dt5 = 0.5 * dt

    !$omp target data &
    !$omp   map(to: h_var, dp, dz, den, denfac, ccn, c_praut, eis) &
    !$omp   map(tofrom: tz, vtr, qv, ql, qr, qi, qs, qg, qa, evap1, m1_rain, w1) &
    !$omp   map(from: r1) &
    !$omp   map(alloc: dl, dm, revap, isubl, qadum, ze, zt, no_fall)

    ! -----------------------------------------------------------------------
    ! terminal speed of rain
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             evap1 (i, j, k) = 0.
             m1_rain (i, j, k) = 0.
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    call check_column_3d (is, ie, js, je, ktop, kbot, qr, no_fall)

    ! -----------------------------------------------------------------------
    ! auto - conversion
    ! assuming linear subgrid vertical distribution of cloud water
    ! following lin et al. 1994, mwr
    ! -----------------------------------------------------------------------

    ! Use In-Cloud condensates
    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             if (.not. do_qa) then
                qadum (i, j, k) = max(qa (i, j, k) ,qcmin)
             else
                qadum (i, j, k) = 1.0
             endif
             ql (i, j, k) = ql (i, j, k) / qadum (i, j, k)
             qi (i, j, k) = qi (i, j, k) / qadum (i, j, k)
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    if (irain_f /= 0) then

       ! -----------------------------------------------------------------------
       ! no subgrid variability
       ! -----------------------------------------------------------------------

       !$omp target teams distribute parallel do simd collapse(3) private(fac_rc, qc0, qv, dq, sink)
       do k = ktop, kbot
          do i = is, ie
             do j = js, je
                fac_rc = min (1.0, eis (i, j) / 10.0) ** 2 ! Estimated inversion strength determine stable regime
                fac_rc = rc * (rthreshs * fac_rc + rthreshu * (1.0 - fac_rc)) ** 3
                qc0 = fac_rc * ccn (i, j, k)
                if (tz (i, j, k) > t_wfr) then
                   qc = qc0 / den (i, j, k)
                   dq = ql (i, j, k) - qc
                   if (dq > 0.) then
                      sink = min (dq, dt * c_praut (i, j, k) * den (i, j, k) * exp (so3 * log (ql (i, j, k))))
                      sink = min (ql0_max / qadum(i, j, k), ql(i, j, k), max(0., sink))
                      ql (i, j, k) = ql (i, j, k) - sink
                      qr (i, j, k) = qr (i, j, k) + sink * qadum(i, j, k)
                   end if
                end if
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

    else ! irain_f

       ! -----------------------------------------------------------------------
       ! with subgrid variability
       ! -----------------------------------------------------------------------

       call linear_prof_3d (ie - is + 1, je - js + 1, kbot - ktop + 1, ql, dl, z_slope_liq, h_var)

       !$omp target teams distribute parallel do simd collapse(3) private(fac_rc, qc0, qc, dq, sink)
       do k = ktop, kbot
          do i = is, ie
             do j = js, je
                fac_rc = min (1.0, eis (i, j) / 10.0) ** 2 ! Estimated inversion strength determine stable regime
                fac_rc = rc * (rthreshs * fac_rc + rthreshu * (1.0 - fac_rc)) ** 3
                qc0 = fac_rc * ccn (i, j, k)
                if (tz (i, j, k) > t_wfr + dt_fr) then
                   dl (i, j, k) = min (max (qcmin, dl (i, j, k)), 0.5 * ql (i, j, k))
                   ! --------------------------------------------------------------------
                   ! as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                   ! --------------------------------------------------------------------
                   qc = qc0 / den (i, j, k)
                   dq = 0.5 * (ql (i, j, k) + dl (i, j, k) - qc)
                   ! --------------------------------------------------------------------
                   ! dq = dl if qc == q_minus = ql - dl
                   ! dq = 0 if qc == q_plus = ql + dl
                   ! --------------------------------------------------------------------
                   if (dq > 0.) then ! q_plus > qc
                      ! --------------------------------------------------------------------
                      ! revised continuous form: linearly decays (with subgrid dl) to zero at qc == ql + dl
                      ! --------------------------------------------------------------------
                      sink = min (1., dq / dl (i, j, k)) * dt * c_praut (i, j, k) * den (i, j, k) * exp (so3 * log (ql (i, j, k)))
                      sink = min(ql0_max / qadum(i, j, k), ql(i, j, k), max(0., sink))
                      ql (i, j, k) = ql (i, j, k) - sink
                      qr (i, j, k) = qr (i, j, k) + sink*qadum(i, j, k)
                   end if
                end if
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

    end if ! irain_f

    ! Revert In-Cloud condensate

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             ql (i, j, k) = ql (i, j, k) * qadum (i, j, k)
             qi (i, j, k) = qi (i, j, k) * qadum (i, j, k)
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! fall speed of rain
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3) private(qden)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             if (no_fall (i, j)) then
                vtr (i, j, k) = vf_min
             else if (const_vr) then
                vtr (i, j, k) = vr_fac ! ifs_2016: 4.0
             else
                qden = qr (i, j, k) * den (i, j, k)
                if (qr (i, j, k) < thr) then
                   vtr (i, j, k) = vr_min
                else
                   vtr (i, j, k) = &
                        vr_fac * &
                        vconr * &
                        sqrt (min (10., sfcrho / den (i, j, k))) * &
                        exp (0.2 * log (qden / normr))
                   vtr (i, j, k) = min (vr_max, max (vr_min, vtr (i, j, k)))
                end if
             end if
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          ze (i, j, kbot + 1) = zs
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! !$omp ordered
    do k = kbot, ktop, - 1
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             ze (i, j, k) = ze (i, j, k + 1) - dz (i, j, k) ! dz < 0
          end do
       end do
       !$omp end target teams distribute parallel do simd
    end do
    ! !$omp end ordered

    ! -----------------------------------------------------------------------
    ! evaporation and accretion of rain for the first 1 / 2 time step
    ! -----------------------------------------------------------------------
    call revap_racc_3d (is, ie, js, je, ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)

    !$omp target teams distribute parallel do simd collapse(3) private(qden)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             evap1 (i, j, k) = revap (i, j, k)
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    if (do_sedi_w) then
       !$omp target teams distribute parallel do simd collapse(3)
       do k = ktop, kbot
          do i = is, ie
             do j = js, je
                dm (i, j, k) = dp (i, j, k) * ( &
                     1. + &
                     qv (i, j, k) + &
                     ql (i, j, k) + &
                     qr (i, j, k) + &
                     qi (i, j, k) + &
                     qs (i, j, k) + &
                     qg (i, j, k))
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd
    endif

    ! -----------------------------------------------------------------------
    ! mass flux induced by falling rain
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(2)
    do j = js, je
       do i = is, ie
          if (no_fall (i, j)) then
             r1 (i, j) = 0.0
          end if
       end do
    end do
    !$omp end target teams distribute parallel do simd

    if (use_ppm) then

       !$omp target teams distribute collapse(2)
       do i = is, ie
          do j = js, je

             if (.not. no_fall (i, j)) then
                zt (i, j, ktop) = ze (i, j, ktop)
                !$omp parallel do simd
                do k = ktop + 1, kbot
                   zt (i, j, k) = ze (i, j, k) - dt * (vtr (i, j, k - 1) + vtr (i, j, k)) / 2.0
                enddo
                !$omp end parallel do simd
                zt (i, j, kbot + 1) = zs - dt * vtr (i, j, kbot)
                ! !$omp ordered
                do k = ktop, kbot
                   if (zt (i, j, k + 1) >= zt (i, j, k)) zt (i, j, k + 1) = zt (i, j, k) - dz_min
                enddo
                ! !$omp end ordered
                call lagrangian_fall_ppm ( &
                     ktop, kbot, zs, ze (i, j, :), zt (i, j, :), &
                     dp (i, j, :), qr (i, j, :), r1 (i, j), m1_rain (i, j, :), mono_prof)
             end if ! not no_fall

          end do
       end do
       !$omp end target teams distribute

    else

       call implicit_fall_3d (dt, is, ie, js, je, ktop, kbot, ze, vtr, dp, qr, r1, m1_rain, no_fall)

    endif

    ! -----------------------------------------------------------------------
    ! vertical velocity transportation during sedimentation
    ! -----------------------------------------------------------------------

    if (do_sedi_w) then

       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             w1 (i, j, ktop) = ( &
                  dm (i, j, ktop) * w1 (i, j, ktop) + &
                  m1_rain (i, j, ktop) * vtr (i, j, ktop) &
                  ) / (dm (i, j, ktop) - m1_rain (i, j, ktop))
          end do
       end do
       !$omp end target teams distribute parallel do simd

       !$omp target teams distribute parallel do simd collapse(3)
       do k = ktop + 1, kbot
          do i = is, ie
             do j = js, je
                w1 (i, j, k) = ( &
                     dm (i, j, k) * w1 (i, j, k) - &
                     m1_rain (i, j, k - 1) * vtr (i, j, k - 1) + &
                     m1_rain (i, j, k) * vtr (i, j, k) &
                     ) / (dm (i, j, k) + m1_rain (i, j, k - 1) - m1_rain (i, j, k))
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

    end if

    ! -----------------------------------------------------------------------
    ! heat transportation during sedimentation
    ! -----------------------------------------------------------------------

    if (do_sedi_heat) then
       call sedi_heat_3d (is, ie, js, je, ktop, kbot, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)
    end if

    ! -----------------------------------------------------------------------
    ! evaporation and accretion of rain for the remaing 1 / 2 time step
    ! -----------------------------------------------------------------------

    call revap_racc_3d (is, ie, js, je, ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop + 1, kbot
       do i = is, ie
          do j = js, je
             evap1 (i, j, k) = evap1 (i, j, k) + revap (i, j, k)
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp end target data

  end subroutine warm_rain_3d

  ! -----------------------------------------------------------------------
  !> evaporation of rain
  ! -----------------------------------------------------------------------

  subroutine revap_racc_3d (is, ie, js, je, ktop, kbot, dt, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in) :: dt ! time step (s)

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: h_var
    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: den, denfac

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: tz, qv, qr, ql, qi, qs, qg, qa

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: revap

    real :: lhl, cvm, q_liq, q_sol, lcpk
    real :: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin
    real :: fac_revp
    real :: TOT_PREC_LS, AREA_LS_PRC, AREA_LS_PRC_K
    integer :: i, j, k

    TOT_PREC_LS = 0.
    AREA_LS_PRC = 0.

    !$omp target data &
    !$omp   map(to: h_var, den, denfac) &
    !$omp   map(tofrom: tz, qv, qr, ql, qi, qs, qg, qa, revap)

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             revap(i, j, k) = 0.
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    !$omp target teams distribute collapse(2)
    do k = ktop, kbot
       do i = is, ie

          !$omp parallel do simd &
          !$omp   firstprivate( &
          !$omp     TOT_PREC_LS, AREA_LS_PRC) &
          !$omp   private( &
          !$omp     fac_revp, lhl, q_liq, q_sol, cvm, lcpk, tin, qpz, &
          !$omp     qsat, dqh, dqv, q_minus, q_plus, dq, qden, t2, evap, &
          !$omp     sink, dqsdt)
          do j = js, je

             TOT_PREC_LS = TOT_PREC_LS  + ( ( qr (i, j, k) + qs (i, j, k) + qg (i, j, k) ) * den (i, j, k) )
             AREA_LS_PRC = AREA_LS_PRC  + ( qa (i, j, k) * ( qr (i, j, k) + qs (i, j, k) + qg (i, j, k) ) * den (i, j, k) )

             if (tz (i, j, k) > t_wfr .and. qr (i, j, k) > qpmin) then

                !! area and timescale efficiency on revap
                !                       AREA_LS_PRC_K = 0.0
                !if (TOT_PREC_LS > 0.0) AREA_LS_PRC_K = MAX( AREA_LS_PRC/TOT_PREC_LS, 1.E-6 )
                !fac_revp = 1. - exp (- AREA_LS_PRC_K * dt / tau_revp)
                fac_revp = 1. - exp (- dt / tau_revp)

                ! -----------------------------------------------------------------------
                ! define heat capacity and latent heat coefficient
                ! -----------------------------------------------------------------------

                lhl = lv00 + d0_vap * tz (i, j, k)
                q_liq = ql (i, j, k) + qr (i, j, k)
                q_sol = qi (i, j, k) + qs (i, j, k) + qg (i, j, k)
                cvm = c_air + qv (i, j, k) * c_vap + q_liq * c_liq + q_sol * c_ice
                lcpk = lhl / cvm

                tin = tz (i, j, k) - lcpk * ql (i, j, k) ! presence of clouds suppresses the rain evap
                qpz = qv (i, j, k) + ql (i, j, k)
                qsat = wqs2 (tin, den (i, j, k), dqsdt)
                dqh = max (ql (i, j, k), h_var(i, j, k) * max (qpz, qcmin))
                dqh = min (dqh, 0.2 * qpz) ! new limiter
                dqv = qsat - qv (i, j, k) ! use this to prevent super - sat the gird box
                q_minus = qpz - dqh
                q_plus = qpz + dqh

                ! -----------------------------------------------------------------------
                ! qsat must be > q_minus to activate evaporation
                ! qsat must be < q_plus to activate accretion
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! rain evaporation
                ! -----------------------------------------------------------------------

                if (dqv > qvmin .and. qsat > q_minus) then
                   if (qsat > q_plus) then
                      dq = qsat - qpz
                   else
                      ! -----------------------------------------------------------------------
                      ! q_minus < qsat < q_plus
                      ! dq == dqh if qsat == q_minus
                      ! -----------------------------------------------------------------------
                      dq = 0.25 * (q_minus - qsat) ** 2 / dqh
                   endif
                   qden = qr (i, j, k) * den (i, j, k)
                   t2 = tin * tin
                   evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                        exp (0.725 * log (qden))) / (crevp (4) * t2 + crevp (5) * qsat * den (i, j, k))
                   evap = min (qr (i, j, k), dt * fac_revp * evap, dqv / (1. + lcpk * dqsdt))
                   qr (i, j, k) = qr (i, j, k) - evap
                   qv (i, j, k) = qv (i, j, k) + evap
                   q_liq = q_liq - evap
                   cvm = c_air + qv (i, j, k) * c_vap + q_liq * c_liq + q_sol * c_ice
                   tz (i, j, k) = tz (i, j, k) - evap * lhl / cvm
                   revap (i, j, k) = evap / dt
                endif

                ! -----------------------------------------------------------------------
                ! accretion: pracc
                ! -----------------------------------------------------------------------

                if (qr (i, j, k) > qpmin .and. ql (i, j, k) > qcmin .and. qsat < q_minus) then
                   sink = dt * denfac (i, j, k) * cracw * exp (0.95 * log (qr (i, j, k) * den (i, j, k)))
                   sink = sink / (1. + sink) * ql (i, j, k)
                   ql (i, j, k) = ql (i, j, k) - sink
                   qr (i, j, k) = qr (i, j, k) + sink
                endif

             end if ! warm - rain

          end do
          !$omp end parallel do

       end do
    end do
    !$omp end target teams distribute

    !$omp end target data

  end subroutine revap_racc_3d

  ! -----------------------------------------------------------------------
  !> definition of vertical subgrid variability
  !! used for cloud ice and cloud water autoconversion
  !! qi -- > ql & ql -- > qr
  !! edges: qe == qbar + / - dm
  ! -----------------------------------------------------------------------

  subroutine linear_prof_3d (im, jm, km, q, dm, z_var, h_var)

    implicit none

    integer, intent (in) :: im, jm, km

    real, intent (in) :: q (im, jm, km), h_var(im, jm, km)

    real, intent (out) :: dm (im, jm, km)

    logical, intent (in) :: z_var

    real :: dq, dq_p1

    integer :: i, j, k

    !$omp target data map(to: q, h_var) map(from: dm)

    if (z_var) then

       !$omp target teams distribute parallel do simd collapse(2)
       do i = 1, im
          do j = 1, jm
             dm (i, j, 1) = 0.
          end do
       end do
       !$omp end target teams distribute parallel do simd

       ! -----------------------------------------------------------------------
       ! use twice the strength of the positive definiteness limiter (lin et al 1994)
       ! -----------------------------------------------------------------------

       !$omp target teams distribute parallel do simd collapse(3) private(dq, dq_p1)
       do k = 2, km - 1
          do i = 1, im
             do j = 1, jm
                dq = 0.5 * (q (i, j, k) - q (i, j, k - 1))
                dq_p1 = 0.5 * (q (i, j, k + 1) - q (i, j, k))
                dm (i, j, k) = 0.5 * min (abs (dq + dq_p1), 0.5 * q (i, j, k))
                if (dq * dq_p1 <= 0.) then
                   if (dq > 0.) then ! local max
                      dm (i, j, k) = min (dm (i, j, k), dq, - dq_p1)
                   else
                      dm (i, j, k) = 0.
                   endif
                endif
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

       !$omp target teams distribute parallel do simd collapse(2)
       do i = 1, im
          do j = 1, jm
             dm (i, j, km) = 0.
          end do
       end do
       !$omp end target teams distribute parallel do simd

       ! -----------------------------------------------------------------------
       ! impose a presumed background horizontal variability that is proportional to the value itself
       ! -----------------------------------------------------------------------

       !$omp target teams distribute parallel do simd collapse(3)
       do k = 1, km
          do i = 1, im
             do j = 1, jm
                dm (i, j, k) = max (dm (i, j, k), qvmin, h_var(i, j, k) * q (i, j, k))
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

    else ! z_var


       !$omp target teams distribute parallel do simd collapse(3)
       do k = 1, km
          do i = 1, im
             do j = 1, jm
                dm (i, j, k) = max (qvmin, h_var(i, j, k) * q (i, j, k))
             end do
          end do
       end do
       !$omp end target teams distribute parallel do simd

    endif

    !$omp end target data

  end subroutine linear_prof_3d

  ! =======================================================================
  !> ice cloud microphysics processes
  !! bulk cloud micro - physics; processes splitting
  !! with some un - split sub - grouping
  !! time implicit (when possible) accretion and autoconversion
  !>@author: Shian-Jiann lin, gfdl
  ! =======================================================================

  subroutine icloud_3d ( &
       is, ie, js, je, ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, &
       den, denfac, vts, vtg, vtr, qak, dts, subl1, h_var, ccn, cnv_fraction, srf_type)

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: p1, dp1, den, denfac, vts, vtg, vtr

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: tzk, qvk, qlk, qrk, qik, qsk, qgk, qak

    real, intent (out), dimension (is:ie, js:je, ktop:kbot) :: subl1

    real, intent (in) :: dts

    real, intent(in) , dimension (is:ie, js:je) :: cnv_fraction, srf_type

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: h_var, ccn

    real, dimension (is:ie, js:je, ktop:kbot) :: di
    real, dimension (is:ie, js:je, ktop:kbot) :: cvm, q_liq, q_sol

    real :: icpk, tcpk, lhl, lhi
    real :: rdts, fac_g2v, fac_i2s, fac_imlt, fac_frz
    real :: tz, qv, ql, qr, qi, qs, qg, melt, ifrac, newqi, newql
    real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, praci, psaci
    real :: pgmlt, psmlt, pgfr, pgaut, psaut, pgsub
    real :: tc, tsq, dqs0, qden, qim, qsm
    real :: factor, sink, qi_crt
    real :: tmp, qsw, qsi, dqsdt, dq
    real :: dtmp, qc, q_plus, q_minus

    integer :: i, j, k, it

    !$omp target data &
    !$omp   map(to: p1, dp1, den, denfac, vts, vtg, vtr, cnv_fraction, srf_type, h_var, ccn) &
    !$omp   map(tofrom: tzk, qvk, qlk, qrk, qik, qsk, qgk, qak) &
    !$omp   map(from: subl1) &
    !$omp   map(alloc: di, cvm, q_liq, q_sol)

    rdts = 1. / dts

    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------

    fac_i2s = 1. - exp (- dts / tau_i2s)
    fac_g2v = 1. - exp (- dts / tau_g2v)

    fac_imlt = 1. - exp (- dts / tau_imlt)
    fac_frz  = 1. - exp (- dts / tau_frz)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             q_liq (i, j, k) = qlk (i, j, k) + qrk (i, j, k)
             q_sol (i, j, k) = qik (i, j, k) + qsk (i, j, k) + qgk (i, j, k)
             cvm (i, j, k) = c_air + qvk (i, j, k) * c_vap + q_liq (i, j, k) * c_liq + q_sol (i, j, k) * c_ice
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! sources of cloud ice: pihom, cold rain, and the sat_adj
    ! (initiation plus deposition)
    ! sources of snow: cold rain, auto conversion + accretion (from cloud ice)
    ! sat_adj (deposition; requires pre - existing snow) ; initial snow comes from auto conversion
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3) private(lhi, icpk, melt, tmp, sink, qi_crt)
    do k = ktop, kbot
       do i = is, ie
          ! !$omp parallel do simd - gfortran?
          do j = js, je

             lhi = li00 + dc_ice * tzk (i, j, k)

             if (tzk (i, j, k) > tice .and. qik (i, j, k) > qcmin) then

                icpk = lhi / cvm (i, j, k)
                ! -----------------------------------------------------------------------
                ! pimlt: instant melting of cloud ice
                ! -----------------------------------------------------------------------

                melt = min (qik (i, j, k), fac_imlt * (tzk (i, j, k) - tice) / icpk)
                tmp = min (melt, dim (ql_mlt, qlk (i, j, k))) ! max ql amount

                ! new total condensate / old condensate
                qak(i, j, k) = max (0.0, min (1., qak (i, j, k) * max (qik (i, j, k) + qlk (i, j, k) - melt + tmp, 0.0) / &
                     max (qik (i, j, k) + qlk (i, j, k), qcmin)))

                qlk (i, j, k) = qlk (i, j, k) + tmp
                qrk (i, j, k) = qrk (i, j, k) + melt - tmp
                qik (i, j, k) = qik (i, j, k) - melt
                q_liq (i, j, k) = q_liq (i, j, k) + melt
                q_sol (i, j, k) = q_sol (i, j, k) - melt
                cvm (i, j, k) = c_air + qvk (i, j, k) * c_vap + q_liq (i, j, k) * c_liq + q_sol (i, j, k) * c_ice
                tzk (i, j, k) = tzk (i, j, k) - melt * lhi / cvm (i, j, k)

             elseif (tzk (i, j, k) <= tice .and. qlk (i, j, k) > qcmin) then

                ! -----------------------------------------------------------------------
                ! pihom: homogeneous freezing of cloud water into cloud ice
                ! this is the 1st occurance of liquid water freezing in the split mp process
                ! -----------------------------------------------------------------------

                sink = fac_frz * &
                     new_ice_condensate (tzk (i, j, k), qlk (i, j, k), qik (i, j, k), cnv_fraction (i, j), srf_type (i, j))
                qi_crt = qi0_crt / den (i, j, k)
                tmp = min (sink, dim (qi_crt, qik (i, j, k)))

                ! new total condensate / old condensate
                qak (i, j, k) = max (0.0, min (1., qak(i, j, k) * max (qik (i, j, k) + qlk (i, j, k) - sink + tmp, 0.0) / &
                     max (qik (i, j, k) + qlk (i, j, k), qcmin)))

                qlk (i, j, k) = qlk (i, j, k) - sink
                qsk (i, j, k) = qsk (i, j, k) + sink - tmp
                qik (i, j, k) = qik (i, j, k) + tmp
                q_liq (i, j, k) = q_liq (i, j, k) - sink
                q_sol (i, j, k) = q_sol (i, j, k) + sink
                cvm (i, j, k) = c_air + qvk (i, j, k) * c_vap + q_liq (i, j, k) * c_liq + q_sol (i, j, k) * c_ice
                tzk (i, j, k) = tzk (i, j, k) + sink * lhi / cvm (i, j, k)

             end if

          end do
          ! !$omp end parallel do - gfortran?
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! ! ! -----------------------------------------------------------------------
    ! ! ! vertical subgrid variability
    ! ! ! -----------------------------------------------------------------------

    ! ! call linear_prof (kbot - ktop + 1, qik (ktop), di (ktop), z_slope_ice, h_var)

    ! ! -----------------------------------------------------------------------
    ! ! update capacity heat and latend heat coefficient
    ! ! -----------------------------------------------------------------------

    ! !!$acc loop seq
    ! do k = ktop, kbot
    !    lhl = lv00 + d0_vap * tzk (k)
    !    lhi = li00 + dc_ice * tzk (k)
    !    icpk = lhi / cvm (k)
    !    tcpk = lhl / cvm (k) + icpk

    !    ! -----------------------------------------------------------------------
    !    ! do nothing above p_min
    !    ! -----------------------------------------------------------------------

    !    if (p1 (k) < p_min) cycle

    !    tz = tzk (k)
    !    qv = qvk (k)
    !    ql = qlk (k)
    !    qi = qik (k)
    !    qr = qrk (k)
    !    qs = qsk (k)
    !    qg = qgk (k)

    !    pgacr = 0.
    !    pgacw = 0.
    !    tc = tz - tice

    !    if (tc .ge. 0.) then

    !       ! -----------------------------------------------------------------------
    !       ! melting of snow
    !       ! -----------------------------------------------------------------------

    !       dqs0 = ces0 / p1 (k) - qv

    !       if (qs > qpmin) then

    !          ! -----------------------------------------------------------------------
    !          ! psacw: accretion of cloud water by snow
    !          ! only rate is used (for snow melt) since tc > 0.
    !          ! -----------------------------------------------------------------------

    !          if (ql > qcmin) then
    !             factor = denfac (k) * csacw * exp (0.8125 * log (qs * den (k)))
    !             psacw = factor / (1. + dts * factor) * ql ! rate
    !          else
    !             psacw = 0.
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! psacr: accretion of rain by melted snow
    !          ! pracs: accretion of snow by rain
    !          ! -----------------------------------------------------------------------

    !          if (qr > qpmin) then
    !             psacr = min (acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), &
    !                  den (k)), qr * rdts)
    !             pracs = acr3d (vtr (k), vts (k), qs, qr, cracs, acco (1, 1), den (k))
    !          else
    !             psacr = 0.
    !             pracs = 0.
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! total snow sink:
    !          ! psmlt: snow melt (due to rain accretion)
    !          ! -----------------------------------------------------------------------

    !          psmlt = max (0., smlt (tc, dqs0, qs * den (k), psacw, psacr, csmlt, &
    !               den (k), denfac (k)))
    !          sink = min (qs, dts * (psmlt + pracs), tc / icpk)
    !          qs = qs - sink
    !          ! sjl, 20170321:
    !          tmp = min (sink, dim (qs_mlt, ql)) ! max ql due to snow melt

    !          ! new total condensate / old condensate
    !          qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql+tmp,0.0  ) / &
    !               max(qi+ql    ,qcmin) ) )

    !          ql = ql + tmp
    !          qr = qr + sink - tmp
    !          ! qr = qr + sink
    !          ! sjl, 20170321:
    !          q_liq (k) = q_liq (k) + sink
    !          q_sol (k) = q_sol (k) - sink
    !          cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
    !          tz = tz - sink * lhi / cvm (k)
    !          tc = tz - tice

    !       endif

    !       ! -----------------------------------------------------------------------
    !       ! update capacity heat and latend heat coefficient
    !       ! -----------------------------------------------------------------------

    !       lhi = li00 + dc_ice * tz
    !       icpk = lhi / cvm (k)

    !       ! -----------------------------------------------------------------------
    !       ! melting of graupel
    !       ! -----------------------------------------------------------------------

    !       if (qg > qpmin .and. tc > 0.) then

    !          ! -----------------------------------------------------------------------
    !          ! pgacr: accretion of rain by graupel
    !          ! -----------------------------------------------------------------------

    !          if (qr > qpmin) &
    !               pgacr = min (acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
    !               den (k)), rdts * qr)

    !          ! -----------------------------------------------------------------------
    !          ! pgacw: accretion of cloud water by graupel
    !          ! -----------------------------------------------------------------------

    !          qden = qg * den (k)
    !          if (ql > qcmin) then
    !             factor = cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
    !             pgacw = factor / (1. + dts * factor) * ql ! rate
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! pgmlt: graupel melt
    !          ! -----------------------------------------------------------------------

    !          pgmlt = dts * gmlt (tc, dqs0, qden, pgacw, pgacr, cgmlt, den (k))
    !          pgmlt = min (max (0., pgmlt), qg, tc / icpk)
    !          qg = qg - pgmlt
    !          qr = qr + pgmlt
    !          q_liq (k) = q_liq (k) + pgmlt
    !          q_sol (k) = q_sol (k) - pgmlt
    !          cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
    !          tz = tz - pgmlt * lhi / cvm (k)

    !       endif

    !    else

    !       ! -----------------------------------------------------------------------
    !       ! cloud ice proc:
    !       ! -----------------------------------------------------------------------

    !       ! -----------------------------------------------------------------------
    !       ! psaci: accretion of cloud ice by snow
    !       ! -----------------------------------------------------------------------

    !       if (qi > 3.e-7) then ! cloud ice sink terms

    !          if (qs > qpmin) then
    !             ! -----------------------------------------------------------------------
    !             ! sjl added (following lin eq. 23) the temperature dependency
    !             ! to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
    !             ! -----------------------------------------------------------------------
    !             factor = dts * denfac (k) * csaci * exp (0.05 * tc + 0.8125 * log (qs * den (k)))
    !             psaci = factor / (1. + factor) * qi
    !          else
    !             psaci = 0.
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! psaut: autoconversion: cloud ice -- > snow
    !          ! -----------------------------------------------------------------------

    !          ! -----------------------------------------------------------------------
    !          ! similar to lfo 1983: eq. 21 solved implicitly
    !          ! threshold from wsm6 scheme, hong et al 2004, eq (13) : qi0_crt ~0.8e-4
    !          ! -----------------------------------------------------------------------

    !          ! !!!!!!!        qim = qi0_crt / den (k)
    !          ! GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C on qi0_crt
    !          qim = ice_fraction(tz,cnv_fraction,srf_type) * qi0_crt / den (k)

    !          ! -----------------------------------------------------------------------
    !          ! assuming linear subgrid vertical distribution of cloud ice
    !          ! the mismatch computation following lin et al. 1994, mwr
    !          ! -----------------------------------------------------------------------

    !          if (const_vi) then
    !             tmp = fac_i2s
    !          else
    !             tmp = fac_i2s * exp (0.025 * tc)
    !          endif

    !          di (k) = max (di (k), qcmin)
    !          q_plus = qi + di (k)
    !          if (q_plus > (qim + qcmin)) then
    !             if (qim > (qi - di (k))) then
    !                dq = (0.25 * (q_plus - qim) ** 2) / di (k)
    !             else
    !                dq = qi - qim
    !             endif
    !             psaut = tmp * dq
    !          else
    !             psaut = 0.
    !          endif
    !          sink = min (qi, psaci + psaut)

    !          ! new total condensate / old condensate
    !          qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql-sink+tmp,0.0  ) / &
    !               max(qi+ql         ,qcmin) ) )

    !          qi = qi - sink
    !          qs = qs + sink

    !          ! -----------------------------------------------------------------------
    !          ! pgaci: accretion of cloud ice by graupel
    !          ! -----------------------------------------------------------------------

    !          if (qg > qpmin) then
    !             ! -----------------------------------------------------------------------
    !             ! factor = dts * cgaci / sqrt (den (k)) * exp (0.05 * tc + 0.875 * log (qg * den (k)))
    !             ! simplified form: remove temp dependency & set the exponent "0.875" -- > 1
    !             ! -----------------------------------------------------------------------
    !             factor = dts * cgaci * sqrt (den (k)) * qg
    !             pgaci = factor / (1. + factor) * qi
    !             qi = qi - pgaci
    !             qg = qg + pgaci
    !          endif

    !       endif

    !       ! -----------------------------------------------------------------------
    !       ! cold - rain proc:
    !       ! -----------------------------------------------------------------------

    !       ! -----------------------------------------------------------------------
    !       ! rain to ice, snow, graupel processes:
    !       ! -----------------------------------------------------------------------

    !       tc = tz - tice

    !       if (qr > qpmin .and. tc < 0.) then

    !          ! -----------------------------------------------------------------------
    !          ! * sink * terms to qr: psacr + pgfr
    !          ! source terms to qs: psacr
    !          ! source terms to qg: pgfr
    !          ! -----------------------------------------------------------------------

    !          ! -----------------------------------------------------------------------
    !          ! psacr accretion of rain by snow
    !          ! -----------------------------------------------------------------------

    !          if (qs > qpmin) then ! if snow exists
    !             psacr = dts * acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), den (k))
    !          else
    !             psacr = 0.
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! pgfr: rain freezing -- > graupel
    !          ! -----------------------------------------------------------------------

    !          pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
    !               exp (1.75 * log (qr * den (k)))

    !          ! -----------------------------------------------------------------------
    !          ! total sink to qr
    !          ! -----------------------------------------------------------------------

    !          sink = psacr + pgfr
    !          factor = min (sink, qr, - tc / icpk) / max (sink, qpmin)

    !          psacr = factor * psacr
    !          pgfr = factor * pgfr

    !          sink = psacr + pgfr
    !          qr = qr - sink
    !          qs = qs + psacr
    !          qg = qg + pgfr
    !          q_liq (k) = q_liq (k) - sink
    !          q_sol (k) = q_sol (k) + sink
    !          cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
    !          tz = tz + sink * lhi / cvm (k)

    !       endif

    !       ! -----------------------------------------------------------------------
    !       ! update capacity heat and latend heat coefficient
    !       ! -----------------------------------------------------------------------

    !       lhi = li00 + dc_ice * tz
    !       icpk = lhi / cvm (k)

    !       ! -----------------------------------------------------------------------
    !       ! graupel production terms:
    !       ! -----------------------------------------------------------------------

    !       if (qs > qpmin) then

    !          ! -----------------------------------------------------------------------
    !          ! accretion: snow -- > graupel
    !          ! -----------------------------------------------------------------------

    !          if (qg > qpmin) then
    !             sink = dts * acr3d (vtg (k), vts (k), qs, qg, cgacs, acco (1, 4), den (k))
    !          else
    !             sink = 0.
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! autoconversion snow -- > graupel
    !          ! -----------------------------------------------------------------------

    !          qsm = qs0_crt / den (k)
    !          if (qs > qsm) then
    !             factor = dts * 1.e-3 * exp (0.09 * (tz - tice))
    !             sink = sink + factor / (1. + factor) * (qs - qsm)
    !          endif
    !          sink = min (qs, sink)
    !          qs = qs - sink
    !          qg = qg + sink

    !       endif ! snow existed

    !       if (qg > qpmin .and. tz < tice0) then

    !          ! -----------------------------------------------------------------------
    !          ! pgacw: accretion of cloud water by graupel
    !          ! -----------------------------------------------------------------------

    !          if (ql > qcmin) then
    !             qden = qg * den (k)
    !             factor = dts * cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
    !             pgacw = factor / (1. + factor) * ql
    !          else
    !             pgacw = 0.
    !          endif

    !          ! -----------------------------------------------------------------------
    !          ! pgacr: accretion of rain by graupel
    !          ! -----------------------------------------------------------------------

    !          if (qr > qpmin) then
    !             pgacr = min (dts * acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
    !                  den (k)), qr)
    !          else
    !             pgacr = 0.
    !          endif

    !          sink = pgacr + pgacw
    !          factor = min (sink, dim (tice, tz) / icpk) / max (sink, qpmin)
    !          pgacr = factor * pgacr
    !          pgacw = factor * pgacw

    !          sink = pgacr + pgacw
    !          qg = qg + sink
    !          qr = qr - pgacr
    !          ql = ql - pgacw
    !          q_liq (k) = q_liq (k) - sink
    !          q_sol (k) = q_sol (k) + sink
    !          cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
    !          tz = tz + sink * lhi / cvm (k)

    !       endif

    !    endif

    !    tzk (k) = tz
    !    qvk (k) = qv
    !    qlk (k) = ql
    !    qik (k) = qi
    !    qrk (k) = qr
    !    qsk (k) = qs
    !    qgk (k) = qg

    ! enddo

    ! ! -----------------------------------------------------------------------
    ! ! subgrid cloud microphysics
    ! ! -----------------------------------------------------------------------

    ! call subgrid_z_proc (ktop, kbot, p1, den, denfac, dts, tzk, qvk, &
    !      qlk, qrk, qik, qsk, qgk, qak, subl1, h_var, ccn, cnv_fraction, srf_type)

    !$omp end target data

  end subroutine icloud_3d

  ! =======================================================================
  !>temperature sensitive high vertical resolution processes
  ! =======================================================================

  subroutine subgrid_z_proc (ktop, kbot, p1, den, denfac, dts, tz, qv, &
       ql, qr, qi, qs, qg, qa, subl1, h_var, ccn, cnv_fraction, srf_type)

    implicit none
    !$omp declare target

    integer, intent (in) :: ktop, kbot

    real, intent (in), dimension (ktop:kbot) :: p1, den, denfac

    real, intent (in) :: dts, cnv_fraction, srf_type

    real, intent (in), dimension (ktop:kbot) :: h_var, ccn

    real, intent (inout), dimension (ktop:kbot) :: tz, qv, ql, qr, qi, qs, qg, qa

    real, intent (out), dimension (ktop:kbot) :: subl1

    real, dimension (ktop:kbot) :: lcpk, icpk, tcpk, tcp3, lhl, lhi
    real, dimension (ktop:kbot) :: cvm, q_liq, q_sol, q_cond

    real :: fac_v2l, fac_l2v, fac_i2v

    real :: pidep, qi_crt

    ! -----------------------------------------------------------------------
    ! qstar over water may be accurate only down to - 80 deg c with ~10% uncertainty
    ! must not be too large to allow psc
    ! -----------------------------------------------------------------------

    real :: rh, rqi, tin, qsw, qsi, qpz, qstar
    real :: dqsdt, dwsdt, dq, dq0, factor, tmp, oldqa
    real :: dqh, q_plus, q_minus, dt_evap
    real :: evap, subl, sink, tc, pisub, q_adj, dtmp
    real :: pssub, pgsub, tsq, qden, fac_g2v, fac_v2g, fac_s2v, fac_v2s
    real :: ifrac, newqi, fac_frz
    real :: rh_adj, rh_rain

    integer :: k

    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------

    fac_v2l = 1. - exp (- dts / tau_v2l)
    fac_l2v = 1. - exp (- dts / tau_l2v)
    fac_i2v = 1. - exp (- dts / tau_i2v)
    fac_s2v = 1. - exp (- dts / tau_s2v)
    fac_v2s = 1. - exp (- dts / tau_v2s)
    fac_g2v = 1. - exp (- dts / tau_g2v)
    fac_v2g = 1. - exp (- dts / tau_v2g)
    fac_frz = 1. - exp (- dts / tau_frz)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
       lhl (k) = lv00 + d0_vap * tz (k)
       lhi (k) = li00 + dc_ice * tz (k)
       q_liq (k) = ql (k) + qr (k)
       q_sol (k) = qi (k) + qs (k) + qg (k)
       cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
       lcpk (k) = lhl (k) / cvm (k)
       icpk (k) = lhi (k) / cvm (k)
       tcpk (k) = lcpk (k) + icpk (k)
       tcp3 (k) = lcpk (k) + icpk (k) * min (1., dim (tice, tz (k)) / (tice - t_wfr))
    enddo

    do k = ktop, kbot

       rh_adj = 1. - h_var(k) - rh_inc
       rh_rain = max (0.35, rh_adj - rh_inr)

       subl1(k) = 0.0

       if (p1 (k) < p_min) cycle

       ! -----------------------------------------------------------------------
       ! instant deposit all water vapor to cloud ice when temperature is super low
       ! -----------------------------------------------------------------------

       if (tz (k) < t_min) then
          sink = dim (qv (k), qvmin)
          qv (k) = qv (k) - sink
          qi (k) = qi (k) + sink
          q_sol (k) = q_sol (k) + sink
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) + sink * (lhl (k) + lhi (k)) / cvm (k)
          if (do_qa) qa (k) = 1. ! air fully saturated; 100 % cloud cover
          cycle
       endif

       ! -----------------------------------------------------------------------
       ! update heat capacity and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhl (k) = lv00 + d0_vap * tz (k)
       lhi (k) = li00 + dc_ice * tz (k)
       lcpk (k) = lhl (k) / cvm (k)
       icpk (k) = lhi (k) / cvm (k)
       tcpk (k) = lcpk (k) + icpk (k)
       tcp3 (k) = lcpk (k) + icpk (k) * min (1., dim (tice, tz (k)) / (tice - t_wfr))

       ! -----------------------------------------------------------------------
       ! instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
       ! -----------------------------------------------------------------------
       qpz = qv (k) + ql (k) + qi (k)
       tin = tz (k) - (lhl (k) * (ql (k) + qi (k)) + lhi (k) * qi (k)) / (c_air + &
            qpz * c_vap + qr (k) * c_liq + (qs (k) + qg (k)) * c_ice)
       if (tin > t_sub + 6.) then
          rh = qpz / iqs1 (tin, den (k))
          if (rh < rh_adj) then ! qpz / rh_adj < qs
             tz (k) = tin
             qv (k) = qpz
             ql (k) = 0.
             qi (k) = 0.
             if (do_qa) qa (k) = 0.
             cycle ! cloud free
          endif
       endif

       ! -----------------------------------------------------------------------
       ! cloud water < -- > vapor adjustment: LS evaporation
       ! -----------------------------------------------------------------------
       if (do_evap) then
          qsw = wqs2 (tz (k), den (k), dwsdt)
          dq0 = qsw - qv (k)
          if (dq0 > qvmin) then
             factor = min (1., fac_l2v * (10. * dq0 / qsw))
             evap = min (ql (k), factor * ql(k) / (1. + tcp3 (k) * dwsdt))
          else
             evap = 0.0
          endif
          qv (k) = qv (k) + evap
          ql (k) = ql (k) - evap
          q_liq (k) = q_liq (k) - evap
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) - evap * lhl (k) / cvm (k)
       endif

       ! -----------------------------------------------------------------------
       ! update heat capacity and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhi (k) = li00 + dc_ice * tz (k)
       icpk (k) = lhi (k) / cvm (k)

       ! -----------------------------------------------------------------------
       ! enforce complete freezing below - t_wfr
       ! -----------------------------------------------------------------------

       dtmp = t_wfr - tz (k)
       if (dtmp > 0. .and. ql (k) > qcmin) then
          sink = min (ql (k), fac_frz * dtmp / icpk (k))
          ql (k) = ql (k) - sink
          qi (k) = qi (k) + sink
          q_liq (k) = q_liq (k) - sink
          q_sol (k) = q_sol (k) + sink
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) + sink * lhi (k) / cvm (k)
       endif

       ! -----------------------------------------------------------------------
       ! update heat capacity and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhi (k) = li00 + dc_ice * tz (k)
       icpk (k) = lhi (k) / cvm (k)

       ! -----------------------------------------------------------------------
       ! bigg mechanism heterogeneous freezing on existing cloud nuclei
       ! -----------------------------------------------------------------------
       tc = tice - tz (k)
       if (do_bigg .and. ql (k) > qcmin .and. tc > 0.) then
          sink = fac_frz * (100.0/rhor/ccn(k)) * dts * (exp (0.66 * tc) - 1.) * den (k) * ql (k) * ql (k)
          sink = min (ql (k), tc / icpk (k), sink)
          ql (k) = ql (k) - sink
          qi (k) = qi (k) + sink
          q_liq (k) = q_liq (k) - sink
          q_sol (k) = q_sol (k) + sink
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) + sink * lhi (k) / cvm (k)
       endif ! significant ql existed

       ! -----------------------------------------------------------------------
       ! update capacity heat and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhl (k) = lv00 + d0_vap * tz (k)
       lhi (k) = li00 + dc_ice * tz (k)
       lcpk (k) = lhl (k) / cvm (k)
       icpk (k) = lhi (k) / cvm (k)
       tcpk (k) = lcpk (k) + icpk (k)

       ! -----------------------------------------------------------------------
       ! sublimation / deposition of LS ice
       ! -----------------------------------------------------------------------

       if (tz (k) < tice) then
          qsi = iqs2 (tz (k), den (k), dqsdt)
          dq = (qv (k) - qsi)
          sink = min(qi(k), dq / (1. + tcpk (k) * dqsdt))
          if (qi (k) > qcmin) then
             ! eq 9, hong et al. 2004, mwr
             ! for a and b, see dudhia 1989: page 3103 eq (b7) and (b8)
             pidep = dts * dq * 349138.78 * exp (0.875 * log (qi (k) * den (k))) &
                  / (qsi * den (k) * lat2 / (0.0243 * rvgas * tz (k) ** 2) + 4.42478e4)
          else
             pidep = 0.
          endif
          if (dq > 0.) then ! vapor - > ice
             ! deposition
             ifrac = ice_fraction(tz (k),cnv_fraction,srf_type)
             tmp = tice - tz (k)
             qi_crt = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
             qi_crt = max (qi_crt, 1.82e-6) * qi_lim * ifrac / den (k)
             sink = min (sink, max (qi_crt - qi (k), pidep), tmp / tcpk (k))
          else ! ice -- > vapor
             ! sublimation
             if (do_subl) then
                pidep = pidep * min (1., dim (tz (k), t_sub) * 0.2)
                sink = fac_i2v * max (pidep, sink, - qi (k))
             else
                sink = 0.
             endif
          endif
          qv (k) = qv (k) - sink
          qi (k) = qi (k) + sink
          q_sol (k) = q_sol (k) + sink
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) + sink * (lhl (k) + lhi (k)) / cvm (k)
       endif

       ! -----------------------------------------------------------------------
       ! update capacity heat and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhl (k) = lv00 + d0_vap * tz (k)
       lhi (k) = li00 + dc_ice * tz (k)
       lcpk (k) = lhl (k) / cvm (k)
       icpk (k) = lhi (k) / cvm (k)
       tcpk (k) = lcpk (k) + icpk (k)

       ! -----------------------------------------------------------------------
       ! sublimation / deposition of snow
       ! this process happens for all temp rage
       ! -----------------------------------------------------------------------

       if (qs (k) > qpmin) then
          qsi = iqs2 (tz (k), den (k), dqsdt)
          qden = qs (k) * den (k)
          tmp = exp (0.65625 * log (qden))
          tsq = tz (k) * tz (k)
          dq = (qsi - qv (k)) / (1. + tcpk (k) * dqsdt)
          pssub = cssub (1) * tsq * (cssub (2) * sqrt (qden) + cssub (3) * tmp * &
               sqrt (denfac (k))) / (cssub (4) * tsq + cssub (5) * qsi * den (k))
          pssub = (qsi - qv (k)) * dts * pssub
          if (pssub > 0.) then ! qs -- > qv, sublimation
             pssub = min (fac_s2v * pssub * min (1., dim (tz (k), t_sub) * 0.2), qs (k))
             subl1(k) = pssub / dts
          else
             if (tz (k) > tice) then
                pssub = 0. ! no deposition
             else
                pssub = max (fac_v2s * pssub, dq, (tz (k) - tice) / tcpk (k))
             endif
          endif
          qs (k) = qs (k) - pssub
          qv (k) = qv (k) + pssub
          q_sol (k) = q_sol (k) - pssub
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) - pssub * (lhl (k) + lhi (k)) / cvm (k)
       endif

       ! -----------------------------------------------------------------------
       ! update capacity heat and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhl (k) = lv00 + d0_vap * tz (k)
       lhi (k) = li00 + dc_ice * tz (k)
       lcpk (k) = lhl (k) / cvm (k)
       icpk (k) = lhi (k) / cvm (k)
       tcpk (k) = lcpk (k) + icpk (k)

       ! -----------------------------------------------------------------------
       ! simplified 2 - way grapuel sublimation - deposition mechanism
       ! -----------------------------------------------------------------------

       if (qg (k) > qpmin) then
          qsi = iqs2 (tz (k), den (k), dqsdt)
          dq = (qv (k) - qsi) / (1. + tcpk (k) * dqsdt)
          pgsub = (qv (k) / qsi - 1.) * qg (k)
          if (pgsub > 0.) then ! deposition
             if (tz (k) > tice) then
                pgsub = 0. ! no deposition
             else
                pgsub = min (fac_v2g * pgsub, 0.2 * dq, ql (k) + qr (k), &
                     (tice - tz (k)) / tcpk (k))
             endif
          else ! submilation
             pgsub = max (fac_g2v * pgsub, dq) * min (1., dim (tz (k), t_sub) * 0.1)
             subl1(k) = subl1(k) + pgsub / dts
          endif
          qg (k) = qg (k) + pgsub
          qv (k) = qv (k) - pgsub
          q_sol (k) = q_sol (k) + pgsub
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) + pgsub * (lhl (k) + lhi (k)) / cvm (k)
       endif

#ifdef USE_MIN_EVAP
       ! -----------------------------------------------------------------------
       ! update capacity heat and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhl (k) = lv00 + d0_vap * tz (k)
       lcpk (k) = lhl (k) / cvm (k)

       ! -----------------------------------------------------------------------
       ! * minimum evap of rain in dry environmental air
       ! -----------------------------------------------------------------------

       if (qr (k) > qpmin) then
          qsw = wqs2 (tz (k), den (k), dqsdt)
          sink = min (qr (k), dim (rh_rain * qsw, qv (k)) / (1. + lcpk (k) * dqsdt))
          qv (k) = qv (k) + sink
          qr (k) = qr (k) - sink
          q_liq (k) = q_liq (k) - sink
          cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
          tz (k) = tz (k) - sink * lhl (k) / cvm (k)
       endif
#endif

       ! -----------------------------------------------------------------------
       ! update capacity heat and latend heat coefficient
       ! -----------------------------------------------------------------------

       lhl (k) = lv00 + d0_vap * tz (k)
       cvm (k) = c_air + (qv (k) + q_liq (k) + q_sol (k)) * c_vap
       lcpk (k) = lhl (k) / cvm (k)

       ! -----------------------------------------------------------------------
       ! compute cloud fraction
       ! -----------------------------------------------------------------------
       if (.not. do_qa) cycle

       ! -----------------------------------------------------------------------
       ! combine water species
       ! -----------------------------------------------------------------------
       if (preciprad) then
          q_sol (k) = qi (k) + qs (k) + qg (k)
          q_liq (k) = ql (k) + qr (k)
       else
          q_sol (k) = qi (k)
          q_liq (k) = ql (k)
       endif
       q_cond (k) = q_liq (k) + q_sol (k)

       qpz = qv (k) + q_cond (k) ! qpz is conserved

       ! -----------------------------------------------------------------------
       ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
       ! -----------------------------------------------------------------------

       tin = tz (k) - (lcpk (k) * q_cond (k) + icpk (k) * q_sol (k)) ! minimum temperature
       ! tin = tz (k) - ((lv00 + d0_vap * tz (k)) * q_cond (k) + &
       ! (li00 + dc_ice * tz (k)) * q_sol (k)) / (c_air + qpz * c_vap)

       ! -----------------------------------------------------------------------
       ! determine saturated specific humidity
       ! -----------------------------------------------------------------------

       if (tin <= t_wfr) then
          ! ice phase:
          qstar = iqs1 (tin, den (k))
       elseif (tin >= tice) then
          ! liquid phase:
          qstar = wqs1 (tin, den (k))
       else
          ! mixed phase:
          qsi = iqs1 (tin, den (k))
          qsw = wqs1 (tin, den (k))
          if (q_cond (k) > 3.e-6) then
             rqi = q_sol (k) / q_cond (k)
          else
             ! WMP impose CALIPSO ice polynomial from 0 C to -40 C
             rqi = ice_fraction(tin,cnv_fraction,srf_type)
          endif
          qstar = rqi * qsi + (1. - rqi) * qsw
       endif

       ! -----------------------------------------------------------------------
       ! assuming subgrid linear distribution in horizontal; this is effectively a smoother for the
       ! binary cloud scheme
       ! -----------------------------------------------------------------------
       if (qpz > qcmin) then
          ! partial cloudiness by pdf:
          dq = max (qcmin, h_var(k) * qpz)
          q_plus = qpz + dq ! cloud free if qstar > q_plus
          q_minus = qpz - dq
          if (icloud_f == 3) then
             ! triangular
             if(q_plus.le.qstar) then
                ! little/no cloud cover
             elseif ( (qpz.le.qstar).and.(qstar.lt.q_plus) ) then ! partial cloud cover
                qa (k) = max(qcmin, min(1., qa (k) + (q_plus-qstar)*(q_plus-qstar) / ( (q_plus-q_minus)*(q_plus-qpz) )))
             elseif ( (q_minus.le.qstar).and.(qstar.lt.qpz) ) then ! partial cloud cover
                qa (k) = max(qcmin, min(1., qa (k) + 1. - ( (qstar-q_minus)*(qstar-q_minus) / ( (q_plus-q_minus)*(qpz-q_minus) ))))
             elseif ( qstar.le.q_minus ) then
                qa (k) = 1.0 ! air fully saturated; 100 % cloud cover
             endif
          else
             ! top-hat
             if(q_plus.le.qstar) then
                ! little/no cloud cover
             elseif (qstar < q_plus .and. q_cond (k) > qc_crt) then
                qa (k) = max(qcmin, min(1., qa (k) + (q_plus - qstar) / (dq + dq) )) ! partial cloud cover
             elseif (qstar .le. q_minus) then
                qa (k) = 1.0 ! air fully saturated; 100 % cloud cover
             endif
          endif
       endif

    enddo

  end subroutine subgrid_z_proc

  ! =======================================================================
  !>@brief The subroutine 'terminal_fall' computes terminal fall speed.
  !>@details It considers cloud ice, snow, and graupel's melting during fall.
  ! =======================================================================

  subroutine terminal_fall_3d ( &
       dtm, is, ie, js, je, ktop, kbot, & ! in
       tz, qv, ql, qr, qg, qs, qi, & !in/out
       dz, dp, den, vtg, vts, vti, & ! in
       r1, g1, s1, i1, & ! out
       m1_sol, w1) ! in/out

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in) :: dtm ! time step (s)

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: vtg, vts, vti, den, dp, dz

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: qv, ql, qr, qg, qs, qi, tz, m1_sol, w1

    real, intent (out), dimension(is:ie, js:je) :: r1, g1, s1, i1

    real, dimension (is:ie, js:je, ktop:kbot + 1) :: ze, zt

    real :: qsat, dqsdt, evap, dtime
    real :: factor, frac
    real :: tmp, precip, tc, sink

    real, dimension (is:ie, js:je, ktop:kbot) :: icpk, cvm
    real, dimension (is:ie, js:je, ktop:kbot) :: m1, dm

    real :: q_liq, q_sol, lcpk, lhl, lhi

    real :: zs
    real :: fac_imlt

    integer :: i, j, k, m
    integer, dimension (is:ie, js:je) :: k0

    logical :: exit_flag
    logical, dimension (is:ie, js:je) :: no_fall

    zs = 0.

    fac_imlt = 1. - exp (- dtm / tau_imlt)

    !$omp target data &
    !$omp   map(to: vtg, vts, vti, den, dp, dz) &
    !$omp   map(tofrom: qv, ql, qr, qg, qs, qi, tz, m1_sol, w1) &
    !$omp   map(from: r1, g1, s1, i1) &
    !$omp   map(alloc: ze, zt, icpk, cvm, m1, dm, k0, no_fall)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3) private(lhi, q_liq, q_sol)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             m1_sol (i, j, k) = 0.
             ! lhl (k) = lv00 + d0_vap * tz (k)
             lhi = li00 + dc_ice * tz (i, j, k)
             q_liq = ql (i, j, k) + qr (i, j, k)
             q_sol = qi (i, j, k) + qs (i, j, k) + qg (i, j, k)
             cvm (i, j, k) = c_air + qv (i, j, k) * c_vap + q_liq * c_liq + q_sol * c_ice
             ! lcpk = lhl (k) / cvm (k)
             icpk (i, j, k) = lhi / cvm (i, j, k)
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! find significant melting level
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          k0 (i, j) = kbot
       end do
    end do
    !$omp end target teams distribute parallel do simd

    do k = ktop, kbot - 1
       do i = is, ie
          do j = js, je
             if (tz (i, j, k) > tice .and. k0 (i, j) == kbot) then
                k0 (i, j) = k
             endif
          end do
       end do
    enddo
    !$omp target update to(k0)

    ! -----------------------------------------------------------------------
    ! melting of cloud_ice (before fall) :
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3) private(tc, q_liq, q_sol, lhi, sink, tmp)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             if (k > k0 (i, j)) then
                tc = tz (i, j, k) - tice
                if (qi (i, j, k) > qcmin .and. tc > 0.) then
                   q_liq = ql (i, j, k) + qr (i, j, k)
                   q_sol = qi (i, j, k) + qs (i, j, k) + qg (i, j, k)
                   lhi = li00 + dc_ice * tz (i, j, k)
                   sink = min (qi (i, j, k), fac_imlt * tc / icpk (i, j, k))
                   tmp = min (sink, dim (ql_mlt, ql (i, j, k)))
                   ql (i, j, k) = ql (i, j, k) + tmp
                   qr (i, j, k) = qr (i, j, k) + sink - tmp
                   qi (i, j, k) = qi (i, j, k) - sink
                   q_liq = q_liq + sink
                   q_sol = q_sol - sink

                   cvm (i, j, k) = c_air + qv (i, j, k) * c_vap + q_liq * c_liq + q_sol * c_ice
                   tz (i, j, k) = tz (i, j, k) - sink * lhi / cvm (i, j, k)
                   tc = tz (i, j, k) - tice
                end if
             end if
          end do
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! turn off melting when cloud microphysics time step is small
    ! -----------------------------------------------------------------------

    ! TODO: Do we really need this block? k0 is being set in the next block for all dtm!
    if (dtm < 60.) then
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             k0 (i, j) = kbot
          end do
       end do
       !$omp end target teams distribute parallel do simd
    end if

    ! sjl, turn off melting of falling cloud ice, snow and graupel
    !$omp target teams distribute parallel do collapse(2)
    do i = is, ie
       do j = js, je
          k0 (i, j) = kbot
       end do
    end do
    !$omp end target teams distribute parallel do simd
    ! sjl, turn off melting of falling cloud ice, snow and graupel

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          ze (i, j, kbot + 1) = zs
       end do
    end do
    !$omp end target teams distribute parallel do simd
    ! r1 = sum(ze, dim=3)

    ! TODO: Figure out how to do this
    do k = kbot, ktop, - 1
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             ze (i, j, k) = ze (i, j, k + 1) - dz (i, j, k) ! dz < 0
          end do
       end do
       !$omp end target teams distribute parallel do simd
    enddo

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          zt (i, j, ktop) = ze (i, j, ktop)
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! update capacity heat and latend heat coefficient
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do collapse(3) private(lhi)
    ! do k = k0, kbot
    do k = kbot, kbot ! TODO: Is this correct? Should it be k0 to kbot?
       do i = is, ie
          do j = js, je
             lhi = li00 + dc_ice * tz (i, j, k)
             icpk (i, j, k) = lhi / cvm (i, j, k)
          end do
       end do
    enddo

    ! -----------------------------------------------------------------------
    ! melting of falling cloud ice into rain
    ! -----------------------------------------------------------------------

    call check_column_3d (is, ie, js, je, ktop, kbot, qi, no_fall)

    !$omp target teams distribute collapse(2)
    do i = is, ie
       do j = js, je

          if (vi_fac < 1.e-5 .or. no_fall (i, j)) then

             i1 (i, j) = 0.

          else

             !$omp parallel do simd
             do k = ktop + 1, kbot
                zt (i, j, k) = ze (i, j, k) - dtm * (vti (i, j, k - 1) + vti (i, j, k))/2.0
             enddo
             !$omp end parallel do simd

             zt (i, j, kbot + 1) = zs - dtm * vti (i, j, kbot)

             ! !$omp ordered
             do k = ktop, kbot
                if (zt (i, j, k + 1) >= zt (i, j, k)) then
                   zt (i, j, k + 1) = zt (i, j, k) - dz_min
                end if
             enddo
             ! !$omp end ordered

             if (k0 (i, j) < kbot) then
                ! !$omp ordered
                do k = kbot - 1, k0 (i, j), - 1
                   if (qi (i, j, k) > qcmin) then
                      exit_flag = .true.
                      !$omp parallel do simd
                      do m = k + 1, kbot
                         if (zt (i, j, k + 1) >= ze (i, j, m) .and. exit_flag) exit_flag = .false.
                         if (zt (i, j, k) < ze (i, j, m + 1) .and. tz (i, j, m) > tice .and. exit_flag) then
                            dtime = min (1.0, (ze (i, j, m) - ze (i, j, m + 1)) / (max (vr_min, vti (i, j, k)) * tau_imlt))
                            sink = min (qi (i, j, k) * dp (i, j, k) / dp (i, j, m), dtime * (tz (i, j, m) - tice) / icpk (i, j, m))
                            tmp = min (sink, dim (ql_mlt, ql (i, j, m)))
                            ql (i, j, m) = ql (i, j, m) + tmp
                            qr (i, j, m) = qr (i, j, m) - tmp + sink
                            tz (i, j, m) = tz (i, j, m) - sink * icpk (i, j, m)
                            qi (i, j, k) = qi (i, j, k) - sink * dp (i, j, m) / dp (i, j, k)
                         endif
                      enddo
                      !$omp end parallel do
                   endif
                enddo
                ! !$omp end ordered
             endif

             if (do_sedi_w) then
                !$omp parallel do
                do k = ktop, kbot
                   dm (i, j, k) = &
                        dp (i, j, k) * ( &
                        1. + &
                        qv (i, j, k) + &
                        ql (i, j, k) + &
                        qr (i, j, k) + &
                        qi (i, j, k) + &
                        qs (i, j, k) + &
                        qg (i, j, k) &
                        )
                enddo
             endif ! do_sedi_w

          end if ! vi_fac

       end do
    end do
    !$omp end target teams distribute

    if (use_ppm) then
       !$omp target teams distribute parallel do collapse(2)
       do i = is, ie
          do j = js, je
             if (.not. no_fall(i, j)) then
                call lagrangian_fall_ppm ( &
                     ktop, kbot, zs, ze (i, j, :), zt (i, j, :), dp (i, j, :), &
                     qi (i, j, :), i1 (i, j), m1_sol (i, j, :), mono_prof)
             end if
          end do
       end do
    else
       call implicit_fall_3d (dtm, is, ie, js, je, ktop, kbot, ze, vti, dp, qi, i1, m1_sol, no_fall)
    endif

    if (do_sedi_w) then
       !$omp target teams distribute collapse(2)
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                w1 (i, j, ktop) = &
                     ( &
                     dm (i, j, ktop) * w1 (i, j, ktop) + &
                     m1_sol (i, j, ktop) * vti (i, j, ktop) &
                     ) / (dm (i, j, ktop) - m1_sol (i, j, ktop))
                !$omp parallel do simd
                do k = ktop + 1, kbot
                   w1 (i, j, k) = &
                        ( &
                        dm (i, j, k) * w1 (i, j, k) - &
                        m1_sol (i, j, k - 1) * vti (i, j, k - 1) + &
                        m1_sol (i, j, k) * vti (i, j, k) &
                        ) / (dm (i, j, k) + m1_sol (i, j, k - 1) - m1_sol (i, j, k))
                enddo
                !$omp end parallel do simd
             end if ! no_fall
          end do
       end do
       !$omp end target teams distribute
    endif

    ! -----------------------------------------------------------------------
    ! melting of falling snow into rain
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          r1 (i, j) = 0.
       end do
    end do
    !$omp end target teams distribute parallel do simd

    call check_column_3d (is, ie, js, je, ktop, kbot, qs, no_fall)

    !$omp target teams distribute collapse(2)
    do i = is, ie
       do j = js, je

          if (no_fall (i, j)) then

             s1 (i, j) = 0.

          else

             !$omp parallel do
             do k = ktop + 1, kbot
                zt (i, j, k) = ze (i, j, k) - dtm * (vts (i, j, k - 1) + vts (i, j, k))/2.0
             enddo
             !$omp end parallel do
             zt (i, j, kbot + 1) = zs - dtm * vts (i, j, kbot)

             ! !$omp ordered
             do k = ktop, kbot
                if (zt (i, j, k + 1) >= zt (i, j, k)) then
                   zt (i, j, k + 1) = zt (i, j, k) - dz_min
                end if
             enddo
             ! !$omp end ordered

             if (k0 (i, j) < kbot) then
                ! !$omp ordered
                do k = kbot - 1, k0 (i, j), - 1
                   if (qs (i, j, k) > qpmin) then
                      exit_flag = .true.
                      ! !$omp ordered
                      do m = k + 1, kbot
                         if (zt (i, j, k + 1) >= ze (i, j, m) .and. exit_flag) exit_flag = .false.
                         if (exit_flag) then
                            dtime = min (dtm, (ze (i, j, m) - ze (i, j, m + 1)) / (vr_min + vts (i, j, k)))
                            if (zt (i, j, k) < ze (i, j, m + 1) .and. tz (i, j, m) > tice) then
                               dtime = min (1.0, dtime / tau_smlt)
                               sink = min ( &
                                    qs (i, j, k) * dp (i, j, k) / dp (i, j, m), &
                                    dtime * (tz (i, j, m) - tice) / icpk (i, j, m) &
                                    )
                               tz (i, j, m) = tz (i, j, m) - sink * icpk (i, j, m)
                               qs (i, j, k) = qs (i, j, k) - sink * dp (i, j, m) / dp (i, j, k)
                               if (zt (i, j, k) < zs) then
                                  r1 (i, j) = r1 (i,j) + sink * dp (i, j, m) ! precip as rain
                               else
                                  ! qr source here will fall next time step (therefore, can evap)
                                  qr (i, j, m) = qr (i, j, m) + sink
                               endif
                            endif
                         endif ! exit_flag
                         if (qs (i, j, k) < qpmin .and. exit_flag) exit_flag = .false.
                      enddo
                      ! !$omp end ordered
                   endif ! qs (i, j, k) > qpmin
                enddo ! k
                ! !$omp end ordered
             endif ! k0 < kbot

             if (do_sedi_w) then
                !$omp parallel do
                do k = ktop, kbot
                   dm (i, j, k) = &
                        dp (i, j, k) * ( &
                        1. + &
                        qv (i, j, k) + &
                        ql (i, j, k) + &
                        qr (i, j, k) + &
                        qi (i, j, k) + &
                        qs (i, j, k) + &
                        qg (i, j, k) &
                        )
                enddo
                !$omp end parallel do
             endif ! do_sedi_w

          end if ! no_fall

       end do
    end do
    !$omp end target teams distribute

    if (use_ppm) then
       !$omp target teams distribute parallel do collapse(2)
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                call lagrangian_fall_ppm ( &
                     ktop, kbot, zs, ze (i, j, :), zt (i, j, :), dp (i, j, :), &
                     qs (i, j, :), s1 (i, j), m1 (i, j, :), mono_prof)
             end if
          end do
       end do
       !$omp end target teams distribute parallel do
    else
       call implicit_fall_3d (dtm, is, ie, js, je, ktop, kbot, ze, vts, dp, qs, s1, m1, no_fall)
    endif

    !$omp target teams distribute collapse(2)
    do i = is, ie
       do j = js, je

          if (.not. no_fall (i, j)) then

             !$omp parallel do
             do k = ktop, kbot
                m1_sol (i, j, k) = m1_sol (i, j, k) + m1 (i, j, k)
             enddo
             !$omp end parallel do

             if (do_sedi_w) then
                w1 (i, j, ktop) = &
                     ( &
                     dm (i, j, ktop) * w1 (i, j, ktop) + &
                     m1 (i, j, ktop) * vts (i, j, ktop) &
                     ) / (dm (i, j, ktop) - m1 (i, j, ktop))
                !$omp parallel do simd
                do k = ktop + 1, kbot
                   w1 (i, j, k) = &
                        ( &
                        dm (i, j, k) * w1 (i, j, k) - &
                        m1 (i, j, k - 1) * vts (i, j, k - 1) + &
                        m1 (i, j, k) * vts (i, j, k) &
                        ) / (dm (i, j, k) + m1 (i, j, k - 1) - m1 (i, j, k))
                end do
                !$omp end parallel do simd
             end if ! do_sedi_w

          end if ! no_fall

       end do
    end do
    !$omp end target teams distribute

    ! ----------------------------------------------
    ! melting of falling graupel into rain
    ! ----------------------------------------------

    call check_column_3d (is, ie, js, je, ktop, kbot, qg, no_fall)

    !$omp target teams distribute collapse(2)
    do i = is, ie
       do j = js, je

          if (no_fall (i, j)) then

             g1 (i, j) = 0.

          else

             !$omp parallel do simd
             do k = ktop + 1, kbot
                zt (i, j, k) = ze (i, j, k) - dtm * (vtg (i, j, k - 1) + vtg (i, j, k))/2.0
             enddo
             !$omp end parallel do simd
             zt (i, j, kbot + 1) = zs - dtm * vtg (i, j, kbot)

             ! !$omp ordered
             do k = ktop, kbot
                if (zt (i, j, k + 1) >= zt (i, j, k)) then
                   zt (i, j, k + 1) = zt (i, j, k) - dz_min
                end if
             enddo
             ! !$omp ordered

             if (k0 (i, j) < kbot) then
                ! !$omp ordered
                do k = kbot - 1, k0 (i, j), - 1
                   if (qg (i, j, k) > qpmin) then
                      exit_flag = .true.
                      ! !$omp ordered
                      do m = k + 1, kbot
                         if (zt (i, j, k + 1) >= ze (i, j, m) .and. exit_flag) exit_flag = .false.
                         if (exit_flag) then
                            dtime = min (dtm, (ze (i, j, m) - ze (i, j, m + 1)) / vtg (i, j, k))
                            if (zt (i, j, k) < ze (i, j, m + 1) .and. tz (i, j, m) > tice) then
                               dtime = min (1., dtime / tau_g2r)
                               sink = min ( &
                                    qg (i, j, k) * dp (i, j, k) / dp (i, j, m), &
                                    dtime * (tz (i, j, m) - tice) / icpk (i, j, m))
                               tz (i, j, m) = tz (i, j, m) - sink * icpk (i, j, m)
                               qg (i, j, k) = qg (i, j, k) - sink * dp (i, j, m) / dp (i, j, k)
                               if (zt (i, j, k) < zs) then
                                  r1 (i, j) = r1 (i, j) + sink * dp (i, j, m)
                               else
                                  qr (i, j, m) = qr (i, j, m) + sink
                               endif
                            endif
                         endif ! exit_flag
                         if (qg (i, j, k) < qpmin .and. exit_flag) exit_flag = .false.
                      enddo ! m
                      ! !$omp end ordered
                   endif ! qg > qpmin
                enddo
                ! !$omp end ordered
             endif ! k0 < kbot

             if (do_sedi_w) then
                !$omp parallel do simd
                do k = ktop, kbot
                   dm (i, j, k) = dp (i, j, k) * ( &
                        1. + &
                        qv (i, j, k) + &
                        ql (i, j, k) + &
                        qr (i, j, k) + &
                        qi (i, j, k) + &
                        qs (i, j, k) + &
                        qg (i, j, k))
                end do
                !$omp end parallel do simd
             end if ! do_sedi_w

          end if ! no_fall

       end do
    end do
    !$omp end target teams distribute

    if (use_ppm) then
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                call lagrangian_fall_ppm ( &
                     ktop, kbot, zs, ze (i, j, :), zt (i, j, :), dp (i, j, :), &
                     qg (i, j, :), g1 (i, j), m1 (i, j, :), mono_prof)
             end if
          end do
       end do
       !$omp end target teams distribute parallel do simd
    else
       call implicit_fall_3d (dtm, is, ie, js, je, ktop, kbot, ze, vtg, dp, qg, g1, m1, no_fall)
    endif

    !$omp target teams distribute collapse(2)
    do i = is, ie
       do j = js, je

          if (.not. no_fall (i, j)) then
             !$omp parallel do simd
             do k = ktop, kbot
                m1_sol (i, j, k) = m1_sol (i, j, k) + m1 (i, j, k)
             enddo
             !$omp end parallel do simd

             if (do_sedi_w) then
                w1 (i, j, ktop) = ( &
                     dm (i, j, ktop) * w1 (i, j, ktop) + &
                     m1 (i, j, ktop) * vtg (i, j, ktop) &
                     ) / (dm (i, j, ktop) - m1 (i, j, ktop))
                !$omp parallel do
                do k = ktop + 1, kbot
                   w1 (i, j, k) = ( &
                        dm (i, j, k) * w1 (i, j, k) - &
                        m1 (i, j, k - 1) * vtg (i, j, k - 1) + &
                        m1 (i, j, k) * vtg (i, j, k) &
                        ) / (dm (i, j, k) + m1 (i, j, k - 1) - m1 (i, j, k))
                enddo
                !$omp end parallel do
             endif ! do_sedi_w

          endif ! no_fall

       end do
    end do
    !$omp end target teams distribute

    !$omp end target data

  end subroutine terminal_fall_3d

  ! =======================================================================
  !>@brief The subroutine 'check_column' checks
  !!       if the water species is large enough to fall.
  ! =======================================================================

  subroutine check_column (ktop, kbot, q, no_fall)

    implicit none
    !$omp declare target

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: q (ktop:kbot)

    logical, intent (out) :: no_fall

    integer :: k

    no_fall = .true.

    do k = ktop, kbot
       if (q (k) > qpmin) then
          no_fall = .false.
          exit
       endif
    enddo

  end subroutine check_column

  ! =======================================================================
  !>@brief The subroutine 'check_column' checks
  !!       if the water species is large enough to fall.
  ! =======================================================================

  subroutine check_column_3d (is, ie, js, je, ktop, kbot, q, no_fall)

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in) :: q (is:ie, js:je, ktop:kbot)

    logical, intent (out) :: no_fall (is:ie, js:je)

    integer :: i, j, k

    !$omp target data map(to: q) map(from: no_fall)

    !$omp target teams distribute parallel do collapse(2)
    do i = is, ie
       do j = js, je
          no_fall (i, j) = .true.
       end do
    end do

    do k = ktop, kbot
       !$omp target teams distribute parallel do collapse(2)
       do i = is, ie
          do j = js, je
             if (q (i, j, k) > qpmin .and. no_fall(i, j)) then
                no_fall (i, j) = .false.
             end if
          end do
       end do
    end do

    !$omp end target data

  end subroutine check_column_3d

  ! =======================================================================
  !>@brief The subroutine 'implicit_fall' computes the time-implicit monotonic
  !! scheme.
  !>@author Shian-Jiann Lin, 2016
  ! =======================================================================

  subroutine implicit_fall_3d (dt, is, ie, js, je, ktop, kbot, ze, vt, dp, q, precip, m1, no_fall)

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in) :: dt

    real, intent (in), dimension (is:ie, js:je, ktop:kbot + 1) :: ze

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: vt, dp

    logical, intent (in), dimension (is:ie, js:je) :: no_fall

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: q

    real, intent (out), dimension (is:ie, js:je, ktop:kbot) :: m1

    real, intent (out), dimension (is:ie, js:je) :: precip

    real, dimension (is:ie, js:je, ktop:kbot) :: qm

    integer :: i, j, k

    !$omp target data &
    !$omp   map(to: ze, vt, dp, no_fall) &
    !$omp   map(tofrom: q) &
    !$omp   map(from: m1, precip) &
    !$omp   map(alloc: qm)

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i,j)) q (i, j, k) = q (i, j, k) * dp (i, j, k)
          end do
       end do
    enddo
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! sedimentation: non - vectorizable loop
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          if (.not. no_fall (i, j)) then
             qm (i, j, ktop) = q (i, j, ktop) / ((ze (i, j, ktop) - ze (i, j, ktop + 1)) + (dt * vt (i, j, ktop)))
          end if
       end do
    end do
    !$omp end target teams distribute parallel do simd

    do k = ktop + 1, kbot
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                qm (i, j, k) = &
                     (q (i, j, k) + (dt * vt (i, j, k-1)) * qm (i, j, k - 1)) / &
                     ((ze (i, j, k) - ze (i, j, k + 1)) + (dt * vt (i, j, k)))
             end if
          end do
       end do
       !$omp end target teams distribute parallel do simd
    enddo

    ! -----------------------------------------------------------------------
    ! qm is density at this stage
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                qm (i, j, k) = qm (i, j, k) * (ze (i, j, k) - ze (i, j, k + 1))
             end if
          end do
       end do
    enddo
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! output mass fluxes: non - vectorizable loop
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          if (.not. no_fall (i, j)) then
             m1 (i, j, ktop) = q (i, j, ktop) - qm (i, j, ktop)
          end if
       end do
    end do
    !$omp end target teams distribute parallel do simd

    do k = ktop + 1, kbot
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                m1 (i, j, k) = m1 (i, j, k - 1) + q (i, j, k) - qm (i, j, k)
             end if
          end do
       end do
       !$omp end target teams distribute parallel do simd
    enddo

    !$omp target teams distribute parallel do simd collapse(2)
    do i = is, ie
       do j = js, je
          if (.not. no_fall (i, j)) then
             precip (i, j) = m1 (i, j, kbot)
          end if
       end do
    end do
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! update:
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je
             if (.not. no_fall (i, j)) then
                q (i, j, k) = qm (i, j, k) / dp (i, j, k)
             end if
          end do
       end do
    enddo
    !$omp end target teams distribute parallel do simd

    !$omp end target data

  end subroutine implicit_fall_3d

  ! =======================================================================
  !> lagrangian scheme
  !  developed by sj lin, ????
  ! =======================================================================

  subroutine lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, q, precip, m1, mono)

    implicit none

    ! TODO: This routine has not been touched since we run with use_ppm = .false.
    !$omp declare target

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: zs

    logical, intent (in) :: mono

    real, intent (in), dimension (ktop:kbot + 1) :: ze, zt

    real, intent (in), dimension (ktop:kbot) :: dp

    ! m1: flux
    real, intent (inout), dimension (ktop:kbot) :: q, m1

    real, intent (out) :: precip

    real, dimension (ktop:kbot) :: qm, dz

    real :: a4 (4, ktop:kbot)

    real :: pl, pr, delz, esl

    integer :: k, k0, n, m

    real, parameter :: r3 = 1. / 3., r23 = 2. / 3.

    ! -----------------------------------------------------------------------
    ! density:
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
       dz (k) = zt (k) - zt (k + 1) ! note: dz is positive
       q (k) = q (k) * dp (k)
       a4 (1, k) = q (k) / dz (k)
       qm (k) = 0.
    enddo

    ! -----------------------------------------------------------------------
    ! construct vertical profile with zt as coordinate
    ! -----------------------------------------------------------------------

    call cs_profile (a4 (1, ktop), dz (ktop), kbot - ktop + 1, mono)

    k0 = ktop
    do k = ktop, kbot
       do n = k0, kbot
          if (ze (k) <= zt (n) .and. ze (k) >= zt (n + 1)) then
             pl = (zt (n) - ze (k)) / dz (n)
             if (zt (n + 1) <= ze (k + 1)) then
                ! entire new grid is within the original grid
                pr = (zt (n) - ze (k + 1)) / dz (n)
                qm (k) = a4 (2, n) + 0.5 * (a4 (4, n) + a4 (3, n) - a4 (2, n)) * (pr + pl) - &
                     a4 (4, n) * r3 * (pr * (pr + pl) + pl ** 2)
                qm (k) = qm (k) * (ze (k) - ze (k + 1))
                k0 = n
                goto 555
             else
                qm (k) = (ze (k) - zt (n + 1)) * (a4 (2, n) + 0.5 * (a4 (4, n) + &
                     a4 (3, n) - a4 (2, n)) * (1. + pl) - a4 (4, n) * (r3 * (1. + pl * (1. + pl))))
                if (n < kbot) then
                   do m = n + 1, kbot
                      ! locate the bottom edge: ze (k + 1)
                      if (ze (k + 1) < zt (m + 1)) then
                         qm (k) = qm (k) + q (m)
                      else
                         delz = zt (m) - ze (k + 1)
                         esl = delz / dz (m)
                         qm (k) = qm (k) + delz * (a4 (2, m) + 0.5 * esl * &
                              (a4 (3, m) - a4 (2, m) + a4 (4, m) * (1. - r23 * esl)))
                         k0 = m
                         goto 555
                      endif
                   enddo
                endif
                goto 555
             endif
          endif
       enddo
555    continue
    enddo

    m1 (ktop) = q (ktop) - qm (ktop)
    do k = ktop + 1, kbot
       m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (kbot)

    ! convert back to * dry * mixing ratio:
    ! dp must be dry air_mass (because moist air mass will be changed due to terminal fall) .

    do k = ktop, kbot
       q (k) = qm (k) / dp (k)
    enddo

  end subroutine lagrangian_fall_ppm

  subroutine cs_profile (a4, del, km, do_mono)

    implicit none
    !$omp declare target

    integer, intent (in) :: km !< vertical dimension

    real, intent (in) :: del (km)

    logical, intent (in) :: do_mono

    real, intent (inout) :: a4 (4, km)

    real, parameter :: qp_min = 1.e-6

    real :: gam (km)
    real :: q (km + 1)
    real :: d4, bet, a_bot, grat, pmp, lac
    real :: pmp_1, lac_1, pmp_2, lac_2
    real :: da1, da2, a6da

    integer :: k

    logical extm (km)

    grat = del (2) / del (1) ! grid ratio
    bet = grat * (grat + 0.5)
    q (1) = (2. * grat * (grat + 1.) * a4 (1, 1) + a4 (1, 2)) / bet
    gam (1) = (1. + grat * (grat + 1.5)) / bet

    do k = 2, km
       d4 = del (k - 1) / del (k)
       bet = 2. + 2. * d4 - gam (k - 1)
       q (k) = (3. * (a4 (1, k - 1) + d4 * a4 (1, k)) - q (k - 1)) / bet
       gam (k) = d4 / bet
    enddo

    a_bot = 1. + d4 * (d4 + 1.5)
    q (km + 1) = (2. * d4 * (d4 + 1.) * a4 (1, km) + a4 (1, km - 1) - a_bot * q (km)) &
         / (d4 * (d4 + 0.5) - a_bot * gam (km))

    do k = km, 1, - 1
       q (k) = q (k) - gam (k) * q (k + 1)
    enddo

    ! -----------------------------------------------------------------------
    ! apply constraints
    ! -----------------------------------------------------------------------

    do k = 2, km
       gam (k) = a4 (1, k) - a4 (1, k - 1)
    enddo

    ! -----------------------------------------------------------------------
    ! apply large - scale constraints to all fields if not local max / min
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! top:
    ! -----------------------------------------------------------------------

    q (1) = max (q (1), 0.)
    q (2) = min (q (2), max (a4 (1, 1), a4 (1, 2)))
    q (2) = max (q (2), min (a4 (1, 1), a4 (1, 2)), 0.)

    ! -----------------------------------------------------------------------
    ! interior:
    ! -----------------------------------------------------------------------

    do k = 3, km - 1
       if (gam (k - 1) * gam (k + 1) > 0.) then
          q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
          q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
       else
          if (gam (k - 1) > 0.) then
             ! there exists a local max
             q (k) = max (q (k), min (a4 (1, k - 1), a4 (1, k)))
          else
             ! there exists a local min
             q (k) = min (q (k), max (a4 (1, k - 1), a4 (1, k)))
             q (k) = max (q (k), 0.0)
          endif
       endif
    enddo

    ! -----------------------------------------------------------------------
    ! bottom :
    ! -----------------------------------------------------------------------

    q (km) = min (q (km), max (a4 (1, km - 1), a4 (1, km)))
    q (km) = max (q (km), min (a4 (1, km - 1), a4 (1, km)), 0.)
    ! q (km + 1) = max (q (km + 1), 0.)

    ! -----------------------------------------------------------------------
    ! f (s) = al + s * [ (ar - al) + a6 * (1 - s) ] (0 <= s <= 1)
    ! -----------------------------------------------------------------------

    do k = 1, km - 1
       a4 (2, k) = q (k)
       a4 (3, k) = q (k + 1)
    enddo

    do k = 2, km - 1
       if (gam (k) * gam (k + 1) > 0.0) then
          extm (k) = .false.
       else
          extm (k) = .true.
       endif
    enddo

    if (do_mono) then
       do k = 3, km - 2
          if (extm (k)) then
             ! positive definite constraint only if true local extrema
             if (a4 (1, k) < qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                a4 (2, k) = a4 (1, k)
                a4 (3, k) = a4 (1, k)
             endif
          else
             a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
             if (abs (a4 (4, k)) > abs (a4 (2, k) - a4 (3, k))) then
                ! check within the smooth region if subgrid profile is non - monotonic
                pmp_1 = a4 (1, k) - 2.0 * gam (k + 1)
                lac_1 = pmp_1 + 1.5 * gam (k + 2)
                a4 (2, k) = min (max (a4 (2, k), min (a4 (1, k), pmp_1, lac_1)), &
                     max (a4 (1, k), pmp_1, lac_1))
                pmp_2 = a4 (1, k) + 2.0 * gam (k)
                lac_2 = pmp_2 - 1.5 * gam (k - 1)
                a4 (3, k) = min (max (a4 (3, k), min (a4 (1, k), pmp_2, lac_2)), &
                     max (a4 (1, k), pmp_2, lac_2))
             endif
          endif
       enddo
    else
       do k = 3, km - 2
          if (extm (k)) then
             if (a4 (1, k) < qp_min .or. extm (k - 1) .or. extm (k + 1)) then
                a4 (2, k) = a4 (1, k)
                a4 (3, k) = a4 (1, k)
             endif
          endif
       enddo
    endif

    do k = 1, km - 1
       a4 (4, k) = 6. * a4 (1, k) - 3. * (a4 (2, k) + a4 (3, k))
    enddo

    k = km - 1
    if (extm (k)) then
       a4 (2, k) = a4 (1, k)
       a4 (3, k) = a4 (1, k)
       a4 (4, k) = 0.
    else
       da1 = a4 (3, k) - a4 (2, k)
       da2 = da1 ** 2
       a6da = a4 (4, k) * da1
       if (a6da < - da2) then
          a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
          a4 (3, k) = a4 (2, k) - a4 (4, k)
       elseif (a6da > da2) then
          a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
          a4 (2, k) = a4 (3, k) - a4 (4, k)
       endif
    endif

    call cs_limiters (km - 1, a4)

    ! -----------------------------------------------------------------------
    ! bottom layer:
    ! -----------------------------------------------------------------------

    a4 (2, km) = a4 (1, km)
    a4 (3, km) = a4 (1, km)
    a4 (4, km) = 0.

  end subroutine cs_profile

  subroutine cs_limiters (km, a4)

    implicit none
    !$omp declare target

    integer, intent (in) :: km

    real, intent (inout) :: a4 (4, km) !< ppm array

    real, parameter :: r12 = 1. / 12.

    integer :: k

    ! -----------------------------------------------------------------------
    ! positive definite constraint
    ! -----------------------------------------------------------------------

    do k = 1, km
       if (abs (a4 (3, k) - a4 (2, k)) < - a4 (4, k)) then
          if ((a4 (1, k) + 0.25 * (a4 (3, k) - a4 (2, k)) ** 2 / a4 (4, k) + a4 (4, k) * r12) < 0.) then
             if (a4 (1, k) < a4 (3, k) .and. a4 (1, k) < a4 (2, k)) then
                a4 (3, k) = a4 (1, k)
                a4 (2, k) = a4 (1, k)
                a4 (4, k) = 0.
             elseif (a4 (3, k) > a4 (2, k)) then
                a4 (4, k) = 3. * (a4 (2, k) - a4 (1, k))
                a4 (3, k) = a4 (2, k) - a4 (4, k)
             else
                a4 (4, k) = 3. * (a4 (3, k) - a4 (1, k))
                a4 (2, k) = a4 (3, k) - a4 (4, k)
             endif
          endif
       endif
    enddo

  end subroutine cs_limiters

  ! =======================================================================
  !>@brief The subroutine 'fall_speed' calculates vertical fall speed.
  ! =======================================================================

  subroutine fall_speed (ktop, kbot, pl, cnv_fraction, anv_icefall, lsc_icefall, &
       den, qs, qi, qg, ql, tk, vts, vti, vtg)

    implicit none
    !$omp declare target

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: cnv_fraction, anv_icefall, lsc_icefall
    real, intent (in) :: pl, den, qs, qi, qg, ql, tk
    real, intent (out) :: vts, vti, vtg

    ! fall velocity constants:

    real, parameter :: thi = 1.0e-8 !< cloud ice threshold for terminal fall
    real, parameter :: thg = 1.0e-8
    real, parameter :: ths = 1.0e-8

    real, parameter :: aaC = - 4.18334e-5
    real, parameter :: bbC = - 0.00525867
    real, parameter :: ccC = - 0.0486519
    real, parameter :: ddC = 0.00251197
    real, parameter :: eeC = 1.91523

    real, parameter :: aaL = - 1.70704e-5
    real, parameter :: bbL = - 0.00319109
    real, parameter :: ccL = - 0.0169876
    real, parameter :: ddL = 0.00410839
    real, parameter :: eeL = 1.93644

    ! marshall - palmer constants

    real, parameter :: vcons = 6.6280504
    real, parameter :: vcong = 87.2382675
    real, parameter :: norms = 942477796.076938
    real, parameter :: normg = 5026548245.74367

    real :: tc, rhof

    real :: vi1, viCNV, viLSC, IWC
    real :: rBB, C0, C1, DIAM, lnP
    integer :: k

    ! -----------------------------------------------------------------------
    ! marshall - palmer formula
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! try the local air density -- for global model; the true value could be
    ! much smaller than sfcrho over high mountains
    ! -----------------------------------------------------------------------

    rhof = sqrt (min (10., sfcrho / den))

    ! -----------------------------------------------------------------------
    ! ice:
    ! -----------------------------------------------------------------------

    if (const_vi) then
       vti = vi_fac
    else
       vi1 = 0.01 * vi_fac
       if (qi < thi) then ! this is needed as the fall - speed maybe problematic for small qi
          vti = vf_min
       else
          tc = tk - tice ! deg C
          IWC = qi * den * 1.e3 ! Units are g/m3
          ! -----------------------------------------------------------------------
          ! use deng and mace (2008, grl)
          ! https://doi.org/10.1029/2008GL035054
          ! -----------------------------------------------------------------------
          viLSC   = lsc_icefall*10.0**(log10(IWC) * (tc * (aaL * tc + bbL) + ccL) + ddL * tc + eeL)
          viCNV   = anv_icefall*10.0**(log10(IWC) * (tc * (aaC * tc + bbC) + ccC) + ddC * tc + eeC)
          ! -----------------------------------------------------------------------
          ! use Mishra et al (2014, JGR) 'Parameterization of ice fall speeds in
          !                               ice clouds: Results from SPartICus'
          ! -----------------------------------------------------------------------
          !viLSC  = MAX(10.0,lsc_icefall*(1.411*tc(k) + 11.71*log10(IWC*1.e3) + 82.35))
          !viCNV  = MAX(10.0,anv_icefall*(1.119*tc(k) + 14.21*log10(IWC*1.e3) + 68.85))
          ! Combine
          vti = viLSC*(1.0-cnv_fraction) + viCNV*(cnv_fraction)
          ! Update units from cm/s to m/s
          vti = vi1 * vti
          ! Include pressure sensitivity (eq 14 in https://doi.org/10.1175/JAS-D-12-0124.1)
          !------ice cloud effective radius ----- [klaus wyser, 1998]
          !if(tk(k)>t_ice) then
          !   rBB  = -2.
          !else
          !   rBB  = -2. + log10(IWC/50.)*(1.e-3*(t_ice-tk(k))**1.5)
          !endif
          !rBB   = MIN((MAX(rBB,-6.)),-2.)
          !DIAM  = 2.0*(377.4 + 203.3 * rBB+ 37.91 * rBB **2 + 2.3696 * rBB **3)
          !lnP   = log(pl(k)/100.0)
          !C0    = -1.04 + 0.298*lnP
          !C1    =  0.67 - 0.097*lnP
          ! apply pressure scaling
          !vti = vti * (C0 + C1*log(DIAM))
          ! Limits
          vti = min (vi_max, max (vf_min, vti))
       endif
    endif

    ! -----------------------------------------------------------------------
    ! snow:
    ! -----------------------------------------------------------------------

    if (const_vs) then
       vts = vs_fac ! 1. ifs_2016
    else
       if (qs < ths) then
          vts = vf_min
       else
          vts = vs_fac * vcons * rhof * exp (0.0625 * log (qs * den / norms))
          vts = min (vs_max, max (vf_min, vts))
       endif
    endif

    ! -----------------------------------------------------------------------
    ! graupel:
    ! -----------------------------------------------------------------------

    if (const_vg) then
       vtg = vg_fac ! 2.
    else
       if (qg < thg) then
          vtg = vf_min
       else
          vtg = vg_fac * vcong * rhof * sqrt (sqrt (sqrt (qg * den / normg)))
          vtg = min (vg_max, max (vf_min, vtg))
       endif
    endif

  end subroutine fall_speed

  ! =======================================================================
  !>@brief The subroutine 'setup'm' sets up
  !! gfdl cloud microphysics parameters.
  ! =======================================================================

  subroutine setupm

    implicit none

    real :: gcon, cd, scm3, pisq, act (8)
    real :: vdifu, tcond
    real :: visk
    real :: ch2o, hltf
    real :: hlts, hltc, ri50

    real, parameter :: gam263 = 1.456943, gam275 = 1.608355, gam290 = 1.827363, &
         gam325 = 2.54925, gam350 = 3.323363, gam380 = 4.694155, &
         gam425 = 8.285063, gam450 = 11.631769, gam480 = 17.837789, &
         gam625 = 184.860962, gam680 = 496.604067

    ! intercept parameters

    real, parameter :: rnzr = 8.0e6 ! lin83
    real, parameter :: rnzs = 3.0e6 ! lin83
    real, parameter :: rnzg = 4.0e6 ! rh84

    ! density parameters

    real, parameter :: rhos = 0.1e3 !< lin83 (snow density; 1 / 10 of water)
    real, parameter :: rhog = 0.4e3 !< rh84 (graupel density)
    real, parameter :: acc (3) = (/ 5.0, 2.0, 0.5 /)

    integer :: i, k

    pie = 4. * atan (1.0)

    vdifu = 2.11e-5
    tcond = 2.36e-2

    visk = 1.259e-5
    hlts = 2.8336e6
    hltc = 2.5e6
    hltf = 3.336e5

    ch2o = 4.1855e3
    ri50 = 1.e-4

    pisq = pie * pie
    scm3 = (visk / vdifu) ** (1. / 3.)

    cracs = pisq * rnzr * rnzs * rhos
    csacr = pisq * rnzr * rnzs * rhor
    cgacr = pisq * rnzr * rnzg * rhor
    cgacs = pisq * rnzg * rnzs * rhos
    cgacs = cgacs * c_pgacs

    ! act: 1 - 2:racs (s - r) ; 3 - 4:sacr (r - s) ;
    ! 5 - 6:gacr (r - g) ; 7 - 8:gacs (s - g)

    act (1) = pie * rnzs * rhos
    act (2) = pie * rnzr * rhor
    act (6) = pie * rnzg * rhog
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    act (7) = act (1)
    act (8) = act (6)

    do i = 1, 3
       do k = 1, 4
          acco (i, k) = acc (i) / (act (2 * k - 1) ** ((7 - i) * 0.25) * act (2 * k) ** (i * 0.25))
       enddo
    enddo

    gcon = 40.74 * sqrt (sfcrho) ! 44.628

    csacw = pie * rnzs * clin * gam325 / (4. * act (1) ** 0.8125)
    ! decreasing csacw to reduce cloud water --- > snow

    craci = pie * rnzr * alin * gam380 / (4. * act (2) ** 0.95)
    csaci = csacw * c_psaci

    cgacw = pie * rnzg * gam350 * gcon / (4. * act (6) ** 0.875)

    cgaci = cgacw * c_pgaci

    cracw = craci ! cracw = 3.27206196043822
    cracw = c_cracw * cracw

    ! subl and revp: five constants for three separate processes

    cssub (1) = 2. * pie * vdifu * tcond * rvgas * rnzs
    cgsub (1) = 2. * pie * vdifu * tcond * rvgas * rnzg
    crevp (1) = 2. * pie * vdifu * tcond * rvgas * rnzr
    cssub (2) = 0.78 / sqrt (act (1))
    cgsub (2) = 0.78 / sqrt (act (6))
    crevp (2) = 0.78 / sqrt (act (2))
    cssub (3) = 0.31 * scm3 * gam263 * sqrt (clin / visk) / act (1) ** 0.65625
    cgsub (3) = 0.31 * scm3 * gam275 * sqrt (gcon / visk) / act (6) ** 0.6875
    crevp (3) = 0.31 * scm3 * gam290 * sqrt (alin / visk) / act (2) ** 0.725
    cssub (4) = tcond * rvgas
    cssub (5) = hlts ** 2 * vdifu
    cgsub (4) = cssub (4)
    crevp (4) = cssub (4)
    cgsub (5) = cssub (5)
    crevp (5) = hltc ** 2 * vdifu

    cgfr (1) = 20.e2 * pisq * rnzr * rhor / act (2) ** 1.75
    cgfr (2) = 0.66

    ! smlt: five constants (lin et al. 1983)

    csmlt (1) = 2. * pie * tcond * rnzs / hltf
    csmlt (2) = 2. * pie * vdifu * rnzs * hltc / hltf
    csmlt (3) = cssub (2)
    csmlt (4) = cssub (3)
    csmlt (5) = ch2o / hltf

    ! gmlt: five constants

    cgmlt (1) = 2. * pie * tcond * rnzg / hltf
    cgmlt (2) = 2. * pie * vdifu * rnzg * hltc / hltf
    cgmlt (3) = cgsub (2)
    cgmlt (4) = cgsub (3)
    cgmlt (5) = ch2o / hltf

    es0 = 6.107799961e2 ! ~6.1 mb
    ces0 = eps * es0

    !$omp target update to( &
    !$omp     ces0, cracs, cracw, &
    !$omp     csaci, csacr, csacw, cgaci, cgacr, cgacs, cgacw, &
    !$omp     cssub(:), crevp(:), csmlt(:), cgmlt(:), cgfr(:), acco(:,:))

  end subroutine setupm

  ! =======================================================================
  ! initialization of gfdl cloud microphysics
  !>@brief The subroutine 'gfdl_cloud_microphys_init' initializes the GFDL
  !! cloud microphysics.
  ! =======================================================================

  subroutine gfdl_cloud_microphys_init ()
    implicit none
    integer :: file_handle
    character (len = 64) :: file_name = 'input-data/input.nml'

    integer :: rc
    logical :: exists

    inquire (file = trim (file_name), exist = exists)
    if (.not. exists) then
       write (6, *) 'gfdl - mp :: namelist file: ', trim (file_name), ' does not exist'
       stop
    else
       open(newunit = file_handle, file = file_name, status = 'old')
       read(nml = gfdl_cloud_microphysics_nml, unit = file_handle, iostat=rc)
       if (rc /=0) error stop "Could not read input namelist file"
       close(file_handle)
    endif

    !$omp target update to( &
    !$omp     tau_revp, tau_v2l, tau_l2v, tau_i2v, tau_s2v, tau_v2s, tau_g2v, &
    !$omp     tau_v2g, tau_frz, tau_imlt, tau_smlt, tau_i2s, tau_g2r, &
    !$omp     tice, tice0, rh_inc, rh_inr, t_min, do_qa, t_sub, do_evap, &
    !$omp     do_bigg, qi_lim, do_subl, preciprad, icloud_f, qc_crt, z_slope_ice, &
    !$omp     c_paut, prog_ccn, fix_negative, sedi_transport, ql_mlt, qs_mlt, qi0_crt, qs0_crt, &
    !$omp     const_vi, vi_fac, vi_max, const_vs, vs_fac, vs_max, const_vg, vg_fac, vg_max, const_vr, vr_fac, vr_max, &
    !$omp     use_ppm, mono_prof, rthreshs, rthreshu, irain_f, z_slope_liq, do_sedi_heat, &
    !$omp     ql0_max, dt_fr, sat_adj0, dw_land, dw_ocean, c_psaci, c_pgacs, &
    !$omp     ccn_l, ccn_o, c_cracw, use_ccn, de_ice, mp_time)

    if (do_setup) then
       call setup_con
       call setupm
       do_setup = .false.
    endif

    ! if (root_proc) write (logunit, nml = gfdl_cloud_microphys_nml)
    !
    ! id_vtr = register_diag_field (mod_name, 'vt_r', axes (1:3), time, &
    ! 'rain fall speed', 'm / s', missing_value = missing_value)
    ! id_vts = register_diag_field (mod_name, 'vt_s', axes (1:3), time, &
    ! 'snow fall speed', 'm / s', missing_value = missing_value)
    ! id_vtg = register_diag_field (mod_name, 'vt_g', axes (1:3), time, &
    ! 'graupel fall speed', 'm / s', missing_value = missing_value)
    ! id_vti = register_diag_field (mod_name, 'vt_i', axes (1:3), time, &
    ! 'ice fall speed', 'm / s', missing_value = missing_value)

    ! id_droplets = register_diag_field (mod_name, 'droplets', axes (1:3), time, &
    ! 'droplet number concentration', '# / m3', missing_value = missing_value)
    ! id_rh = register_diag_field (mod_name, 'rh_lin', axes (1:2), time, &
    ! 'relative humidity', 'n / a', missing_value = missing_value)

    ! id_rain = register_diag_field (mod_name, 'rain_lin', axes (1:2), time, &
    ! 'rain_lin', 'mm / day', missing_value = missing_value)
    ! id_snow = register_diag_field (mod_name, 'snow_lin', axes (1:2), time, &
    ! 'snow_lin', 'mm / day', missing_value = missing_value)
    ! id_graupel = register_diag_field (mod_name, 'graupel_lin', axes (1:2), time, &
    ! 'graupel_lin', 'mm / day', missing_value = missing_value)
    ! id_ice = register_diag_field (mod_name, 'ice_lin', axes (1:2), time, &
    ! 'ice_lin', 'mm / day', missing_value = missing_value)
    ! id_prec = register_diag_field (mod_name, 'prec_lin', axes (1:2), time, &
    ! 'prec_lin', 'mm / day', missing_value = missing_value)

    ! if (root_proc) write (*, *) 'prec_lin diagnostics initialized.', id_prec

    ! id_cond = register_diag_field (mod_name, 'cond_lin', axes (1:2), time, &
    ! 'total condensate', 'kg / m ** 2', missing_value = missing_value)
    ! id_var = register_diag_field (mod_name, 'var_lin', axes (1:2), time, &
    ! 'subgrid variance', 'n / a', missing_value = missing_value)

    ! call qsmith_init

    ! testing the water vapor tables

    ! if (mp_debug .and. root_proc) then
    ! write (*, *) 'testing water vapor tables in gfdl_cloud_microphys'
    ! tmp = tice - 90.
    ! do k = 1, 25
    ! q1 = wqsat_moist (tmp, 0., 1.e5)
    ! q2 = qs1d_m (tmp, 0., 1.e5)
    ! write (*, *) nint (tmp - tice), q1, q2, 'dq = ', q1 - q2
    ! tmp = tmp + 5.
    ! enddo
    ! endif

    ! if (root_proc) write (*, *) 'gfdl_cloud_micrphys diagnostics initialized.'

    ! gfdl_mp_clock = mpp_clock_id ('gfdl_cloud_microphys', grain = clock_routine)

    module_is_initialized = .true.

  end subroutine gfdl_cloud_microphys_init

  ! =======================================================================
  ! qsmith table initialization
  !>@brief The subroutine 'setup_con' sets up constants and calls 'qsmith_init'.
  ! =======================================================================

  subroutine setup_con

    implicit none

    ! root_proc = (mpp_pe () .eq.mpp_root_pe ())

    rgrav = 1. / grav

    if (.not. qsmith_tables_initialized) call qsmith_init

    qsmith_tables_initialized = .true.

  end subroutine setup_con

  ! =======================================================================
  !>@brief The function 'acr3d' is an accretion function (lin et al. 1983)
  ! =======================================================================

  real function acr3d (v1, v2, q1, q2, c, cac, rho)

    implicit none
    !$omp declare target

    real, intent (in) :: v1, v2, c, rho
    real, intent (in) :: q1, q2 ! mixing ratio!!!
    real, intent (in) :: cac (3)

    real :: t1, s1, s2

    ! integer :: k
    !
    ! real :: a
    !
    ! a = 0.0
    ! do k = 1, 3
    ! a = a + cac (k) * ((q1 * rho) ** ((7 - k) * 0.25) * (q2 * rho) ** (k * 0.25))
    ! enddo
    ! acr3d = c * abs (v1 - v2) * a / rho

    ! optimized

    t1 = sqrt (q1 * rho)
    s1 = sqrt (q2 * rho)
    s2 = sqrt (s1) ! s1 = s2 ** 2
    acr3d = c * abs (v1 - v2) * q1 * s2 * (cac (1) * t1 + cac (2) * sqrt (t1) * s2 + cac (3) * s1)

  end function acr3d

  ! =======================================================================
  !> melting of snow function (lin et al. 1983)
  !  note: psacw and psacr must be calc before smlt is called
  ! =======================================================================

  real function smlt (tc, dqs, qsrho, psacw, psacr, c, rho, rhofac)

    implicit none
    !$omp declare target

    real, intent (in) :: tc, dqs, qsrho, psacw, psacr, c (5), rho, rhofac

    smlt = (c (1) * tc / rho - c (2) * dqs) * (c (3) * sqrt (qsrho) + &
         c (4) * qsrho ** 0.65625 * sqrt (rhofac)) + c (5) * tc * (psacw + psacr)

  end function smlt

  ! =======================================================================
  !> melting of graupel function (lin et al. 1983)
  !  note: pgacw and pgacr must be calc before gmlt is called
  ! =======================================================================

  real function gmlt (tc, dqs, qgrho, pgacw, pgacr, c, rho)

    implicit none
    !$omp declare target

    real, intent (in) :: tc, dqs, qgrho, pgacw, pgacr, c (5), rho

    gmlt = (c (1) * tc / rho - c (2) * dqs) * (c (3) * sqrt (qgrho) + &
         c (4) * qgrho ** 0.6875 / rho ** 0.25) + c (5) * tc * (pgacw + pgacr)

  end function gmlt

  ! =======================================================================
  ! initialization
  ! prepare saturation water vapor pressure tables
  ! =======================================================================
  !>@brief The subroutine 'qsmith_init' initializes lookup tables for saturation
  !! water vapor pressure for the following utility routines that are designed
  !! to return qs consistent with the assumptions in FV3.
  !>@details The calculations are highly accurate values based on the Clausius-Clapeyron
  !! equation.
  ! =======================================================================
  subroutine qsmith_init

    implicit none

    integer, parameter :: length = TABLE_LENGTH

    integer :: i

    if (.not. tables_are_initialized) then

       ! root_proc = (mpp_pe () .eq. mpp_root_pe ())
       ! if (root_proc) print *, ' gfdl mp: initializing qs tables'

       ! debug code
       ! print *, mpp_pe (), allocated (table), allocated (table2), &
       ! allocated (table3), allocated (tablew), allocated (des), &
       ! allocated (des2), allocated (des3), allocated (desw)
       ! end debug code

       ! generate es table (dt = 0.1 deg. c)

       call qs_table (length)
       call qs_table2 (length)
       call qs_table3 (length)
       call qs_tablew (length)

       do i = 1, length - 1
          des (i) = max (0., table (i + 1) - table (i))
          des2 (i) = max (0., table2 (i + 1) - table2 (i))
          des3 (i) = max (0., table3 (i + 1) - table3 (i))
          desw (i) = max (0., tablew (i + 1) - tablew (i))
       enddo
       des (length) = des (length - 1)
       des2 (length) = des2 (length - 1)
       des3 (length) = des3 (length - 1)
       desw (length) = desw (length - 1)

       tables_are_initialized = .true.

    !$omp target update to(table2(:), des2(:), tablew(:), desw(:))

    endif

  end subroutine qsmith_init

  ! =======================================================================
  ! compute the saturated specific humidity for table ii
  !>@brief The function 'wqs1' returns the saturation vapor pressure over pure
  !! liquid water for a given temperature and air density.
  ! =======================================================================

  real function wqs1 (ta, den)

    implicit none
    !$omp declare target

    !> pure water phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min(2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs1 = es / (rvgas * ta * den)

  end function wqs1

  ! =======================================================================
  ! compute the gradient of saturated specific humidity for table ii
  !>@brief The function 'wqs2' returns the saturation vapor pressure over pure
  !! liquid water for a given temperature and air density, as well as the
  !! analytic dqs/dT: rate of change of saturation vapor pressure WRT temperature.
  ! =======================================================================

  real function wqs2 (ta, den, dqdt)

    implicit none
    !$omp declare target

    !> pure water phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real, intent (out) :: dqdt

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.

    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = 10. * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)

  end function wqs2

  ! =======================================================================
  ! compute wet buld temperature
  !>@brief The function 'wet_bulb' uses 'wqs2' to compute the wet-bulb temperature
  !! from the mixing ratio and the temperature.
  ! =======================================================================

  real function wet_bulb (q, t, den)

    implicit none

    real, intent (in) :: t, q, den

    real :: qs, tp, dqdt

    wet_bulb = t
    qs = wqs2 (wet_bulb, den, dqdt)
    tp = 0.5 * (qs - q) / (1. + lcp * dqdt) * lcp
    wet_bulb = wet_bulb - tp

    ! tp is negative if super - saturated
    if (tp > 0.01) then
       qs = wqs2 (wet_bulb, den, dqdt)
       tp = (qs - q) / (1. + lcp * dqdt) * lcp
       wet_bulb = wet_bulb - tp
    endif

  end function wet_bulb

  ! =======================================================================
  !>@brief The function 'iqs1' computes the saturated specific humidity
  !! for table iii
  ! =======================================================================

  real function iqs1 (ta, den)

    implicit none
    !$omp declare target

    !> water - ice phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs1 = es / (rvgas * ta * den)

  end function iqs1

  ! =======================================================================
  !>@brief The function 'iqs2' computes the gradient of saturated specific
  !! humidity for table iii
  ! =======================================================================

  real function iqs2 (ta, den, dqdt)

    implicit none
    !$omp declare target

    !> water - ice phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real, intent (out) :: dqdt

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min(2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    dqdt = 10. * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)

  end function iqs2

  ! =======================================================================
  !>@brief The function 'qs1d_moist' computes the gradient of saturated
  !! specific humidity for table iii.
  ! =======================================================================

  real function qs1d_moist (ta, qv, pa, dqdt)

    implicit none

    real, intent (in) :: ta, pa, qv

    real, intent (out) :: dqdt

    real :: es, ap1, tmin, eps10

    integer :: it

    tmin = table_ice - 160.
    eps10 = 10. * eps
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    qs1d_moist = eps * es * (1. + zvir * qv) / pa
    it = ap1 - 0.5
    dqdt = eps10 * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) * (1. + zvir * qv) / pa

  end function qs1d_moist

  ! =======================================================================
  ! compute the gradient of saturated specific humidity for table ii
  !>@brief The function 'wqsat2_moist' computes the saturated specific humidity
  !! for pure liquid water , as well as des/dT.
  ! =======================================================================

  real function wqsat2_moist (ta, qv, pa, dqdt)

    implicit none

    real, intent (in) :: ta, pa, qv

    real, intent (out) :: dqdt

    real :: es, ap1, tmin, eps10

    integer :: it

    tmin = table_ice - 160.
    eps10 = 10. * eps
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqsat2_moist = eps * es * (1. + zvir * qv) / pa
    it = ap1 - 0.5
    dqdt = eps10 * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) * (1. + zvir * qv) / pa

  end function wqsat2_moist

  ! =======================================================================
  ! compute the saturated specific humidity for table ii
  !>@brief The function 'wqsat_moist' computes the saturated specific humidity
  !! for pure liquid water.
  ! =======================================================================

  real function wqsat_moist (ta, qv, pa)

    implicit none

    real, intent (in) :: ta, pa, qv

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min(2621., ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqsat_moist = eps * es * (1. + zvir * qv) / pa

  end function wqsat_moist

  ! =======================================================================
  !>@brief The function 'qs1d_m' computes the saturated specific humidity
  !! for table iii
  ! =======================================================================

  real function qs1d_m (ta, qv, pa)

    implicit none

    real, intent (in) :: ta, pa, qv

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    qs1d_m = eps * es * (1. + zvir * qv) / pa

  end function qs1d_m

  ! =======================================================================
  !>@brief The function 'd_sat' computes the difference in saturation
  !! vapor * density * between water and ice
  ! =======================================================================

  real function d_sat (ta, den)

    implicit none

    real, intent (in) :: ta, den

    real :: es_w, es_i, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es_w = tablew (it) + (ap1 - it) * desw (it)
    es_i = table2 (it) + (ap1 - it) * des2 (it)
    d_sat = dim (es_w, es_i) / (rvgas * ta * den) ! take positive difference

  end function d_sat

  ! =======================================================================
  !>@brief The function 'esw_table' computes the saturated water vapor
  !! pressure for table ii
  ! =======================================================================

  real function esw_table (ta)

    implicit none

    real, intent (in) :: ta

    real :: ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    esw_table = tablew (it) + (ap1 - it) * desw (it)

  end function esw_table

  ! =======================================================================
  !>@brief The function 'es2_table' computes the saturated water
  !! vapor pressure for table iii
  ! =======================================================================

  real function es2_table (ta)

    implicit none

    real, intent (in) :: ta

    real :: ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (ta, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es2_table = table2 (it) + (ap1 - it) * des2 (it)

  end function es2_table

  ! =======================================================================
  !>@brief The subroutine 'esw_table1d' computes the saturated water vapor
  !! pressure for table ii.
  ! =======================================================================

  subroutine esw_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
       ap1 = 10. * dim (ta (i), tmin) + 1.
       ap1 = min (2621., ap1)
       it = ap1
       es (i) = tablew (it) + (ap1 - it) * desw (it)
    enddo

  end subroutine esw_table1d

  ! =======================================================================
  !>@brief The subroutine 'es3_table1d' computes the saturated water vapor
  !! pressure for table iii.
  ! =======================================================================

  subroutine es2_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
       ap1 = 10. * dim (ta (i), tmin) + 1.
       ap1 = min (2621., ap1)
       it = ap1
       es (i) = table2 (it) + (ap1 - it) * des2 (it)
    enddo

  end subroutine es2_table1d

  ! =======================================================================
  !>@brief The subroutine 'es3_table1d' computes the saturated water vapor
  !! pressure for table iv.
  ! =======================================================================

  subroutine es3_table1d (ta, es, n)

    implicit none

    integer, intent (in) :: n

    real, intent (in) :: ta (n)

    real, intent (out) :: es (n)

    real :: ap1, tmin

    integer :: i, it

    tmin = table_ice - 160.

    do i = 1, n
       ap1 = 10. * dim (ta (i), tmin) + 1.
       ap1 = min (2621., ap1)
       it = ap1
       es (i) = table3 (it) + (ap1 - it) * des3 (it)
    enddo

  end subroutine es3_table1d

  ! =======================================================================
  !>@brief saturation water vapor pressure table ii
  ! 1 - phase table
  ! =======================================================================

  subroutine qs_tablew (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: tmin, tem, fac0, fac1, fac2

    integer :: i

    tmin = table_ice - 160.

    ! -----------------------------------------------------------------------
    ! compute es over water
    ! -----------------------------------------------------------------------

    do i = 1, n
       tem = tmin + delt * real (i - 1)
       fac0 = (tem - t_ice) / (tem * t_ice)
       fac1 = fac0 * lv0
       fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
       tablew (i) = e00 * exp (fac2)
    enddo

  end subroutine qs_tablew

  ! =======================================================================
  !>@brief saturation water vapor pressure table iii
  ! 2 - phase table
  ! =======================================================================

  subroutine qs_table2 (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: tmin, tem0, tem1, fac0, fac1, fac2

    integer :: i, i0, i1

    tmin = table_ice - 160.

    do i = 1, n
       tem0 = tmin + delt * real (i - 1)
       fac0 = (tem0 - t_ice) / (tem0 * t_ice)
       if (i <= 1600) then
          ! -----------------------------------------------------------------------
          ! compute es over ice between - 160 deg c and 0 deg c.
          ! -----------------------------------------------------------------------
          fac1 = fac0 * li2
          fac2 = (d2ice * log (tem0 / t_ice) + fac1) / rvgas
       else
          ! -----------------------------------------------------------------------
          ! compute es over water between 0 deg c and 102 deg c.
          ! -----------------------------------------------------------------------
          fac1 = fac0 * lv0
          fac2 = (dc_vap * log (tem0 / t_ice) + fac1) / rvgas
       endif
       table2 (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! smoother around 0 deg c
    ! -----------------------------------------------------------------------

    i0 = 1600
    i1 = 1601
    tem0 = 0.25 * (table2 (i0 - 1) + 2. * table (i0) + table2 (i0 + 1))
    tem1 = 0.25 * (table2 (i1 - 1) + 2. * table (i1) + table2 (i1 + 1))
    table2 (i0) = tem0
    table2 (i1) = tem1

  end subroutine qs_table2

  ! =======================================================================
  !>@brief saturation water vapor pressure table iv
  ! 2 - phase table with " - 2 c" as the transition point
  ! =======================================================================

  subroutine qs_table3 (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: esbasw, tbasw, esbasi, tmin, tem, aa, b, c, d, e
    real :: tem0, tem1

    integer :: i, i0, i1

    esbasw = 1013246.0
    tbasw = table_ice + 100.
    esbasi = 6107.1
    tmin = table_ice - 160.

    do i = 1, n
       tem = tmin + delt * real (i - 1)
       ! if (i <= 1600) then
       if (i <= 1580) then ! change to - 2 c
          ! -----------------------------------------------------------------------
          ! compute es over ice between - 160 deg c and 0 deg c.
          ! see smithsonian meteorological tables page 350.
          ! -----------------------------------------------------------------------
          aa = - 9.09718 * (table_ice / tem - 1.)
          b = - 3.56654 * alog10 (table_ice / tem)
          c = 0.876793 * (1. - tem / table_ice)
          e = alog10 (esbasi)
          table3 (i) = 0.1 * 10 ** (aa + b + c + e)
       else
          ! -----------------------------------------------------------------------
          ! compute es over water between - 2 deg c and 102 deg c.
          ! see smithsonian meteorological tables page 350.
          ! -----------------------------------------------------------------------
          aa = - 7.90298 * (tbasw / tem - 1.)
          b = 5.02808 * alog10 (tbasw / tem)
          c = - 1.3816e-7 * (10 ** ((1. - tem / tbasw) * 11.344) - 1.)
          d = 8.1328e-3 * (10 ** ((tbasw / tem - 1.) * (- 3.49149)) - 1.)
          e = alog10 (esbasw)
          table3 (i) = 0.1 * 10 ** (aa + b + c + d + e)
       endif
    enddo

    ! -----------------------------------------------------------------------
    ! smoother around - 2 deg c
    ! -----------------------------------------------------------------------

    i0 = 1580
    i1 = 1581
    tem0 = 0.25 * (table3 (i0 - 1) + 2. * table (i0) + table3 (i0 + 1))
    tem1 = 0.25 * (table3 (i1 - 1) + 2. * table (i1) + table3 (i1 + 1))
    table3 (i0) = tem0
    table3 (i1) = tem1

  end subroutine qs_table3

  ! =======================================================================
  ! compute the saturated specific humidity for table
  ! note: this routine is based on "moist" mixing ratio
  !>@brief The function 'qs_blend' computes the saturated specific humidity
  !! with a blend of water and ice depending on the temperature.
  ! =======================================================================

  real function qs_blend (t, p, q)

    implicit none

    real, intent (in) :: t, p, q

    real :: es, ap1, tmin

    integer :: it

    tmin = table_ice - 160.
    ap1 = 10. * dim (t, tmin) + 1.
    ap1 = min (2621., ap1)
    it = ap1
    es = table (it) + (ap1 - it) * des (it)
    qs_blend = eps * es * (1. + zvir * q) / p

  end function qs_blend

  ! =======================================================================
  !>@brief saturation water vapor pressure table i
  ! 3 - phase table
  ! =======================================================================

  subroutine qs_table (n)

    implicit none

    integer, intent (in) :: n

    real :: delt = 0.1
    real :: tmin, tem, esh40
    real :: wice, wh2o, fac0, fac1, fac2
    real :: esupc (400)

    integer :: i
    real :: tc

    tmin = table_ice - 160.

    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg c and 0 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1600
       tem = tmin + delt * real (i - 1)
       fac0 = (tem - t_ice) / (tem * t_ice)
       fac1 = fac0 * li2
       fac2 = (d2ice * log (tem / t_ice) + fac1) / rvgas
       table (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! compute es over water between - 40 deg c and 102 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1421
       tem = 233.16 + delt * real (i - 1)
       fac0 = (tem - t_ice) / (tem * t_ice)
       fac1 = fac0 * lv0
       fac2 = (dc_vap * log (tem / t_ice) + fac1) / rvgas
       esh40 = e00 * exp (fac2)
       if (i <= 400) then
          esupc (i) = esh40
       else
          table (i + 1200) = esh40
       endif
    enddo

    ! -----------------------------------------------------------------------
    ! derive blended es over ice and supercooled water between - 40 deg c and 0 deg c
    ! -----------------------------------------------------------------------

    do i = 1, 400
       tem = 233.16 + delt * real (i - 1)
       ! wice = 0.05 * (table_ice - tem)
       ! wh2o = 0.05 * (tem - 253.16)
       ! GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C
       wice = ice_fraction(tem,0.0,0.0)
       wh2o = 1.0 - wice
       table (i + 1200) = wice * table (i + 1200) + wh2o * esupc (i)
    enddo

  end subroutine qs_table

  ! =======================================================================
  !>@brief The subroutine 'neg_adj' fixes negative water species.
  !>@details This is designed for 6-class micro-physics schemes.
  ! =======================================================================

  subroutine neg_adj (is, ie, js, je, ktop, kbot, pt, dp, qv, ql, qr, qi, qs, qg)

    implicit none

    integer, intent (in) :: is, ie, js, je, ktop, kbot

    real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: dp

    real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: pt, qv, ql, qr, qi, qs, qg

    real :: lcpk, icpk

    real :: dq, cvm

    integer :: i, j, k

    !$omp target data map(to: dp) map(tofrom: pt, qv, ql, qr, qi, qs, qg)

    ! -----------------------------------------------------------------------
    ! define heat capacity and latent heat coefficient
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(3)  private(cvm, lcpk, icpk)
    do k = ktop, kbot
       do i = is, ie
          do j = js, je

             cvm = c_air + &
                  qv (i, j, k) * c_vap + &
                  (qr (i, j, k) + &
                  ql (i, j, k)) * c_liq + &
                  (qi (i, j, k) + qs (i, j, k) + qg (i, j, k)) * c_ice
             lcpk = (lv00 + d0_vap * pt (i, j, k)) / cvm
             icpk = (li00 + dc_ice * pt (i, j, k)) / cvm

             ! -----------------------------------------------------------------------
             ! ice phase:
             ! -----------------------------------------------------------------------

             ! if cloud ice < 0, borrow from snow
             if (qi (i, j, k) < 0.) then
                qs (i, j, k) = qs (i, j, k) + qi (i, j, k)
                qi (i, j, k) = 0.
             endif

             ! if snow < 0, borrow from graupel
             if (qs (i, j, k) < 0.) then
                qg (i, j, k) = qg (i, j, k) + qs (i, j, k)
                qs (i, j, k) = 0.
             endif
             ! if graupel < 0, borrow from rain
             if (qg (i, j, k) < 0.) then
                qr (i, j, k) = qr (i, j, k) + qg (i, j, k)
                pt (i, j, k) = pt (i, j, k) - qg (i, j, k) * icpk ! heating
                qg (i, j, k) = 0.
             endif

             ! -----------------------------------------------------------------------
             ! liquid phase:
             ! -----------------------------------------------------------------------

             ! if rain < 0, borrow from cloud water
             if (qr (i, j, k) < 0.) then
                ql (i, j, k) = ql (i, j, k) + qr (i, j, k)
                qr (i, j, k) = 0.
             endif
             ! if cloud water < 0, borrow from water vapor
             if (ql (i, j, k) < 0.) then
                qv (i, j, k) = qv (i, j, k) + ql (i, j, k)
                pt (i, j, k) = pt (i, j, k) - ql (i, j, k) * lcpk ! heating
                ql (i, j, k) = 0.
             endif

          enddo
       enddo
    enddo
    !$omp end target teams distribute parallel do simd

    ! -----------------------------------------------------------------------
    ! fix water vapor; borrow from below
    ! -----------------------------------------------------------------------

    do k = ktop, kbot - 1
       !$omp target teams distribute parallel do simd collapse(2)
       do i = is, ie
          do j = js, je
             if (qv (i, j, k) < 0.) then
                qv (i, j, k + 1) = qv (i, j, k + 1) + qv (i, j, k) * dp (i, j, k) / dp (i, j, k + 1)
                qv (i, j, k) = 0.
             endif
          enddo
       enddo
       !$omp end target teams distribute parallel do simd
    enddo

    ! -----------------------------------------------------------------------
    ! bottom layer; borrow from above
    ! -----------------------------------------------------------------------

    !$omp target teams distribute parallel do simd collapse(2) private(dq)
    do i = is, ie
       do j = js, je
          if (qv (i, j, kbot) < 0. .and. qv (i, j, kbot - 1) > 0.) then
             dq = min (- qv (i, j, kbot) * dp (i, j, kbot), qv (i, j, kbot - 1) * dp (i, j, kbot - 1))
             qv (i, j, kbot - 1) = qv (i, j, kbot - 1) - dq / dp (i, j, kbot - 1)
             qv (i, j, kbot) = qv (i, j, kbot) + dq / dp (i, j, kbot)
          endif
       enddo
    enddo
    !$omp end target teams distribute parallel do simd

    !$omp end target data

  end subroutine neg_adj

  ! =======================================================================
  ! compute global sum
  !>@brief quick local sum algorithm
  ! =======================================================================

  !real function g_sum (p, ifirst, ilast, jfirst, jlast, area, mode)
  !
  ! use mpp_mod, only: mpp_sum
  !
  ! implicit none
  !
  ! integer, intent (in) :: ifirst, ilast, jfirst, jlast
  ! integer, intent (in) :: mode ! if == 1 divided by area
  !
  ! real, intent (in), dimension (ifirst:ilast, jfirst:jlast) :: p, area
  !
  ! integer :: i, j
  !
  ! real :: gsum
  !
  ! if (global_area < 0.) then
  ! global_area = 0.
  ! do j = jfirst, jlast
  ! do i = ifirst, ilast
  ! global_area = global_area + area (i, j)
  ! enddo
  ! enddo
  ! call mpp_sum (global_area)
  ! endif
  !
  ! gsum = 0.
  ! do j = jfirst, jlast
  ! do i = ifirst, ilast
  ! gsum = gsum + p (i, j) * area (i, j)
  ! enddo
  ! enddo
  ! call mpp_sum (gsum)
  !
  ! if (mode == 1) then
  ! g_sum = gsum / global_area
  ! else
  ! g_sum = gsum
  ! endif
  !
  !end function g_sum

  ! ==========================================================================
  !>@brief The subroutine 'interpolate_z' interpolates to a prescribed height.
  ! ==========================================================================

  subroutine interpolate_z (is, ie, js, je, km, zl, hgt, a3, a2)

    implicit none

    integer, intent (in) :: is, ie, js, je, km

    real, intent (in), dimension (is:ie, js:je, km) :: a3

    real, intent (in), dimension (is:ie, js:je, km + 1) :: hgt !< hgt (k) > hgt (k + 1)

    real, intent (in) :: zl

    real, intent (out), dimension (is:ie, js:je) :: a2

    real, dimension (km) :: zm !< middle layer height

    integer :: i, j, k

    ! !$omp parallel do default (none) shared (is, ie, js, je, km, hgt, zl, a2, a3) private (zm)

    do j = js, je
       do i = is, ie
          do k = 1, km
             zm (k) = 0.5 * (hgt (i, j, k) + hgt (i, j, k + 1))
          enddo
          if (zl >= zm (1)) then
             a2 (i, j) = a3 (i, j, 1)
          elseif (zl <= zm (km)) then
             a2 (i, j) = a3 (i, j, km)
          else
             do k = 1, km - 1
                if (zl <= zm (k) .and. zl >= zm (k + 1)) then
                   a2 (i, j) = a3 (i, j, k) + (a3 (i, j, k + 1) - a3 (i, j, k)) * (zm (k) - zl) / (zm (k) - zm (k + 1))
                   exit
                endif
             enddo
          endif
       enddo
    enddo

  end subroutine interpolate_z

  ! =======================================================================
  !>@brief The subroutine 'cloud_diagnosis' diagnoses the radius of cloud
  !! species.
  ! =======================================================================

  subroutine cloud_diagnosis (is, ie, js, je, den, qw, qi, qr, qs, qg, t, &
       qcw, qci, qcr, qcs, qcg, rew, rei, rer, res, reg)

    implicit none

    integer, intent (in) :: is, ie, js, je

    real, intent (in), dimension (is:ie, js:je) :: den, t
    real, intent (in), dimension (is:ie, js:je) :: qw, qi, qr, qs, qg !< units: kg / kg

    real, intent (out), dimension (is:ie, js:je) :: qcw, qci, qcr, qcs, qcg !< units: kg / m^3
    real, intent (out), dimension (is:ie, js:je) :: rew, rei, rer, res, reg !< units: micron

    integer :: i, j

    real :: lambdar, lambdas, lambdag

    real :: rhow = 1.0e3, rhor = 1.0e3, rhos = 1.0e2, rhog = 4.0e2
    real :: n0r = 8.0e6, n0s = 3.0e6, n0g = 4.0e6
    real :: alphar = 0.8, alphas = 0.25, alphag = 0.5
    real :: gammar = 17.837789, gammas = 8.2850630, gammag = 11.631769
    real :: qmin = 1.0e-5, ccn = 1.0e8, beta = 1.22

    ! real :: rewmin = 1.0, rewmax = 25.0
    ! real :: reimin = 10.0, reimax = 300.0
    ! real :: rermin = 25.0, rermax = 225.0
    ! real :: resmin = 300, resmax = 1000.0
    ! real :: regmin = 1000.0, regmax = 1.0e5
    real :: rewmin = 5.0, rewmax = 10.0
    real :: reimin = 10.0, reimax = 150.0
    real :: rermin = 0.0, rermax = 10000.0
    real :: resmin = 0.0, resmax = 10000.0
    real :: regmin = 0.0, regmax = 10000.0

    do j = js, je
       do i = is, ie

          ! -----------------------------------------------------------------------
          ! cloud water (martin et al., 1994)
          ! -----------------------------------------------------------------------

          if (qw (i, j) .gt. qmin) then
             qcw (i, j) = den (i, j) * qw (i, j)
             rew (i, j) = exp (1.0 / 3.0 * log ((3 * qcw (i, j)) / (4 * pi * rhow * ccn))) * 1.0e6
             rew (i, j) = max (rewmin, min (rewmax, rew (i, j)))
          else
             qcw (i, j) = 0.0
             rew (i, j) = rewmin
          endif

          ! -----------------------------------------------------------------------
          ! cloud ice (heymsfield and mcfarquhar, 1996)
          ! -----------------------------------------------------------------------

          if (qi (i, j) .gt. qmin) then
             qci (i, j) = den (i, j) * qi (i, j)
             if (t (i, j) - tice .lt. - 50) then
                rei (i, j) = beta / 9.917 * exp ((1 - 0.891) * log (1.0e3 * qci (i, j))) * 1.0e3
             elseif (t (i, j) - tice .lt. - 40) then
                rei (i, j) = beta / 9.337 * exp ((1 - 0.920) * log (1.0e3 * qci (i, j))) * 1.0e3
             elseif (t (i, j) - tice .lt. - 30) then
                rei (i, j) = beta / 9.208 * exp ((1 - 0.945) * log (1.0e3 * qci (i, j))) * 1.0e3
             else
                rei (i, j) = beta / 9.387 * exp ((1 - 0.969) * log (1.0e3 * qci (i, j))) * 1.0e3
             endif
             rei (i, j) = max (reimin, min (reimax, rei (i, j)))
          else
             qci (i, j) = 0.0
             rei (i, j) = reimin
          endif

          ! -----------------------------------------------------------------------
          ! rain (lin et al., 1983)
          ! -----------------------------------------------------------------------

          if (qr (i, j) .gt. qmin) then
             qcr (i, j) = den (i, j) * qr (i, j)
             lambdar = exp (0.25 * log (pi * rhor * n0r / qcr (i, j)))
             rer (i, j) = 0.5 * exp (log (gammar / 6) / alphar) / lambdar * 1.0e6
             rer (i, j) = max (rermin, min (rermax, rer (i, j)))
          else
             qcr (i, j) = 0.0
             rer (i, j) = rermin
          endif

          ! -----------------------------------------------------------------------
          ! snow (lin et al., 1983)
          ! -----------------------------------------------------------------------

          if (qs (i, j) .gt. qmin) then
             qcs (i, j) = den (i, j) * qs (i, j)
             lambdas = exp (0.25 * log (pi * rhos * n0s / qcs (i, j)))
             res (i, j) = 0.5 * exp (log (gammas / 6) / alphas) / lambdas * 1.0e6
             res (i, j) = max (resmin, min (resmax, res (i, j)))
          else
             qcs (i, j) = 0.0
             res (i, j) = resmin
          endif

          ! -----------------------------------------------------------------------
          ! graupel (lin et al., 1983)
          ! -----------------------------------------------------------------------

          if (qg (i, j) .gt. qmin) then
             qcg (i, j) = den (i, j) * qg (i, j)
             lambdag = exp (0.25 * log (pi * rhog * n0g / qcg (i, j)))
             reg (i, j) = 0.5 * exp (log (gammag / 6) / alphag) / lambdag * 1.0e6
             reg (i, j) = max (regmin, min (regmax, reg (i, j)))
          else
             qcg (i, j) = 0.0
             reg (i, j) = regmin
          endif

       enddo
    enddo

  end subroutine cloud_diagnosis

  real function new_ice_condensate(tk, qlk, qik, cnv_fraction, srf_type)
    !$omp declare target

    real, intent(in) :: tk, qlk, qik, cnv_fraction, srf_type
    real :: ptc, ifrac

    ifrac = ice_fraction(tk,cnv_fraction, srf_type)
    if (qlk > qcmin) then
       new_ice_condensate = max(0.0,min(ifrac*(qlk+qik) - qik, qlk))
    else
       new_ice_condensate = 0.0
    endif

  end function new_ice_condensate

  real function new_liq_condensate(tk, qlk, qik, cnv_fraction, srf_type)

    real, intent(in) :: tk, qlk, qik, cnv_fraction, srf_type
    real :: lfrac

    lfrac = 1.0 - ice_fraction(tk,cnv_fraction, srf_type)
    if (qik > qcmin) then
       new_liq_condensate = max(0.0,min(lfrac*(qlk+qik) - qlk, qik))
    else
       new_liq_condensate = 0.0
    endif

  end function new_liq_condensate

  function ICE_FRACTION (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
    !$acc routine seq
    !$omp declare target
    real, intent(in) :: TEMP,CNV_FRACTION,SRF_TYPE
    real             :: ICEFRCT
    real             :: tc, ptc
    real             :: ICEFRCT_C, ICEFRCT_M

    ! In anvil/convective clouds
    real, parameter :: aT_ICE_ALL = 252.16
    real, parameter :: aT_ICE_MAX = 268.16
    real, parameter :: aICEFRPWR  = 2.0
    ! Over snow/ice SRF_TYPE = 2
    real, parameter :: iT_ICE_ALL = 236.16
    real, parameter :: iT_ICE_MAX = 261.16
    real, parameter :: iICEFRPWR  = 6.0
    ! Over Land     SRF_TYPE = 1
    real, parameter :: lT_ICE_ALL = 239.16
    real, parameter :: lT_ICE_MAX = 261.16
    real, parameter :: lICEFRPWR  = 2.0
    ! Over Oceans   SRF_TYPE = 0
    real, parameter :: oT_ICE_ALL = 238.16
    real, parameter :: oT_ICE_MAX = 263.16
    real, parameter :: oICEFRPWR  = 4.0

    ! Anvil clouds
    ! Anvil-Convective sigmoidal function like figure 6(right)
    ! Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
    ICEFRCT_C  = 0.00
    if ( TEMP <= aT_ICE_ALL ) then
       ICEFRCT_C = 1.000
    else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
       ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
    end if
    ICEFRCT_C = MIN(ICEFRCT_C,1.00)
    ICEFRCT_C = MAX(ICEFRCT_C,0.00)
    ICEFRCT_C = ICEFRCT_C**aICEFRPWR
#ifdef MODIS_ICE_POLY
    ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384)
    tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
    ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
    ICEFRCT_M = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
#else
    ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
    if (SRF_TYPE == 2.0) then
       ! Over snow/ice
       ICEFRCT_M  = 0.00
       if ( TEMP <= iT_ICE_ALL ) then
          ICEFRCT_M = 1.000
       else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
          ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL ) ) )
       end if
       ICEFRCT_M = MIN(ICEFRCT_M,1.00)
       ICEFRCT_M = MAX(ICEFRCT_M,0.00)
       ICEFRCT_M = ICEFRCT_M**iICEFRPWR
    else if (SRF_TYPE > 1.0) then
       ! Over Land
       ICEFRCT_M  = 0.00
       if ( TEMP <= lT_ICE_ALL ) then
          ICEFRCT_M = 1.000
       else if ( (TEMP > lT_ICE_ALL) .AND. (TEMP <= lT_ICE_MAX) ) then
          ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - lT_ICE_ALL ) / ( lT_ICE_MAX - lT_ICE_ALL ) ) )
       end if
       ICEFRCT_M = MIN(ICEFRCT_M,1.00)
       ICEFRCT_M = MAX(ICEFRCT_M,0.00)
       ICEFRCT_M = ICEFRCT_M**lICEFRPWR
    else
       ! Over Oceans
       ICEFRCT_M  = 0.00
       if ( TEMP <= oT_ICE_ALL ) then
          ICEFRCT_M = 1.000
       else if ( (TEMP > oT_ICE_ALL) .AND. (TEMP <= oT_ICE_MAX) ) then
          ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - oT_ICE_ALL ) / ( oT_ICE_MAX - oT_ICE_ALL ) ) )
       end if
       ICEFRCT_M = MIN(ICEFRCT_M,1.00)
       ICEFRCT_M = MAX(ICEFRCT_M,0.00)
       ICEFRCT_M = ICEFRCT_M**oICEFRPWR
    endif
#endif
    ! Combine the Convective and MODIS functions
    ICEFRCT  = ICEFRCT_M*(1.0-CNV_FRACTION) + ICEFRCT_C*(CNV_FRACTION)

  end function ICE_FRACTION

end module gfdl2_cloud_microphys_mod
