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

  ! use mpp_mod, only: mpp_pe, mpp_root_pe

  ! use fms_mod, only: write_version_number, open_namelist_file, &
  !                    check_nml_error, close_file, file_exist,  &
  !                    fms_init
    use GEOSmoist_Process_Library, only: sigma, ice_fraction, LDRADIUS4
    use MAPL, only: MAPL_AM_I_ROOT

    implicit none

    private

    public gfdl_cloud_microphys_driver, gfdl_cloud_microphys_init, gfdl_cloud_microphys_end
    public cloud_diagnosis
    public ICE_LSC_VFALL_PARAM, ICE_CNV_VFALL_PARAM

    integer :: ICE_LSC_VFALL_PARAM = 1
    integer :: ICE_CNV_VFALL_PARAM = 1

    real :: missing_value = - 1.e10

    logical :: module_is_initialized = .false.
    logical :: qsmith_tables_initialized = .false.

    character (len = 17) :: mod_name = 'gfdl_cloud_microphys'

    real, parameter :: grav = 9.80665 !< gfs: acceleration due to gravity
    real, parameter :: rdgas = 287.05 !< gfs: gas constant for dry air
    real, parameter :: rvgas = 461.50 !< gfs: gas constant for water vapor
    real, parameter :: cp_air = 1004.6 !< gfs: heat capacity of dry air at constant pressure
    real, parameter :: hlv = 2.5e6 !< gfs: latent heat of evaporation
    real, parameter :: hlf = 3.3358e5 !< gfs: latent heat of fusion
    real, parameter :: pi = 3.1415926535897931 !< gfs: ratio of circle circumference to diameter
    real, parameter :: cp_vap = 4.0 * rvgas !< 1846.0, heat capacity of water vapore at constnat pressure
    real, parameter :: cv_air = cp_air - rdgas !< 717.55, heat capacity of dry air at constant volume
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

    integer, parameter :: es_table_length = 2621
    real   , parameter :: es_table_tmin = table_ice - 160.
    real   , parameter :: delt = 0.1
    real   , parameter :: rdelt = 1.0/delt

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

    real, parameter :: dz_min = 1.e-2 ! use for correcting flipped height

    real, parameter :: sfcrho = 1.2 !< surface air density
    real, parameter :: rhor = 1.e3 !< density of rain water, lin83
    
    real, parameter :: rc = (4. / 3.) * pi * rhor

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

    integer :: icloud_f = 3 !< cloud scheme
    integer :: irain_f = 0 !< cloud water to rain auto conversion scheme

    logical :: de_ice = .false. !< to prevent excessive build - up of cloud ice from external sources
    logical :: sedi_transport = .true. !< transport of momentum in sedimentation
    logical :: do_sedi_w = .false. !< transport of vertical motion in sedimentation
    logical :: do_sedi_heat = .false. !< transport of heat in sedimentation
    logical :: do_bigg = .false. !< do bigg mechanism freezing of supercooled liquid on aerosol nuclei
    logical :: do_evap = .true. !< do evaporation
    logical :: do_subl = .true. !< do sublimation
    logical :: in_cloud = .true. !< use in-cloud autoconversion
    logical :: do_qa = .false. !< do inline cloud fraction (WMP: in FV3 dynamics)
    logical :: do_icepsettle = .true. ! include ice pressure settling function
    logical :: preciprad = .true. !< consider precipitates in cloud fraciton calculation
    logical :: fix_negative = .true. !< fix negative water species
    logical :: do_setup = .true. !< setup constants and parameters
    logical :: p_nonhydro = .false. !< perform hydrosatic adjustment on air density

    real, allocatable :: table (:), table2 (:), table3 (:), tablew (:)
    real, allocatable :: des (:), des2 (:), des3 (:), desw (:)

    logical :: tables_are_initialized = .false.

    real, parameter :: dt_fr = 8. !< epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
    ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
    ! dt_fr can be considered as the error bar

    real :: p_min = 100. !< minimum pressure (pascal) for mp to operate

    ! -----------------------------------------------------------------------
    !> namelist parameters
    ! -----------------------------------------------------------------------

    real :: cld_min = 0.05 !< minimum cloud fraction
    real :: tice = 273.16 !< set tice = 165. to trun off ice - phase phys (kessler emulator)

    real :: log_10 = log (10.)
    real :: tice0 = 273.16 - 0.01
    real :: t_wfr = 273.16 - 40.0 ! supercooled water can exist down to - 40 c, which is the "absolute"

    real :: t_min = 273.16 - 95.16 !< min temp to freeze all water vapor
    real :: t_sub = 273.16 - 89.16 !< min temp for sublimation of cloud ice
    real :: mp_time = 150. !< maximum micro - physics time step (sec)

    ! relative humidity increment

    real :: rh_inc = 0.30 !< rh increment for complete evaporation of cloud water and cloud ice
    real :: rh_inr = 0.30 !< rh increment for minimum evaporation of rain
    real :: rh_ins = -99. !< rh increment for sublimation of snow [wmp: not used]

    ! conversion time scale

    real :: tau_r2g  = -9999. !< rain freezing during fast_sat [wmp: not used]
    real :: tau_l2r  = -9999. !< cloud water to rain auto - conversion [wmp: not used]
    real :: tau_v2l  = -9999. !< water vapor to cloud water (condensation) [wmp: not used]
    real :: tau_l2v  =   300. !< cloud water to water vapor (evaporation)
    real :: tau_i2v  =   300. !< cloud ice to water vapor (sublimation)
    real :: tau_s2v  =   600. !< snow sublimation
    real :: tau_g2v  =   900. !< graupel sublimation
    real :: tau_v2s  = 21600. !< snow deposition -- make it a slow process
    real :: tau_v2g  = 21600. !< graupel deposition -- make it a slow process
    real :: tau_i2s  =  1000. !< cloud ice to snow auto - conversion
    real :: tau_revp =   600. !< rain re-evaporation
    real :: tau_frz  =   600. !< timescale for liquid-ice freezing
    real :: tau_imlt =   600. !< cloud ice melting
    real :: tau_smlt =   900. !< snow melting
    real :: tau_gmlt =  1200. !< graupel melting to rain

    real :: rthreshu =  7.0e-6 !< unstable critical cloud drop radius (micro m)
    real :: rthreshs = 10.0e-6 !<   stable critical cloud drop radius (micro m)
    
    real :: sat_adj0 = 0.90 !< adjustment factor (0: no, 1: full) during fast_sat_adj

    real :: qc_crt = 5.0e-8 !< mini condensate mixing ratio to allow partial cloudiness

    real :: qi_lim = 1. !< cloud ice limiter to prevent large ice build up by deposition

    real :: ql_mlt = 2.0e-3 !< max value of cloud water allowed from melted cloud ice
    real :: qs_mlt = 1.0e-6 !< max cloud water due to snow melt
    
    real :: qi_gen = 9.82679e-5 !< max cloud ice generation at -40 C

    ! cloud condensate upper bounds: "safety valves" for ql & qi
    real :: ql0_max = 2.0e-3 !< max cloud water value (auto converted to rain)

    ! critical autoconverion parameters
    real :: qi0_crt = 1.0e-4 !< cloud ice to snow autoconversion threshold
                             !! qi0_crt can be dependent on horizontal resolution
                             !! this sensitivity could be handled with onemsig later in the code
    real :: qs0_crt = 0.6e-3 !< snow to graupel density threshold (0.6e-3 in purdue lin scheme)

    real :: c_paut  = 1.00  !< autoconversion cloud water to rain (use 0.5 to reduce autoconversion)

    ! collection efficiencies for accretion
    !   Dry processes (frozen to/from frozen)
    real :: c_psaci = 0.05  !< accretion: cloud ice to snow
    real :: c_pgacs = 0.01  !< accretion: snow to graupel
    real :: c_pgaci = 0.01  !< accretion: cloud ice to graupel
    !   Wet processes (liquid to/from frozen)
    real :: c_cracw = 1.00  !< accretion: cloud water to rain

    ! accretion efficiencies
    real :: alin = 842.0 !< "a" in lin 1983, [Rain] (increase to ehance ql/qi -- > qr)
    real :: clin = 4.8   !< "c" in lin 1983, [Snow] (increase to ehance ql/qi -- > qs)
    real :: gcon = 40.74 * sqrt (sfcrho) ! [Graupel] (increase to ehance ql/qi -- > qg)

    ! fall velocity tuning constants:

    logical :: const_vi = .false. !< if .t. the constants are specified by v * _fac
    logical :: const_vs = .false. !< if .t. the constants are specified by v * _fac
    logical :: const_vg = .false. !< if .t. the constants are specified by v * _fac
    logical :: const_vr = .false. !< if .t. the constants are specified by v * _fac

    ! ice fall speed ranges based on https://doi.org/10.1002/2013JD020602 fig. 9
    ! bounds of fall speed (with variable speed option) for precip base on
    ! https://www.atmos.albany.edu/facstaff/rfovell/ATM562/lin-etal-1983.pdf fig. 2

    real :: vi_fac = 1.
    real :: vs_fac = 1.
    real :: vg_fac = 1.
    real :: vr_fac = 1.
    real :: vh_fac = 1.

    real :: vi_min = 0.01 !< minimum fall speed or constant fall speed
    real :: vs_min = 1.   !< minimum fall speed or constant fall speed
    real :: vg_min = 3.   !< minimum fall speed or constant fall speed
    real :: vr_min = 4.   !< minimum fall speed or constant fall speed
    real :: vh_min = 9.   !< minimum fall speed or constant fall speed
    
    real :: vi_max = 1.0 !< max fall speed for ice
    real :: vs_max = 3.0 !< max fall speed for snow
    real :: vg_max = 12.0 !< max fall speed for graupel
    real :: vr_max = 12.0 !< max fall speed for rain
    real :: vh_max = 19.0 !< max fall speed for hail

    ! cloud microphysics switchers

    logical :: fast_sat_adj = .false. !< has fast saturation adjustments
    logical :: z_slope_liq  = .true.  !< use linear mono slope for autocconversions
    logical :: z_slope_ice  = .true.  !< use linear mono slope for autocconversions
    logical :: use_ppm      = .false. !< use ppm fall scheme
    logical :: mono_prof    = .true.  !< perform terminal fall with mono ppm scheme
    logical :: mp_print     = .false. !< cloud microphysics debugging printout

    ! real :: global_area = - 1.

    ! -----------------------------------------------------------------------
    ! namelist
    ! -----------------------------------------------------------------------

    namelist / gfdl_cloud_microphysics_nml /                                  &
        mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt,                   &
        vi_fac, vr_fac, vs_fac, vg_fac,                                       &
        vi_min, vr_min, vs_min, vg_min, ql_mlt, do_qa, fix_negative, vi_max,  &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, &
        qi0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
        const_vs, const_vg, const_vr, rthreshu, rthreshs, qc_crt, &
        tau_g2v, tau_v2g, tau_s2v, tau_v2s, &
        tau_revp, tau_frz, do_bigg, do_evap, do_subl, &
        sat_adj0, tau_imlt, tau_v2l, tau_l2v, tau_i2v, &
        tau_i2s, tau_l2r, qi_lim, c_paut, c_psaci, c_pgacs, c_pgaci,  &
        z_slope_liq, z_slope_ice, c_cracw, alin, clin,              &
        preciprad, cld_min, use_ppm, mono_prof, in_cloud,         &
        do_icepsettle, &
        do_sedi_heat, sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f, mp_print

    public                                                                    &
        mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_gmlt,                   &
        vi_min, vr_min, vs_min, vg_min, ql_mlt, do_qa, fix_negative, vi_max,  &
        vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, &
        qi0_crt, fast_sat_adj, rh_inc, rh_ins, rh_inr, const_vi,     &
        const_vs, const_vg, const_vr, rthreshu, rthreshs, qc_crt, &
        tau_g2v, tau_v2g, tau_s2v, tau_v2s, &
        tau_revp, tau_frz, do_bigg, do_evap, do_subl, &
        sat_adj0, tau_imlt, tau_v2l, tau_l2v, tau_i2v, &
        tau_i2s, tau_l2r, qi_lim, c_paut, c_psaci, c_pgacs, c_pgaci,  &
        z_slope_liq, z_slope_ice, c_cracw, alin, clin,              &
        preciprad, cld_min, use_ppm, mono_prof, in_cloud,         &
        do_icepsettle, &
        do_sedi_heat, sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f, mp_print

contains

! -----------------------------------------------------------------------
! the driver of the gfdl cloud microphysics
! -----------------------------------------------------------------------

!>@brief The subroutine 'gfdl_cloud_microphys_driver' executes the full GFDL
!! cloud microphysics.
subroutine gfdl_cloud_microphys_driver (qv, ql, qr, qi, qs, qg, qa, qnl, qni,   &
        qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, pt_dt, pt, w,    &
        uin, vin, udt, vdt, dz, delp, area, dt_in,                        &
        land, cnv_fraction, srf_type, eis,                                &
        rhcrit, anv_icefall, lsc_icefall,                                 &
        revap, isubl, vti, vts, vtg, vtr,                                 &
        rain, snow, ice,                                                  &
        graupel, m2_rain, m2_sol, hydrostatic, phys_hydrostatic,          &
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
    real, intent (in), dimension (:, :, :) :: pt, qv, ql, qr, qg, qa, qnl, qni

    real, intent (inout), dimension (:, :, :) :: qi, qs
    real, intent (inout), dimension (:, :, :) :: pt_dt, qa_dt, udt, vdt, w
    real, intent (inout), dimension (:, :, :) :: qv_dt, ql_dt, qr_dt
    real, intent (inout), dimension (:, :, :) :: qi_dt, qs_dt, qg_dt

    real, intent (out), dimension (:, :) :: rain, snow, ice, graupel
    real, intent (out), dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
    real, intent (out), dimension (:, :, :) :: revap ! Rain evaporation
    real, intent (out), dimension (:, :, :) :: isubl ! Ice sublimation
    real, intent (out), dimension (:, :, :) :: vti, vts, vtg, vtr ! Fall speed exports

    ! logical :: used

    real :: mpdt, rdt, dts, convt, tot_prec

    integer :: i, j, k
    integer :: is, ie, js, je !< physics window
    integer :: ks, ke !< vertical dimension
    integer :: days, ntimes

    real, dimension (iie - iis + 1, jje - jjs + 1) :: prec_mp, prec1, cond, w_var, rh0

    real, dimension (iie - iis + 1, jje - jjs + 1, kke - kks + 1) :: qn2

    real :: allmax

    is = 1
    js = 1
    ks = 1
    ie = iie - iis + 1
    je = jje - jjs + 1
    ke = kke - kks + 1

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

    do j = js, je
        call mpdrv (hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, qg,&
            qa, qnl, qni, dz, is, ie, js, je, ks, ke, ktop, kbot, j, dt_in, ntimes,  &
            rain (:, j), snow (:, j), graupel (:, j), ice (:, j), m2_rain,     &
            m2_sol, cond (:, j), area (:, j),                                  &
            land (:, j), cnv_fraction(:, j), srf_type(:, j), eis(:,j),         &
            rhcrit, anv_icefall, lsc_icefall,                                  &
            revap, isubl,                                                      &
            udt, vdt, pt_dt,                                                   &
            qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, w_var,            &
            vtr, vts, vtg, vti, qn2)
    enddo

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
subroutine mpdrv (hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs,     &
        qg, qa, qnl, qni, dz, is, ie, js, je, ks, ke, ktop, kbot, j, dt_in, ntimes, &
        rain, snow, graupel, ice, m2_rain, m2_sol, cond, area1, land, &
        cnv_fraction, srf_type, eis, rhcrit, anv_icefall, lsc_icefall, revap, isubl, &
        u_dt, v_dt, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,   &
        w_var, vt_r, vt_s, vt_g, vt_i, qn2)

    implicit none

    logical, intent (in) :: hydrostatic

    integer, intent (in) :: j, is, ie, js, je, ks, ke
    integer, intent (in) :: ntimes, ktop, kbot

    real, intent (in) :: dt_in

    real, intent (in), dimension (is:) :: area1, land
    real, intent (in), dimension (is:) :: cnv_fraction
    real, intent (in), dimension (is:) :: srf_type
    real, intent (in), dimension (is:) :: eis

    real, intent (in), dimension (is:, js:, ks:) :: rhcrit
   
    real, intent (in) :: anv_icefall, lsc_icefall

    real, intent (in), dimension (is:, js:, ks:) :: uin, vin, delp, pt, dz
    real, intent (in), dimension (is:, js:, ks:) :: qv, qi, ql, qr, qs, qg, qa, qnl, qni
    
    real, intent (inout), dimension (is:, js:, ks:) :: u_dt, v_dt, w, pt_dt, qa_dt
    real, intent (inout), dimension (is:, js:, ks:) :: qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
    real, intent (  out), dimension (is:, js:, ks:) :: revap, isubl

    real, intent (inout), dimension (is:) :: rain, snow, ice, graupel, cond

    real, intent (out), dimension (is:, js:) :: w_var

    real, intent (out), dimension (is:, js:, ks:) :: vt_r, vt_s, vt_g, vt_i, qn2

    real, intent (out), dimension (is:, js:, ks:) :: m2_rain, m2_sol
    
    real, dimension (ktop:kbot) :: h_var1d 
    real, dimension (ktop:kbot) :: qvz, qlz, qrz, qiz, qsz, qgz, qaz
    real, dimension (ktop:kbot) :: vtiz, vtsz, vtgz, vtrz
    real, dimension (ktop:kbot) :: dp0, dp1, dz0, dz1
    real, dimension (ktop:kbot) :: qv0, ql0, qr0, qi0, qs0, qg0, qa0
    real, dimension (ktop:kbot) :: t0, den, den0, tz, p1, denfac
    real, dimension (ktop:kbot) :: ccn_l, ccn_i, c_praut, m1_rain, m1_sol, m1, evap1, subl1
    real, dimension (ktop:kbot) :: u0, v0, u1, v1, w1

    real :: onemsig, fac_eis
    real :: cpaut, rh_adj, rh_rain
    real :: r1, s1, i1, g1, rdt
    real :: dts
    real :: s_leng, t_land, t_ocean, h_var
    real :: cvm, tmp, omq
    real :: dqi, qio, qin

    integer :: i, k, n

    dts = dt_in / real (ntimes)
    rdt = 1. / dt_in

    ! -----------------------------------------------------------------------
    ! use local variables
    ! -----------------------------------------------------------------------

    do i = is, ie

        do k = ktop, kbot

            t0 (k) = pt (i, j, k)
            tz (k) = t0 (k)
            dp1 (k) = delp (i, j, k)
            dp0 (k) = dp1 (k) ! moist air mass * grav
            
            ! -----------------------------------------------------------------------
           ! import horizontal subgrid variability with pressure dependence
           ! total water subgrid deviation in horizontal direction
           ! default area dependent form: use dx ~ 100 km as the base
           ! -----------------------------------------------------------------------
            h_var1d(k) = min(0.30,1.0 - rhcrit(i,j,k)) ! restricted to 70%

            ! -----------------------------------------------------------------------
            ! convert moist mixing ratios to dry mixing ratios
            ! -----------------------------------------------------------------------
            
            qvz (k) = qv (i, j, k)
            qlz (k) = ql (i, j, k)
            qiz (k) = qi (i, j, k)
            qrz (k) = qr (i, j, k)
            qsz (k) = qs (i, j, k)
            qgz (k) = qg (i, j, k)

            ! dp1: dry air_mass
            ! dp1 (k) = dp1 (k) * (1. - (qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)))
            dp1 (k) = dp1 (k) * (1. - qvz (k)) ! gfs
            omq = dp0 (k) / dp1 (k)

            qvz (k) = qvz (k) * omq
            qlz (k) = qlz (k) * omq
            qrz (k) = qrz (k) * omq
            qiz (k) = qiz (k) * omq
            qsz (k) = qsz (k) * omq
            qgz (k) = qgz (k) * omq
            
            qa0 (k) = qa (i, j, k)
            qaz (k) = qa (i, j, k)
            dz0 (k) = dz (i, j, k)

            den0 (k) = - dp1 (k) / (grav * dz0 (k)) ! density of dry air
            p1 (k) = den0 (k) * rdgas * t0 (k) ! dry air pressure

            ! -----------------------------------------------------------------------
            ! save a copy of old value for computing tendencies
            ! -----------------------------------------------------------------------

            qv0 (k) = qvz (k)
            ql0 (k) = qlz (k)
            qr0 (k) = qrz (k)
            qi0 (k) = qiz (k)
            qs0 (k) = qsz (k)
            qg0 (k) = qgz (k)

            ! -----------------------------------------------------------------------
            ! for sedi_momentum
            ! -----------------------------------------------------------------------

            m1 (k) = 0.
            u0 (k) = uin (i, j, k)
            v0 (k) = vin (i, j, k)
            u1 (k) = u0 (k)
            v1 (k) = v0 (k)

        enddo

        if (do_sedi_w) then
            do k = ktop, kbot
                w1 (k) = w (i, j, k)
            enddo
        endif

        ! 1 minus sigma used to control resolution sensitive parameters
        onemsig = 1.0 - sigma(sqrt(area1(i)))

        ! -----------------------------------------------------------------------
        ! calculate cloud condensation nuclei (ccn)
        ! the following is based on klein eq. 15
        ! -----------------------------------------------------------------------
        if (srf_type(i) < 2.0) then ! exclude snow/ice covered regions
          fac_eis = min(1.0,eis(i)/15.0)**2 ! Estimated inversion strength determine stable regime
          cpaut = c_paut * (0.75*fac_eis + (1.0-fac_eis)) ! scaling autoconversion for stable->unstable
        else
          fac_eis = 0.0
          cpaut = c_paut
        endif
 
        ! ccn needs units #/m^3 
            do k = ktop, kbot
            ! qnl has units # / m^3
            ccn_i (k) = qni (i, j, k)
            ccn_l (k) = qnl (i, j, k)
            c_praut (k) = cpaut * (ccn_l (k) * rhor) ** (- 1. / 3.)
            enddo

        ! -----------------------------------------------------------------------
        ! fix all negative water species
        ! -----------------------------------------------------------------------

        if (fix_negative) &
            call neg_adj (ktop, kbot, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

        m2_rain (i, j, :) = 0.
        m2_sol (i, j, :) = 0.
        revap (i, j, :) = 0.
        isubl (i, j, :) = 0.

        do n = 1, ntimes

            ! -----------------------------------------------------------------------
            ! dry air density
            ! -----------------------------------------------------------------------

            if (p_nonhydro) then
                do k = ktop, kbot
                    dz1 (k) = dz0 (k)
                    den (k) = den0 (k) ! dry air density remains the same
                    denfac (k) = sqrt (sfcrho / den (k))
                enddo
            else
                do k = ktop, kbot
                    dz1 (k) = dz0 (k) * tz (k) / t0 (k) ! hydrostatic balance
                    den (k) = den0 (k) * dz0 (k) / dz1 (k)
                    denfac (k) = sqrt (sfcrho / den (k))
                enddo
            endif

            ! -----------------------------------------------------------------------
            ! sedimentation of cloud ice, snow, and graupel
            ! -----------------------------------------------------------------------

            call fall_speed (ktop, kbot, p1, cnv_fraction(i), anv_icefall, lsc_icefall, &
                             onemsig, den, qsz, qiz, qgz, qlz, tz, vtsz, vtiz, vtgz)

            call terminal_fall (dts, ktop, kbot, tz, qvz, qlz, qrz, qgz, qsz, qiz, &
                dz1, dp1, den, vtgz, vtsz, vtiz, r1, g1, s1, i1, m1_sol, w1)

            rain (i) = rain (i) + r1 ! from melted snow & ice that reached the ground
            snow (i) = snow (i) + s1
            graupel (i) = graupel (i) + g1
            ice (i) = ice (i) + i1

            ! -----------------------------------------------------------------------
            ! heat transportation during sedimentation
            ! -----------------------------------------------------------------------

            if (do_sedi_heat) &
                call sedi_heat (ktop, kbot, dp1, m1_sol, dz1, tz, qvz, qlz, qrz, qiz, &
                    qsz, qgz, c_ice)

            ! -----------------------------------------------------------------------
            ! warm rain processes
            ! -----------------------------------------------------------------------
            
            call warm_rain (dts, ktop, kbot, dp1, dz1, tz, qvz, qlz, qrz, qiz, qsz, &
                qgz, qaz, fac_eis, onemsig, den, denfac, ccn_l, c_praut, vtrz,   &
                r1, evap1, m1_rain, w1, h_var1d)

            rain (i) = rain (i) + r1

            do k = ktop, kbot
                revap (i,j,k) = revap (i,j,k) + evap1(k)
                m2_rain (i, j, k) = m2_rain (i, j, k) + m1_rain (k)
                m2_sol (i, j, k) = m2_sol (i, j, k) + m1_sol (k)
                m1 (k) = m1 (k) + m1_rain (k) + m1_sol (k)
            enddo

            ! -----------------------------------------------------------------------
            ! ice - phase microphysics
            ! -----------------------------------------------------------------------

            call icloud (ktop, kbot, tz, p1, qvz, qlz, qrz, qiz, qsz, qgz, dp1, den, &
                denfac, vtsz, vtgz, vtrz, qaz, dts, subl1, h_var1d, &
                ccn_i, cnv_fraction(i), srf_type(i), onemsig)

            do k = ktop, kbot
                isubl (i,j,k) = isubl (i,j,k) + subl1(k)
            enddo


        enddo ! ntimes

        ! -----------------------------------------------------------------------
        ! momentum transportation during sedimentation
        ! note: dp1 is dry mass; dp0 is the old moist (total) mass
        ! -----------------------------------------------------------------------

        if (sedi_transport) then
            do k = ktop + 1, kbot
                u1 (k) = (dp0 (k) * u1 (k) + m1 (k - 1) * u1 (k - 1)) / (dp0 (k) + m1 (k - 1))
                v1 (k) = (dp0 (k) * v1 (k) + m1 (k - 1) * v1 (k - 1)) / (dp0 (k) + m1 (k - 1))
                u_dt (i, j, k) = u_dt (i, j, k) + (u1 (k) - u0 (k)) * rdt
                v_dt (i, j, k) = v_dt (i, j, k) + (v1 (k) - v0 (k)) * rdt
            enddo
        endif

        if (do_sedi_w) then
            do k = ktop, kbot
                w (i, j, k) = w1 (k)
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! fix all negative water species
        ! -----------------------------------------------------------------------

        if (fix_negative) &
            call neg_adj (ktop, kbot, tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

        ! -----------------------------------------------------------------------
        ! update moist air mass (actually hydrostatic pressure)
        ! convert to dry mixing ratios
        ! -----------------------------------------------------------------------

        do k = ktop, kbot
            omq = dp1 (k) / dp0 (k)
            qv_dt (i, j, k) = qv_dt (i, j, k) + rdt * (qvz (k) - qv0 (k)) * omq
            ql_dt (i, j, k) = ql_dt (i, j, k) + rdt * (qlz (k) - ql0 (k)) * omq
            qr_dt (i, j, k) = qr_dt (i, j, k) + rdt * (qrz (k) - qr0 (k)) * omq
            qi_dt (i, j, k) = qi_dt (i, j, k) + rdt * (qiz (k) - qi0 (k)) * omq
            qs_dt (i, j, k) = qs_dt (i, j, k) + rdt * (qsz (k) - qs0 (k)) * omq
            qg_dt (i, j, k) = qg_dt (i, j, k) + rdt * (qgz (k) - qg0 (k)) * omq
            cvm = c_air + qvz (k) * c_vap + (qrz (k) + qlz (k)) * c_liq + (qiz (k) + qsz (k) + qgz (k)) * c_ice
            pt_dt (i, j, k) = pt_dt (i, j, k) + rdt * (tz (k) - t0 (k)) * cvm / cp_air
        enddo

        ! -----------------------------------------------------------------------
        ! update cloud fraction tendency
        ! -----------------------------------------------------------------------
        if (.not. do_qa) then
        do k = ktop, kbot
             qa_dt (i, j, k) = qa_dt (i, j, k) + rdt * (                          &
                      qa0(k)*SQRT( (qiz(k)+qlz(k)) / max(qi0(k)+ql0(k),qcmin) ) - & ! New Cloud -
                      qa0(k) )                                                      ! Old Cloud
        enddo
        endif

        ! -----------------------------------------------------------------------
        ! fms diagnostics:
        ! -----------------------------------------------------------------------

        ! if (id_cond > 0) then
        ! do k = ktop, kbot ! total condensate
        ! cond (i) = cond (i) + dp1 (k) * (qlz (k) + qrz (k) + qsz (k) + qiz (k) + qgz (k))
        ! enddo
        ! endif
        !
        ! if (id_vtr > 0) then
          do k = ktop, kbot
          vt_r (i, j, k) = vtrz (k)
          enddo
        ! endif
        !
        ! if (id_vts > 0) then
          do k = ktop, kbot
          vt_s (i, j, k) = vtsz (k)
          enddo
        ! endif
        !
        ! if (id_vtg > 0) then
          do k = ktop, kbot
          vt_g (i, j, k) = vtgz (k)
          enddo
        ! endif
        !
        ! if (id_vts > 0) then
          do k = ktop, kbot
          vt_i (i, j, k) = vtiz (k)
          enddo
        ! endif
        !
        ! if (id_droplets > 0) then
        ! do k = ktop, kbot
        ! qn2 (i, j, k) = ccn (k)
        ! enddo
        ! endif

    enddo

end subroutine mpdrv

! -----------------------------------------------------------------------
!> sedimentation of heat
! -----------------------------------------------------------------------

subroutine sedi_heat (ktop, kbot, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)

    implicit none

    ! input q fields are dry mixing ratios, and dm is dry air mass

    integer, intent (in) :: ktop, kbot

    real, intent (in), dimension (ktop:kbot) :: dm, m1, dz, qv, ql, qr, qi, qs, qg

    real, intent (inout), dimension (ktop:kbot) :: tz

    real, intent (in) :: cw ! heat capacity

    real, dimension (ktop:kbot) :: dgz, cvn

    real :: tmp

    integer :: k

    do k = ktop, kbot
        dgz (k) = - 0.5 * grav * dz (k) ! > 0
        cvn (k) = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * &
            c_liq + (qi (k) + qs (k) + qg (k)) * c_ice)
    enddo

    ! -----------------------------------------------------------------------
    ! sjl, july 2014
    ! assumption: the ke in the falling condensates is negligible compared to the potential energy
    ! that was unaccounted for. local thermal equilibrium is assumed, and the loss in pe is transformed
    ! into internal energy (to heat the whole grid box)
    ! backward time - implicit upwind transport scheme:
    ! dm here is dry air mass
    ! -----------------------------------------------------------------------

    k = ktop
    tmp = cvn (k) + m1 (k) * cw
    tz (k) = (tmp * tz (k) + m1 (k) * dgz (k)) / tmp

    ! -----------------------------------------------------------------------
    ! implicit algorithm: can't be vectorized
    ! needs an inner i - loop for vectorization
    ! -----------------------------------------------------------------------

    do k = ktop + 1, kbot
        tz (k) = ((cvn (k) + cw * (m1 (k) - m1 (k - 1))) * tz (k) + m1 (k - 1) * &
            cw * tz (k - 1) + dgz (k) * (m1 (k - 1) + m1 (k))) / (cvn (k) + cw * m1 (k))
    enddo

end subroutine sedi_heat

! -----------------------------------------------------------------------
!> warm rain cloud microphysics
! -----------------------------------------------------------------------

subroutine warm_rain (dt, ktop, kbot, dp, dz, tz, qv, ql, qr, qi, qs, qg, qa, &
        fac_eis, onemsig, &
        den, denfac, ccn, c_praut, vtr, r1, evap1, m1_rain, w1, h_var)

    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: dt !< time step (s)

    real, intent (in), dimension (ktop:kbot) :: h_var

    real, intent (in), dimension (ktop:kbot) :: dp, dz, den
    real, intent (in), dimension (ktop:kbot) :: denfac, ccn, c_praut

    real, intent (in) :: onemsig
    real, intent (in) :: fac_eis !< estimated inversion strength
    
    real, intent (inout), dimension (ktop:kbot) :: tz, vtr
    real, intent (inout), dimension (ktop:kbot) :: qv, ql, qr, qi, qs, qg, qa
    real, intent (inout), dimension (ktop:kbot) :: evap1, m1_rain, w1

    real, intent (out) :: r1

    real, parameter :: so3 = 7. / 3.

    real, dimension (ktop:kbot) :: dl, dm, revap, isubl, qadum
    real, dimension (ktop:kbot + 1) :: ze, zt

    real :: sink, dq, qc
    real :: fac_rc, qden
    real :: zs = 0.
    real :: dt5

    integer :: k

    ! fall velocity constants:

    real, parameter :: vconr = 2503.23638966667
    real, parameter :: normr = 25132741228.7183
    real, parameter :: thr = 1.e-8

    logical :: no_fall

    dt5 = 0.5 * dt

    ! -----------------------------------------------------------------------
    ! terminal speed of rain
    ! -----------------------------------------------------------------------

    evap1 (:) = 0.
    m1_rain (:) = 0.

    call check_column (ktop, kbot, qr, no_fall)
   

    ! -----------------------------------------------------------------------
    ! auto - conversion
    ! assuming linear subgrid vertical distribution of cloud water
    ! following lin et al. 1994, mwr
    ! -----------------------------------------------------------------------

    ! Use In-Cloud condensates
    if (in_cloud) then
      qadum = max(qa,max(qcmin,onemsig))
    else
      qadum = 1.0
    endif
    ql = ql/qadum
    qi = qi/qadum

    fac_rc = rc * (rthreshs*fac_eis + rthreshu*(1.0-fac_eis)) ** 3
 
    if (irain_f /= 0) then
        
        ! -----------------------------------------------------------------------
        ! no subgrid varaibility
        ! -----------------------------------------------------------------------
        
        do k = ktop, kbot
            if (qadum(k) >= onemsig) then
            if (tz (k) > t_wfr) then
                qc = fac_rc * ccn (k) / den (k)
                dq = ql (k) - qc
                if (dq > 0.) then
                    sink = min (dq, dt * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                    sink = min(ql0_max/qadum(k), ql(k), max(0.,sink))
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink*qadum(k)
                    ! new total condensate / old condensate
                    qa (k) = max(0.0,min(1.,qa (k) * max(qadum(k)*(qi (k)+ql (k)-sink),0.0  ) / &
                                                     max(qadum(k)*(qi (k)+ql (k)     ),qcmin) ) )
                endif
            endif
            endif
        enddo
        
    else


        ! -----------------------------------------------------------------------
        ! with subgrid variability
        ! -----------------------------------------------------------------------
        call linear_prof (kbot - ktop + 1, ql (ktop), dl (ktop), z_slope_liq, h_var)

        do k = ktop, kbot
            if (qadum(k) >= onemsig) then
            if (tz (k) > t_wfr + dt_fr) then
                dl (k) = min (max (qcmin, dl (k)), 0.5 * ql (k))
                ! --------------------------------------------------------------------
                ! as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                ! --------------------------------------------------------------------
                qc = fac_rc * ccn (k) / den (k)
                dq = 0.5 * (ql (k) + dl (k) - qc)
                ! --------------------------------------------------------------------
                ! dq = dl if qc == q_minus = ql - dl
                ! dq = 0 if qc == q_plus = ql + dl
                ! --------------------------------------------------------------------
                if (dq > 0.) then ! q_plus > qc
                    ! --------------------------------------------------------------------
                    ! revised continuous form: linearly decays (with subgrid dl) to zero at qc == ql + dl
                    ! --------------------------------------------------------------------
                    sink = min (1., dq / dl (k)) * dt * c_praut (k) * den (k) * exp (so3 * log (ql (k)))
                    sink = min(ql0_max/qadum(k), ql(k), max(0.,sink))
                    ql (k) = ql (k) - sink
                    qr (k) = qr (k) + sink*qadum(k)
                    ! new total condensate / old condensate
                    qa (k) = max(0.0,min(1.,qa (k) * max(qadum(k)*(qi (k)+ql (k)     ),0.0  ) / &
                                                     max(qadum(k)*(qi (k)+ql (k)+sink),qcmin) ) )
                endif
            endif
            endif
        enddo

    endif

    ! Revert In-Cloud condensate
    ql = ql*qadum
    qi = qi*qadum

        ! -----------------------------------------------------------------------
        ! fall speed of rain
        ! -----------------------------------------------------------------------
        
        if (no_fall) then
            vtr (:) = vr_min
        elseif (const_vr) then
            vtr (:) = 0.5*(vr_min+vr_max)
        else
            do k = ktop, kbot
                qden = qr (k) * den (k)
                if (qr (k) < thr) then
                    vtr (k) = vr_min
                else
                    vtr (k) = vr_fac * vconr * sqrt (min (10., sfcrho / den (k))) * &
                        exp (0.2 * log (qden / normr))
                    vtr (k) = min (vr_max, max (vr_min, vtr (k)))
                endif
            enddo
        endif

        ze (kbot + 1) = zs
        do k = kbot, ktop, - 1
            ze (k) = ze (k + 1) - dz (k) ! dz < 0
        enddo

        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the first 1 / 2 time step
        ! -----------------------------------------------------------------------

        call revap_racc (ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)
        evap1 = revap

        if (do_sedi_w) then
            do k = ktop, kbot
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! mass flux induced by falling rain
        ! -----------------------------------------------------------------------
        
        if (no_fall) then
            r1 = 0.0
        elseif (use_ppm) then
            zt (ktop) = ze (ktop)
            do k = ktop + 1, kbot
                zt (k) = ze (k) - dt * (vtr (k - 1) + vtr (k))/2.0
            enddo
            zt (kbot + 1) = zs - dt * vtr (kbot)

            do k = ktop, kbot
                if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            call lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, qr, r1, m1_rain, mono_prof)
        else
            call implicit_fall (dt, ktop, kbot, ze, vtr, dp, qr, r1, m1_rain)
        endif

        ! -----------------------------------------------------------------------
        ! vertical velocity transportation during sedimentation
        ! -----------------------------------------------------------------------

        if (do_sedi_w) then
            w1 (ktop) = (dm (ktop) * w1 (ktop) + m1_rain (ktop) * vtr (ktop)) / (dm (ktop) - m1_rain (ktop))
            do k = ktop + 1, kbot
                w1 (k) = (dm (k) * w1 (k) - m1_rain (k - 1) * vtr (k - 1) + m1_rain (k) * vtr (k)) &
                     / (dm (k) + m1_rain (k - 1) - m1_rain (k))
            enddo
        endif

        ! -----------------------------------------------------------------------
        ! heat transportation during sedimentation
        ! -----------------------------------------------------------------------

        if (do_sedi_heat) &
            call sedi_heat (ktop, kbot, dp, m1_rain, dz, tz, qv, ql, qr, qi, qs, qg, c_liq)

        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the remaing 1 / 2 time step
        ! -----------------------------------------------------------------------

        call revap_racc (ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)
        evap1 = evap1 + revap

end subroutine warm_rain

! -----------------------------------------------------------------------
!> evaporation of rain
! -----------------------------------------------------------------------

subroutine revap_racc (ktop, kbot, dt, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)

    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: dt ! time step (s)
    
    real, intent (in), dimension (ktop:kbot) :: h_var
    real, intent (in), dimension (ktop:kbot) :: den, denfac

    real, intent (inout), dimension (ktop:kbot) :: tz, qv, qr, ql, qi, qs, qg, qa

    real, intent (inout), dimension (ktop:kbot) :: revap

    real, dimension (ktop:kbot) :: lhl, cvm, q_liq, q_sol, lcpk

    real :: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
    real :: qpz, dq, dqh, tin
    real :: fac_revp 
    integer :: k

    revap(:) = 0.

    do k = ktop, kbot
        
        if (tz (k) > t_wfr .and. qr (k) > qpmin) then

           ! timescale efficiency on revap
            fac_revp = 1. - exp (- dt / tau_revp) * qa (k) ! Increase revp when no clouds present (WMP)
            
            ! -----------------------------------------------------------------------
            ! define heat capacity and latent heat coefficient
            ! -----------------------------------------------------------------------

            lhl (k) = lv00 + d0_vap * tz (k)
            q_liq (k) = ql (k) + qr (k)
            q_sol (k) = qi (k) + qs (k) + qg (k)
            cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
            lcpk (k) = lhl (k) / cvm (k)

            tin = tz (k) - lcpk (k) * ql (k) ! presence of clouds suppresses the rain evap
            qpz = qv (k) + ql (k)
            qsat = wqs2 (tin, den (k), dqsdt)
            dqh = max (ql (k), h_var(k) * max (qpz, qcmin))
            dqh = min (dqh, 0.2 * qpz) ! new limiter
            dqv = qsat - qv (k) ! use this to prevent super - sat the grid box
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
                qden = qr (k) * den (k)
                t2 = tin * tin
                evap = crevp (1) * t2 * dq * (crevp (2) * sqrt (qden) + crevp (3) * &
                    exp (0.725 * log (qden))) / (crevp (4) * t2 + crevp (5) * qsat * den (k))
                evap = min (qr (k), dt * fac_revp * evap, dqv / (1. + lcpk (k) * dqsdt))
                qr (k) = qr (k) - evap
                qv (k) = qv (k) + evap
                q_liq (k) = q_liq (k) - evap
                cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tz (k) = tz (k) - evap * lhl (k) / cvm (k)
                revap(k) = evap / dt
            endif

            ! -----------------------------------------------------------------------
            ! accretion: pracc
            ! -----------------------------------------------------------------------
            
            if (qr (k) > qpmin .and. ql (k) > qcmin .and. qsat < q_minus) then
                sink = dt * denfac (k) * cracw * exp (0.95 * log (qr (k) * den (k)))
                sink = ql (k) * sink / (1. + sink)

                ! new total condensate / old condensate
                qa (k) = max(0.0,min(1.,qa (k) * max(qi (k)+ql (k)-sink,0.0  ) / &
                                                 max(qi (k)+ql (k)     ,qcmin) ) )

                ql (k) = ql (k) - sink
                qr (k) = qr (k) + sink
            endif

        endif ! warm - rain
    enddo

end subroutine revap_racc

! -----------------------------------------------------------------------
!> definition of vertical subgrid variability
!! used for cloud ice and cloud water autoconversion
!! qi -- > ql & ql -- > qr
!! edges: qe == qbar + / - dm
! -----------------------------------------------------------------------

subroutine linear_prof (km, q, dm, z_var, h_var)

    implicit none

    integer, intent (in) :: km

    real, intent (in) :: q (km), h_var(km)

    real, intent (out) :: dm (km)

    logical, intent (in) :: z_var

    real :: dq (km)

    integer :: k

    if (z_var) then
        do k = 2, km
            dq (k) = 0.5 * (q (k) - q (k - 1))
        enddo
        dm (1) = 0.

        ! -----------------------------------------------------------------------
        ! use twice the strength of the positive definiteness limiter (lin et al 1994)
        ! -----------------------------------------------------------------------

        do k = 2, km - 1
            dm (k) = 0.5 * min (abs (dq (k) + dq (k + 1)), 0.5 * q (k))
            if (dq (k) * dq (k + 1) <= 0.) then
                if (dq (k) > 0.) then ! local max
                    dm (k) = min (dm (k), dq (k), - dq (k + 1))
                else
                    dm (k) = 0.
                endif
            endif
        enddo
        dm (km) = 0.

        ! -----------------------------------------------------------------------
        ! impose a presumed background horizontal variability that is proportional to the value itself
        ! -----------------------------------------------------------------------

        do k = 1, km
            dm (k) = max (dm (k), qvmin, h_var(k) * q (k))
        enddo
    else
        do k = 1, km
            dm (k) = max (qvmin, h_var(k) * q (k))
        enddo
    endif

end subroutine linear_prof

! =======================================================================
!> ice cloud microphysics processes
!! bulk cloud micro - physics; processes splitting
!! with some un - split sub - grouping
!! time implicit (when possible) accretion and autoconversion
!>@author: Shian-Jiann lin, gfdl
! =======================================================================

subroutine icloud (ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, &
        den, denfac, vts, vtg, vtr, qak, dts, subl1, h_var, ccn, cnv_fraction, srf_type, onemsig)
    
    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in), dimension (ktop:kbot) :: p1, dp1, den, denfac, vts, vtg, vtr
    
    real, intent (inout), dimension (ktop:kbot) :: tzk, qvk, qlk, qrk, qik, qsk, qgk, qak

    real, intent (out), dimension (ktop:kbot) :: subl1

    real, intent (in) :: dts, cnv_fraction, srf_type, onemsig

    real, intent (in), dimension (ktop:kbot) :: h_var, ccn
    
    real, dimension (ktop:kbot) :: lcpk, icpk, tcpk, di, lhl, lhi
    real, dimension (ktop:kbot) :: cvm, q_liq, q_sol

    real :: rdts, fac_i2s, fac_imlt, fac_frz, newice, newliq
    real :: tz, qv, ql, qr, qi, qs, qg, melt, frez, ifrac, newqi, newql
    real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, praci, psaci
    real :: pgmlt, psmlt, pgfr, pgaut, psaut, pgsub
    real :: tc, tsq, dqs0, qden, qim, qsm
    real :: factor, sink, qi_crt
    real :: tmp, qsw, qsi, dqsdt, dq
    real :: dtmp, qc, q_plus, q_minus
    real :: qadum
    real :: critical_qi_factor

    integer :: k, it
    
    rdts = 1. / dts
   
    ! -----------------------------------------------------------------------
    ! define conversion scalar / factor
    ! -----------------------------------------------------------------------

    fac_i2s  = 1. - exp (- dts / tau_i2s)
    fac_imlt = 1. - exp (- dts / tau_imlt)
    fac_frz  = 1. - exp (- dts / tau_frz)
    
    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        lhi (k) = li00 + dc_ice * tzk (k)
        q_liq (k) = qlk (k) + qrk (k)
        q_sol (k) = qik (k) + qsk (k) + qgk (k)
        cvm (k) = c_air + qvk (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
        icpk (k) = lhi (k) / cvm (k)
    enddo

    ! -----------------------------------------------------------------------
    ! sources of cloud ice: pihom, cold rain, and the sat_adj
    ! (initiation plus deposition)
    ! sources of snow: cold rain, auto conversion + accretion (from cloud ice)
    ! sat_adj (deposition; requires pre - existing snow) ; initial snow comes from auto conversion
    ! -----------------------------------------------------------------------

    do k = ktop, kbot

      ! Use In-Cloud condensates
      if (in_cloud) then
        qadum = max(qak (k),max(qcmin,onemsig))
      else    
        qadum = 1.0
      endif
      if (qadum >= onemsig) then

        ql = qlk (k)/qadum
        qi = qik (k)/qadum

        if (tzk (k) > tice .and. qi > qcmin) then
            ! -----------------------------------------------------------------------
            ! pimlt: melting of cloud ice
            ! -----------------------------------------------------------------------
            newliq = new_liq_condensate(tzk (k), ql, qi, cnv_fraction, srf_type)
            melt = fac_imlt * min (qi, newliq, (tzk (k) - tice) / icpk (k))
            tmp = min (melt, dim (ql_mlt/qadum, ql))

            ! new total condensate / old condensate
            qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql-melt+tmp,0.0  ) / &
                                             max(qi+ql         ,qcmin) ) )

            ql      = ql      +         tmp
            qrk (k) = qrk (k) + (melt - tmp)*qadum
            qi      = qi      -  melt
            q_liq (k) = q_liq (k) + melt*qadum
            q_sol (k) = q_sol (k) - melt*qadum
            cvm (k) = c_air + qvk (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
            tzk (k) = tzk (k) - melt*qadum * lhi (k) / cvm (k)
        elseif (tzk (k) <= tice .and. ql > qcmin) then
            ! -----------------------------------------------------------------------
            ! pihom: homogeneous freezing of cloud water into cloud ice
            ! this is the 1st occurance of liquid water freezing in the split mp process
            ! -----------------------------------------------------------------------
            newice = new_ice_condensate(tzk (k), ql, qi, cnv_fraction, srf_type)
            frez = fac_frz * min(ql, newice, ql * (tice - tzk (k)) / icpk (k))
            qi_crt = qi_gen / den (k)
            tmp = min (frez, dim (qi_crt/qadum, qi))

            ! new total condensate / old condensate
            qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql-frez+tmp,0.0  ) / &
                                             max(qi+ql         ,qcmin) ) )

            ql      = ql      -  frez
            qsk (k) = qsk (k) + (frez - tmp)*qadum
            qi      = qi      +         tmp
            q_liq (k) = q_liq (k) - frez*qadum
            q_sol (k) = q_sol (k) + frez*qadum
            cvm (k) = c_air + qvk (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
            tzk (k) = tzk (k) + frez*qadum * lhi (k) / cvm (k)
        endif

        ! Revert In-Cloud condensate
        qlk (k) = ql*qadum
        qik (k) = qi*qadum

      endif

    enddo

    ! -----------------------------------------------------------------------
    ! vertical subgrid variability
    ! -----------------------------------------------------------------------

    call linear_prof (kbot - ktop + 1, qik (ktop), di (ktop), z_slope_ice, h_var)

    ! -----------------------------------------------------------------------
    ! update capacity heat and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        lhl (k) = lv00 + d0_vap * tzk (k)
        lhi (k) = li00 + dc_ice * tzk (k)
        lcpk (k) = lhl (k) / cvm (k)
        icpk (k) = lhi (k) / cvm (k)
        tcpk (k) = lcpk (k) + icpk (k)
    enddo

    do k = ktop, kbot

        ! -----------------------------------------------------------------------
        ! do nothing above p_min
        ! -----------------------------------------------------------------------

        if (p1 (k) < p_min) cycle

        tz = tzk (k)
        qv = qvk (k)
        ql = qlk (k)
        qi = qik (k)
        qr = qrk (k)
        qs = qsk (k)
        qg = qgk (k)

        pgacr = 0.
        pgacw = 0.

        if (tz > tice) then

            ! -----------------------------------------------------------------------
            ! melting of snow
            ! -----------------------------------------------------------------------

            dqs0 = ces0 / p1 (k) - qv

            if (qs > qpmin) then

                ! -----------------------------------------------------------------------
                ! psacw: accretion of cloud water by snow
                ! only rate is used (for snow melt) since tz > tice
                ! -----------------------------------------------------------------------

                if (ql > qcmin) then
                    factor = denfac (k) * csacw * exp (0.8125 * log (qs * den (k)))
                    psacw = factor / (1. + dts * factor) * ql ! rate
                else
                    psacw = 0.
                endif

                ! -----------------------------------------------------------------------
                ! psacr: accretion of rain by melted snow
                ! pracs: accretion of snow by rain
                ! -----------------------------------------------------------------------

                if (qr > qpmin) then
                    psacr = min (acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), &
                        den (k)), qr * rdts)
                    pracs = acr3d (vtr (k), vts (k), qs, qr, cracs, acco (1, 1), den (k))
                else
                    psacr = 0.
                    pracs = 0.
                endif

                ! -----------------------------------------------------------------------
                ! total snow sink:
                ! psmlt: snow melt (due to rain accretion)
                ! -----------------------------------------------------------------------

                tc = tz - tice
                psmlt = max (0., smlt (tc, dqs0, qs * den (k), psacw, psacr, csmlt, &
                    den (k), denfac (k)))
                sink = min (qs, dts * (psmlt + pracs), tc / icpk (k))
                qs = qs - sink
                tmp = min (sink, dim (qs_mlt, ql)) ! max ql due to snow melt

                ! new total condensate / old condensate
                qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql+tmp,0.0  ) / &
                                                 max(qi+ql    ,qcmin) ) )

                ql = ql + tmp
                qr = qr + sink - tmp
                q_liq (k) = q_liq (k) + sink
                q_sol (k) = q_sol (k) - sink
                cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tz = tz - sink * lhi (k) / cvm (k)

            endif

            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------

            lhi (k) = li00 + dc_ice * tz
            icpk (k) = lhi (k) / cvm (k)

            ! -----------------------------------------------------------------------
            ! melting of graupel
            ! -----------------------------------------------------------------------

            if (qg > qpmin .and. tz > tice) then

                ! -----------------------------------------------------------------------
                ! pgacr: accretion of rain by graupel
                ! -----------------------------------------------------------------------

                if (qr > qpmin) &
                    pgacr = min (acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
                        den (k)), rdts * qr)

                ! -----------------------------------------------------------------------
                ! pgacw: accretion of cloud water by graupel
                ! -----------------------------------------------------------------------

                qden = qg * den (k)
                if (ql > qcmin) then
                    factor = cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                    pgacw = factor / (1. + dts * factor) * ql ! rate
                endif

                ! -----------------------------------------------------------------------
                ! pgmlt: graupel melt
                ! -----------------------------------------------------------------------

                tc = tz - tice
                pgmlt = dts * gmlt (tc, dqs0, qg * den (k), pgacw, pgacr, cgmlt, den (k))
                pgmlt = min (max (0., pgmlt), qg, tc / icpk (k))
                qg = qg - pgmlt
                qr = qr + pgmlt
                q_liq (k) = q_liq (k) + pgmlt
                q_sol (k) = q_sol (k) - pgmlt
                cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tz = tz - pgmlt * lhi (k) / cvm (k)

            endif

        else

            ! -----------------------------------------------------------------------
            ! cloud ice proc:
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! psaci: accretion of cloud ice by snow
            ! -----------------------------------------------------------------------

            if (qi > 3.e-7 .and. tz < tice) then ! cloud ice sink terms

                tc = tz - tice

                if (qs > qpmin) then
                    ! -----------------------------------------------------------------------
                    ! sjl added (following lin eq. 23) the temperature dependency
                    ! to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
                    ! -----------------------------------------------------------------------
                    factor = dts * denfac (k) * csaci * exp (0.05 * tc + 0.8125 * log (qs * den (k)))
                    psaci = qi * factor / (1. + factor)
                else
                    psaci = 0.
                endif

                ! -----------------------------------------------------------------------
                ! psaut: autoconversion: cloud ice -- > snow
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! similar to lfo 1983: eq. 21 solved implicitly
                ! threshold from wsm6 scheme, hong et al 2004, eq (13)
                ! slight decrease in critical_qi_factor at warmer temps
                ! -----------------------------------------------------------------------

                ! qi0_crt (ice to snow conversion) has strong resolution dependence
                !    account for this using onemsig to convert more ice to snow at coarser resolutions
                critical_qi_factor = qi0_crt*(onemsig + 1.e-1*(1.0-onemsig))
             
                ! use ice_fraction to convert more ice to snow closer to 0-cel 
                qim = ice_fraction(tz, cnv_fraction, srf_type) * critical_qi_factor / den (k)
 
                ! -----------------------------------------------------------------------
                ! assuming linear subgrid vertical distribution of cloud ice
                ! the mismatch computation following lin et al. 1994, mwr
                ! -----------------------------------------------------------------------

                di (k) = max (di (k), qcmin)
                q_plus = qi + di (k)
                if (q_plus > (qim + qcmin)) then
                    if (qim > (qi - di (k))) then
                        dq = (0.25 * (q_plus - qim) ** 2) / di (k)
                    else
                        dq = qi - qim
                    endif
                    psaut = fac_i2s * exp (0.025 * tc) * dq
                else
                    psaut = 0.
                endif
                sink = min (qi, psaci + psaut)

                ! new total condensate / old condensate
                qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql-sink,0.0  ) / &
                                                 max(qi+ql     ,qcmin) ) )

                qi = qi - sink
                qs = qs + sink

                ! -----------------------------------------------------------------------
                ! pgaci: accretion of cloud ice by graupel
                ! -----------------------------------------------------------------------

                if (qg > qpmin) then
                    ! -----------------------------------------------------------------------
                    ! factor = dts * cgaci / sqrt (den (k)) * exp (0.05 * tc + 0.875 * log (qg * den (k)))
                    ! simplified form: remove temp dependency & set the exponent "0.875" -- > 1
                    ! -----------------------------------------------------------------------
                    factor = dts * cgaci * sqrt (den (k)) * qg
                    pgaci = factor / (1. + factor) * qi

                    ! new total condensate / old condensate
                    qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql-pgaci,0.0  ) / &
                                                     max(qi+ql      ,qcmin) ) )

                    qi = qi - pgaci
                    qg = qg + pgaci
                endif

            endif

            ! -----------------------------------------------------------------------
            ! cold - rain proc:
            ! -----------------------------------------------------------------------

            ! -----------------------------------------------------------------------
            ! rain to ice, snow, graupel processes:
            ! -----------------------------------------------------------------------

            if (qr > qpmin .and. tz < tice) then

                tc = tz - tice

                ! -----------------------------------------------------------------------
                ! * sink * terms to qr: psacr + pgfr
                ! source terms to qs: psacr
                ! source terms to qg: pgfr
                ! -----------------------------------------------------------------------

                ! -----------------------------------------------------------------------
                ! psacr accretion of rain by snow
                ! -----------------------------------------------------------------------

                if (qs > qpmin) then ! if snow exists
                    psacr = dts * acr3d (vts (k), vtr (k), qr, qs, csacr, acco (1, 2), den (k))
                else
                    psacr = 0.
                endif

                ! -----------------------------------------------------------------------
                ! pgfr: rain freezing -- > graupel
                ! -----------------------------------------------------------------------

                pgfr = dts * cgfr (1) / den (k) * (exp (- cgfr (2) * tc) - 1.) * &
                    exp (1.75 * log (qr * den (k)))

                ! -----------------------------------------------------------------------
                ! total sink to qr
                ! -----------------------------------------------------------------------

                sink = psacr + pgfr
                factor = min (sink, qr, - tc / icpk (k)) / max (sink, qpmin)

                psacr = factor * psacr
                pgfr = factor * pgfr

                sink = psacr + pgfr
                qr = qr - sink
                qs = qs + psacr
                qg = qg + pgfr
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink
                cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tz = tz + sink * lhi (k) / cvm (k)

            endif

            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------

            lhi (k) = li00 + dc_ice * tz
            icpk (k) = lhi (k) / cvm (k)

            ! -----------------------------------------------------------------------
            ! graupel production terms:
            ! -----------------------------------------------------------------------

            if (qs > qpmin .and. tz < tice0) then

                ! -----------------------------------------------------------------------
                ! accretion: snow -- > graupel
                ! -----------------------------------------------------------------------

                if (qg > qpmin) then
                    sink = dts * acr3d (vtg (k), vts (k), qs, qg, cgacs, acco (1, 4), den (k))
                else
                    sink = 0.
                endif

                ! -----------------------------------------------------------------------
                ! autoconversion snow -- > graupel
                ! -----------------------------------------------------------------------

                qsm = qs0_crt / den (k)
                if (qs > qsm) then
                    factor = dts * 1.e-3 * exp (0.09 * (tz - tice))
                    sink = sink + (qs - qsm) * factor / (1. + factor)
                endif
                sink = min (qs, sink)
                qs = qs - sink
                qg = qg + sink

            endif ! snow existed

            if (qg > qpmin .and. tz < tice0) then

                ! -----------------------------------------------------------------------
                ! pgacw: accretion of cloud water by graupel
                ! -----------------------------------------------------------------------

                if (ql > qcmin) then
                    qden = qg * den (k)
                    factor = dts * cgacw * qden / sqrt (den (k) * sqrt (sqrt (qden)))
                    pgacw = factor / (1. + factor) * ql
                else
                    pgacw = 0.
                endif

                ! -----------------------------------------------------------------------
                ! pgacr: accretion of rain by graupel
                ! -----------------------------------------------------------------------

                if (qr > qpmin) then
                    pgacr = min (dts * acr3d (vtg (k), vtr (k), qr, qg, cgacr, acco (1, 3), &
                        den (k)), qr)
                else
                    pgacr = 0.
                endif

                sink = pgacr + pgacw
                factor = min (sink, dim (tice, tz) / icpk (k)) / max (sink, qpmin)
                pgacr = factor * pgacr
                pgacw = factor * pgacw

                sink = pgacr + pgacw
                qg = qg + sink
                qr = qr - pgacr
                ql = ql - pgacw
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink
                cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tz = tz + sink * lhi (k) / cvm (k)

            endif

        endif

        tzk (k) = tz
        qvk (k) = qv
        qlk (k) = ql
        qik (k) = qi
        qrk (k) = qr
        qsk (k) = qs
        qgk (k) = qg

    enddo

    ! -----------------------------------------------------------------------
    ! subgrid cloud microphysics
    ! -----------------------------------------------------------------------

    call subgrid_z_proc (ktop, kbot, p1, den, denfac, dts, tzk, qvk, &
        qlk, qrk, qik, qsk, qgk, qak, subl1, h_var, ccn, cnv_fraction, srf_type, onemsig)

end subroutine icloud

! =======================================================================
!>temperature sensitive high vertical resolution processes
! =======================================================================

subroutine subgrid_z_proc (ktop, kbot, p1, den, denfac, dts, tz, qv, &
    ql, qr, qi, qs, qg, qa, subl1, h_var, ccn, cnv_fraction, srf_type, onemsig)
    
    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in), dimension (ktop:kbot) :: p1, den, denfac

    real, intent (in) :: dts, cnv_fraction, srf_type, onemsig

    real, intent (in), dimension (ktop:kbot) :: h_var, ccn

    real, intent (inout), dimension (ktop:kbot) :: tz, qv, ql, qr, qi, qs, qg, qa

    real, intent (out), dimension (ktop:kbot) :: subl1

    real, dimension (ktop:kbot) :: lcpk, icpk, tcpk, tcp3, lhl, lhi
    real, dimension (ktop:kbot) :: cvm, q_liq, q_sol, q_cond

    real :: fac_l2v, fac_i2v

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
        
      rh_adj =             1. - h_var(k) - rh_inc
      rh_rain = max (0.35, 1. - h_var(k) - rh_inr)

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
            qa (k) = 1. ! air fully saturated; 100 % cloud cover
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
        ! cloud water < -- > vapor adjustment: LS evaporation
        ! -----------------------------------------------------------------------
        if (do_evap) then
           evap = 0.0
           subl = 0.0
           tin = tz (k) - (lhl (k) * (ql (k) + qi (k)) + lhi (k) * qi (k)) / (c_air + &
                 qpz * c_vap + qr (k) * c_liq + (qs (k) + qg (k)) * c_ice)
           if (tin > t_sub + 6.) then
             qpz = qv (k) + ql (k) + qi (k)
             rh = qpz / iqs1 (tin, den (k))
             if (rh < rh_adj) then
                ! instant evap of all liquid
                evap = ql(k)
                ! instant subl of all ice
                subl = ql(k)
             else
                ! partial evap of liquid
                qsw = wqs2 (tz (k), den (k), dwsdt)
                dq0 = qsw - qv (k)
                if (dq0 > qvmin) then
                   factor = min (1., fac_l2v * (10. * dq0 / qsw))
                   evap = min (ql (k), factor * ql(k) / (1. + tcp3 (k) * dwsdt))
                endif
              endif
           endif
 
           evap = evap*onemsig ! resolution dependent evap 0:1 coarse:fine
           subl = subl*onemsig ! resolution dependent subl 0:1 coarse:fine

           ! new total condensate / old condensate
           qa(k) = max(0.0,min(1.,qa(k) * max(qi(k)+ql(k)-evap-subl,0.0  ) / &
                                          max(qi(k)+ql(k)          ,qcmin) ) )
           qv (k) = qv (k) + evap + subl
           ql (k) = ql (k) - evap
           qi (k) = qi (k) - subl
           q_liq (k) = q_liq (k) - evap
           q_sol (k) = q_sol (k) - subl
           cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
           tz (k) = tz (k) - evap * lhl (k) / cvm (k) - subl * lhi (k) / cvm (k)
        endif

        ! -----------------------------------------------------------------------
        ! update heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------

        lhi (k) = li00 + dc_ice * tz (k)
        icpk (k) = lhi (k) / cvm (k)

        ! -----------------------------------------------------------------------
        ! enforce complete freezing when ice_fraction==1
        ! -----------------------------------------------------------------------

        ifrac = ice_fraction(tz (k),cnv_fraction,srf_type)
        if (ifrac == 1. .and. ql (k) > qcmin) then
            sink = ql (k)
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
                 tmp = tice - tz (k)
              qi_crt = 4.92e-11 * exp (1.33 * log (1.e3 * exp (0.1 * tmp)))
              qi_crt = max (qi_crt, 1.82e-6) * min (qi_lim, 0.1 * tmp) / den (k)
                sink = min (sink, max (qi_crt - qi (k), pidep), tmp / tcpk (k))
            else ! ice -- > vapor
              ! sublimation
              if (do_subl) then
                pidep = pidep * min (1., dim (tz (k), t_sub) * 0.2)
                sink = fac_i2v * max (pidep, sink, - qi (k))
                subl1(k) = subl1(k) - sink / dts
              else
                sink = 0.
              endif
              sink = sink*onemsig ! resolution dependent subl 0:1 coarse:fine
            endif
            ! new total condensate / old condensate
            qa(k) = max(0.0,min(1.,qa(k) * max(qi(k)+ql(k)+sink,0.0  ) / &
                                           max(qi(k)+ql(k)     ,qcmin) ) )
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
                subl1(k) = subl1(k) + pssub / dts
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

subroutine terminal_fall (dtm, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
        den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1)

    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: dtm ! time step (s)

    real, intent (in), dimension (ktop:kbot) :: vtg, vts, vti, den, dp, dz

    real, intent (inout), dimension (ktop:kbot) :: qv, ql, qr, qg, qs, qi, tz, m1_sol, w1

    real, intent (out) :: r1, g1, s1, i1

    real, dimension (ktop:kbot + 1) :: ze, zt
    
    real :: qsat, dqsdt, evap, dtime
    real :: factor, frac
    real :: tmp, precip, tc, sink

    real, dimension (ktop:kbot) :: lcpk, icpk, cvm, q_liq, q_sol, lhl, lhi
    real, dimension (ktop:kbot) :: m1, dm

    real :: zs = 0.
    real :: fac_imlt

    integer :: k, k0, m

    logical :: no_fall
    
    fac_imlt = 1. - exp (- dtm / tau_imlt)
    
    ! -----------------------------------------------------------------------
    ! define heat capacity and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        m1_sol (k) = 0.
        lhl (k) = lv00 + d0_vap * tz (k)
        lhi (k) = li00 + dc_ice * tz (k)
        q_liq (k) = ql (k) + qr (k)
        q_sol (k) = qi (k) + qs (k) + qg (k)
        cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
        lcpk (k) = lhl (k) / cvm (k)
        icpk (k) = lhi (k) / cvm (k)
    enddo

    ! -----------------------------------------------------------------------
    ! find significant melting level
    ! -----------------------------------------------------------------------

    k0 = kbot
    do k = ktop, kbot - 1
        if (tz (k) > tice) then
            k0 = k
            exit
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! melting of cloud_ice (before fall) :
    ! -----------------------------------------------------------------------

    do k = k0, kbot
        tc = tz (k) - tice
        if (qi (k) > qcmin .and. tc > 0.) then
            sink = min (qi (k), fac_imlt * tc / icpk (k))
            tmp = min (sink, dim (ql_mlt, ql (k)))
            ql (k) = ql (k) + tmp
            qr (k) = qr (k) + sink - tmp
            qi (k) = qi (k) - sink
            q_liq (k) = q_liq (k) + sink
            q_sol (k) = q_sol (k) - sink
            cvm (k) = c_air + qv (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
            tz (k) = tz (k) - sink * lhi (k) / cvm (k)
            tc = tz (k) - tice
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! turn off melting when cloud microphysics time step is small
    ! -----------------------------------------------------------------------

    if (dtm < 300.) k0 = kbot

    ze (kbot + 1) = zs
    do k = kbot, ktop, - 1
        ze (k) = ze (k + 1) - dz (k) ! dz < 0
    enddo

    zt (ktop) = ze (ktop)

    ! -----------------------------------------------------------------------
    ! update capacity heat and latend heat coefficient
    ! -----------------------------------------------------------------------

    do k = k0, kbot
        lhi (k) = li00 + dc_ice * tz (k)
        icpk (k) = lhi (k) / cvm (k)
    enddo

    ! -----------------------------------------------------------------------
    ! melting of falling cloud ice into rain
    ! -----------------------------------------------------------------------

    call check_column (ktop, kbot, qi, no_fall)

    if (vi_min < 1.e-5 .or. no_fall) then
        i1 = 0.
    else

        do k = ktop + 1, kbot
            zt (k) = ze (k) - dtm * (vti (k - 1) + vti (k))/2.0
        enddo
        zt (kbot + 1) = zs - dtm * vti (kbot)

        do k = ktop, kbot
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < kbot) then
            do k = kbot - 1, k0, - 1
                if (qi (k) > qcmin) then
                    do m = k + 1, kbot
                        if (zt (k + 1) >= ze (m)) exit
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (dtm, (ze (m) - ze (m + 1)) / vti (k))
                            dtime = min (1.0, dtime / tau_imlt)
                            sink = min (qi (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tmp = min (sink, dim (ql_mlt, ql (m)))
                            ql (m) = ql (m) + tmp
                            qr (m) = qr (m) - tmp + sink
                            tz (m) = tz (m) - sink * icpk (m)
                            qi (k) = qi (k) - sink * dp (m) / dp (k)
                        endif
                    enddo
                endif
            enddo
        endif

        if (do_sedi_w) then
            do k = ktop, kbot
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, qi, i1, m1_sol, mono_prof)
        else
            call implicit_fall (dtm, ktop, kbot, ze, vti, dp, qi, i1, m1_sol)
        endif

        if (do_sedi_w) then
            w1 (ktop) = (dm (ktop) * w1 (ktop) + m1_sol (ktop) * vti (ktop)) / (dm (ktop) - m1_sol (ktop))
            do k = ktop + 1, kbot
                w1 (k) = (dm (k) * w1 (k) - m1_sol (k - 1) * vti (k - 1) + m1_sol (k) * vti (k)) &
                     / (dm (k) + m1_sol (k - 1) - m1_sol (k))
            enddo
        endif

    endif

    ! -----------------------------------------------------------------------
    ! melting of falling snow into rain
    ! -----------------------------------------------------------------------

    r1 = 0.

    call check_column (ktop, kbot, qs, no_fall)

    if (no_fall) then
        s1 = 0.
    else

        do k = ktop + 1, kbot
            zt (k) = ze (k) - dtm * (vts (k - 1) + vts (k))/2.0
        enddo
        zt (kbot + 1) = zs - dtm * vts (kbot)

        do k = ktop, kbot
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < kbot) then
            do k = kbot - 1, k0, - 1
                if (qs (k) > qpmin) then
                    do m = k + 1, kbot
                        if (zt (k + 1) >= ze (m)) exit
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (dtm, (ze (m) - ze (m + 1)) / vts (k))
                            dtime = min (1.0, dtime / tau_smlt)
                            sink = min (qs (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qs (k) = qs (k) - sink * dp (m) / dp (k)
                            if (zt (k) < zs) then
                                r1 = r1 + sink * dp (m) ! precip as rain
                            else
                                ! qr source here will fall next time step (therefore, can evap)
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qs (k) < qpmin) exit
                    enddo
                endif
            enddo
        endif

        if (do_sedi_w) then
            do k = ktop, kbot
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, qs, s1, m1, mono_prof)
        else
            call implicit_fall (dtm, ktop, kbot, ze, vts, dp, qs, s1, m1)
        endif

        do k = ktop, kbot
            m1_sol (k) = m1_sol (k) + m1 (k)
        enddo

        if (do_sedi_w) then
            w1 (ktop) = (dm (ktop) * w1 (ktop) + m1 (ktop) * vts (ktop)) / (dm (ktop) - m1 (ktop))
            do k = ktop + 1, kbot
                w1 (k) = (dm (k) * w1 (k) - m1 (k - 1) * vts (k - 1) + m1 (k) * vts (k)) &
                     / (dm (k) + m1 (k - 1) - m1 (k))
            enddo
        endif

    endif

    ! ----------------------------------------------
    ! melting of falling graupel into rain
    ! ----------------------------------------------

    call check_column (ktop, kbot, qg, no_fall)

    if (no_fall) then
        g1 = 0.
    else

        do k = ktop + 1, kbot
            zt (k) = ze (k) - dtm * (vtg (k - 1) + vtg (k))/2.0
        enddo
        zt (kbot + 1) = zs - dtm * vtg (kbot)

        do k = ktop, kbot
            if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
        enddo

        if (k0 < kbot) then
            do k = kbot - 1, k0, - 1
                if (qg (k) > qpmin) then
                    do m = k + 1, kbot
                        if (zt (k + 1) >= ze (m)) exit
                        if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                            dtime = min (dtm, (ze (m) - ze (m + 1)) / vtg (k))
                            dtime = min (1.0, dtime / tau_gmlt)
                            sink = min (qg (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                            tz (m) = tz (m) - sink * icpk (m)
                            qg (k) = qg (k) - sink * dp (m) / dp (k)
                            if (zt (k) < zs) then
                                r1 = r1 + sink * dp (m)
                            else
                                qr (m) = qr (m) + sink
                            endif
                        endif
                        if (qg (k) < qpmin) exit
                    enddo
                endif
            enddo
        endif

        if (do_sedi_w) then
            do k = ktop, kbot
                dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
            enddo
        endif

        if (use_ppm) then
            call lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, qg, g1, m1, mono_prof)
        else
            call implicit_fall (dtm, ktop, kbot, ze, vtg, dp, qg, g1, m1)
        endif

        do k = ktop, kbot
            m1_sol (k) = m1_sol (k) + m1 (k)
        enddo

        if (do_sedi_w) then
            w1 (ktop) = (dm (ktop) * w1 (ktop) + m1 (ktop) * vtg (ktop)) / (dm (ktop) - m1 (ktop))
            do k = ktop + 1, kbot
                w1 (k) = (dm (k) * w1 (k) - m1 (k - 1) * vtg (k - 1) + m1 (k) * vtg (k)) &
                     / (dm (k) + m1 (k - 1) - m1 (k))
            enddo
        endif

    endif

end subroutine terminal_fall

! =======================================================================
!>@brief The subroutine 'check_column' checks
!!       if the water species is large enough to fall.
! =======================================================================

subroutine check_column (ktop, kbot, q, no_fall)

    implicit none

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
!>@brief The subroutine 'implicit_fall' computes the time-implicit monotonic
!! scheme.
!>@author Shian-Jiann Lin, 2016
! =======================================================================

subroutine implicit_fall (dt, ktop, kbot, ze, vt, dp, q, precip, m1)

    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: dt

    real, intent (in), dimension (ktop:kbot + 1) :: ze

    real, intent (in), dimension (ktop:kbot) :: vt, dp

    real, intent (inout), dimension (ktop:kbot) :: q

    real, intent (out), dimension (ktop:kbot) :: m1

    real, intent (out) :: precip

    real, dimension (ktop:kbot) :: dz, qm, dd

    integer :: k

    do k = ktop, kbot
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dt * vt (k)
        q (k) = q (k) * dp (k)
    enddo

    ! -----------------------------------------------------------------------
    ! sedimentation: non - vectorizable loop
    ! -----------------------------------------------------------------------

    qm (ktop) = q (ktop) / (dz (ktop) + dd (ktop))
    do k = ktop + 1, kbot
        qm (k) = (q (k) + dd (k - 1) * qm (k - 1)) / (dz (k) + dd (k))
    enddo

    ! -----------------------------------------------------------------------
    ! qm is density at this stage
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        qm (k) = qm (k) * dz (k)
    enddo

    ! -----------------------------------------------------------------------
    ! output mass fluxes: non - vectorizable loop
    ! -----------------------------------------------------------------------

    m1 (ktop) = q (ktop) - qm (ktop)
    do k = ktop + 1, kbot
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (kbot)

    ! -----------------------------------------------------------------------
    ! update:
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        q (k) = qm (k) / dp (k)
    enddo

end subroutine implicit_fall

! =======================================================================
!> lagrangian scheme
!  developed by sj lin, ????
! =======================================================================

subroutine lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, q, precip, m1, mono)

    implicit none

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
        555 continue
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
                       onemsig, den, qs, qi, qg, ql, tk, vts, vti, vtg)

    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in) :: cnv_fraction, anv_icefall, lsc_icefall, onemsig
    real, intent (in), dimension (ktop:kbot) :: pl, den, qs, qi, qg, ql, tk
    real, intent (out), dimension (ktop:kbot) :: vts, vti, vtg

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

    real, dimension (ktop:kbot) :: rhof

    real :: tc
    real :: zero=0.0
    real :: viCNV, viLSC, IWC
    real :: rBB, C0, C1, DIAM, lnP
    integer :: k

    ! -----------------------------------------------------------------------
    ! marshall - palmer formula
    ! -----------------------------------------------------------------------

    ! -----------------------------------------------------------------------
    ! try the local air density -- for global model; the true value could be
    ! much smaller than sfcrho over high mountains
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        rhof (k) = sqrt (min (10., sfcrho / den (k)))
    enddo

    ! -----------------------------------------------------------------------
    ! ice:
    ! -----------------------------------------------------------------------

    if (const_vi) then
        vti (:) = 0.5*(vi_min+vi_max)
    else
        do k = ktop, kbot
            if (qi (k) < thi) then
                vti (k) = vi_min
            else
                tc     = tk (k) - tice ! deg C
                IWC    = qi (k) * den (k) * 1.e3 ! Units are g/m3

               if (ICE_LSC_VFALL_PARAM == 1) then
               ! -----------------------------------------------------------------------
               ! use deng and mace (2008, grl)
               ! https://doi.org/10.1029/2008GL035054
               ! -----------------------------------------------------------------------
                viLSC   = lsc_icefall*10.0**(log10(IWC) * (tc * (aaL * tc + bbL) + ccL) + ddL * tc + eeL)
               else
               ! -----------------------------------------------------------------------
               ! use Mishra et al (2014, JGR) 'Parameterization of ice fall speeds in
               !                               ice clouds: Results from SPartICus'
               ! -----------------------------------------------------------------------
                viLSC  = MAX(10.0,lsc_icefall*(1.411*tc + 11.71*log10(IWC*1.e3) + 82.35))
               endif
                
               if (ICE_CNV_VFALL_PARAM == 1) then
               ! -----------------------------------------------------------------------
               ! use deng and mace (2008, grl)
               ! https://doi.org/10.1029/2008GL035054
               ! -----------------------------------------------------------------------
                viCNV   = anv_icefall*10.0**(log10(IWC) * (tc * (aaC * tc + bbC) + ccC) + ddC * tc + eeC)
               else
               ! -----------------------------------------------------------------------
               ! use Mishra et al (2014, JGR) 'Parameterization of ice fall speeds in
               !                               ice clouds: Results from SPartICus'
               ! -----------------------------------------------------------------------
                viCNV  = MAX(10.0,anv_icefall*(1.119*tc + 14.21*log10(IWC*1.e3) + 68.85))
               endif

               ! Combine
                vti (k) = viLSC*(1.0-cnv_fraction) + viCNV*(cnv_fraction)

                if (do_icepsettle) then
                  ! Include pressure sensitivity (eq 14 in https://doi.org/10.1175/JAS-D-12-0124.1)
                  DIAM = 2.0*LDRADIUS4(pl(k)/100.0,tk(k),qi(k),zero,zero,2)*1.e6 ! microns
                  lnP   = log(pl(k)/100.0)
                  C0    = -1.04 + 0.298*lnP
                  C1    =  0.67 - 0.097*lnP
                  ! apply pressure scaling
                  vti (k) = vti (k) * (C0 + C1*log(DIAM))
                endif

               ! Update units from cm/s to m/s
                vti (k) = 0.01 * vi_fac * vti (k)

               ! Limits
                vti (k) = min (vi_max, max (vi_min, vti (k)))
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! snow:
    ! -----------------------------------------------------------------------

    if (const_vs) then
        vts (:) = 0.5*(vs_min+vs_max)
    else
        do k = ktop, kbot
            if (qs (k) < ths) then
                vts (k) = vs_min
            else
                vts (k) = vs_fac * vcons * rhof (k) * exp (0.0625 * log (qs (k) * den (k) / norms))
                vts (k) = min (vs_max, max (vs_min, vts (k)))
            endif
        enddo
    endif

    ! -----------------------------------------------------------------------
    ! graupel:
    ! -----------------------------------------------------------------------

    if (const_vg) then
        vtg (:) = 0.5*(vg_min+vg_max)
    else
        do k = ktop, kbot
            if (qg (k) < thg) then
                vtg (k) = vg_min
            else
                vtg (k) = vg_fac * vcong * rhof (k) * sqrt (sqrt (sqrt (qg (k) * den (k) / normg)))
                vtg (k) = min (vg_max, max (vg_min, vtg (k)))
            endif
        enddo
    endif

end subroutine fall_speed

! =======================================================================
!>@brief The subroutine 'setup'm' sets up
!! gfdl cloud microphysics parameters.
! =======================================================================

subroutine setupm

    implicit none

    real :: cd, scm3, pisq, act (8)
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

    c_paut = c_paut * 0.104 * grav / 1.717e-5

    cracs = pisq * rnzr * rnzs * rhos
    csacr = pisq * rnzr * rnzs * rhor
    cgacr = pisq * rnzr * rnzg * rhor
    cgacs = pisq * rnzg * rnzs * rhos
    cgacs = cgacs * c_pgacs

  ! act: 1 - 2:racs (s - r)
  !      3 - 4:sacr (r - s)
  !      5 - 6:gacr (r - g)
  !      7 - 8:gacs (s - g)

    act (1) = pie * rnzs * rhos
    act (2) = pie * rnzr * rhor
    act (3) = act (2)
    act (4) = act (1)
    act (5) = act (2)
    act (6) = pie * rnzg * rhog
    act (7) = act (1)
    act (8) = act (6)

    do i = 1, 3
        do k = 1, 4
            acco (i, k) = acc (i) / (act (2 * k - 1) ** ((7 - i) * 0.25) * act (2 * k) ** (i * 0.25))
        enddo
    enddo

    ! decreasing alin will reduce accretion of rain from cloud ice/water
    craci = pie * rnzr * alin * gam380 / (4. * act (2) ** 0.95)
    cracw = c_cracw * craci

    ! decreasing clin will reduce accretion of snow from cloud water/ice
    csacw = pie * rnzs * clin * gam325 / (4. * act (1) ** 0.8125)
    csaci = c_psaci * csacw

    ! decreasing gcon will reduce accretion of graupel from cloud ice/water
    cgacw = pie * rnzg * gcon * gam350 / (4. * act (6) ** 0.875)
    cgaci = c_pgaci * cgacw

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
    cgsub (4) = cssub (4)
    crevp (4) = cssub (4)
    cssub (5) = hlts ** 2 * vdifu
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

end subroutine setupm

! =======================================================================
! initialization of gfdl cloud microphysics
!>@brief The subroutine 'gfdl_cloud_microphys_init' initializes the GFDL
!! cloud microphysics.
! =======================================================================

subroutine gfdl_cloud_microphys_init ()
    implicit none
    integer :: nlunit
    character (len = 64) :: fn_nml = 'input.nml'

    integer :: ios, ierr
    logical :: exists

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml = gfdl_cloud_microphysics_nml)
#else
    inquire (file = trim (fn_nml), exist = exists)
    if (.not. exists) then
        write (6, *) 'gfdl - mp :: namelist file: ', trim (fn_nml), ' does not exist'
        stop
    else
        !nlunit=open_namelist_file()
        !rewind (nlunit)
        open(NEWUNIT=nlunit,file=trim(fn_nml), form='formatted',access='sequential',iostat=ios)
        if(ios /= 0) stop 'open namelist file gfdl_cloud_microphys_init failed, bailing out...'
        rewind (nlunit, iostat=ios)
        if(ios /= 0) stop 'rewind namelist file gfdl_cloud_microphys_init failed, bailing out...'
     ! Read Main namelist
        read (nlunit,gfdl_cloud_microphysics_nml,iostat=ios)
        if(ios /= 0) stop 'read namelist gfdl_cloud_microphys_init failed, bailing out...'
        !ierr = check_nml_error(ios,'gfdl_cloud_microphysics_nml')
        !call close_file(nlunit)
        close(nlunit, iostat=ios)
        if(ios /= 0) stop 'close namelist file gfdl_cloud_microphys_init failed, bailing out...'
    endif
#endif

    if (MAPL_AM_I_ROOT()) then
        write (*, *) " ================================================================== "
        write (*, *) "gfdl_cloud_microphys_mod"
        write (*, nml = gfdl_cloud_microphysics_nml)
        write (*, *) " ================================================================== "
    endif

    if (do_setup) then
        call setup_con
        call setupm
        do_setup = .false.
    endif

    module_is_initialized = .true.

end subroutine gfdl_cloud_microphys_init

! =======================================================================
! end of gfdl cloud microphysics
!>@brief The subroutine 'gfdl_cloud_microphys_init' terminates the GFDL
!! cloud microphysics.
! =======================================================================

subroutine gfdl_cloud_microphys_end

    implicit none

    deallocate (table)
    deallocate (table2)
    deallocate (table3)
    deallocate (tablew)
    deallocate (des)
    deallocate (des2)
    deallocate (des3)
    deallocate (desw)

    tables_are_initialized = .false.

end subroutine gfdl_cloud_microphys_end

! =======================================================================
! qsmith table initialization
!>@brief The subroutine 'setup_con' sets up constants and calls 'qsmith_init'.
! =======================================================================

subroutine setup_con

    implicit none

    rgrav = 1. / grav

    if (.not. qsmith_tables_initialized) call qsmith_init

    qsmith_tables_initialized = .true.

end subroutine setup_con

! =======================================================================
!>@brief The function 'acr3d' is an accretion function (lin et al. 1983)
! =======================================================================

real function acr3d (v1, v2, q1, q2, c, cac, rho)

    implicit none

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

    integer :: i

    if (.not. tables_are_initialized) then

        ! generate es table (dt = 0.1 deg. c)

        allocate (table (es_table_length))
        allocate (table2 (es_table_length))
        allocate (table3 (es_table_length))
        allocate (tablew (es_table_length))
        allocate (des (es_table_length))
        allocate (des2 (es_table_length))
        allocate (des3 (es_table_length))
        allocate (desw (es_table_length))

        call qs_table (es_table_length)
        call qs_table2 (es_table_length)
        call qs_table3 (es_table_length)
        call qs_tablew (es_table_length)

        do i = 1, es_table_length - 1
            des (i) = max (0., table (i + 1) - table (i))
            des2 (i) = max (0., table2 (i + 1) - table2 (i))
            des3 (i) = max (0., table3 (i + 1) - table3 (i))
            desw (i) = max (0., tablew (i + 1) - tablew (i))
        enddo
        des (es_table_length) = des (es_table_length - 1)
        des2 (es_table_length) = des2 (es_table_length - 1)
        des3 (es_table_length) = des3 (es_table_length - 1)
        desw (es_table_length) = desw (es_table_length - 1)

        tables_are_initialized = .true.

    endif

end subroutine qsmith_init

! =======================================================================
! compute the saturated specific humidity for table ii
!>@brief The function 'wqs1' returns the saturation vapor pressure over pure
!! liquid water for a given temperature and air density.
! =======================================================================

real function wqs1 (ta, den)

    implicit none

    !> pure water phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es

    integer :: it, ap1

    ap1 = rdelt * max(0.0,dim (ta, es_table_tmin)) + 1.
    ap1 = min(es_table_length, ap1)
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

    !> pure water phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real, intent (out) :: dqdt

    real :: es

    integer :: it, ap1


    if (.not. tables_are_initialized) call qsmith_init

    ap1 = rdelt * max(0.0,dim (ta, es_table_tmin)) + 1.
    ap1 = min (es_table_length, ap1)
    it = ap1
    es = tablew (it) + (ap1 - it) * desw (it)
    wqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    ! finite diff, del_t = 0.1:
    dqdt = rdelt * (desw (it) + (ap1 - it) * (desw (it + 1) - desw (it))) / (rvgas * ta * den)

end function wqs2

! =======================================================================
!>@brief The function 'iqs1' computes the saturated specific humidity
!! for table iii
! =======================================================================

real function iqs1 (ta, den)

    implicit none

    !> water - ice phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real :: es

    integer :: it, ap1

    ap1 = rdelt * max(0.0,dim (ta, es_table_tmin))+ 1.
    ap1 = min (es_table_length, ap1)
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

    !> water - ice phase; universal dry / moist formular using air density
    !> input "den" can be either dry or moist air density

    real, intent (in) :: ta, den

    real, intent (out) :: dqdt

    real :: es

    integer :: it, ap1

    ap1 = rdelt * max(0.0,dim (ta, es_table_tmin)) + 1.
    ap1 = min(es_table_length, ap1)
    it = ap1
    es = table2 (it) + (ap1 - it) * des2 (it)
    iqs2 = es / (rvgas * ta * den)
    it = ap1 - 0.5
    dqdt = rdelt * (des2 (it) + (ap1 - it) * (des2 (it + 1) - des2 (it))) / (rvgas * ta * den)

end function iqs2

! =======================================================================
!>@brief saturation water vapor pressure table ii
! 1 - phase table
! =======================================================================

subroutine qs_tablew (n)

    implicit none

    integer, intent (in) :: n

    real :: tem, fac0, fac1, fac2

    integer :: i

    ! -----------------------------------------------------------------------
    ! compute es over water
    ! -----------------------------------------------------------------------

    do i = 1, n
        tem = es_table_tmin + delt * real (i - 1)
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

    real :: tem0, tem1, fac0, fac1, fac2

    integer :: i, i0, i1

    do i = 1, n
        tem0 = es_table_tmin + delt * real (i - 1)
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

    real :: esbasw, tbasw, esbasi, tem, aa, b, c, d, e
    real :: tem0, tem1

    integer :: i, i0, i1

    esbasw = 1013246.0
    tbasw = table_ice + 100.
    esbasi = 6107.1

    do i = 1, n
        tem = es_table_tmin + delt * real (i - 1)
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

    real :: es

    integer :: it, ap1

    ap1 = rdelt * max(0.0,dim (t, es_table_tmin)) + 1.
    ap1 = min (es_table_length, ap1)
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

    real :: tem, esh40
    real :: wice, wh2o, fac0, fac1, fac2
    real :: esupc (400)

    integer :: i
    real :: tc

    ! -----------------------------------------------------------------------
    ! compute es over ice between - 160 deg c and 0 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, 1600
        tem = es_table_tmin + delt * real (i - 1)
        fac0 = (tem - t_ice) / (tem * t_ice)
        fac1 = fac0 * li2
        fac2 = (d2ice * log (tem / t_ice) + fac1) / rvgas
        table (i) = e00 * exp (fac2)
    enddo

    ! -----------------------------------------------------------------------
    ! compute es over water between - 40 deg c and 102 deg c.
    ! -----------------------------------------------------------------------

    do i = 1, es_table_length-1200
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
! GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C
        wice = ice_fraction(tem,0.0,0.0)
        wh2o = 1.0 - wice
        table (i + 1200) = wice * table (i + 1200) + wh2o * esupc (i)
    enddo

end subroutine qs_table

! =======================================================================
! compute the saturated specific humidity and the gradient of saturated specific humidity
! input t in deg k, p in pa; p = rho rdry tv, moist pressure
!>@brief The function 'qsmith' computes the saturated specific humidity
!! with a blend of water and ice depending on the temperature in 3D.
!@details It als oincludes the option for computing des/dT.
! =======================================================================

subroutine qsmith (im, km, ks, t, p, q, qs, dqdt)

    implicit none

    integer, intent (in) :: im, km, ks

    real, intent (in), dimension (im, km) :: t, p, q

    real, intent (out), dimension (im, km) :: qs

    real, intent (out), dimension (im, km), optional :: dqdt

    real :: eps10

    real, dimension (im, km) :: es

    integer :: i, k, it, ap1

    eps10 = 10. * eps

    if (.not. tables_are_initialized) then
        call qsmith_init
    endif

    do k = ks, km
        do i = 1, im
            ap1 = rdelt * max(0.0,dim (t (i, k), es_table_tmin)) + 1.
            ap1 = min (es_table_length, ap1)
            it = ap1
            es (i, k) = table (it) + (ap1 - it) * des (it)
            qs (i, k) = eps * es (i, k) * (1. + zvir * q (i, k)) / p (i, k)
        enddo
    enddo

    if (present (dqdt)) then
        do k = ks, km
            do i = 1, im
                ap1 = rdelt * max(0.0,dim (t (i, k), es_table_tmin)) + 1. 
                ap1 = min (es_table_length, ap1) - 0.5
                it = ap1
                dqdt (i, k) = eps10 * (des (it) + (ap1 - it) * (des (it + 1) - des (it))) * (1. + zvir * q (i, k)) / p (i, k)
            enddo
        enddo
    endif

end subroutine qsmith

! =======================================================================
!>@brief The subroutine 'neg_adj' fixes negative water species.
!>@details This is designed for 6-class micro-physics schemes.
! =======================================================================

subroutine neg_adj (ktop, kbot, pt, dp, qv, ql, qr, qi, qs, qg)

    implicit none

    integer, intent (in) :: ktop, kbot

    real, intent (in), dimension (ktop:kbot) :: dp

    real, intent (inout), dimension (ktop:kbot) :: pt, qv, ql, qr, qi, qs, qg

    real, dimension (ktop:kbot) :: lcpk, icpk

    real :: dq, cvm

    integer :: k

    ! -----------------------------------------------------------------------
    ! define heat capacity and latent heat coefficient
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        cvm = c_air + qv (k) * c_vap + (qr (k) + ql (k)) * c_liq + (qi (k) + qs (k) + qg (k)) * c_ice
        lcpk (k) = (lv00 + d0_vap * pt (k)) / cvm
        icpk (k) = (li00 + dc_ice * pt (k)) / cvm
    enddo

    do k = ktop, kbot

        ! -----------------------------------------------------------------------
        ! ice phase:
        ! -----------------------------------------------------------------------

        ! if cloud ice < 0, borrow from snow
        if (qi (k) < 0.) then
            qs (k) = qs (k) + qi (k)
            qi (k) = 0.
        endif
        ! if snow < 0, borrow from graupel
        if (qs (k) < 0.) then
            qg (k) = qg (k) + qs (k)
            qs (k) = 0.
        endif
        ! if graupel < 0, borrow from rain
        if (qg (k) < 0.) then
            qr (k) = qr (k) + qg (k)
            pt (k) = pt (k) - qg (k) * icpk (k) ! heating
            qg (k) = 0.
        endif

        ! -----------------------------------------------------------------------
        ! liquid phase:
        ! -----------------------------------------------------------------------

        ! if rain < 0, borrow from cloud water
        if (qr (k) < 0.) then
            ql (k) = ql (k) + qr (k)
            qr (k) = 0.
        endif
        ! if cloud water < 0, borrow from water vapor
        if (ql (k) < 0.) then
            qv (k) = qv (k) + ql (k)
            pt (k) = pt (k) - ql (k) * lcpk (k) ! heating
            ql (k) = 0.
        endif

    enddo

    ! -----------------------------------------------------------------------
    ! fix water vapor; borrow from below
    ! -----------------------------------------------------------------------

    do k = ktop, kbot - 1
        if (qv (k) < 0.) then
            qv (k + 1) = qv (k + 1) + qv (k) * dp (k) / dp (k + 1)
            qv (k) = 0.
        endif
    enddo

    ! -----------------------------------------------------------------------
    ! bottom layer; borrow from above
    ! -----------------------------------------------------------------------

    if (qv (kbot) < 0. .and. qv (kbot - 1) > 0.) then
        dq = min (- qv (kbot) * dp (kbot), qv (kbot - 1) * dp (kbot - 1))
        qv (kbot - 1) = qv (kbot - 1) - dq / dp (kbot - 1)
        qv (kbot) = qv (kbot) + dq / dp (kbot)
    endif

end subroutine neg_adj

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

real function new_liq_condensate(tk, qlk, qik, cnv_fraction, srf_type)
        
     real, intent(in) :: tk, qlk, qik, cnv_fraction, srf_type
     real :: ptc, ifrac

     ifrac = ice_fraction(tk,cnv_fraction, srf_type)
     new_liq_condensate = min(max(0.0,(1.0-ifrac)*(qlk+qik) - qlk),qik)
    
end function new_liq_condensate

real function new_ice_condensate(tk, qlk, qik, cnv_fraction, srf_type)

     real, intent(in) :: tk, qlk, qik, cnv_fraction, srf_type
     real :: ptc, ifrac

     ifrac = ice_fraction(tk,cnv_fraction, srf_type)
     new_ice_condensate = min(max(0.0,ifrac*(qlk+qik) - qik),qlk)

end function new_ice_condensate

end module gfdl2_cloud_microphys_mod
