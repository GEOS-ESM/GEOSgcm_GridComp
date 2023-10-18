module moist_subroutines_cloud_microphys

    use Process_Library_standalone

    implicit none

    private

    public :: gfdl_cloud_microphys_driver, update_microphys_constants

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
    
    real, allocatable :: table (:), table2 (:), table3 (:), tablew (:)
    real, allocatable :: des (:), des2 (:), des3 (:), desw (:)
    
    logical :: tables_are_initialized = .true.
    
    ! logical :: root_proc
    ! integer :: id_rh, id_vtr, id_vts, id_vtg, id_vti, id_rain, id_snow, id_graupel, &
    ! id_ice, id_prec, id_cond, id_var, id_droplets
    ! integer :: gfdl_mp_clock ! clock for timing of driver routine
    
    real :: dt_fr = 8. !< epsilon on homogeneous freezing of cloud water at t_wfr + dt_fr
    ! minimum temperature water can exist (moore & molinero nov. 2011, nature)
    ! dt_fr can be considered as the error bar
    
    real :: p_min = 100. !< minimum pressure (pascal) for mp to operate
    
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
    real :: t_wfr = 273.16 - 40.0 ! supercooled water can exist down to - 40 c, which is the "absolute"
 
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

!$acc declare create(missing_value,  &
!$acc               cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw, &
!$acc               acco, cssub, cgsub, crevp, cgfr, csmlt, cgmlt, es0, ces0, pie, rgrav, &
!$acc               c_air, c_vap, lati, latv, lats, lat2, lcp, icp, tcp, d0_vap, &
!$acc               lv00, icloud_f, irain_f, de_ice, sedi_transport, do_sedi_w, do_sedi_heat, &
!$acc               prog_ccn, do_bigg, do_evap, do_subl, do_qa, preciprad, fix_negative, do_setup, p_nonhydro, table, &
!$acc               table2, table3, tablew, des, des2, des3, desw, dt_fr, p_min, cld_min, &
!$acc               tice, log_10, tice0, t_wfr, t_min, t_sub, mp_time, rh_inc, rh_inr, rh_ins, &
!$acc               tau_r2g, tau_smlt, tau_g2r, tau_imlt, tau_i2s, tau_l2r, tau_v2l, tau_l2v, &
!$acc               tau_i2v, tau_s2v, tau_v2s, tau_g2v, tau_v2g, tau_revp, tau_frz, dw_land, dw_ocean, ccn_o, &
!$acc               ccn_l, rthreshu, rthreshs, sat_adj0, qc_crt, qi_lim, ql_mlt, qs_mlt, ql_gen, qi_gen, &
!$acc               ql0_max, qi0_max, qi0_crt, qr0_crt, qs0_crt, c_paut, c_psaci, c_piacr, &
!$acc               c_cracw, c_pgacs, c_pgaci, alin, clin, const_vi, const_vs, const_vg, const_vr, &
!$acc               vi_fac, vs_fac, vg_fac, vr_fac, vi_max, vs_max, vg_max, vr_max, fast_sat_adj, &
!$acc               z_slope_liq, z_slope_ice, use_ccn, use_ppm, mono_prof, mp_print, tables_are_initialized)

    contains

    ! -----------------------------------------------------------------------
    !> sedimentation of heat
    ! -----------------------------------------------------------------------

    subroutine sedi_heat (ktop, kbot, dm, m1, dz, tz, qv, ql, qr, qi, qs, qg, cw)
    !$acc routine seq
        implicit none
    
        ! input q fields are dry mixing ratios, and dm is dry air mass
        
        integer, intent (in) :: ktop, kbot
        
        real, intent (in), dimension (ktop:kbot) :: dm, m1, dz, qv, ql, qr, qi, qs, qg
        
        real, intent (inout), dimension (ktop:kbot) :: tz
        
        real, intent (in) :: cw ! heat capacity
        
        real :: dgz, cvn
        
        real :: tmp
        
        integer :: k
        
        ! -----------------------------------------------------------------------
        ! sjl, july 2014
        ! assumption: the ke in the falling condensates is negligible compared to the potential energy
        ! that was unaccounted for. local thermal equilibrium is assumed, and the loss in pe is transformed
        ! into internal energy (to heat the whole grid box)
        ! backward time - implicit upwind transport scheme:
        ! dm here is dry air mass
        ! -----------------------------------------------------------------------
        
        k = ktop
        cvn = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * &
            c_liq + (qi (k) + qs (k) + qg (k)) * c_ice)
        tz (k) = ((cvn + m1 (k) * cw) * tz (k) + m1 (k) * (- 0.5 * grav * dz (k))) / (cvn + m1 (k) * cw)
        
        ! -----------------------------------------------------------------------
        ! implicit algorithm: can't be vectorized
        ! needs an inner i - loop for vectorization
        ! -----------------------------------------------------------------------
        
        do k = ktop + 1, kbot
            dgz = - 0.5 * grav * dz (k) ! > 0
            cvn = dm (k) * (cv_air + qv (k) * cv_vap + (qr (k) + ql (k)) * &
                c_liq + (qi (k) + qs (k) + qg (k)) * c_ice)
            tz (k) = ((cvn + cw * (m1 (k) - m1 (k - 1))) * tz (k) + m1 (k - 1) * &
                cw * tz (k - 1) + dgz * (m1 (k - 1) + m1 (k))) / (cvn + cw * m1 (k))
        enddo
        
    end subroutine sedi_heat

    ! -----------------------------------------------------------------------
    !> warm rain cloud microphysics
    ! -----------------------------------------------------------------------

    subroutine warm_rain (dt, ktop, kbot, dp, dz, tz, qv, ql, qr, qi, qs, qg, qa, &
        eis, &
        den, denfac, ccn, c_praut, vtr, r1, evap1, m1_rain, w1, h_var)
    !$acc routine vector
        implicit none
    
        integer, intent (in) :: ktop, kbot
        
        real, intent (in) :: dt !< time step (s)

        real, intent (in), dimension (ktop:kbot) :: h_var
        
        real, intent (in), dimension (ktop:kbot) :: dp, dz, den
        real, intent (in), dimension (ktop:kbot) :: denfac, ccn, c_praut

        real, intent (in) :: eis !< estimated inversion strength
        
        real, intent (inout), dimension (ktop:kbot) :: tz, vtr
        real, intent (inout), dimension (ktop:kbot) :: qv, ql, qr, qi, qs, qg, qa
        real, intent (inout), dimension (ktop:kbot) :: evap1, m1_rain, w1

        real, intent (out) :: r1
        
        real, parameter :: so3 = 7. / 3.
        
        real, dimension (ktop:kbot) :: dl, dm, revap, isubl, qadum
        real, dimension (ktop:kbot + 1) :: ze, zt
        
        real :: sink, dq, qc0, qc
        real :: fac_rc, qden
        real :: zs
        real :: dt5
        
        integer :: k
        
        ! fall velocity constants:
        
        real, parameter :: vconr = 2503.23638966667
        real, parameter :: normr = 25132741228.7183
        real, parameter :: thr = 1.e-8
        
        logical :: no_fall

        zs = 0.
    
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
        if (.not. do_qa) then
            qadum = max(qa,qcmin)
        else
            qadum = 1.0
        endif
        ql = ql/qadum
        qi = qi/qadum

        fac_rc = min(1.0,eis/10.0)**2 ! Estimated inversion strength determine stable regime
        fac_rc = rc * (rthreshs*fac_rc + rthreshu*(1.0-fac_rc)) ** 3
    
        if (irain_f /= 0) then
            
            ! -----------------------------------------------------------------------
            ! no subgrid varaibility
            ! -----------------------------------------------------------------------
!$acc loop vector private(qc0, qc, dq, sink)
            do k = ktop, kbot
                qc0 = fac_rc * ccn (k)
                if (tz (k) > t_wfr) then
                    qc = qc0 / den (k)
                    dq = ql (k) - qc
                    if (dq > 0.) then
                        sink = min (dq, dt * c_praut (k) * den (k) * exp (so3 * log (ql (k))))
                        sink = min(ql0_max/qadum(k), ql(k), max(0.,sink))
                        ql (k) = ql (k) - sink
                        qr (k) = qr (k) + sink*qadum(k)
                    endif
                endif
            enddo
        
        else
            
            ! -----------------------------------------------------------------------
            ! with subgrid variability
            ! -----------------------------------------------------------------------

            call linear_prof (kbot - ktop + 1, ql (ktop), dl (ktop), z_slope_liq, h_var)
!$acc loop vector private(qc0, qc, dq, sink)
            do k = ktop, kbot
                qc0 = fac_rc * ccn (k)
                if (tz (k) > t_wfr + dt_fr) then
                    dl (k) = min (max (qcmin, dl (k)), 0.5 * ql (k))
                    ! --------------------------------------------------------------------
                    ! as in klein's gfdl am2 stratiform scheme (with subgrid variations)
                    ! --------------------------------------------------------------------
                    qc = qc0 / den (k)
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
            vtr (:) = vf_min
        elseif (const_vr) then
            vtr (:) = vr_fac ! ifs_2016: 4.0
        else
!$acc loop vector private(qden)
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
!$acc loop seq
        do k = kbot, ktop, - 1
            ze (k) = ze (k + 1) - dz (k) ! dz < 0
        enddo
        
        ! -----------------------------------------------------------------------
        ! evaporation and accretion of rain for the first 1 / 2 time step
        ! -----------------------------------------------------------------------
        
        call revap_racc (ktop, kbot, dt5, tz, qv, ql, qr, qi, qs, qg, qa, revap, den, denfac, h_var)
        evap1 = revap
 
        if (do_sedi_w) then
!$acc loop vector
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
!$acc loop vector
            do k = ktop + 1, kbot
                zt (k) = ze (k) - dt * (vtr (k - 1) + vtr (k))/2.0
            enddo
            zt (kbot + 1) = zs - dt * vtr (kbot)
!$acc loop seq
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
!$acc loop vector
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
        !$acc routine seq
        implicit none
    
        integer, intent (in) :: ktop, kbot
        
        real, intent (in) :: dt ! time step (s)
        
        real, intent (in), dimension (ktop:kbot) :: h_var
        real, intent (in), dimension (ktop:kbot) :: den, denfac
        
        real, intent (inout), dimension (ktop:kbot) :: tz, qv, qr, ql, qi, qs, qg, qa

        real, intent (inout), dimension (ktop:kbot) :: revap

        real :: lhl, cvm, q_liq, q_sol, lcpk  
        real :: dqv, qsat, dqsdt, evap, t2, qden, q_plus, q_minus, sink
        real :: qpz, dq, dqh, tin
        real :: fac_revp 
        real :: TOT_PREC_LS, AREA_LS_PRC, AREA_LS_PRC_K 
        integer :: k
    
        revap(:) = 0.
    
        TOT_PREC_LS = 0.
        AREA_LS_PRC = 0.
        do k = ktop, kbot
            
            TOT_PREC_LS = TOT_PREC_LS  + (          ( qr (k) + qs (k) + qg (k) ) * den (k) )
            AREA_LS_PRC = AREA_LS_PRC  + ( qa (k) * ( qr (k) + qs (k) + qg (k) ) * den (k) )

            if (tz (k) > t_wfr .and. qr (k) > qpmin) then

            !! area and timescale efficiency on revap
            !                       AREA_LS_PRC_K = 0.0
            !if (TOT_PREC_LS > 0.0) AREA_LS_PRC_K = MAX( AREA_LS_PRC/TOT_PREC_LS, 1.E-6 )
            !fac_revp = 1. - exp (- AREA_LS_PRC_K * dt / tau_revp)
                fac_revp = 1. - exp (- dt / tau_revp)
                
                ! -----------------------------------------------------------------------
                ! define heat capacity and latent heat coefficient
                ! -----------------------------------------------------------------------
                
                lhl = lv00 + d0_vap * tz (k)
                q_liq = ql (k) + qr (k)
                q_sol = qi (k) + qs (k) + qg (k)
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                lcpk = lhl / cvm
                
                tin = tz (k) - lcpk * ql (k) ! presence of clouds suppresses the rain evap
                qpz = qv (k) + ql (k)
                qsat = wqs2 (tin, den (k), dqsdt)
                dqh = max (ql (k), h_var(k) * max (qpz, qcmin))
                dqh = min (dqh, 0.2 * qpz) ! new limiter
                dqv = qsat - qv (k) ! use this to prevent super - sat the gird box
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
                    evap = min (qr (k), dt * fac_revp * evap, dqv / (1. + lcpk * dqsdt))
                    qr (k) = qr (k) - evap
                    qv (k) = qv (k) + evap
                    q_liq = q_liq - evap
                    cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                    tz (k) = tz (k) - evap * lhl / cvm
                    revap(k) = evap / dt
                endif
                
                ! ! -----------------------------------------------------------------------
                ! ! accretion: pracc
                ! ! -----------------------------------------------------------------------
                
                if (qr (k) > qpmin .and. ql (k) > qcmin .and. qsat < q_minus) then
                    sink = dt * denfac (k) * cracw * exp (0.95 * log (qr (k) * den (k)))
                    sink = sink / (1. + sink) * ql (k)
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
    !$acc routine vector
        implicit none
    
        integer, intent (in) :: km
        
        real, intent (in) :: q (km), h_var(km)
        
        real, intent (out) :: dm (km)
        
        logical, intent (in) :: z_var
        
        real :: dq, dq_p1
        
        integer :: k

        if (z_var) then
            dm (1) = 0.
            
            ! -----------------------------------------------------------------------
            ! use twice the strength of the positive definiteness limiter (lin et al 1994)
            ! -----------------------------------------------------------------------
!$acc loop vector private(dq, dq_p1)
            do k = 2, km - 1
                dq = 0.5 * (q (k) - q (k - 1))
                dq_p1 = 0.5 * (q (k + 1) - q (k))
                dm (k) = 0.5 * min (abs (dq + dq_p1), 0.5 * q (k))
                if (dq * dq_p1 <= 0.) then
                    if (dq > 0.) then ! local max
                        dm (k) = min (dm (k), dq, - dq_p1)
                    else
                        dm (k) = 0.
                    endif
                endif
            enddo
            dm (km) = 0.
            
            ! -----------------------------------------------------------------------
            ! impose a presumed background horizontal variability that is proportional to the value itself
            ! -----------------------------------------------------------------------
!$acc loop vector
            do k = 1, km
                dm (k) = max (dm (k), qvmin, h_var(k) * q (k))
            enddo
        else
!$acc loop vector
            do k = 1, km
                                dm (k) = max (qvmin, h_var(k) * q (k))
            enddo
        endif
        
    end subroutine linear_prof

    subroutine icloud (ktop, kbot, tzk, p1, qvk, qlk, qrk, qik, qsk, qgk, dp1, &
            den, denfac, vts, vtg, vtr, qak, dts, subl1, h_var, ccn, cnv_fraction, srf_type)
    !$acc routine vector
        implicit none

        integer, intent (in) :: ktop, kbot
        
        real, intent (in), dimension (ktop:kbot) :: p1, dp1, den, denfac, vts, vtg, vtr
        
        real, intent (inout), dimension (ktop:kbot) :: tzk, qvk, qlk, qrk, qik, qsk, qgk, qak

        real, intent (out), dimension (ktop:kbot) :: subl1

        real, intent (in) :: dts, cnv_fraction, srf_type
    
        real, intent (in), dimension (ktop:kbot) :: h_var, ccn
                
        real, dimension (ktop:kbot) :: di
        real, dimension (ktop:kbot) :: cvm, q_liq, q_sol
        
        real :: icpk, tcpk, lhl, lhi
        real :: rdts, fac_g2v, fac_i2s, fac_imlt, fac_frz
        real :: tz, qv, ql, qr, qi, qs, qg, melt, ifrac, newqi, newql
        real :: pracs, psacw, pgacw, psacr, pgacr, pgaci, praci, psaci
        real :: pgmlt, psmlt, pgfr, pgaut, psaut, pgsub
        real :: tc, tsq, dqs0, qden, qim, qsm
        real :: factor, sink, qi_crt
        real :: tmp, qsw, qsi, dqsdt, dq
        real :: dtmp, qc, q_plus, q_minus
        
        integer :: k, it
    
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
!$acc loop vector
        do k = ktop, kbot
            q_liq (k) = qlk (k) + qrk (k)
            q_sol (k) = qik (k) + qsk (k) + qgk (k)
            cvm (k) = c_air + qvk (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
        enddo
    
        ! -----------------------------------------------------------------------
        ! sources of cloud ice: pihom, cold rain, and the sat_adj
        ! (initiation plus deposition)
        ! sources of snow: cold rain, auto conversion + accretion (from cloud ice)
        ! sat_adj (deposition; requires pre - existing snow) ; initial snow comes from auto conversion
        ! -----------------------------------------------------------------------
!*****
! Note : If 'sink' gets added as a private variable, the code will not verify
!*****
!$acc loop vector private(lhi, icpk, melt, tmp)
        do k = ktop, kbot
            lhi = li00 + dc_ice * tzk (k)
            if (tzk (k) > tice .and. qik (k) > qcmin) then
                icpk = lhi / cvm (k)
                ! -----------------------------------------------------------------------
                ! pimlt: instant melting of cloud ice
                ! -----------------------------------------------------------------------
                
                melt = min (qik (k), fac_imlt * (tzk (k) - tice) / icpk)
                tmp = min (melt, dim (ql_mlt, qlk (k))) ! max ql amount

                ! new total condensate / old condensate 
                qak(k) = max(0.0,min(1.,qak(k) * max(qik(k)+qlk(k)-melt+tmp,0.0  ) / &
                                                max(qik(k)+qlk(k)         ,qcmin) ) )

                qlk (k) = qlk (k) + tmp
                qrk (k) = qrk (k) + melt - tmp
                qik (k) = qik (k) - melt
                q_liq (k) = q_liq (k) + melt
                q_sol (k) = q_sol (k) - melt
                cvm (k) = c_air + qvk (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tzk (k) = tzk (k) - melt * lhi / cvm (k)

            elseif (tzk (k) <= tice .and. qlk (k) > qcmin) then

                ! -----------------------------------------------------------------------
                ! pihom: homogeneous freezing of cloud water into cloud ice
                ! this is the 1st occurance of liquid water freezing in the split mp process
                ! -----------------------------------------------------------------------

                sink = fac_frz * new_ice_condensate(tzk (k), qlk (k), qik (k), cnv_fraction, srf_type)
                qi_crt = qi0_crt / den (k)
                tmp = min (sink, dim (qi_crt, qik (k)))

                ! new total condensate / old condensate 
                qak(k) = max(0.0,min(1.,qak(k) * max(qik(k)+qlk(k)-sink+tmp,0.0  ) / &
                                                max(qik(k)+qlk(k)         ,qcmin) ) )

                qlk (k) = qlk (k) - sink
                qsk (k) = qsk (k) + sink - tmp
                qik (k) = qik (k) + tmp
                q_liq (k) = q_liq (k) - sink
                q_sol (k) = q_sol (k) + sink
                cvm (k) = c_air + qvk (k) * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                tzk (k) = tzk (k) + sink * lhi / cvm (k)
                
            endif
        enddo
  
        ! -----------------------------------------------------------------------
        ! vertical subgrid variability
        ! -----------------------------------------------------------------------
        
        call linear_prof (kbot - ktop + 1, qik (ktop), di (ktop), z_slope_ice, h_var)
        
        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------
!$acc loop seq
        do k = ktop, kbot
            lhi = li00 + dc_ice * tzk (k)
            lhl = lv00 + d0_vap * tzk (k)
            icpk = lhi / cvm (k)
            tcpk = lhl / cvm (k) + icpk
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
            tc = tz - tice
            
            if (tc .ge. 0.) then
            
                ! -----------------------------------------------------------------------
                ! melting of snow
                ! -----------------------------------------------------------------------
                
                dqs0 = ces0 / p1 (k) - qv
                
                if (qs > qpmin) then
                
                    ! -----------------------------------------------------------------------
                    ! psacw: accretion of cloud water by snow
                    ! only rate is used (for snow melt) since tc > 0.
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
                    
                    psmlt = max (0., smlt (tc, dqs0, qs * den (k), psacw, psacr, csmlt, &
                        den (k), denfac (k)))
                    sink = min (qs, dts * (psmlt + pracs), tc / icpk)
                    qs = qs - sink
                    ! sjl, 20170321:
                    tmp = min (sink, dim (qs_mlt, ql)) ! max ql due to snow melt

                    ! new total condensate / old condensate 
                    qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql+tmp,0.0  ) / &
                                                        max(qi+ql    ,qcmin) ) )

                    ql = ql + tmp
                    qr = qr + sink - tmp
                    ! qr = qr + sink
                    ! sjl, 20170321:
                    q_liq (k) = q_liq (k) + sink
                    q_sol (k) = q_sol (k) - sink
                    cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    tz = tz - sink * lhi / cvm (k)
                    tc = tz - tice
                    
                endif
            
                ! -----------------------------------------------------------------------
                ! update capacity heat and latend heat coefficient
                ! -----------------------------------------------------------------------
                
                lhi = li00 + dc_ice * tz
                icpk = lhi / cvm (k)
            
                ! -----------------------------------------------------------------------
                ! melting of graupel
                ! -----------------------------------------------------------------------
                
                if (qg > qpmin .and. tc > 0.) then
                    
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
                    
                    pgmlt = dts * gmlt (tc, dqs0, qden, pgacw, pgacr, cgmlt, den (k))
                    pgmlt = min (max (0., pgmlt), qg, tc / icpk)
                    qg = qg - pgmlt
                    qr = qr + pgmlt
                    q_liq (k) = q_liq (k) + pgmlt
                    q_sol (k) = q_sol (k) - pgmlt
                    cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    tz = tz - pgmlt * lhi / cvm (k)
                    
                endif
            
            else
            
                ! -----------------------------------------------------------------------
                ! cloud ice proc:
                ! -----------------------------------------------------------------------
                
                ! -----------------------------------------------------------------------
                ! psaci: accretion of cloud ice by snow
                ! -----------------------------------------------------------------------
                
                if (qi > 3.e-7) then ! cloud ice sink terms
                
                    if (qs > qpmin) then
                        ! -----------------------------------------------------------------------
                        ! sjl added (following lin eq. 23) the temperature dependency
                        ! to reduce accretion, use esi = exp (0.05 * tc) as in hong et al 2004
                        ! -----------------------------------------------------------------------
                        factor = dts * denfac (k) * csaci * exp (0.05 * tc + 0.8125 * log (qs * den (k)))
                        psaci = factor / (1. + factor) * qi
                    else
                        psaci = 0.
                    endif
                
                    ! -----------------------------------------------------------------------
                    ! psaut: autoconversion: cloud ice -- > snow
                    ! -----------------------------------------------------------------------
                    
                    ! -----------------------------------------------------------------------
                    ! similar to lfo 1983: eq. 21 solved implicitly
                    ! threshold from wsm6 scheme, hong et al 2004, eq (13) : qi0_crt ~0.8e-4
                    ! -----------------------------------------------------------------------
                
    !!!!!!!!        qim = qi0_crt / den (k)
    ! GEOS ! WMP impose CALIPSO ice polynomial from 0 C to -40 C on qi0_crt  
                    qim = ice_fraction(tz,cnv_fraction,srf_type) * qi0_crt / den (k)
 
                    ! -----------------------------------------------------------------------
                    ! assuming linear subgrid vertical distribution of cloud ice
                    ! the mismatch computation following lin et al. 1994, mwr
                    ! -----------------------------------------------------------------------
                    
                    if (const_vi) then
                        tmp = fac_i2s
                    else
                        tmp = fac_i2s * exp (0.025 * tc)
                    endif
                
                    di (k) = max (di (k), qcmin)
                    q_plus = qi + di (k)
                    if (q_plus > (qim + qcmin)) then
                        if (qim > (qi - di (k))) then
                            dq = (0.25 * (q_plus - qim) ** 2) / di (k)
                        else
                            dq = qi - qim
                        endif
                        psaut = tmp * dq
                    else
                        psaut = 0.
                    endif
                    sink = min (qi, psaci + psaut)

                    ! new total condensate / old condensate 
                    qak(k) = max(0.0,min(1.,qak(k) * max(qi+ql-sink+tmp,0.0  ) / &
                                                    max(qi+ql         ,qcmin) ) )

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
                
                tc = tz - tice
                
                if (qr > qpmin .and. tc < 0.) then
                
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
                    factor = min (sink, qr, - tc / icpk) / max (sink, qpmin)
                    
                    psacr = factor * psacr
                    pgfr = factor * pgfr
                    
                    sink = psacr + pgfr
                    qr = qr - sink
                    qs = qs + psacr
                    qg = qg + pgfr
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink
                    cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    tz = tz + sink * lhi / cvm (k)
                    
                endif
            
                ! -----------------------------------------------------------------------
                ! update capacity heat and latend heat coefficient
                ! -----------------------------------------------------------------------
                
                lhi = li00 + dc_ice * tz
                icpk = lhi / cvm (k)
                
                ! -----------------------------------------------------------------------
                ! graupel production terms:
                ! -----------------------------------------------------------------------
                
                if (qs > qpmin) then
                
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
                        sink = sink + factor / (1. + factor) * (qs - qsm)
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
                    factor = min (sink, dim (tice, tz) / icpk) / max (sink, qpmin)
                    pgacr = factor * pgacr
                    pgacw = factor * pgacw
                    
                    sink = pgacr + pgacw
                    qg = qg + sink
                    qr = qr - pgacr
                    ql = ql - pgacw
                    q_liq (k) = q_liq (k) - sink
                    q_sol (k) = q_sol (k) + sink
                    cvm (k) = c_air + qv * c_vap + q_liq (k) * c_liq + q_sol (k) * c_ice
                    tz = tz + sink * lhi / cvm (k)
                    
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
            qlk, qrk, qik, qsk, qgk, qak, subl1, h_var, ccn, cnv_fraction, srf_type)

    end subroutine icloud

    ! =======================================================================
    !>temperature sensitive high vertical resolution processes
    ! =======================================================================

    subroutine subgrid_z_proc (ktop, kbot, p1, den, denfac, dts, tz, qv, &
        ql, qr, qi, qs, qg, qa, subl1, h_var, ccn, cnv_fraction, srf_type)
    !$acc routine seq
        implicit none
    
        integer, intent (in) :: ktop, kbot
        
        real, intent (in), dimension (ktop:kbot) :: p1, den, denfac
        
        real, intent (in) :: dts, cnv_fraction, srf_type

        real, intent (in), dimension (ktop:kbot) :: h_var, ccn

        real, intent (inout), dimension (ktop:kbot) :: tz, qv, ql, qr, qi, qs, qg, qa

        real, intent (out), dimension (ktop:kbot) :: subl1
        
        real :: lcpk, icpk, tcpk, tcp3, lhl, lhi
        real :: cvm, q_liq, q_sol, q_cond
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
        
        ! do k = ktop, kbot
        !     lhl = lv00 + d0_vap * tz (k)
        !     lhi = li00 + dc_ice * tz (k)
        !     q_liq = ql (k) + qr (k)
        !     q_sol = qi (k) + qs (k) + qg (k)
        !     cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
        !     lcpk = lhl / cvm
        !     icpk = lhi / cvm
        !     tcpk = lcpk + icpk
        !     tcp3 = lcpk + icpk * min (1., dim (tice, tz (k)) / (tice - t_wfr))
        ! enddo
        
        do k = ktop, kbot
            lhl = lv00 + d0_vap * tz (k)
            lhi = li00 + dc_ice * tz (k)
            q_liq = ql (k) + qr (k)
            q_sol = qi (k) + qs (k) + qg (k)
            cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice

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
                q_sol = q_sol + sink
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) + sink * (lhl + lhi) / cvm
                if (do_qa) qa (k) = 1. ! air fully saturated; 100 % cloud cover
                cycle
            endif
            
            ! -----------------------------------------------------------------------
            ! update heat capacity and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhl = lv00 + d0_vap * tz (k)
            lhi = li00 + dc_ice * tz (k)
            lcpk = lhl / cvm
            icpk = lhi / cvm
            tcpk = lcpk + icpk
            tcp3 = lcpk + icpk * min (1., dim (tice, tz (k)) / (tice - t_wfr))
        
            ! -----------------------------------------------------------------------
            ! instant evaporation / sublimation of all clouds if rh < rh_adj -- > cloud free
            ! -----------------------------------------------------------------------
            qpz = qv (k) + ql (k) + qi (k)
            tin = tz (k) - (lhl * (ql (k) + qi (k)) + lhi * qi (k)) / (c_air + &
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
                    evap = min (ql (k), factor * ql(k) / (1. + tcp3 * dwsdt))
                else
                    evap = 0.0
                endif
                qv (k) = qv (k) + evap
                ql (k) = ql (k) - evap
                q_liq = q_liq - evap
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) - evap * lhl / cvm
            endif

            ! -----------------------------------------------------------------------
            ! update heat capacity and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhi = li00 + dc_ice * tz (k)
            icpk = lhi / cvm
        
            ! -----------------------------------------------------------------------
            ! enforce complete freezing below - t_wfr
            ! -----------------------------------------------------------------------
            
            dtmp = t_wfr - tz (k)
            if (dtmp > 0. .and. ql (k) > qcmin) then
                sink = min (ql (k), fac_frz * dtmp / icpk)
                ql (k) = ql (k) - sink
                qi (k) = qi (k) + sink
                q_liq = q_liq - sink
                q_sol = q_sol + sink
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) + sink * lhi / cvm
            endif
            
            ! -----------------------------------------------------------------------
            ! update heat capacity and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhi = li00 + dc_ice * tz (k)
            icpk = lhi / cvm
            
            ! -----------------------------------------------------------------------
            ! bigg mechanism heterogeneous freezing on existing cloud nuclei
            ! -----------------------------------------------------------------------
            tc = tice - tz (k)
            if (do_bigg .and. ql (k) > qcmin .and. tc > 0.) then
                sink = fac_frz * (100.0/rhor/ccn(k)) * dts * (exp (0.66 * tc) - 1.) * den (k) * ql (k) * ql (k)
                sink = min (ql (k), tc / icpk, sink)
                ql (k) = ql (k) - sink
                qi (k) = qi (k) + sink
                q_liq = q_liq - sink
                q_sol = q_sol + sink
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) + sink * lhi / cvm
            endif ! significant ql existed
        
            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhl = lv00 + d0_vap * tz (k)
            lhi = li00 + dc_ice * tz (k)
            lcpk = lhl / cvm
            icpk = lhi / cvm
            tcpk = lcpk + icpk
        
            ! -----------------------------------------------------------------------
            ! sublimation / deposition of LS ice
            ! -----------------------------------------------------------------------
            
            if (tz (k) < tice) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                dq = (qv (k) - qsi)
                sink = min(qi(k), dq / (1. + tcpk * dqsdt))
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
                        sink = min (sink, max (qi_crt - qi (k), pidep), tmp / tcpk)
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
                q_sol = q_sol + sink
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) + sink * (lhl + lhi) / cvm
            endif

            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhl = lv00 + d0_vap * tz (k)
            lhi = li00 + dc_ice * tz (k)
            lcpk = lhl / cvm
            icpk = lhi / cvm
            tcpk = lcpk + icpk
        
            ! -----------------------------------------------------------------------
            ! sublimation / deposition of snow
            ! this process happens for all temp rage
            ! -----------------------------------------------------------------------
            
            if (qs (k) > qpmin) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                qden = qs (k) * den (k)
                tmp = exp (0.65625 * log (qden))
                tsq = tz (k) * tz (k)
                dq = (qsi - qv (k)) / (1. + tcpk * dqsdt)
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
                        pssub = max (fac_v2s * pssub, dq, (tz (k) - tice) / tcpk)
                    endif
                endif
                qs (k) = qs (k) - pssub
                qv (k) = qv (k) + pssub
                q_sol = q_sol - pssub
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) - pssub * (lhl + lhi) / cvm
            endif
            
            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhl = lv00 + d0_vap * tz (k)
            lhi = li00 + dc_ice * tz (k)
            lcpk = lhl / cvm
            icpk = lhi / cvm
            tcpk = lcpk + icpk
        
            ! -----------------------------------------------------------------------
            ! simplified 2 - way grapuel sublimation - deposition mechanism
            ! -----------------------------------------------------------------------
            
            if (qg (k) > qpmin) then
                qsi = iqs2 (tz (k), den (k), dqsdt)
                dq = (qv (k) - qsi) / (1. + tcpk * dqsdt)
                pgsub = (qv (k) / qsi - 1.) * qg (k)
                if (pgsub > 0.) then ! deposition
                    if (tz (k) > tice) then
                        pgsub = 0. ! no deposition
                    else
                        pgsub = min (fac_v2g * pgsub, 0.2 * dq, ql (k) + qr (k), &
                            (tice - tz (k)) / tcpk)
                    endif
                else ! submilation
                    pgsub = max (fac_g2v * pgsub, dq) * min (1., dim (tz (k), t_sub) * 0.1)
                    subl1(k) = subl1(k) + pgsub / dts
                endif
                qg (k) = qg (k) + pgsub
                qv (k) = qv (k) - pgsub
                q_sol = q_sol + pgsub
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) + pgsub * (lhl + lhi) / cvm
            endif
        
#ifdef USE_MIN_EVAP
            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhl = lv00 + d0_vap * tz (k)
            lcpk = lhl / cvm
            
            ! -----------------------------------------------------------------------
            ! * minimum evap of rain in dry environmental air
            ! -----------------------------------------------------------------------
            
            if (qr (k) > qpmin) then
                qsw = wqs2 (tz (k), den (k), dqsdt)
                sink = min (qr (k), dim (rh_rain * qsw, qv (k)) / (1. + lcpk * dqsdt))
                qv (k) = qv (k) + sink
                qr (k) = qr (k) - sink
                q_liq = q_liq - sink
                cvm = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) - sink * lhl / cvm
            endif
#endif
        
            ! -----------------------------------------------------------------------
            ! update capacity heat and latend heat coefficient
            ! -----------------------------------------------------------------------
            
            lhl = lv00 + d0_vap * tz (k)
            cvm = c_air + (qv (k) + q_liq + q_sol) * c_vap
            lcpk = lhl / cvm
        
            ! -----------------------------------------------------------------------
            ! compute cloud fraction
            ! -----------------------------------------------------------------------
            if (.not. do_qa) cycle

            ! -----------------------------------------------------------------------
            ! combine water species
            ! -----------------------------------------------------------------------
            if (preciprad) then
                q_sol = qi (k) + qs (k) + qg (k)
                q_liq = ql (k) + qr (k)
            else
                q_sol = qi (k)
                q_liq = ql (k)
            endif
            q_cond = q_liq + q_sol
            
            qpz = qv (k) + q_cond ! qpz is conserved
        
            ! -----------------------------------------------------------------------
            ! use the "liquid - frozen water temperature" (tin) to compute saturated specific humidity
            ! -----------------------------------------------------------------------
            
            tin = tz (k) - (lcpk * q_cond + icpk * q_sol) ! minimum temperature
            ! tin = tz (k) - ((lv00 + d0_vap * tz (k)) * q_cond + &
            ! (li00 + dc_ice * tz (k)) * q_sol) / (c_air + qpz * c_vap)
            
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
                if (q_cond > 3.e-6) then
                    rqi = q_sol / q_cond
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
                    elseif (qstar < q_plus .and. q_cond > qc_crt) then
                        qa (k) = max(qcmin, min(1., qa (k) + (q_plus - qstar) / (dq + dq) )) ! partial cloud cover
                    elseif (qstar .le. q_minus) then
                        qa (k) = 1.0 ! air fully saturated; 100 % cloud cover
                    endif
                endif
            endif
            
        enddo
        
    end subroutine subgrid_z_proc

    subroutine terminal_fall (dtm, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
            den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1)
    !$acc routine vector
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
        
        real, dimension (ktop:kbot) :: icpk, cvm
        real, dimension (ktop:kbot) :: m1, dm
        
        real :: q_liq, q_sol, lcpk, lhl, lhi

        real :: zs
        real :: fac_imlt
        
        integer :: k, k0, m
        
        logical :: no_fall, exit_flag
    
        zs = 0.

        fac_imlt = 1. - exp (- dtm / tau_imlt)
        
        ! -----------------------------------------------------------------------
        ! define heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------
!$acc loop vector private(lhi, q_liq, q_sol)
        do k = ktop, kbot
            m1_sol (k) = 0.
            ! lhl (k) = lv00 + d0_vap * tz (k)
            lhi = li00 + dc_ice * tz (k)
            q_liq = ql (k) + qr (k)
            q_sol = qi (k) + qs (k) + qg (k)
            cvm (k) = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
            ! lcpk = lhl (k) / cvm (k)
            icpk (k) = lhi / cvm (k)
        enddo
    
        ! -----------------------------------------------------------------------
        ! find significant melting level
        ! -----------------------------------------------------------------------
        
        k0 = kbot
        exit_flag = .true.
!$acc loop seq
        do k = ktop, kbot - 1
            if (tz (k) > tice .and. exit_flag) then
                k0 = k
                exit_flag = .false.
            endif
        enddo
    
        ! -----------------------------------------------------------------------
        ! melting of cloud_ice (before fall) :
        ! -----------------------------------------------------------------------
!$acc loop vector private(tc, q_liq, q_sol, lhi, sink, tmp)
        do k = k0, kbot
            tc = tz (k) - tice
            if (qi (k) > qcmin .and. tc > 0.) then
                q_liq = ql (k) + qr (k)
                q_sol = qi (k) + qs (k) + qg (k)
                lhi = li00 + dc_ice * tz (k)

                sink = min (qi (k), fac_imlt * tc / icpk (k))
                tmp = min (sink, dim (ql_mlt, ql (k)))
                ql (k) = ql (k) + tmp
                qr (k) = qr (k) + sink - tmp
                qi (k) = qi (k) - sink
                q_liq = q_liq + sink
                q_sol = q_sol - sink
                cvm (k) = c_air + qv (k) * c_vap + q_liq * c_liq + q_sol * c_ice
                tz (k) = tz (k) - sink * lhi / cvm (k)
                tc = tz (k) - tice
            endif
        enddo
    
        ! -----------------------------------------------------------------------
        ! turn off melting when cloud microphysics time step is small
        ! -----------------------------------------------------------------------
        
        if (dtm < 60.) k0 = kbot
        
        ! sjl, turn off melting of falling cloud ice, snow and graupel
        k0 = kbot
        ! sjl, turn off melting of falling cloud ice, snow and graupel
        
        ze (kbot + 1) = zs
!$acc loop seq
        do k = kbot, ktop, - 1
            ze (k) = ze (k + 1) - dz (k) ! dz < 0
        enddo
    
        zt (ktop) = ze (ktop)
        
        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------
!$acc loop vector private(lhi)
        do k = k0, kbot
            lhi = li00 + dc_ice * tz (k)
            icpk (k) = lhi / cvm (k)
        enddo
    
        ! -----------------------------------------------------------------------
        ! melting of falling cloud ice into rain
        ! -----------------------------------------------------------------------
        
        call check_column (ktop, kbot, qi, no_fall)
        
        if (vi_fac < 1.e-5 .or. no_fall) then
            i1 = 0.
        else
!$acc loop vector
            do k = ktop + 1, kbot
                zt (k) = ze (k) - dtm * (vti (k - 1) + vti (k))/2.0
            enddo
            zt (kbot + 1) = zs - dtm * vti (kbot)
!$acc loop seq
            do k = ktop, kbot
                if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            
            if (k0 < kbot) then
!$acc loop seq
                do k = kbot - 1, k0, - 1
                    if (qi (k) > qcmin) then
                        exit_flag = .true.
!$acc loop seq
                        do m = k + 1, kbot
                            if (zt (k + 1) >= ze (m) .and. exit_flag) exit_flag = .false.
                            if (zt (k) < ze (m + 1) .and. tz (m) > tice .and. exit_flag) then
                                dtime = min (1.0, (ze (m) - ze (m + 1)) / (max (vr_min, vti (k)) * tau_imlt))
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
!$acc loop vector
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
!$acc loop vector
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
!$acc loop vector
            do k = ktop + 1, kbot
                zt (k) = ze (k) - dtm * (vts (k - 1) + vts (k))/2.0
            enddo
            zt (kbot + 1) = zs - dtm * vts (kbot)
!$acc loop seq
            do k = ktop, kbot
                if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            
            if (k0 < kbot) then
!$acc loop seq
                do k = kbot - 1, k0, - 1
                    if (qs (k) > qpmin) then
                        exit_flag = .true.
!$acc loop seq
                        do m = k + 1, kbot
                            if (zt (k + 1) >= ze (m) .and. exit_flag) exit_flag = .false.
                            if(exit_flag) then
                                dtime = min (dtm, (ze (m) - ze (m + 1)) / (vr_min + vts (k)))
                                if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
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
                            endif
                            if (qs (k) < qpmin .and. exit_flag) exit_flag = .false.
                        enddo
                    endif
                enddo
            endif
            
            if (do_sedi_w) then
!$acc loop vector
                do k = ktop, kbot
                    dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
                enddo
            endif
            
            if (use_ppm) then
                call lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, qs, s1, m1, mono_prof)
            else
                call implicit_fall (dtm, ktop, kbot, ze, vts, dp, qs, s1, m1)
            endif
!$acc loop vector
            do k = ktop, kbot
                m1_sol (k) = m1_sol (k) + m1 (k)
            enddo
            
            if (do_sedi_w) then
                w1 (ktop) = (dm (ktop) * w1 (ktop) + m1 (ktop) * vts (ktop)) / (dm (ktop) - m1 (ktop))
!$acc loop vector
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
!$acc loop vector
            do k = ktop + 1, kbot
                zt (k) = ze (k) - dtm * (vtg (k - 1) + vtg (k))/2.0
            enddo
            zt (kbot + 1) = zs - dtm * vtg (kbot)
!$acc loop seq
            do k = ktop, kbot
                if (zt (k + 1) >= zt (k)) zt (k + 1) = zt (k) - dz_min
            enddo
            
            if (k0 < kbot) then
!$acc loop seq
                do k = kbot - 1, k0, - 1
                    if (qg (k) > qpmin) then
                        exit_flag = .true.
!$acc loop seq
                        do m = k + 1, kbot
                            if (zt (k + 1) >= ze (m) .and. exit_flag) exit_flag = .false.
                            if(exit_flag) then
                                dtime = min (dtm, (ze (m) - ze (m + 1)) / vtg (k))
                                if (zt (k) < ze (m + 1) .and. tz (m) > tice) then
                                    dtime = min (1., dtime / tau_g2r)
                                    sink = min (qg (k) * dp (k) / dp (m), dtime * (tz (m) - tice) / icpk (m))
                                    tz (m) = tz (m) - sink * icpk (m)
                                    qg (k) = qg (k) - sink * dp (m) / dp (k)
                                    if (zt (k) < zs) then
                                        r1 = r1 + sink * dp (m)
                                    else
                                        qr (m) = qr (m) + sink
                                    endif
                                endif
                            endif
                            if (qg (k) < qpmin .and. exit_flag) exit_flag = .false.
                        enddo
                    endif
                enddo
            endif
            
            if (do_sedi_w) then
!$acc loop vector
                do k = ktop, kbot
                    dm (k) = dp (k) * (1. + qv (k) + ql (k) + qr (k) + qi (k) + qs (k) + qg (k))
                enddo
            endif
        
            if (use_ppm) then
                call lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, qg, g1, m1, mono_prof)
            else
                call implicit_fall (dtm, ktop, kbot, ze, vtg, dp, qg, g1, m1)
            endif
!$acc loop vector 
            do k = ktop, kbot
                m1_sol (k) = m1_sol (k) + m1 (k)
            enddo
            
            if (do_sedi_w) then
                w1 (ktop) = (dm (ktop) * w1 (ktop) + m1 (ktop) * vtg (ktop)) / (dm (ktop) - m1 (ktop))
!$acc loop vector
                do k = ktop + 1, kbot
                    w1 (k) = (dm (k) * w1 (k) - m1 (k - 1) * vtg (k - 1) + m1 (k) * vtg (k)) &
                        / (dm (k) + m1 (k - 1) - m1 (k))
                enddo
            endif
            
        endif
        
    end subroutine terminal_fall

    subroutine terminal_fall_3d (dtm, is, ie, js, je, ktop, kbot, tz, qv, ql, qr, qg, qs, qi, dz, dp, &
        den, vtg, vts, vti, r1, g1, s1, i1, m1_sol, w1)
    !$acc routine vector
        implicit none

        integer, intent (in) :: is, ie, js, je, ktop, kbot
        
        real, intent (in) :: dtm ! time step (s)
        
        real, intent (in), dimension (is:ie, js:je, ktop:kbot) :: vtg, vts, vti, den, dp, dz
        
        real, intent (inout), dimension (is:ie, js:je, ktop:kbot) :: qv, ql, qr, qg, qs, qi, tz, m1_sol, w1
        
        real, dimension(is:ie, js:je), intent(out) :: r1, g1, s1, i1

        real, dimension (is:ie, js:je,ktop:kbot + 1) :: ze, zt
        
        real :: qsat, dqsdt, evap, dtime
        real :: factor, frac
        real :: tmp, precip, tc, sink
        
        real, dimension (is:ie, js:je,ktop:kbot) :: icpk, cvm
        real, dimension (is:ie, js:je,ktop:kbot) :: m1, dm
        
        real :: q_liq, q_sol, lcpk, lhl, lhi

        real :: zs
        real :: fac_imlt
        
        integer :: i, j, k, m

        integer, dimension(is:ie, js:je) :: k0
        
        logical :: exit_flag

        logical, dimension(is:ie, js:je) :: no_fall

        zs = 0.

        fac_imlt = 1. - exp (- dtm / tau_imlt)
        
        ! -----------------------------------------------------------------------
        ! define heat capacity and latend heat coefficient
        ! -----------------------------------------------------------------------
        !$acc loop vector private(lhi, q_liq, q_sol)
        do k = ktop, kbot
            do j = js, je
                do i = is, ie
                    m1_sol (i,j,k) = 0.
                    ! lhl (k) = lv00 + d0_vap * tz (k)
                    lhi = li00 + dc_ice * tz (i,j,k)
                    q_liq = ql (i,j,k) + qr (i,j,k)
                    q_sol = qi (i,j,k) + qs (i,j,k) + qg (i,j,k)
                    cvm (i,j,k) = c_air + qv (i,j,k) * c_vap + q_liq * c_liq + q_sol * c_ice
                    ! lcpk = lhl (k) / cvm (k)
                    icpk (i,j,k) = lhi / cvm (i,j,k)
                enddo
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! find significant melting level
        ! -----------------------------------------------------------------------
        
        k0 = kbot

        !$acc loop seq
        do k = ktop, kbot - 1
            do j = js, je
                do i = is, ie
                    if (tz (i,j,k) > tice .and. k0(i,j).eq.kbot) then
                        k0(i,j) = k
                    endif
                enddo
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! melting of cloud_ice (before fall) :
        ! -----------------------------------------------------------------------
        !$acc loop vector private(tc, q_liq, q_sol, lhi, sink, tmp)
        do k = ktop, kbot
            do j = js, je
                do i = is, ie
                    if(k.ge.k0(i,j)) then
                        tc = tz (i,j,k) - tice
                        if (qi (i,j,k) > qcmin .and. tc > 0.) then
                            q_liq = ql (i,j,k) + qr (i,j,k)
                            q_sol = qi (i,j,k) + qs (i,j,k) + qg (i,j,k)
                            lhi = li00 + dc_ice * tz (i,j,k)

                            sink = min (qi (i,j,k), fac_imlt * tc / icpk (i,j,k))
                            tmp = min (sink, dim (ql_mlt, ql (i,j,k)))
                            ql (i,j,k) = ql (i,j,k) + tmp
                            qr (i,j,k) = qr (i,j,k) + sink - tmp
                            qi (i,j,k) = qi (i,j,k) - sink
                            q_liq = q_liq + sink
                            q_sol = q_sol - sink
                            cvm (i,j,k) = c_air + qv (i,j,k) * c_vap + q_liq * c_liq + q_sol * c_ice
                            tz (i,j,k) = tz (i,j,k) - sink * lhi / cvm (i,j,k)
                            tc = tz (i,j,k) - tice
                        endif
                    endif
                enddo
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! turn off melting when cloud microphysics time step is small
        ! -----------------------------------------------------------------------
        
        if (dtm < 60.) k0 = kbot
        
        ! sjl, turn off melting of falling cloud ice, snow and graupel
        k0 = kbot
        ! sjl, turn off melting of falling cloud ice, snow and graupel
        
        ze (:,:,kbot + 1) = zs
        !$acc loop seq
        do k = kbot, ktop, - 1
            do j = js, je
                do i = is, ie
                    ze (i,j,k) = ze (i,j,k + 1) - dz (i,j,k) ! dz < 0
                enddo
            enddo
        enddo

        zt (:,:,ktop) = ze (:,:,ktop)
        
        ! -----------------------------------------------------------------------
        ! update capacity heat and latend heat coefficient
        ! -----------------------------------------------------------------------
        !$acc loop vector private(lhi)
        ! do k = k0, kbot
        do k = kbot, kbot
            do j = js, je
                do i = is, ie
                    lhi = li00 + dc_ice * tz (i,j,k)
                    icpk (i,j,k) = lhi / cvm (i,j,k)
                enddo
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! melting of falling cloud ice into rain
        ! -----------------------------------------------------------------------
        
        call check_column_3d (qi, no_fall)
        
        do j = js, je
            do i = is, ie
                if (vi_fac < 1.e-5 .or. no_fall(i,j)) then
                    i1(i,j) = 0.
                else
                !$acc loop vector
                    do k = ktop + 1, kbot
                        zt (i,j,k) = ze (i,j,k) - dtm * (vti (i,j,k - 1) + vti (i,j,k))/2.0
                    enddo
                    zt (i,j,kbot + 1) = zs - dtm * vti (i,j,kbot)
                !$acc loop seq
                    do k = ktop, kbot
                        if (zt (i,j,k + 1) >= zt (i,j,k)) zt (i,j,k + 1) = zt (i,j,k) - dz_min
                    enddo
                    
                    if (k0(i,j) < kbot) then
                !$acc loop seq
                        do k = kbot - 1, k0(i,j), - 1
                            if (qi (i,j,k) > qcmin) then
                                exit_flag = .true.
                                !$acc loop seq
                                do m = k + 1, kbot
                                    if (zt (i,j,k + 1) >= ze (i,j,m) .and. exit_flag) exit_flag = .false.
                                    if (zt (i,j,k) < ze (i,j,m + 1) .and. tz (i,j,m) > tice .and. exit_flag) then
                                        dtime = min (1.0, (ze (i,j,m) - ze (i,j,m + 1)) / (max (vr_min, vti (i,j,k)) * tau_imlt))
                                        sink = min (qi (i,j,k) * dp (i,j,k) / dp (i,j,m), dtime * (tz (i,j,m) - tice) / icpk (i,j,m))
                                        tmp = min (sink, dim (ql_mlt, ql (i,j,m)))
                                        ql (i,j,m) = ql (i,j,m) + tmp
                                        qr (i,j,m) = qr (i,j,m) - tmp + sink
                                        tz (i,j,m) = tz (i,j,m) - sink * icpk (i,j,m)
                                        qi (i,j,k) = qi (i,j,k) - sink * dp (i,j,m) / dp (i,j,k)
                                    endif
                                enddo
                            endif
                        enddo
                    endif
                    
                    if (do_sedi_w) then
                        !$acc loop vector
                        do k = ktop, kbot
                            dm (i,j,k) = dp (i,j,k) * (1. + qv (i,j,k) + ql (i,j,k) + qr (i,j,k) + qi (i,j,k) + qs (i,j,k) + qg (i,j,k))
                        enddo
                    endif
                    
                    if (use_ppm) then
                        call lagrangian_fall_ppm (ktop, kbot, zs, ze(i,j,:), zt(i,j,:), dp(i,j,:), qi(i,j,:), i1(i,j), m1_sol(i,j,:), mono_prof)
                    else
                        call implicit_fall (dtm, ktop, kbot, ze(i,j,:), vti(i,j,:), dp(i,j,:), qi(i,j,:), i1(i,j), m1_sol(i,j,:))
                    endif
                    
                    if (do_sedi_w) then
                        w1 (i,j,ktop) = (dm (i,j,ktop) * w1 (i,j,ktop) + m1_sol (i,j,ktop) * vti (i,j,ktop)) / (dm (i,j,ktop) - m1_sol (i,j,ktop))
                        !$acc loop vector
                        do k = ktop + 1, kbot
                            w1 (i,j,k) = (dm (i,j,k) * w1 (i,j,k) - m1_sol (i,j,k - 1) * vti (i,j,k - 1) + m1_sol (i,j,k) * vti (i,j,k)) &
                                / (dm (i,j,k) + m1_sol (i,j,k - 1) - m1_sol (i,j,k))
                        enddo
                    endif
                    
                endif
            enddo
        enddo

        ! -----------------------------------------------------------------------
        ! melting of falling snow into rain
        ! -----------------------------------------------------------------------
        
        r1 = 0.
        
        call check_column_3d (qs, no_fall)
        
        do j = js, je
            do i = is, ie
                if (no_fall(i,j)) then
                    s1(i,j) = 0.
                else
                    !$acc loop vector
                    do k = ktop + 1, kbot
                        zt (i,j,k) = ze (i,j,k) - dtm * (vts (i,j,k - 1) + vts (i,j,k))/2.0
                    enddo
                    zt (i,j,kbot + 1) = zs - dtm * vts (i,j,kbot)
                    !$acc loop seq
                    do k = ktop, kbot
                        if (zt (i,j,k + 1) >= zt (i,j,k)) zt (i,j,k + 1) = zt (i,j,k) - dz_min
                    enddo
                    
                    if (k0(i,j) < kbot) then
                        !$acc loop seq
                        do k = kbot - 1, k0(i,j), - 1
                            if (qs (i,j,k) > qpmin) then
                                exit_flag = .true.
                                !$acc loop seq
                                do m = k + 1, kbot
                                    if (zt (i,j,k + 1) >= ze (i,j,m) .and. exit_flag) exit_flag = .false.
                                    if(exit_flag) then
                                        dtime = min (dtm, (ze (i,j,m) - ze (i,j,m + 1)) / (vr_min + vts (i,j,k)))
                                        if (zt (i,j,k) < ze (i,j,m + 1) .and. tz (i,j,m) > tice) then
                                            dtime = min (1.0, dtime / tau_smlt)
                                            sink = min (qs (i,j,k) * dp (i,j,k) / dp (i,j,m), dtime * (tz (i,j,m) - tice) / icpk (i,j,m))
                                            tz (i,j,m) = tz (i,j,m) - sink * icpk (i,j,m)
                                            qs (i,j,k) = qs (i,j,k) - sink * dp (i,j,m) / dp (i,j,k)
                                            if (zt (i,j,k) < zs) then
                                                r1(i,j) = r1(i,j) + sink * dp (i,j,m) ! precip as rain
                                            else
                                                ! qr source here will fall next time step (therefore, can evap)
                                                qr (i,j,m) = qr (i,j,m) + sink
                                            endif
                                        endif
                                    endif
                                    if (qs (i,j,k) < qpmin .and. exit_flag) exit_flag = .false.
                                enddo
                            endif
                        enddo
                    endif
                    
                    if (do_sedi_w) then
                        !$acc loop vector
                        do k = ktop, kbot
                            dm (i,j,k) = dp (i,j,k) * (1. + qv (i,j,k) + ql (i,j,k) + qr (i,j,k) + qi (i,j,k) + qs (i,j,k) + qg (i,j,k))
                        enddo
                    endif
                    
                    if (use_ppm) then
                        call lagrangian_fall_ppm (ktop, kbot, zs, ze(i,j,:), zt(i,j,:), dp(i,j,:), qs(i,j,:), s1(i,j), m1(i,j,:), mono_prof)
                    else
                        call implicit_fall (dtm, ktop, kbot, ze(i,j,:), vts(i,j,:), dp(i,j,:), qs(i,j,:), s1(i,j), m1(i,j,:))
                    endif
            !$acc loop vector
                    do k = ktop, kbot
                        m1_sol (i,j,k) = m1_sol (i,j,k) + m1 (i,j,k)
                    enddo
                    
                    if (do_sedi_w) then
                        w1 (i,j,ktop) = (dm (i,j,ktop) * w1 (i,j,ktop) + m1 (i,j,ktop) * vts (i,j,ktop)) / (dm (i,j,ktop) - m1 (i,j,ktop))
            !$acc loop vector
                        do k = ktop + 1, kbot
                            w1 (i,j,k) = (dm (i,j,k) * w1 (i,j,k) - m1 (i,j,k - 1) * vts (i,j,k - 1) + m1 (i,j,k) * vts (i,j,k)) &
                                / (dm (i,j,k) + m1 (i,j,k - 1) - m1 (i,j,k))
                        enddo
                    endif
                    
                endif
            enddo
        enddo

        ! ----------------------------------------------
        ! melting of falling graupel into rain
        ! ----------------------------------------------
        
        call check_column_3d (qg, no_fall)
        
        do j = js, je
            do i = is, ie
                if (no_fall(i,j)) then
                    g1(i,j) = 0.
                else
            !$acc loop vector
                    do k = ktop + 1, kbot
                        zt (i,j,k) = ze (i,j,k) - dtm * (vtg (i,j,k - 1) + vtg (i,j,k))/2.0
                    enddo
                    zt (i,j,kbot + 1) = zs - dtm * vtg (i,j,kbot)
            !$acc loop seq
                    do k = ktop, kbot
                        if (zt (i,j,k + 1) >= zt (i,j,k)) zt (i,j,k + 1) = zt (i,j,k) - dz_min
                    enddo
                    
                    if (k0(i,j) < kbot) then
            !$acc loop seq
                        do k = kbot - 1, k0(i,j), - 1
                            if (qg (i,j,k) > qpmin) then
                                exit_flag = .true.
            !$acc loop seq
                                do m = k + 1, kbot
                                    if (zt (i,j,k + 1) >= ze (i,j,m) .and. exit_flag) exit_flag = .false.
                                    if(exit_flag) then
                                        dtime = min (dtm, (ze (i,j,m) - ze (i,j,m + 1)) / vtg (i,j,k))
                                        if (zt (i,j,k) < ze (i,j,m + 1) .and. tz (i,j,m) > tice) then
                                            dtime = min (1., dtime / tau_g2r)
                                            sink = min (qg (i,j,k) * dp (i,j,k) / dp (i,j,m), dtime * (tz (i,j,m) - tice) / icpk (i,j,m))
                                            tz (i,j,m) = tz (i,j,m) - sink * icpk (i,j,m)
                                            qg (i,j,k) = qg (i,j,k) - sink * dp (i,j,m) / dp (i,j,k)
                                            if (zt (i,j,k) < zs) then
                                                r1(i,j) = r1(i,j) + sink * dp (i,j,m)
                                            else
                                                qr (i,j,m) = qr (i,j,m) + sink
                                            endif
                                        endif
                                    endif
                                    if (qg (i,j,k) < qpmin .and. exit_flag) exit_flag = .false.
                                enddo
                            endif
                        enddo
                    endif
                    
                    if (do_sedi_w) then
            !$acc loop vector
                        do k = ktop, kbot
                            dm (i,j,k) = dp (i,j,k) * (1. + qv (i,j,k) + ql (i,j,k) + qr (i,j,k) + qi (i,j,k) + qs (i,j,k) + qg (i,j,k))
                        enddo
                    endif
                
                    if (use_ppm) then
                        call lagrangian_fall_ppm (ktop, kbot, zs, ze(i,j,:), zt(i,j,:), dp(i,j,:), qg(i,j,:), g1(i,j), m1(i,j,:), mono_prof)
                    else
                        call implicit_fall (dtm, ktop, kbot, ze(i,j,:), vtg(i,j,:), dp(i,j,:), qg(i,j,:), g1(i,j), m1(i,j,:))
                    endif
            !$acc loop vector 
                    do k = ktop, kbot
                        m1_sol (i,j,k) = m1_sol (i,j,k) + m1 (i,j,k)
                    enddo
                    
                    if (do_sedi_w) then
                        w1 (i,j,ktop) = (dm (i,j,ktop) * w1 (i,j,ktop) + m1 (i,j,ktop) * vtg (i,j,ktop)) / (dm (i,j,ktop) - m1 (i,j,ktop))
            !$acc loop vector
                        do k = ktop + 1, kbot
                            w1 (i,j,k) = (dm (i,j,k) * w1 (i,j,k) - m1 (i,j,k - 1) * vtg (i,j,k - 1) + m1 (i,j,k) * vtg (i,j,k)) &
                                / (dm (i,j,k) + m1 (i,j,k - 1) - m1 (i,j,k))
                        enddo
                    endif
                    
                endif
            enddo
        enddo
        
    end subroutine terminal_fall_3d

    ! =======================================================================
    !>@brief The subroutine 'check_column' checks
    !!       if the water species is large enough to fall.
    ! =======================================================================

    subroutine check_column (ktop, kbot, q, no_fall)
    !$acc routine seq
        implicit none
    
        integer, intent (in) :: ktop, kbot
        
        real, intent (in) :: q (ktop:kbot)
        
        logical, intent (out) :: no_fall
        
        integer :: k
        
        no_fall = .true.
        
        do k = ktop, kbot
            if (q (k) > qpmin .and. no_fall) then
                no_fall = .false.
                exit
            endif
        enddo
        
    end subroutine check_column

    subroutine check_column_3d (q, no_fall)
        !$acc routine seq
        implicit none
        
        real, intent (in) :: q (:,:,:)
        
        logical, intent (out) :: no_fall(:,:)
        
        integer :: i,j,k, is, ie, js, je, ktop, kbot

        is = lbound(q,1)
        ie = ubound(q,1)
        js = lbound(q,2)
        je = ubound(q,2)
        ktop = lbound(q,3)
        kbot = ubound(q,3)
        
        no_fall = .true.
        
        do k = ktop, kbot
            do j = js,je
                do i = is, ie
                    if (q (i,j,k) > qpmin .and. no_fall(i,j)) then
                        no_fall(i,j) = .false.
                    endif
                enddo
            enddo
        enddo
        
    end subroutine check_column_3d

    ! =======================================================================
    !>@brief The subroutine 'implicit_fall' computes the time-implicit monotonic 
    !! scheme.
    !>@author Shian-Jiann Lin, 2016
    ! =======================================================================

    subroutine implicit_fall (dt, ktop, kbot, ze, vt, dp, q, precip, m1)
    !$acc routine vector
        implicit none
    
        integer, intent (in) :: ktop, kbot
        
        real, intent (in) :: dt
        
        real, intent (in), dimension (ktop:kbot + 1) :: ze
        
        real, intent (in), dimension (ktop:kbot) :: vt, dp
        
        real, intent (inout), dimension (ktop:kbot) :: q
        
        real, intent (out), dimension (ktop:kbot) :: m1
        
        real, intent (out) :: precip
        
        real, dimension (ktop:kbot) :: qm

        integer :: k
!$acc loop vector
        do k = ktop, kbot
            q (k) = q (k) * dp (k)
        enddo
        
        ! -----------------------------------------------------------------------
        ! sedimentation: non - vectorizable loop
        ! -----------------------------------------------------------------------
        
        qm (ktop) = q (ktop) / ((ze (ktop) - ze (ktop + 1)) + (dt * vt (ktop)))
!$acc loop seq
        do k = ktop + 1, kbot
            qm (k) = (q (k) + (dt * vt (k-1)) * qm (k - 1)) / ((ze (k) - ze (k + 1)) + (dt * vt (k)))
        enddo
        
        ! -----------------------------------------------------------------------
        ! qm is density at this stage
        ! -----------------------------------------------------------------------
!$acc loop vector
        do k = ktop, kbot
            qm (k) = qm (k) * (ze (k) - ze (k + 1))
        enddo
        
        ! -----------------------------------------------------------------------
        ! output mass fluxes: non - vectorizable loop
        ! -----------------------------------------------------------------------
        
        m1 (ktop) = q (ktop) - qm (ktop)
!$acc loop seq
        do k = ktop + 1, kbot
            m1 (k) = m1 (k - 1) + q (k) - qm (k)
        enddo
        precip = m1 (kbot)
        
        ! -----------------------------------------------------------------------
        ! update:
        ! -----------------------------------------------------------------------
!$acc loop vector
        do k = ktop, kbot
            q (k) = qm (k) / dp (k)
        enddo
    end subroutine implicit_fall

    ! =======================================================================
    !> lagrangian scheme
    !  developed by sj lin, ????
    ! =======================================================================

    subroutine lagrangian_fall_ppm (ktop, kbot, zs, ze, zt, dp, q, precip, m1, mono)
        !$acc routine seq
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
    !$acc routine seq
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
    !$acc routine seq
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
                        den, qs, qi, qg, ql, tk, vts, vti, vtg)
    !$acc routine seq
        implicit none

        integer, intent (in) :: ktop, kbot
        
        real, intent (in) :: cnv_fraction, anv_icefall, lsc_icefall 
        ! real, intent (in), dimension (ktop:kbot) :: pl, den, qs, qi, qg, ql, tk
        ! real, intent (out), dimension (ktop:kbot) :: vts, vti, vtg
        real, intent (in)  :: pl, den, qs, qi, qg, ql, tk
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
        
        ! do k = ktop, kbot
        !     rhof (k) = sqrt (min (10., sfcrho / den (k)))
        ! enddo
        
        ! -----------------------------------------------------------------------
        ! ice:
        ! -----------------------------------------------------------------------
        
        if (const_vi) then
            vti = vi_fac
        else
            vi1 = 0.01 * vi_fac
!!$acc loop vector private(tc, IWC, viLSC, viCNV)
            ! do k = ktop, kbot
                if (qi < thi) then ! this is needed as the fall - speed maybe problematic for small qi
                    vti = vf_min
                else
                    tc = tk - tice ! deg C
                    IWC    = qi * den * 1.e3 ! Units are g/m3
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
                    !vti (k) = vti (k) * (C0 + C1*log(DIAM))
                    ! Limits
                    vti = min (vi_max, max (vf_min, vti ))
                endif
            ! enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! snow:
        ! -----------------------------------------------------------------------
        
        if (const_vs) then
            vts = vs_fac ! 1. ifs_2016
        else
!!$acc loop vector private(rhof)
            ! do k = ktop, kbot
                if (qs < ths) then
                    vts = vf_min
                else
                    rhof = sqrt (min (10., sfcrho / den ))
                    vts  = vs_fac * vcons * rhof * exp (0.0625 * log (qs  * den  / norms))
                    vts  = min (vs_max, max (vf_min, vts ))
                endif
            ! enddo
        endif
        
        ! -----------------------------------------------------------------------
        ! graupel:
        ! -----------------------------------------------------------------------
        
        if (const_vg) then
            vtg = vg_fac ! 2.
        else
!!$acc loop vector private(rhof)
            ! do k = ktop, kbot
                if (qg  < thg) then
                    vtg  = vf_min
                else
                    rhof = sqrt (min (10., sfcrho / den ))
                    vtg  = vg_fac * vcong * rhof * sqrt (sqrt (sqrt (qg  * den  / normg)))
                    vtg  = min (vg_max, max (vf_min, vtg ))
                endif
            ! enddo
        endif
        
    end subroutine fall_speed

    ! =======================================================================
    !>@brief The function 'acr3d' is an accretion function (lin et al. 1983)
    ! =======================================================================

    real function acr3d (v1, v2, q1, q2, c, cac, rho)
    !$acc routine seq
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
    !$acc routine seq
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
    !$acc routine seq
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
!$acc routine seq
        implicit none
        
        integer, parameter :: length = 2621
        
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
            
            ! allocate (table (length))
            ! allocate (table2 (length))
            ! allocate (table3 (length))
            ! allocate (tablew (length))
            ! allocate (des (length))
            ! allocate (des2 (length))
            ! allocate (des3 (length))
            ! allocate (desw (length))
            
            ! call qs_table (length)
            ! call qs_table2 (length)
            ! call qs_table3 (length)
            ! call qs_tablew (length)
            
            ! do i = 1, length - 1
            !     des (i) = max (0., table (i + 1) - table (i))
            !     des2 (i) = max (0., table2 (i + 1) - table2 (i))
            !     des3 (i) = max (0., table3 (i + 1) - table3 (i))
            !     desw (i) = max (0., tablew (i + 1) - tablew (i))
            ! enddo
            ! des (length) = des (length - 1)
            ! des2 (length) = des2 (length - 1)
            ! des3 (length) = des3 (length - 1)
            ! desw (length) = desw (length - 1)
            
            ! tables_are_initialized = .true.
            
        endif
        
    end subroutine qsmith_init

    ! =======================================================================
    ! compute the saturated specific humidity for table ii
    !>@brief The function 'wqs1' returns the saturation vapor pressure over pure
    !! liquid water for a given temperature and air density.
    ! =======================================================================

    real function wqs1 (ta, den)
    !$acc routine seq
        implicit none
        
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
    !$acc routine seq
        implicit none
        
        !> pure water phase; universal dry / moist formular using air density
        !> input "den" can be either dry or moist air density
        
        real, intent (in) :: ta, den
        
        real, intent (out) :: dqdt
        
        real :: es, ap1, tmin
        
        integer :: it
        
        tmin = table_ice - 160.
        
        if (.not. tables_are_initialized) call qsmith_init
        
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
    !>@brief The function 'iqs1' computes the saturated specific humidity
    !! for table iii
    ! =======================================================================

    real function iqs1 (ta, den)
    !$acc routine seq
        implicit none
        
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
!$acc routine seq
        implicit none
        
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

    subroutine neg_adj (pt, dp, qv, ql, qr, qi, qs, qg)
    !!$acc routine vector
        implicit none
    
        ! integer, intent (in) :: ktop, kbot
        
        real, intent (in), dimension (:,:,:) :: dp
        
        real, intent (inout), dimension (:,:,:) :: pt, qv, ql, qr, qi, qs, qg
        
        !real, dimension (ktop:kbot) :: lcpk, icpk
        
        integer :: is, ie, js, je, ktop, kbot

        real :: lcpk, icpk

        real :: dq, cvm
        
        integer :: i, j, k

        is = lbound(pt,1)
        ie = ubound(pt,1)

        js = lbound(pt,2)
        je = ubound(pt,2)

        ktop = lbound(pt,3)
        kbot = ubound(pt,3)
        
        ! -----------------------------------------------------------------------
        ! define heat capacity and latent heat coefficient
        ! -----------------------------------------------------------------------
        !$acc parallel loop vector private(cvm, lcpk, icpk)

        do k = ktop, kbot
            do j = js, je
                do i = is, ie
                    cvm = c_air + qv (i,j,k) * c_vap + (qr (i,j,k) + ql (i,j,k)) * c_liq + (qi (i,j,k) + qs (i,j,k) + qg (i,j,k)) * c_ice
                    lcpk  = (lv00 + d0_vap * pt (i,j,k)) / cvm
                    icpk  = (li00 + dc_ice * pt (i,j,k)) / cvm
                    ! -----------------------------------------------------------------------
                    ! ice phase:
                    ! -----------------------------------------------------------------------
            
                    ! if cloud ice < 0, borrow from snow
                    if (qi (i,j,k) < 0.) then
                        qs (i,j,k) = qs (i,j,k) + qi (i,j,k)
                        qi (i,j,k) = 0.
                    endif
                    ! if snow < 0, borrow from graupel
                    if (qs (i,j,k) < 0.) then
                        qg (i,j,k) = qg (i,j,k) + qs (i,j,k)
                        qs (i,j,k) = 0.
                    endif
                    ! if graupel < 0, borrow from rain
                    if (qg (i,j,k) < 0.) then
                        qr (i,j,k) = qr (i,j,k) + qg (i,j,k)
                        pt (i,j,k) = pt (i,j,k) - qg (i,j,k) * icpk ! heating
                        qg (i,j,k) = 0.
                    endif
            
                    ! -----------------------------------------------------------------------
                    ! liquid phase:
                    ! -----------------------------------------------------------------------
                    
                    ! if rain < 0, borrow from cloud water
                    if (qr (i,j,k) < 0.) then
                        ql (i,j,k) = ql (i,j,k) + qr (i,j,k)
                        qr (i,j,k) = 0.
                    endif
                    ! if cloud water < 0, borrow from water vapor
                    if (ql (i,j,k) < 0.) then
                        qv (i,j,k) = qv (i,j,k) + ql (i,j,k)
                        pt (i,j,k) = pt (i,j,k) - ql (i,j,k) * lcpk ! heating
                        ql (i,j,k) = 0.
                    endif
            
                enddo
            enddo
        enddo
        !$acc end parallel loop
        
        ! -----------------------------------------------------------------------
        ! fix water vapor; borrow from below
        ! -----------------------------------------------------------------------
        !$acc parallel loop seq
        do k = ktop, kbot - 1
            do j = js, je
                do i = is, ie
                    if (qv (i,j,k) < 0.) then
                        qv (i,j,k + 1) = qv (i,j,k + 1) + qv (i,j,k) * dp (i,j,k) / dp (i,j,k + 1)
                        qv (i,j,k) = 0.
                    endif
                enddo
            enddo
        enddo
        !$acc end parallel loop
        
        ! -----------------------------------------------------------------------
        ! bottom layer; borrow from above
        ! -----------------------------------------------------------------------
        !$acc parallel loop collapse(2) private(dq)
        do j = js, je
            do i = is, ie
                if (qv (i,j,kbot) < 0. .and. qv (i,j,kbot - 1) > 0.) then
                    dq = min (- qv (i,j,kbot) * dp (i,j,kbot), qv (i,j,kbot - 1) * dp (i,j,kbot - 1))
                    qv (i,j,kbot - 1) = qv (i,j,kbot - 1) - dq / dp (i,j,kbot - 1)
                    qv (i,j,kbot) = qv (i,j,kbot) + dq / dp (i,j,kbot)
                endif
            enddo
        enddo
        !$acc end parallel loop
        
    end subroutine neg_adj

    real function new_ice_condensate(tk, qlk, qik, cnv_fraction, srf_type)
!$acc routine seq
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
!$acc routine seq
        real, intent(in) :: tk, qlk, qik, cnv_fraction, srf_type
        real :: lfrac

        lfrac = 1.0 - ice_fraction(tk,cnv_fraction, srf_type)
        if (qik > qcmin) then
            new_liq_condensate = max(0.0,min(lfrac*(qlk+qik) - qlk, qik))
        else
            new_liq_condensate = 0.0
        endif

    end function new_liq_condensate

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
        qg, qa, qn, dz, is, ie, js, je, ks, ke, ktop, kbot, dt_in, ntimes, &
        rain, snow, graupel, ice, m2_rain, m2_sol, cond, area1, land, &
        cnv_fraction, srf_type, eis, rhcrit, anv_icefall, lsc_icefall, revap, isubl,                 &
        u_dt, v_dt, pt_dt, qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt,   &
        w_var, vt_r, vt_s, vt_g, vt_i, qn2)

        implicit none
    
        logical, intent (in) :: hydrostatic
        
        integer, intent (in) :: is, ie, js, je, ks, ke
        integer, intent (in) :: ntimes, ktop, kbot
        
        real, intent (in) :: dt_in

        real, intent (in), dimension(is:, js:) :: area1, land, cnv_fraction, srf_type, eis

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

        real, dimension(is:ie, js:je, ktop:kbot) :: tz, h_var1d, qvz, qlz, qrz, qiz, qsz, qgz, qaz
        real, dimension(is:ie, js:je, ktop:kbot) :: dp1, dz1
        real, dimension(is:ie, js:je, ktop:kbot) :: w1
        real, dimension(is:ie, js:je, ktop:kbot) :: vtiz, vtsz, vtgz, vtrz
        real, dimension(is:ie, js:je, ktop:kbot) :: qv0, ql0, qr0, qi0, qs0, qg0
        real, dimension(is:ie, js:je, ktop:kbot) :: den, p1, denfac
        real, dimension(is:ie, js:je, ktop:kbot) :: ccn, c_praut, m1_rain, m1_sol, m1, evap1, subl1
        
        real, dimension(is:ie, js:ie) :: r1, s1, i1, g1

        real :: cpaut, rh_adj, rh_rain, t0, den0
        real :: rdt, ccn0
        real :: dts
        real :: s_leng, t_land, t_ocean, h_var
        real :: cvm, tmp, omq
        real :: dqi, qio, qin
        real :: u1_k, u1_km1, v1_k, v1_km1
        
        integer :: i, j, k, n
    
        dts = dt_in / real (ntimes)
        rdt = 1. / dt_in

        ! -----------------------------------------------------------------------
        ! calculate cloud condensation nuclei (ccn)
        ! the following is based on klein eq. 15
        ! -----------------------------------------------------------------------
                
        cpaut = c_paut * 0.104 * grav / 1.717e-5

!$acc data copyin(dts, rdt, cpaut, anv_icefall, lsc_icefall) &
!$acc      present(area1, land, cnv_fraction, srf_type, eis, &
!$acc              rhcrit, uin, vin, delp, pt, dz, &
!$acc              qv, qi, ql, qr, qs, qg, qa, qn, &
!$acc              u_dt, v_dt, w, pt_dt, qa_dt, qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, &
!$acc              revap, isubl, rain, snow, ice, graupel, cond, w_var, &
!$acc              vt_r, vt_s, vt_g, vt_i, qn2, m2_rain, m2_sol) &
!$acc       create(h_var1d, qvz, qlz, qrz, qiz, qsz, qgz, qaz, &
!$acc         vtiz, vtsz, vtgz, vtrz, dp1, dz1, &
!$acc         qv0, ql0, qr0, qi0, qs0, qg0, &
!$acc         den, tz, p1, denfac, &
!$acc         ccn, c_praut, m1_rain, m1_sol, m1, evap1, subl1, &
!$acc         w1, r1, s1, i1, g1, v1_km1, v1_k)
        

        ! -----------------------------------------------------------------------
        ! use local variables
        ! -----------------------------------------------------------------------
        !$acc parallel loop gang vector collapse(3) private(t0, omq, den0)
        do k = ktop, kbot
            do j = js, je
                do i = is, ie

                    m2_rain(i,j,k) = 0.
                    m2_sol(i,j,k) = 0.
                    revap(i,j,k) = 0.
                    isubl(i,j,k) = 0.

                    t0  = pt (i, j, k)
                    tz (i,j,k) = t0
                    dp1 (i,j,k) = delp (i, j, k)
                    ! dp0 (k) = dp1 (k) ! moist air mass * grav
                    
                    ! -----------------------------------------------------------------------
                    ! import horizontal subgrid variability with pressure dependence
                    ! total water subgrid deviation in horizontal direction
                    ! default area dependent form: use dx ~ 100 km as the base
                    ! -----------------------------------------------------------------------
                    h_var1d(i,j,k) = min(0.30,1.0 - rhcrit(i,j,k)) ! restricted to 70%
        
                    ! -----------------------------------------------------------------------
                    ! convert moist mixing ratios to dry mixing ratios
                    ! -----------------------------------------------------------------------
                    
                    qvz (i,j,k) = qv (i, j, k)
                    qlz (i,j,k) = ql (i, j, k)
                    qiz (i,j,k) = qi (i, j, k)
                    qrz (i,j,k) = qr (i, j, k)
                    qsz (i,j,k) = qs (i, j, k)
                    qgz (i,j,k) = qg (i, j, k)
                    
                    ! dp1: dry air_mass
                    ! dp1 (k) = dp1 (k) * (1. - (qvz (k) + qlz (k) + qrz (k) + qiz (k) + qsz (k) + qgz (k)))
                    dp1 (i,j,k) = dp1 (i,j,k) * (1. - qvz (i,j,k)) ! gfs
                    omq = delp (i, j, k) / dp1 (i,j,k)
                    
                    qvz (i,j,k) = qvz (i,j,k) * omq
                    qlz (i,j,k) = qlz (i,j,k) * omq
                    qrz (i,j,k) = qrz (i,j,k) * omq
                    qiz (i,j,k) = qiz (i,j,k) * omq
                    qsz (i,j,k) = qsz (i,j,k) * omq
                    qgz (i,j,k) = qgz (i,j,k) * omq
                    
                    qaz (i,j,k) = qa (i, j, k)

                    den0 = - dp1 (i,j,k) / (grav * dz (i, j, k)) ! density of dry air
                    p1 (i,j,k) = den0 * rdgas * t0 ! dry air pressure
                    
                    ! -----------------------------------------------------------------------
                    ! save a copy of old value for computing tendencies
                    ! -----------------------------------------------------------------------
                    
                    qv0 (i,j,k) = qvz (i,j,k)
                    ql0 (i,j,k) = qlz (i,j,k)
                    qr0 (i,j,k) = qrz (i,j,k)
                    qi0 (i,j,k) = qiz (i,j,k)
                    qs0 (i,j,k) = qsz (i,j,k)
                    qg0 (i,j,k) = qgz (i,j,k)
                    
                    ! -----------------------------------------------------------------------
                    ! for sedi_momentum
                    ! -----------------------------------------------------------------------
                    
                    m1 (i,j,k) = 0.
                    
                    if(do_sedi_w) w1(i,j,k) = w(i, j, k)

                    if (prog_ccn) then
                        ccn (i,j,k) = qn (i, j, k)
                        c_praut (i,j,k) = cpaut * (ccn (i,j,k) * rhor) ** (- 1. / 3.)
                    else
                        ccn (i,j,k) = qn (i, j, k)
                        !!! use GEOS ccn: ccn (k) = (ccn_l * land (i) + ccn_o * (1. - land (i))) * 1.e6
                        c_praut (i,j,k) = cpaut * (ccn (i,j,k) * rhor) ** (- 1. / 3.)
                    endif

                enddo
            enddo
        enddo
        !$acc end parallel loop
        ! -----------------------------------------------------------------------
        ! fix all negative water species
        ! -----------------------------------------------------------------------
        
        if (fix_negative) &
            call neg_adj (tz, dp1, qvz, qlz, qrz, qiz, qsz, qgz)

        do n = 1, ntimes
      
            ! -----------------------------------------------------------------------
            ! dry air density
            ! -----------------------------------------------------------------------
            !$acc parallel loop collapse(3) private(t0)
            do k = ktop, kbot
                do j = js, je
                    do i = is, ie  
                        if (p_nonhydro) then
                            dz1 (i,j,k) = dz (i, j, k)
                            ! den (k) = den0 (k) ! dry air density remains the same
                            den(i,j,k) = - dp1 (i,j,k) / (grav * dz (i, j, k)) ! density of dry air
                            denfac (i,j,k) = sqrt (sfcrho / den (i,j,k))
                        else
                            t0 = pt (i, j, k)
                            dz1 (i,j,k) = dz (i, j, k) * tz (i,j,k) / t0 ! hydrostatic balance
                            ! den (k) = den0 (k) * dz (i, j, k) / dz1 (k)
                            den(i,j,k) = (- dp1 (i,j,k) / (grav * dz (i, j, k))) * dz(i,j,k) / dz1(i,j,k) ! density of dry air
                            denfac (i,j,k) = sqrt (sfcrho / den (i,j,k))
                        endif

                        call fall_speed (ktop, kbot, p1(i,j,k), cnv_fraction(i,j), anv_icefall, lsc_icefall, &
                                    den(i,j,k), qsz(i,j,k), qiz(i,j,k), qgz(i,j,k), qlz(i,j,k), tz(i,j,k), vtsz(i,j,k), vtiz(i,j,k), vtgz(i,j,k))
                    enddo
                enddo
            enddo
            !$acc end parallel loop

            !$acc parallel loop
            do j = js, je
                do i = is, ie  

                    call terminal_fall (dts, ktop, kbot, tz(i,j,:), qvz(i,j,:), qlz(i,j,:), qrz(i,j,:), qgz(i,j,:), qsz(i,j,:), qiz(i,j,:), &
                        dz1(i,j,:), dp1(i,j,:), den(i,j,:), vtgz(i,j,:), vtsz(i,j,:), vtiz(i,j,:), r1(i,j), g1(i,j), s1(i,j), i1(i,j), m1_sol(i,j,:), w1(i,j,:))
                    
                    rain (i,j) = rain (i,j) + r1(i,j) ! from melted snow & ice that reached the ground
                    snow (i,j) = snow (i,j) + s1(i,j)
                    graupel (i,j) = graupel (i,j) + g1(i,j)
                    ice (i,j) = ice (i,j) + i1(i,j)
                    
!                     ! -----------------------------------------------------------------------
                    ! heat transportation during sedimentation
                    ! -----------------------------------------------------------------------
                    
                    if (do_sedi_heat) &
                        call sedi_heat (ktop, kbot, dp1(i,j,:), m1_sol(i,j,:), dz1(i,j,:), tz(i,j,:), qvz(i,j,:), qlz(i,j,:), qrz(i,j,:), qiz(i,j,:), &
                            qsz(i,j,:), qgz(i,j,:), c_ice)
                    
                    ! -----------------------------------------------------------------------
                    ! warm rain processes
                    ! -----------------------------------------------------------------------
                    
                    call warm_rain (dts, ktop, kbot, dp1(i,j,:), dz1(i,j,:), tz(i,j,:), qvz(i,j,:), qlz(i,j,:), qrz(i,j,:), qiz(i,j,:), qsz(i,j,:), &
                        qgz(i,j,:), qaz(i,j,:), eis(i,j), den(i,j,:), denfac(i,j,:), ccn(i,j,:), c_praut(i,j,:), vtrz(i,j,:),   &
                        r1(i,j), evap1(i,j,:), m1_rain(i,j,:), w1(i,j,:), h_var1d(i,j,:))

                    rain (i,j) = rain (i,j) + r1(i,j)
                enddo
            enddo
            !$acc end parallel loop

            !$acc parallel loop collapse(3)
            do k = ktop, kbot
                do j = js, je
                    do i = is, ie  
                        revap (i,j,k) = revap (i,j,k) + evap1(i,j,k)
                        m2_rain (i, j, k) = m2_rain (i, j, k) + m1_rain (i,j,k)
                        m2_sol (i, j, k) = m2_sol (i, j, k) + m1_sol (i,j,k)
                        m1 (i,j,k) = m1 (i,j,k) + m1_rain (i,j,k) + m1_sol (i,j,k)
                    enddo
                enddo
            enddo    
            !$acc end parallel loop             
                    ! -----------------------------------------------------------------------
                    ! ice - phase microphysics
                    ! -----------------------------------------------------------------------
            !$acc parallel loop collapse(2)
            do j = js, je
                do i = is, ie 
                    call icloud (ktop, kbot, tz(i,j,:), p1(i,j,:), qvz(i,j,:), qlz(i,j,:), qrz(i,j,:), qiz(i,j,:), qsz(i,j,:), qgz(i,j,:), dp1(i,j,:), den(i,j,:), &
                        denfac(i,j,:), vtsz(i,j,:), vtgz(i,j,:), vtrz(i,j,:), qaz(i,j,:), dts, subl1(i,j,:), h_var1d(i,j,:), &
                        ccn(i,j,:), cnv_fraction(i,j), srf_type(i,j))
                enddo
            enddo
            !$acc end parallel loop

            !$acc parallel loop collapse(3)
            do k = ktop, kbot
                do j = js, je
                    do i = is, ie
                        isubl (i,j,k) = isubl (i,j,k) + subl1(i,j,k)
                    enddo
                enddo
            enddo
            !$acc end parallel loop
                                
        enddo ! ntimes

        !$acc parallel loop collapse(2)
        do j = js, je
            do i = is, ie
            
                ! -----------------------------------------------------------------------
                ! momentum transportation during sedimentation
                ! note: dp1 is dry mass; dp0 is the old moist (total) mass
                ! -----------------------------------------------------------------------
                
                if (sedi_transport) then
                    v1_km1 = vin (i, j, ktop)
                    u1_km1 = uin (i, j, ktop)
                    !$acc loop seq
                    do k = ktop + 1, kbot
                        u1_k = uin(i, j, k)
                        v1_k = vin(i, j, k)
                        u1_k = (delp (i, j, k) * u1_k + m1 (i,j,k - 1) * u1_km1) / (delp (i, j, k) + m1 (i,j,k - 1))
                        v1_k = (delp (i, j, k) * v1_k + m1 (i,j,k - 1) * v1_km1) / (delp (i, j, k) + m1 (i,j,k - 1))
                        u_dt (i, j, k) = u_dt (i, j, k) + (u1_k - uin (i, j, k)) * rdt
                        v_dt (i, j, k) = v_dt (i, j, k) + (v1_k - vin (i, j, k)) * rdt
                        u1_km1 = u1_k
                        v1_km1 = v1_k
                    enddo
                endif
            enddo
        enddo
        !$acc end parallel loop

        !$acc parallel loop collapse(3) private(t0, omq, cvm)
        do k = ktop,kbot
            do j = js, je
                do i = is, ie
            
                    if (do_sedi_w) then
                            w (i, j, k) = w1 (i,j,k)
                    endif

                    ! -----------------------------------------------------------------------
                    ! update moist air mass (actually hydrostatic pressure)
                    ! convert to dry mixing ratios
                    ! -----------------------------------------------------------------------
                    t0 = pt (i, j, k)
                    omq = dp1 (i,j,k) / delp (i, j, k)
                    qv_dt (i, j, k) = qv_dt (i, j, k) + rdt * (qvz (i,j,k) - qv0 (i,j,k)) * omq
                    ql_dt (i, j, k) = ql_dt (i, j, k) + rdt * (qlz (i,j,k) - ql0 (i,j,k)) * omq
                    qr_dt (i, j, k) = qr_dt (i, j, k) + rdt * (qrz (i,j,k) - qr0 (i,j,k)) * omq
                    qi_dt (i, j, k) = qi_dt (i, j, k) + rdt * (qiz (i,j,k) - qi0 (i,j,k)) * omq
                    qs_dt (i, j, k) = qs_dt (i, j, k) + rdt * (qsz (i,j,k) - qs0 (i,j,k)) * omq
                    qg_dt (i, j, k) = qg_dt (i, j, k) + rdt * (qgz (i,j,k) - qg0 (i,j,k)) * omq
                    cvm = c_air + qvz (i,j,k) * c_vap + (qrz (i,j,k) + qlz (i,j,k)) * c_liq + (qiz (i,j,k) + qsz (i,j,k) + qgz (i,j,k)) * c_ice
                    pt_dt (i, j, k) = pt_dt (i, j, k) + rdt * (tz (i,j,k) - t0) * cvm / cp_air

            
                    ! -----------------------------------------------------------------------
                    ! update cloud fraction tendency
                    ! -----------------------------------------------------------------------
                    if (.not. do_qa) then
                        qa_dt (i, j, k) = qa_dt (i, j, k) + rdt * (                          &
                            qa (i, j, k)*SQRT( (qiz(i,j,k)+qlz(i,j,k)) / max(qi0(i,j,k)+ql0(i,j,k),qcmin) ) - & ! New Cloud -
                            qa (i, j, k) )                                                      ! Old Cloud
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
                ! do k = ktop, kbot
                ! vt_r (i, j, k) = vtrz (k)
                ! enddo
                ! endif
                !
                ! if (id_vts > 0) then
                ! do k = ktop, kbot
                ! vt_s (i, j, k) = vtsz (k)
                ! enddo
                ! endif
                !
                ! if (id_vtg > 0) then
                ! do k = ktop, kbot
                ! vt_g (i, j, k) = vtgz (k)
                ! enddo
                ! endif
                !
                ! if (id_vts > 0) then
                ! do k = ktop, kbot
                ! vt_i (i, j, k) = vtiz (k)
                ! enddo
                ! endif
                !
                ! if (id_droplets > 0) then
                ! do k = ktop, kbot
                ! qn2 (i, j, k) = ccn (k)
                ! enddo
                ! endif
            enddo
            enddo
        enddo
        !$acc end parallel loop
!$acc end data
    end subroutine mpdrv

    subroutine gfdl_cloud_microphys_driver (qv, ql, qr, qi, qs, qg, qa, qn,   &
        qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, pt_dt, pt, w,    &
        uin, vin, udt, vdt, dz, delp, area, dt_in,                        &
        land, cnv_fraction, srf_type, eis,                                &
        rhcrit, anv_icefall, lsc_icefall,                                 &
        revap, isubl,                                                     &
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

        lcp = latv / cp_air
        icp = lati / cp_air
        tcp = (latv + lati) / cp_air
!$acc update device(c_air, c_vap, p_nonhydro, d0_vap, lv00, do_sedi_w, latv, lati, &
!$acc               lats, lat2, lcp, icp, tcp)

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

!$acc data copyin(uin, vin, delp, pt, qv, ql, qr, qi, qs, qg, &
!$acc             qa, qn, dz, area, land, cnv_fraction, srf_type, eis, &
!$acc             rhcrit, ntimes) &
!$acc      copyout(m2_rain, m2_sol, revap, isubl, w_var, vt_r, vt_s, vt_g, &
!$acc              vt_i, qn2) &
!$acc      copy(w, rain, snow, graupel, ice, cond, udt, vdt, pt_dt, qv_dt, &
!$acc           ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt)

        call mpdrv (hydrostatic, uin, vin, w, delp, pt, qv, ql, qr, qi, qs, qg,&
            qa, qn, dz, is, ie, js, je, ks, ke, ktop, kbot, dt_in, ntimes,  &
            rain, snow, graupel, ice, m2_rain,     &
            m2_sol, cond, area,                                  &
            land, cnv_fraction, srf_type, eis,         &
            rhcrit, anv_icefall, lsc_icefall,                                  &
            revap, isubl,                                                      &
            udt, vdt, pt_dt,                                                   &
            qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, w_var, vt_r,      &
            vt_s, vt_g, vt_i, qn2)

!$acc end data
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

    subroutine update_microphys_constants(dirName, rank_str)
    
        implicit none

        integer :: fileID, N

        character*100, intent(in) :: dirName, rank_str

        N = 2621

        allocate(table(N))
        allocate(table2(N))
        allocate(table3(N))
        allocate(tablew(N))
        allocate(des(N))
        allocate(des2(N))
        allocate(des3(N))
        allocate(desw(N))

        open(newunit=fileID, file=trim(dirName) // '/cracs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cracs
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cracs) = ', sum(cracs)

        open(newunit=fileID, file=trim(dirName) // '/csacr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) csacr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(csacr) = ', sum(csacr)

        open(newunit=fileID, file=trim(dirName) // '/cgacr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgacr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgacr) = ', sum(cgacr)

        open(newunit=fileID, file=trim(dirName) // '/cgacs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgacs
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgacs) = ', sum(cgacs)

        open(newunit=fileID, file=trim(dirName) // '/csacw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) csacw
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(csacw) = ', sum(csacw)

        open(newunit=fileID, file=trim(dirName) // '/craci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) craci
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(craci) = ', sum(craci)

        open(newunit=fileID, file=trim(dirName) // '/csaci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) csaci
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(csaci) = ', sum(csaci)

        open(newunit=fileID, file=trim(dirName) // '/cgacw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgacw
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgacw) = ', sum(cgacw)

        open(newunit=fileID, file=trim(dirName) // '/cgaci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgaci
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgaci) = ', sum(cgaci)

        open(newunit=fileID, file=trim(dirName) // '/cracw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cracw
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cracw) = ', sum(cracw)

        open(newunit=fileID, file=trim(dirName) // '/acco_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) acco
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(acco) = ', sum(acco)

        open(newunit=fileID, file=trim(dirName) // '/cssub_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cssub
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cssub) = ', sum(cssub)

        open(newunit=fileID, file=trim(dirName) // '/cgsub_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgsub
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgsub) = ', sum(cgsub)

        open(newunit=fileID, file=trim(dirName) // '/crevp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) crevp
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(crevp) = ', sum(crevp)

        open(newunit=fileID, file=trim(dirName) // '/cgfr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgfr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgfr) = ', sum(cgfr)

        open(newunit=fileID, file=trim(dirName) // '/csmlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) csmlt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(csmlt) = ', sum(csmlt)

        open(newunit=fileID, file=trim(dirName) // '/cgmlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cgmlt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cgmlt) = ', sum(cgmlt)

        open(newunit=fileID, file=trim(dirName) // '/es0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) es0
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(es0) = ', sum(es0)

        open(newunit=fileID, file=trim(dirName) // '/ces0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ces0
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(ces0) = ', sum(ces0)

        open(newunit=fileID, file=trim(dirName) // '/pie_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) pie
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(pie) = ', sum(pie)

        open(newunit=fileID, file=trim(dirName) // '/rgrav_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rgrav
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(rgrav) = ', sum(rgrav)

        open(newunit=fileID, file=trim(dirName) // '/c_air_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_air
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In c_air = ', c_air

        open(newunit=fileID, file=trim(dirName) // '/c_vap_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_vap
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In c_vap = ', c_vap

        open(newunit=fileID, file=trim(dirName) // '/lati_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) lati
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(lati) = ', sum(lati)

        open(newunit=fileID, file=trim(dirName) // '/latv_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) latv
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(latv) = ', sum(latv)

        open(newunit=fileID, file=trim(dirName) // '/lats_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) lats
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(lats) = ', sum(lats)

        open(newunit=fileID, file=trim(dirName) // '/lat2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) lat2
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(lat2) = ', sum(lat2)

        open(newunit=fileID, file=trim(dirName) // '/lcp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) lcp
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(lcp) = ', sum(lcp)

        open(newunit=fileID, file=trim(dirName) // '/icp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) icp
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(icp) = ', sum(icp)

        open(newunit=fileID, file=trim(dirName) // '/tcp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tcp
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tcp) = ', sum(tcp)

        open(newunit=fileID, file=trim(dirName) // '/d0_vap_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) d0_vap
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(d0_vap) = ', sum(d0_vap)

        open(newunit=fileID, file=trim(dirName) // '/lv00_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) lv00
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(lv00) = ', sum(lv00)

        open(newunit=fileID, file=trim(dirName) // '/icloud_f_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) icloud_f
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(icloud_f) = ', sum(icloud_f)

        open(newunit=fileID, file=trim(dirName) // '/irain_f_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) irain_f
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(irain_f) = ', sum(irain_f)

        open(newunit=fileID, file=trim(dirName) // '/de_ice_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) de_ice
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(de_ice) = ', sum(de_ice)

        open(newunit=fileID, file=trim(dirName) // '/sedi_transport_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) sedi_transport
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(sedi_transport) = ', sum(sedi_transport)

        open(newunit=fileID, file=trim(dirName) // '/do_sedi_w_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_sedi_w
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_sedi_w) = ', sum(do_sedi_w)

        open(newunit=fileID, file=trim(dirName) // '/do_sedi_heat_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_sedi_heat
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_sedi_heat) = ', sum(do_sedi_heat)

        open(newunit=fileID, file=trim(dirName) // '/prog_ccn_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) prog_ccn
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(prog_ccn) = ', sum(prog_ccn)

        open(newunit=fileID, file=trim(dirName) // '/do_bigg_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_bigg
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_bigg) = ', sum(do_bigg)

        open(newunit=fileID, file=trim(dirName) // '/do_evap_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_evap
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_evap) = ', sum(do_evap)

        open(newunit=fileID, file=trim(dirName) // '/do_subl_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_subl
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_subl) = ', sum(do_subl)

        open(newunit=fileID, file=trim(dirName) // '/do_qa_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_qa
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_qa) = ', sum(do_qa)

        open(newunit=fileID, file=trim(dirName) // '/preciprad_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) preciprad
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(preciprad) = ', sum(preciprad)

        open(newunit=fileID, file=trim(dirName) // '/fix_negative_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) fix_negative
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(fix_negative) = ', sum(fix_negative)

        open(newunit=fileID, file=trim(dirName) // '/do_setup_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_setup
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(do_setup) = ', sum(do_setup)

        open(newunit=fileID, file=trim(dirName) // '/p_nonhydro_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) p_nonhydro
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(p_nonhydro) = ', sum(p_nonhydro)

        open(newunit=fileID, file=trim(dirName) // '/table_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) table
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(table) = ', sum(table)

        open(newunit=fileID, file=trim(dirName) // '/table2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) table2
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(table2) = ', sum(table2)

        open(newunit=fileID, file=trim(dirName) // '/table3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) table3
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(table3) = ', sum(table3)

        open(newunit=fileID, file=trim(dirName) // '/tablew_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tablew
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tablew) = ', sum(tablew)

        open(newunit=fileID, file=trim(dirName) // '/des_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) des
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(des) = ', sum(des)

        open(newunit=fileID, file=trim(dirName) // '/des2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) des2
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(des2) = ', sum(des2)

        open(newunit=fileID, file=trim(dirName) // '/des3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) des3
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(des3) = ', sum(des3)

        open(newunit=fileID, file=trim(dirName) // '/desw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) desw
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(desw) = ', sum(desw)

        open(newunit=fileID, file=trim(dirName) // '/dt_fr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dt_fr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(dt_fr) = ', sum(tables_are_idt_frnitialized)

        open(newunit=fileID, file=trim(dirName) // '/p_min_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) p_min
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(p_min) = ', sum(p_min)

        open(newunit=fileID, file=trim(dirName) // '/cld_min_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cld_min
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(cld_min) = ', sum(cld_min)

        open(newunit=fileID, file=trim(dirName) // '/tice_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tice
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tice) = ', sum(tice)

        open(newunit=fileID, file=trim(dirName) // '/log_10_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) log_10
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(log_10) = ', sum(log_10)

        open(newunit=fileID, file=trim(dirName) // '/tice0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tice0
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tice0) = ', sum(tice0)

        open(newunit=fileID, file=trim(dirName) // '/t_wfr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) t_wfr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(t_wfr) = ', sum(t_wfr)

        open(newunit=fileID, file=trim(dirName) // '/t_min_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) t_min
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(t_min) = ', sum(t_min)

        open(newunit=fileID, file=trim(dirName) // '/t_sub_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) t_sub
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(t_sub) = ', sum(t_sub)

        open(newunit=fileID, file=trim(dirName) // '/mp_time_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mp_time
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(mp_time) = ', sum(mp_time)

        open(newunit=fileID, file=trim(dirName) // '/rh_inc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rh_inc
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(rh_inc) = ', sum(rh_inc)

        open(newunit=fileID, file=trim(dirName) // '/rh_inr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rh_inr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(rh_inr) = ', sum(rh_inr)

        open(newunit=fileID, file=trim(dirName) // '/rh_ins_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rh_ins
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(rh_ins) = ', sum(rh_ins)

        open(newunit=fileID, file=trim(dirName) // '/tau_r2g_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_r2g
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_r2g) = ', sum(tau_r2g)

        open(newunit=fileID, file=trim(dirName) // '/tau_smlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_smlt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_smlt) = ', sum(tau_smlt)

        open(newunit=fileID, file=trim(dirName) // '/tau_g2r_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_g2r
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_g2r) = ', sum(tau_g2r)

        open(newunit=fileID, file=trim(dirName) // '/tau_imlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_imlt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_imlt) = ', sum(tau_imlt)

        open(newunit=fileID, file=trim(dirName) // '/tau_i2s_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_i2s
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_i2s) = ', sum(tau_i2s)

        open(newunit=fileID, file=trim(dirName) // '/tau_l2r_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_l2r
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_l2r) = ', sum(tau_l2r)

        open(newunit=fileID, file=trim(dirName) // '/tau_v2l_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_v2l
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_v2l) = ', sum(tau_v2l)

        open(newunit=fileID, file=trim(dirName) // '/tau_l2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_l2v
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_l2v) = ', sum(tau_l2v)

        open(newunit=fileID, file=trim(dirName) // '/tau_i2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_i2v
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_i2v) = ', sum(tau_i2v)

        open(newunit=fileID, file=trim(dirName) // '/tau_s2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_s2v
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_s2v) = ', sum(tau_s2v)

        open(newunit=fileID, file=trim(dirName) // '/tau_v2s_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_v2s
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_v2s) = ', sum(tau_v2s)

        open(newunit=fileID, file=trim(dirName) // '/tau_g2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_g2v
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_g2v) = ', sum(tau_g2v)

        open(newunit=fileID, file=trim(dirName) // '/tau_v2g_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_v2g
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_v2g) = ', sum(tau_v2g)

        open(newunit=fileID, file=trim(dirName) // '/tau_revp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_revp
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_revp) = ', sum(tau_revp)

        open(newunit=fileID, file=trim(dirName) // '/tau_frz_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_frz
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(tau_frz) = ', sum(tau_frz)

        open(newunit=fileID, file=trim(dirName) // '/dw_land_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dw_land
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(dw_land) = ', sum(dw_land)

        open(newunit=fileID, file=trim(dirName) // '/dw_ocean_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dw_ocean
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(dw_ocean) = ', sum(dw_ocean)

        open(newunit=fileID, file=trim(dirName) // '/ccn_o_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ccn_o
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(ccn_o) = ', sum(ccn_o)

        open(newunit=fileID, file=trim(dirName) // '/ccn_l_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ccn_l
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(ccn_l) = ', sum(ccn_l)

        open(newunit=fileID, file=trim(dirName) // '/rthreshu_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rthreshu
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(rthreshu) = ', sum(rthreshu)

        open(newunit=fileID, file=trim(dirName) // '/rthreshs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rthreshs
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(rthreshs) = ', sum(rthreshs)

        open(newunit=fileID, file=trim(dirName) // '/sat_adj0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) sat_adj0
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(sat_adj0) = ', sum(sat_adj0)

        open(newunit=fileID, file=trim(dirName) // '/qc_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qc_crt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qc_crt) = ', sum(qc_crt)

        open(newunit=fileID, file=trim(dirName) // '/qi_lim_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi_lim
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qi_lim) = ', sum(qi_lim)

        open(newunit=fileID, file=trim(dirName) // '/ql_mlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ql_mlt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(ql_mlt) = ', sum(ql_mlt)

        open(newunit=fileID, file=trim(dirName) // '/qs_mlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qs_mlt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qs_mlt) = ', sum(qs_mlt)

        open(newunit=fileID, file=trim(dirName) // '/ql_gen_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ql_gen
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(ql_gen) = ', sum(ql_gen)

        open(newunit=fileID, file=trim(dirName) // '/qi_gen_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi_gen
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qi_gen) = ', sum(qi_gen)

        open(newunit=fileID, file=trim(dirName) // '/ql0_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ql0_max
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(ql0_max) = ', sum(ql0_max)

        open(newunit=fileID, file=trim(dirName) // '/qi0_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi0_max
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qi0_max) = ', sum(qi0_max)

        open(newunit=fileID, file=trim(dirName) // '/qi0_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi0_crt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qi0_crt) = ', sum(qi0_crt)

        open(newunit=fileID, file=trim(dirName) // '/qr0_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qr0_crt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qr0_crt) = ', sum(qr0_crt)

        open(newunit=fileID, file=trim(dirName) // '/qs0_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qs0_crt
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(qs0_crt) = ', sum(qs0_crt)

        open(newunit=fileID, file=trim(dirName) // '/c_paut_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_paut
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(c_paut) = ', sum(c_paut)

        open(newunit=fileID, file=trim(dirName) // '/c_psaci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_psaci
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(c_psaci) = ', sum(c_psaci)

        open(newunit=fileID, file=trim(dirName) // '/c_piacr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_piacr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(c_piacr) = ', sum(c_piacr)

        open(newunit=fileID, file=trim(dirName) // '/c_cracw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_cracw
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(c_cracw) = ', sum(c_cracw)

        open(newunit=fileID, file=trim(dirName) // '/c_pgacs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_pgacs
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(c_pgacs) = ', sum(c_pgacs)

        open(newunit=fileID, file=trim(dirName) // '/c_pgaci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_pgaci
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(c_pgaci) = ', sum(c_pgaci)

        open(newunit=fileID, file=trim(dirName) // '/alin_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) alin
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(alin) = ', sum(alin)

        open(newunit=fileID, file=trim(dirName) // '/clin_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) clin
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(clin) = ', sum(clin)

        open(newunit=fileID, file=trim(dirName) // '/const_vi_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vi
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(const_vi) = ', sum(const_vi)

        open(newunit=fileID, file=trim(dirName) // '/const_vs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vs
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(const_vs) = ', sum(const_vs)

        open(newunit=fileID, file=trim(dirName) // '/const_vg_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vg
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(const_vg) = ', sum(const_vg)

        open(newunit=fileID, file=trim(dirName) // '/const_vr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vr
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(const_vr) = ', sum(const_vr)

        open(newunit=fileID, file=trim(dirName) // '/vi_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vi_fac
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vi_fac) = ', sum(vi_fac)

        open(newunit=fileID, file=trim(dirName) // '/vs_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vs_fac
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vs_fac) = ', sum(vs_fac)

        open(newunit=fileID, file=trim(dirName) // '/vg_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vg_fac
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vg_fac) = ', sum(vg_fac)

        open(newunit=fileID, file=trim(dirName) // '/vr_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vr_fac
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vr_fac) = ', sum(vr_fac)

        open(newunit=fileID, file=trim(dirName) // '/vi_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vi_max
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vi_max) = ', sum(vi_max)

        open(newunit=fileID, file=trim(dirName) // '/vs_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vs_max
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vs_max) = ', sum(vs_max)

        open(newunit=fileID, file=trim(dirName) // '/vg_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vg_max
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vg_max) = ', sum(vg_max)

        open(newunit=fileID, file=trim(dirName) // '/vr_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vr_max
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(vr_max) = ', sum(vr_max)

        open(newunit=fileID, file=trim(dirName) // '/fast_sat_adj_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) fast_sat_adj
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(fast_sat_adj) = ', sum(fast_sat_adj)

        open(newunit=fileID, file=trim(dirName) // '/z_slope_liq_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) z_slope_liq
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(z_slope_liq) = ', sum(z_slope_liq)

        open(newunit=fileID, file=trim(dirName) // '/z_slope_ice_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) z_slope_ice
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(z_slope_ice) = ', sum(z_slope_ice)

        open(newunit=fileID, file=trim(dirName) // '/use_ccn_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) use_ccn
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(use_ccn) = ', sum(use_ccn)

        open(newunit=fileID, file=trim(dirName) // '/use_ppm_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) use_ppm
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(use_ppm) = ', sum(use_ppm)

        open(newunit=fileID, file=trim(dirName) // '/mono_prof_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mono_prof
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(mono_prof) = ', sum(mono_prof)

        open(newunit=fileID, file=trim(dirName) // '/mp_print_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mp_print
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In sum(mp_print) = ', sum(mp_print)

!$acc update device(cracs, csacr, cgacr, cgacs, csacw, craci, csaci, cgacw, cgaci, cracw, &
!$acc               acco, cssub, cgsub, crevp, cgfr, csmlt, cgmlt, es0, ces0, pie, rgrav, &
!$acc               c_air, c_vap, lati, latv, lats, lat2, lcp, icp, tcp, d0_vap, &
!$acc               lv00, icloud_f, irain_f, de_ice, sedi_transport, do_sedi_w, do_sedi_heat, &
!$acc               prog_ccn, do_bigg, do_evap, do_subl, do_qa, preciprad, fix_negative, do_setup, p_nonhydro, table, &
!$acc               table2, table3, tablew, des, des2, des3, desw, dt_fr, p_min, cld_min, &
!$acc               tice, log_10, tice0, t_wfr, t_min, t_sub, mp_time, rh_inc, rh_inr, rh_ins, &
!$acc               tau_r2g, tau_smlt, tau_g2r, tau_imlt, tau_i2s, tau_l2r, tau_v2l, tau_l2v, &
!$acc               tau_i2v, tau_s2v, tau_v2s, tau_g2v, tau_v2g, tau_revp, tau_frz, dw_land, dw_ocean, ccn_o, &
!$acc               ccn_l, rthreshu, rthreshs, sat_adj0, qc_crt, qi_lim, ql_mlt, qs_mlt, ql_gen, qi_gen, &
!$acc               ql0_max, qi0_max, qi0_crt, qr0_crt, qs0_crt, c_paut, c_psaci, c_piacr, &
!$acc               c_cracw, c_pgacs, c_pgaci, alin, clin, const_vi, const_vs, const_vg, const_vr, &
!$acc               vi_fac, vs_fac, vg_fac, vr_fac, vi_max, vs_max, vg_max, vr_max, fast_sat_adj, &
!$acc               z_slope_liq, z_slope_ice, use_ccn, use_ppm, mono_prof, mp_print)

    end subroutine

end module

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.