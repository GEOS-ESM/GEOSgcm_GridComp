      real(kind=kind_phys), dimension(IMJM) :: cnvprcp
      real(kind=kind_phys), dimension(IMJM) :: cnvprcpb
      real(kind=kind_phys), dimension(IMJM) :: totprcp
      real(kind=kind_phys), dimension(IMJM) :: totice
      real(kind=kind_phys), dimension(IMJM) :: totsnw
      real(kind=kind_phys), dimension(IMJM) :: totgrp
      real(kind=kind_phys), dimension(IMJM) :: totprcpb
      real(kind=kind_phys), dimension(IMJM) :: toticeb
      real(kind=kind_phys), dimension(IMJM) :: totsnwb
      real(kind=kind_phys), dimension(IMJM) :: totgrpb
      real(kind=kind_phys), dimension(IMJM) :: tdomr
      real(kind=kind_phys), dimension(IMJM) :: tdomzr
      real(kind=kind_phys), dimension(IMJM) :: tdomip
      real(kind=kind_phys), dimension(IMJM) :: tdoms

      real(kind=kind_phys), dimension(IMJM) ::  cldwrk
      ! dt3dt, dq3dt, du3dt, dv3dt upd_mf, dwn_mf, det_mf only allocated if ldiag3d == .true.
      real(kind=kind_phys), dimension(IMJM,LM,8) :: du3dt
      real(kind=kind_phys), dimension(IMJM,LM,8) :: dv3dt
      real(kind=kind_phys), dimension(IMJM,LM,11) :: dt3dt
      real(kind=kind_phys), dimension(IMJM,LM,13) :: dq3dt

      real(kind=kind_phys), dimension(IMJM)  :: rainc

      integer, parameter    :: ntrac=9, nn=ntrac-2
      integer    :: nncl
      integer, parameter    :: ntk=4, nsamftrac=2, nrcm=2, tracers_total=2
      integer, parameter    :: itc=0, ntchm=0
      integer, parameter    :: ntqv=1, ntcw=2, ntiw=3, ntclamt=8, ntrw=5, ntsw=6, & 
                               ntrnc=-99, ntsnc=-99, ntgl=7, ntgnc=-99, ntinc=-99, &
                               ntlnc=-99, ntke=9, ntwa=-99, ntia=-99, ntchs=-99, &
                               ntoz=4
      integer, parameter    :: ntot3d=6, ncnvw=6, ncnvc=-999, ncnvcld3d=1
      integer, parameter    :: nleffr=1, nieffr=2, nreffr=3, nseffr=4, ngeffr=5
      integer, parameter    :: npdf3d=0, num_p3d=5
      integer, parameter    :: imp_physics_mg=10, imp_physics_zhao_carr=99, &
                               imp_physics_zhao_carr_pdf=98,  &
                               imp_physics_gfdl=11, imp_physics_thompson=8, &
                               imp_physics_wsm6=6,  imp_physics_fer_hires=15
      logical, parameter  :: flag_for_dcnv_generic_tend = .true.
      logical, parameter  :: flag_for_scnv_generic_tend = .true.
      integer, parameter    :: lsm_noah=1, lsm_noahmp=2, lsm_ruc=3, lsm_noah_wrfv4=4
      integer, parameter  :: imfshalcnv_samf=2, imfshalcnv_sas=1
      !real(kind=kind_phys)  :: tem
      !real(kind=kind_phys), save  :: dxmin, dxmax, dxinv
      real(kind=kind_phys), dimension(IMJM) :: pgr

      integer             , dimension(IMJM) :: islmsk
      real(kind=kind_phys), dimension(IMJM) :: work1, work2, psurf
      real(kind=kind_phys), dimension(IMJM,LM) :: dudt, dvdt, dtdt, dtdtc
      real(kind=kind_phys), dimension(IMJM,LM,ntrac) ::  dqdt
      real(kind=kind_phys), parameter   :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      !character(len=ESMF_MAXSTR)      :: errmsg
      integer               :: errflg

      integer, dimension(IMJM)  :: kpbl, kinver

      real(kind=kind_phys)      :: rhcbot, rhcpbl, rhctop
      real(kind=kind_phys), dimension(IMJM, LM)         :: prsl, prslk
      real(kind=kind_phys), dimension(IMJM, LM+1)       :: prsi
      !real(kind=kind_phys), dimension(IMJM)               :: xlon, xlat
      real(kind=kind_phys), dimension(IMJM, LM)         :: gt0
      real(kind=kind_phys), dimension(IMJM, LM, ntrac)  :: gq0

      real(kind=kind_phys), dimension(IMJM, LM)       :: rhc !, save_qc
      ! save_qi is not allocated for Zhao-Carr MP
      !real(kind=kind_phys), dimension(IMJM, LM)          :: save_qi
      real(kind=kind_phys), dimension(IMJM, LM)          :: save_tcp ! ONLY ALLOCATE FOR THOMPSON! TODO
      real(kind=kind_phys), dimension(IMJM, LM, nn)   :: clw

      !logical  :: ldiag3d, qdiag3d, do_cnvgwd, cplchm
      !real(kind=kind_phys), dimension(IMJM,LM)  :: gq0_water_vapor
      real(kind=kind_phys), dimension(IMJM,LM)  :: gu0, gv0
      real(kind=kind_phys), dimension(IMJM,LM)  :: save_u, save_v, save_t !, save_qv
      ! dqdti only allocated if cplchm is .true.
      real(kind=kind_phys), dimension(IMJM,LM)      :: dqdti

      !logical  :: lssav

      !real(kind=kind_phys)  :: frain
      real(kind=kind_phys), dimension(IMJM)  :: rain1, cld1d
      real(kind=kind_phys), dimension(IMJM,LM)  :: ud_mf, dd_mf, dt_mf

      real(kind=kind_phys), dimension(IMJM,LM)  :: upd_mf, dwn_mf, det_mf
      real(kind=kind_phys), dimension(IMJM,LM)  :: cnvw, cnvc
      ! The following arrays may not be allocated, depending on certain flags and microphysics schemes.
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated (it is still invalid Fortran code, though).
      real(kind=kind_phys), dimension(IMJM,LM,ntot3d)  :: phy_f3d


      integer   :: islimsk(IMJM)
     ! real(kind=kind_phys)  :: con_cliq, con_cp, con_cvap, con_eps, con_epsm1,   &
     !&   con_fvirt, con_g, con_hvap, con_rd, con_rv, con_t0c
      real(kind=kind_phys)  ::  delt
      real(kind=kind_phys)  :: delp(IMJM,LM),         &
     &   phil(IMJM,LM), dot(IMJM,LM)
     !&   garea(IMJM), dot(IMJM,LM), phil(IMJM,LM)
      real(kind=kind_phys), dimension(ntchm)  :: fscav
      !logical   :: hwrf_samfdeep
      !real(kind=kind_phys)  :: nthresh
      real(kind=kind_phys), dimension(IMJM)   :: ca_deep
      real(kind=kind_phys), dimension(IMJM)  :: condition
      !logical   :: do_ca,ca_closure,ca_entr,ca_trigger

      integer   :: kcnv(IMJM)
      ! DH* TODO - check dimensions of qtr, ntr+2 correct?  *DH
      !real(kind=kind_phys)  ::   qtr(IMJM,LM,ntr+2),        &
     !&   q1(IMJM,LM), t1(IMJM,LM),   u1(IMJM,LM), v1(IMJM,LM)
      !real(kind=kind_phys)  ::   t1(IMJM,LM),   u1(IMJM,LM), v1(IMJM,LM)

      integer  :: kbot(IMJM), ktop(IMJM)
      real(kind=kind_phys)  :: raincd(IMJM), raincs(IMJM)

      ! GJF* These variables are conditionally allocated depending on whether the
      !     Morrison-Gettelman microphysics is used, so they must be declared
      !     using assumed shape.
      real(kind=kind_phys), dimension(IMJM,LM)  ::            &
     &   qlcn_mg, qicn_mg, w_upi, cnv_mfd, cnv_dqldt, clcn_mg                    &
     &,  cnv_fice, cnv_ndrop, cnv_nice, cf_upi
      ! *GJF
      integer :: mp_phys, mp_phys_mg

     ! real(kind=kind_phys)  :: clam_deep,  c0s_deep,  c1_deep,              &
     !&                     betal_deep,   betas_deep,   asolfac_deep,        &
     !&                     evfact_deep,  evfactl_deep, pgcon_deep


      ! Post code for SAS/SAMF
      !logical  :: shcnvcw
      ! The following arrays may not be allocated, depending on certain flags and microphysics schemes.
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated.


      real(kind=kind_phys)  :: hpbl(IMJM)
!
      !logical            :: hwrf_samfshal

      !integer                                   :: tracers_total, &
      !  ntlnc, ntinc

      !real(kind=kind_phys)                      :: con_pi

      real(kind=kind_phys), dimension(IMJM,LM,ntrac)  :: qgrs

      !logical                                         :: do_aw

      real(kind=kind_phys), dimension(IMJM,LM,ntrac)  :: save_q


      !integer  :: kdt
      !logical  :: cal_pre, cplflx

      real(kind=kind_phys), dimension(IMJM)                :: tsfc
      real(kind=kind_phys), dimension(IMJM)             :: ice, snow, graupel
      real(kind=kind_phys), dimension(IMJM)                 :: rainmp, icemp, snowmp, graupelmp ! conditionally allocated in GFS_typedefs (imp_physics == GFDL or Thompson)
      real(kind=kind_phys), dimension(IMJM,nrcm)           :: rann
      real(kind=kind_phys), dimension(IMJM,LM)           :: del
      real(kind=kind_phys), dimension(IMJM,LM+1)         :: phii

      real(kind=kind_phys), dimension(IMJM)       :: sr
      real(kind=kind_phys), dimension(IMJM)       :: rain, tprcp, srflag, pwat
      real(kind=kind_phys), dimension(IMJM)        :: rain_cpl, rainc_cpl, snow_cpl

      ! Stochastic physics / surface perturbations
      !logical  :: do_sppt, ca_global
      real(kind=kind_phys), dimension(IMJM,LM)      :: dtdtr
      real(kind=kind_phys), dimension(IMJM)        :: drain_cpl, dsnow_cpl

      ! Rainfall variables previous time step
      !integer  :: lsm, lsm_ruc, lsm_noahmp
      real(kind=kind_phys), dimension(IMJM)        :: raincprv
      real(kind=kind_phys), dimension(IMJM)        :: rainncprv
      real(kind=kind_phys), dimension(IMJM)        :: iceprv
      real(kind=kind_phys), dimension(IMJM)        :: snowprv
      real(kind=kind_phys), dimension(IMJM)        :: graupelprv
      real(kind=kind_phys), dimension(IMJM)        :: draincprv
      real(kind=kind_phys), dimension(IMJM)        :: drainncprv
      real(kind=kind_phys), dimension(IMJM)        :: diceprv
      real(kind=kind_phys), dimension(IMJM)        :: dsnowprv
      real(kind=kind_phys), dimension(IMJM)        :: dgraupelprv

      real(kind=kind_phys), parameter :: con_p001= 0.001d0
      real(kind=kind_phys), parameter :: con_day = 86400.d0
      real(kind=kind_phys), parameter :: rainmin = 1.0d-13
      ! *DH

      ! interface variables
      !real(kind=kind_phys), dimension(1:IMJM)          :: frland
      !real(kind=kind_phys), dimension(1:IMJM,1:LM)   :: gq0, gq0_ntcw, gq0_ntrw, gq0_ntiw, &
      !                                                                gq0_ntsw, gq0_ntgl, gq0_ntclamt
      real(kind=kind_phys), dimension(1:IMJM,1:LM)   :: vvl

      ! rain/snow/ice/graupel/precip amounts, fraction of frozen precip
      real(kind_phys),       dimension(1:IMJM) :: prcpmp

      !logical  :: hydrostatic, phys_hydrostatic

      !logical  :: lradar
      real(kind=kind_phys), dimension(1:IMJM,1:LM) :: refl_10cm
      !logical  :: reset !, effr_in
      real(kind=kind_phys), dimension(1:IMJM,1:LM) :: rew, rei, rer, res, reg

       !real(kind=kind_phys)  :: du3dt_pbl(IMJM,LM)
       !real(kind=kind_phys)  :: du3dt_orogwd(IMJM,LM)
       !real(kind=kind_phys)  :: du3dt_deepcnv(IMJM,LM)
       !real(kind=kind_phys)  :: du3dt_congwd(IMJM,LM)
       !real(kind=kind_phys)  :: du3dt_rdamp(IMJM,LM)
       !real(kind=kind_phys)  :: du3dt_shalcnv(IMJM,LM)
       !real(kind=kind_phys)  :: du3dt_phys(IMJM,LM)
       !real(kind=kind_phys) :: dv3dt_pbl(IMJM,LM)
       !real(kind=kind_phys)  :: dv3dt_orogwd(IMJM,LM)
       !real(kind=kind_phys)  :: dv3dt_deepcnv(IMJM,LM)
       !real(kind=kind_phys)  :: dv3dt_congwd(IMJM,LM)
       !real(kind=kind_phys)  :: dv3dt_rdamp(IMJM,LM)
       !real(kind=kind_phys)  :: dv3dt_shalcnv(IMJM,LM)
       !real(kind=kind_phys)  :: dv3dt_phys(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_lw(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_sw(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_pbl(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_deepcnv(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_shalcnv(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_mp(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_orogwd(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_rdamp(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_congwd(IMJM,LM)
       !real(kind=kind_phys)  :: dt3dt_phys(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_pbl(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_deepcnv(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_shalcnv(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_mp(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_o3pbl(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_o3prodloss(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_o3mix(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_o3tmp(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_o3column(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_phys(IMJM,LM)
       !real(kind=kind_phys)  :: dq3dt_o3phys(IMJM,LM)
