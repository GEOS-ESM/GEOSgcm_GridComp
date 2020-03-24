module shoc

! Implementation of the Simplified High Order Closure (SHOC) scheme 
! of Bogenschutz and Krueger (2013), J. Adv. Model. Earth Syst, 5, 195-211,
! doi: 10.1002/jame.200118. (further referred to as BK13)
! in a single column form suitable for use in a GCM physics package. 
! Alex Belochitski, heavily based on the code of Peter Bogenschutz.
! S Moorthi - optimization, cleanup, improve and customize for gsm
!
! N Arnold - implemented in GEOS
!

 use MAPL_ConstantsMod, only: ggr    => MAPL_GRAV,   &
                              cp     => MAPL_CP,     &
                              rgas   => MAPL_RGAS,   &
                              rv     => MAPL_RVAP,   &
                              lcond  => MAPL_ALHL,   &
                              lfus   => MAPL_ALHS,   &
                              pi     => MAPL_PI,     &
                              MAPL_H2OMW, MAPL_AIRMW

 use MAPL,              only: MAPL_UNDEF

 use MAPL_SatVaporMod,  only: MAPL_EQsat 

 implicit none

 private

 public run_shoc

 contains

 subroutine run_shoc( nx, ny, nzm, nz, dtn, dm_inv,              &  ! in
                 prsl_inv, phii_inv, phil_inv, u_inv, v_inv,     &  ! in
                 omega_inv, hflx, evap, buoy_mf_inv,             &  ! in
                 tabs_inv, qwv_inv, qi_inv, qc_inv, qpi_inv,     &  ! in 
                 qpl_inv, cld_sgs_inv, wthv_sec_inv, prnum,      &  ! in
                 tke_inv, tkh_inv,                               &  ! inout
                 isotropy_inv, w3_canuto_inv,                    &  ! out
                 tkesbdiss_inv, tkesbbuoy_inv,                   &  ! out
                 tkesbshear_inv,tkesbtrans_inv,                  &  ! out
                 smixt_inv,smixt_oc_inv,smixt_ic_inv,            &  ! out
                 smixt1_inv,smixt2_inv,smixt3_inv,               &  ! out
                 brunt_inv,shear_inv,                            &  ! out
                 LAM, TSCL, VON, CKVAL, CEFAC, CESFAC,           &  ! tuning param in
                 THLTUN, QWTUN, QWTHLTUN, DO_TRANS, DO_CLDLEN,   &
                 USE_MF_PDF, USE_MF_BUOY, BUOY_OPTION )                 


  real, parameter :: lsub = lcond+lfus,         &
                     fac_cond = lcond/cp,       &
                     fac_fus = lfus/cp,         &
                     fac_sub = lsub/cp,         &
                     ggri = 1.0/ggr,            &
                     kapa = rgas/cp,            &
                     gocp = ggr/cp,             &
                     rog = rgas*ggri,           &
                     sqrt2 = sqrt(2.0),         &
                     sqrtpii = 1.0/sqrt(pi+pi), &
                     epsterm = rgas/rv,         &
                     twoby3 = 2.0/3.0,          &
                     onebeps = 1.0/epsterm,     &
                     twoby15 = 2.0 / 15.0,      &
                     onebrvcp= 1.0/(rv*cp),     &
                     tkef1=0.5,                 &
                     tkef2=1.0-tkef1,           &
                     tkhmax=1000.0,             & 
                     epsv=MAPL_H2OMW/MAPL_AIRMW

  integer, intent(in) :: nx    ! Number of points in the physics window in the x
  integer, intent(in) :: ny    ! and y directions

  integer, intent(in) :: nzm   ! Number of vertical layers
  integer, intent(in) :: nz    ! Number of layer interfaces  (= nzm + 1)   
  integer, intent(in) :: DO_TRANS     ! TKE transport term on/off
  integer, intent(in) :: DO_CLDLEN    ! TKE transport term on/off
  integer, intent(in) :: USE_MF_PDF   ! Use mass flux contrib to PDF
  integer, intent(in) :: USE_MF_BUOY  ! Explicitly add mass flux contribution to buoyancy
  integer, intent(in) :: BUOY_OPTION  ! choose source of TKE buoyancy term
                                      ! 0=local stability
                                      ! 1=single-gaussian PDF,
                                      ! 2=double-gaussian PDF (MoistGridComp)
  real, intent(in   ) :: dtn          ! Physics time step, s 
  real, intent(in   ) :: hflx(nx,ny)  ! Surface sensible heat flux
  real, intent(in   ) :: evap(nx,ny)  ! Surface evaporation

  real, intent(in   ) :: dm_inv   (nx,ny,nzm) ! mean layer mass   
  real, intent(in   ) :: prsl_inv (nx,ny,nzm) ! mean layer presure   
  real, intent(in   ) :: phii_inv (nx,ny,nz ) ! interface geopotential height
  real, intent(in   ) :: phil_inv (nx,ny,nzm) ! layer geopotential height  
  real, intent(in   ) :: u_inv    (nx,ny,nzm) ! u-wind, m/s
  real, intent(in   ) :: v_inv    (nx,ny,nzm) ! v-wind, m/s
  real, intent(in   ) :: omega_inv(nx,ny,nzm) ! omega, Pa/s
  real, intent(in   ) :: buoy_mf_inv(nx,ny,nzm) ! buoyancy flux MF, K*m/s
  real, intent(in   ) :: wthv_sec_inv(nx,ny,nzm) ! Buoyancy flux, K*m/s

  real, intent(in   ) :: tabs_inv   (nx,ny,nzm) ! temperature, K
  real, intent(in   ) :: qwv_inv    (nx,ny,nzm) ! water vapor mixing ratio, kg/kg
  real, intent(in   ) :: qc_inv     (nx,ny,nzm) ! cloud water mixing ratio, kg/kg
  real, intent(in   ) :: qi_inv     (nx,ny,nzm) ! cloud ice   mixing ratio, kg/kg
  real, intent(in   ) :: qpl_inv    (nx,ny,nzm) ! rain mixing ratio, kg/kg
  real, intent(in   ) :: qpi_inv    (nx,ny,nzm) ! snow mixing ratio, kg/kg
  real, intent(in   ) :: cld_sgs_inv(nx,ny,nzm) ! sgs cloud fraction
  real, intent(inout) :: tke_inv    (nx,ny,nzm) ! turbulent kinetic energy. m**2/s**2
  real, intent(inout) :: tkh_inv    (nx,ny,nzm) ! eddy diffusivity
  real, intent(inout) :: prnum      (nx,ny,nzm) ! turbulent Prandtl number

  real, intent(  out) :: isotropy_inv(nx,ny,nzm) ! return to isotropy timescale
  real, intent(  out) :: w3_canuto_inv(nx,ny,nz) ! canuto w3 estimate

  real, dimension(:,:,:), pointer :: tkesbdiss_inv  ! dissipation
  real, dimension(:,:,:), pointer :: tkesbbuoy_inv  ! buoyancy production 
  real, dimension(:,:,:), pointer :: tkesbshear_inv ! shear production 
  real, dimension(:,:,:), pointer :: tkesbtrans_inv ! tke transport 

  real, dimension(:,:,:), pointer :: smixt_inv    ! dissipation length scale
  real, dimension(:,:,:), pointer :: smixt1_inv   ! length scale, term 1
  real, dimension(:,:,:), pointer :: smixt2_inv   ! length scale, term 2
  real, dimension(:,:,:), pointer :: smixt3_inv   ! length scale, term 3
  real, dimension(:,:,:), pointer :: smixt_oc_inv ! dissipation length scale
  real, dimension(:,:,:), pointer :: smixt_ic_inv ! dissipation length scale
  real, dimension(:,:,:), pointer :: brunt_inv    ! Brunt vaisala frequency
  real, dimension(:,:,:), pointer :: shear_inv    ! squared shear diagnostic

  real, intent(in   ) :: LAM, TSCL, VON, CKVAL, CEFAC,   &  ! tuning param
                         CESFAC, THLTUN, QWTUN, QWTHLTUN                 

! SHOC tunable parameters

  real :: lambda
  real, parameter :: min_tke = 1e-6  ! Minumum TKE value, m**2/s**2 
  real, parameter :: max_tke = 10.    ! Maximum TKE value, m**2/s**2 

! Maximum turbulent eddy length scale, m
  real, parameter :: max_eddy_length_scale  = 2000. 

! Maximum "return-to-isotropy" time scale, s
  real, parameter :: max_eddy_dissipation_time_scale = 2000.  

! Constants for the TKE dissipation term based on Deardorff (1980)
  real, parameter :: pt19=0.19, pt51=0.51
  real, parameter :: Cs  = 0.15
  real :: Ck, Ce, Ces

! Use double gaussian PDF-based buoyancy
!  integer, parameter :: use_DG = 1

! real, parameter :: Ce  = Ck**3/(0.7*Cs**4) 
! real, parameter :: Ces = Ce/0.7*3.0

! real, parameter :: vonk=0.35  ! Von Karman constant
  real :: vonk, tscale
  real, parameter :: w_tol_sqd = 4.0e-04   ! Min vlaue of second moment of w
  real, parameter :: w_thresh  = 0.0, thresh = 0.0


! These parameters are a tie-in with a microphysical scheme
! Double check their values for the Zhao-Carr scheme.
  real, parameter :: tbgmin = 233.16    ! Minimum temperature for cloud water., K
  real, parameter :: tbgmax = 273.16    ! Maximum temperature for cloud ice, K
  real, parameter :: a_bg   = 1.0/(tbgmax-tbgmin)

! Parameters to tune the second order moments-  No tuning is performed currently
  real :: thl2tune, qw2tune, qwthl2tune
  real, parameter :: thl_tol  = 1.e-2, rt_tol = 1.e-4

! Number iterations for TKE solution
  integer, parameter :: nitr=6

! Local variables. Note that pressure is in millibars in the SHOC code.

  real zl      (nx,ny,nzm)  ! height of the pressure levels above surface, m
  real zi      (nx,ny,nz)   ! height of the interface levels, m
  real adzl    (nx,ny,nzm)  ! layer thickness i.e. zi(k+1)-zi(k) - defined at levels
  real adzi    (nx,ny,nz)   ! level thickness i.e. zl(k)-zl(k-1) - defined at interface
 
  real hl      (nx,ny,nzm)  ! liquid/ice water static energy , K
  real hvl     (nx,ny,nzm)  ! liquid/ice water static energy , K
  real qv      (nx,ny,nzm)  ! water vapor, kg/kg
  real qcl     (nx,ny,nzm)  ! liquid water  (condensate), kg/kg
  real qci     (nx,ny,nzm)  ! ice water  (condensate), kg/kg
  real w       (nx,ny,nzm)  ! z-wind, m/s
  real bet     (nx,ny,nzm)  ! ggr/tv0
  real gamaz   (nx,ny,nzm)  ! ggr/cp*z
  real dm      (nx,ny,nzm)  
  real prsl    (nx,ny,nzm)  
  real u       (nx,ny,nzm)  
  real v       (nx,ny,nzm)  
  real omega   (nx,ny,nzm)  
  real tabs    (nx,ny,nzm)
  real qwv     (nx,ny,nzm)
  real qpl     (nx,ny,nzm)
  real qpi     (nx,ny,nzm)
  real cld_sgs (nx,ny,nzm)
  real tke     (nx,ny,nzm)
  real tkh     (nx,ny,nzm)
  real wthv_sec(nx,ny,nzm)
  real buoy_mf (nx,ny,nzm)
  real tkesbdiss(nx,ny,nzm)
  real tkesbbuoy(nx,ny,nzm)
  real tkesbshear(nx,ny,nzm)
  real tkesbtrans(nx,ny,nzm)

! Moments of the trivariate double Gaussian PDF for the SGS total water mixing ratio
! SGS liquid/ice static energy, and vertical velocity

  real qw_sec   (nx,ny,nzm) ! Second moment total water mixing ratio, kg^2/kg^2
  real thl_sec  (nx,ny,nzm) ! Second moment liquid/ice static energy, K^2
  real qwthl_sec(nx,ny,nzm) ! Covariance tot. wat. mix. ratio and static energy, K*kg/kg
  real wqw_sec  (nx,ny,nzm) ! Turbulent flux of tot. wat. mix., kg/kg*m/s
  real wthl_sec (nx,ny,nzm) ! Turbulent flux of liquid/ice static energy, K*m/s
  real w_sec    (nx,ny,nzm) ! Second moment of vertical velocity, m**2/s**2
  real w3       (nx,ny,nzm) ! Third moment of vertical velocity, m**3/s**3
  real wqp_sec  (nx,ny,nzm) ! Turbulent flux of precipitation, kg/kg*m/s

! Eddy length formulation 
  real smixt    (nx,ny,nzm) ! Turbulent length scale, m
  real smixt1    (nx,ny,nzm) ! Turbulent length scale, m
  real smixt2    (nx,ny,nzm) ! Turbulent length scale, m
  real smixt3    (nx,ny,nzm) ! Turbulent length scale, m
  real smixt_incld(nx,ny,nzm) ! Turbulent length scale, m
  real smixt_outcld(nx,ny,nzm) ! Turbulent length scale, m
  real isotropy (nx,ny,nzm) ! "Return-to-isotropy" eddy dissipation time scale, s
! real isotropy_debug (nx,ny,nzm) ! Return to isotropy scale, s without artificial limits
  real brunt    (nx,ny,nzm) ! Moist Brunt-Vaisalla frequency, s^-1
  real conv_vel2(nx,ny,nzm) ! Convective velocity scale cubed, m^3/s^3


! Local variables

  real, dimension(nx,ny,nzm) :: total_water, brunt2, def2, thv

  real, dimension(nx,ny)     :: denom, numer, l_inf, cldarr

  real lstarn,    depth,    omn,      betdz,    bbb,        &
       term,      qsatt,    dqsat,    thedz,    conv_var,   &  
       tkes,      pval,     pkap,     thlsec,   qwsec,      &
       qwthlsec,  wqwsec,   wthlsec,  dum,      sm,         &
       prespot,   wrk,      wrk1,     wrk2,     wrk3,       & 
       tkeavg,    dtqw,     dtqi

  integer i,j,k,km1,ku,kd,ka,kb,kinv,strt,fnsh,cnvl

! set parameter values
  lambda  = LAM  ! used in return-to-isotropy timescale
  Ck  = CKVAL     ! Coeff in the eddy diffusivity - TKE relationship, see Eq. 7 in BK13 
  Ce  = CEFAC*Ck**3/Cs**4   ! diss ~= Ce * sqrt(tke)
  Ces = CESFAC*Ce
  vonk = VON     ! Von Karman constant Moorthi - as in GFS
  tscale = TSCL  ! time scale set based off of similarity results of BK13, s
  thl2tune = THLTUN
  qw2tune = QWTUN
  qwthl2tune = QWTHLTUN


  wthl_sec = 0.

! Map GEOS variables to those of SHOC

  do k=1,nz
    do j=1,ny
      do i=1,nx
        kinv = nz-k+1
        zi(i,j,kinv)    = phii_inv(i,j,k)-phii_inv(i,j,nz)
      enddo
    enddo
  enddo

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        kinv = nzm-k+1
        zl(i,j,kinv)       = phil_inv(i,j,k)-phii_inv(i,j,nz)
        tkh(i,j,kinv)      = tkh_inv(i,j,k)
        dm(i,j,kinv)       = dm_inv(i,j,k)
        prsl(i,j,kinv)     = prsl_inv(i,j,k)
        u(i,j,kinv)        = u_inv(i,j,k)
        v(i,j,kinv)        = v_inv(i,j,k)
        omega(i,j,kinv)    = omega_inv(i,j,k)
        tabs(i,j,kinv)     = tabs_inv(i,j,k)
        qwv(i,j,kinv)      = qwv_inv(i,j,k)
        qcl(i,j,kinv)       = qc_inv(i,j,k)
        qci(i,j,kinv)       = qi_inv(i,j,k)
        cld_sgs(i,j,kinv)  = cld_sgs_inv(i,j,k)
        tke(i,j,kinv)      = tke_inv(i,j,k)
        wthv_sec(i,j,kinv) = wthv_sec_inv(i,j,k)
        buoy_mf(i,j,kinv)  = USE_MF_BUOY*buoy_mf_inv(i,j,k)
      enddo
    enddo
  enddo
  if (BUOY_OPTION==2 .and. USE_MF_PDF==1) buoy_mf = 0.0

           
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        wrk          = 1.0 / prsl(i,j,k)
        qv(i,j,k)    = max(qwv(i,j,k), 0.0)
        thv(i,j,k)   = tabs(i,j,k) * (1.0+epsv*qv(i,j,k))
        w(i,j,k)     = - rog * omega(i,j,k) * thv(i,j,k) * wrk
        qpl(i,j,k)     = 0.0  ! comment or remove when using with prognostic rain/snow
        qpi(i,j,k)     = 0.0  ! comment or remove when using with prognostic rain/snow
        wqp_sec(i,j,k) = 0.0  ! Turbulent flux of precipiation
        total_water(i,j,k) = qcl(i,j,k) + qci(i,j,k) + qv(i,j,k)

        prespot        = (100000.0*wrk) ** kapa ! Exner function
        bet(i,j,k)     = ggr/(tabs(i,j,k)*prespot)     ! Moorthi
        thv(i,j,k)     = thv(i,j,k)*prespot            ! Moorthi
!
! Lapse rate * height = reference temperature
        gamaz(i,j,k) = gocp * zl(i,j,k)

! Liquid/ice water static energy - ! Note the the units are degrees K
        hl(i,j,k) = tabs(i,j,k) + gamaz(i,j,k) - fac_cond*(qcl(i,j,k)+qpl(i,j,k)) &
                                               - fac_fus *(qci(i,j,k)+qpi(i,j,k))
        hvl(i,j,k) = tabs(i,j,k)*(1.0+epsv*qv(i,j,k)-qcl(i,j,k))+gamaz(i,j,k) &
                     - fac_cond*qcl(i,j,k) - fac_fus *qci(i,j,k)
        w3(i,j,k) = 0.0
      enddo
    enddo
  enddo

   
! Define vertical grid increments for later use in the vertical differentiation

  do k=2,nzm
    km1 = k - 1
    do j=1,ny
      do i=1,nx
        adzi(i,j,k)   = (zl(i,j,k) - zl(i,j,km1))
        adzl(i,j,km1) = (zi(i,j,k) - zi(i,j,km1))   ! from 1 to nzm-1
      enddo
    enddo
  enddo
  do j=1,ny
    do i=1,nx
      adzi(i,j,1)   = (zl(i,j,1)-zi(i,j,1)) !  unused in the code
      adzi(i,j,nz)  = zi(i,j,nz)-zl(i,j,nzm) ! at the top - probably unused
      adzl(i,j,nzm) = adzi(i,j,nzm)
!
      wthl_sec(i,j,1) = hflx(i,j)/cp
      wqw_sec(i,j,1)  = evap(i,j)
    enddo
  enddo

  call tke_shoc()        ! Integrate prognostic TKE equation forward in time

!  print *,'new tke(1,1,:)=',tke(1,1,:)

  
! diagnose second order moments of the subgrid PDF following
! Redelsperger J.L., and G. Sommeria, 1986, JAS, 43, 2619-2635 sans the use of stabilty
! weighting functions - Result is in global variables w_sec, thl_sec, qw_sec, and qwthl_sec
    

! Second moment of vertical velocity.
! Note that Eq 6 in BK13 gives a different expression that is dependent on 
! vertical gradient of grid scale vertical velocity 

    do k=1,nzm
      ku = k+1
      kd = k-1
      ka = ku
      kb = k
      if (k == 1) then
        kd = k
        kb = ka
      elseif (k == nzm) then
        ku = k
        ka = kb
      endif
      do j=1,ny
        do i=1,nx
          if (tke(i,j,k) > 0.0) then
! JDC Modify
!            wrk  = 0.5*(tkh(i,j,ka)+tkh(i,j,kb))*(w(i,j,ku) - w(i,j,kd)) &
!                 / (sqrt(tke(i,j,k)) * (zl(i,j,ku) - zl(i,j,kd)))
            wrk = 0.0
            w_sec(i,j,k) = max(twoby3 * tke(i,j,k) - twoby15 * wrk, 0.0)
          else
            w_sec(i,j,k) = 0.0
          endif
        enddo
      enddo
    enddo
    
    do k=2,nzm
       
      km1 = k-1
      do j=1,ny
        do i=1,nx

! Use backward difference in the vertical, use averaged values of "return-to-isotropy"
! time scale and diffusion coefficient

          wrk1 = 1.0 / adzi(i,j,k)        ! adzi(k) = (zl(k)-zl(km1))
          wrk3 = tkh(i,j,k) * wrk1

          sm   = 0.5*(isotropy(i,j,k)+isotropy(i,j,km1))*wrk1*wrk3 ! Tau*Kh/dz^2
             
! SGS vertical flux liquid/ice water static energy. Eq 1 in BK13
             
          wrk1            = hl(i,j,k) - hl(i,j,km1)
          wthl_sec(i,j,k) = - wrk3 * wrk1

! SGS vertical flux of total water. Eq 2 in BK13

          wrk2           = total_water(i,j,k) - total_water(i,j,km1)
          wqw_sec(i,j,k) = - wrk3 * wrk2

! Second moment of liquid/ice water static energy. Eq 4 in BK13

          thl_sec(i,j,k) = thl2tune * sm * wrk1 * wrk1

! Second moment of total water mixing ratio.  Eq 3 in BK13
             
          qw_sec(i,j,k) = qw2tune * sm * wrk2 * wrk2
             
! Covariance of total water mixing ratio and liquid/ice water static energy.
! Eq 5 in BK13
             
          qwthl_sec(i,j,k) = qwthl2tune * sm * wrk1 * wrk2

        enddo ! i  loop
      enddo   ! j  loop
    enddo     ! k  loop

    do j=1,ny
      do i=1,nx
!       wthl_sec(i,j,1)  = wthl_sec(i,j,2)  ! these are set to surface fluxes above
!       wqw_sec(i,j,1)   = wqw_sec(i,j,2)
        thl_sec(i,j,1)   = thl_sec(i,j,2)
        qw_sec(i,j,1)    = qw_sec(i,j,2)
        qwthl_sec(i,j,1) = qwthl_sec(i,j,2)
      enddo
    enddo


! Diagnose the third moment of SGS vertical velocity

  call canuto()


! If using single-gaussian PDF for buoyancy flux,
! this will overwrite wthv_sec

if (BUOY_OPTION==1) then
  call buoyancy_single_gaussian()
end if


!=== Assign exports and flip vertical ===!

  tkh_inv(:,:,1:nzm)       = tkh(:,:,nzm:1:-1)
  isotropy_inv(:,:,1:nzm)  = isotropy(:,:,nzm:1:-1)
  tke_inv(:,:,1:nzm)       = tke(:,:,nzm:1:-1)

  w3_canuto_inv(:,:,1:nzm) = w3(:,:,nzm:1:-1)

  if (associated(tkesbdiss_inv))  tkesbdiss_inv(:,:,1:nzm)  = tkesbdiss(:,:,nzm:1:-1)
  if (associated(tkesbbuoy_inv))  tkesbbuoy_inv(:,:,1:nzm)  = tkesbbuoy(:,:,nzm:1:-1)
  if (associated(tkesbshear_inv)) tkesbshear_inv(:,:,1:nzm) = tkesbshear(:,:,nzm:1:-1)
  if (associated(tkesbtrans_inv)) tkesbtrans_inv(:,:,1:nzm) = tkesbtrans(:,:,nzm:1:-1)

  if (associated(smixt_inv))    smixt_inv(:,:,1:nzm)    = smixt(:,:,nzm:1:-1)
  if (associated(smixt1_inv))   smixt1_inv(:,:,1:nzm)   = smixt1(:,:,nzm:1:-1)
  if (associated(smixt2_inv))   smixt2_inv(:,:,1:nzm)   = smixt2(:,:,nzm:1:-1)
  if (associated(smixt3_inv))   smixt3_inv(:,:,1:nzm)   = smixt3(:,:,nzm:1:-1)
  if (associated(smixt_oc_inv)) smixt_oc_inv(:,:,1:nzm) = smixt_outcld(:,:,nzm:1:-1)
  if (associated(smixt_ic_inv)) smixt_ic_inv(:,:,1:nzm) = smixt_incld(:,:,nzm:1:-1)
  
  if (associated(brunt_inv))    brunt_inv(:,:,1:nzm)    = brunt2(:,:,nzm:1:-1)
  if (associated(shear_inv))    shear_inv(:,:,1:nzm)    = def2(:,:,nzm:1:-1)

!========================================!


contains

  subroutine tke_shoc()

! This subroutine solves the TKE equation, 
! Heavily based on SAM's tke_full.f90 by Marat Khairoutdinov

    real grd,betdz,Cek,Cee,lstarn, lstarp, bbb, omn, omp,qsatt,dqsat, smix,         &
         buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss,a_prod_bu_debug, buoy_sgs_debug, &
         tscale1, wrk, wrk1, wtke, wtk2, rdtn
    integer i,j,k,ku,kd,itr

    rdtn = 1.0 / dtn

    call tke_shear_prod(def2)   ! Calculate shear production of TKE

!    print *,'def2=',sqrt(def2(20,20,:))

! Ensure values of TKE are reasonable

    do k=1,nzm
      do j=1,ny
        do i=1,nx
          tke(i,j,k)        = max(min_tke,tke(i,j,k))
          tkesbdiss(i,j,k)  = 0.
          tkesbshear(i,j,k) = 0.
          tkesbbuoy(i,j,k)  = 0.
          tkesbtrans(i,j,k) = 0.
        enddo
      enddo
    enddo

    call eddy_length()   ! Find turbulent mixing length
    call check_eddy()    ! Make sure it's reasonable

!---------------------------------------------------------------
! UWMT-style estimate of TKE transport within convective layers
!---------------------------------------------------------------

    if (DO_TRANS/=0) then
    cnvl = 0
    do i=1,nx
     do j=1,ny
      fnsh=0
      do k=1,nzm
        if (wthv_sec(i,j,k)+buoy_mf(i,j,k) .ge.0.01 .and. cnvl.eq.0) then ! if bot of conv layer
!        if (mf(i,j,k).ge.0.001 .and. cnvl.eq.0) then ! if bot of conv layer
          cnvl = 1
          strt = k
        end if        

        ! if top of conv layer, relax tke within conv layer to
        ! pressure-weighted mean tke
        if (cnvl.eq.1 .and. wthv_sec(i,j,k)+buoy_mf(i,j,k).lt.0.01) then
!        if (cnvl.eq.1 .and. mf(i,j,k).lt.0.001) then
          fnsh = k-1
          if (strt<fnsh) then
             tkeavg = sum(tke(i,j,strt:fnsh)*dm(i,j,strt:fnsh))/sum(dm(i,j,strt:fnsh))
             tscale1 = sqrt(tkeavg)*(fnsh-strt+1.)/max(sum(smixt(i,j,strt:fnsh)),10.)
             tkesbtrans(i,j,strt:fnsh) = (tkeavg-tke(i,j,strt:fnsh))*tscale1
          end if

          ! check above unstable layer for tkh>0
          do while(tkh(i,j,fnsh+1).gt.0.05*maxval(tkh(i,j,strt:fnsh)) .and. fnsh.lt.nzm-2)
            fnsh=fnsh+1
          end do
          ! fnsh is index of last interface with tkh>5% of max
          ! but is index of layer above that interface
          fnsh = fnsh-1

          cnvl = 0
        end if ! if top of conv layer
      end do ! k loop
     end do ! j loop
    end do ! i loop
    end if ! do_transport


    do k=1,nzm      
      ku = k+1
      kd = k
      
     Cek = Ce/0.7

      if(k == 1) then
        ku = 2
        kd = 2
        Cek = Ces
      elseif(k == nzm) then
        ku = k
        kd = k
        Cek = Ces
      end if

      
      do j=1,ny
        do i=1,nx
          grd = adzl(i,j,k)             !  adzl(k) = zi(k+1)-zi(k)

!         wrk = zl(i,j,k) / grd + 1.5
!         cek = 1.0 + 2.0 / (wrk*wrk -3.3)

! TKE boyancy production term. wthv_sec (buoyancy flux) is calculated in
! Moist GridComp. The value used here is from the previous time step

!         a_prod_bu = bet(i,j,k)*wthv_sec(i,j,k)
          a_prod_bu = (ggr / thv(i,j,k)) * (wthv_sec(i,j,k)+buoy_mf(i,j,k))

! If wthv_sec from subgrid PDF is not available use Brunt-Vaisalla frequency from eddy_length()
!         wrk  = (0.5*ck)  * (tkh(i,j,ku)+tkh(i,j,kd))
          wrk  = 0.5 * (tkh(i,j,ku)+tkh(i,j,kd))

!Obtain Brunt-Vaisalla frequency from diagnosed SGS buoyancy flux
!Presumably it is more precise than BV freq. calculated in  eddy_length()?
          
          buoy_sgs = brunt(i,j,k)
!          buoy_sgs = - a_prod_bu / (wrk + 0.0001)   ! tkh is eddy thermal diffussivity
!          buoy_sgs = - a_prod_bu / (prnum*wrk + 0.0001)   ! tk is eddy viscosity

!Compute $c_k$ (variable Cee) for the TKE dissipation term following Deardorff (1980)

          if (buoy_sgs <= 0.0) then
            smix = grd
          else
!           smix = min(grd,max(0.1*grd, sqrt(0.76*wrk/sqrt(buoy_sgs+1.e-10))))
!           smix = min(grd,max(0.1*grd, sqrt(0.76*wrk/(Ck*sqrt(buoy_sgs+1.e-10)))))
            smix = min(grd,max(0.1*grd, 0.76*sqrt(tke(i,j,k)/(buoy_sgs+1.e-10))))
          end if

          ratio     = smix/grd

          Cee = Cek* (pt19 + pt51*ratio)

          wrk   = 0.5 * wrk * (prnum(i,j,ku) + prnum(i,j,kd))
          a_prod_sh = min(min(tkhmax,(wrk+0.0001))*def2(i,j,k),0.001)    ! TKE shear production term

! smixt (turb. mixing lenght) is calculated in eddy_length() 
! Explicitly integrate TKE equation forward in time
!         a_diss     = Cee/smixt(i,j,k)*tke(i,j,k)**1.5 ! TKE dissipation term
!         tke(i,j,k) = max(0.,tke(i,j,k)+dtn*(max(0.,a_prod_sh+a_prod_bu)-a_diss))

! Semi-implicitly integrate TKE equation forward in time

          wtke = tke(i,j,k)
          wtk2 = wtke
          wrk  = (dtn*Cee)/smixt(i,j,k)
          wrk1 = wtke + dtn*(a_prod_sh+a_prod_bu+tkesbtrans(i,j,k))

          do itr=1,nitr                        ! iterate for implicit solution
            wtke   = min(max(min_tke, wtke), max_tke)
            a_diss = wrk*sqrt(wtke)            ! Coefficient in the TKE dissipation term
            if (a_diss.ne.-1.) then
              wtke   = wrk1 / (1.+a_diss)
            else
              wtke   = wrk1 / (1.01+a_diss)
            end if
            wtke   = tkef1*wtke + tkef2*wtk2   ! tkef1+tkef2 = 1.0
            wtk2   = wtke

          enddo

          tke(i,j,k) = min(max(min_tke, wtke), max_tke)

          tscale1    = (dtn+dtn) / a_diss        ! See Eq 8 in BK13

          a_diss     = (a_diss/dtn)*tke(i,j,k)    ! TKE dissipation term, epsilon


! Calculate "return-to-isotropy" eddy dissipation time scale, see Eq. 8 in BK13

          if (buoy_sgs <= 0.0) then
            isotropy(i,j,k) = min(max_eddy_dissipation_time_scale,tscale1)
          else
            isotropy(i,j,k) = min(max_eddy_dissipation_time_scale,          &
                             tscale1/(1.0+lambda*buoy_sgs*tscale1*tscale1))
          endif


! TKE budget terms

          tkesbdiss(i,j,k)       = -a_diss
          tkesbshear(i,j,k)      = a_prod_sh
          tkesbbuoy(i,j,k)       = a_prod_bu
!          tkesbtrans(i,j,k)      = a_prod_tr
!         tkesbbuoy_debug(i,j,k) = a_prod_bu_debug
!         tkebuoy_sgs(i,j,k)     = buoy_sgs

        end do ! i loop
      end do   ! j loop
    end do     ! k
!
    wrk = 0.5 * ck
    do k=2,nzm
      do j=1,ny
        do i=1,nx

          wrk1 = wrk / (prnum(i,j,k) + prnum(i,j,k-1))

          tkh(i,j,k) = wrk1 * (isotropy(i,j,k) + isotropy(i,j,k-1))     &
                            * (tke(i,j,k)      + tke(i,j,k-1)) ! Eddy thermal diffusivity
          tkh(i,j,k) = min(tkh(i,j,k),tkhmax)
        end do ! i
      end do ! j
    end do ! k



  end subroutine tke_shoc

 
  subroutine tke_shear_prod(def2)

! Calculate TKE shear production term 

    real, intent(out):: def2(nx,ny,nzm)

    real    rdzw_up, rdzw_dn, wrku(2), wrkv(2), wrkw(2)
    real    txd(nx,ny)
    integer i,j,k,kb,kc
    
   do k=1,nzm
     do j=1,ny
       do i=1,nx
         def2(i,j,k) = 0.0
       enddo
     enddo
   enddo

! Calculate TKE shear production term 

!    print *,'adzl=',adzl(1,1,:)
!    print *,'adzi=',adzi(1,1,:)

    do k=1,nzm  
       
      kb = k-1
      kc = k+1
    
      if (k == 1) then           
          
        do j=1,ny
          do i=1,nx
            rdzw_up     = 1./adzi(i,j,kc)
            wrku(1)     = (u(i,j,kc)-u(i,j,k))*rdzw_up
            wrkv(1)     = (v(i,j,kc)-v(i,j,k))*rdzw_up
!           wrkw(1)     = (w(i,j,kc)-w(i,j,k))*rdzw_up
            def2(i,j,1) = wrku(1)*wrku(1) + wrkv(1)*wrkv(1) !+ 2*wrkw(1) * wrkw(1)
            txd(i,j)    = rdzw_up
          enddo
        enddo
                       
      elseif (k < nzm ) then
        do j=1,ny
          do i=1,nx
!            print *,'adzi(k)=',adzi(i,j,kc)
            rdzw_up     = 1./adzi(i,j,kc)
            rdzw_dn     = txd(i,j)
            wrku(1)     = (u(i,j,kc)-u(i,j,k))*rdzw_up
            wrku(2)     = (u(i,j,k)-u(i,j,kb))*rdzw_dn
            wrkv(1)     = (v(i,j,kc)-v(i,j,k))*rdzw_up
            wrkv(2)     = (v(i,j,k)-v(i,j,kb))*rdzw_dn
!           wrkw(1)     = (w(i,j,kc)-w(i,j,k))*rdzw_up
!           wrkw(2)     = (w(i,j,k)-w(i,j,kb))*rdzw_dn

            def2(i,j,k) = 0.5 * (wrku(1)*wrku(1) + wrku(2)*wrku(2)     &
                               + wrkv(1)*wrkv(1) + wrkv(2)*wrkv(2))  ! &
!                              + wrkw(1)*wrkw(1) + wrkw(2)*wrkw(2)
            txd(i,j)    = rdzw_up
          enddo
        enddo
      else
        do j=1,ny
          do i=1,nx
            rdzw_dn     = txd(i,j)
            wrku(2)     = (u(i,j,k)-u(i,j,kb))*rdzw_dn
            wrkv(2)     = (v(i,j,k)-v(i,j,kb))*rdzw_dn
!           wrkw(2)     = (w(i,j,k)-w(i,j,kb))*rdzw_dn
            def2(i,j,k) = wrku(2)*wrku(2) + wrkv(2)*wrkv(2)  !+ 2*wrkw(2) * wrkw(2)
          enddo
        enddo
      endif

    end do     ! k  loop


  end subroutine tke_shear_prod

  subroutine eddy_length()

! This subroutine computes the turbulent length scale based on a new
! formulation described in BK13

! Local variables
    real    wrk, wrk1, wrk2, wrk3
    integer i, j, k, kk, kl, ku, kb, kc, kli, kui
    
    do j=1,ny
      do i=1,nx
        cldarr(i,j) = 0.0
        numer(i,j)  = 0.0
        denom(i,j)  = 0.0
      enddo
    enddo
    
! Find the length scale outside of clouds, that includes boundary layers.
    
    do k=1,nzm
      do j=1,ny
        do i=1,nx
             
! Reinitialize the mixing length related arrays to zero
          smixt(i,j,k)    = 1.0   ! shoc_mod module variable smixt
          brunt(i,j,k)    = 0.0

!Eq. 11 in BK13 (Eq. 4.13 in Pete's dissertation)
!Outside of cloud, integrate from the surface to the cloud base
!Should the 'if' below check if the cloud liquid < a small constant instead?

          if (qcl(i,j,k)+qci(i,j,k) <= 1e-6 .and. cldarr(i,j).eq.0.0) then 
!          if (qcl(i,j,k)+qci(i,j,k) <= 0.) then 
            tkes       = sqrt(tke(i,j,k)) * adzl(i,j,k)
            numer(i,j) = numer(i,j) + tkes*zl(i,j,k) ! Numerator in Eq. 11 in BK13
            denom(i,j) = denom(i,j) + tkes           ! Denominator in Eq. 11 in BK13
          else
            cldarr(i,j) = 1.0   ! Take note of columns containing cloud.
          endif
        enddo
      enddo
    enddo

! Calculate the measure of PBL depth,  Eq. 11 in BK13 (Is this really PBL depth?)
    do j=1,ny
      do i=1,nx
        if (denom(i,j) >  0.0 .and. numer(i,j) > 0.0) then
          l_inf(i,j) = max(min(0.1 * (numer(i,j)/denom(i,j)),300.),10.)
        else
          l_inf(i,j) = 100.
        endif
      enddo
    enddo
    
    
!Calculate length scale outside of cloud, Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)
    do k=1,nzm

      kb = k-1
      kc = k+1
       
      do j=1,ny
        do i=1,nx

!  vars module variable bet (=ggr/tv0) ; grid module variable  adzi
       
          if (k == 1) then
            kb = 1
            kc = 2
            thedz = adzi(i,j,kc)
          elseif (k == nzm) then
            kb = nzm-1
            kc = nzm
            thedz = adzi(i,j,k)
          else
            thedz = (adzi(i,j,kc)+adzi(i,j,k)) !  = (z(k+1)-z(k-1))
          endif
          betdz = bet(i,j,k) / thedz
           
             
          tkes = sqrt(tke(i,j,k))
             
! Compute local Brunt-Vaisalla frequency
             
          wrk = qcl(i,j,k) + qci(i,j,k)
          if (wrk > 0.0) then            ! If in the cloud
            
! Find the in-cloud Brunt-Vaisalla frequency
                
             omn = qcl(i,j,k) / (wrk+1.e-20) ! Ratio of liquid water to total water

! Latent heat of phase transformation based on relative water phase content
! fac_cond = lcond/cp, fac_fus = lfus/cp

             lstarn = fac_cond + (1.-omn)*fac_fus

! Saturation mixing ratio over water/ice wrt temp  based on relative water phase content

!             qsatt =     omn  * qsatw(tabs(i,j,k),prsl(i,j,k))                &
!                   + (1.-omn) * qsati(tabs(i,j,k),prsl(i,j,k))
             qsatt =     omn  * MAPL_EQsat(tabs(i,j,k),prsl(i,j,k),dtqw)                &
                   + (1.-omn) * MAPL_EQsat(tabs(i,j,k),prsl(i,j,k),dtqi,OverIce=.TRUE.)

! Derivative of saturation mixing ratio over water/ice wrt temp. based on relative water phase content
!             dqsat =     omn  * dtqsatw(tabs(i,j,k),prsl(i,j,k))             &
!                   + (1.-omn) * dtqsati(tabs(i,j,k),prsl(i,j,k))
              dqsat =  omn * dtqw + (1.-omn) * dtqi

! liquid/ice moist static energy static energy divided by cp?

             bbb = (1. + epsv*qsatt-wrk-qpl(i,j,k)-qpi(i,j,k)                &
                 + 1.61*tabs(i,j,k)*dqsat) / (1.+lstarn*dqsat)

! Calculate Brunt-Vaisalla frequency using centered differences in the vertical

             brunt(i,j,k) = betdz*(bbb*(hl(i,j,kc)-hl(i,j,kb))               &
                          + (bbb*lstarn - (1.+lstarn*dqsat)*tabs(i,j,k))     &
                          * (total_water(i,j,kc)-total_water(i,j,kb))        & 
                          + (bbb*fac_cond - (1.+fac_cond*dqsat)*tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))  &
                          + (bbb*fac_sub  - (1.+fac_sub*dqsat)*tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )
                
          else                       ! outside of cloud
                
! Find outside-of-cloud Brunt-Vaisalla frequency
! Only unsaturated air, rain and snow contribute to virt. pot. temp. 
! liquid/ice moist static energy divided by cp?

             bbb = 1. + epsv*qv(i,j,k) - qpl(i,j,k) - qpi(i,j,k)
             brunt(i,j,k) = betdz*( bbb*(hl(i,j,kc)-hl(i,j,kb))                        &
                          + epsv*tabs(i,j,k)*(total_water(i,j,kc)-total_water(i,j,kb)) &
                          + (bbb*fac_cond-tabs(i,j,k))*(qpl(i,j,kc)-qpl(i,j,kb))       &
                          + (bbb*fac_sub -tabs(i,j,k))*(qpi(i,j,kc)-qpi(i,j,kb)) )

          end if

             
! Reduction of mixing length in the stable regions (where B.-V. freq. > 0) is required.
! Here we find regions of Brunt-Vaisalla freq. > 0 for later use. 

            if (brunt(i,j,k) >= 1e-5) then
              brunt2(i,j,k) = brunt(i,j,k)
            else
              brunt2(i,j,k) = 1e-5
            endif
             
! Calculate turbulent length scale in the boundary layer.
! See Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)

! Keep the length scale adequately small near the surface following Blackadar (1984)
! Note that this is not documented in BK13 and was added later for SP-CAM runs

!         if (k == 1) then
!           term = 600.*tkes
!           smixt(i,j,k) = term + (0.4*zl(i,j,k)-term)*exp(-zl(i,j,k)*0.01)
!         else

! tscale is the eddy turnover time scale in the boundary layer and is 
! an empirically derived constant 

            if (tkes > 0.0 .and. l_inf(i,j) > 0.0) then
              wrk1 = (tscale*tkes*vonk*zl(i,j,k))
              wrk2 = (tscale*tkes*l_inf(i,j))
              wrk3 = tke(i,j,k) /(0.01 * brunt2(i,j,k))
!              wrk1 = 1.0 / (tscale*vonk*zl(i,j,k))
!              wrk2 = 1.0 / (tscale*l_inf(i,j))
              smixt1(i,j,k) = wrk1
              smixt2(i,j,k) = wrk2
              smixt3(i,j,k) = wrk3
!              smixt3(i,j,k) = sqrt(brunt2(i,j,k)) / (0.7*tkes)
              wrk1 = 1.0 / (1./wrk1 + 1./wrk2 + 1./wrk3)
              smixt(i,j,k) = min(max_eddy_length_scale, 2.8284*sqrt(wrk1)/0.3)
!              smixt(i,j,k) = min(max_eddy_length_scale,        sqrt(wrk1)/0.3)

           endif
           
!         end if  ! not k=1

         smixt_outcld(i,j,k) = smixt(i,j,k)
  
           
!====== Length scale testing ======

         if (zl(i,j,k).gt.5000.) smixt(i,j,k)=min(smixt(i,j,k),50.)


!==================================

        end do
          
      end do
    end do
    
    
! Now find the in-cloud turbulence length scale
! See Eq. 13 in BK13 (Eq. 4.18 in Pete's disseration)
    
!======== Buoyancy flux from local stability ======!
if (BUOY_OPTION==0) then
   wthv_sec = -1.0*(thv/ggr)*brunt*tkh
endif
!==================================================!
    
! determine cubed convective velocity scale (conv_vel2) inside the cloud

!   call conv_scale()         ! inlining the relevant code

!!! npa, commented out
!    do j=1,ny
!      do i=1,nx
!        conv_vel2(i,j,1) = 0. ! Convective velocity scale cubed
!      enddo
!    enddo
                              ! Integrate velocity scale in the vertical
!    do k=2,nzm
!      do j=1,ny
!        do i=1,nx
!          conv_vel2(i,j,k) = conv_vel2(i,j,k-1)                               &
!                           + 2.5*adzi(i,j,k)*bet(i,j,k)*wthv_sec(i,j,k)
!        enddo
!      enddo
!    enddo
!!! npa
    
    do j=1,ny
      do i=1,nx
          
        if (cldarr(i,j) == 1) then ! If there's a cloud in this column 
             
          kl = 0
          ku = 0
          do k=2,nzm-3
                
! Look for the cloud base in this column  
! thresh (=0) is a  variable local to eddy_length(). Should be a module constant.
            wrk = qcl(i,j,k) + qci(i,j,k)
            if (wrk > thresh .and. kl == 0) then
              kl = k
            endif
                
! Look for the cloud top in this column
            if (wrk > thresh .and. qcl(i,j,k+1)+qci(i,j,k+1) <= thresh) then
              ku = k
!!! npa
! conv_vel2 (Cubed convective velocity scale) is calculated in conv_scale()
! Use the value of conv_vel2 at the top of the cloud. 
!              conv_var = conv_vel2(i,j,k)**(1./3.)
!!! npa
            endif
                
! Compute the mixing length scale for the cloud layer that we just found
!!! npa:  no need to exclude clouds 1 layer thick
!            if (kl > 0 .and. ku > 0 .and. ku-kl > 1) then
            if (kl > 0 .and. ku > 0 .and. ku-kl > 0) then

              conv_var=0.
              do kk=kl,ku
                conv_var = conv_var+2.5*adzi(i,j,kk)*bet(i,j,kk)*wthv_sec(i,j,kk)
              end do
              conv_var = conv_var**(1./3.)

!!! npa                
              if (conv_var > 0) then ! If convective vertical velocity scale > 0
                 
                depth = min((zl(i,j,ku)-zl(i,j,kl)) + adzi(i,j,kl), 1000.)
                      
                     
                do kk=kl,ku
! in-cloud turbulence length scale, Eq. 13 in BK13 (Eq. 4.18)

!                  wrk = conv_var/(depth*sqrt(tke(i,j,kk)))
!                  wrk = wrk * wrk + 0.01*brunt2(i,j,kk)/tke(i,j,kk)
! JDC Modify (based on Eq 13 in BK13)
                  wrk = conv_var/(depth*depth*sqrt(tke(i,j,kk)))
                  smixt_incld(i,j,kk) = 1./wrk
                  wrk = wrk + 0.01*brunt2(i,j,kk)/tke(i,j,kk)
                  
! npa modify, weight in-cloud length scale by local cloud fraction
                  if (DO_CLDLEN/=0) then
                    smixt(i,j,kk) = (1.-cld_sgs(i,j,kk))*smixt(i,j,kk) + cld_sgs(i,j,kk)*3.3*sqrt(1.0/wrk)
                    smixt(i,j,kk) = min(max_eddy_length_scale, smixt(i,j,kk))
                  end if

                enddo
                      
              endif ! If convective vertical velocity scale > 0
              kl = 0.
              ku = 0.
            endif ! if inside the cloud layer
                
          enddo   ! k=2,nzm-3
        endif     ! if in the cloudy column
      enddo       ! i=1,nx
    enddo         ! j=1,ny
    
    
  end subroutine eddy_length


  subroutine conv_scale()

! This subroutine calculates the cubed convective velocity scale needed 
! for the definition of the length scale in clouds
! See Eq. 16 in BK13 (Eq. 4.21 in Pete's dissertation)

    integer i, j, k

!!!!!!!!!
!! A bug in formulation of conv_vel
!  Obtain it by averaging conv_vel2 in the horizontal
!!!!!!!!!!

!   conv_vel(1)=0.      ! Horizontally averaged convective velocity scale cubed 
    do j=1,ny
      do i=1,nx
        conv_vel2(i,j,1) = 0. ! Convective velocity scale cubed
      enddo
    enddo
! Integrate velocity scale in the vertical
    do k=2,nzm
!     conv_vel(k)=conv_vel(k-1)
      do j=1,ny
        do i=1,nx
!**********************************************************************
!Do not include grid-scale contribution to convective velocity scale in GCM applications 
!         conv_vel(k)=conv_vel(k-1)+2.5*adzi(k)*bet(k)*(tvwle(k)+tvws(k))
!         conv_vel(k)=conv_vel(k)+2.5*adzi(i,j,k)*bet(i,j,k)*(tvws(k))
!Do not include grid-scale contribution to convective velocity scale in GCM applications 
!         conv_vel2(i,j,k)=conv_vel2(i,j,k-1)+2.5*adzi(k)*bet(k)*(tvwle(k)+wthv_sec(i,j,k))
!**********************************************************************

          conv_vel2(i,j,k) = conv_vel2(i,j,k-1)                               &
                           + 2.5*adzi(i,j,k)*bet(i,j,k)*wthv_sec(i,j,k)
        enddo
      enddo
    enddo

  end subroutine conv_scale


  subroutine check_eddy()

! This subroutine checks eddy length values 

    integer i, j, k, kb, ks, zend
    real    wrk

    do k=1,nzm

      if (k == nzm) then
        kb = k
      else
        kb = k+1
      endif

      do j=1,ny
        do i=1,nx

          wrk = 0.1*adzl(i,j,k)
                                                            ! Minimum 0.1 of local dz
          smixt(i,j,k) = max(wrk, min(max_eddy_length_scale,smixt(i,j,k))) 

          if (qcl(i,j,kb) == 0 .and. qcl(i,j,k) > 0 .and. brunt(i,j,k) > 1.e-4) then
!If just above the cloud top and atmosphere is stable, set to  0.1 of local dz
            smixt(i,j,k) = wrk
          endif

        end do ! i
      end do   ! j
    end do     ! k

  end subroutine check_eddy

  subroutine canuto()

! Subroutine impements an analytic expression for the third moment of SGS vertical velocity
! based on Canuto et at, 2001, JAS, 58, 1169-1172 (further referred to as C01)
! This allows to avoid having a prognostic equation for the third moment.
! Result is returned in a global variable w3 defined at the interface levels.
    
! Local variables
    integer i, j, k, kb, kc

    real bet2,   f0,     f1,  f2,    f3,   f4,  f5,  iso, isosqr,             &
         omega0,  omega1, omega2, X0,  Y0,    X1,   Y1,  AA0, AA1, buoy_sgs2, &
         thedz,   thedz2, cond,   wrk, wrk1,  wrk2, wrk3, avew
!
! See Eq. 7 in C01 (B.7 in Pete's dissertation)
    real, parameter :: c=7.0, a0=0.52/(c*c*(c-2.)), a1=0.87/(c*c),      &
                       a2=0.5/c, a3=0.6/(c*(c-2.)), a4=2.4/(3.*c+5.),   &
                       a5=0.6/(c*(3.*c+5))
!Moorthi               a5=0.6/(c*(3.+5.*c))
    
!   do k=1,nzm
    do k=2,nzm

      kb = k-1
      kc = k+1
       
      do j=1,ny
        do i=1,nx
             
          if(k == 1) then
            kb = 1
            kc = 2
            thedz  = adzl(i,j,kc)
            thedz2 = thedz
          elseif(k == nzm) then
            kb = nzm-1
            kc = nzm
            thedz  = adzl(i,j,k)
            thedz2 = thedz
!            thedz  = 1.0 / adzi(i,j,k)
!            thedz2 = 1.0 / adzl(i,j,k-1)
          else
             thedz  = adzl(i,j,k)
             thedz2 = adzl(i,j,kc)+adzl(i,j,k)
!             thedz  = 1.0 / adzi(i,j,k)                   ! Moorthi jul08
!             thedz2 = 1.0 / (adzl(i,j,k)+adzl(i,j,kb))      ! Moorthi jul08
          endif

!          if (abs(thedz).lt.1e-20) print *,'thedz=',thedz
!          if (abs(thedz2).lt.1e-20) print *,'thedz2=',thedz2
          if (abs(thedz).le.1e-10) thedz = sign(1e-10,thedz) 
          if (abs(thedz).eq.1e-10) print *,'thedz'
          if (abs(thedz2).le.1e-10) thedz2 = sign(1e-10,thedz2) 
          if (abs(thedz2).eq.1e-10) print *,'thedz2'

          iso       = 0.5*(isotropy(i,j,k)+isotropy(i,j,kb))
          isosqr    = iso*iso ! Two-level average of "return-to-isotropy" time scale squared
          buoy_sgs2 = isosqr*0.5*(brunt(i,j,k)+brunt(i,j,kb))
          bet2      = 0.5*(bet(i,j,k)+bet(i,j,kb))  !Two-level average of BV frequency squared

        
! Compute functions f0-f5, see Eq, 8 in C01 (B.8 in Pete's dissertation)
        

          avew = 0.5*(w_sec(i,j,k)+w_sec(i,j,kb))
          if (abs(avew).ge.1e10) avew = sign(1e10,avew) 
          if (abs(avew).eq.1e10) print *,'avew'
          cond = 1.2*sqrt(max(1.0e-16,2.*avew*avew*avew))
          wrk1 = bet2*iso
          wrk2 = thedz2*wrk1*wrk1*iso
          wrk3 = thl_sec(i,j,kc) - thl_sec(i,j,kb)

          f0   = wrk2 * wrk1 * wthl_sec(i,j,k) * wrk3

          wrk  = wthl_sec(i,j,kc) - wthl_sec(i,j,kb)
             
          f1   = wrk2 * (wrk*wthl_sec(i,j,k) + 0.5*avew*wrk3)
             
          wrk1 = bet2*isosqr
! JDC Modify:  different from SHOC and Canuto et al 2001
!          f2   = thedz*wrk1*wthl_sec(i,j,k)*(w_sec(i,j,k)-w_sec(i,j,kb))     &
!               + (thedz2+thedz2)*bet(i,j,k)*isosqr*wrk
          f2   = thedz*wrk1*wthl_sec(i,j,k)*(w_sec(i,j,k)-w_sec(i,j,kb))     &
               + (thedz2+thedz2)*bet(i,j,k)*isosqr*avew*wrk
! JDC modify: different from SHOC and Canuto et al 2001 
!          f3   = thedz2*wrk1*wrk + thedz*bet2*isosqr*(wthl_sec(i,j,k)*(tke(i,j,k)-tke(i,j,kb)))
          f3 = thedz2*wrk1*avew*wrk + thedz*wrk1*(wthl_sec(i,j,k)*(tke(i,j,k)-tke(i,j,kb)))

          wrk1 = thedz*iso*avew
          f4   = wrk1*(w_sec(i,j,k)-w_sec(i,j,kb) + tke(i,j,k)-tke(i,j,kb))
 
          f5   = wrk1*(w_sec(i,j,k)-w_sec(i,j,kb))
             
       
! Compute the "omega" terms, see Eq. 6 in C01 (B.6 in Pete's dissertation)

          dum = 1.-a5*buoy_sgs2
          if (abs(dum).le.1e-16) dum = sign(1e-16,dum) 
          if (abs(dum).eq.1e-16) print *,'1.-a5*buoy_sgs2'
          omega0 = a4 / dum
          omega1 = omega0 / (c+c)
          omega2 = omega1*f3+(5./4.)*omega0*f4
 
! Compute the X0, Y0, X1, Y1 terms,  see Eq. 5 a-b in C01  (B.5 in Pete's dissertation)

          dum = 1.-(a1+a3)*buoy_sgs2
          if (abs(dum).le.1e-16) dum = sign(1e-16,dum) 
          if (abs(dum).eq.1e-16) print *,'1.-(a1+a3)*buoy_sgs2'
          wrk1 = 1.0 / dum
          dum = 1.-a3*buoy_sgs2
          if (abs(dum).le.1e-16) dum = sign(1e-16,dum) 
          if (abs(dum).eq.1e-16) print *,'1.-a3*buoy_sgs2'
          wrk2 = 1.0 / dum
          X0   = wrk1 * (a2*buoy_sgs2*(1.-a3*buoy_sgs2))
          Y0   = wrk2 * (2.*a2*buoy_sgs2*X0)
          X1   = wrk1 * (a0*f0+a1*f1+a2*(1.-a3*buoy_sgs2)*f2)
          Y1   = wrk2 * (2.*a2*(buoy_sgs2*X1+(a0/a1)*f0+f1))

! Compute the A0, A1 terms,  see Eq. 5d in C01 (B.5 in Pete's dissertation)

          AA0 = omega0*X0 + omega1*Y0
          AA1 = omega0*X1 + omega1*Y1 + omega2

! Finally, we have the third moment of w, see Eq. 4c in C01 (B.4 in Pete's dissertation)
! cond is an estimate of third moment from second oment - If the third moment is larger
! than the estimate - limit w3.

           
           dum = c-1.2*X0+AA0
           if (abs(dum).le.1e-16) dum = sign(1e-16,dum)
           if (abs(dum).eq.1e-16) print *,'c-1.2*X0+AA0=',dum
!!! aab
!           w3(i,j,k) = max(-cond, min(cond, (AA1-1.2*X1-1.5*f5)/dum))
           w3(i,j,k) = (AA1-1.2*X1-1.5*f5)/(c-1.2*X0+AA0)
!!! aab


! Implemetation of the C01 approach in this subroutine is nearly complete
! (the missing part are Eqs. 5c and 5e which are very simple)
! therefore it's easy to diagnose other third order moments obtained in C01 using this code. 

        end do
      end do
    end do
    do j=1,ny
      do i=1,nx
        w3(i,j,1) = w3(i,j,2)
      enddo
    enddo
    
  end subroutine canuto


  subroutine buoyancy_single_gaussian()

! Compute SGS buoyancy flux and SGS cloud fraction, using analytic
! SINGLE-gaussian PDF for moisture and liquid water static energy.

! Local variables

    integer i,j,k,ku,kd
    real a,b,c,alpha,beta,cc,s,stds


! Initialize for statistics
!    do k=1,nzm
!      wqlsb(k) = 0.
!      wqisb(k) = 0.
!    enddo

    DO k=1,nzm
      
      kd = k
      ku = k + 1
      if (k == nzm) ku = k
      
      DO j=1,ny
        DO i=1,nx

! Initialize cloud variables to zero  
!          diag_qn   = 0.0
!          diag_frac = 0.0
!          diag_ql   = 0.0
!          diag_qi   = 0.0

          pval = prsl(i,j,k)
          pkap = (pval/100000.0) ** kapa
             
! Read in liquid/ice static energy, total water mixing ratio, 
! and vertical velocity to variables PDF needs
!          thl_first = hl(i,j,k)
!          qw_first  = total_water(i,j,k)
!          w_first   = w(i,j,k)

! GET ALL INPUT VARIABLES ON THE SAME GRID
! Points to be computed with relation to thermo point
! Read in points that need to be averaged

          thlsec   = max(0., 0.5*(thl_sec(i,j,kd)+thl_sec(i,j,ku)) )
          qwsec    = max(0., 0.5*(qw_sec(i,j,kd)+qw_sec(i,j,ku)) )
          qwthlsec = 0.5 * (qwthl_sec(i,j,kd) + qwthl_sec(i,j,ku))
          wqwsec   = 0.5 * (wqw_sec(i,j,kd)   + wqw_sec(i,j,ku))
          wthlsec  = 0.5 * (wthl_sec(i,j,kd)  + wthl_sec(i,j,ku))   


! Compute square roots of some variables so we don't have to do it again
!          if (w_sec(i,j,k) > 0.0) then
!            sqrtw2   = sqrt(w_sec(i,j,k))
!          else
!            sqrtw2   = 0.0
!          endif
!          if (thlsec > 0.0) then
!            sqrtthl  = sqrt(thlsec)
!          else
!            sqrtthl  = 0.0
!          endif
!          if (qwsec > 0.0) then
!            sqrtqt   = sqrt(qwsec)
!          else
!            sqrtqt   = 0.0
!          endif
             
          ! following Bechtold et al 1995
          b = tabs(i,j,k) - fac_cond*qcl(i,j,k)  ! Bechtold eqn 4
          c = MAPL_EQsat(b,prsl(i,j,k),dtqw)      
          s = total_water(i,j,k) - MAPL_EQsat(tabs(i,j,k),prsl(i,j,k))          

          a = 1. / (1.+fac_cond*dtqw)            ! Bechtold eqn 6
          b = a * (1./pkap) * dtqw
          c = a*(total_water(i,j,k) - c)

          alpha = 0.61 * (tabs(i,j,k)/pkap)      ! Bechtold eqn 12
          beta = (1./pkap) * fac_cond - 1.61*(tabs(i,j,k)/pkap)

          ! Bechtold eqn 
          stds = alpha**2 *qwsec + beta**2 * thlsec &
                 - 2.*alpha*beta*qwthlsec

          ! Cloud fraction, assuming single Gaussian qt distribution
          cc = 0.5 + 0.5*erf(s/(sqrt2*stds))  

          wthv_sec(i,j,k) = wthlsec*( 1.+0.61*total_water(i,j,k) &
                                     - beta*b*cc )               &
                           + wqwsec*( alpha + beta*a*cc )

        ENDDO
      ENDDO
    ENDDO

  end subroutine buoyancy_single_gaussian






! Saturation vapor pressure and mixing ratio subroutines
! Based on Flatau et al (1992), J. App. Met., 31, 1507-1513
! Code by Marat Khairoutdinov
 

  real function esatw(t)   !!! returns e in mb
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                       &
         6.11239921,       0.443987641,     0.142986287e-1, &
         0.264847430e-3,   0.302950461e-5,  0.206739458e-7, &
         0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real dt
    dt    = max(-80.,t-273.16)
    esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  end function esatw

  real function qsatw(t,p)
!    implicit none
    real t	! temperature (K)
    real p	! pressure    (Pa)
    real esat
!   esat  = fpvs(t)
!!    esat  = fpvsl(t)
    esat = MAPL_EQsat(t)
    qsatw = 0.622 * esat/max(esat,p-0.378*esat) 
!   esat  = esatw(t)
!   qsatw = 0.622 * esat/max(esat,p-esat) 
  end function qsatw
  
  
  real function esati(t)    !!! returns e in mb
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                     &
         6.11147274,     0.503160820,     0.188439774e-1, &
         0.420895665e-3, 0.615021634e-5,  0.602588177e-7, &
         0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
    real dt
!    real esatw
    if(t > 273.15) then
       esati = esatw(t)
    else if(t.gt.185.) then
       dt    = t-273.16
       esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
    else   ! use some additional interpolation below 184K
       dt    = max(-100.,t-273.16)
       esati = 0.00763685 + dt*(0.000151069+dt*7.48215e-07)
    end if
  end function esati
        
  real function qsati(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
    real esat !,esati
!   esat  = fpvs(t)
!!    esat  = fpvsi(t)
    esat = MAPL_EQsat(t,OverIce=.TRUE.)
    qsati = 0.622 * esat/max(esat,p-0.378*esat)
!   esat  = esati(t)
!   qsati = 0.622 * esat/max(esat,p-esat)
  end function qsati
  
  real function dtesatw(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                        &
         0.443956472,      0.285976452e-1,   0.794747212e-3, &
         0.121167162e-4,   0.103167413e-6,   0.385208005e-9, &
        -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
    real dt
    dt      = max(-80.,t-273.16)
    dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  end function dtesatw
        
  real function dtqsatw(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
!    real dtesatw
    dtqsatw = 100.0*0.622*dtesatw(t)/p
  end function dtqsatw
  
  real function dtesati(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                      &
         0.503223089,     0.377174432e-1,  0.126710138e-2, &
         0.249065913e-4,  0.312668753e-6,  0.255653718e-8, &
         0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
    real dt
!    real dtesatw
    if(t > 273.15) then
       dtesati = dtesatw(t)
    else if(t > 185.) then
       dt      = t-273.16
       dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
    else  ! use additional interpolation below 185K
       dt      = max(-100.,t-273.16)
       dtesati = 0.0013186 + dt*(2.60269e-05+dt*1.28676e-07)
    end if
  end function dtesati
  
  
  real function dtqsati(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
!    real dtesati
    dtqsati = 100.0*0.622*dtesati(t)/p
  end function dtqsati
  
 end subroutine run_shoc

end module shoc
