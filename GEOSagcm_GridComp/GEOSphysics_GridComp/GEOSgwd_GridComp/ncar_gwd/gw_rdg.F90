module gw_rdg

!
! This module handles gravity waves from orographic sources, and was
! extracted from gw_drag in May 2013.
!
!use shr_const_mod, only: pi => shr_const_pi
use gw_common, only: gw_drag_prof, GWBand, pi,cpair,rair, &
                     calc_taucd, momentum_flux, momentum_fixer, &
                     energy_momentum_adjust, energy_change, energy_fixer
use gw_utils, only:  GW_PRC, dot_2d, midpoint_interp


implicit none
private
save


! Public interface(s)
public :: gw_rdg_init
public :: gw_rdg_ifc

! Tunable Parameters
!--------------------
logical            :: do_divstream

!===========================================
! Parameters for DS2017 (do_divstream=.T.)
!===========================================
! Amplification factor - 1.0 for
! high-drag/windstorm regime
real, protected :: C_BetaMax_DS

! Max Ratio Fr2:Fr1 - 1.0
real, protected :: C_GammaMax

! Normalized limits  for Fr2(Frx) function
real, protected :: Frx0
real, protected :: Frx1


!===========================================
! Parameters for SM2000
!===========================================
! Amplification factor - 1.0 for
! high-drag/windstorm regime
real, protected :: C_BetaMax_SM



! NOTE: Critical inverse Froude number Fr_c is 
! 1./(SQRT(2.)~0.707 in SM2000
! (should be <= 1)
real, protected :: Fr_c


logical :: do_smooth_regimes
logical :: do_adjust_tauoro
logical :: do_backward_compat


! Limiters (min/max values)
! min surface displacement height for orographic waves (m)
real, protected :: orohmin
! min wind speed for orographic waves
real, protected :: orovmin
! min stratification allowing wave behavior
real, protected :: orostratmin
! min stratification allowing wave behavior
real, protected :: orom2min

real, parameter :: PINTADJ_0 = 0.9/19.0
real, parameter :: PINTADJ_1 = TAN(20.*pi/21.-0.5*pi)
real, parameter :: PINTADJ_2 = 0.5*pi
real, parameter :: PINTADJ_3 = 21./pi

! Some description of GW spectrum
type(GWBand)   :: band         ! I hate this variable  ... it just hides information from view


!==========================================================================
contains
!==========================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CCPP Interface routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!------------------------
!------------------------------------
!> \section arg_table_gw_rdg_init  Argument Table
!! \htmlinclude gw_rdg_init.html
subroutine gw_rdg_init (gw_dc, fcrit2, wavelength, pgwv)
#include <netcdf.inc>

  real, intent(in) :: gw_dc,fcrit2,wavelength
  integer, intent(in)  :: pgwv

  !==============================================
  !  Create "Band" structure
  !----------------------------------------------

  band  = GWBand(pgwv, gw_dc, fcrit2, wavelength )
 
  ! Set the local variables
  do_divstream        = .TRUE.
  C_BetaMax_DS        = 0.
  C_GammaMax          = 2.
  Frx0                = 2.
  Frx1                = 3.
  C_BetaMax_SM        = 2.
  Fr_c                = fcrit2
  do_smooth_regimes   = .FALSE.
  do_adjust_tauoro    = .TRUE.
  do_backward_compat  = .FALSE.
  orohmin             = 0.01
  orovmin             = 1.0e-3
  orostratmin         = 0.002
  orom2min            = 0.1
 
end subroutine gw_rdg_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!------------------------------------
!> \section arg_table_gw_rdg_ifc  Argument Table
!! \htmlinclude gw_rdg_ifc.html
subroutine gw_rdg_ifc( &
   ncol, pver, pverp, pcnst, n_rdg, dt, &
   u, v, t, pint, pmid, delp, rdelp, &
   piln, zm, zi, z, &
   ni, nm, rhoi, &
   kvtt,         &
   kwvrdg, effrdg, &
   hwdth, clngt, gbxar, &
   mxdis, angll, anixy, &
   rdg_cd_llb, trpd_leewv, &
   utrdg, vtrdg, ttrdg, &
   flx_heat )

   !!character(len=5), intent(in) :: type         ! BETA or GAMMA
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: pverp        ! Layer Vertical dimension
   integer,          intent(in) :: pver         ! Intfc Vertical dimension
   integer,          intent(in) :: pcnst        ! constituent dimension
   integer,          intent(in) :: n_rdg
   real,         intent(in) :: dt           ! Time step.

   real,         intent(in) :: u(ncol,pver)    ! Midpoint zonal winds. ( m s-1)
   real,         intent(in) :: v(ncol,pver)    ! Midpoint meridional winds. ( m s-1)
   real,         intent(in) :: t(ncol,pver)    ! Midpoint temperatures. (K)
   real,         intent(in) :: delp(ncol,pver) ! Delta(interface pressures).
   real,         intent(in) :: rdelp(ncol,pver) ! inverse Delta.
   real,         intent(in) :: pmid(ncol,pver)   ! midpoint pressures.
   real,         intent(in) :: pint(ncol,pverp)  ! interface pressures.
   real,         intent(in) :: piln(ncol,pverp)  ! Log of interface pressures.
   real,         intent(in) :: zm(ncol,pver)   ! Midpoint altitudes above ground (m).
   real,         intent(in) :: zi(ncol,pverp) ! Interface altitudes above ground (m).
   real,         intent(in) :: z(ncol)          ! Surface elevation (m).
   real, intent(in) :: kvtt(ncol,pverp) ! Molecular thermal diffusivity.
   !!real,         intent(in) :: q(ncol,pver,pcnst) ! Constituent array.
   !!real,         intent(in) :: dse(ncol,pver)  ! Dry static energy.

   real, intent(in) :: nm(ncol,pver)   ! Midpoint altitudes above ground (m).
   real, intent(in) :: ni(ncol,pverp) ! Interface altitudes above ground (m).
   real, intent(in) :: rhoi(ncol,pverp) ! Interface altitudes above ground (m).

   real,         intent(in) :: kwvrdg(ncol,n_rdg) ! horiz wavenumber.
   real,         intent(in) :: effrdg(ncol,n_rdg) ! efficiency factor of ridge scheme.
   real,         intent(in) :: hwdth(ncol,n_rdg) ! width of ridges.
   real,         intent(in) :: clngt(ncol,n_rdg) ! length of ridges.
   real,         intent(in) :: gbxar(ncol)      ! gridbox area

   real,         intent(in) :: mxdis(ncol,n_rdg) ! Height estimate for ridge (m).
   real,         intent(in) :: angll(ncol,n_rdg) ! orientation of ridges.
   real,         intent(in) :: anixy(ncol,n_rdg) ! Anisotropy parameter.

   real,         intent(in) :: rdg_cd_llb      ! Drag coefficient for low-level flow
   logical,      intent(in) :: trpd_leewv


   ! OUTPUTS
   real, intent(out) :: utrdg(ncol,pver)       ! Cum. zonal wind tendency
   real, intent(out) :: vtrdg(ncol,pver)       ! Cum. meridional wind tendency
   real, intent(out) :: ttrdg(ncol,pver)       ! Cum. temperature tendency
   !!real,       intent(out) :: qtrdg(ncol,pver,pcnst) ! Cum. consituent tendencies
   real, intent(inout) :: flx_heat(ncol)       ! Energy change

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn, icnst

   real(GW_PRC), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(GW_PRC), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(GW_PRC), allocatable :: c(:,:)

   ! Isotropic source flag [anisotropic orography].
   integer  :: isoflag(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)
   integer :: bwv_level(ncol)
   integer :: tlb_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real :: ubm(ncol,pver)
   real :: ubi(ncol,pverp)

   ! Unit vectors of source wind (zonal and meridional components).
   real :: xv(ncol)
   real :: yv(ncol)

   ! Averages over source region.
   real :: ubmsrc(ncol) ! On-ridge wind.
   real :: usrc(ncol)   ! Zonal wind.
   real :: vsrc(ncol)   ! Meridional wind.
   real :: nsrc(ncol)   ! B-V frequency.
   real :: rsrc(ncol)   ! Density.

   ! normalized wavenumber
   real :: m2src(ncol)

   ! Top of low-level flow layer.
   real :: tlb(ncol)

   ! Bottom of linear wave region.
   real :: bwv(ncol)

   ! Froude numbers for flow/drag regimes
   real :: Fr1(ncol)
   real :: Fr2(ncol)
   real :: Frx(ncol)

   ! Wave Reynolds stresses at source level
   real :: tauoro(ncol)
   real :: taudsw(ncol)

   ! Surface streamline displacement height for linear waves.
   real :: hdspwv(ncol)

   ! Surface streamline displacement height for downslope wind regime.
   real :: hdspdw(ncol)

   ! Wave breaking level
   real :: wbr(ncol)

   ! Momentum fluxes used by fixer.
   real :: um_flux(ncol), vm_flux(ncol)

   ! Energy change used by fixer.
   real :: de(ncol)

   ! Reynolds stress for waves propagating in each cardinal direction.
   real :: taucd(ncol,pver+1,4)

   real :: utgw(ncol,pver)       ! zonal wind tendency
   real :: vtgw(ncol,pver)       ! meridional wind tendency
   real :: ttgw(ncol,pver)       ! temperature tendency
#ifdef CAM
   real :: qtgw(ncol,pver,pcnst) ! constituents tendencies
#endif

   real :: pint_adj(ncol,pver+1)
   real :: zfac_layer

   character(len=4) :: type         ! BETA or GAMMA (just BETA for now)
   character(len=1) :: cn
   character(len=9) :: fname(4)
   !----------------------------------------------------------------------------

   ! Allocate wavenumber fields.
   allocate(tau(ncol,  -band%ngwv:band%ngwv  , pverp))
   allocate(gwut(ncol,pver,-band%ngwv:band%ngwv  ))
   allocate(c(ncol,-band%ngwv:band%ngwv))

   type='BETA'

   ! initialize accumulated momentum fluxes and tendencies
   utrdg = 0.
   vtrdg = 0.
   ttrdg = 0.
  
!GEOS pressure scaling near model top
!  zfac_layer = 1000.0 ! 10mb
!  pint_adj = (pint/zfac_layer)**3 
!  pint_adj = 0.5*(1+TANH(((2.0*pint/zfac_layer)-1)/0.25))
   pint_adj = 1.0
!  pint_adj = 0.1 + (0.9/19.0) * &
!                   ((ATAN((2.*(pint-10000.0)/(75000.0-10000.0)-1.) * TAN(20.*pi/21.-0.5*pi)) + 0.5*pi) * 21./pi - 1.)
!  adjust strength from surface (1.0) to 10,000m (0.1)
!  do k=1,pver+1 
!    pint_adj(:,k)= 0.1 + PINTADJ_0 * ((ATAN((2.*(z+zi(:,k)-2500.0)/(-2500.0)-1.) * PINTADJ_1) + PINTADJ_2) * PINTADJ_3 - 1.)
!  enddo

   isoflag = 0
 
   do nn = 1, n_rdg
  
    call gw_rdg_src(ncol, pver, pint, pmid, delp, &
         u, v, t, mxdis(:,nn), angll(:,nn), anixy(:,nn), kwvrdg(:,nn), isoflag, zi, nm, &
         src_level, tend_level, bwv_level, tlb_level, tau, ubm, ubi, xv, yv,  & 
         ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, c)

    call gw_rdg_belowpeak(ncol, pver, rdg_cd_llb, &
         t, mxdis(:,nn), anixy(:,nn), kwvrdg(:,nn), & 
         zi, nm, ni, rhoi, &
         src_level, tau, & 
         ubmsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, & 
         tauoro, taudsw, hdspwv, hdspdw)

    call gw_rdg_break_trap(ncol, pver, &
         zi, nm, ni, ubm, ubi, rhoi, kwvrdg(:,nn) , bwv, tlb, wbr, & 
         src_level, tlb_level, hdspwv, hdspdw,  mxdis(:,nn), & 
         tauoro, taudsw, tau, & 
         ldo_trapped_waves=trpd_leewv)

     call gw_drag_prof(ncol, pver, band, pint, delp, rdelp, & 
          src_level, tend_level,dt, t, &
          piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
          effrdg(:,nn), c, kvtt, tau, utgw, vtgw, &
          ttgw, gwut, &
          kwvrdg=kwvrdg(:,nn), satfac_in=1.0, tau_adjust=pint_adj, &
          tndmax_in=40.0/86400.0)

#ifdef NCAR_ADJUST
   ! ! Project stress into directional components.
   ! taucd = calc_taucd(ncol, pver, band%ngwv, tend_level, tau, c, xv, yv, ubi)

   ! ! Find momentum flux, and use it to fix the wind tendencies below
   ! ! the gravity wave region.
   ! call momentum_flux(tend_level, taucd, um_flux, vm_flux)
   ! call momentum_fixer(ncol, pver, tend_level, pint, um_flux, vm_flux, utgw, vtgw)

   ! ! Find energy change in the current state, and use fixer to apply
   ! ! the difference in lower levels.
   ! call energy_change(ncol, pver, dt, delp, u, v, utgw, vtgw, ttgw, de)
   ! call energy_fixer(ncol, pver, tend_level, pint, de-flx_heat, ttgw)
   ! flx_heat=de
#else
   ! call energy_momentum_adjust(ncol, pver, src_level, band, pint, delp, c, tau, &
   !                    effrdg(:,nn), t, ubm, ubi, xv, yv, utgw, vtgw, ttgw)
#endif

     ! Add the tendencies from each ridge to the totals.
     do k = 1, pver
        utrdg(:,k) = utrdg(:,k) + utgw(:,k)
        vtrdg(:,k) = vtrdg(:,k) + vtgw(:,k)
        ttrdg(:,k) = ttrdg(:,k) + ttgw(:,k)
     end do

#ifdef CAM
! disable tracer mixing in GW for now.
      do icnst = 1, pcnst
      do k = 1, pver
         qtrdg(:,k,icnst) = qtrdg(:,k,icnst) + qtgw(:,k,icnst)
      end do
      end do
      if (nn <= 6) then
         write(cn, '(i1)') nn
      end if
#endif

   end do ! end of loop over multiple ridges

   if (trim(type) == 'BETA') then
      fname(1) = 'TAUGWX'
      fname(2) = 'TAUGWY'
      fname(3) = 'UTGWORO'
      fname(4) = 'VTGWORO'
   else if (trim(type) == 'GAMMA') then
      fname(1) = 'TAURDGGMX'
      fname(2) = 'TAURDGGMY'
      fname(3) = 'UTRDGGM'
      fname(4) = 'VTRDGGM'
   else
      call endrun('gw_rdg_calc: FATAL: type must be either BETA or GAMMA'&
                  //' type= '//type)
   end if

   deallocate(tau, gwut, c)

 end subroutine gw_rdg_ifc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Non - interface subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------
subroutine gw_rdg_src(ncol, pver , pint, pmid, delp, &
     u, v, t, mxdis, angxy, anixy, kwvrdg, iso, zi, nm, &
     src_level, tend_level, bwv_level ,tlb_level , tau, ubm, ubi, xv, yv,  & 
     ubmsrc, usrc, vsrc, nsrc, rsrc, m2src, tlb, bwv, Fr1, Fr2, Frx, c)

  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Vertical dimension.
  integer, intent(in) :: pver

  ! Interface pressures. (Pa)
  real, intent(in) :: pint(ncol,pver+1)
  ! Midpoint pressures. (Pa)
  real, intent(in) :: pmid(ncol,pver)
  ! Delta Interface pressures. (Pa)
  real, intent(in) :: delp(ncol,pver)


  ! Midpoint zonal/meridional winds. ( m s-1)
  real, intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures. (K)
  real, intent(in) :: t(ncol,pver)
  ! Height estimate for ridge (m) [anisotropic orography].
  real, intent(in) :: mxdis(ncol)
  ! Angle of ridge axis w/resp to north (degrees) [anisotropic orography].
  real, intent(in) :: angxy(ncol)
  ! Anisotropy parameter [anisotropic orography].
  real, intent(in) :: anixy(ncol)
  ! horiz wavenumber [anisotropic orography].
  real, intent(in) :: kwvrdg(ncol)
  ! Isotropic source flag [anisotropic orography].
  integer, intent(in)  :: iso(ncol)
  ! Interface altitudes above ground (m).
  real, intent(in) :: zi(ncol,pver+1)
  ! Midpoint Brunt-Vaisalla frequencies (s-1).
  real, intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)
  integer, intent(out) :: bwv_level(ncol),tlb_level(ncol)

  ! Averages over source region.
  real, intent(out) :: nsrc(ncol) ! B-V frequency.
  real, intent(out) :: rsrc(ncol) ! Density.
  real, intent(out) :: usrc(ncol) ! Zonal wind.
  real, intent(out) :: vsrc(ncol) ! Meridional wind.
  real, intent(out) :: ubmsrc(ncol) ! On-ridge wind.
  ! Top of low-level flow layer.
  real, intent(out) :: tlb(ncol)
  ! Bottom of linear wave region.
  real, intent(out) :: bwv(ncol)
  ! normalized wavenumber
  real, intent(out) :: m2src(ncol)


  ! Wave Reynolds stress.
  real(GW_PRC), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real, intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(GW_PRC), intent(out) :: c(ncol,-band%ngwv:band%ngwv)
  ! Froude numbers for flow/drag regimes
  real, intent(out) :: Fr1(ncol), Fr2(ncol), Frx(ncol)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Surface streamline displacement height (2*sgh).
  real :: hdsp(ncol)

  ! Difference in interface pressure across source region.
  real :: dpsrc(ncol)
  ! Thickness of downslope wind region.
  real :: ddw(ncol)
  ! Thickness of linear wave region.
  real :: dwv(ncol)
  ! Wind speed in source region.
  real :: wmsrc(ncol)

  real :: ragl(ncol) 
  
!--------------------------------------------------------------------------
! Check that ngwav is equal to zero, otherwise end the job
!--------------------------------------------------------------------------
  !!  if (band%ngwv /= 0) call endrun(' gw_rdg_src :: ERROR - band%ngwv must be zero and it is not')

!--------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the appropiate
! values of wind, stability, etc. for determining the wave source are
! averages over the depth of the atmosphere penterated by the typical
! mountain.
! Reduces to the bottom midpoint values when mxdis=0, such as over ocean.
!--------------------------------------------------------------------------

  hdsp      = mxdis ! no longer multipied by 2
  src_level = pver+1
  bwv_level = -1
  tlb_level = -1

  tau(:,0,:) = 0.0

  ! Find depth of "source layer" for mountain waves
  ! i.e., between ground and mountain top
  do k = pver, 1, -1
     do i = 1, ncol
        ! Need to have h >= z(k+1) here or code will bomb when h=0.
        if ( (hdsp(i) >= zi(i,k+1)) .and. (hdsp(i) < zi(i,k))   ) then
           src_level(i) = k  
        end if
     end do
  end do

  rsrc = 0.
  usrc = 0. 
  vsrc = 0.
  nsrc = 0.
  do i = 1, ncol
      do k = pver, src_level(i), -1
           rsrc(i) = rsrc(i) + pmid(i,k) / ( rair * t(i,k))* delp(i,k)
           usrc(i) = usrc(i) + u(i,k) * delp(i,k)
           vsrc(i) = vsrc(i) + v(i,k) * delp(i,k)
           nsrc(i) = nsrc(i) + nm(i,k)* delp(i,k)
     end do
  end do


  do i = 1, ncol
     dpsrc(i) = pint(i,pver+1) - pint(i,src_level(i))
  end do

  rsrc = rsrc / dpsrc
  usrc = usrc / dpsrc
  vsrc = vsrc / dpsrc
  nsrc = nsrc / dpsrc

  wmsrc = sqrt( usrc**2 + vsrc**2 )


  ! Get the unit vector components
  ! Want agl=0 with U>0 to give xv=1

  ragl = angxy * pi/180.

  ! protect from wierd "bad" angles 
  ! that may occur if hdsp is zero
  where( hdsp <= orohmin )
     ragl = 0.
  end where

  yv   =-sin( ragl )
  xv   = cos( ragl )


  ! Kluge in possible "isotropic" obstacle.
  where( ( iso == 1 ) .and. (wmsrc > orovmin) )
       xv = usrc/wmsrc    
       yv = vsrc/wmsrc
  end where


  ! Project the local wind at midpoints into the on-ridge direction
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do
  ubmsrc = dot_2d(usrc , vsrc , xv, yv)

  ! Ensure on-ridge wind is positive at source level
  do k = 1, pver
     ubm(:,k) = sign( ubmsrc*0.+1. , ubmsrc ) *  ubm(:,k)
  end do

                  ! Sean says just use 1. as 
                  ! first argument
  xv  = sign( ubmsrc*0.+1. , ubmsrc ) *  xv
  yv  = sign( ubmsrc*0.+1. , ubmsrc ) *  yv

  ! Now make ubmsrc positive and protect
  ! against zero
  ubmsrc = abs(ubmsrc)
  ubmsrc = max( 0.01 , ubmsrc )
  

  ! The minimum stratification allowing GW behavior
  ! should really depend on horizontal scale since
  !
  !      m^2 ~ (N/U)^2 - k^2
  !
  ! Should also think about parameterizing
  ! trapped lee-waves.  

  
  ! This needs to be made constistent with later
  ! treatment of nonhydrostatic effects.
  m2src = ( (nsrc/(ubmsrc+0.01))**2 - kwvrdg**2 ) /((nsrc/(ubmsrc+0.01))**2)


  !-------------------------------------------------------------
  ! Calculate provisional limits (in Z [m]) for 3 regimes. This
  ! will modified later if wave breaking or trapping are
  ! diagnosed
  !
  !                                            ^ 
  !                                            | *** linear propagation ***
  !  (H) -------- mountain top -------------   | *** or wave breaking  ****     
  !                                            | *** regimes  *************
  ! (BWV)------ bottom of linear waves ----    |
  !                    :                       |
  !                 *******                    |
  !                    :                       |
  ! (TLB)--- top of flow diversion layer---    '
  !                   :
  !        **** flow diversion *****  
  !                    :
  !============================================

  !============================================
  ! For Dividing streamline para (DS2017)
  !--------------------------------------------
  ! High-drag downslope wind regime exists
  ! between bottom of linear waves and top of
  ! flow diversion. Linear waves can only 
  ! attain vertical displacment of f1*U/N. So,
  ! bottom of linear waves is given by
  !
  !        BWV = H - Fr1*U/N 
  !
  ! Downslope wind layer begins at BWV and 
  ! extends below it until some maximum high
  ! drag obstacle height Fr2*U/N is attained
  ! (where Fr2 >= f1).  Below downslope wind
  ! there is flow diversion, so top of 
  ! diversion layer (TLB) is equivalent to
  ! bottom of downslope wind layer and is;
  !
  !       TLB = H - Fr2*U/N
  !
  !-----------------------------------------

  ! Critical inverse Froude number
  !-----------------------------------------------
  Fr1(:) = Fr_c * 1.00
  Frx(:) = hdsp(:)*nsrc(:)/abs( ubmsrc(:) ) / Fr_c

  if ( do_divstream ) then
     !------------------------------------------------
     ! Calculate Fr2(Frx) for DS2017   
     !------------------------------------------------
     where(Frx <= Frx0)
          Fr2(:) = Fr1(:) + Fr1(:)* C_GammaMax * anixy(:)
     elsewhere((Frx > Frx0).and.(Frx <= Frx1) )
          Fr2(:) = Fr1(:) + Fr1(:)* C_GammaMax * anixy(:) &
                   * (Frx1 - Frx(:))/(Frx1-Frx0)    
     elsewhere(Frx > Frx1) 
          Fr2(:)=Fr1(:)
     endwhere
  else
  !------------------------------------------   
  ! Regime distinctions entirely carried by
  ! amplification of taudsw (next subr)
  !------------------------------------------
     Fr2(:)=Fr1(:)
  end if   


  
  where( m2src > orom2min ) 
     ddw  = Fr2 * ( abs(ubmsrc) )/nsrc
  elsewhere
     ddw  = 0.
  endwhere


  ! If TLB is less than zero then obstacle is not
  ! high enough to produce an low-level diversion layer
  tlb = mxdis - ddw
  where( tlb < 0.)
     tlb = 0.
  endwhere
  do k = pver, pver/2, -1
     do i = 1, ncol
         if ( (tlb(i) > zi(i,k+1)) .and. (tlb(i) <= zi(i,k))   ) then
           tlb_level(i) = k
        end if
     end do
  end do


  ! Find *BOTTOM* of linear wave layer (BWV)
  !where ( nsrc > orostratmin )
  where( m2src > orom2min ) 
      dwv  = Fr1 * ( abs(ubmsrc) )/nsrc
  elsewhere
     dwv  = -9.999e9 ! if weak strat - no waves
  endwhere

  bwv = mxdis - dwv
  where(( bwv < 0.) .or. (dwv < 0.) )
     bwv = 0.
  endwhere
  do k = pver,1, -1
     do i = 1, ncol
        if ( (bwv(i) > zi(i,k+1)) .and. (bwv(i) <= zi(i,k))   ) then
           bwv_level(i) = k+1
        end if
     end do
  end do



  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)
  ubi(:,2:pver) = midpoint_interp(ubm)
  ubi(:,pver+1) = ubm(:,pver)

  ! Allow wind tendencies all the way to the model bottom.
  tend_level = pver

  ! No spectrum; phase speed is just 0.
  c = 0.

  where( m2src < orom2min ) 
     tlb = mxdis
     tlb_level = src_level
  endwhere


end subroutine gw_rdg_src


!==========================================================================

subroutine gw_rdg_belowpeak(ncol, pver, rdg_cd_llb, &
     t, mxdis, anixy, kwvrdg, zi, nm, ni, rhoi, &
     src_level , tau,  & 
     ubmsrc, nsrc, rsrc, m2src,tlb,bwv,Fr1,Fr2,Frx, & 
     tauoro,taudsw, hdspwv,hdspdw  )

  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Vertical dimension.
  integer, intent(in) :: pver
  ! Drag coefficient for low-level flow
  real, intent(in) :: rdg_cd_llb


  ! Midpoint temperatures. (K)
  real, intent(in) :: t(ncol,pver)
  ! Height estimate for ridge (m) [anisotropic orography].
  real, intent(in) :: mxdis(ncol)
  ! Anisotropy parameter [0-1] [anisotropic orography].
  real, intent(in) :: anixy(ncol)
  ! Inverse cross-ridge lengthscale (m-1) [anisotropic orography].
  real, intent(in) :: kwvrdg(ncol)
  ! Interface altitudes above ground (m).
  real, intent(in) :: zi(ncol,pver+1)
  ! Midpoint Brunt-Vaisalla frequencies (s-1).
  real, intent(in) :: nm(ncol,pver)
  ! Interface Brunt-Vaisalla frequencies (s-1).
  real, intent(in) :: ni(ncol,pver+1)
  ! Interface density (kg m-3).
  real, intent(in) :: rhoi(ncol,pver+1)

  ! Indices of top gravity wave source level
  integer, intent(inout) :: src_level(ncol)

  ! Wave Reynolds stress.
  real(GW_PRC), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Top of low-level flow layer.
  real, intent(inout) :: tlb(ncol)
  ! Bottom of linear wave region.
  real, intent(inout) :: bwv(ncol)
  ! surface stress from linear waves.
  real, intent(out) :: tauoro(ncol)
  ! surface stress for downslope wind regime.
  real, intent(out) :: taudsw(ncol)

  ! Surface streamline displacement height for linear waves.
  real, intent(out) :: hdspwv(ncol)
  ! Surface streamline displacement height for downslope wind regime.
  real, intent(out) :: hdspdw(ncol)



  ! Froude numbers for flow/drag regimes
  real, intent(in) :: Fr1(ncol), Fr2(ncol),Frx(ncol)

  ! Averages over source region.
  real, intent(in) :: m2src(ncol) ! normalized non-hydro wavenumber
  real, intent(in) :: nsrc(ncol)  ! B-V frequency.
  real, intent(in) :: rsrc(ncol)  ! Density.
  real, intent(in) :: ubmsrc(ncol) ! On-ridge wind.


  !logical, intent(in), optional :: forcetlb

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  real :: Coeff_LB(ncol),tausat,ubsrcx(ncol),dswamp
  real :: taulin(ncol),BetaMax

  ! ubsrcx introduced to account for situations with high shear, strong strat.
  do i = 1, ncol
        ubsrcx(i)    = max( ubmsrc(i)  , 0. )
  end do

  do i = 1, ncol
     if ( m2src(i) > orom2min )   then 
        hdspwv(i) = min( mxdis(i) , Fr1(i) * ubsrcx(i) / nsrc(i) )
     else
        hdspwv(i) = 0.
     end if
  end do
  
  if (do_divstream) then
     do i = 1, ncol
        if ( m2src(i) > orom2min )   then 
           hdspdw(i) = min( mxdis(i) , Fr2(i) * ubsrcx(i) / nsrc(i) )
        else
           hdspdw(i) = 0.
        end if
     end do
  else
     do i = 1, ncol
        ! Needed only to mark where a DSW occurs
        if ( m2src(i) > orom2min )   then 
           hdspdw(i) = mxdis(i) 
        else
           hdspdw(i) = 0.
        end if
     end do
  end if

  ! Calculate form drag coefficient ("CD")
  !--------------------------------------
  Coeff_LB = rdg_cd_llb*anixy

  ! Determine the orographic c=0 source term following McFarlane (1987).
  ! Set the source top interface index to pver, if the orographic term is
  ! zero.
  ! 
  ! This formula is basically from
  !
  !      tau(src) = rho * u' * w'
  ! where 
  !      u' ~ N*h'  and w' ~ U*h'/b  (b="breite")
  !
  ! and 1/b has been replaced with k (kwvrdg) 
  !
  do i = 1, ncol
     if ( ( src_level(i) > 0 ) .and. ( m2src(i) > orom2min ) ) then
        tauoro(i) = kwvrdg(i) * ( hdspwv(i)**2 ) * rsrc(i) * nsrc(i) &
             * ubsrcx(i)
        taudsw(i) = kwvrdg(i) * ( hdspdw(i)**2 ) * rsrc(i) * nsrc(i) &
             * ubsrcx(i)
     else
        tauoro(i) = 0.
        taudsw(i) = 0.
     end if
  end do

  if (do_divstream) then
     do i = 1, ncol
           taulin(i) = 0.
     end do
  !---------------------------------------
  ! Need linear drag when divstream is not used
  !---------------------------------------
  else
     do i = 1, ncol
        if ( ( src_level(i) > 0 ) .and. ( m2src(i) > orom2min ) ) then
           taulin(i) = kwvrdg(i) * ( mxdis(i)**2 ) * rsrc(i) * nsrc(i) &
                * ubsrcx(i)
        else
           taulin(i) = 0.
        end if
     end do
  end if

  if ( do_divstream ) then
  ! Amplify DSW between Frx=1. and Frx=Frx1
     do i = 1,ncol
        dswamp=0.
        BetaMax   = C_BetaMax_DS * anixy(i)      
        if ( (Frx(i)>1.).and.(Frx(i)<=Frx1)) then
           dswamp = (Frx(i)-1.)*(Frx1-Frx(i))/(0.25*(Frx1-1.)**2)
        end if
        taudsw(i) = (1. + BetaMax*dswamp)*taudsw(i)
     end do
  else
  !-------------------
  ! Scinocca&McFarlane
  !--------------------
     do i = 1, ncol
        BetaMax   = C_BetaMax_SM * anixy(i)      
        if ( (Frx(i) >=1.) .and. (Frx(i) < 1.5) ) then
           dswamp = 2. * BetaMax * (Frx(i) -1.)
        else if ( ( Frx(i) >= 1.5 ) .and. (Frx(i) < 3. ) ) then
           dswamp = ( 1. + BetaMax - (0.666**2) ) * ( 0.666*(3. - Frx(i) ))**2  & 
                      + ( 1. / Frx(i) )**2  -1.
        else
           dswamp    = 0.      
        end if
        if ( (Frx(i) >=1.) .and. (Frx(i) < 3.) ) then
          taudsw(i) = (1. + dswamp )*taulin(i) - tauoro(i)
        else
          taudsw(i) = 0.   
        endif
        ! This code defines "taudsw" as SUM of freely-propagating
        ! waves +DSW enhancement. Different than in SM2000
        taudsw(i) = taudsw(i) + tauoro(i) 
     end do
 !----------------------------------------------------
  end if

  
  do i = 1, ncol
     if ( m2src(i) > orom2min )   then 
        where ( ( zi(i,:) < mxdis(i) ) .and. ( zi(i,:) >= bwv(i) ) )
             tau(i,0,:) =  tauoro(i)
        else where ( ( zi(i,:) < bwv(i) ) .and. ( zi(i,:) >= tlb(i) ) )
             tau(i,0,:) =  tauoro(i) +( taudsw(i)-tauoro(i) )* &
                                         ( bwv(i) - zi(i,:) ) / &
                                         ( bwv(i) - tlb(i) )
        endwhere
        ! low-level form drag on obstacle. Quantity kwvrdg (~1/b) appears for consistency
        ! with tauoro and taudsw forms. Should be weighted by L*b/A_g before applied to flow.
        where ( ( zi(i,:) < tlb(i) ) .and. ( zi(i,:) >= 0. ) )
             tau(i,0,:) =  taudsw(i) +  &
                           Coeff_LB(i) * kwvrdg(i) * rsrc(i) * 0.5 * (ubsrcx(i)**2) * ( tlb(i) - zi(i,:) )
        endwhere
 
        if (do_smooth_regimes) then
        !  This blocks accounts for case where both mxdis and tlb fall
        !  between adjacent edges
           do k=1,pver
              if ( (zi(i,k) >= tlb(i)).and.(zi(i,k+1) < tlb(i)).and. &
                   (zi(i,k) >= mxdis(i)).and.(zi(i,k+1) < mxdis(i)) ) then
                 src_level(i) = src_level(i)-1
                 tau(i,0,k) = tauoro(i)
              end if
           end do
        end if 

     else     !----------------------------------------------
             ! This block allows low-level dynamics to occur
             ! even if m2 is less than orom2min
        where ( ( zi(i,:) < tlb(i) ) .and. ( zi(i,:) >= 0. ) )
               tau(i,0,:) =  taudsw(i) +  &
                   Coeff_LB(i) * kwvrdg(i) * rsrc(i) * 0.5 * &
                   (ubsrcx(i)**2) * ( tlb(i) - zi(i,:) )
        endwhere
     endif
  end do

  ! This may be redundant with newest version of gw_drag_prof.
  ! That code reaches down to level k=src_level+1. (jtb 1/5/16)
  do i = 1, ncol
     k=src_level(i)
     if ( ni(i,k) > orostratmin ) then
         tausat    =  (Fr_c**2) * kwvrdg(i) * rhoi(i,k) * ubsrcx(i)**3 / &
              (1.*ni(i,k)) 
     else
         tausat = 0.
     endif 
     tau(i,0,src_level(i)) = min( tauoro(i), tausat ) 
  end do



  ! Final clean-up. Do nothing if obstacle less than orohmin
  do i = 1, ncol
     if ( mxdis(i) < orohmin ) then
        tau(i,0,:) = 0. 
        tauoro(i)  = 0.
        taudsw(i)  = 0.
     endif 
  end do

          ! Disable vertical propagation if Scorer param is 
          ! too small.
  do i = 1, ncol
     if ( m2src(i) <= orom2min ) then
        src_level(i)=1
     endif 
  end do



end subroutine gw_rdg_belowpeak

!==========================================================================
subroutine gw_rdg_break_trap(ncol, pver, &
     zi, nm, ni, ubm, ubi, rhoi, kwvrdg, bwv, tlb, wbr, & 
     src_level, tlb_level, & 
     hdspwv, hdspdw, mxdis, &
     tauoro, taudsw,  tau, & 
     ldo_trapped_waves, wdth_kwv_scale_in )
  !-----------------------------------------------------------------------
  ! Parameterization of high-drag regimes and trapped lee-waves for CAM
  !
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol
  ! Vertical dimension.
  integer, intent(in) :: pver

  ! Horz wavenumber for ridge (1/m) [anisotropic orography].
  real, intent(in) :: kwvrdg(ncol)
  ! Interface altitudes above ground (m).
  real, intent(in) :: zi(ncol,pver+1)
  ! Midpoint Brunt-Vaisalla frequencies (s-1).
  real, intent(in) :: nm(ncol,pver)
  ! Interface Brunt-Vaisalla frequencies (s-1).
  real, intent(in) :: ni(ncol,pver+1)

  ! Indices of gravity wave sources.
  integer, intent(inout) :: src_level(ncol), tlb_level(ncol)

  ! Wave Reynolds stress.
  real(GW_PRC), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Wave Reynolds stresses at source.
  real, intent(inout) :: taudsw(ncol),tauoro(ncol)
  ! Projection of wind at midpoints and interfaces.
  real, intent(in) :: ubm(ncol,pver)
  real, intent(in) :: ubi(ncol,pver+1)
  ! Interface density (kg m-3).
  real, intent(in) :: rhoi(ncol,pver+1)

  ! Top of low-level flow layer.
  real, intent(in) :: tlb(ncol)
  ! Bottom of linear wave region.
  real, intent(in) :: bwv(ncol)

  ! Surface streamline displacement height for linear waves.
  real, intent(in) :: hdspwv(ncol)
  ! Surface streamline displacement height for downslope wind regime.
  real, intent(in) :: hdspdw(ncol)
  ! Ridge height.
  real, intent(in) :: mxdis(ncol)


  ! Wave breaking level
  real, intent(out) :: wbr(ncol)

  logical, intent(in), optional :: ldo_trapped_waves
  real, intent(in), optional :: wdth_kwv_scale_in

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k, kp1, non_hydro
  real:: m2(ncol,pver),delz(ncol),tausat(ncol),trn(ncol)
  real:: wbrx(ncol)
  real:: phswkb(ncol,pver+1)
  logical :: lldo_trapped_waves
  real:: wdth_kwv_scale
  ! Indices of important levels.
  integer :: trn_level(ncol)

  if (present(ldo_trapped_waves)) then
     lldo_trapped_waves = ldo_trapped_waves
     if(lldo_trapped_waves) then
       non_hydro = 1
     else
       non_hydro = 0
     endif
  else
     lldo_trapped_waves = .false.
     non_hydro = 0
  endif

  if (present(wdth_kwv_scale_in)) then
     wdth_kwv_scale = wdth_kwv_scale_in
  else
     wdth_kwv_scale = 1.
  endif

  ! Calculate vertical wavenumber**2
  !---------------------------------
  m2 = (nm  / (abs(ubm)+.01))**2
  do k=pver,1,-1
     m2(:,k) = m2(:,k) - non_hydro*(wdth_kwv_scale*kwvrdg)**2
     ! sweeping up, zero out m2 above first occurence
     ! of m2(:,k)<=0
     kp1=min( k+1, pver )
     where( (m2(:,k) <= 0.0 ).or.(m2(:,kp1) <= 0.0 ) )
        m2(:,k) = 0.
     endwhere
  end do

  ! Take square root of m**2 and 
  ! do vertical integral to find
  ! WKB phase.
  !-----------------------------
  m2 = SQRT( m2 )
  phswkb(:,:)=0
  do k=pver,1,-1
     where( zi(:,k) > tlb(:) )
        delz(:) = min( zi(:,k)-zi(:,k+1) , zi(:,k)-tlb(:) ) 
        phswkb(:,k) = phswkb(:,k+1) + m2(:,k)*delz(:) 
     endwhere
  end do

  ! Identify top edge of layer in which phswkb reaches 3*pi/2
  ! - approximately the "breaking level"
  !----------------------------------------------------------
  wbr(:)=0.
  wbrx(:)=0.
  if (do_smooth_regimes) then
     do k=pver,1,-1
     where( (phswkb(:,k+1)<1.5*pi).and.(phswkb(:,k)>=1.5*pi) & 
            .and.(hdspdw(:)>hdspwv(:)) )
        wbr(:)  = zi(:,k)  
        ! Extrapolation to make regime
        ! transitions smoother
        wbrx(:) = zi(:,k)   - ( phswkb(:,k) -  1.5*pi ) &
                            / ( m2(:,k) + 1.e-6 )
        src_level(:) = k-1
     endwhere
     end do
  else
     do k=pver,1,-1
     where( (phswkb(:,k+1)<1.5*pi).and.(phswkb(:,k)>=1.5*pi) & 
            .and.(hdspdw(:)>hdspwv(:)) )
        wbr(:)  = zi(:,k)
        src_level(:) = k
     endwhere
     end do
  end if

  ! Adjust tauoro at new source levels if needed.
  ! This is problematic if Fr_c<1.0. Not sure why.
  !----------------------------------------------------------
  if (do_adjust_tauoro) then 
     do i = 1,ncol
        if (wbr(i) > 0. ) then
            tausat(i) = (Fr_c**2) * kwvrdg(i)  * rhoi( i, src_level(i) ) & 
                      * abs(ubi(i , src_level(i) ))**3  &
                      / ni( i , src_level(i) ) 
            tauoro(i) = min( tauoro(i), tausat(i) )
        end if
     end do
  end if

  if (do_smooth_regimes) then
     do i = 1, ncol
     do k=1,pver+1
        if ( ( zi(i,k) <= wbr(i) ) .and. ( zi(i,k) > tlb(i) ) ) then
           tau(i,0,k) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                          ( wbrx(i) - zi(i,k) ) / &
                          ( wbrx(i) - tlb(i)  )
           tau(i,0,k) = max( tau(i,0,k), tauoro(i) ) 
        endif
     end do   
     end do
  else
  ! Following is for backwards B4B compatibility with earlier versions
  ! ("N1" and "N5" -- Note: "N5" used do_backward_compat=.true.)
     if (.not.do_backward_compat) then
        do i = 1, ncol
        do k=1,pver+1
           if ( ( zi(i,k) <  wbr(i) ) .and. ( zi(i,k) >= tlb(i) ) ) then
              tau(i,0,k) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                            ( wbr(i) - zi(i,k) ) / &
                            ( wbr(i) - tlb(i)  )
           endif
        end do   
        end do
     else
        do i = 1, ncol
        do k=1,pver+1
           if ( ( zi(i,k) <= wbr(i) ) .and. ( zi(i,k) > tlb(i) ) ) then
              tau(i,0,k) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                            ( wbr(i) - zi(i,k) ) / &
                            ( wbr(i) - tlb(i)  )
           endif
        end do   
        end do
     end if
  end if
  
  if (lldo_trapped_waves) then 
     
  ! Identify top edge of layer in which Scorer param drops below 0
  ! - approximately the "turning level"
  !----------------------------------------------------------
     trn(:)=1.e8
     trn_level(:) = 0 ! pver+1
     where( m2(:,pver)<= 0. )
         trn(:) = zi(:,pver)
         trn_level(:) = pver
     endwhere
     do k=pver-1,1,-1
        where( (m2(:,k+1)> 0.).and.(m2(:,k)<= 0.) )
           trn(:) = zi(:,k)
           trn_level(:) = k
        endwhere
     end do

     do i = 1,ncol
     ! Case: Turning below mountain top
        if ( (trn(i) < mxdis(i)).and.(trn_level(i)>=1) ) then
            tau(i,0,:) =  tau(i,0,:) - max( tauoro(i),taudsw(i) )
            tau(i,0,:) =  max( tau(i,0,:) , 0. )
            tau(i,0,1:tlb_level(i))=0.
            src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning but no breaking
        if ( (wbr(i) == 0. ).and.(trn(i)>mxdis(i)).and.(trn_level(i)>=1) ) then
           where ( ( zi(i,:) <= trn(i) ) .and. ( zi(i,:) >= bwv(i) ) )
               tau(i,0,:) =  tauoro(i) * &
                             ( trn(i) - zi(i,:) ) / &
                             ( trn(i) - bwv(i)  )
           end where
           src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning AND breaking. Turning ABOVE breaking
        if ( (wbr(i) > 0. ).and.(trn(i) >= wbr(i)).and.(trn_level(i)>=1) ) then
           where ( ( zi(i,:) <= trn(i) ) .and. ( zi(i,:) >= wbr(i) ) )
               tau(i,0,:) =   tauoro(i) * &
                             ( trn(i) - zi(i,:) ) / &
                             ( trn(i) - wbr(i)  )
           endwhere
           src_level(i) = 1 ! disable any more tau calculation
        end if
        ! Case: Turning AND breaking. Turning BELOW breaking
        if ( (wbr(i) > 0. ).and.(trn(i) < wbr(i)).and.(trn_level(i)>=1) ) then
           tauoro(i) = 0.
           where ( ( zi(i,:) < wbr(i) ) .and. ( zi(i,:) >= tlb(i) ) )
               tau(i,0,:) =  tauoro(i) + (taudsw(i)-tauoro(i)) * &
                             ( wbr(i) - zi(i,:) ) / &
                             ( wbr(i) - tlb(i)  )
           endwhere
           src_level(i) = 1 ! disable any more tau calculation
        end if
     end do
  end if

end subroutine gw_rdg_break_trap


!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine endrun(msg)

   integer :: iulog

   character(len=*), intent(in), optional :: msg    ! string to be printed

    iulog=6

   if (present (msg)) then
      write(iulog,*)'ENDRUN:', msg
   else
      write(iulog,*)'ENDRUN: called without a message string'
   end if

   stop

end subroutine endrun


end module gw_rdg
