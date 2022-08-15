module gw_common

use gw_utils, only: GW_PRC, midpoint_interp
!
! This module contains code common to different gravity wave
! parameterizations.
!
implicit none
private
save

! Public interface.

public :: GWBand

public :: gw_common_init
public :: gw_newtonian_set
public :: gw_prof
public :: gw_drag_prof
public :: qbo_hdepth_scaling
public :: calc_taucd, momentum_flux, momentum_fixer
public :: energy_momentum_adjust, energy_change, energy_fixer
public :: hr_cf

public :: west, east, north, south
public :: pi
public :: gravit
public :: rair
public :: cpair

! Index the cardinal directions.
integer, parameter :: west = 1
integer, parameter :: east = 2
integer, parameter :: south = 3
integer, parameter :: north = 4

!++jtb (03/2020)
! Some physical constants used by GW codes
!-----------------------------------------
! 3.14159...
real(GW_PRC), parameter :: pi = acos(-1._GW_PRC)
! Acceleration due to gravity.
real(GW_PRC), protected :: gravit = huge(1._GW_PRC)
! Gas constant for dry air.
real(GW_PRC), protected :: rair = huge(1._GW_PRC)
! Specific heat for dry air.
real(GW_PRC), protected :: cpair = huge(1._GW_PRC)
! rair/gravit
real(GW_PRC) :: rog = huge(1._GW_PRC)


! Scaling factor for generating QBO
real(GW_PRC), protected :: qbo_hdepth_scaling
! Whether or not to enforce an upper boundary condition of tau = 0.
logical :: tau_0_ubc = .false.
! Inverse Prandtl number.
real(GW_PRC) :: prndl
! Heating rate conversion factor
real(GW_PRC), protected :: hr_cf


!
! Private variables
!

! Interface levels for gravity wave sources.
integer :: ktop = huge(1)

! Background diffusivity.
real(GW_PRC), parameter :: dback = 0.05_GW_PRC


! Newtonian cooling coefficients.
real, allocatable :: alpha(:)

!
! Limits to keep values reasonable.
!

! Minimum non-zero stress.
real(GW_PRC), parameter :: taumin = 1.e-10_GW_PRC
! Maximum wind tendency from stress divergence (before efficiency applied).
! 400 m/s/day
real(GW_PRC), parameter :: tndmax = 400._GW_PRC / 86400._GW_PRC
! Maximum allowed change in u-c (before efficiency applied).
real(GW_PRC), parameter :: umcfac = 0.5_GW_PRC
! Minimum value of (u-c)**2.
real(GW_PRC), parameter :: ubmc2mn = 0.01_GW_PRC

! Type describing a band of wavelengths into which gravity waves can be
! emitted.
! Currently this has to have uniform spacing (i.e. adjacent elements of
! cref are exactly dc apart).
type :: GWBand
   ! Dimension of the spectrum.
   integer :: ngwv
   ! Delta between nearest phase speeds [m/s].
   real :: dc
   ! Reference speeds [m/s].
   real, allocatable :: cref(:)
   ! Critical Froude number, squared (usually 1, but CAM3 used 0.5).
   real :: fcrit2
   ! Horizontal wave number [1/m].
   real :: kwv
   ! Effective horizontal wave number [1/m] (fcrit2*kwv).
   real :: effkwv
end type GWBand

interface GWBand
   module procedure new_GWBand
end interface

contains

!==========================================================================

! Constructor for a GWBand that calculates derived components.
function new_GWBand(ngwv, dc, fcrit2, wavelength) result(band)
  ! Used directly to set the type's components.
  integer, intent(in) :: ngwv
  real, intent(in) :: dc
  real, intent(in) :: fcrit2
  ! Wavelength in meters.
  real, intent(in) :: wavelength

  ! Output.
  type(GWBand) :: band

  ! Wavenumber index.
  integer :: l

  ! Simple assignments.
  band%ngwv = ngwv
  band%dc = dc
  band%fcrit2 = fcrit2

  ! Uniform phase speed reference grid.
  allocate(band%cref(-ngwv:ngwv))
  band%cref = [( dc * l, l = -ngwv, ngwv )]

  ! Wavenumber and effective wavenumber come from the wavelength.
  band%kwv = 2._GW_PRC*pi / wavelength
  band%effkwv = band%fcrit2 * band%kwv

end function new_GWBand

!==========================================================================

subroutine gw_common_init(   &
     tau_0_ubc_in, ktop_in, gravit_in, rair_in, cpair_in, & 
     prndl_in, qbo_hdepth_scaling_in, hr_cf_in, errstring)

  logical,  intent(in) :: tau_0_ubc_in
  integer,  intent(in) :: ktop_in
  real, intent(in) :: gravit_in
  real, intent(in) :: rair_in       ! Gas constant for dry air (J kg-1 K-1)
  real, intent(in) :: cpair_in      ! Heat cap. for dry air (J kg-1 K-1)
  real, intent(in) :: prndl_in
  real, intent(in) :: qbo_hdepth_scaling_in
  real, intent(in) :: hr_cf_in
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  integer :: pver
  integer :: ierr

  errstring = ""

  tau_0_ubc = tau_0_ubc_in
  ktop   = ktop_in
  gravit = gravit_in
  rair   = rair_in
  cpair  = cpair_in
  prndl  = prndl_in
  qbo_hdepth_scaling = qbo_hdepth_scaling_in
  hr_cf = hr_cf_in

  rog = rair/gravit

end subroutine gw_common_init

!==========================================================================

!!subroutine gw_prof (ncol, pver, p, t, rhoi, nm, ni)
subroutine gw_prof (ncol, pver, pint, pmid , t, rhoi, nm, ni)
  !-----------------------------------------------------------------------
  ! Compute profiles of background state quantities for the multiple
  ! gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------
  !------------------------------Arguments--------------------------------
  ! Column and vertical dimensions.
  integer, intent(in) :: ncol,pver
  ! Pressure coordinates.
  real, intent(in) :: pmid(ncol,pver) 
  real, intent(in) :: pint(ncol,pver+1) 

  ! Midpoint temperatures.
  real, intent(in) :: t(ncol,pver)

  ! Interface density.
  real, intent(out) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real, intent(out) :: nm(ncol,pver), ni(ncol,pver+1)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i,k

  ! dt/dp
  real(GW_PRC) :: dtdp
  ! Brunt-Vaisalla frequency squared.
  real(GW_PRC) :: n2

  ! Interface temperature.
  real(GW_PRC) :: ti(ncol,pver+1)

  ! Minimum value of Brunt-Vaisalla frequency squared.
  real(GW_PRC), parameter :: n2min = 5.e-5_GW_PRC

  !------------------------------------------------------------------------
  ! Determine the interface densities and Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------



  ! The top interface values are calculated assuming an isothermal
  ! atmosphere above the top level.
  k = 1
  do i = 1, ncol
     ti(i,k) = t(i,k)
     rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
     ni(i,k) = sqrt(gravit*gravit / (cpair*ti(i,k)))
  end do

  ! Interior points use centered differences.
  ti(:,2:pver) = midpoint_interp(real(t,GW_PRC))
  do k = 2, pver
     do i = 1, ncol
        rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
        dtdp = (t(i,k)-t(i,k-1)) / ( pmid(i,k)-pmid(i,k-1) ) ! * p%rdst(i,k-1)
        n2 = gravit*gravit/ti(i,k) * (1._GW_PRC/cpair - rhoi(i,k)*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     end do
  end do

  ! Bottom interface uses bottom level temperature, density; next interface
  ! B-V frequency.
  k = pver+1
  do i = 1, ncol
     ti(i,k) = t(i,k-1)
     rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
     ni(i,k) = ni(i,k-1)
  end do

  !------------------------------------------------------------------------
  ! Determine the midpoint Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------
  nm = midpoint_interp(ni)

end subroutine gw_prof

!==========================================================================
subroutine gw_drag_prof(ncol, pver, band, pint, delp, rdelp, & 
     src_level, tend_level, dt, t,    &
     piln, rhoi,    nm,   ni, ubm,  ubi,  xv,    yv,   &
     effgw,      c, kvtt, tau,  utgw,  vtgw, &
     ttgw,  gwut, ro_adjust, tau_adjust, &
     kwvrdg, satfac_in, lapply_effgw_in )

  !-----------------------------------------------------------------------
  ! Solve for the drag profile from the multiple gravity wave drag
  ! parameterization.
  ! 1. scan up from the wave source to determine the stress profile
  ! 2. scan down the stress profile to determine the tendencies
  !     => apply bounds to the tendency
  !          a. from wkb solution
  !          b. from computational stability constraints
  !     => adjust stress on interface below to reflect actual bounded
  !        tendency
  !-----------------------------------------------------------------------

  use linear_1d_operators, only: TriDiagDecomp

  !------------------------------Arguments--------------------------------
  ! Column and vertical dimension.
  integer, intent(in) :: ncol,pver
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  ! Interface pressures.
  real, intent(in) :: pint(ncol,pver+1)
  ! Delta Interface pressures.
  real, intent(in) :: delp(ncol,pver)
  ! Inverse of Delta Interface pressures.
  real, intent(in) :: rdelp(ncol,pver)
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Time step.
  real, intent(in) :: dt

  ! Midpoint and interface temperatures.
  real, intent(in) :: t(ncol,pver)
  ! Log of interface pressures.
  real, intent(in) :: piln(ncol,pver+1)
  ! Interface densities.
  real, intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real, intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real, intent(in) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real, intent(in) :: effgw(ncol)
  ! Wave phase speeds for each column.
  real(GW_PRC), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Molecular thermal diffusivity.
  real, intent(in) :: kvtt(ncol,pver+1)

!++jtb 
! remove q and dse for now (3/26/20)
  ! Constituent array.
  !real(GW_PRC), intent(in) :: q(:,:,:)
  ! Dry static energy.
  !real(GW_PRC), intent(in) :: dse(ncol,pver)
!--jtb

  ! Wave Reynolds stress.
  real(GW_PRC), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Zonal/meridional wind tendencies.
  real, intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real, intent(out) :: ttgw(ncol,pver)
!++jtb 3/2020
  ! Gravity wave constituent tendency.
  !real(GW_PRC), intent(out) :: qtgw(:,:,:)
!--jtb

  ! Gravity wave wind tendency for each wave.
  real(GW_PRC), intent(out) :: gwut(ncol,pver,-band%ngwv:band%ngwv)

  ! Adjustment parameter for IGWs.
  real, intent(in), optional :: &
       ro_adjust(ncol,-band%ngwv:band%ngwv,pver+1)

  ! Adjustment parameter for TAU.
  real, intent(in), optional :: &
       tau_adjust(ncol,pver+1)

  ! Diagnosed horizontal wavenumber for ridges.
  real, intent(in), optional :: &
       kwvrdg(ncol)

  ! Factor for saturation calculation. Here backwards 
  ! compatibility. I believe it should be 1.0 (jtb). 
  ! Looks like it has been 2.0 for a while in CAM.
  real, intent(in), optional :: &
       satfac_in

  logical, intent(in), optional :: lapply_effgw_in

  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  integer :: k, l, m, i, kmax

  ! Lowest tendency and source levels.
  integer :: kbot_tend, kbot_src

  ! "Total" and saturation diffusivity.
  real(GW_PRC) :: d(ncol)
  ! Imaginary part of vertical wavenumber.
  real(GW_PRC) :: mi(ncol)
  ! Stress after damping.
  real(GW_PRC) :: taudmp(ncol)
  ! Saturation stress.
  real(GW_PRC) :: tausat(ncol)
  ! (ub-c) and (ub-c)**2
  real(GW_PRC) :: ubmc(ncol), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, and at wave l).
  real(GW_PRC) :: ubt(ncol,pver), ubtl(ncol)
  real(GW_PRC) :: wrk(ncol)
  ! Ratio used for ubt tndmax limiting.
  real(GW_PRC) :: ubt_lim_ratio(ncol)
  ! Temporary effkwv
  real(GW_PRC) :: effkwv(ncol)

  ! saturation factor. Defaults to 2.0
  ! unless overidden by satfac_in
  real(GW_PRC) :: satfac

  real(GW_PRC) :: near_zero = tiny(1.0_GW_PRC)

  logical :: lapply_effgw

  ! LU decomposition.
  type(TriDiagDecomp) :: decomp

  if (present(satfac_in)) then
     satfac = satfac_in
  else
     satfac = 2.0
  endif

  ! Default behavior is to apply effgw and
  ! tendency limiters as designed by Sean
  ! Santos (lapply_effgw=.TRUE.). However,
  ! WACCM non-oro GW need to be retuned before
  ! this can done to them. --jtb 03/02/16
  if (present(lapply_effgw_in)) then
      lapply_effgw = lapply_effgw_in
  else
      lapply_effgw = .TRUE.
  endif

  
  ! Lowest levels that loops need to iterate over.
  kbot_tend = maxval(tend_level)
  kbot_src = maxval(src_level)

  ! Initialize gravity wave drag tendencies to zero.

  utgw    = 0.0
  vtgw    = 0.0
  gwut    = 0.0
  ttgw    = 0.0

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi     = 0.0
  taudmp = 0.0
  tausat = 0.0
  ubmc   = 0.0
  ubmc2  = 0.0
  wrk    = 0.0

  if (present(kwvrdg)) then
    effkwv = kwvrdg
  else
    effkwv = band%effkwv
  endif

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
!$OMP parallel do default(none) shared(kbot_src,ktop,kvtt,band,ubi,c,effkwv,rhoi,ni,satfac, &
!$OMP                                  ro_adjust,ncol,alpha,piln,t,rog,src_level,tau_adjust,tau) &
!$OMP                          private(k,d,l,i,tausat,taudmp,ubmc,ubmc2,wrk,mi)
  do k = kbot_src, ktop, -1  !++ but this is in model now 
     
     ! Determine the diffusivity for each column.

     d = dback + kvtt(:,k)

     do l = -band%ngwv, band%ngwv

        do i=1,ncol

        tausat(i) = 0.0
        taudmp(i) = 0.0

        if (src_level(i) >= k) then

          ! Determine the absolute value of the saturation stress.
          ! Define critical levels where the sign of (u-c) changes between
          ! interfaces.
          ubmc(i) = ubi(i,k) - c(i,l)

          ! Test to see if u-c has the same sign here as the level below.
          if (ubmc(i) > 0.0 .eqv. ubi(i,k+1) > c(i,l)) then
             if (ni(i,k) /= 0.0) & 
                 tausat(i) = abs( effkwv(i) * rhoi(i,k) * ubmc(i)**3 / &
                                  (satfac*ni(i,k)) )
             if (present(ro_adjust)) &
                 tausat(i) = tausat(i) * sqrt(ro_adjust(i,l,k))
             if (present(tau_adjust)) &
                 tausat(i) = tausat(i) * tau_adjust(i,k)
          endif

          ! Compute stress for each wave. The stress at this level is the
          ! min of the saturation stress and the stress at the level below
          ! reduced by damping. The sign of the stress must be the same as
          ! at the level below.
          ubmc2(i) = max(ubmc(i)**2, ubmc2mn)
          mi(i) = ni(i,k) / (2.0 * effkwv(i) * ubmc2(i)) * &  ! Is this 2.0 related to satfac?
                 (alpha(k) + ni(i,k)**2/ubmc2(i) * d(i))
          wrk(i) = -2.0*mi(i)*rog*t(i,k)*(piln(i,k+1) - piln(i,k))
          wrk(i) = max( wrk(i), -200.0 ) * exp(wrk(i))
          taudmp(i) = tau(i,l,k+1) * exp(wrk(i))
          ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
          ! we limit tau, so instead limit the arrays used to set it.
          if (tausat(i) <= taumin) tausat(i) = 0.0
          if (taudmp(i) <= taumin) taudmp(i) = 0.0
          tau(i,l,k) = min(taudmp(i), tausat(i))

        endif

        end do
     end do
  end do
  
  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0.0

  ! Apply efficiency to completed stress profile.
  if (lapply_effgw) then
!$OMP parallel do default(none) shared(kbot_tend,ktop,band,tau,tend_level,ncol,effgw) &
!$OMP                          private(k,l,i)
    do k = ktop, kbot_tend+1
        do l = -band%ngwv, band%ngwv
           do i=1,ncol
            if (k-1 <= tend_level(i)) then
              tau(i,l,k) = tau(i,l,k) * effgw(i)
            end if
           end do
        end do
     end do
  end if

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------


  ! Loop over levels from top to bottom
!$OMP parallel do default(none) shared(kbot_tend,ktop,band,ncol,tau,delp,rdelp,c,ubm,dt,gravit,utgw,vtgw, &
!$OMP                                  gwut,ubt,xv,yv,lapply_effgw,ubt_lim_ratio,tend_level,near_zero) &
!$OMP                          private(k,l,i,ubtl)
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0

     do l = -band%ngwv, band%ngwv    ! loop over wave

        do i=1,ncol

          ! Determine the wind tendency, including excess stress carried down
          ! from above.
          ubtl(i) = gravit * (tau(i,l,k+1)-tau(i,l,k)) * rdelp(i,k) ! p%rdel(i,k)  !/1/D_pint


          ! Apply first tendency limit to maintain numerical stability.
          ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
          !    (u^n+1 = u^n + du/dt * dt)
          ! The limiter is somewhat stricter, so that we don't come anywhere
          ! near reversing c-u.
          ubtl(i) = min(ubtl(i), umcfac * abs(c(i,l)-ubm(i,k)) / dt)

          ! Note: Here the limiter is being applied to each component wave
          ! seperately; BEFORE adding spectrum (conv., frontal) and BEFORE 
          ! applying effgw_{} (all GW)
          if (.not. lapply_effgw) ubtl(i) = min(ubtl(i), tndmax)
        
          if (k <= tend_level(i)) then

           ! Save tendency for each wave (for later computation of kzz).
           ! sign function returns magnitude of ubtl with sign of c-ubm 
           ! Renders ubt/ubm check for mountain waves unecessary
           gwut(i,k,l) = sign(ubtl(i), c(i,l)-ubm(i,k))
           ubt(i,k) = ubt(i,k) + gwut(i,k,l)

          end if

        end do

     end do

     if (lapply_effgw) then
        ! Apply second tendency limit to maintain numerical stability.
        ! Enforce du/dt < tndmax so that ridicuously large tendencies are not
        ! permitted.
        ! This can only happen above tend_level, so don't bother checking the
        ! level explicitly.
        do i=1,ncol
          if (abs(ubt(i,k)) > tndmax) then
           ubt_lim_ratio(i) = tndmax/abs(ubt(i,k))
           ubt(i,k) = ubt_lim_ratio(i) * ubt(i,k)
          else
           ubt_lim_ratio(i) = 1.0
          end if
        end do
     else
        ubt_lim_ratio = 1.0
     end if

     
     do l = -band%ngwv, band%ngwv
        gwut(:,k,l) = ubt_lim_ratio*gwut(:,k,l)
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.
        ! Include a protection on SMALL gwut to prevent floating point
        ! issues.
        !--------------------------------------------------
        do i=1,ncol
         if ( abs(gwut(i,k,l)) < near_zero ) then
           gwut(i,k,l) = 0.0
         end if   
         if (k <= tend_level(i)) then
           tau(i,l,k+1) = tau(i,l,k) + & 
                abs(gwut(i,k,l)) * delp(i,k) / gravit 
                !!! abs(gwut(i,k,l)) * p%del(i,k) / gravit 
         end if
        end do
     end do

     ! Project the mean wind tendency onto the components.
     do i=1,ncol
       if (k <= tend_level(i)) then
        utgw(i,k) = ubt(i,k) * xv(i)
        vtgw(i,k) = ubt(i,k) * yv(i)
       end if
     end do

     ! End of level loop.
  end do


  ! Block to undo Sean Santos mods to effgw and limiters.
  ! Here because non-oro GW in WACCM need extensive re-tuning
  ! before Sean's mods can be adopted. --jtb 03/02/16
  !==========================================
  if (.not.(lapply_effgw)) then
!$OMP parallel do default(none) shared(kbot_tend,ktop,band,ncol,tend_level,tau,effgw) &
!$OMP                          private(k,l)
     do k = ktop, kbot_tend+1
        do l = -band%ngwv, band%ngwv
           do i=1,ncol
             if (k-1 <= tend_level(i)) then
              tau(i,l,k) = tau(i,l,k) * effgw(i)
             end if
           end do
        end do
     end do
!$OMP parallel do default(none) shared(kbot_tend,ktop,band,gwut,utgw,vtgw,effgw) &
!$OMP                          private(k,l)
     do k = ktop, kbot_tend
        do l = -band%ngwv, band%ngwv
           gwut(:,k,l) = gwut(:,k,l) * effgw
        end do
        utgw(:,k) = utgw(:,k) * effgw
        vtgw(:,k) = vtgw(:,k) * effgw
     end do
  end if
  !===========================================

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
!$OMP parallel do default(none) shared(kbot_tend,ktop,band,ttgw,ubm,c,gwut) &
!$OMP                          private(k,l)
  do k = ktop, kbot_tend
     do l = -band%ngwv, band%ngwv
        ttgw(:,k) = ttgw(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
     end do
  end do
  ttgw = ttgw / cpair

end subroutine gw_drag_prof


!==========================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gw_newtonian_set( pver, pref )

use interpolate_data, only: lininterp

  integer,  intent(in)  :: pver
  real, intent(in)  :: pref( pver+1 )
  
  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer :: k
  integer, parameter :: nalph = 71

  real :: alpha0(nalph) = [ &
       0.1,         0.1,         0.1,         0.1,         &    
       0.1,         0.1,         0.1,         0.1,         &    
       0.1,         0.1,         0.10133333,  0.104,       &
       0.108,       0.112,       0.116,       0.12066667,  &    
       0.126,       0.132,       0.138,       0.144,       &    
       0.15133333,  0.16,        0.17,        0.18,        &    
       0.19,        0.19933333,  0.208,       0.216,       &    
       0.224,       0.232,       0.23466667,  0.232,       &
       0.224,       0.216,       0.208,       0.20133333,  &    
       0.196,       0.192,       0.188,       0.184,       &    
       0.18266667,  0.184,       0.188,       0.192,       &    
       0.196,       0.19333333,  0.184,       0.168,       &    
       0.152,       0.136,       0.12133333,  0.108,       &
       0.096,       0.084,       0.072,       0.061,       &    
       0.051,       0.042,       0.033,       0.024,       &    
       0.017666667, 0.014,       0.013,       0.012,       &    
       0.011,       0.010333333, 0.01,        0.01,        &    
       0.01,        0.01,        0.01                         &                 
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real :: palph(nalph) = [ &
       2.06115E-06, 2.74280E-06, 3.64988E-06, 4.85694E-06, &
       6.46319E-06, 8.60065E-06, 1.14450E-05, 1.52300E-05, &
       2.02667E-05, 2.69692E-05, 3.58882E-05, 4.77568E-05, &
       6.35507E-05, 8.45676E-05, 0.000112535, 0.000149752, &
       0.000199277, 0.000265180, 0.000352878, 0.000469579, &
       0.000624875, 0.000831529, 0.00110653,  0.00147247,  &
       0.00195943,  0.00260744,  0.00346975,  0.00461724,  &
       0.00614421,  0.00817618,  0.0108801,   0.0144783,   &
       0.0192665,   0.0256382,   0.0341170,   0.0453999,   &
       0.0604142,   0.0803939,   0.106981,    0.142361,    &    
       0.189442,    0.252093,    0.335463,    0.446404,    &    
       0.594036,    0.790490,    1.05192,     1.39980,     &    
       1.86273,     2.47875,     3.29851,     4.38936,     &    
       5.84098,     7.77266,     10.3432,     13.7638,     &    
       18.3156,     24.3728,     32.4332,     43.1593,     &    
       57.4326,     76.4263,     101.701,     135.335,     &    
       180.092,     239.651,     318.907,     424.373,     &    
       564.718,     751.477,     1000.                        &                 
       ]

  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa
  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400.0
     alpha0(k) = max(alpha0(k), 1.e-6)
     palph(k) = palph(k)*1.e2
  end do

  allocate (alpha(pver+1))
  ! interpolate to current vertical grid and obtain alpha
  call lininterp (alpha0, palph, nalph , alpha, pref, pver+1)

end subroutine gw_newtonian_set

!==========================================================================

! Calculate Reynolds stress for waves propagating in each cardinal
! direction.

function calc_taucd(ncol, pver, ngwv, tend_level, tau, c, xv, yv, ubi) &
     result(taucd)

  ! Column and gravity wave wavenumber dimensions.
  integer, intent(in) :: ncol, pver, ngwv
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(:)
  ! Wave Reynolds stress.
  real(GW_PRC), intent(in) :: tau(:,-ngwv:,:)
  ! Wave phase speeds for each column.
  real(GW_PRC), intent(in) :: c(:,-ngwv:)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(in) :: xv(:), yv(:)
  ! Projection of wind at interfaces.
  real, intent(in) :: ubi(:,:)

  real :: taucd(ncol,pver+1,4)

  ! Indices.
  integer :: i, k, l

  ! ubi at tend_level.
  real :: ubi_tend(ncol)

  ! Signed wave Reynolds stress.
  real :: tausg(ncol)

  ! Reynolds stress for waves propagating behind and forward of the wind.
  real :: taub(ncol)
  real :: tauf(ncol)

  taucd = 0.
  tausg = 0.

  ubi_tend = (/ (ubi(i,tend_level(i)+1), i = 1, ncol) /)

  do k = ktop, maxval(tend_level)+1

     taub = 0.
     tauf = 0.

     do l = -ngwv, ngwv
        where (k-1 <= tend_level)

           tausg = sign(tau(:,l,k), c(:,l)-ubi(:,k))

           where ( c(:,l) < ubi_tend )
              taub = taub + tausg
           elsewhere
              tauf = tauf + tausg
           end where

        end where
     end do

     where (k-1 <= tend_level)
        where (xv > 0.)
           taucd(:,k,east) = tauf * xv
           taucd(:,k,west) = taub * xv
        elsewhere
           taucd(:,k,east) = taub * xv
           taucd(:,k,west) = tauf * xv
        end where

        where ( yv > 0.)
           taucd(:,k,north) = tauf * yv
           taucd(:,k,south) = taub * yv
        elsewhere
           taucd(:,k,north) = taub * yv
           taucd(:,k,south) = tauf * yv
        end where
     end where

  end do

end function calc_taucd

!==========================================================================

! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Projected stresses.
  real, intent(in) :: taucd(:,:,:)
  ! Components of momentum change sourced from the bottom.
  real, intent(out) :: um_flux(:), vm_flux(:)

  integer :: i

  ! Tendency for U & V below source level.
  do i = 1, size(tend_level)
     um_flux(i) = taucd(i,tend_level(i)+1, east) + &
                  taucd(i,tend_level(i)+1, west)
     vm_flux(i) = taucd(i,tend_level(i)+1,north) + &
                  taucd(i,tend_level(i)+1,south)
  end do

end subroutine momentum_flux

!==========================================================================

! Subtracts a change in momentum in the gravity wave levels from wind
! tendencies in lower levels, ensuring momentum conservation.
subroutine momentum_fixer(ncol, pver, tend_level, p, um_flux, vm_flux, utgw, vtgw)

  integer, intent(in) :: ncol, pver
  ! Bottom stress level.
  integer, intent(in) :: tend_level(ncol)
  ! Pressure coordinates.
  real, intent(in) :: p(ncol,pver+1)
  ! Components of momentum change sourced from the bottom.
  real, intent(in) :: um_flux(ncol), vm_flux(ncol)
  ! Wind tendencies.
  real, intent(inout) :: utgw(ncol,pver), vtgw(ncol,pver)

  ! Indices.
  integer :: i, k
  ! Reciprocal of total mass.
  real :: rdm(ncol)

  ! Total mass from ground to source level: rho*dz = dp/gravit
  do i = 1, ncol
     rdm(i) = gravit/(p(i,pver+1)-p(i,tend_level(i)+1))
  end do

  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        utgw(:,k) = utgw(:,k) + -um_flux*rdm
        vtgw(:,k) = vtgw(:,k) + -vm_flux*rdm
     end where
  end do
  
end subroutine momentum_fixer

!==========================================================================

! Calculate the change in total energy from tendencies up to this point.
subroutine energy_change(ncol,pver, dt, pdel, u, v, dudt, dvdt, dsdt, de)

  integer, intent(in) :: ncol, pver
  ! Time step.
  real, intent(in) :: dt
  ! Pressure coordinates.
  real, intent(in) :: pdel(ncol,pver+1)
  ! Winds at start of time step.
  real, intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Wind tendencies.
  real, intent(in) :: dudt(ncol,pver), dvdt(ncol,pver)
  ! Heating tendency.
  real, intent(in) :: dsdt(ncol,pver)
  ! Change in energy.
  real, intent(out) :: de(ncol)

  ! Level index.
  integer :: k

  ! Net gain/loss of total energy in the column.
  de = 0.0
  do k = 1, pver
     de = de + pdel(:,k)/gravit * (dsdt(:,k) + &
          dudt(:,k)*(u(:,k)+dudt(:,k)*0.5*dt) + &
          dvdt(:,k)*(v(:,k)+dvdt(:,k)*0.5*dt) )
  end do

end subroutine energy_change

!==========================================================================

! Subtract change in energy from the heating tendency in the levels below
! the gravity wave region.
subroutine energy_fixer(ncol, pver, tend_level, pint, de, ttgw)

  integer, intent(in) :: ncol, pver
  ! Bottom stress level.
  integer, intent(in) :: tend_level(ncol)
  ! Pressure coordinates.
  real, intent(in) :: pint(ncol,pver+1)
  ! Change in energy.
  real, intent(in) :: de(ncol)
  ! Heating tendency.
  real, intent(inout) :: ttgw(ncol,pver)

  ! Column/level indices.
  integer :: i, k
  ! Energy change to apply divided by all the mass it is spread across.
  real :: de_dm(ncol)

  do i = 1, ncol
     de_dm(i) = -de(i)*gravit/(pint(i,pver+1)-pint(i,tend_level(i)+1))
  end do

  ! Subtract net gain/loss of total energy below tend_level.
  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        ttgw(:,k) = ttgw(:,k) + de_dm
     end where
  end do

end subroutine energy_fixer

!==========================================================================

subroutine energy_momentum_adjust(ncol, pver, kbot, band, pint, delp, c, tau, &
                               effgw, t, ubm, ubi, xv, yv, utgw, vtgw, ttgw)

  integer, intent(in) :: ncol, pver, kbot
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure interfaces.
  real, intent(in) :: pint(ncol,pver+1)
  ! Pressure thickness.
  real, intent(in) :: delp(ncol,pver)
  ! Wave phase speeds for each column.
  real(GW_PRC), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Wave Reynolds stress.
  real(GW_PRC), intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Efficiency
  real, intent(in) :: effgw(ncol)
  ! Temperature and winds
  real, intent(in) :: t(ncol,pver), ubm(ncol,pver), ubi(ncol,pver)
  ! projected winds.
  real, intent(in) :: xv(ncol), yv(ncol)
  ! Tendencies.
  real, intent(inout) :: utgw(ncol,pver)
  real, intent(inout) :: vtgw(ncol,pver)
  real, intent(inout) :: ttgw(ncol,pver)

  real :: zlb,pm,rhom,cmu,fpmx,fpmy,fe,fpe,fpml,fpel,fpmt,fpet,dusrcl,dvsrcl,dtsrcl
  real :: utfac,uhtmax
  integer  :: ktop

  ! Level index.
  integer :: i,k,l

! GEOS efficiency and energy/momentum adjustments
  ktop=1
  do i=1,ncol
! Calculate launch level height
    zlb = 0.
    do k = ktop+1, pver
       if (k >= kbot+1) then
! Define layer pressure and density
          pm   = (pint(i,k-1)+pint(i,k))*0.5
          rhom = pm/(rair*t(i,k))
          zlb  = zlb + delp(i,k)/rair/rhom
       end if
    end do
   !-----------------------------------------------------------------------
   ! Calculates energy and momentum flux profiles
   !-----------------------------------------------------------------------
    do l = -band%ngwv, band%ngwv
       do k = ktop, pver
          if ( k <= kbot ) then
             cmu  = c(i,l)-ubi(i,k)
             fpmx =        sign(1.0,cmu)*tau(i,l,k)*xv(i)*effgw(i)
             fpmy =        sign(1.0,cmu)*tau(i,l,k)*yv(i)*effgw(i)
             fe   =    cmu*sign(1.0,cmu)*tau(i,l,k)      *effgw(i)
             fpe  = c(i,l)*sign(1.0,cmu)*tau(i,l,k)      *effgw(i)
             if (k == kbot) then
                fpml = fpmx*xv(i)+fpmy*yv(i)
                fpel = fpe
             end if
             if (k == ktop) then
                fpmt = fpmx*xv(i)+fpmy*yv(i)
                fpet = fpe
             end if
          end if
          if (k >= kbot+1) then
! Define layer pressure and density
             pm   = (pint(i,k-1)+pint(i,k))*0.5
             rhom = pm/(rair*t(i,k))
! Compute sub-source tendencies
             dusrcl = - (fpml-fpmt)/(rhom*zlb)*xv(i)
             dvsrcl = - (fpml-fpmt)/(rhom*zlb)*yv(i)
             dtsrcl = -((fpel-fpet)-ubm(i,k)*(fpml-fpmt))/  &
                              (rhom*zlb*cpair)
! Add sub-source wind and temperature tendencies
             utgw(i,k) = utgw(i,k) + dusrcl
             vtgw(i,k) = vtgw(i,k) + dvsrcl
             ttgw(i,k) = ttgw(i,k) + dtsrcl
          end if
       end do
    end do
   !-----------------------------------------------------------------------
   ! Adjust efficiency factor to prevent unrealistically strong forcing
   !-----------------------------------------------------------------------
    uhtmax = 0.0
    utfac  = 1.0
    do k = ktop, pver
       uhtmax = max(sqrt(utgw(i,k)**2 + vtgw(i,k)**2), uhtmax)
    end do
    if (uhtmax > tndmax) utfac = tndmax/uhtmax
    do k = ktop, pver
       utgw(i,k) = utgw(i,k)*utfac
       vtgw(i,k) = vtgw(i,k)*utfac
       ttgw(i,k) = ttgw(i,k)*utfac
    end do

  end do  ! i=1,ncol

end subroutine energy_momentum_adjust 

!==========================================================================
end module gw_common
