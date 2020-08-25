module gw_common

!
! This module contains code common to different gravity wave
! parameterizations.
  !

!!use shr_kind_mod,   only: r8=>shr_kind_r8
!!use gw_utils, only: r8
use coords_1d, only: Coords1D

implicit none
private
save

! Public interface.

public :: GWBand

public :: gw_newtonian_set
public :: gw_common_init
public :: gw_prof
public :: gw_drag_prof
public :: qbo_hdepth_scaling

!++jtb 
!  These go away for now (3/26/20)
!public :: calc_taucd, momentum_flux, momentum_fixer
!public :: energy_change, energy_fixer
!public :: coriolis_speed, adjust_inertial
!--

!public :: pver

public :: west, east, north, south
public :: pi
public :: gravit
public :: rair
public :: cpair

integer,parameter :: r8 = selected_real_kind(12) ! 8 byte real

! Index the cardinal directions.
integer, parameter :: west = 1
integer, parameter :: east = 2
integer, parameter :: south = 3
integer, parameter :: north = 4


!++jtb (03/2020)
! Some physical constants used by GW codes
!-----------------------------------------
! 3.14159...
real(r8), parameter :: pi = acos(-1._r8)
! Acceleration due to gravity.
real(r8), protected :: gravit = huge(1._r8)
! Gas constant for dry air.
real(r8), protected :: rair = huge(1._r8)
! Specific heat for dry air.
real(r8), protected :: cpair = huge(1._r8)
! rair/gravit
real(r8) :: rog = huge(1._r8)


! Scaling factor for generating QBO
real(r8), protected :: qbo_hdepth_scaling
! Whether or not to enforce an upper boundary condition of tau = 0.
logical :: tau_0_ubc = .false.
! Inverse Prandtl number.
real(r8) :: prndl



!
! Private variables
!

! Interface levels for gravity wave sources.
integer :: ktop = huge(1)

! Background diffusivity.
real(r8), parameter :: dback = 0.05_r8


! Newtonian cooling coefficients.
real(r8), allocatable :: alpha(:)

!
! Limits to keep values reasonable.
!

! Minimum non-zero stress.
real(r8), parameter :: taumin = 1.e-10_r8
! Maximum wind tendency from stress divergence (before efficiency applied).
! 400 m/s/day
real(r8), parameter :: tndmax = 400._r8 / 86400._r8
! Maximum allowed change in u-c (before efficiency applied).
real(r8), parameter :: umcfac = 0.5_r8
! Minimum value of (u-c)**2.
real(r8), parameter :: ubmc2mn = 0.01_r8

! Type describing a band of wavelengths into which gravity waves can be
! emitted.
! Currently this has to have uniform spacing (i.e. adjacent elements of
! cref are exactly dc apart).
type :: GWBand
   ! Dimension of the spectrum.
   integer :: ngwv
   ! Delta between nearest phase speeds [m/s].
   real(r8) :: dc
   ! Reference speeds [m/s].
   real(r8), allocatable :: cref(:)
   ! Critical Froude number, squared (usually 1, but CAM3 used 0.5).
   real(r8) :: fcrit2
   ! Horizontal wave number [1/m].
   real(r8) :: kwv
   ! Effective horizontal wave number [1/m] (fcrit2*kwv).
   real(r8) :: effkwv
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
  real(r8), intent(in) :: dc
  real(r8), intent(in) :: fcrit2
  ! Wavelength in meters.
  real(r8), intent(in) :: wavelength

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
  band%kwv = 2._r8*pi / wavelength
  band%effkwv = band%fcrit2 * band%kwv

end function new_GWBand

!==========================================================================

subroutine gw_common_init(   &
     tau_0_ubc_in, ktop_in, gravit_in, rair_in, cpair_in, & 
     prndl_in, qbo_hdepth_scaling_in, errstring)

  logical,  intent(in) :: tau_0_ubc_in
  integer,  intent(in) :: ktop_in
  real(r8), intent(in) :: gravit_in
  real(r8), intent(in) :: rair_in       ! Gas constant for dry air (J kg-1 K-1)
  real(r8), intent(in) :: cpair_in      ! Heat cap. for dry air (J kg-1 K-1)
  real(r8), intent(in) :: prndl_in
  real(r8), intent(in) :: qbo_hdepth_scaling_in
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
  use gw_utils, only: midpoint_interp
  !------------------------------Arguments--------------------------------
  ! Column and vertical dimensions.
  integer, intent(in) :: ncol,pver
  ! Pressure coordinates.
  !type(Coords1D), intent(in) :: p
  real(r8), intent(in) :: pmid(ncol,pver) 
  real(r8), intent(in) :: pint(ncol,pver+1) 

  ! Specific heat of dry air, constant pressure.
  !real(r8), intent(in) :: cpair
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)

  ! Interface density.
  real(r8), intent(out) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(out) :: nm(ncol,pver), ni(ncol,pver+1)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i,k

  ! dt/dp
  real(r8) :: dtdp
  ! Brunt-Vaisalla frequency squared.
  real(r8) :: n2

  ! Interface temperature.
  real(r8) :: ti(ncol,pver+1)

  ! Minimum value of Brunt-Vaisalla frequency squared.
  real(r8), parameter :: n2min = 5.e-5_r8

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
  ti(:,2:pver) = midpoint_interp(t)
  do k = 2, pver
     do i = 1, ncol
        rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
        dtdp = (t(i,k)-t(i,k-1)) / ( pmid(i,k)-pmid(i,k-1) ) ! * p%rdst(i,k-1)
        n2 = gravit*gravit/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
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
     ttgw,  egwdffi,   gwut, dttdf, dttke, ro_adjust, &
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

  use gw_diffusion, only: gw_ediff, gw_diff_tend
  use linear_1d_operators, only: TriDiagDecomp

  !------------------------------Arguments--------------------------------
  ! Column and vertical dimension.
  integer, intent(in) :: ncol,pver
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  ! Interface pressures.
  real(r8), intent(in) :: pint(ncol,pver+1)
  ! Delta Interface pressures.
  real(r8), intent(in) :: delp(ncol,pver)
  ! Inverse of Delta Interface pressures.
  real(r8), intent(in) :: rdelp(ncol,pver)
      !type(Coords1D), intent(in) :: p
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Time step.
  real(r8), intent(in) :: dt

  ! Midpoint and interface temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Log of interface pressures.
  real(r8), intent(in) :: piln(ncol,pver+1)
  ! Interface densities.
  real(r8), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(in) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real(r8), intent(in) :: effgw(ncol)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Molecular thermal diffusivity.
  real(r8), intent(in) :: kvtt(ncol,pver+1)

!++jtb 
! remove q and dse for now (3/26/20)
  ! Constituent array.
  !real(r8), intent(in) :: q(:,:,:)
  ! Dry static energy.
  !real(r8), intent(in) :: dse(ncol,pver)
!--jtb

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Zonal/meridional wind tendencies.
  real(r8), intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real(r8), intent(out) :: ttgw(ncol,pver)
!++jtb 3/2020
  ! Gravity wave constituent tendency.
  !real(r8), intent(out) :: qtgw(:,:,:)
!--jtb

  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,pver+1)

  ! Gravity wave wind tendency for each wave.
  real(r8), intent(out) :: gwut(ncol,pver,-band%ngwv:band%ngwv)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8), intent(out) :: dttdf(ncol,pver)
  real(r8), intent(out) :: dttke(ncol,pver)

  ! Adjustment parameter for IGWs.
  real(r8), intent(in), optional :: &
       ro_adjust(ncol,-band%ngwv:band%ngwv,pver+1)

  ! Diagnosed horizontal wavenumber for ridges.
  real(r8), intent(in), optional :: &
       kwvrdg(ncol)

  ! Factor for saturation calculation. Here backwards 
  ! compatibility. I believe it should be 1.0 (jtb). 
  ! Looks like it has been 2.0 for a while in CAM.
  real(r8), intent(in), optional :: &
       satfac_in

  logical, intent(in), optional :: lapply_effgw_in

  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  integer :: k, l, m, i, kmax

  ! Lowest tendency and source levels.
  integer :: kbot_tend, kbot_src

  ! "Total" and saturation diffusivity.
  real(r8) :: d(ncol)
  ! Imaginary part of vertical wavenumber.
  real(r8) :: mi(ncol)
  ! Stress after damping.
  real(r8) :: taudmp(ncol)
  ! Saturation stress.
  real(r8) :: tausat(ncol)
  ! (ub-c) and (ub-c)**2
  real(r8) :: ubmc(ncol), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, and at wave l).
  real(r8) :: ubt(ncol,pver), ubtl(ncol)
  real(r8) :: wrk(ncol)
  ! Ratio used for ubt tndmax limiting.
  real(r8) :: ubt_lim_ratio(ncol)



  ! saturation factor. Defaults to 2.0
  ! unless overidden by satfac_in
  real(r8) :: satfac

  logical :: lapply_effgw

  ! LU decomposition.
  type(TriDiagDecomp) :: decomp


  allocate( alpha(pver+1) )
  alpha(:)=0._r8

  if (present(satfac_in)) then
     satfac = satfac_in
  else
     satfac = 2._r8
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

  utgw = 0._r8
  vtgw = 0._r8

  gwut = 0._r8

  dttke = 0._r8
  ttgw = 0._r8
  egwdffi=0._r8
  dttdf=0._r8

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi = 0._r8
  taudmp = 0._r8
  tausat = 0._r8
  ubmc = 0._r8
  ubmc2 = 0._r8
  wrk = 0._r8

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
  ! do k = kbot_src-1, ktop, -1 !++jtb I think this is right 
  do k = kbot_src, ktop, -1  !++ but this is in model now 
     
     ! Determine the diffusivity for each column.

     d = dback + kvtt(:,k)

     do l = -band%ngwv, band%ngwv

        ! Determine the absolute value of the saturation stress.
        ! Define critical levels where the sign of (u-c) changes between
        ! interfaces.
        ubmc = ubi(:,k) - c(:,l)

        tausat = 0.0_r8

        if (present(kwvrdg)) then
           where (src_level >= k)
              ! Test to see if u-c has the same sign here as the level below.
              where (ubmc > 0.0_r8 .eqv. ubi(:,k+1) > c(:,l))
                 tausat = abs(  kwvrdg  * rhoi(:,k) * ubmc**3 / &
                    (satfac*ni(:,k)))
              end where
           end where
        else
           where (src_level >= k)
              ! Test to see if u-c has the same sign here as the level below.
              where (ubmc > 0.0_r8 .eqv. ubi(:,k+1) > c(:,l))
                 tausat = abs(band%effkwv * rhoi(:,k) * ubmc**3 / &
                    (satfac*ni(:,k)))
              end where
           end where
        end if

        if (present(ro_adjust)) then
           where (src_level >= k)
              tausat = tausat * sqrt(ro_adjust(:,l,k))
           end where
        end if

        if (present(kwvrdg)) then
           where (src_level >= k)
              ! Compute stress for each wave. The stress at this level is the
              ! min of the saturation stress and the stress at the level below
              ! reduced by damping. The sign of the stress must be the same as
              ! at the level below.

              ubmc2 = max(ubmc**2, ubmc2mn)
              mi = ni(:,k) / (2._r8 *   kwvrdg * ubmc2) * &  ! Is this 2._r8 related to satfac?
                 (alpha(k) + ni(:,k)**2/ubmc2 * d)
              wrk = -2._r8*mi*rog*t(:,k)*(piln(:,k+1) - piln(:,k))
                wrk = max( wrk, -200._r8 )

              taudmp = tau(:,l,k+1)

              ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
              ! we limit tau, so instead limit the arrays used to set it.
              where (tausat <= taumin) tausat = 0._r8
              where (taudmp <= taumin) taudmp = 0._r8

              tau(:,l,k) = min(taudmp, tausat)
           end where

        else
           where (src_level >= k)

              ! Compute stress for each wave. The stress at this level is the
              ! min of the saturation stress and the stress at the level below
              ! reduced by damping. The sign of the stress must be the same as
              ! at the level below.
              ubmc2 = max(ubmc**2, ubmc2mn)
              mi = ni(:,k) / (2._r8 * band%kwv * ubmc2) * &
                 (alpha(k) + ni(:,k)**2/ubmc2 * d)
              wrk = -2._r8*mi*rog*t(:,k)*(piln(:,k+1) - piln(:,k))
                wrk = max( wrk, -200._r8 )

              taudmp = tau(:,l,k+1) * exp(wrk)

              ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
              ! we limit tau, so instead limit the arrays used to set it.
              where (tausat <= taumin) tausat = 0._r8
              where (taudmp <= taumin) taudmp = 0._r8

              tau(:,l,k) = min(taudmp, tausat)
           end where
        endif

     end do
  end do
  


  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0._r8

  ! Apply efficiency to completed stress profile.
  if (lapply_effgw) then
     do k = ktop, kbot_tend+1
        do l = -band%ngwv, band%ngwv
           where (k-1 <= tend_level)
              tau(:,l,k) = tau(:,l,k) * effgw
           end where
        end do
     end do
  end if

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------


  ! Loop over levels from top to bottom
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0_r8

     do l = -band%ngwv, band%ngwv    ! loop over wave

        ! Determine the wind tendency, including excess stress carried down
        ! from above.
        ubtl = gravit * (tau(:,l,k+1)-tau(:,l,k)) * rdelp(:,k) ! p%rdel(:,k)  !/1/D_pint


        ! Apply first tendency limit to maintain numerical stability.
        ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
        !    (u^n+1 = u^n + du/dt * dt)
        ! The limiter is somewhat stricter, so that we don't come anywhere
        ! near reversing c-u.
        ubtl = min(ubtl, umcfac * abs(c(:,l)-ubm(:,k)) / dt)

        ! Note: Here the limiter is being applied to each component wave
        ! seperately; BEFORE adding spectrum (conv., frontal) and BEFORE 
        ! applying effgw_{} (all GW)
        if (.not. lapply_effgw) ubtl = min(ubtl, tndmax)
        

        where (k <= tend_level)

           ! Save tendency for each wave (for later computation of kzz).
           ! sign function returns magnitude of ubtl with sign of c-ubm 
           ! Renders ubt/ubm check for mountain waves unecessary
           gwut(:,k,l) = sign(ubtl, c(:,l)-ubm(:,k))
           ubt(:,k) = ubt(:,k) + gwut(:,k,l)

        end where

     end do

     if (lapply_effgw) then
        ! Apply second tendency limit to maintain numerical stability.
        ! Enforce du/dt < tndmax so that ridicuously large tendencies are not
        ! permitted.
        ! This can only happen above tend_level, so don't bother checking the
        ! level explicitly.
        where (abs(ubt(:,k)) > tndmax)
           ubt_lim_ratio = tndmax/abs(ubt(:,k))
           ubt(:,k) = ubt_lim_ratio * ubt(:,k)
        elsewhere
           ubt_lim_ratio = 1._r8
        end where
     else
        ubt_lim_ratio = 1._r8
     end if

     
     do l = -band%ngwv, band%ngwv
        gwut(:,k,l) = ubt_lim_ratio*gwut(:,k,l)
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.
        where (k <= tend_level)
           tau(:,l,k+1) = tau(:,l,k) + & 
                abs(gwut(:,k,l)) * delp(:,k) / gravit 
                !!! abs(gwut(:,k,l)) * p%del(:,k) / gravit 
        end where
     end do

     ! Project the mean wind tendency onto the components.
     where (k <= tend_level)
        utgw(:,k) = ubt(:,k) * xv
        vtgw(:,k) = ubt(:,k) * yv
     end where

     ! End of level loop.
  end do


  ! Block to undo Sean Santos mods to effgw and limiters.
  ! Here because non-oro GW in WACCM need extensive re-tuning
  ! before Sean's mods can be adopted. --jtb 03/02/16
  !==========================================
  if (.not.(lapply_effgw)) then
     do k = ktop, kbot_tend+1
        do l = -band%ngwv, band%ngwv
           where (k-1 <= tend_level)
              tau(:,l,k) = tau(:,l,k) * effgw
           end where
        end do
     end do
     do k = ktop, kbot_tend
        do l = -band%ngwv, band%ngwv
           gwut(:,k,l) = gwut(:,k,l) * effgw
        end do
        utgw(:,k) = utgw(:,k) * effgw
        vtgw(:,k) = vtgw(:,k) * effgw
     end do
  end if
  !===========================================


!++jtb
! Let's disable vertical diffusion calculations in GW
! until we have better physical ideas. (JTB 3/26/20).
! Also avoids having to pass tracers into GW code
!=====================================
#if 0
  ! Calculate effective diffusivity and LU decomposition for the
  ! vertical diffusion solver.
  call gw_ediff (ncol, pver, band%ngwv, kbot_tend, ktop, tend_level, &
       gwut, ubm, nm, rhoi, dt, prndl, gravit, p, c, &
       egwdffi, decomp, ro_adjust=ro_adjust)

  ! Calculate tendency on each constituent.
  do m = 1, size(q,3)

     call gw_diff_tend(ncol, pver, kbot_tend, ktop, q(:,:,m), &
          dt, decomp, qtgw(:,:,m))

  enddo

  ! Calculate tendency from diffusing dry static energy (dttdf).
  call gw_diff_tend(ncol, pver, kbot_tend, ktop, dse, dt, decomp, dttdf)

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
  do l = -band%ngwv, band%ngwv
     do k = ktop, kbot_tend
        dttke(:,k) = dttke(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
     end do
  end do

  ttgw = dttke + dttdf

  ! Deallocate decomp.
  call decomp%finalize()
#else

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
  do l = -band%ngwv, band%ngwv
     do k = ktop, kbot_tend
        dttke(:,k) = dttke(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
     end do
  end do

  ttgw = dttke / cpair

#endif

deallocate( alpha )

end subroutine gw_drag_prof


!==========================================================================


!==========================================================================
! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Projected stresses.
  real(r8), intent(in) :: taucd(:,:,:)
  ! Components of momentum change sourced from the bottom.
  real(r8), intent(out) :: um_flux(:), vm_flux(:)

  integer :: i

  ! Tendency for U & V below source level.
  do i = 1, size(tend_level)
     um_flux(i) = taucd(i,tend_level(i)+1, east) + &
                  taucd(i,tend_level(i)+1, west)
     vm_flux(i) = taucd(i,tend_level(i)+1,north) + &
                  taucd(i,tend_level(i)+1,south)
  end do

end subroutine momentum_flux


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gw_newtonian_set( pver, pref, alpha )

use interpolate_data, only: lininterp

  integer,  intent(in)  :: pver
  real(r8), intent(in)  :: pref( pver+1 )
  
  ! Interpolated Newtonian cooling coefficients.
  real(r8), intent(out) :: alpha(pver+1)

  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer :: k
  integer, parameter :: nalph = 71
  real(r8) :: alpha0(nalph) = [ &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.1_r8,         0.1_r8,         &
       0.1_r8,         0.1_r8,         0.10133333_r8,  0.104_r8,       &
       0.108_r8,       0.112_r8,       0.116_r8,       0.12066667_r8,  &
       0.126_r8,       0.132_r8,       0.138_r8,       0.144_r8,       &
       0.15133333_r8,  0.16_r8,        0.17_r8,        0.18_r8,        &
       0.19_r8,        0.19933333_r8,  0.208_r8,       0.216_r8,       &
       0.224_r8,       0.232_r8,       0.23466667_r8,  0.232_r8,       &
       0.224_r8,       0.216_r8,       0.208_r8,       0.20133333_r8,  &
       0.196_r8,       0.192_r8,       0.188_r8,       0.184_r8,       &
       0.18266667_r8,  0.184_r8,       0.188_r8,       0.192_r8,       &
       0.196_r8,       0.19333333_r8,  0.184_r8,       0.168_r8,       &
       0.152_r8,       0.136_r8,       0.12133333_r8,  0.108_r8,       &
       0.096_r8,       0.084_r8,       0.072_r8,       0.061_r8,       &
       0.051_r8,       0.042_r8,       0.033_r8,       0.024_r8,       &
       0.017666667_r8, 0.014_r8,       0.013_r8,       0.012_r8,       &
       0.011_r8,       0.010333333_r8, 0.01_r8,        0.01_r8,        &
       0.01_r8,        0.01_r8,        0.01_r8                         &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(r8) :: palph(nalph) = [ &
       2.06115E-06_r8, 2.74280E-06_r8, 3.64988E-06_r8, 4.85694E-06_r8, &
       6.46319E-06_r8, 8.60065E-06_r8, 1.14450E-05_r8, 1.52300E-05_r8, &
       2.02667E-05_r8, 2.69692E-05_r8, 3.58882E-05_r8, 4.77568E-05_r8, &
       6.35507E-05_r8, 8.45676E-05_r8, 0.000112535_r8, 0.000149752_r8, &
       0.000199277_r8, 0.000265180_r8, 0.000352878_r8, 0.000469579_r8, &
       0.000624875_r8, 0.000831529_r8, 0.00110653_r8,  0.00147247_r8,  &
       0.00195943_r8,  0.00260744_r8,  0.00346975_r8,  0.00461724_r8,  &
       0.00614421_r8,  0.00817618_r8,  0.0108801_r8,   0.0144783_r8,   &
       0.0192665_r8,   0.0256382_r8,   0.0341170_r8,   0.0453999_r8,   &
       0.0604142_r8,   0.0803939_r8,   0.106981_r8,    0.142361_r8,    &
       0.189442_r8,    0.252093_r8,    0.335463_r8,    0.446404_r8,    &
       0.594036_r8,    0.790490_r8,    1.05192_r8,     1.39980_r8,     &
       1.86273_r8,     2.47875_r8,     3.29851_r8,     4.38936_r8,     &
       5.84098_r8,     7.77266_r8,     10.3432_r8,     13.7638_r8,     &
       18.3156_r8,     24.3728_r8,     32.4332_r8,     43.1593_r8,     &
       57.4326_r8,     76.4263_r8,     101.701_r8,     135.335_r8,     &
       180.092_r8,     239.651_r8,     318.907_r8,     424.373_r8,     &
       564.718_r8,     751.477_r8,     1000._r8                        &
       ]

  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa
  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._r8
     alpha0(k) = max(alpha0(k), 1.e-6_r8)
     palph(k) = palph(k)*1.e2_r8
  end do

  ! interpolate to current vertical grid and obtain alpha
  call lininterp (alpha0  ,palph, nalph , alpha  , pref , pver+1)


end subroutine gw_newtonian_set

end module gw_common
