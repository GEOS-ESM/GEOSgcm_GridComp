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
public :: gw_prof
public :: gw_drag_prof
public :: qbo_hdepth_scaling
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
real(GW_PRC), allocatable :: alpha(:)

!
! Limits to keep values reasonable.
!

! Minimum non-zero stress.
real(GW_PRC), parameter :: taumin = 1.e-10_GW_PRC
! Maximum wind tendency from stress divergence (before efficiency applied).
!! 400 m/s/day
!real(GW_PRC), parameter :: tndmax = 400._GW_PRC / 86400._GW_PRC
! 80 m/s/day
real(GW_PRC), parameter :: tndmax = 80._GW_PRC / 86400._GW_PRC
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
   real(GW_PRC) :: dc
   ! Reference speeds [m/s].
   real(GW_PRC), allocatable :: cref(:)
   ! Critical Froude number, squared (usually 1, but CAM3 used 0.5).
   real(GW_PRC) :: fcrit2
   ! Horizontal wave number [1/m].
   real(GW_PRC) :: kwv
   ! Effective horizontal wave number [1/m] (fcrit2*kwv).
   real(GW_PRC) :: effkwv
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

  ! Specific heat of dry air, constant pressure.
  !real(GW_PRC), intent(in) :: cpair
  ! Midpoint temperatures.
  real, intent(in) :: t(ncol,pver)

  ! Interface density.
  real(GW_PRC), intent(out) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(GW_PRC), intent(out) :: nm(ncol,pver), ni(ncol,pver+1)

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
     ttgw,  egwdffi,   gwut, dttdf, dttke, ro_adjust, tau_adjust, &
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
  real(GW_PRC), intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(GW_PRC), intent(in) :: nm(ncol,pver), ni(ncol,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(GW_PRC), intent(in) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(GW_PRC), intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real(GW_PRC), intent(in) :: effgw(ncol)
  ! Wave phase speeds for each column.
  real(GW_PRC), intent(in) :: c(ncol,-band%ngwv:band%ngwv)
  ! Molecular thermal diffusivity.
  real(GW_PRC), intent(in) :: kvtt(ncol,pver+1)

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
  real(GW_PRC), intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real(GW_PRC), intent(out) :: ttgw(ncol,pver)
!++jtb 3/2020
  ! Gravity wave constituent tendency.
  !real(GW_PRC), intent(out) :: qtgw(:,:,:)
!--jtb

  ! Effective gravity wave diffusivity at interfaces.
  real(GW_PRC), intent(out) :: egwdffi(ncol,pver+1)

  ! Gravity wave wind tendency for each wave.
  real(GW_PRC), intent(out) :: gwut(ncol,pver,-band%ngwv:band%ngwv)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(GW_PRC), intent(out) :: dttdf(ncol,pver)
  real(GW_PRC), intent(out) :: dttke(ncol,pver)

  ! Adjustment parameter for IGWs.
  real(GW_PRC), intent(in), optional :: &
       ro_adjust(ncol,-band%ngwv:band%ngwv,pver+1)

  ! Adjustment parameter for TAU.
  real(GW_PRC), intent(in), optional :: &
       tau_adjust(ncol,pver+1)

  ! Diagnosed horizontal wavenumber for ridges.
  real(GW_PRC), intent(in), optional :: &
       kwvrdg(ncol)

  ! Factor for saturation calculation. Here backwards 
  ! compatibility. I believe it should be 1.0 (jtb). 
  ! Looks like it has been 2.0 for a while in CAM.
  real(GW_PRC), intent(in), optional :: &
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

  logical :: lapply_effgw

  ! LU decomposition.
  type(TriDiagDecomp) :: decomp


  allocate( alpha(pver+1) )
  alpha(:)=0._GW_PRC

  if (present(satfac_in)) then
     satfac = satfac_in
  else
     satfac = 2._GW_PRC
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

  utgw = 0._GW_PRC
  vtgw = 0._GW_PRC

  gwut = 0._GW_PRC

  dttke = 0._GW_PRC
  ttgw = 0._GW_PRC
  egwdffi=0._GW_PRC
  dttdf=0._GW_PRC

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi = 0._GW_PRC
  taudmp = 0._GW_PRC
  tausat = 0._GW_PRC
  ubmc = 0._GW_PRC
  ubmc2 = 0._GW_PRC
  wrk = 0._GW_PRC

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

        tausat(i) = 0.0_GW_PRC
        taudmp(i) = 0.0_GW_PRC

        if (src_level(i) >= k) then

          ! Determine the absolute value of the saturation stress.
          ! Define critical levels where the sign of (u-c) changes between
          ! interfaces.
          ubmc(i) = ubi(i,k) - c(i,l)

          ! Test to see if u-c has the same sign here as the level below.
          if (ubmc(i) > 0.0_GW_PRC .eqv. ubi(i,k+1) > c(i,l)) then
              tausat(i) = abs(effkwv(i) * rhoi(i,k) * ubmc(i)**3 / &
                 (satfac*ni(i,k)))
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
          mi(i) = ni(i,k) / (2._GW_PRC * effkwv(i) * ubmc2(i)) * &  ! Is this 2._GW_PRC related to satfac?
                 (alpha(k) + ni(i,k)**2/ubmc2(i) * d(i))
          wrk(i) = -2._GW_PRC*mi(i)*rog*t(i,k)*(piln(i,k+1) - piln(i,k))
          wrk(i) = max( wrk(i), -200._GW_PRC ) * exp(wrk(i))
          taudmp(i) = tau(i,l,k+1) * exp(wrk(i))
          ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
          ! we limit tau, so instead limit the arrays used to set it.
          if (tausat(i) <= taumin) tausat(i) = 0._GW_PRC
          if (taudmp(i) <= taumin) taudmp(i) = 0._GW_PRC
          tau(i,l,k) = min(taudmp(i), tausat(i))

        endif

        end do
     end do
  end do
  
  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0._GW_PRC

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
!$OMP                                  gwut,ubt,xv,yv,lapply_effgw,ubt_lim_ratio,tend_level) &
!$OMP                          private(k,l,i,ubtl)
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0_GW_PRC

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
           ubt_lim_ratio(i) = 1._GW_PRC
          end if
        end do
     else
        ubt_lim_ratio = 1._GW_PRC
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
         if ( abs(gwut(i,k,l)) < 1.e-15_GW_PRC ) then
           gwut(i,k,l) = 0._GW_PRC
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
!$OMP parallel do default(none) shared(kbot_tend,ktop,band,dttke,ubm,c,gwut) &
!$OMP                          private(k,l)
  do k = ktop, kbot_tend
     do l = -band%ngwv, band%ngwv
        dttke(:,k) = dttke(:,k) - (ubm(:,k) - c(:,l)) * gwut(:,k,l)
     end do
  end do

  ttgw = dttke / cpair

deallocate( alpha )

end subroutine gw_drag_prof


!==========================================================================

#ifdef NOT_USED
!==========================================================================
! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Projected stresses.
  real(GW_PRC), intent(in) :: taucd(:,:,:)
  ! Components of momentum change sourced from the bottom.
  real(GW_PRC), intent(out) :: um_flux(:), vm_flux(:)

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
  real(GW_PRC), intent(in)  :: pref( pver+1 )
  
  ! Interpolated Newtonian cooling coefficients.
  real(GW_PRC), intent(out) :: alpha(pver+1)

  ! Levels of pre-calculated Newtonian cooling (1/day).
  ! The following profile is digitized from:
  ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5

  integer :: k
  integer, parameter :: nalph = 71
  real(GW_PRC) :: alpha0(nalph) = [ &
       0.1_GW_PRC,         0.1_GW_PRC,         0.1_GW_PRC,         0.1_GW_PRC,         &
       0.1_GW_PRC,         0.1_GW_PRC,         0.1_GW_PRC,         0.1_GW_PRC,         &
       0.1_GW_PRC,         0.1_GW_PRC,         0.10133333_GW_PRC,  0.104_GW_PRC,       &
       0.108_GW_PRC,       0.112_GW_PRC,       0.116_GW_PRC,       0.12066667_GW_PRC,  &
       0.126_GW_PRC,       0.132_GW_PRC,       0.138_GW_PRC,       0.144_GW_PRC,       &
       0.15133333_GW_PRC,  0.16_GW_PRC,        0.17_GW_PRC,        0.18_GW_PRC,        &
       0.19_GW_PRC,        0.19933333_GW_PRC,  0.208_GW_PRC,       0.216_GW_PRC,       &
       0.224_GW_PRC,       0.232_GW_PRC,       0.23466667_GW_PRC,  0.232_GW_PRC,       &
       0.224_GW_PRC,       0.216_GW_PRC,       0.208_GW_PRC,       0.20133333_GW_PRC,  &
       0.196_GW_PRC,       0.192_GW_PRC,       0.188_GW_PRC,       0.184_GW_PRC,       &
       0.18266667_GW_PRC,  0.184_GW_PRC,       0.188_GW_PRC,       0.192_GW_PRC,       &
       0.196_GW_PRC,       0.19333333_GW_PRC,  0.184_GW_PRC,       0.168_GW_PRC,       &
       0.152_GW_PRC,       0.136_GW_PRC,       0.12133333_GW_PRC,  0.108_GW_PRC,       &
       0.096_GW_PRC,       0.084_GW_PRC,       0.072_GW_PRC,       0.061_GW_PRC,       &
       0.051_GW_PRC,       0.042_GW_PRC,       0.033_GW_PRC,       0.024_GW_PRC,       &
       0.017666667_GW_PRC, 0.014_GW_PRC,       0.013_GW_PRC,       0.012_GW_PRC,       &
       0.011_GW_PRC,       0.010333333_GW_PRC, 0.01_GW_PRC,        0.01_GW_PRC,        &
       0.01_GW_PRC,        0.01_GW_PRC,        0.01_GW_PRC                         &
       ]

  ! Pressure levels that were used to calculate alpha0 (hPa).
  real(GW_PRC) :: palph(nalph) = [ &
       2.06115E-06_GW_PRC, 2.74280E-06_GW_PRC, 3.64988E-06_GW_PRC, 4.85694E-06_GW_PRC, &
       6.46319E-06_GW_PRC, 8.60065E-06_GW_PRC, 1.14450E-05_GW_PRC, 1.52300E-05_GW_PRC, &
       2.02667E-05_GW_PRC, 2.69692E-05_GW_PRC, 3.58882E-05_GW_PRC, 4.77568E-05_GW_PRC, &
       6.35507E-05_GW_PRC, 8.45676E-05_GW_PRC, 0.000112535_GW_PRC, 0.000149752_GW_PRC, &
       0.000199277_GW_PRC, 0.000265180_GW_PRC, 0.000352878_GW_PRC, 0.000469579_GW_PRC, &
       0.000624875_GW_PRC, 0.000831529_GW_PRC, 0.00110653_GW_PRC,  0.00147247_GW_PRC,  &
       0.00195943_GW_PRC,  0.00260744_GW_PRC,  0.00346975_GW_PRC,  0.00461724_GW_PRC,  &
       0.00614421_GW_PRC,  0.00817618_GW_PRC,  0.0108801_GW_PRC,   0.0144783_GW_PRC,   &
       0.0192665_GW_PRC,   0.0256382_GW_PRC,   0.0341170_GW_PRC,   0.0453999_GW_PRC,   &
       0.0604142_GW_PRC,   0.0803939_GW_PRC,   0.106981_GW_PRC,    0.142361_GW_PRC,    &
       0.189442_GW_PRC,    0.252093_GW_PRC,    0.335463_GW_PRC,    0.446404_GW_PRC,    &
       0.594036_GW_PRC,    0.790490_GW_PRC,    1.05192_GW_PRC,     1.39980_GW_PRC,     &
       1.86273_GW_PRC,     2.47875_GW_PRC,     3.29851_GW_PRC,     4.38936_GW_PRC,     &
       5.84098_GW_PRC,     7.77266_GW_PRC,     10.3432_GW_PRC,     13.7638_GW_PRC,     &
       18.3156_GW_PRC,     24.3728_GW_PRC,     32.4332_GW_PRC,     43.1593_GW_PRC,     &
       57.4326_GW_PRC,     76.4263_GW_PRC,     101.701_GW_PRC,     135.335_GW_PRC,     &
       180.092_GW_PRC,     239.651_GW_PRC,     318.907_GW_PRC,     424.373_GW_PRC,     &
       564.718_GW_PRC,     751.477_GW_PRC,     1000._GW_PRC                        &
       ]

  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa
  do k=1,nalph
     alpha0(k) = alpha0(k) / 86400._GW_PRC
     alpha0(k) = max(alpha0(k), 1.e-6_GW_PRC)
     palph(k) = palph(k)*1.e2_GW_PRC
  end do

  ! interpolate to current vertical grid and obtain alpha
  call lininterp (alpha0  ,palph, nalph , alpha  , pref , pver+1)


end subroutine gw_newtonian_set
#endif

end module gw_common
