module gw_oro

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

  use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  use gw_common, only: GWBand, rair, gw_drag_prof 

implicit none
private
save

public :: gw_oro_ifc
public :: gw_oro_src
public :: gw_oro_init

integer,parameter :: r8 = selected_real_kind(12) ! 8 byte real

real(R8),parameter :: PI      = 3.14159265358979323846_R8  ! pi

real(r8) :: gw_oro_south_fac

contains

!==========================================================================

!------------------------------------
subroutine gw_oro_init (band )
#include <netcdf.inc>

  type(GWBand), intent(inout) :: band
  real(r8) :: gw_dc,wavelength
  integer  :: pgwv

  

! Need to call GWBand for oro waves

  ! Hardwire for now
  gw_dc = 2.5_r8
  pgwv  = 0
  wavelength = 1.e5_r8
  band  = GWBand(pgwv, gw_dc, 1.0_r8, wavelength )

  gw_oro_south_fac = 2.0_r8
  
end subroutine gw_oro_init

!------------------------------------
subroutine gw_oro_src(ncol,pver, band, &
     pint, pmid, delp, &
     u, v, t, sgh, zm, nm, &
     src_level, tend_level, tau, ubm, ubi, xv, yv, c)
  !-----------------------------------------------------------------------
  ! Orographic source for multiple gravity wave drag parameterization.
  !
  ! The stress is returned for a single wave with c=0, over orography.
  ! For points where the orographic variance is small (including ocean),
  ! the returned stress is zero.
  !-----------------------------------------------------------------------
  !!! use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !!! use gw_common, only: GWBand,rair

  ! Column dimension.
  integer, intent(in) :: ncol, pver
  ! Band to emit orographic waves in.
  ! Regardless, we will only ever emit into l = 0.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  !type(Coords1D), intent(in) :: p

  ! Interface pressures. (Pa)
  real(r8), intent(in) :: pint(ncol,pver+1)
  ! Midpoint pressures. (Pa)
  real(r8), intent(in) :: pmid(ncol,pver)
  ! Delta Interface pressures. (Pa)
  real(r8), intent(in) :: delp(ncol,pver)


  ! Midpoint zonal/meridional winds.
  real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Standard deviation of orography.
  real(r8), intent(in) :: sgh(ncol)
  ! Midpoint altitudes.
  real(r8), intent(in) :: zm(ncol,pver)
  ! Midpoint Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Surface streamline displacement height (2*sgh).
  real(r8) :: hdsp(ncol)
  ! Max orographic standard deviation to use.
  real(r8) :: sghmax
  ! c=0 stress from orography.
  real(r8) :: tauoro(ncol)
  ! Averages over source region.
  real(r8) :: nsrc(ncol) ! B-V frequency.
  real(r8) :: rsrc(ncol) ! Density.
  real(r8) :: usrc(ncol) ! Zonal wind.
  real(r8) :: vsrc(ncol) ! Meridional wind.

  ! Difference in interface pressure across source region.
  real(r8) :: dpsrc(ncol)

  ! Limiters (min/max values)
  ! min surface displacement height for orographic waves
  real(r8), parameter :: orohmin = 10._r8
  ! min wind speed for orographic waves
  real(r8), parameter :: orovmin = 2._r8

!--------------------------------------------------------------------------
! Average the basic state variables for the wave source over the depth of
! the orographic standard deviation. Here we assume that the appropiate
! values of wind, stability, etc. for determining the wave source are
! averages over the depth of the atmosphere penterated by the typical
! mountain.
! Reduces to the bottom midpoint values when sgh=0, such as over ocean.
!--------------------------------------------------------------------------

  hdsp = 2.0_r8 * sgh

  k = pver
  src_level = k-1
  rsrc = pmid(:,k)/(rair*t(:,k)) * delp(:,k)
  usrc = u(:,k) * delp(:,k)
  vsrc = v(:,k) * delp(:,k)
  nsrc = nm(:,k)* delp(:,k)

  do k = pver-1, 1, -1
     do i = 1, ncol
        if (hdsp(i) > sqrt(zm(i,k)*zm(i,k+1))) then
           src_level(i) = k-1
           rsrc(i) = rsrc(i) + &
                pmid(i,k) / (rair*t(i,k)) * delp(i,k)
           usrc(i) = usrc(i) + u(i,k) * delp(i,k)
           vsrc(i) = vsrc(i) + v(i,k) * delp(i,k)
           nsrc(i) = nsrc(i) + nm(i,k)* delp(i,k)
        end if
     end do
     ! Break the loop when all source levels found.
     if (all(src_level >= k)) exit
  end do


  do i = 1, ncol
     dpsrc(i) = pint(i,pver+1) - pint(i,src_level(i)+1)
  end do

  
  rsrc = rsrc / dpsrc
  usrc = usrc / dpsrc
  vsrc = vsrc / dpsrc
  nsrc = nsrc / dpsrc


  ! Get the unit vector components and magnitude at the surface.
  call get_unit_vector(usrc, vsrc, xv, yv, ubi(:,pver+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  ! Determine the orographic c=0 source term following McFarlane (1987).
  ! Set the source top interface index to pver, if the orographic term is
  ! zero.
  do i = 1, ncol
     if ((ubi(i,pver+1) > orovmin) .and. (hdsp(i) > orohmin)) then
        sghmax = band%fcrit2 * (ubi(i,pver+1) / nsrc(i))**2
        tauoro(i) = 0.5_r8 * band%kwv * min(hdsp(i)**2, sghmax) * &
             rsrc(i) * nsrc(i) * ubi(i,pver+1)
     else
        tauoro(i) = 0._r8
        src_level(i) = pver
     end if
  end do

  ! Set the phase speeds and wave numbers in the direction of the source
  ! wind. Set the source stress magnitude (positive only, note that the
  ! sign of the stress is the same as (c-u).
  tau = 0._r8
  do k = pver, minval(src_level), -1
     where (src_level <= k) tau(:,0,k+1) = tauoro
  end do

  ! Allow wind tendencies all the way to the model bottom.
  tend_level = pver

  ! No spectrum; phase speed is just 0.
  c = 0._r8

end subroutine gw_oro_src


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Main Interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gw_oro_ifc( band, &
   ncol, pver, dt, effgw_oro,  &
   u, v, t, pint, pmid, &
   delp, rdelp, piln, &
   zm, zi, nm, ni, rhoi, kvtt,  &
   sgh, lats, &
   utgw,vtgw,ttgw, flx_heat)

   !!!use coords_1d,  only: Coords1D
   !!!use gw_convect,     only: gw_beres_src

!!use cesm_const_mod,   only: pi=>shr_const_pi    !, cl=>shr_kind_cl
!!! use gw_common,  only: gw_drag_prof 
!use coords_1d,  only: Coords1D
!++ jtb 3/2020
!! use cesm_physics_types,  only: physics_ptend
!! use cesm_constituent, only: pcnst

   type(GWBand),     intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,          intent(in) :: ncol         ! number of atmospheric columns
   integer,          intent(in) :: pver         ! number of vertical layers
   real(r8),         intent(in) :: dt           ! Time step.
   real(r8),         intent(in) :: effgw_oro

   real(r8),         intent(in) :: u(ncol,pver)      ! Midpoint zonal winds. ( m s-1)
   real(r8),         intent(in) :: v(ncol,pver)      ! Midpoint meridional winds. ( m s-1)
   real(r8),         intent(in) :: t(ncol,pver)      ! Midpoint temperatures. (K)
   real(r8),         intent(in) :: piln(ncol,pver+1) ! Log of interface pressures.
   real(r8),         intent(in) :: pmid(ncol,pver)   ! Midpoint pressures. (Pa)
   real(r8),         intent(in) :: pint(ncol,pver+1) ! Interface pressures. (Pa)
   real(r8),         intent(in) :: delp(ncol,pver)   ! Layer pressures thickness. (Pa)
   real(r8),         intent(in) :: rdelp(ncol,pver)  ! Inverse pressure thickness. (Pa-1)
   real(r8),         intent(in) :: zm(ncol,pver)     ! Midpoint altitudes above ground (m).
   real(r8),         intent(in) :: zi(ncol,pver+1)   ! Interface altitudes above ground (m).
   real(r8),         intent(in) :: nm(ncol,pver)     ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: ni(ncol,pver+1)   ! Interface Brunt-Vaisalla frequencies (s-1).
   real(r8),         intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(r8),         intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.
!++jtb 3/2020
   !!! real(r8),         intent(in) :: q(:,:,:)          ! Constituent array.
   !!! real(r8),         intent(in) :: dse(ncol,pver)    ! Dry static energy.

   real(r8),         intent(in) :: sgh(ncol)       ! subgrid orographic std dev (m)
   real(r8),         intent(in) :: lats(ncol)      ! latitudes


   !! type(physics_ptend), intent(inout):: ptend   ! Parameterization net tendencies.

   real(r8),        intent(out) :: flx_heat(ncol)
   real(r8),        intent(out) :: utgw(ncol,pver)       ! zonal wind tendency
   real(r8),        intent(out) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(r8),        intent(out) :: ttgw(ncol,pver)       ! temperature tendency

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn

   real(r8), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(r8), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(r8), allocatable :: c(:,:)

   ! Efficiency for a gravity wave source.
   real(r8) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real(r8) :: ubm(ncol,pver)
   real(r8) :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real(r8) :: xv(ncol)
   real(r8) :: yv(ncol)

   ! Averages over source region.
   real(r8) :: ubmsrc(ncol) ! On-ridge wind.
   real(r8) :: usrc(ncol)   ! Zonal wind.
   real(r8) :: vsrc(ncol)   ! Meridional wind.
   real(r8) :: nsrc(ncol)   ! B-V frequency.
   real(r8) :: rsrc(ncol)   ! Density.


   ! Wave Reynolds stresses at source level
   real(r8) :: tauoro(ncol)
   real(r8) :: taudsw(ncol)

   ! Wave breaking level
   real(r8) :: wbr(ncol)

   !!! real(r8) :: qtgw(ncol,pver,pcnst) ! constituents tendencies

   ! Heating depth [m] and maximum heating in each column.
   real(r8) :: hdepth(ncol), maxq0(ncol)

   ! Effective gravity wave diffusivity at interfaces.
   real(r8) :: egwdffi(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(r8) :: dttdf(ncol,pver)
   real(r8) :: dttke(ncol,pver)

   ! Wave stress in zonal/meridional direction
   real(r8) :: taurx(ncol,pver+1)
   real(r8) :: taurx0(ncol,pver+1)
   real(r8) :: taury(ncol,pver+1)
   real(r8) :: taury0(ncol,pver+1)


   ! Energy change used by fixer.
   real(r8) :: de(ncol)
   logical  :: gw_apply_tndmax  	!- default .TRUE. for Anisotropic: "Sean" limiters

   character(len=1) :: cn
   character(len=9) :: fname(4)

   integer :: i,j

   !----------------------------------------------------------------------------


   ! Allocate wavenumber fields.
   allocate(tau(ncol,band%ngwv:band%ngwv,pver+1))
   allocate(gwut(ncol,pver,band%ngwv:band%ngwv))
   allocate(c(ncol,band%ngwv:band%ngwv))


     gw_apply_tndmax  = .FALSE.


! Efficiency of gravity wave momentum transfer.
     effgw(:) = effgw_oro


! Determine the orographic wave source
        call gw_oro_src(ncol, pver, band, pint, pmid, delp, &
             u, v, t, sgh, zm, nm, &
             src_level, tend_level, tau, ubm, ubi, xv, yv, c)


     do i = 1, ncol
        if (lats(i) < 0._r8) then
           tau(i,:,:) = tau(i,:,:) * gw_oro_south_fac
        end if
     end do
     
     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, pver, band, pint, delp, rdelp, & 
          src_level, tend_level,   dt, t,    &
          piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
          effgw,c,          kvtt,  tau,  utgw,  vtgw, &
          ttgw, egwdffi,  gwut, dttdf, dttke,            &
          satfac_in = 1._r8,                                   &
          lapply_effgw_in=gw_apply_tndmax)

     flx_heat(:ncol) = 0._r8

     
end subroutine gw_oro_ifc


!************************************************************************
!!handle_err
!************************************************************************
!
!!ROUTINE:      handle_err
!!DESCRIPTION:  error handler
!--------------------------------------------------------------------------

subroutine handle_err(status)
  
  implicit         none
  
#include <netcdf.inc>
  
  integer          status
  
  if (status .ne. nf_noerr) then
    print *, nf_strerror(status)
    stop 'Stopped'
  endif
  
end subroutine handle_err



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





! Short routine to get the indices of a set of values rounded to their
! nearest points on a grid.
function index_of_nearest(x, grid) result(idx)
  real(r8), intent(in) :: x(:)
  real(r8), intent(in) :: grid(:)

  integer :: idx(size(x))

  real(r8) :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2._r8

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

end module gw_oro
