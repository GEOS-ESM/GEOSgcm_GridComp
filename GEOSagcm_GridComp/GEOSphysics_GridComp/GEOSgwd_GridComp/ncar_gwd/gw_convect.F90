module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

  use gw_utils, only: GW_PRC, GW_R8, get_unit_vector, dot_2d, midpoint_interp
  use gw_common, only: GWBand, qbo_hdepth_scaling, gw_drag_prof, hr_cf, &
                       calc_taucd, momentum_flux, momentum_fixer, &
                       energy_momentum_adjust, energy_change, energy_fixer 

  use MAPL_ConstantsMod, only: MAPL_RGAS, MAPL_CP, MAPL_GRAV

implicit none
private
save

public :: BeresSourceDesc
public :: gw_beres_ifc
public :: gw_beres_src
public :: gw_beres_init

real, parameter :: PI      = 3.14159265358979323846 ! pi
real, parameter :: rad2deg = 180./PI

type :: BeresSourceDesc
   logical :: active
   ! Whether wind speeds are shifted to be relative to storm cells.
   logical :: storm_shift
   ! Index for level where wind speed is used as the source speed.
   integer :: k
   ! Heating depths below this value [m] will be ignored.
   real :: min_hdepth
   ! Source for wave spectrum
   real :: spectrum_source
   ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
   integer :: maxh
   integer :: maxuh
   ! Heating depths [m].
   real, allocatable :: hd(:)
   ! Table of source spectra.
   real, allocatable :: mfcc(:,:,:)
   ! Forced background for extratropics
   real, allocatable :: taubck(:,:)
end type BeresSourceDesc


contains

!==========================================================================

!------------------------------------
subroutine gw_beres_init (file_name, band, desc, pgwv, gw_dc, fcrit2, wavelength, &
                          spectrum_source, min_hdepth, storm_shift, taubgnd, active, ncol, lats)
#include <netcdf.inc>

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(inout) :: band

  type(BeresSourceDesc), intent(inout) :: desc

  integer, intent(in) :: pgwv, ncol
  real, intent(in) :: gw_dc, fcrit2, wavelength
  real, intent(in) :: spectrum_source, min_hdepth, taubgnd
  logical, intent(in) :: storm_shift, active
  real, intent(in) :: lats(ncol)

  ! Stuff for Beres convective gravity wave source.
  real(GW_R8), allocatable :: mfcc(:,:,:), hdcc(:)
  integer  :: hd_mfcc , mw_mfcc, ps_mfcc, ngwv_file, ps_mfcc_mid

  ! For forced background extratropical wave speed
  real    :: c4, latdeg, flat_gw
  real, allocatable :: c0(:), cw4(:)
  integer :: i, kc

  ! Vars needed by NetCDF operators
  integer  :: ncid, dimid, varid, status
  
  status = nf_open(file_name , 0, ncid)

  status = NF_INQ_DIMID(ncid, 'PS', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, ps_mfcc )

  status = NF_INQ_DIMID(ncid, 'MW', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, mw_mfcc )

  status = NF_INQ_DIMID(ncid, 'HD', dimid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_INQ_DIMLEN(ncid, dimid, hd_mfcc )

  allocate( mfcc(hd_mfcc , mw_mfcc, ps_mfcc) )
  allocate( hdcc(hd_mfcc) )
   
  status = NF_INQ_VARID(ncid, 'HD', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, hdcc )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = NF_INQ_VARID(ncid, 'mfcc', varid)
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
  status = NF_GET_VAR_DOUBLE(ncid, varid, mfcc )
  IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

  status = nf_close (ncid)

  band  = GWBand(pgwv, gw_dc, fcrit2, wavelength )

  ! These dimensions; {HD,MW,PS}_MFCC, came from Beres forcing file.

  ! Get HD (heating depth) dimension.
  desc%maxh = HD_MFCC

  ! Get MW (mean wind) dimension.
  desc%maxuh = MW_MFCC

  ! Get PS (phase speed) dimension.
  ngwv_file = PS_MFCC

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2
  ! midpoint of spectrum in netcdf file is ps_mfcc (odd number) -1 divided by 2, plus 1
  ! E.g., ps_mfcc = 5. So, ps_mfcc_mid = 3
  !       1  2  3  4  5
  !      -2 -1  0 +1 +2
  ps_mfcc_mid= (ngwv_file-1)/2

  desc%active = active
  if (active) then

    allocate(desc%hd(desc%maxh) , stat=status )

    allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,-band%ngwv:band%ngwv), stat=status )

    desc%mfcc( : , -desc%maxuh:desc%maxuh , -band%ngwv            :band%ngwv             ) & 
       = mfcc( :,             :           , -band%ngwv+ps_mfcc_mid:band%ngwv+ps_mfcc_mid )
  
    ! While not currently documented in the file, it uses kilometers. Convert
    ! to meters.
    desc%hd = hdcc * 1000.0

    desc%spectrum_source = spectrum_source

    desc%min_hdepth = min_hdepth

    desc%storm_shift=storm_shift

    ! Intialize forced background wave speeds
    allocate(desc%taubck(ncol,-band%ngwv:band%ngwv))
    allocate(c0(-band%ngwv:band%ngwv))
    allocate(cw4(-band%ngwv:band%ngwv))
    desc%taubck = 0.0
    c0  = 0.0
    cw4 = 0.0
    do kc = -4,4
        c4 =  10.0*kc
       cw4(kc) =  exp(-(c4/30.)**2)
    enddo
    do kc = -band%ngwv,band%ngwv
       c0(kc) =  10.0*(4.0/real(band%ngwv))*kc
       desc%taubck(:,kc) =  exp(-(c0(kc)/30.)**2)
    enddo
    do i=1,ncol
      ! include forced background stress in extra tropics
       ! Determine the background stress at c=0
       ! Include dependence on latitude:
       latdeg = lats(i)*rad2deg
       if (-15.3 < latdeg .and. latdeg < 15.3) then
         flat_gw =  0.10
       else if (latdeg > -31. .and. latdeg <= -15.3) then
         flat_gw =  0.10
       else if (latdeg <  31. .and. latdeg >=  15.3) then
         flat_gw =  0.10
       else if (latdeg > -60. .and. latdeg <= -31.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
       else if (latdeg <  60. .and. latdeg >=  31.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
       else if (latdeg <= -60.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
       else if (latdeg >=  60.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
       end if
       desc%taubck(i,:) = taubgnd*0.001*flat_gw*desc%taubck(i,:)*(sum(cw4)/sum(desc%taubck(i,:)))
    enddo
    deallocate( c0, cw4 )
  end if
    
end subroutine gw_beres_init

!------------------------------------
subroutine gw_beres_src(ncol, pver, band, desc, u, v, &
     dqcdt, netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     c, hdepth, maxq0, lats)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
  !!!use gw_utils, only: get_unit_vector, dot_2d, midpoint_interp
  !!!use gw_common, only: GWBand, qbo_hdepth_scaling

!------------------------------Arguments--------------------------------
  ! Column and vertical dimensions.
  integer, intent(in) :: ncol, pver

  ! Wavelengths triggered by convection.
  type(GWBand), intent(in) :: band

  ! Settings for convection type (e.g. deep vs shallow).
  type(BeresSourceDesc), intent(in) :: desc

  ! Midpoint zonal/meridional winds.
  real, intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Condensate tendency due to large-scale (kg kg-1 s-1)
  real, intent(in) :: dqcdt(ncol,pver)
  ! Heating rate due to convection.
  real, intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real, intent(in) :: zm(ncol,pver)
  ! latitudes.
  real, intent(in) :: lats(ncol)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(GW_PRC), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real, intent(out) :: ubm(ncol,pver)
  real, intent(out) :: ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(GW_PRC), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Heating depth [m] and maximum heating in each column.
  real, intent(out) :: hdepth(ncol), maxq0(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real :: uconv(ncol), vconv(ncol)

  ! Maximum heating rate.
  real(GW_PRC) :: q0(ncol)

  ! Bottom/top heating range index.
  integer  :: boti(ncol), topi(ncol)
  ! Index for looking up heating depth dimension in the table.
  integer  :: hd_idx(ncol)
  ! Mean wind in heating region.
  real(GW_PRC) :: uh(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer :: Umini(ncol), Umaxi(ncol)
  ! Source level tau for a column.
  real(GW_PRC) :: tau0(-band%ngwv:band%ngwv)
  ! Speed of convective cells relative to storm.
  real(GW_PRC) :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Averaging length.
  real, parameter :: AL = 1.0e5

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------
  tau = 0.0
  hdepth = 0.0
  q0 = 0.0
  tau0 = 0.0

  !------------------------------------------------------------------------
  ! Determine wind and unit vectors approximately at the source level, then
  ! project winds.
  !------------------------------------------------------------------------

  ! Source wind speed and direction.
  uconv = u(:,desc%k)
  vconv = v(:,desc%k)

  ! Get the unit vector components and magnitude at the source level.
  call get_unit_vector(uconv, vconv, xv, yv, ubi(:,desc%k+1))

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = midpoint_interp(ubm)

  !-----------------------------------------------------------------------
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.
  boti = 0
  topi = 0
  do k = pver, 1, -1
     do i = 1, ncol
        if (boti(i) == 0) then
           ! Detect if we are outside the maximum range (where z = 20 km).
           if (zm(i,k) >= 20000.0) then
              boti(i) = k
              topi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0) boti(i) = k
           end if
        else if (topi(i) == 0) then
           ! Detect if we are outside the maximum range (z = 20 km).
           if (zm(i,k) >= 20000.0) then
              topi(i) = k
           else
              ! First spot where heating rate is no longer positive.
              if (netdt(i,k) <= 0.0) topi(i) = k
           end if
        end if
     end do
     ! When all done, exit.
     if (all(topi /= 0)) exit
  end do

  ! Heating depth in m.
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
  hdepth = hdepth*qbo_hdepth_scaling

  hd_idx = index_of_nearest(hdepth, desc%hd)

  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though.
  where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0

  ! Maximum heating rate.
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        q0 = max(q0, netdt(:,k))
     end where
  end do

  !output max heating rate in K/day
  maxq0 = q0*86400.0

  ! Multipy by conversion factor
  q0 = q0 * hr_cf

  if (desc%storm_shift) then

     ! Find the cell speed where the storm speed is > 10 m/s.
     ! Storm speed is taken to be the source wind speed.
     CS = sign(max(abs(ubm(:,desc%k))-10.0, 0.0), ubm(:,desc%k))

     ! Average wind in heating region, relative to storm cells.
     uh = 0.0
     do k = minval(topi), maxval(boti)
        where (k >= topi .and. k <= boti)
           uh = uh + ubm(:,k)/(boti-topi+1)
        end where
     end do

     uh = uh - CS

  else

     ! For shallow convection, wind is relative to ground, and "heating
     ! region" wind is just the source level wind.
     uh = ubm(:,desc%k)

  end if

  ! Limit uh to table range.
  uh = min(uh,  real(desc%maxuh))
  uh = max(uh, -real(desc%maxuh))

  ! Speeds for critical level filtering.
  Umini =  band%ngwv
  Umaxi = -band%ngwv
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:,k)/band%dc))
        Umaxi = max(Umaxi, nint(ubm(:,k)/band%dc))
     end where
  end do

  Umini = max(Umini, -band%ngwv)
  Umaxi = min(Umaxi, band%ngwv)

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if the heating depth is large enough, else set
     ! tau0 = 0.
     !---------------------------------------------------------------------

     if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = desc%mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc%storm_shift) then
           ! For deep convection, the wind was relative to storm cells, so
           ! shift the spectrum so that it is now relative to the ground.
           shift = -nint(CS(i)/band%dc)
           tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude.
        tau0 = tau0*(q0(i)**2)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0
 
        tau(i,:,topi(i)+1) = tau0

     elseif (dqcdt(i,desc%k) > 1.e-8) then ! frontal region (large-scale forcing)

      ! include forced background stress in extra tropical large-scale systems
       ! Set the phase speeds and wave numbers in the direction of the source wind.
       ! Set the source stress magnitude (positive only, note that the sign of the 
       ! stress is the same as (c-u).
       tau(i,:,desc%k+1) = desc%taubck(i,:)
       topi(i) = desc%k

    endif

  enddo
  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = topi
  tend_level = topi

  ! Set phase speeds; just use reference speeds.
  c = spread(band%cref, 1, ncol)

end subroutine gw_beres_src



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Main Interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gw_beres_ifc( band, &
   ncol, pver, dt, effgw_dp,  &
   u, v, t, pref, pint, delp, rdelp, piln, &
   zm, zi, nm, ni, rhoi, kvtt,  &
   dqcdt, &
   netdt,desc,lats, &
   utgw,vtgw,ttgw,flx_heat)

   type(BeresSourceDesc), intent(inout) :: desc
   type(GWBand), intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,      intent(in) :: ncol         ! number of atmospheric columns
   integer,      intent(in) :: pver         ! number of vertical layers
   real,         intent(in) :: dt           ! Time step.
   real,         intent(in) :: effgw_dp

   real,         intent(in) :: u(ncol,pver)      ! Midpoint zonal winds. ( m s-1)
   real,         intent(in) :: v(ncol,pver)      ! Midpoint meridional winds. ( m s-1)
   real,         intent(in) :: t(ncol,pver)      ! Midpoint temperatures. (K)
   real,         intent(in) :: dqcdt(ncol,pver)  ! Condensate tendency due to large-scale (kg kg-1 s-1)
   real,         intent(in) :: netdt(ncol,pver)  ! Convective heating rate (K s-1)
   real,         intent(in) :: pref(pver+1)      ! Reference pressure at interfaces (Pa !!! )
   real,         intent(in) :: piln(ncol,pver+1) ! Log of interface pressures.
   real,         intent(in) :: pint(ncol,pver+1) ! Interface pressures. (Pa)
   real,         intent(in) :: delp(ncol,pver)   ! Layer pressures thickness. (Pa)
   real,         intent(in) :: rdelp(ncol,pver)  ! Inverse pressure thickness. (Pa-1)
   real,         intent(in) :: zm(ncol,pver)     ! Midpoint altitudes above ground (m).
   real,         intent(in) :: zi(ncol,pver+1)   ! Interface altitudes above ground (m).
   real, intent(in) :: nm(ncol,pver)     ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real, intent(in) :: ni(ncol,pver+1)   ! Interface Brunt-Vaisalla frequencies (s-1).
   real, intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real, intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.

   real,         intent(in) :: lats(ncol)      ! latitudes

   real, intent(out) :: utgw(ncol,pver)       ! zonal wind tendency
   real, intent(out) :: vtgw(ncol,pver)       ! meridional wind tendency
   real, intent(out) :: ttgw(ncol,pver)       ! temperature tendency
   real, intent(inout) :: flx_heat(ncol)        ! Energy change

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn

   real(GW_PRC), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(GW_PRC), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(GW_PRC), allocatable :: c(:,:)

   ! Efficiency for a gravity wave source.
   real :: effgw(ncol)

   ! Momentum fluxes used by fixer.
   real :: um_flux(ncol), vm_flux(ncol)

   ! Energy change used by fixer.
   real :: de(ncol)

   ! Reynolds stress for waves propagating in each cardinal direction.
   real :: taucd(ncol,pver+1,4)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real :: ubm(ncol,pver)
   real :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real :: xv(ncol)
   real :: yv(ncol)

   ! Heating depth [m] and maximum heating in each column.
   real :: hdepth(ncol), maxq0(ncol)

   real :: pint_adj(ncol,pver+1)
   real :: zfac_layer

   character(len=1) :: cn
   character(len=9) :: fname(4)

   integer :: i,j,l

   !----------------------------------------------------------------------------


   ! Allocate wavenumber fields.
   allocate(tau(ncol,-band%ngwv:band%ngwv,pver+1))
   allocate(gwut(ncol,pver,-band%ngwv:band%ngwv))
   allocate(c(ncol,-band%ngwv:band%ngwv))

     ! Efficiency of gravity wave momentum transfer.
     ! This is really only to remove the pole points.
     where (pi/2.0 - abs(lats(:ncol)) >= 1.e-4 )  !-4*epsilon(1.0))
        effgw = effgw_dp
     elsewhere
        effgw = 0.0
     end where

     do k = 0, pver
        ! spectrum source index
        if (pref(k+1) < desc%spectrum_source) desc%k = k+1
     end do

     ! Determine wave sources for Beres deep scheme
     call gw_beres_src(ncol, pver, band, desc, &
          u, v, dqcdt, netdt, zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, c, hdepth, maxq0, lats)

!WMP pressure scaling near model top
!    zfac_layer = 15.0 ! 0.15mb
!    pint_adj = 0.5*(1+TANH(((2.0*pint/zfac_layer)-1)/0.25))
     pint_adj = 1.0

      ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, pver, band, pint, delp, rdelp, & 
          src_level, tend_level, dt, t,    &
          piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
          effgw, c, kvtt, tau, utgw, vtgw, &
          ttgw, gwut, tau_adjust=pint_adj)

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
   ! call energy_momentum_adjust(ncol, pver, desc%k, band, pint, delp, c, tau, &
   !                             effgw, t, ubm, ubi, xv, yv, utgw, vtgw, ttgw)
#endif
 
   deallocate(tau, gwut, c)

end subroutine gw_beres_ifc


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
  real,     intent(in) :: x(:)
  real, intent(in) :: grid(:)

  integer :: idx(size(x))

  real :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2.d0

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

end module gw_convect
