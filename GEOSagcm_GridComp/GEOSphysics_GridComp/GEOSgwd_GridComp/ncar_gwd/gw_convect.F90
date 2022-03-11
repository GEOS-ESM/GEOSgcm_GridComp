module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

  use gw_utils, only: GW_PRC, get_unit_vector, dot_2d, midpoint_interp
  use gw_common, only: GWBand, qbo_hdepth_scaling, gw_drag_prof, hr_cf 

  use MAPL_ConstantsMod, only: MAPL_RGAS, MAPL_CP, MAPL_GRAV

implicit none
private
save

public :: BeresSourceDesc
public :: gw_beres_ifc
public :: gw_beres_src
public :: gw_beres_init

real,parameter :: PI      = 3.14159265358979323846 ! pi
real,    parameter :: TNDMAX  = 500. / 86400.  ! maximum wind tendency

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
end type BeresSourceDesc



contains

!==========================================================================

!------------------------------------
subroutine gw_beres_init (file_name, band, desc, pgwv, gw_dc, fcrit2, wavelength, &
                          spectrum_source, min_hdepth, storm_shift, active)
#include <netcdf.inc>

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(inout) :: band

  type(BeresSourceDesc), intent(inout) :: desc

  integer, intent(in) :: pgwv
  real, intent(in) :: gw_dc, fcrit2, wavelength
  real, intent(in) :: spectrum_source, min_hdepth
  logical, intent(in) :: storm_shift, active

  ! Stuff for Beres convective gravity wave source.
  real(kind(1.d0)), allocatable :: mfcc(:,:,:), hdcc(:)
  integer  :: hd_mfcc , mw_mfcc, ps_mfcc, ngwv_file, ps_mfcc_mid

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
    desc%hd = hdcc *1000.

    desc%spectrum_source = spectrum_source

    desc%min_hdepth = min_hdepth

    desc%storm_shift=storm_shift

  end if
    
end subroutine gw_beres_init

!------------------------------------
subroutine gw_beres_src(ncol, pver, band, desc, u, v, &
     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     c, hdepth, maxq0)
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
  ! Heating rate due to convection.
  real, intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real, intent(in) :: zm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real(GW_PRC), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real(GW_PRC), intent(out) :: ubm(ncol,pver)
  real(GW_PRC), intent(out) :: ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real(GW_PRC), intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(GW_PRC), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Heating depth [m] and maximum heating in each column.
  real(GW_PRC), intent(out) :: hdepth(ncol), maxq0(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real(GW_PRC) :: uconv(ncol), vconv(ncol)

  ! Maximum heating rate.
  real :: q0(ncol)

  ! Bottom/top heating range index.
  integer  :: boti(ncol), topi(ncol)
  ! Index for looking up heating depth dimension in the table.
  integer  :: hd_idx(ncol)
  ! Mean wind in heating region.
  real :: uh(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer :: Umini(ncol), Umaxi(ncol)
  ! Source level tau for a column.
  real(kind(1.d0)) :: tau0(-band%ngwv:band%ngwv)
  ! Speed of convective cells relative to storm.
  real :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Averaging length.
  real(kind(1.d0)), parameter :: AL = 1.0e5

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
     ubm(:,k) = dot_2d(real(u(:,k),GW_PRC), real(v(:,k),GW_PRC), xv, yv)
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
              if (.not. (netdt(i,k) > 0.0)) topi(i) = k
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
        tau0 = tau0*DBLE(q0(i)**2)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0
 
        tau(i,:,topi(i)+1) = tau0

     end if ! heating depth above min and not at the pole

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
   netdt,desc,lats, &
   utgw,vtgw,ttgw, flx_heat)

   type(BeresSourceDesc), intent(inout) :: desc
   type(GWBand), intent(in) :: band         ! I hate this variable  ... it just hides information from view
   integer,      intent(in) :: ncol         ! number of atmospheric columns
   integer,      intent(in) :: pver         ! number of vertical layers
   real,         intent(in) :: dt           ! Time step.
   real,         intent(in) :: effgw_dp

   real,         intent(in) :: u(ncol,pver)      ! Midpoint zonal winds. ( m s-1)
   real,         intent(in) :: v(ncol,pver)      ! Midpoint meridional winds. ( m s-1)
   real,         intent(in) :: t(ncol,pver)      ! Midpoint temperatures. (K)
   real,         intent(in) :: netdt(ncol,pver)  ! Convective heating rate (K s-1)
   real,         intent(in) :: pref(pver+1)      ! Reference pressure at interfaces (Pa !!! )
   real,         intent(in) :: piln(ncol,pver+1) ! Log of interface pressures.
   real,         intent(in) :: pint(ncol,pver+1) ! Interface pressures. (Pa)
   real,         intent(in) :: delp(ncol,pver)   ! Layer pressures thickness. (Pa)
   real,         intent(in) :: rdelp(ncol,pver)  ! Inverse pressure thickness. (Pa-1)
   real,         intent(in) :: zm(ncol,pver)     ! Midpoint altitudes above ground (m).
   real,         intent(in) :: zi(ncol,pver+1)   ! Interface altitudes above ground (m).
   real(GW_PRC), intent(in) :: nm(ncol,pver)     ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real(GW_PRC), intent(in) :: ni(ncol,pver+1)   ! Interface Brunt-Vaisalla frequencies (s-1).
   real(GW_PRC), intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real(GW_PRC), intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.

   real,         intent(in) :: lats(ncol)      ! latitudes

   real,         intent(out) :: flx_heat(ncol)
   real(GW_PRC), intent(out) :: utgw(ncol,pver)       ! zonal wind tendency
   real(GW_PRC), intent(out) :: vtgw(ncol,pver)       ! meridional wind tendency
   real(GW_PRC), intent(out) :: ttgw(ncol,pver)       ! temperature tendency

   !---------------------------Local storage-------------------------------

   integer :: k, m, nn

   real(GW_PRC), allocatable :: tau(:,:,:)  ! wave Reynolds stress
   ! gravity wave wind tendency for each wave
   real(GW_PRC), allocatable :: gwut(:,:,:)
   ! Wave phase speeds for each column
   real(GW_PRC), allocatable :: c(:,:)

   ! Efficiency for a gravity wave source.
   real(GW_PRC) :: effgw(ncol)

   ! Indices of top gravity wave source level and lowest level where wind
   ! tendencies are allowed.
   integer :: src_level(ncol)
   integer :: tend_level(ncol)

   ! Projection of wind at midpoints and interfaces.
   real(GW_PRC) :: ubm(ncol,pver)
   real(GW_PRC) :: ubi(ncol,pver+1)

   ! Unit vectors of source wind (zonal and meridional components).
   real(GW_PRC) :: xv(ncol)
   real(GW_PRC) :: yv(ncol)

   ! Averages over source region.
   real(GW_PRC) :: ubmsrc(ncol) ! On-ridge wind.
   real(GW_PRC) :: usrc(ncol)   ! Zonal wind.
   real(GW_PRC) :: vsrc(ncol)   ! Meridional wind.
   real(GW_PRC) :: nsrc(ncol)   ! B-V frequency.
   real(GW_PRC) :: rsrc(ncol)   ! Density.

   ! Heating depth [m] and maximum heating in each column.
   real(GW_PRC) :: hdepth(ncol), maxq0(ncol)

   ! Effective gravity wave diffusivity at interfaces.
   real(GW_PRC) :: egwdffi(ncol,pver+1)

   ! Temperature tendencies from diffusion and kinetic energy.
   real(GW_PRC) :: dttdf(ncol,pver)
   real(GW_PRC) :: dttke(ncol,pver)

   real(GW_PRC) :: pint_adj(ncol,pver+1)
   real(GW_PRC) :: zfac_layer

   logical, parameter :: gw_apply_tndmax = .TRUE.!- default .TRUE. for Anisotropic: "Sean" limiters

   real(GW_PRC) :: zlb,pm,rhom,cmu,fpmx,fpmy,fe,fpe,fpml,fpel,fpmt,fpet,dusrcl,dvsrcl,dtsrcl
   real(GW_PRC) :: utfac,uhtmax
   integer  :: ktop

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
        ! sprectrum source index
        if (pref(k+1) < desc%spectrum_source) desc%k = k+1
     end do

     ! Determine wave sources for Beres deep scheme
     call gw_beres_src(ncol, pver, band , desc, &
          u, v, netdt, zm, src_level, tend_level, tau, &
          ubm, ubi, xv, yv, c, hdepth, maxq0)

!WMP pressure scaling from GEOS top 0.01mb to zfac_layer
     pint_adj = 1.0
     zfac_layer = 100.0 ! 1.0mb
     where (pint < zfac_layer)
       pint_adj = 1./19. * &
                  ((atan( (2.*(pint-1.0)/(zfac_layer-1.0)-1.) * &
                  tan(20.*PI/21.-0.5*PI) ) + 0.5*PI) * 21./PI - 1.)
     endwhere
!WMP pressure scaling from GEOS

     ! Solve for the drag profile with orographic sources.
     call gw_drag_prof(ncol, pver, band, pint, delp, rdelp, & 
          src_level, tend_level, dt, t,    &
          piln, rhoi, nm, ni, ubm, ubi, xv, yv,   &
          effgw, c, kvtt, tau, utgw, vtgw, &
          ttgw, egwdffi, gwut, dttdf, dttke, tau_adjust=pint_adj)

! GEOS efficiency and energy/momentum adjustments
  ktop=1
  do i=1,ncol
! Calculate launch level height
    zlb = 0.
    do k = ktop+1, pver
       if (k >= desc%k+1) then
! Define layer pressure and density
          pm   = (pint(i,k-1)+pint(i,k))*0.5
          rhom = pm/(MAPL_RGAS*t(i,k))
          zlb  = zlb + delp(i,k)/MAPL_GRAV/rhom
       end if
    end do
   !-----------------------------------------------------------------------
   ! Calculates energy and momentum flux profiles
   !-----------------------------------------------------------------------
    do l = -band%ngwv, band%ngwv
       do k = ktop, pver
          if ( k <= desc%k ) then
             cmu  = c(i,l)-ubi(i,k)
             fpmx =        sign(1.0,cmu)*tau(i,l,k)*xv(i)*effgw(i)
             fpmy =        sign(1.0,cmu)*tau(i,l,k)*yv(i)*effgw(i)
             fe   =    cmu*sign(1.0,cmu)*tau(i,l,k)      *effgw(i)
             fpe  = c(i,l)*sign(1.0,cmu)*tau(i,l,k)      *effgw(i)
             if (k == desc%k) then
                fpml = fpmx*xv(i)+fpmy*yv(i)
                fpel = fpe
             end if
             if (k == ktop) then
                fpmt = fpmx*xv(i)+fpmy*yv(i)
                fpet = fpe
             end if
          end if
          if (k >= desc%k+1) then
! Define layer pressure and density
             pm   = (pint(i,k-1)+pint(i,k))*0.5
             rhom = pm/(MAPL_RGAS*t(i,k))
! Compute sub-source tendencies
             dusrcl = - (fpml-fpmt)/(rhom*zlb)*xv(i)
             dvsrcl = - (fpml-fpmt)/(rhom*zlb)*yv(i)
             dtsrcl = -((fpel-fpet)-ubm(i,k)*(fpml-fpmt))/  &
                              (rhom*zlb*MAPL_CP)
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
    if (uhtmax > TNDMAX) utfac = TNDMAX/uhtmax
    do k = ktop, pver
       utgw(i,k) = utgw(i,k)*utfac
       vtgw(i,k) = vtgw(i,k)*utfac
       ttgw(i,k) = ttgw(i,k)*utfac
    end do

  end do  ! i=1,ncol

   flx_heat(:ncol) = 0.0

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
  real(GW_PRC), intent(in) :: x(:)
  real        , intent(in) :: grid(:)

  integer :: idx(size(x))

  real(GW_PRC) :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2.0_GW_PRC

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

end module gw_convect
