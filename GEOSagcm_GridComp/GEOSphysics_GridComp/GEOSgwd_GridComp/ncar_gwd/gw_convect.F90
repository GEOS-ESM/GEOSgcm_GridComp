module gw_convect

!
! This module handles gravity waves from convection, and was extracted from
! gw_drag in May 2013.
!

  use gw_utils, only: GW_PRC, GW_R8, get_unit_vector, dot_2d, midpoint_interp
  use gw_common, only: GWBand, gw_drag_prof, tau_0_ubc_cnv, tau_0_ubc_frt, &
                       energy_momentum_adjust, gw_flux_diagnostics
  use MAPL_Constants, only: MAPL_RGAS, MAPL_CP, MAPL_GRAV

implicit none
private

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
   ! Heating depths below this value [m] will be ignored.
   real :: min_hdepth
   ! Source for wave spectrum
   real :: spectrum_source
   ! Index for level where wind speed is used as the source speed.
   integer, allocatable :: k(:)
   ! tendency limiter
   real :: tndmax
   ! Table bounds, for convenience. (Could be inferred from shape(mfcc).)
   integer :: maxh
   integer :: maxuh
   ! Heating depths [m].
   real, allocatable :: hd(:)
   ! Convective heating rate conversion factor
   real :: hr_cf
   ! Scaling factor for generating QBO
   real :: qbo_hdepth_scaling
   ! Table of source spectra.
   real, allocatable :: mfcc(:,:,:)
   ! Forced background for extratropics
   real, allocatable :: taubck(:,:)
   ! Efficiency TR:ET function
   real, allocatable :: effbck(:)
   real :: et_bkg_dtdtm_forcing
   real :: et_bkg_speed_forcing
end type BeresSourceDesc


contains

!==========================================================================

!------------------------------------
subroutine gw_beres_init (file_name, band, desc, pgwv, gw_dc, ew_crit_thresh, ww_crit_thresh, fcrit2, wavelength, &
                          spectrum_source, hr_cf, qbo_hdepth_scaling, min_hdepth, storm_shift, eff_tr, eff_et, &
                          tau_bkg, et_fac_dtdtm, et_fac_speed, tndmax, &
                          active, ncol, lats)
#include <netcdf.inc>

  character(len=*), intent(in) :: file_name
  type(GWBand), intent(inout) :: band

  type(BeresSourceDesc), intent(inout) :: desc

  integer, intent(in) :: pgwv, ncol
  real, intent(in) :: gw_dc, ew_crit_thresh, ww_crit_thresh, fcrit2, wavelength
  real, intent(in) :: spectrum_source, hr_cf, qbo_hdepth_scaling, min_hdepth, eff_tr, eff_et, tau_bkg, tndmax
  logical, intent(in) :: storm_shift, active
  real, intent(in) :: et_fac_dtdtm, et_fac_speed, lats(ncol)

  ! Stuff for Beres convective gravity wave source.
  real(GW_R8), allocatable :: mfcc(:,:,:), hdcc(:)
  integer  :: hd_mfcc , mw_mfcc, ps_mfcc, ngwv_file, ps_mfcc_mid

  ! For forced background extratropical wave speed
  real    :: latdeg, flat_gw
  real, allocatable :: cw(:)
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

  band  = GWBand(pgwv, gw_dc, ew_crit_thresh, ww_crit_thresh, fcrit2, wavelength )

  ! These dimensions; {HD,MW,PS}_MFCC, came from Beres forcing file.

  ! Get HD (heating depth) dimension.
  desc%maxh = HD_MFCC

  ! Get MW (mean wind) dimension.
  desc%maxuh = MW_MFCC

  ! Get PS (phase speed) dimension.
  ngwv_file = PS_MFCC

  ! Number in each direction is half of total (and minus phase speed of 0).
  desc%maxuh = (desc%maxuh-1)/2

  ! midpoint of spectrum in netcdf file is ps_mfcc (odd number) divided by 2, plus 1
  ! E.g., ps_mfcc = 81. So, ps_mfcc_mid = 41
  !       1   11  21  31 32 33 34 35 36 37 38 39 40 41 42 43 ... 
  !      -40 -30 -20 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1  0 +1 +2 ...
  ps_mfcc_mid= INT(ngwv_file/2) + 1

  desc%active = active
  if (active) then

    allocate(desc%hd(desc%maxh) , stat=status )

    allocate(desc%mfcc(desc%maxh,-desc%maxuh:desc%maxuh,-band%ngwv:band%ngwv), stat=status )

    desc%mfcc( : , -desc%maxuh:desc%maxuh , -band%ngwv            :band%ngwv             ) & 
       = mfcc( :,             :           , -band%ngwv+ps_mfcc_mid:band%ngwv+ps_mfcc_mid )
  
    ! While not currently documented in the file, it uses kilometers. Convert
    ! to meters.
    desc%hd = hdcc * 1000.0

    ! Source level index allocated, filled later
    desc%spectrum_source = spectrum_source
    allocate(desc%k(ncol))

    desc%hr_cf = hr_cf

    desc%qbo_hdepth_scaling = qbo_hdepth_scaling

    desc%min_hdepth = min_hdepth

    desc%storm_shift = storm_shift

    desc%tndmax = tndmax

    ! Intialize forced background wave speeds
    allocate(desc%effbck(ncol))
    allocate(desc%taubck(ncol,-band%ngwv:band%ngwv))
    allocate(cw(-band%ngwv:band%ngwv))
    desc%effbck = 1.0
    desc%taubck = 0.0
    cw  = 0.0
    ! Create Gaussian weights using actual band%cref values
    do kc = -band%ngwv,band%ngwv
       cw(kc) =  exp(-(band%cref(kc)/25.)**2)
    enddo
    desc%et_bkg_dtdtm_forcing = et_fac_dtdtm
    desc%et_bkg_speed_forcing = et_fac_speed
    do i=1,ncol
      ! include forced background stress in extra tropics
      ! Determine the background stress at c=0
       if (desc%et_bkg_dtdtm_forcing /= 0.0 .or. desc%et_bkg_speed_forcing /= 0.0) then
          flat_gw = 0.05 ! weak background forcing
          desc%taubck(i,:) = tau_bkg*0.001*flat_gw*cw
         ! efficiency function
          desc%effbck(i) = eff_tr*cos(lats(i))**2 + &
                           eff_et*sin(lats(i))**2
       else
         ! Include dependence on latitude:
          latdeg = lats(i)*rad2deg
          if (ABS(latdeg) <  60.) then
            flat_gw =  max(0.15,0.50*exp(-((abs(latdeg)-60.)/23.)**2))
          elseif (ABS(latdeg) >= 60.) then
            flat_gw =           0.50*exp(-((abs(latdeg)-60.)/70.)**2)
          endif
          desc%taubck(i,:) = tau_bkg*0.001*flat_gw*cw
         ! efficiency function
          desc%effbck(i) = eff_tr*cos(lats(i))**2 + &
                           eff_et*sin(lats(i))**2
       endif
    enddo
    deallocate( cw )
  end if
    
end subroutine gw_beres_init

!------------------------------------
subroutine gw_beres_src(ncol, pver, band, desc, pint, u, v, &
     netdt, zm, src_level, tend_level, tau, tau_0_ubc, ubm, ubi, xv, yv, &
     bkg_tau, bkg_tau_cnv, bkg_tau_dry, bkg_tau_mst, &
     c, dtdtm, speed)
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

!------------------------------Arguments--------------------------------
  ! Column and vertical dimensions.
  integer, intent(in) :: ncol, pver

  ! Wavelengths triggered by convection.
  type(GWBand), intent(in) :: band

  ! Settings for convection type (e.g. deep vs shallow).
  type(BeresSourceDesc), intent(inout) :: desc

  ! Edge pressures
  real, intent(in) :: pint(ncol,pver+1)
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

  ! background wave stress forcings
  real, intent(out) :: bkg_tau(ncol)
  real, intent(out) :: bkg_tau_cnv(ncol)
  real, intent(out) :: bkg_tau_dry(ncol)
  real, intent(out) :: bkg_tau_mst(ncol)

  ! Wave Reynolds stress.
  real(GW_PRC), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real, intent(out) :: ubm(ncol,pver)
  real, intent(out) :: ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(out) :: xv(ncol), yv(ncol)
  ! Phase speeds.
  real(GW_PRC), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Frontal and Jet proxy inputs
  real, intent(in) :: dtdtm(ncol,pver)  ! Microphysics temperature tendency / latent heating (K s-1)
  real, intent(in) :: speed(ncol)       ! Katabatic proxy: Max wind speed in lowest 300m stable layer (m s-1)

  ! tau_0_ubc column dependence
  real(GW_PRC), intent(out) :: tau_0_ubc(ncol)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real :: uconv(ncol), vconv(ncol), ubi1d(ncol)

  ! Heating depth [m] and maximum heating in each column.
  real(GW_PRC) :: hdepth(ncol)

  ! Maximum heating rate.
  real(GW_PRC) :: q0(ncol)
  real(GW_PRC) :: moist_mult(ncol)
  real(GW_PRC) :: dry_mult, phys_mult

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
  integer :: thread
  integer :: k_ceiling, k_max, k_300m
  real :: target_press
  ! Define a critical heating rate threshold (e.g., 1 K/day = 1.157e-5 K/s)
  real, parameter :: dtdtm_critical = 1.157e-5  ! 1 K/day in K/s

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------
  tau = 0.0
  hdepth = 0.0
  q0 = 0.0
  tau0 = 0.0
  ubi = 0.0
  bkg_tau = 0.0
  bkg_tau_cnv = 0.0
  bkg_tau_dry = 0.0
  bkg_tau_mst = 0.0

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
  hdepth = hdepth*desc%qbo_hdepth_scaling

  hd_idx = index_of_nearest(hdepth, desc%hd)

  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though.
  where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0

  ! =========================================================================
  ! NEW CONFIGURATION: DYNAMIC SOURCE & RESOLUTION DECOUPLING BLOCK
  ! =========================================================================

  do i = 1, ncol
     if (hd_idx(i) > 0) then
        ! -------------------------------------------------------------------
        ! A. DEEP CONVECTIVE SOURCE REGIME (Preserves spectrum_source Baseline)
        ! -------------------------------------------------------------------
        q0(i) = 0.0
        do k = topi(i), boti(i)
           if (netdt(i,k) > q0(i)) q0(i) = netdt(i,k)
        end do

        ! Apply the legacy resolution scale awareness factor (hr_cf)
        q0(i) = q0(i) * desc%hr_cf

        ! Rigidly anchor convective launch level to spectrum_source
        desc%k(i) = 1
        do k = 0, pver-2
           if (pint(i,k+1) < desc%spectrum_source) desc%k(i) = k+1
        end do

        tau_0_ubc(i) = tau_0_ubc_cnv
     else
        ! -------------------------------------------------------------------
        ! B. SHALLOW/FRONTAL SOURCE REGIME (Frontal Masking & Launch Level)
        ! -------------------------------------------------------------------

        ! Identify an index for a safe physical ceiling for frontal scheme
        k_ceiling = 1
        do k = 1, pver
           if (pint(i,k+1) >= desc%spectrum_source) then
               k_ceiling = k
              exit
           end if
        end do

        ! Find 300m launch level (dry GWD default)
        k_300m = pver
        target_press = pint(i,pver+1) - 3000.0  ! 30 hPa above the surface, ~300m depth
        do k = pver, k_ceiling, -1
           if (pint(i,k+1) <= target_press) then
               k_300m = k
               exit
           endif
        end do

        q0(i) = dtdtm_critical
        k_max = k_300m ! Default to 300m for dry GWD

        if (desc%et_bkg_dtdtm_forcing /= 0.0) then
           ! Scan upward from surface to  desc%k(i) to find shallow high-lat peaks
           do k = pver,  k_ceiling, -1  
              if (dtdtm(i,k) > q0(i)) then
                 q0(i) = dtdtm(i,k)
                 k_max = k ! Capture dynamic launch layer
              endif
           end do
           ! Scale moist multiplier based on instantaneous tropospheric heating rates
           ! Cap at 1-10x, then scale by tuning factor
           moist_mult(i) = MAX(1.0, MIN(10.0, q0(i) * 86400.0) * desc%et_bkg_dtdtm_forcing)
        else
           moist_mult(i) = 1.0
        endif

        ! Update dynamic launch allocations for frontal waves
        desc%k(i) = k_max
        topi(i)   = k_max
        boti(i)   = pver ! Anchor bottom array to surface

        tau_0_ubc(i) = tau_0_ubc_frt
     endif
  end do

  ! =========================================================================
  ! DOWNSTREAM WIND CALCULATIONS (Now fully responsive to dynamic desc%k)
  ! =========================================================================
  ! Source wind speed and direction.
  do i=1,ncol
   uconv(i) = u(i,desc%k(i))
   vconv(i) = v(i,desc%k(i))
  enddo

  ! Get the unit vector components and magnitude at the source level.
  ubi1d = 0.0
  call get_unit_vector(uconv, vconv, xv, yv, ubi1d)
  do i=1,ncol
   ! SAFEGUARD CAP: Ensure interface indexing does not overshoot pver bounds
   ubi(i, min(pver, desc%k(i)+1)) = ubi1d(i)
  enddo

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)
  ubi(:,2:pver) = midpoint_interp(ubm)

  uh = 0.0
  do i=1,ncol
     if (desc%storm_shift .and. (hd_idx(i) > 0)) then
         ! Average wind in heating region, relative to storm cells.
          do k = topi(i), boti(i)
             uh(i) = uh(i) + ubm(i,k)
          end do
          uh(i) = uh(i)/(boti(i)-topi(i)+1)
         ! Find the cell speed where the storm speed is > 10 m/s.
         ! Storm speed is taken to be the source wind speed.
          CS(i) = sign(max(abs(ubm(i,desc%k(i)))-10.0, 0.0), ubm(i,desc%k(i)))
          uh(i) = uh(i) - CS(i)
     else
         ! For shallow convection, wind is relative to ground, and "heating
         ! region" wind is just the source level wind.
          uh(i) = ubm(i,desc%k(i))
     endif
  enddo

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

     if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------
        tau0 = desc%mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc%storm_shift) then
           shift = -nint(CS(i)/band%dc)
           tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude. (q0 already accounts for hr_cf here)
        tau0 = tau0*q0(i)*q0(i)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0
 
        tau(i,:,topi(i)+1) = tau0

        ! export background tau
        bkg_tau(i) = maxval(tau0)
        bkg_tau_cnv(i) = maxval(tau0)

     else

        tau(i,:,:) = 0.0
        if (desc%et_bkg_dtdtm_forcing /= 0.0 .or. desc%et_bkg_speed_forcing /= 0.0) then
          ! -----------------------------------------------------------------
          ! Frontal Detection via Physical Multipliers (Calculated Above)
          ! -----------------------------------------------------------------
          ! Proxy 2: The Dry Wind (Katabatic winds - evaluated locally)
           if (desc%et_bkg_speed_forcing /= 0.0) then
               ! A baseline 5 m/s wind yields a 1.0x multiplier (no extra drag).
               ! A linear 5 to 25 m/s ramp (1 - 10)x capped
               ! Scale by tuning factor: 0.5 => 1:5x, 1.0 => 1:10x, 2.0 => 1:20x
               dry_mult = MAX(1.0, MIN(10.0, (1.0 + (speed(i) - 5.0) * (9.0 / 20.0)) * desc%et_bkg_speed_forcing))
           else
               dry_mult = 1.0
           endif
           ! export dry bkg tau
           bkg_tau_dry(i) = maxval(desc%taubck(i,:)) * dry_mult
           ! export background tau
           bkg_tau_mst(i) = maxval(desc%taubck(i,:)) * moist_mult(i)
           ! Find the dominant physical forcing mechanism
           phys_mult = MAX(moist_mult(i), dry_mult)
           ! Apply physical forcing
           tau(i,:,desc%k(i)+1) = desc%taubck(i,:) * phys_mult
           ! export background tau
           bkg_tau(i) = maxval(desc%taubck(i,:)) * phys_mult
        else
           ! Fallback background stress assignment
           tau(i,:,desc%k(i)+1) = desc%taubck(i,:)
           ! export background tau
           bkg_tau(i) = maxval(desc%taubck(i,:))
        endif

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

subroutine gw_beres_ifc(band, ncol, pver, dt, effgw_dp, &
     u, v, t, pref, pint, delp, rdelp, piln, &
     zm, zi, nm, ni, rhoi, kvtt, &
     netdt, desc, alpha, &
     bkg_tau, bkg_tau_cnv, bkg_tau_dry, bkg_tau_mst, &
     utgw, vtgw, ttgw, flx_heat, dtdtm, speed, &
     taugwx, taugwy, fegw, fepgw, &
     taugwx_east, taugwx_west, taugwy_east, taugwy_west, &
     fegw_east, fegw_west, fepgw_east, fepgw_west, &
     taugwx_sfc, taugwy_sfc)

  !-----------------------------------------------------------------------
  ! Interface routine for Beres gravity wave drag parameterization
  !
  ! This routine orchestrates the complete GWD calculation:
  ! 1. Determine wave sources from convection (gw_beres_src)
  ! 2. Propagate waves and compute drag profiles (gw_drag_prof)
  ! 3. Apply efficiency and stability constraints (energy_momentum_adjust)
  ! 4. Compute diagnostic fluxes (gw_flux_diagnostics)
  !
  ! References:
  ! Beres et al. (2004): "A method of specifying the gravity wave spectrum
  ! above convection based on latent heating properties and background wind"
  ! J. Atmos. Sci., Vol 61, No. 3, pp. 324-337.
  !-----------------------------------------------------------------------

  !------------------------------Arguments--------------------------------
  ! Configuration and dimensions
  type(GWBand), intent(in) :: band
  type(BeresSourceDesc), intent(inout) :: desc
  integer, intent(in) :: ncol                ! Number of atmospheric columns
  integer, intent(in) :: pver                ! Number of vertical layers
  real, intent(in) :: dt                     ! Time step (s)
  real, intent(in) :: effgw_dp               ! Deep convection GW efficiency

  ! Wind fields
  real, intent(in) :: u(ncol,pver)           ! Zonal wind (m/s)
  real, intent(in) :: v(ncol,pver)           ! Meridional wind (m/s)

  ! Thermodynamic fields
  real, intent(in) :: t(ncol,pver)           ! Temperature (K)
  real, intent(in) :: netdt(ncol,pver)       ! Convective heating rate (K/s)
  real, intent(in) :: dtdtm(ncol,pver)       ! Microphysics heating rate (K/s)

  ! Pressure coordinates
  real, intent(in) :: pref(pver+1)           ! Reference pressure at interfaces (Pa)
  real, intent(in) :: pint(ncol,pver+1)      ! Interface pressures (Pa)
  real, intent(in) :: piln(ncol,pver+1)      ! Log of interface pressures
  real, intent(in) :: delp(ncol,pver)        ! Layer pressure thickness (Pa)
  real, intent(in) :: rdelp(ncol,pver)       ! Inverse pressure thickness (Pa^-1)

  ! Altitude coordinates
  real, intent(in) :: zm(ncol,pver)          ! Midpoint altitudes (m)
  real, intent(in) :: zi(ncol,pver+1)        ! Interface altitudes (m)

  ! Stability parameters
  real, intent(in) :: nm(ncol,pver)          ! Brunt-Vaisala frequency at midpoints (s^-1)
  real, intent(in) :: ni(ncol,pver+1)        ! Brunt-Vaisala frequency at interfaces (s^-1)

  ! Density and diffusivity
  real, intent(in) :: rhoi(ncol,pver+1)      ! Interface density (kg/m^3)
  real, intent(in) :: kvtt(ncol,pver+1)      ! Molecular thermal diffusivity (m^2/s)

  ! Vertical structure
  real, intent(in) :: alpha(:)                ! Vertical damping profile

  ! Surface forcing proxy
  real, intent(in) :: speed(ncol)             ! Max wind in stable surface layer (m/s)

  !------ Output: Wind and temperature tendencies ------
  real, intent(out) :: utgw(ncol,pver)       ! Zonal wind tendency (m/s^2)
  real, intent(out) :: vtgw(ncol,pver)       ! Meridional wind tendency (m/s^2)
  real, intent(out) :: ttgw(ncol,pver)       ! Temperature tendency (K/s)

  !------ Output: Background stress diagnostics ------
  real, intent(out) :: bkg_tau(ncol)         ! Total background stress (Pa)
  real, intent(out) :: bkg_tau_cnv(ncol)     ! Convective source stress (Pa)
  real, intent(out) :: bkg_tau_dry(ncol)     ! Dry/katabatic source stress (Pa)
  real, intent(out) :: bkg_tau_mst(ncol)     ! Moist source stress (Pa)

  !------ Output: Energy diagnostics ------
  real, intent(inout) :: flx_heat(ncol)      ! Energy change (J/m^2)

  !------ Output: Momentum and energy flux diagnostics ------
  ! Total fluxes
  real, intent(out) :: taugwx(ncol,pver)     ! Zonal momentum flux (Pa)
  real, intent(out) :: taugwy(ncol,pver)     ! Meridional momentum flux (Pa)
  real, intent(out) :: fegw(ncol,pver)       ! Kinetic energy flux (W/m^2)
  real, intent(out) :: fepgw(ncol,pver)      ! Phase-speed weighted energy flux (W/m^2)

  ! Eastward-propagating waves (c > 0)
  real, intent(out) :: taugwx_east(ncol,pver)  ! Zonal momentum flux (Pa)
  real, intent(out) :: taugwy_east(ncol,pver)  ! Meridional momentum flux (Pa)
  real, intent(out) :: fegw_east(ncol,pver)    ! Kinetic energy flux (W/m^2)
  real, intent(out) :: fepgw_east(ncol,pver)   ! Phase-speed weighted energy flux (W/m^2)

  ! Westward-propagating waves (c < 0)
  real, intent(out) :: taugwx_west(ncol,pver)  ! Zonal momentum flux (Pa)
  real, intent(out) :: taugwy_west(ncol,pver)  ! Meridional momentum flux (Pa)
  real, intent(out) :: fegw_west(ncol,pver)    ! Kinetic energy flux (W/m^2)
  real, intent(out) :: fepgw_west(ncol,pver)   ! Phase-speed weighted energy flux (W/m^2)

  !------ Output: Surface stress ------
  real, intent(out) :: taugwx_sfc(ncol)          ! Zonal surface stress (Pa)
  real, intent(out) :: taugwy_sfc(ncol)          ! Meridional surface stress (Pa)

  !---------------------------Local Storage-------------------------------
  ! Wave stress and phase speeds
  real(GW_PRC), allocatable :: tau(:,:,:)    ! Wave Reynolds stress (Pa)
  real(GW_PRC), allocatable :: c(:,:)        ! Phase speeds (m/s)
  real(GW_PRC), allocatable :: gwut(:,:,:)   ! Wind tendency per wave band

  ! Wind projections
  real :: ubm(ncol,pver)                     ! Wind projection at midpoints
  real :: ubi(ncol,pver+1)                   ! Wind projection at interfaces
  real :: xv(ncol), yv(ncol)                 ! Unit vectors of source wind

  ! Source and tendency levels
  integer :: src_level(ncol)                 ! Top gravity wave source level
  integer :: tend_level(ncol)                ! Lowest level for wind tendencies

  ! Efficiency and stress
  real :: effgw(ncol)                        ! GW momentum transfer efficiency
  real :: tau_0_ubc(ncol)                    ! Upper boundary condition for stress

  integer :: i, l

  !-----------------------------------------------------------------------
  ! Allocate working arrays
  !-----------------------------------------------------------------------
  allocate(tau(ncol,-band%ngwv:band%ngwv,pver+1))
  allocate(c(ncol,-band%ngwv:band%ngwv))
  allocate(gwut(ncol,pver,-band%ngwv:band%ngwv))

  !-----------------------------------------------------------------------
  ! Step 1: Determine wave sources from convection
  !-----------------------------------------------------------------------
  call gw_beres_src(ncol, pver, band, desc, pint, &
       u, v, netdt, zm, src_level, tend_level, tau, &
       tau_0_ubc, ubm, ubi, xv, yv, &
       bkg_tau, bkg_tau_cnv, bkg_tau_dry, bkg_tau_mst, &
       c, dtdtm=dtdtm, speed=speed)

  !-----------------------------------------------------------------------
  ! Step 2: Propagate waves and compute drag profiles
  !-----------------------------------------------------------------------
  call gw_drag_prof(ncol, pver, band, pint, delp, rdelp, &
       src_level, tend_level, dt, t, &
       piln, rhoi, nm, ni, ubm, ubi, xv, yv, &
       c, kvtt, tau, tau_0_ubc, utgw, vtgw, ttgw, gwut, alpha)

  !-----------------------------------------------------------------------
  ! Step 3: Apply efficiency scaling and stability constraints
  !-----------------------------------------------------------------------
  effgw = effgw_dp * desc%effbck

  call energy_momentum_adjust(ncol, pver, band, pint, delp, u, v, dt, c, tau, &
       effgw, t, ubm, ubi, xv, yv, utgw, vtgw, ttgw, &
       tend_level, tndmax_in=desc%tndmax)

  !-----------------------------------------------------------------------
  ! Step 4: Compute diagnostic fluxes (after efficiency applied)
  !-----------------------------------------------------------------------
  call gw_flux_diagnostics(ncol, pver, band, c, ubi, tau, xv, yv, &
       src_level, tend_level, pint, &
       taugwx, taugwy, fegw, fepgw, &
       taugwx_east, taugwx_west, taugwy_east, taugwy_west, &
       fegw_east, fegw_west, fepgw_east, fepgw_west)

  !-----------------------------------------------------------------------
  ! Compute surface stress (at bottom interface, k=pver+1)
  !-----------------------------------------------------------------------
  taugwx_sfc = 0.0
  taugwy_sfc = 0.0
  do i = 1, ncol
     do l = -band%ngwv, band%ngwv
        ! Surface stress from this wave
        taugwx_sfc(i) = taugwx_sfc(i) + real(tau(i,l,pver+1)) * xv(i)
        taugwy_sfc(i) = taugwy_sfc(i) + real(tau(i,l,pver+1)) * yv(i)
     enddo
  enddo

  !-----------------------------------------------------------------------
  ! Clean up
  !-----------------------------------------------------------------------
  deallocate(tau, c, gwut)

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
