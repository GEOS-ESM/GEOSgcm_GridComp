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

public :: BeresSourceDesc
public :: gw_beres_ifc
! public :: gw_beres_src
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
   real, allocatable :: k(:)
   ! tendency limiter
   real :: tndmax
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
  spectrum_source, min_hdepth, storm_shift, taubgnd, tndmax, active, ncol, lats)
#include <netcdf.inc>

   character(len=*), intent(in) :: file_name
   type(GWBand), intent(inout) :: band
 
   type(BeresSourceDesc), intent(inout) :: desc
 
   integer, intent(in) :: pgwv, ncol
   real, intent(in) :: gw_dc, fcrit2, wavelength
   real, intent(in) :: spectrum_source, min_hdepth, taubgnd, tndmax
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
 
     ! Source level index allocated, filled later
     desc%spectrum_source = spectrum_source
     allocate(desc%k(ncol))
 
     desc%min_hdepth = min_hdepth
 
     desc%storm_shift = storm_shift
 
     desc%tndmax = tndmax
 
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
! subroutine gw_beres_src(ncol, pver, band, desc, u, v, &
!      netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
!      c, hdepth, maxq0)
! !-----------------------------------------------------------------------
! ! Driver for multiple gravity wave drag parameterization.
! !
! ! The parameterization is assumed to operate only where water vapor
! ! concentrations are negligible in determining the density.
! !
! ! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! ! specifying the gravity wave spectrum above convection based on latent
! ! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! ! pp. 324-337.
! !
! !-----------------------------------------------------------------------
!   !!!use gw_common, only: GWBand, qbo_hdepth_scaling

! !------------------------------Arguments--------------------------------
!   ! Column and vertical dimensions.
!   integer, intent(in) :: ncol, pver

!   ! Wavelengths triggered by convection.
!   type(GWBand), intent(in) :: band

!   ! Settings for convection type (e.g. deep vs shallow).
!   type(BeresSourceDesc), intent(in) :: desc

!   ! Midpoint zonal/meridional winds.
!   real(r8), intent(in) :: u(ncol,pver), v(ncol,pver)
!   ! Heating rate due to convection.
!   real(r8), intent(in) :: netdt(:,:)
!   ! Midpoint altitudes.
!   real(r8), intent(in) :: zm(ncol,pver)

!   ! Indices of top gravity wave source level and lowest level where wind
!   ! tendencies are allowed.
!   integer, intent(out) :: src_level(ncol)
!   integer, intent(out) :: tend_level(ncol)

!   ! Wave Reynolds stress.
!   real(r8), intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
!   ! Projection of wind at midpoints and interfaces.
!   real(r8), intent(out) :: ubm(ncol,pver), ubi(ncol,pver+1)
!   ! Unit vectors of source wind (zonal and meridional components).
!   real(r8), intent(out) :: xv(ncol), yv(ncol)
!   ! Phase speeds.
!   real(r8), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

!   ! Heating depth [m] and maximum heating in each column.
!   real(r8), intent(out) :: hdepth(ncol), maxq0(ncol)

! !---------------------------Local Storage-------------------------------
!   ! Column and level indices.
!   integer :: i, k

!   integer :: l, topi_sum, topi_index, boti_index

!   ! Zonal/meridional wind at roughly the level where the convection occurs.
!   real(r8) :: uconv(ncol), vconv(ncol)

!   ! Maximum heating rate.
!   real(r8) :: q0(ncol)

!   ! Bottom/top heating range index.
!   integer  :: boti(ncol), topi(ncol)
!   ! Index for looking up heating depth dimension in the table.
!   integer  :: hd_idx(ncol)
!   ! Mean wind in heating region.
!   real(r8) :: uh(ncol)
!   ! Min/max wavenumber for critical level filtering.
!   integer :: Umini(ncol), Umaxi(ncol)
!   ! Source level tau for a column.
!   real(r8) :: tau0(-band%ngwv:band%ngwv, ncol)
!   ! Speed of convective cells relative to storm.
!   real(r8) :: CS(ncol)
!   ! Index to shift spectra relative to ground.
!   integer :: shift

!   ! Heating rate conversion factor.
!   real(r8), parameter :: CF = 20._r8
!   ! Averaging length.
!   real(r8), parameter :: AL = 1.0e5_r8

! !!$acc data present(desc, desc%hd, desc%mfcc, desc%k, u, v, tau, ubi, xv, yv, hdepth,  band, band%ngwv, &
! !$acc data present(desc, u, v, tau, ubi, xv, yv, hdepth, band, &
! !$acc              ubm, zm, netdt, maxq0, src_level, tend_level, c) &
! !$acc       create(uconv, vconv, q0, boti, topi, hd_idx, uh, Umini, Umaxi, tau0, CS)
! !$acc data present(desc%hd, desc%mfcc, desc%k, band%ngwv, band%cref, band%dc)

! !$acc parallel
!   !----------------------------------------------------------------------
!   ! Initialize tau array
!   !----------------------------------------------------------------------

! !$acc loop gang vector
!    do i = 1,ncol
!       hdepth(i) = 0.0_r8
!       q0(i) = 0.0_r8

!       ! Source wind speed and direction.
!       ! uconv(i) = u(i,desc%k)
!       ! vconv(i) = v(i,desc%k)
!    enddo
! !$acc end parallel

! !$acc parallel
! !$acc loop gang vector collapse(3)
!    do k = 1,pver+1
!       do l = -band%ngwv,band%ngwv
!          do i = 1,ncol
!             tau(i,l,k) = 0.0_r8
!          enddo
!       enddo
!    enddo
! !$acc end parallel

! !$acc parallel
! !$acc loop gang vector collapse(2)
!    do i = 1, ncol
!       do l = -band%ngwv,band%ngwv
!          tau0(l,i) = 0.0_r8
!       enddo
!    enddo
! !$acc end parallel

!   !------------------------------------------------------------------------
!   ! Determine wind and unit vectors approximately at the source level, then
!   ! project winds.
!   !------------------------------------------------------------------------

! !$acc parallel
!   ! Get the unit vector components and magnitude at the source level.
! !$acc loop gang vector
!    do i = 1,ncol
!       ubi(i,desc%k+1) = sqrt(u(i,desc%k)*u(i,desc%k) + v(i,desc%k)*v(i,desc%k))

!       if(ubi(i,desc%k+1) > 0._r8) then
!          xv(i) = u(i,desc%k)/ubi(i,desc%k+1)
!          yv(i) = v(i,desc%k)/ubi(i,desc%k+1)
!       else
!          xv(i) = 0._r8
!          yv(i) = 0._r8
!       endif
!    enddo
! !$acc end parallel

! !$acc parallel

!   ! Project the local wind at midpoints onto the source wind.
! !$acc loop gang vector collapse(2)
!   do k = 1, pver
!       do i = 1,ncol
!          ubm(i,k) = u(i,k) * xv(i) + v(i,k) * yv(i)
!       enddo
!   end do

! !$acc end parallel

! !$acc parallel

!   ! Compute the interface wind projection by averaging the midpoint winds.
!   ! Use the top level wind at the top interface.
! !$acc loop gang vector
!    do i = 1,ncol
!       ubi(i,1) = ubm(i,1)
!    enddo
! !$acc end parallel

! !$acc parallel

! !$acc loop gang vector collapse(2)
!    do k = 2,pver
!       do i = 1,ncol
!          ubi(i,k) = 0.5_r8 * (ubm(i,k-1)+ubm(i,k))
!       enddo
!    enddo
! !$acc end parallel
!   !-----------------------------------------------------------------------
!   ! Calculate heating depth.
!   !
!   ! Heating depth is defined as the first height range from the bottom in
!   ! which heating rate is continuously positive.
!   !-----------------------------------------------------------------------

!   ! First find the indices for the top and bottom of the heating range.
! !$acc parallel
! !$acc loop gang vector
!    do i = 1,ncol
!       boti(i) = 0
!       topi(i) = 0
!    enddo
! !$acc loop seq
!    do k = pver, 1, -1
!       topi_sum = 0
! !$acc loop vector reduction(+ : topi_sum)
!       do i = 1,ncol
!          if (topi(i).eq.0) then
!             topi_sum = topi_sum + 1
!          endif
!       enddo
!       ! if (all(topi /= 0).eq..FALSE.) then
!       if (topi_sum.gt.0) then
! !$acc loop gang vector
!          do i = 1, ncol
!             if (boti(i) == 0) then
!                ! Detect if we are outside the maximum range (where z = 20 km).
!                if (zm(i,k) >= 20000._r8) then
!                   boti(i) = k
!                   topi(i) = k
!                else
!                   ! First spot where heating rate is positive.
!                   if (netdt(i,k) > 0.0_r8) boti(i) = k
!                end if
!             else if (topi(i) == 0) then
!                ! Detect if we are outside the maximum range (z = 20 km).
!                if (zm(i,k) >= 20000._r8) then
!                   topi(i) = k
!                else
!                   ! First spot where heating rate is no longer positive.
!                   if (.not. (netdt(i,k) > 0.0_r8)) topi(i) = k
!                end if
!             end if
!          end do
!       endif
!   end do
! !$acc end parallel

! !$acc parallel
!   ! Heating depth in m.
! !$acc loop gang vector
!    do i = 1,ncol
!       hdepth(i) = zm(i,topi(i))-zm(i,boti(i))*qbo_hdepth_scaling
!    enddo

!   ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
! !   hdepth = hdepth*qbo_hdepth_scaling
!   ! CK : Incorporated qbo_hdepth_scaling into above loop

! !   hd_idx = index_of_nearest(hdepth, desc%hd)
!    call index_of_nearest_s(hdepth, desc%hd, hd_idx)

!   ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
!   ! either not big enough for the lowest table entry, or it is below the
!   ! minimum allowed for this convection type.
!   ! Values above the max in the table still get the highest value, though.
! !   where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0

! !$acc end parallel

! !$acc kernels
!       topi_index = minval(topi)
!       boti_index = maxval(boti)
! !$acc end kernels

! !$acc parallel

! !$acc loop gang vector
!   do i = 1, ncol
!       if(hdepth(i) < max(desc%min_hdepth, desc%hd(1))) then
!          hd_idx(i) = 0
!       endif
!   enddo
! !$acc end parallel
!   ! Maximum heating rate.

! !$acc parallel
! ! NOTE : Check with NVIDIA about appropriate use of reduction in OpenACC
! !        and whether this is an approprate use case to apply reduction to q0

! !$acc loop gang vector collapse(2)
!   do k = topi_index, boti_index
!       do i = 1, ncol
!          ! if (k >= topi(i) .and. k <= boti(i)) then
! !$acc atomic update         
!             q0(i) = max(q0(i), netdt(i,k))
!          ! endif
!       enddo
!   end do
! !$acc end parallel

! !$acc parallel
! !$acc loop gang vector
!    do i = 1,ncol
!       ! output max heating rate in K/day
!       maxq0(i) = q0(i)*24._r8*3600._r8
!       ! Multipy by conversion factor
!       q0(i) = q0(i) * CF
!    enddo
! !$acc end parallel

!   if (desc%storm_shift) then
! !$acc parallel
!      ! Find the cell speed where the storm speed is > 10 m/s.
!      ! Storm speed is taken to be the source wind speed.
! !$acc loop gang vector
!       do i = 1,ncol
!          CS(i) = sign(max(abs(ubm(i,desc%k))-10._r8, 0._r8), ubm(i,desc%k))
!          ! Average wind in heating region, relative to storm cells.
!          uh(i) = 0._r8
!       enddo
! !$acc end parallel
! ! NOTE : This is a similar section to check whether reduction is used properly

! !$acc parallel
! !$acc loop gang vector collapse(2)
!      do k = topi_index, boti_index
!       do i = 1,ncol
!       !   where (k >= topi .and. k <= boti)
! !$acc atomic update
!            uh(i) = uh(i) + ubm(i,k)/(boti(i)-topi(i)+1)
!       !   end where
!         enddo
!      end do
! !$acc end parallel

! !$acc parallel
! !$acc loop gang vector
!       do i = 1,ncol
!          uh(i) = uh(i) - CS(i)
!       enddo
! !$acc end parallel
!   else
!      ! For shallow convection, wind is relative to ground, and "heating
!      ! region" wind is just the source level wind.
! !$acc parallel
! !$acc loop gang vector
!       do i = 1,ncol
!          uh(i) = ubm(i,desc%k)
!       enddo
! !$acc end parallel
!   end if

! !$acc parallel
! !$acc loop gang vector
!    do i = 1,ncol
!   ! Limit uh to table range.
!       uh(i) = max(min(uh(i), real(desc%maxuh, r8)), -real(desc%maxuh, r8))

!   ! Speeds for critical level filtering.
!       Umini(i) =  band%ngwv
!       Umaxi(i) = -band%ngwv   
!    enddo
! !$acc end parallel

! !$acc parallel
! !$acc loop gang vector collapse(2)
!   do k = topi_index, boti_index
!       do i = 1,ncol
!          ! if(k >= topi(i) .and. k <= boti(i)) then
! !$acc atomic update
!             Umini(i) = min(Umini(i), nint(ubm(i,k)/band%dc))
! !$acc atomic update
!             Umaxi(i) = max(Umaxi(i), nint(ubm(i,k)/band%dc))
!          ! endif
!      enddo
!   end do
! !$acc end parallel

! !$acc parallel
! !$acc loop gang vector
!    do i = 1,ncol
!       Umini(i) = max(Umini(i), -band%ngwv)
!       Umaxi(i) = min(Umaxi(i), band%ngwv)
!    enddo
! !$acc end parallel

! !$acc parallel
!   !-----------------------------------------------------------------------
!   ! Gravity wave sources
!   !-----------------------------------------------------------------------
!   ! Start loop over all columns.
!   !-----------------------------------------------------------------------
! !$acc loop gang private(shift)
!   do i=1,ncol

!      !---------------------------------------------------------------------
!      ! Look up spectrum only if the heating depth is large enough, else set
!      ! tau0 = 0.
!      !---------------------------------------------------------------------

!      if (hd_idx(i) > 0) then
!         !------------------------------------------------------------------
!         ! Look up the spectrum using depth and uh while adjusting magnitude.
!         !------------------------------------------------------------------
! !$acc loop vector
!         do l = -band%ngwv,band%ngwv
!          tau0(l, i) = desc%mfcc(hd_idx(i),nint(uh(i)),l) * q0(i)*q0(i)/AL
!         enddo

!         if (desc%storm_shift) then
!            ! For deep convection, the wind was relative to storm cells, so
!            ! shift the spectrum so that it is now relative to the ground.
!            shift = -nint(CS(i)/band%dc)
!            tau0(:,i) = eoshift(tau0(:,i), shift)
!         end if

!         ! Adjust magnitude.
!       !   tau0(:,i) = tau0(:,i)*q0(i)*q0(i)/AL

!         ! Adjust for critical level filtering.
!         tau0(Umini(i):Umaxi(i), i) = 0.0_r8
 
!         tau(i,:,topi(i)+1) = tau0(:,i)

!      end if ! heating depth above min and not at the pole

!   enddo
! !$acc end parallel
!   !-----------------------------------------------------------------------
!   ! End loop over all columns.
!   !-----------------------------------------------------------------------
!   ! Output the source level.
! !$acc parallel
! !$acc loop gang vector
!    do i = 1,ncol
!       src_level(i) = topi(i)
!       tend_level(i) = topi(i)
!    enddo

! !$acc end parallel

! !$acc kernels
!   ! Set phase speeds; just use reference speeds.
!   c = spread(band%cref, 1, ncol)
! !$acc end kernels
! !$acc end data
! !$acc end data

! end subroutine gw_beres_src



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Main Interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gw_beres_ifc( band, &
   ncol, pver, dt, effgw_dp,  &
   u, v, t, pref, pint, delp, rdelp, piln, &
   zm, zi, nm, ni, rhoi, kvtt,  &
   netdt,desc,lats, alpha, &
   utgw,vtgw,ttgw,flx_heat,dqcdt)

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
   real,         intent(in) :: nm(ncol,pver)     ! Midpoint Brunt-Vaisalla frequencies (s-1).
   real,         intent(in) :: ni(ncol,pver+1)   ! Interface Brunt-Vaisalla frequencies (s-1).
   real,         intent(in) :: rhoi(ncol,pver+1) ! Interface density (kg m-3).
   real,         intent(in) :: kvtt(ncol,pver+1) ! Molecular thermal diffusivity.

   real,         intent(in) :: lats(ncol)      ! latitudes
   real,         intent(in) :: alpha(:)

   real,         intent(out) :: utgw(ncol,pver)       ! zonal wind tendency
   real,         intent(out) :: vtgw(ncol,pver)       ! meridional wind tendency
   real,         intent(out) :: ttgw(ncol,pver)       ! temperature tendency
   real,         intent(inout) :: flx_heat(ncol)        ! Energy change

   real, optional, intent(in) :: dqcdt(ncol,pver)  ! Condensate tendency due to large-scale (kg kg-1 s-1)

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

!    ! Allocate wavenumber fields.
!    allocate(tau(ncol,-band%ngwv:band%ngwv,pver+1))
!    allocate(gwut(ncol,pver,-band%ngwv:band%ngwv))
!    allocate(c(ncol,-band%ngwv:band%ngwv))

! !!$acc data present(effgw, pref, desc, desc%k, lats, flx_heat) &
! !$acc data present(pref, desc, lats, flx_heat) &
! !$acc       create(tau, gwut, c, effgw, src_level, tend_level, ubm, ubi, &
! !$acc              xv, yv, ubmsrc, usrc, vsrc, nsrc, rsrc, tauoro, taudsw, &
! !$acc              wbr, hdepth, maxq0, egwdffi, dttdf, dttke, &
! !$acc              taurx, taurx0, taury, taury0, de)
! !$acc data present(desc%k)

!      gw_apply_tndmax  = .FALSE.

! !$acc parallel
!      ! Efficiency of gravity wave momentum transfer.
!      ! This is really only to remove the pole points.
!    !   where (pi/2._r8 - abs(lats(:ncol)) >= 1.e-4 )  !-4*epsilon(1._r8))
!    !      effgw = effgw_dp
!    !   elsewhere
!    !      effgw = 0._r8
!    !   end where
! !$acc loop gang vector
!       do i = 1,ncol
!          if((pi/2._r8 - abs(lats(i))) >= 1.e-4) then
!             effgw(i) = effgw_dp
!          else
!             effgw(i) = 0._r8
!          endif
!       enddo

! !$acc loop seq
!      do k = 0, pver
!         ! 700 hPa index
!         if (pref(k+1) < 70000._r8) desc%k = k+1
!      end do
! !$acc end parallel

!      ! Determine wave sources for Beres deep scheme
!      call gw_beres_src(ncol, pver, band , desc, &
!           u, v, netdt, zm, src_level, tend_level, tau, &
!           ubm, ubi, xv, yv, c, hdepth, maxq0)

!      ! satfac_in is 2 by default for CAM5

!      ! Solve for the drag profile with orographic sources.
!      call gw_drag_prof(ncol, pver, band, pint, delp, rdelp, & 
!           src_level, tend_level,   dt, t,    &
!           piln, rhoi,       nm,   ni, ubm,  ubi,  xv,    yv,   &
!           effgw,c,          kvtt,  tau,  utgw,  vtgw, &
!           ttgw, egwdffi,  gwut, dttdf, dttke,            &
!           satfac_in = 1._r8,                                   &
!           lapply_effgw_in=gw_apply_tndmax)

!      ! For orographic waves, don't bother with taucd, since there are no
!      ! momentum conservation routines or directional diagnostics.

!      !  add the diffusion coefficients
!      !do k = 1, pver+1
!      !   egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
!      !end do

!      ! Add the orographic tendencies to the spectrum tendencies.
!      ! Don't calculate fixers, since we are too close to the ground to
!      ! spread momentum/energy differences across low layers.
! !++jtb 3/2020
!      !do k = 1, pver
!      !   ptend%u(:ncol,k) = ptend%u(:ncol,k) + utgw(:,k)
!      !   ptend%v(:ncol,k) = ptend%v(:ncol,k) + vtgw(:,k)
!      !   ptend%s(:ncol,k) = ptend%s(:ncol,k) + ttgw(:,k)
!      !end do

!      ! Calculate energy change for output to CAM's energy checker.
!      ! This is sort of cheating; we don't have a good a priori idea of the
!      ! energy coming from surface stress, so we just integrate what we and
!      ! actually have so far and overwrite flx_heat with that.
! !++jtb 3/2020
!      ! call energy_change(dt, p, u, v, ptend%u(:ncol,:), &
!      !     ptend%v(:ncol,:), ptend%s(:ncol,:), de)
!      !flx_heat(:ncol) = de
! !$acc parallel
! !$acc loop gang vector
!       do i = 1,ncol
!          flx_heat(i) = 0._r8
!       enddo
! !$acc end parallel

!      !do m = 1, pcnst
!      !   do k = 1, pver
!      !      ptend%q(:ncol,k,m) = ptend%q(:ncol,k,m) + qtgw(:,k,m)
!      !   end do
!      !end do
! !$acc end data
! !$acc end data

!    deallocate(tau, gwut, c)

end subroutine gw_beres_ifc


! !************************************************************************
! !!handle_err
! !************************************************************************
! !
! !!ROUTINE:      handle_err
! !!DESCRIPTION:  error handler
! !--------------------------------------------------------------------------

! subroutine handle_err(status)
  
!   implicit         none
  
! #include <netcdf.inc>
  
!   integer          status
  
!   if (status .ne. nf_noerr) then
!     print *, nf_strerror(status)
!     stop 'Stopped'
!   endif
  
! end subroutine handle_err



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine endrun(msg)

!    integer :: iulog

!    character(len=*), intent(in), optional :: msg    ! string to be printed

!     iulog=6

!    if (present (msg)) then
!       write(iulog,*)'ENDRUN:', msg
!    else
!       write(iulog,*)'ENDRUN: called without a message string'
!    end if

!    stop

! end subroutine endrun

! ! Short routine to get the indices of a set of values rounded to their
! ! nearest points on a grid.
! function index_of_nearest(x, grid) result(idx)
! !$acc routine seq
!   real(r8), intent(in) :: x(:)
!   real(r8), intent(in) :: grid(:)

!   integer :: idx(size(x))

!   real(r8) :: interfaces(size(grid)-1)
!   integer :: i, n

!   n = size(grid)
!   interfaces = (grid(:n-1) + grid(2:))/2._r8

!   idx = 1
!   do i = 1, n-1
!      where (x > interfaces(i)) idx = i + 1
!   end do

! end function index_of_nearest

! subroutine index_of_nearest_s(x, grid, idx)
! !$acc routine vector
!   real(r8), intent(in) :: x(:)
!   real(r8), intent(in) :: grid(:)
!   integer,  intent(out) :: idx(:)

!   integer :: i
! !$acc loop vector
!   do i = 1, size(grid)-1
!    idx(i) = 1
!    if (x(i) > ((grid(i) + grid(i+1))/2._r8)) then
!       idx(i) = i+1
!    endif
!   enddo
! !$acc end loop

! end subroutine

subroutine handle_err(status)
  
   implicit         none
   
#include <netcdf.inc>
   
   integer          status
   
   if (status .ne. nf_noerr) then
     print *, nf_strerror(status)
     stop 'Stopped'
   endif
   
 end subroutine handle_err

end module gw_convect

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”

! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at

! http://www.apache.org/licenses/LICENSE-2.0

! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.
