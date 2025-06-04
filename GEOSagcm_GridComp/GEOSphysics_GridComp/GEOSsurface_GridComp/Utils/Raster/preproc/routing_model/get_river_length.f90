program main
!Main purpose: Determines main river channel lengths for each catchment by using HydroSHEDS data of distance to sink.

use river_read
use constant, only : nc, nlon, nlat, nlonh, nlath, cur_avg, cur_min, cur_max

implicit none

real*8, allocatable                :: lon(:), lat(:)
real, allocatable                  :: ldn1m(:,:), elev1m(:,:)
integer, allocatable               :: catid(:,:), flag_slp(:)

real*8, allocatable                :: lonh(:), lath(:)
real, allocatable                  :: ldnh(:,:), elev_15s(:,:)

! Declare arrays to hold routing and catchment characteristics:
real, allocatable, dimension(:)    :: lon_dn, lat_dn, lon_up, lat_up, dist_ref, dist_ref2, ldn_min, ldn_max, riv_len, str_len, slp
real, allocatable, dimension(:)    :: lon_min, lon_max, lat_min, lat_max, area, elevdiff_ref, elevdiff
integer, allocatable, dimension(:) :: xi_min, yi_min, xi_max, yi_max
integer, allocatable, dimension(:) :: downid

! Loop indices and temporary variables
integer xi, yi
integer                            :: num, i, j, cid, did, k
integer                            :: data1, data12
real*8                             :: data2
real                               :: data7, data9, data10
real                               :: elev_temp

character(len=100)                 :: file_pfafmap !input/SRTM_PfafData.nc
character(len=100)                 :: file_ldn !input/hyd_glo_ldn_15s.nc
character(len=100)                 :: file_hyelev !input/hyd_glo_dem_15s.nc
character(len=100)                 :: file_pfafrout !input/Pfafcatch-routing.dat

  if (command_argument_count() /= 4) then
      print *, "no <file_pfafmap> <file_ldn> <file_hyelev> <file_pfafrout> found"
      stop
  endif
  call get_command_argument(1, file_pfafmap)
  call get_command_argument(2, file_ldn)
  call get_command_argument(3, file_hyelev)
  call get_command_argument(4, file_pfafrout)

!-----------------------------------------------------------------------
! Regrid LDN (length to sink) from HydroSHEDS data

allocate(ldn1m(nlon, nlat), catid(nlon, nlat))
allocate(lon(nlon), lat(nlat))
! Read longitude, latitude, and catchment index data from SRTM Pfaf data
call read_ncfile_double1d(trim(file_pfafmap), "longitude", lon, nlon)
call read_ncfile_double1d(trim(file_pfafmap), "latitude", lat, nlat)
call read_ncfile_int2d(trim(file_pfafmap), "CatchIndex", catid, nlon, nlat)
ldn1m = -1.
where(catid == -9999) ldn1m = -9999.

! Allocate high-resolution LDN array and read data from HydroSHEDS 15s file
allocate(ldnh(nlonh, nlath))
call read_ncfile_real2d(trim(file_ldn), "Band1", ldnh, nlonh, nlath)
where(ldnh.lt.4.e9) ldnh = ldnh / 1.e3  ! Convert from meters to kilometers

! Regrid: For each grid cell in the M09 grid, assign the minimum LDN value from the corresponding high-res block.
do xi = 1, nlon
  do yi = 2041, 10440
    if (ldn1m(xi, yi) .ne. -9999.) then
      ldn1m(xi, yi) = minval(ldnh(4*xi-3:4*xi, 4*yi-3-8160:4*yi-8160))
      if (ldn1m(xi, yi) .gt. 4.e9) ldn1m(xi, yi) = -1.
    end if
  end do
end do
print *, maxval(ldn1m)

! Allocate arrays to store minimum and maximum LDN for each catchment and their corresponding grid indices
allocate(ldn_min(nc), ldn_max(nc), xi_min(nc), yi_min(nc), xi_max(nc), yi_max(nc))
ldn_min = 1.e20
ldn_max = -9999.
xi_min = -9999; yi_min = -9999; xi_max = -9999; yi_max = -9999
do i = 1, nlon
  do j = 1, nlat
    if (catid(i, j) >= 1) then
      cid = catid(i, j)
      if (ldn1m(i, j) > 0. .and. ldn1m(i, j) < ldn_min(cid)) then
        ldn_min(cid) = ldn1m(i, j)
        xi_min(cid) = i
        yi_min(cid) = j
      endif
      if (ldn1m(i, j) > 0. .and. ldn1m(i, j) > ldn_max(cid)) then
        ldn_max(cid) = ldn1m(i, j)
        xi_max(cid) = i
        yi_max(cid) = j        
      endif      
    endif
  end do
end do
where(ldn_min == 1.e20) ldn_min = -9999

!-----------------------------------------------------------------------
! Compute elevation at 1-minute resolution from high-resolution DEM (15s)
allocate(elev_15s(nlonh, nlath), elev1m(nlon, nlat))
call read_ncfile_real2d(trim(file_hyelev), "Band1", elev_15s, nlonh, nlath)
where(elev_15s > 30000.) elev_15s = 0.
elev1m = 0.
do xi = 1, nlon
  do yi = 2041, 10440
      elev1m(xi, yi) = sum(elev_15s(4*xi-3:4*xi, 4*yi-3-8160:4*yi-8160)) / 16.
  end do
end do


deallocate(ldnh, elev_15s)
!-----------------------------------------------------------------------
! Get reference distances using routing data

open(77, file=trim(file_pfafrout), form="formatted", status="old")
read(77, *) num
allocate(lon_dn(nc), lat_dn(nc), lon_up(nc), lat_up(nc), dist_ref(nc), dist_ref2(nc))
allocate(lon_min(nc), lon_max(nc), lat_min(nc), lat_max(nc), area(nc), elevdiff_ref(nc), elevdiff(nc))

! Read routing and catchment geometry data from the Pfafcatch routing file
do i = 1, nc
  read(77, *) data1, data2, lon_min(i), lon_max(i), lat_min(i), lat_max(i), data7, area(i), data9, data10, elevdiff_ref(i), data12, lon_dn(i), lat_dn(i), lon_up(i), lat_up(i)
end do

! Compute spherical distances reference
do i = 1, nc
  dist_ref(i) = spherical_distance(lon_dn(i), lat_dn(i), lon_up(i), lat_up(i))
  dist_ref2(i) = spherical_distance(lon_min(i), lat_min(i), lon_max(i), lat_max(i))
end do
where(dist_ref > dist_ref2 .or. dist_ref == 0.) dist_ref = 0.5 * dist_ref2

!--------------------------------------------------------------------
! Get initial guess of river length (riv_len) based on LDN differences and elevation differences

allocate(riv_len(nc), downid(nc), flag_slp(nc))
open(77, file="output/downstream_1D_new_noadj.txt")
read(77, *) downid

flag_slp = 1

riv_len = -9999.
elevdiff = -9999.
do i = 1, nc
  if (downid(i) >= 1) then
    did = downid(i)
    if (.not. (riv_len(did) >= cur_min * dist_ref(did) .and. riv_len(did) <= cur_max * dist_ref(did))) then
      riv_len(did) = ldn_min(i) - ldn_min(did)
      if (xi_min(i) > 0 .and. xi_min(did) > 0) then
        elevdiff(did) = max(0., elev1m(xi_min(i), yi_min(i)) - elev1m(xi_min(did), yi_min(did)))
        flag_slp(did) = 1
      else
        elevdiff(did) = elevdiff_ref(did)
        flag_slp(did) = 0
      endif
    else if (flag_slp(did) == 0 .or. elevdiff(did) == 0.) then
      riv_len(did) = ldn_min(i) - ldn_min(did)
      if (xi_min(i) > 0 .and. xi_min(did) > 0) then
        elevdiff(did) = max(0., elev1m(xi_min(i), yi_min(i)) - elev1m(xi_min(did), yi_min(did)))
        flag_slp(did) = 1
      else
        elevdiff(did) = elevdiff_ref(did)
        flag_slp(did) = 0
      endif
    endif
  endif
end do

do i = 1, nc
  if (riv_len(i) == -9999.) then
    riv_len(i) = (ldn_max(i) - ldn_min(i)) * 0.5
    if (xi_min(i) > 0) then
      elevdiff(i) = max(0., 0.5 * elev1m(xi_max(i), yi_max(i)) - 0.5 * elev1m(xi_min(i), yi_min(i)))
    else
      elevdiff(i) = elevdiff_ref(i)
      flag_slp(i) = 0
    endif
  endif
end do

k = 0
do i = 1, nc
  if (.not. (riv_len(i) >= cur_min * dist_ref(i) .and. riv_len(i) <= cur_max * dist_ref(i))) then
    riv_len(i) = cur_avg * dist_ref(i)
    elevdiff(i) = elevdiff_ref(i)
    flag_slp(i) = 0
    k = k + 1
  endif
end do
open(88, file="output/Pfaf_lriv_PR.txt")
do i = 1, nc
  write(88, *) riv_len(i)
end do

!--------------------------------------------------------------------
! Calculate the length scale of local streams based on catchment area and river length.
allocate(str_len(nc))
str_len = area / riv_len / 4. * cur_avg
open(88, file="output/Pfaf_lstr_PR.txt")
do i = 1, nc
  write(88, *) str_len(i)
end do
!--------------------------------------------------------------------
! Calculate the catchment slope from elevation difference and river length.
allocate(slp(nc))
slp = elevdiff * 1.e-3 / riv_len
where(slp.lt.1.e-5) flag_slp = 0
where(slp.lt.1.e-5) slp = 1.e-5
print *, sum(flag_slp)
open(88, file="temp/Pfaf_slope.txt")
do i = 1, nc
  write(88, *) slp(i)
end do
print *, minval(slp)
open(88, file="temp/Pfaf_slope_flag.txt")
do i = 1, nc
  write(88, *) flag_slp(i)
end do

!--------------------------------------------------------------------
contains

function spherical_distance(lon_dn, lat_dn, lon_up, lat_up) result(distance)
    implicit none
    !------------------------------------------------------------
    ! Function: spherical_distance
    ! Purpose : Calculates the great-circle distance between two geographic
    !           points using the Haversine formula.
    !
    ! Input:
    !   lon_dn, lat_dn - Longitude and latitude of the first point (degrees)
    !   lon_up, lat_up - Longitude and latitude of the second point (degrees)
    !
    ! Output:
    !   distance       - Great-circle distance between the two points (kilometers)
    !------------------------------------------------------------
    real, intent(in) :: lon_dn, lat_dn   ! Coordinates of downstream point
    real, intent(in) :: lon_up, lat_up     ! Coordinates of upstream point
    real :: distance                     ! Computed distance (km)
    real :: R, dlon, dlat, a, c          ! Intermediate variables

    R = 6371.0                         ! Earth's radius in kilometers
    dlon = (lon_up - lon_dn) * (acos(-1.0) / 180.0)  ! Delta longitude (radians)
    dlat = (lat_up - lat_dn) * (acos(-1.0) / 180.0)  ! Delta latitude (radians)

    a = sin(dlat / 2.0)**2 + cos(lat_dn * (acos(-1.0) / 180.0)) * &
        cos(lat_up * (acos(-1.0) / 180.0)) * sin(dlon / 2.0)**2
    c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
    distance = R * c

end function spherical_distance

end program main