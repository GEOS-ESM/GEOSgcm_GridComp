program main
!Main purpose: Assigns a catchment‐tile index from catchment definition files to each model grid cell for M36 grid.

use river_read
use constant, only : nmax=>nmax36, nc, nlon, nlat

implicit none


! Variable declarations:
integer              :: id, xi, yi, i, flag, subi
integer              :: nsub(nc)                   ! Array to store the number of sub-areas for each catchment
integer, allocatable :: xsub(:,:), ysub(:,:), subi_global(:,:)
! xsub and ysub: 2D arrays to store mapped x and y indices for sub-catchments (not using subi_global in this code)
real, allocatable    :: asub(:,:)        ! 2D array to store aggregated area for each sub-catchment

real*8, allocatable  :: lon(:), lat(:)  ! Arrays to hold longitude and latitude values from the NetCDF file
integer, allocatable :: loni(:), lati(:)
! loni and lati: Arrays holding mapping indices from 1-minute resolution data files
integer, allocatable :: catchind(:,:) ! 2D array holding catchment indices for each grid cell
real, allocatable    :: cellarea(:,:)    ! 2D array containing the area of each grid cell

! Define file path for input routing data:
character(len=900)   :: file_path !"input/CatchIndex.nc"   

if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
endif
call get_command_argument(1, file_path)

! Allocate arrays with the specified dimensions:
allocate(xsub(nmax, nc), ysub(nmax, nc), asub(nmax, nc))
allocate(catchind(nlon, nlat), cellarea(nlon, nlat))
allocate(lon(nlon), lat(nlat))
allocate(loni(nlon), lati(nlat))

! Read grid information from the NetCDF file "CatchIndex.nc":
call read_ncfile_double1d(trim(file_path), "longitude", lon, nlon)
call read_ncfile_double1d(trim(file_path), "latitude", lat, nlat)
call read_ncfile_int2d(trim(file_path), "CatchIndex", catchind, nlon, nlat)
! Read cell area data from the NetCDF file "cellarea.nc":
call read_ncfile_real2d("temp/cellarea.nc", "data", cellarea, nlon, nlat)
cellarea = cellarea / 1.e6  ! Convert cell area (e.g., from m^2 to km^2)

! Read mapping indices for the 1-minute resolution grid from text files:
open(10, file="temp/lati_1m_M36.txt")
read(10, *) lati
open(11, file="temp/loni_1m_M36.txt")
read(11, *) loni

! Initialize aggregation arrays to zero:
nsub = 0
xsub = 0
ysub = 0
asub = 0.

! Loop over all grid cells to aggregate cell areas by catchment and sub-area:
do xi = 1, nlon
  do yi = 1, nlat
    if (catchind(xi, yi) >= 1) then
      ! The grid cell belongs to a catchment
      id = catchind(xi, yi)  ! Get the catchment id for the current cell
      flag = 0              ! Reset flag to indicate whether a matching sub-area is found
      
      ! If the catchment already has one or more sub-areas, check for a matching sub-area:
      if (nsub(id) >= 1) then
        do i = 1, nsub(id)
          if (loni(xi) == xsub(i, id) .and. lati(yi) == ysub(i, id)) then
            flag = 1
            ! If a match is found, accumulate the cell area into the existing sub-area:
            asub(i, id) = asub(i, id) + cellarea(xi, yi)
            exit  ! Exit the inner loop since a matching sub-area has been found
          endif
        end do
      endif
      
      ! If no matching sub-area was found, create a new sub-area:
      if (flag == 0) then
        nsub(id) = nsub(id) + 1
        xsub(nsub(id), id) = loni(xi)
        ysub(nsub(id), id) = lati(yi)
        asub(nsub(id), id) = cellarea(xi, yi)
      endif

    endif
  end do
end do

! Open output files to write the aggregated sub-catchment information:
open(50, file="output/Pfaf_nsub_M36.txt")
open(51, file="output/Pfaf_xsub_M36.txt")
open(52, file="output/Pfaf_ysub_M36.txt")
open(53, file="output/Pfaf_asub_M36.txt")
! Loop over all catchments and write:
do i = 1, nc
  write(50, *) nsub(i)  ! Write the number of sub-areas for catchment i
  write(51, '(150(1x,i3))') xsub(:, i)  ! Write the x indices for all sub-areas (formatted as 3-digit integers)
  write(52, '(150(1x,i3))') ysub(:, i)  ! Write the y indices for all sub-areas (formatted as 3-digit integers)
  write(53, '(150(1x,f10.4))') asub(:, i)  ! Write the aggregated areas for all sub-areas (formatted as floating-point numbers)
end do

! Print the maximum number of sub-areas found for any catchment and its location:
print *, maxval(nsub)
print *, maxloc(nsub)

end
