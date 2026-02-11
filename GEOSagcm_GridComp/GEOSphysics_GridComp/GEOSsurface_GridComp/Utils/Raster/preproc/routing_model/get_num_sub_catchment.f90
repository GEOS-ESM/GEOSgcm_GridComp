program main
!Main purpose: Assigns a catchment‐tile index from catchment definition files to each model grid cell for M09 grid.

use river_read
use constant, only : nmax=>nmax09, nc, nlon, nlat

use EASE_pfaf_fracMod, only : EASE_find_subs

implicit none

! Variable declarations:
integer              :: id, xi, yi, i, flag, subi
integer              :: nsub(nc)       ! Array storing the number of sub-areas for each catchment

! Allocatable arrays for sub-catchment information:
integer, allocatable :: xsub(:,:), ysub(:,:)
real, allocatable    :: asub(:,:)   ! Aggregated area for each sub-catchment

! Arrays for grid data and mapping:
real*8, allocatable  :: lon(:), lat(:)      ! Longitude and latitude arrays from NetCDF file
integer, allocatable :: loni(:), lati(:)     ! Mapped integer indices from 1-minute resolution files
integer, allocatable :: catchind(:,:)        ! 2D array of catchment indices for each grid cell


! Define file path for input routing data:
character(len=900)   :: file_path !"input/CatchIndex.nc"   

if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
endif
call get_command_argument(1, file_path)

! Allocate arrays based on the defined dimensions:
allocate(xsub(nmax, nc), ysub(nmax, nc), asub(nmax, nc))
allocate(catchind(nlon, nlat))
allocate(lon(nlon), lat(nlat))
allocate(loni(nlon), lati(nlat))

! Read grid longitude, latitude, catchment index, and cell area data from NetCDF files:
call read_ncfile_double1d(trim(file_path), "longitude", lon, nlon)
call read_ncfile_double1d(trim(file_path), "latitude", lat, nlat)
call read_ncfile_int2d(trim(file_path), "CatchIndex", catchind, nlon, nlat)

! Read mapped grid indices for 1-minute resolution from text files:
open(10, file="temp/lati_1m_M09.txt")
read(10, *) lati
open(11, file="temp/loni_1m_M09.txt")
read(11, *) loni

! Initialize aggregation arrays:
nsub = 0         ! Set number of sub-areas per catchment to zero
xsub = 0         ! Initialize x-coordinate array for sub-catchments to zero
ysub = 0         ! Initialize y-coordinate array for sub-catchments to zero
asub = 0.        ! Initialize aggregated area values to zero

! Loop over each 1m grid cell to accumulate cell areas into sub-catchments:
call EASE_Find_subs(catchind, loni, lati, lat, nsub, xsub, ysub, asub)

! Open output files to write the aggregated sub-catchment information:
open(50, file="temp/Pfaf_nsub_M09.txt")
open(51, file="temp/Pfaf_xsub_M09.txt")
open(52, file="temp/Pfaf_ysub_M09.txt")
open(53, file="temp/Pfaf_asub_M09.txt")

! For each catchment, write:
!   - Number of sub-areas
!   - X indices of sub-areas (formatted in groups of 458 integers)
!   - Y indices of sub-areas (formatted similarly)
!   - Aggregated area values of sub-areas (formatted as floating-point numbers)
do i = 1, nc
  write(50, *) nsub(i)
  write(51, '(458(1x,i4))') xsub(:, i)
  write(52, '(458(1x,i4))') ysub(:, i)
  write(53, '(458(1x,f10.4))') asub(:, i)
end do

! Print the maximum number of sub-areas found for any catchment and its location:
print *, maxval(nsub)
print *, maxloc(nsub)

end
