program main
!Main purpose: Gets the area for each catchment-tile for M09 grid.

use river_read
use constant, only: nmax=>nmax09,nc,nlon,nlat,nlat09,nlon09,nt_global=>nt_global09

implicit none
! Require explicit declaration of all variables

! Declare variables for indices, flags, and temporary storage
integer             :: id, xi, yi, i, j, flag, subi, x_m09, y_m09, it
! Allocatable arrays to hold sub-catchment coordinate indices and global sub-catchment information
integer,allocatable :: xsub(:,:), ysub(:,:), subi_global(:,:)

! Allocatable array to store sub-catchment area data
real,allocatable    :: asub(:,:)

! Allocatable double precision arrays for storing longitude and latitude values from file
real*8,allocatable  :: lon(:), lat(:)
! Allocatable integer arrays for mapping longitude and latitude indices
integer,allocatable :: loni(:), lati(:)
! 2D arrays: catchind holds catchment index for each grid cell; map_tile maps M09 grid cells to global indices
integer,allocatable :: catchind(:,:), map_tile(:,:)
! Arrays for cell areas from the original grid, aggregated area on the M grid, and area per global tile
real,allocatable    :: cellarea(:,:), area_m09(:,:), area_tile(:)

! Define file path for input routing data:
character(len=900)  :: file_path !"input/CatchIndex.nc"   

if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
endif
call get_command_argument(1, file_path)

! Allocate arrays for sub-catchment data
! Allocate 2D arrays with dimensions (nmax, nc) for sub-catchment coordinate indices and areas
allocate(xsub(nmax,nc), ysub(nmax,nc), asub(nmax,nc))

! Allocate arrays for the catchment index grid and cell area data
allocate(catchind(nlon,nlat), cellarea(nlon,nlat))
! Allocate 1D arrays for longitude and latitude values
allocate(lon(nlon), lat(nlat))
! Allocate arrays for integer mappings of longitude and latitude indices
allocate(loni(nlon), lati(nlat))

! Read longitude and latitude data from the NetCDF file
call read_ncfile_double1d(trim(file_path), "longitude", lon, nlon)
call read_ncfile_double1d(trim(file_path), "latitude", lat, nlat)
! Read 2D catchment index data from the same file
call read_ncfile_int2d(trim(file_path), "CatchIndex", catchind, nlon, nlat)
! Read cell area data 
call read_ncfile_real2d("temp/cellarea.nc", "data", cellarea, nlon, nlat)
! Scale cell area values (from m^2 to km^2)
cellarea = cellarea/1.e6

! Read mapping indices from text files for the M09 grid conversion
! Read integer latitude indices mapping for each original grid 
open(10, file="temp/lati_1m_M09.txt")
read(10, *) lati
! Read integer longitude indices mapping for each original grid 
open(11, file="temp/loni_1m_M09.txt")
read(11, *) loni

! Allocate and initialize the aggregated area array for the M09 grid
allocate(area_m09(nlon09, nlat09))
area_m09 = 0.
! Loop over the original grid and accumulate cell areas into the M09 grid using mapping indices
do xi = 1, nlon
  do yi = 1, nlat
    if (catchind(xi,yi) >= 1) then
      x_m09 = loni(xi)
      y_m09 = lati(yi)
      ! For grid cells with a valid catchment index, add their cell area to the corresponding M09 grid cell      
      area_m09(x_m09, y_m09) = area_m09(x_m09, y_m09) + cellarea(xi,yi)
    endif
  enddo
enddo

! Allocate the map_tile array and read its data from a NetCDF file
allocate(map_tile(nlon09, nlat09))
call read_ncfile_int2d("temp/map_tile_M09.nc", "data", map_tile, nlon09, nlat09)
! Allocate the global area array to hold area data for each tile
allocate(area_tile(nt_global))
area_tile = -9999.

! Map the aggregated M09 grid areas to the global tile indices using the map_tile array
do i = 1, nlon09
  do j = 1, nlat09
    it = map_tile(i, j)
    if (it > 0) then    
      area_tile(it) = area_m09(i, j)
    endif
  enddo
enddo

! Write the global tile area data to an output text file
open(88, file="output/area_M09_1d.txt")
do i = 1, nt_global
  write(88, *) area_tile(i)
enddo

end
