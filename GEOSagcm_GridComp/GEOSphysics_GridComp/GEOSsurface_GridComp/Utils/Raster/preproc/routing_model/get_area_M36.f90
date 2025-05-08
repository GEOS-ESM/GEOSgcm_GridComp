program main
!Main purpose: Gets the area for each catchment-tile for M36 grid.

use river_read ! Use the module "river_read" to access functions for reading NetCDF files
use constant, only: nmax=>nmax36,nc,nlon,nlat,nlat36,nlon36,nt_global=>nt_global36

implicit none
! Require explicit declaration of all variables

! Declare variables for indices and temporary storage
integer             :: id, xi, yi, i, j, flag, subi, x_m36, y_m36, it

! Declare allocatable arrays for sub-catchment information
integer,allocatable :: xsub(:,:), ysub(:,:), subi_global(:,:)
! Arrays for storing sub-catchment coordinate indices and sub-catchment areas
real,allocatable    :: asub(:,:)

! Declare arrays for grid and mapping information
! Arrays for longitude and latitude values from the NetCDF file
real*8,allocatable  :: lon(:), lat(:)
! Arrays for integer mappings of longitude and latitude indices
integer,allocatable :: loni(:), lati(:)
! 2D arrays: "catchind" holds catchment index for each original grid cell; "map_tile" maps aggregated grid cells to global indices
integer,allocatable :: catchind(:,:), map_tile(:,:)
! "cellarea" holds the area of each grid cell; "area_m36" is the aggregated area on the M36 grid;
! "area_tile" will store the area for each global tile based on the aggregated grid
real,allocatable    :: cellarea(:,:), area_m36(:,:), area_tile(:)

! Define file path for input routing data:
character(len=900)  :: file_path !"input/CatchIndex.nc"   

if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
endif
call get_command_argument(1, file_path)

! Allocate arrays for sub-catchment information with dimensions (nmax, nc)
allocate(xsub(nmax,nc), ysub(nmax,nc), asub(nmax,nc))
! Allocate arrays for the original grid: catchment indices and cell areas
allocate(catchind(nlon,nlat), cellarea(nlon,nlat))
! Allocate arrays for longitude and latitude values
allocate(lon(nlon), lat(nlat))
! Allocate arrays for the mapping of longitude and latitude indices
allocate(loni(nlon), lati(nlat))

! Read longitude and latitude data from the NetCDF file
call read_ncfile_double1d(trim(file_path), "longitude", lon, nlon)
call read_ncfile_double1d(trim(file_path), "latitude", lat, nlat)
! Read the 2D catchment index data 
call read_ncfile_int2d(trim(file_path), "CatchIndex", catchind, nlon, nlat)
! Read the cell area data 
call read_ncfile_real2d("temp/cellarea.nc", "data", cellarea, nlon, nlat)
! Convert cell areas (from m^2 to km^2) by scaling with 1.e6
cellarea = cellarea / 1.e6

! Open text files to read the mapping indices for the aggregated grid (M36)
! Read the latitude mapping: converts original latitude indices to M36 grid indices
open(10, file="temp/lati_1m_M36.txt")
read(10, *) lati
! Read the longitude mapping: converts original longitude indices to M36 grid indices
open(11, file="temp/loni_1m_M36.txt")
read(11, *) loni

! Allocate and initialize the aggregated area array for the M36 grid
allocate(area_m36(nlon36, nlat36))
area_m36 = 0.
! Loop over each grid cell in the original grid
do xi = 1, nlon
  do yi = 1, nlat
    if (catchind(xi,yi) >= 1) then
      ! For cells that belong to a catchment (valid catchment index)
      x_m36 = loni(xi)
      y_m36 = lati(yi)
      ! Accumulate the cell area into the corresponding aggregated grid cell
      area_m36(x_m36, y_m36) = area_m36(x_m36, y_m36) + cellarea(xi,yi)
    endif
  enddo
enddo

! Allocate the map_tile array and read its data from a NetCDF file
allocate(map_tile(nlon36, nlat36))
call read_ncfile_int2d("temp/map_tile_M36.nc", "data", map_tile, nlon36, nlat36)
! Allocate the global area array to hold the area for each global tile
allocate(area_tile(nt_global))
area_tile = -9999.

! Map the aggregated grid areas to the global tile indices using the map_tile mapping
do i = 1, nlon36
  do j = 1, nlat36
    it = map_tile(i, j)
    if (it > 0) then
      area_tile(it) = area_m36(i, j)
    endif
  enddo
enddo

! Write the global tile area data to an output text file
open(88, file="output/area_M36_1d.txt")
do i = 1, nt_global
  write(88, *) area_tile(i)
enddo

end