program main
!Main purpose: Assigns a catchment‐tile index from maptile files to each sub-catchment for M09 grid.

use river_read        ! Use custom module for reading NetCDF files
use constant, only: nlat=>nlat09, nlon=>nlon09, nmax=>nmax09, nc
implicit none

! Declare allocatable arrays for grid mapping and sub-catchment indices
integer,allocatable :: map_tile(:,:), subx(:,:), suby(:,:), subi(:,:)

! Declare integer variables for looping and temporary storage
integer             :: i, x, y, j, it

! Allocate the map_tile array for the aggregated grid (M09) with dimensions (nlon, nlat)
allocate(map_tile(nlon,nlat))
! Read the mapping data from a NetCDF file into the map_tile array
call read_ncfile_int2d("temp/map_tile_M09.nc", "data", map_tile, nlon, nlat)

! Allocate subx, suby, and subi arrays to store sub-catchment coordinate data and indices
allocate(subx(nmax,nc), suby(nmax,nc), subi(nmax,nc))

! Open and read the x-coordinates of sub-catchments from a text file into subx
open(77, file="output/Pfaf_xsub_M09.txt")
read(77, *) subx

! Open and read the y-coordinates of sub-catchments from a text file into suby
open(77, file="output/Pfaf_ysub_M09.txt")
read(77, *) suby

! Initialize the subi array to zero
subi = 0

! Loop over each catchment
do i = 1, nc
  ! Loop over each possible sub-area within a catchment
  do j = 1, nmax
    x = subx(j, i)
    y = suby(j, i)
    ! If the x-coordinate is non-zero, then the sub-area exists
    if (x /= 0) then
      ! If x exists but y is zero, then there is an error and the program stops
      if (y == 0) stop
      ! Map the sub-area indices from the aggregated grid using the map_tile array
      subi(j, i) = map_tile(x, y)
    endif
  enddo
enddo

! Open an output file to write the computed sub-catchment tile indices
open(88, file="output/Pfaf_isub_M09.txt")
do i = 1, nc
  write(88, '(150(i8))') subi(:, i)
enddo

end program main
