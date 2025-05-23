program main
!Main purpose: Assigns a catchment‐tile index from maptile files to each sub-catchment for M36 grid.

use river_read         ! Use custom module for reading NetCDF files
use constant, only: nlat=>nlat36, nlon=>nlon36, nmax=>nmax36, nc
implicit none

! Declare allocatable arrays for grid mapping and sub-catchment indices
integer,allocatable :: map_tile(:,:), subx(:,:), suby(:,:), subi(:,:)

! Declare integer variables for loop indices and temporary storage
integer             :: i, x, y, j, it

! Allocate the mapping array for the aggregated grid with dimensions (nlon, nlat)
allocate(map_tile(nlon, nlat))
! Read mapping data from the NetCDF file into the map_tile array
call read_ncfile_int2d("temp/map_tile_M36.nc", "data", map_tile, nlon, nlat)

! Allocate arrays to store sub-catchment x and y coordinates and their mapped indices
allocate(subx(nmax, nc), suby(nmax, nc), subi(nmax, nc))

! Open and read the sub-catchment x-coordinates from a text file into the subx array
open(77, file="output/Pfaf_xsub_M36.txt")
read(77, *) subx

! Open and read the sub-catchment y-coordinates from a text file into the suby array
open(77, file="output/Pfaf_ysub_M36.txt")
read(77, *) suby

! Initialize the sub-area index array to zero
subi = 0

! Loop over each catchment
do i = 1, nc
  ! Loop over each potential sub-area within the current catchment
  do j = 1, nmax
    x = subx(j, i)  ! Retrieve the x-coordinate for the sub-area
    y = suby(j, i)  ! Retrieve the y-coordinate for the sub-area
    if (x /= 0) then  ! Check if a valid sub-area exists (non-zero x-coordinate)
      if (y == 0) stop  ! If x is valid but y is zero, there is an error, so stop the program
      subi(j, i) = map_tile(x, y)  ! Map the sub-area coordinates to a global tile index using map_tile
    endif
  enddo
enddo

! Open an output file to write the computed sub-catchment tile indices
open(88, file="output/Pfaf_isub_M36.txt")
do i = 1, nc
  write(88, '(150(i7))') subi(:, i)
enddo

end program main

