program main
!Main purpose: Extracts the latitude/longitude boundaries of each catchment-tile from catchment definition files for M36 grid.

use constant, only : nt=>nt36

implicit none

! Declare allocatable arrays for catchment ID, parent catchment ID, and boundary coordinates
integer, allocatable, dimension(:) :: id, catid
real, allocatable, dimension(:)    :: lon_left, lon_right, lat_bottom, lat_top

integer :: i, ntot  ! 'i' is the loop counter; 'ntot' holds the total number of catchments read from the file

! Define file path for input routing data:
character(len=900) :: file_path !"input/catchment_M36.def"   

if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
endif
call get_command_argument(1, file_path)

! Allocate arrays with size nt
allocate(id(nt), catid(nt), lon_left(nt), lon_right(nt), lat_bottom(nt), lat_top(nt))

! Open the catchment definition file for the M36 grid and read the total number of catchments (header)
open(77, file=trim(file_path))
read(77, *) ntot

! Loop over each catchment and read: id, catchment id, left/right longitudes, bottom/top latitudes
do i = 1, nt
  read(77, *) id(i), catid(i), lon_left(i), lon_right(i), lat_bottom(i), lat_top(i)
end do

! Write the left boundary longitudes to an output file
open(88, file="temp/lon_left_M36.txt")
do i = 1, nt
  write(88, *) lon_left(i)
end do

! Write the right boundary longitudes to an output file
open(88, file="temp/lon_right_M36.txt")
do i = 1, nt
  write(88, *) lon_right(i)
end do

! Write the bottom boundary latitudes to an output file
open(88, file="temp/lat_bottom_M36.txt")
do i = 1, nt
  write(88, *) lat_bottom(i)
end do

! Write the top boundary latitudes to an output file
open(88, file="temp/lat_upper_M36.txt")
do i = 1, nt
  write(88, *) lat_top(i)
end do

end
