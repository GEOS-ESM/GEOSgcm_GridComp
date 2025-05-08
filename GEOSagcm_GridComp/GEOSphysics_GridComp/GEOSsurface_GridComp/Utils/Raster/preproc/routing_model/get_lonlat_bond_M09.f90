program main
!Main purpose: Extracts the latitude/longitude boundaries of each catchment-tile from catchment definition files for M09 grid.

use constant, only : nt=>nt09

implicit none

! Declare allocatable arrays for catchment ID, parent catchment ID,
! and geographical boundaries (longitude and latitude extents)
integer, allocatable, dimension(:) :: id, catid
real, allocatable, dimension(:)    :: lon_left, lon_right, lat_bottom, lat_top

integer :: i, ntot  ! Loop counter and total number of catchments read from file

! Define file path for input routing data:
character(len=900) :: file_path !"input/catchment_M09.def"   

if (command_argument_count() /= 1) then
    print *, "no <file_path> found"
    stop
endif
call get_command_argument(1, file_path)

! Allocate arrays with size nt
allocate(id(nt), catid(nt), lon_left(nt), lon_right(nt), lat_bottom(nt), lat_top(nt))

! Open input file that contains catchment definitions
open(77, file=trim(file_path))
! Read total number of catchments (ntot) from the file header
read(77, *) ntot

! Loop over each catchment and read the definitions:
! id, catchment id, left and right longitudes, bottom and top latitudes
do i = 1, nt
  read(77, *) id(i), catid(i), lon_left(i), lon_right(i), lat_bottom(i), lat_top(i)
end do

! Write the left longitude values to a temporary output file
open(88, file="temp/lon_left_M09.txt")
do i = 1, nt
  write(88, *) lon_left(i)
end do

! Write the right longitude values to a temporary output file
open(88, file="temp/lon_right_M09.txt")
do i = 1, nt
  write(88, *) lon_right(i)
end do

! Write the bottom latitude values to a temporary output file
open(88, file="temp/lat_bottom_M09.txt")
do i = 1, nt
  write(88, *) lat_bottom(i)
end do

! Write the upper (top) latitude values to a temporary output file
open(88, file="temp/lat_upper_M09.txt")
do i = 1, nt
  write(88, *) lat_top(i)
end do

end