program main

implicit none
integer,parameter :: nt=1684725

integer,allocatable,dimension(:) :: id, catid
real,allocatable,dimension(:) :: lon_left,lon_right,lat_bottom,lat_top

integer :: i,ntot

allocate(id(nt),catid(nt),lon_left(nt),lon_right(nt),lat_bottom(nt),lat_top(nt))
open(77,file="input/catchment_M09.def")
read(77,*)ntot
do i=1,nt
  read(77,*)id(i),catid(i),lon_left(i),lon_right(i),lat_bottom(i),lat_top(i)
enddo

open(88,file="temp/lon_left_M09.txt")
do i=1,nt
  write(88,*)lon_left(i)
enddo
open(88,file="temp/lon_right_M09.txt")
do i=1,nt
  write(88,*)lon_right(i)
enddo
open(88,file="temp/lat_bottom_M09.txt")
do i=1,nt
  write(88,*)lat_bottom(i)
enddo
open(88,file="temp/lat_upper_M09.txt")
do i=1,nt
  write(88,*)lat_top(i)
enddo

end