program main

use river_read        ! Use custom module for reading NetCDF files
implicit none

integer,parameter :: nlat=1624,nlon=3856
integer,parameter :: nmax=458           ! Maximum number of sub-areas per catchment
integer,parameter :: nc=291284          ! Total number of catchments

integer,allocatable :: map_tile(:,:),subx(:,:),suby(:,:),subi(:,:)

integer :: i,x,y,j,it

allocate(map_tile(nlon,nlat))
call read_ncfile_int2d("temp/map_tile_M09.nc", "data", map_tile, nlon, nlat)
allocate(subx(nmax,nc),suby(nmax,nc),subi(nmax,nc))
open(77,file="output/Pfaf_xsub_M09.txt"); read(77,*)subx
open(77,file="output/Pfaf_ysub_M09.txt"); read(77,*)suby
subi=0
do i=1,nc
  do j=1,nmax
    x=subx(j,i)
    y=suby(j,i)
    if(x/=0)then
      if(y==0)stop
      subi(j,i)=map_tile(x,y)
    endif
  enddo
enddo

open(88,file="output/Pfaf_isub_M09.txt")
do i=1,nc
	write(88,'(150(i8))') subi(:,i)
enddo




end program main
