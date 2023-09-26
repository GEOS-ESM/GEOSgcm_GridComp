program main

use omp_lib
use rwncfile
implicit none

integer,parameter :: nlon=43200
integer,parameter :: nlat=21600
integer,parameter :: nt=2592

real*8,allocatable,dimension(:) :: lon,lat
integer,allocatable,dimension(:,:) :: landocean,mask
integer,allocatable,dimension(:) :: mask1d

integer :: i,j,xi,yi,tid

allocate(landocean(nlon,nlat))
allocate(lon(nlon),lat(nlat))
call read_ncfile_double1d("inputs/TM0072xTM0036.nc","lon",lon,nlon)
call read_ncfile_double1d("inputs/TM0072xTM0036.nc","lat",lat,nlat)
call read_ncfile_int2d("inputs/TM0072xTM0036.nc","data",landocean,nlon,nlat)


allocate(mask(nlon,nlat),mask1d(nt))
open(77,file="outputs/mask_MAPL_1d_TM0072xTM0036.txt")
read(77,*)mask1d
do i=1,nlon
  do j=1,nlat
  	tid=landocean(i,j)
    mask(i,j)=mask1d(tid)
  enddo
enddo
call create_ncfile_int2d("outputs/mask_MAPL_TM0072xTM0036.nc","data",mask,lon,lat,nlon,nlat)


end
