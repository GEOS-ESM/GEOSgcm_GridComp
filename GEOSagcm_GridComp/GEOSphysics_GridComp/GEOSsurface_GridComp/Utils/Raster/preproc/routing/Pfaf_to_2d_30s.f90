program main

use omp_lib
use rwncfile
implicit none

character(len=100) :: var1="outlet_sinky_allcat_TM0072xTM0036_mask"
character(len=100) :: var2="outlet_sinkx_allcat_TM0072xTM0036_mask"
character(len=100) :: map="TM0072xTM0036-Pfafstetter_Greenland_real.nc"
integer,parameter :: nc=291809
integer,parameter :: nlon=43200
integer,parameter :: nlat=21600

real*8,allocatable :: lon(:),lat(:)
integer,allocatable :: catchind(:,:)
integer,allocatable :: data2d(:,:)
integer,allocatable :: data_Pfaf(:)

integer :: xi,yi,id


allocate(catchind(nlon,nlat),data2d(nlon,nlat))
allocate(lon(nlon),lat(nlat))
call read_ncfile_double1d("inputs/"//trim(map),"lon",lon,nlon)
call read_ncfile_double1d("inputs/"//trim(map),"lat",lat,nlat)
call read_ncfile_int2d("inputs/"//trim(map),"data",catchind,nlon,nlat)

allocate(data_Pfaf(nc))

open(77,file="outputs/"//trim(var1)//".txt")
read(77,*)data_Pfaf
data2d=-999
do xi=1,nlon
  do yi=1,nlat
    if(catchind(xi,yi)>=1.and.catchind(xi,yi)<=nc)then
      id=catchind(xi,yi)
      data2d(xi,yi)=data_Pfaf(id)
    endif 
  enddo
enddo
call create_ncfile_int2d_fill("outputs/"//trim(var1)//"_2d.nc","data",data2d,lon,lat,nlon,nlat,-999.)

open(77,file="outputs/"//trim(var2)//".txt")
read(77,*)data_Pfaf
data2d=-999
do xi=1,nlon
  do yi=1,nlat
    if(catchind(xi,yi)>=1.and.catchind(xi,yi)<=nc)then
      id=catchind(xi,yi)
      data2d(xi,yi)=data_Pfaf(id)
    endif 
  enddo
enddo
call create_ncfile_int2d_fill("outputs/"//trim(var2)//"_2d.nc","data",data2d,lon,lat,nlon,nlat,-999.)




end
