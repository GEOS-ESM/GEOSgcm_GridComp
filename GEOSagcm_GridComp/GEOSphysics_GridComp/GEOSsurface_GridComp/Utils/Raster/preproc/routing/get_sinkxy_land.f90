program main

use rwncfile
implicit none

integer,parameter :: ns=22612
integer,parameter :: nlat=21600
integer,parameter :: nlon=43200

real*8,allocatable,dimension(:) :: lats,lons,lat30s,lon30s,lat_dis,lon_dis
integer,allocatable,dimension(:) :: lati,loni

integer :: i,temp(1)

allocate(lats(ns),lons(ns),lati(ns),loni(ns))
allocate(lat30s(nlat),lon30s(nlon),lat_dis(nlat),lon_dis(nlon))
open(77,file="outputs/outlet_sinklat.txt")
read(77,*)lats
open(77,file="outputs/outlet_sinklon.txt")
read(77,*)lons
open(77,file="inputs/lat_30s.txt")
read(77,*)lat30s
open(77,file="inputs/lon_30s.txt")
read(77,*)lon30s

do i=1,ns
  lat_dis=abs(lat30s-lats(i))
  temp=minloc(lat_dis)
  lati(i)=temp(1)
enddo
do i=1,ns
  lon_dis=abs(lon30s-lons(i))
  temp=minloc(lon_dis)
  loni(i)=temp(1)
enddo

open(88,file="outputs/outlet_sinky.txt")
do i=1,ns
  write(88,*)lati(i)
enddo
open(88,file="outputs/outlet_sinkx.txt")
do i=1,ns
  write(88,*)loni(i)
enddo

end