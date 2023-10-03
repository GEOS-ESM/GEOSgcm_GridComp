program main

use omp_lib
use rwncfile
implicit none

integer,parameter :: nc=291809
integer,parameter :: nlon=43200
integer,parameter :: nlat=21600
integer,parameter :: nlon_G=8400
integer,parameter :: nlat_G=4800
integer,parameter :: loni_min=12001
integer,parameter :: loni_max=20400
integer,parameter :: lati_min=16801
integer,parameter :: lati_max=21600

integer,parameter :: id_glac=286926
integer,parameter :: id_lake=286925
integer,parameter :: id_landend=284954

real*8,allocatable,dimension(:) :: lon,lat,lon_G,lat_G
integer,allocatable,dimension(:,:) :: landocean,Greenland
integer,allocatable,dimension(:) :: Pfaf_real, countc

integer :: i,j

allocate(landocean(nlon,nlat))
allocate(lon(nlon),lat(nlat))
call read_ncfile_double1d("inputs/TM0072xTM0036-Pfafstetter.nc","lon",lon,nlon)
call read_ncfile_double1d("inputs/TM0072xTM0036-Pfafstetter.nc","lat",lat,nlat)
call read_ncfile_int2d("inputs/TM0072xTM0036-Pfafstetter.nc","data",landocean,nlon,nlat)

allocate(Greenland(nlon_G,nlat_G))
allocate(lon_G(nlon_G),lat_G(nlat_G))
call read_ncfile_double1d("inputs/GreenlandID_30s.nc","lon",lon_G,nlon_G)
call read_ncfile_double1d("inputs/GreenlandID_30s.nc","lat",lat_G,nlat_G)
call read_ncfile_int2d("inputs/GreenlandID_30s.nc","data",Greenland,nlon_G,nlat_G)

where(Greenland/=-9999.and.(landocean(loni_min:loni_max,lati_min:lati_max)<=id_landend.or.&
	landocean(loni_min:loni_max,lati_min:lati_max)==id_glac ))&
 landocean(loni_min:loni_max,lati_min:lati_max)=Greenland


where(landocean>id_landend.and.landocean<id_lake) landocean=-9999
where(landocean==id_lake.or.landocean==id_glac) landocean=0

allocate(Pfaf_real(id_landend))
open(77,file="inputs/TM0072xTM0036-Pfaf_real.txt")
read(77,*)Pfaf_real

do i=1,nlon
  do j=1,nlat
    if(landocean(i,j)<=id_landend.and.landocean(i,j)>=1)then
      landocean(i,j)=Pfaf_real(landocean(i,j))
    else if(landocean(i,j)>=700000000)then
      landocean(i,j)=landocean(i,j)-700000000+291284
    endif
  enddo
enddo



call create_ncfile_int2d("outputs/TM0072xTM0036-Pfafstetter_Greenland_real.nc","data",landocean,lon,lat,nlon,nlat)


end