program main

use omp_lib
use rwncfile
implicit none

integer,parameter :: nc=291284
integer,parameter :: nl=22087
integer,parameter :: ng=525
integer,parameter :: nlon=21600
integer,parameter :: nlat=10800

real*8,allocatable :: lon(:),lat(:),long(:),latg(:),lons(:),lats(:)
integer,allocatable :: catchind(:,:)
real,allocatable :: acah(:,:)
integer,allocatable :: down(:),sx(:),sy(:),msk(:)
real,allocatable :: acas(:)

integer :: id,xi,yi,i,k,xis,yis,ntot

ntot=nl+ng
allocate(catchind(nlon,nlat),acah(nlon,nlat))
allocate(lon(nlon),lat(nlat))
allocate(sx(nc),sy(nc),acas(nc),down(nc),msk(nc))
allocate(long(ng),latg(ng),lons(ntot),lats(ntot))

call read_ncfile_double1d("inputs/CatchIndex.nc","lon",lon,nlon)
call read_ncfile_double1d("inputs/CatchIndex.nc","lat",lat,nlat)
call read_ncfile_int2d("inputs/CatchIndex.nc","data",catchind,nlon,nlat)
call read_ncfile_real2d("inputs/HydroSHEDS_drainage_area.nc","data",acah,nlon,nlat)


open(77,file="inputs/downstream_1D_new_noadj.txt")
read(77,*)down
open(77,file="inputs/Pfaf_msk.txt")
read(77,*)msk

acas=-9999.
sx=0
sy=0
do xi=1,nlon
  do yi=1,nlat
    if(catchind(xi,yi)>=1)then
      id=catchind(xi,yi)
      if(down(id)==-1.and.acah(xi,yi)>=acas(id))then
        acas(id)=acah(xi,yi)
        sx(id)=xi
        sy(id)=yi
      endif
    endif
  enddo
enddo

where(down/=-1)sx=-1
where(down/=-1)sy=-1
k=0
do i=1,nc
  if(msk(i)==2)then
    k=k+1
    lons(k)=lon(sx(i))
    lats(k)=lat(sy(i))
  endif
enddo
!print *,k

open(77,file="inputs/Greenland_outlets_lat.txt")
read(77,*)latg
open(77,file="inputs/Greenland_outlets_lon.txt")
read(77,*)long

lons(k+1:ntot)=long
lats(k+1:ntot)=latg


open(88,file="outputs/outlet_sinklat.txt")
do i=1,ntot
    write(88,*)lats(i)
enddo
open(88,file="outputs/outlet_sinklon.txt")
do i=1,ntot
    write(88,*)lons(i)
enddo




end
