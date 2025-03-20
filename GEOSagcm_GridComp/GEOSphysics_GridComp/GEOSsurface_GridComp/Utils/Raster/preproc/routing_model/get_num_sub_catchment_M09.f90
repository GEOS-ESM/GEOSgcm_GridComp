program main

use river_read
implicit none

integer,parameter :: nmax=458
integer,parameter :: nc=291284
integer,parameter :: nlon=21600
integer,parameter :: nlat=10800

integer :: id,xi,yi,i,flag,subi
integer :: nsub(nc)
integer,allocatable :: xsub(:,:),ysub(:,:)
real,allocatable :: asub(:,:)

real*8,allocatable :: lon(:),lat(:)
integer,allocatable :: loni(:),lati(:)
integer,allocatable :: catchind(:,:)
real,allocatable :: cellarea(:,:)

allocate(xsub(nmax,nc),ysub(nmax,nc),asub(nmax,nc))
allocate(catchind(nlon,nlat),cellarea(nlon,nlat))
allocate(lon(nlon),lat(nlat))
allocate(loni(nlon),lati(nlat))


call read_ncfile_double1d("input/CatchIndex.nc","lon",lon,nlon)
call read_ncfile_double1d("input/CatchIndex.nc","lat",lat,nlat)
call read_ncfile_int2d("input/CatchIndex.nc","data",catchind,nlon,nlat)
call read_ncfile_real2d("temp/cellarea.nc","data",cellarea,nlon,nlat)
cellarea=cellarea/1.e6


open(10,file="temp/lati_1m_M09.txt")
read(10,*)lati
open(11,file="temp/loni_1m_M09.txt")
read(11,*)loni

nsub=0
xsub=0
ysub=0
asub=0.
do xi=1,nlon
  do yi=1,nlat
    if(catchind(xi,yi)>=1)then

      id=catchind(xi,yi)
      flag=0
      if(nsub(id)>=1)then
        do i=1,nsub(id)
          if(loni(xi)==xsub(i,id).and.lati(yi)==ysub(i,id))then
            flag=1
            asub(i,id)=asub(i,id)+cellarea(xi,yi)            
            exit
          endif
        enddo
      endif
      if(flag==0)then
        nsub(id)=nsub(id)+1
        xsub(nsub(id),id)=loni(xi)
        ysub(nsub(id),id)=lati(yi)
        asub(nsub(id),id)=cellarea(xi,yi)
      endif

    endif
  enddo
enddo

open(50,file="output/Pfaf_nsub_M09.txt")
open(51,file="output/Pfaf_xsub_M09.txt")
open(52,file="output/Pfaf_ysub_M09.txt")
open(53,file="output/Pfaf_asub_M09.txt")
do i=1,nc
  write(50,*)nsub(i)
  write(51,'(458(1x,i4))')xsub(:,i)
  write(52,'(458(1x,i4))')ysub(:,i)
  write(53,'(458(1x,f10.4))')asub(:,i)
enddo

print *,maxval(nsub)
print *,maxloc(nsub)



end
