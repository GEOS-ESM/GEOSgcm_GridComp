program main

use river_read
implicit none

integer,parameter :: nmax=150
integer,parameter :: nc=291284
integer,parameter :: nlon=21600
integer,parameter :: nlat=10800
integer,parameter :: nlat36=406,nlon36=964
integer,parameter :: nt_global=112573

integer :: id,xi,yi,i,j,flag,subi,x_m36,y_m36,it
integer :: nsub(nc)
integer,allocatable :: xsub(:,:),ysub(:,:),subi_global(:,:)
real,allocatable :: asub(:,:)

real*8,allocatable :: lon(:),lat(:)
integer,allocatable :: loni(:),lati(:)
integer,allocatable :: catchind(:,:),map_tile(:,:)
real,allocatable :: cellarea(:,:),area_m36(:,:),area_tile(:)
real*8,allocatable :: lat36(:),lon36(:)


!allocate(subi_global(nmax,nc))
!open(77,file="Pfaf_isub_M36.txt",status="old",action="read"); read(77,*)subi_global; close(77)
!open(90,file="subi.txt",action="write")
!do i=1,nc
!  write(90,'(150(i7))')subi_global(:,i)
!end do
!print *,"successful"
!stop

allocate(xsub(nmax,nc),ysub(nmax,nc),asub(nmax,nc))
allocate(catchind(nlon,nlat),cellarea(nlon,nlat))
allocate(lon(nlon),lat(nlat))
allocate(loni(nlon),lati(nlat))


call read_ncfile_double1d("input/CatchIndex.nc","lon",lon,nlon)
call read_ncfile_double1d("input/CatchIndex.nc","lat",lat,nlat)
call read_ncfile_int2d("input/CatchIndex.nc","data",catchind,nlon,nlat)
call read_ncfile_real2d("temp/cellarea.nc","data",cellarea,nlon,nlat)
cellarea=cellarea/1.e6


open(10,file="temp/lati_1m_M36.txt")
read(10,*)lati
open(11,file="temp/loni_1m_M36.txt")
read(11,*)loni


allocate(area_m36(nlon36,nlat36))
area_m36=0.
do xi=1,nlon
  do yi=1,nlat
    if(catchind(xi,yi)>=1)then
      x_m36=loni(xi)
      y_m36=lati(yi)
      area_m36(x_m36,y_m36)=area_m36(x_m36,y_m36)+cellarea(xi,yi)
    endif
  enddo
enddo

allocate(map_tile(nlon36,nlat36))
call read_ncfile_int2d("temp/map_tile_M36.nc","data",map_tile,nlon36,nlat36)
allocate(area_tile(nt_global))
area_tile=-9999.
do i=1,nlon36
  do j=1,nlat36
    it=map_tile(i,j)
    if(it>0)then
      area_tile(it)=area_m36(i,j)
    endif
  enddo
enddo

open(88,file="output/area_M36_1d.txt")
do i=1,nt_global
  write(88,*)area_tile(i)
enddo


end
