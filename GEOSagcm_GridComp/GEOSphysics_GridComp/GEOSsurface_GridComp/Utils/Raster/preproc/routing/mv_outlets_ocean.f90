program main

use omp_lib
use rwncfile
implicit none

integer,parameter :: nsh=788292
integer,parameter :: ns=22612
integer,parameter :: nlon=43200
integer,parameter :: nlat=21600

integer,allocatable :: mask(:,:)

real*8 :: lons(ns),lats(ns)
integer :: lonsi(ns),latsi(ns)
real*8 :: lons_adj(ns),lats_adj(ns)
real*8,allocatable :: lonsh(:),latsh(:)
integer :: catid(ns),flag(ns)

real :: dist(ns)

integer :: i,j
real :: dy,dy2,dx,dx2,dxA,dxB,dist_temp 

allocate(lonsh(nsh),latsh(nsh))

open(77,file="outputs/outlet_sinklon.txt")
read(77,*)lons
open(77,file="outputs/outlet_sinklat.txt")
read(77,*)lats
open(77,file="outputs/outlet_sinkx.txt")
read(77,*)lonsi
open(77,file="outputs/outlet_sinky.txt")
read(77,*)latsi

open(77,file="outputs/lon_oceanbond_list_TM0072xTM0036_mask.txt")
read(77,*)lonsh
open(77,file="outputs/lat_oceanbond_list_TM0072xTM0036_mask.txt")
read(77,*)latsh

allocate(mask(nlon,nlat))
call read_ncfile_int2d("outputs/TM0072xTM0036_mask.nc","data",mask,nlon,nlat)

!$OMP PARALLEL default(shared) private(i,j,dy,dy2,dx,dx2,dxA,dxB,dist_temp) 
!$OMP DO
do i=1,ns
  !if(mod(i,100)==0) print *,i
  IF(mask(lonsi(i),latsi(i))==0)THEN
  dist(i)=1.e12
  do j=1,nsh
    dy=abs(lats(i)-latsh(j))
    dy2=dy*dy
    dxA=abs(lons(i)-lonsh(j)) 
    dxB=360.-dxA
    dx=min(dxA,dxB)
    dx2=dx*dx
    dist_temp=sqrt(dx2+dy2)
    if(dist_temp<dist(i))then
      dist(i)=dist_temp
      lons_adj(i)=lonsh(j)
      lats_adj(i)=latsh(j)
    endif
  enddo
  ELSE
    lons_adj(i)=lons(i)
    lats_adj(i)=lats(i)
    dist(i)=0.
  ENDIF
enddo 
!$OMP END DO
!$OMP END PARALLEL

open(88,file="outputs/outlet_sinklat_TM0072xTM0036_mask.txt")
do i=1,ns
  write(88,*)lats_adj(i)
enddo
open(88,file="outputs/outlet_sinklon_TM0072xTM0036_mask.txt")
do i=1,ns
  write(88,*)lons_adj(i)
enddo


end
