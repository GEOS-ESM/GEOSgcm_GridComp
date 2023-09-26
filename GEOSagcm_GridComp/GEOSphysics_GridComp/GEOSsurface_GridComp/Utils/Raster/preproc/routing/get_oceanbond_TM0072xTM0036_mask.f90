program main

use omp_lib
use rwncfile
implicit none

integer,parameter :: nlon=43200
integer,parameter :: nlat=21600
real*8,allocatable :: lon(:),lat(:)

integer,allocatable :: catchind(:,:)
integer,allocatable :: boundary(:,:)

integer :: xi,yi,id
integer :: xp1,xm1,yp1,ym1

allocate(catchind(nlon,nlat),boundary(nlon,nlat),lon(nlon),lat(nlat))
call read_ncfile_double1d("outputs/TM0072xTM0036_mask.nc","lon",lon,nlon)
call read_ncfile_double1d("outputs/TM0072xTM0036_mask.nc","lat",lat,nlat)
call read_ncfile_int2d("outputs/TM0072xTM0036_mask.nc","data",catchind,nlon,nlat)

boundary=catchind
boundary=-9999

!$OMP PARALLEL default(shared) private(xi,yi,id) 
!$OMP DO
do xi=2,nlon-1
  !if(mod(xi,100)==0)then
  !  print *,xi
  !endif
  do yi=2,nlat-1
    id=catchind(xi,yi)
    if(id==1)then
      boundary(xi,yi)=0 
      if(catchind(xi+1,yi)==1.and.&
         catchind(xi+1,yi-1)==1.and.&
         catchind(xi  ,yi-1)==1.and.&
         catchind(xi-1,yi-1)==1.and.&
         catchind(xi-1,yi)==1.and.&
         catchind(xi-1,yi+1)==1.and.&
         catchind(xi  ,yi+1)==1.and.&
         catchind(xi+1,yi+1)==1)then
         boundary(xi,yi)=-9999
      endif
    endif
  enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

call create_ncfile_int2d("outputs/TM0072xTM0036_mask_oceanboundary.nc","data",boundary,lon,lat,nlon,nlat)

end
