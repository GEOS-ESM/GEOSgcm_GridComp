program main

use omp_lib
use rwncfile
implicit none

integer,parameter :: nsh=788292 !304483  !1877262
integer,parameter :: nlonh=43200
integer,parameter :: nlath=21600
real*8,allocatable :: lonh(:),lath(:)
integer,allocatable :: mskh(:,:)
real*8,allocatable :: lonsh(:),latsh(:)

integer i,xi,yi,k


allocate(mskh(nlonh,nlath))
allocate(lonh(nlonh),lath(nlath))
call read_ncfile_double1d("outputs/TM0072xTM0036_mask_oceanboundary.nc","lon",lonh,nlonh)
call read_ncfile_double1d("outputs/TM0072xTM0036_mask_oceanboundary.nc","lat",lath,nlath)
call read_ncfile_int2d("outputs/TM0072xTM0036_mask_oceanboundary.nc","data",mskh,nlonh,nlath)

allocate(lonsh(nsh),latsh(nsh))

k=0
!!$OMP PARALLEL default(shared) shared(k) private(xi,yi) 
!!$OMP DO
do xi=1,nlonh
  do yi=1,nlath
    if(mskh(xi,yi)==0)then
      k=k+1
      lonsh(k)=lonh(xi)
      latsh(k)=lath(yi)
    endif
  enddo
enddo
!!$OMP END DO
!!$OMP END PARALLEL
!print *,k


open(88,file="outputs/lon_oceanbond_list_TM0072xTM0036_mask.txt")
do i=1,nsh
  write(88,*)lonsh(i)
enddo
open(88,file="outputs/lat_oceanbond_list_TM0072xTM0036_mask.txt")
do i=1,nsh
  write(88,*)latsh(i)
enddo





end
