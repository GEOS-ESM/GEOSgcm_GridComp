program main

use omp_lib
use rwncfile

implicit none

character(len=100) :: lonfile="outlet_sinkx_allcat_TM0072xTM0036_mask_2d.nc"
character(len=100) :: latfile="outlet_sinky_allcat_TM0072xTM0036_mask_2d.nc"
integer, parameter :: nx=43200, ny=21600
integer, allocatable  :: lats(:,:), lons(:,:)
integer i,j

allocate(lats(nx,ny), lons(nx,ny))

call read_ncfile_int2d("outputs/"//trim(latfile),"data",lats,nx,ny)
call read_ncfile_int2d("outputs/"//trim(lonfile),"data",lons,nx,ny)


open(30,file="Outlet_latlon.43200x21600",form="unformatted")


do j = 1, ny
  write (30) lons(:,j)
  write (30) lats(:,j)
end do


end
