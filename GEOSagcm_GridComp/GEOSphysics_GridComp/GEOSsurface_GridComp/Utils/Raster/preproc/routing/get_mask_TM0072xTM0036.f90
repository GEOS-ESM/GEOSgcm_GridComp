program main

use omp_lib
use rwncfile

implicit none

integer,parameter :: nlon=43200
integer,parameter :: nlat=21600
real*8,allocatable :: lon(:),lat(:)

integer,allocatable,dimension(:,:) :: mask_mapl,mask_rst,mask

allocate(mask_mapl(nlon,nlat),mask_rst(nlon,nlat),lon(nlon),lat(nlat),mask(nlon,nlat))

call read_ncfile_double1d("inputs/TM0072xTM0036-Pfafstetter_Greenland_real.nc","lon",lon,nlon)
call read_ncfile_double1d("inputs/TM0072xTM0036-Pfafstetter_Greenland_real.nc","lat",lat,nlat)
call read_ncfile_int2d("inputs/TM0072xTM0036-Pfafstetter_Greenland_real.nc","data",mask_rst,nlon,nlat)
call read_ncfile_int2d("outputs/mask_MAPL_TM0072xTM0036.nc","data",mask_mapl,nlon,nlat)


mask=0
where(mask_rst==-9999..and.mask_mapl==1)mask=1

call create_ncfile_int2d("outputs/TM0072xTM0036_mask.nc","data",mask,lon,lat,nlon,nlat)



end