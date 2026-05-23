program main
!Main purpose: Assigns a catchment‐tile index from catchment definition files to each model grid cell for EASE grid.

use MAPL
use river_ncfile_helper
use routing_model_constants, only : nc, nlon, nlat, np=>np_tot, nlatE

implicit none

! Variable declarations:
integer, parameter   :: nsub_max_init=9999
integer              :: id, xi, yi, i, j, flag, subi, nt09, nmax, cid, k

! Allocatable arrays for sub-catchment information:
integer, allocatable :: xsub0(:,:), ysub0(:,:), xsub(:,:), ysub(:,:), nsub(:)
real, allocatable    :: asub0(:,:), asub(:,:)  ! Aggregated area for each sub-catchment
real, allocatable    :: pfaf_area(:)
integer, allocatable :: type_tile(:),xi_tile(:),yi_tile(:),catid_tile(:)
real,    allocatable :: area_tile(:)

integer,           allocatable :: iTable(:,:)
real(kind=REAL64), allocatable :: rTable(:,:)  

! Arrays for grid data and mapping:
real*8, allocatable  :: lon(:), lat(:)         ! Longitude and latitude arrays from NetCDF file
integer, allocatable :: loni(:), lati(:)       ! Mapped integer indices from 1-minute resolution files
integer, allocatable :: catchind(:,:)          ! 2D array of catchment indices for each grid cell


! Define file path for input routing data:
character(len=900)   :: file_path1 

if (command_argument_count() /= 1) then
    print *, "no appropriate files found"
    stop
endif
call get_command_argument(1, file_path1)

! Allocate arrays based on the defined dimensions:
allocate(xsub(nsub_max_init, nc), ysub(nsub_max_init, nc), asub(nsub_max_init, nc), nsub(nc))

! Initialize aggregation arrays:
nsub  = 0         ! Set number of sub-areas per catchment to zero
xsub  = 0         ! Initialize x-coordinate array for sub-catchments to zero
ysub  = 0         ! Initialize y-coordinate array for sub-catchments to zero
asub  = 0.        ! Initialize aggregated area values to zero

allocate(type_tile(np),xi_tile(np),yi_tile(np),catid_tile(np),area_tile(np))

call read_ncfile_int1d(trim(file_path1),"typ",type_tile,np)
call read_ncfile_int1d(trim(file_path1),"i_indg",xi_tile,np)
xi_tile = xi_tile + 1
call read_ncfile_int1d(trim(file_path1),"j_indg",yi_tile,np)
yi_tile = nlatE - yi_tile
call read_ncfile_int1d(trim(file_path1),"pfaf_index",catid_tile,np) 
call read_ncfile_real1d(trim(file_path1),"area",area_tile,np) 
area_tile = area_tile * MAPL_RADIUS**2 / 1.e6 !km^2

do i=1,np
  if(type_tile(i)==100)then
    cid=catid_tile(i)
    nsub(cid)=nsub(cid)+1
    xsub(nsub(cid),cid)=xi_tile(i)
    ysub(nsub(cid),cid)=yi_tile(i)
    asub(nsub(cid),cid)=area_tile(i)
  endif
enddo
deallocate(type_tile,xi_tile,yi_tile,catid_tile,area_tile)
nmax = maxval(nsub)
print *,"total sub-catchment number: ",sum(nsub)
! Open output files to write the aggregated sub-catchment information:
! For each catchment, write:
!   - Number of sub-areas
!   - X indices of sub-areas 
!   - Y indices of sub-areas 
!   - Aggregated area values of sub-areas
!   - Catchment Index of sub-areas 
open(50, file="temp/Pfaf_nsub_EASE.txt")
do i=1,nc
  write(50, *) nsub(i)  
end do
open(51, file="temp/Pfaf_xsub_EASE.txt")
open(52, file="temp/Pfaf_ysub_EASE.txt")
open(53, file="temp/Pfaf_asub_EASE.txt")
open(54, file="temp/Pfaf_csub_EASE.txt")
do i=1,nc
  if(nsub(i)>0)then
    do j=1,nsub(i)
      write(51, *)xsub(j, i) 
      write(52, *)ysub(j, i)
      write(53, *)asub(j, i)
      write(54, *)i     
    enddo
  endif
end do

allocate(pfaf_area(nc))
pfaf_area = 0.
do i = 1, nc
  pfaf_area(i) = sum(asub(:,i))
end do

! Write catchment areas to an output file:
open(88, file="temp/Pfaf_area.txt")
do i = 1, nc
  write(88, *) pfaf_area(i)
end do


end program main
