program main
!Main purpose: Assigns a catchment‐tile index from catchment definition files to each model grid cell for M09 grid.

use MAPL
use river_ncfile_helper
use routing_model_constants, only : nc, nlon, nlat

implicit none

! Variable declarations:
integer, parameter   :: nsub_max_init=9999
integer              :: id, xi, yi, i, flag, subi, nt09, nmax
integer              :: nsub(nc)               ! Array storing the number of sub-areas for each catchment

! Allocatable arrays for sub-catchment information:
integer, allocatable :: xsub0(:,:), ysub0(:,:), xsub(:,:), ysub(:,:)
real, allocatable    :: asub0(:,:), asub(:,:)  ! Aggregated area for each sub-catchment
real, allocatable    :: pfaf_area(:)

! Arrays for grid data and mapping:
real*8, allocatable  :: lon(:), lat(:)         ! Longitude and latitude arrays from NetCDF file
integer, allocatable :: loni(:), lati(:)       ! Mapped integer indices from 1-minute resolution files
integer, allocatable :: catchind(:,:)          ! 2D array of catchment indices for each grid cell


! Define file path for input routing data:
character(len=900)   :: file_path1 !"input/CatchIndex.nc" 
character(len=900)   :: file_path2 !"/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/v12/geometry/EASEv2_M09/rst/EASEv2_M09_3856x1624.rst"  
character(len=900)   :: file_path3 !""

if (command_argument_count() /= 3) then
    print *, "no appropriate files found"
    stop
endif
call get_command_argument(1, file_path1)
call get_command_argument(2, file_path2)
call get_command_argument(3, file_path3)

! Allocate arrays based on the defined dimensions:
allocate(xsub0(nsub_max_init, nc), ysub0(nsub_max_init, nc), asub0(nsub_max_init, nc))
allocate(catchind(nlon, nlat))
allocate(lon(nlon), lat(nlat))
allocate(loni(nlon), lati(nlat))

! Read grid longitude, latitude, catchment index, and cell area data from NetCDF files:
call read_ncfile_double1d(trim(file_path1), "longitude", lon, nlon)
call read_ncfile_double1d(trim(file_path1), "latitude", lat, nlat)
call read_ncfile_int2d(trim(file_path1), "CatchIndex", catchind, nlon, nlat)

! Read mapped grid indices for 1-minute resolution from text files:
open(10, file="temp/lati_1m_M09.txt")
read(10, *) lati
open(11, file="temp/loni_1m_M09.txt")
read(11, *) loni

! Initialize aggregation arrays:
nsub = 0         ! Set number of sub-areas per catchment to zero
xsub = 0         ! Initialize x-coordinate array for sub-catchments to zero
ysub = 0         ! Initialize y-coordinate array for sub-catchments to zero
asub = 0.        ! Initialize aggregated area values to zero

open(77, file=trim(file_path3));read(77, *) nt09
! Loop over each 1m grid cell to accumulate cell areas into sub-catchments:
call EASE_Find_subs(catchind, loni, lati, lat, nt09, 2*nlon, 2*nlat, trim(file_path2), &
	                nsub, xsub0, ysub0, asub0)
nmax = maxval(nsub)
print *,nmax
allocate(xsub(nmax, nc), ysub(nmax, nc), asub(nmax, nc))
xsub=xsub0(1:nmax,:)
ysub=ysub0(1:nmax,:)
asub=asub0(1:nmax,:)
deallocate(xsub0,ysub0,asub0) 
! Open output files to write the aggregated sub-catchment information:
open(50, file="temp/Pfaf_nsub_M09.txt")
open(51, file="temp/Pfaf_xsub_M09.txt")
open(52, file="temp/Pfaf_ysub_M09.txt")
open(53, file="temp/Pfaf_asub_M09.txt")

! For each catchment, write:
!   - Number of sub-areas
!   - X indices of sub-areas (formatted in groups of nmax integers)
!   - Y indices of sub-areas (formatted similarly)
!   - Aggregated area values of sub-areas (formatted as floating-point numbers)
do i = 1, nc
  write(50, *) nsub(i)
  write(51, '(*(1x,i4))') xsub(:, i)
  write(52, '(*(1x,i4))') ysub(:, i)
  write(53, '(*(1x,f10.4))') asub(:, i)
end do

! Print the maximum number of sub-areas found for any catchment and its location:
print *, maxval(nsub)
print *, maxloc(nsub)

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



contains

  ! ----------------------------------------------------------------------------------
  
  subroutine EASE_Find_subs(catchind, loni, lati, lat, nland, nlon30s, nlat30s, rstfile, nsub, xsub0, ysub0, asub0)
    
    integer,           intent(in)  :: catchind(:,:)
    integer,           intent(in)  :: loni(:)          ! lon index of overlying EASE grid cell;                 size = nlon
    integer,           intent(in)  :: lati(:)          ! lat index of overlying EASE grid cell;                 size = nlat
    real(kind=REAL64), intent(in)  :: lat(:)           ! lat of raster grid cell;                               size = nlat
    integer,           intent(in)  :: nland            ! number of land tiles
    integer,           intent(in)  :: nlon30s          ! number of lon coordinate in 30-s grid
    integer,           intent(in)  :: nlat30s          ! number of lat coordinate in 30-s grid
    character(*),      intent(in)  :: rstfile          ! path of *.rst file
    integer,           intent(out) :: nsub(:)          ! # raster grid cells that contribute to Pfaf catchment; size = nPfaf
    integer,           intent(out) :: xsub0(:,:)       ! lon index of overlying EASE grid cell;                 size = nsub_max_init x nPfaf
    integer,           intent(out) :: ysub0(:,:)       ! lat index of overlying EASE grid cell;                 size = nsub_max_init x nPfaf
    real,              intent(out) :: asub0(:,:)       ! area of raster grid cell;                              size = nsub_max_init x nPfaf

    integer :: xi, yi, flag, id, i, j, nlon, nlat
    real(kind=REAL64)             :: cellarea, delta, area30s(2,2)
    real(kind=REAL64),allocatable :: lat30s(:), cellarea30s(:)
    integer,allocatable           :: rst(:,:)

    nlon  = size(loni)
    nlat  = size(lati)
    nsub  = 0 
    xsub0 = 0
    ysub0 = 0
    asub0 = 0.

    allocate(rst(nlon30s,nlat30s),lat30s(nlat30s),cellarea30s(nlat30s)) 
    delta=180./nlat30s
    ! Calculate the lat of 30-s grid.
    do j=1,nlat  
      if(lat(1)<0.d0)then
        lat30s(2*j-1)=lat(j)-delta/2.d0
        lat30s(2*j)=lat(j)+delta/2.d0
      else
        lat30s(2*j-1)=lat(j)+delta/2.d0
        lat30s(2*j)=lat(j)-delta/2.d0
      endif      
    enddo
    delta = 2*MAPL_PI_R8/nlon30s * MAPL_PI_R8/nlat30s
    open(20,file=trim(rstfile),form="unformatted",status="old")
    ! Read rst file and calculate the area of 30-s gridcell.
    do j=1,nlat30s
      read(20) rst(:,j)
      cellarea30s(j) = cos(lat30s(j) * MAPL_DEGREES_TO_RADIANS_R8)*delta * MAPL_RADIUS/1.e3*MAPL_RADIUS/1.e3 !km^2
    enddo
    where (rst<1.or.rst>nland) rst=0 !For no-land tiles in the rst file, set rst to 0.

    delta = 2*MAPL_PI_R8/nlon * MAPL_PI_R8/nlat                    ! pre-compute term for raster grid cell area
    ! Loop through all raster grid cells to aggregate cell areas by catchment and sub-area:
    do yi = 1, nlat
      cellarea = cos(lat(yi) * MAPL_DEGREES_TO_RADIANS_R8)*delta   ! area of raster grid cell for latitude yi [radians^2]
      cellarea = cellarea*MAPL_RADIUS/1.e3*MAPL_RADIUS/1.e3        ! convert to km^2
      do xi = 1, nlon
        if (catchind(xi, yi) >= 1) then
          area30s(:,1)=cellarea30s(2*yi-1) !area30s(2,2) is the area of 30s cells corresponding to the (xi,yi) of 1-min cell.
          area30s(:,2)=cellarea30s(2*yi)
          !set land area to 0 in the area30s(2,2), so the sum(area30s(2,2)) is the no-land area in the (xi,yi) of 1-min cell.
          where(rst(2*xi-1:2*xi,2*yi-1:2*yi)/=0) area30s=0.d0          
          ! The raster grid cell belongs to a catchment
          id = catchind(xi, yi)  ! Get the catchment id for the current cell
          flag = 0               ! Reset flag to indicate whether a matching sub-area is found      
          ! If the catchment already has one or more sub-areas, check for a matching sub-area:
          if (nsub(id) >= 1) then
            do i = 1, nsub(id)
              if (loni(xi) == xsub0(i, id) .and. lati(yi) == ysub0(i, id)) then
                flag = 1
                ! If a match is found, accumulate the cell area into the existing sub-area:
                ! subtract the no-land area sum(area30s) here, so only land area is summed to sub-area
                asub0(i, id) = asub0(i, id) + max(0.,cellarea - sum(area30s)) 
                exit  ! Exit the inner loop since a matching sub-area has been found
              endif
            end do
          endif
          ! If no matching sub-area was found, create a new sub-area:
          if (flag == 0) then
            nsub(           id) = nsub(id) + 1
            xsub0(nsub(id), id) = loni(xi)
            ysub0(nsub(id), id) = lati(yi)
            ! subtract the no-land area sum(area30s) here, so only land area is summed to sub-area
            asub0(nsub(id), id) = max(0.,cellarea - sum(area30s))
          endif
        endif
      end do
    end do

    deallocate(rst,lat30s,cellarea30s)

  end subroutine EASE_Find_subs


end program main
