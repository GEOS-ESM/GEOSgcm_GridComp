program main
  
  use routing_constant,only : nc,ng,nlon,nlat
  implicit none
  
  character(len=100)  :: var1="outlet_sinky_allcat"
  character(len=100)  :: var2="outlet_sinkx_allcat"
  character(len=100)  :: map="Pfafstetter_Greenland_real"
  
  real*8,allocatable  :: lon(:),lat(:)
  integer,allocatable :: catchind(:,:)
  integer,allocatable :: lons(:,:),lats(:,:)
  integer,allocatable :: data_Pfaf(:)
  
  integer :: i,j,xi,yi,id,nall
  
! Transform the 1d list to the unformatted Fortran binary file "Outlet_latlon.43200x21600" that can be read directly by "mk_runofftbl.F90" of makebcs.
! nc= number of land catchments (excluding Greenland), 291284 in our case.
! ng= number of Greenland catchments, 525 in our case
! nall =number of the total catchments (including Greenland)
! catchind = catchment index:  1-291284 for land catchments, and 291285-291809 for Greenland catchments
! data_Pfaf= lat (or lon) index (on a 30s map) of the final sink point for each catchment,
! lats= a map (with a resolution of 30s) of lat values. For each pixel, the value is the latitude of its final outlet point
! lons= a map (with a resolution of 30s) of lon values. For each pixel, the value is the longitude of its final outlet point

  nall=nc+ng
  allocate(catchind(nlon,nlat),lons(nlon,nlat),lats(nlon,nlat))
  allocate(lon(nlon),lat(nlat))

! Read the raster array of catchment indices
  open(30,file="outputs/"//trim(map),form="unformatted")
  do j = 1,nlat
     read (30) catchind(:,j)
  end do
  
  allocate(data_Pfaf(nall))
! Read in the latitudes associated with each catchment   
  open(77,file="outputs/"//trim(var1)//".txt")
  read(77,*)data_Pfaf
  lats=-999
! For each raster point, find the catchment index and use that, in conjunction with data_Pfaf, to compute a raster array of latitudes  
  do xi=1,nlon
     do yi=1,nlat
        if(catchind(xi,yi)>=1.and.catchind(xi,yi)<=nall)then
           id=catchind(xi,yi)
           lats(xi,yi)=data_Pfaf(id)
        endif
     enddo
  enddo

! Read in the longitudes associated with each catchment   
  open(77,file="outputs/"//trim(var2)//".txt")
  read(77,*)data_Pfaf
  lons=-999
 ! For each raster point, find the catchment index and use that, in conjunction with data_Pfaf, to compute a raster array of longitudes  
  do xi=1,nlon
     do yi=1,nlat
        if(catchind(xi,yi)>=1.and.catchind(xi,yi)<=nall)then
           id=catchind(xi,yi)
           lons(xi,yi)=data_Pfaf(id)
        endif
     enddo
  enddo
  
  open(30,file="Outlet_latlon.43200x21600",form="unformatted")
  do j = 1, nlat
     write (30) lons(:,j)
     write (30) lats(:,j)
  end do

end program main
