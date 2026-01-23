program main
  
  use routing_constant,only : nl,ng,nlon,nlat
  implicit none
  
  real*8,allocatable,dimension(:)  :: lats,lons,lat30s,lon30s,lat_dis,lon_dis
  integer,allocatable,dimension(:) :: lati,loni
  
  integer :: i,temp(1),ns
  real*8  :: dlat,dlon

! Convert outlet locations in degree lat/lon to indices on the 30 arc-sec raster grid.
! lats = lat degree of the outlets
! lons = lon degree of the outlets
! lati = lat index of the outlets on the 30s map
! loni =lon index of the outlets on the 30s map 
! lat30s = lat coordinate  of the 30s map
! lon30s =lon coordinate of the 30s map
! lat_dis = distance of lat between the center of each 30s pixel and the outlet point.
! lon_dis = distance of lon between the center of each 30s pixel and the outlet point.
! nlat = number of latitude indices. For 30s map, nlat = 21600
! nlon = number of longitude indices.  For 30s map, nlon = 43200

  ns=nl+ng
  allocate(lats(ns),lons(ns),lati(ns),loni(ns))
  allocate(lat30s(nlat),lon30s(nlon),lat_dis(nlat),lon_dis(nlon))
  open(77,file="outputs/outlet_sinklat.txt")
  read(77,*)lats
  open(77,file="outputs/outlet_sinklon.txt")
  read(77,*)lons
  
  dlat=180.D0/nlat
  dlon=360.D0/nlon
  lat30s(1)=-90.D0+dlat/2.D0
  lon30s(1)=-180.D0+dlon/2.D0
  do i=2,nlat
    lat30s(i)=lat30s(i-1)+dlat
  enddo
  do i=2,nlon
    lon30s(i)=lon30s(i-1)+dlon
  enddo  

! For each sink catchment, find the 30s-latitude and 30s-longitude closest to its outlet point.  
  do i=1,ns
     lat_dis=abs(lat30s-lats(i))
     temp=minloc(lat_dis)
     lati(i)=temp(1)
  enddo
  do i=1,ns
     lon_dis=abs(lon30s-lons(i))
     temp=minloc(lon_dis)
     loni(i)=temp(1)
  enddo
  
  open(88,file="outputs/outlet_sinky.txt")
  do i=1,ns
     write(88,*)lati(i)
  enddo
  open(88,file="outputs/outlet_sinkx.txt")
  do i=1,ns
     write(88,*)loni(i)
  enddo
  
end program main
