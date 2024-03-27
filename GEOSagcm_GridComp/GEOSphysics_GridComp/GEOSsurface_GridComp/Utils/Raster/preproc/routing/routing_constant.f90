module routing_constant
  
  ! hardwired constants for GEOS river routing scheme

  implicit none
  public
  
  integer,parameter :: nlon       =  43200  ! number of lat of world grid with 30 sec resolution
  integer,parameter :: nlat       =  21600  ! number of lon of world grid with 30 sec resolution
  integer,parameter :: nlon1m     =  21600  ! number of lat of world grid with  1 min resolution
  integer,parameter :: nlat1m     =  10800  ! number of lon of world grid with  1 min resolution
  integer,parameter :: nlon_G     =   8400  ! number of lat of Greenland grid (30 sec resolution)
  integer,parameter :: nlat_G     =   4800  ! number of lon of Greenland grid (30 sec resolution)
  integer,parameter :: loni_min   =  12001  !                          ind of lon start of Greenland grid in world grid (30 sec resolution)
  integer,parameter :: loni_max   =  20400  ! = loni_min + nlon_G - 1, ind of lon end   of Greenland grid in world grid (30 sec resolution)
  integer,parameter :: lati_min   =  16801  !                          ind of lat start of Greenland grid in world grid (30 sec resolution)
  integer,parameter :: lati_max   =  21600  ! = lati_min + nlat_G - 1, ind of lat end   of Greenland grid in world grid (30 sec resolution)
  integer,parameter :: nc         = 291284  ! number of catchments in land
  integer,parameter :: ng         =    525  ! number of catchments in Greenland
  integer,parameter :: nl         =  22116  ! number of outlets to ocean in land (not including Greenland)

end module routing_constant
