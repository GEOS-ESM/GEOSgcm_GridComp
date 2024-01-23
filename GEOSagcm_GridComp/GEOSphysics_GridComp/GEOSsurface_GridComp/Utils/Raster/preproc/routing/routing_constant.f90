module routing_constant
  
  implicit none
  public
  
  integer,parameter :: nlon       =  43200  ! number of lat of the world grid with 30 sec resolution
  integer,parameter :: nlat       =  21600  ! number of lon of the world grid with 30 sec resolution
  integer,parameter :: nlon1m     =  21600  ! number of lat of the world grid with  1 min resolution
  integer,parameter :: nlat1m     =  10800  ! number of lon of the world grid with  1 min resolution
  integer,parameter :: nlon_G     =   8400  ! number of lat of the Greenland grid (30 sec resolution)
  integer,parameter :: nlat_G     =   4800  ! number of lon of the Greenland grid (30 sec resolution)
  integer,parameter :: loni_min   =  12001  ! index of the lon start of the Greenland grid in the world grid (30 sec resolution)
  integer,parameter :: loni_max   =  20400  ! index of the lon end of the Greenland in the world grid (30 sec resolution)
  integer,parameter :: lati_min   =  16801  ! index of the lat start of the Greenland grid in the world grid (30 sec resolution)
  integer,parameter :: lati_max   =  21600  ! index of the lat end of the Greenland in the world grid (30 sec resolution)
  integer,parameter :: nc         = 291284  ! number of catchments in land
  integer,parameter :: ng         =    525  ! number of catchments in Greenland
  integer,parameter :: nl         =  22116  ! number of outlets to ocean in land (not including Greenland)

end module routing_constant
