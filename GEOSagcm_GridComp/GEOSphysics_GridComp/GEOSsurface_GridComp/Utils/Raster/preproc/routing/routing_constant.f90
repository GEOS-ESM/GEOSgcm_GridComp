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
  integer,parameter :: id_glac    = 290191  ! index of glacier tiles in the Pfafstetter.rst
  integer,parameter :: id_lake    = 290190  ! index of lake tiles in the Pfafstetter.rst
  integer,parameter :: id_landend = 290188  ! index of the last land tile in the Pfafstetter.rst
  integer,parameter :: nc         = 291284  ! number of catchments in land
  integer,parameter :: ns         =  22612  ! number of outlets to ocean
  integer,parameter :: ng         =    525  ! number of catchments in Greenland
  integer,parameter :: nl         =  22087  ! number of outlets to ocean in land (not including Greenland)
  integer,parameter :: nall       = 291809  ! total number of catchments in land and Greenland

end module routing_constant
