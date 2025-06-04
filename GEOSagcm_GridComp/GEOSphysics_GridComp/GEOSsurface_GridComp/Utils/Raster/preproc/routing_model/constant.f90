module constant
!module for constants used in the river routing pre-processing package  

  implicit none
  public

  ! Define constant parameters 
  integer,parameter :: nmax09=458         ! Maximum number of sub-catchments per catchment for M09
  integer,parameter :: nmax36=150         ! Maximum number of sub-catchments per catchment for M36
  integer,parameter :: nc=291284          ! Total number of catchments
  integer,parameter :: nlon=21600         ! Number of longitude grid points in the original grid
  integer,parameter :: nlat=10800         ! Number of latitude grid points in the original grid
  integer,parameter :: nlat09=1624, nlon09=3856  ! Dimensions for the M09 grid
  integer,parameter :: nlat36=406, nlon36=964  ! Dimensions for the aggregated M36 grid  
  integer,parameter :: nt_global09=1684725  ! Total number of global tiles for area mapping for M09
  integer,parameter :: nt_global36=112573   ! Total number of global tiles for area mapping for M36
  ! Define grid dimensions for 15-second resolution data (HydroSHEDS high-res grid)
  integer,parameter :: nlonh = 86400
  integer,parameter :: nlath = 33600  

  integer,parameter :: nl_USGS = 3352492           ! Total number of USGS data records
  integer,parameter :: nt09=1684725, nt36=112573    !Total number of catchment gridcell in M09 and M36
  integer,parameter :: nupmax = 34    ! Maximum number of upstream catchments to record

  !river curve parameters
  real,parameter    :: cur_avg = 1.4        
  real,parameter    :: cur_min = 0.5        
  real,parameter    :: cur_max = 5.0 

  integer,parameter :: nga = 9067  ! Number of GAGE-II records

  !lake input parameters:
  integer,parameter :: no = 1459201  ! Total number of outlet records in the outlet files
  integer,parameter :: nvl = 3409    ! Number of lakes that pass the filtering criteria (area >= 50)
  integer,parameter :: nvo = 3917    ! Number of outlet records after matching with lakes
  integer,parameter :: nl_lake = 1426967  ! Total number of lake records in the input files  


end module constant