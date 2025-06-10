program main
!Main purpose: Calculates the K parameter used in the river routing model.

  use k_module   ! Import custom module "k_module" which contains necessary subroutines and functions
  use constant, only: nl=>nl_USGS,nlat,nlon,nc
  implicit none

  ! Declare variables and allocatable arrays
  integer, allocatable            :: lati(:), loni(:)       ! Arrays to store grid indices (latitude and longitude) for stations
  real, allocatable               :: data(:, :)  ! 2D data array to store USGS data
  integer, allocatable            :: catid_full(:), catid(:)  ! Arrays to store full and filtered catchment ids for stations
  real, allocatable, dimension(:) :: vel, dis  ! Arrays to store velocity and distance data from the USGS dataset
  integer, allocatable            :: nv(:), flag_gageii(:)  ! Arrays for number of values per station and GAGE-II validation flags
  real, allocatable               :: Qclmt_full(:), slp_full(:), KKobs_full(:), KImodel_full(:)  
  ! climatology discharge (Qclmt), slope (slp), observed K factor (KKobs), and initial appromaxtion of modeled K factor (KImodel)
  real, allocatable               :: Qclmt(:), slp(:), KKobs(:), KImodel(:)  
  ! Arrays for filtered parameters after station selection
  real, allocatable               :: KImodel_all(:)  ! Array to store all modeled K values from the dataset
  real, allocatable               :: lats_full(:), lons_full(:)  ! Arrays to store the full latitude and longitude values for the grid
  real*8, allocatable             :: MU_axis(:), slp_axis(:), clmt_axis(:), p_axis(:)  
  ! Arrays for parameter axes: M factor (MU), slope exponent (slp_axis), runoff exponent (clmt_axis), and an additional parameter axis (p_axis)

  ! Set model parameters and scaling factors
  real    :: mm = 0.35, MU = 0.45, exp_slp = 0.2, exp_clmt = -0.2  ! Base model parameters (mm, MU, and exponents for slope and climatology discharge)
  !real :: mm=0.4, MU=0.1, exp_slp=0.5, exp_clmt=0.2  ! Alternative parameter set (commented out)
  real    :: fac_str = 1.  ! Scaling factor for stream K

  ! Declare additional integer and real variables for looping and statistical calculations
  integer :: nt, ns, np, i, j, k, p, count
  real    :: ccr(10,10,10), rms(10,10,10)  ! 3D arrays to hold correlation coefficients and RMS errors over a parameter space
  !real :: ccr(20,10), rms(20,10)  ! Alternative array dimensions (commented out)
  real    :: ccrp, rmsp  ! Variables to store computed correlation coefficient and RMS error for a given parameter set
 
  character(len=900) :: file_vel
  character(len=900) :: file_dis
  character(len=900) :: file_usid
  character(len=900) :: file_lats
  character(len=900) :: file_lons
  character(len=900) :: file_lat1m
  character(len=900) :: file_lon1m
  character(len=900) :: file_pfafmap
  character(len=900) :: file_gage_id 
  character(len=900) :: file_gage_acar

  if (command_argument_count() /= 10) then
      print *, "no appropriate files found"
      stop
  endif
  call get_command_argument(1, file_vel)
  call get_command_argument(2, file_dis)
  call get_command_argument(3, file_usid)
  call get_command_argument(4, file_lats)
  call get_command_argument(5, file_lons)
  call get_command_argument(6, file_lat1m)
  call get_command_argument(7, file_lon1m)
  call get_command_argument(8, file_pfafmap)  
  call get_command_argument(9, file_gage_id)
  call get_command_argument(10, file_gage_acar)  

  ! Read USGS data and process it
  call read_usgs_data(file_vel, file_dis, nl, data)  ! Read USGS data (nl records) into the 2D array "data"
  call process_usgs_data(file_usid, nl, ns, data, nv, nt, vel, dis)  
  ! Process the USGS data to extract the number of stations (ns), velocity, distance

  ! Determine the nearest grid coordinates for each station based on the full grid latitude and longitude arrays
  call find_nearest_coords(file_lats, file_lons, file_lat1m, file_lon1m, ns, nlat, nlon, lats_full, lons_full, lati, loni)

  ! Allocate arrays for parameter axes (each with 10 discrete values)
  allocate(MU_axis(10), slp_axis(10), clmt_axis(10))
  ! Initialize the correlation and RMS error arrays with a default invalid value
  ccr = -9999.
  rms = -9999.
  count = 0
  ! Set up the parameter axis for MU (M factor): values from 0 to 0.45 in increments of 0.05
  do k = 1, 10
    MU_axis(k) = (k - 1) * 0.05
  enddo
  ! Set up the parameter axis for slope exponent: values from 0 to 0.9 in increments of 0.1
  do i = 1, 10
    slp_axis(i) = (i - 1) * 0.1 
  enddo
  ! Set up the parameter axis for climate exponent: values from -0.8 to 1.2 in increments of 0.2
  do j = 1, 10
    clmt_axis(j) = (j - 1) * 0.2 - 0.8
  enddo

  !do k=1,10
  !do i=1,10
  !do j=1,10 

 ! count = count + 1  ! Increment the count of parameter combinations (currently only one iteration)

  !MU = MU_axis(k)
  !exp_slp = slp_axis(i)
  !exp_clmt = clmt_axis(j)

!  print *, "count=", count
  print *, "M=", MU, ", exp_slp=", exp_slp, ", exp_clmt=", exp_clmt

  ! Retrieve station information and associated parameter data based on grid indices and model parameters
  call get_station_inf(file_pfafmap, ns, nc, nlat, nlon, lati, loni, catid_full, Qclmt_full, slp_full, KImodel_all, exp_slp, exp_clmt, fac_str)
  ! filtering stations using the GAGE-II dataset criteria
  call get_valide_stations_gageii(file_gage_id, file_gage_acar, ns, nc, catid_full, flag_gageii)
  ! Perform regression analysis using the USGS data
  call regression(nt, vel, dis, nv, ns, Qclmt_full, slp_full, KKobs_full, KImodel_full, exp_slp, exp_clmt, mm, MU)
  ! Filter stations based on predefined criteria 
  call filter_station(nc, ns, np, lats_full, lons_full, Qclmt_full, slp_full, catid_full, KKobs_full, KImodel_full, Qclmt, slp, catid, KKobs, KImodel, flag_gageii)
  ! Calculate the modeled K parameter for each station
  !call cal_Kmodel(ns, np, nc, MU, exp_slp, exp_clmt, Qclmt, slp, KKobs, KImodel, KImodel_all, catid, catid_full, ccr(k,i,j), rms(k,i,j))
  call cal_Kmodel(ns, np, nc, MU, exp_slp, exp_clmt, Qclmt, slp, KKobs, KImodel, KImodel_all, catid, catid_full, ccrp, rmsp)
  
  ! Print the computed correlation coefficient and RMS error
  print *, "ccr=", ccrp
  print *, "rms=", rmsp

  !enddo
  !enddo
  !enddo

  ! The following calls would write the 3D parameter space results to NetCDF files (currently commented out)
  !call create_ncfile_real3d("ccr_clmtxslpxMU_10x10x10_mm0p35.nc", "data", ccr, MU_axis, slp_axis, clmt_axis, 10, 10, 10)
  !call create_ncfile_real3d("rms_clmtxslpxMU_10x10x10_mm0p35.nc", "data", rms, MU_axis, slp_axis, clmt_axis, 10, 10, 10)

end program main