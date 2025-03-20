program main

  use k_module

  implicit none

  integer, parameter :: nl = 3352492
  integer, parameter :: nlat = 10800, nlon = 21600
  integer, parameter :: nc = 291284

  ! Declare variables
  integer, allocatable :: lati(:), loni(:)
  real, allocatable :: data(:, :)  ! 2D data array
  integer,allocatable :: catid_full(:),catid(:)
  real,allocatable, dimension(:) :: vel, dis
  integer, allocatable :: nv(:),flag_gageii(:)
  real,allocatable :: Qclmt_full(:),slp_full(:),KKobs_full(:),KImodel_full(:)
  real,allocatable :: Qclmt(:),slp(:),KKobs(:),KImodel(:)
  real,allocatable :: KImodel_all(:)
  real,allocatable :: lats_full(:), lons_full(:)  
  real*8,allocatable :: MU_axis(:),slp_axis(:),clmt_axis(:),p_axis(:)

  real :: mm=0.35, MU=0.45, exp_slp=0.2, exp_clmt=-0.2 !MU ~(-0.6)
  !real :: mm=0.4, MU=0.1, exp_slp=0.5, exp_clmt=0.2
  real :: fac_str=1.  

  integer :: nt,ns,np,i,j,k,p,count
  real :: ccr(10,10,10),rms(10,10,10)
  !real :: ccr(20,10),rms(20,10)
  real :: ccrp, rmsp
 

  call read_usgs_data(nl, data)
  call process_usgs_data(nl, ns, data, nv, nt, vel, dis) 
  !stop
  call find_nearest_coords(ns, nlat, nlon, lats_full, lons_full, lati, loni)

  allocate(MU_axis(10),slp_axis(10),clmt_axis(10))

  ccr=-9999.
  rms=-9999.
  count=0

  do k=1,10
    MU_axis(k)=(k-1)*0.05
  enddo
  do i=1,10
    slp_axis(i)=(i-1)*0.1 
  enddo
  do j=1,10
    clmt_axis(j)=(j-1)*0.2-0.8
  enddo

  !do k=1,10
  !do i=1,10
  !do j=1,10 

  count=count+1

  !MU=MU_axis(k)
  !exp_slp=slp_axis(i)
  !exp_clmt=clmt_axis(j)


  print *,"count=",count
  print *,"M=",MU,", exp_slp=",exp_slp,", exp_clmt=",exp_clmt

  call get_station_inf(ns, nc, nlat, nlon, lati, loni, catid_full, Qclmt_full, slp_full, KImodel_all,exp_slp,exp_clmt,fac_str)
  call get_valide_stations_gageii(ns,nc,catid_full,flag_gageii)
  call regression(nt,vel,dis,nv,ns,Qclmt_full,slp_full,KKobs_full,KImodel_full,exp_slp,exp_clmt,mm,MU) 
  call filter_station(nc,ns,np,lats_full,lons_full,Qclmt_full,slp_full,catid_full,KKobs_full,KImodel_full,Qclmt,slp,catid,KKobs,KImodel,flag_gageii)
  !call cal_Kmodel(ns,np,nc,MU,exp_slp,exp_clmt,Qclmt,slp,KKobs,KImodel,KImodel_all,catid,catid_full,ccr(k,i,j),rms(k,i,j))
  call cal_Kmodel(ns,np,nc,MU,exp_slp,exp_clmt,Qclmt,slp,KKobs,KImodel,KImodel_all,catid,catid_full,ccrp,rmsp)
  
  print *,"ccr=",ccrp
  print *,"rms=",rmsp

  !enddo
  !enddo
  !enddo



  !call create_ncfile_real3d("ccr_clmtxslpxMU_10x10x10_mm0p35.nc","data",ccr,MU_axis,slp_axis,clmt_axis,10,10,10)
  !call create_ncfile_real3d("rms_clmtxslpxMU_10x10x10_mm0p35.nc","data",rms,MU_axis,slp_axis,clmt_axis,10,10,10) 



end program main