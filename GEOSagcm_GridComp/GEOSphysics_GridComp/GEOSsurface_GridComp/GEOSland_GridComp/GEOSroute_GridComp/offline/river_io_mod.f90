module river_io

use interp
use rwncfile

implicit none
private

public :: read_input,read_restart,read_runoff,write_output

real*8, parameter :: rho = 1.D3      ! Water density in kg/m^3
character(len=500) :: input_dir="/discover/nobackup/yzeng3/work/river_routing_model_offline/input/" ! Directory for input files
character(len=500) :: output_dir="/discover/nobackup/yzeng3/river_output/" ! Directory for output files
character(len=500) :: runoff_dir="/discover/nobackup/yzeng3/GEOldas_output/" ! Directory for runoff files

integer :: nlon=964 !for M36, change to 3856 for M09
integer :: nlat=406 !for M36, change to 1624 for M09

contains

!------------------------------
subroutine read_input(nc,ny,upmax,days_in_year,fac_kstr,qstr_clmt,qri_clmt,nts,upID,nup,llc_ori,lstr,qin_clmt,K,Kstr,days_acc_year,days_acc_noleap,days_acc_leap,inputdir)
  ! Input parameters:
  integer,intent(in) :: nc, ny, upmax  ! nc: number of catchments, ny: number of years, upmax: max number of upstream catchments
  integer,intent(in) :: days_in_year(ny)  ! Array of days in each year
  real*8,intent(in) :: fac_kstr  ! Scaling factor for streamflow
  real*8,intent(out) :: qstr_clmt(nc), qri_clmt(nc)  ! Climate streamflow (qstr_clmt) and routing inflow (qri_clmt) in kg/s
  integer,intent(out) :: nts(nc), upID(upmax,nc), nup(nc)  ! Number of time steps, upstream IDs, and number of upstream catchments
  real*8,intent(out) :: llc_ori(nc), lstr(nc), qin_clmt(nc), K(nc), Kstr(nc)  ! Original stream length (llc_ori), stream length (lstr), climate inflow (qin_clmt), and hydraulic parameters (K, Kstr)
  integer,intent(out) :: days_acc_year(ny), days_acc_noleap(12), days_acc_leap(12)  ! Accumulated days in regular and leap years
  character(len=500),intent(out) :: inputdir

  ! Days in each month for no-leap and leap years
  integer,dimension(12) :: days_in_mon_noleap=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  integer,dimension(12) :: days_in_mon_leap=(/31,29,31,30,31,30,31,31,30,31,30,31/)
  integer :: i

  inputdir=input_dir
  ! Read input data from files
  open(77,file=trim(input_dir)//"/Pfaf_qstr.txt")
  read(77,*)qstr_clmt  ! Read streamflow climatology (m3/s)
  qstr_clmt=qstr_clmt*rho  ! Convert to kg/s

  open(77,file=trim(input_dir)//"/Pfaf_qri.txt")
  read(77,*)qri_clmt  ! Read routing inflow (m3/s)
  qri_clmt=qri_clmt*rho  ! Convert to kg/s

  open(77,file=trim(input_dir)//"/Pfaf_qin.txt")
  read(77,*)qin_clmt  ! Read climate inflow (m3/s)
  qin_clmt=qin_clmt*rho  ! Convert to kg/s

  open(77,file=trim(input_dir)//"/Pfaf_tosink.txt")
  read(77,*)nts  ! Read number of steps to endpoint

  open(77,file=trim(input_dir)//"/upstream_1D.txt")
  read(77,*)upID  ! Read upstream IDs

  open(77,file=trim(input_dir)//"/Pfaf_upnum.txt")
  read(77,*)nup  ! Read number of upstream catchments

  open(77,file=trim(input_dir)//"/Pfaf_lriv_PR.txt")
  read(77,*)llc_ori  ! Read original stream length (km)
  llc_ori=llc_ori*1.D3  ! Convert km to meters

  open(77,file=trim(input_dir)//"/Pfaf_lstr_PR.txt")
  read(77,*)lstr  ! Read stream length (km)
  lstr=lstr*1.D3  ! Convert km to meters

  open(77,file=trim(input_dir)//"Pfaf_Kv_PR_0p35_0p45_0p2_n0p2.txt")
  read(77,*)K  ! Read hydraulic parameter K

  open(77,file=trim(input_dir)//"Pfaf_Kstr_PR_fac1_0p35_0p45_0p2_n0p2.txt")
  read(77,*)Kstr  ! Read hydraulic parameter Kstr
  Kstr=fac_kstr*Kstr  ! Apply scaling factor to Kstr

  ! Calculate accumulated days for regular years
  days_acc_year(1)=0
  do i=2,ny
    days_acc_year(i)=days_acc_year(i-1)+days_in_year(i-1)
  end do

  ! Calculate accumulated days for no-leap and leap years
  days_acc_noleap(1)=0
  days_acc_leap(1)=0
  do i=2,12
    days_acc_noleap(i)=days_acc_noleap(i-1)+days_in_mon_noleap(i-1)
    days_acc_leap(i)=days_acc_leap(i-1)+days_in_mon_leap(i-1)
  end do

end subroutine read_input
!------------------------------
subroutine read_restart(iter,is_coldstart,ny,nc,days_acc_year,days_acc_noleap,days_acc_leap,Ws,Wr,Wr_res,Wr_lake)
  ! Input parameters:
  integer,intent(in) :: iter  ! Current iteration
  logical,intent(inout) :: is_coldstart  ! Flag for cold start condition
  integer,intent(in) :: ny, nc  ! ny: number of years, nc: number of catchments
  integer,intent(in) :: days_acc_year(ny), days_acc_noleap(12), days_acc_leap(12)  ! Accumulated days for each year and for no-leap/leap years
  real*8,intent(inout) :: Ws(nc), Wr(nc), Wr_res(nc), Wr_lake(nc)  ! Water storage in soil (Ws), routing (Wr), reservoir (Wr_res), and lake (Wr_lake)

  ! Local variables:
  character(len=50) :: iter_s, yr_s, mon_s, day_s  ! Strings for iteration, year, month, and day
  integer :: step_prev, i, yr_cur, mon_cur, day_cur, d_res  ! Step count, loop index, current year, month, day, and day residual
  integer :: days_acc_mon(12)  ! Accumulated days per month

  ! Convert iteration number to string format
  write(iter_s,'(I5.5)')iter
  print *,trim(iter_s)

  ! If first iteration or cold start, read initial data
  if(iter==1.or.is_coldstart)then
    ! Read initial water storage data from files for cold start
    open(77,file=trim(input_dir)//"/Pfaf_Ws_Kv_M0.10_mm0.40_20170330_OL7000.txt")
    read(77,*)Ws  ! Read soil water storage (Ws)
    
    open(77,file=trim(input_dir)//"/Pfaf_Wr_Kv_M0.10_mm0.40_20170330_OL7000.txt")
    read(77,*)Wr  ! Read routing water storage (Wr)
    
    !----reservoir module-------
    open(77,file=trim(input_dir)//"/Pfaf_Wr_res_Kv_M0.10_mm0.40_20170330_OL7000.txt")
    read(77,*)Wr_res  ! Read reservoir water storage (Wr_res)
    
    !----lake module------------
    open(77,file=trim(input_dir)//"/Pfaf_Wr_lake_Kv_M0.10_mm0.40_20170330_OL7000.txt")
    read(77,*)Wr_lake  ! Read lake water storage (Wr_lake)

    ! Set cold start flag to False after initialization
    is_coldstart=.False.

  else
    ! For non-cold start, calculate the current year and day from the previous iteration
    step_prev = iter - 1
    do i = ny, 1, -1
      if(step_prev > days_acc_year(i))then
        yr_cur = 1989 + i  ! Calculate the current year
        d_res = step_prev - days_acc_year(i)  ! Calculate residual days
        exit
      endif
    enddo

    ! Determine whether the current year is a leap year
    if(mod(yr_cur,4) == 0)then
      days_acc_mon = days_acc_leap  ! Use leap year days if it is a leap year
    else
      days_acc_mon = days_acc_noleap  ! Use no-leap year days if it is not a leap year
    endif

    ! Determine the current month and day from the residual days
    do i = 12, 1, -1
      if(d_res > days_acc_mon(i))then
        mon_cur = i  ! Current month
        day_cur = d_res - days_acc_mon(i)  ! Current day
        exit
      endif
    enddo

    ! Convert year, month, and day to string format
    write(yr_s,'(I4)')yr_cur
    write(mon_s,'(I2.2)')mon_cur
    write(day_s,'(I2.2)')day_cur

    ! Read water storage data for the specific date (year, month, day)
    open(77,file=trim(output_dir)//"/Pfaf_Ws_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
    read(77,*)Ws  ! Read soil water storage (Ws)

    open(77,file=trim(output_dir)//"/Pfaf_Wr_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
    read(77,*)Wr  ! Read routing water storage (Wr)

    !----reservoir module-------
    open(77,file=trim(output_dir)//"Pfaf_Wr_res_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
    read(77,*)Wr_res  ! Read reservoir water storage (Wr_res)

    !----lake module------------
    open(77,file=trim(output_dir)//"Pfaf_Wr_lake_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
    read(77,*)Wr_lake  ! Read lake water storage (Wr_lake)

    ! Optionally scale the water storage values (commented out)
    ! Ws = Ws * 1.D9
    ! Wr = Wr * 1.D9
  endif

end subroutine read_restart
!------------------------------
subroutine read_runoff(nc,ny,iter,days_acc_year,days_acc_noleap,days_acc_leap,Qrunf,yr_s,mon_s,day_s,d_res,mon_cur)
  integer,intent(in) :: nc,ny,iter
  integer,intent(in) :: days_acc_year(ny),days_acc_noleap(12),days_acc_leap(12)
  real*8,intent(inout) :: Qrunf(nc)
  character(len=50),intent(inout) :: yr_s,mon_s,day_s
  integer,intent(out) :: d_res,mon_cur

  real*8,allocatable,dimension(:,:,:) :: runoff,runoffr,baseflow  ! Declare 3D arrays for runoff and baseflow

  integer :: i,yr_cur,day_cur
  integer :: days_acc_mon(12)  ! Array to store accumulated days for current month


  ! Determine current year based on iteration days
  do i=ny,1,-1
    if(iter>days_acc_year(i))then
      yr_cur=1989+i  ! Set current year
      d_res=iter-days_acc_year(i)  ! Calculate residual days
      exit
    endif
  enddo

  ! Set days_acc_mon based on whether the current year is a leap year
  if(mod(yr_cur,4)==0)then
    days_acc_mon=days_acc_leap  ! Use leap year days
  else
    days_acc_mon=days_acc_noleap  ! Use non-leap year days
  endif

  ! Determine current month and day based on residual days
  do i=12,1,-1
    if(d_res>days_acc_mon(i))then
      mon_cur=i  ! Set current month
      day_cur=d_res-days_acc_mon(i)  ! Set current day
      exit
    endif
  enddo

  ! Write current year, month, and day as strings
  write(yr_s,'(I4)')yr_cur
  write(mon_s,'(I2.2)')mon_cur
  write(day_s,'(I2.2)')day_cur
  print *,trim(yr_s)," ",trim(mon_s)," ",trim(day_s)

  ! Allocate memory for runoff, runoffr, and baseflow arrays
  allocate(runoff(nlon,nlat,1),runoffr(nlon,nlat,1),baseflow(nlon,nlat,1))

  ! Read runoff and baseflow data from NetCDF files
  call read_ncfile_double3d(trim(runoff_dir)//"/Y"//trim(yr_s)//"/M"//trim(mon_s)//"/SMAP_Nature_v10.0_M36.tavg24_2d_lnd_Nx."//trim(yr_s)//trim(mon_s)//trim(day_s)//"_1200z.nc4","RUNOFF",runoff,nlon,nlat,1)
  call read_ncfile_double3d(trim(runoff_dir)//"/Y"//trim(yr_s)//"/M"//trim(mon_s)//"/SMAP_Nature_v10.0_M36.tavg24_2d_lnd_Nx."//trim(yr_s)//trim(mon_s)//trim(day_s)//"_1200z.nc4","BASEFLOW",baseflow,nlon,nlat,1)

  ! Combine runoff and baseflow, and convert to daily values
  runoff=runoff+baseflow
  runoff=runoff*86400.D0  ! Convert to mm/day

  ! Reverse the y-direction of the runoff array
  do i=1,406
    runoffr(:,i,:)=runoff(:,407-i,:)
  enddo
  runoff=runoffr

  ! Convert from mm/day to kg/s and store in Qrunf
  Qrunf=M36_to_cat(runoff(:,:,1),nlon,nlat,nc,input_dir)

  ! Deallocate the arrays to free memory
  deallocate(runoff,runoffr,baseflow)

  ! The following lines are commented out, but they suggest reading runoff from a text file instead of NetCDF
  !open(77,file="/Users/zsp/Desktop/work/river/OL7000_Pfaf/runoff_"//trim(yr_s)//trim(mon_s)//trim(day_s)//".txt")
  !read(77,*)Qrunf
  !Qrunf=Qrunf*rho !m3/s -> kg/s

end subroutine read_runoff
!------------------------------
subroutine write_output(nc,yr_s,mon_s,day_s,Qout,Ws,Wr,Q_res,Wr_res,Q_lake,Wr_lake)
  integer,intent(in) :: nc
  character(len=50),intent(in) :: yr_s,mon_s,day_s
  real*8,intent(in) :: Qout(nc),Ws(nc),Wr(nc),Q_res(nc),Wr_res(nc),Q_lake(nc),Wr_lake(nc)

  integer :: i

  ! Open file to write Qout (discharge) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Qr_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Qout(i)/1.D3 ! Convert from m^3/s to km^3/s
  enddo

  ! Open file to write Ws (soil water storage) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Ws_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Ws(i) ! Write Ws values, unit in kg
  enddo

  ! Open file to write Wr (river water storage) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Wr_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Wr(i) ! Write Wr values, unit in kg
  enddo

  !-----------reservoir module----------------
  ! Open file to write Q_res (reservoir discharge) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Q_res_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Q_res(i)/1.D3 ! Convert from m^3/s to km^3/s
  enddo

  ! Open file to write Wr_res (reservoir water storage) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Wr_res_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Wr_res(i) ! Write Wr_res values, unit in kg
  enddo
  !-------------------------------------------

  !-----------lake module---------------------
  ! Open file to write Q_lake (lake discharge) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Q_lake_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Q_lake(i)/1.D3 ! Convert from m^3/s to km^3/s
  enddo

  ! Open file to write Wr_lake (lake water storage) values and write to the file
  open(88,file=trim(output_dir)//"/Pfaf_Wr_lake_Kv_"//trim(yr_s)//trim(mon_s)//trim(day_s)//"_OL7000.txt")
  do i=1,nc
    write(88,*)Wr_lake(i) ! Write Wr_lake values, unit in kg
  enddo
  !-------------------------------------------

  ! Print out the sum of Wr (river water storage) in petagrams (10^12 kg)
  print *,"sum of Wr is ", sum(Wr)/1.D12
  ! Print out the sum of Wr_lake (lake water storage) in petagrams (10^12 kg)
  print *,"sum of Wr_lake is ", sum(Wr_lake)/1.D12
  ! Print out the sum of Wr_res (reservoir water storage) in petagrams (10^12 kg)
  print *,"sum of Wr_res is ", sum(Wr_res)/1.D12

end subroutine write_output

end module river_io
