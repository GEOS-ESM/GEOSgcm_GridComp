module reservoir

use rwncfile

implicit none
private
public :: res_init, res_cal

!----Reservoir module constants----------

real*8, parameter  :: fac_elec_a = 0.30D0 ! Coefficient for hydropower calculation
real*8, parameter  :: fac_elec_b = 2.00D0 ! Exponent for hydropower calculation
real*8, parameter  :: fac_irr_a = 0.01D0 ! Coefficient for irrigation calculation (arid areas)
real*8, parameter  :: fac_irr_b = 3.00D0  ! Scaling factor for irrigation (arid areas)
real*8, parameter  :: fac_sup_a = 0.03D0  ! Coefficient for water supply calculation
real*8, parameter  :: fac_sup_b = 2.00D0  ! Exponent for water supply calculation
real*8, parameter  :: fac_other_a = 0.20D0 ! Coefficient for other reservoir types
real*8, parameter  :: fac_other_b = 2.00D0 ! Exponent for other reservoir types
integer, parameter :: fac_fld = 1         ! Flood control parameter

real*8, parameter :: dt = 86400.D0        ! Time step in seconds (1 day)

real*8, parameter :: ai_thres = 0.5D0     ! Aridity index threshold for irrigation reservoirs
real*8, parameter :: rho = 1.D3           ! Water density (kg/m^3)

!-----------------------------------------

contains

!------------------------------------------
! Initialization subroutine for reservoirs
subroutine res_init(input_dir,nres,nc,use_res,active_res,Wr_res,Q_res,type_res,cap_res,Qavg_res,ai_res,fld_res,Qfld_thres,irr_sea_frac,cat2res,wid_res)
  character(len=500),intent(in) :: input_dir
  ! Define the number of reservoirs (nres) and the number of catchments (nc)
  integer,intent(in) :: nres,nc
  ! Logical variable to check if reservoirs are used
  logical,intent(in) :: use_res
  ! Input/output arrays for reservoir attributes: active reservoirs, types, capacities, etc.
  integer,intent(inout),allocatable :: active_res(:),type_res(:),fld_res(:),cat2res(:)
  real*8,intent(inout),allocatable :: Wr_res(:),Q_res(:),cap_res(:),Qavg_res(:),ai_res(:),Qfld_thres(:),irr_sea_frac(:,:)
  real*8,intent(inout),allocatable :: wid_res(:)

  ! Internal arrays for various reservoir-related data
  integer,allocatable,dimension(:) :: flag_grand,catid_grand,elec_grand,irrsup_grand,fld_grand,supply_grand,irr_grand,realuse_grand
  integer,allocatable,dimension(:) :: nav_grand,rec_grand,other_grand
  real*8,allocatable,dimension(:) :: cap_grand,area_max_res,Qavg_grand,ai_grand,area_grand,power_grand,area_res
  real*8,allocatable,dimension(:,:) :: Wres_tar

  ! Define the flood threshold variable and a counter variable
  character(len=2) :: fld_thres  
  integer :: i,cid,rid

!----------reservoir module--------------
  ! Allocate memory for each array
  allocate(flag_grand(nres),catid_grand(nres),active_res(nc))
  allocate(Wr_res(nc),Q_res(nc))
  allocate(elec_grand(nres),type_res(nc),cap_grand(nres),cap_res(nc),area_grand(nres))
  allocate(area_res(nc),area_max_res(nc))
  allocate(irrsup_grand(nres))
  allocate(fld_grand(nres),fld_res(nc),Qfld_thres(nc),supply_grand(nres))
  allocate(irr_grand(nres))
  allocate(cat2res(nc))
  allocate(nav_grand(nres),rec_grand(nres))
  allocate(other_grand(nres))
  allocate(wid_res(nc))
  allocate(realuse_grand(nres))

  ! Open reservoir-related data files and read the corresponding arrays
  open(77,file=trim(input_dir)//"/catid_dam_corr_aca_grand5000.txt")
  read(77,*)catid_grand
  open(77,file=trim(input_dir)//"/flag_all_res.txt")
  read(77,*)flag_grand
  open(77,file=trim(input_dir)//"/cap_max_grand.txt")
  read(77,*)cap_grand
  cap_grand=cap_grand*1.D6*rho ! Convert capacity from million cubic meters (MCM) to kilograms (kg)
  open(77,file=trim(input_dir)//"/hydroelec_grand.txt")
  read(77,*)elec_grand
  !open(77,file=trim(input_dir)//"/Qavg_res_2016_2020_OL7000.txt")
  !read(77,*)Qavg_grand
  !Qavg_grand=Qavg_grand*rho ! Convert flow rate from cubic meters per second (m3/s) to kilograms per second (kg/s)
  !open(77,file=trim(input_dir)//"/ai_grand.txt")
  !read(77,*)ai_grand
  open(77,file=trim(input_dir)//"/irrmainsec_noelec_grand.txt")
  read(77,*)irrsup_grand
  open(77,file=trim(input_dir)//"/fldmainsec_grand.txt")
  read(77,*)fld_grand
  write(fld_thres,'(I2.2)')fac_fld
  open(77,file=trim(input_dir)//"/Pfaf_flood_qr_thres"//trim(fld_thres)//".txt")
  read(77,*)Qfld_thres ! Read flood thresholds in cubic meters per second (m3/s)
  Qfld_thres=Qfld_thres*rho ! Convert threshold from cubic meters per second to kilograms per second (kg/s)
  open(77,file=trim(input_dir)//"/watersupply_grand.txt")
  read(77,*)supply_grand
  open(77,file=trim(input_dir)//"/irr_grand.txt")
  read(77,*)irr_grand
  open(77,file=trim(input_dir)//"/nav_grand.txt")
  read(77,*)nav_grand
  open(77,file=trim(input_dir)//"/rec_grand.txt")
  read(77,*)rec_grand
  open(77,file=trim(input_dir)//"/other_grand.txt")
  read(77,*)other_grand
  open(77,file=trim(input_dir)//"/area_skm_grand.txt")
  read(77,*)area_grand
  area_grand=area_grand*1.D6 ! Convert area from square kilometers (km2) to square meters (m2)
  !open(77,file=trim(input_dir)//"/power_grand.txt")
  !read(77,*)power_grand

  ! Set initial reservoir ID mapping
  cat2res=0
  do i=1,nres
    if(flag_grand(i)==1)then  
      cid=catid_grand(i)
      cat2res(cid)=i ! Link reservoirs with catchments: multiple reservoirs in a catchment share attributes that can be accessed via cat2res
    endif
  enddo

  ! Initialize reservoir properties
  cap_res = 0.D0           ! Set reservoir capacity to zero
  area_res = 0.D0          ! Set reservoir area to zero
  area_max_res = 0.D0      ! Set max reservoir area to zero
  type_res = 0             ! Set reservoir type to zero
  !Qavg_res = 0.D0          ! Set average reservoir flow rate to zero
  !ai_res = 0.D0            ! Set irrigation index to zero
  fld_res = 0              ! Set flood status to zero
  active_res = 0           ! Set active reservoirs to zero
  realuse_grand = 0        ! Initialize real use for each reservoir to zero

  ! Loop over all reservoirs
  do i = 1, nres
    if(flag_grand(i) == 1) then     ! If the reservoir is flagged as active
      cid = catid_grand(i)          ! Get the catchment ID for the reservoir
      cap_res(cid) = cap_res(cid) + cap_grand(i)  ! Sum up the capacities for reservoirs in the same catchment
      area_res(cid) = area_res(cid) + area_grand(i) ! Sum up the areas for reservoirs in the same catchment
      !Qavg_res(cid) = Qavg_grand(i)               ! Assign average flow rate to the catchment
      if(fld_grand(i) == 1) fld_res(cid) = 1      ! Mark the catchment if it has flood control
    endif
  enddo

  ! Compute reservoir width from area (square root of the area)
  wid_res = sqrt(area_res)

  ! Assign reservoir type 7 (Other use) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(other_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 7             ! Type 7 for other uses
        cat2res(cid) = i              ! Map the catchment to the reservoir
        area_max_res(cid) = area_grand(i) ! Update the maximum area for the catchment
      endif
    endif
  enddo

  ! Assign reservoir type 6 (Recreational use) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(rec_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 6
        cat2res(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 5 (Navigational use) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(nav_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 5
        cat2res(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 4 (Water supply) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(supply_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 4
        cat2res(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 3 (Irrigation) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(irr_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 3
        cat2res(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 2 (Electricity generation) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(elec_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 2
        cat2res(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 1 (Irrigation supply) with specific conditions
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(irrsup_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res(cid) = 1            ! Assign type 1 for irrigation supply
        !ai_res(cid) = ai_grand(i)     ! Assign irrigation index to the catchment
        cat2res(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Mark active reservoirs based on type or flood control status
  do i = 1, nc
    if(type_res(i) /= 0 .or. fld_res(i) == 1) then
      active_res(i) = 1
    endif
  enddo

  ! Assign real reservoir usage based on type, with error checking
  do i = 1, nres
    if(flag_grand(i) == 1) then 
      cid = catid_grand(i)
      rid = cat2res(cid)
      if(rid > 0) then
        if(type_res(cid) == 0 .and. fld_res(cid) == 0) then
          print *, "type_res(cid) == 0"
          stop
        endif
        if(type_res(cid) == 0) then
          realuse_grand(i) = -1    ! Invalid reservoir use type
        else
          realuse_grand(i) = type_res(cid)  ! Assign the actual use type
        endif
      else
        print *, "rid == 0"
        stop
      endif
    endif
  enddo

  ! Read irrigation and reservoir target data from NetCDF files
  ! call read_ncfile_double2d(trim(input_dir)//"/irr_grand_frac.nc", "data", irr_sea_frac, nres, 12)
  ! call read_ncfile_double2d(trim(input_dir)//"/Wr_tar_Dang.nc", "data", Wres_tar, 365, nres)

  ! Wres_tar = Wres_tar * 1.D6 * rho  ! Convert from million cubic meters (MCM) to kilograms (kg)

  ! Deactivate reservoirs if the use_res flag is set to False
  if(use_res == .False.) active_res = 0

end subroutine res_init

!-----------------------
! Reservoir calculation subroutine
subroutine res_cal(active_res,active_lake,Qout,Q_lake,type_res,cat2res,Q_res,wid_res,fld_res,Wr_res,Qfld_thres,cap_res,B1,B2)
  integer, intent(in) :: active_res, type_res, active_lake, cat2res, fld_res
  real*8, intent(in) :: Qout, Q_lake, wid_res, Qfld_thres, cap_res
  real*8, intent(inout) :: Q_res, Wr_res, B1, B2

  integer :: rid  ! Reservoir ID
  real*8 :: Qin_res, coe, irrfac, alp_res  ! Variables for inflow, coefficients, and factors

  ! If the reservoir is active
  if (active_res == 1) then

    ! Determine the inflow to the reservoir (from river or lake)
    if (active_lake == 0) then
      Qin_res = Qout  ! Inflow from river
    else
      Qin_res = Q_lake  ! Inflow from lake
    endif

    ! Irrigation reservoir 
    if (type_res == 1 .or. type_res == 3) then 
      alp_res = fac_irr_a * ((1.D0 / (wid_res / 1.D3)) ** fac_irr_b) / 3600.D0   ! irrigation coefficient
      Q_res = alp_res * Wr_res  ! Outflow based on water storage

    ! Hydropower reservoir
    else if (type_res == 2) then 
      alp_res = fac_elec_a * ((1.D0 / (wid_res / 1.D3)) ** fac_elec_b) / 3600.D0  ! Hydropower coefficient
      Q_res = alp_res * Wr_res  ! Outflow based on water storage

    ! Water supply reservoir
    else if (type_res == 4) then 
      alp_res = fac_sup_a * ((1.D0 / (wid_res / 1.D3)) ** fac_sup_b) / 3600.D0  ! Supply coefficient
      Q_res = alp_res * Wr_res  ! Outflow based on water storage

    ! Other reservoir types
    else if (type_res == 5 .or. type_res == 6 .or. type_res == 7 .or. type_res == 0) then 
      alp_res = fac_other_a * ((1.D0 / (wid_res / 1.D3)) ** fac_other_b) / 3600.D0  ! Generic reservoir coefficient
      Q_res = alp_res * Wr_res  ! Outflow based on water storage
    endif

    ! Ensure outflow is within reasonable bounds
    Q_res = max(0.D0, Q_res)  ! Ensure non-negative outflow
    Q_res = min(Q_res, Wr_res / dt + Qin_res)  ! Limit outflow to prevent exceeding inflow and storage
    if (fld_res == 1) Q_res = min(Q_res, Qfld_thres)  ! Limit outflow for flood control
    Wr_res = Wr_res + dt * (Qin_res - Q_res)  ! Update water storage in the reservoir
    Wr_res = max(0.D0, Wr_res)  ! Ensure non-negative storage

    ! If the storage exceeds capacity, adjust outflow and storage
    if (Wr_res > cap_res) then
      Q_res = Q_res + (Wr_res - cap_res) / dt  ! Adjust outflow for overflow
      Wr_res = cap_res  ! Limit storage to reservoir capacity
    endif

    ! Output the calculated outflow and zero out the second output variable (B2)
    B1 = Q_res
    B2 = 0.D0

  endif

end subroutine res_cal

end module reservoir
