module reservoir


implicit none
private
public :: res_init, res_cal

!----Reservoir module constants----------
integer,parameter :: nres=7250
integer,parameter :: nlake=3917

real, parameter  :: fac_elec_a = 0.30 ! Coefficient for hydropower calculation
real, parameter  :: fac_elec_b = 2.00 ! Exponent for hydropower calculation
real, parameter  :: fac_irr_a = 0.01 ! Coefficient for irrigation calculation (arid areas)
real, parameter  :: fac_irr_b = 3.00  ! Scaling factor for irrigation (arid areas)
real, parameter  :: fac_sup_a = 0.03  ! Coefficient for water supply calculation
real, parameter  :: fac_sup_b = 2.00  ! Exponent for water supply calculation
real, parameter  :: fac_other_a = 0.20 ! Coefficient for other reservoir types
real, parameter  :: fac_other_b = 2.00 ! Exponent for other reservoir types
integer, parameter :: fac_fld = 1         ! Flood control parameter

real, parameter  :: fac_a_slake = 0.003  ! Factor for small lakes
real, parameter  :: fac_b_slake = 0.40   ! Exponent for small lakes
real, parameter  :: fac_a_llake = 0.01   ! Factor for large lakes
real, parameter  :: fac_b_llake = 0.60   ! Exponent for large lakes
real, parameter  :: thr_wid_lake = 1.e5    ! Threshold lake width (in m)

!real, parameter :: dt = 86400.       ! Time step in seconds (1 day)
real, parameter :: rho = 1.e3         ! Water density (kg/m^3)

!-----------------------------------------

contains

!------------------------------------------
! Initialization subroutine for reservoirs
subroutine res_init(input_dir,nall,nc,minCatch,maxCatch,use_res,active_res,type_res,cap_res,fld_res,Qfld_thres,cat2res,wid_res)
  character(len=500),intent(in) :: input_dir
  ! Define the number of reservoirs (nres) and the number of catchments (nc)
  integer,intent(in) :: nall,nc,minCatch,maxCatch
  ! Logical variable to check if reservoirs are used
  logical,intent(in) :: use_res
  ! Input/output arrays for reservoir attributes: active reservoirs, types, capacities, etc.
  integer,intent(inout),pointer :: active_res(:),type_res(:),fld_res(:),cat2res(:)
  real,intent(inout),pointer :: cap_res(:),Qfld_thres(:)
  real,intent(inout),pointer :: wid_res(:)

  ! Internal arrays for various reservoir-related data
  integer,allocatable,dimension(:) :: flag_grand,catid_grand,elec_grand,fld_grand,supply_grand,irr_grand,realuse_grand
  integer,allocatable,dimension(:) :: nav_grand,rec_grand,other_grand
  integer,allocatable,dimension(:) :: type_res_all,cat2res_all
  real,allocatable,dimension(:) :: cap_grand,area_max_res,Qavg_grand,ai_grand,area_grand,power_grand,area_res
  real,allocatable,dimension(:,:) :: Wres_tar
  real,pointer :: buff_global(:)=>NULL(),area_all(:)=>NULL()
  integer,pointer :: fld_all(:)=>NULL() !buff_global_int(:)=>NULL()
  real :: value_max

  integer,allocatable,dimension(:) :: flag_lake,catid_lake
  real,allocatable,dimension(:) :: area_lake

  ! Define the flood threshold variable and a counter variable
  character(len=2) :: fld_thres  
  integer :: i,cid,rid

!----------reservoir module--------------
  ! Allocate memory for each array
  allocate(flag_grand(nres),catid_grand(nres),active_res(nc),Qfld_thres(nc))
  allocate(elec_grand(nres),type_res(nc),type_res_all(nall),cap_grand(nres),cap_res(nc),area_grand(nres))
  allocate(area_res(nc),area_max_res(nall))
  allocate(fld_grand(nres),fld_res(nc),supply_grand(nres))
  allocate(irr_grand(nres))
  allocate(cat2res(nc),cat2res_all(nall))
  allocate(nav_grand(nres),rec_grand(nres))
  allocate(other_grand(nres))
  allocate(wid_res(nc))
  allocate(realuse_grand(nres))

  allocate(flag_lake(nlake),catid_lake(nlake),area_lake(nlake))

  ! Open reservoir-related data files and read the corresponding arrays
  open(77,file=trim(input_dir)//"/catid_dam_corr_aca_grand5000.txt",status="old",action="read")
  read(77,*)catid_grand;close(77)
  open(77,file=trim(input_dir)//"/flag_all_res.txt",status="old",action="read")
  read(77,*)flag_grand;close(77)
  open(77,file=trim(input_dir)//"/cap_max_grand.txt",status="old",action="read")
  read(77,*)cap_grand;close(77)
  cap_grand=cap_grand*1.e6! Convert capacity from million cubic meters (MCM) to m3
  open(77,file=trim(input_dir)//"/hydroelec_grand.txt",status="old",action="read")
  read(77,*)elec_grand;close(77)
  open(77,file=trim(input_dir)//"/fldmainsec_grand.txt",status="old",action="read")
  read(77,*)fld_grand;close(77)
  write(fld_thres,'(I2.2)')fac_fld

  open(77,file=trim(input_dir)//"/watersupply_grand.txt",status="old",action="read")
  read(77,*)supply_grand;close(77)
  open(77,file=trim(input_dir)//"/irr_grand.txt",status="old",action="read")
  read(77,*)irr_grand;close(77)
  open(77,file=trim(input_dir)//"/nav_grand.txt",status="old",action="read")
  read(77,*)nav_grand;close(77)
  open(77,file=trim(input_dir)//"/rec_grand.txt",status="old",action="read")
  read(77,*)rec_grand;close(77)
  open(77,file=trim(input_dir)//"/other_grand.txt",status="old",action="read")
  read(77,*)other_grand;close(77)
  open(77,file=trim(input_dir)//"/area_skm_grand.txt",status="old",action="read")
  read(77,*)area_grand;close(77)
  area_grand=area_grand*1.e6 ! Convert area from square kilometers (km2) to square meters (m2)

  allocate(buff_global(nall))
  !open(77,file=trim(input_dir)//"/Pfaf_flood_qr_thres"//trim(fld_thres)//".txt");read(77,*)buff_global;close(77)
  Qfld_thres=0.!buff_global(minCatch:maxCatch)
  deallocate(buff_global) 

  !lake input
  open(77, file = trim(input_dir)//"/lake_outlet_flag_valid_2097.txt")
  read(77, *) flag_lake;close(77)
  open(77, file = trim(input_dir)//"/lake_outlet_catid.txt")
  read(77, *) catid_lake;close(77)
  open(77, file = trim(input_dir)//"/lake_outlet_lakearea.txt")
  read(77, *) area_lake;close(77) ! km^2  
  area_lake=area_lake*1.e6 

  ! Set initial reservoir ID mapping 
  cat2res_all=0
  do i=1,nres
    if(flag_grand(i)==1)then  
      cid=catid_grand(i)
      cat2res_all(cid)=i ! Link reservoirs with catchments: multiple reservoirs in a catchment share attributes that can be accessed via cat2res
    endif
  enddo

  ! Initialize reservoir properties
  cap_res = 0.0           ! Set reservoir capacity to zero
  area_res = 0.0          ! Set reservoir area to zero
  area_max_res = 0.0      ! Set max reservoir area to zero
  type_res_all = 0             ! Set reservoir type to zero
  fld_res = 0              ! Set flood status to zero
  active_res = 0           ! Set active reservoirs to zero
  realuse_grand = 0        ! Initialize real use for each reservoir to zero

  ! Loop over all reservoirs
  allocate(buff_global(nall),fld_all(nall),area_all(nall))
  buff_global=0.
  area_all=0.
  fld_all=0
  do i = 1, nres
    if(flag_grand(i) == 1) then     ! If the reservoir is flagged as active
      cid = catid_grand(i)          ! Get the catchment ID for the reservoir
      buff_global(cid) = buff_global(cid) + cap_grand(i)  ! Sum up the capacities for reservoirs in the same catchment
      area_all(cid) = area_all(cid) + area_grand(i) ! Sum up the areas for reservoirs in the same catchment
      !Qavg_res(cid) = Qavg_grand(i)               ! Assign average flow rate to the catchment
      if(fld_grand(i) == 1) fld_all(cid) = 1      ! Mark the catchment if it has flood control
    endif
  enddo
  cap_res=buff_global(minCatch:maxCatch)
  value_max=huge(value_max)
  where(cap_res==0.) cap_res=value_max
  !area_res=buff_global2(minCatch:maxCatch)
  fld_res=fld_all(minCatch:maxCatch)
  deallocate(buff_global)

  ! Assign reservoir type 6 (Other use) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(other_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res_all(cid) = 6             
        cat2res_all(cid) = i              ! Map the catchment to the reservoir
        area_max_res(cid) = area_grand(i) ! Update the maximum area for the catchment
      endif
    endif
  enddo

  ! Assign reservoir type 5 (Recreational use) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(rec_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res_all(cid) = 5
        cat2res_all(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 4 (Navigational use) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(nav_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res_all(cid) = 4
        cat2res_all(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 3 (Water supply) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(supply_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res_all(cid) = 3
        cat2res_all(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 2 (Electricity generation) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(elec_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res_all(cid) = 2
        cat2res_all(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Assign reservoir type 1 (Irrigation) to the largest reservoir in a catchment
  do i = 1, nres
    if(flag_grand(i) == 1) then
      cid = catid_grand(i)
      if(irr_grand(i) == 1 .and. area_grand(i) >= area_max_res(cid)) then
        type_res_all(cid) = 1
        cat2res_all(cid) = i
        area_max_res(cid) = area_grand(i)
      endif
    endif
  enddo

  ! Set up natural lakes
  do i = 1, nlake
    if(flag_lake(i) == 1 .and. catid_lake(i) > 0) then
      cid = catid_lake(i)
      if(type_res_all(cid)==0.and.fld_all(cid)==0)then
        type_res_all(cid) = -1 !for lake
        cat2res_all(cid) = i
        area_all(cid) = area_lake(i)
      endif
    endif
  enddo

  type_res=type_res_all(minCatch:maxCatch)
  cat2res=cat2res_all(minCatch:maxCatch)
  area_res=area_all(minCatch:maxCatch)
  ! Compute reservoir width from area (square root of the area)
  wid_res = sqrt(area_res)!m

  ! Mark active reservoirs based on type or flood control status
  do i = 1, nc
    if(type_res(i) /= 0 .or. fld_res(i) == 1) then
      active_res(i) = 1
    endif
  enddo

  ! Deactivate reservoirs if the use_res flag is set to False
  if(use_res == .False.) active_res = 0

  deallocate(flag_grand,catid_grand,elec_grand,type_res_all,cap_grand,area_grand)
  deallocate(area_res,area_max_res,fld_grand,supply_grand,irr_grand)
  deallocate(cat2res_all,nav_grand,rec_grand,other_grand,realuse_grand)
  deallocate(flag_lake,catid_lake,area_lake,area_all,fld_all)

end subroutine res_init

!-----------------------
! Reservoir calculation subroutine
subroutine res_cal(active_res,Qout,type_res,cat2res,Q_res,wid_res,fld_res,Wr_res,Qfld_thres,cap_res,dt)
  integer, intent(in) :: active_res, type_res, cat2res, fld_res
  real, intent(in) :: Qout, wid_res, Qfld_thres, cap_res
  real, intent(inout) :: Q_res, Wr_res
  real, intent(in) :: dt

  integer :: rid  ! Reservoir ID
  real :: Qin_res, coe, irrfac, alp_res  ! Variables for inflow, coefficients, and factors

  ! If the reservoir is active
  if (active_res == 1) then

    ! Determine the inflow to the reservoir 
    Qin_res = Qout  ! Inflow from river

    ! Irrigation reservoir 
    if (type_res == 1) then 
      alp_res = fac_irr_a * ((1.0 / (wid_res / 1.e3)) ** fac_irr_b) / 3600.0   ! irrigation coefficient

    ! Hydropower reservoir
    else if (type_res == 2) then 
      alp_res = fac_elec_a * ((1.0 / (wid_res / 1.e3)) ** fac_elec_b) / 3600.0  ! Hydropower coefficient

    ! Water supply reservoir
    else if (type_res == 3) then 
      alp_res = fac_sup_a * ((1.0 / (wid_res / 1.e3)) ** fac_sup_b) / 3600.0  ! Supply coefficient

    ! Other reservoir types
    else if (type_res == 4 .or. type_res == 5 .or. type_res == 6 .or. type_res == 0) then 
      alp_res = fac_other_a * ((1.0 / (wid_res / 1.e3)) ** fac_other_b) / 3600.0  ! Generic reservoir coefficient

    ! Natural lake 
    else if (type_res == -1) then  
      ! Determine lake type based on area and calculate alpha
      if (wid_res >= thr_wid_lake) then
        alp_res = fac_a_llake * ( (1. / (wid_res / 1.e3)) ** fac_b_llake ) / 3600.
      else
        alp_res = fac_a_slake * ( (1./ (wid_res / 1.e3)) ** fac_b_slake ) / 3600.
      endif

    endif

    Q_res = alp_res * Wr_res

    ! Ensure outflow is within reasonable bounds
    Q_res = max(0.0, Q_res)  ! Ensure non-negative outflow
    Q_res = min(Q_res, Wr_res / dt + Qin_res)  ! Limit outflow to prevent exceeding inflow and storage
    !if (fld_res == 1) Q_res = min(Q_res, Qfld_thres)  ! Limit outflow for flood control
    Wr_res = Wr_res + dt * (Qin_res - Q_res)  ! Update water storage in the reservoir
    Wr_res = max(0.0, Wr_res)  ! Ensure non-negative storage

    ! If the storage exceeds capacity, adjust outflow and storage
    if (Wr_res > cap_res) then
      Q_res = Q_res + (Wr_res - cap_res) / dt  ! Adjust outflow for overflow
      Wr_res = cap_res  ! Limit storage to reservoir capacity
    endif

  endif

end subroutine res_cal

end module reservoir