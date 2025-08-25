module lake


implicit none
private
public :: lake_init, lake_cal

! Define parameters for small and large lakes
real*8, parameter  :: fac_a_slake = 0.003D0  ! Factor for small lakes
real*8, parameter  :: fac_b_slake = 0.40D0   ! Exponent for small lakes
real*8, parameter  :: fac_a_llake = 0.01D0   ! Factor for large lakes
real*8, parameter  :: fac_b_llake = 0.60D0   ! Exponent for large lakes
real*8, parameter  :: thr_area_lake = 1D4    ! Threshold lake area (in km^2)

! Define constants
real*8, parameter :: dt = 86400.D0   ! Time step in seconds (1 day)
real*8, parameter :: rho = 1.D3      ! Water density in kg/m^3

contains

!------------------------------
! Initialization subroutine for lakes
subroutine lake_init(input_dir, use_lake, nc, nlake, nres, active_res, active_lake, area_lake, Wr_lake, Q_lake)
  character(len=500),intent(in) :: input_dir
  logical, intent(in) :: use_lake         ! Flag to use lake module
  integer, intent(in) :: nc, nlake, nres  ! Number of catchments, lakes, reservoirs
  integer, intent(in) :: active_res(nres) ! Active reservoirs
  integer, allocatable, intent(inout) :: active_lake(:) ! Active lakes (output)
  real*8, allocatable, intent(inout) :: area_lake(:), Wr_lake(:), Q_lake(:)  ! Lake areas, water storage, outflow

  integer, allocatable :: flag_valid_laked(:), catid_laked(:)
  real*8, allocatable :: area_laked(:)

  integer :: i, cid

  ! Allocate arrays for lake attributes
  allocate(flag_valid_laked(nlake), catid_laked(nlake), area_laked(nlake))
  allocate(active_lake(nc), area_lake(nc))
  allocate(Wr_lake(nc), Q_lake(nc))

  ! Read lake outlet and area data from external files
  open(77, file = trim(input_dir)//"/lake_outlet_flag_valid_2097.txt")
  read(77, *) flag_valid_laked
  open(77, file = trim(input_dir)//"/lake_outlet_catid.txt")
  read(77, *) catid_laked
  open(77, file = trim(input_dir)//"/lake_outlet_lakearea.txt")
  read(77, *) area_laked  ! km^2

  ! Initialize lake attributes to zero
  area_lake = 0.D0
  active_lake = 0

  ! Assign active lakes and their areas based on data
  do i = 1, nlake
    if (flag_valid_laked(i) == 1) then
      cid = catid_laked(i)
      active_lake(cid) = 1
      area_lake(cid) = area_laked(i)
    endif
  enddo

  ! Deactivate lakes where reservoirs are active
  where (active_res == 1) active_lake = 0

  ! If lakes are not being used, set active lakes to zero
  if (use_lake .eqv. .False.) active_lake = 0

end subroutine lake_init

!------------------------------
! Calculation subroutine for lakes
subroutine lake_cal(active_lake, area_lake, Q_lake, Wr_lake, Qout, B1, B2)
  integer, intent(in) :: active_lake         ! Flag indicating if lake is active
  real*8, intent(in) :: area_lake, Qout      ! Lake area, outlet flow rate
  real*8, intent(inout) :: Q_lake, Wr_lake   ! Lake inflow, water storage
  real*8, intent(inout) :: B1, B2            ! Output variables (Q_lake, some other parameter)

  real*8 :: alp_lake                         ! Alpha parameter for lake flow calculation

  ! Process only active lakes
  if (active_lake == 1) then

    ! Determine lake type based on area and calculate alpha
    if (area_lake >= thr_area_lake) then
      alp_lake = fac_a_llake * ( (1.D0 / sqrt(area_lake)) ** fac_b_llake ) / 3600.D0
    else
      alp_lake = fac_a_slake * ( (1.D0 / sqrt(area_lake)) ** fac_b_slake ) / 3600.D0
    endif

    ! Compute lake outflow based on alpha and water storage
    Q_lake = alp_lake * Wr_lake

    ! Ensure that outflow is non-negative and does not exceed available water
    Q_lake = max(0.D0, Q_lake)
    Q_lake = min(Q_lake, Wr_lake / dt + Qout)

    ! Update water storage in lake
    Wr_lake = Wr_lake + dt * (Qout - Q_lake)
    Wr_lake = max(0.D0, Wr_lake)

    ! Assign output values
    B1 = Q_lake
    B2 = 0.D0

  endif

end subroutine lake_cal

end module lake