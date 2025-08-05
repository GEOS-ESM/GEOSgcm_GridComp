module GridcellType

  use MAPL_Constants   , ONLY : MAPL_PI
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use clm_varcon       , only : ispval, max_lunit
  use clm_varpar       , only : numpft, num_zon, num_veg, var_pft

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: gridcell_type

     ! topological mapping functionality, local 1d gdc arrays
     integer , pointer :: gindex           (:) ! global index
     real(r8), pointer :: area             (:) ! total land area, gridcell (km^2)
     real(r8), pointer :: lat              (:) ! latitude (radians)
     real(r8), pointer :: lon              (:) ! longitude (radians)
     real(r8), pointer :: latdeg           (:) ! latitude (degrees)
     real(r8), pointer :: londeg           (:) ! longitude (degrees)
     logical , pointer :: active           (:) ! just needed for symmetry with other subgrid types

     integer,  pointer :: nbedrock         (:) ! index of uppermost bedrock layer

     ! Daylength
     real(r8) , pointer :: max_dayl        (:) ! maximum daylength for this grid cell (s)
     real(r8) , pointer :: dayl            (:) ! daylength (seconds)
     real(r8) , pointer :: prev_dayl       (:) ! daylength from previous timestep (seconds)

     ! indices into landunit-level arrays for landunits in this grid cell (ispval implies
     ! this landunit doesn't exist on this grid cell) [1:max_lunit, begg:endg]
     ! (note that the spatial dimension is last here, in contrast to most 2-d variables;
     ! this is for efficiency, since most loops will go over g in the outer loop, and
     ! landunit type in the inner loop)
     integer , pointer :: landunit_indices (:,:)

    contains

     procedure, public :: Init

  end type gridcell_type
  type(gridcell_type), public, target :: grc

  contains

!-----------------------------------------------
  subroutine Init(this, bounds, nch, cnpft, lats, lons)

  ! !DESCRIPTION:
! Initialize CTSM gridcell type needed for calling CTSM routines                                 
! jk Apr 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
!                                                                                                                       
! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type),                            intent(in) :: bounds
    integer,                                      intent(in) :: nch    ! number of Catchment tiles
    real, dimension(nch,num_zon,num_veg,var_pft), intent(in) :: cnpft  ! pft-level (patch-level) restart variable array
    real, dimension(nch),                         intent(in) :: lats   ! Catchment tile latitudes in radians
    real, dimension(nch),                         intent(in) :: lons   ! Catchment tile longitudes in radians
    class(gridcell_type)                                     :: this

    !LOCAL
    integer :: begg, endg
    integer :: nc
    !----------------------------

    begg = bounds%begg;  endg = bounds%endg


    ! The following is set in InitGridCells
    allocate(this%gindex    (begg:endg)) ; this%gindex    (:) = ispval
    allocate(this%area      (begg:endg)) ; this%area      (:) = nan
    allocate(this%lat       (begg:endg)) ; this%lat       (:) = nan
    allocate(this%lon       (begg:endg)) ; this%lon       (:) = nan
    allocate(this%latdeg    (begg:endg)) ; this%latdeg    (:) = nan
    allocate(this%londeg    (begg:endg)) ; this%londeg    (:) = nan
    allocate(this%active    (begg:endg)) ; this%active    (:) = .true.
    allocate(this%nbedrock  (begg:endg)) ; this%nbedrock  (:) = ispval

    ! This is initiailized in module DayLength
    allocate(this%max_dayl  (begg:endg)) ; this%max_dayl  (:) = nan
    allocate(this%dayl      (begg:endg)) ; this%dayl      (:) = nan
    allocate(this%prev_dayl (begg:endg)) ; this%prev_dayl (:) = nan

    allocate(this%landunit_indices(1:max_lunit, begg:endg)); this%landunit_indices(:,:) = ispval

    ! initialize variables from restart file or set to cold start value

    do nc = 1,nch        ! catchment tile loop

       this%lat    (nc) = lats(nc)
       this%lon    (nc) = lons(nc)
       this%latdeg (nc) = lats(nc) / MAPL_PI * 180.
       this%londeg (nc) = lons(nc) / MAPL_PI * 180.
       this%londeg (nc) = this%londeg(nc)+180 ! convert from [-180 180] to [0 360]
       this%dayl   (nc) = cnpft (nc,1,1, 28) ! variable used to be patch level and is now gridcell level; assume all patches in gridcell have same day length
       
       this%prev_dayl(nc) = this%dayl(nc) ! following previous Catchment-CN versions, daylength of previous day is initialized as daylength of current day; changed for subsequent time steps in CN_DriverMod

    end do ! nc
  end subroutine Init
end module GridcellType
