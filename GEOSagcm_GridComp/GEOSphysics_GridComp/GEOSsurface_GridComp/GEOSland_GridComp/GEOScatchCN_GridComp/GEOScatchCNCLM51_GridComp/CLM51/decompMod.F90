module decompMod

  use shr_kind_mod     , only: r8 => shr_kind_r8
  use clm_varpar       , only: NUM_ZON, NUM_VEG, numpft

  ! !PUBLIC TYPES:
  implicit none
  save
!

  ! Define possible bounds subgrid levels
  integer, parameter, public :: BOUNDS_SUBGRID_GRIDCELL = 1
  integer, parameter, public :: BOUNDS_SUBGRID_LANDUNIT = 2
  integer, parameter, public :: BOUNDS_SUBGRID_COLUMN   = 3
  integer, parameter, public :: BOUNDS_SUBGRID_PATCH    = 4
  integer, parameter, public :: BOUNDS_SUBGRID_COHORT   = 5

  ! !PUBLIC MEMBER FUNCTIONS:

  public get_beg            ! get beg bound for a given subgrid level
  public get_end            ! get end bound for a given subgrid level

  type bounds_type
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending patch index
     integer :: begCohort, endCohort ! beginning and ending cohort indices

     integer :: level            ! whether defined on the proc or clump level
     integer :: clump_index      ! if defined on the clump level, this gives the clump index

   contains

     procedure, public :: Init

  end type bounds_type
  type(bounds_type), public, target, save :: bounds

 contains

!----------------------------------------------------
  subroutine Init(this, nch)

  ! !ARGUMENTS:                                                                                                           
    implicit none

  ! INPUT:
    integer, intent(in) :: nch         ! number of Catchment tiles
    class(bounds_type)  :: this
  !----------------------------------

  this%begg = 1 ; this%endg = nch
  this%begl = 1 ; this%endl = nch
  this%begc = 1 ; this%endc = nch*NUM_ZON
  this%begp = 1 ; this%endp = nch*NUM_ZON*(numpft+1)

  end subroutine Init


  !-----------------------------------------------------------------------
  pure function get_beg(bounds, subgrid_level) result(beg_index)
    !
    ! !DESCRIPTION:
    ! Get beginning bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! BOUNDS_SUBGRID_GRIDCELL, BOUNDS_SUBGRID_LANDUNIT, etc.
    !
    ! Returns -1 for invalid subgrid_level (does not abort in this case, in order to keep
    ! this function pure).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: beg_index  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: subgrid_level
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_beg'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case (BOUNDS_SUBGRID_GRIDCELL)
       beg_index = bounds%begg
    case (BOUNDS_SUBGRID_LANDUNIT)
       beg_index = bounds%begl
    case (BOUNDS_SUBGRID_COLUMN)
       beg_index = bounds%begc
    case (BOUNDS_SUBGRID_PATCH)
       beg_index = bounds%begp
    case (BOUNDS_SUBGRID_COHORT)
       beg_index = bounds%begCohort
    case default
       beg_index = -1
    end select

  end function get_beg

  !-----------------------------------------------------------------------
  pure function get_end(bounds, subgrid_level) result(end_index)
    !
    ! !DESCRIPTION:
    ! Get end bounds for a given subgrid level
    !
    ! subgrid_level should be one of the constants defined in this module:
    ! BOUNDS_SUBGRID_GRIDCELL, BOUNDS_SUBGRID_LANDUNIT, etc.
    !
    ! Returns -1 for invalid subgrid_level (does not abort in this case, in order to keep
    ! this function pure).
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: end_index  ! function result
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: subgrid_level
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'get_end'
    !-----------------------------------------------------------------------

    select case (subgrid_level)
    case (BOUNDS_SUBGRID_GRIDCELL)
       end_index = bounds%endg
    case (BOUNDS_SUBGRID_LANDUNIT)
       end_index = bounds%endl
    case (BOUNDS_SUBGRID_COLUMN)
       end_index = bounds%endc
    case (BOUNDS_SUBGRID_PATCH)
       end_index = bounds%endp
    case (BOUNDS_SUBGRID_COHORT)
       end_index = bounds%endCohort
    case default
       end_index = -1
    end select

  end function get_end
end module decompMod
