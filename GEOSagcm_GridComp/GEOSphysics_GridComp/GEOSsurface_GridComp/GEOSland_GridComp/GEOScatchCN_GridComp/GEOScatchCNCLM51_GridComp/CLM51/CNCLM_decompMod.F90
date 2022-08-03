module CNCLM_decompMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use clm_varpar       , only: NUM_ZON, NUM_VEG, numpft

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  public :: init_bounds

  type bounds_type
     integer :: begg, endg       ! beginning and ending gridcell index
     integer :: begl, endl       ! beginning and ending landunit index
     integer :: begc, endc       ! beginning and ending column index
     integer :: begp, endp       ! beginning and ending patch index
     integer :: begCohort, endCohort ! beginning and ending cohort indices

     integer :: level            ! whether defined on the proc or clump level
     integer :: clump_index      ! if defined on the clump level, this gives the clump index
  end type bounds_type
  type(bounds_type), public, target, save :: bounds

 contains

!----------------------------------------------------
  subroutine init_bounds(nch, this)

  ! !ARGUMENTS:                                                                                                           
    implicit none

  ! INPUT:
    integer, intent(in) :: nch         ! number of Catchment tiles
    type(bounds_type), intent(inout) :: this
  !----------------------------------

  this%begg = 1 ; this%endg = nch
  this%begl = 1 ; this%endl = nch
  this%begc = 1 ; this%endc = nch*NUM_ZON
  this%begp = 1 ; this%endp = nch*NUM_ZON*(numpft+1)

  end subroutine init_bounds
end module CNCLM_decompMod
