module clm_varpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:

  use clm_varpar_shared, only :    &
       VAR_COL => VAR_COL_40,      &
       VAR_PFT => VAR_PFT_40,      &
       NUM_PFT => NUM_PFT_CN_40,   &
       NUM_ZON => NUM_ZON_CN,      &
       NUM_VEG => NUM_VEG_CN_40
  
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Define number of levels

  integer, parameter :: nlevsoi     =   1     ! number of hydrologically active soil layers
  integer, parameter :: nlevgrnd    =   1     ! number of ground layers (includes lower layers that are hydrologically inactive)
  integer, parameter :: nlevsno     =   0     ! maximum number of snow layers

! Define indices used in surface file read
! maxpatch_pft     = max number of plant functional types in naturally vegetated landunit

  integer            :: maxpatch_pft

! clm_varpar_init seems to do something similar; less prone to error to move
! these three lines there? (slevis)
  
  integer, parameter :: max_pft_per_col   = NUM_PFT + 1

! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters

! !REVISION HISTORY:
! Created by Mariana Vertenstein

  ! CatchCN parameters
  ! ------------------

  real, parameter, PUBLIC, dimension(NUM_ZON) :: CN_zone_weight = (/0.10,0.45,0.45/) ! gkw: tunable; must sum to 1   
  integer, parameter, PUBLIC :: map_cat(0:NUM_PFT) = (/4,3,3,3,1,1,2,2,2,5,5,5,6,4,4,4,4,4,4,4/) ! gkw: 0 -> 6, since 8 now gone

  real, parameter, PUBLIC    :: firefac = 0.1
  ! gkw: fire tuning factor. 0: threshold WPWET; 1: threshold 1; <0: no fires

!EOP
!-----------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varpar_init
!
! !INTERFACE:
  subroutine clm_varpar_init()
!
! !DESCRIPTION:
! This subroutine initializes parameters in clm_varpar
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  maxpatch_pft   = NUM_PFT + 1 ! MAXPATCH_PFT ! gkw: this was set via compiler directive

  end subroutine clm_varpar_init

!------------------------------------------------------------------------------
end module clm_varpar
