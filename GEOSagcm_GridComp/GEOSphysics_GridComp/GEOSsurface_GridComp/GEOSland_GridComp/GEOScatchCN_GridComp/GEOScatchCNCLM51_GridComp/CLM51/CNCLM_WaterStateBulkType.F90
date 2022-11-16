module WaterStateBulkType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water state variables that just apply to bulk
  ! water. Note that this type extends the base waterstate_type, so the full
  ! waterstatebulk_type contains the union of the fields defined here and the fields
  ! defined in waterstate_type.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varpar     , only : nlevmaxurbgrnd, nlevsno
  use clm_varcon     , only : spval
  use WaterStateType , only : waterstate_type
  !
  implicit none
  save
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: init_waterstatebulk_type
  !
  ! !PUBLIC TYPES:
  type, extends(waterstate_type), public :: waterstatebulk_type

     real(r8), pointer :: snow_persistence_col   (:)   ! col length of time that ground has had non-zero snow thickness (sec)
     real(r8), pointer :: int_snow_col           (:)   ! col integrated snowfall (mm H2O)

  end type waterstatebulk_type
  type(waterstatebulk_type), public, target, save :: waterstatebulk_inst

contains

!---------------------------------------------
  subroutine init_waterstatebulk_type(bounds, this)

  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    type(waterstatebulk_type), intent(inout):: this

    !LOCAL
    integer :: begp, endp
    integer :: begc, endc
    integer :: begl, endl
    integer :: begg, endg
    !--------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc
    begl = bounds%begl ; endl= bounds%endl
    begg = bounds%begg ; endg = bounds%endg

    allocate(this%snow_persistence_col   (begc:endc))                     ; this%snow_persistence_col   (:)   = nan
    allocate(this%int_snow_col           (begc:endc))                     ; this%int_snow_col           (:)   = nan

  end subroutine init_waterstatebulk_type

end module WaterStateBulkType
