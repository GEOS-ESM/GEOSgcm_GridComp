module SaturatedExcessRunoffMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Type and associated routines for calculating surface runoff due to saturated surface
  !
  ! This also includes calculations of fsat (fraction of each column that is saturated)
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use decompMod    , only : bounds_type
  use abortutils   , only : endrun
  use clm_varcon   , only : spval
  use nanMod       , only : nan

  implicit none
  save
  private

  ! !PUBLIC TYPES:

  type, public :: saturated_excess_runoff_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: fsat_col(:) ! fractional area with water table at surface

     ! Private data members
     integer :: fsat_method
     real(r8), pointer :: fcov_col(:) ! fractional impermeable area

    contains 

     procedure, public :: Init

  end type saturated_excess_runoff_type

  type, private :: params_type
     real(r8) :: fff  ! Decay factor for fractional saturated area (1/m)
  end type params_type
  type(params_type), private ::  params_inst

contains

!--------------------------------------------------------------
  subroutine Init(this, bounds)

    ! !USES:
  !                                                                                                    
  ! !ARGUMENTS:                                                         
    implicit none
    ! INPUT/OUTPUT
    type(bounds_type),                  intent(in) :: bounds
    class(saturated_excess_runoff_type)            :: this

    ! LOCAL
    integer :: begc, endc
  !-------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate(this%fsat_col(begc:endc))                 ; this%fsat_col(:)                 = nan
    allocate(this%fcov_col(begc:endc))                 ; this%fcov_col(:)                 = nan

 end subroutine Init

end module SaturatedExcessRunoffMod
