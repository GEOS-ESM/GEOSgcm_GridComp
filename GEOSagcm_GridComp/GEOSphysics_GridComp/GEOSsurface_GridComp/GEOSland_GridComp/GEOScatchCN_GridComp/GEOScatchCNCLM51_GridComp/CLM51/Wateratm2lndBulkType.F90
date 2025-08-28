module Wateratm2lndBulkType


  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water atm2lnd variables that just apply to bulk
  ! water. Note that this type extends the base wateratm2lnd_type, so the full
  ! wateratm2lndbulk_type contains the union of the fields defined here and the fields
  ! defined in wateratm2lnd_type.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use decompMod      , only : bounds_type
  use abortutils     , only : endrun
  use PatchType      , only : patch
  use clm_varctl     , only : iulog, use_fates, use_cn, use_cndv
  use clm_varcon     , only : spval
  use WaterAtm2lndType , only : wateratm2lnd_type

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC TYPES:
  type, extends(wateratm2lnd_type), public :: wateratm2lndbulk_type

     real(r8), pointer :: volrmch_grc                   (:)   ! rof volr main channel (m3)
     real(r8), pointer :: volr_grc                      (:)   ! rof volr total volume (m3)
     real(r8), pointer :: forc_rh_grc                   (:)   ! atmospheric relative humidity (%)
     real(r8) , pointer :: prec365_col                  (:)   ! col 365-day running mean of tot. precipitation (see comment in UpdateAccVars regarding why this is col-level despite other prec accumulators being patch-level)
     real(r8) , pointer :: prec60_patch                 (:)   ! patch 60-day running mean of tot. precipitation (mm/s)
     real(r8) , pointer :: prec10_patch                 (:)   ! patch 10-day running mean of tot. precipitation (mm/s)
     real(r8) , pointer :: rh30_patch                   (:)   ! patch 30-day running mean of relative humidity
     real(r8) , pointer :: prec24_patch                 (:)   ! patch 24-hour running mean of tot. precipitation (mm/s)
     real(r8) , pointer :: rh24_patch                   (:)   ! patch 24-hour running mean of relative humidity

    contains

     procedure, public :: InitBulk

  end type wateratm2lndbulk_type

  contains

  !------------------------------------------------------------------------
  subroutine InitBulk(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
      use nanMod      , only : nan
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    class(wateratm2lndbulk_type), intent(inout)  :: this

    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival  = 0.0_r8  ! initial value
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    call this%Init(bounds)

    allocate(this%volr_grc                      (begg:endg))        ; this%volr_grc    (:)   = ival
    allocate(this%volrmch_grc                   (begg:endg))        ; this%volrmch_grc (:)   = ival
    allocate(this%forc_rh_grc                   (begg:endg))        ; this%forc_rh_grc (:)   = ival
    allocate(this%prec365_col                   (begc:endc))        ; this%prec365_col (:)   = nan
    allocate(this%prec60_patch                  (begp:endp))        ; this%prec60_patch(:)   = spval
    allocate(this%prec10_patch                  (begp:endp))        ; this%prec10_patch(:)   = spval
    allocate(this%rh30_patch                    (begp:endp))        ; this%rh30_patch  (:)   = spval
    if (use_fates) then
       allocate(this%prec24_patch               (begp:endp))        ; this%prec24_patch(:)   = nan
       allocate(this%rh24_patch                 (begp:endp))        ; this%rh24_patch  (:)   = nan
    end if


  end subroutine InitBulk
end module Wateratm2lndBulkType
