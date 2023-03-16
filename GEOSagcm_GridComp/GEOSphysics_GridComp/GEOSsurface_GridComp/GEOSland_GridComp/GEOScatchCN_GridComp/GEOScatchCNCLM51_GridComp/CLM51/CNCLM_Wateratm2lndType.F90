module Wateratm2lndType

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water atm2lnd variables that apply to both bulk water
  ! and water tracers.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varcon     , only : spval

  implicit none
  save
  private

  ! !PUBLIC MEMBER FUNCTIONS:
  !
  ! !PUBLIC TYPES:
  type, public :: wateratm2lnd_type

     real(r8), pointer :: forc_q_not_downscaled_grc     (:)   ! not downscaled atm specific humidity (kg/kg)
     real(r8), pointer :: forc_rain_not_downscaled_grc  (:)   ! not downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_not_downscaled_grc  (:)   ! not downscaled atm snow rate [mm/s]
     real(r8), pointer :: forc_q_downscaled_col         (:)   ! downscaled atm specific humidity (kg/kg)
     real(r8), pointer :: forc_flood_grc                (:)   ! rof flood (mm/s)
     real(r8), pointer :: forc_rain_downscaled_col      (:)   ! downscaled atm rain rate [mm/s]
     real(r8), pointer :: forc_snow_downscaled_col      (:)   ! downscaled atm snow rate [mm/s]

     real(r8), pointer :: rain_to_snow_conversion_col   (:)   ! amount of rain converted to snow via precipitation repartitioning (mm/s)
     real(r8), pointer :: snow_to_rain_conversion_col   (:)   ! amount of snow converted to rain via precipitation repartitioning (mm/s)

   contains

     procedure, public :: Init

  end type wateratm2lnd_type
!  type(wateratm2lnd_type), public, target, save :: wateratm2lnd_inst

  contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
      use nanMod      , only : nan
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    class(wateratm2lnd_type), intent(inout)      :: this
    !
    ! !LOCAL VARIABLES:
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%forc_q_not_downscaled_grc      (begg:endg))
    allocate(this%forc_rain_not_downscaled_grc   (begg:endg))
    allocate(this%forc_snow_not_downscaled_grc   (begg:endg))
    allocate(this%forc_q_downscaled_col          (begc:endc))
    allocate(this%forc_flood_grc                 (begg:endg))
    allocate(this%forc_rain_downscaled_col       (begc:endc))
    allocate(this%forc_snow_downscaled_col       (begc:endc))
    allocate(this%rain_to_snow_conversion_col    (begc:endc))
    allocate(this%snow_to_rain_conversion_col    (begc:endc))

    this%forc_rain_not_downscaled_grc(begg:endg) = spval
    this%forc_snow_not_downscaled_grc(begg:endg) = spval
    this%forc_q_downscaled_col(begc:endc) = spval
    this%forc_flood_grc(begg:endg) = spval
    this%forc_rain_downscaled_col(begc:endc) = spval
    this%forc_snow_downscaled_col(begc:endc) = spval

  end subroutine Init
end module Wateratm2lndType
