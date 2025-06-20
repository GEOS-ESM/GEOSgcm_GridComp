module WaterStateType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water state variables that apply to both bulk water
  ! and water tracers.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun
  use decompMod      , only : bounds_type
  use clm_varpar     , only : nlevgrnd, nlevsoi, nlevurb, nlevmaxurbgrnd, nlevsno
  use clm_varcon     , only : spval
  use LandunitType   , only : lun
  use ColumnType     , only : col
  use nanMod         , only : nan

  implicit none
  save
  private
!
! !PUBLIC MEMBER FUNCTIONS:

  !
  ! !PUBLIC TYPES:
  type, public :: waterstate_type

     real(r8), pointer :: h2osno_no_layers_col   (:)   ! col snow that is not resolved into layers; this is non-zero only if there is too little snow for there to be explicit snow layers (mm H2O)
     real(r8), pointer :: h2osoi_liq_col         (:,:) ! col liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_ice_col         (:,:) ! col ice lens (kg/m2) (new) (-nlevsno+1:nlevgrnd)    
     real(r8), pointer :: h2osoi_vol_col         (:,:) ! col volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osoi_vol_prs_grc     (:,:) ! grc volumetric soil water prescribed (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevgrnd)
     real(r8), pointer :: h2osfc_col             (:)   ! col surface water (mm H2O)
     real(r8), pointer :: snocan_patch           (:)   ! patch canopy snow water (mm H2O)
     real(r8), pointer :: liqcan_patch           (:)   ! patch canopy liquid water (mm H2O)

     real(r8), pointer :: wa_col                 (:)   ! col water in the unconfined aquifer (mm)

     ! For the following dynbal baseline variables: positive values are subtracted to
     ! avoid counting liquid water content of "virtual" states; negative values are added
     ! to account for missing states in the model.
     real(r8), pointer :: dynbal_baseline_liq_col(:)   ! baseline liquid water content subtracted from each column's total liquid water calculation (mm H2O)
     real(r8), pointer :: dynbal_baseline_ice_col(:)   ! baseline ice content subtracted from each column's total ice calculation (mm H2O)

     real(r8) :: aquifer_water_baseline                ! baseline value for water in the unconfined aquifer (wa_col) for this bulk / tracer (mm)

   contains

    procedure , public :: Init

  end type waterstate_type
 ! type(waterstate_type), public, target, save :: waterstate_inst

contains

!---------------------------------------------
  subroutine Init(this, bounds)

  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(waterstate_type), intent(inout)        :: this

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

    allocate( this%h2osfc_col (begc:endc))   ;  this%h2osfc_col(begc:endc) = 0._r8
    allocate( this%snocan_patch (begp:endp)) ;  this%snocan_patch(begp:endp) = 0._r8
    allocate( this%liqcan_patch (begp:endp)) ;  this%liqcan_patch(begp:endp) = 0._r8

    allocate(this%h2osoi_vol_col(begc:endc,1:nlevmaxurbgrnd))         ;  this%h2osoi_vol_col(begc:endc,         1:) = spval
    allocate(this%h2osoi_vol_prs_grc(begg:endg,1:nlevgrnd))           ;  this%h2osoi_vol_prs_grc(begg:endg,     1:) = spval
    allocate(this%h2osoi_liq_col(begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ;  this%h2osoi_liq_col(begc:endc,-nlevsno+1:) = spval
    allocate(this%h2osoi_ice_col(begc:endc,-nlevsno+1:nlevmaxurbgrnd)) ;  this%h2osoi_ice_col(bounds%begc:bounds%endc,-nlevsno+1:) = spval
   
    allocate( this%wa_col (begc:endc))                 ; this%wa_col(begc:endc) = spval
    allocate( this%h2osno_no_layers_col (begc:endc))   ; this%h2osno_no_layers_col(begc:endc) = nan
    allocate( this%dynbal_baseline_liq_col (begc:endc)); this%dynbal_baseline_liq_col(begc:endc) = nan
    allocate( this%dynbal_baseline_ice_col (begc:endc)); this%dynbal_baseline_ice_col(begc:endc) = nan

  end subroutine Init
end module WaterStateType
