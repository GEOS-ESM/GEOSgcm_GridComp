module WaterDiagnosticType

#include "shr_assert.h"

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Defines a derived type containing water diagnostic variables that apply to both bulk
  ! water and water tracers. Diagnostic variables are neither fundamental state variables
  ! nor fluxes between those fundamental states, but are typically derived from those
  ! states and/or fluxes.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use clm_varcon     , only : spval
  use LandunitType   , only : lun                
  use WaterStateType, only : waterstate_type
  use WaterFluxType, only : waterflux_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  type, public :: waterdiagnostic_type

     real(r8), pointer :: snowice_col            (:)   ! col average snow ice lens
     real(r8), pointer :: snowliq_col            (:)   ! col average snow liquid water

     real(r8), pointer :: h2ocan_patch           (:)   ! patch total canopy water (liq+ice) (mm H2O)
     real(r8), pointer :: total_plant_stored_h2o_col(:) ! col water that is bound in plants, including roots, sapwood, leaves, etc
                                                        ! in most cases, the vegetation scheme does not have a dynamic
                                                        ! water storage in plants, and thus 0.0 is a suitable for the trivial case.
                                                        ! When FATES is coupled in with plant hydraulics turned on, this storage
                                                        ! term is set to non-zero. (kg/m2 H2O)

     real(r8), pointer :: h2osoi_liqice_10cm_col (:)   ! col liquid water + ice lens in top 10cm of soil (kg/m2)
     real(r8), pointer :: tws_grc                (:)   ! grc total water storage (mm H2O)
     real(r8), pointer :: q_ref2m_patch          (:)   ! patch 2 m height surface specific humidity (kg/kg)
     real(r8), pointer :: qg_snow_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_soil_col            (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_h2osfc_col          (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qg_col                 (:)   ! col ground specific humidity [kg/kg]
     real(r8), pointer :: qaf_lun                (:)   ! lun urban canopy air specific humidity (kg/kg)

   contains

     procedure, public          :: Init

  end type waterdiagnostic_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
 !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    ! !ARGUMENTS:
    class(waterdiagnostic_type), intent(inout) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp, endp
    integer           :: begc, endc
    integer           :: begg, endg
    !------------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    allocate(this%h2ocan_patch     (begp:endp))
    this%h2ocan_patch(begp:endp) = spval 

    allocate(this%h2osoi_liqice_10cm_col     (begc:endc))
    this%h2osoi_liqice_10cm_col(begc:endc) = spval

    allocate(this%tws_grc     (begg:endg))
    this%tws_grc(begg:endg) = spval

    allocate(this%q_ref2m_patch     (begp:endp))
    this%q_ref2m_patch(begp:endp) = spval

    ! Snow properties - these will be vertically averaged over the snow profile

    allocate(this%snowliq_col     (begc:endc))
    this%snowliq_col(begc:endc) = spval

    allocate(this%snowice_col     (begc:endc))
    this%snowice_col(begc:endc) = spval

  end subroutine Init

end module WaterDiagnosticType
