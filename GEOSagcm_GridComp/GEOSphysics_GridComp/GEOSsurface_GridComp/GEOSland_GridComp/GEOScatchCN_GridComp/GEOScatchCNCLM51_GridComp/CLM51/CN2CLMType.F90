#include "MAPL_Generic.h"

module CN2CLMType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use MAPL_ExceptionHandling

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: cn2clm_type

      real(r8), pointer :: forc_hdm_cn2clm(:)         ! Human population density
      real(r8), pointer :: forc_lnfm_cn2clm(:)        ! Lightning frequency
      real(r8), pointer :: btran2_patch_cn2clm(:)     ! patch root zone soil wetness factor (0 to 1)
   contains 

    procedure, public :: Init
  
  end type cn2clm_type
  type(cn2clm_type), public, target, save :: cn2clm_inst

contains

!--------------------------------------------------------------
  subroutine Init(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM canopy state type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    ! INPUT/OUTPUT
    type(bounds_type),                                intent(in) :: bounds
    class(cn2clm_type)                                           :: this

    ! LOCAL
    integer :: begp, endp
    integer :: begg, endg
 
    !---------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begg = bounds%begg ; endg = bounds%endg


    allocate(this%forc_hdm_cn2clm     (begg:endg))           ; this%forc_hdm_cn2clm     (:)   = nan
    allocate(this%forc_lnfm_cn2clm    (begg:endg))           ; this%forc_lnfm_cn2clm    (:)   = nan
    allocate(this%btran2_patch_cn2clm (begp:endp))           ; this%btran2_patch_cn2clm (:)   = nan

  end subroutine Init

end module CN2CLMType
