module OzoneBaseMod

  use shr_kind_mod    , only : r8 => shr_kind_r8
  use nanMod          , only : nan
  use decompMod       , only : bounds_type

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: ozone_base_type

     ! Public data members
     ! These should be treated as read-only by other modules (except that they can be
     ! modified by extensions of the ozone_base_type)
     real(r8), pointer, public :: o3coefvsha_patch(:)  ! ozone coefficient for photosynthesis, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefvsun_patch(:)  ! ozone coefficient for photosynthesis, sunlit leaves (0 - 1)
     real(r8), pointer, public :: o3coefgsha_patch(:)  ! ozone coefficient for conductance, shaded leaves (0 - 1)
     real(r8), pointer, public :: o3coefgsun_patch(:)  ! ozone coefficient for conductance, sunlit leaves (0 - 1)

   contains

    procedure, public :: Init

  end type ozone_base_type
  type(ozone_base_type), public, target, save :: ozone_inst

contains

!------------------------------------------------
  subroutine Init(this, bounds)

  ! !DESCRIPTION:
  ! Initialize CTSM ozone base type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(ozone_base_type)        :: this

    ! LOCAL
    integer :: begp, endp
    !-----------------------
   
    begp = bounds%begp ; endp = bounds%endp

    allocate(this%o3coefvsha_patch(begp:endp))  ; this%o3coefvsha_patch(:) = nan
    allocate(this%o3coefvsun_patch(begp:endp))  ; this%o3coefvsun_patch(:) = nan
    allocate(this%o3coefgsha_patch(begp:endp))  ; this%o3coefgsha_patch(:) = nan
    allocate(this%o3coefgsun_patch(begp:endp))  ; this%o3coefgsun_patch(:) = nan

    this%o3coefvsha_patch = 1.
    this%o3coefvsun_patch = 1.
    this%o3coefgsha_patch = 1.
    this%o3coefgsun_patch = 1.

  end subroutine Init

end module OzoneBaseMod
