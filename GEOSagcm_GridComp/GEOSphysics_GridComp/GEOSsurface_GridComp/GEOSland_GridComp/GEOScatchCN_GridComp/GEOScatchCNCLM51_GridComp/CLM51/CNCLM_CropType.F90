module CropType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use clm_varcon       , only : spval
  use nanMod           , only : nan
  use CNCLM_decompMod  , only : bounds_type
  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_crop_type

  type, public :: crop_type

     ! Note that cropplant and harvdate could be 2D to facilitate rotation
     integer , pointer :: nyrs_crop_active_patch  (:)   ! number of years this crop patch has been active (0 for non-crop patches)
     logical , pointer :: croplive_patch          (:)   ! patch Flag, true if planted, not harvested
     logical , pointer :: cropplant_patch         (:)   ! patch Flag, true if planted
     integer , pointer :: harvdate_patch          (:)   ! patch harvest date
     real(r8), pointer :: fertnitro_patch         (:)   ! patch fertilizer nitrogen
     real(r8), pointer :: gddplant_patch          (:)   ! patch accum gdd past planting date for crop       (ddays)
     real(r8), pointer :: gddtsoi_patch           (:)   ! patch growing degree-days from planting (top two soil layers) (ddays)
     real(r8), pointer :: vf_patch                (:)   ! patch vernalization factor for cereal
     real(r8), pointer :: cphase_patch            (:)   ! phenology phase
     real(r8), pointer :: latbaset_patch          (:)   ! Latitude vary baset for gddplant (degree C)
     character(len=20) :: baset_mapping
     real(r8) :: baset_latvary_intercept
     real(r8) :: baset_latvary_slope

 end type crop_type
 type(crop_type), public, target, save :: crop_inst

contains

!------------------------------------------------------
  subroutine init_crop_type(bounds, this)

  ! !DESCRIPTION:
  ! Initialize CTSM crop type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    type(crop_type), intent(inout):: this

    !LOCAL
    integer, intent(in) :: begp, endp

    !---------------------------------

    begp = bounds%begp ; endp = bounds%endp

    allocate(this%nyrs_crop_active_patch(begp:endp)) ; this%nyrs_crop_active_patch(:) = 0
    allocate(this%croplive_patch (begp:endp)) ; this%croplive_patch (:) = .false.
    allocate(this%cropplant_patch(begp:endp)) ; this%cropplant_patch(:) = .false.
    allocate(this%harvdate_patch (begp:endp)) ; this%harvdate_patch (:) = huge(1)
    allocate(this%fertnitro_patch (begp:endp)) ; this%fertnitro_patch (:) = spval
    allocate(this%gddplant_patch (begp:endp)) ; this%gddplant_patch (:) = spval
    allocate(this%gddtsoi_patch  (begp:endp)) ; this%gddtsoi_patch  (:) = spval
    allocate(this%vf_patch       (begp:endp)) ; this%vf_patch       (:) = 0.0_r8
    allocate(this%cphase_patch   (begp:endp)) ; this%cphase_patch   (:) = 0.0_r8
    allocate(this%latbaset_patch (begp:endp)) ; this%latbaset_patch (:) = spval

   end subroutine init_crop_type

end module CropType
