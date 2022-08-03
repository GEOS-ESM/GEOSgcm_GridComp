module CNCLM_CNDVType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing routines to drive the annual dynamic vegetation
  ! that works with CN, reset related variables,
  ! and initialize/reset time invariant variables
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_log_mod  , only : errMsg => shr_log_errMsg
  use abortutils   , only : endrun
  use decompMod    , only : bounds_type
  use clm_varctl   , only : use_cndv, iulog

  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_dgvs_type

  ! DGVM state variables structure
  type, public :: dgvs_type
     real(r8), pointer, public :: agdd_patch        (:) ! patch accumulated growing degree days above 5
     real(r8), pointer, public :: agddtw_patch      (:) ! patch accumulated growing degree days above twmax
     real(r8), pointer, public :: agdd20_patch      (:) ! patch 20-yr running mean of agdd
     real(r8), pointer, public :: tmomin20_patch    (:) ! patch 20-yr running mean of tmomin
     logical , pointer, public :: present_patch     (:) ! patch whether PATCH present in patch
     logical , pointer, public :: pftmayexist_patch (:) ! patch if .false. then exclude seasonal decid patches from tropics
     real(r8), pointer, public :: nind_patch        (:) ! patch number of individuals (#/m**2)
     real(r8), pointer, public :: lm_ind_patch      (:) ! patch individual leaf mass
     real(r8), pointer, public :: lai_ind_patch     (:) ! patch LAI per individual
     real(r8), pointer, public :: fpcinc_patch      (:) ! patch foliar projective cover increment (fraction) 
     real(r8), pointer, public :: fpcgrid_patch     (:) ! patch foliar projective cover on gridcell (fraction)
     real(r8), pointer, public :: fpcgridold_patch  (:) ! patch last yr's fpcgrid
     real(r8), pointer, public :: crownarea_patch   (:) ! patch area that each individual tree takes up (m^2)
     real(r8), pointer, public :: greffic_patch     (:)
     real(r8), pointer, public :: heatstress_patch  (:)

 end type dgvs_type
 type(dgvs_type), public, target, save :: dgvs_inst

contains

!------------------------------------------------------
  subroutine init_dgvs_type(bounds, this)

  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    type(solarabs_type), intent(inout):: this

    !LOCAL
    integer, intent(in) :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%agdd_patch        (begp:endp)) ;     this%agdd_patch        (:) = nan
    allocate(this%agddtw_patch      (begp:endp)) ;     this%agddtw_patch      (:) = nan
    allocate(this%agdd20_patch      (begp:endp)) ;     this%agdd20_patch      (:) = nan
    allocate(this%tmomin20_patch    (begp:endp)) ;     this%tmomin20_patch    (:) = nan
    allocate(this%present_patch     (begp:endp)) ;     this%present_patch     (:) = .false.
    allocate(this%pftmayexist_patch (begp:endp)) ;     this%pftmayexist_patch (:) = .true.
    allocate(this%nind_patch        (begp:endp)) ;     this%nind_patch        (:) = nan
    allocate(this%lm_ind_patch      (begp:endp)) ;     this%lm_ind_patch      (:) = nan
    allocate(this%lai_ind_patch     (begp:endp)) ;     this%lai_ind_patch     (:) = nan
    allocate(this%fpcinc_patch      (begp:endp)) ;     this%fpcinc_patch      (:) = nan
    allocate(this%fpcgrid_patch     (begp:endp)) ;     this%fpcgrid_patch     (:) = nan
    allocate(this%fpcgridold_patch  (begp:endp)) ;     this%fpcgridold_patch  (:) = nan
    allocate(this%crownarea_patch   (begp:endp)) ;     this%crownarea_patch   (:) = nan
    allocate(this%greffic_patch     (begp:endp)) ;     this%greffic_patch     (:) = nan
    allocate(this%heatstress_patch  (begp:endp)) ;     this%heatstress_patch  (:) = nan

  end subroutine init_dgvs_type

end module CNCLM_CNDVType
