module CNDVType

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

  ! !PUBLIC DATA TYPES:
  !
  ! DGVM-specific ecophysiological constants structure (patch-level)
  type, public :: dgv_ecophyscon_type
     real(r8), pointer :: crownarea_max(:)   ! patch tree maximum crown area [m2]
     real(r8), pointer :: tcmin(:)           ! patch minimum coldest monthly mean temperature [units?]
     real(r8), pointer :: tcmax(:)           ! patch maximum coldest monthly mean temperature [units?]
     real(r8), pointer :: gddmin(:)          ! patch minimum growing degree days (at or above 5 C)
     real(r8), pointer :: twmax(:)           ! patch upper limit of temperature of the warmest month [units?]
     real(r8), pointer :: reinickerp(:)      ! patch parameter in allometric equation
     real(r8), pointer :: allom1(:)          ! patch parameter in allometric
     real(r8), pointer :: allom2(:)          ! patch parameter in allometric
     real(r8), pointer :: allom3(:)          ! patch parameter in allometric
  end type dgv_ecophyscon_type
  type(dgv_ecophyscon_type), public :: dgv_ecophyscon

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

   contains
  
     procedure , public :: Init

 end type dgvs_type
 type(dgvs_type), public, target, save :: dgvs_inst

contains

!------------------------------------------------------
  subroutine Init(this, bounds)

      use nanMod      , only : nan
    use clm_varpar     , only : maxveg
    use pftconMod      , only : allom1s, allom2s, allom1, allom2, allom3, reinickerp
    use pftconMod      , only : nbrdlf_dcd_brl_shrub
    use pftconMod      , only : pftcon

  ! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(bounds_type), intent(in) :: bounds
    class(dgvs_type)              :: this

    !LOCAL
    integer :: begp, endp
    integer :: m
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


    allocate(dgv_ecophyscon%crownarea_max (0:maxveg))
    allocate(dgv_ecophyscon%tcmin         (0:maxveg))
    allocate(dgv_ecophyscon%tcmax         (0:maxveg))
    allocate(dgv_ecophyscon%gddmin        (0:maxveg))
    allocate(dgv_ecophyscon%twmax         (0:maxveg))
    allocate(dgv_ecophyscon%reinickerp    (0:maxveg))
    allocate(dgv_ecophyscon%allom1        (0:maxveg))
    allocate(dgv_ecophyscon%allom2        (0:maxveg))
    allocate(dgv_ecophyscon%allom3        (0:maxveg))

    do m = 0,maxveg
       dgv_ecophyscon%crownarea_max(m) = pftcon%pftpar20(m)
       dgv_ecophyscon%tcmin(m)         = pftcon%pftpar28(m)
       dgv_ecophyscon%tcmax(m)         = pftcon%pftpar29(m)
       dgv_ecophyscon%gddmin(m)        = pftcon%pftpar30(m)
       dgv_ecophyscon%twmax(m)         = pftcon%pftpar31(m)
       dgv_ecophyscon%reinickerp(m)    = reinickerp
       dgv_ecophyscon%allom1(m)        = allom1
       dgv_ecophyscon%allom2(m)        = allom2
       dgv_ecophyscon%allom3(m)        = allom3
       ! modification for shrubs by X.D.Z
       if (pftcon%is_shrub(m)) then
          dgv_ecophyscon%allom1(m) = allom1s
          dgv_ecophyscon%allom2(m) = allom2s
       end if
    end do
  end subroutine Init

end module CNDVType
