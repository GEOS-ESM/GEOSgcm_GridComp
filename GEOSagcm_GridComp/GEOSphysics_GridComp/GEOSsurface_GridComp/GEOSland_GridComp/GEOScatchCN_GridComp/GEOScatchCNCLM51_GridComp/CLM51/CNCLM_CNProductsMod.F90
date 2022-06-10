module CNCLM_CNProductsMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use MAPL_ExceptionHandling
  use nanMod           , only : nan
  use CNCLM_decompMod  , only : bounds_type
  use clm_varpar       , only : num_zon, var_col, cn_zone_weight

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_cn_products_type

  ! !PUBLIC TYPES:
  type, public :: cn_products_type
     private
     ! ------------------------------------------------------------------------
     ! Public instance variables
     ! ------------------------------------------------------------------------

     real(r8), pointer, public :: product_loss_grc(:)   ! (g[C or N]/m2/s) total decomposition loss from ALL product pools
     real(r8), pointer, public :: cropprod1_grc(:)  ! (g[C or N]/m2) crop product pool (grain + biofuel), 1-year lifespan
     real(r8), pointer, public :: tot_woodprod_grc(:)  ! (g[C or N]/m2) total wood product pool

     ! ------------------------------------------------------------------------
     ! Private instance variables
     ! ------------------------------------------------------------------------

  !   class(species_base_type), allocatable :: species    ! C, N, C13, C14, etc.

     ! States
     real(r8), pointer :: prod10_grc(:)       ! (g[C or N]/m2) wood product pool, 10-year lifespan
     real(r8), pointer :: prod100_grc(:)      ! (g[C or N]/m2) wood product pool, 100-year lifespan

     ! Fluxes: gains
     real(r8), pointer :: dwt_prod10_gain_grc(:)  ! (g[C or N]/m2/s) dynamic landcover addition to 10-year wood product pool
     real(r8), pointer :: dwt_prod100_gain_grc(:) ! (g[C or N]/m2/s) dynamic landcover addition to 100-year wood product pool
     real(r8), pointer :: dwt_woodprod_gain_grc(:) ! (g[C or N]/m2/s) dynamic landcover addition to wood product pools
     real(r8), pointer :: dwt_cropprod1_gain_grc(:) ! (g[C or N]/m2/s) dynamic landcover addition to 1-year crop product pool
     real(r8), pointer :: hrv_deadstem_to_prod10_patch(:)  ! (g[C or N]/m2/s) dead stem harvest to 10-year wood product pool
     real(r8), pointer :: hrv_deadstem_to_prod10_grc(:)  ! (g[C or N]/m2/s) dead stem harvest to 10-year wood product pool
     real(r8), pointer :: hrv_deadstem_to_prod100_patch(:) ! (g[C or N]/m2/s) dead stem harvest to 100-year wood product pool
     real(r8), pointer :: hrv_deadstem_to_prod100_grc(:) ! (g[C or N]/m2/s) dead stem harvest to 100-year wood product pool
     real(r8), pointer :: grain_to_cropprod1_patch(:) ! (g[C or N]/m2/s) grain to 1-year crop product pool
     real(r8), pointer :: grain_to_cropprod1_grc(:) ! (g[C or N]/m2/s) grain to 1-year crop product pool

     ! Fluxes: losses
     real(r8), pointer :: cropprod1_loss_grc(:)    ! (g[C or N]/m2/s) decomposition loss from 1-yr crop product pool
     real(r8), pointer :: prod10_loss_grc(:)       ! (g[C or N]/m2/s) decomposition loss from 10-yr wood product pool
     real(r8), pointer :: prod100_loss_grc(:)      ! (g[C or N]/m2/s) decomposition loss from 100-yr wood product pool
     real(r8), pointer :: tot_woodprod_loss_grc(:) ! (g[C or N]/m2/s) decompomposition loss from all wood product pools

  end type cn_products_type
  type(cn_products_type), public, target, save :: cn_products_inst

contains

!--------------------------------------------------------------
  subroutine init_cn_products_type(bounds, nch, cncol, species,  this)

  ! !DESCRIPTION:
  ! Initialize CTSM wood products type  needed for calling CTSM routines                                 
  ! jk Oct 2021: type is allocated and initialized to NaN; values are assigned from Catchment states before calls to CLM subroutines are made
  ! this type is only used to be able to pass Catchment states and fluxes to CLM subroutines in the format they expect         
  !                                                                                                                       
  ! !ARGUMENTS:                                                                                                           
    implicit none
    ! INPUT/OUTPUT
    type(bounds_type),                    intent(in) :: bounds
    integer,                              intent(in) :: nch       ! number of Catchment tiles
    real, dimension(nch,num_zon,var_col), intent(in) :: cncol     ! column-level restart variable array 
    character(*),                         intent(in) :: species   ! C or N
    type(cn_products_type),               intent(inout):: this

    ! LOCAL
    integer :: begp, endp
    integer :: begg, endg
    integer :: nc, nz
    !---------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begg = bounds%begg ; endg = bounds%endg

    allocate(this%cropprod1_grc(begg:endg)) ; this%cropprod1_grc(:) = nan
    allocate(this%prod10_grc(begg:endg)) ; this%prod10_grc(:) = nan
    allocate(this%prod100_grc(begg:endg)) ; this%prod100_grc(:) = nan
    allocate(this%tot_woodprod_grc(begg:endg)) ; this%tot_woodprod_grc(:) = nan

    allocate(this%dwt_prod10_gain_grc(begg:endg)) ; this%dwt_prod10_gain_grc(:) = nan
    allocate(this%dwt_prod100_gain_grc(begg:endg)) ; this%dwt_prod100_gain_grc(:) = nan
    allocate(this%dwt_woodprod_gain_grc(begg:endg)) ; this%dwt_woodprod_gain_grc(:) = nan

    allocate(this%dwt_cropprod1_gain_grc(begg:endg)) ; this%dwt_cropprod1_gain_grc(:) = nan

    allocate(this%hrv_deadstem_to_prod10_patch(begp:endp)) ; this%hrv_deadstem_to_prod10_patch(:) = nan
    allocate(this%hrv_deadstem_to_prod10_grc(begg:endg)) ; this%hrv_deadstem_to_prod10_grc(:) = nan

    allocate(this%hrv_deadstem_to_prod100_patch(begp:endp)) ; this%hrv_deadstem_to_prod100_patch(:) = nan
    allocate(this%hrv_deadstem_to_prod100_grc(begg:endg)) ; this%hrv_deadstem_to_prod100_grc(:) = nan

    allocate(this%grain_to_cropprod1_patch(begp:endp)) ; this%grain_to_cropprod1_patch(:) = nan
    allocate(this%grain_to_cropprod1_grc(begg:endg)) ; this%grain_to_cropprod1_grc(:) = nan

    allocate(this%cropprod1_loss_grc(begg:endg)) ; this%cropprod1_loss_grc(:) = nan
    allocate(this%prod10_loss_grc(begg:endg)) ; this%prod10_loss_grc(:) = nan
    allocate(this%prod100_loss_grc(begg:endg)) ; this%prod100_loss_grc(:) = nan
    allocate(this%tot_woodprod_loss_grc(begg:endg)) ; this%tot_woodprod_loss_grc(:) = nan
    allocate(this%product_loss_grc(begg:endg)) ; this%product_loss_grc(:) = nan


    ! initialize variables from restart file or set to cold start value

    do nc = 1,nch        ! catchment tile loop
    
       this%prod100_grc(nc) = 0
       this%prod10_grc(nc)  = 0

       do nz = 1,num_zon    ! CN zone loop

          if (trim(species) == 'C') then
             this%prod100_grc(nc) = this%prod100_grc(nc) + cncol(nc,nz,7)*CN_zone_weight(nz)
             this%prod10_grc(nc)  = this%prod10_grc(nc) + cncol(nc,nz,8)*CN_zone_weight(nz)
          elseif (trim(species) == 'N') then
             this%prod100_grc(nc) = this%prod100_grc(nc) + cncol(nc,nz,21)*CN_zone_weight(nz)
             this%prod10_grc(nc)  = this%prod10_grc(nc) + cncol(nc,nz,22)*CN_zone_weight(nz)
          else
             _ASSERT(.FALSE.,'unknown species')
          end if 

       end do ! nz
    end do ! nc
  end subroutine init_cn_products_type

end module CNCLM_CNProductsMod
