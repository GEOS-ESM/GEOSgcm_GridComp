module CNProductsMod

#include "MAPL_Generic.h"
#include "shr_assert.h"

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use MAPL_ExceptionHandling
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use clm_varpar       , only : num_zon, var_col, cn_zone_weight, numpft
  use clm_time_manager , only : get_step_size_real
  use PatchType        , only : patch
  use clm_varcon       , only : spval

  ! !PUBLIC TYPES:
  implicit none
  save

!
! !PUBLIC MEMBER FUNCTIONS:

  ! !PUBLIC TYPES:
  type, public :: cn_products_type
    
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

  contains

     ! Science routines
     procedure, public  :: UpdateProducts
     procedure, private :: PartitionWoodFluxes
     procedure, private :: PartitionGrainFluxes
     procedure, private :: ComputeSummaryVars
     procedure, public  :: Init
 
  end type cn_products_type
  type(cn_products_type), public, target, save :: cn_products_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

!--------------------------------------------------------------
  subroutine Init(this, bounds, nch, cncol, species, rc)

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
    class(cn_products_type)                          :: this
    integer, optional,                    intent(out) :: rc

    ! LOCAL
    integer :: begp, endp
    integer :: begg, endg
    integer :: nc, nz, p, np
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

    allocate(this%hrv_deadstem_to_prod10_patch(begp:endp)) ; this%hrv_deadstem_to_prod10_patch(:) = spval
    allocate(this%hrv_deadstem_to_prod10_grc(begg:endg)) ; this%hrv_deadstem_to_prod10_grc(:) = nan

    allocate(this%hrv_deadstem_to_prod100_patch(begp:endp)) ; this%hrv_deadstem_to_prod100_patch(:) = spval
    allocate(this%hrv_deadstem_to_prod100_grc(begg:endg)) ; this%hrv_deadstem_to_prod100_grc(:) = nan

    allocate(this%grain_to_cropprod1_patch(begp:endp)) ; this%grain_to_cropprod1_patch(:) = nan
    allocate(this%grain_to_cropprod1_grc(begg:endg)) ; this%grain_to_cropprod1_grc(:) = nan

    allocate(this%cropprod1_loss_grc(begg:endg)) ; this%cropprod1_loss_grc(:) = nan
    allocate(this%prod10_loss_grc(begg:endg)) ; this%prod10_loss_grc(:) = nan
    allocate(this%prod100_loss_grc(begg:endg)) ; this%prod100_loss_grc(:) = nan
    allocate(this%tot_woodprod_loss_grc(begg:endg)) ; this%tot_woodprod_loss_grc(:) = nan
    allocate(this%product_loss_grc(begg:endg)) ; this%product_loss_grc(:) = nan

    this%dwt_cropprod1_gain_grc(begg:endg) = 0._r8
    this%dwt_prod10_gain_grc(begg:endg) = 0._r8
    this%dwt_prod100_gain_grc(begg:endg) = 0._r8
    this%grain_to_cropprod1_grc(begg:endg) = 0._r8
  

    ! initialize variables from restart file or set to cold start value

    np = 0
    do nc = 1,nch        ! catchment tile loop
    
       this%prod100_grc(nc) = 0._r8
       this%prod10_grc(nc)  = 0._r8
       this%cropprod1_grc(nc) = 0._r8
       this%tot_woodprod_grc(nc) = 0._r8

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

          do p = 0,numpft  ! PFT index loop
             np = np + 1
             this%hrv_deadstem_to_prod10_patch(np) = 0._r8
             this%hrv_deadstem_to_prod100_patch(np) = 0._r8
             this%grain_to_cropprod1_patch(np) = 0._r8
          end do ! p
       end do ! nz
    end do ! nc
  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine UpdateProducts(this, bounds, &
       num_soilp, filter_soilp, &
       dwt_wood_product_gain_patch, &
       wood_harvest_patch, &
       dwt_crop_product_gain_patch, &
       grain_to_cropprod_patch)
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from wood and grain product pools, and update product pool
    ! state variables for both loss and gain terms
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches

    ! dynamic landcover addition to wood product pools (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_wood_product_gain_patch( bounds%begp: )

    ! wood harvest addition to wood product pools (g/m2/s) [patch]
    real(r8), intent(in) :: wood_harvest_patch( bounds%begp: )

    ! dynamic landcover addition to crop product pools (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_crop_product_gain_patch( bounds%begp: )

    ! grain to crop product pool (g/m2/s) [patch]
    real(r8), intent(in) :: grain_to_cropprod_patch( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer  :: g        ! indices
    real(r8) :: dt       ! time step (seconds)
    real(r8) :: kprod1   ! decay constant for 1-year product pool
    real(r8) :: kprod10  ! decay constant for 10-year product pool
    real(r8) :: kprod100 ! decay constant for 100-year product pool
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(dwt_wood_product_gain_patch) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(wood_harvest_patch) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(dwt_crop_product_gain_patch) == (/bounds%endp/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(grain_to_cropprod_patch) == (/bounds%endp/)), sourcefile, __LINE__)

    call this%PartitionWoodFluxes(bounds, &
         num_soilp, filter_soilp, &
         dwt_wood_product_gain_patch(bounds%begp:bounds%endp), &
         wood_harvest_patch(bounds%begp:bounds%endp))

    call this%PartitionGrainFluxes(bounds, &
         num_soilp, filter_soilp, &
         dwt_crop_product_gain_patch(bounds%begp:bounds%endp), &
         grain_to_cropprod_patch(bounds%begp:bounds%endp))

    ! calculate losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial state over 1, 10 and 100 years,
    ! respectively, using a discrete-time fractional decay algorithm.
    kprod1  = 7.2e-8
    kprod10 = 7.2e-9
    kprod100 = 7.2e-10

    do g = bounds%begg, bounds%endg
       ! calculate fluxes out of product pools (1/sec)
       this%cropprod1_loss_grc(g) = this%cropprod1_grc(g) * kprod1
       this%prod10_loss_grc(g)    = this%prod10_grc(g)    * kprod10
       this%prod100_loss_grc(g)   = this%prod100_grc(g)   * kprod100
    end do

    ! set time steps
    dt = get_step_size_real()

    ! update product state variables
    do g = bounds%begg, bounds%endg

       ! fluxes into wood & crop product pools, from landcover change
       this%cropprod1_grc(g) = this%cropprod1_grc(g) + this%dwt_cropprod1_gain_grc(g)*dt
       this%prod10_grc(g)    = this%prod10_grc(g)    + this%dwt_prod10_gain_grc(g)*dt
       this%prod100_grc(g)   = this%prod100_grc(g)   + this%dwt_prod100_gain_grc(g)*dt

       ! fluxes into wood & crop product pools, from harvest
       this%cropprod1_grc(g) = this%cropprod1_grc(g) + this%grain_to_cropprod1_grc(g)*dt
       this%prod10_grc(g)    = this%prod10_grc(g)    + this%hrv_deadstem_to_prod10_grc(g)*dt
       this%prod100_grc(g)   = this%prod100_grc(g)   + this%hrv_deadstem_to_prod100_grc(g)*dt

       ! fluxes out of wood & crop product pools, from decomposition
       this%cropprod1_grc(g) = this%cropprod1_grc(g) - this%cropprod1_loss_grc(g)*dt
       this%prod10_grc(g)    = this%prod10_grc(g)    - this%prod10_loss_grc(g)*dt
       this%prod100_grc(g)   = this%prod100_grc(g)   - this%prod100_loss_grc(g)*dt

    end do

    call this%ComputeSummaryVars(bounds)

  end subroutine UpdateProducts

  !-----------------------------------------------------------------------
  subroutine PartitionWoodFluxes(this, bounds, &
       num_soilp, filter_soilp, &
       dwt_wood_product_gain_patch, &
       wood_harvest_patch)
    !
    ! !DESCRIPTION:
    ! Partition input wood fluxes into 10 and 100 year product pools
    !
    ! !USES:
    use pftconMod    , only : pftcon
    use subgridAveMod, only : p2g
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches

    ! dynamic landcover addition to wood product pools (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_wood_product_gain_patch( bounds%begp: )

    ! wood harvest addition to wood product pools (g/m2/s) [patch]
    real(r8), intent(in) :: wood_harvest_patch( bounds%begp: )

    !
    ! !LOCAL VARIABLES:
    integer :: fp
    integer :: p
    integer :: g
    real(r8) :: pprod10       ! PFT proportion of deadstem to 10-year product pool
    real(r8) :: pprod100      ! PFT proportion of deadstem to 100-year product pool
    real(r8) :: pprod_tot     ! PFT proportion of deadstem to any product pool
    real(r8) :: pprod10_frac  ! PFT fraction of deadstem to product pool that goes to 10-year product pool
    real(r8) :: pprod100_frac ! PFT fraction of deadstem to product pool that goes to 100-year product pool

    character(len=*), parameter :: subname = 'PartitionWoodFluxes'
    !-----------------------------------------------------------------------

    ! Partition patch-level harvest fluxes to 10 and 100-year product pools
    do fp = 1, num_soilp
       p = filter_soilp(fp)
       this%hrv_deadstem_to_prod10_patch(p)  = &
            wood_harvest_patch(p) * pftcon%pprodharv10(patch%itype(p))
       this%hrv_deadstem_to_prod100_patch(p) = &
            wood_harvest_patch(p) * (1.0_r8 - pftcon%pprodharv10(patch%itype(p)))
    end do

    ! Average harvest fluxes from patch to gridcell
    call p2g(bounds, &
         this%hrv_deadstem_to_prod10_patch(bounds%begp:bounds%endp), &
         this%hrv_deadstem_to_prod10_grc(bounds%begg:bounds%endg), &
         p2c_scale_type = 'unity', &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    call p2g(bounds, &
         this%hrv_deadstem_to_prod100_patch(bounds%begp:bounds%endp), &
         this%hrv_deadstem_to_prod100_grc(bounds%begg:bounds%endg), &
         p2c_scale_type = 'unity', &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    ! Zero the dwt gains
    do g = bounds%begg, bounds%endg
       this%dwt_prod10_gain_grc(g) = 0._r8
       this%dwt_prod100_gain_grc(g) = 0._r8
    end do


    ! Partition dynamic land cover fluxes to 10 and 100-year product pools.
    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       ! Note that pprod10 + pprod100 do NOT sum to 1: some fraction of the dwt changes
       ! was lost to other fluxes. dwt_wood_product_gain_patch gives the amount that goes
       ! to all product pools, so we need to determine the fraction of that flux that
       ! goes to each pool.
       pprod10 = pftcon%pprod10(patch%itype(p))
       pprod100 = pftcon%pprod100(patch%itype(p))
       pprod_tot = pprod10 + pprod100
       if (pprod_tot > 0) then
          pprod10_frac = pprod10 / pprod_tot
          pprod100_frac = pprod100 / pprod_tot
       else
          ! Avoid divide by 0
          pprod10_frac = 0._r8
          pprod100_frac = 0._r8
       end if

       ! Note that the patch-level fluxes are expressed per unit gridcell area. So, to go
       ! from patch-level fluxes to gridcell-level fluxes, we simply add up the various
       ! patch contributions, without having to multiply by any area weightings.
       this%dwt_prod10_gain_grc(g) = this%dwt_prod10_gain_grc(g) + &
            dwt_wood_product_gain_patch(p) * pprod10_frac
       this%dwt_prod100_gain_grc(g) = this%dwt_prod100_gain_grc(g) + &
            dwt_wood_product_gain_patch(p) * pprod100_frac
    end do

  end subroutine PartitionWoodFluxes

  !-----------------------------------------------------------------------
  subroutine PartitionGrainFluxes(this, bounds, &
       num_soilp, filter_soilp, &
       dwt_crop_product_gain_patch, &
       grain_to_cropprod_patch)
    !
    ! !DESCRIPTION:
    ! Partition input grain fluxes into crop product pools
    !
    ! For now this doesn't do much, since there is just a single (1-year) crop product
    ! pool. But this provides the capability to add different crop product pools in the
    ! future, without requiring any changes to code outside of this class. It also gives
    ! symmetry with the wood fluxes.
    !
    ! !USES:
    use subgridAveMod, only : p2g
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches

    ! dynamic landcover addition to crop product pool (g/m2/s) [patch]; although this is
    ! a patch-level flux, it is expressed per unit GRIDCELL area
    real(r8), intent(in) :: dwt_crop_product_gain_patch( bounds%begp: )

    ! grain to crop product pool(s) (g/m2/s) [patch]
    real(r8)                , intent(in)    :: grain_to_cropprod_patch( bounds%begp: )
    !
    ! !LOCAL VARIABLES:
    integer :: fp
    integer :: p
    integer :: g

    character(len=*), parameter :: subname = 'PartitionGrainFluxes'
    !-----------------------------------------------------------------------

    ! Determine gains from crop harvest

    do fp = 1, num_soilp
       p = filter_soilp(fp)

       ! For now all crop product is put in the 1-year crop product pool
       this%grain_to_cropprod1_patch(p) = grain_to_cropprod_patch(p)
    end do

    call p2g(bounds, &
         this%grain_to_cropprod1_patch(bounds%begp:bounds%endp), &
         this%grain_to_cropprod1_grc(bounds%begg:bounds%endg), &
         p2c_scale_type = 'unity', &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    ! Determine gains from dynamic landcover

    do g = bounds%begg, bounds%endg
       this%dwt_cropprod1_gain_grc(g) = 0._r8
    end do

    do p = bounds%begp, bounds%endp
       g = patch%gridcell(p)

       ! Note that the patch-level fluxes are expressed per unit gridcell area. So, to go
       ! from patch-level fluxes to gridcell-level fluxes, we simply add up the various
       ! patch contributions, without having to multiply by any area weightings.
       this%dwt_cropprod1_gain_grc(g) = this%dwt_cropprod1_gain_grc(g) + &
            dwt_crop_product_gain_patch(p)
    end do

  end subroutine PartitionGrainFluxes

  !-----------------------------------------------------------------------
  subroutine ComputeSummaryVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Compute summary variables in this object: sums across multiple product pools
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(cn_products_type) , intent(inout) :: this
    type(bounds_type)       , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: g        ! indices

    character(len=*), parameter :: subname = 'ComputeSummaryVars'
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg

       ! total wood products
       this%tot_woodprod_grc(g) = &
            this%prod10_grc(g) + &
            this%prod100_grc(g)

       ! total loss from wood products
       this%tot_woodprod_loss_grc(g) = &
            this%prod10_loss_grc(g) + &
            this%prod100_loss_grc(g)

       ! total loss from ALL products
       this%product_loss_grc(g) = &
            this%cropprod1_loss_grc(g) + &
            this%prod10_loss_grc(g) + &
            this%prod100_loss_grc(g)

       this%dwt_woodprod_gain_grc(g) = &
            this%dwt_prod100_gain_grc(g) + &
            this%dwt_prod10_gain_grc(g)
    end do

  end subroutine ComputeSummaryVars

end module CNProductsMod
