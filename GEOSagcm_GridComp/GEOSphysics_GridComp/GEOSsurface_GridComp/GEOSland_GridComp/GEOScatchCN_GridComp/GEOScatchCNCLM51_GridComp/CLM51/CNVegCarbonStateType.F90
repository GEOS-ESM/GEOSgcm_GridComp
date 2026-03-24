module CNVegCarbonStateType

#include "shr_assert.h"

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use clm_varctl       , only : iulog, use_cndv, use_crop, use_matrixcn
  use clm_varpar       , only : numpft, num_zon, num_veg, &
                                var_col, var_pft, CN_zone_weight
  use clm_varcon       , only : spval
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use pftconMod        , only : noveg, npcropmin, pftcon
  use PatchType        , only : patch
  use abortutils       , only : endrun
  use shr_log_mod      , only : errMsg => shr_log_errMsg

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: cnveg_carbonstate_type

     integer :: species  ! c12, c13, c14

     real(r8), pointer :: grainc_patch                        (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_patch                (:) ! (gC/m2) grain C storage (crop model)     
     real(r8), pointer :: grainc_xfer_patch                   (:) ! (gC/m2) grain C transfer (crop model)    
     real(r8), pointer :: matrix_cap_grainc_patch             (:) ! (gC/m2) Capacity of grain C
     real(r8), pointer :: matrix_cap_grainc_storage_patch     (:) ! (gC/m2) Capacity of grain storage C      
     real(r8), pointer :: matrix_cap_grainc_xfer_patch        (:) ! (gC/m2) Capacity of grain transfer C     
     real(r8), pointer :: leafc_patch                         (:) ! (gC/m2) leaf C            
     real(r8), pointer :: leafc_storage_patch                 (:) ! (gC/m2) leaf C storage    
     real(r8), pointer :: leafc_xfer_patch                    (:) ! (gC/m2) leaf C transfer   
     real(r8), pointer :: matrix_cap_leafc_patch              (:) ! (gC/m2) Capacity of leaf C
     real(r8), pointer :: matrix_cap_leafc_storage_patch      (:) ! (gC/m2) Capacity of leaf C storage       
     real(r8), pointer :: matrix_cap_leafc_xfer_patch         (:) ! (gC/m2) Capacity of leaf C transfer      
     real(r8), pointer :: leafc_storage_xfer_acc_patch        (:) ! (gC/m2) Accmulated leaf C transfer       
     real(r8), pointer :: storage_cdemand_patch               (:) ! (gC/m2)       C use from the C storage pool 
     real(r8), pointer :: frootc_patch                        (:) ! (gC/m2) fine root C       
     real(r8), pointer :: frootc_storage_patch                (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_patch                   (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: matrix_cap_frootc_patch             (:) ! (gC/m2) Capacity of fine root C          
     real(r8), pointer :: matrix_cap_frootc_storage_patch     (:) ! (gC/m2) Capacity of fine root C storage  
     real(r8), pointer :: matrix_cap_frootc_xfer_patch        (:) ! (gC/m2) Capacity of fine root C transfer 
     real(r8), pointer :: livestemc_patch                     (:) ! (gC/m2) live stem C       
     real(r8), pointer :: livestemc_storage_patch             (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemc_xfer_patch                (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: matrix_cap_livestemc_patch          (:) ! (gC/m2) Capacity of live stem C          
     real(r8), pointer :: matrix_cap_livestemc_storage_patch  (:) ! (gC/m2) Capacity of live stem C storage  
     real(r8), pointer :: matrix_cap_livestemc_xfer_patch     (:) ! (gC/m2) Capacity of live stem C transfer 
     real(r8), pointer :: deadstemc_patch                     (:) ! (gC/m2) dead stem C       
     real(r8), pointer :: deadstemc_storage_patch             (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemc_xfer_patch                (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: matrix_cap_deadstemc_patch          (:) ! (gC/m2) Capacity of dead stem C
     real(r8), pointer :: matrix_cap_deadstemc_storage_patch  (:) ! (gC/m2) Capacity of dead stem C storage
     real(r8), pointer :: matrix_cap_deadstemc_xfer_patch     (:) ! (gC/m2) Capacity of dead stem C transfer
     real(r8), pointer :: livecrootc_patch                    (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootc_storage_patch            (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootc_xfer_patch               (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: matrix_cap_livecrootc_patch         (:) ! (gC/m2) Capacity of live coarse root C
     real(r8), pointer :: matrix_cap_livecrootc_storage_patch (:) ! (gC/m2) Capacity of live coarse root C storage
     real(r8), pointer :: matrix_cap_livecrootc_xfer_patch    (:) ! (gC/m2) Capacity of live coarse root C transfer
     real(r8), pointer :: deadcrootc_patch                    (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootc_storage_patch            (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootc_xfer_patch               (:) ! (gC/m2) dead coarse root C transfer
     real(r8), pointer :: matrix_cap_deadcrootc_patch         (:) ! (gC/m2) Capacity of dead coarse root C
     real(r8), pointer :: matrix_cap_deadcrootc_storage_patch (:) ! (gC/m2) Capacity of dead coarse root C storage
     real(r8), pointer :: matrix_cap_deadcrootc_xfer_patch    (:) ! (gC/m2) Capacity of dead coarse root C transfer
     real(r8), pointer :: gresp_storage_patch                 (:) ! (gC/m2) growth respiration storage
     real(r8), pointer :: gresp_xfer_patch                    (:) ! (gC/m2) growth respiration transfer
     real(r8), pointer :: cpool_patch                         (:) ! (gC/m2) temporary photosynthate C pool
     real(r8), pointer :: xsmrpool_patch                      (:) ! (gC/m2) abstract C pool to meet excess MR demand
     real(r8), pointer :: xsmrpool_loss_patch                 (:) ! (gC/m2) abstract C pool to meet excess MR demand loss
     real(r8), pointer :: ctrunc_patch                        (:) ! (gC/m2) patch-level sink for C truncation
     real(r8), pointer :: woodc_patch                         (:) ! (gC/m2) wood C
     real(r8), pointer :: leafcmax_patch                      (:) ! (gC/m2) ann max leaf C
     real(r8), pointer :: totc_patch                          (:) ! (gC/m2) total patch-level carbon, including cpool
     real(r8), pointer :: rootc_col                           (:) ! (gC/m2) root carbon at column level (fire)
     real(r8), pointer :: leafc_col                           (:) ! (gC/m2) column-level leafc (fire)
     real(r8), pointer :: deadstemc_col                       (:) ! (gC/m2) column-level deadstemc (fire)
     real(r8), pointer :: fuelc_col                           (:) ! fuel load outside cropland
     real(r8), pointer :: fuelc_crop_col                      (:) ! fuel load for cropland
     real(r8), pointer :: cropseedc_deficit_patch             (:) ! (gC/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
! initial pool size of year for matrix
     real(r8), pointer :: leafc0_patch                        (:) ! (gC/m2) Initial value of leaf C for SASU
     real(r8), pointer :: leafc0_storage_patch                (:) ! (gC/m2) Initial value of leaf C storage for SASU         
     real(r8), pointer :: leafc0_xfer_patch                   (:) ! (gC/m2) Initial value of leaf C transfer for SASU        
     real(r8), pointer :: frootc0_patch                       (:) ! (gC/m2) Initial value of fine root C for SASU            
     real(r8), pointer :: frootc0_storage_patch               (:) ! (gC/m2) Initial value of fine root C storage for SASU    
     real(r8), pointer :: frootc0_xfer_patch                  (:) ! (gC/m2) Initial value of fine root C transfer for SASU   
     real(r8), pointer :: livestemc0_patch                    (:) ! (gC/m2) Initial value of live stem C for SASU            
     real(r8), pointer :: livestemc0_storage_patch            (:) ! (gC/m2) Initial value of live stem C storage for SASU    
     real(r8), pointer :: livestemc0_xfer_patch               (:) ! (gC/m2) Initial value of live stem C transfer for SASU   
     real(r8), pointer :: deadstemc0_patch                    (:) ! (gC/m2) Initial value of dead stem C for SASU            
     real(r8), pointer :: deadstemc0_storage_patch            (:) ! (gC/m2) Initial value of dead stem C storage for SASU    
     real(r8), pointer :: deadstemc0_xfer_patch               (:) ! (gC/m2) Initial value of dead stem C transfer for SASU   
     real(r8), pointer :: livecrootc0_patch                   (:) ! (gC/m2) Initial value of live coarse root C for SASU     
     real(r8), pointer :: livecrootc0_storage_patch           (:) ! (gC/m2) Initial value of live coarse root C storage for SASU
     real(r8), pointer :: livecrootc0_xfer_patch              (:) ! (gC/m2) Initial value of live coarse root C transfer for SASU
     real(r8), pointer :: deadcrootc0_patch                   (:) ! (gC/m2) Initial value of dead coarse root C for SASU     
     real(r8), pointer :: deadcrootc0_storage_patch           (:) ! (gC/m2) Initial value of dead coarse root C storage for SASU
     real(r8), pointer :: deadcrootc0_xfer_patch              (:) ! (gC/m2) Initial value of dead coarse root C transfer for SASU
     real(r8), pointer :: grainc0_patch                       (:) ! (gC/m2) Initial value of fine grain C for SASU           
     real(r8), pointer :: grainc0_storage_patch               (:) ! (gC/m2) Initial value of fine grain C storage for SASU   
     real(r8), pointer :: grainc0_xfer_patch                  (:) ! (gC/m2) Initial value of fine grain C transfer for SASU  
                                                                                                                             
     ! pools for dynamic landcover                                                                                           
     real(r8), pointer :: seedc_grc                           (:) ! (gC/m2) gridcell-level pool for seeding new PFTs via dynamic landcover  
                                                                                                                             
     ! summary (diagnostic) state variables, not involved in mass balance                                                    
     real(r8), pointer :: dispvegc_patch                      (:) ! (gC/m2) displayed veg carbon, excluding storage and cpool
     real(r8), pointer :: storvegc_patch                      (:) ! (gC/m2) stored vegetation carbon, excluding cpool        
     real(r8), pointer :: totvegc_patch                       (:) ! (gC/m2) total vegetation carbon, excluding cpool         
     real(r8), pointer :: totvegc_col                         (:) ! (gC/m2) total vegetation carbon, excluding cpool averaged to column (p2c)

     ! Total C pools                                                                                                         
     real(r8), pointer :: totc_p2c_col                        (:) ! (gC/m2) totc_patch averaged to col 
     real(r8), pointer :: totc_col                            (:) ! (gC/m2) total column carbon, incl veg and cpool
     real(r8), pointer :: totecosysc_col                      (:) ! (gC/m2) total ecosystem carbon, incl veg but excl cpool 
     real(r8), pointer :: totc_grc                            (:) ! (gC/m2) total gridcell carbon                            
                                                                                                                             
! Accumulation variables are accumulated for a whole year. They are used for matrix spinup and calculation of diagnostic variables
     real(r8), pointer :: matrix_calloc_leaf_acc_patch        (:) ! (gC/m2/year) Input C allocated to leaf during this year  
     real(r8), pointer :: matrix_calloc_leafst_acc_patch      (:) ! (gC/m2/year) Input C allocated to leaf storage during this year         
     real(r8), pointer :: matrix_calloc_froot_acc_patch       (:) ! (gC/m2/year) Input C allocated to fine root during this year
     real(r8), pointer :: matrix_calloc_frootst_acc_patch     (:) ! (gC/m2/year) Input C allocated to fine root storage during this year    
     real(r8), pointer :: matrix_calloc_livestem_acc_patch    (:) ! (gC/m2/year) Input C allocated to live stem during this year
     real(r8), pointer :: matrix_calloc_livestemst_acc_patch  (:) ! (gC/m2/year) Input C allocated to live stem storage during this year    
     real(r8), pointer :: matrix_calloc_deadstem_acc_patch    (:) ! (gC/m2/year) Input C allocated to dead stem during this year
     real(r8), pointer :: matrix_calloc_deadstemst_acc_patch  (:) ! (gC/m2/year) Input C allocated to dead stem storage during this year    
     real(r8), pointer :: matrix_calloc_livecroot_acc_patch   (:) ! (gC/m2/year) Input C allocated to live coarse root during this year     
     real(r8), pointer :: matrix_calloc_livecrootst_acc_patch (:) ! (gC/m2/year) Input C allocated to live coarse root storage during this year
     real(r8), pointer :: matrix_calloc_deadcroot_acc_patch   (:) ! (gC/m2/year) Input C allocated to dead coarse root during this year     
     real(r8), pointer :: matrix_calloc_deadcrootst_acc_patch (:) ! (gC/m2/year) Input C allocated to dead coarse root storage during this year
     real(r8), pointer :: matrix_calloc_grain_acc_patch       (:) ! (gC/m2/year) Input C allocated to grain during this year 
     real(r8), pointer :: matrix_calloc_grainst_acc_patch     (:) ! (gC/m2/year) Input C allocated to grain storage during this year        
                                                                                                                             
     real(r8), pointer :: matrix_ctransfer_leafst_to_leafxf_acc_patch           (:) ! (gC/m2/year) C transfer from leaf storage to leaf transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_leafxf_to_leaf_acc_patch             (:) ! (gC/m2/year) C transfer from leaf transfer to leaf pool during this year
     real(r8), pointer :: matrix_ctransfer_frootst_to_frootxf_acc_patch         (:) ! (gC/m2/year) C transfer from fine root storage to fine root transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_frootxf_to_froot_acc_patch           (:) ! (gC/m2/year) C transfer from fine root transfer to fine root pool during this year
     real(r8), pointer :: matrix_ctransfer_livestemst_to_livestemxf_acc_patch   (:) ! (gC/m2/year) C transfer from live stem storage to live stem transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_livestemxf_to_livestem_acc_patch     (:) ! (gC/m2/year) C transfer from live stem transfer to live stem pool during this year                                           
     real(r8), pointer :: matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch   (:) ! (gC/m2/year) C transfer from dead stem storage to dead stem transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_deadstemxf_to_deadstem_acc_patch     (:) ! (gC/m2/year) C transfer from dead stem transfer to dead stem pool during this year
     real(r8), pointer :: matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch (:) ! (gC/m2/year) C transfer from live coarse root storage to live coarse root transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   (:) ! (gC/m2/year) C transfer from live coarse root transfer to live coarse root pool during this year
     real(r8), pointer :: matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch (:) ! (gC/m2/year) C transfer from dead coarse root storage to dead coarse root transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   (:) ! (gC/m2/year) C transfer from dead coarse root transfer to dead coarse root pool during this year
     real(r8), pointer :: matrix_ctransfer_grainst_to_grainxf_acc_patch         (:) ! (gC/m2/year) C transfer from grain storage to grain transfer pool during this year
     real(r8), pointer :: matrix_ctransfer_grainxf_to_grain_acc_patch           (:) ! (gC/m2/year) C transfer from grain transfer to grain pool during this year
     real(r8), pointer :: matrix_ctransfer_livestem_to_deadstem_acc_patch       (:) ! (gC/m2/year) C transfer from live stem to dead stem pool during this year
     real(r8), pointer :: matrix_ctransfer_livecroot_to_deadcroot_acc_patch     (:) ! (gC/m2/year) C transfer from live coarse root to dead coarse root pool during this year
                                                                                                                             
     real(r8), pointer :: matrix_cturnover_leaf_acc_patch             (:) ! (gC/m2/year) C turnover from leaf
     real(r8), pointer :: matrix_cturnover_leafst_acc_patch           (:) ! (gC/m2/year) C turnover from leaf storage        
     real(r8), pointer :: matrix_cturnover_leafxf_acc_patch           (:) ! (gC/m2/year) C turnover from leaf transfer
     real(r8), pointer :: matrix_cturnover_froot_acc_patch            (:) ! (gC/m2/year) C turnover from fine root           
     real(r8), pointer :: matrix_cturnover_frootst_acc_patch          (:) ! (gC/m2/year) C turnover from fine root storage
     real(r8), pointer :: matrix_cturnover_frootxf_acc_patch          (:) ! (gC/m2/year) C turnover from fine root transfer  
     real(r8), pointer :: matrix_cturnover_livestem_acc_patch         (:) ! (gC/m2/year) C turnover from live stem
     real(r8), pointer :: matrix_cturnover_livestemst_acc_patch       (:) ! (gC/m2/year) C turnover from live stem storage   
     real(r8), pointer :: matrix_cturnover_livestemxf_acc_patch       (:) ! (gC/m2/year) C turnover from live stem transfer
     real(r8), pointer :: matrix_cturnover_deadstem_acc_patch         (:) ! (gC/m2/year) C turnover from dead stem           
     real(r8), pointer :: matrix_cturnover_deadstemst_acc_patch       (:) ! (gC/m2/year) C turnover from dead stem storage
     real(r8), pointer :: matrix_cturnover_deadstemxf_acc_patch       (:) ! (gC/m2/year) C turnover from dead stem transfer  
     real(r8), pointer :: matrix_cturnover_livecroot_acc_patch        (:) ! (gC/m2/year) C turnover from live coarse root
     real(r8), pointer :: matrix_cturnover_livecrootst_acc_patch      (:) ! (gC/m2/year) C turnover from live coarse root storage
     real(r8), pointer :: matrix_cturnover_livecrootxf_acc_patch      (:) ! (gC/m2/year) C turnover from live coarse root transfer
     real(r8), pointer :: matrix_cturnover_deadcroot_acc_patch        (:) ! (gC/m2/year) C turnover from dead coarse root    
     real(r8), pointer :: matrix_cturnover_deadcrootst_acc_patch      (:) ! (gC/m2/year) C turnover from dead coarse root storage
     real(r8), pointer :: matrix_cturnover_deadcrootxf_acc_patch      (:) ! (gC/m2/year) C turnover from dead coarse root transfer          
     real(r8), pointer :: matrix_cturnover_grain_acc_patch            (:) ! (gC/m2/year) C turnover from grain 
     real(r8), pointer :: matrix_cturnover_grainst_acc_patch          (:) ! (gC/m2/year) C turnover from grain storage       
     real(r8), pointer :: matrix_cturnover_grainxf_acc_patch          (:) ! (gC/m2/year) C turnover from grain transfer
                                                                                                                             
     real(r8), pointer :: grainc_SASUsave_patch               (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainc_storage_SASUsave_patch       (:) ! (gC/m2) grain C storage (crop model)                     
     real(r8), pointer :: leafc_SASUsave_patch                (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafc_storage_SASUsave_patch        (:) ! (gC/m2) leaf C storage                                   
     real(r8), pointer :: leafc_xfer_SASUsave_patch           (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: frootc_SASUsave_patch               (:) ! (gC/m2) fine root C                                      
     real(r8), pointer :: frootc_storage_SASUsave_patch       (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootc_xfer_SASUsave_patch          (:) ! (gC/m2) fine root C transfer                             
     real(r8), pointer :: livestemc_SASUsave_patch            (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemc_storage_SASUsave_patch    (:) ! (gC/m2) live stem C storage                              
     real(r8), pointer :: livestemc_xfer_SASUsave_patch       (:) ! (gC/m2) live stem C transfer                             
     real(r8), pointer :: deadstemc_SASUsave_patch            (:) ! (gC/m2) dead stem C                                      
     real(r8), pointer :: deadstemc_storage_SASUsave_patch    (:) ! (gC/m2) dead stem C storage                              
     real(r8), pointer :: deadstemc_xfer_SASUsave_patch       (:) ! (gC/m2) dead stem C transfer                             
     real(r8), pointer :: livecrootc_SASUsave_patch           (:) ! (gC/m2) live coarse root C                               
     real(r8), pointer :: livecrootc_storage_SASUsave_patch   (:) ! (gC/m2) live coarse root C storage                       
     real(r8), pointer :: livecrootc_xfer_SASUsave_patch      (:) ! (gC/m2) live coarse root C transfer                      
     real(r8), pointer :: deadcrootc_SASUsave_patch           (:) ! (gC/m2) dead coarse root C                               
     real(r8), pointer :: deadcrootc_storage_SASUsave_patch   (:) ! (gC/m2) dead coarse root C storage                       
     real(r8), pointer :: deadcrootc_xfer_SASUsave_patch      (:) ! (gC/m2) dead coarse root C transfer   
     logical, private  :: dribble_crophrv_xsmrpool_2atm

  contains

     procedure , public  :: Summary => Summary_carbonstate
     procedure , public  :: ZeroDWT
     procedure , public  :: Init
     procedure , private :: InitReadNML     ! Read in namelist

 end type cnveg_carbonstate_type

type(cnveg_carbonstate_type), public, target, save :: cnveg_carbonstate_inst

  real(r8), public  :: spinup_factor_deadwood = 1.0_r8        ! Spinup factor used for this simulation
  real(r8), public  :: spinup_factor_AD       = 10.0_r8       ! Spinup factor used when in Accelerated Decomposition mode

  ! !PRIVATE DATA:

  type, private :: cnvegcarbonstate_const_type
      ! !PRIVATE MEMBER DATA:
      real(r8) :: initial_vegC = 20._r8    ! Initial vegetation carbon for leafc/frootc and storage
  end type
  type(cnvegcarbonstate_const_type), private :: cnvegcstate_const    ! Constants used here

  character(len=*), parameter :: sourcefile = &
       __FILE__

contains

!----------------------------------------------
  subroutine Init(this, bounds, NLFilename, nch, ityp, fveg, cncol, cnpft)

! !DESCRIPTION:
! Initialize CTSM carbon states
! jk Apr 2021: type is allocated and initialized to NaN;
! if data arrays from restart file are passed (cncol and cnpft), the type is then initialized with these values
!
! !ARGUMENTS:
    implicit none

  ! INPUT
    type(bounds_type),                            intent(in) :: bounds
    character(len=*) ,                            intent(in) :: NLFilename                 ! Namelist filename
    integer,                                      intent(in) :: nch ! number of tiles
    integer, dimension(nch,NUM_VEG,NUM_ZON),      intent(in) :: ityp ! PFT index
    real, dimension(nch,NUM_VEG,NUM_ZON),         intent(in) :: fveg    ! PFT fraction
    real, dimension(nch,NUM_ZON,VAR_COL),         intent(in) :: cncol ! gkw: column CN restart
    real, dimension(nch,NUM_ZON,NUM_VEG,VAR_PFT), intent(in) :: cnpft ! gkw: PFT CN restart
    class(cnveg_carbonstate_type)                             :: this

    ! LOCAL
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: begg, endg
    integer  :: np, nc, nz, p, nv, n
    !--------------------------------------------------------

    begp = bounds%begp  ; endp = bounds%endp
    begg = bounds%begg  ; endg = bounds%endg
    begc = bounds%begc  ; endc = bounds%endc

    allocate(this%leafc_patch                            (begp:endp)) ; this%leafc_patch                        (:) = spval    
    allocate(this%leafc_storage_patch                    (begp:endp)) ; this%leafc_storage_patch                (:) = nan    
    allocate(this%leafc_xfer_patch                       (begp:endp)) ; this%leafc_xfer_patch                   (:) = nan    
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_leafc_patch              (begp:endp)) ; this%matrix_cap_leafc_patch             (:) = nan    
       allocate(this%matrix_cap_leafc_storage_patch      (begp:endp)) ; this%matrix_cap_leafc_storage_patch     (:) = nan    
       allocate(this%matrix_cap_leafc_xfer_patch         (begp:endp)) ; this%matrix_cap_leafc_xfer_patch        (:) = nan    
    end if                                                                                                                   
    allocate(this%leafc_storage_xfer_acc_patch           (begp:endp)) ; this%leafc_storage_xfer_acc_patch       (:) = nan    
    allocate(this%storage_cdemand_patch                  (begp:endp)) ; this%storage_cdemand_patch              (:) = nan    
    allocate(this%frootc_patch                           (begp:endp)) ; this%frootc_patch                       (:) = nan    
    allocate(this%frootc_storage_patch                   (begp:endp)) ; this%frootc_storage_patch               (:) = nan    
    allocate(this%frootc_xfer_patch                      (begp:endp)) ; this%frootc_xfer_patch                  (:) = nan    
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_frootc_patch             (begp:endp)) ; this%matrix_cap_frootc_patch            (:) = nan    
       allocate(this%matrix_cap_frootc_storage_patch     (begp:endp)) ; this%matrix_cap_frootc_storage_patch    (:) = nan    
       allocate(this%matrix_cap_frootc_xfer_patch        (begp:endp)) ; this%matrix_cap_frootc_xfer_patch       (:) = nan    
    end if                                                                                                                   
    allocate(this%livestemc_patch                        (begp:endp)) ; this%livestemc_patch                    (:) = nan    
    allocate(this%livestemc_storage_patch                (begp:endp)) ; this%livestemc_storage_patch            (:) = nan    
    allocate(this%livestemc_xfer_patch                   (begp:endp)) ; this%livestemc_xfer_patch               (:) = nan    
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_livestemc_patch          (begp:endp)) ; this%matrix_cap_livestemc_patch         (:) = nan    
       allocate(this%matrix_cap_livestemc_storage_patch  (begp:endp)) ; this%matrix_cap_livestemc_storage_patch (:) = nan    
       allocate(this%matrix_cap_livestemc_xfer_patch     (begp:endp)) ; this%matrix_cap_livestemc_xfer_patch    (:) = nan    
    end if                                                                                                                   
    allocate(this%deadstemc_patch                        (begp:endp)) ; this%deadstemc_patch                    (:) = spval
    allocate(this%deadstemc_storage_patch                (begp:endp)) ; this%deadstemc_storage_patch            (:) = nan    
    allocate(this%deadstemc_xfer_patch                   (begp:endp)) ; this%deadstemc_xfer_patch               (:) = nan    
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_deadstemc_patch          (begp:endp)) ; this%matrix_cap_deadstemc_patch         (:) = nan
       allocate(this%matrix_cap_deadstemc_storage_patch  (begp:endp)) ; this%matrix_cap_deadstemc_storage_patch (:) = nan    
       allocate(this%matrix_cap_deadstemc_xfer_patch     (begp:endp)) ; this%matrix_cap_deadstemc_xfer_patch    (:) = nan    
    end if                                                                                                                   
    allocate(this%livecrootc_patch                       (begp:endp)) ; this%livecrootc_patch                   (:) = nan    
    allocate(this%livecrootc_storage_patch               (begp:endp)) ; this%livecrootc_storage_patch           (:) = nan    
    allocate(this%livecrootc_xfer_patch                  (begp:endp)) ; this%livecrootc_xfer_patch              (:) = nan    
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_livecrootc_patch         (begp:endp)) ; this%matrix_cap_livecrootc_patch        (:) = nan    
       allocate(this%matrix_cap_livecrootc_storage_patch (begp:endp)) ; this%matrix_cap_livecrootc_storage_patch(:) = nan    
       allocate(this%matrix_cap_livecrootc_xfer_patch    (begp:endp)) ; this%matrix_cap_livecrootc_xfer_patch   (:) = nan    
    end if                                                                                                                   
    allocate(this%deadcrootc_patch                       (begp:endp)) ; this%deadcrootc_patch                   (:) = nan    
    allocate(this%deadcrootc_storage_patch               (begp:endp)) ; this%deadcrootc_storage_patch           (:) = nan    
    allocate(this%deadcrootc_xfer_patch                  (begp:endp)) ; this%deadcrootc_xfer_patch              (:) = nan    
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_deadcrootc_patch         (begp:endp)) ; this%matrix_cap_deadcrootc_patch        (:) = nan    
       allocate(this%matrix_cap_deadcrootc_storage_patch (begp:endp)) ; this%matrix_cap_deadcrootc_storage_patch(:) = nan    
       allocate(this%matrix_cap_deadcrootc_xfer_patch    (begp:endp)) ; this%matrix_cap_deadcrootc_xfer_patch   (:) = nan    
    end if                                                                                                                   
    allocate(this%gresp_storage_patch                    (begp:endp)) ; this%gresp_storage_patch                (:) = nan    
    allocate(this%gresp_xfer_patch                       (begp:endp)) ; this%gresp_xfer_patch                   (:) = nan    
    allocate(this%cpool_patch                            (begp:endp)) ; this%cpool_patch                        (:) = nan    
    allocate(this%xsmrpool_patch                         (begp:endp)) ; this%xsmrpool_patch                     (:) = nan    
    allocate(this%xsmrpool_loss_patch                    (begp:endp)) ; this%xsmrpool_loss_patch                (:) = nan    
    allocate(this%ctrunc_patch                           (begp:endp)) ; this%ctrunc_patch                       (:) = nan    
    allocate(this%dispvegc_patch                         (begp:endp)) ; this%dispvegc_patch                     (:) = nan    
    allocate(this%storvegc_patch                         (begp:endp)) ; this%storvegc_patch                     (:) = nan    
    allocate(this%leafcmax_patch                         (begp:endp)) ; this%leafcmax_patch                     (:) = nan    
    allocate(this%totc_patch                             (begp:endp))  ; this%totc_patch                        (:) = spval
    allocate(this%grainc_patch                           (begp:endp)) ; this%grainc_patch                       (:) = nan    
    allocate(this%grainc_storage_patch                   (begp:endp)) ; this%grainc_storage_patch               (:) = nan    
    allocate(this%grainc_xfer_patch                      (begp:endp)) ; this%grainc_xfer_patch                  (:) = nan  
    if(use_matrixcn)then                                                                                                     
       allocate(this%matrix_cap_grainc_patch             (begp:endp)) ; this%matrix_cap_grainc_patch            (:) = nan    
       allocate(this%matrix_cap_grainc_storage_patch     (begp:endp)) ; this%matrix_cap_grainc_storage_patch    (:) = nan    
       allocate(this%matrix_cap_grainc_xfer_patch        (begp:endp)) ; this%matrix_cap_grainc_xfer_patch       (:) = nan    
    end if                                                                                                                   
    allocate(this%woodc_patch                            (begp:endp)) ; this%woodc_patch                        (:) = nan     
!initial pool size of year for matrix                                                                                        
    if(use_matrixcn)then                                                                                                     
       allocate(this%leafc0_patch                        (begp:endp)) ; this%leafc0_patch                       (:) = nan    
       allocate(this%leafc0_storage_patch                (begp:endp)) ; this%leafc0_storage_patch               (:) = nan    
       allocate(this%leafc0_xfer_patch                   (begp:endp)) ; this%leafc0_xfer_patch                  (:) = nan    
       allocate(this%frootc0_patch                       (begp:endp)) ; this%frootc0_patch                      (:) = nan    
       allocate(this%frootc0_storage_patch               (begp:endp)) ; this%frootc0_storage_patch              (:) = nan    
       allocate(this%frootc0_xfer_patch                  (begp:endp)) ; this%frootc0_xfer_patch                 (:) = nan    
       allocate(this%livestemc0_patch                    (begp:endp)) ; this%livestemc0_patch                   (:) = nan    
       allocate(this%livestemc0_storage_patch            (begp:endp)) ; this%livestemc0_storage_patch           (:) = nan    
       allocate(this%livestemc0_xfer_patch               (begp:endp)) ; this%livestemc0_xfer_patch              (:) = nan    
       allocate(this%deadstemc0_patch                    (begp:endp)) ; this%deadstemc0_patch                   (:) = nan    
       allocate(this%deadstemc0_storage_patch            (begp:endp)) ; this%deadstemc0_storage_patch           (:) = nan    
       allocate(this%deadstemc0_xfer_patch               (begp:endp)) ; this%deadstemc0_xfer_patch              (:) = nan    
       allocate(this%livecrootc0_patch                   (begp:endp)) ; this%livecrootc0_patch                  (:) = nan    
       allocate(this%livecrootc0_storage_patch           (begp:endp)) ; this%livecrootc0_storage_patch          (:) = nan    
       allocate(this%livecrootc0_xfer_patch              (begp:endp)) ; this%livecrootc0_xfer_patch             (:) = nan    
       allocate(this%deadcrootc0_patch                   (begp:endp)) ; this%deadcrootc0_patch                  (:) = nan    
       allocate(this%deadcrootc0_storage_patch           (begp:endp)) ; this%deadcrootc0_storage_patch          (:) = nan    
       allocate(this%deadcrootc0_xfer_patch              (begp:endp)) ; this%deadcrootc0_xfer_patch             (:) = nan    
       allocate(this%grainc0_patch                       (begp:endp)) ; this%grainc0_patch                      (:) = nan    
       allocate(this%grainc0_storage_patch               (begp:endp)) ; this%grainc0_storage_patch              (:) = nan    
       allocate(this%grainc0_xfer_patch                  (begp:endp)) ; this%grainc0_xfer_patch                 (:) = nan    
                                                                                                                             
       allocate(this%leafc_SASUsave_patch                (begp:endp)) ; this%leafc_SASUsave_patch               (:) = nan    
       allocate(this%leafc_storage_SASUsave_patch        (begp:endp)) ; this%leafc_storage_SASUsave_patch       (:) = nan 
       allocate(this%leafc_xfer_SASUsave_patch           (begp:endp)) ; this%leafc_xfer_SASUsave_patch          (:) = nan    
       allocate(this%frootc_SASUsave_patch               (begp:endp)) ; this%frootc_SASUsave_patch              (:) = nan    
       allocate(this%frootc_storage_SASUsave_patch       (begp:endp)) ; this%frootc_storage_SASUsave_patch      (:) = nan    
       allocate(this%frootc_xfer_SASUsave_patch          (begp:endp)) ; this%frootc_xfer_SASUsave_patch         (:) = nan    
       allocate(this%livestemc_SASUsave_patch            (begp:endp)) ; this%livestemc_SASUsave_patch           (:) = nan    
       allocate(this%livestemc_storage_SASUsave_patch    (begp:endp)) ; this%livestemc_storage_SASUsave_patch   (:) = nan    
       allocate(this%livestemc_xfer_SASUsave_patch       (begp:endp)) ; this%livestemc_xfer_SASUsave_patch      (:) = nan    
       allocate(this%deadstemc_SASUsave_patch            (begp:endp)) ; this%deadstemc_SASUsave_patch           (:) = nan    
       allocate(this%deadstemc_storage_SASUsave_patch    (begp:endp)) ; this%deadstemc_storage_SASUsave_patch   (:) = nan    
       allocate(this%deadstemc_xfer_SASUsave_patch       (begp:endp)) ; this%deadstemc_xfer_SASUsave_patch      (:) = nan    
       allocate(this%livecrootc_SASUsave_patch           (begp:endp)) ; this%livecrootc_SASUsave_patch          (:) = nan    
       allocate(this%livecrootc_storage_SASUsave_patch   (begp:endp)) ; this%livecrootc_storage_SASUsave_patch  (:) = nan    
       allocate(this%livecrootc_xfer_SASUsave_patch      (begp:endp)) ; this%livecrootc_xfer_SASUsave_patch     (:) = nan    
       allocate(this%deadcrootc_SASUsave_patch           (begp:endp)) ; this%deadcrootc_SASUsave_patch          (:) = nan    
       allocate(this%deadcrootc_storage_SASUsave_patch   (begp:endp)) ; this%deadcrootc_storage_SASUsave_patch  (:) = nan    
       allocate(this%deadcrootc_xfer_SASUsave_patch      (begp:endp)) ; this%deadcrootc_xfer_SASUsave_patch     (:) = nan    
       allocate(this%grainc_SASUsave_patch               (begp:endp)) ; this%grainc_SASUsave_patch              (:) = nan    
       allocate(this%grainc_storage_SASUsave_patch       (begp:endp)) ; this%grainc_storage_SASUsave_patch      (:) = nan    
                                                                                                                             
       allocate(this%matrix_calloc_leaf_acc_patch        (begp:endp)); this%matrix_calloc_leaf_acc_patch        (:) = nan    
       allocate(this%matrix_calloc_leafst_acc_patch      (begp:endp)); this%matrix_calloc_leafst_acc_patch      (:) = nan    
       allocate(this%matrix_calloc_froot_acc_patch       (begp:endp)); this%matrix_calloc_froot_acc_patch       (:) = nan    
       allocate(this%matrix_calloc_frootst_acc_patch     (begp:endp)); this%matrix_calloc_frootst_acc_patch     (:) = nan    
       allocate(this%matrix_calloc_livestem_acc_patch    (begp:endp)); this%matrix_calloc_livestem_acc_patch    (:) = nan    
       allocate(this%matrix_calloc_livestemst_acc_patch  (begp:endp)); this%matrix_calloc_livestemst_acc_patch  (:) = nan    
       allocate(this%matrix_calloc_deadstem_acc_patch    (begp:endp)); this%matrix_calloc_deadstem_acc_patch    (:) = nan    
       allocate(this%matrix_calloc_deadstemst_acc_patch  (begp:endp)); this%matrix_calloc_deadstemst_acc_patch  (:) = nan    
       allocate(this%matrix_calloc_livecroot_acc_patch   (begp:endp)); this%matrix_calloc_livecroot_acc_patch   (:) = nan    
       allocate(this%matrix_calloc_livecrootst_acc_patch (begp:endp)); this%matrix_calloc_livecrootst_acc_patch (:) = nan    
       allocate(this%matrix_calloc_deadcroot_acc_patch   (begp:endp)); this%matrix_calloc_deadcroot_acc_patch   (:) = nan    
       allocate(this%matrix_calloc_deadcrootst_acc_patch (begp:endp)); this%matrix_calloc_deadcrootst_acc_patch (:) = nan    
       allocate(this%matrix_calloc_grain_acc_patch       (begp:endp)); this%matrix_calloc_grain_acc_patch       (:) = nan
       allocate(this%matrix_calloc_grainst_acc_patch     (begp:endp)); this%matrix_calloc_grainst_acc_patch     (:) = nan    
                                                                                                                             
       allocate(this%matrix_ctransfer_leafst_to_leafxf_acc_patch           (begp:endp))                                      
       this%matrix_ctransfer_leafst_to_leafxf_acc_patch                    (:) = nan                                         
       allocate(this%matrix_ctransfer_leafxf_to_leaf_acc_patch             (begp:endp))                                      
       this%matrix_ctransfer_leafxf_to_leaf_acc_patch                      (:) = nan                                         
       allocate(this%matrix_ctransfer_frootst_to_frootxf_acc_patch         (begp:endp))                                      
       this%matrix_ctransfer_frootst_to_frootxf_acc_patch                  (:) = nan                                         
       allocate(this%matrix_ctransfer_frootxf_to_froot_acc_patch           (begp:endp))                                      
       this%matrix_ctransfer_frootxf_to_froot_acc_patch                    (:) = nan                                         
       allocate(this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch   (begp:endp))                                      
       this%matrix_ctransfer_livestemst_to_livestemxf_acc_patch            (:) = nan                                         
       allocate(this%matrix_ctransfer_livestemxf_to_livestem_acc_patch     (begp:endp))                                      
       this%matrix_ctransfer_livestemxf_to_livestem_acc_patch              (:) = nan                                         
       allocate(this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch   (begp:endp))                                      
       this%matrix_ctransfer_deadstemst_to_deadstemxf_acc_patch            (:) = nan                                         
       allocate(this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch     (begp:endp))                                      
       this%matrix_ctransfer_deadstemxf_to_deadstem_acc_patch              (:) = nan                                         
       allocate(this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch (begp:endp))                                      
       this%matrix_ctransfer_livecrootst_to_livecrootxf_acc_patch          (:) = nan                                         
       allocate(this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch   (begp:endp))                                      
       this%matrix_ctransfer_livecrootxf_to_livecroot_acc_patch            (:) = nan                                         
       allocate(this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch (begp:endp))                                      
       this%matrix_ctransfer_deadcrootst_to_deadcrootxf_acc_patch          (:) = nan                                         
       allocate(this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch   (begp:endp))                                      
       this%matrix_ctransfer_deadcrootxf_to_deadcroot_acc_patch            (:) = nan                                         
       allocate(this%matrix_ctransfer_grainst_to_grainxf_acc_patch         (begp:endp))                                      
       this%matrix_ctransfer_grainst_to_grainxf_acc_patch                  (:) = nan                                         
       allocate(this%matrix_ctransfer_grainxf_to_grain_acc_patch           (begp:endp))                                      
       this%matrix_ctransfer_grainxf_to_grain_acc_patch                    (:) = nan                                         
       allocate(this%matrix_ctransfer_livestem_to_deadstem_acc_patch       (begp:endp))                                      
       this%matrix_ctransfer_livestem_to_deadstem_acc_patch                (:) = nan 
       allocate(this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch     (begp:endp))                                      
       this%matrix_ctransfer_livecroot_to_deadcroot_acc_patch              (:) = nan                                         
                                                                                                                             
       allocate(this%matrix_cturnover_leaf_acc_patch        (begp:endp)) ; this%matrix_cturnover_leaf_acc_patch        (:) = nan
       allocate(this%matrix_cturnover_leafst_acc_patch      (begp:endp)) ; this%matrix_cturnover_leafst_acc_patch      (:) = nan
       allocate(this%matrix_cturnover_leafxf_acc_patch      (begp:endp)) ; this%matrix_cturnover_leafxf_acc_patch      (:) = nan
       allocate(this%matrix_cturnover_froot_acc_patch       (begp:endp)) ; this%matrix_cturnover_froot_acc_patch       (:) = nan
       allocate(this%matrix_cturnover_frootst_acc_patch     (begp:endp)) ; this%matrix_cturnover_frootst_acc_patch     (:) = nan
       allocate(this%matrix_cturnover_frootxf_acc_patch     (begp:endp)) ; this%matrix_cturnover_frootxf_acc_patch     (:) = nan
       allocate(this%matrix_cturnover_livestem_acc_patch    (begp:endp)) ; this%matrix_cturnover_livestem_acc_patch    (:) = nan
       allocate(this%matrix_cturnover_livestemst_acc_patch  (begp:endp)) ; this%matrix_cturnover_livestemst_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_livestemxf_acc_patch  (begp:endp)) ; this%matrix_cturnover_livestemxf_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_deadstem_acc_patch    (begp:endp)) ; this%matrix_cturnover_deadstem_acc_patch    (:) = nan
       allocate(this%matrix_cturnover_deadstemst_acc_patch  (begp:endp)) ; this%matrix_cturnover_deadstemst_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_deadstemxf_acc_patch  (begp:endp)) ; this%matrix_cturnover_deadstemxf_acc_patch  (:) = nan
       allocate(this%matrix_cturnover_livecroot_acc_patch   (begp:endp)) ; this%matrix_cturnover_livecroot_acc_patch   (:) = nan
       allocate(this%matrix_cturnover_livecrootst_acc_patch (begp:endp)) ; this%matrix_cturnover_livecrootst_acc_patch (:) = nan
       allocate(this%matrix_cturnover_livecrootxf_acc_patch (begp:endp)) ; this%matrix_cturnover_livecrootxf_acc_patch (:) = nan
       allocate(this%matrix_cturnover_deadcroot_acc_patch   (begp:endp)) ; this%matrix_cturnover_deadcroot_acc_patch   (:) = nan
       allocate(this%matrix_cturnover_deadcrootst_acc_patch (begp:endp)) ; this%matrix_cturnover_deadcrootst_acc_patch (:) = nan
       allocate(this%matrix_cturnover_deadcrootxf_acc_patch (begp:endp)) ; this%matrix_cturnover_deadcrootxf_acc_patch (:) = nan
       allocate(this%matrix_cturnover_grain_acc_patch       (begp:endp)) ; this%matrix_cturnover_grain_acc_patch       (:) = nan
       allocate(this%matrix_cturnover_grainst_acc_patch     (begp:endp)) ; this%matrix_cturnover_grainst_acc_patch     (:) = nan
       allocate(this%matrix_cturnover_grainxf_acc_patch     (begp:endp)) ; this%matrix_cturnover_grainxf_acc_patch     (:) = nan
    end if                                                                                                                   
                                                                                                                             
    allocate(this%cropseedc_deficit_patch  (begp:endp)) ; this%cropseedc_deficit_patch  (:) = nan                            
    allocate(this%seedc_grc                (begg:endg)) ; this%seedc_grc                (:) = nan                            
    allocate(this%rootc_col                (begc:endc)) ; this%rootc_col                (:) = nan                            
    allocate(this%leafc_col                (begc:endc)) ; this%leafc_col                (:) = nan                            
    allocate(this%deadstemc_col            (begc:endc)) ; this%deadstemc_col            (:) = nan                            
    allocate(this%fuelc_col                (begc:endc)) ; this%fuelc_col                (:) = nan 
    allocate(this%fuelc_crop_col           (begc:endc)) ; this%fuelc_crop_col           (:) = nan                            
                                                                                                                             
    allocate(this%totvegc_patch            (begp:endp)) ; this%totvegc_patch            (:) = spval                            
    allocate(this%totvegc_col              (begc:endc)) ; this%totvegc_col              (:) = nan                            
                                                                                                                             
    allocate(this%totc_p2c_col             (begc:endc)) ; this%totc_p2c_col             (:) = nan                            
    allocate(this%totc_col                 (begc:endc)) ; this%totc_col                 (:) = spval
    allocate(this%totecosysc_col           (begc:endc)) ; this%totecosysc_col           (:) = nan                            
    allocate(this%totc_grc                 (begg:endg)) ; this%totc_grc                 (:) = nan                            
  
 ! initialize variables from restart file or set to cold start value
 n = 0
 np = 0
    do nc = 1,nch        ! catchment tile loop

       this%seedc_grc(nc) = 0.

       do nz = 1,num_zon    ! CN zone loop
          n = n + 1

          this%totvegc_col(n)  = cncol(nc,nz, 6)
          this%seedc_grc  (nc) = this%seedc_grc(nc) + cncol(nc,nz,9)*CN_zone_weight(nz) 
          this%totc_col   (n)  = cncol(nc,nz,14)  
 
          do p = 0,numpft  ! PFT index loop
             np = np + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
                
                     ! "old" variables: CNCLM45 and before
                     this%cpool_patch             (np) = cnpft(nc,nz,nv,  1)
                     this%deadcrootc_patch        (np) = cnpft(nc,nz,nv,  2)
                     this%deadcrootc_storage_patch(np) = cnpft(nc,nz,nv,  3)
                     this%deadcrootc_xfer_patch   (np) = cnpft(nc,nz,nv,  4)
                     this%deadstemc_patch         (np) = cnpft(nc,nz,nv,  5)
                     this%deadstemc_storage_patch (np) = cnpft(nc,nz,nv,  6)
                     this%deadstemc_xfer_patch    (np) = cnpft(nc,nz,nv,  7)
                     this%frootc_patch            (np) = cnpft(nc,nz,nv,  8)
                     this%frootc_storage_patch    (np) = cnpft(nc,nz,nv,  9)
                     this%frootc_xfer_patch       (np) = cnpft(nc,nz,nv,  10)
                     this%gresp_storage_patch     (np) = cnpft(nc,nz,nv,  11)
                     this%gresp_xfer_patch        (np) = cnpft(nc,nz,nv,  12)
                     this%leafc_patch             (np) = cnpft(nc,nz,nv,  13)
                     this%leafc_storage_patch     (np) = cnpft(nc,nz,nv,  14)
                     this%leafc_xfer_patch        (np) = cnpft(nc,nz,nv,  15)
                     this%livecrootc_patch        (np) = cnpft(nc,nz,nv,  16)
                     this%livecrootc_storage_patch(np) = cnpft(nc,nz,nv,  17)
                     this%livecrootc_xfer_patch   (np) = cnpft(nc,nz,nv,  18)
                     this%livestemc_patch         (np) = cnpft(nc,nz,nv,  19)
                     this%livestemc_storage_patch (np) = cnpft(nc,nz,nv,  20)
                     this%livestemc_xfer_patch    (np) = cnpft(nc,nz,nv,  21)
                     this%ctrunc_patch            (np) = cnpft(nc,nz,nv,  22)
                     this%xsmrpool_patch          (np) = cnpft(nc,nz,nv,  23)
            
                     this%totvegc_patch           (np) = &
                           this%leafc_patch(np)              + &  
                           this%leafc_storage_patch(np)      + &  
                           this%leafc_xfer_patch(np)         + &  
                           this%frootc_patch(np)             + &  
                           this%frootc_storage_patch(np)     + &  
                           this%frootc_xfer_patch(np)        + &  
                           this%livestemc_patch(np)          + &  
                           this%livestemc_storage_patch(np)  + &  
                           this%livestemc_xfer_patch(np)     + &  
                           this%deadstemc_patch(np)          + &  
                           this%deadstemc_storage_patch(np)  + &  
                           this%deadstemc_xfer_patch(np)     + &  
                           this%livecrootc_patch(np)         + &  
                           this%livecrootc_storage_patch(np) + &
                           this%livecrootc_xfer_patch(np)    + &  
                           this%deadcrootc_patch(np)         + &  
                           this%deadcrootc_storage_patch(np) + &
                           this%deadcrootc_xfer_patch(np)    + &  
                           this%gresp_storage_patch(np)      + &  
                           this%gresp_xfer_patch(np)         + &  
                           this%cpool_patch(np)
         
                 end if                                                                          
            end do !nv
       end do ! p
     end do ! nz
  end do ! nc                                                                                                         

  call this%InitReadNML  ( NLFilename )

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitReadNML(this, NLFilename)
    !
    ! !DESCRIPTION:
    ! Read the namelist for CNVegCarbonState
    !
    !USES:
    use fileutils      , only : getavu, relavu, opnfil
    use shr_nl_mod     , only : shr_nl_find_group_name
    use spmdMod        , only : masterproc, mpicom
    use shr_mpi_mod    , only : shr_mpi_bcast
    use clm_varctl     , only : iulog
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)                       :: this
    character(len=*)             , intent(in)           :: NLFilename                 ! Namelist filename
    !
    ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file

    character(len=*), parameter :: subname = 'InitReadNML'
    character(len=*), parameter :: nmlname = 'cnvegcarbonstate'   ! MUST match what is in namelist below
    !-----------------------------------------------------------------------
    real(r8) :: initial_vegC
    namelist /cnvegcarbonstate/ initial_vegC

    initial_vegC = cnvegcstate_const%initial_vegC

    if (masterproc) then
       unitn = getavu()
       write(iulog,*) 'Read in '//nmlname//'  namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, nmlname, status=ierr)
       if (ierr == 0) then
          read(unitn, nml=cnvegcarbonstate, iostat=ierr)
          if (ierr /= 0) then
             call endrun(msg="ERROR reading "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
          end if
       else
          call endrun(msg="ERROR could NOT find "//nmlname//"namelist"//errmsg(sourcefile, __LINE__))
       end if
       call relavu( unitn )
    end if

    call shr_mpi_bcast (initial_vegC            , mpicom)

    cnvegcstate_const%initial_vegC = initial_vegC

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) nmlname//' settings:'
       write(iulog,nml=cnvegcarbonstate)    ! Name here MUST be the same as in nmlname above!
       write(iulog,*) ' '
    end if

    !-----------------------------------------------------------------------
  end subroutine InitReadNML

  !-----------------------------------------------------------------------
  subroutine Summary_carbonstate(this, bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       soilbiogeochem_cwdc_col, soilbiogeochem_totlitc_col, soilbiogeochem_totsomc_col, &
       soilbiogeochem_ctrunc_col)
    !
    ! !USES:
    use subgridAveMod, only : p2c
    use clm_time_manager , only : get_nstep

    !
    ! !DESCRIPTION:
    ! Perform patch and column-level carbon summary calculations
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type)  :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_allc        ! number of columns in allc filter
    integer           , intent(in) :: filter_allc(:)  ! filter for all active columns
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    real(r8)          , intent(in) :: soilbiogeochem_cwdc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_totlitc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_totsomc_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_ctrunc_col(bounds%begc:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l       ! indices
    integer  :: fp,fc           ! lake filter indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(soilbiogeochem_cwdc_col)    == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(soilbiogeochem_totlitc_col) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(soilbiogeochem_totsomc_col) == (/bounds%endc/)), sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(soilbiogeochem_ctrunc_col)  == (/bounds%endc/)), sourcefile, __LINE__)

    ! calculate patch -level summary of carbon state

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
       this%dispvegc_patch(p) =        &
            this%leafc_patch(p)      + &
            this%frootc_patch(p)     + &
            this%livestemc_patch(p)  + &
            this%deadstemc_patch(p)  + &
            this%livecrootc_patch(p) + &
            this%deadcrootc_patch(p)

       ! stored vegetation carbon, excluding cpool (STORVEGC)
       this%storvegc_patch(p) =                &
            this%cpool_patch(p)              + &
            this%leafc_storage_patch(p)      + &
            this%frootc_storage_patch(p)     + &
            this%livestemc_storage_patch(p)  + &
            this%deadstemc_storage_patch(p)  + &
            this%livecrootc_storage_patch(p) + &
            this%deadcrootc_storage_patch(p) + &
            this%leafc_xfer_patch(p)         + &
            this%frootc_xfer_patch(p)        + &
            this%livestemc_xfer_patch(p)     + &
            this%deadstemc_xfer_patch(p)     + &
            this%livecrootc_xfer_patch(p)    + &
            this%deadcrootc_xfer_patch(p)    + &
            this%gresp_storage_patch(p)      + &
            this%gresp_xfer_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%storvegc_patch(p) =            &
               this%storvegc_patch(p)       + &
               this%grainc_storage_patch(p) + &
               this%grainc_xfer_patch(p)

          this%dispvegc_patch(p) =            &
               this%dispvegc_patch(p)       + &
               this%grainc_patch(p)
       end if

       ! total vegetation carbon, excluding cpool (TOTVEGC)
       this%totvegc_patch(p) = &
            this%dispvegc_patch(p) + &
            this%storvegc_patch(p)

       ! total patch-level carbon, including xsmrpool, ctrunc
       this%totc_patch(p) = &
            this%totvegc_patch(p) + &
            this%xsmrpool_patch(p) + &
            this%ctrunc_patch(p)

       if (use_crop) then
          this%totc_patch(p) = this%totc_patch(p) + this%cropseedc_deficit_patch(p) + &
               this%xsmrpool_loss_patch(p)
       end if

       ! (WOODC) - wood C
       this%woodc_patch(p) = &
            this%deadstemc_patch(p)    + &
            this%livestemc_patch(p)    + &
            this%deadcrootc_patch(p)   + &
            this%livecrootc_patch(p)

    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegc_patch(bounds%begp:bounds%endp), &
         this%totvegc_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totc_patch(bounds%begp:bounds%endp), &
         this%totc_p2c_col(bounds%begc:bounds%endc))

    do fc = 1,num_allc
       c = filter_allc(fc)

       ! total ecosystem carbon, including veg but excluding cpool (TOTECOSYSC)
       this%totecosysc_col(c) =    &
            soilbiogeochem_cwdc_col(c)    + &
            soilbiogeochem_totlitc_col(c) + &
            soilbiogeochem_totsomc_col(c) + &
            this%totvegc_col(c)

       ! total column carbon, including veg and cpool (TOTCOLC)
       this%totc_col(c) =  this%totc_p2c_col(c) + &
            soilbiogeochem_cwdc_col(c)      + &
            soilbiogeochem_totlitc_col(c)   + &
            soilbiogeochem_totsomc_col(c)   + &
            soilbiogeochem_ctrunc_col(c)

    end do

  end subroutine Summary_carbonstate

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_carbonstate_type) :: this
    type(bounds_type), intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegc_patch(p)   = 0._r8
       this%storvegc_patch(p)   = 0._r8
       this%totc_patch(p)       = 0._r8
    end do

  end subroutine ZeroDwt

end module CNVegCarbonStateType
