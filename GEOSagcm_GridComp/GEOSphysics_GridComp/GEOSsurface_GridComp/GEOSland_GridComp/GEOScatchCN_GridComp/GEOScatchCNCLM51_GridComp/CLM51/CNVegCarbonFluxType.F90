module CNVegCarbonFluxType

#include "MAPL_Generic.h"
#include "shr_assert.h"

  use shr_kind_mod     , only : r8 => shr_kind_r8
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools,&
                                nvegcpool,ncphtrans,ncgmtrans,ncfitrans,&
                                ncphouttrans,ncgmouttrans,ncfiouttrans
  use clm_varpar       , only : nlevdecomp_full, nlevgrnd,nlevdecomp
  use clm_varpar       , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                ilivestem,ilivestem_st,ilivestem_xf,&
                                ideadstem,ideadstem_st,ideadstem_xf,&
                                ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                igrain,igrain_st,igrain_xf,ioutc
  use clm_varpar       , only : numpft, num_zon, num_veg, &
                                var_col, var_pft, CN_zone_weight
  use clm_varctl       , only : use_crop, use_matrixcn, use_cndv, use_grainproduct, iulog
  use clm_varcon       , only : dzsoi_decomp
  use pftconMod        , only : npcropmin
  use clm_varcon       , only : spval
  use ColumnType       , only : col
  use PatchType        , only : patch
  use AnnualFluxDribbler  , only : annual_flux_dribbler_type, annual_flux_dribbler_gridcell
  use MAPL_ExceptionHandling
  use abortutils       , only : endrun
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
!  use SPMMod                             , only : sparse_matrix_type, diag_matrix_type, vector_type


  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: cnveg_carbonflux_type

     ! gap mortality fluxes
     real(r8), pointer :: m_leafc_to_litter_patch                   (:)     ! leaf C mortality (gC/m2/s)
     real(r8), pointer :: m_leafc_storage_to_litter_patch           (:)     ! leaf C storage mortality (gC/m2/s)
     real(r8), pointer :: m_leafc_xfer_to_litter_patch              (:)     ! leaf C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_frootc_to_litter_patch                  (:)     ! fine root C mortality (gC/m2/s)
     real(r8), pointer :: m_frootc_storage_to_litter_patch          (:)     ! fine root C storage mortality (gC/m2/s)
     real(r8), pointer :: m_frootc_xfer_to_litter_patch             (:)     ! fine root C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_livestemc_to_litter_patch               (:)     ! live stem C mortality (gC/m2/s)
     real(r8), pointer :: m_livestemc_storage_to_litter_patch       (:)     ! live stem C storage mortality (gC/m2/s)
     real(r8), pointer :: m_livestemc_xfer_to_litter_patch          (:)     ! live stem C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_deadstemc_to_litter_patch               (:)     ! dead stem C mortality (gC/m2/s)
     real(r8), pointer :: m_deadstemc_storage_to_litter_patch       (:)     ! dead stem C storage mortality (gC/m2/s)
     real(r8), pointer :: m_deadstemc_xfer_to_litter_patch          (:)     ! dead stem C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_livecrootc_to_litter_patch              (:)     ! live coarse root C mortality (gC/m2/s)
     real(r8), pointer :: m_livecrootc_storage_to_litter_patch      (:)     ! live coarse root C storage mortality (gC/m2/s)
     real(r8), pointer :: m_livecrootc_xfer_to_litter_patch         (:)     ! live coarse root C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_deadcrootc_to_litter_patch              (:)     ! dead coarse root C mortality (gC/m2/s)
     real(r8), pointer :: m_deadcrootc_storage_to_litter_patch      (:)     ! dead coarse root C storage mortality (gC/m2/s)
     real(r8), pointer :: m_deadcrootc_xfer_to_litter_patch         (:)     ! dead coarse root C transfer mortality (gC/m2/s)
     real(r8), pointer :: m_gresp_storage_to_litter_patch           (:)     ! growth respiration storage mortality (gC/m2/s)
     real(r8), pointer :: m_gresp_xfer_to_litter_patch              (:)     ! growth respiration transfer mortality (gC/m2/s)

     ! harvest mortality fluxes
     real(r8), pointer :: hrv_leafc_to_litter_patch                 (:)     ! leaf C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_leafc_storage_to_litter_patch         (:)     ! leaf C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_leafc_xfer_to_litter_patch            (:)     ! leaf C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_frootc_to_litter_patch                (:)     ! fine root C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_frootc_storage_to_litter_patch        (:)     ! fine root C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_frootc_xfer_to_litter_patch           (:)     ! fine root C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_to_litter_patch             (:)     ! live stem C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_storage_to_litter_patch     (:)     ! live stem C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livestemc_xfer_to_litter_patch        (:)     ! live stem C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadstemc_storage_to_litter_patch     (:)     ! dead stem C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadstemc_xfer_to_litter_patch        (:)     ! dead stem C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livecrootc_to_litter_patch            (:)     ! live coarse root C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livecrootc_storage_to_litter_patch    (:)     ! live coarse root C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_livecrootc_xfer_to_litter_patch       (:)     ! live coarse root C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadcrootc_to_litter_patch            (:)     ! dead coarse root C harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadcrootc_storage_to_litter_patch    (:)     ! dead coarse root C storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_deadcrootc_xfer_to_litter_patch       (:)     ! dead coarse root C transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_gresp_storage_to_litter_patch         (:)     ! growth respiration storage harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_gresp_xfer_to_litter_patch            (:)     ! growth respiration transfer harvest mortality (gC/m2/s)
     real(r8), pointer :: hrv_xsmrpool_to_atm_patch                 (:)     ! excess MR pool harvest mortality (gC/m2/s)

     ! fire fluxes 
     real(r8), pointer :: m_leafc_to_fire_patch                     (:)     ! (gC/m2/s) fire C emissions from leafc 
     real(r8), pointer :: m_leafc_storage_to_fire_patch             (:)     ! (gC/m2/s) fire C emissions from leafc_storage             
     real(r8), pointer :: m_leafc_xfer_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from leafc_xfer
     real(r8), pointer :: m_livestemc_to_fire_patch                 (:)     ! (gC/m2/s) fire C emissions from livestemc
     real(r8), pointer :: m_livestemc_storage_to_fire_patch         (:)     ! (gC/m2/s) fire C emissions from livestemc_storage       
     real(r8), pointer :: m_livestemc_xfer_to_fire_patch            (:)     ! (gC/m2/s) fire C emissions from livestemc_xfer
     real(r8), pointer :: m_deadstemc_to_fire_patch                 (:)     ! (gC/m2/s) fire C emissions from deadstemc_xfer
     real(r8), pointer :: m_deadstemc_storage_to_fire_patch         (:)     ! (gC/m2/s) fire C emissions from deadstemc_storage         
     real(r8), pointer :: m_deadstemc_xfer_to_fire_patch            (:)     ! (gC/m2/s) fire C emissions from deadstemc_xfer
     real(r8), pointer :: m_frootc_to_fire_patch                    (:)     ! (gC/m2/s) fire C emissions from frootc
     real(r8), pointer :: m_frootc_storage_to_fire_patch            (:)     ! (gC/m2/s) fire C emissions from frootc_storage
     real(r8), pointer :: m_frootc_xfer_to_fire_patch               (:)     ! (gC/m2/s) fire C emissions from frootc_xfer
     real(r8), pointer :: m_livecrootc_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from livecrootc
     real(r8), pointer :: m_livecrootc_storage_to_fire_patch        (:)     ! (gC/m2/s) fire C emissions from livecrootc_storage     
     real(r8), pointer :: m_livecrootc_xfer_to_fire_patch           (:)     ! (gC/m2/s) fire C emissions from livecrootc_xfer
     real(r8), pointer :: m_deadcrootc_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from deadcrootc
     real(r8), pointer :: m_deadcrootc_storage_to_fire_patch        (:)     ! (gC/m2/s) fire C emissions from deadcrootc_storage 
     real(r8), pointer :: m_deadcrootc_xfer_to_fire_patch           (:)     ! (gC/m2/s) fire C emissions from deadcrootc_xfer
     real(r8), pointer :: m_gresp_storage_to_fire_patch             (:)     ! (gC/m2/s) fire C emissions from gresp_storage 
     real(r8), pointer :: m_gresp_xfer_to_fire_patch                (:)     ! (gC/m2/s) fire C emissions from gresp_xfer
     real(r8), pointer :: m_leafc_to_litter_fire_patch              (:)     ! (gC/m2/s) from leafc to litter c due to fire
     real(r8), pointer :: m_leafc_storage_to_litter_fire_patch      (:)     ! (gC/m2/s) from leafc_storage to litter C  due to fire               
     real(r8), pointer :: m_leafc_xfer_to_litter_fire_patch         (:)     ! (gC/m2/s) from leafc_xfer to litter C  due to fire               
     real(r8), pointer :: m_livestemc_to_litter_fire_patch          (:)     ! (gC/m2/s) from livestemc to litter C  due to fire               
     real(r8), pointer :: m_livestemc_storage_to_litter_fire_patch  (:)     ! (gC/m2/s) from livestemc_storage to litter C due to fire      
     real(r8), pointer :: m_livestemc_xfer_to_litter_fire_patch     (:)     ! (gC/m2/s) from livestemc_xfer to litter C due to fire      
     real(r8), pointer :: m_livestemc_to_deadstemc_fire_patch       (:)     ! (gC/m2/s) from livestemc to deadstemc due to fire       
     real(r8), pointer :: m_deadstemc_to_litter_fire_patch          (:)     ! (gC/m2/s) from deadstemc to litter C due to fire      
     real(r8), pointer :: m_deadstemc_storage_to_litter_fire_patch  (:)     ! (gC/m2/s) from deadstemc_storage to litter C due to fire               
     real(r8), pointer :: m_deadstemc_xfer_to_litter_fire_patch     (:)     ! (gC/m2/s) from deadstemc_xfer to litter C due to fire               
     real(r8), pointer :: m_frootc_to_litter_fire_patch             (:)     ! (gC/m2/s) from frootc to litter C due to fire               
     real(r8), pointer :: m_frootc_storage_to_litter_fire_patch     (:)     ! (gC/m2/s) from frootc_storage to litter C due to fire               
     real(r8), pointer :: m_frootc_xfer_to_litter_fire_patch        (:)     ! (gC/m2/s) from frootc_xfer to litter C due to fire               
     real(r8), pointer :: m_livecrootc_to_litter_fire_patch         (:)     ! (gC/m2/s) from livecrootc to litter C due to fire                     
     real(r8), pointer :: m_livecrootc_storage_to_litter_fire_patch (:)     ! (gC/m2/s) from livecrootc_storage to litter C due to fire                     
     real(r8), pointer :: m_livecrootc_xfer_to_litter_fire_patch    (:)     ! (gC/m2/s) from livecrootc_xfer to litter C due to fire                     
     real(r8), pointer :: m_livecrootc_to_deadcrootc_fire_patch     (:)     ! (gC/m2/s) from livecrootc to deadstemc due to fire        
     real(r8), pointer :: m_deadcrootc_to_litter_fire_patch         (:)     ! (gC/m2/s) from deadcrootc to litter C due to fire                       
     real(r8), pointer :: m_deadcrootc_storage_to_litter_fire_patch (:)     ! (gC/m2/s) from deadcrootc_storage to litter C due to fire                       
     real(r8), pointer :: m_deadcrootc_xfer_to_litter_fire_patch    (:)     ! (gC/m2/s) from deadcrootc_xfer to litter C due to fire                       
     real(r8), pointer :: m_gresp_storage_to_litter_fire_patch      (:)     ! (gC/m2/s) from gresp_storage to litter C due to fire                       
     real(r8), pointer :: m_gresp_xfer_to_litter_fire_patch         (:)     ! (gC/m2/s) from gresp_xfer to litter C due to fire                       

     ! phenology fluxes from transfer pools                     
     real(r8), pointer :: grainc_xfer_to_grainc_patch               (:)     ! grain C growth from storage for prognostic crop(gC/m2/s)
     real(r8), pointer :: leafc_xfer_to_leafc_patch                 (:)     ! leaf C growth from storage (gC/m2/s)
     real(r8), pointer :: frootc_xfer_to_frootc_patch               (:)     ! fine root C growth from storage (gC/m2/s)
     real(r8), pointer :: livestemc_xfer_to_livestemc_patch         (:)     ! live stem C growth from storage (gC/m2/s)
     real(r8), pointer :: deadstemc_xfer_to_deadstemc_patch         (:)     ! dead stem C growth from storage (gC/m2/s)
     real(r8), pointer :: livecrootc_xfer_to_livecrootc_patch       (:)     ! live coarse root C growth from storage (gC/m2/s)
     real(r8), pointer :: deadcrootc_xfer_to_deadcrootc_patch       (:)     ! dead coarse root C growth from storage (gC/m2/s)

     ! leaf and fine root litterfall fluxes                          
     real(r8), pointer :: leafc_to_litter_patch                     (:)     ! leaf C litterfall (gC/m2/s)
     real(r8), pointer :: leafc_to_litter_fun_patch                 (:)     ! leaf C litterfall used by FUN (gC/m2/s)
     real(r8), pointer :: frootc_to_litter_patch                    (:)     ! fine root C litterfall (gC/m2/s)
     real(r8), pointer :: livestemc_to_litter_patch                 (:)     ! live stem C litterfall (gC/m2/s)
     real(r8), pointer :: grainc_to_food_patch                      (:)     ! grain C to food for prognostic crop(gC/m2/s)

     real(r8), pointer :: leafc_to_biofuelc_patch                   (:)     ! leaf C to biofuel C (gC/m2/s)
     real(r8), pointer :: livestemc_to_biofuelc_patch               (:)     ! livestem C to biofuel C (gC/m2/s)
     real(r8), pointer :: grainc_to_seed_patch                      (:)     ! grain C to seed for prognostic crop(gC/m2/s)

     ! maintenance respiration fluxes     
     real(r8), pointer :: cpool_to_resp_patch                       (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_leafc_resp_patch                 (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_leafc_storage_resp_patch         (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_frootc_resp_patch                (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_frootc_storage_resp_patch        (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_livecrootc_resp_patch            (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_livecrootc_storage_resp_patch    (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_livestemc_resp_patch             (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: cpool_to_livestemc_storage_resp_patch     (:)     ! CNflex excess C maintenance respiration (gC/m2/s)
     real(r8), pointer :: leaf_mr_patch                             (:)     ! leaf maintenance respiration (gC/m2/s)
     real(r8), pointer :: froot_mr_patch                            (:)     ! fine root maintenance respiration (gC/m2/s)
     real(r8), pointer :: livestem_mr_patch                         (:)     ! live stem maintenance respiration (gC/m2/s)
     real(r8), pointer :: livecroot_mr_patch                        (:)     ! live coarse root maintenance respiration (gC/m2/s)
     real(r8), pointer :: grain_mr_patch                            (:)     ! crop grain or organs maint. respiration (gC/m2/s)
     real(r8), pointer :: leaf_curmr_patch                          (:)     ! leaf maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: froot_curmr_patch                         (:)     ! fine root maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: livestem_curmr_patch                      (:)     ! live stem maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: livecroot_curmr_patch                     (:)     ! live coarse root maintenance respiration from current GPP (gC/m2/s)
     real(r8), pointer :: grain_curmr_patch                         (:)     ! crop grain or organs maint. respiration from current GPP (gC/m2/s)
     real(r8), pointer :: leaf_xsmr_patch                           (:)     ! leaf maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: froot_xsmr_patch                          (:)     ! fine root maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: livestem_xsmr_patch                       (:)     ! live stem maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: livecroot_xsmr_patch                      (:)     ! live coarse root maintenance respiration from storage (gC/m2/s)
     real(r8), pointer :: grain_xsmr_patch                          (:)     ! crop grain or organs maint. respiration from storage (gC/m2/s)

     ! photosynthesis fluxes                                   
     real(r8), pointer :: psnsun_to_cpool_patch                     (:)     ! C fixation from sunlit canopy (gC/m2/s)
     real(r8), pointer :: psnshade_to_cpool_patch                   (:)     ! C fixation from shaded canopy (gC/m2/s)

     ! allocation fluxes, from current GPP                     
     real(r8), pointer :: cpool_to_xsmrpool_patch                   (:)     ! allocation to maintenance respiration storage pool (gC/m2/s)
     real(r8), pointer :: cpool_to_grainc_patch                     (:)     ! allocation to grain C for prognostic crop(gC/m2/s)
     real(r8), pointer :: cpool_to_grainc_storage_patch             (:)     ! allocation to grain C storage for prognostic crop(gC/m2/s)
     real(r8), pointer :: cpool_to_leafc_patch                      (:)     ! allocation to leaf C (gC/m2/s)
     real(r8), pointer :: cpool_to_leafc_storage_patch              (:)     ! allocation to leaf C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_frootc_patch                     (:)     ! allocation to fine root C (gC/m2/s)
     real(r8), pointer :: cpool_to_frootc_storage_patch             (:)     ! allocation to fine root C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_livestemc_patch                  (:)     ! allocation to live stem C (gC/m2/s)
     real(r8), pointer :: cpool_to_livestemc_storage_patch          (:)     ! allocation to live stem C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_deadstemc_patch                  (:)     ! allocation to dead stem C (gC/m2/s)
     real(r8), pointer :: cpool_to_deadstemc_storage_patch          (:)     ! allocation to dead stem C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_livecrootc_patch                 (:)     ! allocation to live coarse root C (gC/m2/s)
     real(r8), pointer :: cpool_to_livecrootc_storage_patch         (:)     ! allocation to live coarse root C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_deadcrootc_patch                 (:)     ! allocation to dead coarse root C (gC/m2/s)
     real(r8), pointer :: cpool_to_deadcrootc_storage_patch         (:)     ! allocation to dead coarse root C storage (gC/m2/s)
     real(r8), pointer :: cpool_to_gresp_storage_patch              (:)     ! allocation to growth respiration storage (gC/m2/s)


     ! growth respiration fluxes                               
     real(r8), pointer :: xsmrpool_to_atm_patch                     (:)     ! excess MR pool harvest mortality (gC/m2/s)
     real(r8), pointer :: xsmrpool_to_atm_col                       (:)     ! excess MR pool harvest mortality (gC/m2/s) (p2c)
     real(r8), pointer :: xsmrpool_to_atm_grc                       (:)     ! excess MR pool harvest mortality (gC/m2/s) (p2g)
     real(r8), pointer :: cpool_leaf_gr_patch                       (:)     ! leaf growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_leaf_storage_gr_patch               (:)     ! leaf growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_leaf_gr_patch                    (:)     ! leaf growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_froot_gr_patch                      (:)     ! fine root growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_froot_storage_gr_patch              (:)     ! fine root  growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_froot_gr_patch                   (:)     ! fine root  growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_livestem_gr_patch                   (:)     ! live stem growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_livestem_storage_gr_patch           (:)     ! live stem growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_livestem_gr_patch                (:)     ! live stem growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_deadstem_gr_patch                   (:)     ! dead stem growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_deadstem_storage_gr_patch           (:)     ! dead stem growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_deadstem_gr_patch                (:)     ! dead stem growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_livecroot_gr_patch                  (:)     ! live coarse root growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_livecroot_storage_gr_patch          (:)     ! live coarse root growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_livecroot_gr_patch               (:)     ! live coarse root growth respiration from storage (gC/m2/s)
     real(r8), pointer :: cpool_deadcroot_gr_patch                  (:)     ! dead coarse root growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_deadcroot_storage_gr_patch          (:)     ! dead coarse root growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_deadcroot_gr_patch               (:)     ! dead coarse root growth respiration from storage (gC/m2/s)

     ! growth respiration for prognostic crop model
     real(r8), pointer :: cpool_grain_gr_patch                      (:)     ! grain growth respiration (gC/m2/s)
     real(r8), pointer :: cpool_grain_storage_gr_patch              (:)     ! grain growth respiration to storage (gC/m2/s)
     real(r8), pointer :: transfer_grain_gr_patch                   (:)     ! grain growth respiration from storage (gC/m2/s)

     ! annual turnover of storage to transfer pools            
     real(r8), pointer :: grainc_storage_to_xfer_patch              (:)     ! grain C shift storage to transfer for prognostic crop model (gC/m2/s)
     real(r8), pointer :: leafc_storage_to_xfer_patch               (:)     ! leaf C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: frootc_storage_to_xfer_patch              (:)     ! fine root C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: livestemc_storage_to_xfer_patch           (:)     ! live stem C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: deadstemc_storage_to_xfer_patch           (:)     ! dead stem C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: livecrootc_storage_to_xfer_patch          (:)     ! live coarse root C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: deadcrootc_storage_to_xfer_patch          (:)     ! dead coarse root C shift storage to transfer (gC/m2/s)
     real(r8), pointer :: gresp_storage_to_xfer_patch               (:)     ! growth respiration shift storage to transfer (gC/m2/s)

     ! turnover of livewood to deadwood
     real(r8), pointer :: livestemc_to_deadstemc_patch              (:)     ! live stem C turnover (gC/m2/s)
     real(r8), pointer :: livecrootc_to_deadcrootc_patch            (:)     ! live coarse root C turnover (gC/m2/s)

     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_c_to_litr_met_c_col             (:,:)   ! C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
     real(r8), pointer :: phenology_c_to_litr_cel_c_col             (:,:)   ! C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
     real(r8), pointer :: phenology_c_to_litr_lig_c_col             (:,:)   ! C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)

     ! gap mortality
     real(r8), pointer :: gap_mortality_c_to_litr_met_c_col         (:,:)   ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
     real(r8), pointer :: gap_mortality_c_to_litr_cel_c_col         (:,:)   ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
     real(r8), pointer :: gap_mortality_c_to_litr_lig_c_col         (:,:)   ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
     real(r8), pointer :: gap_mortality_c_to_cwdc_col               (:,:)   ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)

     ! fire
     real(r8), pointer :: fire_mortality_c_to_cwdc_col              (:,:)   ! C fluxes associated with fire mortality to CWD pool (gC/m3/s)


     ! harvest
     real(r8), pointer :: harvest_c_to_litr_met_c_col               (:,:)   ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
     real(r8), pointer :: harvest_c_to_litr_cel_c_col               (:,:)   ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
     real(r8), pointer :: harvest_c_to_litr_lig_c_col               (:,:)   ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
     real(r8), pointer :: harvest_c_to_cwdc_col                     (:,:)   ! C fluxes associated with harvest to CWD pool (gC/m3/s)
     real(r8), pointer :: grainc_to_cropprodc_patch                 (:)     ! grain C to crop product pool (gC/m2/s)
     real(r8), pointer :: grainc_to_cropprodc_col                   (:)     ! grain C to crop product pool (gC/m2/s)

     ! fire fluxes
     real(r8), pointer :: m_decomp_cpools_to_fire_vr_col            (:,:,:) ! vertically-resolved decomposing C fire loss (gC/m3/s)
     real(r8), pointer :: m_decomp_cpools_to_fire_col               (:,:)   ! vertically-integrated (diagnostic) decomposing C fire loss (gC/m2/s)
     real(r8), pointer :: m_c_to_litr_met_fire_col                  (:,:)   ! C from leaf, froot, xfer and storage C to litter labile C by fire (gC/m3/s) 
     real(r8), pointer :: m_c_to_litr_cel_fire_col                  (:,:)   ! C from leaf, froot, xfer and storage C to litter cellulose C by fire (gC/m3/s) 
     real(r8), pointer :: m_c_to_litr_lig_fire_col                  (:,:)   ! C from leaf, froot, xfer and storage C to litter lignin C by fire (gC/m3/s) 

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedc_to_leaf_patch                   (:)     ! (gC/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedc_to_leaf_grc                     (:)     ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_seedc_to_deadstem_patch               (:)     ! (gC/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedc_to_deadstem_grc                 (:)     ! (gC/m2/s) dwt_seedc_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_cflux_patch                      (:)     ! (gC/m2/s) conversion C flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_conv_cflux_grc                        (:)     ! (gC/m2/s) dwt_conv_cflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_cflux_dribbled_grc               (:)     ! (gC/m2/s) dwt_conv_cflux_grc dribbled evenly throughout the year
     real(r8), pointer :: dwt_wood_productc_gain_patch              (:)     ! (gC/m2/s) addition to wood product pools from landcover change; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_crop_productc_gain_patch              (:)     ! (gC/m2/s) addition to crop product pools from landcover change; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_slash_cflux_patch                     (:)     ! (gC/m2/s) conversion slash flux due to landcover change
     real(r8), pointer :: dwt_slash_cflux_grc                       (:)     ! (gC/m2/s) dwt_slash_cflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_frootc_to_litr_met_c_col              (:,:)   ! (gC/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootc_to_litr_cel_c_col              (:,:)   ! (gC/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootc_to_litr_lig_c_col              (:,:)   ! (gC/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootc_to_cwdc_col                (:,:)   ! (gC/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootc_to_cwdc_col                (:,:)   ! (gC/m3/s) dead coarse root to CWD due to landcover change

     ! crop fluxes
     real(r8), pointer :: crop_seedc_to_leaf_patch                  (:)     ! (gC/m2/s) seed source to leaf, for crops

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: gpp_before_downreg_patch                  (:)     ! (gC/m2/s) gross primary production before down regulation
     real(r8), pointer :: current_gr_patch                          (:)     ! (gC/m2/s) growth resp for new growth displayed in this timestep
     real(r8), pointer :: transfer_gr_patch                         (:)     ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
     real(r8), pointer :: storage_gr_patch                          (:)     ! (gC/m2/s) growth resp for growth sent to storage for later display
     real(r8), pointer :: plant_calloc_patch                        (:)     ! (gC/m2/s) total allocated C flux 
     real(r8), pointer :: excess_cflux_patch                        (:)     ! (gC/m2/s) C flux not allocated due to downregulation 
     real(r8), pointer :: prev_leafc_to_litter_patch                (:)     ! (gC/m2/s) previous timestep leaf C litterfall flux 
     real(r8), pointer :: prev_frootc_to_litter_patch               (:)     ! (gC/m2/s) previous timestep froot C litterfall flux 
     real(r8), pointer :: availc_patch                              (:)     ! (gC/m2/s) C flux available for allocation 
     real(r8), pointer :: xsmrpool_recover_patch                    (:)     ! (gC/m2/s) C flux assigned to recovery of negative cpool
     real(r8), pointer :: xsmrpool_c13ratio_patch                   (:)     ! C13/C(12+13) ratio for xsmrpool (proportion)

     real(r8), pointer :: cwdc_hr_col                               (:)     ! (gC/m2/s) col-level coarse woody debris C heterotrophic respiration
     real(r8), pointer :: cwdc_loss_col                             (:)     ! (gC/m2/s) col-level coarse woody debris C loss
     real(r8), pointer :: litterc_loss_col                          (:)     ! (gC/m2/s) col-level litter C loss
     real(r8), pointer :: frootc_alloc_patch                        (:)     ! (gC/m2/s) patch-level fine root C alloc
     real(r8), pointer :: frootc_loss_patch                         (:)     ! (gC/m2/s) patch-level fine root C loss
     real(r8), pointer :: leafc_alloc_patch                         (:)     ! (gC/m2/s) patch-level leaf C alloc
     real(r8), pointer :: leafc_loss_patch                          (:)     ! (gC/m2/s) patch-level leaf C loss
     real(r8), pointer :: woodc_alloc_patch                         (:)     ! (gC/m2/s) patch-level wood C alloc
     real(r8), pointer :: woodc_loss_patch                          (:)     ! (gC/m2/s


     real(r8), pointer :: gpp_patch                                 (:)     ! (gC/m2/s) patch gross primary production 
     real(r8), pointer :: gpp_col                                   (:)     ! (gC/m2/s) column GPP flux before downregulation  (p2c)         
     real(r8), pointer :: rr_patch                                  (:)     ! (gC/m2/s) root respiration (fine root MR + total root GR)
     real(r8), pointer :: rr_col                                    (:)     ! (gC/m2/s) root respiration (fine root MR + total root GR) (p2c)
     real(r8), pointer :: mr_patch                                  (:)     ! (gC/m2/s) maintenance respiration
     real(r8), pointer :: gr_patch                                  (:)     ! (gC/m2/s) total growth respiration
     real(r8), pointer :: ar_patch                                  (:)     ! (gC/m2/s) patch autotrophic respiration (MR + GR)
     real(r8), pointer :: ar_col                                    (:)     ! (gC/m2/s) column autotrophic respiration (MR + GR) (p2c)      
     real(r8), pointer :: npp_patch                                 (:)     ! (gC/m2/s) patch net primary production
     real(r8), pointer :: npp_col                                   (:)     ! (gC/m2/s) column net primary production (p2c)                  
     real(r8), pointer :: agnpp_patch                               (:)     ! (gC/m2/s) aboveground NPP
     real(r8), pointer :: bgnpp_patch                               (:)     ! (gC/m2/s) belowground NPP
     real(r8), pointer :: litfall_patch                             (:)     ! (gC/m2/s) patch litterfall (leaves and fine roots)
     real(r8), pointer :: wood_harvestc_patch                       (:)     ! (gC/m2/s) patch-level wood harvest (to product pools)
     real(r8), pointer :: wood_harvestc_col                         (:)     ! (gC/m2/s) column-level wood harvest (to product pools) (p2c)
     real(r8), pointer :: slash_harvestc_patch                      (:)     ! (gC/m2/s) patch-level slash from harvest (to litter)
     real(r8), pointer :: cinputs_patch                             (:)     ! (gC/m2/s) patch-level carbon inputs (for balance checking)
     real(r8), pointer :: coutputs_patch                            (:)     ! (gC/m2/s) patch-level carbon outputs (for balance checking)
     real(r8), pointer :: sr_col                                    (:)     ! (gC/m2/s) total soil respiration (HR + root resp)
     real(r8), pointer :: er_col                                    (:)     ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
     real(r8), pointer :: litfire_col                               (:)     ! (gC/m2/s) litter fire losses
     real(r8), pointer :: somfire_col                               (:)     ! (gC/m2/s) soil organic matter fire losses
     real(r8), pointer :: totfire_col                               (:)     ! (gC/m2/s) total ecosystem fire losses
     real(r8), pointer :: hrv_xsmrpool_to_atm_col                   (:)     ! (gC/m2/s) excess MR pool harvest mortality (p2c)

     ! fire code
     real(r8), pointer :: fire_closs_patch                          (:)     ! (gC/m2/s) total fire C loss 
     real(r8), pointer :: fire_closs_p2c_col                        (:)     ! (gC/m2/s) patch2col averaged column-level fire C loss (p2c)
     real(r8), pointer :: fire_closs_col                            (:)     ! (gC/m2/s) total patch-level fire C loss 


     ! temporary and annual sums
     real(r8), pointer :: tempsum_litfall_patch                     (:)     ! (gC/m2/yr) temporary annual sum of litfall (CNDV only for now)
     real(r8), pointer :: annsum_litfall_patch                      (:)     ! (gC/m2/yr) annual sum of litfall (CNDV only for now)
     real(r8), pointer :: tempsum_npp_patch                         (:)     ! (gC/m2/yr) temporary annual sum of NPP 
     real(r8), pointer :: annsum_npp_patch                          (:)     ! (gC/m2/yr) annual sum of NPP 
     real(r8), pointer :: annsum_npp_col                            (:)     ! (gC/m2/yr) annual sum of NPP, averaged from patch-level
     real(r8), pointer :: lag_npp_col                               (:)     ! (gC/m2/yr) lagged net primary production

     ! Summary C fluxes. 
     real(r8), pointer :: nep_col        (:) ! (gC/m2/s) net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink
     real(r8), pointer :: nbp_grc        (:) ! (gC/m2/s) net biome production, includes fire, landuse, harvest and hrv_xsmrpool flux, positive for sink (same as net carbon exchange between land and atmosphere)
     real(r8), pointer :: nee_grc        (:) ! (gC/m2/s) net ecosystem exchange of carbon, includes fire and hrv_xsmrpool, excludes landuse and harvest flux, positive for source 

     ! Dynamic landcover fluxnes
     real(r8), pointer :: landuseflux_grc(:) ! (gC/m2/s) dwt_conv_cflux+product_closs
     real(r8), pointer :: npp_Nactive_patch                         (:)     ! C used by mycorrhizal uptake    (gC/m2/s)
     real(r8), pointer :: npp_burnedoff_patch                       (:)     ! C that cannot be used for N uptake   (gC/m2/s)
     real(r8), pointer :: npp_Nnonmyc_patch                         (:)     ! C used by non-myc uptake        (gC/m2/s)
     real(r8), pointer :: npp_Nam_patch                             (:)     ! C used by AM plant              (gC/m2/s)
     real(r8), pointer :: npp_Necm_patch                            (:)     ! C used by ECM plant             (gC/m2/s)
     real(r8), pointer :: npp_Nactive_no3_patch                     (:)     ! C used by mycorrhizal uptake    (gC/m2/s)
     real(r8), pointer :: npp_Nactive_nh4_patch                     (:)     ! C used by mycorrhizal uptake    (gC/m2/s)
     real(r8), pointer :: npp_Nnonmyc_no3_patch                     (:)     ! C used by non-myc               (gC/m2/s)
     real(r8), pointer :: npp_Nnonmyc_nh4_patch                     (:)     ! C used by non-myc               (gC/m2/s)
     real(r8), pointer :: npp_Nam_no3_patch                         (:)     ! C used by AM plant              (gC/m2/s)
     real(r8), pointer :: npp_Nam_nh4_patch                         (:)     ! C used by AM plant              (gC/m2/s)
     real(r8), pointer :: npp_Necm_no3_patch                        (:)     ! C used by ECM plant             (gC/m2/s)
     real(r8), pointer :: npp_Necm_nh4_patch                        (:)     ! C used by ECM plant             (gC/m2/s)
     real(r8), pointer :: npp_Nfix_patch                            (:)     ! C used by Symbiotic BNF         (gC/m2/s)
     real(r8), pointer :: npp_Nretrans_patch                        (:)     ! C used by retranslocation       (gC/m2/s)
     real(r8), pointer :: npp_Nuptake_patch                         (:)     ! Total C used by N uptake in FUN (gC/m2/s)
     real(r8), pointer :: npp_growth_patch                          (:)     ! Total C u for growth in FUN      (gC/m2/s)   
     real(r8), pointer :: leafc_change_patch                        (:)     ! Total used C from leaves        (gC/m2/s)
     real(r8), pointer :: soilc_change_patch                        (:)     ! Total used C from soil          (gC/m2/s)

     ! Matrix for C flux index
     real(r8), pointer :: matrix_Cinput_patch                       (:)      ! I-matrix for carbon input
     real(r8), pointer :: matrix_C13input_patch                     (:)      ! I-matrix for C13 input
     real(r8), pointer :: matrix_C14input_patch                     (:)      ! I-matrix for C14 input
     real(r8), pointer :: matrix_alloc_patch                        (:,:)    ! B-matrix for carbon allocation

     real(r8), pointer :: matrix_phtransfer_patch                   (:,:)    ! A-matrix_phenology
     real(r8), pointer :: matrix_phturnover_patch                   (:,:)    ! K-matrix_phenology
     integer,  pointer :: matrix_phtransfer_doner_patch             (:)      ! A-matrix_phenology non-zero indices (column indices)
     integer,  pointer :: matrix_phtransfer_receiver_patch          (:)      ! A-matrix_phenology non-zero indices (row indices)
     integer,  pointer :: actpatch_fire                             (:)      ! Patch indices with fire in current time step
     integer           :: num_actpatch_fire                                  ! Number of patches with fire in current time step

     real(r8), pointer :: matrix_gmtransfer_patch                   (:,:)    ! A-matrix_gap mortality
     real(r8), pointer :: matrix_gmturnover_patch                   (:,:)    ! K-matrix_gap mortality
     integer,  pointer :: matrix_gmtransfer_doner_patch             (:)      ! A-matrix_gap mortality non-zero indices (column indices)
     integer,  pointer :: matrix_gmtransfer_receiver_patch          (:)      ! A-matrix_gap mortality non-zero indices (row indices)

     real(r8), pointer :: matrix_fitransfer_patch                   (:,:)    ! A-matrix_fire
     real(r8), pointer :: matrix_fiturnover_patch                   (:,:)    ! K-matrix_fire     
     integer,  pointer :: matrix_fitransfer_doner_patch             (:)      ! A-matrix_fire non-zero indices (column indices)
     integer,  pointer :: matrix_fitransfer_receiver_patch          (:)      ! A-matrix_fire non-zero indices (row indices)

!     real(r8), pointer :: soilc_change_col                          (:)     ! Total used C from soil          (gC/m2/s)
!   matrix variables 
     integer ileafst_to_ileafxf_ph                    ! Index of phenology related C transfer from leaf storage pool to leaf transfer pool
     integer ileafxf_to_ileaf_ph                      ! Index of phenology related C transfer from leaf transfer pool to leaf pool
     integer ifrootst_to_ifrootxf_ph                  ! Index of phenology related C transfer from fine root storage pool to fine root transfer pool
     integer ifrootxf_to_ifroot_ph                    ! Index of phenology related C transfer from fine root transfer pool to fine root pool
     integer ilivestemst_to_ilivestemxf_ph            ! Index of phenology related C transfer from live stem storage pool to live stem transfer pool
     integer ilivestemxf_to_ilivestem_ph              ! Index of phenology related C transfer from live stem transfer pool to live stem pool
     integer ideadstemst_to_ideadstemxf_ph            ! Index of phenology related C transfer from dead stem storage pool to dead stem transfer pool
     integer ideadstemxf_to_ideadstem_ph              ! Index of phenology related C transfer from dead stem transfer pool to dead stem pool
     integer ilivecrootst_to_ilivecrootxf_ph          ! Index of phenology related C transfer from live coarse root storage pool to live coarse root transfer pool
     integer ilivecrootxf_to_ilivecroot_ph            ! Index of phenology related C transfer from live coarse root transfer pool to live coarse root pool
     integer ideadcrootst_to_ideadcrootxf_ph          ! Index of phenology related C transfer from dead coarse root storage pool to dead coarse root transfer pool
     integer ideadcrootxf_to_ideadcroot_ph            ! Index of phenology related C transfer from dead coarse root transfer pool to dead coarse root pool
     integer ilivestem_to_ideadstem_ph                ! Index of phenology related C transfer from live stem pool to dead stem pool
     integer ilivecroot_to_ideadcroot_ph              ! Index of phenology related C transfer from live coarse root pool to dead coarse root pool
     integer ileaf_to_iout_ph                         ! Index of phenology related C transfer from leaf pool to outside of vegetation pools
     integer ifroot_to_iout_ph                        ! Index of phenology related C transfer from fine root pool to outside of vegetation pools
     integer ilivestem_to_iout_ph                     ! Index of phenology related C transfer from live stem pool to outside of vegetation pools 
     integer igrain_to_iout_ph                        ! Index of phenology related C transfer from grain pool to outside of vegetation pools
     integer ileaf_to_iout_gm                         ! Index of gap mortality related C transfer from leaf pool to outside of vegetation pools
     integer ileafst_to_iout_gm                       ! Index of gap mortality related C transfer from leaf storage pool to outside of vegetation pools               
     integer ileafxf_to_iout_gm                       ! Index of gap mortality related C transfer from leaf transfer pool to outside of vegetation pools
     integer ifroot_to_iout_gm                        ! Index of gap mortality related C transfer from fine root pool to outside of vegetation pools
     integer ifrootst_to_iout_gm                      ! Index of gap mortality related C transfer from fine root storage pool to outside of vegetation pools
     integer ifrootxf_to_iout_gm                      ! Index of gap mortality related C transfer from fine root transfer pool to outside of vegetation pools
     integer ilivestem_to_iout_gm                     ! Index of gap mortality related C transfer from live stem pool to outside of vegetation pools
     integer ilivestemst_to_iout_gm                   ! Index of gap mortality related C transfer from live stem storage pool to outside of vegetation pools
     integer ilivestemxf_to_iout_gm                   ! Index of gap mortality related C transfer from live stem transfer pool to outside of vegetation pools
     integer ideadstem_to_iout_gm                     ! Index of gap mortality related C transfer from dead stem pool to outside of vegetation pools
     integer ideadstemst_to_iout_gm                   ! Index of gap mortality related C transfer from dead stem storage pool to outside of vegetation pools
     integer ideadstemxf_to_iout_gm                   ! Index of gap mortality related C transfer from dead stem transfer pool to outside of vegetation pools
     integer ilivecroot_to_iout_gm                    ! Index of gap mortality related C transfer from live coarse root pool to outside of vegetation pools
     integer ilivecrootst_to_iout_gm                  ! Index of gap mortality related C transfer from live coarse root storage pool to outside of vegetation pools
     integer ilivecrootxf_to_iout_gm                  ! Index of gap mortality related C transfer from live coarse root transfer pool to outside of vegetation pools
     integer ideadcroot_to_iout_gm                    ! Index of gap mortality related C transfer from dead coarse root pool to outside of vegetation pools
     integer ideadcrootst_to_iout_gm                  ! Index of gap mortality related C transfer from dead coarse root storage pool to outside of vegetation pools
     integer ideadcrootxf_to_iout_gm                  ! Index of gap mortality related C transfer from dead coarse root transfer pool to outside of vegetation pools
     integer ileaf_to_iout_fi                         ! Index of fire related C transfer from leaf pool to outside of vegetation pools
     integer ileafst_to_iout_fi                       ! Index of fire related C transfer from leaf storage pool to outside of vegetation pools
     integer ileafxf_to_iout_fi                       ! Index of fire related C transfer from leaf transfer pool to outside of vegetation pools
     integer ifroot_to_iout_fi                        ! Index of fire related C transfer from fine root pool to outside of vegetation pools
     integer ifrootst_to_iout_fi                      ! Index of fire related C transfer from fine root storage pool to outside of vegetation pools
     integer ifrootxf_to_iout_fi                      ! Index of fire related C transfer from fine root transfer pool to outside of vegetation pools
     integer ilivestem_to_iout_fi                     ! Index of fire related C transfer from live stem pool to outside of vegetation pools
     integer ilivestemst_to_iout_fi                   ! Index of fire related C transfer from live stem storage pool to outside of vegetation pools
     integer ilivestemxf_to_iout_fi                   ! Index of fire related C transfer from live stem transfer pool to outside of vegetation pools
     integer ideadstem_to_iout_fi                     ! Index of fire related C transfer from dead stem pool to outside of vegetation pools
     integer ideadstemst_to_iout_fi                   ! Index of fire related C transfer from dead stem storage pool to outside of vegetation pools
     integer ideadstemxf_to_iout_fi                   ! Index of fire related C transfer from dead stem transfer pool to outside of vegetation pools
     integer ilivecroot_to_iout_fi                    ! Index of fire related C transfer from live coarse root pool to outside of vegetation pools
     integer ilivecrootst_to_iout_fi                  ! Index of fire related C transfer from live coarse root storage pool to outside of vegetation pools
     integer ilivecrootxf_to_iout_fi                  ! Index of fire related C transfer from live coarse root transfer pool to outside of vegetation pools
     integer ideadcroot_to_iout_fi                    ! Index of fire related C transfer from dead coarse root pool to outside of vegetation pools
     integer ideadcrootst_to_iout_fi                  ! Index of fire related C transfer from dead coarse root storage pool to outside of vegetation pools
     integer ideadcrootxf_to_iout_fi                  ! Index of fire related C transfer from dead coarse root transfer pool to outside of vegetation pools
     integer ilivestem_to_ideadstem_fi                ! Index of fire related C transfer from live stem pool to dead stem pools
     integer ilivecroot_to_ideadcroot_fi              ! Index of fire related C transfer from live coarse root pool to dead coarse root pools

     integer,pointer :: list_phc_phgmc     (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKphc to AKphc+AKgmc
     integer,pointer :: list_gmc_phgmc     (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKgmc to AKphc+AKgmc
     integer,pointer :: list_phc_phgmfic   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKphc to AKphc+AKgmc+AKfic
     integer,pointer :: list_gmc_phgmfic   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKgmc to AKphc+AKgmc+AKfic
     integer,pointer :: list_fic_phgmfic   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKfic to AKphc+AKgmc+AKfic
     integer,pointer :: list_aphc          (:)        ! Indices of non-diagnoal entries in full sparse matrix Aph for C cycle
     integer,pointer :: list_agmc          (:)        ! Indices of non-diagnoal entries in full sparse matrix Agm for C cycle
     integer,pointer :: list_afic          (:)        ! Indices of non-diagnoal entries in full sparse matrix Afi for C cycle

!     type(sparse_matrix_type)      :: AKphvegc        ! Aph*Kph for C cycle in sparse matrix format
!     type(sparse_matrix_type)      :: AKgmvegc        ! Agm*Kgm for C cycle in sparse matrix format
!     type(sparse_matrix_type)      :: AKfivegc        ! Afi*Kfi for C cycle in sparse matrix format
!     type(sparse_matrix_type)      :: AKallvegc       ! Aph*Kph + Agm*Kgm + Afi*Kfi for C cycle in sparse matrix format 
!
!     type(vector_type)             :: Xvegc           ! Vegetation C of each compartment in a vector format
!     type(vector_type)             :: Xveg13c         ! Vegetation C13 of each compartment in a vector format
!     type(vector_type)             :: Xveg14c         ! Vegetation C14 of each compartment in a vector format

     ! Objects that help convert once-per-year dynamic land cover changes into fluxes
     ! that are dribbled throughout the year
     type(annual_flux_dribbler_type) :: dwt_conv_cflux_dribbler
     type(annual_flux_dribbler_type) :: hrv_xsmrpool_to_atm_dribbler
     logical, private  :: dribble_crophrv_xsmrpool_2atm

 contains

     procedure , public  :: SetValues
     procedure , public  :: Summary => Summary_carbonflux
     procedure , public  :: ZeroDWT
     procedure , public  :: Init

 end type cnveg_carbonflux_type

type(cnveg_carbonflux_type), public, target, save :: cnveg_carbonflux_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

!---------------------------------------
 subroutine Init(this, bounds, nch, ityp, fveg, cncol, cnpft, carbon_type, cn5_cold_start, rc)

! !DESCRIPTION:
! Initialize CTSM carbon fluxes
! jk Apr 2021: type is allocated and initialized to NaN;
! if data arrays from restart file are passed (cncol and cnpft), the type is then initialized with these values
!
! !ARGUMENTS:
    implicit none

  ! INPUT
    type(bounds_type),                            intent(in) :: bounds
    integer,                                      intent(in) :: nch ! number of tiles
    integer, dimension(nch,NUM_VEG,NUM_ZON),      intent(in) :: ityp ! PFT index
    real, dimension(nch,NUM_VEG,NUM_ZON),         intent(in) :: fveg    ! PFT fraction
    real, dimension(nch,NUM_ZON,VAR_COL),         intent(in) :: cncol ! gkw: column CN restart
    real, dimension(nch,NUM_ZON,NUM_VEG,VAR_PFT), intent(in) :: cnpft ! gkw: PFT CN restart
    character(len=*) , intent(in) :: carbon_type ! one of ['c12', c13','c14']
    logical, optional,                            intent(in) :: cn5_cold_start
    class(cnveg_carbonflux_type)                              :: this
    integer, optional,                            intent(out) :: rc

    ! LOCAL
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: begg, endg
    integer  :: np, nc, nz, p, nv, n, nl
    logical  :: cold_start = .false.
    logical           :: allows_non_annual_delta
    character(len=:), allocatable :: carbon_type_suffix
    !--------------------------------------------------------

    ! check whether a cn5_cold_start option was set and change cold_start accordingly
    if (present(cn5_cold_start) .and. (cn5_cold_start.eqv..true.)) then
       cold_start = .true.
    end if

    ! jkolassa: if cold_start is false, check that both CNCOL and CNPFT have the expected size for CNCLM50, else abort 
    if ((cold_start.eqv..false.) .and. ((size(cncol,3).ne.var_col) .or. &
       (size(cnpft,4).ne.var_pft))) then
       _ASSERT(.FALSE.,'option CNCLM50_cold_start = .FALSE. requires a CNCLM50 restart file')
    end if

    allocate(this%matrix_phtransfer_doner_patch(1:18))
    allocate(this%matrix_phtransfer_receiver_patch(1:18))

    this%ileafst_to_ileafxf_ph           = 1
    this%matrix_phtransfer_doner_patch(this%ileafst_to_ileafxf_ph)              = ileaf_st
    this%matrix_phtransfer_receiver_patch(this%ileafst_to_ileafxf_ph)           = ileaf_xf

    this%ileafxf_to_ileaf_ph             = 2
    this%matrix_phtransfer_doner_patch(this%ileafxf_to_ileaf_ph)                = ileaf_xf
    this%matrix_phtransfer_receiver_patch(this%ileafxf_to_ileaf_ph)             = ileaf

    this%ifrootst_to_ifrootxf_ph         = 3
    this%matrix_phtransfer_doner_patch(this%ifrootst_to_ifrootxf_ph)            = ifroot_st
    this%matrix_phtransfer_receiver_patch(this%ifrootst_to_ifrootxf_ph)         = ifroot_xf

    this%ifrootxf_to_ifroot_ph           = 4
    this%matrix_phtransfer_doner_patch(this%ifrootxf_to_ifroot_ph)              = ifroot_xf
    this%matrix_phtransfer_receiver_patch(this%ifrootxf_to_ifroot_ph)           = ifroot

    this%ilivestem_to_ideadstem_ph       = 5
    this%matrix_phtransfer_doner_patch(this%ilivestem_to_ideadstem_ph)          = ilivestem
    this%matrix_phtransfer_receiver_patch(this%ilivestem_to_ideadstem_ph)       = ideadstem

    this%ilivestemst_to_ilivestemxf_ph   = 6
    this%matrix_phtransfer_doner_patch(this%ilivestemst_to_ilivestemxf_ph)      = ilivestem_st
    this%matrix_phtransfer_receiver_patch(this%ilivestemst_to_ilivestemxf_ph)   = ilivestem_xf

    this%ilivestemxf_to_ilivestem_ph     = 7
    this%matrix_phtransfer_doner_patch(this%ilivestemxf_to_ilivestem_ph)        = ilivestem_xf
    this%matrix_phtransfer_receiver_patch(this%ilivestemxf_to_ilivestem_ph)     = ilivestem

    this%ideadstemst_to_ideadstemxf_ph   = 8
    this%matrix_phtransfer_doner_patch(this%ideadstemst_to_ideadstemxf_ph)      = ideadstem_st
    this%matrix_phtransfer_receiver_patch(this%ideadstemst_to_ideadstemxf_ph)   = ideadstem_xf

    this%ideadstemxf_to_ideadstem_ph     = 9
    this%matrix_phtransfer_doner_patch(this%ideadstemxf_to_ideadstem_ph)        = ideadstem_xf
    this%matrix_phtransfer_receiver_patch(this%ideadstemxf_to_ideadstem_ph)     = ideadstem

    this%ilivecroot_to_ideadcroot_ph     = 10
    this%matrix_phtransfer_doner_patch(this%ilivecroot_to_ideadcroot_ph)        = ilivecroot
    this%matrix_phtransfer_receiver_patch(this%ilivecroot_to_ideadcroot_ph)     = ideadcroot

    this%ilivecrootst_to_ilivecrootxf_ph = 11
    this%matrix_phtransfer_doner_patch(this%ilivecrootst_to_ilivecrootxf_ph)    = ilivecroot_st
    this%matrix_phtransfer_receiver_patch(this%ilivecrootst_to_ilivecrootxf_ph) = ilivecroot_xf

    this%ilivecrootxf_to_ilivecroot_ph   = 12
    this%matrix_phtransfer_doner_patch(this%ilivecrootxf_to_ilivecroot_ph)      = ilivecroot_xf
    this%matrix_phtransfer_receiver_patch(this%ilivecrootxf_to_ilivecroot_ph)   = ilivecroot

    this%ideadcrootst_to_ideadcrootxf_ph = 13
    this%matrix_phtransfer_doner_patch(this%ideadcrootst_to_ideadcrootxf_ph)    = ideadcroot_st
    this%matrix_phtransfer_receiver_patch(this%ideadcrootst_to_ideadcrootxf_ph) = ideadcroot_xf

    this%ideadcrootxf_to_ideadcroot_ph   = 14
    this%matrix_phtransfer_doner_patch(this%ideadcrootxf_to_ideadcroot_ph)      = ideadcroot_xf
    this%matrix_phtransfer_receiver_patch(this%ideadcrootxf_to_ideadcroot_ph)   = ideadcroot

    this%ileaf_to_iout_ph                = 15
    this%matrix_phtransfer_doner_patch(this%ileaf_to_iout_ph)                   = ileaf
    this%matrix_phtransfer_receiver_patch(this%ileaf_to_iout_ph)                = ioutc

    this%ifroot_to_iout_ph               = 16
    this%matrix_phtransfer_doner_patch(this%ifroot_to_iout_ph)                  = ifroot
    this%matrix_phtransfer_receiver_patch(this%ifroot_to_iout_ph)               = ioutc

    this%ilivestem_to_iout_ph            = 17
    this%matrix_phtransfer_doner_patch(this%ilivestem_to_iout_ph)               = ilivestem
    this%matrix_phtransfer_receiver_patch(this%ilivestem_to_iout_ph)            = ioutc

    if(use_crop)then
       this%igrain_to_iout_ph            = 18
       this%matrix_phtransfer_doner_patch(this%igrain_to_iout_ph)               = igrain
       this%matrix_phtransfer_receiver_patch(this%igrain_to_iout_ph)            = ioutc
    end if

    allocate(this%matrix_gmtransfer_doner_patch(1:18))
    allocate(this%matrix_gmtransfer_receiver_patch(1:18))

    this%ileaf_to_iout_gm                = 1
    this%matrix_gmtransfer_doner_patch(this%ileaf_to_iout_gm)                   = ileaf
    this%matrix_gmtransfer_receiver_patch(this%ileaf_to_iout_gm)                = ioutc

    this%ileafst_to_iout_gm              = 2
    this%matrix_gmtransfer_doner_patch(this%ileafst_to_iout_gm)                 = ileaf_st
    this%matrix_gmtransfer_receiver_patch(this%ileafst_to_iout_gm)              = ioutc

    this%ileafxf_to_iout_gm              = 3
    this%matrix_gmtransfer_doner_patch(this%ileafxf_to_iout_gm)                 = ileaf_xf
    this%matrix_gmtransfer_receiver_patch(this%ileafxf_to_iout_gm)              = ioutc

    this%ifroot_to_iout_gm               = 4
    this%matrix_gmtransfer_doner_patch(this%ifroot_to_iout_gm)                  = ifroot
    this%matrix_gmtransfer_receiver_patch(this%ifroot_to_iout_gm)               = ioutc

    this%ifrootst_to_iout_gm             = 5
    this%matrix_gmtransfer_doner_patch(this%ifrootst_to_iout_gm)                = ifroot_st
    this%matrix_gmtransfer_receiver_patch(this%ifrootst_to_iout_gm)             = ioutc

    this%ifrootxf_to_iout_gm             = 6
    this%matrix_gmtransfer_doner_patch(this%ifrootxf_to_iout_gm)                = ifroot_xf
    this%matrix_gmtransfer_receiver_patch(this%ifrootxf_to_iout_gm)             = ioutc

    this%ilivestem_to_iout_gm            = 7
    this%matrix_gmtransfer_doner_patch(this%ilivestem_to_iout_gm)               = ilivestem
    this%matrix_gmtransfer_receiver_patch(this%ilivestem_to_iout_gm)            = ioutc

    this%ilivestemst_to_iout_gm          = 8
    this%matrix_gmtransfer_doner_patch(this%ilivestemst_to_iout_gm)             = ilivestem_st
    this%matrix_gmtransfer_receiver_patch(this%ilivestemst_to_iout_gm)          = ioutc

    this%ilivestemxf_to_iout_gm          = 9
    this%matrix_gmtransfer_doner_patch(this%ilivestemxf_to_iout_gm)             = ilivestem_xf
    this%matrix_gmtransfer_receiver_patch(this%ilivestemxf_to_iout_gm)          = ioutc

    this%ideadstem_to_iout_gm            = 10
    this%matrix_gmtransfer_doner_patch(this%ideadstem_to_iout_gm)               = ideadstem
    this%matrix_gmtransfer_receiver_patch(this%ideadstem_to_iout_gm)            = ioutc

    this%ideadstemst_to_iout_gm          = 11
    this%matrix_gmtransfer_doner_patch(this%ideadstemst_to_iout_gm)             = ideadstem_st
    this%matrix_gmtransfer_receiver_patch(this%ideadstemst_to_iout_gm)          = ioutc

    this%ideadstemxf_to_iout_gm          = 12
    this%matrix_gmtransfer_doner_patch(this%ideadstemxf_to_iout_gm)             = ideadstem_xf
    this%matrix_gmtransfer_receiver_patch(this%ideadstemxf_to_iout_gm)          = ioutc

    this%ilivecroot_to_iout_gm           = 13
    this%matrix_gmtransfer_doner_patch(this%ilivecroot_to_iout_gm)              = ilivecroot
    this%matrix_gmtransfer_receiver_patch(this%ilivecroot_to_iout_gm)           = ioutc

    this%ilivecrootst_to_iout_gm         = 14
    this%matrix_gmtransfer_doner_patch(this%ilivecrootst_to_iout_gm)            = ilivecroot_st
    this%matrix_gmtransfer_receiver_patch(this%ilivecrootst_to_iout_gm)         = ioutc

    this%ilivecrootxf_to_iout_gm         = 15
    this%matrix_gmtransfer_doner_patch(this%ilivecrootxf_to_iout_gm)            = ilivecroot_xf
    this%matrix_gmtransfer_receiver_patch(this%ilivecrootxf_to_iout_gm)         = ioutc

    this%ideadcroot_to_iout_gm           = 16
    this%matrix_gmtransfer_doner_patch(this%ideadcroot_to_iout_gm)              = ideadcroot
    this%matrix_gmtransfer_receiver_patch(this%ideadcroot_to_iout_gm)           = ioutc

    this%ideadcrootst_to_iout_gm         = 17
    this%matrix_gmtransfer_doner_patch(this%ideadcrootst_to_iout_gm)            = ideadcroot_st
    this%matrix_gmtransfer_receiver_patch(this%ideadcrootst_to_iout_gm)         = ioutc

    this%ideadcrootxf_to_iout_gm         = 18
    this%matrix_gmtransfer_doner_patch(this%ideadcrootxf_to_iout_gm)            = ideadcroot_xf
    this%matrix_gmtransfer_receiver_patch(this%ideadcrootxf_to_iout_gm)         = ioutc

    allocate(this%matrix_fitransfer_doner_patch(1:20))
    allocate(this%matrix_fitransfer_receiver_patch(1:20))

    this%ilivestem_to_ideadstem_fi       = 1
    this%matrix_fitransfer_doner_patch(this%ilivestem_to_ideadstem_fi)          = ilivestem
    this%matrix_fitransfer_receiver_patch(this%ilivestem_to_ideadstem_fi)       = ideadstem

    this%ilivecroot_to_ideadcroot_fi     = 2
    this%matrix_fitransfer_doner_patch(this%ilivecroot_to_ideadcroot_fi)        = ilivecroot
    this%matrix_fitransfer_receiver_patch(this%ilivecroot_to_ideadcroot_fi)     = ideadcroot

    this%ileaf_to_iout_fi                = 3
    this%matrix_fitransfer_doner_patch(this%ileaf_to_iout_fi)                   = ileaf
    this%matrix_fitransfer_receiver_patch(this%ileaf_to_iout_fi)                = ioutc

    this%ileafst_to_iout_fi              = 4
    this%matrix_fitransfer_doner_patch(this%ileafst_to_iout_fi)                 = ileaf_st
    this%matrix_fitransfer_receiver_patch(this%ileafst_to_iout_fi)              = ioutc

    this%ileafxf_to_iout_fi              = 5
    this%matrix_fitransfer_doner_patch(this%ileafxf_to_iout_fi)                 = ileaf_xf
    this%matrix_fitransfer_receiver_patch(this%ileafxf_to_iout_fi)              = ioutc

    this%ifroot_to_iout_fi               = 6
    this%matrix_fitransfer_doner_patch(this%ifroot_to_iout_fi)                  = ifroot
    this%matrix_fitransfer_receiver_patch(this%ifroot_to_iout_fi)               = ioutc

    this%ifrootst_to_iout_fi             = 7
    this%matrix_fitransfer_doner_patch(this%ifrootst_to_iout_fi)                = ifroot_st
    this%matrix_fitransfer_receiver_patch(this%ifrootst_to_iout_fi)             = ioutc

    this%ifrootxf_to_iout_fi             = 8
    this%matrix_fitransfer_doner_patch(this%ifrootxf_to_iout_fi)                = ifroot_xf
    this%matrix_fitransfer_receiver_patch(this%ifrootxf_to_iout_fi)             = ioutc

    this%ilivestem_to_iout_fi            = 9
    this%matrix_fitransfer_doner_patch(this%ilivestem_to_iout_fi)               = ilivestem
    this%matrix_fitransfer_receiver_patch(this%ilivestem_to_iout_fi)            = ioutc

    this%ilivestemst_to_iout_fi          = 10
    this%matrix_fitransfer_doner_patch(this%ilivestemst_to_iout_fi)             = ilivestem_st
    this%matrix_fitransfer_receiver_patch(this%ilivestemst_to_iout_fi)          = ioutc

    this%ilivestemxf_to_iout_fi          = 11
    this%matrix_fitransfer_doner_patch(this%ilivestemxf_to_iout_fi)             = ilivestem_xf
    this%matrix_fitransfer_receiver_patch(this%ilivestemxf_to_iout_fi)          = ioutc

    this%ideadstem_to_iout_fi            = 12
    this%matrix_fitransfer_doner_patch(this%ideadstem_to_iout_fi)               = ideadstem
    this%matrix_fitransfer_receiver_patch(this%ideadstem_to_iout_fi)            = ioutc

    this%ideadstemst_to_iout_fi          = 13
    this%matrix_fitransfer_doner_patch(this%ideadstemst_to_iout_fi)             = ideadstem_st
    this%matrix_fitransfer_receiver_patch(this%ideadstemst_to_iout_fi)          = ioutc

    this%ideadstemxf_to_iout_fi          = 14
    this%matrix_fitransfer_doner_patch(this%ideadstemxf_to_iout_fi)             = ideadstem_xf
    this%matrix_fitransfer_receiver_patch(this%ideadstemxf_to_iout_fi)          = ioutc

    this%ilivecroot_to_iout_fi           = 15
    this%matrix_fitransfer_doner_patch(this%ilivecroot_to_iout_fi)              = ilivecroot
    this%matrix_fitransfer_receiver_patch(this%ilivecroot_to_iout_fi)           = ioutc

    this%ilivecrootst_to_iout_fi         = 16
    this%matrix_fitransfer_doner_patch(this%ilivecrootst_to_iout_fi)            = ilivecroot_st
    this%matrix_fitransfer_receiver_patch(this%ilivecrootst_to_iout_fi)         = ioutc

    this%ilivecrootxf_to_iout_fi         = 17
    this%matrix_fitransfer_doner_patch(this%ilivecrootxf_to_iout_fi)            = ilivecroot_xf
    this%matrix_fitransfer_receiver_patch(this%ilivecrootxf_to_iout_fi)         = ioutc

    this%ideadcroot_to_iout_fi           = 18
    this%matrix_fitransfer_doner_patch(this%ideadcroot_to_iout_fi)              = ideadcroot
    this%matrix_fitransfer_receiver_patch(this%ideadcroot_to_iout_fi)           = ioutc


    this%ideadcrootst_to_iout_fi         = 19
    this%matrix_fitransfer_doner_patch(this%ideadcrootst_to_iout_fi)            = ideadcroot_st
    this%matrix_fitransfer_receiver_patch(this%ideadcrootst_to_iout_fi)         = ioutc

    this%ideadcrootxf_to_iout_fi         = 20
    this%matrix_fitransfer_doner_patch(this%ideadcrootxf_to_iout_fi)            = ideadcroot_xf
    this%matrix_fitransfer_receiver_patch(this%ideadcrootxf_to_iout_fi)         = ioutc

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%m_leafc_to_litter_patch                   (begp:endp)) ; this%m_leafc_to_litter_patch                   (:) = nan
    allocate(this%m_frootc_to_litter_patch                  (begp:endp)) ; this%m_frootc_to_litter_patch                  (:) = nan
    allocate(this%m_leafc_storage_to_litter_patch           (begp:endp)) ; this%m_leafc_storage_to_litter_patch           (:) = nan
    allocate(this%m_frootc_storage_to_litter_patch          (begp:endp)) ; this%m_frootc_storage_to_litter_patch          (:) = nan
    allocate(this%m_livestemc_storage_to_litter_patch       (begp:endp)) ; this%m_livestemc_storage_to_litter_patch       (:) = nan
    allocate(this%m_deadstemc_storage_to_litter_patch       (begp:endp)) ; this%m_deadstemc_storage_to_litter_patch       (:) = nan
    allocate(this%m_livecrootc_storage_to_litter_patch      (begp:endp)) ; this%m_livecrootc_storage_to_litter_patch      (:) = nan
    allocate(this%m_deadcrootc_storage_to_litter_patch      (begp:endp)) ; this%m_deadcrootc_storage_to_litter_patch      (:) = nan
    allocate(this%m_leafc_xfer_to_litter_patch              (begp:endp)) ; this%m_leafc_xfer_to_litter_patch              (:) = nan
    allocate(this%m_frootc_xfer_to_litter_patch             (begp:endp)) ; this%m_frootc_xfer_to_litter_patch             (:) = nan
    allocate(this%m_livestemc_xfer_to_litter_patch          (begp:endp)) ; this%m_livestemc_xfer_to_litter_patch          (:) = nan
    allocate(this%m_deadstemc_xfer_to_litter_patch          (begp:endp)) ; this%m_deadstemc_xfer_to_litter_patch          (:) = nan
    allocate(this%m_livecrootc_xfer_to_litter_patch         (begp:endp)) ; this%m_livecrootc_xfer_to_litter_patch         (:) = nan
    allocate(this%m_deadcrootc_xfer_to_litter_patch         (begp:endp)) ; this%m_deadcrootc_xfer_to_litter_patch         (:) = nan
    allocate(this%m_livestemc_to_litter_patch               (begp:endp)) ; this%m_livestemc_to_litter_patch               (:) = nan
    allocate(this%m_deadstemc_to_litter_patch               (begp:endp)) ; this%m_deadstemc_to_litter_patch               (:) = nan
    allocate(this%m_livecrootc_to_litter_patch              (begp:endp)) ; this%m_livecrootc_to_litter_patch              (:) = nan
    allocate(this%m_deadcrootc_to_litter_patch              (begp:endp)) ; this%m_deadcrootc_to_litter_patch              (:) = nan
    allocate(this%m_gresp_storage_to_litter_patch           (begp:endp)) ; this%m_gresp_storage_to_litter_patch           (:) = nan
    allocate(this%m_gresp_xfer_to_litter_patch              (begp:endp)) ; this%m_gresp_xfer_to_litter_patch              (:) = nan
    allocate(this%hrv_leafc_to_litter_patch                 (begp:endp)) ; this%hrv_leafc_to_litter_patch                 (:) = nan
    allocate(this%hrv_leafc_storage_to_litter_patch         (begp:endp)) ; this%hrv_leafc_storage_to_litter_patch         (:) = nan
    allocate(this%hrv_leafc_xfer_to_litter_patch            (begp:endp)) ; this%hrv_leafc_xfer_to_litter_patch            (:) = nan
    allocate(this%hrv_frootc_to_litter_patch                (begp:endp)) ; this%hrv_frootc_to_litter_patch                (:) = nan
    allocate(this%hrv_frootc_storage_to_litter_patch        (begp:endp)) ; this%hrv_frootc_storage_to_litter_patch        (:) = nan
    allocate(this%hrv_frootc_xfer_to_litter_patch           (begp:endp)) ; this%hrv_frootc_xfer_to_litter_patch           (:) = nan
    allocate(this%hrv_livestemc_to_litter_patch             (begp:endp)) ; this%hrv_livestemc_to_litter_patch             (:) = nan
    allocate(this%hrv_livestemc_storage_to_litter_patch     (begp:endp)) ; this%hrv_livestemc_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_livestemc_xfer_to_litter_patch        (begp:endp)) ; this%hrv_livestemc_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_deadstemc_storage_to_litter_patch     (begp:endp)) ; this%hrv_deadstemc_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_deadstemc_xfer_to_litter_patch        (begp:endp)) ; this%hrv_deadstemc_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_livecrootc_to_litter_patch            (begp:endp)) ; this%hrv_livecrootc_to_litter_patch            (:) = nan
    allocate(this%hrv_livecrootc_storage_to_litter_patch    (begp:endp)) ; this%hrv_livecrootc_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_livecrootc_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livecrootc_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_deadcrootc_to_litter_patch            (begp:endp)) ; this%hrv_deadcrootc_to_litter_patch            (:) = nan
    allocate(this%hrv_deadcrootc_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadcrootc_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_deadcrootc_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadcrootc_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_gresp_storage_to_litter_patch         (begp:endp)) ; this%hrv_gresp_storage_to_litter_patch         (:) = nan
    allocate(this%hrv_gresp_xfer_to_litter_patch            (begp:endp)) ; this%hrv_gresp_xfer_to_litter_patch            (:) = nan
    allocate(this%hrv_xsmrpool_to_atm_patch                 (begp:endp)) ; this%hrv_xsmrpool_to_atm_patch                 (:) = 0.0_r8
    allocate(this%m_leafc_to_fire_patch                     (begp:endp)) ; this%m_leafc_to_fire_patch                     (:) = nan
    allocate(this%m_leafc_storage_to_fire_patch             (begp:endp)) ; this%m_leafc_storage_to_fire_patch             (:) = nan
    allocate(this%m_leafc_xfer_to_fire_patch                (begp:endp)) ; this%m_leafc_xfer_to_fire_patch                (:) = nan
    allocate(this%m_livestemc_to_fire_patch                 (begp:endp)) ; this%m_livestemc_to_fire_patch                 (:) = nan
    allocate(this%m_livestemc_storage_to_fire_patch         (begp:endp)) ; this%m_livestemc_storage_to_fire_patch         (:) = nan
    allocate(this%m_livestemc_xfer_to_fire_patch            (begp:endp)) ; this%m_livestemc_xfer_to_fire_patch            (:) = nan
    allocate(this%m_deadstemc_to_fire_patch                 (begp:endp)) ; this%m_deadstemc_to_fire_patch                 (:) = nan
    allocate(this%m_deadstemc_storage_to_fire_patch         (begp:endp)) ; this%m_deadstemc_storage_to_fire_patch         (:) = nan
    allocate(this%m_deadstemc_xfer_to_fire_patch            (begp:endp)) ; this%m_deadstemc_xfer_to_fire_patch            (:) = nan
    allocate(this%m_frootc_to_fire_patch                    (begp:endp)) ; this%m_frootc_to_fire_patch                    (:) = nan
    allocate(this%m_frootc_storage_to_fire_patch            (begp:endp)) ; this%m_frootc_storage_to_fire_patch            (:) = nan
    allocate(this%m_frootc_xfer_to_fire_patch               (begp:endp)) ; this%m_frootc_xfer_to_fire_patch               (:) = nan
    allocate(this%m_livecrootc_to_fire_patch                (begp:endp)) ; this%m_livecrootc_to_fire_patch                (:) = nan
    allocate(this%m_livecrootc_storage_to_fire_patch        (begp:endp)) ; this%m_livecrootc_storage_to_fire_patch        (:) = nan
    allocate(this%m_livecrootc_xfer_to_fire_patch           (begp:endp)) ; this%m_livecrootc_xfer_to_fire_patch           (:) = nan
    allocate(this%m_deadcrootc_to_fire_patch                (begp:endp)) ; this%m_deadcrootc_to_fire_patch                (:) = nan
    allocate(this%m_deadcrootc_storage_to_fire_patch        (begp:endp)) ; this%m_deadcrootc_storage_to_fire_patch        (:) = nan
    allocate(this%m_deadcrootc_xfer_to_fire_patch           (begp:endp)) ; this%m_deadcrootc_xfer_to_fire_patch           (:) = nan
    allocate(this%m_gresp_storage_to_fire_patch             (begp:endp)) ; this%m_gresp_storage_to_fire_patch             (:) = nan
    allocate(this%m_gresp_xfer_to_fire_patch                (begp:endp)) ; this%m_gresp_xfer_to_fire_patch                (:) = nan
    allocate(this%m_leafc_to_litter_fire_patch              (begp:endp)) ; this%m_leafc_to_litter_fire_patch              (:) = nan
    allocate(this%m_leafc_storage_to_litter_fire_patch      (begp:endp)) ; this%m_leafc_storage_to_litter_fire_patch      (:) = nan
    allocate(this%m_leafc_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_leafc_xfer_to_litter_fire_patch         (:) = nan
    allocate(this%m_livestemc_to_litter_fire_patch          (begp:endp)) ; this%m_livestemc_to_litter_fire_patch          (:) = nan
    allocate(this%m_livestemc_storage_to_litter_fire_patch  (begp:endp)) ; this%m_livestemc_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_livestemc_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_livestemc_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_livestemc_to_deadstemc_fire_patch       (begp:endp)) ; this%m_livestemc_to_deadstemc_fire_patch       (:) = nan
    allocate(this%m_deadstemc_to_litter_fire_patch          (begp:endp)) ; this%m_deadstemc_to_litter_fire_patch          (:) = nan
    allocate(this%m_deadstemc_storage_to_litter_fire_patch  (begp:endp)) ; this%m_deadstemc_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_deadstemc_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_deadstemc_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootc_to_litter_fire_patch             (begp:endp)) ; this%m_frootc_to_litter_fire_patch             (:) = nan
    allocate(this%m_frootc_storage_to_litter_fire_patch     (begp:endp)) ; this%m_frootc_storage_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootc_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_frootc_xfer_to_litter_fire_patch        (:) = nan
    allocate(this%m_livecrootc_to_litter_fire_patch         (begp:endp)) ; this%m_livecrootc_to_litter_fire_patch         (:) = nan
    allocate(this%m_livecrootc_storage_to_litter_fire_patch (begp:endp)) ; this%m_livecrootc_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_livecrootc_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livecrootc_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_livecrootc_to_deadcrootc_fire_patch     (begp:endp)) ; this%m_livecrootc_to_deadcrootc_fire_patch     (:) = nan
    allocate(this%m_deadcrootc_to_litter_fire_patch         (begp:endp)) ; this%m_deadcrootc_to_litter_fire_patch         (:) = nan
    allocate(this%m_deadcrootc_storage_to_litter_fire_patch (begp:endp)) ; this%m_deadcrootc_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_deadcrootc_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadcrootc_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_gresp_storage_to_litter_fire_patch      (begp:endp)) ; this%m_gresp_storage_to_litter_fire_patch      (:) = nan
    allocate(this%m_gresp_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_gresp_xfer_to_litter_fire_patch         (:) = nan
    allocate(this%leafc_xfer_to_leafc_patch                 (begp:endp)) ; this%leafc_xfer_to_leafc_patch                 (:) = nan
    allocate(this%frootc_xfer_to_frootc_patch               (begp:endp)) ; this%frootc_xfer_to_frootc_patch               (:) = nan
    allocate(this%livestemc_xfer_to_livestemc_patch         (begp:endp)) ; this%livestemc_xfer_to_livestemc_patch         (:) = nan
    allocate(this%deadstemc_xfer_to_deadstemc_patch         (begp:endp)) ; this%deadstemc_xfer_to_deadstemc_patch         (:) = nan
    allocate(this%livecrootc_xfer_to_livecrootc_patch       (begp:endp)) ; this%livecrootc_xfer_to_livecrootc_patch       (:) = nan
    allocate(this%deadcrootc_xfer_to_deadcrootc_patch       (begp:endp)) ; this%deadcrootc_xfer_to_deadcrootc_patch       (:) = nan
    allocate(this%leafc_to_litter_patch                     (begp:endp)) ; this%leafc_to_litter_patch                     (:) = nan
    allocate(this%leafc_to_litter_fun_patch                 (begp:endp)) ; this%leafc_to_litter_fun_patch                 (:) = nan
    allocate(this%frootc_to_litter_patch                    (begp:endp)) ; this%frootc_to_litter_patch                    (:) = nan
    allocate(this%cpool_to_resp_patch                       (begp:endp)) ; this%cpool_to_resp_patch                       (:) = nan
    allocate(this%cpool_to_leafc_resp_patch                 (begp:endp)) ; this%cpool_to_leafc_resp_patch                 (:) = nan
    allocate(this%cpool_to_leafc_storage_resp_patch         (begp:endp)) ; this%cpool_to_leafc_storage_resp_patch         (:) = nan
    allocate(this%cpool_to_frootc_resp_patch                (begp:endp)) ; this%cpool_to_frootc_resp_patch                (:) = nan
    allocate(this%cpool_to_frootc_storage_resp_patch        (begp:endp)) ; this%cpool_to_frootc_storage_resp_patch        (:) = nan
    allocate(this%cpool_to_livecrootc_resp_patch            (begp:endp)) ; this%cpool_to_livecrootc_resp_patch            (:) = nan
    allocate(this%cpool_to_livecrootc_storage_resp_patch    (begp:endp)) ; this%cpool_to_livecrootc_storage_resp_patch    (:) = nan
    allocate(this%cpool_to_livestemc_resp_patch             (begp:endp)) ; this%cpool_to_livestemc_resp_patch             (:) = nan
    allocate(this%cpool_to_livestemc_storage_resp_patch     (begp:endp)) ; this%cpool_to_livestemc_storage_resp_patch     (:) = nan
    allocate(this%leaf_mr_patch                             (begp:endp)) ; this%leaf_mr_patch                             (:) = nan
    allocate(this%froot_mr_patch                            (begp:endp)) ; this%froot_mr_patch                            (:) = nan
    allocate(this%livestem_mr_patch                         (begp:endp)) ; this%livestem_mr_patch                         (:) = nan
    allocate(this%livecroot_mr_patch                        (begp:endp)) ; this%livecroot_mr_patch                        (:) = nan
    allocate(this%grain_mr_patch                            (begp:endp)) ; this%grain_mr_patch                            (:) = nan
    allocate(this%leaf_curmr_patch                          (begp:endp)) ; this%leaf_curmr_patch                          (:) = nan
    allocate(this%froot_curmr_patch                         (begp:endp)) ; this%froot_curmr_patch                         (:) = nan
    allocate(this%livestem_curmr_patch                      (begp:endp)) ; this%livestem_curmr_patch                      (:) = nan
    allocate(this%livecroot_curmr_patch                     (begp:endp)) ; this%livecroot_curmr_patch                     (:) = nan
    allocate(this%grain_curmr_patch                         (begp:endp)) ; this%grain_curmr_patch                         (:) = nan
    allocate(this%leaf_xsmr_patch                           (begp:endp)) ; this%leaf_xsmr_patch                           (:) = nan
    allocate(this%froot_xsmr_patch                          (begp:endp)) ; this%froot_xsmr_patch                          (:) = nan
    allocate(this%livestem_xsmr_patch                       (begp:endp)) ; this%livestem_xsmr_patch                       (:) = nan
    allocate(this%livecroot_xsmr_patch                      (begp:endp)) ; this%livecroot_xsmr_patch                      (:) = nan
    allocate(this%grain_xsmr_patch                          (begp:endp)) ; this%grain_xsmr_patch                          (:) = nan
    allocate(this%psnsun_to_cpool_patch                     (begp:endp)) ; this%psnsun_to_cpool_patch                     (:) = nan
    allocate(this%psnshade_to_cpool_patch                   (begp:endp)) ; this%psnshade_to_cpool_patch                   (:) = nan
    allocate(this%cpool_to_xsmrpool_patch                   (begp:endp)) ; this%cpool_to_xsmrpool_patch                   (:) = nan
    allocate(this%cpool_to_leafc_patch                      (begp:endp)) ; this%cpool_to_leafc_patch                      (:) = nan
    allocate(this%cpool_to_leafc_storage_patch              (begp:endp)) ; this%cpool_to_leafc_storage_patch              (:) = nan
    allocate(this%cpool_to_frootc_patch                     (begp:endp)) ; this%cpool_to_frootc_patch                     (:) = nan
    allocate(this%cpool_to_frootc_storage_patch             (begp:endp)) ; this%cpool_to_frootc_storage_patch             (:) = nan
    allocate(this%cpool_to_livestemc_patch                  (begp:endp)) ; this%cpool_to_livestemc_patch                  (:) = nan
    allocate(this%cpool_to_livestemc_storage_patch          (begp:endp)) ; this%cpool_to_livestemc_storage_patch          (:) = nan
    allocate(this%cpool_to_deadstemc_patch                  (begp:endp)) ; this%cpool_to_deadstemc_patch                  (:) = nan
    allocate(this%cpool_to_deadstemc_storage_patch          (begp:endp)) ; this%cpool_to_deadstemc_storage_patch          (:) = nan
    allocate(this%cpool_to_livecrootc_patch                 (begp:endp)) ; this%cpool_to_livecrootc_patch                 (:) = nan
    allocate(this%cpool_to_livecrootc_storage_patch         (begp:endp)) ; this%cpool_to_livecrootc_storage_patch         (:) = nan
    allocate(this%cpool_to_deadcrootc_patch                 (begp:endp)) ; this%cpool_to_deadcrootc_patch                 (:) = nan
    allocate(this%cpool_to_deadcrootc_storage_patch         (begp:endp)) ; this%cpool_to_deadcrootc_storage_patch         (:) = nan
    allocate(this%cpool_to_gresp_storage_patch              (begp:endp)) ; this%cpool_to_gresp_storage_patch              (:) = nan
    allocate(this%cpool_leaf_gr_patch                       (begp:endp)) ; this%cpool_leaf_gr_patch                       (:) = nan
    allocate(this%cpool_leaf_storage_gr_patch               (begp:endp)) ; this%cpool_leaf_storage_gr_patch               (:) = nan
    allocate(this%transfer_leaf_gr_patch                    (begp:endp)) ; this%transfer_leaf_gr_patch                    (:) = nan
    allocate(this%cpool_froot_gr_patch                      (begp:endp)) ; this%cpool_froot_gr_patch                      (:) = nan
    allocate(this%cpool_froot_storage_gr_patch              (begp:endp)) ; this%cpool_froot_storage_gr_patch              (:) = nan
    allocate(this%transfer_froot_gr_patch                   (begp:endp)) ; this%transfer_froot_gr_patch                   (:) = nan
    allocate(this%cpool_livestem_gr_patch                   (begp:endp)) ; this%cpool_livestem_gr_patch                   (:) = nan
    allocate(this%cpool_livestem_storage_gr_patch           (begp:endp)) ; this%cpool_livestem_storage_gr_patch           (:) = nan
    allocate(this%transfer_livestem_gr_patch                (begp:endp)) ; this%transfer_livestem_gr_patch                (:) = nan
    allocate(this%cpool_deadstem_gr_patch                   (begp:endp)) ; this%cpool_deadstem_gr_patch                   (:) = nan
    allocate(this%cpool_deadstem_storage_gr_patch           (begp:endp)) ; this%cpool_deadstem_storage_gr_patch           (:) = nan
    allocate(this%transfer_deadstem_gr_patch                (begp:endp)) ; this%transfer_deadstem_gr_patch                (:) = nan
    allocate(this%cpool_livecroot_gr_patch                  (begp:endp)) ; this%cpool_livecroot_gr_patch                  (:) = nan
    allocate(this%cpool_livecroot_storage_gr_patch          (begp:endp)) ; this%cpool_livecroot_storage_gr_patch          (:) = nan
    allocate(this%transfer_livecroot_gr_patch               (begp:endp)) ; this%transfer_livecroot_gr_patch               (:) = nan
    allocate(this%cpool_deadcroot_gr_patch                  (begp:endp)) ; this%cpool_deadcroot_gr_patch                  (:) = nan
    allocate(this%cpool_deadcroot_storage_gr_patch          (begp:endp)) ; this%cpool_deadcroot_storage_gr_patch          (:) = nan
    allocate(this%transfer_deadcroot_gr_patch               (begp:endp)) ; this%transfer_deadcroot_gr_patch               (:) = nan
    allocate(this%leafc_storage_to_xfer_patch               (begp:endp)) ; this%leafc_storage_to_xfer_patch               (:) = nan
    allocate(this%frootc_storage_to_xfer_patch              (begp:endp)) ; this%frootc_storage_to_xfer_patch              (:) = nan
    allocate(this%livestemc_storage_to_xfer_patch           (begp:endp)) ; this%livestemc_storage_to_xfer_patch           (:) = nan
    allocate(this%deadstemc_storage_to_xfer_patch           (begp:endp)) ; this%deadstemc_storage_to_xfer_patch           (:) = nan
    allocate(this%livecrootc_storage_to_xfer_patch          (begp:endp)) ; this%livecrootc_storage_to_xfer_patch          (:) = nan
    allocate(this%deadcrootc_storage_to_xfer_patch          (begp:endp)) ; this%deadcrootc_storage_to_xfer_patch          (:) = nan
    allocate(this%gresp_storage_to_xfer_patch               (begp:endp)) ; this%gresp_storage_to_xfer_patch               (:) = nan
    allocate(this%livestemc_to_deadstemc_patch              (begp:endp)) ; this%livestemc_to_deadstemc_patch              (:) = nan
    allocate(this%livecrootc_to_deadcrootc_patch            (begp:endp)) ; this%livecrootc_to_deadcrootc_patch            (:) = nan
    allocate(this%current_gr_patch                          (begp:endp)) ; this%current_gr_patch                          (:) = nan
    allocate(this%transfer_gr_patch                         (begp:endp)) ; this%transfer_gr_patch                         (:) = nan
    allocate(this%storage_gr_patch                          (begp:endp)) ; this%storage_gr_patch                          (:) = nan
    allocate(this%plant_calloc_patch                        (begp:endp)) ; this%plant_calloc_patch                        (:) = nan
    allocate(this%excess_cflux_patch                        (begp:endp)) ; this%excess_cflux_patch                        (:) = nan
    allocate(this%prev_leafc_to_litter_patch                (begp:endp)) ; this%prev_leafc_to_litter_patch                (:) = nan
    allocate(this%prev_frootc_to_litter_patch               (begp:endp)) ; this%prev_frootc_to_litter_patch               (:) = nan
    allocate(this%gpp_before_downreg_patch                  (begp:endp)) ; this%gpp_before_downreg_patch                  (:) = nan
    allocate(this%availc_patch                              (begp:endp)) ; this%availc_patch                              (:) = nan
    allocate(this%xsmrpool_recover_patch                    (begp:endp)) ; this%xsmrpool_recover_patch                    (:) = nan
    allocate(this%xsmrpool_c13ratio_patch                   (begp:endp)) ; this%xsmrpool_c13ratio_patch                   (:) = nan

    allocate(this%cpool_to_grainc_patch                     (begp:endp)) ; this%cpool_to_grainc_patch                     (:) = nan
    allocate(this%cpool_to_grainc_storage_patch             (begp:endp)) ; this%cpool_to_grainc_storage_patch             (:) = nan
    allocate(this%livestemc_to_litter_patch                 (begp:endp)) ; this%livestemc_to_litter_patch                 (:) = nan
    allocate(this%grainc_to_food_patch                      (begp:endp)) ; this%grainc_to_food_patch                      (:) = nan
    allocate(this%leafc_to_biofuelc_patch                   (begp:endp)) ; this%leafc_to_biofuelc_patch                   (:) = nan
    allocate(this%livestemc_to_biofuelc_patch               (begp:endp)) ; this%livestemc_to_biofuelc_patch               (:) = nan
    allocate(this%grainc_to_seed_patch                      (begp:endp)) ; this%grainc_to_seed_patch                      (:) = nan
    allocate(this%grainc_xfer_to_grainc_patch               (begp:endp)) ; this%grainc_xfer_to_grainc_patch               (:) = nan
    allocate(this%cpool_grain_gr_patch                      (begp:endp)) ; this%cpool_grain_gr_patch                      (:) = nan
    allocate(this%cpool_grain_storage_gr_patch              (begp:endp)) ; this%cpool_grain_storage_gr_patch              (:) = nan
    allocate(this%transfer_grain_gr_patch                   (begp:endp)) ; this%transfer_grain_gr_patch                   (:) = nan
    allocate(this%xsmrpool_to_atm_patch                     (begp:endp)) ; this%xsmrpool_to_atm_patch                     (:) = 0.0_r8
    allocate(this%xsmrpool_to_atm_col                       (begc:endc)) ; this%xsmrpool_to_atm_col                       (:) = 0.0_r8
    allocate(this%xsmrpool_to_atm_grc                       (begg:endg)) ; this%xsmrpool_to_atm_grc                       (:) = 0.0_r8
    allocate(this%grainc_storage_to_xfer_patch              (begp:endp)) ; this%grainc_storage_to_xfer_patch              (:) = nan
    allocate(this%frootc_alloc_patch                        (begp:endp)) ; this%frootc_alloc_patch                        (:) = nan
    allocate(this%frootc_loss_patch                         (begp:endp)) ; this%frootc_loss_patch                         (:) = nan
    allocate(this%leafc_alloc_patch                         (begp:endp)) ; this%leafc_alloc_patch                         (:) = nan
    allocate(this%leafc_loss_patch                          (begp:endp)) ; this%leafc_loss_patch                          (:) = nan
    allocate(this%woodc_alloc_patch                         (begp:endp)) ; this%woodc_alloc_patch                         (:) = nan
    allocate(this%woodc_loss_patch                          (begp:endp)) ; this%woodc_loss_patch                          (:) = nan

    allocate(this%phenology_c_to_litr_met_c_col     (begc:endc,1:nlevdecomp_full));
    this%phenology_c_to_litr_met_c_col (:,:)=nan

    allocate(this%phenology_c_to_litr_cel_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_cel_c_col (:,:)=nan
    allocate(this%phenology_c_to_litr_lig_c_col     (begc:endc,1:nlevdecomp_full)); this%phenology_c_to_litr_lig_c_col (:,:)=nan

    allocate(this%gap_mortality_c_to_litr_met_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_met_c_col(:,:)=nan
    allocate(this%gap_mortality_c_to_litr_cel_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_cel_c_col(:,:)=nan
    allocate(this%gap_mortality_c_to_litr_lig_c_col (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_litr_lig_c_col(:,:)=nan

    allocate(this%gap_mortality_c_to_cwdc_col       (begc:endc,1:nlevdecomp_full)); this%gap_mortality_c_to_cwdc_col  (:,:)=nan
    allocate(this%fire_mortality_c_to_cwdc_col      (begc:endc,1:nlevdecomp_full)); this%fire_mortality_c_to_cwdc_col (:,:)=nan
    allocate(this%m_c_to_litr_met_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_met_fire_col     (:,:)=nan
    allocate(this%m_c_to_litr_cel_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_cel_fire_col     (:,:)=nan
    allocate(this%m_c_to_litr_lig_fire_col          (begc:endc,1:nlevdecomp_full)); this%m_c_to_litr_lig_fire_col     (:,:)=nan
    allocate(this%harvest_c_to_litr_met_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_met_c_col  (:,:)=nan
    allocate(this%harvest_c_to_litr_cel_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_cel_c_col  (:,:)=nan
    allocate(this%harvest_c_to_litr_lig_c_col       (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_litr_lig_c_col  (:,:)=nan
    allocate(this%harvest_c_to_cwdc_col             (begc:endc,1:nlevdecomp_full)); this%harvest_c_to_cwdc_col        (:,:)=nan

    allocate(this%dwt_slash_cflux_patch             (begp:endp))                  ; this%dwt_slash_cflux_patch        (:) =nan
    allocate(this%dwt_slash_cflux_grc               (begg:endg))                  ; this%dwt_slash_cflux_grc          (:) =nan
    allocate(this%dwt_frootc_to_litr_met_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_met_c_col (:,:)=nan
    allocate(this%dwt_frootc_to_litr_cel_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_cel_c_col (:,:)=nan
    allocate(this%dwt_frootc_to_litr_lig_c_col      (begc:endc,1:nlevdecomp_full)); this%dwt_frootc_to_litr_lig_c_col (:,:)=nan
    allocate(this%dwt_livecrootc_to_cwdc_col        (begc:endc,1:nlevdecomp_full)); this%dwt_livecrootc_to_cwdc_col   (:,:)=nan
    allocate(this%dwt_deadcrootc_to_cwdc_col        (begc:endc,1:nlevdecomp_full)); this%dwt_deadcrootc_to_cwdc_col   (:,:)=nan

    allocate(this%dwt_seedc_to_leaf_patch           (begp:endp))                  ; this%dwt_seedc_to_leaf_patch   (:)  =nan
    allocate(this%dwt_seedc_to_leaf_grc             (begg:endg))                  ; this%dwt_seedc_to_leaf_grc     (:)  =nan
    allocate(this%dwt_seedc_to_deadstem_patch       (begp:endp))                  ; this%dwt_seedc_to_deadstem_patch(:)  =nan
    allocate(this%dwt_seedc_to_deadstem_grc         (begg:endg))                  ; this%dwt_seedc_to_deadstem_grc (:)  =nan
    allocate(this%dwt_conv_cflux_patch              (begp:endp))                  ; this%dwt_conv_cflux_patch      (:)  =nan
    allocate(this%dwt_conv_cflux_grc                (begg:endg))                  ; this%dwt_conv_cflux_grc        (:)  =nan
    allocate(this%dwt_conv_cflux_dribbled_grc       (begg:endg))                  ; this%dwt_conv_cflux_dribbled_grc(:)  =nan
    allocate(this%dwt_wood_productc_gain_patch      (begp:endp))                  ; this%dwt_wood_productc_gain_patch(:)  =nan
    allocate(this%dwt_crop_productc_gain_patch      (begp:endp))                  ; this%dwt_crop_productc_gain_patch(:) =nan

    allocate(this%crop_seedc_to_leaf_patch          (begp:endp))                  ; this%crop_seedc_to_leaf_patch  (:)  =nan

    allocate(this%cwdc_hr_col                       (begc:endc))                  ; this%cwdc_hr_col               (:)  =nan
    allocate(this%cwdc_loss_col                     (begc:endc))                  ; this%cwdc_loss_col             (:)  =nan
    allocate(this%litterc_loss_col                  (begc:endc))                  ; this%litterc_loss_col          (:)  =nan

    allocate(this%grainc_to_cropprodc_patch(begp:endp))
    this%grainc_to_cropprodc_patch(:) = spval

    allocate(this%grainc_to_cropprodc_col(begc:endc))
    this%grainc_to_cropprodc_col(:) = nan

    allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%m_decomp_cpools_to_fire_vr_col(:,:,:)= nan

    allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))
    this%m_decomp_cpools_to_fire_col(:,:)= nan

    allocate(this%m_decomp_cpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%m_decomp_cpools_to_fire_vr_col(:,:,:)= nan

    allocate(this%m_decomp_cpools_to_fire_col(begc:endc,1:ndecomp_pools))
    this%m_decomp_cpools_to_fire_col(:,:)= nan

    allocate(this%rr_patch                (begp:endp)) ; this%rr_patch                (:) = spval
    allocate(this%mr_patch                (begp:endp)) ; this%mr_patch                (:) = nan
    allocate(this%gr_patch                (begp:endp)) ; this%gr_patch                (:) = nan
    allocate(this%ar_patch                (begp:endp)) ; this%ar_patch                (:) = spval
    allocate(this%npp_patch               (begp:endp)) ; this%npp_patch               (:) = spval
    allocate(this%agnpp_patch             (begp:endp)) ; this%agnpp_patch             (:) = nan
    allocate(this%bgnpp_patch             (begp:endp)) ; this%bgnpp_patch             (:) = nan
    allocate(this%litfall_patch           (begp:endp)) ; this%litfall_patch           (:) = nan
    allocate(this%wood_harvestc_patch     (begp:endp)) ; this%wood_harvestc_patch     (:) = nan
    allocate(this%slash_harvestc_patch    (begp:endp)) ; this%slash_harvestc_patch    (:) = nan
    allocate(this%cinputs_patch           (begp:endp)) ; this%cinputs_patch           (:) = nan
    allocate(this%coutputs_patch          (begp:endp)) ; this%coutputs_patch          (:) = nan
    allocate(this%gpp_patch               (begp:endp)) ; this%gpp_patch               (:) = spval
    allocate(this%fire_closs_patch        (begp:endp)) ; this%fire_closs_patch        (:) = spval
    allocate(this%sr_col                  (begc:endc)) ; this%sr_col                  (:) = nan
    allocate(this%er_col                  (begc:endc)) ; this%er_col                  (:) = nan
    allocate(this%litfire_col             (begc:endc)) ; this%litfire_col             (:) = nan
    allocate(this%somfire_col             (begc:endc)) ; this%somfire_col             (:) = nan
    allocate(this%totfire_col             (begc:endc)) ; this%totfire_col             (:) = nan
    allocate(this%rr_col                  (begc:endc)) ; this%rr_col                  (:) = nan
    allocate(this%ar_col                  (begc:endc)) ; this%ar_col                  (:) = nan
    allocate(this%gpp_col                 (begc:endc)) ; this%gpp_col                 (:) = nan
    allocate(this%npp_col                 (begc:endc)) ; this%npp_col                 (:) = nan
    allocate(this%fire_closs_p2c_col      (begc:endc)) ; this%fire_closs_p2c_col      (:) = nan
    allocate(this%fire_closs_col          (begc:endc)) ; this%fire_closs_col          (:) = spval
    allocate(this%wood_harvestc_col       (begc:endc)) ; this%wood_harvestc_col       (:) = nan
    allocate(this%hrv_xsmrpool_to_atm_col (begc:endc)) ; this%hrv_xsmrpool_to_atm_col (:) = 0.0_r8
    allocate(this%tempsum_npp_patch       (begp:endp)) ; this%tempsum_npp_patch       (:) = nan
    allocate(this%annsum_npp_patch        (begp:endp)) ; this%annsum_npp_patch        (:) = spval
    allocate(this%tempsum_litfall_patch   (begp:endp)) ; this%tempsum_litfall_patch   (:) = nan
    allocate(this%annsum_litfall_patch    (begp:endp)) ; this%annsum_litfall_patch    (:) = nan
    allocate(this%annsum_npp_col          (begc:endc)) ; this%annsum_npp_col          (:) = nan
    allocate(this%lag_npp_col             (begc:endc)) ; this%lag_npp_col             (:) = spval


    allocate(this%nep_col                 (begc:endc)) ; this%nep_col                 (:) = spval
    allocate(this%nbp_grc                 (begg:endg)) ; this%nbp_grc                 (:) = nan
    allocate(this%nee_grc                 (begg:endg)) ; this%nee_grc                 (:) = nan
    allocate(this%landuseflux_grc         (begg:endg)) ; this%landuseflux_grc         (:) = nan
    allocate(this%npp_Nactive_patch       (begp:endp)) ; this%npp_Nactive_patch       (:) = nan
    allocate(this%npp_burnedoff_patch     (begp:endp)) ; this%npp_burnedoff_patch     (:) = nan
    allocate(this%npp_Nnonmyc_patch       (begp:endp)) ; this%npp_Nnonmyc_patch       (:) = nan
    allocate(this%npp_Nam_patch           (begp:endp)) ; this%npp_Nam_patch           (:) = nan
    allocate(this%npp_Necm_patch          (begp:endp)) ; this%npp_Necm_patch          (:) = nan
    allocate(this%npp_Nactive_no3_patch   (begp:endp)) ; this%npp_Nactive_no3_patch   (:) = nan
    allocate(this%npp_Nactive_nh4_patch   (begp:endp)) ; this%npp_Nactive_nh4_patch   (:) = nan
    allocate(this%npp_Nnonmyc_no3_patch   (begp:endp)) ; this%npp_Nnonmyc_no3_patch   (:) = nan
    allocate(this%npp_Nnonmyc_nh4_patch   (begp:endp)) ; this%npp_Nnonmyc_nh4_patch   (:) = nan
    allocate(this%npp_Nam_no3_patch       (begp:endp)) ; this%npp_Nam_no3_patch       (:) = nan
    allocate(this%npp_Nam_nh4_patch       (begp:endp)) ; this%npp_Nam_nh4_patch       (:) = nan
    allocate(this%npp_Necm_no3_patch      (begp:endp)) ; this%npp_Necm_no3_patch      (:) = nan
    allocate(this%npp_Necm_nh4_patch      (begp:endp)) ; this%npp_Necm_nh4_patch      (:) = nan
    allocate(this%npp_Nfix_patch          (begp:endp)) ; this%npp_Nfix_patch          (:) = nan
    allocate(this%npp_Nretrans_patch      (begp:endp)) ; this%npp_Nretrans_patch      (:) = nan
    allocate(this%npp_Nuptake_patch       (begp:endp)) ; this%npp_Nuptake_patch       (:) = nan
    allocate(this%npp_growth_patch        (begp:endp)) ; this%npp_growth_patch       (:) = nan
    allocate(this%leafc_change_patch      (begp:endp)) ; this%leafc_change_patch      (:) = nan
    allocate(this%soilc_change_patch      (begp:endp)) ; this%soilc_change_patch      (:) = nan


 ! initialize variables from restart file or set to cold start value

 this%dwt_conv_cflux_dribbled_grc(begg:endg) = 0._r8
 this%dwt_conv_cflux_grc(begg:endg) = 0._r8


 n = 0
 np = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          n = n + 1

          this%annsum_npp_col (n) = cncol(nc,nz, 33)


          do nl = 1, nlevdecomp_full
             this%dwt_frootc_to_litr_met_c_col(n,nl) = 0._r8
             this%dwt_frootc_to_litr_cel_c_col(n,nl) = 0._r8
             this%dwt_frootc_to_litr_lig_c_col(n,nl) = 0._r8
             this%dwt_livecrootc_to_cwdc_col(n,nl)   = 0._r8
             this%dwt_deadcrootc_to_cwdc_col(n,nl)   = 0._r8
          end do


          do p = 0,numpft  ! PFT index loop
             np = np + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then

                   this%gpp_patch(p)                   = 0._r8
                   this%excess_cflux_patch(p)          = 0._r8
                   this%leafc_to_litter_fun_patch(p)   = 0._r8
                   this%plant_calloc_patch(p)          = 0._r8

                  ! "old" variables: CNCLM45 and before
                  this%annsum_npp_patch            (np) = cnpft(nc,nz,nv, 26)
                  this%prev_frootc_to_litter_patch (np) = cnpft(nc,nz,nv, 41)
                  this%prev_leafc_to_litter_patch  (np) = cnpft(nc,nz,nv, 42)
                  this%tempsum_npp_patch           (np) = cnpft(nc,nz,nv, 45)
                  this%xsmrpool_recover_patch      (np) = cnpft(nc,nz,nv, 47)
                  

                  ! "new" variables: introduced in CNCLM50
                  if (cold_start.eqv..false.) then
                      this%annsum_litfall_patch(np)    = cnpft(nc,nz,nv, 82)
                      this%tempsum_litfall_patch(np)   = cnpft(nc,nz,nv, 83)  
                  elseif (cold_start) then
                      this%annsum_litfall_patch(np)    = 0._r8      
                      this%tempsum_litfall_patch(np)   = 0._r8 
                  else
                     _ASSERT(.FALSE.,'missing CNCLM50_cold_start setting')
                  end if 

                 end if
            end do !nv

            this%excess_cflux_patch(np)          = 0._r8
            this%leafc_to_litter_fun_patch(np)   = 0._r8
            this%plant_calloc_patch(np)          = 0._r8
            this%dwt_wood_productc_gain_patch(np) = 0._r8   ! following CNCLM45 setting
            this%dwt_crop_productc_gain_patch(np) = 0._r8   ! following CNCLM45 setting


       end do ! p
     end do ! nz
  end do ! nc     

    ! Construct restart field names consistently to what is done in SpeciesNonIsotope &
    ! SpeciesIsotope, to aid future migration to that infrastructure
    if (carbon_type == 'c12') then
       carbon_type_suffix = 'c'
    else if (carbon_type == 'c13') then
       carbon_type_suffix = 'c_13'
    else if (carbon_type == 'c14') then
       carbon_type_suffix = 'c_14'
    else
       write(iulog,*) 'CNVegCarbonFluxType InitAllocate: Unknown carbon_type: ', trim(carbon_type)
       call endrun(msg='CNVegCarbonFluxType InitAllocate: Unknown carbon_type: ' // &
            errMsg(sourcefile, __LINE__))
    end if

    if (use_cndv) then
       allows_non_annual_delta = .true.  
    else                                          
       allows_non_annual_delta = .false.          
    end if
    this%dwt_conv_cflux_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &                       
         name = 'dwt_conv_flux_' // carbon_type_suffix, &
         units = 'gC/m^2', &                      
         allows_non_annual_delta = allows_non_annual_delta)
    this%hrv_xsmrpool_to_atm_dribbler = annual_flux_dribbler_gridcell( &
         bounds = bounds, &              
         name = 'hrv_xsmrpool_to_atm_' // carbon_type_suffix, &
         units = 'gC/m^2', &
         allows_non_annual_delta = .false.)

  end subroutine Init

!-----------------------------------------

  subroutine SetValues ( this, nvegcpool, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set carbon state fluxes
    !
    ! !ARGUMENTS:
    class (cnveg_carbonflux_type) :: this
    integer , intent(in) :: num_patch,nvegcpool
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k,l    ! indices
    !------------------------------------------------------------------------


    do fi = 1,num_patch
       i = filter_patch(fi)

       this%m_leafc_to_litter_patch(i)                   = value_patch
       this%m_frootc_to_litter_patch(i)                  = value_patch
       this%m_leafc_storage_to_litter_patch(i)           = value_patch
       this%m_frootc_storage_to_litter_patch(i)          = value_patch
       this%m_livestemc_storage_to_litter_patch(i)       = value_patch
       this%m_deadstemc_storage_to_litter_patch(i)       = value_patch
       this%m_livecrootc_storage_to_litter_patch(i)      = value_patch
       this%m_deadcrootc_storage_to_litter_patch(i)      = value_patch
       this%m_leafc_xfer_to_litter_patch(i)              = value_patch
       this%m_frootc_xfer_to_litter_patch(i)             = value_patch
       this%m_livestemc_xfer_to_litter_patch(i)          = value_patch
       this%m_deadstemc_xfer_to_litter_patch(i)          = value_patch
       this%m_livecrootc_xfer_to_litter_patch(i)         = value_patch
       this%m_deadcrootc_xfer_to_litter_patch(i)         = value_patch
       this%m_livestemc_to_litter_patch(i)               = value_patch
       this%m_deadstemc_to_litter_patch(i)               = value_patch
       this%m_livecrootc_to_litter_patch(i)              = value_patch
       this%m_deadcrootc_to_litter_patch(i)              = value_patch
       this%m_gresp_storage_to_litter_patch(i)           = value_patch
       this%m_gresp_xfer_to_litter_patch(i)              = value_patch
       this%hrv_leafc_to_litter_patch(i)                 = value_patch
       this%hrv_leafc_storage_to_litter_patch(i)         = value_patch
       this%hrv_leafc_xfer_to_litter_patch(i)            = value_patch
       this%hrv_frootc_to_litter_patch(i)                = value_patch
       this%hrv_frootc_storage_to_litter_patch(i)        = value_patch
       this%hrv_frootc_xfer_to_litter_patch(i)           = value_patch
       this%hrv_livestemc_to_litter_patch(i)             = value_patch
       this%hrv_livestemc_storage_to_litter_patch(i)     = value_patch
       this%hrv_livestemc_xfer_to_litter_patch(i)        = value_patch
       this%hrv_deadstemc_storage_to_litter_patch(i)     = value_patch
       this%hrv_deadstemc_xfer_to_litter_patch(i)        = value_patch
       this%hrv_livecrootc_to_litter_patch(i)            = value_patch
       this%hrv_livecrootc_storage_to_litter_patch(i)    = value_patch
       this%hrv_livecrootc_xfer_to_litter_patch(i)       = value_patch
       this%hrv_deadcrootc_to_litter_patch(i)            = value_patch
       this%hrv_deadcrootc_storage_to_litter_patch(i)    = value_patch
       this%hrv_deadcrootc_xfer_to_litter_patch(i)       = value_patch
       this%hrv_gresp_storage_to_litter_patch(i)         = value_patch
       this%hrv_gresp_xfer_to_litter_patch(i)            = value_patch
       this%hrv_xsmrpool_to_atm_patch(i)                 = value_patch


       this%m_leafc_to_fire_patch(i)                     = value_patch
       this%m_leafc_storage_to_fire_patch(i)             = value_patch
       this%m_leafc_xfer_to_fire_patch(i)                = value_patch
       this%m_livestemc_to_fire_patch(i)                 = value_patch
       this%m_livestemc_storage_to_fire_patch(i)         = value_patch
       this%m_livestemc_xfer_to_fire_patch(i)            = value_patch
       this%m_deadstemc_to_fire_patch(i)                 = value_patch
       this%m_deadstemc_storage_to_fire_patch(i)         = value_patch
       this%m_deadstemc_xfer_to_fire_patch(i)            = value_patch
       this%m_frootc_to_fire_patch(i)                    = value_patch
       this%m_frootc_storage_to_fire_patch(i)            = value_patch
       this%m_frootc_xfer_to_fire_patch(i)               = value_patch
       this%m_livecrootc_to_fire_patch(i)                = value_patch
       this%m_livecrootc_storage_to_fire_patch(i)        = value_patch
       this%m_livecrootc_xfer_to_fire_patch(i)           = value_patch
       this%m_deadcrootc_to_fire_patch(i)                = value_patch
       this%m_deadcrootc_storage_to_fire_patch(i)        = value_patch
       this%m_deadcrootc_xfer_to_fire_patch(i)           = value_patch
       this%m_gresp_storage_to_fire_patch(i)             = value_patch
       this%m_gresp_xfer_to_fire_patch(i)                = value_patch

       this%m_leafc_to_litter_fire_patch(i)              = value_patch
       this%m_leafc_storage_to_litter_fire_patch(i)      = value_patch
       this%m_leafc_xfer_to_litter_fire_patch(i)         = value_patch
       this%m_livestemc_to_litter_fire_patch(i)          = value_patch
       this%m_livestemc_storage_to_litter_fire_patch(i)  = value_patch
       this%m_livestemc_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_livestemc_to_deadstemc_fire_patch(i)       = value_patch
       this%m_deadstemc_to_litter_fire_patch(i)          = value_patch
       this%m_deadstemc_storage_to_litter_fire_patch(i)  = value_patch
       this%m_deadstemc_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_frootc_to_litter_fire_patch(i)             = value_patch
       this%m_frootc_storage_to_litter_fire_patch(i)     = value_patch
       this%m_frootc_xfer_to_litter_fire_patch(i)        = value_patch
       this%m_livecrootc_to_litter_fire_patch(i)         = value_patch
       this%m_livecrootc_storage_to_litter_fire_patch(i) = value_patch
       this%m_livecrootc_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_livecrootc_to_deadcrootc_fire_patch(i)     = value_patch
       this%m_deadcrootc_to_litter_fire_patch(i)         = value_patch
       this%m_deadcrootc_storage_to_litter_fire_patch(i) = value_patch
       this%m_deadcrootc_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_gresp_storage_to_litter_fire_patch(i)      = value_patch
       this%m_gresp_xfer_to_litter_fire_patch(i)         = value_patch

       this%leafc_xfer_to_leafc_patch(i)                 = value_patch
       this%frootc_xfer_to_frootc_patch(i)               = value_patch
       this%livestemc_xfer_to_livestemc_patch(i)         = value_patch
       this%deadstemc_xfer_to_deadstemc_patch(i)         = value_patch
       this%livecrootc_xfer_to_livecrootc_patch(i)       = value_patch
       this%deadcrootc_xfer_to_deadcrootc_patch(i)       = value_patch
       this%leafc_to_litter_patch(i)                     = value_patch
       this%frootc_to_litter_patch(i)                    = value_patch
       this%cpool_to_resp_patch(i)                       = value_patch
       this%cpool_to_leafc_resp_patch(i)                 = value_patch
       this%cpool_to_leafc_storage_resp_patch(i)         = value_patch
       this%cpool_to_frootc_resp_patch(i)                = value_patch
       this%cpool_to_frootc_storage_resp_patch(i)        = value_patch
       this%cpool_to_livecrootc_resp_patch(i)            = value_patch
       this%cpool_to_livecrootc_storage_resp_patch(i)    = value_patch
       this%cpool_to_livestemc_resp_patch(i)             = value_patch
       this%cpool_to_livestemc_storage_resp_patch(i)     = value_patch
       this%leaf_mr_patch(i)                             = value_patch
       this%froot_mr_patch(i)                            = value_patch
       this%livestem_mr_patch(i)                         = value_patch
       this%livecroot_mr_patch(i)                        = value_patch
       this%grain_mr_patch(i)                            = value_patch
       this%leaf_curmr_patch(i)                          = value_patch
       this%froot_curmr_patch(i)                         = value_patch
       this%livestem_curmr_patch(i)                      = value_patch
       this%livecroot_curmr_patch(i)                     = value_patch
       this%grain_curmr_patch(i)                         = value_patch
       this%leaf_xsmr_patch(i)                           = value_patch
       this%froot_xsmr_patch(i)                          = value_patch
       this%livestem_xsmr_patch(i)                       = value_patch
       this%livecroot_xsmr_patch(i)                      = value_patch
       this%grain_xsmr_patch(i)                          = value_patch
       this%psnsun_to_cpool_patch(i)                     = value_patch
       this%psnshade_to_cpool_patch(i)                   = value_patch
       this%cpool_to_xsmrpool_patch(i)                   = value_patch
       this%cpool_to_leafc_patch(i)                      = value_patch
       this%cpool_to_leafc_storage_patch(i)              = value_patch
       this%cpool_to_frootc_patch(i)                     = value_patch
       this%cpool_to_frootc_storage_patch(i)             = value_patch
       this%cpool_to_livestemc_patch(i)                  = value_patch
       this%cpool_to_livestemc_storage_patch(i)          = value_patch
       this%cpool_to_deadstemc_patch(i)                  = value_patch
       this%cpool_to_deadstemc_storage_patch(i)          = value_patch
       this%cpool_to_livecrootc_patch(i)                 = value_patch
       this%cpool_to_livecrootc_storage_patch(i)         = value_patch
       this%cpool_to_deadcrootc_patch(i)                 = value_patch
       this%cpool_to_deadcrootc_storage_patch(i)         = value_patch
       this%cpool_to_gresp_storage_patch(i)              = value_patch
       this%cpool_leaf_gr_patch(i)                       = value_patch
       this%cpool_leaf_storage_gr_patch(i)               = value_patch
       this%transfer_leaf_gr_patch(i)                    = value_patch
       this%cpool_froot_gr_patch(i)                      = value_patch
       this%cpool_froot_storage_gr_patch(i)              = value_patch
       this%transfer_froot_gr_patch(i)                   = value_patch
       this%cpool_livestem_gr_patch(i)                   = value_patch
       this%cpool_livestem_storage_gr_patch(i)           = value_patch
       this%transfer_livestem_gr_patch(i)                = value_patch
       this%cpool_deadstem_gr_patch(i)                   = value_patch
       this%cpool_deadstem_storage_gr_patch(i)           = value_patch
       this%transfer_deadstem_gr_patch(i)                = value_patch
       this%cpool_livecroot_gr_patch(i)                  = value_patch
       this%cpool_livecroot_storage_gr_patch(i)          = value_patch
       this%transfer_livecroot_gr_patch(i)               = value_patch
       this%cpool_deadcroot_gr_patch(i)                  = value_patch
       this%cpool_deadcroot_storage_gr_patch(i)          = value_patch
       this%transfer_deadcroot_gr_patch(i)               = value_patch
       this%leafc_storage_to_xfer_patch(i)               = value_patch
       this%frootc_storage_to_xfer_patch(i)              = value_patch
       this%livestemc_storage_to_xfer_patch(i)           = value_patch
       this%deadstemc_storage_to_xfer_patch(i)           = value_patch
       this%livecrootc_storage_to_xfer_patch(i)          = value_patch
       this%deadcrootc_storage_to_xfer_patch(i)          = value_patch
       this%gresp_storage_to_xfer_patch(i)               = value_patch
       this%livestemc_to_deadstemc_patch(i)              = value_patch
       this%livecrootc_to_deadcrootc_patch(i)            = value_patch

       this%current_gr_patch(i)                          = value_patch
       this%transfer_gr_patch(i)                         = value_patch
       this%storage_gr_patch(i)                          = value_patch
       this%frootc_alloc_patch(i)                        = value_patch
       this%frootc_loss_patch(i)                         = value_patch
       this%leafc_alloc_patch(i)                         = value_patch
       this%leafc_loss_patch(i)                          = value_patch
       this%woodc_alloc_patch(i)                         = value_patch
       this%woodc_loss_patch(i)                          = value_patch

       this%crop_seedc_to_leaf_patch(i)                  = value_patch
       this%grainc_to_cropprodc_patch(i)                 = value_patch
!   Matrix
       if(use_matrixcn)then
          this%matrix_Cinput_patch(i)                       = value_patch
          this%matrix_C13input_patch(i)                       = value_patch
          this%matrix_C14input_patch(i)                       = value_patch
       end if
    end do
    if(use_matrixcn)then
       do j = 1, nvegcpool
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_alloc_patch(i,j)       = value_patch
             this%matrix_phturnover_patch (i,j) = value_patch
             this%matrix_gmturnover_patch (i,j) = value_patch
             this%matrix_fiturnover_patch (i,j) = value_patch
          end do
       end do

       do j = 1, ncphtrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_phtransfer_patch (i,j) = value_patch
          end do
       end do

       do j = 1, ncgmtrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_gmtransfer_patch (i,j) = value_patch
          end do
       end do

       do j = 1, ncfitrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_fitransfer_patch (i,j) = value_patch
          end do
       end do
    end if


    if ( use_crop )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%xsmrpool_to_atm_patch(i)         = value_patch
          this%livestemc_to_litter_patch(i)     = value_patch
          this%grainc_to_food_patch(i)          = value_patch

          this%leafc_to_biofuelc_patch(i)       = value_patch
          this%livestemc_to_biofuelc_patch(i)   = value_patch

          this%grainc_to_seed_patch(i)          = value_patch
          this%grainc_xfer_to_grainc_patch(i)   = value_patch
          this%cpool_to_grainc_patch(i)         = value_patch
          this%cpool_to_grainc_storage_patch(i) = value_patch
          this%cpool_grain_gr_patch(i)          = value_patch
          this%cpool_grain_storage_gr_patch(i)  = value_patch
          this%transfer_grain_gr_patch(i)       = value_patch
          this%grainc_storage_to_xfer_patch(i)  = value_patch
       end do
    end if


    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          this%phenology_c_to_litr_met_c_col(i,j)     = value_column
          this%phenology_c_to_litr_cel_c_col(i,j)     = value_column
          this%phenology_c_to_litr_lig_c_col(i,j)     = value_column

          this%gap_mortality_c_to_litr_met_c_col(i,j) = value_column
          this%gap_mortality_c_to_litr_cel_c_col(i,j) = value_column
          this%gap_mortality_c_to_litr_lig_c_col(i,j) = value_column
          this%gap_mortality_c_to_cwdc_col(i,j)       = value_column

          this%fire_mortality_c_to_cwdc_col(i,j)      = value_column
          this%m_c_to_litr_met_fire_col(i,j)          = value_column
          this%m_c_to_litr_cel_fire_col(i,j)          = value_column
          this%m_c_to_litr_lig_fire_col(i,j)          = value_column

          this%harvest_c_to_litr_met_c_col(i,j)       = value_column
          this%harvest_c_to_litr_cel_c_col(i,j)       = value_column
          this%harvest_c_to_litr_lig_c_col(i,j)       = value_column
          this%harvest_c_to_cwdc_col(i,j)             = value_column

       end do
    end do


    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%m_decomp_cpools_to_fire_vr_col(i,j,k) = value_column
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%m_decomp_cpools_to_fire_col(i,k) = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%grainc_to_cropprodc_col(i)       = value_column
       this%cwdc_hr_col(i)                   = value_column
       this%cwdc_loss_col(i)                 = value_column
       this%litterc_loss_col(i)              = value_column

    end do


    do fi = 1,num_patch
       i = filter_patch(fi)

       this%gpp_patch(i)           = value_patch
       this%mr_patch(i)            = value_patch
       this%gr_patch(i)            = value_patch
       this%ar_patch(i)            = value_patch
       this%rr_patch(i)            = value_patch
       this%npp_patch(i)           = value_patch
       this%agnpp_patch(i)         = value_patch
       this%bgnpp_patch(i)         = value_patch
       this%litfall_patch(i)       = value_patch
       this%wood_harvestc_patch(i) = value_patch
       this%slash_harvestc_patch(i) = value_patch
       this%cinputs_patch(i)       = value_patch
       this%coutputs_patch(i)      = value_patch
       this%fire_closs_patch(i)    = value_patch
       this%npp_Nactive_patch(i)     = value_patch
       this%npp_burnedoff_patch(i)     = value_patch
       this%npp_Nnonmyc_patch(i)     = value_patch
       this%npp_Nam_patch(i)         = value_patch
       this%npp_Necm_patch(i)        = value_patch
       this%npp_Nactive_no3_patch(i) = value_patch
       this%npp_Nactive_nh4_patch(i) = value_patch
       this%npp_Nnonmyc_no3_patch(i) = value_patch
       this%npp_Nnonmyc_nh4_patch(i) = value_patch
       this%npp_Nam_no3_patch(i)     = value_patch
       this%npp_Nam_nh4_patch(i)     = value_patch
       this%npp_Necm_no3_patch(i)    = value_patch
       this%npp_Necm_nh4_patch(i)    = value_patch
       this%npp_Nfix_patch(i)        = value_patch
       this%npp_Nretrans_patch(i)    = value_patch
       this%npp_Nuptake_patch(i)     = value_patch
       this%npp_growth_patch(i)      = value_patch
       this%leafc_change_patch(i)    = value_patch
       this%soilc_change_patch(i)    = value_patch
    end do

    do fi = 1,num_column
       i  = filter_column(fi)

       this%sr_col(i)                  = value_column
       this%er_col(i)                  = value_column
       this%litfire_col(i)             = value_column
       this%somfire_col(i)             = value_column
       this%totfire_col(i)             = value_column

       ! Zero p2c column fluxes
       this%rr_col(i)                  = value_column
       this%ar_col(i)                  = value_column
       this%gpp_col(i)                 = value_column
       this%npp_col(i)                 = value_column
       this%fire_closs_col(i)          = value_column
       this%wood_harvestc_col(i)       = value_column
       this%hrv_xsmrpool_to_atm_col(i) = value_column
       this%nep_col(i)                 = value_column
       if ( use_crop )then
          this%xsmrpool_to_atm_col(i)  = value_column
       end if

    end do

  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary_carbonflux(this, &
       bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       isotope, soilbiogeochem_hr_col, soilbiogeochem_lithr_col, &
       soilbiogeochem_decomp_cascade_ctransfer_col, &
       product_closs_grc)
    !
    ! !DESCRIPTION:
    ! Perform patch and column-level carbon summary calculations
    !
    ! !USES:
    use clm_time_manager                   , only: get_step_size_real
    use clm_varcon                         , only: secspday
    use clm_varctl                         , only: nfix_timeconst, carbon_resp_opt
    use subgridAveMod                      , only: p2c, c2g
    use SoilBiogeochemDecompCascadeConType , only: decomp_cascade_con
    use CNSharedParamsMod                  , only: use_fun
    !
    ! !ARGUMENTS:
    class(cnveg_carbonflux_type)   :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    character(len=*)  , intent(in) :: isotope
    real(r8)          , intent(in) :: soilbiogeochem_hr_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_lithr_col(bounds%begc:)
    real(r8)          , intent(in) :: soilbiogeochem_decomp_cascade_ctransfer_col(bounds%begc:,1:)
    real(r8)          , intent(in) :: product_closs_grc(bounds%begg:)
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l,g     ! indices
    integer  :: fp,fc           ! lake filter indices
    real(r8) :: nfixlags, dtime ! temp variables for making lagged npp
    real(r8) :: maxdepth        ! depth to integrate soil variables
    real(r8) :: nep_grc(bounds%begg:bounds%endg)        ! nep_col averaged to gridcell
    real(r8) :: fire_closs_grc(bounds%begg:bounds%endg) ! fire_closs_col averaged to gridcell
    real(r8) :: hrv_xsmrpool_to_atm_grc(bounds%begg:bounds%endg) ! hrv_xsmrpool_to_atm_col averaged to gridcell (gC/m2/s)
    real(r8) :: hrv_xsmrpool_to_atm_delta_grc(bounds%begg:bounds%endg) ! hrv_xsmrpool_to_atm_col averaged to gridcell, expressed as a delta (not a flux) (gC/m2)
    real(r8) :: hrv_xsmrpool_to_atm_dribbled_grc(bounds%begg:bounds%endg) ! hrv_xsmrpool_to_atm, dribbled over the year (gC/m2/s)
    real(r8) :: dwt_conv_cflux_delta_grc(bounds%begg:bounds%endg)    ! dwt_conv_cflux_grc expressed as a total delta (not a flux) (gC/m2)
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL_FL((ubound(product_closs_grc) == (/bounds%endg/)), sourcefile, __LINE__)

    ! calculate patch-level summary carbon fluxes and states

    dtime = get_step_size_real()

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! maintenance respiration (MR)
       if ( trim(isotope) == 'c13' .or. trim(isotope) == 'c14') then
          this%leaf_mr_patch(p)      = this%leaf_curmr_patch(p)      + this%leaf_xsmr_patch(p)
          this%froot_mr_patch(p)     = this%froot_curmr_patch(p)     + this%froot_xsmr_patch(p)
          this%livestem_mr_patch(p)  = this%livestem_curmr_patch(p)  + this%livestem_xsmr_patch(p)
          this%livecroot_mr_patch(p) = this%livecroot_curmr_patch(p) + this%livecroot_xsmr_patch(p)
       endif

       this%mr_patch(p)  = &
            this%leaf_mr_patch(p)     + &
            this%froot_mr_patch(p)    + &
            this%livestem_mr_patch(p) + &
            this%livecroot_mr_patch(p)

       if (carbon_resp_opt == 1) then
          this%mr_patch(p)  = &
               this%cpool_to_resp_patch(p)     + &
               this%leaf_mr_patch(p)     + &
               this%froot_mr_patch(p)    + &
               this%livestem_mr_patch(p) + &
               this%livecroot_mr_patch(p)
       end if
       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%mr_patch(p) = &
               this%mr_patch(p) + &
               this%grain_mr_patch(p)
       end if

       ! growth respiration (GR)

       ! current GR is respired this time step for new growth displayed in this timestep
       this%current_gr_patch(p) = &
            this%cpool_leaf_gr_patch(p)      + &
            this%cpool_froot_gr_patch(p)     + &
            this%cpool_livestem_gr_patch(p)  + &
            this%cpool_deadstem_gr_patch(p)  + &
            this%cpool_livecroot_gr_patch(p) + &
            this%cpool_deadcroot_gr_patch(p)
       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%current_gr_patch(p) = this%current_gr_patch(p) + &
               this%cpool_grain_gr_patch(p)
       end if


       ! transfer GR is respired this time step for transfer growth displayed in this timestep
       this%transfer_gr_patch(p) = &
            this%transfer_leaf_gr_patch(p)      + &
            this%transfer_froot_gr_patch(p)     + &
            this%transfer_livestem_gr_patch(p)  + &
            this%transfer_deadstem_gr_patch(p)  + &
            this%transfer_livecroot_gr_patch(p) + &
            this%transfer_deadcroot_gr_patch(p)
       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%transfer_gr_patch(p) = this%transfer_gr_patch(p) + &
               this%transfer_grain_gr_patch(p)
       end if

       ! storage GR is respired this time step for growth sent to storage for later display
       this%storage_gr_patch(p) = &
            this%cpool_leaf_storage_gr_patch(p)      + &
            this%cpool_froot_storage_gr_patch(p)     + &
            this%cpool_livestem_storage_gr_patch(p)  + &
            this%cpool_deadstem_storage_gr_patch(p)  + &
            this%cpool_livecroot_storage_gr_patch(p) + &
            this%cpool_deadcroot_storage_gr_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%storage_gr_patch(p) = this%storage_gr_patch(p) + &
               this%cpool_grain_storage_gr_patch(p)
       end if

       ! GR is the sum of current + transfer + storage GR
       this%gr_patch(p) = &
            this%current_gr_patch(p)  + &
            this%transfer_gr_patch(p) + &
            this%storage_gr_patch(p)


       ! autotrophic respiration (AR) adn 
       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%ar_patch(p) =           &
               this%mr_patch(p)      + &
               this%gr_patch(p)
          if ( .not. this%dribble_crophrv_xsmrpool_2atm ) this%ar_patch(p) = this%ar_patch(p) + &
                                             this%xsmrpool_to_atm_patch(p) ! xsmr... is -ve (slevis)
       else
             this%ar_patch(p) =           &
                  this%mr_patch(p)      + &
                  this%gr_patch(p)
       end if

       if (use_fun) then
          this%ar_patch(p) = this%ar_patch(p) + this%soilc_change_patch(p)
       end if

       ! gross primary production (GPP)
       this%gpp_patch(p) = &
            this%psnsun_to_cpool_patch(p) + &
            this%psnshade_to_cpool_patch(p)

       ! net primary production (NPP)      
       this%npp_patch(p) =      &
               this%gpp_patch(p) - &
               this%ar_patch(p)


       ! root respiration (RR)
       this%rr_patch(p) =         &
            this%froot_mr_patch(p)                   + &
            this%livecroot_mr_patch(p)               + &
            this%cpool_froot_gr_patch(p)             + &
            this%cpool_livecroot_gr_patch(p)         + &
            this%cpool_deadcroot_gr_patch(p)         + &
            this%transfer_froot_gr_patch(p)          + &
            this%transfer_livecroot_gr_patch(p)      + &
            this%transfer_deadcroot_gr_patch(p)      + &
            this%cpool_froot_storage_gr_patch(p)     + &
            this%cpool_livecroot_storage_gr_patch(p) + &
            this%cpool_deadcroot_storage_gr_patch(p)

       ! update the annual NPP accumulator, for use in allocation code 
       if (trim(isotope) == 'bulk') then
          this%tempsum_npp_patch(p) = &
               this%tempsum_npp_patch(p) + &
               this%npp_patch(p)
       end if

       ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
       ! This is supposed to correspond as closely as possible to
       ! field measurements of AGNPP, so it ignores the storage pools
       ! and only treats the fluxes into displayed pools.

       this%agnpp_patch(p) = &
            this%cpool_to_leafc_patch(p)                  + &
            this%leafc_xfer_to_leafc_patch(p)             + &
            this%cpool_to_livestemc_patch(p)              + &
            this%livestemc_xfer_to_livestemc_patch(p)     + &
            this%cpool_to_deadstemc_patch(p)              + &
            this%deadstemc_xfer_to_deadstemc_patch(p)


       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%agnpp_patch(p) =      &
               this%agnpp_patch(p) + &
               this%cpool_to_grainc_patch(p)            + &
               this%grainc_xfer_to_grainc_patch(p)
       end if

       ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
       ! This is supposed to correspond as closely as possible to
       ! field measurements of BGNPP, so it ignores the storage pools
       ! and only treats the fluxes into displayed pools.

       this%bgnpp_patch(p) = &
            this%cpool_to_frootc_patch(p)                   + &
            this%frootc_xfer_to_frootc_patch(p)             + &
            this%cpool_to_livecrootc_patch(p)               + &
            this%livecrootc_xfer_to_livecrootc_patch(p)     + &
            this%cpool_to_deadcrootc_patch(p)               + &
            this%deadcrootc_xfer_to_deadcrootc_patch(p)

       ! litterfall (LITFALL)

       this%litfall_patch(p) = &
            this%leafc_to_litter_patch(p)                     + &
            this%frootc_to_litter_patch(p)                    + &
            this%m_leafc_to_litter_patch(p)                   + &
            this%m_leafc_storage_to_litter_patch(p)           + &
            this%m_leafc_xfer_to_litter_patch(p)              + &
            this%m_frootc_to_litter_patch(p)                  + &
            this%m_frootc_storage_to_litter_patch(p)          + &
            this%m_frootc_xfer_to_litter_patch(p)             + &
            this%m_livestemc_to_litter_patch(p)               + &
            this%m_livestemc_storage_to_litter_patch(p)       + &
            this%m_livestemc_xfer_to_litter_patch(p)          + &
            this%m_deadstemc_to_litter_patch(p)               + &
            this%m_deadstemc_storage_to_litter_patch(p)       + &
            this%m_deadstemc_xfer_to_litter_patch(p)          + &
            this%m_livecrootc_to_litter_patch(p)              + &
            this%m_livecrootc_storage_to_litter_patch(p)      + &
            this%m_livecrootc_xfer_to_litter_patch(p)         + &
            this%m_deadcrootc_to_litter_patch(p)              + &
            this%m_deadcrootc_storage_to_litter_patch(p)      + &
            this%m_deadcrootc_xfer_to_litter_patch(p)         + &
            this%m_gresp_storage_to_litter_patch(p)           + &
            this%m_gresp_xfer_to_litter_patch(p)              + &

            this%m_leafc_to_litter_fire_patch(p)              + &
            this%m_leafc_storage_to_litter_fire_patch(p)      + &
            this%m_leafc_xfer_to_litter_fire_patch(p)         + &
            this%m_livestemc_to_litter_fire_patch(p)          + &
            this%m_livestemc_storage_to_litter_fire_patch(p)  + &
            this%m_livestemc_xfer_to_litter_fire_patch(p)     + &
            this%m_deadstemc_to_litter_fire_patch(p)          + &
            this%m_deadstemc_storage_to_litter_fire_patch(p)  + &
            this%m_deadstemc_xfer_to_litter_fire_patch(p)     + &
            this%m_frootc_to_litter_fire_patch(p)             + &
            this%m_frootc_storage_to_litter_fire_patch(p)     + &
            this%m_frootc_xfer_to_litter_fire_patch(p)        + &
            this%m_livecrootc_to_litter_fire_patch(p)         + &
            this%m_livecrootc_storage_to_litter_fire_patch(p) + &
            this%m_livecrootc_xfer_to_litter_fire_patch(p)    + &
            this%m_deadcrootc_to_litter_fire_patch(p)         + &
            this%m_deadcrootc_storage_to_litter_fire_patch(p) + &
            this%m_deadcrootc_xfer_to_litter_fire_patch(p)    + &
            this%m_gresp_storage_to_litter_fire_patch(p)      + &
            this%m_gresp_xfer_to_litter_fire_patch(p)         + &

            this%hrv_leafc_to_litter_patch(p)                 + &
            this%hrv_leafc_storage_to_litter_patch(p)         + &
            this%hrv_leafc_xfer_to_litter_patch(p)            + &
            this%hrv_frootc_to_litter_patch(p)                + &
            this%hrv_frootc_storage_to_litter_patch(p)        + &
            this%hrv_frootc_xfer_to_litter_patch(p)           + &
            this%hrv_livestemc_to_litter_patch(p)             + &
            this%hrv_livestemc_storage_to_litter_patch(p)     + &
            this%hrv_livestemc_xfer_to_litter_patch(p)        + &
            this%hrv_deadstemc_storage_to_litter_patch(p)     + &
            this%hrv_deadstemc_xfer_to_litter_patch(p)        + &
            this%hrv_livecrootc_to_litter_patch(p)            + &
            this%hrv_livecrootc_storage_to_litter_patch(p)    + &
            this%hrv_livecrootc_xfer_to_litter_patch(p)       + &
            this%hrv_deadcrootc_to_litter_patch(p)            + &
            this%hrv_deadcrootc_storage_to_litter_patch(p)    + &
            this%hrv_deadcrootc_xfer_to_litter_patch(p)       + &
            this%hrv_gresp_storage_to_litter_patch(p)         + &
            this%hrv_gresp_xfer_to_litter_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%litfall_patch(p) =      &
               this%litfall_patch(p) + &
               this%livestemc_to_litter_patch(p)

          if (.not. use_grainproduct) then
             this%litfall_patch(p) = &
                  this%litfall_patch(p) + &
                  this%grainc_to_food_patch(p)
          end if
       end if

       ! update the annual litfall accumulator, for use in mortality code

       if (use_cndv) then
          this%tempsum_litfall_patch(p) = &
               this%tempsum_litfall_patch(p) + &
               this%leafc_to_litter_patch(p) + &
               this%frootc_to_litter_patch(p)
       end if

       ! patch-level carbon losses to fire changed by F. Li and S. Levis

       this%fire_closs_patch(p) = &
            this%m_leafc_to_fire_patch(p)                + &
            this%m_leafc_storage_to_fire_patch(p)        + &
            this%m_leafc_xfer_to_fire_patch(p)           + &
            this%m_frootc_to_fire_patch(p)               + &
            this%m_frootc_storage_to_fire_patch(p)       + &
            this%m_frootc_xfer_to_fire_patch(p)          + &
            this%m_livestemc_to_fire_patch(p)            + &
            this%m_livestemc_storage_to_fire_patch(p)    + &
            this%m_livestemc_xfer_to_fire_patch(p)       + &
            this%m_deadstemc_to_fire_patch(p)            + &
            this%m_deadstemc_storage_to_fire_patch(p)    + &
            this%m_deadstemc_xfer_to_fire_patch(p)       + &
            this%m_livecrootc_to_fire_patch(p)           + &
            this%m_livecrootc_storage_to_fire_patch(p)   + &
            this%m_livecrootc_xfer_to_fire_patch(p)      + &
            this%m_deadcrootc_to_fire_patch(p)           + &
            this%m_deadcrootc_storage_to_fire_patch(p)   + &
            this%m_deadcrootc_xfer_to_fire_patch(p)      + &
            this%m_gresp_storage_to_fire_patch(p)        + &
            this%m_gresp_xfer_to_fire_patch(p)

       ! new summary variables for CLAMP

       ! (FROOTC_ALLOC) - fine root C allocation
       this%frootc_alloc_patch(p) = &
            this%frootc_xfer_to_frootc_patch(p)    + &
            this%cpool_to_frootc_patch(p)

       ! (FROOTC_LOSS) - fine root C loss changed by F. Li and S. Levis
       this%frootc_loss_patch(p) = &
            this%m_frootc_to_litter_patch(p)       + &
            this%m_frootc_to_fire_patch(p)         + &
            this%m_frootc_to_litter_fire_patch(p)  + &
            this%hrv_frootc_to_litter_patch(p)     + &
            this%frootc_to_litter_patch(p)

       ! (LEAFC_ALLOC) - leaf C allocation
       this%leafc_alloc_patch(p) = &
            this%leafc_xfer_to_leafc_patch(p)    + &
            this%cpool_to_leafc_patch(p)

       ! (LEAFC_LOSS) - leaf C loss changed by F. Li and S. Levis
       this%leafc_loss_patch(p) = &
            this%m_leafc_to_litter_patch(p)      + &
            this%m_leafc_to_fire_patch(p)        + &
            this%m_leafc_to_litter_fire_patch(p) + &
            this%hrv_leafc_to_litter_patch(p)    + &
            this%leafc_to_litter_patch(p)

       ! (WOODC_ALLOC) - wood C allocation
       this%woodc_alloc_patch(p) = &
            this%livestemc_xfer_to_livestemc_patch(p)   + &
            this%deadstemc_xfer_to_deadstemc_patch(p)   + &
            this%livecrootc_xfer_to_livecrootc_patch(p) + &
            this%deadcrootc_xfer_to_deadcrootc_patch(p) + &
            this%cpool_to_livestemc_patch(p)            + &
            this%cpool_to_deadstemc_patch(p)            + &
            this%cpool_to_livecrootc_patch(p)           + &
            this%cpool_to_deadcrootc_patch(p)



       ! (WOODC_LOSS) - wood C loss
       this%woodc_loss_patch(p) = &
            this%m_livestemc_to_litter_patch(p)            + &
            this%m_deadstemc_to_litter_patch(p)            + &
            this%m_livecrootc_to_litter_patch(p)           + &
            this%m_deadcrootc_to_litter_patch(p)           + &
            this%m_livestemc_to_fire_patch(p)              + &
            this%m_deadstemc_to_fire_patch(p)              + &
            this%m_livecrootc_to_fire_patch(p)             + &
            this%m_deadcrootc_to_fire_patch(p)             + &
            this%hrv_livestemc_to_litter_patch(p)          + &
            this%hrv_livestemc_storage_to_litter_patch(p)  + &
            this%hrv_livestemc_xfer_to_litter_patch(p)     + &
            this%wood_harvestc_patch(p)                    + &
            this%hrv_deadstemc_storage_to_litter_patch(p)  + &
            this%hrv_deadstemc_xfer_to_litter_patch(p)     + &
            this%hrv_livecrootc_to_litter_patch(p)         + &
            this%hrv_livecrootc_storage_to_litter_patch(p) + &
            this%hrv_livecrootc_xfer_to_litter_patch(p)    + &
            this%hrv_deadcrootc_to_litter_patch(p)         + &
            this%hrv_deadcrootc_storage_to_litter_patch(p) + &
            this%hrv_deadcrootc_xfer_to_litter_patch(p)


       ! (Slash Harvest Flux) - Additional Wood Harvest Veg C Losses
       this%slash_harvestc_patch(p) = &
            this%hrv_leafc_to_litter_patch(p)              + &
            this%hrv_leafc_storage_to_litter_patch(p)      + &
            this%hrv_leafc_xfer_to_litter_patch(p)         + &
            this%hrv_frootc_to_litter_patch(p)             + &
            this%hrv_frootc_storage_to_litter_patch(p)     + &
            this%hrv_frootc_xfer_to_litter_patch(p)        + &
            this%hrv_livestemc_to_litter_patch(p)          + &
            this%hrv_livestemc_storage_to_litter_patch(p)  + &
            this%hrv_livestemc_xfer_to_litter_patch(p)     + &
            this%hrv_deadstemc_storage_to_litter_patch(p)  + &
            this%hrv_deadstemc_xfer_to_litter_patch(p)     + &
            this%hrv_livecrootc_to_litter_patch(p)         + &
            this%hrv_livecrootc_storage_to_litter_patch(p) + &
            this%hrv_livecrootc_xfer_to_litter_patch(p)    + &
            this%hrv_deadcrootc_to_litter_patch(p)         + &
            this%hrv_deadcrootc_storage_to_litter_patch(p) + &
            this%hrv_deadcrootc_xfer_to_litter_patch(p)    + &
            this%hrv_xsmrpool_to_atm_patch(p)              + &
            this%hrv_gresp_storage_to_litter_patch(p)      + &
            this%hrv_gresp_xfer_to_litter_patch(p)

    end do  ! end of patches loop


    !------------------------------------------------
    ! column variables
    !------------------------------------------------

    ! use p2c routine to get selected column-average patch-level fluxes and states

    call p2c(bounds, num_soilc, filter_soilc, &
         this%hrv_xsmrpool_to_atm_patch(bounds%begp:bounds%endp), &
         this%hrv_xsmrpool_to_atm_col(bounds%begc:bounds%endc))

    if (use_crop .and. this%dribble_crophrv_xsmrpool_2atm) then
       call p2c(bounds, num_soilc, filter_soilc, &
            this%xsmrpool_to_atm_patch(bounds%begp:bounds%endp), &
            this%xsmrpool_to_atm_col(bounds%begc:bounds%endc))

       call c2g( bounds = bounds, &
            carr = this%xsmrpool_to_atm_col(bounds%begc:bounds%endc), &
            garr = this%xsmrpool_to_atm_grc(bounds%begg:bounds%endg), &
            c2l_scale_type = 'unity', &
            l2g_scale_type = 'unity')
    end if

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_closs_patch(bounds%begp:bounds%endp), &
         this%fire_closs_p2c_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%npp_patch(bounds%begp:bounds%endp), &
         this%npp_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%rr_patch(bounds%begp:bounds%endp), &
         this%rr_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%ar_patch(bounds%begp:bounds%endp), &
         this%ar_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%gpp_patch(bounds%begp:bounds%endp), &
         this%gpp_col(bounds%begc:bounds%endc))


    ! this code is to calculate an exponentially-relaxed npp value for use in NDynamics code

    if ( trim(isotope) == 'bulk') then
       if (nfix_timeconst > 0._r8 .and. nfix_timeconst < 500._r8 ) then
          nfixlags = nfix_timeconst * secspday
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             if ( this%lag_npp_col(c) /= spval ) then
                this%lag_npp_col(c) = &
                     this%lag_npp_col(c) * exp(-dtime/nfixlags) + &
                     this%npp_col(c) * (1._r8 - exp(-dtime/nfixlags))
             else
                ! first timestep
                this%lag_npp_col(c) = this%npp_col(c)
             endif
          end do
       endif
    endif


    ! vertically integrate column-level carbon fire losses
!    if(.not. use_soil_matrixcn)then
       do l = 1, ndecomp_pools
          do j = 1,nlevdecomp
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                this%m_decomp_cpools_to_fire_col(c,l) = &
                     this%m_decomp_cpools_to_fire_col(c,l) + &
                     this%m_decomp_cpools_to_fire_vr_col(c,j,l)*dzsoi_decomp(j)
             end do
          end do
       end do
!    end if !not use_soil_matrixcn


    do fc = 1,num_soilc
       c = filter_soilc(fc)

       g = col%gridcell(c)

       ! litter fire losses (LITFIRE)
       this%litfire_col(c) = 0._r8

       ! soil organic matter fire losses (SOMFIRE)
       this%somfire_col(c) = 0._r8

       ! total ecosystem fire losses (TOTFIRE)
       this%totfire_col(c) = &
            this%litfire_col(c) + &
            this%somfire_col(c)

       ! carbon losses to fire, including patch losses
       this%fire_closs_col(c) = this%fire_closs_p2c_col(c)
       do l = 1, ndecomp_pools
          this%fire_closs_col(c) = &
               this%fire_closs_col(c) + &
               this%m_decomp_cpools_to_fire_col(c,l)
       end do

       ! total soil respiration, heterotrophic + root respiration (SR)
       this%sr_col(c) = &
            this%rr_col(c) + &
            soilbiogeochem_hr_col(c)

       ! total ecosystem respiration, autotrophic + heterotrophic (ER)
       this%er_col(c) = &
            this%ar_col(c) + &
            soilbiogeochem_hr_col(c)

       ! coarse woody debris heterotrophic respiration
       this%cwdc_hr_col(c) = 0._r8
       ! net ecosystem production, excludes fire flux, landcover change, 
       ! and loss from wood products, positive for sink (NEP)
       this%nep_col(c) = &
            this%gpp_col(c) - &
            this%er_col(c)

    end do

    call c2g( bounds = bounds, &
         carr = this%nep_col(bounds%begc:bounds%endc), &
         garr = nep_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    call c2g( bounds = bounds, &
         carr = this%fire_closs_col(bounds%begc:bounds%endc), &
         garr = fire_closs_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')

    call c2g( bounds = bounds, &
         carr = this%hrv_xsmrpool_to_atm_col(bounds%begc:bounds%endc), &
         garr = hrv_xsmrpool_to_atm_grc(bounds%begg:bounds%endg), &
         c2l_scale_type = 'unity', &
         l2g_scale_type = 'unity')
    hrv_xsmrpool_to_atm_delta_grc(bounds%begg:bounds%endg) = &
         hrv_xsmrpool_to_atm_grc(bounds%begg:bounds%endg) * dtime
    call this%hrv_xsmrpool_to_atm_dribbler%set_curr_delta(bounds, &
         hrv_xsmrpool_to_atm_delta_grc(bounds%begg:bounds%endg))
    call this%hrv_xsmrpool_to_atm_dribbler%get_curr_flux(bounds, &
         hrv_xsmrpool_to_atm_dribbled_grc(bounds%begg:bounds%endg))


    dwt_conv_cflux_delta_grc(bounds%begg:bounds%endg) = &
         this%dwt_conv_cflux_grc(bounds%begg:bounds%endg) * dtime
    call this%dwt_conv_cflux_dribbler%set_curr_delta(bounds, &
         dwt_conv_cflux_delta_grc(bounds%begg:bounds%endg))
    call this%dwt_conv_cflux_dribbler%get_curr_flux(bounds, &
         this%dwt_conv_cflux_dribbled_grc(bounds%begg:bounds%endg))

    do g = bounds%begg, bounds%endg
       ! net ecosystem exchange of carbon, includes fire flux and hrv_xsmrpool flux,
       ! positive for source (NEE)
       this%nee_grc(g) = &
            -nep_grc(g)       + &
            fire_closs_grc(g) + &
            hrv_xsmrpool_to_atm_dribbled_grc(g)

       this%landuseflux_grc(g) = &
            this%dwt_conv_cflux_dribbled_grc(g)   + &
            product_closs_grc(g)

       ! net biome production of carbon, positive for sink
       this%nbp_grc(g) = &
            -this%nee_grc(g)        - &
            this%landuseflux_grc(g)
       if ( this%dribble_crophrv_xsmrpool_2atm ) this%nbp_grc(g) = this%nbp_grc(g) - this%xsmrpool_to_atm_grc(g)
    end do

    ! coarse woody debris C loss
    do fc = 1,num_soilc
       c  = filter_soilc(fc)
       this%cwdc_loss_col(c)  = 0._r8
    end do
    associate(is_cwd    => decomp_cascade_con%is_cwd) ! TRUE => pool is a cwd pool   
      do l = 1, ndecomp_pools
         if ( is_cwd(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%cwdc_loss_col(c) = &
                    this%cwdc_loss_col(c) + &
                    this%m_decomp_cpools_to_fire_col(c,l)
            end do
         end if
      end do
      do k = 1, ndecomp_cascade_transitions
         if ( is_cwd(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%cwdc_loss_col(c) = &
                    this%cwdc_loss_col(c) + &
                    soilbiogeochem_decomp_cascade_ctransfer_col(c,k)
            end do
         end if
      end do
    end associate



    ! litter C loss      
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%litterc_loss_col(c) = soilbiogeochem_lithr_col(c)
    end do
    associate(is_litter => decomp_cascade_con%is_litter) ! TRUE => pool is a litter pool
      do l = 1, ndecomp_pools
         if ( is_litter(l) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%litterc_loss_col(c) = &
                    this%litterc_loss_col(c) + &
                    this%m_decomp_cpools_to_fire_col(c,l)
            end do
         end if
      end do
      do k = 1, ndecomp_cascade_transitions
         if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               this%litterc_loss_col(c) = &
                    this%litterc_loss_col(c) + &
                    soilbiogeochem_decomp_cascade_ctransfer_col(c,k)
            end do
         end if
      end do
    end associate

  end subroutine Summary_carbonflux
!-----------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_carbonflux_type) :: this
    type(bounds_type), intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c, g, j          ! indices
    !-----------------------------------------------------------------------

    ! set conversion and product pool fluxes to 0 at the beginning of every timestep

    do g = bounds%begg, bounds%endg
       this%dwt_seedc_to_leaf_grc(g)        = 0._r8
       this%dwt_seedc_to_deadstem_grc(g)    = 0._r8
       this%dwt_conv_cflux_grc(g)           = 0._r8
       this%dwt_slash_cflux_grc(g)           = 0._r8
    end do

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          this%dwt_frootc_to_litr_met_c_col(c,j)    = 0._r8
          this%dwt_frootc_to_litr_cel_c_col(c,j)    = 0._r8
          this%dwt_frootc_to_litr_lig_c_col(c,j)    = 0._r8
          this%dwt_livecrootc_to_cwdc_col(c,j)      = 0._r8
          this%dwt_deadcrootc_to_cwdc_col(c,j)      = 0._r8
       end do
    end do

  end subroutine ZeroDwt

  !-----------------------------------------------------------------------
end module CNVegCarbonFluxType


