module CNVegNitrogenFluxType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar       , only : nlevdecomp_full, nlevgrnd,nlevdecomp
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp,nvegnpool,&
                                nnphtrans,nngmtrans,nnfitrans,nnphouttrans,&
                                nngmouttrans,nnfiouttrans
  use clm_varpar        , only : ileaf,ileaf_st,ileaf_xf,ifroot,ifroot_st,ifroot_xf,&
                                 ilivestem,ilivestem_st,ilivestem_xf,&
                                 ideadstem,ideadstem_st,ideadstem_xf,&
                                 ilivecroot,ilivecroot_st,ilivecroot_xf,&
                                 ideadcroot,ideadcroot_st,ideadcroot_xf,&
                                 igrain,igrain_st,igrain_xf,iretransn,ioutn
  use clm_varpar       , only : numpft, num_zon, num_veg, &
                                var_col, var_pft, CN_zone_weight
  use clm_varcon       , only : spval, ispval, dzsoi_decomp
  use clm_varctl       , only : use_nitrif_denitrif, use_vertsoilc, use_crop, use_matrixcn
  use PatchType        , only : patch
  use CNSharedParamsMod , only : use_fun
  use LandunitType      , only : lun
  use landunit_varcon  , only : istsoil, istcrop


  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

  type, public :: cnveg_nitrogenflux_type

     ! gap mortality fluxes
     real(r8), pointer :: m_leafn_to_litter_patch                   (:)     ! patch leaf N mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_to_litter_patch                  (:)     ! patch fine root N mortality (gN/m2/s)
     real(r8), pointer :: m_leafn_storage_to_litter_patch           (:)     ! patch leaf N storage mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_storage_to_litter_patch          (:)     ! patch fine root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_storage_to_litter_patch       (:)     ! patch live stem N storage mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_storage_to_litter_patch       (:)     ! patch dead stem N storage mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_storage_to_litter_patch      (:)     ! patch live coarse root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_storage_to_litter_patch      (:)     ! patch dead coarse root N storage mortality (gN/m2/s)
     real(r8), pointer :: m_leafn_xfer_to_litter_patch              (:)     ! patch leaf N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_frootn_xfer_to_litter_patch             (:)     ! patch fine root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_xfer_to_litter_patch          (:)     ! patch live stem N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_xfer_to_litter_patch          (:)     ! patch dead stem N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_xfer_to_litter_patch         (:)     ! patch live coarse root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_xfer_to_litter_patch         (:)     ! patch dead coarse root N transfer mortality (gN/m2/s)
     real(r8), pointer :: m_livestemn_to_litter_patch               (:)     ! patch live stem N mortality (gN/m2/s)
     real(r8), pointer :: m_deadstemn_to_litter_patch               (:)     ! patch dead stem N mortality (gN/m2/s)
     real(r8), pointer :: m_livecrootn_to_litter_patch              (:)     ! patch live coarse root N mortality (gN/m2/s)
     real(r8), pointer :: m_deadcrootn_to_litter_patch              (:)     ! patch dead coarse root N mortality (gN/m2/s)
     real(r8), pointer :: m_retransn_to_litter_patch                (:)     ! patch retranslocated N pool mortality (gN/m2/s)

     ! harvest fluxes
     real(r8), pointer :: hrv_leafn_to_litter_patch                 (:)     ! patch leaf N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_to_litter_patch                (:)     ! patch fine root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_leafn_storage_to_litter_patch         (:)     ! patch leaf N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_storage_to_litter_patch        (:)     ! patch fine root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_storage_to_litter_patch     (:)     ! patch live stem N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_storage_to_litter_patch     (:)     ! patch dead stem N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_storage_to_litter_patch    (:)     ! patch live coarse root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_storage_to_litter_patch    (:)     ! patch dead coarse root N storage harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_leafn_xfer_to_litter_patch            (:)     ! patch leaf N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_frootn_xfer_to_litter_patch           (:)     ! patch fine root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_xfer_to_litter_patch        (:)     ! patch live stem N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadstemn_xfer_to_litter_patch        (:)     ! patch dead stem N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_xfer_to_litter_patch       (:)     ! patch live coarse root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_xfer_to_litter_patch       (:)     ! patch dead coarse root N transfer harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livestemn_to_litter_patch             (:)     ! patch live stem N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_livecrootn_to_litter_patch            (:)     ! patch live coarse root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_deadcrootn_to_litter_patch            (:)     ! patch dead coarse root N harvest mortality (gN/m2/s)
     real(r8), pointer :: hrv_retransn_to_litter_patch              (:)     ! patch retranslocated N pool harvest mortality (gN/m2/s)
     real(r8), pointer :: grainn_to_cropprodn_patch                 (:)     ! patch grain N to crop product pool (gN/m2/s)
     real(r8), pointer :: grainn_to_cropprodn_col                   (:)     ! col grain N to crop product pool (gN/m2/s)
     real(r8), pointer :: m_n_to_litr_met_fire_col                  (:,:)   ! col N from leaf, froot, xfer and storage N to litter labile N by fire (gN/m3/s)
     real(r8), pointer :: m_n_to_litr_cel_fire_col                  (:,:)   ! col N from leaf, froot, xfer and storage N to litter cellulose N by fire (gN/m3/s) 
     real(r8), pointer :: m_n_to_litr_lig_fire_col                  (:,:)   ! col N from leaf, froot, xfer and storage N to litter lignin N by fire (gN/m3/s) 
     real(r8), pointer :: harvest_n_to_litr_met_n_col               (:,:)   ! col N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
     real(r8), pointer :: harvest_n_to_litr_cel_n_col               (:,:)   ! col N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
     real(r8), pointer :: harvest_n_to_litr_lig_n_col               (:,:)   ! col N fluxes associated with harvest to litter lignin pool (gN/m3/s)
     real(r8), pointer :: harvest_n_to_cwdn_col                     (:,:)   ! col N fluxes associated with harvest to CWD pool (gN/m3/s)

     ! fire N fluxes 
     real(r8), pointer :: m_decomp_npools_to_fire_vr_col            (:,:,:) ! col vertically-resolved decomposing N fire loss (gN/m3/s)
     real(r8), pointer :: m_decomp_npools_to_fire_col               (:,:)   ! col vertically-integrated (diagnostic) decomposing N fire loss (gN/m2/s)
     real(r8), pointer :: m_leafn_to_fire_patch                     (:)     ! patch (gN/m2/s) fire N emissions from leafn 
     real(r8), pointer :: m_leafn_storage_to_fire_patch             (:)     ! patch (gN/m2/s) fire N emissions from leafn_storage            
     real(r8), pointer :: m_leafn_xfer_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from leafn_xfer     
     real(r8), pointer :: m_livestemn_to_fire_patch                 (:)     ! patch (gN/m2/s) fire N emissions from livestemn 
     real(r8), pointer :: m_livestemn_storage_to_fire_patch         (:)     ! patch (gN/m2/s) fire N emissions from livestemn_storage      
     real(r8), pointer :: m_livestemn_xfer_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from livestemn_xfer
     real(r8), pointer :: m_deadstemn_to_fire_patch                 (:)     ! patch (gN/m2/s) fire N emissions from deadstemn
     real(r8), pointer :: m_deadstemn_storage_to_fire_patch         (:)     ! patch (gN/m2/s) fire N emissions from deadstemn_storage         
     real(r8), pointer :: m_deadstemn_xfer_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from deadstemn_xfer
     real(r8), pointer :: m_frootn_to_fire_patch                    (:)     ! patch (gN/m2/s) fire N emissions from frootn
     real(r8), pointer :: m_frootn_storage_to_fire_patch            (:)     ! patch (gN/m2/s) fire N emissions from frootn_storage
     real(r8), pointer :: m_frootn_xfer_to_fire_patch               (:)     ! patch (gN/m2/s) fire N emissions from frootn_xfer
     real(r8), pointer :: m_livecrootn_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from m_livecrootn_to_fire
     real(r8), pointer :: m_livecrootn_storage_to_fire_patch        (:)     ! patch (gN/m2/s) fire N emissions from livecrootn_storage     
     real(r8), pointer :: m_livecrootn_xfer_to_fire_patch           (:)     ! patch (gN/m2/s) fire N emissions from livecrootn_xfer
     real(r8), pointer :: m_deadcrootn_to_fire_patch                (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn
     real(r8), pointer :: m_deadcrootn_storage_to_fire_patch        (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn_storage  
     real(r8), pointer :: m_deadcrootn_xfer_to_fire_patch           (:)     ! patch (gN/m2/s) fire N emissions from deadcrootn_xfer
     real(r8), pointer :: m_retransn_to_fire_patch                  (:)     ! patch (gN/m2/s) fire N emissions from retransn
     real(r8), pointer :: m_leafn_to_litter_fire_patch              (:)     ! patch (gN/m2/s) from leafn to litter N  due to fire               
     real(r8), pointer :: m_leafn_storage_to_litter_fire_patch      (:)     ! patch (gN/m2/s) from leafn_storage to litter N  due to fire                              
     real(r8), pointer :: m_leafn_xfer_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from leafn_xfer to litter N  due to fire                              
     real(r8), pointer :: m_livestemn_to_litter_fire_patch          (:)     ! patch (gN/m2/s) from livestemn to litter N  due to fire                              
     real(r8), pointer :: m_livestemn_storage_to_litter_fire_patch  (:)     ! patch (gN/m2/s) from livestemn_storage to litter N  due to fire                                     
     real(r8), pointer :: m_livestemn_xfer_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from livestemn_xfer to litter N  due to fire                                     
     real(r8), pointer :: m_livestemn_to_deadstemn_fire_patch       (:)     ! patch (gN/m2/s) from livestemn to deadstemn N  due to fire                                     
     real(r8), pointer :: m_deadstemn_to_litter_fire_patch          (:)     ! patch (gN/m2/s) from deadstemn to litter N  due to fire                                     
     real(r8), pointer :: m_deadstemn_storage_to_litter_fire_patch  (:)     ! patch (gN/m2/s) from deadstemn_storage to litter N  due to fire                                               
     real(r8), pointer :: m_deadstemn_xfer_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from deadstemn_xfer to litter N  due to fire                                               
     real(r8), pointer :: m_frootn_to_litter_fire_patch             (:)     ! patch (gN/m2/s) from frootn to litter N  due to fire                                               
     real(r8), pointer :: m_frootn_storage_to_litter_fire_patch     (:)     ! patch (gN/m2/s) from frootn_storage to litter N  due to fire                                               
     real(r8), pointer :: m_frootn_xfer_to_litter_fire_patch        (:)     ! patch (gN/m2/s) from frootn_xfer to litter N  due to fire        
     real(r8), pointer :: m_livecrootn_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from livecrootn to litter N  due to fire                                               
     real(r8), pointer :: m_livecrootn_storage_to_litter_fire_patch (:)     ! patch (gN/m2/s) from livecrootn_storage to litter N  due to fire                                                     
     real(r8), pointer :: m_livecrootn_xfer_to_litter_fire_patch    (:)     ! patch (gN/m2/s) from livecrootn_xfer to litter N  due to fire                                                     
     real(r8), pointer :: m_livecrootn_to_deadcrootn_fire_patch     (:)     ! patch (gN/m2/s) from livecrootn_xfer to deadcrootn due to fire                                                     
     real(r8), pointer :: m_deadcrootn_to_litter_fire_patch         (:)     ! patch (gN/m2/s) from deadcrootn to deadcrootn due to fire                                                       
     real(r8), pointer :: m_deadcrootn_storage_to_litter_fire_patch (:)     ! patch (gN/m2/s) from deadcrootn_storage to deadcrootn due to fire                                                        
     real(r8), pointer :: m_deadcrootn_xfer_to_litter_fire_patch    (:)     ! patch (gN/m2/s) from deadcrootn_xfer to deadcrootn due to fire                                                         
     real(r8), pointer :: m_retransn_to_litter_fire_patch           (:)     ! patch (gN/m2/s) from retransn to deadcrootn due to fire                                                         
     real(r8), pointer :: fire_nloss_patch                          (:)     ! patch total patch-level fire N loss (gN/m2/s) 
     real(r8), pointer :: fire_nloss_col                            (:)     ! col total column-level fire N loss (gN/m2/s)
     real(r8), pointer :: fire_nloss_p2c_col                        (:)     ! col patch2col column-level fire N loss (gN/m2/s) (p2c)
     real(r8), pointer :: fire_mortality_n_to_cwdn_col              (:,:)   ! col N fluxes associated with fire mortality to CWD pool (gN/m3/s)

     ! phenology fluxes from transfer pool
     real(r8), pointer :: grainn_xfer_to_grainn_patch               (:)     ! patch grain N growth from storage for prognostic crop model (gN/m2/s)
     real(r8), pointer :: leafn_xfer_to_leafn_patch                 (:)     ! patch leaf N growth from storage (gN/m2/s)
     real(r8), pointer :: frootn_xfer_to_frootn_patch               (:)     ! patch fine root N growth from storage (gN/m2/s)
     real(r8), pointer :: livestemn_xfer_to_livestemn_patch         (:)     ! patch live stem N growth from storage (gN/m2/s)
     real(r8), pointer :: deadstemn_xfer_to_deadstemn_patch         (:)     ! patch dead stem N growth from storage (gN/m2/s)
     real(r8), pointer :: livecrootn_xfer_to_livecrootn_patch       (:)     ! patch live coarse root N growth from storage (gN/m2/s)
     real(r8), pointer :: deadcrootn_xfer_to_deadcrootn_patch       (:)     ! patch dead coarse root N growth from storage (gN/m2/s)

     ! litterfall fluxes
     real(r8), pointer :: livestemn_to_litter_patch                 (:)     ! patch livestem N to litter (gN/m2/s)
     real(r8), pointer :: grainn_to_food_patch                      (:)     ! patch grain N to food for prognostic crop (gN/m2/s)
     real(r8), pointer :: leafn_to_biofueln_patch                   (:)     ! patch leaf N to biofuel N (gN/m2/s)
     real(r8), pointer :: livestemn_to_biofueln_patch               (:)     ! patch livestem N to biofuel N (gN/m2/s)
     real(r8), pointer :: grainn_to_seed_patch                      (:)     ! patch grain N to seed for prognostic crop (gN/m2/s)
     real(r8), pointer :: leafn_to_litter_patch                     (:)     ! patch leaf N litterfall (gN/m2/s)
     real(r8), pointer :: leafn_to_retransn_patch                   (:)     ! patch leaf N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: frootn_to_retransn_patch                  (:)     ! patch fine root N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: frootn_to_litter_patch                    (:)     ! patch fine root N litterfall (gN/m2/s)

     ! allocation fluxes
     real(r8), pointer :: retransn_to_npool_patch                   (:)     ! patch deployment of retranslocated N (gN/m2/s)  
     real(r8), pointer :: free_retransn_to_npool_patch              (:)     ! patch deployment of free retranslocated N (gN/m2/s)           
     real(r8), pointer :: sminn_to_npool_patch                      (:)     ! patch deployment of soil mineral N uptake (gN/m2/s)
     real(r8), pointer :: npool_to_grainn_patch                     (:)     ! patch allocation to grain N for prognostic crop (gN/m2/s)
     real(r8), pointer :: npool_to_grainn_storage_patch             (:)     ! patch allocation to grain N storage for prognostic crop (gN/m2/s)
     real(r8), pointer :: npool_to_leafn_patch                      (:)     ! patch allocation to leaf N (gN/m2/s)
     real(r8), pointer :: npool_to_leafn_storage_patch              (:)     ! patch allocation to leaf N storage (gN/m2/s)
     real(r8), pointer :: npool_to_frootn_patch                     (:)     ! patch allocation to fine root N (gN/m2/s)
     real(r8), pointer :: npool_to_frootn_storage_patch             (:)     ! patch allocation to fine root N storage (gN/m2/s)
     real(r8), pointer :: npool_to_livestemn_patch                  (:)     ! patch allocation to live stem N (gN/m2/s)
     real(r8), pointer :: npool_to_livestemn_storage_patch          (:)     ! patch allocation to live stem N storage (gN/m2/s)
     real(r8), pointer :: npool_to_deadstemn_patch                  (:)     ! patch allocation to dead stem N (gN/m2/s)
     real(r8), pointer :: npool_to_deadstemn_storage_patch          (:)     ! patch allocation to dead stem N storage (gN/m2/s)
     real(r8), pointer :: npool_to_livecrootn_patch                 (:)     ! patch allocation to live coarse root N (gN/m2/s)
     real(r8), pointer :: npool_to_livecrootn_storage_patch         (:)     ! patch allocation to live coarse root N storage (gN/m2/s)
     real(r8), pointer :: npool_to_deadcrootn_patch                 (:)     ! patch allocation to dead coarse root N (gN/m2/s)
     real(r8), pointer :: npool_to_deadcrootn_storage_patch         (:)     ! patch allocation to dead coarse root N storage (gN/m2/s)

     ! annual turnover of storage to transfer pools           
     real(r8), pointer :: grainn_storage_to_xfer_patch              (:)     ! patch grain N shift storage to transfer for prognostic crop (gN/m2/s)
     real(r8), pointer :: leafn_storage_to_xfer_patch               (:)     ! patch leaf N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: frootn_storage_to_xfer_patch              (:)     ! patch fine root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: livestemn_storage_to_xfer_patch           (:)     ! patch live stem N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: deadstemn_storage_to_xfer_patch           (:)     ! patch dead stem N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: livecrootn_storage_to_xfer_patch          (:)     ! patch live coarse root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: deadcrootn_storage_to_xfer_patch          (:)     ! patch dead coarse root N shift storage to transfer (gN/m2/s)
     real(r8), pointer :: fert_patch                                (:)     ! patch applied fertilizer (gN/m2/s)
     real(r8), pointer :: fert_counter_patch                        (:)     ! patch >0 fertilize; <=0 not
     real(r8), pointer :: soyfixn_patch                             (:)     ! patch soybean fixed N (gN/m2/s)

     ! turnover of livewood to deadwood, with retranslocation 
     real(r8), pointer :: livestemn_to_deadstemn_patch              (:)     ! patch live stem N turnover (gN/m2/s)
     real(r8), pointer :: livestemn_to_retransn_patch               (:)     ! patch live stem N to retranslocated N pool (gN/m2/s)
     real(r8), pointer :: livecrootn_to_deadcrootn_patch            (:)     ! patch live coarse root N turnover (gN/m2/s)
     real(r8), pointer :: livecrootn_to_retransn_patch              (:)     ! patch live coarse root N to retranslocated N pool (gN/m2/s)

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: ndeploy_patch                             (:)     ! patch total N deployed to growth and storage (gN/m2/s)
     real(r8), pointer :: wood_harvestn_patch                       (:)     ! patch total N losses to wood product pools (gN/m2/s)
     real(r8), pointer :: wood_harvestn_col                         (:)     ! col total N losses to wood product pools (gN/m2/s) (p2c)
     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_n_to_litr_met_n_col             (:,:)   ! col N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
     real(r8), pointer :: phenology_n_to_litr_cel_n_col             (:,:)   ! col N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
     real(r8), pointer :: phenology_n_to_litr_lig_n_col             (:,:)   ! col N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)

     ! gap mortality fluxes
     real(r8), pointer :: gap_mortality_n_to_litr_met_n_col         (:,:)   ! col N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_litr_cel_n_col         (:,:)   ! col N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_litr_lig_n_col         (:,:)   ! col N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
     real(r8), pointer :: gap_mortality_n_to_cwdn_col               (:,:)   ! col N fluxes associated with gap mortality to CWD pool (gN/m3/s)

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedn_to_leaf_patch                   (:)     ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedn_to_leaf_grc                     (:)     ! (gN/m2/s) dwt_seedn_to_leaf_patch summed to the gridcell-level
     real(r8), pointer :: dwt_seedn_to_deadstem_patch               (:)     ! (gN/m2/s) seed source to patch-level; although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_seedn_to_deadstem_grc                 (:)     ! (gN/m2/s) dwt_seedn_to_deadstem_patch summed to the gridcell-level
     real(r8), pointer :: dwt_conv_nflux_patch                      (:)     ! (gN/m2/s) conversion N flux (immediate loss to atm); although this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_conv_nflux_grc                        (:)     ! (gN/m2/s) dwt_conv_nflux_patch summed to the gridcell-level
     real(r8), pointer :: dwt_wood_productn_gain_patch              (:)     ! patch (gN/m2/s) addition to wood product pools from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_crop_productn_gain_patch              (:)     ! patch (gN/m2/s) addition to crop product pool from landcover change; even though this is a patch-level flux, it is expressed per unit GRIDCELL area
     real(r8), pointer :: dwt_frootn_to_litr_met_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootn_to_litr_cel_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootn_to_litr_lig_n_col              (:,:)   ! col (gN/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootn_to_cwdn_col                (:,:)   ! col (gN/m3/s) dead coarse root to CWD due to landcover change

     ! crop fluxes
     real(r8), pointer :: crop_seedn_to_leaf_patch                  (:)     ! patch (gN/m2/s) seed source to leaf, for crops

     ! Misc
     real(r8), pointer :: plant_ndemand_patch                       (:)     ! N flux required to support initial GPP (gN/m2/s)
     real(r8), pointer :: avail_retransn_patch                      (:)     ! N flux available from retranslocation pool (gN/m2/s)
     real(r8), pointer :: plant_nalloc_patch                        (:)     ! total allocated N flux (gN/m2/s)
     real(r8), pointer :: plant_ndemand_retrans_patch               (:)     ! The N demand pool generated for FUN2.0; mainly used for deciduous trees (gN/m2/s)
     real(r8), pointer :: plant_ndemand_season_patch                (:)     ! The N demand pool for seasonal deciduous (gN/m2/s)
     real(r8), pointer :: plant_ndemand_stress_patch                (:)     ! The N demand pool for stress deciduous   (gN/m2/s)
     real(r8), pointer :: Nactive_patch                             (:)     ! N acquired by mycorrhizal uptake  (gN/m2/s)
     real(r8), pointer :: Nnonmyc_patch                             (:)     ! N acquired by non-myc uptake      (gN/m2/s)
     real(r8), pointer :: Nam_patch                                 (:)     ! N acquired by AM plant            (gN/m2/s)
     real(r8), pointer :: Necm_patch                                (:)     ! N acquired by ECM plant           (gN/m2/s)
     real(r8), pointer :: Nactive_no3_patch                         (:)     ! N acquired by mycorrhizal uptake  (gN/m2/s)
     real(r8), pointer :: Nactive_nh4_patch                         (:)     ! N acquired by mycorrhizal uptake  (gN/m2/s)
     real(r8), pointer :: Nnonmyc_no3_patch                         (:)     ! N acquired by non-myc             (gN/m2/s)
     real(r8), pointer :: Nnonmyc_nh4_patch                         (:)     ! N acquired by non-myc             (gN/m2/s)
     real(r8), pointer :: Nam_no3_patch                             (:)     ! N acquired by AM plant            (gN/m2/s)
     real(r8), pointer :: Nam_nh4_patch                             (:)     ! N acquired by AM plant            (gN/m2/s)
     real(r8), pointer :: Necm_no3_patch                            (:)     ! N acquired by ECM plant           (gN/m2/s)
     real(r8), pointer :: Necm_nh4_patch                            (:)     ! N acquired by ECM plant           (gN/m2/s)
     real(r8), pointer :: Nfix_patch                                (:)     ! N acquired by Symbiotic BNF       (gN/m2/s)
     real(r8), pointer :: Npassive_patch                            (:)     ! N acquired by passive uptake      (gN/m2/s)
     real(r8), pointer :: Nretrans_patch                            (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nretrans_org_patch                        (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nretrans_season_patch                     (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nretrans_stress_patch                     (:)     ! N acquired by retranslocation     (gN/m2/s)
     real(r8), pointer :: Nuptake_patch                             (:)     ! Total N uptake of FUN             (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_patch                  (:)     ! Total soil N uptake of FUN        (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_vr_patch               (:,:)   ! Total layer soil N uptake of FUN  (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_no3_vr_patch           (:,:)   ! Total layer no3 uptake of FUN     (gN/m2/s)
     real(r8), pointer :: sminn_to_plant_fun_nh4_vr_patch           (:,:)   ! Total layer nh4 uptake of FUN     (gN/m2/s)
     real(r8), pointer :: cost_nfix_patch                           (:)     ! Average cost of fixation          (gN/m2/s)
     real(r8), pointer :: cost_nactive_patch                        (:)     ! Average cost of active uptake     (gN/m2/s)
     real(r8), pointer :: cost_nretrans_patch                       (:)     ! Average cost of retranslocation   (gN/m2/s)
     real(r8), pointer :: nuptake_npp_fraction_patch                (:)     ! frac of npp spent on N acquisition   (gN/m2/s)
         ! Matrix
     real(r8), pointer :: matrix_nalloc_patch                       (:,:)   ! B-matrix for nitrogen allocation
     real(r8), pointer :: matrix_Ninput_patch                       (:)     ! I-matrix for nitrogen input
        
     real(r8), pointer :: matrix_nphtransfer_patch                  (:,:)   ! A-matrix_phenologh for nitrogen
     real(r8), pointer :: matrix_nphturnover_patch                  (:,:)   ! K-matrix_phenologh for nitrogen
     integer,  pointer :: matrix_nphtransfer_doner_patch            (:)     ! A-matrix_phenology non-zero indices (column indices) for nitrogen
     integer,  pointer :: matrix_nphtransfer_receiver_patch         (:)     ! A-matrix_phenology non-zero indices (row indices) for nitrogen

     real(r8), pointer :: matrix_ngmtransfer_patch                  (:,:)   ! A-matrix_gap mortality for nitrogen
     real(r8), pointer :: matrix_ngmturnover_patch                  (:,:)   ! K-matrix_gap mortality for nitrogen 
     integer,  pointer :: matrix_ngmtransfer_doner_patch            (:)     ! A-matrix_gap mortality non-zero indices (column indices) for nitrogen
     integer,  pointer :: matrix_ngmtransfer_receiver_patch         (:)     ! A-matrix_gap mortality non-zero indices (row indices) for nitrogen

     real(r8), pointer :: matrix_nfitransfer_patch                  (:,:)   ! A-matrix_fire for nitrogen
     real(r8), pointer :: matrix_nfiturnover_patch                  (:,:)   ! K-matrix_fire for nitrogen
     integer,  pointer :: matrix_nfitransfer_doner_patch            (:)     ! A-matrix_fire non-zero indices (column indices) for nitrogen
     integer,  pointer :: matrix_nfitransfer_receiver_patch         (:)     ! A-matrix_fire non-zero indices (row indices) for nitrogen

     integer ileafst_to_ileafxf_ph                    ! Index of phenology related N transfer from leaf storage pool to leaf transfer pool
     integer ileafxf_to_ileaf_ph                      ! Index of phenology related N transfer from leaf transfer pool to leaf pool  
     integer ifrootst_to_ifrootxf_ph                  ! Index of phenology related N transfer from fine root storage pool to fine root transfer pool
     integer ifrootxf_to_ifroot_ph                    ! Index of phenology related N transfer from fine root transfer pool to fine root pool  
     integer ilivestemst_to_ilivestemxf_ph            ! Index of phenology related N transfer from live stem storage pool to live stem transfer pool
     integer ilivestemxf_to_ilivestem_ph              ! Index of phenology related N transfer from live stem transfer pool to live stem pool  
     integer ideadstemst_to_ideadstemxf_ph            ! Index of phenology related N transfer from dead stem storage pool to dead stem transfer pool
     integer ideadstemxf_to_ideadstem_ph              ! Index of phenology related N transfer from dead stem transfer pool to dead stem pool  
     integer ilivecrootst_to_ilivecrootxf_ph          ! Index of phenology related N transfer from live coarse root storage pool to live coarse root transfer pool
     integer ilivecrootxf_to_ilivecroot_ph            ! Index of phenology related N transfer from live coarse root transfer pool to live coarse root pool  
     integer ideadcrootst_to_ideadcrootxf_ph          ! Index of phenology related N transfer from dead coarse root storage pool to dead coarse root transfer pool
     integer ideadcrootxf_to_ideadcroot_ph            ! Index of phenology related N transfer from dead coarse root transfer pool to dead coarse root pool  
     integer ilivestem_to_ideadstem_ph                ! Index of phenology related N transfer from live stem pool to dead stem pool  
     integer ilivecroot_to_ideadcroot_ph              ! Index of phenology related N transfer from live coarse root pool to dead coarse root pool  
     integer iretransn_to_ileaf_ph                    ! Index of phenology related N transfer from retranslocation pool to leaf pool
     integer iretransn_to_ileafst_ph                  ! Index of phenology related N transfer from retranslocation pool to leaf storage pool
     integer iretransn_to_ifroot_ph                   ! Index of phenology related N transfer from retranslocation pool to fine root pool
     integer iretransn_to_ifrootst_ph                 ! Index of phenology related N transfer from retranslocation pool to fine root storage pool
     integer iretransn_to_ilivestem_ph                ! Index of phenology related N transfer from retranslocation pool to live stem pool
     integer iretransn_to_ilivestemst_ph              ! Index of phenology related N transfer from retranslocation pool to live stem storage pool
     integer iretransn_to_ideadstem_ph                ! Index of phenology related N transfer from retranslocation pool to dead stem pool
     integer iretransn_to_ideadstemst_ph              ! Index of phenology related N transfer from retranslocation pool to dead stem storage pool
     integer iretransn_to_ilivecroot_ph               ! Index of phenology related N transfer from retranslocation pool to live coarse root pool
     integer iretransn_to_ilivecrootst_ph             ! Index of phenology related N transfer from retranslocation pool to live coarse root storage pool
     integer iretransn_to_ideadcroot_ph               ! Index of phenology related N transfer from retranslocation pool to dead coarse root pool
     integer iretransn_to_ideadcrootst_ph             ! Index of phenology related N transfer from retranslocation pool to dead coarse root storage pool
     integer iretransn_to_igrain_ph                   ! Index of phenology related N transfer from retranslocation pool to grain pool
     integer iretransn_to_igrainst_ph                 ! Index of phenology related N transfer from retranslocation pool to grain storage pool
     integer ileaf_to_iout_ph                         ! Index of phenology related N transfer from leaf pool to outside of vegetation pools  
     integer ifroot_to_iout_ph                        ! Index of phenology related N transfer from fine root pool to outside of vegetation pools  
     integer ilivestem_to_iout_ph                     ! Index of phenology related N transfer from live stem pool to outside of vegetation pools  
     integer ileaf_to_iretransn_ph                    ! Index of phenology related N transfer from leaf pool to retranslocation pools
     integer ifroot_to_iretransn_ph                   ! Index of phenology related N transfer from fine root pool to retranslocation pools
     integer ilivestem_to_iretransn_ph                ! Index of phenology related N transfer from live stem pool to retranslocation pools
     integer ilivecroot_to_iretransn_ph               ! Index of phenology related N transfer from live coarse root pool to retranslocation pools
     integer igrain_to_iout_ph                        ! Index of phenology related N transfer from grain pool to outside of vegetation pools  
     integer iretransn_to_iout_ph                     ! Index of phenology related N transfer from retranslocation pool to outside of vegetation pools
     integer ileaf_to_iout_gm                         ! Index of gap mortality related N transfer from leaf pool to outside of vegetation pools
     integer ileafst_to_iout_gm                       ! Index of gap mortality related N transfer from leaf storage pool to outside of vegetation pools
     integer ileafxf_to_iout_gm                       ! Index of gap mortality related N transfer from leaf transfer pool to outside of vegetation pools
     integer ifroot_to_iout_gm                        ! Index of gap mortality related N transfer from fine root pool to outside of vegetation pools
     integer ifrootst_to_iout_gm                      ! Index of gap mortality related N transfer from fine root storage pool to outside of vegetation pools
     integer ifrootxf_to_iout_gm                      ! Index of gap mortality related N transfer from fine root transfer pool to outside of vegetation pools
     integer ilivestem_to_iout_gm                     ! Index of gap mortality related N transfer from live stem pool to outside of vegetation pools
     integer ilivestemst_to_iout_gm                   ! Index of gap mortality related N transfer from live stem storage pool to outside of vegetation pools
     integer ilivestemxf_to_iout_gm                   ! Index of gap mortality related N transfer from live stem transfer pool to outside of vegetation pools
     integer ideadstem_to_iout_gm                     ! Index of gap mortality related N transfer from dead stem pool to outside of vegetation pools
     integer ideadstemst_to_iout_gm                   ! Index of gap mortality related N transfer from dead stem storage pool to outside of vegetation pools
     integer ideadstemxf_to_iout_gm                   ! Index of gap mortality related N transfer from dead stem transfer pool to outside of vegetation pools
     integer ilivecroot_to_iout_gm                    ! Index of gap mortality related N transfer from live coarse root pool to outside of vegetation pools
     integer ilivecrootst_to_iout_gm                  ! Index of gap mortality related N transfer from live coarse root storage pool to outside of vegetation pools
     integer ilivecrootxf_to_iout_gm                  ! Index of gap mortality related N transfer from live coarse root transfer pool to outside of vegetation pools
     integer ideadcroot_to_iout_gm                    ! Index of gap mortality related N transfer from dead coarse root pool to outside of vegetation pools
     integer ideadcrootst_to_iout_gm                  ! Index of gap mortality related N transfer from dead coarse root storage pool to outside of vegetation pools
     integer ideadcrootxf_to_iout_gm                  ! Index of gap mortality related N transfer from dead coarse root transfer pool to outside of vegetation pools
     integer iretransn_to_iout_gm                     ! Index of gap mortality related N transfer from retranslocation to outside of vegetation pools
     integer ileaf_to_iout_fi                         ! Index of fire related N transfer from leaf pool to outside of vegetation pools
     integer ileafst_to_iout_fi                       ! Index of fire related N transfer from leaf storage pool to outside of vegetation pools
     integer ileafxf_to_iout_fi                       ! Index of fire related N transfer from leaf transfer pool to outside of vegetation pools
     integer ifroot_to_iout_fi                        ! Index of fire related N transfer from fine root pool to outside of vegetation pools
     integer ifrootst_to_iout_fi                      ! Index of fire related N transfer from fine root storage pool to outside of vegetation pools
     integer ifrootxf_to_iout_fi                      ! Index of fire related N transfer from fine root transfer pool to outside of vegetation pools
     integer ilivestem_to_iout_fi                     ! Index of fire related N transfer from live stem pool to outside of vegetation pools
     integer ilivestemst_to_iout_fi                   ! Index of fire related N transfer from live stem storage pool to outside of vegetation pools
     integer ilivestemxf_to_iout_fi                   ! Index of fire related N transfer from live stem transfer pool to outside of vegetation pools
     integer ideadstem_to_iout_fi                     ! Index of fire related N transfer from dead stem pool to outside of vegetation pools
     integer ideadstemst_to_iout_fi                   ! Index of fire related N transfer from dead stem storage pool to outside of vegetation pools
     integer ideadstemxf_to_iout_fi                   ! Index of fire related N transfer from dead stem transfer pool to outside of vegetation pools
     integer ilivecroot_to_iout_fi                    ! Index of fire related N transfer from live coarse root pool to outside of vegetation pools
     integer ilivecrootst_to_iout_fi                  ! Index of fire related N transfer from live coarse root storage pool to outside of vegetation pools
     integer ilivecrootxf_to_iout_fi                  ! Index of fire related N transfer from live coarse root transfer pool to outside of vegetation pools
     integer ideadcroot_to_iout_fi                    ! Index of fire related N transfer from dead coarse root pool to outside of vegetation pools
     integer ideadcrootst_to_iout_fi                  ! Index of fire related N transfer from dead coarse root storage pool to outside of vegetation pools
     integer ideadcrootxf_to_iout_fi                  ! Index of fire related N transfer from dead coarse root transfer pool to outside of vegetation pools
     integer iretransn_to_iout_fi                     ! Index of fire related N transfer from retranslocation transfer pool to outside of vegetation pools
     integer ilivestem_to_ideadstem_fi                ! Index of fire related N transfer from live stem pool to dead stem pools
     integer ilivecroot_to_ideadcroot_fi              ! Index of fire related N transfer from live coarse root pool to dead coarse root pools

     integer,pointer :: list_phn_phgmn     (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKphn to AKphn+AKgmn
     integer,pointer :: list_gmn_phgmn     (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKgmn to AKphn+AKgmn
     integer,pointer :: list_phn_phgmfin   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKphn to AKphn+AKgmn+AKfin
     integer,pointer :: list_gmn_phgmfin   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKgmn to AKphn+AKgmn+AKfin
     integer,pointer :: list_fin_phgmfin   (:)        ! Index mapping for sparse matrix addition (save to reduce computational cost): from AKfin to AKphn+AKgmn+AKfin
     integer,pointer :: list_aphn          (:)        ! Indices of non-diagnoal entries in full sparse matrix Aph for N cycle
     integer,pointer :: list_agmn          (:)        ! Indices of non-diagnoal entries in full sparse matrix Agm for N cycle
     integer,pointer :: list_afin          (:)        ! Indices of non-diagnoal entries in full sparse matrix Afi for N cycle

 contains

     procedure , public  :: SetValues
     procedure , public  :: Summary => Summary_nitrogenflux
     procedure , public  :: ZeroDWT
     procedure , public  :: Init

 end type cnveg_nitrogenflux_type

type(cnveg_nitrogenflux_type), public, target, save :: cnveg_nitrogenflux_inst

contains

!---------------------------------------
 subroutine Init(this, bounds, nch, ityp, fveg, cncol, cnpft)

! !DESCRIPTION:
! Initialize CTSM nitrogen fluxes
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
    class(cnveg_nitrogenflux_type)                           :: this

    ! LOCAL
    integer  :: begp, endp
    integer  :: begc, endc
    integer  :: begg, endg
    integer  :: np, nc, nz, p, nv, n, l, j
    !--------------------------------

    allocate(this%matrix_nphtransfer_doner_patch(1:37)) 
    allocate(this%matrix_nphtransfer_receiver_patch(1:37))   

    this%ileaf_to_iretransn_ph           = 1
    this%matrix_nphtransfer_doner_patch(this%ileaf_to_iretransn_ph)              = ileaf
    this%matrix_nphtransfer_receiver_patch(this%ileaf_to_iretransn_ph)           = iretransn

    this%ileafst_to_ileafxf_ph           = 2
    this%matrix_nphtransfer_doner_patch(this%ileafst_to_ileafxf_ph)              = ileaf_st
    this%matrix_nphtransfer_receiver_patch(this%ileafst_to_ileafxf_ph)           = ileaf_xf

    this%ileafxf_to_ileaf_ph             = 3
    this%matrix_nphtransfer_doner_patch(this%ileafxf_to_ileaf_ph)                = ileaf_xf
    this%matrix_nphtransfer_receiver_patch(this%ileafxf_to_ileaf_ph)             = ileaf

    this%ifroot_to_iretransn_ph          = 4
    this%matrix_nphtransfer_doner_patch(this%ifroot_to_iretransn_ph)             = ifroot
    this%matrix_nphtransfer_receiver_patch(this%ifroot_to_iretransn_ph)          = iretransn

    this%ifrootst_to_ifrootxf_ph         = 5
    this%matrix_nphtransfer_doner_patch(this%ifrootst_to_ifrootxf_ph)            = ifroot_st
    this%matrix_nphtransfer_receiver_patch(this%ifrootst_to_ifrootxf_ph)         = ifroot_xf

    this%ifrootxf_to_ifroot_ph           = 6
    this%matrix_nphtransfer_doner_patch(this%ifrootxf_to_ifroot_ph)              = ifroot_xf
    this%matrix_nphtransfer_receiver_patch(this%ifrootxf_to_ifroot_ph)           = ifroot

    this%ilivestem_to_ideadstem_ph       = 7
    this%matrix_nphtransfer_doner_patch(this%ilivestem_to_ideadstem_ph)          = ilivestem
    this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_ideadstem_ph)       = ideadstem

    this%ilivestem_to_iretransn_ph       = 8
    this%matrix_nphtransfer_doner_patch(this%ilivestem_to_iretransn_ph)          = ilivestem
    this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_iretransn_ph)       = iretransn

    this%ilivestemst_to_ilivestemxf_ph   = 9
    this%matrix_nphtransfer_doner_patch(this%ilivestemst_to_ilivestemxf_ph)      = ilivestem_st
    this%matrix_nphtransfer_receiver_patch(this%ilivestemst_to_ilivestemxf_ph)   = ilivestem_xf

    this%ilivestemxf_to_ilivestem_ph     = 10
    this%matrix_nphtransfer_doner_patch(this%ilivestemxf_to_ilivestem_ph)        = ilivestem_xf
    this%matrix_nphtransfer_receiver_patch(this%ilivestemxf_to_ilivestem_ph)     = ilivestem

    this%ideadstemst_to_ideadstemxf_ph   = 11
    this%matrix_nphtransfer_doner_patch(this%ideadstemst_to_ideadstemxf_ph)      = ideadstem_st
    this%matrix_nphtransfer_receiver_patch(this%ideadstemst_to_ideadstemxf_ph)   = ideadstem_xf

    this%ideadstemxf_to_ideadstem_ph     = 12
    this%matrix_nphtransfer_doner_patch(this%ideadstemxf_to_ideadstem_ph)        = ideadstem_xf
    this%matrix_nphtransfer_receiver_patch(this%ideadstemxf_to_ideadstem_ph)     = ideadstem

    this%ilivecroot_to_ideadcroot_ph     = 13
    this%matrix_nphtransfer_doner_patch(this%ilivecroot_to_ideadcroot_ph)        = ilivecroot
    this%matrix_nphtransfer_receiver_patch(this%ilivecroot_to_ideadcroot_ph)     = ideadcroot

    this%ilivecroot_to_iretransn_ph      = 14
    this%matrix_nphtransfer_doner_patch(this%ilivecroot_to_iretransn_ph)         = ilivecroot
    this%matrix_nphtransfer_receiver_patch(this%ilivecroot_to_iretransn_ph)      = iretransn

    this%ilivecrootst_to_ilivecrootxf_ph = 15
    this%matrix_nphtransfer_doner_patch(this%ilivecrootst_to_ilivecrootxf_ph)    = ilivecroot_st
    this%matrix_nphtransfer_receiver_patch(this%ilivecrootst_to_ilivecrootxf_ph) = ilivecroot_xf

    this%ilivecrootxf_to_ilivecroot_ph   = 16
    this%matrix_nphtransfer_doner_patch(this%ilivecrootxf_to_ilivecroot_ph)      = ilivecroot_xf
    this%matrix_nphtransfer_receiver_patch(this%ilivecrootxf_to_ilivecroot_ph)   = ilivecroot

    this%ideadcrootst_to_ideadcrootxf_ph = 17
    this%matrix_nphtransfer_doner_patch(this%ideadcrootst_to_ideadcrootxf_ph)    = ideadcroot_st
    this%matrix_nphtransfer_receiver_patch(this%ideadcrootst_to_ideadcrootxf_ph) = ideadcroot_xf

    this%ideadcrootxf_to_ideadcroot_ph   = 18
    this%matrix_nphtransfer_doner_patch(this%ideadcrootxf_to_ideadcroot_ph)      = ideadcroot_xf
    this%matrix_nphtransfer_receiver_patch(this%ideadcrootxf_to_ideadcroot_ph)   = ideadcroot

    this%iretransn_to_ileaf_ph           = 19
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ileaf_ph)              = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ileaf_ph)           = ileaf

    this%iretransn_to_ileafst_ph         = 20
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ileafst_ph)            = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ileafst_ph)         = ileaf_st

    this%iretransn_to_ifroot_ph          = 21
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ifroot_ph)             = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ifroot_ph)          = ifroot

    this%iretransn_to_ifrootst_ph        = 22
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ifrootst_ph)           = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ifrootst_ph)        = ifroot_st
    this%iretransn_to_ilivestem_ph       = 23
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivestem_ph)          = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivestem_ph)       = ilivestem

    this%iretransn_to_ilivestemst_ph     = 24
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivestemst_ph)        = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivestemst_ph)     = ilivestem_st

    this%iretransn_to_ideadstem_ph       = 25
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadstem_ph)          = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadstem_ph)       = ideadstem

    this%iretransn_to_ideadstemst_ph     = 26
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadstemst_ph)        = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadstemst_ph)     = ideadstem_st

    this%iretransn_to_ilivecroot_ph      = 27
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivecroot_ph)         = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivecroot_ph)      = ilivecroot

    this%iretransn_to_ilivecrootst_ph    = 28
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ilivecrootst_ph)       = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ilivecrootst_ph)    = ilivecroot_st

    this%iretransn_to_ideadcroot_ph      = 29
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadcroot_ph)         = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadcroot_ph)      = ideadcroot

    this%iretransn_to_ideadcrootst_ph    = 30
    this%matrix_nphtransfer_doner_patch(this%iretransn_to_ideadcrootst_ph)       = iretransn
    this%matrix_nphtransfer_receiver_patch(this%iretransn_to_ideadcrootst_ph)    = ideadcroot_st

    if(.not. use_crop)then
       this%ileaf_to_iout_ph             = 31
       this%matrix_nphtransfer_doner_patch(this%ileaf_to_iout_ph)                = ileaf
       this%matrix_nphtransfer_receiver_patch(this%ileaf_to_iout_ph)             = ioutn

       this%ifroot_to_iout_ph            = 32
       this%matrix_nphtransfer_doner_patch(this%ifroot_to_iout_ph)               = ifroot
       this%matrix_nphtransfer_receiver_patch(this%ifroot_to_iout_ph)            = ioutn

       this%ilivestem_to_iout_ph         = 33
       this%matrix_nphtransfer_doner_patch(this%ilivestem_to_iout_ph)            = ilivestem
       this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_iout_ph)         = ioutn

       this%iretransn_to_iout_ph         = 34
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_iout_ph)            = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_iout_ph)         = ioutn
    else
       this%iretransn_to_igrain_ph       = 31
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_igrain_ph)          = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_igrain_ph)       = igrain

       this%iretransn_to_igrainst_ph     = 32
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_igrainst_ph)        = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_igrainst_ph)     = igrain_st

       this%ileaf_to_iout_ph             = 33
       this%matrix_nphtransfer_doner_patch(this%ileaf_to_iout_ph)                = ileaf
       this%matrix_nphtransfer_receiver_patch(this%ileaf_to_iout_ph)             = ioutn

       this%ifroot_to_iout_ph            = 34
       this%matrix_nphtransfer_doner_patch(this%ifroot_to_iout_ph)               = ifroot
       this%matrix_nphtransfer_receiver_patch(this%ifroot_to_iout_ph)            = ioutn

       this%ilivestem_to_iout_ph         = 35
       this%matrix_nphtransfer_doner_patch(this%ilivestem_to_iout_ph)            = ilivestem
       this%matrix_nphtransfer_receiver_patch(this%ilivestem_to_iout_ph)         = ioutn

       this%igrain_to_iout_ph            = 36
       this%matrix_nphtransfer_doner_patch(this%igrain_to_iout_ph)               = igrain
       this%matrix_nphtransfer_receiver_patch(this%igrain_to_iout_ph)            = ioutn

       this%iretransn_to_iout_ph         = 37
       this%matrix_nphtransfer_doner_patch(this%iretransn_to_iout_ph)            = iretransn
       this%matrix_nphtransfer_receiver_patch(this%iretransn_to_iout_ph)         = ioutn
    end if

    allocate(this%matrix_ngmtransfer_doner_patch(1:19))
    allocate(this%matrix_ngmtransfer_receiver_patch(1:19))

    this%ileaf_to_iout_gm                = 1
    this%matrix_ngmtransfer_doner_patch(this%ileaf_to_iout_gm)                   = ileaf
    this%matrix_ngmtransfer_receiver_patch(this%ileaf_to_iout_gm)                = ioutn

    this%ileafst_to_iout_gm              = 2
    this%matrix_ngmtransfer_doner_patch(this%ileafst_to_iout_gm)                 = ileaf_st
    this%matrix_ngmtransfer_receiver_patch(this%ileafst_to_iout_gm)              = ioutn

    this%ileafxf_to_iout_gm              = 3
    this%matrix_ngmtransfer_doner_patch(this%ileafxf_to_iout_gm)                 = ileaf_xf
    this%matrix_ngmtransfer_receiver_patch(this%ileafxf_to_iout_gm)              = ioutn

    this%ifroot_to_iout_gm               = 4
    this%matrix_ngmtransfer_doner_patch(this%ifroot_to_iout_gm)                  = ifroot
    this%matrix_ngmtransfer_receiver_patch(this%ifroot_to_iout_gm)               = ioutn

    this%ifrootst_to_iout_gm             = 5
    this%matrix_ngmtransfer_doner_patch(this%ifrootst_to_iout_gm)                = ifroot_st
    this%matrix_ngmtransfer_receiver_patch(this%ifrootst_to_iout_gm)             = ioutn

    this%ifrootxf_to_iout_gm             = 6
    this%matrix_ngmtransfer_doner_patch(this%ifrootxf_to_iout_gm)                = ifroot_xf
    this%matrix_ngmtransfer_receiver_patch(this%ifrootxf_to_iout_gm)             = ioutn

    this%ilivestem_to_iout_gm            = 7
    this%matrix_ngmtransfer_doner_patch(this%ilivestem_to_iout_gm)               = ilivestem
    this%matrix_ngmtransfer_receiver_patch(this%ilivestem_to_iout_gm)            = ioutn

    this%ilivestemst_to_iout_gm          = 8
    this%matrix_ngmtransfer_doner_patch(this%ilivestemst_to_iout_gm)             = ilivestem_st
    this%matrix_ngmtransfer_receiver_patch(this%ilivestemst_to_iout_gm)          = ioutn

    this%ilivestemxf_to_iout_gm          = 9
    this%matrix_ngmtransfer_doner_patch(this%ilivestemxf_to_iout_gm)             = ilivestem_xf
    this%matrix_ngmtransfer_receiver_patch(this%ilivestemxf_to_iout_gm)          = ioutn

    this%ideadstem_to_iout_gm            = 10
    this%matrix_ngmtransfer_doner_patch(this%ideadstem_to_iout_gm)               = ideadstem
    this%matrix_ngmtransfer_receiver_patch(this%ideadstem_to_iout_gm)            = ioutn

    this%ideadstemst_to_iout_gm          = 11
    this%matrix_ngmtransfer_doner_patch(this%ideadstemst_to_iout_gm)             = ideadstem_st
    this%matrix_ngmtransfer_receiver_patch(this%ideadstemst_to_iout_gm)          = ioutn

    this%ideadstemxf_to_iout_gm          = 12
    this%matrix_ngmtransfer_doner_patch(this%ideadstemxf_to_iout_gm)             = ideadstem_xf
    this%matrix_ngmtransfer_receiver_patch(this%ideadstemxf_to_iout_gm)          = ioutn

    this%ilivecroot_to_iout_gm           = 13
    this%matrix_ngmtransfer_doner_patch(this%ilivecroot_to_iout_gm)              = ilivecroot
    this%matrix_ngmtransfer_receiver_patch(this%ilivecroot_to_iout_gm)           = ioutn

    this%ilivecrootst_to_iout_gm         = 14
    this%matrix_ngmtransfer_doner_patch(this%ilivecrootst_to_iout_gm)            = ilivecroot_st
    this%matrix_ngmtransfer_receiver_patch(this%ilivecrootst_to_iout_gm)         = ioutn


    this%ilivecrootxf_to_iout_gm         = 15
    this%matrix_ngmtransfer_doner_patch(this%ilivecrootxf_to_iout_gm)            = ilivecroot_xf
    this%matrix_ngmtransfer_receiver_patch(this%ilivecrootxf_to_iout_gm)         = ioutn

    this%ideadcroot_to_iout_gm           = 16
    this%matrix_ngmtransfer_doner_patch(this%ideadcroot_to_iout_gm)              = ideadcroot
    this%matrix_ngmtransfer_receiver_patch(this%ideadcroot_to_iout_gm)           = ioutn

    this%ideadcrootst_to_iout_gm         = 17
    this%matrix_ngmtransfer_doner_patch(this%ideadcrootst_to_iout_gm)            = ideadcroot_st
    this%matrix_ngmtransfer_receiver_patch(this%ideadcrootst_to_iout_gm)         = ioutn

    this%ideadcrootxf_to_iout_gm         = 18
    this%matrix_ngmtransfer_doner_patch(this%ideadcrootxf_to_iout_gm)            = ideadcroot_xf
    this%matrix_ngmtransfer_receiver_patch(this%ideadcrootxf_to_iout_gm)         = ioutn

    this%iretransn_to_iout_gm            = 19
    this%matrix_ngmtransfer_doner_patch(this%iretransn_to_iout_gm)               = iretransn
    this%matrix_ngmtransfer_receiver_patch(this%iretransn_to_iout_gm)            = ioutn

    allocate(this%matrix_nfitransfer_doner_patch(1:21))
    allocate(this%matrix_nfitransfer_receiver_patch(1:21))

    this%ilivestem_to_ideadstem_fi       = 1
    this%matrix_nfitransfer_doner_patch(this%ilivestem_to_ideadstem_fi)          = ilivestem
    this%matrix_nfitransfer_receiver_patch(this%ilivestem_to_ideadstem_fi)       = ideadstem

    this%ilivecroot_to_ideadcroot_fi     = 2
    this%matrix_nfitransfer_doner_patch(this%ilivecroot_to_ideadcroot_fi)        = ilivecroot
    this%matrix_nfitransfer_receiver_patch(this%ilivecroot_to_ideadcroot_fi)     = ideadcroot

    this%ileaf_to_iout_fi                = 3
    this%matrix_nfitransfer_doner_patch(this%ileaf_to_iout_fi)                   = ileaf
    this%matrix_nfitransfer_receiver_patch(this%ileaf_to_iout_fi)                = ioutn

    this%ileafst_to_iout_fi              = 4
    this%matrix_nfitransfer_doner_patch(this%ileafst_to_iout_fi)                 = ileaf_st
    this%matrix_nfitransfer_receiver_patch(this%ileafst_to_iout_fi)              = ioutn

    this%ileafxf_to_iout_fi              = 5
    this%matrix_nfitransfer_doner_patch(this%ileafxf_to_iout_fi)                 = ileaf_xf
    this%matrix_nfitransfer_receiver_patch(this%ileafxf_to_iout_fi)              = ioutn

    this%ifroot_to_iout_fi               = 6
    this%matrix_nfitransfer_doner_patch(this%ifroot_to_iout_fi)                  = ifroot
    this%matrix_nfitransfer_receiver_patch(this%ifroot_to_iout_fi)               = ioutn

    this%ifrootst_to_iout_fi             = 7
    this%matrix_nfitransfer_doner_patch(this%ifrootst_to_iout_fi)                = ifroot_st
    this%matrix_nfitransfer_receiver_patch(this%ifrootst_to_iout_fi)             = ioutn

    this%ifrootxf_to_iout_fi             = 8
    this%matrix_nfitransfer_doner_patch(this%ifrootxf_to_iout_fi)                = ifroot_xf
    this%matrix_nfitransfer_receiver_patch(this%ifrootxf_to_iout_fi)             = ioutn

    this%ilivestem_to_iout_fi            = 9
    this%matrix_nfitransfer_doner_patch(this%ilivestem_to_iout_fi)               = ilivestem
    this%matrix_nfitransfer_receiver_patch(this%ilivestem_to_iout_fi)            = ioutn

    this%ilivestemst_to_iout_fi          = 10
    this%matrix_nfitransfer_doner_patch(this%ilivestemst_to_iout_fi)             = ilivestem_st
    this%matrix_nfitransfer_receiver_patch(this%ilivestemst_to_iout_fi)          = ioutn

    this%ilivestemxf_to_iout_fi          = 11
    this%matrix_nfitransfer_doner_patch(this%ilivestemxf_to_iout_fi)             = ilivestem_xf
    this%matrix_nfitransfer_receiver_patch(this%ilivestemxf_to_iout_fi)          = ioutn

    this%ideadstem_to_iout_fi            = 12
    this%matrix_nfitransfer_doner_patch(this%ideadstem_to_iout_fi)               = ideadstem
    this%matrix_nfitransfer_receiver_patch(this%ideadstem_to_iout_fi)            = ioutn

    this%ideadstemst_to_iout_fi          = 13
    this%matrix_nfitransfer_doner_patch(this%ideadstemst_to_iout_fi)             = ideadstem_st
    this%matrix_nfitransfer_receiver_patch(this%ideadstemst_to_iout_fi)          = ioutn

    this%ideadstemxf_to_iout_fi          = 14
    this%matrix_nfitransfer_doner_patch(this%ideadstemxf_to_iout_fi)             = ideadstem_xf
    this%matrix_nfitransfer_receiver_patch(this%ideadstemxf_to_iout_fi)          = ioutn

    this%ilivecroot_to_iout_fi           = 15
    this%matrix_nfitransfer_doner_patch(this%ilivecroot_to_iout_fi)              = ilivecroot
    this%matrix_nfitransfer_receiver_patch(this%ilivecroot_to_iout_fi)           = ioutn

    this%ilivecrootst_to_iout_fi         = 16
    this%matrix_nfitransfer_doner_patch(this%ilivecrootst_to_iout_fi)            = ilivecroot_st
    this%matrix_nfitransfer_receiver_patch(this%ilivecrootst_to_iout_fi)         = ioutn

    this%ilivecrootxf_to_iout_fi         = 17
    this%matrix_nfitransfer_doner_patch(this%ilivecrootxf_to_iout_fi)            = ilivecroot_xf
    this%matrix_nfitransfer_receiver_patch(this%ilivecrootxf_to_iout_fi)         = ioutn


    this%ideadcroot_to_iout_fi           = 18
    this%matrix_nfitransfer_doner_patch(this%ideadcroot_to_iout_fi)              = ideadcroot
    this%matrix_nfitransfer_receiver_patch(this%ideadcroot_to_iout_fi)           = ioutn

    this%ideadcrootst_to_iout_fi         = 19
    this%matrix_nfitransfer_doner_patch(this%ideadcrootst_to_iout_fi)            = ideadcroot_st
    this%matrix_nfitransfer_receiver_patch(this%ideadcrootst_to_iout_fi)         = ioutn

    this%ideadcrootxf_to_iout_fi         = 20
    this%matrix_nfitransfer_doner_patch(this%ideadcrootxf_to_iout_fi)            = ideadcroot_xf
    this%matrix_nfitransfer_receiver_patch(this%ideadcrootxf_to_iout_fi)         = ioutn

    this%iretransn_to_iout_fi            = 21
    this%matrix_nfitransfer_doner_patch(this%iretransn_to_iout_fi)               = iretransn
    this%matrix_nfitransfer_receiver_patch(this%iretransn_to_iout_fi)            = ioutn

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg = bounds%endg

    allocate(this%m_leafn_to_litter_patch                   (begp:endp)) ; this%m_leafn_to_litter_patch                   (:) = nan
    allocate(this%m_frootn_to_litter_patch                  (begp:endp)) ; this%m_frootn_to_litter_patch                  (:) = nan
    allocate(this%m_leafn_storage_to_litter_patch           (begp:endp)) ; this%m_leafn_storage_to_litter_patch           (:) = nan
    allocate(this%m_frootn_storage_to_litter_patch          (begp:endp)) ; this%m_frootn_storage_to_litter_patch          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_patch       (begp:endp)) ; this%m_livestemn_storage_to_litter_patch       (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_patch       (begp:endp)) ; this%m_deadstemn_storage_to_litter_patch       (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_patch      (begp:endp)) ; this%m_livecrootn_storage_to_litter_patch      (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_patch      (begp:endp)) ; this%m_deadcrootn_storage_to_litter_patch      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_patch              (begp:endp)) ; this%m_leafn_xfer_to_litter_patch              (:) = nan
    allocate(this%m_frootn_xfer_to_litter_patch             (begp:endp)) ; this%m_frootn_xfer_to_litter_patch             (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_patch          (begp:endp)) ; this%m_livestemn_xfer_to_litter_patch          (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_patch          (begp:endp)) ; this%m_deadstemn_xfer_to_litter_patch          (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_patch         (begp:endp)) ; this%m_livecrootn_xfer_to_litter_patch         (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_patch         (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_patch         (:) = nan
    allocate(this%m_livestemn_to_litter_patch               (begp:endp)) ; this%m_livestemn_to_litter_patch               (:) = nan
    allocate(this%m_deadstemn_to_litter_patch               (begp:endp)) ; this%m_deadstemn_to_litter_patch               (:) = nan
    allocate(this%m_livecrootn_to_litter_patch              (begp:endp)) ; this%m_livecrootn_to_litter_patch              (:) = nan
    allocate(this%m_deadcrootn_to_litter_patch              (begp:endp)) ; this%m_deadcrootn_to_litter_patch              (:) = nan
    allocate(this%m_retransn_to_litter_patch                (begp:endp)) ; this%m_retransn_to_litter_patch                (:) = nan
    allocate(this%hrv_leafn_to_litter_patch                 (begp:endp)) ; this%hrv_leafn_to_litter_patch                 (:) = nan
    allocate(this%hrv_frootn_to_litter_patch                (begp:endp)) ; this%hrv_frootn_to_litter_patch                (:) = nan
    allocate(this%hrv_leafn_storage_to_litter_patch         (begp:endp)) ; this%hrv_leafn_storage_to_litter_patch         (:) = nan
    allocate(this%hrv_frootn_storage_to_litter_patch        (begp:endp)) ; this%hrv_frootn_storage_to_litter_patch        (:) = nan
    allocate(this%hrv_livestemn_storage_to_litter_patch     (begp:endp)) ; this%hrv_livestemn_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_deadstemn_storage_to_litter_patch     (begp:endp)) ; this%hrv_deadstemn_storage_to_litter_patch     (:) = nan
    allocate(this%hrv_livecrootn_storage_to_litter_patch    (begp:endp)) ; this%hrv_livecrootn_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_deadcrootn_storage_to_litter_patch    (begp:endp)) ; this%hrv_deadcrootn_storage_to_litter_patch    (:) = nan
    allocate(this%hrv_leafn_xfer_to_litter_patch            (begp:endp)) ; this%hrv_leafn_xfer_to_litter_patch            (:) = nan
    allocate(this%hrv_frootn_xfer_to_litter_patch           (begp:endp)) ; this%hrv_frootn_xfer_to_litter_patch           (:) = nan
    allocate(this%hrv_livestemn_xfer_to_litter_patch        (begp:endp)) ; this%hrv_livestemn_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_deadstemn_xfer_to_litter_patch        (begp:endp)) ; this%hrv_deadstemn_xfer_to_litter_patch        (:) = nan
    allocate(this%hrv_livecrootn_xfer_to_litter_patch       (begp:endp)) ; this%hrv_livecrootn_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_deadcrootn_xfer_to_litter_patch       (begp:endp)) ; this%hrv_deadcrootn_xfer_to_litter_patch       (:) = nan
    allocate(this%hrv_livestemn_to_litter_patch             (begp:endp)) ; this%hrv_livestemn_to_litter_patch             (:) = nan
    allocate(this%hrv_livecrootn_to_litter_patch            (begp:endp)) ; this%hrv_livecrootn_to_litter_patch            (:) = nan
    allocate(this%hrv_deadcrootn_to_litter_patch            (begp:endp)) ; this%hrv_deadcrootn_to_litter_patch            (:) = nan
    allocate(this%hrv_retransn_to_litter_patch              (begp:endp)) ; this%hrv_retransn_to_litter_patch              (:) = nan

    allocate(this%m_leafn_to_fire_patch                     (begp:endp)) ; this%m_leafn_to_fire_patch                     (:) = nan
    allocate(this%m_leafn_storage_to_fire_patch             (begp:endp)) ; this%m_leafn_storage_to_fire_patch             (:) = nan
    allocate(this%m_leafn_xfer_to_fire_patch                (begp:endp)) ; this%m_leafn_xfer_to_fire_patch                (:) = nan
    allocate(this%m_livestemn_to_fire_patch                 (begp:endp)) ; this%m_livestemn_to_fire_patch                 (:) = nan
    allocate(this%m_livestemn_storage_to_fire_patch         (begp:endp)) ; this%m_livestemn_storage_to_fire_patch         (:) = nan
    allocate(this%m_livestemn_xfer_to_fire_patch            (begp:endp)) ; this%m_livestemn_xfer_to_fire_patch            (:) = nan
    allocate(this%m_deadstemn_to_fire_patch                 (begp:endp)) ; this%m_deadstemn_to_fire_patch                 (:) = nan
    allocate(this%m_deadstemn_storage_to_fire_patch         (begp:endp)) ; this%m_deadstemn_storage_to_fire_patch         (:) = nan
    allocate(this%m_deadstemn_xfer_to_fire_patch            (begp:endp)) ; this%m_deadstemn_xfer_to_fire_patch            (:) = nan
    allocate(this%m_frootn_to_fire_patch                    (begp:endp)) ; this%m_frootn_to_fire_patch                    (:) = nan
    allocate(this%m_frootn_storage_to_fire_patch            (begp:endp)) ; this%m_frootn_storage_to_fire_patch            (:) = nan
    allocate(this%m_frootn_xfer_to_fire_patch               (begp:endp)) ; this%m_frootn_xfer_to_fire_patch               (:) = nan
    allocate(this%m_livecrootn_to_fire_patch                (begp:endp)) ;
    allocate(this%m_livecrootn_storage_to_fire_patch        (begp:endp)) ; this%m_livecrootn_storage_to_fire_patch        (:) = nan
    allocate(this%m_livecrootn_xfer_to_fire_patch           (begp:endp)) ; this%m_livecrootn_xfer_to_fire_patch           (:) = nan
    allocate(this%m_deadcrootn_to_fire_patch                (begp:endp)) ; this%m_deadcrootn_to_fire_patch                (:) = nan
    allocate(this%m_deadcrootn_storage_to_fire_patch        (begp:endp)) ; this%m_deadcrootn_storage_to_fire_patch        (:) = nan
    allocate(this%m_deadcrootn_xfer_to_fire_patch           (begp:endp)) ; this%m_deadcrootn_xfer_to_fire_patch           (:) = nan
    allocate(this%m_retransn_to_fire_patch                  (begp:endp)) ; this%m_retransn_to_fire_patch                  (:) = nan

    allocate(this%m_leafn_to_litter_fire_patch              (begp:endp)) ; this%m_leafn_to_litter_fire_patch              (:) = nan
    allocate(this%m_leafn_storage_to_litter_fire_patch      (begp:endp)) ; this%m_leafn_storage_to_litter_fire_patch      (:) = nan
    allocate(this%m_leafn_xfer_to_litter_fire_patch         (begp:endp)) ; this%m_leafn_xfer_to_litter_fire_patch         (:) = nan
    allocate(this%m_livestemn_to_litter_fire_patch          (begp:endp)) ; this%m_livestemn_to_litter_fire_patch          (:) = nan
    allocate(this%m_livestemn_storage_to_litter_fire_patch  (begp:endp)) ; this%m_livestemn_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_livestemn_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_livestemn_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_livestemn_to_deadstemn_fire_patch       (begp:endp)) ; this%m_livestemn_to_deadstemn_fire_patch       (:) = nan
    allocate(this%m_deadstemn_to_litter_fire_patch          (begp:endp)) ; this%m_deadstemn_to_litter_fire_patch          (:) = nan
    allocate(this%m_deadstemn_storage_to_litter_fire_patch  (begp:endp)) ; this%m_deadstemn_storage_to_litter_fire_patch  (:) = nan
    allocate(this%m_deadstemn_xfer_to_litter_fire_patch     (begp:endp)) ; this%m_deadstemn_xfer_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootn_to_litter_fire_patch             (begp:endp)) ; this%m_frootn_to_litter_fire_patch             (:) = nan
    allocate(this%m_frootn_storage_to_litter_fire_patch     (begp:endp)) ; this%m_frootn_storage_to_litter_fire_patch     (:) = nan
    allocate(this%m_frootn_xfer_to_litter_fire_patch        (begp:endp)) ; this%m_frootn_xfer_to_litter_fire_patch        (:) = nan
    allocate(this%m_livecrootn_to_litter_fire_patch         (begp:endp)) ; this%m_livecrootn_to_litter_fire_patch         (:) = nan
    allocate(this%m_livecrootn_storage_to_litter_fire_patch (begp:endp)) ; this%m_livecrootn_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_livecrootn_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_livecrootn_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_livecrootn_to_deadcrootn_fire_patch     (begp:endp)) ; this%m_livecrootn_to_deadcrootn_fire_patch     (:) = nan
    allocate(this%m_deadcrootn_to_litter_fire_patch         (begp:endp)) ; this%m_deadcrootn_to_litter_fire_patch         (:) = nan
    allocate(this%m_deadcrootn_storage_to_litter_fire_patch (begp:endp)) ; this%m_deadcrootn_storage_to_litter_fire_patch (:) = nan
    allocate(this%m_deadcrootn_xfer_to_litter_fire_patch    (begp:endp)) ; this%m_deadcrootn_xfer_to_litter_fire_patch    (:) = nan
    allocate(this%m_retransn_to_litter_fire_patch           (begp:endp)) ; this%m_retransn_to_litter_fire_patch           (:) = nan


    allocate(this%leafn_xfer_to_leafn_patch                 (begp:endp)) ; this%leafn_xfer_to_leafn_patch                 (:) = nan
    allocate(this%frootn_xfer_to_frootn_patch               (begp:endp)) ; this%frootn_xfer_to_frootn_patch               (:) = nan
    allocate(this%livestemn_xfer_to_livestemn_patch         (begp:endp)) ; this%livestemn_xfer_to_livestemn_patch         (:) = nan
    allocate(this%deadstemn_xfer_to_deadstemn_patch         (begp:endp)) ; this%deadstemn_xfer_to_deadstemn_patch         (:) = nan
    allocate(this%livecrootn_xfer_to_livecrootn_patch       (begp:endp)) ; this%livecrootn_xfer_to_livecrootn_patch       (:) = nan
    allocate(this%deadcrootn_xfer_to_deadcrootn_patch       (begp:endp)) ; this%deadcrootn_xfer_to_deadcrootn_patch       (:) = nan
    allocate(this%leafn_to_litter_patch                     (begp:endp)) ; this%leafn_to_litter_patch                     (:) = nan
    allocate(this%leafn_to_retransn_patch                   (begp:endp)) ; this%leafn_to_retransn_patch                   (:) = nan
    allocate(this%frootn_to_retransn_patch                  (begp:endp)) ; this%frootn_to_retransn_patch                  (:) = nan
    allocate(this%frootn_to_litter_patch                    (begp:endp)) ; this%frootn_to_litter_patch                    (:) = nan
    allocate(this%retransn_to_npool_patch                   (begp:endp)) ; this%retransn_to_npool_patch                   (:) = nan
    allocate(this%free_retransn_to_npool_patch              (begp:endp)) ; this%free_retransn_to_npool_patch              (:) = nan
    allocate(this%sminn_to_npool_patch                      (begp:endp)) ; this%sminn_to_npool_patch                      (:) = nan

    allocate(this%npool_to_leafn_patch                      (begp:endp)) ; this%npool_to_leafn_patch                      (:) = nan
    allocate(this%npool_to_leafn_storage_patch              (begp:endp)) ; this%npool_to_leafn_storage_patch              (:) = nan
    allocate(this%npool_to_frootn_patch                     (begp:endp)) ; this%npool_to_frootn_patch                     (:) = nan
    allocate(this%npool_to_frootn_storage_patch             (begp:endp)) ; this%npool_to_frootn_storage_patch             (:) = nan
    allocate(this%npool_to_livestemn_patch                  (begp:endp)) ; this%npool_to_livestemn_patch                  (:) = nan
    allocate(this%npool_to_livestemn_storage_patch          (begp:endp)) ; this%npool_to_livestemn_storage_patch          (:) = nan
    allocate(this%npool_to_deadstemn_patch                  (begp:endp)) ; this%npool_to_deadstemn_patch                  (:) = nan
    allocate(this%npool_to_deadstemn_storage_patch          (begp:endp)) ; this%npool_to_deadstemn_storage_patch          (:) = nan
    allocate(this%npool_to_livecrootn_patch                 (begp:endp)) ; this%npool_to_livecrootn_patch                 (:) = nan
    allocate(this%npool_to_livecrootn_storage_patch         (begp:endp)) ; this%npool_to_livecrootn_storage_patch         (:) = nan
    allocate(this%npool_to_deadcrootn_patch                 (begp:endp)) ; this%npool_to_deadcrootn_patch                 (:) = nan
    allocate(this%npool_to_deadcrootn_storage_patch         (begp:endp)) ; this%npool_to_deadcrootn_storage_patch         (:) = nan
    allocate(this%leafn_storage_to_xfer_patch               (begp:endp)) ; this%leafn_storage_to_xfer_patch               (:) = nan
    allocate(this%frootn_storage_to_xfer_patch              (begp:endp)) ; this%frootn_storage_to_xfer_patch              (:) = nan
    allocate(this%livestemn_storage_to_xfer_patch           (begp:endp)) ; this%livestemn_storage_to_xfer_patch           (:) = nan
    allocate(this%deadstemn_storage_to_xfer_patch           (begp:endp)) ; this%deadstemn_storage_to_xfer_patch           (:) = nan
    allocate(this%livecrootn_storage_to_xfer_patch          (begp:endp)) ; this%livecrootn_storage_to_xfer_patch          (:) = nan
    allocate(this%deadcrootn_storage_to_xfer_patch          (begp:endp)) ; this%deadcrootn_storage_to_xfer_patch          (:) = nan
    allocate(this%livestemn_to_deadstemn_patch              (begp:endp)) ; this%livestemn_to_deadstemn_patch              (:) = nan
    allocate(this%livestemn_to_retransn_patch               (begp:endp)) ; this%livestemn_to_retransn_patch               (:) = nan
    allocate(this%livecrootn_to_deadcrootn_patch            (begp:endp)) ; this%livecrootn_to_deadcrootn_patch            (:) = nan
    allocate(this%livecrootn_to_retransn_patch              (begp:endp)) ; this%livecrootn_to_retransn_patch              (:) = nan
    allocate(this%ndeploy_patch                             (begp:endp)) ; this%ndeploy_patch                             (:) = nan
    allocate(this%wood_harvestn_patch                       (begp:endp)) ; this%wood_harvestn_patch                       (:) = nan
    allocate(this%fire_nloss_patch                          (begp:endp)) ; this%fire_nloss_patch                          (:) = spval
    allocate(this%npool_to_grainn_patch                     (begp:endp)) ; this%npool_to_grainn_patch                     (:) = nan
    allocate(this%npool_to_grainn_storage_patch             (begp:endp)) ; this%npool_to_grainn_storage_patch             (:) = nan
    allocate(this%livestemn_to_litter_patch                 (begp:endp)) ; this%livestemn_to_litter_patch                 (:) = nan
    allocate(this%grainn_to_food_patch                      (begp:endp)) ; this%grainn_to_food_patch                      (:) = nan
    allocate(this%leafn_to_biofueln_patch                   (begp:endp)) ; this%leafn_to_biofueln_patch                   (:) = nan
    allocate(this%livestemn_to_biofueln_patch               (begp:endp)) ; this%livestemn_to_biofueln_patch               (:) = nan
    allocate(this%grainn_to_seed_patch                      (begp:endp)) ; this%grainn_to_seed_patch                      (:) = nan
    allocate(this%grainn_xfer_to_grainn_patch               (begp:endp)) ; this%grainn_xfer_to_grainn_patch               (:) = nan
    allocate(this%grainn_storage_to_xfer_patch              (begp:endp)) ; this%grainn_storage_to_xfer_patch              (:) = nan
    allocate(this%fert_patch                                (begp:endp)) ; this%fert_patch                                (:) = nan
    allocate(this%fert_counter_patch                        (begp:endp)) ; this%fert_counter_patch                        (:) = nan
    allocate(this%soyfixn_patch                             (begp:endp)) ; this%soyfixn_patch                             (:) = nan

    allocate(this%grainn_to_cropprodn_patch                 (begp:endp)) ; this%grainn_to_cropprodn_patch                 (:) = spval
    allocate(this%grainn_to_cropprodn_col                   (begc:endc)) ; this%grainn_to_cropprodn_col                   (:) = nan

    allocate(this%fire_nloss_col                            (begc:endc)) ; this%fire_nloss_col                            (:) = nan
    allocate(this%fire_nloss_p2c_col                        (begc:endc)) ; this%fire_nloss_p2c_col                        (:) = nan

    allocate(this%m_n_to_litr_met_fire_col     (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_met_fire_col     (:,:) = nan
    allocate(this%m_n_to_litr_cel_fire_col     (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_cel_fire_col     (:,:) = nan
    allocate(this%m_n_to_litr_lig_fire_col     (begc:endc,1:nlevdecomp_full)) ; this%m_n_to_litr_lig_fire_col     (:,:) = nan

    allocate(this%dwt_seedn_to_leaf_patch      (begp:endp))                   ; this%dwt_seedn_to_leaf_patch      (:)   = nan
    allocate(this%dwt_seedn_to_leaf_grc        (begg:endg))                   ; this%dwt_seedn_to_leaf_grc        (:)   = nan
    allocate(this%dwt_seedn_to_deadstem_patch  (begp:endp))                   ; this%dwt_seedn_to_deadstem_patch  (:)   = nan
    allocate(this%dwt_seedn_to_deadstem_grc    (begg:endg))                   ; this%dwt_seedn_to_deadstem_grc    (:)   = nan
    allocate(this%dwt_conv_nflux_patch         (begp:endp))                   ; this%dwt_conv_nflux_patch         (:)   = nan
    allocate(this%dwt_conv_nflux_grc           (begg:endg))                   ; this%dwt_conv_nflux_grc           (:)   = nan
    allocate(this%dwt_wood_productn_gain_patch (begp:endp))                   ; this%dwt_wood_productn_gain_patch (:)   = nan
    allocate(this%dwt_crop_productn_gain_patch (begp:endp))                   ; this%dwt_crop_productn_gain_patch (:)   = nan
    allocate(this%wood_harvestn_col            (begc:endc))                   ; this%wood_harvestn_col            (:)   = nan

    allocate(this%dwt_frootn_to_litr_met_n_col (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_met_n_col (:,:) = nan
    allocate(this%dwt_frootn_to_litr_cel_n_col (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_cel_n_col (:,:) = nan
    allocate(this%dwt_frootn_to_litr_lig_n_col (begc:endc,1:nlevdecomp_full)) ; this%dwt_frootn_to_litr_lig_n_col (:,:) = nan
    allocate(this%dwt_livecrootn_to_cwdn_col   (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootn_to_cwdn_col   (:,:) = nan
    allocate(this%dwt_deadcrootn_to_cwdn_col   (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootn_to_cwdn_col   (:,:) = nan

    allocate(this%crop_seedn_to_leaf_patch     (begp:endp))                   ; this%crop_seedn_to_leaf_patch     (:)   = nan

    allocate(this%m_decomp_npools_to_fire_vr_col    (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(this%m_decomp_npools_to_fire_col       (begc:endc,1:ndecomp_pools                  ))

    this%m_decomp_npools_to_fire_vr_col   (:,:,:) = nan
    this%m_decomp_npools_to_fire_col      (:,:)   = nan

    allocate(this%phenology_n_to_litr_met_n_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_n_to_litr_cel_n_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_n_to_litr_lig_n_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_litr_met_n_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_litr_cel_n_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_litr_lig_n_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_n_to_cwdn_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%fire_mortality_n_to_cwdn_col      (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_met_n_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_cel_n_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_litr_lig_n_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_n_to_cwdn_col             (begc:endc, 1:nlevdecomp_full))

    this%phenology_n_to_litr_met_n_col     (:,:) = nan
    this%phenology_n_to_litr_cel_n_col     (:,:) = nan
    this%phenology_n_to_litr_lig_n_col     (:,:) = nan
    this%gap_mortality_n_to_litr_met_n_col (:,:) = nan
    this%gap_mortality_n_to_litr_cel_n_col (:,:) = nan
    this%gap_mortality_n_to_litr_lig_n_col (:,:) = nan
    this%gap_mortality_n_to_cwdn_col       (:,:) = nan
    this%fire_mortality_n_to_cwdn_col      (:,:) = nan
    this%harvest_n_to_litr_met_n_col       (:,:) = nan
    this%harvest_n_to_litr_cel_n_col       (:,:) = nan
    this%harvest_n_to_litr_lig_n_col       (:,:) = nan
    this%harvest_n_to_cwdn_col             (:,:) = nan

    allocate(this%plant_ndemand_patch         (begp:endp)) ;    this%plant_ndemand_patch         (:) = spval
    allocate(this%avail_retransn_patch        (begp:endp)) ;    this%avail_retransn_patch        (:) = nan
    allocate(this%plant_nalloc_patch          (begp:endp)) ;    this%plant_nalloc_patch          (:) = nan

    allocate(this%plant_ndemand_retrans_patch (begp:endp)) ;    this%plant_ndemand_retrans_patch (:) = nan
    allocate(this%plant_ndemand_season_patch  (begp:endp)) ;    this%plant_ndemand_season_patch  (:) = nan
    allocate(this%plant_ndemand_stress_patch  (begp:endp)) ;    this%plant_ndemand_stress_patch  (:) = nan
    allocate(this%Nactive_patch               (begp:endp)) ;    this%Nactive_patch               (:) = nan
    allocate(this%Nnonmyc_patch               (begp:endp)) ;    this%Nnonmyc_patch               (:) = nan
    allocate(this%Nam_patch                   (begp:endp)) ;    this%Nam_patch                   (:) = nan
    allocate(this%Necm_patch                  (begp:endp)) ;    this%Necm_patch                  (:) = nan
    allocate(this%Nactive_no3_patch           (begp:endp)) ;    this%Nactive_no3_patch           (:) = nan
    allocate(this%Nactive_nh4_patch           (begp:endp)) ;    this%Nactive_nh4_patch           (:) = nan
    allocate(this%Nnonmyc_no3_patch           (begp:endp)) ;    this%Nnonmyc_no3_patch           (:) = nan
    allocate(this%Nnonmyc_nh4_patch           (begp:endp)) ;    this%Nnonmyc_nh4_patch           (:) = nan
    allocate(this%Nam_no3_patch               (begp:endp)) ;    this%Nam_no3_patch               (:) = nan
    allocate(this%Nam_nh4_patch               (begp:endp)) ;    this%Nam_nh4_patch               (:) = nan
    allocate(this%Necm_no3_patch              (begp:endp)) ;    this%Necm_no3_patch              (:) = nan
    allocate(this%Necm_nh4_patch              (begp:endp)) ;    this%Necm_nh4_patch              (:) = nan
    allocate(this%Npassive_patch              (begp:endp)) ;    this%Npassive_patch              (:) = nan
    allocate(this%Nfix_patch                  (begp:endp)) ;    this%Nfix_patch                  (:) = nan
    allocate(this%Nretrans_patch              (begp:endp)) ;    this%Nretrans_patch              (:) = nan
    allocate(this%Nretrans_org_patch          (begp:endp)) ;    this%Nretrans_org_patch          (:) = nan
    allocate(this%Nretrans_season_patch       (begp:endp)) ;    this%Nretrans_season_patch       (:) = nan
    allocate(this%Nretrans_stress_patch       (begp:endp)) ;    this%Nretrans_stress_patch       (:) = nan
    allocate(this%Nuptake_patch               (begp:endp)) ;    this%Nuptake_patch               (:) = nan
    allocate(this%sminn_to_plant_fun_patch    (begp:endp)) ;    this%sminn_to_plant_fun_patch    (:) = nan
    allocate(this%sminn_to_plant_fun_vr_patch (begp:endp,1:nlevdecomp_full))
    this%sminn_to_plant_fun_vr_patch          (:,:) = nan
    allocate(this%sminn_to_plant_fun_no3_vr_patch (begp:endp,1:nlevdecomp_full))
    this%sminn_to_plant_fun_no3_vr_patch      (:,:) = nan
    allocate(this%sminn_to_plant_fun_nh4_vr_patch (begp:endp,1:nlevdecomp_full))
    this%sminn_to_plant_fun_nh4_vr_patch      (:,:) = nan
    allocate(this%cost_nfix_patch              (begp:endp)) ;    this%cost_nfix_patch            (:) = nan
    allocate(this%cost_nactive_patch           (begp:endp)) ;    this%cost_nactive_patch         (:) = nan
    allocate(this%cost_nretrans_patch          (begp:endp)) ;    this%cost_nretrans_patch        (:) = nan
    allocate(this%nuptake_npp_fraction_patch   (begp:endp)) ;    this%nuptake_npp_fraction_patch            (:) = nan

 ! initialize variables from restart file or set to cold start value
 n = 0
 np = 0
    do nc = 1,nch        ! catchment tile loop
       do nz = 1,num_zon    ! CN zone loop
          n = n + 1

          

          do p = 0,numpft  ! PFT index loop
             np = np + 1
             do nv = 1,num_veg ! defined veg loop
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then

                  this%plant_ndemand_patch (np) = cnpft(nc,nz,nv, 75)
                 
                 end if
            end do !nv
            this%dwt_wood_productn_gain_patch(np) = 0.   ! following CNCLM45 setting
            this%dwt_crop_productn_gain_patch(np) = 0.   ! following CNCLM45 setting
       end do ! p
     end do ! nz
  end do ! nc     

    do p = begp,endp
       l = patch%landunit(p)

       if ( use_crop )then
          this%fert_counter_patch(p)  = spval
          this%fert_patch(p)          = 0._r8
          this%soyfixn_patch(p)       = 0._r8
       end if

       if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
          this%fert_counter_patch(p)  = 0._r8
       end if
       if ( use_fun ) then      !previously set to spval for special land units
          if (lun%ifspecial(l)) then
             this%plant_ndemand_patch(p)        = 0._r8
             this%avail_retransn_patch(p)       = 0._r8
             this%plant_nalloc_patch(p)         = 0._r8
             this%Npassive_patch(p)             = 0._r8
             this%Nactive_patch(p)              = 0._r8
             this%Nnonmyc_patch(p)              = 0._r8
             this%Nam_patch(p)                  = 0._r8
             this%Necm_patch(p)                 = 0._r8
             if (use_nitrif_denitrif) then
                this%Nactive_no3_patch(p)       = 0._r8
                this%Nactive_nh4_patch(p)       = 0._r8
                this%Nnonmyc_no3_patch(p)       = 0._r8
                this%Nnonmyc_nh4_patch(p)       = 0._r8
                this%Nam_no3_patch(p)           = 0._r8
                this%Nam_nh4_patch(p)           = 0._r8
                this%Necm_no3_patch(p)          = 0._r8
                this%Necm_nh4_patch(p)          = 0._r8
             end if
             this%Nfix_patch(p)                 = 0._r8
             this%Nretrans_patch(p)             = 0._r8
             this%Nretrans_org_patch(p)         = 0._r8
             this%Nretrans_season_patch(p)      = 0._r8
             this%Nretrans_stress_patch(p)      = 0._r8
             this%Nuptake_patch(p)              = 0._r8
             this%sminn_to_plant_fun_patch(p)   = 0._r8
             this%cost_nfix_patch               = 0._r8
             this%cost_nactive_patch            = 0._r8
             this%cost_nretrans_patch           = 0._r8
             this%nuptake_npp_fraction_patch    = 0._r8

             do j = 1, nlevdecomp
                this%sminn_to_plant_fun_vr_patch(p,j)       = 0._r8
                this%sminn_to_plant_fun_no3_vr_patch(p,j)   = 0._r8
                this%sminn_to_plant_fun_nh4_vr_patch(p,j)   = 0._r8
             end do
          end if
       end if
    end do


  end subroutine Init

!------------------------------------------
  subroutine SetValues ( this,nvegnpool, &
       num_patch, filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    integer , intent(in) :: num_patch,nvegnpool
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------

    do fi = 1,num_patch
       i=filter_patch(fi)

       this%m_leafn_to_litter_patch(i)                   = value_patch
       this%m_frootn_to_litter_patch(i)                  = value_patch
       this%m_leafn_storage_to_litter_patch(i)           = value_patch
       this%m_frootn_storage_to_litter_patch(i)          = value_patch
       this%m_livestemn_storage_to_litter_patch(i)       = value_patch
       this%m_deadstemn_storage_to_litter_patch(i)       = value_patch
       this%m_livecrootn_storage_to_litter_patch(i)      = value_patch
       this%m_deadcrootn_storage_to_litter_patch(i)      = value_patch
       this%m_leafn_xfer_to_litter_patch(i)              = value_patch
       this%m_frootn_xfer_to_litter_patch(i)             = value_patch
       this%m_livestemn_xfer_to_litter_patch(i)          = value_patch
       this%m_deadstemn_xfer_to_litter_patch(i)          = value_patch
       this%m_livecrootn_xfer_to_litter_patch(i)         = value_patch
       this%m_deadcrootn_xfer_to_litter_patch(i)         = value_patch
       this%m_livestemn_to_litter_patch(i)               = value_patch
       this%m_deadstemn_to_litter_patch(i)               = value_patch
       this%m_livecrootn_to_litter_patch(i)              = value_patch
       this%m_deadcrootn_to_litter_patch(i)              = value_patch
       this%m_retransn_to_litter_patch(i)                = value_patch
       this%hrv_leafn_to_litter_patch(i)                 = value_patch
       this%hrv_frootn_to_litter_patch(i)                = value_patch
       this%hrv_leafn_storage_to_litter_patch(i)         = value_patch
       this%hrv_frootn_storage_to_litter_patch(i)        = value_patch
       this%hrv_livestemn_storage_to_litter_patch(i)     = value_patch
       this%hrv_deadstemn_storage_to_litter_patch(i)     = value_patch
       this%hrv_livecrootn_storage_to_litter_patch(i)    = value_patch
       this%hrv_deadcrootn_storage_to_litter_patch(i)    = value_patch
       this%hrv_leafn_xfer_to_litter_patch(i)            = value_patch
       this%hrv_frootn_xfer_to_litter_patch(i)           = value_patch
       this%hrv_livestemn_xfer_to_litter_patch(i)        = value_patch
       this%hrv_deadstemn_xfer_to_litter_patch(i)        = value_patch
       this%hrv_livecrootn_xfer_to_litter_patch(i)       = value_patch
       this%hrv_deadcrootn_xfer_to_litter_patch(i)       = value_patch
       this%hrv_livestemn_to_litter_patch(i)             = value_patch
       this%hrv_livecrootn_to_litter_patch(i)            = value_patch
       this%hrv_deadcrootn_to_litter_patch(i)            = value_patch
       this%hrv_retransn_to_litter_patch(i)              = value_patch

       this%m_leafn_to_fire_patch(i)                     = value_patch
       this%m_leafn_storage_to_fire_patch(i)             = value_patch
       this%m_leafn_xfer_to_fire_patch(i)                = value_patch
       this%m_livestemn_to_fire_patch(i)                 = value_patch
       this%m_livestemn_storage_to_fire_patch(i)         = value_patch
       this%m_livestemn_xfer_to_fire_patch(i)            = value_patch
       this%m_deadstemn_to_fire_patch(i)                 = value_patch
       this%m_deadstemn_storage_to_fire_patch(i)         = value_patch
       this%m_deadstemn_xfer_to_fire_patch(i)            = value_patch
       this%m_frootn_to_fire_patch(i)                    = value_patch
       this%m_frootn_storage_to_fire_patch(i)            = value_patch
       this%m_frootn_xfer_to_fire_patch(i)               = value_patch
       this%m_livecrootn_to_fire_patch(i)                = value_patch
       this%m_livecrootn_storage_to_fire_patch(i)        = value_patch
       this%m_livecrootn_xfer_to_fire_patch(i)           = value_patch
       this%m_deadcrootn_to_fire_patch(i)                = value_patch
       this%m_deadcrootn_storage_to_fire_patch(i)        = value_patch
       this%m_deadcrootn_xfer_to_fire_patch(i)           = value_patch
       this%m_retransn_to_fire_patch(i)                  = value_patch


       this%m_leafn_to_litter_fire_patch(i)              = value_patch
       this%m_leafn_storage_to_litter_fire_patch(i)      = value_patch
       this%m_leafn_xfer_to_litter_fire_patch(i)         = value_patch
       this%m_livestemn_to_litter_fire_patch(i)          = value_patch
       this%m_livestemn_storage_to_litter_fire_patch(i)  = value_patch
       this%m_livestemn_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_livestemn_to_deadstemn_fire_patch(i)       = value_patch
       this%m_deadstemn_to_litter_fire_patch(i)          = value_patch
       this%m_deadstemn_storage_to_litter_fire_patch(i)  = value_patch
       this%m_deadstemn_xfer_to_litter_fire_patch(i)     = value_patch
       this%m_frootn_to_litter_fire_patch(i)             = value_patch
       this%m_frootn_storage_to_litter_fire_patch(i)     = value_patch
       this%m_frootn_xfer_to_litter_fire_patch(i)        = value_patch
       this%m_livecrootn_to_litter_fire_patch(i)         = value_patch
       this%m_livecrootn_storage_to_litter_fire_patch(i) = value_patch
       this%m_livecrootn_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_livecrootn_to_deadcrootn_fire_patch(i)     = value_patch
       this%m_deadcrootn_to_litter_fire_patch(i)         = value_patch
       this%m_deadcrootn_storage_to_litter_fire_patch(i) = value_patch
       this%m_deadcrootn_xfer_to_litter_fire_patch(i)    = value_patch
       this%m_retransn_to_litter_fire_patch(i)           = value_patch


       this%leafn_xfer_to_leafn_patch(i)                 = value_patch
       this%frootn_xfer_to_frootn_patch(i)               = value_patch
       this%livestemn_xfer_to_livestemn_patch(i)         = value_patch
       this%deadstemn_xfer_to_deadstemn_patch(i)         = value_patch
       this%livecrootn_xfer_to_livecrootn_patch(i)       = value_patch
       this%deadcrootn_xfer_to_deadcrootn_patch(i)       = value_patch
       this%leafn_to_litter_patch(i)                     = value_patch
       this%leafn_to_retransn_patch(i)                   = value_patch
       this%frootn_to_litter_patch(i)                    = value_patch
       this%retransn_to_npool_patch(i)                   = value_patch
       this%free_retransn_to_npool_patch(i)              = value_patch
       this%sminn_to_npool_patch(i)                      = value_patch
       this%npool_to_leafn_patch(i)                      = value_patch
       this%npool_to_leafn_storage_patch(i)              = value_patch
       this%npool_to_frootn_patch(i)                     = value_patch
       this%npool_to_frootn_storage_patch(i)             = value_patch
       this%npool_to_livestemn_patch(i)                  = value_patch
       this%npool_to_livestemn_storage_patch(i)          = value_patch
       this%npool_to_deadstemn_patch(i)                  = value_patch
       this%npool_to_deadstemn_storage_patch(i)          = value_patch
       this%npool_to_livecrootn_patch(i)                 = value_patch
       this%npool_to_livecrootn_storage_patch(i)         = value_patch
       this%npool_to_deadcrootn_patch(i)                 = value_patch
       this%npool_to_deadcrootn_storage_patch(i)         = value_patch
       this%leafn_storage_to_xfer_patch(i)               = value_patch
       this%frootn_storage_to_xfer_patch(i)              = value_patch
       this%livestemn_storage_to_xfer_patch(i)           = value_patch
       this%deadstemn_storage_to_xfer_patch(i)           = value_patch
       this%livecrootn_storage_to_xfer_patch(i)          = value_patch
       this%deadcrootn_storage_to_xfer_patch(i)          = value_patch
       this%livestemn_to_deadstemn_patch(i)              = value_patch
       this%livestemn_to_retransn_patch(i)               = value_patch
       this%livecrootn_to_deadcrootn_patch(i)            = value_patch
       this%livecrootn_to_retransn_patch(i)              = value_patch
       this%ndeploy_patch(i)                             = value_patch
       this%wood_harvestn_patch(i)                       = value_patch
       this%fire_nloss_patch(i)                          = value_patch

       this%crop_seedn_to_leaf_patch(i)                  = value_patch
       this%grainn_to_cropprodn_patch(i)                 = value_patch
    end do

    if ( use_crop )then
       do fi = 1,num_patch
          i = filter_patch(fi)
          this%livestemn_to_litter_patch(i)              = value_patch
          this%grainn_to_food_patch(i)                   = value_patch
          this%leafn_to_biofueln_patch(i)                = value_patch
          this%livestemn_to_biofueln_patch(i)            = value_patch
          this%grainn_to_seed_patch(i)                   = value_patch
          this%grainn_xfer_to_grainn_patch(i)            = value_patch
          this%npool_to_grainn_patch(i)                  = value_patch
          this%npool_to_grainn_storage_patch(i)          = value_patch
          this%grainn_storage_to_xfer_patch(i)           = value_patch
          this%soyfixn_patch(i)                          = value_patch
          this%frootn_to_retransn_patch(i)               = value_patch
       end do
    end if

    do j = 1, nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)

          ! phenology: litterfall and crop fluxes associated wit
          this%phenology_n_to_litr_met_n_col(i,j)        = value_column
          this%phenology_n_to_litr_cel_n_col(i,j)        = value_column
          this%phenology_n_to_litr_lig_n_col(i,j)        = value_column

          ! gap mortality
          this%gap_mortality_n_to_litr_met_n_col(i,j)    = value_column
          this%gap_mortality_n_to_litr_cel_n_col(i,j)    = value_column
          this%gap_mortality_n_to_litr_lig_n_col(i,j)    = value_column
          this%gap_mortality_n_to_cwdn_col(i,j)          = value_column

          ! fire
          this%fire_mortality_n_to_cwdn_col(i,j)         = value_column
          this%m_n_to_litr_met_fire_col(i,j)             = value_column
          this%m_n_to_litr_cel_fire_col(i,j)             = value_column
          this%m_n_to_litr_lig_fire_col(i,j)             = value_column

          ! harvest
          this%harvest_n_to_litr_met_n_col(i,j)          = value_column
          this%harvest_n_to_litr_cel_n_col(i,j)          = value_column
          this%harvest_n_to_litr_lig_n_col(i,j)          = value_column
          this%harvest_n_to_cwdn_col(i,j)                = value_column
       end do
    end do

    do fi = 1,num_column
       i = filter_column(fi)

       this%grainn_to_cropprodn_col(i)       = value_column
       this%fire_nloss_col(i)                = value_column

       ! Zero p2c column fluxes
       this%fire_nloss_col(i) = value_column
       this%wood_harvestn_col(i) = value_column
    end do

    do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%m_decomp_npools_to_fire_col(i,k) = value_column
       end do
    end do
! Matrix
    if(use_matrixcn)then
       do j = 1, nvegnpool
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_nalloc_patch(i,j)       = value_patch
             this%matrix_nphturnover_patch (i,j) = value_patch
             this%matrix_ngmturnover_patch (i,j) = value_patch
             this%matrix_nfiturnover_patch (i,j) = value_patch
          end do
       end do

       do j = 1, nnphtrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_nphtransfer_patch (i,j) = value_patch
          end do
       end do

       do j = 1, nngmtrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_ngmtransfer_patch (i,j) = value_patch
          end do
       end do

       do j = 1, nnfitrans
          do fi = 1,num_patch
             i = filter_patch(fi)
             this%matrix_nfitransfer_patch (i,j) = value_patch
          end do
       end do

    end if
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%m_decomp_npools_to_fire_vr_col(i,j,k) = value_column
          end do
       end do
    end do

  end subroutine SetValues

 !-----------------------------------------------------------------------
  subroutine Summary_nitrogenflux(this, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use clm_varpar    , only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl    , only: use_nitrif_denitrif
    use subgridAveMod , only: p2c
    !
    ! !ARGUMENTS:
    class (cnveg_nitrogenflux_type) :: this
    type(bounds_type) , intent(in) :: bounds
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer           , intent(in) :: num_soilp       ! number of soil patches in filter
    integer           , intent(in) :: filter_soilp(:) ! filter for soil patches
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l   ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
       this%ndeploy_patch(p) = &
            this%sminn_to_npool_patch(p) + &
            this%retransn_to_npool_patch(p) + &
            this%free_retransn_to_npool_patch(p)


       ! total patch-level fire N losses
       this%fire_nloss_patch(p) = &
            this%m_leafn_to_fire_patch(p)               + &
            this%m_leafn_storage_to_fire_patch(p)       + &
            this%m_leafn_xfer_to_fire_patch(p)          + &
            this%m_frootn_to_fire_patch(p)              + &
            this%m_frootn_storage_to_fire_patch(p)      + &
            this%m_frootn_xfer_to_fire_patch(p)         + &
            this%m_livestemn_to_fire_patch(p)           + &
            this%m_livestemn_storage_to_fire_patch(p)   + &
            this%m_livestemn_xfer_to_fire_patch(p)      + &
            this%m_deadstemn_to_fire_patch(p)           + &
            this%m_deadstemn_storage_to_fire_patch(p)   + &
            this%m_deadstemn_xfer_to_fire_patch(p)      + &
            this%m_livecrootn_to_fire_patch(p)          + &
            this%m_livecrootn_storage_to_fire_patch(p)  + &
            this%m_livecrootn_xfer_to_fire_patch(p)     + &
            this%m_deadcrootn_to_fire_patch(p)          + &
            this%m_deadcrootn_storage_to_fire_patch(p)  + &
            this%m_deadcrootn_xfer_to_fire_patch(p)     + &
            this%m_retransn_to_fire_patch(p)

    end do

    call p2c(bounds, num_soilc, filter_soilc, &
         this%fire_nloss_patch(bounds%begp:bounds%endp), &
         this%fire_nloss_p2c_col(bounds%begc:bounds%endc))


    ! vertically integrate column-level fire N losses
    do k = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%m_decomp_npools_to_fire_col(c,k) = &
                  this%m_decomp_npools_to_fire_col(c,k) + &
                  this%m_decomp_npools_to_fire_vr_col(c,j,k) * dzsoi_decomp(j)
          end do
       end do
    end do

    ! total column-level fire N losses
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       this%fire_nloss_col(c) = this%fire_nloss_p2c_col(c)
    end do
    do k = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          this%fire_nloss_col(c) = &
               this%fire_nloss_col(c) + &
               this%m_decomp_npools_to_fire_col(c,k)
       end do
    end do

  end subroutine Summary_nitrogenflux

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize flux variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenflux_type) :: this
    type(bounds_type), intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c, g, j          ! indices
    !-----------------------------------------------------------------------

    do g = bounds%begg, bounds%endg
       this%dwt_seedn_to_leaf_grc(g)     = 0._r8
       this%dwt_seedn_to_deadstem_grc(g) = 0._r8
       this%dwt_conv_nflux_grc(g)        = 0._r8
    end do

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          this%dwt_frootn_to_litr_met_n_col(c,j) = 0._r8
          this%dwt_frootn_to_litr_cel_n_col(c,j) = 0._r8
          this%dwt_frootn_to_litr_lig_n_col(c,j) = 0._r8
          this%dwt_livecrootn_to_cwdn_col(c,j)   = 0._r8
          this%dwt_deadcrootn_to_cwdn_col(c,j)   = 0._r8
       end do
    end do

  end subroutine ZeroDwt

end module CNVegNitrogenFluxType
