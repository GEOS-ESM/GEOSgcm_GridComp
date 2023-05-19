module CNVegNitrogenStateType

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R8
  use MAPL_ExceptionHandling
  use clm_varctl       , only : use_matrixcn, use_crop
  use clm_varctl       , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varpar       , only : NUM_ZON, NUM_VEG, VAR_COL, VAR_PFT, &
                                numpft, CN_zone_weight
  use clm_varpar       , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar       , only : nlevdecomp_full, nlevdecomp
  use clm_varcon       , only : spval, ispval, dzsoi_decomp, zisoi
  use nanMod           , only : nan
  use decompMod        , only : bounds_type
  use pftconMod                          , only : npcropmin
  use PatchType        , only : patch

  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:

 !
  type, public :: cnveg_nitrogenstate_type

     real(r8), pointer :: grainn_patch                        (:) ! (gN/m2) grain N (crop)
     real(r8), pointer :: grainn_storage_patch                (:) ! (gN/m2) grain N storage (crop)
     real(r8), pointer :: grainn_xfer_patch                   (:) ! (gN/m2) grain N transfer (crop)
     real(r8), pointer :: matrix_cap_grainn_patch             (:) ! (gN/m2) Capacity of grain N
     real(r8), pointer :: matrix_cap_grainn_storage_patch     (:) ! (gN/m2) Capacity of grain N storage
     real(r8), pointer :: matrix_cap_grainn_xfer_patch        (:) ! (gN/m2) Capacity of grain N transfer
     real(r8), pointer :: leafn_patch                         (:) ! (gN/m2) leaf N
     real(r8), pointer :: leafn_storage_patch                 (:) ! (gN/m2) leaf N storage
     real(r8), pointer :: leafn_xfer_patch                    (:) ! (gN/m2) leaf N transfer
     real(r8), pointer :: matrix_cap_leafn_patch              (:) ! (gN/m2) Capacity of leaf N
     real(r8), pointer :: matrix_cap_leafn_storage_patch      (:) ! (gN/m2) Capacity of leaf N storage
     real(r8), pointer :: matrix_cap_leafn_xfer_patch         (:) ! (gN/m2) Capacity of leaf N transfer
     real(r8), pointer :: leafn_storage_xfer_acc_patch        (:) ! (gN/m2) Accmulated leaf N transfer
     real(r8), pointer :: storage_ndemand_patch               (:) ! (gN/m2) N demand during the offset period
     real(r8), pointer :: frootn_patch                        (:) ! (gN/m2) fine root N
     real(r8), pointer :: frootn_storage_patch                (:) ! (gN/m2) fine root N storage
     real(r8), pointer :: frootn_xfer_patch                   (:) ! (gN/m2) fine root N transfer
     real(r8), pointer :: matrix_cap_frootn_patch             (:) ! (gN/m2) Capacity of fine root N
     real(r8), pointer :: matrix_cap_frootn_storage_patch     (:) ! (gN/m2) Capacity of fine root N storage
     real(r8), pointer :: matrix_cap_frootn_xfer_patch        (:) ! (gN/m2) Capacity of fine root N transfer
     real(r8), pointer :: livestemn_patch                     (:) ! (gN/m2) live stem N
     real(r8), pointer :: livestemn_storage_patch             (:) ! (gN/m2) live stem N storage
     real(r8), pointer :: livestemn_xfer_patch                (:) ! (gN/m2) live stem N transfer
     real(r8), pointer :: deadstemn_patch                     (:) ! (gN/m2) dead stem N
     real(r8), pointer :: deadstemn_storage_patch             (:) ! (gN/m2) dead stem N storage
     real(r8), pointer :: deadstemn_xfer_patch                (:) ! (gN/m2) dead stem N transfer
     real(r8), pointer :: livecrootn_patch                    (:) ! (gN/m2) live coarse root N
     real(r8), pointer :: livecrootn_storage_patch            (:) ! (gN/m2) live coarse root N storage
     real(r8), pointer :: livecrootn_xfer_patch               (:) ! (gN/m2) live coarse root N transfer
     real(r8), pointer :: deadcrootn_patch                    (:) ! (gN/m2) dead coarse root N
     real(r8), pointer :: deadcrootn_storage_patch            (:) ! (gN/m2) dead coarse root N storage
     real(r8), pointer :: deadcrootn_xfer_patch               (:) ! (gN/m2) dead coarse root N transfer
     real(r8), pointer :: matrix_cap_livestemn_patch          (:) ! (gN/m2) Capacity of live stem N
     real(r8), pointer :: matrix_cap_livestemn_storage_patch  (:) ! (gN/m2) Capacity of live stem N storage
     real(r8), pointer :: matrix_cap_livestemn_xfer_patch     (:) ! (gN/m2) Capacity of live stem N transfer
     real(r8), pointer :: matrix_cap_deadstemn_patch          (:) ! (gN/m2) Capacity of dead stem N
     real(r8), pointer :: matrix_cap_deadstemn_storage_patch  (:) ! (gN/m2) Capacity of dead stem N storage
     real(r8), pointer :: matrix_cap_deadstemn_xfer_patch     (:) ! (gN/m2) Capacity of dead stem N transfer
     real(r8), pointer :: matrix_cap_livecrootn_patch         (:) ! (gN/m2) Capacity of live coarse root N
     real(r8), pointer :: matrix_cap_livecrootn_storage_patch (:) ! (gN/m2) Capacity of live coarse root N storage
     real(r8), pointer :: matrix_cap_livecrootn_xfer_patch    (:) ! (gN/m2) Capacity of live coarse root N transfer
     real(r8), pointer :: matrix_cap_deadcrootn_patch         (:) ! (gN/m2) Capacity of dead coarse root N
     real(r8), pointer :: matrix_cap_deadcrootn_storage_patch (:) ! (gN/m2) Capacity of dead coarse root N storage
     real(r8), pointer :: matrix_cap_deadcrootn_xfer_patch    (:) ! (gN/m2) Capacity of dead coarse root N transfer
     real(r8), pointer :: retransn_patch                      (:) ! (gN/m2) plant pool of retranslocated N
     real(r8), pointer :: npool_patch                         (:) ! (gN/m2) temporary plant N pool
     real(r8), pointer :: ntrunc_patch                        (:) ! (gN/m2) patch-level sink for N truncation
     real(r8), pointer :: cropseedn_deficit_patch             (:) ! (gN/m2) pool for seeding new crop growth; this is a NEGATIVE term, indicating the amount of seed usage that needs to be repaid
     real(r8), pointer :: seedn_grc                           (:) ! (gN/m2) gridcell-level pool for seeding new pFTs via dynamic landcover
! Pool for initial step of year for matrix
     real(r8), pointer :: leafn0_patch                        (:) ! (gN/m2) Initial value of leaf N for SASU
     real(r8), pointer :: leafn0_storage_patch                (:) ! (gN/m2) Initial value of leaf N storage for SASU
     real(r8), pointer :: leafn0_xfer_patch                   (:) ! (gN/m2) Initial value of leaf N transfer for SASU
     real(r8), pointer :: frootn0_patch                       (:) ! (gN/m2) Initial value of fine root N for SASU
     real(r8), pointer :: frootn0_storage_patch               (:) ! (gN/m2) Initial value of fine root N storage for SASU
     real(r8), pointer :: frootn0_xfer_patch                  (:) ! (gN/m2) Initial value of fine root N transfer for SASU
     real(r8), pointer :: livestemn0_patch                    (:) ! (gN/m2) Initial value of live stem N for SASU
     real(r8), pointer :: livestemn0_storage_patch            (:) ! (gN/m2) Initial value of live stem N storage for SASU
     real(r8), pointer :: livestemn0_xfer_patch               (:) ! (gN/m2) Initial value of live stem N transfer for SASU
     real(r8), pointer :: deadstemn0_patch                    (:) ! (gN/m2) Initial value of dead stem N for SASU
     real(r8), pointer :: deadstemn0_storage_patch            (:) ! (gN/m2) Initial value of dead stem N storage for SASU
     real(r8), pointer :: deadstemn0_xfer_patch               (:) ! (gN/m2) Initial value of dead stem N transfer for SASU
     real(r8), pointer :: livecrootn0_patch                   (:) ! (gN/m2) Initial value of live coarse root N for SASU
     real(r8), pointer :: livecrootn0_storage_patch           (:) ! (gN/m2) Initial value of live coarse root N storage for SASU
     real(r8), pointer :: livecrootn0_xfer_patch              (:) ! (gN/m2) Initial value of live coarse root N transfer for SASU
     real(r8), pointer :: deadcrootn0_patch                   (:) ! (gN/m2) Initial value of dead coarse root N for SASU
     real(r8), pointer :: deadcrootn0_storage_patch           (:) ! (gN/m2) Initial value of dead coarse root N storage for SASU
     real(r8), pointer :: deadcrootn0_xfer_patch              (:) ! (gN/m2) Initial value of dead coarse root N transfer for SASU
     real(r8), pointer :: retransn0_patch                     (:) ! (gN/m2) Initial value of dead coarse root N transfer for SASU
     real(r8), pointer :: grainn0_patch                       (:) ! (gN/m2) Initial value of grain N for SASU
     real(r8), pointer :: grainn0_storage_patch               (:) ! (gN/m2) Initial value of grain N storage for SASU
     real(r8), pointer :: grainn0_xfer_patch                  (:) ! (gN/m2) Initial value of grain N transfer for SASU

     ! summary (diagnostic) state variables, not involved in mass balance
     real(r8), pointer :: dispvegn_patch                      (:) ! (gN/m2) displayed veg nitrogen, excluding storage
     real(r8), pointer :: storvegn_patch                      (:) ! (gN/m2) stored vegetation nitrogen
     real(r8), pointer :: totvegn_patch                       (:) ! (gN/m2) total vegetation nitrogen
     real(r8), pointer :: totvegn_col                         (:) ! (gN/m2) total vegetation nitrogen (p2c)
     real(r8), pointer :: totn_patch                          (:) ! (gN/m2) total patch-level nitrogen
     real(r8), pointer :: totn_p2c_col                        (:) ! (gN/m2) totn_patch averaged to col
     real(r8), pointer :: totn_col                            (:) ! (gN/m2) total column nitrogen, incl veg
     real(r8), pointer :: totecosysn_col                      (:) ! (gN/m2) total ecosystem nitrogen, incl veg
     real(r8), pointer :: totn_grc                            (:) ! (gN/m2) total gridcell nitrogen
! acc spinup
     real(r8), pointer :: matrix_nalloc_leaf_acc_patch        (:) ! (gN/m2/year) Input N allocated to leaf during this year
     real(r8), pointer :: matrix_nalloc_leafst_acc_patch      (:) ! (gN/m2/year) Input N allocated to leaf storage during this year
     real(r8), pointer :: matrix_nalloc_froot_acc_patch       (:) ! (gN/m2/year) Input N allocated to fine root during this year
     real(r8), pointer :: matrix_nalloc_frootst_acc_patch     (:) ! (gN/m2/year) Input N allocated to fine root storage during this year
     real(r8), pointer :: matrix_nalloc_livestem_acc_patch    (:) ! (gN/m2/year) Input N allocated to live stem during this year
     real(r8), pointer :: matrix_nalloc_livestemst_acc_patch  (:) ! (gN/m2/year) Input N allocated to live stem storage during this year
     real(r8), pointer :: matrix_nalloc_deadstem_acc_patch    (:) ! (gN/m2/year) Input N allocated to dead stem during this year
     real(r8), pointer :: matrix_nalloc_deadstemst_acc_patch  (:) ! (gN/m2/year) Input N allocated to dead stem storage during this year
     real(r8), pointer :: matrix_nalloc_livecroot_acc_patch   (:) ! (gN/m2/year) Input N allocated to live coarse root during this year
     real(r8), pointer :: matrix_nalloc_livecrootst_acc_patch (:) ! (gN/m2/year) Input N allocated to live coarse root storage during this year
     real(r8), pointer :: matrix_nalloc_deadcroot_acc_patch   (:) ! (gN/m2/year) Input N allocated to dead coarse root during this year
     real(r8), pointer :: matrix_nalloc_deadcrootst_acc_patch (:) ! (gN/m2/year) Input N allocated to dead coarse root storage during this year
     real(r8), pointer :: matrix_nalloc_grain_acc_patch       (:) ! (gN/m2/year) Input N allocated to grain during this year
     real(r8), pointer :: matrix_nalloc_grainst_acc_patch     (:) ! (gN/m2/year) Input N allocated to grain storage during this year

     real(r8), pointer :: matrix_ntransfer_leafst_to_leafxf_acc_patch           (:) ! (gN/m2/year) N transfer from leaf storage to leaf transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_leafxf_to_leaf_acc_patch             (:) ! (gN/m2/year) N transfer from leaf transfer to leaf pool during this year
     real(r8), pointer :: matrix_ntransfer_frootst_to_frootxf_acc_patch         (:) ! (gN/m2/year) N transfer from fine root storage to fine root transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_frootxf_to_froot_acc_patch           (:) ! (gN/m2/year) N transfer from fine root transfer to fine root pool during this year
     real(r8), pointer :: matrix_ntransfer_livestemst_to_livestemxf_acc_patch   (:) ! (gN/m2/year) N transfer from live stem storage to live stem transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_livestemxf_to_livestem_acc_patch     (:) ! (gN/m2/year) N transfer from live stem transfer to live stem pool during this year
     real(r8), pointer :: matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch   (:) ! (gN/m2/year) N transfer from dead stem storage to dead stem transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_deadstemxf_to_deadstem_acc_patch     (:) ! (gN/m2/year) N transfer from dead stem transfer to dead stem pool during this year
     real(r8), pointer :: matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch (:) ! (gN/m2/year) N transfer from live coarse root storage to live coarse root transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_livecrootxf_to_livecroot_acc_patch   (:) ! (gN/m2/year) N transfer from live coarse root transfer to live coarse root pool during this year
     real(r8), pointer :: matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch (:) ! (gN/m2/year) N transfer from dead coarse root storage to dead coarse root transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch   (:) ! (gN/m2/year) N transfer from dead coarse root transfer to dead coarse root pool during this year
     real(r8), pointer :: matrix_ntransfer_grainst_to_grainxf_acc_patch         (:) ! (gN/m2/year) N transfer from grain storage to grain transfer pool during this year
     real(r8), pointer :: matrix_ntransfer_grainxf_to_grain_acc_patch           (:) ! (gN/m2/year) N transfer from grain transfer to grain pool during this year
     real(r8), pointer :: matrix_ntransfer_livestem_to_deadstem_acc_patch       (:) ! (gN/m2/year) N transfer from live stem to dead stem pool during this year
     real(r8), pointer :: matrix_ntransfer_livecroot_to_deadcroot_acc_patch     (:) ! (gN/m2/year) N transfer from live coarse root to dead coarse root pool during this year

     real(r8), pointer :: matrix_ntransfer_retransn_to_leaf_acc_patch           (:) ! (gN/m2/year) N transfer from retranslocation to leaf pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_leafst_acc_patch         (:) ! (gN/m2/year) N transfer from retranslocation to leaf storage pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_froot_acc_patch          (:) ! (gN/m2/year) N transfer from retranslocation to fine root  pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_frootst_acc_patch        (:) ! (gN/m2/year) N transfer from retranslocation to fine root storage pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_livestem_acc_patch       (:) ! (gN/m2/year) N transfer from retranslocation to live stem pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_livestemst_acc_patch     (:) ! (gN/m2/year) N transfer from retranslocation to live stem storage pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_deadstem_acc_patch       (:) ! (gN/m2/year) N transfer from retranslocation to dead stem pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_deadstemst_acc_patch     (:) ! (gN/m2/year) N transfer from retranslocation to dead stem storage pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_livecroot_acc_patch      (:) ! (gN/m2/year) N transfer from retranslocation to live coarse root pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_livecrootst_acc_patch    (:) ! (gN/m2/year) N transfer from retranslocation to live coarse root storage pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_deadcroot_acc_patch      (:) ! (gN/m2/year) N transfer from retranslocation to dead coarse root pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_deadcrootst_acc_patch    (:) ! (gN/m2/year) N transfer from retranslocation to dead coarse root storage pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_grain_acc_patch          (:) ! (gN/m2/year) N transfer from retranslocation to grain pool during this year
     real(r8), pointer :: matrix_ntransfer_retransn_to_grainst_acc_patch        (:) ! (gN/m2/year) N transfer from retranslocation to grain storage pool during this year

     real(r8), pointer :: matrix_ntransfer_leaf_to_retransn_acc_patch           (:) ! (gN/m2/year) N transfer from leaf to retranslocation pool during this year
     real(r8), pointer :: matrix_ntransfer_froot_to_retransn_acc_patch          (:) ! (gN/m2/year) N transfer from fine root to retranslocation pool during this year
     real(r8), pointer :: matrix_ntransfer_livestem_to_retransn_acc_patch       (:) ! (gN/m2/year) N transfer from live stem to retranslocation pool during this year
     real(r8), pointer :: matrix_ntransfer_livecroot_to_retransn_acc_patch      (:) ! (gN/m2/year) N transfer from live coarse root to retranslocation pool during this year

     real(r8), pointer :: matrix_nturnover_leaf_acc_patch                       (:) ! (gN/m2/year) N turnover from leaf
     real(r8), pointer :: matrix_nturnover_leafst_acc_patch                     (:) ! (gN/m2/year) N turnover from leaf storage
     real(r8), pointer :: matrix_nturnover_leafxf_acc_patch                     (:) ! (gN/m2/year) N turnover from leaf transfer
     real(r8), pointer :: matrix_nturnover_froot_acc_patch                      (:) ! (gN/m2/year) N turnover from root
     real(r8), pointer :: matrix_nturnover_frootst_acc_patch                    (:) ! (gN/m2/year) N turnover from root storage
     real(r8), pointer :: matrix_nturnover_frootxf_acc_patch                    (:) ! (gN/m2/year) N turnover from root transfer
     real(r8), pointer :: matrix_nturnover_livestem_acc_patch                   (:) ! (gN/m2/year) N turnover from live stem
     real(r8), pointer :: matrix_nturnover_livestemst_acc_patch                 (:) ! (gN/m2/year) N turnover from live stem storage
     real(r8), pointer :: matrix_nturnover_livestemxf_acc_patch                 (:) ! (gN/m2/year) N turnover from live stem transfer
     real(r8), pointer :: matrix_nturnover_deadstem_acc_patch                   (:) ! (gN/m2/year) N turnover from dead stem
     real(r8), pointer :: matrix_nturnover_deadstemst_acc_patch                 (:) ! (gN/m2/year) N turnover from dead stem storage
     real(r8), pointer :: matrix_nturnover_deadstemxf_acc_patch                 (:) ! (gN/m2/year) N turnover from dead stem transfer
     real(r8), pointer :: matrix_nturnover_livecroot_acc_patch                  (:) ! (gN/m2/year) N turnover from live coarse root
     real(r8), pointer :: matrix_nturnover_livecrootst_acc_patch                (:) ! (gN/m2/year) N turnover from live coarse root storage
     real(r8), pointer :: matrix_nturnover_livecrootxf_acc_patch                (:) ! (gN/m2/year) N turnover from live coarse root transfer
     real(r8), pointer :: matrix_nturnover_deadcroot_acc_patch                  (:) ! (gN/m2/year) N turnover from dead coarse root
     real(r8), pointer :: matrix_nturnover_deadcrootst_acc_patch                (:) ! (gN/m2/year) N turnover from dead coarse root storage
     real(r8), pointer :: matrix_nturnover_deadcrootxf_acc_patch                (:) ! (gN/m2/year) N turnover from dead coarse root transfer
     real(r8), pointer :: matrix_nturnover_grain_acc_patch                      (:) ! (gN/m2/year) N turnover from grain
     real(r8), pointer :: matrix_nturnover_grainst_acc_patch                    (:) ! (gN/m2/year) N turnover from grain storage
     real(r8), pointer :: matrix_nturnover_grainxf_acc_patch                    (:) ! (gN/m2/year) N turnover from grain transfer
     real(r8), pointer :: matrix_nturnover_retransn_acc_patch                   (:) ! (gN/m2/year) N turnover from retranslocation transfer

     real(r8), pointer :: grainn_SASUsave_patch               (:) ! (gC/m2) grain C (crop model)
     real(r8), pointer :: grainn_storage_SASUsave_patch       (:) ! (gC/m2) grain C storage (crop model)
     real(r8), pointer :: leafn_SASUsave_patch                (:) ! (gC/m2) leaf C
     real(r8), pointer :: leafn_storage_SASUsave_patch        (:) ! (gC/m2) leaf C storage
     real(r8), pointer :: leafn_xfer_SASUsave_patch           (:) ! (gC/m2) leaf C transfer
     real(r8), pointer :: frootn_SASUsave_patch               (:) ! (gC/m2) fine root C
     real(r8), pointer :: frootn_storage_SASUsave_patch       (:) ! (gC/m2) fine root C storage
     real(r8), pointer :: frootn_xfer_SASUsave_patch          (:) ! (gC/m2) fine root C transfer
     real(r8), pointer :: livestemn_SASUsave_patch            (:) ! (gC/m2) live stem C
     real(r8), pointer :: livestemn_storage_SASUsave_patch    (:) ! (gC/m2) live stem C storage
     real(r8), pointer :: livestemn_xfer_SASUsave_patch       (:) ! (gC/m2) live stem C transfer
     real(r8), pointer :: deadstemn_SASUsave_patch            (:) ! (gC/m2) dead stem C
     real(r8), pointer :: deadstemn_storage_SASUsave_patch    (:) ! (gC/m2) dead stem C storage
     real(r8), pointer :: deadstemn_xfer_SASUsave_patch       (:) ! (gC/m2) dead stem C transfer
     real(r8), pointer :: livecrootn_SASUsave_patch           (:) ! (gC/m2) live coarse root C
     real(r8), pointer :: livecrootn_storage_SASUsave_patch   (:) ! (gC/m2) live coarse root C storage
     real(r8), pointer :: livecrootn_xfer_SASUsave_patch      (:) ! (gC/m2) live coarse root C transfer
     real(r8), pointer :: deadcrootn_SASUsave_patch           (:) ! (gC/m2) dead coarse root C
     real(r8), pointer :: deadcrootn_storage_SASUsave_patch   (:) ! (gC/m2) dead coarse root C storage
     real(r8), pointer :: deadcrootn_xfer_SASUsave_patch      (:) ! (gC/m2) dead coarse root C transfer:wq

  contains 

     procedure , public  :: Summary => Summary_nitrogenstate
     procedure , public  :: ZeroDWT
     procedure , public  :: Init

end type cnveg_nitrogenstate_type
type(cnveg_nitrogenstate_type), public, target, save :: cnveg_nitrogenstate_inst

contains

!-------------------------------------------------------------
  subroutine Init(this, bounds, nch, ityp, fveg, cncol, cnpft)

! !DESCRIPTION:
! Initialize CTSM nitrogen states
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
    class(cnveg_nitrogenstate_type)                           :: this


  ! LOCAL:

    integer  :: begp, endp, begg, endg, begc, endc
    integer  :: np, nc, nz, p, nv, n
  !---------------------------------------------------------------------

    begp = bounds%begp  ; endp = bounds%endp
    begg = bounds%begg  ; endg = bounds%endg
    begc = bounds%begc  ; endc = bounds%endc

    allocate(this%grainn_patch                           (begp:endp)) ; this%grainn_patch                        (:) = nan
    allocate(this%grainn_storage_patch                   (begp:endp)) ; this%grainn_storage_patch                (:) = nan
    allocate(this%grainn_xfer_patch                      (begp:endp)) ; this%grainn_xfer_patch                   (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_grainn_patch             (begp:endp)) ; this%matrix_cap_grainn_patch             (:) = nan
       allocate(this%matrix_cap_grainn_storage_patch     (begp:endp)) ; this%matrix_cap_grainn_storage_patch     (:) = nan
       allocate(this%matrix_cap_grainn_xfer_patch        (begp:endp)) ; this%matrix_cap_grainn_xfer_patch        (:) = nan
    end if
    allocate(this%leafn_patch                            (begp:endp)) ; this%leafn_patch                         (:) = nan
    allocate(this%leafn_storage_patch                    (begp:endp)) ; this%leafn_storage_patch                 (:) = nan
    allocate(this%leafn_xfer_patch                       (begp:endp)) ; this%leafn_xfer_patch                    (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_leafn_patch              (begp:endp)) ; this%matrix_cap_leafn_patch              (:) = nan
       allocate(this%matrix_cap_leafn_storage_patch      (begp:endp)) ; this%matrix_cap_leafn_storage_patch      (:) = nan
       allocate(this%matrix_cap_leafn_xfer_patch         (begp:endp)) ; this%matrix_cap_leafn_xfer_patch         (:) = nan
    end if
    allocate(this%leafn_storage_xfer_acc_patch           (begp:endp)) ; this%leafn_storage_xfer_acc_patch        (:) = nan
    allocate(this%storage_ndemand_patch                  (begp:endp)) ; this%storage_ndemand_patch               (:) = nan
    allocate(this%frootn_patch                           (begp:endp)) ; this%frootn_patch                        (:) = nan
    allocate(this%frootn_storage_patch                   (begp:endp)) ; this%frootn_storage_patch                (:) = nan
    allocate(this%frootn_xfer_patch                      (begp:endp)) ; this%frootn_xfer_patch                   (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_frootn_patch             (begp:endp)) ; this%matrix_cap_frootn_patch             (:) = nan
       allocate(this%matrix_cap_frootn_storage_patch     (begp:endp)) ; this%matrix_cap_frootn_storage_patch     (:) = nan
       allocate(this%matrix_cap_frootn_xfer_patch        (begp:endp)) ; this%matrix_cap_frootn_xfer_patch        (:) = nan
    end if
    allocate(this%livestemn_patch                        (begp:endp)) ; this%livestemn_patch                     (:) = nan
    allocate(this%livestemn_storage_patch                (begp:endp)) ; this%livestemn_storage_patch             (:) = nan
    allocate(this%livestemn_xfer_patch                   (begp:endp)) ; this%livestemn_xfer_patch                (:) = nan
    allocate(this%deadstemn_patch                        (begp:endp)) ; this%deadstemn_patch                     (:) = nan
    allocate(this%deadstemn_storage_patch                (begp:endp)) ; this%deadstemn_storage_patch             (:) = nan
    allocate(this%deadstemn_xfer_patch                   (begp:endp)) ; this%deadstemn_xfer_patch                (:) = nan
    allocate(this%livecrootn_patch                       (begp:endp)) ; this%livecrootn_patch                    (:) = nan
    allocate(this%livecrootn_storage_patch               (begp:endp)) ; this%livecrootn_storage_patch            (:) = nan
    allocate(this%livecrootn_xfer_patch                  (begp:endp)) ; this%livecrootn_xfer_patch               (:) = nan
    allocate(this%deadcrootn_patch                       (begp:endp)) ; this%deadcrootn_patch                    (:) = nan
    allocate(this%deadcrootn_storage_patch               (begp:endp)) ; this%deadcrootn_storage_patch            (:) = nan
    allocate(this%deadcrootn_xfer_patch                  (begp:endp)) ; this%deadcrootn_xfer_patch               (:) = nan
    if(use_matrixcn)then
       allocate(this%matrix_cap_livestemn_patch          (begp:endp)) ; this%matrix_cap_livestemn_patch          (:) = nan
       allocate(this%matrix_cap_livestemn_storage_patch  (begp:endp)) ; this%matrix_cap_livestemn_storage_patch  (:) = nan
       allocate(this%matrix_cap_livestemn_xfer_patch     (begp:endp)) ; this%matrix_cap_livestemn_xfer_patch     (:) = nan
       allocate(this%matrix_cap_deadstemn_patch          (begp:endp)) ; this%matrix_cap_deadstemn_patch          (:) = nan
       allocate(this%matrix_cap_deadstemn_storage_patch  (begp:endp)) ; this%matrix_cap_deadstemn_storage_patch  (:) = nan
       allocate(this%matrix_cap_deadstemn_xfer_patch     (begp:endp)) ; this%matrix_cap_deadstemn_xfer_patch     (:) = nan
       allocate(this%matrix_cap_livecrootn_patch         (begp:endp)) ; this%matrix_cap_livecrootn_patch         (:) = nan
       allocate(this%matrix_cap_livecrootn_storage_patch (begp:endp)) ; this%matrix_cap_livecrootn_storage_patch (:) = nan
       allocate(this%matrix_cap_livecrootn_xfer_patch    (begp:endp)) ; this%matrix_cap_livecrootn_xfer_patch    (:) = nan
       allocate(this%matrix_cap_deadcrootn_patch         (begp:endp)) ; this%matrix_cap_deadcrootn_patch         (:) = nan
       allocate(this%matrix_cap_deadcrootn_storage_patch (begp:endp)) ; this%matrix_cap_deadcrootn_storage_patch (:) = nan
       allocate(this%matrix_cap_deadcrootn_xfer_patch    (begp:endp)) ; this%matrix_cap_deadcrootn_xfer_patch    (:) = nan
    end if
    allocate(this%retransn_patch                         (begp:endp)) ; this%retransn_patch                      (:) = nan
    allocate(this%npool_patch                            (begp:endp)) ; this%npool_patch                         (:) = nan
    allocate(this%ntrunc_patch                           (begp:endp)) ; this%ntrunc_patch                        (:) = nan
    allocate(this%dispvegn_patch                         (begp:endp)) ; this%dispvegn_patch                      (:) = nan
    allocate(this%storvegn_patch                         (begp:endp)) ; this%storvegn_patch                      (:) = nan
    allocate(this%totvegn_patch                          (begp:endp)) ; this%totvegn_patch                       (:) = spval
    allocate(this%totn_patch                             (begp:endp)) ; this%totn_patch                          (:) = spval

    allocate(this%cropseedn_deficit_patch                (begp:endp)) ; this%cropseedn_deficit_patch             (:) = nan
    allocate(this%seedn_grc                              (begg:endg)) ; this%seedn_grc                           (:) = nan
    allocate(this%totvegn_col                            (begc:endc)) ; this%totvegn_col                         (:) = nan
    allocate(this%totn_p2c_col                           (begc:endc)) ; this%totn_p2c_col                        (:) = nan
    allocate(this%totn_col                               (begc:endc)) ; this%totn_col                            (:) = spval
    allocate(this%totecosysn_col                         (begc:endc)) ; this%totecosysn_col                      (:) = nan
    allocate(this%totn_grc                               (begg:endg)) ; this%totn_grc                            (:) = nan

    if(use_matrixcn)then
       allocate(this%leafn0_patch                        (begp:endp)) ; this%leafn0_patch                        (:) = nan
       allocate(this%leafn0_storage_patch                (begp:endp)) ; this%leafn0_storage_patch                (:) = nan
       allocate(this%leafn0_xfer_patch                   (begp:endp)) ; this%leafn0_xfer_patch                   (:) = nan
       allocate(this%frootn0_patch                       (begp:endp)) ; this%frootn0_patch                       (:) = nan
       allocate(this%frootn0_storage_patch               (begp:endp)) ; this%frootn0_storage_patch               (:) = nan
       allocate(this%frootn0_xfer_patch                  (begp:endp)) ; this%frootn0_xfer_patch                  (:) = nan
       allocate(this%livestemn0_patch                    (begp:endp)) ; this%livestemn0_patch                    (:) = nan
       allocate(this%livestemn0_storage_patch            (begp:endp)) ; this%livestemn0_storage_patch            (:) = nan
       allocate(this%livestemn0_xfer_patch               (begp:endp)) ; this%livestemn0_xfer_patch               (:) = nan
       allocate(this%deadstemn0_patch                    (begp:endp)) ; this%deadstemn0_patch                    (:) = nan
       allocate(this%deadstemn0_storage_patch            (begp:endp)) ; this%deadstemn0_storage_patch            (:) = nan
       allocate(this%deadstemn0_xfer_patch               (begp:endp)) ; this%deadstemn0_xfer_patch               (:) = nan
       allocate(this%livecrootn0_patch                   (begp:endp)) ; this%livecrootn0_patch                   (:) = nan
       allocate(this%livecrootn0_storage_patch           (begp:endp)) ; this%livecrootn0_storage_patch           (:) = nan
       allocate(this%livecrootn0_xfer_patch              (begp:endp)) ; this%livecrootn0_xfer_patch              (:) = nan
       allocate(this%deadcrootn0_patch                   (begp:endp)) ; this%deadcrootn0_patch                   (:) = nan
       allocate(this%deadcrootn0_storage_patch           (begp:endp)) ; this%deadcrootn0_storage_patch           (:) = nan
       allocate(this%deadcrootn0_xfer_patch              (begp:endp)) ; this%deadcrootn0_xfer_patch              (:) = nan
       allocate(this%grainn0_patch                       (begp:endp)) ; this%grainn0_patch                       (:) = nan
       allocate(this%grainn0_storage_patch               (begp:endp)) ; this%grainn0_storage_patch               (:) = nan
       allocate(this%grainn0_xfer_patch                  (begp:endp)) ; this%grainn0_xfer_patch                  (:) = nan
       allocate(this%retransn0_patch                     (begp:endp)) ; this%retransn0_patch                     (:) = nan

       allocate(this%leafn_SASUsave_patch                (begp:endp)) ; this%leafn_SASUsave_patch               (:) = nan
       allocate(this%leafn_storage_SASUsave_patch        (begp:endp)) ; this%leafn_storage_SASUsave_patch       (:) = nan
       allocate(this%leafn_xfer_SASUsave_patch           (begp:endp)) ; this%leafn_xfer_SASUsave_patch          (:) = nan
       allocate(this%frootn_SASUsave_patch               (begp:endp)) ; this%frootn_SASUsave_patch              (:) = nan
       allocate(this%frootn_storage_SASUsave_patch       (begp:endp)) ; this%frootn_storage_SASUsave_patch      (:) = nan
       allocate(this%frootn_xfer_SASUsave_patch          (begp:endp)) ; this%frootn_xfer_SASUsave_patch         (:) = nan
       allocate(this%livestemn_SASUsave_patch            (begp:endp)) ; this%livestemn_SASUsave_patch           (:) = nan
       allocate(this%livestemn_storage_SASUsave_patch    (begp:endp)) ; this%livestemn_storage_SASUsave_patch   (:) = nan
       allocate(this%livestemn_xfer_SASUsave_patch       (begp:endp)) ; this%livestemn_xfer_SASUsave_patch      (:) = nan
       allocate(this%deadstemn_SASUsave_patch            (begp:endp)) ; this%deadstemn_SASUsave_patch           (:) = nan
       allocate(this%deadstemn_storage_SASUsave_patch    (begp:endp)) ; this%deadstemn_storage_SASUsave_patch   (:) = nan
       allocate(this%deadstemn_xfer_SASUsave_patch       (begp:endp)) ; this%deadstemn_xfer_SASUsave_patch      (:) = nan
       allocate(this%livecrootn_SASUsave_patch           (begp:endp)) ; this%livecrootn_SASUsave_patch          (:) = nan
       allocate(this%livecrootn_storage_SASUsave_patch   (begp:endp)) ; this%livecrootn_storage_SASUsave_patch  (:) = nan
       allocate(this%livecrootn_xfer_SASUsave_patch      (begp:endp)) ; this%livecrootn_xfer_SASUsave_patch     (:) = nan
       allocate(this%deadcrootn_SASUsave_patch           (begp:endp)) ; this%deadcrootn_SASUsave_patch          (:) = nan
       allocate(this%deadcrootn_storage_SASUsave_patch   (begp:endp)) ; this%deadcrootn_storage_SASUsave_patch  (:) = nan
       allocate(this%deadcrootn_xfer_SASUsave_patch      (begp:endp)) ; this%deadcrootn_xfer_SASUsave_patch     (:) = nan
       allocate(this%grainn_SASUsave_patch               (begp:endp)) ; this%grainn_SASUsave_patch              (:) = nan
       allocate(this%grainn_storage_SASUsave_patch       (begp:endp)) ; this%grainn_storage_SASUsave_patch      (:) = nan

       allocate(this%matrix_nalloc_leaf_acc_patch        (begp:endp)) ; this%matrix_nalloc_leaf_acc_patch        (:) = nan
       allocate(this%matrix_nalloc_leafst_acc_patch      (begp:endp)) ; this%matrix_nalloc_leafst_acc_patch      (:) = nan
       allocate(this%matrix_nalloc_froot_acc_patch       (begp:endp)) ; this%matrix_nalloc_froot_acc_patch       (:) = nan
       allocate(this%matrix_nalloc_frootst_acc_patch     (begp:endp)) ; this%matrix_nalloc_frootst_acc_patch     (:) = nan
       allocate(this%matrix_nalloc_livestem_acc_patch    (begp:endp)) ; this%matrix_nalloc_livestem_acc_patch    (:) = nan
       allocate(this%matrix_nalloc_livestemst_acc_patch  (begp:endp)) ; this%matrix_nalloc_livestemst_acc_patch  (:) = nan
       allocate(this%matrix_nalloc_deadstem_acc_patch    (begp:endp)) ; this%matrix_nalloc_deadstem_acc_patch    (:) = nan
       allocate(this%matrix_nalloc_deadstemst_acc_patch  (begp:endp)) ; this%matrix_nalloc_deadstemst_acc_patch  (:) = nan
       allocate(this%matrix_nalloc_livecroot_acc_patch   (begp:endp)) ; this%matrix_nalloc_livecroot_acc_patch   (:) = nan
       allocate(this%matrix_nalloc_livecrootst_acc_patch (begp:endp)) ; this%matrix_nalloc_livecrootst_acc_patch (:) = nan
       allocate(this%matrix_nalloc_deadcroot_acc_patch   (begp:endp)) ; this%matrix_nalloc_deadcroot_acc_patch   (:) = nan
       allocate(this%matrix_nalloc_deadcrootst_acc_patch (begp:endp)) ; this%matrix_nalloc_deadcrootst_acc_patch (:) = nan
       allocate(this%matrix_nalloc_grain_acc_patch       (begp:endp)) ; this%matrix_nalloc_grain_acc_patch       (:) = nan
       allocate(this%matrix_nalloc_grainst_acc_patch     (begp:endp)) ; this%matrix_nalloc_grainst_acc_patch     (:) = nan

       allocate(this%matrix_ntransfer_leafst_to_leafxf_acc_patch           (begp:endp)) ; this%matrix_ntransfer_leafst_to_leafxf_acc_patch           (:) = nan
       allocate(this%matrix_ntransfer_leafxf_to_leaf_acc_patch             (begp:endp)) ; this%matrix_ntransfer_leafxf_to_leaf_acc_patch             (:) = nan
       allocate(this%matrix_ntransfer_frootst_to_frootxf_acc_patch         (begp:endp)) ; this%matrix_ntransfer_frootst_to_frootxf_acc_patch         (:) = nan
       allocate(this%matrix_ntransfer_frootxf_to_froot_acc_patch           (begp:endp)) ; this%matrix_ntransfer_frootxf_to_froot_acc_patch           (:) = nan
       allocate(this%matrix_ntransfer_livestemst_to_livestemxf_acc_patch   (begp:endp)) ; this%matrix_ntransfer_livestemst_to_livestemxf_acc_patch   (:) = nan
       allocate(this%matrix_ntransfer_livestemxf_to_livestem_acc_patch     (begp:endp)) ; this%matrix_ntransfer_livestemxf_to_livestem_acc_patch     (:) = nan
       allocate(this%matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch   (begp:endp)) ; this%matrix_ntransfer_deadstemst_to_deadstemxf_acc_patch   (:) = nan
       allocate(this%matrix_ntransfer_deadstemxf_to_deadstem_acc_patch     (begp:endp)) ; this%matrix_ntransfer_deadstemxf_to_deadstem_acc_patch     (:) = nan
       allocate(this%matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch (begp:endp)) ; this%matrix_ntransfer_livecrootst_to_livecrootxf_acc_patch (:) = nan
       allocate(this%matrix_ntransfer_livecrootxf_to_livecroot_acc_patch   (begp:endp)) ; this%matrix_ntransfer_livecrootxf_to_livecroot_acc_patch   (:) = nan
       allocate(this%matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch (begp:endp)) ; this%matrix_ntransfer_deadcrootst_to_deadcrootxf_acc_patch (:) = nan
       allocate(this%matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch   (begp:endp)) ; this%matrix_ntransfer_deadcrootxf_to_deadcroot_acc_patch   (:) = nan
       allocate(this%matrix_ntransfer_grainst_to_grainxf_acc_patch         (begp:endp)) ; this%matrix_ntransfer_grainst_to_grainxf_acc_patch         (:) = nan
       allocate(this%matrix_ntransfer_grainxf_to_grain_acc_patch           (begp:endp)) ; this%matrix_ntransfer_grainxf_to_grain_acc_patch           (:) = nan
       allocate(this%matrix_ntransfer_livestem_to_deadstem_acc_patch       (begp:endp)) ; this%matrix_ntransfer_livestem_to_deadstem_acc_patch       (:) = nan
       allocate(this%matrix_ntransfer_livecroot_to_deadcroot_acc_patch     (begp:endp)) ; this%matrix_ntransfer_livecroot_to_deadcroot_acc_patch     (:) = nan

       allocate(this%matrix_ntransfer_retransn_to_leaf_acc_patch           (begp:endp)) ; this%matrix_ntransfer_retransn_to_leaf_acc_patch           (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_leafst_acc_patch         (begp:endp)) ; this%matrix_ntransfer_retransn_to_leafst_acc_patch         (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_froot_acc_patch          (begp:endp)) ; this%matrix_ntransfer_retransn_to_froot_acc_patch          (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_frootst_acc_patch        (begp:endp)) ; this%matrix_ntransfer_retransn_to_frootst_acc_patch        (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_livestem_acc_patch       (begp:endp)) ; this%matrix_ntransfer_retransn_to_livestem_acc_patch       (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_livestemst_acc_patch     (begp:endp)) ; this%matrix_ntransfer_retransn_to_livestemst_acc_patch     (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_deadstem_acc_patch       (begp:endp)) ; this%matrix_ntransfer_retransn_to_deadstem_acc_patch       (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_deadstemst_acc_patch     (begp:endp)) ; this%matrix_ntransfer_retransn_to_deadstemst_acc_patch     (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_livecroot_acc_patch      (begp:endp)) ; this%matrix_ntransfer_retransn_to_livecroot_acc_patch      (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_livecrootst_acc_patch    (begp:endp)) ; this%matrix_ntransfer_retransn_to_livecrootst_acc_patch    (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_deadcroot_acc_patch      (begp:endp)) ; this%matrix_ntransfer_retransn_to_deadcroot_acc_patch      (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_deadcrootst_acc_patch    (begp:endp)) ; this%matrix_ntransfer_retransn_to_deadcrootst_acc_patch    (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_grain_acc_patch          (begp:endp)) ; this%matrix_ntransfer_retransn_to_grain_acc_patch          (:) = nan
       allocate(this%matrix_ntransfer_retransn_to_grainst_acc_patch        (begp:endp)) ; this%matrix_ntransfer_retransn_to_grainst_acc_patch        (:) = nan

       allocate(this%matrix_ntransfer_leaf_to_retransn_acc_patch           (begp:endp)) ; this%matrix_ntransfer_leaf_to_retransn_acc_patch           (:) = nan
       allocate(this%matrix_ntransfer_froot_to_retransn_acc_patch          (begp:endp)) ; this%matrix_ntransfer_froot_to_retransn_acc_patch          (:) = nan
       allocate(this%matrix_ntransfer_livestem_to_retransn_acc_patch       (begp:endp)) ; this%matrix_ntransfer_livestem_to_retransn_acc_patch       (:) = nan
       allocate(this%matrix_ntransfer_livecroot_to_retransn_acc_patch      (begp:endp)) ; this%matrix_ntransfer_livecroot_to_retransn_acc_patch      (:) = nan

       allocate(this%matrix_nturnover_leaf_acc_patch                       (begp:endp)) ; this%matrix_nturnover_leaf_acc_patch                       (:) = nan
       allocate(this%matrix_nturnover_leafst_acc_patch                     (begp:endp)) ; this%matrix_nturnover_leafst_acc_patch                     (:) = nan
       allocate(this%matrix_nturnover_leafxf_acc_patch                     (begp:endp)) ; this%matrix_nturnover_leafxf_acc_patch                     (:) = nan
       allocate(this%matrix_nturnover_froot_acc_patch                      (begp:endp)) ; this%matrix_nturnover_froot_acc_patch                      (:) = nan
       allocate(this%matrix_nturnover_frootst_acc_patch                    (begp:endp)) ; this%matrix_nturnover_frootst_acc_patch                    (:) = nan
       allocate(this%matrix_nturnover_frootxf_acc_patch                    (begp:endp)) ; this%matrix_nturnover_frootxf_acc_patch                    (:) = nan
       allocate(this%matrix_nturnover_livestem_acc_patch                   (begp:endp)) ; this%matrix_nturnover_livestem_acc_patch                   (:) = nan
       allocate(this%matrix_nturnover_livestemst_acc_patch                 (begp:endp)) ; this%matrix_nturnover_livestemst_acc_patch                 (:) = nan
       allocate(this%matrix_nturnover_livestemxf_acc_patch                 (begp:endp)) ; this%matrix_nturnover_livestemxf_acc_patch                 (:) = nan
       allocate(this%matrix_nturnover_deadstem_acc_patch                   (begp:endp)) ; this%matrix_nturnover_deadstem_acc_patch                   (:) = nan
       allocate(this%matrix_nturnover_deadstemst_acc_patch                 (begp:endp)) ; this%matrix_nturnover_deadstemst_acc_patch                 (:) = nan
       allocate(this%matrix_nturnover_deadstemxf_acc_patch                 (begp:endp)) ; this%matrix_nturnover_deadstemxf_acc_patch                 (:) = nan
       allocate(this%matrix_nturnover_livecroot_acc_patch                  (begp:endp)) ; this%matrix_nturnover_livecroot_acc_patch                  (:) = nan
       allocate(this%matrix_nturnover_livecrootst_acc_patch                (begp:endp)) ; this%matrix_nturnover_livecrootst_acc_patch                (:) = nan
       allocate(this%matrix_nturnover_livecrootxf_acc_patch                (begp:endp)) ; this%matrix_nturnover_livecrootxf_acc_patch                (:) = nan
       allocate(this%matrix_nturnover_deadcroot_acc_patch                  (begp:endp)) ; this%matrix_nturnover_deadcroot_acc_patch                  (:) = nan
       allocate(this%matrix_nturnover_deadcrootst_acc_patch                (begp:endp)) ; this%matrix_nturnover_deadcrootst_acc_patch                (:) = nan
       allocate(this%matrix_nturnover_deadcrootxf_acc_patch                (begp:endp)) ; this%matrix_nturnover_deadcrootxf_acc_patch                (:) = nan
       allocate(this%matrix_nturnover_grain_acc_patch                      (begp:endp)) ; this%matrix_nturnover_grain_acc_patch                      (:) = nan
       allocate(this%matrix_nturnover_grainst_acc_patch                    (begp:endp)) ; this%matrix_nturnover_grainst_acc_patch                    (:) = nan
       allocate(this%matrix_nturnover_grainxf_acc_patch                    (begp:endp)) ; this%matrix_nturnover_grainxf_acc_patch                    (:) = nan
       allocate(this%matrix_nturnover_retransn_acc_patch                   (begp:endp)) ; this%matrix_nturnover_retransn_acc_patch                   (:) = nan
    end if

  ! initialize arrays with values from restarts

    n  = 0
    np = 0
    do nc = 1,nch        ! catchment tile loop                                                                 
       
       this%seedn_grc(nc) = 0

       do nz = 1,NUM_ZON    ! CN zone loop
          n = n + 1

          this%seedn_grc(nc) = this%seedn_grc(nc) + cncol(nc,nz,23)*CN_zone_weight(nz)
          this%totn_col(n) = cncol(nc,nz,29)
                                                                        
          do p = 0,numpft  ! PFT index loop                                                                        
             np = np + 1
             do nv = 1,NUM_VEG ! defined veg loop                                                                      
                if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then

                     this%deadcrootn_patch         (np) = cnpft(nc,nz,nv, 48)
                     this%deadcrootn_storage_patch (np) = cnpft(nc,nz,nv, 49)
                     this%deadcrootn_xfer_patch    (np) = cnpft(nc,nz,nv, 50)
                     this%deadstemn_patch          (np) = cnpft(nc,nz,nv, 51)
                     this%deadstemn_storage_patch  (np) = cnpft(nc,nz,nv, 52)
                     this%deadstemn_xfer_patch     (np) = cnpft(nc,nz,nv, 53)
                     this%frootn_patch             (np) = cnpft(nc,nz,nv, 54)
                     this%frootn_storage_patch     (np) = cnpft(nc,nz,nv, 55)
                     this%frootn_xfer_patch        (np) = cnpft(nc,nz,nv, 56)
                     this%leafn_patch              (np) = cnpft(nc,nz,nv, 57)
                     this%leafn_storage_patch      (np) = cnpft(nc,nz,nv, 58)
                     this%leafn_xfer_patch         (np) = cnpft(nc,nz,nv, 59)
                     this%livecrootn_patch         (np) = cnpft(nc,nz,nv, 60)
                     this%livecrootn_storage_patch (np) = cnpft(nc,nz,nv, 61)
                     this%livecrootn_xfer_patch    (np) = cnpft(nc,nz,nv, 62)
                     this%livestemn_patch          (np) = cnpft(nc,nz,nv, 63)
                     this%livestemn_storage_patch  (np) = cnpft(nc,nz,nv, 64)
                     this%livestemn_xfer_patch     (np) = cnpft(nc,nz,nv, 65)
                     this%npool_patch              (np) = cnpft(nc,nz,nv, 66)
                     this%ntrunc_patch             (np) = cnpft(nc,nz,nv, 67)
                     this%retransn_patch           (np) = cnpft(nc,nz,nv, 68)

                end if
              end do !nv                                                                                            
          end do !p                                                                                               
       end do !nz                                                                                               
    end do ! nc  

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine Summary_nitrogenstate(this, bounds, num_allc, filter_allc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp,&
       soilbiogeochem_nitrogenstate_inst)
    !
    ! !USES:
    use subgridAveMod, only : p2c
    use SoilBiogeochemNitrogenStateType, only : soilbiogeochem_nitrogenstate_type
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type)                      :: this
    type(bounds_type)                       , intent(in) :: bounds
    integer                                 , intent(in) :: num_allc        ! number of columns in allc filter
    integer                                 , intent(in) :: filter_allc(:)  ! filter for all active columns
    integer                                 , intent(in) :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in) :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in) :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in) :: filter_soilp(:) ! filter for soil patches
    type(soilbiogeochem_nitrogenstate_type) , intent(in) :: soilbiogeochem_nitrogenstate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p,j,k,l ! indices
    integer  :: fp,fc       ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

    ! --------------------------------------------
    ! patch level summary
    ! --------------------------------------------

    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
       this%dispvegn_patch(p) = &
            this%leafn_patch(p)      + &
            this%frootn_patch(p)     + &
            this%livestemn_patch(p)  + &
            this%deadstemn_patch(p)  + &
            this%livecrootn_patch(p) + &
            this%deadcrootn_patch(p)

       ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
       this%storvegn_patch(p) = &
            this%leafn_storage_patch(p)      + &
            this%frootn_storage_patch(p)     + &
            this%livestemn_storage_patch(p)  + &
            this%deadstemn_storage_patch(p)  + &
            this%livecrootn_storage_patch(p) + &
            this%deadcrootn_storage_patch(p) + &
            this%leafn_xfer_patch(p)         + &
            this%frootn_xfer_patch(p)        + &
            this%livestemn_xfer_patch(p)     + &
            this%deadstemn_xfer_patch(p)     + &
            this%livecrootn_xfer_patch(p)    + &
            this%deadcrootn_xfer_patch(p)    + &
            this%npool_patch(p)              + &
            this%retransn_patch(p)

       if ( use_crop .and. patch%itype(p) >= npcropmin )then
          this%dispvegn_patch(p) = &
               this%dispvegn_patch(p) + &
               this%grainn_patch(p)

          this%storvegn_patch(p) = &
               this%storvegn_patch(p) + &
               this%grainn_storage_patch(p)     + &
               this%grainn_xfer_patch(p) + &
               this%cropseedn_deficit_patch(p)
       end if

       ! total vegetation nitrogen (TOTVEGN)
       this%totvegn_patch(p) = &
            this%dispvegn_patch(p) + &
            this%storvegn_patch(p)

       ! total patch-level carbon (add ntrunc)
       this%totn_patch(p) = &
            this%totvegn_patch(p) + &
            this%ntrunc_patch(p)

    end do

    ! --------------------------------------------
    ! column level summary
    ! --------------------------------------------

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totvegn_patch(bounds%begp:bounds%endp), &
         this%totvegn_col(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         this%totn_patch(bounds%begp:bounds%endp), &
         this%totn_p2c_col(bounds%begc:bounds%endc))


    do fc = 1,num_allc
       c = filter_allc(fc)

       ! total ecosystem nitrogen, including veg (TOTECOSYSN)
       this%totecosysn_col(c) =    &
            soilbiogeochem_nitrogenstate_inst%cwdn_col(c)    + &
            soilbiogeochem_nitrogenstate_inst%totlitn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%totsomn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%sminn_col(c)   + &
            this%totvegn_col(c)

       ! total column nitrogen, including patch (TOTCOLN)

       this%totn_col(c) = this%totn_p2c_col(c)               + &
            soilbiogeochem_nitrogenstate_inst%cwdn_col(c)    + &
            soilbiogeochem_nitrogenstate_inst%totlitn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%totsomn_col(c) + &
            soilbiogeochem_nitrogenstate_inst%sminn_col(c)   + &
            soilbiogeochem_nitrogenstate_inst%ntrunc_col(c)

    end do

  end subroutine Summary_nitrogenstate

  !-----------------------------------------------------------------------
  subroutine ZeroDwt( this, bounds )
    !
    ! !DESCRIPTION
    ! Initialize variables needed for dynamic land use.
    !
    ! !ARGUMENTS:
    class(cnveg_nitrogenstate_type) :: this
    type(bounds_type), intent(in)  :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: p          ! indices
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%dispvegn_patch(p) = 0._r8
       this%storvegn_patch(p) = 0._r8
       this%totvegn_patch(p)  = 0._r8
       this%totn_patch(p)     = 0._r8
    end do

  end subroutine ZeroDwt

end module CNVegNitrogenStateType

