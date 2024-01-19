module clmtype

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION: 
! Define derived type hierarchy. Includes declaration of
! the clm derived type and 1d mapping arrays. 
!
! -------------------------------------------------------- 
! gridcell types can have values of 
! -------------------------------------------------------- 
!   1 => default
! -------------------------------------------------------- 
! landunits types can have values of (see clm_varcon.F90)
! -------------------------------------------------------- 
!   1  => (istsoil) soil (vegetated or bare soil landunit)
!   2  => (istice)  land ice
!   3  => (istdlak) deep lake
!   4  => (istslak) shall lake (not currently implemented)
!   5  => (istwet)  wetland
!   6  => (isturb)  urban 
! -------------------------------------------------------- 
! column types can have values of
! -------------------------------------------------------- 
!   1  => (istsoil)          soil (vegetated or bare soil)
!   2  => (istice)           land ice
!   3  => (istdlak)          deep lake
!   4  => (istslak)          shallow lake 
!   5  => (istwet)           wetland
! -------------------------------------------------------- 
! pft types can have values of
! -------------------------------------------------------- 
!   0  => not vegetated
!   1  => needleleaf evergreen temperate tree
!   2  => needleleaf evergreen boreal tree
!   3  => needleleaf deciduous boreal tree
!   4  => broadleaf evergreen tropical tree
!   5  => broadleaf evergreen temperate tree
!   6  => broadleaf deciduous tropical tree
!   7  => broadleaf deciduous temperate tree
!   8  => broadleaf deciduous boreal tree
!   9  => broadleaf evergreen shrub
!   10 => broadleaf deciduous temperate shrub
!   11 => broadleaf deciduous boreal shrub
!   12 => c3 arctic grass
!   13 => c3 non-arctic grass
!   14 => c4 grass
!   15 => corn
!   16 => wheat
! -------------------------------------------------------- 
!
! !USES:
!
! !PUBLIC TYPES:
  implicit none

  private
!                              
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!*******************************************************************************
!----------------------------------------------------
! Begin definition of conservation check structures
!----------------------------------------------------

!----------------------------------------------------
! carbon balance structure
!----------------------------------------------------
type, public :: carbon_balance_type
   real, pointer :: begcb(:)         !carbon mass, beginning of time step (gC/m**2)
   real, pointer :: endcb(:)         !carbon mass, end of time step (gC/m**2)
   real, pointer :: errcb(:)         !carbon balance error for the timestep (gC/m**2)
end type carbon_balance_type

!----------------------------------------------------
! nitrogen balance structure
!----------------------------------------------------
type, public :: nitrogen_balance_type
   real, pointer :: begnb(:)         !nitrogen mass, beginning of time step (gN/m**2)
   real, pointer :: endnb(:)         !nitrogen mass, end of time step (gN/m**2)
   real, pointer :: errnb(:)         !nitrogen balance error for the timestep (gN/m**2)
end type nitrogen_balance_type

!----------------------------------------------------
! End definition of conservation check structures
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the pft_type level
!----------------------------------------------------
! pft physical state variables structure
!----------------------------------------------------
type, public :: pft_pstate_type
   integer , pointer :: frac_veg_nosno(:)       !fraction of vegetation not covered by snow (0 OR 1) [-] 
   integer , pointer :: frac_veg_nosno_alb(:)   !fraction of vegetation not covered by snow (0 OR 1) [-] 
   real, pointer :: rootfr(:,:)             !fraction of roots in each soil layer  (nlevgrnd)
   real, pointer :: laisun(:)               !sunlit projected leaf area index
   real, pointer :: laisha(:)               !shaded projected leaf area index
   real, pointer :: fsun(:)                 !sunlit fraction of canopy
   real, pointer :: tlai(:)                 !one-sided leaf area index, no burying by snow
   real, pointer :: tsai(:)                 !one-sided stem area index, no burying by snow
   real, pointer :: elai(:)                 !one-sided leaf area index with burying by snow
   real, pointer :: esai(:)                 !one-sided stem area index with burying by snow
   real, pointer :: fwet(:)                 !fraction of canopy that is wet (0 to 1)
   real, pointer :: fdry(:)                 !fraction of foliage that is green and dry [-] (new)
   real, pointer :: htop(:)                 !canopy top (m)
   real, pointer :: hbot(:)                 !canopy bottom (m)
   real, pointer :: forc_hgt_u_pft(:)       !wind forcing height (10m+z0m+d) (m)
end type pft_pstate_type

!----------------------------------------------------
! pft ecophysiological constants structure
!----------------------------------------------------
type, public :: pft_epc_type
   integer , pointer :: noveg(:)                !value for not vegetated
   integer , pointer :: tree(:)                 !tree or not?
   real, pointer :: fnitr(:)                !foliage nitrogen limitation factor (-)
   real, pointer :: c3psn(:)                !photosynthetic pathway: 0. = c4, 1. = c3
   real, pointer :: vcmx25(:)               !max rate of carboxylation at 25C (umol CO2/m**2/s)
   real, pointer :: mp(:)                   !slope of conductance-to-photosynthesis relationship
   real, pointer :: qe25(:)                 !quantum efficiency at 25C (umol CO2 / umol photon)
   real, pointer :: z0mr(:)                 !ratio of momentum roughness length to canopy top height (-)
   real, pointer :: displar(:)              !ratio of displacement height to canopy top height (-)
   ! new variables for CN code
   real, pointer :: dwood(:)           !wood density (gC/m3)
   real, pointer :: slatop(:)    !specific leaf area at top of canopy, projected area basis [m^2/gC]
   real, pointer :: dsladlai(:)  !dSLA/dLAI, projected area basis [m^2/gC]
   real, pointer :: leafcn(:)    !leaf C:N (gC/gN)
   real, pointer :: flnr(:)      !fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
   real, pointer :: woody(:)     !binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: lflitcn(:)      !leaf litter C:N (gC/gN)
   real, pointer :: frootcn(:)      !fine root C:N (gC/gN)
   real, pointer :: livewdcn(:)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
   real, pointer :: deadwdcn(:)     !dead wood (xylem and heartwood) C:N (gC/gN)
   real, pointer :: froot_leaf(:)   !allocation parameter: new fine root C per new leaf C (gC/gC)
   real, pointer :: stem_leaf(:)    !allocation parameter: new stem c per new leaf C (gC/gC)
   real, pointer :: croot_stem(:)   !allocation parameter: new coarse root C per new stem C (gC/gC)
   real, pointer :: flivewd(:)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
   real, pointer :: fcur(:)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
   real, pointer :: lf_flab(:)      !leaf litter labile fraction
   real, pointer :: lf_fcel(:)      !leaf litter cellulose fraction
   real, pointer :: lf_flig(:)      !leaf litter lignin fraction
   real, pointer :: fr_flab(:)      !fine root litter labile fraction
   real, pointer :: fr_fcel(:)      !fine root litter cellulose fraction
   real, pointer :: fr_flig(:)      !fine root litter lignin fraction
   real, pointer :: dw_fcel(:)      !dead wood cellulose fraction
   real, pointer :: dw_flig(:)      !dead wood lignin fraction
   real, pointer :: leaf_long(:)    !leaf longevity (yrs)
   real, pointer :: evergreen(:)    !binary flag for evergreen leaf habit (0 or 1)
   real, pointer :: stress_decid(:) !binary flag for stress-deciduous leaf habit (0 or 1)
   real, pointer :: season_decid(:) !binary flag for seasonal-deciduous leaf habit (0 or 1)
   real, pointer :: resist(:)       !resistance to fire (no units)
   ! gkw: moved these here for lack of a better place
   real, pointer :: xl(:)           ! leaf/stem orientation index
   real, pointer :: rhol(:)         ! leaf reflectance (visible)
   real, pointer :: rhos(:)         ! stem reflectance (visible)
   real, pointer :: taul(:)         ! leaf transmittance (visible)
   real, pointer :: taus(:)         ! stem transmittance (visible)
end type pft_epc_type

!----------------------------------------------------
! pft ecophysiological variables structure
!----------------------------------------------------
type, public :: pft_epv_type
   real, pointer :: dormant_flag(:)         !dormancy flag
   real, pointer :: days_active(:)          !number of days since last dormancy
   real, pointer :: onset_flag(:)           !onset flag
   real, pointer :: onset_counter(:)        !onset days counter
   real, pointer :: onset_gddflag(:)        !onset flag for growing degree day sum
   real, pointer :: onset_fdd(:)            !onset freezing degree days counter
   real, pointer :: onset_gdd(:)            !onset growing degree days
   real, pointer :: onset_swi(:)            !onset soil water index
   real, pointer :: offset_flag(:)          !offset flag
   real, pointer :: offset_counter(:)       !offset days counter
   real, pointer :: offset_fdd(:)           !offset freezing degree days counter
   real, pointer :: offset_swi(:)           !offset soil water index
   real, pointer :: lgsf(:)                 !long growing season factor [0-1]
   real, pointer :: bglfr(:)                !background litterfall rate (1/s)
   real, pointer :: bgtr(:)                 !background transfer growth rate (1/s)
   real, pointer :: dayl(:)                 !daylength (seconds)
   real, pointer :: prev_dayl(:)            !daylength from previous timestep (seconds)
   real, pointer :: annavg_t2m(:)           !annual average 2m air temperature (K)
   real, pointer :: tempavg_t2m(:)          !temporary average 2m air temperature (K)
   real, pointer :: gpp(:)                  !GPP flux before downregulation (gC/m2/s)
   real, pointer :: availc(:)               !C flux available for allocation (gC/m2/s)
   real, pointer :: xsmrpool_recover(:)     !C flux assigned to recovery of negative cpool (gC/m2/s)
   real, pointer :: alloc_pnow(:)           !fraction of current allocation to display as new growth (DIM)
   real, pointer :: c_allometry(:)          !C allocation index (DIM)
   real, pointer :: n_allometry(:)          !N allocation index (DIM)
   real, pointer :: plant_ndemand(:)        !N flux required to support initial GPP (gN/m2/s)
   real, pointer :: tempsum_potential_gpp(:)!temporary annual sum of potential GPP
   real, pointer :: annsum_potential_gpp(:) !annual sum of potential GPP
   real, pointer :: tempmax_retransn(:)     !temporary annual max of retranslocated N pool (gN/m2)
   real, pointer :: annmax_retransn(:)      !annual max of retranslocated N pool (gN/m2)
   real, pointer :: avail_retransn(:)       !N flux available from retranslocation pool (gN/m2/s)
   real, pointer :: plant_nalloc(:)         !total allocated N flux (gN/m2/s)
   real, pointer :: plant_calloc(:)         !total allocated C flux (gC/m2/s)
   real, pointer :: excess_cflux(:)         !C flux not allocated due to downregulation (gC/m2/s)
   real, pointer :: downreg(:)              !fractional reduction in GPP due to N limitation (DIM)
   real, pointer :: prev_leafc_to_litter(:) !previous timestep leaf C litterfall flux (gC/m2/s)
   real, pointer :: prev_frootc_to_litter(:)!previous timestep froot C litterfall flux (gC/m2/s)
   real, pointer :: tempsum_npp(:)          !temporary annual sum of NPP (gC/m2/yr)
   real, pointer :: annsum_npp(:)           !annual sum of NPP (gC/m2/yr)
end type pft_epv_type                        

!----------------------------------------------------
! pft energy state variables structure
!----------------------------------------------------
type, public :: pft_estate_type
   real, pointer :: t_ref2m(:)            !2 m height surface air temperature (Kelvin)
end type pft_estate_type

!----------------------------------------------------
! pft carbon state variables structure
!----------------------------------------------------
type, public :: pft_cstate_type
   real, pointer :: leafcmax(:)           ! (gC/m2) ann max leaf C
   real, pointer :: leafc(:)              ! (gC/m2) leaf C
   real, pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real, pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real, pointer :: frootc(:)             ! (gC/m2) fine root C
   real, pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real, pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real, pointer :: livestemc(:)          ! (gC/m2) live stem C
   real, pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real, pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real, pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real, pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real, pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real, pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real, pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real, pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real, pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real, pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real, pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real, pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real, pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real, pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real, pointer :: xsmrpool(:)           ! (gC/m2) abstract C pool to meet excess MR demand
   real, pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   ! summary (diagnostic) state variables, not involved in mass balance
   real, pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real, pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real, pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real, pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
end type pft_cstate_type

!----------------------------------------------------
! pft nitrogen state variables structure
!----------------------------------------------------
type, public :: pft_nstate_type
   real, pointer :: leafn(:)              ! (gN/m2) leaf N 
   real, pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real, pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real, pointer :: frootn(:)             ! (gN/m2) fine root N
   real, pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real, pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real, pointer :: livestemn(:)          ! (gN/m2) live stem N
   real, pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real, pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real, pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real, pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real, pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real, pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real, pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real, pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real, pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real, pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real, pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real, pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real, pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real, pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   ! summary (diagnostic) state variables, not involved in mass balance
   real, pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
   real, pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
   real, pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
   real, pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
end type pft_nstate_type

!----------------------------------------------------
! pft carbon flux variables structure
!----------------------------------------------------
type, public :: pft_cflux_type
   real, pointer :: psnsun(:)         !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real, pointer :: psnsha(:)         !shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real, pointer :: fco2(:)           !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
   ! new variables for CN code
   ! gap mortality fluxes
   real, pointer :: m_leafc_to_litter(:)                 ! leaf C mortality (gC/m2/s)
   real, pointer :: m_leafc_storage_to_litter(:)         ! leaf C storage mortality (gC/m2/s)
   real, pointer :: m_leafc_xfer_to_litter(:)            ! leaf C transfer mortality (gC/m2/s)
   real, pointer :: m_frootc_to_litter(:)                ! fine root C mortality (gC/m2/s)
   real, pointer :: m_frootc_storage_to_litter(:)        ! fine root C storage mortality (gC/m2/s)
   real, pointer :: m_frootc_xfer_to_litter(:)           ! fine root C transfer mortality (gC/m2/s)
   real, pointer :: m_livestemc_to_litter(:)             ! live stem C mortality (gC/m2/s)
   real, pointer :: m_livestemc_storage_to_litter(:)     ! live stem C storage mortality (gC/m2/s)
   real, pointer :: m_livestemc_xfer_to_litter(:)        ! live stem C transfer mortality (gC/m2/s)
   real, pointer :: m_deadstemc_to_litter(:)             ! dead stem C mortality (gC/m2/s)
   real, pointer :: m_deadstemc_storage_to_litter(:)     ! dead stem C storage mortality (gC/m2/s)
   real, pointer :: m_deadstemc_xfer_to_litter(:)        ! dead stem C transfer mortality (gC/m2/s)
   real, pointer :: m_livecrootc_to_litter(:)            ! live coarse root C mortality (gC/m2/s)
   real, pointer :: m_livecrootc_storage_to_litter(:)    ! live coarse root C storage mortality (gC/m2/s)
   real, pointer :: m_livecrootc_xfer_to_litter(:)       ! live coarse root C transfer mortality (gC/m2/s)
   real, pointer :: m_deadcrootc_to_litter(:)            ! dead coarse root C mortality (gC/m2/s)
   real, pointer :: m_deadcrootc_storage_to_litter(:)    ! dead coarse root C storage mortality (gC/m2/s)
   real, pointer :: m_deadcrootc_xfer_to_litter(:)       ! dead coarse root C transfer mortality (gC/m2/s)
   real, pointer :: m_gresp_storage_to_litter(:)         ! growth respiration storage mortality (gC/m2/s)
   real, pointer :: m_gresp_xfer_to_litter(:)            ! growth respiration transfer mortality (gC/m2/s)
   ! harvest mortality fluxes
   real, pointer :: hrv_leafc_to_litter(:)               ! leaf C harvest mortality (gC/m2/s)
   real, pointer :: hrv_leafc_storage_to_litter(:)       ! leaf C storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_leafc_xfer_to_litter(:)          ! leaf C transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_frootc_to_litter(:)              ! fine root C harvest mortality (gC/m2/s)
   real, pointer :: hrv_frootc_storage_to_litter(:)      ! fine root C storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_frootc_xfer_to_litter(:)         ! fine root C transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_livestemc_to_litter(:)           ! live stem C harvest mortality (gC/m2/s)
   real, pointer :: hrv_livestemc_storage_to_litter(:)   ! live stem C storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_livestemc_xfer_to_litter(:)      ! live stem C transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_deadstemc_to_prod10c(:)          ! dead stem C harvest to 10-year product pool (gC/m2/s)
   real, pointer :: hrv_deadstemc_to_prod100c(:)         ! dead stem C harvest to 100-year product pool (gC/m2/s)
   real, pointer :: hrv_deadstemc_storage_to_litter(:)   ! dead stem C storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_deadstemc_xfer_to_litter(:)      ! dead stem C transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_livecrootc_to_litter(:)          ! live coarse root C harvest mortality (gC/m2/s)
   real, pointer :: hrv_livecrootc_storage_to_litter(:)  ! live coarse root C storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_livecrootc_xfer_to_litter(:)     ! live coarse root C transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_deadcrootc_to_litter(:)          ! dead coarse root C harvest mortality (gC/m2/s)
   real, pointer :: hrv_deadcrootc_storage_to_litter(:)  ! dead coarse root C storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_deadcrootc_xfer_to_litter(:)     ! dead coarse root C transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_gresp_storage_to_litter(:)       ! growth respiration storage harvest mortality (gC/m2/s)
   real, pointer :: hrv_gresp_xfer_to_litter(:)          ! growth respiration transfer harvest mortality (gC/m2/s)
   real, pointer :: hrv_xsmrpool_to_atm(:)               ! excess MR pool harvest mortality (gC/m2/s)
   ! PFT-level fire fluxes
   real, pointer :: m_leafc_to_fire(:)                   ! leaf C fire loss (gC/m2/s)
   real, pointer :: m_leafc_storage_to_fire(:)           ! leaf C storage fire loss (gC/m2/s)
   real, pointer :: m_leafc_xfer_to_fire(:)              ! leaf C transfer fire loss (gC/m2/s)
   real, pointer :: m_frootc_to_fire(:)                  ! fine root C fire loss (gC/m2/s)
   real, pointer :: m_frootc_storage_to_fire(:)          ! fine root C storage fire loss (gC/m2/s)
   real, pointer :: m_frootc_xfer_to_fire(:)             ! fine root C transfer fire loss (gC/m2/s)
   real, pointer :: m_livestemc_to_fire(:)               ! live stem C fire loss (gC/m2/s)
   real, pointer :: m_livestemc_storage_to_fire(:)       ! live stem C storage fire loss (gC/m2/s)
   real, pointer :: m_livestemc_xfer_to_fire(:)          ! live stem C transfer fire loss (gC/m2/s)
   real, pointer :: m_deadstemc_to_fire(:)               ! dead stem C fire loss (gC/m2/s)
   real, pointer :: m_deadstemc_to_litter_fire(:)        ! dead stem C fire mortality to litter (gC/m2/s)
   real, pointer :: m_deadstemc_storage_to_fire(:)       ! dead stem C storage fire loss (gC/m2/s)
   real, pointer :: m_deadstemc_xfer_to_fire(:)          ! dead stem C transfer fire loss (gC/m2/s)
   real, pointer :: m_livecrootc_to_fire(:)              ! live coarse root C fire loss (gC/m2/s)
   real, pointer :: m_livecrootc_storage_to_fire(:)      ! live coarse root C storage fire loss (gC/m2/s)
   real, pointer :: m_livecrootc_xfer_to_fire(:)         ! live coarse root C transfer fire loss (gC/m2/s)
   real, pointer :: m_deadcrootc_to_fire(:)              ! dead coarse root C fire loss (gC/m2/s)
   real, pointer :: m_deadcrootc_to_litter_fire(:)       ! dead coarse root C fire mortality to litter (gC/m2/s)
   real, pointer :: m_deadcrootc_storage_to_fire(:)      ! dead coarse root C storage fire loss (gC/m2/s)
   real, pointer :: m_deadcrootc_xfer_to_fire(:)         ! dead coarse root C transfer fire loss (gC/m2/s)
   real, pointer :: m_gresp_storage_to_fire(:)           ! growth respiration storage fire loss (gC/m2/s)
   real, pointer :: m_gresp_xfer_to_fire(:)              ! growth respiration transfer fire loss (gC/m2/s)
   ! phenology fluxes from transfer pools                     
   real, pointer :: leafc_xfer_to_leafc(:)               ! leaf C growth from storage (gC/m2/s)
   real, pointer :: frootc_xfer_to_frootc(:)             ! fine root C growth from storage (gC/m2/s)
   real, pointer :: livestemc_xfer_to_livestemc(:)       ! live stem C growth from storage (gC/m2/s)
   real, pointer :: deadstemc_xfer_to_deadstemc(:)       ! dead stem C growth from storage (gC/m2/s)
   real, pointer :: livecrootc_xfer_to_livecrootc(:)     ! live coarse root C growth from storage (gC/m2/s)
   real, pointer :: deadcrootc_xfer_to_deadcrootc(:)     ! dead coarse root C growth from storage (gC/m2/s)
   ! leaf and fine root litterfall                           
   real, pointer :: leafc_to_litter(:)                   ! leaf C litterfall (gC/m2/s)
   real, pointer :: frootc_to_litter(:)                  ! fine root C litterfall (gC/m2/s)
   ! maintenance respiration fluxes                          
   real, pointer :: leaf_mr(:)                           ! leaf maintenance respiration (gC/m2/s)
   real, pointer :: froot_mr(:)                          ! fine root maintenance respiration (gC/m2/s)
   real, pointer :: livestem_mr(:)                       ! live stem maintenance respiration (gC/m2/s)
   real, pointer :: livecroot_mr(:)                      ! live coarse root maintenance respiration (gC/m2/s)
   real, pointer :: leaf_curmr(:)                        ! leaf maintenance respiration from current GPP (gC/m2/s)
   real, pointer :: froot_curmr(:)                       ! fine root maintenance respiration from current GPP (gC/m2/s)
   real, pointer :: livestem_curmr(:)                    ! live stem maintenance respiration from current GPP (gC/m2/s)
   real, pointer :: livecroot_curmr(:)                   ! live coarse root maintenance respiration from current GPP (gC/m2/s)
   real, pointer :: leaf_xsmr(:)                         ! leaf maintenance respiration from storage (gC/m2/s)
   real, pointer :: froot_xsmr(:)                        ! fine root maintenance respiration from storage (gC/m2/s)
   real, pointer :: livestem_xsmr(:)                     ! live stem maintenance respiration from storage (gC/m2/s)
   real, pointer :: livecroot_xsmr(:)                    ! live coarse root maintenance respiration from storage (gC/m2/s)
   ! photosynthesis fluxes                                   
   real, pointer :: psnsun_to_cpool(:)                   ! C fixation from sunlit canopy (gC/m2/s)
   real, pointer :: psnshade_to_cpool(:)                 ! C fixation from shaded canopy (gC/m2/s)
   ! allocation fluxes, from current GPP                     
   real, pointer :: cpool_to_xsmrpool(:)                 ! allocation to maintenance respiration storage pool (gC/m2/s)
   real, pointer :: cpool_to_leafc(:)                    ! allocation to leaf C (gC/m2/s)
   real, pointer :: cpool_to_leafc_storage(:)            ! allocation to leaf C storage (gC/m2/s)
   real, pointer :: cpool_to_frootc(:)                   ! allocation to fine root C (gC/m2/s)
   real, pointer :: cpool_to_frootc_storage(:)           ! allocation to fine root C storage (gC/m2/s)
   real, pointer :: cpool_to_livestemc(:)                ! allocation to live stem C (gC/m2/s)
   real, pointer :: cpool_to_livestemc_storage(:)        ! allocation to live stem C storage (gC/m2/s)
   real, pointer :: cpool_to_deadstemc(:)                ! allocation to dead stem C (gC/m2/s)
   real, pointer :: cpool_to_deadstemc_storage(:)        ! allocation to dead stem C storage (gC/m2/s)
   real, pointer :: cpool_to_livecrootc(:)               ! allocation to live coarse root C (gC/m2/s)
   real, pointer :: cpool_to_livecrootc_storage(:)       ! allocation to live coarse root C storage (gC/m2/s)
   real, pointer :: cpool_to_deadcrootc(:)               ! allocation to dead coarse root C (gC/m2/s)
   real, pointer :: cpool_to_deadcrootc_storage(:)       ! allocation to dead coarse root C storage (gC/m2/s)
   real, pointer :: cpool_to_gresp_storage(:)            ! allocation to growth respiration storage (gC/m2/s)
   ! growth respiration fluxes                               
   real, pointer :: cpool_leaf_gr(:)                     ! leaf growth respiration (gC/m2/s)
   real, pointer :: cpool_leaf_storage_gr(:)             ! leaf growth respiration to storage (gC/m2/s)
   real, pointer :: transfer_leaf_gr(:)                  ! leaf growth respiration from storage (gC/m2/s)
   real, pointer :: cpool_froot_gr(:)                    ! fine root growth respiration (gC/m2/s)
   real, pointer :: cpool_froot_storage_gr(:)            ! fine root  growth respiration to storage (gC/m2/s)
   real, pointer :: transfer_froot_gr(:)                 ! fine root  growth respiration from storage (gC/m2/s)
   real, pointer :: cpool_livestem_gr(:)                 ! live stem growth respiration (gC/m2/s)
   real, pointer :: cpool_livestem_storage_gr(:)         ! live stem growth respiration to storage (gC/m2/s)
   real, pointer :: transfer_livestem_gr(:)              ! live stem growth respiration from storage (gC/m2/s)
   real, pointer :: cpool_deadstem_gr(:)                 ! dead stem growth respiration (gC/m2/s)
   real, pointer :: cpool_deadstem_storage_gr(:)         ! dead stem growth respiration to storage (gC/m2/s)
   real, pointer :: transfer_deadstem_gr(:)              ! dead stem growth respiration from storage (gC/m2/s)
   real, pointer :: cpool_livecroot_gr(:)                ! live coarse root growth respiration (gC/m2/s)
   real, pointer :: cpool_livecroot_storage_gr(:)        ! live coarse root growth respiration to storage (gC/m2/s)
   real, pointer :: transfer_livecroot_gr(:)             ! live coarse root growth respiration from storage (gC/m2/s)
   real, pointer :: cpool_deadcroot_gr(:)                ! dead coarse root growth respiration (gC/m2/s)
   real, pointer :: cpool_deadcroot_storage_gr(:)        ! dead coarse root growth respiration to storage (gC/m2/s)
   real, pointer :: transfer_deadcroot_gr(:)             ! dead coarse root growth respiration from storage (gC/m2/s)
   ! annual turnover of storage to transfer pools            
   real, pointer :: leafc_storage_to_xfer(:)             ! leaf C shift storage to transfer (gC/m2/s)
   real, pointer :: frootc_storage_to_xfer(:)            ! fine root C shift storage to transfer (gC/m2/s)
   real, pointer :: livestemc_storage_to_xfer(:)         ! live stem C shift storage to transfer (gC/m2/s)
   real, pointer :: deadstemc_storage_to_xfer(:)         ! dead stem C shift storage to transfer (gC/m2/s)
   real, pointer :: livecrootc_storage_to_xfer(:)        ! live coarse root C shift storage to transfer (gC/m2/s)
   real, pointer :: deadcrootc_storage_to_xfer(:)        ! dead coarse root C shift storage to transfer (gC/m2/s)
   real, pointer :: gresp_storage_to_xfer(:)             ! growth respiration shift storage to transfer (gC/m2/s)
   ! turnover of livewood to deadwood
   real, pointer :: livestemc_to_deadstemc(:)            ! live stem C turnover (gC/m2/s)
   real, pointer :: livecrootc_to_deadcrootc(:)          ! live coarse root C turnover (gC/m2/s)
   ! summary (diagnostic) flux variables, not involved in mass balance
   real, pointer :: gpp(:)            ! (gC/m2/s) gross primary production 
   real, pointer :: mr(:)             ! (gC/m2/s) maintenance respiration
   real, pointer :: current_gr(:)     ! (gC/m2/s) growth resp for new growth displayed in this timestep
   real, pointer :: transfer_gr(:)    ! (gC/m2/s) growth resp for transfer growth displayed in this timestep
   real, pointer :: storage_gr(:)     ! (gC/m2/s) growth resp for growth sent to storage for later display
   real, pointer :: gr(:)             ! (gC/m2/s) total growth respiration
   real, pointer :: ar(:)             ! (gC/m2/s) autotrophic respiration (MR + GR)
   real, pointer :: rr(:)             ! (gC/m2/s) root respiration (fine root MR + total root GR)
   real, pointer :: npp(:)            ! (gC/m2/s) net primary production
   real, pointer :: agnpp(:)          ! (gC/m2/s) aboveground NPP
   real, pointer :: bgnpp(:)          ! (gC/m2/s) belowground NPP
   real, pointer :: litfall(:)        ! (gC/m2/s) litterfall (leaves and fine roots)
   real, pointer :: vegfire(:)        ! (gC/m2/s) pft-level fire loss (obsolete, mark for removal)
   real, pointer :: wood_harvestc(:)  ! (gC/m2/s) pft-level wood harvest (to product pools)
   real, pointer :: pft_cinputs(:)    ! (gC/m2/s) pft-level carbon inputs (for balance checking)
   real, pointer :: pft_coutputs(:)   ! (gC/m2/s) pft-level carbon outputs (for balance checking)
   ! new variables for fire code
   real, pointer :: pft_fire_closs(:) ! (gC/m2/s) total pft-level fire C loss 
end type pft_cflux_type

!----------------------------------------------------
! pft nitrogen flux variables structure
!----------------------------------------------------
type, public :: pft_nflux_type
   ! new variables for CN code
   ! gap mortality fluxes
   real, pointer :: m_leafn_to_litter(:)                ! leaf N mortality (gN/m2/s)
   real, pointer :: m_frootn_to_litter(:)               ! fine root N mortality (gN/m2/s)
   real, pointer :: m_leafn_storage_to_litter(:)        ! leaf N storage mortality (gN/m2/s)
   real, pointer :: m_frootn_storage_to_litter(:)       ! fine root N storage mortality (gN/m2/s)
   real, pointer :: m_livestemn_storage_to_litter(:)    ! live stem N storage mortality (gN/m2/s)
   real, pointer :: m_deadstemn_storage_to_litter(:)    ! dead stem N storage mortality (gN/m2/s)
   real, pointer :: m_livecrootn_storage_to_litter(:)   ! live coarse root N storage mortality (gN/m2/s)
   real, pointer :: m_deadcrootn_storage_to_litter(:)   ! dead coarse root N storage mortality (gN/m2/s)
   real, pointer :: m_leafn_xfer_to_litter(:)           ! leaf N transfer mortality (gN/m2/s)
   real, pointer :: m_frootn_xfer_to_litter(:)          ! fine root N transfer mortality (gN/m2/s)
   real, pointer :: m_livestemn_xfer_to_litter(:)       ! live stem N transfer mortality (gN/m2/s)
   real, pointer :: m_deadstemn_xfer_to_litter(:)       ! dead stem N transfer mortality (gN/m2/s)
   real, pointer :: m_livecrootn_xfer_to_litter(:)      ! live coarse root N transfer mortality (gN/m2/s)
   real, pointer :: m_deadcrootn_xfer_to_litter(:)      ! dead coarse root N transfer mortality (gN/m2/s)
   real, pointer :: m_livestemn_to_litter(:)            ! live stem N mortality (gN/m2/s)
   real, pointer :: m_deadstemn_to_litter(:)            ! dead stem N mortality (gN/m2/s)
   real, pointer :: m_livecrootn_to_litter(:)           ! live coarse root N mortality (gN/m2/s)
   real, pointer :: m_deadcrootn_to_litter(:)           ! dead coarse root N mortality (gN/m2/s)
   real, pointer :: m_retransn_to_litter(:)             ! retranslocated N pool mortality (gN/m2/s)
   ! harvest mortality fluxes
   real, pointer :: hrv_leafn_to_litter(:)                ! leaf N harvest mortality (gN/m2/s)
   real, pointer :: hrv_frootn_to_litter(:)               ! fine root N harvest mortality (gN/m2/s)
   real, pointer :: hrv_leafn_storage_to_litter(:)        ! leaf N storage harvest mortality (gN/m2/s)
   real, pointer :: hrv_frootn_storage_to_litter(:)       ! fine root N storage harvest mortality (gN/m2/s)
   real, pointer :: hrv_livestemn_storage_to_litter(:)    ! live stem N storage harvest mortality (gN/m2/s)
   real, pointer :: hrv_deadstemn_storage_to_litter(:)    ! dead stem N storage harvest mortality (gN/m2/s)
   real, pointer :: hrv_livecrootn_storage_to_litter(:)   ! live coarse root N storage harvest mortality (gN/m2/s)
   real, pointer :: hrv_deadcrootn_storage_to_litter(:)   ! dead coarse root N storage harvest mortality (gN/m2/s)
   real, pointer :: hrv_leafn_xfer_to_litter(:)           ! leaf N transfer harvest mortality (gN/m2/s)
   real, pointer :: hrv_frootn_xfer_to_litter(:)          ! fine root N transfer harvest mortality (gN/m2/s)
   real, pointer :: hrv_livestemn_xfer_to_litter(:)       ! live stem N transfer harvest mortality (gN/m2/s)
   real, pointer :: hrv_deadstemn_xfer_to_litter(:)       ! dead stem N transfer harvest mortality (gN/m2/s)
   real, pointer :: hrv_livecrootn_xfer_to_litter(:)      ! live coarse root N transfer harvest mortality (gN/m2/s)
   real, pointer :: hrv_deadcrootn_xfer_to_litter(:)      ! dead coarse root N transfer harvest mortality (gN/m2/s)
   real, pointer :: hrv_livestemn_to_litter(:)            ! live stem N harvest mortality (gN/m2/s)
   real, pointer :: hrv_deadstemn_to_prod10n(:)           ! dead stem N harvest to 10-year product pool (gN/m2/s)
   real, pointer :: hrv_deadstemn_to_prod100n(:)          ! dead stem N harvest to 100-year product pool (gN/m2/s)
   real, pointer :: hrv_livecrootn_to_litter(:)           ! live coarse root N harvest mortality (gN/m2/s)
   real, pointer :: hrv_deadcrootn_to_litter(:)           ! dead coarse root N harvest mortality (gN/m2/s)
   real, pointer :: hrv_retransn_to_litter(:)             ! retranslocated N pool harvest mortality (gN/m2/s)
   ! fire mortality fluxes
   real, pointer :: m_leafn_to_fire(:)                  ! leaf N fire loss (gN/m2/s)
   real, pointer :: m_leafn_storage_to_fire(:)          ! leaf N storage fire loss (gN/m2/s)
   real, pointer :: m_leafn_xfer_to_fire(:)             ! leaf N transfer fire loss (gN/m2/s)
   real, pointer :: m_frootn_to_fire(:)                 ! fine root N fire loss (gN/m2/s)
   real, pointer :: m_frootn_storage_to_fire(:)         ! fine root N storage fire loss (gN/m2/s)
   real, pointer :: m_frootn_xfer_to_fire(:)            ! fine root N transfer fire loss (gN/m2/s)
   real, pointer :: m_livestemn_to_fire(:)              ! live stem N fire loss (gN/m2/s)
   real, pointer :: m_livestemn_storage_to_fire(:)      ! live stem N storage fire loss (gN/m2/s)
   real, pointer :: m_livestemn_xfer_to_fire(:)         ! live stem N transfer fire loss (gN/m2/s)
   real, pointer :: m_deadstemn_to_fire(:)              ! dead stem N fire loss (gN/m2/s)
   real, pointer :: m_deadstemn_to_litter_fire(:)       ! dead stem N fire mortality to litter (gN/m2/s)
   real, pointer :: m_deadstemn_storage_to_fire(:)      ! dead stem N storage fire loss (gN/m2/s)
   real, pointer :: m_deadstemn_xfer_to_fire(:)         ! dead stem N transfer fire loss (gN/m2/s)
   real, pointer :: m_livecrootn_to_fire(:)             ! live coarse root N fire loss (gN/m2/s)
   real, pointer :: m_livecrootn_storage_to_fire(:)     ! live coarse root N storage fire loss (gN/m2/s)
   real, pointer :: m_livecrootn_xfer_to_fire(:)        ! live coarse root N transfer fire loss (gN/m2/s)
   real, pointer :: m_deadcrootn_to_fire(:)             ! dead coarse root N fire loss (gN/m2/s)
   real, pointer :: m_deadcrootn_to_litter_fire(:)      ! dead coarse root N fire mortality to litter (gN/m2/s)
   real, pointer :: m_deadcrootn_storage_to_fire(:)     ! dead coarse root N storage fire loss (gN/m2/s)
   real, pointer :: m_deadcrootn_xfer_to_fire(:)        ! dead coarse root N transfer fire loss (gN/m2/s)
   real, pointer :: m_retransn_to_fire(:)               ! retranslocated N pool fire loss (gN/m2/s)
   ! phenology fluxes from transfer pool                     
   real, pointer :: leafn_xfer_to_leafn(:)              ! leaf N growth from storage (gN/m2/s)
   real, pointer :: frootn_xfer_to_frootn(:)            ! fine root N growth from storage (gN/m2/s)
   real, pointer :: livestemn_xfer_to_livestemn(:)      ! live stem N growth from storage (gN/m2/s)
   real, pointer :: deadstemn_xfer_to_deadstemn(:)      ! dead stem N growth from storage (gN/m2/s)
   real, pointer :: livecrootn_xfer_to_livecrootn(:)    ! live coarse root N growth from storage (gN/m2/s)
   real, pointer :: deadcrootn_xfer_to_deadcrootn(:)    ! dead coarse root N growth from storage (gN/m2/s)
   ! litterfall fluxes
   real, pointer :: leafn_to_litter(:)                  ! leaf N litterfall (gN/m2/s)
   real, pointer :: leafn_to_retransn(:)                ! leaf N to retranslocated N pool (gN/m2/s)
   real, pointer :: frootn_to_litter(:)                 ! fine root N litterfall (gN/m2/s)
   ! allocation fluxes
   real, pointer :: retransn_to_npool(:)                ! deployment of retranslocated N (gN/m2/s)       
   real, pointer :: sminn_to_npool(:)                   ! deployment of soil mineral N uptake (gN/m2/s)
   real, pointer :: npool_to_leafn(:)                   ! allocation to leaf N (gN/m2/s)
   real, pointer :: npool_to_leafn_storage(:)           ! allocation to leaf N storage (gN/m2/s)
   real, pointer :: npool_to_frootn(:)                  ! allocation to fine root N (gN/m2/s)
   real, pointer :: npool_to_frootn_storage(:)          ! allocation to fine root N storage (gN/m2/s)
   real, pointer :: npool_to_livestemn(:)               ! allocation to live stem N (gN/m2/s)
   real, pointer :: npool_to_livestemn_storage(:)       ! allocation to live stem N storage (gN/m2/s)
   real, pointer :: npool_to_deadstemn(:)               ! allocation to dead stem N (gN/m2/s)
   real, pointer :: npool_to_deadstemn_storage(:)       ! allocation to dead stem N storage (gN/m2/s)
   real, pointer :: npool_to_livecrootn(:)              ! allocation to live coarse root N (gN/m2/s)
   real, pointer :: npool_to_livecrootn_storage(:)      ! allocation to live coarse root N storage (gN/m2/s)
   real, pointer :: npool_to_deadcrootn(:)              ! allocation to dead coarse root N (gN/m2/s)
   real, pointer :: npool_to_deadcrootn_storage(:)      ! allocation to dead coarse root N storage (gN/m2/s)
   ! annual turnover of storage to transfer pools           
   real, pointer :: leafn_storage_to_xfer(:)            ! leaf N shift storage to transfer (gN/m2/s)
   real, pointer :: frootn_storage_to_xfer(:)           ! fine root N shift storage to transfer (gN/m2/s)
   real, pointer :: livestemn_storage_to_xfer(:)        ! live stem N shift storage to transfer (gN/m2/s)
   real, pointer :: deadstemn_storage_to_xfer(:)        ! dead stem N shift storage to transfer (gN/m2/s)
   real, pointer :: livecrootn_storage_to_xfer(:)       ! live coarse root N shift storage to transfer (gN/m2/s)
   real, pointer :: deadcrootn_storage_to_xfer(:)       ! dead coarse root N shift storage to transfer (gN/m2/s)
   ! turnover of livewood to deadwood, with retranslocation 
   real, pointer :: livestemn_to_deadstemn(:)           ! live stem N turnover (gN/m2/s)
   real, pointer :: livestemn_to_retransn(:)            ! live stem N to retranslocated N pool (gN/m2/s)
   real, pointer :: livecrootn_to_deadcrootn(:)         ! live coarse root N turnover (gN/m2/s)
   real, pointer :: livecrootn_to_retransn(:)           ! live coarse root N to retranslocated N pool (gN/m2/s)
   ! summary (diagnostic) flux variables, not involved in mass balance
   real, pointer :: ndeploy(:)                          ! total N deployed to growth and storage (gN/m2/s)
   real, pointer :: pft_ninputs(:)                      ! total N inputs to pft-level (gN/m2/s)
   real, pointer :: pft_noutputs(:)                     ! total N outputs from pft-level (gN/m2/s)
   real, pointer :: wood_harvestn(:)                    ! total N losses to wood product pools (gN/m2/s)
   ! new variables for fire code 
   real, pointer :: pft_fire_nloss(:)                   ! total pft-level fire N loss (gN/m2/s) 
end type pft_nflux_type

!----------------------------------------------------
! End definition of structures defined at the pft_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the column_type level
!----------------------------------------------------
! column physical state variables structure
!----------------------------------------------------
type, public :: column_pstate_type
   real, pointer :: snowdp(:)             !snow height (m)
   real, pointer :: dz(:,:)               !layer thickness (m)  (-nlevsno+1:nlevgrnd) 
   real, pointer :: wf(:)                 !soil water as frac. of whc for top 0.5 m
   ! new variables for CN code
   real, pointer :: psisat(:,:)        !soil water potential at saturation for CN code (MPa)
   real, pointer :: psiwilt(:)         !root-zone soil water potential at wilting point (MPa)
   real, pointer :: soilpsi(:,:)         !soil water potential in each soil layer (MPa)
   real, pointer :: fpi(:)           !fraction of potential immobilization (no units)
   real, pointer :: fpg(:)           !fraction of potential gpp (no units)
   real, pointer :: annsum_counter(:) !seconds since last annual accumulator turnover
   real, pointer :: cannsum_npp(:)    !annual sum of NPP, averaged from pft-level (gC/m2/yr)
   real, pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
   ! new variables for fire code
   real, pointer :: me(:)                 !moisture of extinction (proportion) 
   real, pointer :: fire_prob(:)          !daily fire probability (0-1) 
   real, pointer :: mean_fire_prob(:)     !e-folding mean of daily fire probability (0-1) 
   real, pointer :: fireseasonl(:)        !annual fire season length (days, <= 365) 
   real, pointer :: farea_burned(:)       !timestep fractional area burned (proportion) 
   real, pointer :: ann_farea_burned(:)   !annual total fractional area burned (proportion)
end type column_pstate_type

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
type, public :: column_estate_type
   real, pointer :: t_grnd(:)             !ground temperature (Kelvin)
   real, pointer :: t_soisno(:,:)         !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd) 
end type column_estate_type

!----------------------------------------------------
! column water state variables structure
!----------------------------------------------------
type, public :: column_wstate_type
   real, pointer :: h2osoi_liq(:)         !column liquid water (kg/m2) (new)   
end type column_wstate_type

!----------------------------------------------------
! column carbon state variables structure
!----------------------------------------------------
type, public :: column_cstate_type
   type(pft_cstate_type):: pcs_a              !pft-level carbon state variables averaged to the column
   ! BGC variables
   real, pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real, pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real, pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real, pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real, pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real, pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real, pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real, pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slowest pool)
   real, pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   ! pools for dynamic landcover
   real, pointer :: seedc(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real, pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
   real, pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
   real, pointer :: totprodc(:)           ! (gC/m2) total wood product C
   ! summary (diagnostic) state variables, not involved in mass balance
   real, pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real, pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon
   real, pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real, pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   
end type column_cstate_type

!----------------------------------------------------
! column nitrogen state variables structure
!----------------------------------------------------
type, public :: column_nstate_type
   type(pft_nstate_type):: pns_a              !pft-level nitrogen state variables averaged to the column
   ! BGC variables
   real, pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real, pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real, pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real, pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real, pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real, pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real, pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real, pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real, pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real, pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   ! wood product pools, for dynamic landcover
   real, pointer :: seedn(:)              ! (gN/m2) column-level pool for seeding new PFTs
   real, pointer :: prod10n(:)            ! (gN/m2) wood product N pool, 10-year lifespan
   real, pointer :: prod100n(:)           ! (gN/m2) wood product N pool, 100-year lifespan
   real, pointer :: totprodn(:)           ! (gN/m2) total wood product N
   ! summary (diagnostic) state variables, not involved in mass balance
   real, pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
   real, pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen
   real, pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg 
   real, pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
end type column_nstate_type

!----------------------------------------------------
! column water flux variables structure
!----------------------------------------------------
type, public :: column_wflux_type
   real, pointer :: qflx_drain(:) 	  ! sub-surface runoff (mm H2O /s)
end type column_wflux_type

!----------------------------------------------------
! column carbon flux variables structure
!----------------------------------------------------
type, public :: column_cflux_type
   type(pft_cflux_type):: pcf_a                           !pft-level carbon flux variables averaged to the column
   ! new variables for CN code
   ! column-level gap mortality fluxes
   real, pointer :: m_leafc_to_litr1c(:)              ! leaf C mortality to litter 1 C (gC/m2/s) 
   real, pointer :: m_leafc_to_litr2c(:)              ! leaf C mortality to litter 2 C (gC/m2/s)
   real, pointer :: m_leafc_to_litr3c(:)              ! leaf C mortality to litter 3 C (gC/m2/s)
   real, pointer :: m_frootc_to_litr1c(:)             ! fine root C mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_frootc_to_litr2c(:)             ! fine root C mortality to litter 2 C (gC/m2/s)
   real, pointer :: m_frootc_to_litr3c(:)             ! fine root C mortality to litter 3 C (gC/m2/s)
   real, pointer :: m_livestemc_to_cwdc(:)            ! live stem C mortality to coarse woody debris C (gC/m2/s)
   real, pointer :: m_deadstemc_to_cwdc(:)            ! dead stem C mortality to coarse woody debris C (gC/m2/s)
   real, pointer :: m_livecrootc_to_cwdc(:)           ! live coarse root C mortality to coarse woody debris C (gC/m2/s)
   real, pointer :: m_deadcrootc_to_cwdc(:)           ! dead coarse root C mortality to coarse woody debris C (gC/m2/s)
   real, pointer :: m_leafc_storage_to_litr1c(:)      ! leaf C storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_frootc_storage_to_litr1c(:)     ! fine root C storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_livestemc_storage_to_litr1c(:)  ! live stem C storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_deadstemc_storage_to_litr1c(:)  ! dead stem C storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_livecrootc_storage_to_litr1c(:) ! live coarse root C storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_deadcrootc_storage_to_litr1c(:) ! dead coarse root C storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_gresp_storage_to_litr1c(:)      ! growth respiration storage mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_leafc_xfer_to_litr1c(:)         ! leaf C transfer mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_frootc_xfer_to_litr1c(:)        ! fine root C transfer mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_livestemc_xfer_to_litr1c(:)     ! live stem C transfer mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_deadstemc_xfer_to_litr1c(:)     ! dead stem C transfer mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_livecrootc_xfer_to_litr1c(:)    ! live coarse root C transfer mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_deadcrootc_xfer_to_litr1c(:)    ! dead coarse root C transfer mortality to litter 1 C (gC/m2/s)
   real, pointer :: m_gresp_xfer_to_litr1c(:)         ! growth respiration transfer mortality to litter 1 C (gC/m2/s)
   ! column-level harvest mortality fluxes
   real, pointer :: hrv_leafc_to_litr1c(:)               ! leaf C harvest mortality to litter 1 C (gC/m2/s)                         
   real, pointer :: hrv_leafc_to_litr2c(:)               ! leaf C harvest mortality to litter 2 C (gC/m2/s)                        
   real, pointer :: hrv_leafc_to_litr3c(:)               ! leaf C harvest mortality to litter 3 C (gC/m2/s)                        
   real, pointer :: hrv_frootc_to_litr1c(:)              ! fine root C harvest mortality to litter 1 C (gC/m2/s)                   
   real, pointer :: hrv_frootc_to_litr2c(:)              ! fine root C harvest mortality to litter 2 C (gC/m2/s)                   
   real, pointer :: hrv_frootc_to_litr3c(:)              ! fine root C harvest mortality to litter 3 C (gC/m2/s)                   
   real, pointer :: hrv_livestemc_to_cwdc(:)             ! live stem C harvest mortality to coarse woody debris C (gC/m2/s)        
   real, pointer :: hrv_deadstemc_to_prod10c(:)          ! dead stem C harvest mortality to 10-year product pool (gC/m2/s)        
   real, pointer :: hrv_deadstemc_to_prod100c(:)         ! dead stem C harvest mortality to 100-year product pool (gC/m2/s)        
   real, pointer :: hrv_livecrootc_to_cwdc(:)            ! live coarse root C harvest mortality to coarse woody debris C (gC/m2/s) 
   real, pointer :: hrv_deadcrootc_to_cwdc(:)            ! dead coarse root C harvest mortality to coarse woody debris C (gC/m2/s) 
   real, pointer :: hrv_leafc_storage_to_litr1c(:)       ! leaf C storage harvest mortality to litter 1 C (gC/m2/s)                
   real, pointer :: hrv_frootc_storage_to_litr1c(:)      ! fine root C storage harvest mortality to litter 1 C (gC/m2/s)           
   real, pointer :: hrv_livestemc_storage_to_litr1c(:)   ! live stem C storage harvest mortality to litter 1 C (gC/m2/s)           
   real, pointer :: hrv_deadstemc_storage_to_litr1c(:)   ! dead stem C storage harvest mortality to litter 1 C (gC/m2/s)           
   real, pointer :: hrv_livecrootc_storage_to_litr1c(:)  ! live coarse root C storage harvest mortality to litter 1 C (gC/m2/s)    
   real, pointer :: hrv_deadcrootc_storage_to_litr1c(:)  ! dead coarse root C storage harvest mortality to litter 1 C (gC/m2/s)    
   real, pointer :: hrv_gresp_storage_to_litr1c(:)       ! growth respiration storage harvest mortality to litter 1 C (gC/m2/s)    
   real, pointer :: hrv_leafc_xfer_to_litr1c(:)          ! leaf C transfer harvest mortality to litter 1 C (gC/m2/s)               
   real, pointer :: hrv_frootc_xfer_to_litr1c(:)         ! fine root C transfer harvest mortality to litter 1 C (gC/m2/s)          
   real, pointer :: hrv_livestemc_xfer_to_litr1c(:)      ! live stem C transfer harvest mortality to litter 1 C (gC/m2/s)          
   real, pointer :: hrv_deadstemc_xfer_to_litr1c(:)      ! dead stem C transfer harvest mortality to litter 1 C (gC/m2/s)          
   real, pointer :: hrv_livecrootc_xfer_to_litr1c(:)     ! live coarse root C transfer harvest mortality to litter 1 C (gC/m2/s)   
   real, pointer :: hrv_deadcrootc_xfer_to_litr1c(:)     ! dead coarse root C transfer harvest mortality to litter 1 C (gC/m2/s)   
   real, pointer :: hrv_gresp_xfer_to_litr1c(:)          ! growth respiration transfer harvest mortality to litter 1 C (gC/m2/s)   
   ! column-level fire fluxes
   real, pointer :: m_deadstemc_to_cwdc_fire(:)       ! dead stem C to coarse woody debris C by fire (gC/m2/s)
   real, pointer :: m_deadcrootc_to_cwdc_fire(:)      ! dead coarse root C to to woody debris C by fire (gC/m2/s)
   real, pointer :: m_litr1c_to_fire(:)               ! litter 1 C fire loss (gC/m2/s)
   real, pointer :: m_litr2c_to_fire(:)               ! litter 2 C fire loss (gC/m2/s)
   real, pointer :: m_litr3c_to_fire(:)               ! litter 3 C fire loss (gC/m2/s)
   real, pointer :: m_cwdc_to_fire(:)                 ! coarse woody debris C fire loss (gC/m2/s)
   ! litterfall fluxes
   real, pointer :: leafc_to_litr1c(:)                ! leaf C litterfall to litter 1 C (gC/m2/s)
   real, pointer :: leafc_to_litr2c(:)                ! leaf C litterfall to litter 2 C (gC/m2/s)
   real, pointer :: leafc_to_litr3c(:)                ! leaf C litterfall to litter 3 C (gC/m2/s)
   real, pointer :: frootc_to_litr1c(:)               ! fine root C litterfall to litter 1 C (gC/m2/s)
   real, pointer :: frootc_to_litr2c(:)               ! fine root C litterfall to litter 2 C (gC/m2/s)
   real, pointer :: frootc_to_litr3c(:)               ! fine root C litterfall to litter 3 C (gC/m2/s)
   ! decomposition fluxes
   real, pointer :: cwdc_to_litr2c(:)     ! decomp. of coarse woody debris C to litter 2 C (gC/m2/s)
   real, pointer :: cwdc_to_litr3c(:)     ! decomp. of coarse woody debris C to litter 3 C (gC/m2/s)
   real, pointer :: litr1_hr(:)           ! het. resp. from litter 1 C (gC/m2/s)
   real, pointer :: litr1c_to_soil1c(:)   ! decomp. of litter 1 C to SOM 1 C (gC/m2/s)
   real, pointer :: litr2_hr(:)           ! het. resp. from litter 2 C (gC/m2/s)
   real, pointer :: litr2c_to_soil2c(:)   ! decomp. of litter 2 C to SOM 2 C (gC/m2/s)
   real, pointer :: litr3_hr(:)           ! het. resp. from litter 3 C (gC/m2/s)
   real, pointer :: litr3c_to_soil3c(:)   ! decomp. of litter 3 C to SOM 3 C (gC/m2/s)
   real, pointer :: soil1_hr(:)           ! het. resp. from SOM 1 C (gC/m2/s)
   real, pointer :: soil1c_to_soil2c(:)   ! decomp. of SOM 1 C to SOM 2 C (gC/m2/s)
   real, pointer :: soil2_hr(:)           ! het. resp. from SOM 2 C (gC/m2/s)
   real, pointer :: soil2c_to_soil3c(:)   ! decomp. of SOM 2 C to SOM 3 C (gC/m2/s)
   real, pointer :: soil3_hr(:)           ! het. resp. from SOM 3 C (gC/m2/s)
   real, pointer :: soil3c_to_soil4c(:)   ! decomp. of SOM 3 C to SOM 4 C (gC/m2/s)
   real, pointer :: soil4_hr(:)           ! het. resp. from SOM 4 C (gC/m2/s)
   ! dynamic landcover fluxes
   real, pointer :: dwt_seedc_to_leaf(:)      ! (gC/m2/s) seed source to PFT-level
   real, pointer :: dwt_seedc_to_deadstem(:)  ! (gC/m2/s) seed source to PFT-level
   real, pointer :: dwt_conv_cflux(:)         ! (gC/m2/s) conversion C flux (immediate loss to atm)
   real, pointer :: dwt_prod10c_gain(:)       ! (gC/m2/s) addition to 10-yr wood product pool
   real, pointer :: dwt_prod100c_gain(:)      ! (gC/m2/s) addition to 100-yr wood product pool
   real, pointer :: dwt_frootc_to_litr1c(:)   ! (gC/m2/s) fine root to litter due to landcover change
   real, pointer :: dwt_frootc_to_litr2c(:)   ! (gC/m2/s) fine root to litter due to landcover change
   real, pointer :: dwt_frootc_to_litr3c(:)   ! (gC/m2/s) fine root to litter due to landcover change
   real, pointer :: dwt_livecrootc_to_cwdc(:) ! (gC/m2/s) live coarse root to CWD due to landcover change
   real, pointer :: dwt_deadcrootc_to_cwdc(:) ! (gC/m2/s) dead coarse root to CWD due to landcover change
   real, pointer :: dwt_closs(:)              ! (gC/m2/s) total carbon loss from product pools and conversion
   real, pointer :: landuseflux(:)            ! (gC/m2/s) dwt_closs+product_closs
   real, pointer :: landuptake(:)             ! (gC/m2/s) nee-landuseflux
   ! wood product pool loss fluxes
   real, pointer :: prod10c_loss(:)           ! (gC/m2/s) decomposition loss from 10-yr wood product pool
   real, pointer :: prod100c_loss(:)          ! (gC/m2/s) decomposition loss from 100-yr wood product pool
   real, pointer :: product_closs(:)          ! (gC/m2/s) total wood product carbon loss
   ! summary (diagnostic) flux variables, not involved in mass balance
   real, pointer :: lithr(:)         ! (gC/m2/s) litter heterotrophic respiration 
   real, pointer :: somhr(:)         ! (gC/m2/s) soil organic matter heterotrophic respiration
   real, pointer :: hr(:)            ! (gC/m2/s) total heterotrophic respiration
   real, pointer :: sr(:)            ! (gC/m2/s) total soil respiration (HR + root resp)
   real, pointer :: er(:)            ! (gC/m2/s) total ecosystem respiration, autotrophic + heterotrophic
   real, pointer :: litfire(:)       ! (gC/m2/s) litter fire losses
   real, pointer :: somfire(:)       ! (gC/m2/s) soil organic matter fire losses
   real, pointer :: totfire(:)       ! (gC/m2/s) total ecosystem fire losses
   real, pointer :: nep(:)           ! (gC/m2/s) net ecosystem production, excludes fire, landuse, and harvest flux, positive for sink
   real, pointer :: nbp(:)           ! (gC/m2/s) net biome production, includes fire, landuse, and harvest flux, positive for sink
   real, pointer :: nee(:)           ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source
   real, pointer :: col_cinputs(:)   ! (gC/m2/s) total column-level carbon inputs (for balance check)
   real, pointer :: col_coutputs(:)  ! (gC/m2/s) total column-level carbon outputs (for balance check) 

   ! new variables for fire
   real, pointer :: col_fire_closs(:) ! (gC/m2/s) total column-level fire C loss
end type column_cflux_type

!----------------------------------------------------
! column nitrogen flux variables structure
!----------------------------------------------------
type, public :: column_nflux_type
   type(pft_nflux_type):: pnf_a        !pft-level nitrogen flux variables averaged to the column
   ! new variables for CN code
   ! deposition fluxes
   real, pointer :: ndep_to_sminn(:)                   ! atmospheric N deposition to soil mineral N (gN/m2/s)
   real, pointer :: nfix_to_sminn(:)                   ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s) 
   ! column-level gap mortality fluxes
   real, pointer :: m_leafn_to_litr1n(:)               ! leaf N mortality to litter 1 N (gC/m2/s)
   real, pointer :: m_leafn_to_litr2n(:)               ! leaf N mortality to litter 2 N (gC/m2/s)
   real, pointer :: m_leafn_to_litr3n(:)               ! leaf N mortality to litter 3 N (gC/m2/s)
   real, pointer :: m_frootn_to_litr1n(:)              ! fine root N mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_frootn_to_litr2n(:)              ! fine root N mortality to litter 2 N (gN/m2/s)
   real, pointer :: m_frootn_to_litr3n(:)              ! fine root N mortality to litter 3 N (gN/m2/s)
   real, pointer :: m_livestemn_to_cwdn(:)             ! live stem N mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: m_deadstemn_to_cwdn(:)             ! dead stem N mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: m_livecrootn_to_cwdn(:)            ! live coarse root N mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: m_deadcrootn_to_cwdn(:)            ! dead coarse root N mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: m_retransn_to_litr1n(:)            ! retranslocated N pool mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_leafn_storage_to_litr1n(:)       ! leaf N storage mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_frootn_storage_to_litr1n(:)      ! fine root N storage mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_livestemn_storage_to_litr1n(:)   ! live stem N storage mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_deadstemn_storage_to_litr1n(:)   ! dead stem N storage mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_livecrootn_storage_to_litr1n(:)  ! live coarse root N storage mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_deadcrootn_storage_to_litr1n(:)  ! dead coarse root N storage mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_leafn_xfer_to_litr1n(:)          ! leaf N transfer mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_frootn_xfer_to_litr1n(:)         ! fine root N transfer mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_livestemn_xfer_to_litr1n(:)      ! live stem N transfer mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_deadstemn_xfer_to_litr1n(:)      ! dead stem N transfer mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_livecrootn_xfer_to_litr1n(:)     ! live coarse root N transfer mortality to litter 1 N (gN/m2/s)
   real, pointer :: m_deadcrootn_xfer_to_litr1n(:)     ! dead coarse root N transfer mortality to litter 1 N (gN/m2/s)
   ! column-level harvest fluxes
   real, pointer :: hrv_leafn_to_litr1n(:)               ! leaf N harvest mortality to litter 1 N (gC/m2/s)
   real, pointer :: hrv_leafn_to_litr2n(:)               ! leaf N harvest mortality to litter 2 N (gC/m2/s)
   real, pointer :: hrv_leafn_to_litr3n(:)               ! leaf N harvest mortality to litter 3 N (gC/m2/s)
   real, pointer :: hrv_frootn_to_litr1n(:)              ! fine root N harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_frootn_to_litr2n(:)              ! fine root N harvest mortality to litter 2 N (gN/m2/s)
   real, pointer :: hrv_frootn_to_litr3n(:)              ! fine root N harvest mortality to litter 3 N (gN/m2/s)
   real, pointer :: hrv_livestemn_to_cwdn(:)             ! live stem N harvest mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: hrv_deadstemn_to_prod10n(:)          ! dead stem N harvest mortality to 10-year product pool (gN/m2/s)
   real, pointer :: hrv_deadstemn_to_prod100n(:)         ! dead stem N harvest mortality to 100-year product pool (gN/m2/s)
   real, pointer :: hrv_livecrootn_to_cwdn(:)            ! live coarse root N harvest mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: hrv_deadcrootn_to_cwdn(:)            ! dead coarse root N harvest mortality to coarse woody debris N (gN/m2/s)
   real, pointer :: hrv_retransn_to_litr1n(:)            ! retranslocated N pool harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_leafn_storage_to_litr1n(:)       ! leaf N storage harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_frootn_storage_to_litr1n(:)      ! fine root N storage harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_livestemn_storage_to_litr1n(:)   ! live stem N storage harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_deadstemn_storage_to_litr1n(:)   ! dead stem N storage harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_livecrootn_storage_to_litr1n(:)  ! live coarse root N storage harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_deadcrootn_storage_to_litr1n(:)  ! dead coarse root N storage harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_leafn_xfer_to_litr1n(:)          ! leaf N transfer harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_frootn_xfer_to_litr1n(:)         ! fine root N transfer harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_livestemn_xfer_to_litr1n(:)      ! live stem N transfer harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_deadstemn_xfer_to_litr1n(:)      ! dead stem N transfer harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_livecrootn_xfer_to_litr1n(:)     ! live coarse root N transfer harvest mortality to litter 1 N (gN/m2/s)
   real, pointer :: hrv_deadcrootn_xfer_to_litr1n(:)     ! dead coarse root N transfer harvest mortality to litter 1 N (gN/m2/s)
   ! column-level fire fluxes
   real, pointer :: m_deadstemn_to_cwdn_fire(:)        ! dead stem N to coarse woody debris N by fire (gN/m2/s)
   real, pointer :: m_deadcrootn_to_cwdn_fire(:)       ! dead coarse root N to to woody debris N by fire (gN/m2/s)
   real, pointer :: m_litr1n_to_fire(:)                ! litter 1 N fire loss (gN/m2/s)
   real, pointer :: m_litr2n_to_fire(:)                ! litter 2 N fire loss (gN/m2/s)
   real, pointer :: m_litr3n_to_fire(:)                ! litter 3 N fire loss (gN/m2/s)
   real, pointer :: m_cwdn_to_fire(:)                  ! coarse woody debris N fire loss (gN/m2/s)
   ! litterfall fluxes
   real, pointer :: leafn_to_litr1n(:)       ! leaf N litterfall to litter 1 N (gN/m2/s)
   real, pointer :: leafn_to_litr2n(:)       ! leaf N litterfall to litter 2 N (gN/m2/s)
   real, pointer :: leafn_to_litr3n(:)       ! leaf N litterfall to litter 3 N (gN/m2/s)
   real, pointer :: frootn_to_litr1n(:)      ! fine root N litterfall to litter 1 N (gN/m2/s)
   real, pointer :: frootn_to_litr2n(:)      ! fine root N litterfall to litter 2 N (gN/m2/s)
   real, pointer :: frootn_to_litr3n(:)      ! fine root N litterfall to litter 3 N (gN/m2/s)
   ! decomposition fluxes
   real, pointer :: cwdn_to_litr2n(:)        ! decomp. of coarse woody debris N to litter 2 N (gN/m2/s)
   real, pointer :: cwdn_to_litr3n(:)        ! decomp. of coarse woody debris N to litter 3 N (gN/m2/s)
   real, pointer :: litr1n_to_soil1n(:)      ! decomp. of litter 1 N to SOM 1 N (gN/m2/s)
   real, pointer :: sminn_to_soil1n_l1(:)    ! mineral N flux for decomp. of litter 1 to SOM 1 (gN/m2/s)
   real, pointer :: litr2n_to_soil2n(:)      ! decomp. of litter 2 N to SOM 2 N (gN/m2/s)
   real, pointer :: sminn_to_soil2n_l2(:)    ! mineral N flux for decomp. of litter 2 to SOM 2 (gN/m2/s)
   real, pointer :: litr3n_to_soil3n(:)      ! decomp. of litter 3 N to SOM 3 N (gN/m2/s)
   real, pointer :: sminn_to_soil3n_l3(:)    ! mineral N flux for decomp. of litter 3 to SOM 3 (gN/m2/s)
   real, pointer :: soil1n_to_soil2n(:)      ! decomp. of SOM 1 N to SOM 2 N (gN/m2/s)
   real, pointer :: sminn_to_soil2n_s1(:)    ! mineral N flux for decomp. of SOM 1 to SOM 2 (gN/m2/s)
   real, pointer :: soil2n_to_soil3n(:)      ! decomp. of SOM 2 N to SOM 3 N (gN/m2/s)
   real, pointer :: sminn_to_soil3n_s2(:)    ! mineral N flux for decomp. of SOM 2 to SOM 3 (gN/m2/s)
   real, pointer :: soil3n_to_soil4n(:)      ! decomp. of SOM 3 N to SOM 4 N (gN/m2/s)
   real, pointer :: sminn_to_soil4n_s3(:)    ! mineral N flux for decomp. of SOM 3 to SOM 4 (gN/m2/s)
   real, pointer :: soil4n_to_sminn(:)       ! N mineralization for decomp. of SOM 4 (gN/m2/s)
   ! denitrification fluxes
   real, pointer :: sminn_to_denit_l1s1(:)   ! denitrification for decomp. of litter 1 to SOM 1 (gN/m2/s) 
   real, pointer :: sminn_to_denit_l2s2(:)   ! denitrification for decomp. of litter 2 to SOM 2 (gN/m2/s)
   real, pointer :: sminn_to_denit_l3s3(:)   ! denitrification for decomp. of litter 3 to SOM 3 (gN/m2/s)
   real, pointer :: sminn_to_denit_s1s2(:)   ! denitrification for decomp. of SOM 1 to SOM 2 (gN/m2/s)
   real, pointer :: sminn_to_denit_s2s3(:)   ! denitrification for decomp. of SOM 2 to SOM 3 (gN/m2/s)
   real, pointer :: sminn_to_denit_s3s4(:)   ! denitrification for decomp. of SOM 3 to SOM 4 (gN/m2/s)
   real, pointer :: sminn_to_denit_s4(:)     ! denitrification for decomp. of SOM 4 (gN/m2/s)
   real, pointer :: sminn_to_denit_excess(:) ! denitrification from excess mineral N pool (gN/m2/s)
   ! leaching fluxes
   real, pointer :: sminn_leached(:)         ! soil mineral N pool loss to leaching (gN/m2/s)
   ! dynamic landcover fluxes
   real, pointer :: dwt_seedn_to_leaf(:)      ! (gN/m2/s) seed source to PFT-level
   real, pointer :: dwt_seedn_to_deadstem(:)  ! (gN/m2/s) seed source to PFT-level
   real, pointer :: dwt_conv_nflux(:)         ! (gN/m2/s) conversion N flux (immediate loss to atm)
   real, pointer :: dwt_prod10n_gain(:)       ! (gN/m2/s) addition to 10-yr wood product pool
   real, pointer :: dwt_prod100n_gain(:)      ! (gN/m2/s) addition to 100-yr wood product pool
   real, pointer :: dwt_frootn_to_litr1n(:)   ! (gN/m2/s) fine root to litter due to landcover change
   real, pointer :: dwt_frootn_to_litr2n(:)   ! (gN/m2/s) fine root to litter due to landcover change
   real, pointer :: dwt_frootn_to_litr3n(:)   ! (gN/m2/s) fine root to litter due to landcover change
   real, pointer :: dwt_livecrootn_to_cwdn(:) ! (gN/m2/s) live coarse root to CWD due to landcover change
   real, pointer :: dwt_deadcrootn_to_cwdn(:) ! (gN/m2/s) dead coarse root to CWD due to landcover change
   real, pointer :: dwt_nloss(:)              ! (gN/m2/s) total nitrogen loss from product pools and conversion
   ! wood product pool loss fluxes
   real, pointer :: prod10n_loss(:)           ! (gN/m2/s) decomposition loss from 10-yr wood product pool
   real, pointer :: prod100n_loss(:)          ! (gN/m2/s) decomposition loss from 100-yr wood product pool
   real, pointer :: product_nloss(:)          ! (gN/m2/s) total wood product nitrogen loss
   ! summary (diagnostic) flux variables, not involved in mass balance
   real, pointer :: potential_immob(:)       ! potential N immobilization (gN/m2/s)
   real, pointer :: actual_immob(:)          ! actual N immobilization (gN/m2/s)
   real, pointer :: sminn_to_plant(:)        ! plant uptake of soil mineral N (gN/m2/s)
   real, pointer :: supplement_to_sminn(:)   ! supplemental N supply (gN/m2/s)
   real, pointer :: gross_nmin(:)            ! gross rate of N mineralization (gN/m2/s)
   real, pointer :: net_nmin(:)              ! net rate of N mineralization (gN/m2/s)
   real, pointer :: denit(:)                 ! total rate of denitrification (gN/m2/s)
   real, pointer :: col_ninputs(:)           ! column-level N inputs (gN/m2/s)
   real, pointer :: col_noutputs(:)          ! column-level N outputs (gN/m2/s)
   ! new variables for fire
   real, pointer :: col_fire_nloss(:)        ! total column-level fire N loss (gN/m2/s)
end type column_nflux_type

!----------------------------------------------------
! End definition of structures defined at the column_type level
!----------------------------------------------------
!*******************************************************************************

! gkw: edited to here

!*******************************************************************************
!----------------------------------------------------
! Begin definition of spatial scaling hierarchy
!----------------------------------------------------

!----------------------------------------------------
! define the pft structure
!----------------------------------------------------

type, public :: pft_type

   ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer, pointer :: column(:)        !index into column level quantities
   real, pointer :: wtcol(:)        !weight (relative to column) 
   integer, pointer :: landunit(:)      !index into landunit level quantities
   integer, pointer :: gridcell(:)      !index into gridcell level quantities
   real, pointer :: wtgcell(:)      !weight (relative to gridcell) 

   ! topological mapping functionality
   integer , pointer :: itype(:)        !pft vegetation 

   ! conservation check structures for the pft level
   type(carbon_balance_type)   :: pcbal !carbon balance structure
   type(nitrogen_balance_type) :: pnbal !nitrogen balance structure
   
   ! CN ecophysiological variables
   type(pft_epv_type)    :: pepv        !pft ecophysiological variables
   
   ! state variables defined at the pft level
   type(pft_pstate_type) :: pps         !physical state variables
   type(pft_estate_type) :: pes         !pft energy state
   type(pft_cstate_type) :: pcs         !pft carbon state
   type(pft_nstate_type) :: pns         !pft nitrogen state

   ! flux variables defined at the pft level
   type(pft_cflux_type)  :: pcf         !pft carbon flux
   type(pft_nflux_type)  :: pnf         !pft nitrogen flux
   
end type pft_type

!----------------------------------------------------
! define the column structure
!----------------------------------------------------

type, public :: column_type

   type(pft_type)   :: p       !plant functional type (pft) data structure 

   ! g/l/c/p hierarchy, local g/l/c/p cells only
   integer , pointer :: landunit(:)     !index into landunit level quantities
   real, pointer :: wtlunit(:)      !weight (relative to landunit)
   integer , pointer :: gridcell(:)     !index into gridcell level quantities
   real, pointer :: wtgcell(:)      !weight (relative to gridcell)
   integer , pointer :: pfti(:)         !beginning pft index for each column
   integer , pointer :: pftf(:)         !ending pft index for each column
   integer , pointer :: npfts(:)        !number of pfts for each column
   
   ! topological mapping functionality
   integer , pointer :: itype(:)        !column type

   ! conservation check structures for the column level
   type(carbon_balance_type)   :: ccbal !carbon balance structure
   type(nitrogen_balance_type) :: cnbal !nitrogen balance structure
   
   ! state variables defined at the column level
   type(column_pstate_type) :: cps      !column physical state variables
   type(column_estate_type) :: ces      !column energy state
   type(column_wstate_type) :: cws      !column water state
   type(column_cstate_type) :: ccs      !column carbon state
   type(column_nstate_type) :: cns      !column nitrogen state
   
   ! flux variables defined at the column level
   type(column_wflux_type) :: cwf       !column water flux
   type(column_cflux_type) :: ccf       !column carbon flux
   type(column_nflux_type) :: cnf       !column nitrogen flux

end type column_type

!----------------------------------------------------
! define the geomorphological land unit structure
!----------------------------------------------------

type, public :: landunit_type
   type(column_type) :: c                 !column data structure (soil/snow/canopy columns)

   ! topological mapping functionality
   integer , pointer :: itype(:)        !landunit type
   logical , pointer :: ifspecial(:)    !BOOL: true=>landunit is not vegetated

   ! conservation check structures for the landunit level
   type(carbon_balance_type)   :: lcbal !carbon balance structure
   type(nitrogen_balance_type) :: lnbal !nitrogen balance structure
   
end type landunit_type

!----------------------------------------------------
! define the gridcell structure
!----------------------------------------------------

type, public :: gridcell_type

   type(landunit_type) :: l             !geomorphological landunits

   ! topological mapping functionality, local 1d gdc arrays
   integer , pointer :: gindex(:)       !global index
   real, pointer :: forc_ndep(:) ! nitrogen deposition rate (gN/m2/s)

   ! conservation check structures for the gridcell level
   type(carbon_balance_type)   :: gcbal !carbon balance structure
   type(nitrogen_balance_type) :: gnbal !nitrogen balance structure

end type gridcell_type

!----------------------------------------------------
! define the top-level (model) structure 
!----------------------------------------------------

type, public :: model_type
   ! lower level in hierarch
   type(gridcell_type) :: g    !gridicell data structure
end type model_type

!----------------------------------------------------
! End definition of spatial scaling hierarchy
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Declare single instance of clmtype
!----------------------------------------------------
type(model_type)    , public, target     , save :: clm3

!----------------------------------------------------
! Declare single instance of array of ecophysiological constant types
!----------------------------------------------------
type(pft_epc_type), public, target, save :: pftcon

!
!EOP
!----------------------------------------------------------------------- 
end module clmtype  
