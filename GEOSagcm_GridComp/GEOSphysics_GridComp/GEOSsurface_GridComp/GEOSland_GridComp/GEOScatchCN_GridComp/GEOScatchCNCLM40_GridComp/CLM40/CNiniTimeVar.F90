!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNiniTimeVar
!
! !INTERFACE:
subroutine CNiniTimeVar(begg, endg, begl, endl, begc, endc, begp, endp)

!
! !DESCRIPTION:
! Initializes time varying variables used only in
! coupled carbon-nitrogen mode (CN):
!
! !USES:
   use clmtype
   use clm_varcon  , only: istsoil
   use pftvarcon   , only: noveg
!  use decompMod   , only: get_proc_bounds
!
! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:
! subroutine iniTimeVar in file iniTimeVar.F90
!
! !REVISION HISTORY:
! 10/21/03: Created by Peter Thornton
!
!
! local pointers to implicit in arguments
!
   real, pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
   real, pointer :: woody(:)     ! binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: leafcn(:)    ! leaf C:N (gC/gN)
   real, pointer :: deadwdcn(:)  ! dead wood (xylem and heartwood) C:N (gC/gN)
   integer , pointer :: ivt(:)       ! pft vegetation type
   integer , pointer :: plandunit(:) ! landunit index associated with each pft
   integer , pointer :: clandunit(:) ! landunit index associated with each column
   integer , pointer :: itypelun(:)  ! landunit type
!
! local pointers to implicit out arguments
!
   real, pointer :: forc_hgt_u_pft(:)    !observational height of wind at pft-level [m]
   real, pointer :: annsum_counter(:) ! seconds since last annual accumulator turnover
   real, pointer :: cannsum_npp(:)    ! annual sum of NPP, averaged from pft-level (gC/m2/yr)
   real, pointer :: cannavg_t2m(:)    !annual average of 2m air temperature, averaged from pft-level (K)
   real, pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real, pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real, pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real, pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real, pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real, pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real, pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real, pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slowest pool)
   real, pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real, pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real, pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real, pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real, pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real, pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real, pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real, pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real, pointer :: sminn(:)              ! (gN/m2) soil mineral N
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
   real, pointer :: psnsun(:)             ! sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real, pointer :: psnsha(:)             ! shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real, pointer :: laisun(:)             ! sunlit projected leaf area index
   real, pointer :: laisha(:)             ! shaded projected leaf area index
   real, pointer :: dormant_flag(:)       ! dormancy flag
   real, pointer :: days_active(:)        ! number of days since last dormancy
   real, pointer :: onset_flag(:)         ! onset flag
   real, pointer :: onset_counter(:)      ! onset days counter
   real, pointer :: onset_gddflag(:)      ! onset flag for growing degree day sum
   real, pointer :: onset_fdd(:)          ! onset freezing degree days counter
   real, pointer :: onset_gdd(:)          ! onset growing degree days
   real, pointer :: onset_swi(:)          ! onset soil water index
   real, pointer :: offset_flag(:)        ! offset flag
   real, pointer :: offset_counter(:)     ! offset days counter
   real, pointer :: offset_fdd(:)         ! offset freezing degree days counter
   real, pointer :: offset_swi(:)         ! offset soil water index
   real, pointer :: lgsf(:)               ! long growing season factor [0-1]
   real, pointer :: bglfr(:)              ! background litterfall rate (1/s)
   real, pointer :: bgtr(:)               ! background transfer rate (1/s)
   real, pointer :: dayl(:)               ! daylength (seconds)
   real, pointer :: prev_dayl(:)          ! daylength from previous timestep (seconds)
   real, pointer :: annavg_t2m(:)         ! annual average 2m air temperature (K)
   real, pointer :: tempavg_t2m(:)        ! temporary average 2m air temperature (K)
   real, pointer :: gpp(:)                ! GPP flux before downregulation (gC/m2/s)
   real, pointer :: availc(:)             ! C flux available for allocation (gC/m2/s)
   real, pointer :: xsmrpool_recover(:)   ! C flux assigned to recovery of negative cpool (gC/m2/s)
   real, pointer :: alloc_pnow(:)         ! fraction of current allocation to display as new growth (DIM)
   real, pointer :: c_allometry(:)        ! C allocation index (DIM)
   real, pointer :: n_allometry(:)        ! N allocation index (DIM)
   real, pointer :: plant_ndemand(:)      ! N flux required to support initial GPP (gN/m2/s)
   real, pointer :: tempsum_potential_gpp(:) ! temporary annual sum of plant_ndemand
   real, pointer :: annsum_potential_gpp(:)  ! annual sum of plant_ndemand
   real, pointer :: tempmax_retransn(:)   ! temporary max of retranslocated N pool (gN/m2)
   real, pointer :: annmax_retransn(:)    ! annual max of retranslocated N pool (gN/m2)
   real, pointer :: avail_retransn(:)     ! N flux available from retranslocation pool (gN/m2/s)
   real, pointer :: plant_nalloc(:)       ! total allocated N flux (gN/m2/s)
   real, pointer :: plant_calloc(:)       ! total allocated C flux (gC/m2/s)
   real, pointer :: excess_cflux(:)       ! C flux not allocated due to downregulation (gC/m2/s)
   real, pointer :: downreg(:)            ! fractional reduction in GPP due to N limitation (DIM)
   real, pointer :: tempsum_npp(:)        ! temporary annual sum of NPP
   real, pointer :: annsum_npp(:)         ! annual sum of NPP
   real, pointer :: qflx_drain(:)         ! sub-surface runoff (mm H2O /s)
   ! new variables for fire
   real, pointer :: wf(:)                 ! soil moisture in top 0.5 m
   real, pointer :: me(:)                 ! moisture of extinction (proportion)
   real, pointer :: fire_prob(:)          ! daily fire probability (0-1)
   real, pointer :: mean_fire_prob(:)     ! e-folding mean of daily fire probability (0-1)
   real, pointer :: fireseasonl(:)        ! annual fire season length (days, <= 365)
   real, pointer :: farea_burned(:)       ! timestep fractional area burned (proportion)
   real, pointer :: ann_farea_burned(:)   ! annual total fractional area burned (proportion)
   real, pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real, pointer :: totcolc(:)            ! (gC/m2) total column carbon, incl veg and cpool
   real, pointer :: totecosysc(:)         ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
   real, pointer :: totlitc(:)            ! (gC/m2) total litter carbon
   real, pointer :: totsomc(:)            ! (gC/m2) total soil organic matter carbon

   real, pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real, pointer :: totcoln(:)            ! (gN/m2) total column nitrogen, incl veg
   real, pointer :: totecosysn(:)         ! (gN/m2) total ecosystem nitrogen, incl veg
   real, pointer :: totlitn(:)            ! (gN/m2) total litter nitrogen
   real, pointer :: totsomn(:)            ! (gN/m2) total soil organic matter nitrogen
   real, pointer :: dispvegc(:)           ! (gC/m2) displayed veg carbon, excluding storage and cpool
   real, pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real, pointer :: storvegc(:)           ! (gC/m2) stored vegetation carbon, excluding cpool
   real, pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
   real, pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
   real, pointer :: prev_frootc_to_litter(:)!previous timestep froot C litterfall flux (gC/m2/s)
   real, pointer :: prev_leafc_to_litter(:) !previous timestep leaf C litterfall flux (gC/m2/s)
   real, pointer :: dispvegn(:)           ! (gN/m2) displayed veg nitrogen, excluding storage
   real, pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real, pointer :: storvegn(:)           ! (gN/m2) stored vegetation nitrogen
   real, pointer :: totpftn(:)            ! (gN/m2) total pft-level nitrogen
   real, pointer :: totvegn(:)            ! (gN/m2) total vegetation nitrogen
   real, pointer :: vcmxsha(:)            ! shaded leaf Vcmax (umolCO2/m^2/s)
   real, pointer :: vcmxsun(:)            ! sunlit leaf Vcmax (umolCO2/m^2/s)
   ! dynamic landuse variables
   real, pointer :: seedc(:)              ! (gC/m2) column-level pool for seeding new PFTs
   real, pointer :: prod10c(:)            ! (gC/m2) wood product C pool, 10-year lifespan
   real, pointer :: prod100c(:)           ! (gC/m2) wood product C pool, 100-year lifespan
   real, pointer :: totprodc(:)           ! (gC/m2) total wood product C
   real, pointer :: seedn(:)              ! (gN/m2) column-level pool for seeding new PFTs
   real, pointer :: prod10n(:)            ! (gN/m2) wood product N pool, 10-year lifespan
   real, pointer :: prod100n(:)           ! (gN/m2) wood product N pool, 100-year lifespan
   real, pointer :: totprodn(:)           ! (gN/m2) total wood product N
!
! !LOCAL VARIABLES:
   integer :: g,l,c,p      ! indices
   integer :: begp, endp   ! per-clump/proc beginning and ending pft indices
   integer :: begc, endc   ! per-clump/proc beginning and ending column indices
   integer :: begl, endl   ! per-clump/proc beginning and ending landunit indices
   integer :: begg, endg   ! per-clump/proc gridcell ending gridcell indices
!EOP
!-----------------------------------------------------------------------

    ! assign local pointers at the gridcell level

    ! assign local pointers at the landunit level
    itypelun                       => clm3%g%l%itype

    ! assign local pointers at the column level
    clandunit                      => clm3%g%l%c%landunit
    annsum_counter                 => clm3%g%l%c%cps%annsum_counter
    cannsum_npp                    => clm3%g%l%c%cps%cannsum_npp
    cannavg_t2m                    => clm3%g%l%c%cps%cannavg_t2m
    wf                             => clm3%g%l%c%cps%wf
    me                             => clm3%g%l%c%cps%me
    fire_prob                      => clm3%g%l%c%cps%fire_prob
    mean_fire_prob                 => clm3%g%l%c%cps%mean_fire_prob
    fireseasonl                    => clm3%g%l%c%cps%fireseasonl
    farea_burned                   => clm3%g%l%c%cps%farea_burned
    ann_farea_burned               => clm3%g%l%c%cps%ann_farea_burned
    qflx_drain                     => clm3%g%l%c%cwf%qflx_drain
    cwdc                           => clm3%g%l%c%ccs%cwdc
    litr1c                         => clm3%g%l%c%ccs%litr1c
    litr2c                         => clm3%g%l%c%ccs%litr2c
    litr3c                         => clm3%g%l%c%ccs%litr3c
    soil1c                         => clm3%g%l%c%ccs%soil1c
    soil2c                         => clm3%g%l%c%ccs%soil2c
    soil3c                         => clm3%g%l%c%ccs%soil3c
    soil4c                         => clm3%g%l%c%ccs%soil4c
    
    ! dynamic landuse variables
    seedc                          => clm3%g%l%c%ccs%seedc
    prod10c                        => clm3%g%l%c%ccs%prod10c
    prod100c                       => clm3%g%l%c%ccs%prod100c
    totprodc                       => clm3%g%l%c%ccs%totprodc
    seedn                          => clm3%g%l%c%cns%seedn
    prod10n                        => clm3%g%l%c%cns%prod10n
    prod100n                       => clm3%g%l%c%cns%prod100n
    totprodn                       => clm3%g%l%c%cns%totprodn
    
    cwdn                           => clm3%g%l%c%cns%cwdn
    litr1n                         => clm3%g%l%c%cns%litr1n
    litr2n                         => clm3%g%l%c%cns%litr2n
    litr3n                         => clm3%g%l%c%cns%litr3n
    soil1n                         => clm3%g%l%c%cns%soil1n
    soil2n                         => clm3%g%l%c%cns%soil2n
    soil3n                         => clm3%g%l%c%cns%soil3n
    soil4n                         => clm3%g%l%c%cns%soil4n
    sminn                          => clm3%g%l%c%cns%sminn
    col_ctrunc                     => clm3%g%l%c%ccs%col_ctrunc
    totcolc                        => clm3%g%l%c%ccs%totcolc
    totecosysc                     => clm3%g%l%c%ccs%totecosysc
    totlitc                        => clm3%g%l%c%ccs%totlitc
    totsomc                        => clm3%g%l%c%ccs%totsomc

    col_ntrunc                     => clm3%g%l%c%cns%col_ntrunc
    totcoln                        => clm3%g%l%c%cns%totcoln
    totecosysn                     => clm3%g%l%c%cns%totecosysn
    totlitn                        => clm3%g%l%c%cns%totlitn
    totsomn                        => clm3%g%l%c%cns%totsomn

    ! assign local pointers at the pft level
    ivt                            => clm3%g%l%c%p%itype
    plandunit                      => clm3%g%l%c%p%landunit
    leafc                          => clm3%g%l%c%p%pcs%leafc
    leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
    leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
    frootc                         => clm3%g%l%c%p%pcs%frootc
    frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
    frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
    livestemc                      => clm3%g%l%c%p%pcs%livestemc
    livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
    livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
    deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
    livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
    livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
    livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
    deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
    deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
    deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
    gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
    gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
    cpool                          => clm3%g%l%c%p%pcs%cpool
    xsmrpool                       => clm3%g%l%c%p%pcs%xsmrpool
    forc_hgt_u_pft                 => clm3%g%l%c%p%pps%forc_hgt_u_pft

    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn
    npool                          => clm3%g%l%c%p%pns%npool
    psnsun                         => clm3%g%l%c%p%pcf%psnsun
    psnsha                         => clm3%g%l%c%p%pcf%psnsha
    laisun                         => clm3%g%l%c%p%pps%laisun
    laisha                         => clm3%g%l%c%p%pps%laisha
    dormant_flag                   => clm3%g%l%c%p%pepv%dormant_flag
    days_active                    => clm3%g%l%c%p%pepv%days_active
    onset_flag                     => clm3%g%l%c%p%pepv%onset_flag
    onset_counter                  => clm3%g%l%c%p%pepv%onset_counter
    onset_gddflag                  => clm3%g%l%c%p%pepv%onset_gddflag
    onset_fdd                      => clm3%g%l%c%p%pepv%onset_fdd
    onset_gdd                      => clm3%g%l%c%p%pepv%onset_gdd
    onset_swi                      => clm3%g%l%c%p%pepv%onset_swi
    offset_flag                    => clm3%g%l%c%p%pepv%offset_flag
    offset_counter                 => clm3%g%l%c%p%pepv%offset_counter
    offset_fdd                     => clm3%g%l%c%p%pepv%offset_fdd
    offset_swi                     => clm3%g%l%c%p%pepv%offset_swi
    lgsf                           => clm3%g%l%c%p%pepv%lgsf
    bglfr                          => clm3%g%l%c%p%pepv%bglfr
    bgtr                           => clm3%g%l%c%p%pepv%bgtr
    dayl                           => clm3%g%l%c%p%pepv%dayl
    prev_dayl                      => clm3%g%l%c%p%pepv%prev_dayl
    annavg_t2m                     => clm3%g%l%c%p%pepv%annavg_t2m
    tempavg_t2m                    => clm3%g%l%c%p%pepv%tempavg_t2m
    gpp                            => clm3%g%l%c%p%pepv%gpp
    availc                         => clm3%g%l%c%p%pepv%availc
    xsmrpool_recover                  => clm3%g%l%c%p%pepv%xsmrpool_recover
    alloc_pnow                     => clm3%g%l%c%p%pepv%alloc_pnow
    c_allometry                    => clm3%g%l%c%p%pepv%c_allometry
    n_allometry                    => clm3%g%l%c%p%pepv%n_allometry
    plant_ndemand                  => clm3%g%l%c%p%pepv%plant_ndemand
    tempsum_potential_gpp          => clm3%g%l%c%p%pepv%tempsum_potential_gpp
    annsum_potential_gpp           => clm3%g%l%c%p%pepv%annsum_potential_gpp
    tempmax_retransn               => clm3%g%l%c%p%pepv%tempmax_retransn
    annmax_retransn                => clm3%g%l%c%p%pepv%annmax_retransn
    avail_retransn                 => clm3%g%l%c%p%pepv%avail_retransn
    plant_nalloc                   => clm3%g%l%c%p%pepv%plant_nalloc
    plant_calloc                   => clm3%g%l%c%p%pepv%plant_calloc
    excess_cflux                   => clm3%g%l%c%p%pepv%excess_cflux
    downreg                        => clm3%g%l%c%p%pepv%downreg
    tempsum_npp                    => clm3%g%l%c%p%pepv%tempsum_npp
    annsum_npp                     => clm3%g%l%c%p%pepv%annsum_npp
    dispvegc                       => clm3%g%l%c%p%pcs%dispvegc
    pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc
    storvegc                       => clm3%g%l%c%p%pcs%storvegc
    totpftc                        => clm3%g%l%c%p%pcs%totpftc
    totvegc                        => clm3%g%l%c%p%pcs%totvegc
    prev_frootc_to_litter          => clm3%g%l%c%p%pepv%prev_frootc_to_litter
    prev_leafc_to_litter           => clm3%g%l%c%p%pepv%prev_leafc_to_litter
    dispvegn                       => clm3%g%l%c%p%pns%dispvegn
    pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc
    storvegn                       => clm3%g%l%c%p%pns%storvegn
    totpftn                        => clm3%g%l%c%p%pns%totpftn
    totvegn                        => clm3%g%l%c%p%pns%totvegn
    
    ! assign local pointers for ecophysiological constants
    evergreen                      => pftcon%evergreen
    woody                          => pftcon%woody
    leafcn                         => pftcon%leafcn
    deadwdcn                       => pftcon%deadwdcn

   ! Determine subgrid bounds on this processor
!  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp) ! gkw

   ! Added 5/4/04, PET: initialize forc_hgt_u (gridcell-level),
   ! since this is not initialized before first call to CNVegStructUpdate,
   ! and it is required to set the upper bound for canopy top height.
   ! Changed 3/21/08, KO: still needed but don't have sufficient information 
   ! to set this properly (e.g., pft-level displacement height and roughness 
   ! length). So leave at 30m.
   do p = begp, endp
      forc_hgt_u_pft(p) = 30.
   end do

   ! initialize column-level variables
   do c = begc, endc
      l = clandunit(c)
      if (itypelun(l) == istsoil) then
         ! column physical state variables
         annsum_counter(c) = 0.
         cannsum_npp(c)    = 0.
         cannavg_t2m(c)    = 280.
         wf(c) = 1.0  ! it needs to be non zero so the first time step has no fires
         me(c) = 0.
         fire_prob(c) = 0.
         mean_fire_prob(c) = 0.
         fireseasonl(c) = 0.
         farea_burned(c) = 0.
         ann_farea_burned(c) = 0.
         
         ! needed for CNNLeaching
         qflx_drain(c) = 0.

         ! column carbon state variable initialization
         cwdc(c)   = 0.
         litr1c(c) = 0.
         litr2c(c) = 0.
         litr3c(c) = 0.
         soil1c(c) = 0.
         soil2c(c) = 0.
         soil3c(c) = 0.
         soil4c(c) = 10.
         col_ctrunc(c) = 0.
         totlitc(c)    = 0.
         totsomc(c)    = 0.
         totecosysc(c) = 0.
         totcolc(c)    = 0.

         ! column nitrogen state variables
         cwdn(c)   = cwdc(c) / 500.
         litr1n(c) = litr1c(c) / 90.
         litr2n(c) = litr2c(c) / 90.
         litr3n(c) = litr3c(c) / 90.
         soil1n(c) = soil1c(c) / 12.
         soil2n(c) = soil2c(c) / 12.
         soil3n(c) = soil3c(c) / 10.
         soil4n(c) = soil4c(c) / 10.
         sminn(c) = 0.
         col_ntrunc(c) = 0.
         totlitn(c)    = 0.
         totsomn(c)    = 0.
         totecosysn(c) = 0.
         totcoln(c)    = 0.

	 ! dynamic landcover state variables
         seedc(c)  = 0.
	 prod10c(c)    = 0.
	 prod100c(c)   = 0.
	 totprodc(c)   = 0.
	 seedn(c)      = 0.
	 prod10n(c)    = 0.
	 prod100n(c)   = 0.
	 totprodn(c)   = 0.
	 
	 ! also initialize dynamic landcover fluxes so that they have
	 ! real values on first timestep, prior to calling pftdyn_cnbal
	 clm3%g%l%c%ccf%dwt_seedc_to_leaf(c) = 0.
	 clm3%g%l%c%ccf%dwt_seedc_to_deadstem(c) = 0.
	 clm3%g%l%c%ccf%dwt_conv_cflux(c) = 0.
	 clm3%g%l%c%ccf%dwt_prod10c_gain(c) = 0.
	 clm3%g%l%c%ccf%prod10c_loss(c) = 0.
	 clm3%g%l%c%ccf%dwt_prod100c_gain(c) = 0.
	 clm3%g%l%c%ccf%prod100c_loss(c) = 0.
	 clm3%g%l%c%ccf%dwt_frootc_to_litr1c(c) = 0.
	 clm3%g%l%c%ccf%dwt_frootc_to_litr2c(c) = 0.
	 clm3%g%l%c%ccf%dwt_frootc_to_litr3c(c) = 0.
	 clm3%g%l%c%ccf%dwt_livecrootc_to_cwdc(c) = 0.
	 clm3%g%l%c%ccf%dwt_deadcrootc_to_cwdc(c) = 0.
	 clm3%g%l%c%ccf%dwt_closs(c) = 0.
	 clm3%g%l%c%cnf%dwt_seedn_to_leaf(c) = 0.
	 clm3%g%l%c%cnf%dwt_seedn_to_deadstem(c) = 0.
	 clm3%g%l%c%cnf%dwt_conv_nflux(c) = 0.
	 clm3%g%l%c%cnf%dwt_prod10n_gain(c) = 0.
	 clm3%g%l%c%cnf%prod10n_loss(c) = 0.
	 clm3%g%l%c%cnf%dwt_prod100n_gain(c) = 0.
	 clm3%g%l%c%cnf%prod100n_loss(c) = 0.
	 clm3%g%l%c%cnf%dwt_frootn_to_litr1n(c) = 0.
	 clm3%g%l%c%cnf%dwt_frootn_to_litr2n(c) = 0.
	 clm3%g%l%c%cnf%dwt_frootn_to_litr3n(c) = 0.
	 clm3%g%l%c%cnf%dwt_livecrootn_to_cwdn(c) = 0.
	 clm3%g%l%c%cnf%dwt_deadcrootn_to_cwdn(c) = 0.
	 clm3%g%l%c%cnf%dwt_nloss(c) = 0.
      end if
   end do

   ! initialize pft-level variables
   do p = begp, endp
      l = plandunit(p)
      if (itypelun(l) == istsoil) then
         
         ! carbon state variables
         if (ivt(p) == noveg) then
            leafc(p) = 0.
            leafc_storage(p) = 0.
         else
            if (evergreen(ivt(p)) == 1.) then
               leafc(p) = 1.
               leafc_storage(p) = 0.
            else
               leafc(p) = 0.
               leafc_storage(p) = 1.
            end if
         end if

         leafc_xfer(p) = 0.
         frootc(p) = 0.
         frootc_storage(p) = 0.
         frootc_xfer(p) = 0.
         livestemc(p) = 0.
         livestemc_storage(p) = 0.
         livestemc_xfer(p) = 0.

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(ivt(p)) == 1.) then
            deadstemc(p) = 0.1
         else
            deadstemc(p) = 0.
         end if

         deadstemc_storage(p) = 0.
         deadstemc_xfer(p) = 0.
         livecrootc(p) = 0.
         livecrootc_storage(p) = 0.
         livecrootc_xfer(p) = 0.
         deadcrootc(p) = 0.
         deadcrootc_storage(p) = 0.
         deadcrootc_xfer(p) = 0.
         gresp_storage(p) = 0.
         gresp_xfer(p) = 0.
         cpool(p) = 0.
         xsmrpool(p) = 0.
         pft_ctrunc(p) = 0.
         dispvegc(p) = 0.
         storvegc(p) = 0.
         totpftc(p)  = 0.
         ! calculate totvegc explicitly so that it is available for the isotope 
         ! code on the first time step.
         totvegc(p)  = leafc(p) + leafc_storage(p) + leafc_xfer(p) + frootc(p) +  &
            frootc_storage(p) + frootc_xfer(p) + livestemc(p) + livestemc_storage(p) +  &
            livestemc_xfer(p) + deadstemc(p) + deadstemc_storage(p) + deadstemc_xfer(p) +  &
            livecrootc(p) + livecrootc_storage(p) + livecrootc_xfer(p) + deadcrootc(p) +  &
            deadcrootc_storage(p) + deadcrootc_xfer(p) + gresp_storage(p) +  &
            gresp_xfer(p) + cpool(p)

         ! nitrogen state variables
         if (ivt(p) == noveg) then
            leafn(p) = 0.
            leafn_storage(p) = 0.
         else
            leafn(p) = leafc(p) / leafcn(ivt(p))
            leafn_storage(p) = leafc_storage(p) / leafcn(ivt(p))
         end if

         leafn_xfer(p) = 0.
         frootn(p) = 0.
         frootn_storage(p) = 0.
         frootn_xfer(p) = 0.
         livestemn(p) = 0.
         livestemn_storage(p) = 0.
         livestemn_xfer(p) = 0.

         ! tree types need to be initialized with some stem mass so that
         ! roughness length is not zero in canopy flux calculation

         if (woody(ivt(p)) == 1.) then
            deadstemn(p) = deadstemc(p) / deadwdcn(ivt(p))
         else
            deadstemn(p) = 0.
         end if

         deadstemn_storage(p) = 0.
         deadstemn_xfer(p) = 0.
         livecrootn(p) = 0.
         livecrootn_storage(p) = 0.
         livecrootn_xfer(p) = 0.
         deadcrootn(p) = 0.
         deadcrootn_storage(p) = 0.
         deadcrootn_xfer(p) = 0.
         retransn(p) = 0.
         npool(p) = 0.
         pft_ntrunc(p) = 0.
         dispvegn(p) = 0.
         storvegn(p) = 0.
         totvegn(p)  = 0.
         totpftn(p)  = 0.

         ! initialization for psnsun and psnsha required for
         ! proper arbitrary initialization of allocation routine
         ! in initial ecosysdyn call

         psnsun(p) = 0.
         psnsha(p) = 0.
         laisun(p) = 0.
         laisha(p) = 0.

         ! ecophysiological variables
         ! phenology variables
         dormant_flag(p) = 1.
         days_active(p) = 0.
         onset_flag(p) = 0.
         onset_counter(p) = 0.
         onset_gddflag(p) = 0.
         onset_fdd(p) = 0.
         onset_gdd(p) = 0.
         onset_swi(p) = 0.0
         offset_flag(p) = 0.
         offset_counter(p) = 0.
         offset_fdd(p) = 0.
         offset_swi(p) = 0.
         lgsf(p) = 0.
         bglfr(p) = 0.
         bgtr(p) = 0.
         annavg_t2m(p) = 280.
         tempavg_t2m(p) = 0.

         ! non-phenology variables
         gpp(p) = 0.
         availc(p) = 0.
         xsmrpool_recover(p) = 0.
         alloc_pnow(p) = 1.
         c_allometry(p) = 0.
         n_allometry(p) = 0.
         plant_ndemand(p) = 0.
         tempsum_potential_gpp(p) = 0.
         annsum_potential_gpp(p) = 0.
         tempmax_retransn(p) = 0.
         annmax_retransn(p) = 0.
         avail_retransn(p) = 0.
         plant_nalloc(p) = 0.
         plant_calloc(p) = 0.
         excess_cflux(p) = 0.
         downreg(p) = 0.
         prev_leafc_to_litter(p) = 0.
         prev_frootc_to_litter(p) = 0.
         tempsum_npp(p) = 0.
         annsum_npp(p) = 0.

      end if   ! end of if-istsoil block
   end do   ! end of loop over pfts  

end subroutine CNiniTimeVar
