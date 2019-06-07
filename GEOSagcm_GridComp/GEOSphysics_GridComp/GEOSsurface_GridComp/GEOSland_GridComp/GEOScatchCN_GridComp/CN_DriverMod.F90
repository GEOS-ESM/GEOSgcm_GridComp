module CN_DriverMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CN_Driver
!
! !DESCRIPTION:
! Driver for CN model. Interface between GEOS5 and CLM4 data structures
!
! !USES:
  use clmtype 
  use clmtypeInitMod  
  use CNEcosystemDynMod
  use nanMod
  use clm_varcon,        only: grav, denh2o
  use clm_varpar,        only: clm_varpar_init, numpft
  use CNSetValueMod,     only: CNZeroFluxes_dwt
  use CNAnnualUpdateMod, only: CNAnnualUpdate
  use CNBalanceCheckMod, only: BeginCBalance, BeginNBalance, &
                               CBalanceCheck, NBalanceCheck
  use pftvarcon,         only: noveg
  use clm_time_manager,  only: get_step_size

  implicit none
  private

! local variables to the whole module

! !PUBLIC MEMBER FUNCTIONS:
  public :: CN_Driver
  public :: CN_init
  public :: CN_exit
  public :: get_CN_LAI
!
! !REVISION HISTORY:
! 2011-2015: Created by Greg Walker NASA/GSFC GMAO
!
!EOP
!-----------------------------------------------------------------------

contains

  subroutine CN_Driver(istep,nch,nveg,nzone,daylength,               &
                       tgw,tp1,tp2,tp3,tp4,tp5,tp6,sfmc,rzm,         &
		       wpwet,psis,bee,poros,vgwmax,bflow,totwat,     &
		       tm,psun,psha,lsun,lsha,                       &
                       ityp,fveg,wtzone,sndzn,asnow,ndep,zlai,zsai,ztai,colc,tile_id,ann_t2m, &
		       nppg,gppg,srg,neeg,root,padd,vegc,xsmr,burn,seal,closs,firefac)

  ! !ARGUMENTS:
  implicit none

  integer*8, intent(in) :: istep ! carbon model time step
  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  real*4, intent(in) :: firefac ! gkw: fire tuning parameter. valid range: 0-1; negative value will prevent fires
  real*4, dimension(nch), intent(in) :: daylength ! daylength (seconds)
  real*4, dimension(nch,nzone), intent(in) :: tgw ! soil surface layer temperature (K)
  real*4, dimension(nch,nzone), intent(in) :: rzm ! weighted root-zone moisture content
  real*4, dimension(nch), intent(in) :: tp1,tp2,tp3,tp4,tp5,tp6 ! soil temperatures (K)
  real*4, dimension(nch), intent(in) :: sfmc,wpwet,psis,bee,poros,vgwmax,bflow,totwat ! soil water & parameters
  real*4, dimension(nch), intent(in) :: tm   ! air temperature
  real*4, dimension(nch), intent(in) :: ndep ! nitrogen deposition
  real*4, dimension(nch), intent(in) :: sndzn ! total snow depth
  real*4, dimension(nch), intent(in) :: asnow ! areal snow coverage [0-1]
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! CLM PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg ! catchment vegetation fractions
  real, dimension(nch,nzone), intent(in) :: wtzone ! zone fractions
  real*4, dimension(nch,nveg,nzone), intent(in) :: psha,psun ! photosynthesis
  real*4, dimension(nch,nveg,nzone), intent(in) :: lsha,lsun ! LAI
  integer, dimension(nch), intent(in) :: tile_id ! tile index for debugging
  real*4, dimension(nch), intent(in) :: ann_t2m  ! gkw: annual mean CONUS T 2m to override CN

  real*4, dimension(nch,nveg,nzone), intent(inout) :: zlai ! leaf-area index for tile (subject to burying by snow)
  real*4, dimension(nch,nveg,nzone), intent(inout) :: zsai ! stem-area index for tile
  real*4, dimension(nch,nveg,nzone), intent(inout) :: ztai ! leaf-area index for tile (not buried by snow)
  real*4, dimension(nch,nzone), intent(out) :: colc ! column total carbon

! PFT carbon fluxes averaged to tile 
  real*4, dimension(nch), intent(out) :: nppg ! (gC/m2/s) net primary production [PFT]
  real*4, dimension(nch), intent(out) :: gppg ! (gC/m2/s) gross primary production [PFT]

! column carbon fluxes averaged to tile 
  real*4, dimension(nch), intent(out) :: srg  ! (gC/m2/s) total soil respiration (HR + root resp) [column]
  real*4, dimension(nch), intent(out) :: neeg ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source [column]

! fire diagnostics
  real*4, dimension(nch), intent(out) :: burn, seal, closs

! column & PFT level carbon added to sustain growth
! -------------------------------------------------
  real*4, dimension(nch), intent(out) :: padd  ! (gC/m2/s)
  real*4, dimension(nch), intent(out) :: root, vegc, xsmr  ! (gC/m2)

  integer :: lbg, ubg        ! grid bounds
  integer :: lbl, ubl        ! land-type bounds
  integer :: lbc, ubc        ! column bounds
  integer :: lbp, ubp        ! pft bounds
  integer :: num_soilc       ! number of soil columns in filter
  integer :: num_soilp       ! number of soil pfts in filter
  logical, save :: doalb = .true.         ! assume surface albedo calculation time step
  logical, save :: spin = .false.         ! true if spinup vegetation gkw: critical!
  logical, save :: exit_spin = .false.    ! true if this is first continuation from a spin=true run

  logical, save :: first = .true.
  integer, parameter :: npft = numpft+1 
  integer :: n, c, p, pft_num, idum, pf, i, j, nv, nc, nz, z, icn
  real bare, leafc_tot

  integer, allocatable, save :: filter_soilc(:),filter_soilp(:)
  integer, allocatable, save :: index_soilp(:),zone_soilp(:),zone_soilc(:)
  real, allocatable, save :: leafc_add(:)

  real :: dt   ! time step delta t (seconds)

  logical, pointer :: ifspecial(:)    !BOOL: true=>landunit is not vegetated
  integer, pointer :: litype(:)       !landunit type
  integer, pointer :: citype(:)       !column type
  integer, pointer :: clandunit(:)    !index into landunit level quantities
  integer, pointer :: cgridcell(:)    !index into gridcell level quantities
  integer, pointer :: npfts(:)        !number of pfts for each column
  integer, pointer :: pfti(:)         !beginning pft index for each column
  integer, pointer :: pftf(:)         !ending pft index for each column
  integer, pointer :: pitype(:)       !pft vegetation
  integer, pointer :: pcolumn(:)      !column type
  integer, pointer :: plandunit(:)    !index into landunit level quantities
  integer, pointer :: pgridcell(:)    !index into gridcell level quantities
  real, pointer    :: pwtcol(:)       !weight (relative to column) 
  real, pointer    :: pwtgcell(:)     !weight (relative to gridcell) 

  real, pointer :: t_ref2m(:)        !2 m height surface air temperature (Kelvin)
  real, pointer :: psnsun(:)         !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
  real, pointer :: psnsha(:)         !shaded leaf photosynthesis (umol CO2 /m**2/ s)
  real, pointer :: laisun(:)         !sunlit projected leaf area index
  real, pointer :: laisha(:)         !shaded projected leaf area index
  real, pointer :: forc_ndep(:)      !nitrogen deposition rate (gN/m2/s)
  real, pointer :: t_soisno(:,:)     !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
  real, pointer :: rootfr(:,:)       !fraction of roots in each soil layer  (nlevgrnd)
  real, pointer :: dz(:,:)           !layer thickness (m)  (-nlevsno+1:nlevgrnd)
  real, pointer :: psisat(:,:)       !soil water potential at saturation for CN code (MPa)
  real, pointer :: psiwilt(:)        !root-zone soil water potential at wilting point (MPa)
  real, pointer :: soilpsi(:,:)      !soil water potential in each soil layer (MPa)
  real, pointer :: h2osoi_liq(:)     !column liquid water (kg/m2) (new)
  real, pointer :: wf(:)             !soil water as frac. of whc for top 0.5 m
  real, pointer :: qflx_drain(:)     !sub-surface runoff (mm H2O /s)
  real, pointer :: snowdp(:)         !snow height (m)
  real, pointer :: t_grnd(:)         !ground temperature (Kelvin)
  real, pointer :: forc_hgt_u_pft(:) !wind forcing height (10m+z0m+d) (m)
  real, pointer :: dayl(:)           !daylength (seconds)
  real, pointer :: prev_dayl(:)      !daylength from previous albedo timestep (seconds)
  real, pointer :: elai(:)           !one-sided leaf area index with burying by snow
  real, pointer :: esai(:)           !one-sided stem area index with burying by snow
  real, pointer :: tlai(:)           !one-sided leaf area index, no burying by snow
  real, pointer :: totcolc(:)        ! (gC/m2) total column carbon, incl veg and cpool
  real, pointer :: annavg_t2m(:)     ! annual average 2m air temperature (K)

  real, pointer :: col_ctrunc(:)	 ! (gC/m2) column-level sink for C truncation
  real, pointer :: totecosysc(:)	 ! (gC/m2) total ecosystem carbon, incl veg but excl cpool
  real, pointer :: totlitc(:)		 ! (gC/m2) total litter carbon
  real, pointer :: totsomc(:)		 ! (gC/m2) total soil organic matter carbon
  real, pointer :: cwdc(:), litr1c(:), litr2c(:), litr3c(:), leafc(:)
  real, pointer :: soil1c(:), soil2c(:), soil3c(:), soil4c(:), leafc_storage(:)
  real, pointer :: leafc_xfer(:)	 ! (gC/m2) leaf C transfer
  real, pointer :: frootc(:)		 ! (gC/m2) fine root C
  real, pointer :: frootc_storage(:)	 ! (gC/m2) fine root C storage
  real, pointer :: frootc_xfer(:)	 ! (gC/m2) fine root C transfer
  real, pointer :: livestemc(:) 	 ! (gC/m2) live stem C
  real, pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
  real, pointer :: livestemc_xfer(:)	 ! (gC/m2) live stem C transfer
  real, pointer :: deadstemc(:) 	 ! (gC/m2) dead stem C
  real, pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
  real, pointer :: deadstemc_xfer(:)	 ! (gC/m2) dead stem C transfer
  real, pointer :: livecrootc(:)	 ! (gC/m2) live coarse root C
  real, pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
  real, pointer :: livecrootc_xfer(:)	 ! (gC/m2) live coarse root C transfer
  real, pointer :: deadcrootc(:)	 ! (gC/m2) dead coarse root C
  real, pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
  real, pointer :: deadcrootc_xfer(:)	 ! (gC/m2) dead coarse root C transfer
  real, pointer :: gresp_storage(:)	 ! (gC/m2) growth respiration storage
  real, pointer :: gresp_xfer(:)	 ! (gC/m2) growth respiration transfer
  real, pointer :: cpool(:)		 ! (gC/m2) temporary photosynthate C pool
  real, pointer :: xsmrpool(:)  	 ! (gC/m2) abstract C pool to meet excess MR demand
  real, pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
  real, pointer :: totvegc(:)            ! (gC/m2) total vegetation carbon, excluding cpool
  real, pointer :: totpftc(:)            ! (gC/m2) total pft-level carbon, including cpool
  real, pointer :: cwdn(:)		 ! (gN/m2) coarse woody debris N
  real, pointer :: litr1n(:)		 ! (gN/m2) litter labile N
  real, pointer :: litr2n(:)		 ! (gN/m2) litter cellulose N
  real, pointer :: litr3n(:)		 ! (gN/m2) litter lignin N
  real, pointer :: soil1n(:)		 ! (gN/m2) soil organic matter N (fast pool)
  real, pointer :: soil2n(:)		 ! (gN/m2) soil organic matter N (medium pool)
  real, pointer :: soil3n(:)		 ! (gN/m2) soil orgainc matter N (slow pool)
  real, pointer :: soil4n(:)		 ! (gN/m2) soil orgainc matter N (slowest pool)
  real, pointer :: sminn(:)		 ! (gN/m2) soil mineral N
  real, pointer :: leafn(:)		 ! (gN/m2) leaf N
  real, pointer :: leafn_storage(:)	 ! (gN/m2) leaf N storage
  real, pointer :: leafn_xfer(:)	 ! (gN/m2) leaf N transfer
  real, pointer :: frootn(:)		 ! (gN/m2) fine root N
  real, pointer :: frootn_storage(:)	 ! (gN/m2) fine root N storage
  real, pointer :: frootn_xfer(:)	 ! (gN/m2) fine root N transfer
  real, pointer :: livestemn(:) 	 ! (gN/m2) live stem N
  real, pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
  real, pointer :: livestemn_xfer(:)	 ! (gN/m2) live stem N transfer
  real, pointer :: deadstemn(:) 	 ! (gN/m2) dead stem N
  real, pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
  real, pointer :: deadstemn_xfer(:)	 ! (gN/m2) dead stem N transfer
  real, pointer :: livecrootn(:)	 ! (gN/m2) live coarse root N
  real, pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
  real, pointer :: livecrootn_xfer(:)	 ! (gN/m2) live coarse root N transfer
  real, pointer :: deadcrootn(:)	 ! (gN/m2) dead coarse root N
  real, pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
  real, pointer :: deadcrootn_xfer(:)	 ! (gN/m2) dead coarse root N transfer
  real, pointer :: retransn(:)  	 ! (gN/m2) plant pool of retranslocated N
  real, pointer :: npool(:)		 ! (gN/m2) temporary plant N pool
  real, pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
  real, pointer :: col_ntrunc(:)	 ! (gN/m2) column-level sink for N truncation
  real, pointer :: totcoln(:)		 ! (gN/m2) total column nitrogen, incl veg
  real, pointer :: totecosysn(:)	 ! (gN/m2) total ecosystem nitrogen, incl veg
  real, pointer :: totlitn(:)		 ! (gN/m2) total litter nitrogen
  real, pointer :: totsomn(:)		 ! (gN/m2) total soil organic matter nitrogen
  real, pointer :: seedc(:)		 ! (gC/m2) column-level pool for seeding new PFTs
  real, pointer :: prod10c(:)		 ! (gC/m2) wood product C pool, 10-year lifespan
  real, pointer :: prod100c(:)  	 ! (gC/m2) wood product C pool, 100-year lifespan
  real, pointer :: totprodc(:)  	 ! (gC/m2) total wood product C
  real, pointer :: seedn(:)		 ! (gN/m2) column-level pool for seeding new PFTs
  real, pointer :: prod10n(:)		 ! (gN/m2) wood product N pool, 10-year lifespan
  real, pointer :: prod100n(:)  	 ! (gN/m2) wood product N pool, 100-year lifespan
  real, pointer :: totprodn(:)  	 ! (gN/m2) total wood product N
  real, pointer :: gpp(:)                ! (gC/m2/s) gross primary production 
  real, pointer :: npp(:)                ! (gC/m2/s) net primary production
  real, pointer :: sr(:)                 ! (gC/m2/s) total soil respiration (HR + root resp)
  real, pointer :: nee(:)                ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source
  real, pointer :: farea_burned(:)       ! fractional area burned in this timestep (proportion)
  real, pointer :: me(:)                 ! moisture of extinction (proportion) 
  real, pointer :: fireseasonl(:)        ! annual fire season length (days, <= 365)
  real, pointer :: col_fire_closs(:)     ! (gC/m2/s) total column-level fire C loss

! define size of grid, landunit, column, and PFT
! ----------------------------------------------
  lbg = 1 ; ubg = nch
  lbl = 1 ; ubl = nch
  lbc = 1 ; ubc = nch*nzone
  lbp = 1 ; ubp = nch*nzone*npft ! potential PFT index (0-16); actual will be set in num_soilp filter

  if(first) then
    allocate (filter_soilc(ubc-lbc+1))
    allocate (filter_soilp(ubp-lbp+1))
    allocate (index_soilp(ubp-lbp+1))
    allocate (leafc_add(ubp-lbp+1))
    allocate (zone_soilp(ubp-lbp+1))
    allocate (zone_soilc(ubc-lbc+1))
  endif

! assign local pointers to derived type arrays
! --------------------------------------------
  litype        => clm3%g%l%itype
  ifspecial     => clm3%g%l%ifspecial
  citype        => clm3%g%l%c%itype
  clandunit     => clm3%g%l%c%landunit
  npfts         => clm3%g%l%c%npfts
  cgridcell     => clm3%g%l%c%gridcell
  pfti          => clm3%g%l%c%pfti
  pftf          => clm3%g%l%c%pftf

  pitype        => clm3%g%l%c%p%itype
  pcolumn       => clm3%g%l%c%p%column
  pwtcol        => clm3%g%l%c%p%wtcol
  t_ref2m       => clm3%g%l%c%p%pes%t_ref2m
  psnsha        => clm3%g%l%c%p%pcf%psnsha
  psnsun        => clm3%g%l%c%p%pcf%psnsun
  laisha        => clm3%g%l%c%p%pps%laisha
  laisun        => clm3%g%l%c%p%pps%laisun
  plandunit     => clm3%g%l%c%p%landunit
  pgridcell     => clm3%g%l%c%p%gridcell
  pwtgcell      => clm3%g%l%c%p%wtgcell
  forc_ndep     => clm3%g%forc_ndep
  t_soisno      => clm3%g%l%c%ces%t_soisno
  rootfr        => clm3%g%l%c%p%pps%rootfr
  dz            => clm3%g%l%c%cps%dz
  psisat        => clm3%g%l%c%cps%psisat
  psiwilt       => clm3%g%l%c%cps%psiwilt
  soilpsi       => clm3%g%l%c%cps%soilpsi
  h2osoi_liq    => clm3%g%l%c%cws%h2osoi_liq
  wf            => clm3%g%l%c%cps%wf
  qflx_drain    => clm3%g%l%c%cwf%qflx_drain
  snowdp        => clm3%g%l%c%cps%snowdp
  t_grnd        => clm3%g%l%c%ces%t_grnd
  forc_hgt_u_pft=> clm3%g%l%c%p%pps%forc_hgt_u_pft
  dayl          => clm3%g%l%c%p%pepv%dayl
  prev_dayl     => clm3%g%l%c%p%pepv%prev_dayl
  elai          => clm3%g%l%c%p%pps%elai
  esai          => clm3%g%l%c%p%pps%esai
  tlai          => clm3%g%l%c%p%pps%tlai
  totcolc       => clm3%g%l%c%ccs%totcolc
  annavg_t2m    => clm3%g%l%c%p%pepv%annavg_t2m
  me            => clm3%g%l%c%cps%me

!  prevent negative carbon
!  -----------------------
  cwdc  			 => clm3%g%l%c%ccs%cwdc
  litr1c			 => clm3%g%l%c%ccs%litr1c
  litr2c			 => clm3%g%l%c%ccs%litr2c
  litr3c			 => clm3%g%l%c%ccs%litr3c
  soil1c			 => clm3%g%l%c%ccs%soil1c
  soil2c			 => clm3%g%l%c%ccs%soil2c
  soil3c			 => clm3%g%l%c%ccs%soil3c
  soil4c			 => clm3%g%l%c%ccs%soil4c

  col_ctrunc			 => clm3%g%l%c%ccs%col_ctrunc
  totecosysc			 => clm3%g%l%c%ccs%totecosysc
  totlitc			 => clm3%g%l%c%ccs%totlitc
  totsomc			 => clm3%g%l%c%ccs%totsomc

  cwdn  			 => clm3%g%l%c%cns%cwdn
  litr1n			 => clm3%g%l%c%cns%litr1n
  litr2n			 => clm3%g%l%c%cns%litr2n
  litr3n			 => clm3%g%l%c%cns%litr3n
  soil1n			 => clm3%g%l%c%cns%soil1n
  soil2n			 => clm3%g%l%c%cns%soil2n
  soil3n			 => clm3%g%l%c%cns%soil3n
  soil4n			 => clm3%g%l%c%cns%soil4n
  sminn 			 => clm3%g%l%c%cns%sminn
  col_ctrunc			 => clm3%g%l%c%ccs%col_ctrunc
  totcolc			 => clm3%g%l%c%ccs%totcolc
  totecosysc			 => clm3%g%l%c%ccs%totecosysc
  totlitc			 => clm3%g%l%c%ccs%totlitc
  totsomc			 => clm3%g%l%c%ccs%totsomc

  col_ntrunc			 => clm3%g%l%c%cns%col_ntrunc
  totcoln			 => clm3%g%l%c%cns%totcoln
  totecosysn			 => clm3%g%l%c%cns%totecosysn
  totlitn			 => clm3%g%l%c%cns%totlitn
  totsomn			 => clm3%g%l%c%cns%totsomn

  seedc 			 => clm3%g%l%c%ccs%seedc
  prod10c			 => clm3%g%l%c%ccs%prod10c
  prod100c			 => clm3%g%l%c%ccs%prod100c
  totprodc			 => clm3%g%l%c%ccs%totprodc
  seedn 			 => clm3%g%l%c%cns%seedn
  prod10n			 => clm3%g%l%c%cns%prod10n
  prod100n			 => clm3%g%l%c%cns%prod100n
  totprodn			 => clm3%g%l%c%cns%totprodn

  leafc 			 => clm3%g%l%c%p%pcs%leafc
  leafc_storage 		 => clm3%g%l%c%p%pcs%leafc_storage
  leafc_xfer			 => clm3%g%l%c%p%pcs%leafc_xfer
  frootc			 => clm3%g%l%c%p%pcs%frootc
  frootc_storage		 => clm3%g%l%c%p%pcs%frootc_storage
  frootc_xfer			 => clm3%g%l%c%p%pcs%frootc_xfer
  livestemc			 => clm3%g%l%c%p%pcs%livestemc
  livestemc_storage		 => clm3%g%l%c%p%pcs%livestemc_storage
  livestemc_xfer		 => clm3%g%l%c%p%pcs%livestemc_xfer
  deadstemc			 => clm3%g%l%c%p%pcs%deadstemc
  deadstemc_storage		 => clm3%g%l%c%p%pcs%deadstemc_storage
  deadstemc_xfer		 => clm3%g%l%c%p%pcs%deadstemc_xfer
  livecrootc			 => clm3%g%l%c%p%pcs%livecrootc
  livecrootc_storage		 => clm3%g%l%c%p%pcs%livecrootc_storage
  livecrootc_xfer		 => clm3%g%l%c%p%pcs%livecrootc_xfer
  deadcrootc			 => clm3%g%l%c%p%pcs%deadcrootc
  deadcrootc_storage		 => clm3%g%l%c%p%pcs%deadcrootc_storage
  deadcrootc_xfer		 => clm3%g%l%c%p%pcs%deadcrootc_xfer
  gresp_storage 		 => clm3%g%l%c%p%pcs%gresp_storage
  gresp_xfer			 => clm3%g%l%c%p%pcs%gresp_xfer
  cpool 			 => clm3%g%l%c%p%pcs%cpool
  xsmrpool			 => clm3%g%l%c%p%pcs%xsmrpool
  pft_ctrunc                     => clm3%g%l%c%p%pcs%pft_ctrunc
  totvegc                        => clm3%g%l%c%p%pcs%totvegc
  totpftc                        => clm3%g%l%c%p%pcs%totpftc

  leafn 			 => clm3%g%l%c%p%pns%leafn
  leafn_storage 		 => clm3%g%l%c%p%pns%leafn_storage
  leafn_xfer			 => clm3%g%l%c%p%pns%leafn_xfer
  frootn			 => clm3%g%l%c%p%pns%frootn
  frootn_storage		 => clm3%g%l%c%p%pns%frootn_storage
  frootn_xfer			 => clm3%g%l%c%p%pns%frootn_xfer
  livestemn			 => clm3%g%l%c%p%pns%livestemn
  livestemn_storage		 => clm3%g%l%c%p%pns%livestemn_storage
  livestemn_xfer		 => clm3%g%l%c%p%pns%livestemn_xfer
  deadstemn			 => clm3%g%l%c%p%pns%deadstemn
  deadstemn_storage		 => clm3%g%l%c%p%pns%deadstemn_storage
  deadstemn_xfer		 => clm3%g%l%c%p%pns%deadstemn_xfer
  livecrootn			 => clm3%g%l%c%p%pns%livecrootn
  livecrootn_storage		 => clm3%g%l%c%p%pns%livecrootn_storage
  livecrootn_xfer		 => clm3%g%l%c%p%pns%livecrootn_xfer
  deadcrootn			 => clm3%g%l%c%p%pns%deadcrootn
  deadcrootn_storage		 => clm3%g%l%c%p%pns%deadcrootn_storage
  deadcrootn_xfer		 => clm3%g%l%c%p%pns%deadcrootn_xfer
  retransn			 => clm3%g%l%c%p%pns%retransn
  npool 			 => clm3%g%l%c%p%pns%npool
  pft_ntrunc                     => clm3%g%l%c%p%pns%pft_ntrunc

! define landunit & column settings
! ---------------------------------
  litype(:)       = 1       ! all land-units are soil
  ifspecial(:)    = .false. ! no special land-units; all are soil

  citype(:)       = 1       ! all columns are soil (vegetated or bare)
  clandunit(:)    = 1       ! all landunits are soil
  plandunit(:)    = 1       ! all landunits are soil
  npfts(:)        = npft    ! max number of PFTs per column

  num_soilc = nch*nzone     ! number of columns = number of catchments*zones

! map vegetation types into closest PFT & map PFT into column; assign PFT weight & filter
! ---------------------------------------------------------------------------------------
  num_soilp = 0             ! initialize PFT filter
  filter_soilp(:) = 0       ! set PFT index to invalid number
  index_soilp(:) = 0        ! set veg index to invalid number
  zone_soilp(:) = 0         ! set zone index to invalid number

! loop over catchment tiles [columns]
! -----------------------------------
  n = 0
  
  do nc = 1,nch

! loop over zones
! ---------------
  do nz = 1,nzone

    n = n + 1

    filter_soilc(n) = n           ! 1:1 mapping catchment to column
    zone_soilc(n) = nz            ! for remapping column to tile
    cgridcell(n) = nc             ! catchment
    pfti(n)     = npft*(n-1) + 1  ! starting PFT index
    pftf(n)     = npft*n          ! ending PFT index

    bare = 1.                     ! bare soil for this tile
    do nv = 1,nveg
      bare = bare - fveg(nc,nv,nz)! subtract vegetated fractions
    end do
    if(bare .lt. 1.e-4) bare = 0. ! don't bother with small bare fractions

    pft_num = 0                   ! PFT loop
    do p = pfti(n),pftf(n)
     pitype(p) = pft_num          ! PFT index
     pcolumn(p) = n               ! column index for PFT 
     pwtcol(p) = 0.               ! weight will be zero unless otherwise set
     t_ref2m(p) = tm(nc)          ! 2m air temperature
     pgridcell(p) = nc            ! PFT map into catchment tile
     pwtgcell(p) = 0.             ! PFT weight in catchment tile

! map bare soil if present
! ------------------------
     if(bare.gt.0. .and. pft_num.eq.0) then
       num_soilp = num_soilp + 1
       filter_soilp(num_soilp) = p
       pwtcol(p) = bare
       psnsha(p) = 0.
       psnsun(p) = 0.
       laisha(p) = 0.
       laisun(p) = 0.
       forc_hgt_u_pft(p) = 30. ! gkw: may need from land model; use this for now
       rootfr(p,1) = 0.0
       pwtgcell(p) = bare*wtzone(nc,nz) ! PFT weight in catchment tile
     endif

! map vegetation type
! -------------------
     do nv = 1,nveg
       if(ityp(nc,nv,nz).eq.pft_num .and. fveg(nc,nv,nz).gt.1.e-4) then
         num_soilp = num_soilp + 1
         filter_soilp(num_soilp) = p
         index_soilp(num_soilp) = nv ! for remapping LAI to tile
         zone_soilp(num_soilp)  = nz ! for remapping LAI to tile
         pwtcol(p) = fveg(nc,nv,nz)
         psnsha(p) = psha(nc,nv,nz)
         psnsun(p) = psun(nc,nv,nz)
         laisha(p) = lsha(nc,nv,nz)
         laisun(p) = lsun(nc,nv,nz)
         forc_hgt_u_pft(p) = 30. ! gkw: may need from land model; use this for now
         rootfr(p,1) = 1.0 ! gkw: affects maint resp. test sensitivity
         pwtgcell(p) = fveg(nc,nv,nz)*wtzone(nc,nz) ! PFT weight in catchment tile

! set daylength here (moved from CNPhenology)

         if(first .and. istep==0) dayl(p) = daylength(nc) ! working fix for cold carbon start gkw: 2015-08-01
         prev_dayl(p) = dayl(p)
         dayl(p) = daylength(nc)
       endif
     end do

     pft_num = pft_num + 1
    end do ! end PFT loop

! set column soil properties
! --------------------------
    dz(n,1) = 1.00          ! hydrologically active soil layer thickness (m)
    t_soisno(n,1) = tp1(nc) ! soil layer temperature (K)
    t_grnd(n) = tgw(nc,nz)  ! ground surface temperature (K)

    psisat(n,1) = 1.e-6*psis(nc)*grav*denh2o ! saturated soil water potential m -> Mpa 
    soilpsi(n,1) = 1.e-6*psis(nc)*grav*denh2o*rzm(nc,nz)**(-bee(nc))
!!! psiwilt(n) = -2. ! gkw: CN default value                      ! root-zone wilting soil water potential (Mpa) 
    psiwilt(n) = 1.e-6*psis(nc)*grav*denh2o*wpwet(nc)**(-bee(nc)) ! root-zone wilting soil water potential (Mpa) 

! soil liquid water
! ------------------
    h2osoi_liq(n) = totwat(nc) ! soil liquid water, kg/m2
    qflx_drain(n) = bflow(nc)  ! sub-surface runoff (mm H2O /s); used for soil nitrogen leaching
    wf(n) = sfmc(nc)/poros(nc) ! soil water as frac of WHC for top 5cm; used in CNFireMod

    if(firefac >= 0.) then
      me(n) = wpwet(nc) + firefac*(1. - wpwet(nc))  ! set fire moisture of extinction as a function of WPWET
     else
      me(n) = -1.e6                                 ! gkw: this will prevent fires
    endif

    snowdp(n) = sndzn(nc)*asnow(nc)      ! snow depth (m)
  end do ! end zone loop

  end do ! end catchment loop

! nitrogen deposition
! -------------------
  forc_ndep(:) = ndep(:)

! initialize column & PFT for cold carbon restart if istep=0 (only if no carbon restart found)
! --------------------------------------------------------------------------------------------
  if(first .and. istep==0) then
    call CNiniTimeVar(lbg,ubg,lbl,ubl,lbc,ubc,lbp,ubp) ! Set arbitrary initial conditions for time varying fields used in coupled carbon-nitrogen code
    print *, 'warning: CN model cold carbon state!'
  endif

! gkw: override annual mean 2m T here, since it exists on restart 
  n = 0  
  do nc = 1,nch   ! loop over catchment tiles [columns]
  do nz = 1,nzone ! loop over zones
    n = n + 1
    do p = pfti(n),pftf(n)
      annavg_t2m(p) = ann_t2m(nc)
    end do        ! end PFT loop
  end do          ! end zone loop
  end do          ! end catchment loop

! initialize balance checks
!  ------------------------
  call BeginCBalance(lbc, ubc, num_soilc, filter_soilc)
  call BeginNBalance(lbc, ubc, num_soilc, filter_soilc)

  call CNZeroFluxes_dwt(lbc, ubc, lbp, ubp)

  call CNEcosystemDyn(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                      num_soilp, filter_soilp, doalb, spin)

  call CNAnnualUpdate(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                      num_soilp, filter_soilp)

! check the carbon and nitrogen balance gkw: don't do balance check on first call of cold start
! -------------------------------------
  if(.not.first) then
  call CBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
  call NBalanceCheck(lbc, ubc, num_soilc, filter_soilc)
  else
  first = .false.
  endif

! keep leaf carbon above some minimum value
! -----------------------------------------
  dt = real( get_step_size() )

  do n = 1,num_soilp
    p = filter_soilp(n)
    leafc_add(p) = 0.

    leafc_tot = leafc(p) + leafc_storage(p) + leafc_xfer(p)
    if(leafc_tot < 0.3333) then ! gkw: arbitrary carbon threshold (g/m2)
      if (pftcon%evergreen(pitype(p)) == 1.) then
        leafc(p) = max(leafc(p),0.3333)
        leafc_storage(p) = max(leafc_storage(p),0.)
       else
        leafc(p) = max(leafc(p),0.)
        leafc_storage(p) = max(leafc_storage(p),0.3333)
      end if
      c = pcolumn(p) ! PFT column index

      leafc_add(p) = leafc(p) + leafc_storage(p) + leafc_xfer(p) - leafc_tot
      totcolc(c) = totcolc(c) + leafc_add(p)*pwtcol(p) ! correct carbon balance
      leafc_add(p) = leafc_add(p)/dt
    endif
  end do

! PFT level diags 
! ---------------
  zlai = 0.    ! elai
  zsai = 0.    ! esai
  ztai = -999. ! tlai
  nppg = 0.
  gppg = 0.
  padd = 0.
  root = 0.
  vegc = 0.
  xsmr = 0.

  npp => clm3%g%l%c%p%pcf%npp
  gpp => clm3%g%l%c%p%pcf%gpp

  do n = 1,num_soilp
    p = filter_soilp(n)
    i = index_soilp(n)
    z = zone_soilp(n)
    c = pgridcell(p)
    if(i .gt. 0) then
      zlai(c,i,z) = elai(p)
      zsai(c,i,z) = esai(p)
      ztai(c,i,z) = tlai(p)
      nppg(c) = nppg(c) + npp(p)*pwtgcell(p)
      gppg(c) = gppg(c) + gpp(p)*pwtgcell(p)
      padd(c) = padd(c) + leafc_add(p)*pwtgcell(p)
      root(c) = root(c) + (frootc(p)+frootc_storage(p)+frootc_xfer(p))*pwtgcell(p)
      vegc(c) = vegc(c) + totvegc(p)*pwtgcell(p)
      xsmr(c) = xsmr(c) + xsmrpool(p)*pwtgcell(p)
    endif
  end do

! column level diags
! ------------------
  srg = 0.
  neeg = 0.
  burn = 0.
  seal = 0.
  closs = 0.

  sr  => clm3%g%l%c%ccf%sr
  nee => clm3%g%l%c%ccf%nee
  farea_burned => clm3%g%l%c%cps%farea_burned
  fireseasonl => clm3%g%l%c%cps%fireseasonl
  col_fire_closs => clm3%g%l%c%ccf%col_fire_closs

  do n = 1,num_soilc
    i = filter_soilc(n)
    z = zone_soilc(n)
    c = cgridcell(i)
    colc(c,z) = totcolc(i)
    srg(c)  = srg(c)  + sr(i) *wtzone(c,z)
    neeg(c) = neeg(c) + nee(i)*wtzone(c,z)
    burn(c) = burn(c) + (farea_burned(i)/dt)*wtzone(c,z)   ! burn rate (fraction per second)
    seal(c) = seal(c) + fireseasonl(i)*wtzone(c,z)
    closs(c) = closs(c) + col_fire_closs(i)*wtzone(c,z)
  end do

  end subroutine CN_Driver

  subroutine CN_init(istep,nch,nveg,nzone,ityp,fveg,cncol,var_col,cnpft,var_pft)

  integer*8, intent(in) :: istep
  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg    ! PFT fraction

  integer, intent(in) :: var_col ! number of CN column restart variables
  real*4, dimension(nch,nzone,var_col), intent(in) :: cncol ! gkw: column CN restart 

  integer, intent(in) :: var_pft ! number of CN PFT restart variables
  real*4, dimension(nch,nzone,nveg,var_pft), intent(in) :: cnpft ! gkw: PFT CN restart 

  integer :: n, p, nv, nc, nz, np

! PFT parameters note: index 0 is "noveg"

  real, save, dimension(0:numpft) :: z0mx,displax,c3psx,vcmx2x,mx,qe2x,slatox,dsladlax,leafcx,flnx,fnitx,woodx,lflitcx
  real, save, dimension(0:numpft) :: frootcx,livewdcx,deadwdcx,dwoox,froot_leax,stem_leax,croot_stex,flivewx,fcux
  real, save, dimension(0:numpft) :: lf_flax,lf_fcex,lf_flix,fr_flax,fr_fcex,fr_flix,leaf_lonx,evergreex,resisx
  real, save, dimension(0:numpft) :: stress_decix,season_decix,xlx,rholx,rhosx,taulx,tausx

  data     z0mx    /   0.,0.055,0.055,0.055,0.075,0.075,0.055,0.055,0.055,0.120,0.120,0.120,0.120,0.120,0.120,0.120,0.120,0.120,0.120,0.120/
  data  displax    /   0.,0.670,0.670,0.670,0.670,0.670,0.670,0.670,0.670,0.680,0.680,0.680,0.680,0.680,0.680,0.680,0.680,0.680,0.680,0.680/
  data    c3psx    /   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,   0.,   1.,   1./
  data   vcmx2x    /   0.,  51.,  43.,  51.,  75.,  69.,  40.,  51.,  51.,  17.,  17.,  17.,  33.,  43.,  43.,  43.,  24.,  24.,  50.,  50./
  data       mx    /   9.,   6.,   6.,   6.,   9.,   9.,   9.,   9.,   9.,   9.,   9.,   9.,   9.,   9.,   9.,   9.,   5.,   5.,   9.,   9./
  data     qe2x    /   0.,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.060,0.040,0.040,0.060,0.060/
  data   slatox    /   0.,0.010,0.008,0.024,0.012,0.012,0.030,0.030,0.030,0.012,0.030,0.030,0.030,0.030,0.030,0.030,0.030,0.030,0.030,0.030/
  data dsladlax    /   0.,0.00125,0.00100,0.00300,0.00150,0.00150,0.00400,0.00400,0.00400,0.,0.,0.,0.,0.,  0.,   0.,   0.,   0.,   0.,   0./
  data   leafcx    /   1., 35.0, 40.0, 25.0, 30.0, 30.0, 25.0, 25.0, 25.0, 30.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0/
  data     flnx    /   0.,0.050,0.040,0.080,0.060,0.060,0.090,0.090,0.090,0.060,0.090,0.090,0.090,0.090,0.090,0.090,0.090,0.090,0.100,0.100/
  data    fnitx    /   0.,0.720,0.780,0.790,0.830,0.710,0.660,0.640,0.700,0.620,0.600,0.600,0.760,0.680,0.610,0.610,0.640,0.640,0.610,0.610/
  data    woodx    /   0.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,   0.,   0.,   0.,   0.,   0.,   0./
  data  lflitcx    /   1., 70.0, 80.0, 50.0, 60.0, 60.0, 50.0, 50.0, 50.0, 60.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0/
  data  frootcx    /   1., 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0/
  data livewdcx    /   1., 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data deadwdcx    /   1.,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data    dwoox    /2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5,2.5e5/
  data froot_leax  /   0.,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  2.0,  2.0,  2.0,  2.0,  2.0,  2.0/
  data stem_leax   /   0., -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 ,-1.0, -1.0,  0.2,  0.2,  0.2,  0.2,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data croot_stex  /   0.,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.3,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data flivewx     /   0.,  0.1,  0.1,  0.1  ,0.1,  0.1,  0.1  ,0.1,  0.1,  0.5,  0.5,  0.5,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
!!data fcux        /   0.,  1.0,  1.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data fcux        /   0.,  1.0,  1.0,  0.5,  1.0,  1.0,  0.5,  0.5,  0.5,  1.0,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5/
  data lf_flax     /   0., 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25/
  data lf_fcex     /   0., 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50/
  data lf_flix     /   0., 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25/
  data fr_flax     /   0., 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25/
  data fr_fcex     /   0., 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50/
  data fr_flix     /   0., 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25/
  data leaf_lonx   /   0.,  3.0,  6.0,  1.0,  1.5,  1.5,  1.0,  1.0,  1.0,  1.5,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0/
  data evergreex   /   0.,  1.0,  1.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data resisx      /   1., 0.12, 0.12, 0.12, 0.12, 0.12, 0.12 ,0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 1.00, 1.00/
  data stress_decix/   0.,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0/
  data season_decix/   0.,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0/
  data xlx         /   0., 0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, 0.25, 0.01, 0.25, 0.25, 0.25,-0.30,-0.30,-0.30,-0.30,-0.30,-0.30,-0.30/
  data rholx       /   0., 0.07, 0.07, 0.07, 0.10, 0.10, 0.10, 0.10, 0.10, 0.07, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11/
  data rhosx       /   0., 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31, 0.31/
  data taulx       /   0., 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05/
  data tausx       /   0.,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.120,0.120,0.120,0.120,0.120,0.120,0.120/

  integer :: lbg, ubg        ! grid bounds
  integer :: lbl, ubl        ! land-type bounds
  integer :: lbc, ubc        ! column bounds
  integer :: lbp, ubp        ! pft bounds

! define size of grid, landunit, column, and PFT
! ----------------------------------------------
  lbg = 1 ; ubg = nch                  ! "grid" (tile)
  lbl = 1 ; ubl = nch                  ! one landunit per tile
  lbc = 1 ; ubc = nch*nzone            ! number of zones
  lbp = 1 ; ubp = nch*nzone*(numpft+1) ! potential PFT index (0-19)

! initialize CN model
! -------------------
  call clm_varpar_init()
  call initClmtype(lbg,ubg,lbl,ubl,lbc,ubc,lbp,ubp) ! allocation & initialization

! initialize PFT parameters
! -------------------------
  pftcon%z0mr         = z0mx	     ! ratio of momentum roughness length to canopy top height (-)
  pftcon%displar      = displax	     ! ratio of displacement height to canopy top height (-)
  pftcon%c3psn        = c3psx        ! photosynthetic pathway: 0. = c4, 1. = c3
  pftcon%vcmx25       = vcmx2x	     ! max rate of carboxylation at 25C (umol CO2/m**2/s)
  pftcon%mp           = mx	     ! slope of conductance-to-photosynthesis relationship
  pftcon%qe25         = qe2x	     ! quantum efficiency at 25C (umol CO2 / umol photon)
  pftcon%slatop       = slatox	     ! specific leaf area at top of canopy, projected area basis [m^2/gC]
  pftcon%dsladlai     = dsladlax     ! dSLA/dLAI, projected area basis [m^2/gC]
  pftcon%leafcn       = leafcx	     ! leaf C:N (gC/gN)
  pftcon%flnr         = flnx	     ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
  pftcon%fnitr        = fnitx	     ! foliage nitrogen limitation factor (-)
  pftcon%woody        = woodx	     ! binary flag for woody lifeform (1=woody, 0=not woody)
  pftcon%lflitcn      = lflitcx	     ! leaf litter C:N (gC/gN)
  pftcon%frootcn      = frootcx	     ! fine root C:N (gC/gN)
  pftcon%livewdcn     = livewdcx     ! live wood (phloem and ray parenchyma) C:N (gC/gN)
  pftcon%deadwdcn     = deadwdcx     ! dead wood (xylem and heartwood) C:N (gC/gN)
  pftcon%dwood        = dwoox  	     ! wood density (gC/m3)
  pftcon%froot_leaf   = froot_leax   ! allocation parameter: new fine root C per new leaf C (gC/gC)
  pftcon%stem_leaf    = stem_leax    ! allocation parameter: new stem c per new leaf C (gC/gC)
  pftcon%croot_stem   = croot_stex   ! allocation parameter: new coarse root C per new stem C (gC/gC)
  pftcon%flivewd      = flivewx	     ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  pftcon%fcur         = fcux	     ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  pftcon%lf_flab      = lf_flax	     ! leaf litter labile fraction
  pftcon%lf_fcel      = lf_fcex	     ! leaf litter cellulose fraction
  pftcon%lf_flig      = lf_flix	     ! leaf litter lignin fraction
  pftcon%fr_flab      = fr_flax	     ! fine root litter labile fraction
  pftcon%fr_fcel      = fr_fcex	     ! fine root litter cellulose fraction
  pftcon%fr_flig      = fr_flix	     ! fine root litter lignin fraction
  pftcon%leaf_long    = leaf_lonx    ! leaf longevity (yrs)
  pftcon%evergreen    = evergreex    ! binary flag for evergreen leaf habit (0 or 1)
  pftcon%resist       = resisx	     ! resistance to fire (no units)
  pftcon%stress_decid = stress_decix ! binary flag for stress-deciduous leaf habit (0 or 1)
  pftcon%season_decid = season_decix ! binary flag for seasonal-deciduous leaf habit (0 or 1)
  pftcon%xl           = xlx          ! leaf/stem orientation index
  pftcon%rhol         = rholx        ! leaf reflectance (visible)
  pftcon%rhos         = rhosx        ! stem reflectance (visible)
  pftcon%taul         = taulx        ! leaf transmittance (visible)
  pftcon%taus         = tausx        ! stem transmittance (visible)

! transfer restart vars from to CLM data structures if restart exists
! -------------------------------------------------------------------
  if(istep /= 0) then

  n = 0
  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,nzone    ! CN zone loop
      n = n + 1
      clm3%g%l%c%ccs%col_ctrunc      (n) = cncol(nc,nz, 1)
      clm3%g%l%c%ccs%cwdc            (n) = cncol(nc,nz, 2)
      clm3%g%l%c%ccs%litr1c          (n) = cncol(nc,nz, 3)
      clm3%g%l%c%ccs%litr2c          (n) = cncol(nc,nz, 4)
      clm3%g%l%c%ccs%litr3c          (n) = cncol(nc,nz, 5)
      clm3%g%l%c%ccs%pcs_a%totvegc   (n) = cncol(nc,nz, 6)
      clm3%g%l%c%ccs%prod100c        (n) = cncol(nc,nz, 7)
      clm3%g%l%c%ccs%prod10c         (n) = cncol(nc,nz, 8)
      clm3%g%l%c%ccs%seedc           (n) = cncol(nc,nz, 9)
      clm3%g%l%c%ccs%soil1c          (n) = cncol(nc,nz,10)
      clm3%g%l%c%ccs%soil2c          (n) = cncol(nc,nz,11)
      clm3%g%l%c%ccs%soil3c          (n) = cncol(nc,nz,12)
      clm3%g%l%c%ccs%soil4c          (n) = cncol(nc,nz,13)
      clm3%g%l%c%ccs%totcolc         (n) = cncol(nc,nz,14)
      clm3%g%l%c%ccs%totlitc         (n) = cncol(nc,nz,15)
      clm3%g%l%c%cns%col_ntrunc      (n) = cncol(nc,nz,16)
      clm3%g%l%c%cns%cwdn            (n) = cncol(nc,nz,17)
      clm3%g%l%c%cns%litr1n          (n) = cncol(nc,nz,18)
      clm3%g%l%c%cns%litr2n          (n) = cncol(nc,nz,19)
      clm3%g%l%c%cns%litr3n          (n) = cncol(nc,nz,20)
      clm3%g%l%c%cns%prod100n        (n) = cncol(nc,nz,21)
      clm3%g%l%c%cns%prod10n         (n) = cncol(nc,nz,22)
      clm3%g%l%c%cns%seedn           (n) = cncol(nc,nz,23)
      clm3%g%l%c%cns%sminn           (n) = cncol(nc,nz,24)
      clm3%g%l%c%cns%soil1n          (n) = cncol(nc,nz,25)
      clm3%g%l%c%cns%soil2n          (n) = cncol(nc,nz,26)
      clm3%g%l%c%cns%soil3n          (n) = cncol(nc,nz,27)
      clm3%g%l%c%cns%soil4n          (n) = cncol(nc,nz,28)
      clm3%g%l%c%cns%totcoln         (n) = cncol(nc,nz,29)
      clm3%g%l%c%cps%ann_farea_burned(n) = cncol(nc,nz,30)
      clm3%g%l%c%cps%annsum_counter  (n) = cncol(nc,nz,31)
      clm3%g%l%c%cps%cannavg_t2m     (n) = cncol(nc,nz,32)
      clm3%g%l%c%cps%cannsum_npp     (n) = cncol(nc,nz,33)
      clm3%g%l%c%cps%farea_burned    (n) = cncol(nc,nz,34)
      clm3%g%l%c%cps%fire_prob       (n) = cncol(nc,nz,35)
      clm3%g%l%c%cps%fireseasonl     (n) = cncol(nc,nz,36)
      clm3%g%l%c%cps%fpg             (n) = cncol(nc,nz,37)
      clm3%g%l%c%cps%fpi             (n) = cncol(nc,nz,38)
      clm3%g%l%c%cps%me              (n) = cncol(nc,nz,39)
      clm3%g%l%c%cps%mean_fire_prob  (n) = cncol(nc,nz,40)

      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,nveg ! defined veg loop

          if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
           clm3%g%l%c%p%pcs%cpool                 (np) = cnpft(nc,nz,nv, 1)
           clm3%g%l%c%p%pcs%deadcrootc            (np) = cnpft(nc,nz,nv, 2)
           clm3%g%l%c%p%pcs%deadcrootc_storage    (np) = cnpft(nc,nz,nv, 3)
           clm3%g%l%c%p%pcs%deadcrootc_xfer       (np) = cnpft(nc,nz,nv, 4)
           clm3%g%l%c%p%pcs%deadstemc             (np) = cnpft(nc,nz,nv, 5)
           clm3%g%l%c%p%pcs%deadstemc_storage     (np) = cnpft(nc,nz,nv, 6)
           clm3%g%l%c%p%pcs%deadstemc_xfer        (np) = cnpft(nc,nz,nv, 7)
           clm3%g%l%c%p%pcs%frootc                (np) = cnpft(nc,nz,nv, 8)
           clm3%g%l%c%p%pcs%frootc_storage        (np) = cnpft(nc,nz,nv, 9)
           clm3%g%l%c%p%pcs%frootc_xfer           (np) = cnpft(nc,nz,nv,10)
           clm3%g%l%c%p%pcs%gresp_storage         (np) = cnpft(nc,nz,nv,11)
           clm3%g%l%c%p%pcs%gresp_xfer            (np) = cnpft(nc,nz,nv,12)
           clm3%g%l%c%p%pcs%leafc                 (np) = cnpft(nc,nz,nv,13)
           clm3%g%l%c%p%pcs%leafc_storage         (np) = cnpft(nc,nz,nv,14) 
           clm3%g%l%c%p%pcs%leafc_xfer            (np) = cnpft(nc,nz,nv,15)
           clm3%g%l%c%p%pcs%livecrootc            (np) = cnpft(nc,nz,nv,16) 
           clm3%g%l%c%p%pcs%livecrootc_storage    (np) = cnpft(nc,nz,nv,17)
           clm3%g%l%c%p%pcs%livecrootc_xfer       (np) = cnpft(nc,nz,nv,18)
           clm3%g%l%c%p%pcs%livestemc             (np) = cnpft(nc,nz,nv,19)
           clm3%g%l%c%p%pcs%livestemc_storage     (np) = cnpft(nc,nz,nv,20)
           clm3%g%l%c%p%pcs%livestemc_xfer        (np) = cnpft(nc,nz,nv,21)
           clm3%g%l%c%p%pcs%pft_ctrunc            (np) = cnpft(nc,nz,nv,22)
           clm3%g%l%c%p%pcs%xsmrpool              (np) = cnpft(nc,nz,nv,23)
           clm3%g%l%c%p%pepv%annavg_t2m           (np) = cnpft(nc,nz,nv,24)
           clm3%g%l%c%p%pepv%annmax_retransn      (np) = cnpft(nc,nz,nv,25)
           clm3%g%l%c%p%pepv%annsum_npp           (np) = cnpft(nc,nz,nv,26)
           clm3%g%l%c%p%pepv%annsum_potential_gpp (np) = cnpft(nc,nz,nv,27)
           clm3%g%l%c%p%pepv%dayl                 (np) = cnpft(nc,nz,nv,28)
           clm3%g%l%c%p%pepv%days_active          (np) = cnpft(nc,nz,nv,29)
           clm3%g%l%c%p%pepv%dormant_flag         (np) = cnpft(nc,nz,nv,30)
           clm3%g%l%c%p%pepv%offset_counter       (np) = cnpft(nc,nz,nv,31)
           clm3%g%l%c%p%pepv%offset_fdd           (np) = cnpft(nc,nz,nv,32)
           clm3%g%l%c%p%pepv%offset_flag          (np) = cnpft(nc,nz,nv,33)
           clm3%g%l%c%p%pepv%offset_swi           (np) = cnpft(nc,nz,nv,34)
           clm3%g%l%c%p%pepv%onset_counter        (np) = cnpft(nc,nz,nv,35)
           clm3%g%l%c%p%pepv%onset_fdd            (np) = cnpft(nc,nz,nv,36)
           clm3%g%l%c%p%pepv%onset_flag           (np) = cnpft(nc,nz,nv,37)
           clm3%g%l%c%p%pepv%onset_gdd            (np) = cnpft(nc,nz,nv,38)
           clm3%g%l%c%p%pepv%onset_gddflag        (np) = cnpft(nc,nz,nv,39)
           clm3%g%l%c%p%pepv%onset_swi            (np) = cnpft(nc,nz,nv,40)
           clm3%g%l%c%p%pepv%prev_frootc_to_litter(np) = cnpft(nc,nz,nv,41)
           clm3%g%l%c%p%pepv%prev_leafc_to_litter (np) = cnpft(nc,nz,nv,42)
           clm3%g%l%c%p%pepv%tempavg_t2m          (np) = cnpft(nc,nz,nv,43)
           clm3%g%l%c%p%pepv%tempmax_retransn     (np) = cnpft(nc,nz,nv,44)
           clm3%g%l%c%p%pepv%tempsum_npp          (np) = cnpft(nc,nz,nv,45)
           clm3%g%l%c%p%pepv%tempsum_potential_gpp(np) = cnpft(nc,nz,nv,46)
           clm3%g%l%c%p%pepv%xsmrpool_recover     (np) = cnpft(nc,nz,nv,47)
           clm3%g%l%c%p%pns%deadcrootn            (np) = cnpft(nc,nz,nv,48)
           clm3%g%l%c%p%pns%deadcrootn_storage    (np) = cnpft(nc,nz,nv,49)
           clm3%g%l%c%p%pns%deadcrootn_xfer       (np) = cnpft(nc,nz,nv,50)
           clm3%g%l%c%p%pns%deadstemn             (np) = cnpft(nc,nz,nv,51)
           clm3%g%l%c%p%pns%deadstemn_storage     (np) = cnpft(nc,nz,nv,52)
           clm3%g%l%c%p%pns%deadstemn_xfer        (np) = cnpft(nc,nz,nv,53)
           clm3%g%l%c%p%pns%frootn                (np) = cnpft(nc,nz,nv,54)
           clm3%g%l%c%p%pns%frootn_storage        (np) = cnpft(nc,nz,nv,55)
           clm3%g%l%c%p%pns%frootn_xfer           (np) = cnpft(nc,nz,nv,56)
           clm3%g%l%c%p%pns%leafn                 (np) = cnpft(nc,nz,nv,57)
           clm3%g%l%c%p%pns%leafn_storage         (np) = cnpft(nc,nz,nv,58)
           clm3%g%l%c%p%pns%leafn_xfer            (np) = cnpft(nc,nz,nv,59)
           clm3%g%l%c%p%pns%livecrootn            (np) = cnpft(nc,nz,nv,60)
           clm3%g%l%c%p%pns%livecrootn_storage    (np) = cnpft(nc,nz,nv,61)
           clm3%g%l%c%p%pns%livecrootn_xfer       (np) = cnpft(nc,nz,nv,62)
           clm3%g%l%c%p%pns%livestemn             (np) = cnpft(nc,nz,nv,63)
           clm3%g%l%c%p%pns%livestemn_storage     (np) = cnpft(nc,nz,nv,64)
           clm3%g%l%c%p%pns%livestemn_xfer        (np) = cnpft(nc,nz,nv,65)
           clm3%g%l%c%p%pns%npool                 (np) = cnpft(nc,nz,nv,66)
           clm3%g%l%c%p%pns%pft_ntrunc            (np) = cnpft(nc,nz,nv,67)
           clm3%g%l%c%p%pns%retransn              (np) = cnpft(nc,nz,nv,68)
           clm3%g%l%c%p%pps%elai                  (np) = cnpft(nc,nz,nv,69)
           clm3%g%l%c%p%pps%esai                  (np) = cnpft(nc,nz,nv,70)
           clm3%g%l%c%p%pps%hbot                  (np) = cnpft(nc,nz,nv,71)
           clm3%g%l%c%p%pps%htop                  (np) = cnpft(nc,nz,nv,72)
           clm3%g%l%c%p%pps%tlai                  (np) = cnpft(nc,nz,nv,73)
           clm3%g%l%c%p%pps%tsai                  (np) = cnpft(nc,nz,nv,74)
          endif

        end do ! defined veg loop
      end do   ! PFT index loop

    end do     ! CN zone loop
  end do       ! catchment tile loop

  endif

  return

  end subroutine CN_init

  subroutine CN_exit(nch,nveg,nzone,ityp,fveg,cncol,var_col,cnpft,var_pft)     

  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg    ! PFT fraction

  integer, intent(in) :: var_col ! number of CN column restart variables
  real*4, dimension(nch,nzone,var_col), intent(out) :: cncol ! gkw: column CN restart 

  integer, intent(in) :: var_pft ! number of CN PFT restart variables
  real*4, dimension(nch,nzone,nveg,var_pft), intent(out) :: cnpft ! gkw: PFT CN restart 

  integer :: n, p, nv, nc, nz, np


! copy CN_restart vars to catch_internal_rst
! ------------------------------------------
  n = 0
  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,nzone    ! CN zone loop
      n = n + 1
      cncol(nc,nz, 1) = clm3%g%l%c%ccs%col_ctrunc      (n)
      cncol(nc,nz, 2) = clm3%g%l%c%ccs%cwdc            (n)
      cncol(nc,nz, 3) = clm3%g%l%c%ccs%litr1c          (n)
      cncol(nc,nz, 4) = clm3%g%l%c%ccs%litr2c          (n)
      cncol(nc,nz, 5) = clm3%g%l%c%ccs%litr3c          (n)
      cncol(nc,nz, 6) = clm3%g%l%c%ccs%pcs_a%totvegc   (n)
      cncol(nc,nz, 7) = clm3%g%l%c%ccs%prod100c        (n)
      cncol(nc,nz, 8) = clm3%g%l%c%ccs%prod10c         (n)
      cncol(nc,nz, 9) = clm3%g%l%c%ccs%seedc           (n)
      cncol(nc,nz,10) = clm3%g%l%c%ccs%soil1c          (n)
      cncol(nc,nz,11) = clm3%g%l%c%ccs%soil2c          (n)
      cncol(nc,nz,12) = clm3%g%l%c%ccs%soil3c          (n)
      cncol(nc,nz,13) = clm3%g%l%c%ccs%soil4c          (n)
      cncol(nc,nz,14) = clm3%g%l%c%ccs%totcolc         (n)
      cncol(nc,nz,15) = clm3%g%l%c%ccs%totlitc         (n)
      cncol(nc,nz,16) = clm3%g%l%c%cns%col_ntrunc      (n)
      cncol(nc,nz,17) = clm3%g%l%c%cns%cwdn            (n)
      cncol(nc,nz,18) = clm3%g%l%c%cns%litr1n          (n)
      cncol(nc,nz,19) = clm3%g%l%c%cns%litr2n          (n)
      cncol(nc,nz,20) = clm3%g%l%c%cns%litr3n          (n)
      cncol(nc,nz,21) = clm3%g%l%c%cns%prod100n        (n)
      cncol(nc,nz,22) = clm3%g%l%c%cns%prod10n         (n)
      cncol(nc,nz,23) = clm3%g%l%c%cns%seedn           (n)
      cncol(nc,nz,24) = clm3%g%l%c%cns%sminn           (n)
      cncol(nc,nz,25) = clm3%g%l%c%cns%soil1n          (n)
      cncol(nc,nz,26) = clm3%g%l%c%cns%soil2n          (n)
      cncol(nc,nz,27) = clm3%g%l%c%cns%soil3n          (n)
      cncol(nc,nz,28) = clm3%g%l%c%cns%soil4n          (n)
      cncol(nc,nz,29) = clm3%g%l%c%cns%totcoln         (n)
      cncol(nc,nz,30) = clm3%g%l%c%cps%ann_farea_burned(n)
      cncol(nc,nz,31) = clm3%g%l%c%cps%annsum_counter  (n)
      cncol(nc,nz,32) = clm3%g%l%c%cps%cannavg_t2m     (n)
      cncol(nc,nz,33) = clm3%g%l%c%cps%cannsum_npp     (n)
      cncol(nc,nz,34) = clm3%g%l%c%cps%farea_burned    (n)
      cncol(nc,nz,35) = clm3%g%l%c%cps%fire_prob       (n)
      cncol(nc,nz,36) = clm3%g%l%c%cps%fireseasonl     (n)
      cncol(nc,nz,37) = clm3%g%l%c%cps%fpg             (n)
      cncol(nc,nz,38) = clm3%g%l%c%cps%fpi             (n)
      cncol(nc,nz,39) = clm3%g%l%c%cps%me              (n)
      cncol(nc,nz,40) = clm3%g%l%c%cps%mean_fire_prob  (n)

      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,nveg ! defined veg loop

          if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
            cnpft(nc,nz,nv, 1) = clm3%g%l%c%p%pcs%cpool                 (np)
            cnpft(nc,nz,nv, 2) = clm3%g%l%c%p%pcs%deadcrootc            (np)
            cnpft(nc,nz,nv, 3) = clm3%g%l%c%p%pcs%deadcrootc_storage    (np)
            cnpft(nc,nz,nv, 4) = clm3%g%l%c%p%pcs%deadcrootc_xfer       (np)
            cnpft(nc,nz,nv, 5) = clm3%g%l%c%p%pcs%deadstemc             (np)
            cnpft(nc,nz,nv, 6) = clm3%g%l%c%p%pcs%deadstemc_storage     (np)
            cnpft(nc,nz,nv, 7) = clm3%g%l%c%p%pcs%deadstemc_xfer        (np)
            cnpft(nc,nz,nv, 8) = clm3%g%l%c%p%pcs%frootc                (np)
            cnpft(nc,nz,nv, 9) = clm3%g%l%c%p%pcs%frootc_storage        (np)
            cnpft(nc,nz,nv,10) = clm3%g%l%c%p%pcs%frootc_xfer           (np)
            cnpft(nc,nz,nv,11) = clm3%g%l%c%p%pcs%gresp_storage         (np)
            cnpft(nc,nz,nv,12) = clm3%g%l%c%p%pcs%gresp_xfer            (np)
            cnpft(nc,nz,nv,13) = clm3%g%l%c%p%pcs%leafc                 (np)
            cnpft(nc,nz,nv,14) = clm3%g%l%c%p%pcs%leafc_storage         (np) 
            cnpft(nc,nz,nv,15) = clm3%g%l%c%p%pcs%leafc_xfer            (np)
            cnpft(nc,nz,nv,16) = clm3%g%l%c%p%pcs%livecrootc            (np) 
            cnpft(nc,nz,nv,17) = clm3%g%l%c%p%pcs%livecrootc_storage    (np)
            cnpft(nc,nz,nv,18) = clm3%g%l%c%p%pcs%livecrootc_xfer       (np)
            cnpft(nc,nz,nv,19) = clm3%g%l%c%p%pcs%livestemc             (np)
            cnpft(nc,nz,nv,20) = clm3%g%l%c%p%pcs%livestemc_storage     (np)
            cnpft(nc,nz,nv,21) = clm3%g%l%c%p%pcs%livestemc_xfer        (np)
            cnpft(nc,nz,nv,22) = clm3%g%l%c%p%pcs%pft_ctrunc            (np)
            cnpft(nc,nz,nv,23) = clm3%g%l%c%p%pcs%xsmrpool              (np)
            cnpft(nc,nz,nv,24) = clm3%g%l%c%p%pepv%annavg_t2m           (np)
            cnpft(nc,nz,nv,25) = clm3%g%l%c%p%pepv%annmax_retransn      (np)
            cnpft(nc,nz,nv,26) = clm3%g%l%c%p%pepv%annsum_npp           (np)
            cnpft(nc,nz,nv,27) = clm3%g%l%c%p%pepv%annsum_potential_gpp (np)
            cnpft(nc,nz,nv,28) = clm3%g%l%c%p%pepv%dayl                 (np)
            cnpft(nc,nz,nv,29) = clm3%g%l%c%p%pepv%days_active          (np)
            cnpft(nc,nz,nv,30) = clm3%g%l%c%p%pepv%dormant_flag         (np)
            cnpft(nc,nz,nv,31) = clm3%g%l%c%p%pepv%offset_counter       (np)
            cnpft(nc,nz,nv,32) = clm3%g%l%c%p%pepv%offset_fdd           (np)
            cnpft(nc,nz,nv,33) = clm3%g%l%c%p%pepv%offset_flag          (np)
            cnpft(nc,nz,nv,34) = clm3%g%l%c%p%pepv%offset_swi           (np)
            cnpft(nc,nz,nv,35) = clm3%g%l%c%p%pepv%onset_counter        (np)
            cnpft(nc,nz,nv,36) = clm3%g%l%c%p%pepv%onset_fdd            (np)
            cnpft(nc,nz,nv,37) = clm3%g%l%c%p%pepv%onset_flag           (np)
            cnpft(nc,nz,nv,38) = clm3%g%l%c%p%pepv%onset_gdd            (np)
            cnpft(nc,nz,nv,39) = clm3%g%l%c%p%pepv%onset_gddflag        (np)
            cnpft(nc,nz,nv,40) = clm3%g%l%c%p%pepv%onset_swi            (np)
            cnpft(nc,nz,nv,41) = clm3%g%l%c%p%pepv%prev_frootc_to_litter(np)
            cnpft(nc,nz,nv,42) = clm3%g%l%c%p%pepv%prev_leafc_to_litter (np)
            cnpft(nc,nz,nv,43) = clm3%g%l%c%p%pepv%tempavg_t2m          (np)
            cnpft(nc,nz,nv,44) = clm3%g%l%c%p%pepv%tempmax_retransn     (np)
            cnpft(nc,nz,nv,45) = clm3%g%l%c%p%pepv%tempsum_npp          (np)
            cnpft(nc,nz,nv,46) = clm3%g%l%c%p%pepv%tempsum_potential_gpp(np)
            cnpft(nc,nz,nv,47) = clm3%g%l%c%p%pepv%xsmrpool_recover     (np)
            cnpft(nc,nz,nv,48) = clm3%g%l%c%p%pns%deadcrootn            (np)
            cnpft(nc,nz,nv,49) = clm3%g%l%c%p%pns%deadcrootn_storage    (np)
            cnpft(nc,nz,nv,50) = clm3%g%l%c%p%pns%deadcrootn_xfer       (np)
            cnpft(nc,nz,nv,51) = clm3%g%l%c%p%pns%deadstemn             (np)
            cnpft(nc,nz,nv,52) = clm3%g%l%c%p%pns%deadstemn_storage     (np)
            cnpft(nc,nz,nv,53) = clm3%g%l%c%p%pns%deadstemn_xfer        (np)
            cnpft(nc,nz,nv,54) = clm3%g%l%c%p%pns%frootn                (np)
            cnpft(nc,nz,nv,55) = clm3%g%l%c%p%pns%frootn_storage        (np)
            cnpft(nc,nz,nv,56) = clm3%g%l%c%p%pns%frootn_xfer           (np)
            cnpft(nc,nz,nv,57) = clm3%g%l%c%p%pns%leafn                 (np)
            cnpft(nc,nz,nv,58) = clm3%g%l%c%p%pns%leafn_storage         (np)
            cnpft(nc,nz,nv,59) = clm3%g%l%c%p%pns%leafn_xfer            (np)
            cnpft(nc,nz,nv,60) = clm3%g%l%c%p%pns%livecrootn            (np)
            cnpft(nc,nz,nv,61) = clm3%g%l%c%p%pns%livecrootn_storage    (np)
            cnpft(nc,nz,nv,62) = clm3%g%l%c%p%pns%livecrootn_xfer       (np)
            cnpft(nc,nz,nv,63) = clm3%g%l%c%p%pns%livestemn             (np)
            cnpft(nc,nz,nv,64) = clm3%g%l%c%p%pns%livestemn_storage     (np)
            cnpft(nc,nz,nv,65) = clm3%g%l%c%p%pns%livestemn_xfer        (np)
            cnpft(nc,nz,nv,66) = clm3%g%l%c%p%pns%npool                 (np)
            cnpft(nc,nz,nv,67) = clm3%g%l%c%p%pns%pft_ntrunc            (np)
            cnpft(nc,nz,nv,68) = clm3%g%l%c%p%pns%retransn              (np)
            cnpft(nc,nz,nv,69) = clm3%g%l%c%p%pps%elai                  (np)
            cnpft(nc,nz,nv,70) = clm3%g%l%c%p%pps%esai                  (np)
            cnpft(nc,nz,nv,71) = clm3%g%l%c%p%pps%hbot                  (np)
            cnpft(nc,nz,nv,72) = clm3%g%l%c%p%pps%htop                  (np)
            cnpft(nc,nz,nv,73) = clm3%g%l%c%p%pps%tlai                  (np)
            cnpft(nc,nz,nv,74) = clm3%g%l%c%p%pps%tsai                  (np)
          endif

        end do ! defined veg loop
      end do   ! PFT index loop
    end do     ! CN zone loop
  end do       ! catchment tile loop

  return

  end subroutine CN_exit

  subroutine get_CN_LAI(nch,nveg,nzone,ityp,fveg,elai,esai)

  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg    ! PFT fraction
  real, dimension(nch,nveg,nzone), intent(out) :: elai, esai ! LAI & SAI

  integer :: n, p, nv, nc, nz, np

  elai = 0.
  esai = 0.

  n = 0
  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,nzone    ! CN zone loop
      n = n + 1
      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,nveg ! defined veg loop

! extract LAI & SAI from CN clmtype
! ---------------------------------
          if(ityp(nc,nv,nz)==p .and. ityp(nc,nv,nz)>0 .and. fveg(nc,nv,nz)>1.e-4) then
            elai(nc,nv,nz) = clm3%g%l%c%p%pps%elai(np)
            esai(nc,nv,nz) = clm3%g%l%c%p%pps%esai(np)
          endif

        end do ! defined veg loop
      end do   ! PFT index loop
    end do     ! CN zone loop
  end do       ! catchment tile loop

  end subroutine get_CN_LAI

end module CN_DriverMod
