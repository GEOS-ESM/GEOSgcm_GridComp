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
  use clm_varcon,        only: clm_varcon_init
  use clm_varpar,        only: clm_varpar_init, numpft
  use CNSetValueMod,     only: CNZeroFluxes_dwt
  use CNAnnualUpdateMod, only: CNAnnualUpdate
  use CNBalanceCheckMod, only: BeginCBalance, BeginNBalance, &
                               CBalanceCheck, NBalanceCheck
  use pftvarcon,         only: noveg
  use clm_time_manager,  only: get_step_size, get_nstep, get_curr_date
  use shr_const_mod,     only: SHR_CONST_TKFRZ, SHR_CONST_CDAY
  use pftvarcon,         only: npcropmin, ntree, nbrdlf_dcd_trp_tree
#ifndef CENTURY_DECOMP
  use CNDecompCascadeMod_BGC, only : init_decompcascade
#else
  use CNDecompCascadeMod_CENTURY, only : init_decompcascade
#endif
  use catch_constants,   only: DZTSURF=>CATCH_DZTSURF, DZGT=>CATCH_DZGT
  use SurfParams,        only: LAND_FIX
!  use update_model_para4cn, only : LocalTileID, upd_tileid   ! useful for debugging

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

  subroutine CN_Driver(istep,nch,nveg,nzone,daylength,                        &
                       tgw,tp1,tp2,tp3,tp4,tp5,tp6,sfm,rzm,                   &
		       wpwet,psis,bee,poros,vgwmax,bflow,totwat,              &
		       runsrf,tm,rhm,windm,rainfm,snowfm,                     &
                       prec10d,prec60d,t10d,                                  &
                       psun,psha,lmrsunm,lmrsham,lsun,lsha,                   &
                       car1m,btran2x_rz,btran2x_sf,latitude,longitude,        &
                       ityp,fveg,wtzone,sndzn,asnow,ndep,abm,peatf,           &
                       gdp,hdm,field_cap,lnfm,zlai,zsai,ztai,colc,tile_id,    &
		       ann_t2m,nppg,gppg,srg,neeg,root,padd,vegc,xsmr,        &
                       burn,closs,nfire,som_closs,                            &
                       ndeployg,denitg,sminn_leachedg,sminng,col_fire_nlossg, &
                       leafng,leafcg,gross_nming,net_nming,nfix_to_sminng,    &
                       actual_immobg,fpgg,fpig,sminn_to_plantg,               &
                       sminn_to_npoolg,ndep_to_sminng,totvegng,totlitng,      &
                       totsomng,retransng,retransn_to_npoolg,fuelcg,totlitcg, &
                       cwdcg,rootcg)

  ! !ARGUMENTS:
  implicit none

  integer*8, intent(in) :: istep ! number of CN time steps run
  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  real*4, dimension(nch), intent(in) :: daylength ! daylength (seconds)
  real*4, dimension(nch,nzone), intent(in) :: tgw ! soil surface layer temperature (K)
  real*4, dimension(nch,nzone), intent(in) :: rzm ! weighted root-zone moisture content as frac of WHC
  real*4, dimension(nch,nzone), intent(in) :: sfm ! weighted surface moisture content as frac of WHC
  real*4, dimension(nch), intent(in) :: tp1,tp2,tp3,tp4,tp5,tp6 ! soil temperatures (K)
  real*4, dimension(nch), intent(in) :: wpwet,psis,bee,poros,vgwmax,bflow,totwat ! soil water & parameters
  real*4, dimension(nch), intent(in) :: tm        ! air temperature (K)
  real*4, dimension(nch), intent(in) :: rhm       ! relative humidity (%)
  real*4, dimension(nch), intent(in) :: windm     ! wind speed (m/s)
  real*4, dimension(nch), intent(in) :: rainfm    ! rainfall (convective + largescale) (kg/m2/s)
  real*4, dimension(nch), intent(in) :: snowfm    ! snowfall (kg/m2/s)  
  real*4, dimension(nch), intent(in) :: prec10d   ! 10-day running mean of total precipitation (mm H2O/s)
  real*4, dimension(nch), intent(in) :: prec60d   ! 60-day running mean of total precipitation (mm H2O/s)
  real*4, dimension(nch), intent(in) :: t10d      ! 10-day running mean of 2-m temperature (K)
  real*4, dimension(nch), intent(in) :: runsrf    ! surface runoff (kg/m2/s)
  real*4, dimension(nch), intent(in) :: ndep      ! nitrogen deposition
  real*4, dimension(nch), intent(in) :: abm       ! Peak month for agricultural fire, unitless
  real*4, dimension(nch), intent(in) :: peatf     ! Fraction of peatland, unitless (0-1)
  real*4, dimension(nch), intent(in) :: gdp       ! Real GDP (K 1995US$/capita)
  real*4, dimension(nch), intent(in) :: hdm       ! Human population density in 2010 (individual/km2)
  real*4, dimension(nch), intent(in) :: field_cap ! Field capacity (m3/m3)
  real*4, dimension(nch), intent(in) :: lnfm      ! Lightning frequency [Flashes/km^2/day]
  real*4, dimension(nch), intent(in) :: sndzn     ! snow height of snow covered area (m)
  real*4, dimension(nch), intent(in) :: asnow     ! areal snow coverage [0-1]
  real*4, dimension(nch), intent(in) :: car1m     ! fraction of tile that is saturated area
  real*4, dimension(nch,nzone), intent(in) :: btran2x_rz ! root zone soil wetness, used to calculate btran2 for CNFireMod
  real*4, dimension(nch,nzone), intent(in) :: btran2x_sf ! surface soil wetness, used to calculate btran2 for CNFireMod
  real*4, dimension(nch), intent(in) :: latitude  ! center-of-mass latitude (degree)
  real*4, dimension(nch), intent(in) :: longitude ! center-of-mass longitude (degree)  
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! CLM PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg ! catchment vegetation fractions
  real, dimension(nch,nzone), intent(in) :: wtzone ! zone fractions
  real*4, dimension(nch,nveg,nzone), intent(in) :: psha,psun ! photosynthesis
  real*4, dimension(nch,nveg,nzone), intent(in) :: lmrsunm,lmrsham ! leaf maintenance respiration rate (umol CO2/m**2/s)
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
  real*4, dimension(nch), intent(out) :: burn, closs, nfire, som_closs

! column & PFT level carbon added to sustain growth
! -------------------------------------------------
  real*4, dimension(nch), intent(out) :: padd  ! (gC/m2/s)
  real*4, dimension(nch), intent(out) :: root, vegc, xsmr  ! (gC/m2)
  
! states and fluxes to understand nitrogen cycle, fzeng, 4 April 2019
  real*4, dimension(nch), intent(out) :: ndeployg        ! total N deployed to growth and storage (gN/m2/s)
  real*4, dimension(nch), intent(out) :: denitg          ! total rate of denitrification (gN/m2/s)
  real*4, dimension(nch), intent(out) :: sminn_leachedg  ! soil mineral N pool loss to leaching (gN/m2/s)
  real*4, dimension(nch), intent(out) :: sminng          ! (gN/m2) soil mineral N
  real*4, dimension(nch), intent(out) :: col_fire_nlossg ! (gN/m2/s) total column-level fire N loss
  real*4, dimension(nch), intent(out) :: leafng          ! (gN/m2) leaf N
  real*4, dimension(nch), intent(out) :: leafcg          ! (gC/m2) leaf C
  real*4, dimension(nch), intent(out) :: gross_nming     ! gross rate of N mineralization (gN/m2/s)
  real*4, dimension(nch), intent(out) :: net_nming       ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
  real*4, dimension(nch), intent(out) :: nfix_to_sminng  ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
  real*4, dimension(nch), intent(out) :: actual_immobg   ! vert-int (diagnostic) actual N immobilization (gN/m2/s) 
  real*4, dimension(nch), intent(out) :: fpgg            ! fraction of potential gpp (no units)
  real*4, dimension(nch), intent(out) :: fpig            ! fraction of potential immobilization (no units)
  real*4, dimension(nch), intent(out) :: sminn_to_plantg ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
  real*4, dimension(nch), intent(out) :: sminn_to_npoolg ! deployment of soil mineral N uptake (gN/m2/s)
  real*4, dimension(nch), intent(out) :: ndep_to_sminng  ! atmospheric N deposition to soil mineral N (gN/m2/s)
  real*4, dimension(nch), intent(out) :: totvegng              ! (gN/m2) total vegetation nitrogen
  real*4, dimension(nch), intent(out) :: totlitng              ! (gN/m2) total litter nitrogen
  real*4, dimension(nch), intent(out) :: totsomng              ! (gN/m2) total soil organic matter nitrogen
  real*4, dimension(nch), intent(out) :: retransng             ! (gN/m2) plant pool of retranslocated N
  real*4, dimension(nch), intent(out) :: retransn_to_npoolg    ! deployment of retranslocated N (gN/m2/s)

! states and fluxes to understand the fire model, fzeng, 30 July 2019
  real*4, dimension(nch), intent(out) :: fuelcg          ! fuel avalability for non-crop areas outside tropical closed broadleaf evergreen closed forests (gC/m2)
  real*4, dimension(nch), intent(out) :: totlitcg        ! (gC/m2) total litter carbon
  real*4, dimension(nch), intent(out) :: cwdcg            ! (gC/m2) coarse woody debris C
  real*4, dimension(nch), intent(out) :: rootcg          ! (gC/m2) total root carbon

  integer :: lbg, ubg        ! grid bounds
  integer :: lbl, ubl        ! land-type bounds
  integer :: lbc, ubc        ! column bounds
  integer :: lbp, ubp        ! pft bounds
  integer :: num_soilc       ! number of soil columns in filter
  integer :: num_soilp       ! number of soil pfts in filter
  integer :: num_pcropp      ! number of prog. crop pfts in filter
  
  ! added from /discover/nobackup/fzeng/clm_orig/cesm1_2_2/models/lnd/clm/src/util_share/decompMod.F90, fzeng, 29 Mar 2017
  integer :: numg     ! total number of gridcells on all procs
  integer :: numl     ! total number of landunits on all procs
  integer :: numc     ! total number of columns on all procs
  integer :: nump     ! total number of pfts on all procs
  
  logical, save :: doalb = .true.         ! assume surface albedo calculation time step
  logical, save :: spin = .false.         ! true if spinup vegetation gkw: critical!
  logical, save :: exit_spin = .false.    ! true if this is first continuation from a spin=true run

  logical, save :: first = .true.
  integer, parameter :: npft = numpft+1 
  integer :: n, c, p, pft_num, idum, pf, i, j, nv, nc, nz, z, icn
  real :: bare, leafc_tot
  integer :: itypveg                      ! vegetation type

  integer, allocatable, save :: filter_soilc(:),filter_soilp(:),filter_pcropp(:)
  integer, allocatable, save :: index_soilp(:),zone_soilp(:),zone_soilc(:)
  real, allocatable, save :: leafc_add(:)
  
  integer, allocatable, save :: tileid_soilp(:)     ! fzeng added for debugging

  real :: dt   ! time step delta t (seconds)
  integer*8 :: nstep_cn   ! number of CN model steps run
  integer :: curr_year, curr_mon, curr_day, curr_tod  ! year, month, day, time of day (seconds past 0z) of the current CN time step

  logical, pointer :: ifspecial(:)    !BOOL: true=>landunit is not vegetated
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
  real, pointer    :: cwtgcell(:)     !weight of columns relative to gridcells

  real, pointer :: t_ref2m(:)        !2 m height surface air temperature (Kelvin)
  real, pointer :: psnsun(:)         !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
  real, pointer :: psnsha(:)         !shaded leaf photosynthesis (umol CO2 /m**2/ s)
  real, pointer :: laisun(:)         !sunlit projected leaf area index
  real, pointer :: laisha(:)         !shaded projected leaf area index
  real, pointer :: forc_ndep(:)      !nitrogen deposition rate (gN/m2/s)
  real, pointer :: forc_t(:)         !air temperature (K)
  real, pointer :: forc_rh(:)        !relative humidity (%)
  real, pointer :: forc_wind(:)      !wind speed (m/s)
  real, pointer :: forc_rain(:)      !rainfall (convective + largescale) (mm/s)
  real, pointer :: forc_snow(:)      !snowfall (mm/s)  
  real, pointer :: t_soisno(:,:)     !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
  real, pointer :: rootfr(:,:)       !fraction of roots in each soil layer  (nlevgrnd)
  real, pointer :: dz(:,:)           !layer thickness (m)  (-nlevsno+1:nlevgrnd)
  real, pointer :: psisat(:,:)       !soil water potential at saturation for CN code (MPa)
  real, pointer :: psiwilt(:)        !root-zone soil water potential at wilting point (MPa)
  real, pointer :: soilpsi(:,:)      !soil water potential in each soil layer (MPa)
  real, pointer :: h2osoi_liq(:,:)   !column liquid water (kg/m2) (new)
  real, pointer :: wf(:)             !soil water as frac. of whc for top 0.05 m
  real, pointer :: wf2(:)            !soil water as frac. of whc for top 0.17 m
  real, pointer :: qflx_drain(:)     !sub-surface runoff (mm H2O /s)
  real, pointer :: qflx_surf(:)      !surface runoff (mm H2O /s)
  real, pointer :: t_grnd(:)         !ground temperature (Kelvin)
  real, pointer :: forc_hgt_u_pft(:) !wind forcing height (10m+z0m+d) (m)
  real, pointer :: dayl(:)           !daylength (seconds)
  real, pointer :: prev_dayl(:)      !daylength from previous albedo timestep (seconds)
  real, pointer :: elai(:)           !one-sided leaf area index with burying by snow
  real, pointer :: esai(:)           !one-sided stem area index with burying by snow
  real, pointer :: tlai(:)           !one-sided leaf area index, no burying by snow
  real, pointer :: totcolc(:)        ! (gC/m2) total column carbon, incl veg and cpool
  real, pointer :: annavg_t2m(:)     ! annual average 2m air temperature (K)
  integer, pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]

  real, pointer :: col_ctrunc(:)	 ! (gC/m2) column-level sink for C truncation
  real, pointer :: totlitc(:)		 ! (gC/m2) total litter carbon
  real, pointer :: totsomc(:)		 ! (gC/m2) total soil organic matter carbon
  real, pointer :: leafc(:)
  real, pointer :: leafc_storage(:)
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
  real, pointer :: cwdn(:)		 ! (gN/m2) coarse woody debris N
  real, pointer :: litr1n(:)		 ! (gN/m2) litter labile N
  real, pointer :: litr2n(:)		 ! (gN/m2) litter cellulose N
  real, pointer :: litr3n(:)		 ! (gN/m2) litter lignin N
  real, pointer :: soil1n(:)		 ! (gN/m2) soil organic matter N (fast pool)
  real, pointer :: soil2n(:)		 ! (gN/m2) soil organic matter N (medium pool)
  real, pointer :: soil3n(:)		 ! (gN/m2) soil orgainc matter N (slow pool)
  real, pointer :: soil4n(:)		 ! (gN/m2) soil orgainc matter N (slowest pool)
  real, pointer :: sminn_vr(:,:)	 ! (gN/m2) soil mineral N
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
  real, pointer :: seedc(:)		 ! (gC/m2) column-level pool for seeding new PFTs
  real, pointer :: prod10c(:)		 ! (gC/m2) wood product C pool, 10-year lifespan
  real, pointer :: prod100c(:)  	 ! (gC/m2) wood product C pool, 100-year lifespan
  real, pointer :: seedn(:)		 ! (gN/m2) column-level pool for seeding new PFTs
  real, pointer :: prod10n(:)		 ! (gN/m2) wood product N pool, 10-year lifespan
  real, pointer :: prod100n(:)  	 ! (gN/m2) wood product N pool, 100-year lifespan
  real, pointer :: gpp(:)                ! (gC/m2/s) gross primary production 
  real, pointer :: npp(:)                ! (gC/m2/s) net primary production
  real, pointer :: sr(:)                 ! (gC/m2/s) total soil respiration (HR + root resp)
  real, pointer :: nee(:)                ! (gC/m2/s) net ecosystem exchange of carbon, includes fire, landuse, harvest, and hrv_xsmrpool flux, positive for source
  real, pointer :: farea_burned(:)       ! fractional area burned (proportion/s)
  real, pointer :: col_fire_closs(:)     ! (gC/m2/s) total column-level fire C loss
  real, pointer :: col_nfire(:)          ! (count/km2/s) column-level fire counts, valid only in Reg. C
  real, pointer :: col_somc_fire(:)      ! (gC/m2/s) column-level carbon emissions due to peat burning

! grid-level
  real, pointer :: latdeg(:)             ! latitude (degrees)
  real, pointer :: londeg(:)             ! longitude (degrees)
  real, pointer :: forc_hdm(:)           ! Human population density (individual/km2)
  real, pointer :: forc_lnfm(:)          ! Lightning frequency [Flashes/km^2/day]

! column level (need to organize the ones above here)
  real, pointer :: fsat(:)               ! fractional area with water table at surface  
  real, pointer :: tsoi17(:)             ! soil temperature in top 17cm of soil (Kelvin)
  real, pointer :: sucsat(:,:)           ! minimum soil suction (mm)
  real, pointer :: bd(:,:)               ! bulk density of dry soil material [kg/m^3]
  real, pointer :: watsat(:,:)           ! volumetric soil water at saturation (porosity) (nlevgrnd)
  real, pointer :: bsw(:,:)              ! Clapp and Hornberger "b" (nlevgrnd)
  integer, pointer :: abm_lf(:)          ! global peak month of crop fire emissions
  real, pointer :: peatf_lf(:)           ! global peatland fraction data (0-1) 
  real, pointer :: gdp_lf(:)             ! global real gdp data (k US$/capita)
  real, pointer :: watfc(:,:)            ! volumetric soil water at field capacity (nlevsoi)
  real, pointer :: snow_depth(:)         ! column averaged snow height (m). fzeng: use column averaged snow height instead of snow height of snow covered area as in the original CLM4.5 to avoid sharp drop of LAI due to light snow 
  
! pft level (need to organize the ones above here)
  real, pointer :: btran2(:)  
  real, pointer :: prec10(:)             ! 10-day running mean of tot. precipitation
  real, pointer :: prec60(:)             ! 60-day running mean of tot. precipitation
  real, pointer :: t10(:)                ! 10-day running mean of 2-m temperature (K)
  real, pointer :: lmrsun(:)             ! sunlit leaf maintenance respiration rate (umol CO2/m**2/s) 
  real, pointer :: lmrsha(:)             ! shaded leaf maintenance respiration rate (umol CO2/m**2/s) 
  logical, pointer :: pactive(:)         ! true=>do computations on this pft (see reweightMod for details), fzeng adopted a modified/simpler way to set pactive here

! PFT parameters
  real, pointer :: mxtmp(:)              ! Max Temperature, parameter used in accFlds (degree C) 
  real, pointer :: baset(:)              ! Base Temperature, parameter used in accFlds (degree C) 
  
! states and fluxes to understand nitrogen cycle, fzeng, 4 April 2019
  real, pointer :: ndeploy(:)        ! total N deployed to growth and storage (gN/m2/s)
  real, pointer :: denit(:)          ! total rate of denitrification (gN/m2/s)
  real, pointer :: sminn_leached(:)  ! soil mineral N pool loss to leaching (gN/m2/s)
  real, pointer :: sminn(:)          ! (gN/m2) soil mineral N
  real, pointer :: col_fire_nloss(:) ! (gN/m2/s) total column-level fire N loss
  real, pointer :: gross_nmin(:)     ! gross rate of N mineralization (gN/m2/s)
  real, pointer :: net_nmin(:)       ! vert-int (diagnostic) net rate of N mineralization (gN/m2/s)
  real, pointer :: nfix_to_sminn(:)  ! symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
  real, pointer :: actual_immob(:)   ! vert-int (diagnostic) actual N immobilization (gN/m2/s) 
  real, pointer :: fpg(:)            ! fraction of potential gpp (no units)
  real, pointer :: fpi(:)            ! fraction of potential immobilization (no units)
  real, pointer :: sminn_to_plant(:) ! vert-int (diagnostic) plant uptake of soil mineral N (gN/m2/s)
  real, pointer :: sminn_to_npool(:) ! deployment of soil mineral N uptake (gN/m2/s)
  real, pointer :: ndep_to_sminn(:)  ! atmospheric N deposition to soil mineral N (gN/m2/s)
  real, pointer :: totvegn(:)              ! (gN/m2) total vegetation nitrogen
  real, pointer :: totlitn(:)              ! (gN/m2) total litter nitrogen
  real, pointer :: totsomn(:)              ! (gN/m2) total soil organic matter nitrogen
  real, pointer :: retransn_to_npool(:)    ! deployment of retranslocated N (gN/m2/s)

! states and fluxes to understand the fire model, fzeng, 30 July 2019
  real, pointer :: fuelc(:)          ! fuel avalability for non-crop areas outside tropical closed broadleaf evergreen closed forests (gC/m2)
  real, pointer :: cwdc(:)           ! (gC/m2) coarse woody debris C

! define size of grid, landunit, column, and PFT
! ----------------------------------------------
  lbg = 1 ; ubg = nch
  lbl = 1 ; ubl = nch
  lbc = 1 ; ubc = nch*nzone
  lbp = 1 ; ubp = nch*nzone*npft ! potential PFT index (0-16); actual will be set in num_soilp filter
  
  ! fzeng, 29 Mar 2017, correct? should move to CN_init?
  numg = ubg
  numl = ubl
  numc = ubc
  nump = ubp

  if(first) then
    allocate (filter_soilc(ubc-lbc+1))
    allocate (filter_soilp(ubp-lbp+1))
    allocate (filter_pcropp(ubp-lbp+1))
    allocate (index_soilp(ubp-lbp+1))
    allocate (leafc_add(ubp-lbp+1))
    allocate (zone_soilp(ubp-lbp+1))
    allocate (zone_soilc(ubc-lbc+1)) 
    
    allocate (tileid_soilp(ubp-lbp+1))    ! fzeng added for debugging
       
  endif 

! assign local pointers to derived type members 

! grid level
  latdeg        		 => grc%latdeg
  londeg        		 => grc%londeg
  forc_ndep     		 => grc%forc_ndep
  forc_rh       		 => grc%forc_rh
  forc_wind     		 => grc%forc_wind
  forc_t        		 => grc%forc_t
  forc_rain     		 => grc%forc_rain
  forc_snow     		 => grc%forc_snow
  forc_hdm      		 => grc%forc_hdm
  forc_lnfm     		 => grc%forc_lnfm
  
! land level
  ifspecial     		 => lun%ifspecial  

! column level 
  clandunit     		 => col%landunit
  npfts         		 => col%npfts
  cgridcell     		 => col%gridcell
  cwtgcell                       => col%wtgcell
  pfti          		 => col%pfti
  pftf          		 => col%pftf
  t_soisno      		 => ces%t_soisno
  dz            		 => cps%dz
  psisat        		 => cps%psisat
  psiwilt       		 => cps%psiwilt
  soilpsi       		 => cps%soilpsi
  h2osoi_liq    		 => cws%h2osoi_liq
  qflx_drain    		 => cwf%qflx_drain
  qflx_surf     		 => cwf%qflx_surf
  snow_depth    		 => cps%snow_depth
  t_grnd        		 => ces%t_grnd 
  wf            		 => cps%wf
  wf2           		 => cps%wf2
  fsat          		 => cws%fsat
  tsoi17        		 => ces%tsoi17
  sucsat        		 => cps%sucsat
  bd            		 => cps%bd
  watsat        		 => cps%watsat
  bsw           		 => cps%bsw
  abm_lf        		 => cps%abm_lf
  peatf_lf      		 => cps%peatf_lf
  gdp_lf        		 => cps%gdp_lf
  watfc         		 => cps%watfc
  totcolc       		 => ccs%totcolc
  col_ctrunc			 => ccs%col_ctrunc
  totlitc			 => ccs%totlitc
  totsomc			 => ccs%totsomc  
  sminn_vr 			 => cns%sminn_vr
  col_ntrunc			 => cns%col_ntrunc
  totcoln			 => cns%totcoln
  seedc 			 => ccs%seedc
  prod10c			 => ccs%prod10c
  prod100c			 => ccs%prod100c
  seedn 			 => cns%seedn
  prod10n			 => cns%prod10n
  prod100n			 => cns%prod100n

! pft level  
  btran2        		 => pps%btran2
  prec10        		 => pps%prec10
  prec60        		 => pps%prec60
  t10           		 => pes%t10
  pactive       		 => pft%active
  pitype        		 => pft%itype
  pcolumn       		 => pft%column
  pwtcol        		 => pft%wtcol
  t_ref2m       		 => pes%t_ref2m
  psnsha        		 => pcf%psnsha
  psnsun        		 => pcf%psnsun
  lmrsha                         => pcf%lmrsha
  lmrsun                         => pcf%lmrsun
  laisha        		 => pps%laisha
  laisun        		 => pps%laisun
  plandunit     		 => pft%landunit
  pgridcell     		 => pft%gridcell
  pwtgcell      		 => pft%wtgcell
  rootfr        		 => pps%rootfr  
  forc_hgt_u_pft		 => pps%forc_hgt_u_pft
  dayl          		 => pepv%dayl
  prev_dayl     		 => pepv%prev_dayl
  elai          		 => pps%elai
  esai          		 => pps%esai
  tlai          		 => pps%tlai
  annavg_t2m    		 => pepv%annavg_t2m  
  leafc 			 => pcs%leafc
  leafc_storage 		 => pcs%leafc_storage
  leafc_xfer			 => pcs%leafc_xfer
  frootc			 => pcs%frootc
  frootc_storage		 => pcs%frootc_storage
  frootc_xfer			 => pcs%frootc_xfer
  livestemc			 => pcs%livestemc
  livestemc_storage		 => pcs%livestemc_storage
  livestemc_xfer		 => pcs%livestemc_xfer
  deadstemc			 => pcs%deadstemc
  deadstemc_storage		 => pcs%deadstemc_storage
  deadstemc_xfer		 => pcs%deadstemc_xfer
  livecrootc			 => pcs%livecrootc
  livecrootc_storage		 => pcs%livecrootc_storage
  livecrootc_xfer		 => pcs%livecrootc_xfer
  deadcrootc			 => pcs%deadcrootc
  deadcrootc_storage		 => pcs%deadcrootc_storage
  deadcrootc_xfer		 => pcs%deadcrootc_xfer
  gresp_storage 		 => pcs%gresp_storage
  gresp_xfer			 => pcs%gresp_xfer
  cpool 			 => pcs%cpool
  xsmrpool			 => pcs%xsmrpool
  pft_ctrunc                     => pcs%pft_ctrunc
  totvegc                        => pcs%totvegc
  frac_veg_nosno                 => pps%frac_veg_nosno

  leafn 			 => pns%leafn
  leafn_storage 		 => pns%leafn_storage
  leafn_xfer			 => pns%leafn_xfer
  frootn			 => pns%frootn
  frootn_storage		 => pns%frootn_storage
  frootn_xfer			 => pns%frootn_xfer
  livestemn			 => pns%livestemn
  livestemn_storage		 => pns%livestemn_storage
  livestemn_xfer		 => pns%livestemn_xfer
  deadstemn			 => pns%deadstemn
  deadstemn_storage		 => pns%deadstemn_storage
  deadstemn_xfer		 => pns%deadstemn_xfer
  livecrootn			 => pns%livecrootn
  livecrootn_storage		 => pns%livecrootn_storage
  livecrootn_xfer		 => pns%livecrootn_xfer
  deadcrootn			 => pns%deadcrootn
  deadcrootn_storage		 => pns%deadcrootn_storage
  deadcrootn_xfer		 => pns%deadcrootn_xfer
  retransn			 => pns%retransn
  npool 			 => pns%npool
  pft_ntrunc                     => pns%pft_ntrunc

! PFT parameters
  mxtmp         		 => pftcon%mxtmp
  baset         		 => pftcon%baset

! states and fluxes to understand nitrogen cycle, fzeng, 4 April 2019
  ndeploy                        => pnf%ndeploy         
  denit                          => cnf%denit           
  sminn_leached                  => cnf%sminn_leached
  sminn                          => cns%sminn           
  col_fire_nloss                 => cnf%col_fire_nloss  
  leafn                          => pns%leafn           
  leafc                          => pcs%leafc
  gross_nmin                     => cnf%gross_nmin
  net_nmin                       => cnf%net_nmin  
  nfix_to_sminn                  => cnf%nfix_to_sminn
  actual_immob                   => cnf%actual_immob
  fpg                            => cps%fpg           
  fpi                            => cps%fpi           
  sminn_to_plant                 => cnf%sminn_to_plant
  sminn_to_npool                 => pnf%sminn_to_npool
  ndep_to_sminn                  => cnf%ndep_to_sminn 
  totvegn                        => pns%totvegn          
  totlitn                        => cns%totlitn          
  totsomn                        => cns%totsomn                   
  retransn_to_npool              => pnf%retransn_to_npool

! states and fluxes to understand the fire model, fzeng, 30 July 2019
  fuelc                          => ccs%fuelc
  cwdc                           => ccs%cwdc

! define landunit & column settings
! ---------------------------------
! litype(:)       = 1       ! all land-units are soil
  ifspecial(:)    = .false. ! no special land-units; all are soil

! citype(:)       = 1       ! all columns are soil (vegetated or bare)
  clandunit(:)    = 1       ! all landunits are soil
  plandunit(:)    = 1       ! all landunits are soil
  npfts(:)        = npft    ! max number of PFTs per column

  num_soilc = nch*nzone     ! number of columns = number of catchments*zones

! map vegetation types into closest PFT & map PFT into column; assign PFT weight & filter
! ---------------------------------------------------------------------------------------
  num_soilp = 0             ! initialize PFT filter
  num_pcropp = 0            ! initialize crop PFT filter
  filter_soilp(:) = 0       ! set PFT index to invalid number
  filter_pcropp(:) = 0      ! set prognostic crop PFT index to invalid number
  index_soilp(:) = 0        ! set veg index to invalid number
  zone_soilp(:) = 0         ! set zone index to invalid number

  latdeg = latitude
  londeg = longitude

  dt = real( get_step_size() )
  
  ! Following CLM4.5's CNiniSpecial.F90
  frac_veg_nosno = 0

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
    cwtgcell(n) = wtzone(nc,nz)   ! weight of columns relative to gridcells

! set column soil properties
! --------------------------
    dz(n,1) = 1.00                                                ! hydrologically active soil layer thickness (m)
    psisat(n,1) = 1.e-6*psis(nc)*grav*denh2o                      ! saturated soil water potential m -> Mpa 
    sucsat(n,1) = psis(nc) * 1e3 * (-1)                           ! minimum soil suction (mm), psis is in (m)
    soilpsi(n,1) = 1.e-6*psis(nc)*grav*denh2o*rzm(nc,nz)**(-bee(nc)) 
    psiwilt(n) = 1.e-6*psis(nc)*grav*denh2o*wpwet(nc)**(-bee(nc)) ! root-zone wilting soil water potential (Mpa)
    bd(n,1) = (1.-poros(nc))*2.7e3                                ! see iniTimeConst.F90 in CLM4.5
    watsat(n,1) = poros(nc)
    bsw(n,1) = bee(nc)
    watfc(n,1) = field_cap(nc)

! soil liquid water
! ------------------
    h2osoi_liq(n,1) = totwat(nc)  ! soil liquid water, kg/m2
    qflx_drain(n) = bflow(nc)     ! sub-surface runoff (mm H2O /s); used for soil nitrogen leaching
    qflx_surf(n)  = runsrf(nc)    ! surface runoff (mm H2O /s); used for soil nitrogen leaching
    wf(n)  = sfm(nc,nz)           ! soil water as frac of WHC for top 5cm; used in CNFireMod
    wf2(n) = rzm(nc,nz)           ! soil water as frac of WHC for top 17cm; used in CNFireMod    

! soil temperature
! ----------------
    t_soisno(n,1) = tp1(nc)       ! soil layer temperature (K)
    t_grnd(n) = tgw(nc,nz)        ! ground surface temperature (K)
    tsoi17(n) = (DZTSURF*tgw(nc,nz)+(DZGT(1)-DZTSURF)*tp1(nc)+(0.17-DZGT(1))*tp2(nc))/0.17        ! soil temperature in top 17cm of soil (Kelvin)
                                                                ! fzeng: tgw is for the top 5cm; tp1 is for the 2nd 5cm; tp2 is for the next 10cm
                                                                ! see Koster et al., 2000, JGR
                                                                ! The depths are hard coded here. Improve this?  
                                  

! new parameters for the fire module in CLM4.5
! --------------------------------------------
    abm_lf(n)   = abm(nc)
    peatf_lf(n) = peatf(nc)
    gdp_lf(n)   = gdp(nc)   

! snow depth 
! ----------
    snow_depth(n) = sndzn(nc)*asnow(nc)      ! column averaged snow height (m)
    
    ! assuming nzone = 3, fzeng, 6 Mar 2017
    if(nz==1) then
      fsat(n) = min(max(0.,car1m(nc)/wtzone(nc,nz)),1.)
     elseif(nz==2) then
      fsat(n) = min(max(0.,(car1m(nc)-wtzone(nc,1))/wtzone(nc,nz)),1.)
     else
      fsat(n) = min(max(0.,(car1m(nc)-wtzone(nc,1)-wtzone(nc,2))/wtzone(nc,nz)),1.)
    endif          

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
     prec10(p)  = prec10d(nc)    
     prec60(p)  = prec60d(nc)    
     t10(p)     = t10d(nc) 
     pgridcell(p) = nc            ! PFT map into catchment tile
     pwtgcell(p) = 0.             ! PFT weight in catchment tile
     pactive(p) = .false.         ! pactive will be .false. unless otherwise set, fzeng

! map bare soil if present
! ------------------------
     if(bare.gt.0. .and. pft_num.eq.0) then
       num_soilp = num_soilp + 1
       filter_soilp(num_soilp) = p
       pwtcol(p) = bare
       psnsha(p) = 0.
       psnsun(p) = 0.
       lmrsha(p) = 0.
       lmrsun(p) = 0.
       laisha(p) = 0.
       laisun(p) = 0.
       forc_hgt_u_pft(p) = 30. ! gkw: may need from land model; use this for now
       rootfr(p,1) = 0.0
       pwtgcell(p) = bare*wtzone(nc,nz) ! PFT weight in catchment tile
       btran2(p) = 0.
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
         lmrsha(p) = lmrsham(nc,nv,nz)
         lmrsun(p) = lmrsunm(nc,nv,nz)
         laisha(p) = lsha(nc,nv,nz)
         laisun(p) = lsun(nc,nv,nz)
         forc_hgt_u_pft(p) = 30. ! gkw: may need from land model; use this for now
         rootfr(p,1) = 1.0 ! gkw: affects maint resp. test sensitivity
         pwtgcell(p) = fveg(nc,nv,nz)*wtzone(nc,nz) ! PFT weight in catchment tile
         
         ! calculate btran2 for CNFireMod
         ! use btran derived from rzmc for trees except broadleaf deciduous tropical tree
         ! use btran derived from sfmc for the other vegetation types
         if(pft_num<=ntree .and. pft_num/=nbrdlf_dcd_trp_tree) then  
            btran2(p) = btran2x_rz(nc,nz)
          else
            btran2(p) = btran2x_sf(nc,nz)
         endif          
         
         pactive(p) = .true.     ! set pactive to .true. if there is vegetation, fzeng
         
         tileid_soilp(p) = tile_id(nc)  ! fzeng added for debugging
         
         ! for prog. crop. correct? fzeng, 20 Mar 2017
         if(ityp(nc,nv,nz) >= npcropmin) then
           num_pcropp = num_pcropp + 1
           filter_pcropp(num_pcropp) = p
         endif

         ! set daylength here (moved from CNPhenology)
         if(first .and. istep==0) dayl(p) = daylength(nc) ! working fix for cold carbon start gkw: 2015-08-01
         prev_dayl(p) = dayl(p)
         dayl(p) = daylength(nc)
       
         ! Fraction of vegetation free of snow, following CNVegStructUpdate and CLM4.5's initSurfAlbMod.F90
         if ((elai(p) + esai(p)) > 0._r8) then
            frac_veg_nosno(p) = 1
         else
            frac_veg_nosno(p) = 0
         end if                
         
       endif       
       
     end do

     pft_num = pft_num + 1
    end do ! end PFT loop

  end do ! end zone loop

  end do ! end catchment loop

! nitrogen deposition
! -------------------
  forc_ndep(:) = ndep(:)

! forcing for CNFireMod
! ---------------------      
  forc_rh(:)   = rhm(:)    
  forc_wind(:) = windm(:)
  forc_t(:)    = tm(:)   
  forc_rain(:) = rainfm(:) 
  forc_snow(:) = snowfm(:)
  forc_hdm(:)  = hdm(:) 
  forc_lnfm(:) = lnfm(:) 

! Initialize CN Ecosystem Dynamics (must be after time-manager initialization)
! ----------------------------------------------------------------------------
  if(first) call CNEcosystemDynInit(lbg,ubg,lbc,ubc,lbp,ubp)

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

! update CN time step number, fzeng, 20 Mar 2017
! ----------------------------------------------
  nstep_cn = get_nstep(istep)          ! nstep_cn is actually not being used here. By doing so, istep is passed to the CN routines whenever "nstep = get_nstep()" is called.

! initialize balance checks
!  ------------------------
  call BeginCBalance(lbc, ubc, num_soilc, filter_soilc)
  call BeginNBalance(lbc, ubc, num_soilc, filter_soilc)

  call CNZeroFluxes_dwt(lbc, ubc, lbp, ubp)

  call CNEcosystemDyn(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
                      num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, tileid_soilp)                 

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
  ndeployg = 0.
  leafng   = 0.
  leafcg   = 0.
  sminn_to_npoolg = 0.
  totvegng = 0.
  retransng = 0.
  retransn_to_npoolg = 0.
  rootcg   = 0.
         
  npp => pcf%npp
  gpp => pcf%gpp

  do n = 1,num_soilp
    p = filter_soilp(n)
    i = index_soilp(n)   ! veg index, 1 to 4
    z = zone_soilp(n)    ! zone index, 1 to 3
    c = pgridcell(p)     ! tile index, 1 to nch
    if(i .gt. 0) then    ! veg exists
      zlai(c,i,z) = elai(p)
      zsai(c,i,z) = esai(p)
      ztai(c,i,z) = tlai(p)
      nppg(c) = nppg(c) + npp(p)*pwtgcell(p)
      gppg(c) = gppg(c) + gpp(p)*pwtgcell(p)
      padd(c) = padd(c) + leafc_add(p)*pwtgcell(p)
      root(c) = root(c) + (frootc(p)+frootc_storage(p)+frootc_xfer(p))*pwtgcell(p)
      vegc(c) = vegc(c) + totvegc(p)*pwtgcell(p)
      xsmr(c) = xsmr(c) + xsmrpool(p)*pwtgcell(p)
      ndeployg(c) = ndeployg(c) + ndeploy(p)*pwtgcell(p)
      leafng(c)   = leafng(c)   + leafn(p)*pwtgcell(p)
      leafcg(c)   = leafcg(c)   + leafc(p)*pwtgcell(p)
      sminn_to_npoolg(c)   = sminn_to_npoolg(c)   + sminn_to_npool(p)*pwtgcell(p)
      totvegng(c)   = totvegng(c)   + totvegn(p)*pwtgcell(p)
      retransng(c)  = retransng(c)  + retransn(p)*pwtgcell(p)
      retransn_to_npoolg(c)   = retransn_to_npoolg(c)   + retransn_to_npool(p)*pwtgcell(p)
      rootcg(c) = rootcg(c) + (frootc(p)+frootc_storage(p)+frootc_xfer(p)+ &
                               livecrootc(p)+livecrootc_storage(p)+livecrootc_xfer(p)+ &
                               deadcrootc(p)+deadcrootc_storage(p)+deadcrootc_xfer(p))*pwtgcell(p)
    endif
  end do

! column level diags
! ------------------
  srg = 0.
  neeg = 0.
  burn = 0.
  closs = 0.
  nfire = 0.
  som_closs = 0.  
  denitg          = 0.                
  sminn_leachedg  = 0.                
  sminng          = 0.                
  col_fire_nlossg = 0.                
  gross_nming     = 0.                
  net_nming       = 0.                
  nfix_to_sminng  = 0.                
  actual_immobg   = 0. 
  fpgg            = 0.                
  fpig            = 0.  
  sminn_to_plantg = 0.
  ndep_to_sminng  = 0.
  totlitng        = 0.
  totsomng        = 0.
  fuelcg          = 0.
  totlitcg        = 0.
  cwdcg           = 0.

  sr  => ccf%sr
  nee => ccf%nee
  farea_burned   => cps%farea_burned
  col_fire_closs => ccf%col_fire_closs
  col_nfire      => cps%nfire
  col_somc_fire  => ccf%somc_fire

  do n = 1,num_soilc
    i = filter_soilc(n)
    z = zone_soilc(n)
    c = cgridcell(i)
    colc(c,z) = totcolc(i)
    srg(c)  = srg(c)  + sr(i) *wtzone(c,z)
    neeg(c) = neeg(c) + nee(i)*wtzone(c,z)
    burn(c) = burn(c) + farea_burned(i)*wtzone(c,z)   ! burn rate (fraction per second)
    closs(c) = closs(c) + col_fire_closs(i)*wtzone(c,z)
    nfire(c) = nfire(c) + col_nfire(i)*wtzone(c,z)    ! fire counts (count/km2/s)
    som_closs(c) = som_closs(c) + col_somc_fire(i)*wtzone(c,z)
    denitg(c)          = denitg(c)          + denit(i) * wtzone(c,z)         
    sminn_leachedg(c)  = sminn_leachedg(c)  + sminn_leached(i) * wtzone(c,z)
    sminng(c)          = sminng(c)          + sminn(i) * wtzone(c,z)         
    col_fire_nlossg(c) = col_fire_nlossg(c) + col_fire_nloss(i) * wtzone(c,z)
    gross_nming(c)     = gross_nming(c)     + gross_nmin(i) * wtzone(c,z)    
    net_nming(c)       = net_nming(c)       + net_nmin(i) * wtzone(c,z)      
    nfix_to_sminng(c)  = nfix_to_sminng(c)  + nfix_to_sminn(i) * wtzone(c,z) 
    actual_immobg(c)   = actual_immobg(c)   + actual_immob(i) * wtzone(c,z)
    fpgg(c)            = fpgg(c)            + fpg(i) * wtzone(c,z)
    fpig(c)            = fpig(c)            + fpi(i) * wtzone(c,z)
    sminn_to_plantg(c) = sminn_to_plantg(c) + sminn_to_plant(i) * wtzone(c,z)
    ndep_to_sminng(c)  = ndep_to_sminng(c)  + ndep_to_sminn(i) * wtzone(c,z)
    totlitng(c)        = totlitng(c)        + totlitn(i) * wtzone(c,z)
    totsomng(c)        = totsomng(c)        + totsomn(i) * wtzone(c,z) 
    fuelcg(c)          = fuelcg(c)          + fuelc(i) * wtzone(c,z)
    totlitcg(c)        = totlitcg(c)        + totlitc(i) * wtzone(c,z)
    cwdcg(c)           = cwdcg(c)           + cwdc(i) * wtzone(c,z)  
  end do
  
  if ( .not. LAND_FIX ) then ! jkolassa Oct 2020: the if-wrapper here is to toggle between the LDASsa version used by Fanwei Zeng and Eunjee Lee and current GEOSldas Catchment-CN; there is likely a better way to control this
     where (zlai > 20.) zlai = 20.
     where (zsai > 20.) zsai = 20.
  end if

  end subroutine CN_Driver

  subroutine CN_init(istep,nch,nveg,nzone,ityp,fveg,var_col,var_pft,cncol,cnpft,skip_initCN)

  use clm_varpar,        only: nlevdecomp
  use clm_varcon,        only: dzsoi_decomp

  integer*8, intent(in) :: istep
  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg    ! PFT fraction

  integer, intent(in) :: var_col ! number of CN column restart variables
  real*4, dimension(nch,nzone,var_col), optional, intent(in) :: cncol ! gkw: column CN restart 

  integer, intent(in) :: var_pft ! number of CN PFT restart variables
  real*4, dimension(nch,nzone,nveg,var_pft), optional, intent(in) :: cnpft ! gkw: PFT CN restart 
  logical,optional, intent(in) :: skip_initCN

  integer :: n, p, nv, nc, nz, np, j

! PFT parameters note: index 0 is "noveg"

  real, save, dimension(0:numpft) :: qe2x,z0mx,displax,c3psx,mx,slatox,dsladlax,leafcx,flnx,fnitx,woodx,lflitcx
  real, save, dimension(0:numpft) :: frootcx,livewdcx,deadwdcx,dwoox,froot_leax,stem_leax,croot_stex,flivewx,fcux
  real, save, dimension(0:numpft) :: lf_flax,lf_fcex,lf_flix,fr_flax,fr_fcex,fr_flix,leaf_lonx,evergreex,resisx
  real, save, dimension(0:numpft) :: stress_decix,season_decix,xlx,rholx,rhosx,taulx,tausx
  real, save, dimension(0:numpft) :: cc_dstex,cc_leax,cc_lstex,cc_othex,fm_leax,fm_lstex,fm_othex,fm_roox,fm_lroox,fm_droox      ! new in CLM4.5
  real, save, dimension(0:numpft) :: fd_pfx,fsr_pfx,grperx,grpnox                                                                ! new in CLM4.5 
  
! real, save, dimension(0:numpft) :: fertnitrx,lfemerx,grnfilx,mxmax,hybgdx                                                      ! new in CLM4.5 
! real, save, dimension(0:numpft) :: laimxx,ztopmxx,gddmix,graincx,fleafcx,ffrootcx,fstemcx                                      ! new in CLM4.5 
! real, save, dimension(0:numpft) :: rootprof_betx,aleafx,allconslx,allconssx,arootfx,arootix,astemx,bfacx,declfacx,fleafx       ! new in CLM4.5 
! real, save, dimension(0:numpft) :: basex,mxtmx,planttemx,minplanttemx                                                          ! new in CLM4.5 
! integer, save, dimension(0:numpft) :: mnNHplantdatx,mxNHplantdatx,mnSHplantdatx,mxSHplantdatx                                  ! new in CLM4.5 
    

! The data below is from 
! https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/pftdata/pft-physiology.c130503.nc
! A copy of this file is /discover/nobackup/fzeng/clm4-to-clm4.5/data/pftdata4.5/pft-physiology.c130503.nc
! In CLM4.5, data in this file is read in pftvarcon.F90.

! From /gpfsm/dnb02/smahanam/bcs/Heracles-4_3/Heracles-4_3_MERRA-3/DE_00720x00360_PE_0720x0360/clsm/README
! It matches that in pft-physiology.c130503.nc with 3 additional split types (i.e. 11, 15 and 17).                                            

! NEW  pftname 
!  0   Bare                                              
!  1   Needleleaf evergreen temperate tree                        
!  2   Needleleaf evergreen boreal tree                          
!  3   Needleleaf deciduous boreal tree                           
!  4   Broadleaf evergreen tropical tree                          
!  5   Broadleaf evergreen temperate tree                         
!  6   Broadleaf deciduous tropical tree                         
!  7   Broadleaf deciduous temperate tree                         
!  8   Broadleaf deciduous boreal tree                            
!  9   Broadleaf evergreen temperate shrub                        
! 10   Broadleaf deciduous temperate shrub [moisture + deciduous]                       
! 11   Broadleaf deciduous temperate shrub [moisture stress only] 
! 12   Broadleaf deciduous boreal shrub                           
! 13   Arctic c3 grass                                            
! 14   Cool c3 grass [moisture + deciduous]                                             
! 15   Cool c3 grass [moisture stress only]                       
! 16   Warm c4 grass [moisture + deciduous]                                             
! 17   Warm c4 grass [moisture stress only]                       
! 18   C3 crop [moisture + deciduous]                                                   
! 19   C3 crop [moisture stress only]                                                                                         

! fzeng: qe25 does not exist in pft-physiology.c130503.nc. Give it the same numbers (0.06 for c3 plants and 0.04 for c4 plants) as the CLM4 version here just to get the code compiling. 
! Need to ask Jung-Eun Lee for the CLM4.5 version of fluorescence calculation!!
! When remove qe25 later from here and compute_rc, also need to remove it from clmtype and clmtypeInit!

! Some of these parameters are only used for prognostic crops. Although we removed prognostic crops, keep them here for now to avoid having to check one by one which is not used.
! fzeng, 10 May 2018

! fzeng, 8 May 2019: updated fsr_pfx following /discover/nobackup/fzeng/clm4-to-clm4.5/data/paramdata4.5/clm_params.c140423.nc (Li et al., BG 2014)
! fzeng, 12 July 2019: modified fsr_pfx for type 6 because this type (i.e. broadleaf deciduous tropical trees in ESA) in Africa is classified as woody savanna and savanna in MODIS land cover which is more consistent with CLM4.5CN tree and grass fractions in Africa
 
! pftname            /   0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17     18     19 /      
  data     qe2x      /   0., 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.060, 0.040, 0.040, 0.060, 0.060/ 
  data     z0mx      /   0., 0.055, 0.055, 0.055, 0.075, 0.075, 0.055, 0.055, 0.055, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120/
  data  displax      /   0., 0.670, 0.670, 0.670, 0.670, 0.670, 0.670, 0.670, 0.670, 0.680, 0.680, 0.680, 0.680, 0.680, 0.680, 0.680, 0.680, 0.680, 0.680, 0.680/
  data    c3psx      /   1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    0.,    0.,    1.,    1./
  data   slatox      /   0., 0.010, 0.008, 0.024, 0.012, 0.012, 0.030, 0.030, 0.030, 0.012, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030, 0.030/
  data dsladlax      /   0.,0.00125,0.001, 0.003,0.0015,0.0015, 0.004, 0.004, 0.004,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
  data   leafcx      /   1.,   35.,   40.,   25.,   30.,   30.,   25.,   25.,   25.,   30.,   25.,   25.,   25.,   25.,   25.,   25.,   25.,   25.,   25.,   25./
  data     flnx      /   0.,0.0509,0.0466,0.0546,0.0461,0.0515,0.0716,0.1007,0.1007,0.0517,0.0943,0.0943,0.0943,0.1365,0.1365,0.1365,0.0900,0.0900,0.1758,0.1758/
  data    fnitx      /   0.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1./ ! largely changed from CLM4!
  data    woodx      /   0.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
  data  lflitcx      /   1.,   70.,   80.,   50.,   60.,   60.,   50.,   50.,   50.,   60.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50./
  data  frootcx      /   1.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42.,   42./
  data livewdcx      /   1.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,   50.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
  data deadwdcx      /   1.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,  500.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
  data froot_leax    /   0.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    2.,    2.,    2.,    2.,    2.,    2./
  data stem_leax     /   0.,   -1.,   -1.,   -1.,   -1.,   -1.,   -1.,   -1.,   -1.,   0.2,   0.2,   0.2,   0.2,    0.,    0.,    0.,    0.,    0.,    0.,    0./
  data croot_stex    /   0.,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,    0.,    0.,    0.,    0.,    0.,    0.,    0./
  data flivewx       /   0.,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.5,   0.5,   0.5,   0.1,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data fcux          /   0.,    1.,    1.,    0.,    1.,    1.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./ ! original CLM4.5 values
  data fcux          /   0.,   1.0,   1.0,   0.5,   1.0,   1.0,   0.5,   0.5,   0.5,   1.0,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/ ! our fix, same as CLM4 Catchment-CN, fzeng, 10 May 2018
  data lf_flax       /   0.,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25/
  data lf_fcex       /   0.,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/
  data lf_flix       /   0.,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25/
  data fr_flax       /   0.,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25/
  data fr_fcex       /   0.,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/
  data fr_flix       /   0.,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25/
  data leaf_lonx     /   0.,    3.,    6.,    1.,   1.5,   1.5,    1.,    1.,    1.,   1.5,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1./
  data evergreex     /   0.,    1.,    1.,    0.,    1.,    1.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./ 
  data stress_decix  /   0.,    0.,    0.,    0.,    0.,    0.,    1.,    0.,    0.,    0.,    1.,    1.,    0.,    0.,    1.,    1.,    1.,    1.,    1.,    1./
  data season_decix  /   0.,    0.,    0.,    1.,    0.,    0.,    0.,    1.,    1.,    0.,    0.,    0.,    1.,    1.,    0.,    0.,    0.,    0.,    0.,    0./    
  data xlx           /   0.,  0.01,  0.01,  0.01,   0.1,   0.1,  0.01,  0.25,  0.25,  0.01,  0.25,  0.25,  0.25,  -0.3,  -0.3,  -0.3,  -0.3,  -0.3,  -0.3,  -0.3/
  data rholx         /   0.,  0.07,  0.07,  0.07,   0.1,   0.1,   0.1,   0.1,   0.1,  0.07,   0.1,   0.1,   0.1,  0.11,  0.11,  0.11,  0.11,  0.11,  0.11,  0.11/
  data rhosx         /   0.,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.16,  0.31,  0.31,  0.31,  0.31,  0.31,  0.31,  0.31/
  data taulx         /   0.,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05/
  data tausx         /   0., 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12/

! data cc_dstex      /   0.,  0.22,  0.25,  0.25,  0.22,  0.22,  0.22,  0.22,  0.22,   0.3,   0.3,   0.3,   0.3,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values  
! data cc_leax       /   0.,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values
! data cc_lstex      /   0.,  0.22,  0.25,  0.25,  0.22,  0.22,  0.22,  0.22,  0.22,   0.3,   0.3,   0.3,   0.3,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values
! data cc_othex      /   0.,  0.45,   0.5,   0.5,  0.45,  0.45,  0.45,  0.45,  0.45,  0.55,  0.55,  0.55,  0.55,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values
! data fm_leax       /   0.,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values
! data fm_lstex      /   0.,  0.45,   0.5,   0.5,  0.45,  0.45,  0.35,  0.35,  0.45,  0.55,  0.55,  0.55,  0.55,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values
! data fm_othex      /   0.,  0.45,   0.5,   0.5,  0.45,  0.45,  0.35,  0.35,  0.45,  0.55,  0.55,  0.55,  0.55,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! original CLM4.5 values
! data fm_roox       /   0.,  0.13,  0.15,  0.15,  0.13,  0.13,   0.1,   0.1,  0.13,  0.17,  0.17,  0.17,  0.17,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2/ ! original CLM4.5 values
! data fm_lroox      /   0.,  0.13,  0.15,  0.15,  0.13,  0.13,   0.1,   0.1,  0.13,  0.17,  0.17,  0.17,  0.17,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2/ ! original CLM4.5 values
! data fm_droox      /   0.,  0.13,  0.15,  0.15,  0.13,  0.13,   0.1,   0.1,  0.13,  0.17,  0.17,  0.17,  0.17,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2/ ! original CLM4.5 values

  data cc_dstex      /   0.,  0.22,  0.25,  0.25,  0.22,  0.22,  0.60,  0.22,  0.22,   0.3,   0.3,   0.3,   0.3,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019   
  data cc_leax       /   0.,   0.8,   0.8,   0.8,   0.8,   0.8,  0.80,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data cc_lstex      /   0.,  0.22,  0.25,  0.25,  0.22,  0.22,  0.60,  0.22,  0.22,   0.3,   0.3,   0.3,   0.3,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data cc_othex      /   0.,  0.45,   0.5,   0.5,  0.45,  0.45,  0.68,  0.45,  0.45,  0.55,  0.55,  0.55,  0.55,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data fm_leax       /   0.,   0.8,   0.8,   0.8,   0.8,   0.8,  0.80,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data fm_lstex      /   0.,  0.45,   0.5,   0.5,  0.45,  0.45,  0.64,  0.35,  0.45,  0.55,  0.55,  0.55,  0.55,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data fm_othex      /   0.,  0.45,   0.5,   0.5,  0.45,  0.45,  0.64,  0.35,  0.45,  0.55,  0.55,  0.55,  0.55,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8,   0.8/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data fm_roox       /   0.,  0.13,  0.15,  0.15,  0.13,  0.13,  0.17,   0.1,  0.13,  0.17,  0.17,  0.17,  0.17,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data fm_lroox      /   0.,  0.13,  0.15,  0.15,  0.13,  0.13,  0.17,   0.1,  0.13,  0.17,  0.17,  0.17,  0.17,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data fm_droox      /   0.,  0.13,  0.15,  0.15,  0.13,  0.13,  0.17,   0.1,  0.13,  0.17,  0.17,  0.17,  0.17,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2,   0.2/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019

  data fd_pfx        /   0.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24.,   24./
! data fsr_pfx       /   0.,   0.4,  0.43,  0.43,   0.4,   0.4,   0.4,   0.4,   0.4,  0.46,  0.46,  0.46,  0.46,  0.55,  0.55,  0.55,  0.55,  0.55,  0.55,  0.55/ ! original CLM4.5 values
! data fsr_pfx       /   0.,   0.4,  0.43,  0.43,   0.4,   0.4,   0.4,   0.4,   0.4,  0.46,  0.46,  0.46,  0.46,  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  0.60/ ! original CLM4.5 values after bug fix in Li et al., 2014
  data fsr_pfx       /   0.,   0.4,  0.43,  0.43,   0.4,   0.4,  0.53,   0.4,   0.4,  0.46,  0.46,  0.46,  0.46,  0.60,  0.60,  0.60,  0.60,  0.60,  0.60,  0.60/ ! assume type 6 is 35% tree and 65% grass, fzeng, 2019
  data grperx        /  0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3/
  data grpnox        /   1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1.,    1./
     
! data fertnitrx     /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data lfemerx       /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data grnfilx       /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data mxmax         /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data hybgdx        /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./  
! data laimxx        /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data ztopmxx       /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data gddmix        /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data graincx       /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./  
! data fleafcx       / 999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999./
! data ffrootcx      / 999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999./
! data fstemcx       / 999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999.,  999./  
! data rootprof_betx /   0., 0.976, 0.943, 0.943, 0.962, 0.966, 0.961, 0.966, 0.943, 0.964, 0.964, 0.964, 0.914, 0.914, 0.943, 0.943, 0.943, 0.943, 0.961, 0.961/
! data aleafx        /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./ 
! data allconslx     /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data allconssx     /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data arootfx       /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data arootix       /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data astemx        /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data bfacx         /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data declfacx      /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data fleafx        /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data basex         /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data mxtmx         /   0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0./
! data planttemx     /1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000./
! data minplanttemx  /1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000./
! data mnNHplantdatx /    0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
! data mxNHplantdatx /    0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
! data mnSHplantdatx /    0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0/
! data mxSHplantdatx /    0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0/

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
  call clm_varcon_init()
  call initClmtype(lbg,ubg,lbl,ubl,lbc,ubc,lbp,ubp) ! allocation & initialization
  
! Initialize time-constant arrays of decomposition constants
! ----------------------------------------------------------
  call init_decompcascade(lbc, ubc)

! initialize PFT parameters
! -------------------------
  pftcon%qe25          = qe2x	       ! quantum efficiency at 25C (umol CO2 / umol photon)
  pftcon%z0mr          = z0mx	       ! ratio of momentum roughness length to canopy top height (-)
  pftcon%displar       = displax       ! ratio of displacement height to canopy top height (-)
  pftcon%c3psn         = c3psx         ! photosynthetic pathway: 0. = c4, 1. = c3
  pftcon%slatop        = slatox	       ! specific leaf area at top of canopy, projected area basis [m^2/gC]
  pftcon%dsladlai      = dsladlax      ! dSLA/dLAI, projected area basis [m^2/gC]
  pftcon%leafcn        = leafcx	       ! leaf C:N (gC/gN)
  pftcon%flnr          = flnx	       ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
  pftcon%fnitr         = fnitx	       ! foliage nitrogen limitation factor (-)
  pftcon%woody         = woodx	       ! binary flag for woody lifeform (1=woody, 0=not woody)
  pftcon%lflitcn       = lflitcx       ! leaf litter C:N (gC/gN)
  pftcon%frootcn       = frootcx       ! fine root C:N (gC/gN)
  pftcon%livewdcn      = livewdcx      ! live wood (phloem and ray parenchyma) C:N (gC/gN)
  pftcon%deadwdcn      = deadwdcx      ! dead wood (xylem and heartwood) C:N (gC/gN)
  pftcon%froot_leaf    = froot_leax    ! allocation parameter: new fine root C per new leaf C (gC/gC)
  pftcon%stem_leaf     = stem_leax     ! allocation parameter: new stem c per new leaf C (gC/gC)
  pftcon%croot_stem    = croot_stex    ! allocation parameter: new coarse root C per new stem C (gC/gC)
  pftcon%flivewd       = flivewx       ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  pftcon%fcur          = fcux	       ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  pftcon%lf_flab       = lf_flax       ! leaf litter labile fraction
  pftcon%lf_fcel       = lf_fcex       ! leaf litter cellulose fraction
  pftcon%lf_flig       = lf_flix       ! leaf litter lignin fraction
  pftcon%fr_flab       = fr_flax       ! fine root litter labile fraction
  pftcon%fr_fcel       = fr_fcex       ! fine root litter cellulose fraction
  pftcon%fr_flig       = fr_flix       ! fine root litter lignin fraction
  pftcon%leaf_long     = leaf_lonx     ! leaf longevity (yrs)
  pftcon%evergreen     = evergreex     ! binary flag for evergreen leaf habit (0 or 1)
  pftcon%stress_decid  = stress_decix  ! binary flag for stress-deciduous leaf habit (0 or 1)
  pftcon%season_decid  = season_decix  ! binary flag for seasonal-deciduous leaf habit (0 or 1)
  pftcon%xl            = xlx           ! leaf/stem orientation index
  pftcon%rhol          = rholx         ! leaf reflectance (visible)
  pftcon%rhos          = rhosx         ! stem reflectance (visible)
  pftcon%taul          = taulx         ! leaf transmittance (visible)
  pftcon%taus          = tausx         ! stem transmittance (visible)
  pftcon%cc_dstem      = cc_dstex      ! Combustion completeness factor for dead stem (0 to 1)
  pftcon%cc_leaf       = cc_leax       ! Combustion completeness factor for leaf (0 to 1) 
  pftcon%cc_lstem      = cc_lstex      ! Combustion completeness factor for live stem (0 to 1)
  pftcon%cc_other      = cc_othex      ! Combustion completeness factor for other plant (0 to 1)
  pftcon%fm_leaf       = fm_leax       ! Fire-related mortality factor for leaf (0 to 1)
  pftcon%fm_lstem      = fm_lstex      ! Fire-related mortality factor for live stem (0 to 1) 
  pftcon%fm_other      = fm_othex      ! Fire-related mortality factor for other plant (0 to 1)
  pftcon%fm_root       = fm_roox       ! Fire-related mortality factor for fine roots (0 to 1)
  pftcon%fm_lroot      = fm_lroox      ! Fire-related mortality factor for live roots (0 to 1) 
  pftcon%fm_droot      = fm_droox      ! Fire-related mortality factor for dead roots (0 to 1)
  pftcon%fd_pft        = fd_pfx        ! Fire duration (hr) 
  pftcon%fsr_pft       = fsr_pfx       ! Fire spread rate (m/s)
  pftcon%grperc        = grperx        ! Growth respiration factor (unitless)
  pftcon%grpnow        = grpnox        ! Growth respiration factor (unitless)
  pftcon%dwood         = 2.5e5         ! cn wood density (gC/m3); lpj:2.0e5; from CLM4.5 pftvarcon.F90. Values are same as in CLM4, fzeng.
  
! pftcon%fertnitro     = fertnitrx     ! Max fertilizer to be applied in total (kg N/m2)
! pftcon%lfemerg       = lfemerx       ! Leaf emergence parameter used in CNPhenology (unitless)
! pftcon%grnfill       = grnfilx       ! Grain fill parameter used in CNPhenology (unitless)
! pftcon%mxmat         = mxmax         ! Maximum number of days to maturity parameter in CNPhenology (days)
! pftcon%hybgdd        = hybgdx        ! Growing Degree Days for maturity used in CNPhenology (unitless)  
! pftcon%laimx         = laimxx        ! Maximum Leaf Area Index used in CNVegStructUpdate  
! pftcon%ztopmx        = ztopmxx       ! Canopy top coefficient used in CNVegStructUpdate (m)
! pftcon%gddmin        = gddmix        ! Minimim growing degree days used in CNPhenology 
! pftcon%graincn       = graincx       ! Grain C:N (gC/gN)
! pftcon%fleafcn       = fleafcx       ! Leaf C:N during organ fill (gC/gN) 
! pftcon%ffrootcn      = ffrootcx      ! Fine root C:N during organ fill (gC/gN) 
! pftcon%fstemcn       = fstemcx       ! Stem C:N during organ fill (gC/gN)
! pftcon%rootprof_beta = rootprof_betx ! Rooting beta parameter, for C and N vertical discretization (unitless)

  ! for CNAllocation
! pftcon%aleaff        = aleafx        ! Leaf Allocation coefficient parameter used in CNAllocation
! pftcon%allconsl      = allconslx     ! Leaf Allocation coefficient parameter power used in CNAllocation
! pftcon%allconss      = allconssx     ! Stem Allocation coefficient parameter power used in CNAllocation
! pftcon%arootf        = arootfx       ! Root Allocation coefficient parameter used in CNAllocation
! pftcon%arooti        = arootix       ! Root Allocation coefficient parameter used in CNAllocation
! pftcon%astemf        = astemx        ! Stem Allocation coefficient parameter used in CNAllocation
! pftcon%bfact         = bfacx         ! Exponential factor used in CNAllocation for fraction allocated to leaf
! pftcon%declfact      = declfacx      ! Decline factor for gddmaturity used in CNAllocation
! pftcon%fleafi        = fleafx        ! Leaf Allocation coefficient parameter fraction used in CNAllocation
  
! pftcon%baset         = basex         ! Base Temperature, parameter used in accFlds (degree C)   
! pftcon%mxtmp         = mxtmx         ! Max Temperature, parameter used in accFlds (degree C) 
! pftcon%planttemp     = planttemx     ! Average 10 day temperature needed for planting (K) 
! pftcon%minplanttemp  = minplanttemx  ! Average 5 day daily minimum temperature needed for planting (K)
! pftcon%mnNHplantdate = mnNHplantdatx ! Minimum planting date for the Northern Hemipsphere (MMDD) 
                                       ! Typical U.S. earliest planting dates according to AgroIBIS: Maize Apr 10th; soybean May 15th; spring wheat early Apr; winter wheat Sep 1st
! pftcon%mxNHplantdate = mxNHplantdatx ! Maximum planting date for the Northern Hemipsphere (MMDD)
                                       ! Typical U.S. latest planting dates according to AgroIBIS: Maize May 10th; soybean Jun 20th; spring wheat mid-May; winter wheat early Nov.  
! pftcon%mnSHplantdate = mnSHplantdatx ! Minimum planting date for the Southern Hemipsphere (MMDD), same as min_NH_planting_date, but offset by six months
! pftcon%mxSHplantdate = mxSHplantdatx ! Maximum planting date for the Southern Hemipsphere (MMDD), same as max_NH_planting_date, but offset by six months

! transfer restart vars from to CLM data structures if restart exists
! -------------------------------------------------------------------
  if(istep /= 0) then

  n = 0
  np = 0
  do nc = 1,nch        ! catchment tile loop
    do nz = 1,nzone    ! CN zone loop
      n = n + 1
      ccs%col_ctrunc_vr   (n,1)   = cncol(nc,nz, 1)
      ccs%decomp_cpools_vr(n,1,4) = cncol(nc,nz, 2) ! cwdc
      ccs%decomp_cpools_vr(n,1,1) = cncol(nc,nz, 3) ! litr1c
      ccs%decomp_cpools_vr(n,1,2) = cncol(nc,nz, 4) ! litr2c
      ccs%decomp_cpools_vr(n,1,3) = cncol(nc,nz, 5) ! litr3c
      ccs%totvegc_col     (n)     = cncol(nc,nz, 6)
      ccs%prod100c        (n)     = cncol(nc,nz, 7)
      ccs%prod10c         (n)     = cncol(nc,nz, 8)
      ccs%seedc           (n)     = cncol(nc,nz, 9)
      ccs%decomp_cpools_vr(n,1,5) = cncol(nc,nz,10) ! soil1c
      ccs%decomp_cpools_vr(n,1,6) = cncol(nc,nz,11) ! soil2c
      ccs%decomp_cpools_vr(n,1,7) = cncol(nc,nz,12) ! soil3c
      ccs%decomp_cpools_vr(n,1,8) = cncol(nc,nz,13) ! soil4c
      ccs%totcolc         (n)     = cncol(nc,nz,14)
      ccs%totlitc         (n)     = cncol(nc,nz,15)
      cns%col_ntrunc_vr   (n,1)   = cncol(nc,nz,16)
      cns%decomp_npools_vr(n,1,4) = cncol(nc,nz,17) ! cwdn
      cns%decomp_npools_vr(n,1,1) = cncol(nc,nz,18) ! litr1n
      cns%decomp_npools_vr(n,1,2) = cncol(nc,nz,19) ! litr2n
      cns%decomp_npools_vr(n,1,3) = cncol(nc,nz,20) ! litr3n
      cns%prod100n        (n)     = cncol(nc,nz,21)
      cns%prod10n         (n)     = cncol(nc,nz,22)
      cns%seedn           (n)     = cncol(nc,nz,23)
      cns%sminn_vr        (n,1)   = cncol(nc,nz,24)
      cns%decomp_npools_vr(n,1,5) = cncol(nc,nz,25) ! soil1n
      cns%decomp_npools_vr(n,1,6) = cncol(nc,nz,26) ! soil2n
      cns%decomp_npools_vr(n,1,7) = cncol(nc,nz,27) ! soil3n
      cns%decomp_npools_vr(n,1,8) = cncol(nc,nz,28) ! soil4n
      cns%totcoln         (n)     = cncol(nc,nz,29)
      cps%fpg             (n)     = cncol(nc,nz,30)
      cps%annsum_counter  (n)     = cncol(nc,nz,31)
      cps%cannavg_t2m     (n)     = cncol(nc,nz,32)
      cps%cannsum_npp     (n)     = cncol(nc,nz,33)
      cps%farea_burned    (n)     = cncol(nc,nz,34)
      cps%fpi_vr          (n,1)   = cncol(nc,nz,35)      

      ccs%decomp_cpools   (n,4)   = cncol(nc,nz, 2) ! cwdc
      ccs%decomp_cpools   (n,1)   = cncol(nc,nz, 3) ! litr1c
      ccs%decomp_cpools   (n,2)   = cncol(nc,nz, 4) ! litr2c
      ccs%decomp_cpools   (n,3)   = cncol(nc,nz, 5) ! litr3c
      ccs%decomp_cpools   (n,5)   = cncol(nc,nz,10) ! soil1c
      ccs%decomp_cpools   (n,6)   = cncol(nc,nz,11) ! soil2c
      ccs%decomp_cpools   (n,7)   = cncol(nc,nz,12) ! soil3c
      ccs%decomp_cpools   (n,8)   = cncol(nc,nz,13) ! soil4c

      ccs%decomp_cpools_1m(n,4)   = cncol(nc,nz, 2) ! cwdc
      ccs%decomp_cpools_1m(n,1)   = cncol(nc,nz, 3) ! litr1c
      ccs%decomp_cpools_1m(n,2)   = cncol(nc,nz, 4) ! litr2c
      ccs%decomp_cpools_1m(n,3)   = cncol(nc,nz, 5) ! litr3c
      ccs%decomp_cpools_1m(n,5)   = cncol(nc,nz,10) ! soil1c
      ccs%decomp_cpools_1m(n,6)   = cncol(nc,nz,11) ! soil2c
      ccs%decomp_cpools_1m(n,7)   = cncol(nc,nz,12) ! soil3c
      ccs%decomp_cpools_1m(n,8)   = cncol(nc,nz,13) ! soil4c
      
      cns%decomp_npools   (n,4)   = cncol(nc,nz,17) ! cwdn
      cns%decomp_npools   (n,1)   = cncol(nc,nz,18) ! litr1n
      cns%decomp_npools   (n,2)   = cncol(nc,nz,19) ! litr2n
      cns%decomp_npools   (n,3)   = cncol(nc,nz,20) ! litr3n    
      cns%decomp_npools   (n,5)   = cncol(nc,nz,25) ! soil1n
      cns%decomp_npools   (n,6)   = cncol(nc,nz,26) ! soil2n
      cns%decomp_npools   (n,7)   = cncol(nc,nz,27) ! soil3n
      cns%decomp_npools   (n,8)   = cncol(nc,nz,28) ! soil4n

      cns%decomp_npools_1m(n,4)   = cncol(nc,nz,17) ! cwdn
      cns%decomp_npools_1m(n,1)   = cncol(nc,nz,18) ! litr1n
      cns%decomp_npools_1m(n,2)   = cncol(nc,nz,19) ! litr2n
      cns%decomp_npools_1m(n,3)   = cncol(nc,nz,20) ! litr3n    
      cns%decomp_npools_1m(n,5)   = cncol(nc,nz,25) ! soil1n
      cns%decomp_npools_1m(n,6)   = cncol(nc,nz,26) ! soil2n
      cns%decomp_npools_1m(n,7)   = cncol(nc,nz,27) ! soil3n
      cns%decomp_npools_1m(n,8)   = cncol(nc,nz,28) ! soil4n

      ! total sminn, see CNSummaryMod.F90
      cns%sminn(n) = 0.
      do j = 1, nlevdecomp
          cns%sminn(n) = cns%sminn(n) + &
              cns%sminn_vr(n,j) * dzsoi_decomp(j)
      end do
      
      ! total soil organic matter carbon (TOTSOMC), see CNSummaryMod.F90
      ! compute the initial value of totsomc for CNFireFluxes, fzeng, 23 Aug 2018
      ccs%totsomc(n) = ccs%decomp_cpools(n,5) + ccs%decomp_cpools(n,6) + &
                       ccs%decomp_cpools(n,7) + ccs%decomp_cpools(n,8)

      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,nveg ! defined veg loop

          if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
           pcs%cpool                 (np) = cnpft(nc,nz,nv,  1)
           pcs%deadcrootc            (np) = cnpft(nc,nz,nv,  2)
           pcs%deadcrootc_storage    (np) = cnpft(nc,nz,nv,  3)
           pcs%deadcrootc_xfer       (np) = cnpft(nc,nz,nv,  4)
           pcs%deadstemc             (np) = cnpft(nc,nz,nv,  5)
           pcs%deadstemc_storage     (np) = cnpft(nc,nz,nv,  6)
           pcs%deadstemc_xfer        (np) = cnpft(nc,nz,nv,  7)
           pcs%frootc                (np) = cnpft(nc,nz,nv,  8)
           pcs%frootc_storage        (np) = cnpft(nc,nz,nv,  9)
           pcs%frootc_xfer           (np) = cnpft(nc,nz,nv, 10)
           pcs%gresp_storage         (np) = cnpft(nc,nz,nv, 11)
           pcs%gresp_xfer            (np) = cnpft(nc,nz,nv, 12)
           pcs%leafc                 (np) = cnpft(nc,nz,nv, 13)
           pcs%leafc_storage         (np) = cnpft(nc,nz,nv, 14) 
           pcs%leafc_xfer            (np) = cnpft(nc,nz,nv, 15)
           pcs%livecrootc            (np) = cnpft(nc,nz,nv, 16) 
           pcs%livecrootc_storage    (np) = cnpft(nc,nz,nv, 17)
           pcs%livecrootc_xfer       (np) = cnpft(nc,nz,nv, 18)
           pcs%livestemc             (np) = cnpft(nc,nz,nv, 19)
           pcs%livestemc_storage     (np) = cnpft(nc,nz,nv, 20)
           pcs%livestemc_xfer        (np) = cnpft(nc,nz,nv, 21)
           pcs%pft_ctrunc            (np) = cnpft(nc,nz,nv, 22)
           pcs%xsmrpool              (np) = cnpft(nc,nz,nv, 23)
           pepv%annavg_t2m           (np) = cnpft(nc,nz,nv, 24)
           pepv%annmax_retransn      (np) = cnpft(nc,nz,nv, 25)
           pepv%annsum_npp           (np) = cnpft(nc,nz,nv, 26)
           pepv%annsum_potential_gpp (np) = cnpft(nc,nz,nv, 27)
           pepv%dayl                 (np) = cnpft(nc,nz,nv, 28)
           pepv%days_active          (np) = cnpft(nc,nz,nv, 29)
           pepv%dormant_flag         (np) = cnpft(nc,nz,nv, 30)
           pepv%offset_counter       (np) = cnpft(nc,nz,nv, 31)
           pepv%offset_fdd           (np) = cnpft(nc,nz,nv, 32)
           pepv%offset_flag          (np) = cnpft(nc,nz,nv, 33)
           pepv%offset_swi           (np) = cnpft(nc,nz,nv, 34)
           pepv%onset_counter        (np) = cnpft(nc,nz,nv, 35)
           pepv%onset_fdd            (np) = cnpft(nc,nz,nv, 36)
           pepv%onset_flag           (np) = cnpft(nc,nz,nv, 37)
           pepv%onset_gdd            (np) = cnpft(nc,nz,nv, 38)
           pepv%onset_gddflag        (np) = cnpft(nc,nz,nv, 39)
           pepv%onset_swi            (np) = cnpft(nc,nz,nv, 40)
           pepv%prev_frootc_to_litter(np) = cnpft(nc,nz,nv, 41)
           pepv%prev_leafc_to_litter (np) = cnpft(nc,nz,nv, 42)
           pepv%tempavg_t2m          (np) = cnpft(nc,nz,nv, 43)
           pepv%tempmax_retransn     (np) = cnpft(nc,nz,nv, 44)
           pepv%tempsum_npp          (np) = cnpft(nc,nz,nv, 45)
           pepv%tempsum_potential_gpp(np) = cnpft(nc,nz,nv, 46)
           pepv%xsmrpool_recover     (np) = cnpft(nc,nz,nv, 47)
           pns%deadcrootn            (np) = cnpft(nc,nz,nv, 48)
           pns%deadcrootn_storage    (np) = cnpft(nc,nz,nv, 49)
           pns%deadcrootn_xfer       (np) = cnpft(nc,nz,nv, 50)
           pns%deadstemn             (np) = cnpft(nc,nz,nv, 51)
           pns%deadstemn_storage     (np) = cnpft(nc,nz,nv, 52)
           pns%deadstemn_xfer        (np) = cnpft(nc,nz,nv, 53)
           pns%frootn                (np) = cnpft(nc,nz,nv, 54)
           pns%frootn_storage        (np) = cnpft(nc,nz,nv, 55)
           pns%frootn_xfer           (np) = cnpft(nc,nz,nv, 56)
           pns%leafn                 (np) = cnpft(nc,nz,nv, 57)
           pns%leafn_storage         (np) = cnpft(nc,nz,nv, 58)
           pns%leafn_xfer            (np) = cnpft(nc,nz,nv, 59)
           pns%livecrootn            (np) = cnpft(nc,nz,nv, 60)
           pns%livecrootn_storage    (np) = cnpft(nc,nz,nv, 61)
           pns%livecrootn_xfer       (np) = cnpft(nc,nz,nv, 62)
           pns%livestemn             (np) = cnpft(nc,nz,nv, 63)
           pns%livestemn_storage     (np) = cnpft(nc,nz,nv, 64)
           pns%livestemn_xfer        (np) = cnpft(nc,nz,nv, 65)
           pns%npool                 (np) = cnpft(nc,nz,nv, 66)
           pns%pft_ntrunc            (np) = cnpft(nc,nz,nv, 67)
           pns%retransn              (np) = cnpft(nc,nz,nv, 68)
           pps%elai                  (np) = cnpft(nc,nz,nv, 69)
           pps%esai                  (np) = cnpft(nc,nz,nv, 70)
           pps%hbot                  (np) = cnpft(nc,nz,nv, 71)
           pps%htop                  (np) = cnpft(nc,nz,nv, 72)
           pps%tlai                  (np) = cnpft(nc,nz,nv, 73)
           pps%tsai                  (np) = cnpft(nc,nz,nv, 74)
           pepv%plant_ndemand        (np) = cnpft(nc,nz,nv, 75)
           
	   pcs%totvegc               (np) = pcs%cpool                 (np) + &
                                            pcs%deadcrootc            (np) + &
                                            pcs%deadcrootc_storage    (np) + &
                                            pcs%deadcrootc_xfer       (np) + &
                                            pcs%deadstemc             (np) + &
                                            pcs%deadstemc_storage     (np) + &
                                            pcs%deadstemc_xfer        (np) + &
                                            pcs%frootc                (np) + &
                                            pcs%frootc_storage        (np) + &
                                            pcs%frootc_xfer           (np) + &
                                            pcs%gresp_storage         (np) + &
                                            pcs%gresp_xfer            (np) + &
                                            pcs%leafc                 (np) + &
                                            pcs%leafc_storage         (np) + &
                                            pcs%leafc_xfer            (np) + &
                                            pcs%livecrootc            (np) + &
                                            pcs%livecrootc_storage    (np) + &
                                            pcs%livecrootc_xfer       (np) + &
                                            pcs%livestemc             (np) + &
                                            pcs%livestemc_storage     (np) + &
                                            pcs%livestemc_xfer        (np) 
          
	   ! Make it a cold start if there is no information from the restart file (derived from CLM4 restart)
	   ! fzeng, 1 Aug 2017
	   if (isnan(pcs%leafc_xfer(np))) pcs%leafc_xfer(np) = 0.                                                     
          
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
      cncol(nc,nz, 1) = ccs%col_ctrunc_vr   (n,1)   
      cncol(nc,nz, 2) = ccs%decomp_cpools_vr(n,1,4)  ! cwdc
      cncol(nc,nz, 3) = ccs%decomp_cpools_vr(n,1,1)  ! litr1c
      cncol(nc,nz, 4) = ccs%decomp_cpools_vr(n,1,2)  ! litr2c
      cncol(nc,nz, 5) = ccs%decomp_cpools_vr(n,1,3)  ! litr3c
      cncol(nc,nz, 6) = ccs%totvegc_col     (n)     
      cncol(nc,nz, 7) = ccs%prod100c        (n)     
      cncol(nc,nz, 8) = ccs%prod10c         (n)     
      cncol(nc,nz, 9) = ccs%seedc           (n)     
      cncol(nc,nz,10) = ccs%decomp_cpools_vr(n,1,5)  ! soil1c
      cncol(nc,nz,11) = ccs%decomp_cpools_vr(n,1,6)  ! soil2c
      cncol(nc,nz,12) = ccs%decomp_cpools_vr(n,1,7)  ! soil3c
      cncol(nc,nz,13) = ccs%decomp_cpools_vr(n,1,8)  ! soil4c
      cncol(nc,nz,14) = ccs%totcolc         (n)     
      cncol(nc,nz,15) = ccs%totlitc         (n)     
      cncol(nc,nz,16) = cns%col_ntrunc_vr   (n,1)   
      cncol(nc,nz,17) = cns%decomp_npools_vr(n,1,4)  ! cwdn
      cncol(nc,nz,18) = cns%decomp_npools_vr(n,1,1)  ! litr1n
      cncol(nc,nz,19) = cns%decomp_npools_vr(n,1,2)  ! litr2n
      cncol(nc,nz,20) = cns%decomp_npools_vr(n,1,3)  ! litr3n
      cncol(nc,nz,21) = cns%prod100n        (n)     
      cncol(nc,nz,22) = cns%prod10n         (n)     
      cncol(nc,nz,23) = cns%seedn           (n)     
      cncol(nc,nz,24) = cns%sminn_vr        (n,1)   
      cncol(nc,nz,25) = cns%decomp_npools_vr(n,1,5)  ! soil1n
      cncol(nc,nz,26) = cns%decomp_npools_vr(n,1,6)  ! soil2n
      cncol(nc,nz,27) = cns%decomp_npools_vr(n,1,7)  ! soil3n
      cncol(nc,nz,28) = cns%decomp_npools_vr(n,1,8)  ! soil4n
      cncol(nc,nz,29) = cns%totcoln         (n)     
      cncol(nc,nz,30) = cps%fpg             (n)     
      cncol(nc,nz,31) = cps%annsum_counter  (n)     
      cncol(nc,nz,32) = cps%cannavg_t2m     (n)     
      cncol(nc,nz,33) = cps%cannsum_npp     (n)     
      cncol(nc,nz,34) = cps%farea_burned    (n)              
      cncol(nc,nz,35) = cps%fpi_vr          (n,1)       

      do p = 0,numpft  ! PFT index loop
        np = np + 1
        do nv = 1,nveg ! defined veg loop

          if(ityp(nc,nv,nz)==p .and. fveg(nc,nv,nz)>1.e-4) then
            cnpft(nc,nz,nv,  1) = pcs%cpool                 (np)
            cnpft(nc,nz,nv,  2) = pcs%deadcrootc            (np)
            cnpft(nc,nz,nv,  3) = pcs%deadcrootc_storage    (np)
            cnpft(nc,nz,nv,  4) = pcs%deadcrootc_xfer       (np)
            cnpft(nc,nz,nv,  5) = pcs%deadstemc             (np)
            cnpft(nc,nz,nv,  6) = pcs%deadstemc_storage     (np)
            cnpft(nc,nz,nv,  7) = pcs%deadstemc_xfer        (np)
            cnpft(nc,nz,nv,  8) = pcs%frootc                (np)
            cnpft(nc,nz,nv,  9) = pcs%frootc_storage        (np)
            cnpft(nc,nz,nv, 10) = pcs%frootc_xfer           (np)
            cnpft(nc,nz,nv, 11) = pcs%gresp_storage         (np)
            cnpft(nc,nz,nv, 12) = pcs%gresp_xfer            (np)
            cnpft(nc,nz,nv, 13) = pcs%leafc                 (np)
            cnpft(nc,nz,nv, 14) = pcs%leafc_storage         (np) 
            cnpft(nc,nz,nv, 15) = pcs%leafc_xfer            (np)
            cnpft(nc,nz,nv, 16) = pcs%livecrootc            (np) 
            cnpft(nc,nz,nv, 17) = pcs%livecrootc_storage    (np)
            cnpft(nc,nz,nv, 18) = pcs%livecrootc_xfer       (np)
            cnpft(nc,nz,nv, 19) = pcs%livestemc             (np)
            cnpft(nc,nz,nv, 20) = pcs%livestemc_storage     (np)
            cnpft(nc,nz,nv, 21) = pcs%livestemc_xfer        (np)
            cnpft(nc,nz,nv, 22) = pcs%pft_ctrunc            (np)
            cnpft(nc,nz,nv, 23) = pcs%xsmrpool              (np)
            cnpft(nc,nz,nv, 24) = pepv%annavg_t2m           (np)
            cnpft(nc,nz,nv, 25) = pepv%annmax_retransn      (np)
            cnpft(nc,nz,nv, 26) = pepv%annsum_npp           (np)
            cnpft(nc,nz,nv, 27) = pepv%annsum_potential_gpp (np)
            cnpft(nc,nz,nv, 28) = pepv%dayl                 (np)
            cnpft(nc,nz,nv, 29) = pepv%days_active          (np)
            cnpft(nc,nz,nv, 30) = pepv%dormant_flag         (np)
            cnpft(nc,nz,nv, 31) = pepv%offset_counter       (np)
            cnpft(nc,nz,nv, 32) = pepv%offset_fdd           (np)
            cnpft(nc,nz,nv, 33) = pepv%offset_flag          (np)
            cnpft(nc,nz,nv, 34) = pepv%offset_swi           (np)
            cnpft(nc,nz,nv, 35) = pepv%onset_counter        (np)
            cnpft(nc,nz,nv, 36) = pepv%onset_fdd            (np)
            cnpft(nc,nz,nv, 37) = pepv%onset_flag           (np)
            cnpft(nc,nz,nv, 38) = pepv%onset_gdd            (np)
            cnpft(nc,nz,nv, 39) = pepv%onset_gddflag        (np)
            cnpft(nc,nz,nv, 40) = pepv%onset_swi            (np)
            cnpft(nc,nz,nv, 41) = pepv%prev_frootc_to_litter(np)
            cnpft(nc,nz,nv, 42) = pepv%prev_leafc_to_litter (np)
            cnpft(nc,nz,nv, 43) = pepv%tempavg_t2m          (np)
            cnpft(nc,nz,nv, 44) = pepv%tempmax_retransn     (np)
            cnpft(nc,nz,nv, 45) = pepv%tempsum_npp          (np)
            cnpft(nc,nz,nv, 46) = pepv%tempsum_potential_gpp(np)
            cnpft(nc,nz,nv, 47) = pepv%xsmrpool_recover     (np)
            cnpft(nc,nz,nv, 48) = pns%deadcrootn            (np)
            cnpft(nc,nz,nv, 49) = pns%deadcrootn_storage    (np)
            cnpft(nc,nz,nv, 50) = pns%deadcrootn_xfer       (np)
            cnpft(nc,nz,nv, 51) = pns%deadstemn             (np)
            cnpft(nc,nz,nv, 52) = pns%deadstemn_storage     (np)
            cnpft(nc,nz,nv, 53) = pns%deadstemn_xfer        (np)
            cnpft(nc,nz,nv, 54) = pns%frootn                (np)
            cnpft(nc,nz,nv, 55) = pns%frootn_storage        (np)
            cnpft(nc,nz,nv, 56) = pns%frootn_xfer           (np)
            cnpft(nc,nz,nv, 57) = pns%leafn                 (np)
            cnpft(nc,nz,nv, 58) = pns%leafn_storage         (np)
            cnpft(nc,nz,nv, 59) = pns%leafn_xfer            (np)
            cnpft(nc,nz,nv, 60) = pns%livecrootn            (np)
            cnpft(nc,nz,nv, 61) = pns%livecrootn_storage    (np)
            cnpft(nc,nz,nv, 62) = pns%livecrootn_xfer       (np)
            cnpft(nc,nz,nv, 63) = pns%livestemn             (np)
            cnpft(nc,nz,nv, 64) = pns%livestemn_storage     (np)
            cnpft(nc,nz,nv, 65) = pns%livestemn_xfer        (np)
            cnpft(nc,nz,nv, 66) = pns%npool                 (np)
            cnpft(nc,nz,nv, 67) = pns%pft_ntrunc            (np)
            cnpft(nc,nz,nv, 68) = pns%retransn              (np)
            cnpft(nc,nz,nv, 69) = pps%elai                  (np)
            cnpft(nc,nz,nv, 70) = pps%esai                  (np)
            cnpft(nc,nz,nv, 71) = pps%hbot                  (np)
            cnpft(nc,nz,nv, 72) = pps%htop                  (np)
            cnpft(nc,nz,nv, 73) = pps%tlai                  (np)
            cnpft(nc,nz,nv, 74) = pps%tsai                  (np)           
            cnpft(nc,nz,nv, 75) = pepv%plant_ndemand        (np)
            
          endif

        end do ! defined veg loop
      end do   ! PFT index loop
    end do     ! CN zone loop
  end do       ! catchment tile loop

  return

  end subroutine CN_exit

  subroutine get_CN_LAI(nch,nveg,nzone,ityp,fveg,elai,esai,tlai,tsai)

  integer, intent(in) :: nch ! number of tiles
  integer, intent(in) :: nveg ! number of vegetation types per zone
  integer, intent(in) :: nzone ! number of stress zones per tile
  integer, dimension(nch,nveg,nzone), intent(in) :: ityp ! PFT index
  real, dimension(nch,nveg,nzone), intent(in) :: fveg    ! PFT fraction
  real, dimension(nch,nveg,nzone), intent(out)           :: elai   ! exposed leaf-area index
  real, dimension(nch,nveg,nzone), intent(out), optional :: esai   ! exposed stem-area index
  real, dimension(nch,nveg,nzone), intent(out), optional :: tlai   ! total leaf-area index
  real, dimension(nch,nveg,nzone), intent(out), optional :: tsai   ! total stem-area index

  integer :: n, p, nv, nc, nz, np

                    elai = 0.
  if(present(esai)) esai = 0.
  if(present(tlai)) tlai = 0.
  if(present(tsai)) tsai = 0.

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
                              elai(nc,nv,nz) = pps%elai(np)
            if(present(esai)) esai(nc,nv,nz) = pps%esai(np)
            if(present(tlai)) tlai(nc,nv,nz) = pps%tlai(np)
            if(present(tsai)) tsai(nc,nv,nz) = pps%tsai(np)
          endif

        end do ! defined veg loop
      end do   ! PFT index loop
    end do     ! CN zone loop
  end do       ! catchment tile loop

  if ( .not. LAND_FIX ) then ! jkolassa Oct 2020: the if-wrapper here is to toggle between the LDASsa version used by Fanwei Zeng and Eunjee Lee and current GEOSldas Catchment-CN; there is likely a better way to control this
     where (elai > 20.) elai = 20.
     where (esai > 20.) esai = 20.
  end if 

  end subroutine get_CN_LAI 

end module CN_DriverMod
