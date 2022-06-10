module CNCLM_pftconMod

  use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
  use nanMod           , only : nan
  use clm_varpar       , only : mxpft, numrad
  use netcdf 
  use MAPL_ExceptionHandling


  ! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_pftcon_type

!
! Vegetation type constants
!
  integer, public :: noveg                 =  0   ! Bare                                 
  integer, public :: ndllf_evr_tmp_tree    =  1   ! Needleleaf evergreen temperate tree           
  integer, public :: ndllf_evr_brl_tree    =  2   ! Needleleaf evergreen boreal tree             
  integer, public :: ndllf_dcd_brl_tree    =  3   ! Needleleaf deciduous boreal tree              
  integer, public :: nbrdlf_evr_trp_tree   =  4   ! Broadleaf evergreen tropical tree            
  integer, public :: nbrdlf_evr_tmp_tree   =  5   ! Broadleaf evergreen temperate tree            
  integer, public :: nbrdlf_dcd_trp_tree   =  6   ! Broadleaf deciduous tropical tree             
  integer, public :: nbrdlf_dcd_tmp_tree   =  7   ! Broadleaf deciduous temperate tree            
  integer, public :: nbrdlf_dcd_brl_tree   =  8   ! Broadleaf deciduous boreal tree               
  integer, public :: nbrdlf_evr_shrub      =  9   ! Broadleaf evergreen temperate shrub                     
  integer, public :: nbrdlf_dcd_tmp_shrub  = 10   ! Broadleaf deciduous temperate shrub [moisture + deciduous]
  integer, public :: nbrdlf_dcd_brl_shrub  = 11   ! Broadleaf deciduous boreal shrub            
  integer, public :: nc3_arctic_grass      = 12   ! Arctic c3 grass                              
  integer, public :: nc3_nonarctic_grass   = 13   ! Cool c3 grass [moisture + deciduous]
  integer, public :: nc4_grass             = 14   ! Warm c4 grass [moisture + deciduous]
  integer, public :: nc3crop               = 15   ! C3_crop [moisture + deciduous]   
  integer, public :: npcropmin = nc3crop          ! value for first crop

 !
  type, public :: pftcon_type

     integer , allocatable :: noveg         (:)   ! value for not vegetated
     logical , allocatable :: is_tree       (:)   ! tree or not?
     logical , allocatable :: is_shrub      (:)   ! shrub or not?
     logical , allocatable :: is_grass      (:)   ! grass or not?

     real(r8), allocatable :: dleaf         (:)   ! characteristic leaf dimension (m)
     real(r8), allocatable :: c3psn         (:)   ! photosynthetic pathway: 0. = c4, 1. = c3
     real(r8), allocatable :: xl            (:)   ! leaf/stem orientation index
     real(r8), allocatable :: rhol          (:,:) ! leaf reflectance: 1=vis, 2=nir
     real(r8), allocatable :: rhos          (:,:) ! stem reflectance: 1=vis, 2=nir
     real(r8), allocatable :: taul          (:,:) ! leaf transmittance: 1=vis, 2=nir
     real(r8), allocatable :: taus          (:,:) ! stem transmittance: 1=vis, 2=nir
     real(r8), allocatable :: z0mr          (:)   ! ratio of momentum roughness length to canopy top height (-)
     real(r8), allocatable :: displar       (:)   ! ratio of displacement height to canopy top height (-)
     real(r8), allocatable :: roota_par     (:)   ! CLM rooting distribution parameter [1/m]
     real(r8), allocatable :: rootb_par     (:)   ! CLM rooting distribution parameter [1/m]
     real(r8), allocatable :: crop          (:)   ! crop pft: 0. = not crop, 1. = crop pft
     real(r8), allocatable :: irrigated     (:)   ! irrigated pft: 0. = not, 1. = irrigated
     real(r8), allocatable :: smpso         (:)   ! soil water potential at full stomatal opening (mm)
     real(r8), allocatable :: smpsc         (:)   ! soil water potential at full stomatal closure (mm)
     real(r8), allocatable :: fnitr         (:)   ! foliage nitrogen limitation factor (-)

     !  CN code
     real(r8), allocatable :: dwood         (:)   ! wood density (gC/m3)
     real(r8), allocatable :: slatop        (:)   ! SLA at top of canopy [m^2/gC]
     real(r8), allocatable :: dsladlai      (:)   ! dSLA/dLAI [m^2/gC]
     real(r8), allocatable :: leafcn        (:)   ! leaf C:N [gC/gN]
     real(r8), allocatable :: biofuel_harvfrac (:) ! fraction of stem and leaf cut for harvest, sent to biofuels [unitless]
     real(r8), allocatable :: flnr          (:)   ! fraction of leaf N in Rubisco [no units]
     real(r8), allocatable :: woody         (:)   ! woody lifeform flag (0 or 1)
     real(r8), allocatable :: lflitcn       (:)   ! leaf litter C:N (gC/gN)
     real(r8), allocatable :: frootcn       (:)   ! fine root C:N (gC/gN)
     real(r8), allocatable :: livewdcn      (:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
     real(r8), allocatable :: deadwdcn      (:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
     real(r8), allocatable :: grperc        (:)   ! growth respiration parameter
     real(r8), allocatable :: grpnow        (:)   ! growth respiration parameter
     real(r8), allocatable :: rootprof_beta (:,:) ! CLM rooting distribution parameter for C and N inputs [unitless]
     real(r8), allocatable :: root_radius   (:)   ! root radius (m)
     real(r8), allocatable :: root_density  (:)   ! root density (gC/m3)

     real(r8), allocatable :: dbh  (:)            ! diameter at breast height (m)
     real(r8), allocatable :: fbw  (:)            ! fraction of biomass that is water
     real(r8), allocatable :: nstem  (:)          ! stem density (#/m2)
     real(r8), allocatable :: taper  (:)          ! tapering ratio of height:radius_breast_height
     real(r8), allocatable :: rstem_per_dbh  (:)  ! stem resistance per dbh (s/m/m)
     real(r8), allocatable :: wood_density  (:)   ! wood density (kg/m3)

     !  crop

     ! These arrays give information about the merge of unused crop types to the types CLM
     ! knows about. mergetoclmpft(m) gives the crop type that CLM uses to simulate input
     ! type m (and mergetoclmpft(m) == m implies that CLM simulates crop type m
     ! directly). is_pft_known_to_model(m) is true if CLM simulates crop type m, and false
     ! otherwise. Note that these do NOT relate to whether irrigation is on or off in a
     ! given simulation - that is handled separately.
     integer , allocatable :: mergetoclmpft         (:)
     logical , allocatable :: is_pft_known_to_model (:)

     real(r8), allocatable :: graincn       (:)   ! grain C:N (gC/gN)
     real(r8), allocatable :: mxtmp         (:)   ! parameter used in accFlds
     real(r8), allocatable :: baset         (:)   ! parameter used in accFlds
     real(r8), allocatable :: declfact      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: bfact         (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: aleaff        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: arootf        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: astemf        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: arooti        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: fleafi        (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: allconsl      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: allconss      (:)   ! parameter used in CNAllocation
     real(r8), allocatable :: ztopmx        (:)   ! parameter used in CNVegStructUpdate
     real(r8), allocatable :: laimx         (:)   ! parameter used in CNVegStructUpdate
     real(r8), allocatable :: gddmin        (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: hybgdd        (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: lfemerg       (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: grnfill       (:)   ! parameter used in CNPhenology
     integer , allocatable :: mxmat         (:)   ! parameter used in CNPhenology
     real(r8), allocatable :: mbbopt        (:)   ! Ball-Berry equation slope used in Photosynthesis
     real(r8), allocatable :: medlynslope   (:)   ! Medlyn equation slope used in Photosynthesis
     real(r8), allocatable :: medlynintercept(:)  ! Medlyn equation intercept used in Photosynthesis
     integer , allocatable :: mnNHplantdate (:)   ! minimum planting date for NorthHemisphere (YYYYMMDD)
     integer , allocatable :: mxNHplantdate (:)   ! maximum planting date for NorthHemisphere (YYYYMMDD)
     integer , allocatable :: mnSHplantdate (:)   ! minimum planting date for SouthHemisphere (YYYYMMDD)
     integer , allocatable :: mxSHplantdate (:)   ! maximum planting date for SouthHemisphere (YYYYMMDD)
     real(r8), allocatable :: planttemp     (:)   ! planting temperature used in CNPhenology (K)
     real(r8), allocatable :: minplanttemp  (:)   ! mininum planting temperature used in CNPhenology (K)
     real(r8), allocatable :: froot_leaf    (:)   ! allocation parameter: new fine root C per new leaf C (gC/gC) 
     real(r8), allocatable :: stem_leaf     (:)   ! allocation parameter: new stem c per new leaf C (gC/gC)
     real(r8), allocatable :: croot_stem    (:)   ! allocation parameter: new coarse root C per new stem C (gC/gC)
     real(r8), allocatable :: flivewd       (:)   ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
     real(r8), allocatable :: fcur          (:)   ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
     real(r8), allocatable :: fcurdv        (:)   ! alternate fcur for use with cndv
     real(r8), allocatable :: lf_flab       (:)   ! leaf litter labile fraction
     real(r8), allocatable :: lf_fcel       (:)   ! leaf litter cellulose fraction
     real(r8), allocatable :: lf_flig       (:)   ! leaf litter lignin fraction
     real(r8), allocatable :: fr_flab       (:)   ! fine root litter labile fraction
     real(r8), allocatable :: fr_fcel       (:)   ! fine root litter cellulose fraction
     real(r8), allocatable :: fr_flig       (:)   ! fine root litter lignin fraction
     real(r8), allocatable :: leaf_long     (:)   ! leaf longevity (yrs)
     real(r8), allocatable :: evergreen     (:)   ! binary flag for evergreen leaf habit (0 or 1)
     real(r8), allocatable :: stress_decid  (:)   ! binary flag for stress-deciduous leaf habit (0 or 1)
     real(r8), allocatable :: season_decid  (:)   ! binary flag for seasonal-deciduous leaf habit (0 or 1)
!KO
     real(r8), allocatable :: season_decid_temperate(:) ! binary flag for seasonal-deciduous temperate leaf habit (0 or 1)
!KO
     real(r8), allocatable :: pconv         (:)   ! proportion of deadstem to conversion flux
     real(r8), allocatable :: pprod10       (:)   ! proportion of deadstem to 10-yr product pool
     real(r8), allocatable :: pprod100      (:)   ! proportion of deadstem to 100-yr product pool
     real(r8), allocatable :: pprodharv10   (:)   ! harvest mortality proportion of deadstem to 10-yr pool

     ! pft paraemeters for fire code
     real(r8), allocatable :: cc_leaf       (:)
     real(r8), allocatable :: cc_lstem      (:)
     real(r8), allocatable :: cc_dstem      (:)
     real(r8), allocatable :: cc_other      (:)
     real(r8), allocatable :: fm_leaf       (:)
     real(r8), allocatable :: fm_lstem      (:)
     real(r8), allocatable :: fm_dstem      (:)
     real(r8), allocatable :: fm_other      (:)
     real(r8), allocatable :: fm_root       (:)
     real(r8), allocatable :: fm_lroot      (:)
     real(r8), allocatable :: fm_droot      (:)
     real(r8), allocatable :: fsr_pft       (:)
     real(r8), allocatable :: fd_pft        (:)
     real(r8), allocatable :: rswf_min      (:)
     real(r8), allocatable :: rswf_max      (:)

     ! pft parameters for crop code
     real(r8), allocatable :: manunitro     (:)   ! manure
     real(r8), allocatable :: fleafcn       (:)   ! C:N during grain fill; leaf
     real(r8), allocatable :: ffrootcn      (:)   ! C:N during grain fill; fine root
     real(r8), allocatable :: fstemcn       (:)   ! C:N during grain fill; stem

     real(r8), allocatable :: i_vcad        (:)
     real(r8), allocatable :: s_vcad        (:)
     real(r8), allocatable :: i_flnr        (:)
     real(r8), allocatable :: s_flnr        (:)

     ! pft parameters for CNDV code (from LPJ subroutine pftparameters)
     real(r8), allocatable :: pftpar20      (:)   ! tree maximum crown area (m2)
     real(r8), allocatable :: pftpar28      (:)   ! min coldest monthly mean temperature
     real(r8), allocatable :: pftpar29      (:)   ! max coldest monthly mean temperature
     real(r8), allocatable :: pftpar30      (:)   ! min growing degree days (>= 5 deg C)
     real(r8), allocatable :: pftpar31      (:)   ! upper limit of temperature of the warmest month (twmax)

     ! pft parameters for FUN
     real(r8), allocatable :: a_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: b_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: c_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: s_fix         (:)   ! A BNF parameter
     real(r8), allocatable :: akc_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: akn_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: ekc_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: ekn_active    (:)   ! A mycorrhizal uptake parameter
     real(r8), allocatable :: kc_nonmyc     (:)   ! A non-mycorrhizal uptake parameter
     real(r8), allocatable :: kn_nonmyc     (:)   ! A non-mycorrhizal uptake parameter
     real(r8), allocatable :: kr_resorb     (:)   ! A retrasnlcation parameter
     real(r8), allocatable :: perecm        (:)   ! The fraction of ECM-associated PFT 
     real(r8), allocatable :: fun_cn_flex_a (:)   ! Parameter a of FUN-flexcn link code (def 5)
     real(r8), allocatable :: fun_cn_flex_b (:)   ! Parameter b of FUN-flexcn link code (def 200)
     real(r8), allocatable :: fun_cn_flex_c (:)   ! Parameter b of FUN-flexcn link code (def 80)         
     real(r8), allocatable :: FUN_fracfixers(:)   ! Fraction of C that can be used for fixation.    


     ! pft parameters for dynamic root code
     real(r8), allocatable :: root_dmx(:)     !maximum root depth

  end type pftcon_type

type(pftcon_type), public, target, save :: pftcon

contains

!--------------------------------
  subroutine init_pftcon_type(this)

  ! !DESCRIPTION:
! Initialize CTSM PFT constants                                
!                                                                                                                       
! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    type(pftcon_type), intent(inout):: this

    !LOCAL
    character(300)     :: paramfile
    integer            :: ierr, clm_varid

    real(r8), allocatable, dimension(:)   :: read_tmp_1
    real(r8), allocatable, dimension(:,:) :: read_tmp_2
    integer , allocatable, dimension(:)   :: read_tmp_3

!---------------------------------------------------------

    allocate( read_tmp_1         (0:78)) 
    allocate( read_tmp_2         (0:78,nvariants))
    allocate( read_tmp_3         (0:78))   

    allocate( this%noveg         (0:mxpft)); this%noveg    (:) = huge(1)
    allocate( this%is_tree       (0:mxpft)); this%is_tree  (:) = .false.
    allocate( this%is_shrub      (0:mxpft)); this%is_shrub (:) = .false.
    allocate( this%is_grass      (0:mxpft)); this%is_grass (:) = .false.

    allocate( this%dleaf         (0:mxpft) ); this%dleaf   (:) = nan  !#
    allocate( this%c3psn         (0:mxpft) ); this%c3psn   (:) = nan
    allocate( this%xl            (0:mxpft) ); this%xl      (:) = nan
    allocate( this%rhol          (0:mxpft,numrad) );  this%rhol (:,:) = nan
    allocate( this%rhos          (0:mxpft,numrad) );  this%rhos (:,:) = nan
    allocate( this%taul          (0:mxpft,numrad) );  this%taul (:,:) = nan
    allocate( this%taus          (0:mxpft,numrad) );  this%taus (:,:) = nan
    allocate( this%z0mr          (0:mxpft) ); this%z0mr     (:) = nan
    allocate( this%displar       (0:mxpft) ); this%displar  (:) = nan 
    allocate( this%roota_par     (0:mxpft) ); this%roota_par(:) = nan
    allocate( this%rootb_par     (0:mxpft) ); this%rootb_par(:) = nan
    allocate( this%crop          (0:mxpft) ); this%crop     (:) = nan  !#
    allocate( this%mergetoclmpft (0:mxpft) ); this%mergetoclmpft (:) = nan !#
    allocate( this%is_pft_known_to_model  (0:mxpft) ); this%is_pft_known_to_model(:) = nan !#
    allocate( this%irrigated     (0:mxpft) ); this%irrigated (:) = nan   !#
    allocate( this%smpso         (0:mxpft) ); this%smpso     (:) = nan   !#
    allocate( this%smpsc         (0:mxpft) ); this%smpsc     (:) = nan   !#
    allocate( this%fnitr         (0:mxpft) ); this%fnitr     (:) = nan   !#
    allocate( this%slatop        (0:mxpft) ); this%slatop    (:) = nan   
    allocate( this%dsladlai      (0:mxpft) ); this%dsladlai  (:) = nan   
    allocate( this%leafcn        (0:mxpft) ); this%leafcn    (:) = nan   
    allocate( this%biofuel_harvfrac (0:mxpft) ); this%biofuel_harvfrac(:) = nan !#
    allocate( this%flnr          (0:mxpft) ); this%flnr      (:) = nan   
    allocate( this%woody         (0:mxpft) ); this%woody     (:) = nan
    allocate( this%lflitcn       (0:mxpft) ); this%lflitcn   (:) = nan
    allocate( this%frootcn       (0:mxpft) ); this%frootcn   (:) = nan
    allocate( this%livewdcn      (0:mxpft) ); this%livewdcn  (:) = nan
    allocate( this%deadwdcn      (0:mxpft) ); this%deadwdcn  (:) = nan
    allocate( this%grperc        (0:mxpft) ); this%grperc    (:) = nan  
    allocate( this%grpnow        (0:mxpft) ); this%grpnow    (:) = nan  
    allocate( this%rootprof_beta (0:mxpft,nvariants) ); this%rootprof_beta(:,:) = nan
    allocate( this%graincn       (0:mxpft) ); this%graincn   (:) = nan  
    allocate( this%mxtmp         (0:mxpft) ); this%mxtmp     (:) = nan  
    allocate( this%baset         (0:mxpft) ); this%baset     (:) = nan  
    allocate( this%declfact      (0:mxpft) ); this%declfact  (:) = nan
    allocate( this%bfact         (0:mxpft) ); this%bfact     (:) = nan
    allocate( this%aleaff        (0:mxpft) ); this%aleaff    (:) = nan 
    allocate( this%arootf        (0:mxpft) ); this%arootf    (:) = nan
    allocate( this%astemf        (0:mxpft) ); this%astemf    (:) = nan 
    allocate( this%arooti        (0:mxpft) ); this%arooti    (:) = nan
    allocate( this%fleafi        (0:mxpft) ); this%fleafi    (:) = nan
    allocate( this%allconsl      (0:mxpft) ); this%allconsl  (:) = nan
    allocate( this%allconss      (0:mxpft) ); this%allconss  (:) = nan
    allocate( this%ztopmx        (0:mxpft) ); this%ztopmx    (:) = nan
    allocate( this%laimx         (0:mxpft) ); this%laimx     (:) = nan
    allocate( this%gddmin        (0:mxpft) ); this%gddmin    (:) = nan
    allocate( this%hybgdd        (0:mxpft) ); this%hybgdd    (:) = nan
    allocate( this%lfemerg       (0:mxpft) ); this%lfemerg   (:) = nan
    allocate( this%grnfill       (0:mxpft) ); this%grnfill   (:) = nan
    allocate( this%mbbopt        (0:mxpft) ); this%mbbopt    (:) = nan !#
    allocate( this%medlynslope   (0:mxpft) ); this%medlynslope (:) = nan !#
    allocate( this%medlynintercept(0:mxpft) ); this%medlynintercept = nan !#
    allocate( this%mxmat         (0:mxpft) ); this%mxmat     (:) = nan
    allocate( this%mnNHplantdate (0:mxpft) ); this%mnNHplantdate (:) = huge(1)
    allocate( this%mxNHplantdate (0:mxpft) ); this%mxNHplantdate (:) = huge(1)
    allocate( this%mnSHplantdate (0:mxpft) ); this%mnSHplantdate (:) = huge(1)
    allocate( this%mxSHplantdate (0:mxpft) ); this%mxSHplantdate (:) = huge(1)
    allocate( this%planttemp     (0:mxpft) ); this%planttemp     (:) = nan
    allocate( this%minplanttemp  (0:mxpft) ); this%minplanttemp  (:) = nan
    allocate( this%froot_leaf    (0:mxpft) ); this%froot_leaf    (:) = nan
    allocate( this%stem_leaf     (0:mxpft) ); this%stem_leaf     (:) = nan
    allocate( this%croot_stem    (0:mxpft) ); this%croot_stem    (:) = nan
    allocate( this%flivewd       (0:mxpft) ); this%flivewd       (:) = nan
    allocate( this%fcur          (0:mxpft) ); this%fcur          (:) = nan
    allocate( this%fcurdv        (0:mxpft) ); this%fcurdv        (:) = nan   !#
    allocate( this%lf_flab       (0:mxpft) ); this%lf_flab       (:) = nan
    allocate( this%lf_fcel       (0:mxpft) ); this%lf_fcel       (:) = nan
    allocate( this%lf_flig       (0:mxpft) ); this%lf_flig       (:) = nan
    allocate( this%fr_flab       (0:mxpft) ); this%fr_flab       (:) = nan
    allocate( this%fr_fcel       (0:mxpft) ); this%fr_fcel       (:) = nan
    allocate( this%fr_flig       (0:mxpft) ); this%fr_flig       (:) = nan
    allocate( this%leaf_long     (0:mxpft) ); this%leaf_long     (:) = nan
    allocate( this%evergreen     (0:mxpft) ); this%evergreen     (:) = nan
    allocate( this%stress_decid  (0:mxpft) ); this%stress_decid  (:) = nan
    allocate( this%season_decid  (0:mxpft) ); this%season_decid  (:) = nan
!KO
    allocate( this%season_decid_temperate (0:mxpft) ); this%season_decid_temperate (:) = nan  !#
!KO
    allocate( this%dwood         (0:mxpft) ); this%dwood          (:) = nan
    allocate( this%root_density  (0:mxpft) ); this%root_density   (:) = nan  !#
    allocate( this%root_radius   (0:mxpft) ); this%root_radius    (:) = nan  !#
    allocate( this%pconv         (0:mxpft) ); this%pconv          (:) = nan  !#
    allocate( this%pprod10       (0:mxpft) ); this%pprod10        (:) = nan  !#
    allocate( this%pprod100      (0:mxpft) ); this%pprod100       (:) = nan  !#
    allocate( this%pprodharv10   (0:mxpft) ); this%pprodharv10    (:) = nan  !#
    allocate( this%cc_leaf       (0:mxpft) ); this%cc_leaf        (:) = nan  
    allocate( this%cc_lstem      (0:mxpft) ); this%cc_lstem       (:) = nan 
    allocate( this%cc_dstem      (0:mxpft) ); this%cc_dstem       (:) = nan 
    allocate( this%cc_other      (0:mxpft) ); this%cc_other       (:) = nan
    allocate( this%fm_leaf       (0:mxpft) ); this%fm_leaf        (:) = nan
    allocate( this%fm_lstem      (0:mxpft) ); this%fm_lstem       (:) = nan
    allocate( this%fm_dstem      (0:mxpft) ); this%fm_dstem       (:) = nan
    allocate( this%fm_other      (0:mxpft) ); this%fm_other       (:) = nan
    allocate( this%fm_root       (0:mxpft) ); this%fm_root        (:) = nan
    allocate( this%fm_lroot      (0:mxpft) ); this%fm_lroot       (:) = nan
    allocate( this%fm_droot      (0:mxpft) ); this%fm_droot       (:) = nan
    allocate( this%fsr_pft       (0:mxpft) ); this%fsr_pft        (:) = nan
    allocate( this%fd_pft        (0:mxpft) ); this%fd_pft         (:) = nan
    allocate( this%rswf_max      (0:mxpft) ); this%rswf_max       (:) = nan  !#
    allocate( this%rswf_min      (0:mxpft) ); this%rswf_min       (:) = nan  !#
    allocate( this%manunitro     (0:mxpft) ); this%manunitro      (:) = nan  !#
    allocate( this%fleafcn       (0:mxpft) ); this%fleafcn        (:) = nan
    allocate( this%ffrootcn      (0:mxpft) ); this%ffrootcn       (:) = nan
    allocate( this%fstemcn       (0:mxpft) ); this%fstemcn        (:) = nan
    allocate( this%i_vcad        (0:mxpft) ); this%i_vcad         (:) = nan  !#
    allocate( this%s_vcad        (0:mxpft) ); this%s_vcad         (:) = nan  !#
    allocate( this%i_flnr        (0:mxpft) ); this%i_flnr         (:) = nan  !#
    allocate( this%s_flnr        (0:mxpft) ); this%s_flnr         (:) = nan  !#
    allocate( this%pftpar20      (0:mxpft) ); this%pftpar20       (:) = nan  !#
    allocate( this%pftpar28      (0:mxpft) ); this%pftpar28       (:) = nan  !#
    allocate( this%pftpar29      (0:mxpft) ); this%pftpar29       (:) = nan  !#
    allocate( this%pftpar30      (0:mxpft) ); this%pftpar30       (:) = nan  !#
    allocate( this%pftpar31      (0:mxpft) ); this%pftpar31       (:) = nan  !#
    allocate( this%a_fix         (0:mxpft) ); this%a_fix          (:) = nan  !#
    allocate( this%b_fix         (0:mxpft) ); this%b_fix          (:) = nan  !#
    allocate( this%c_fix         (0:mxpft) ); this%c_fix          (:) = nan  !#
    allocate( this%s_fix         (0:mxpft) ); this%s_fix          (:) = nan  !#
    allocate( this%akc_active    (0:mxpft) ); this%akc_active     (:) = nan  !#
    allocate( this%akn_active    (0:mxpft) ); this%akn_active     (:) = nan  !#
    allocate( this%ekc_active    (0:mxpft) ); this%ekc_active     (:) = nan  !#
    allocate( this%ekn_active    (0:mxpft) ); this%ekn_active     (:) = nan  !#
    allocate( this%kc_nonmyc     (0:mxpft) ); this%kc_nonmyc      (:) = nan  !#
    allocate( this%kn_nonmyc     (0:mxpft) ); this%kn_nonmyc      (:) = nan  !#
    allocate( this%kr_resorb     (0:mxpft) ); this%kr_resorb      (:) = nan  !#
    allocate( this%perecm        (0:mxpft) ); this%perecm         (:) = nan  !#
    allocate( this%root_dmx      (0:mxpft) ); this%root_dmx       (:) = nan  !#
    allocate( this%fun_cn_flex_a (0:mxpft) ); this%fun_cn_flex_a  (:) = nan  !#
    allocate( this%fun_cn_flex_b (0:mxpft) ); this%fun_cn_flex_b  (:) = nan  !#
    allocate( this%fun_cn_flex_c (0:mxpft) ); this%fun_cn_flex_c  (:) = nan  !#
    allocate( this%FUN_fracfixers(0:mxpft) ); this%FUN_fracfixers (:) = nan  !#
    allocate( this%dbh           (0:mxpft) ); this%dbh            (:) = nan  !#
    allocate( this%fbw           (0:mxpft) ); this%fbw            (:) = nan  !#
    allocate( this%nstem         (0:mxpft) ); this%nstem          (:) = nan  !#
    allocate( this%taper         (0:mxpft) ); this%taper          (:) = nan  !#
    allocate( this%rstem_per_dbh (0:mxpft) ); this%rstem_per_dbh  (:) = nan  !#
    allocate( this%wood_density  (0:mxpft) ); this%wood_density   (:) = nan  !#

    ! jkolassa, Dec 2021: read in parameters from CLM parameter file
    ! TO DO: pass parameter file through rc files rather than hardcoding name here

    paramfile = '/discover/nobackup/jkolassa/CLM/parameter_files/ctsm51_params.c210923.nc'
    ierr = NF90_OPEN(trim(paramfile),NF90_NOWRITE,ncid)
    if (ierr/=0) then
       _ASSERT(.FALSE.,'error opening netcdf file')
    end if

    ierr = NF90_INQ_VARID(ncid,'z0mr',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%z0mr(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'displar',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%displar(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'dleaf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%dleaf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'c3psn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%c3psn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rholvis',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rhol(:,1) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rholnir',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rhol(:,2) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rhosvis',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rhos(:,1) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rhosnir',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rhos(:,2) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'taulvis',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%taul(:,1) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'taulnir',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%taul(:,2) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'tausvis',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%taus(:,1) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'tausnir',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%taus(:,2) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'xl',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%xl(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'roota_par',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%roota_par(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rootb_par',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rootb_par(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'slatop',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%slatop(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'dsladlai',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%dsladlai(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'leafcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%leafcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'biofuel_harvfrac',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%biofuel_harvfrac(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'flnr',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%flnr(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'smpso',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%smpso(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'smpsc',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%smpsc(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fnitr',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fnitr(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'woody',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%woody(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'lflitcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%lflitcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'frootcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%frootcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'livewdcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%livewdcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'deadwdcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%deadwdcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'grperc',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%grperc(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'grpnow',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%grpnow(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'froot_leaf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%froot_leaf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'stem_leaf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%stem_leaf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'croot_stem',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%croot_stem(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'flivewd',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%flivewd(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fcur',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fcur(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fcurdv',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fcurdv(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'lf_flab',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%lf_flab(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'lf_fcel',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%lf_fcel(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'lf_flig',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%lf_flig(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fr_flab',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fr_flab(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fr_fcel',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fr_fcel(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fr_flig',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fr_flig(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'leaf_long',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%leaf_long(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'evergreen',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%evergreen(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'stress_decid',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%stress_decid(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'season_decid',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%season_decid(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'season_decid_temperate',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%z0mr(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pftpar20',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pftpar20(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pftpar28',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pftpar28(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pftpar29',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pftpar29(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pftpar30',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pftpar30(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pftpar31',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pftpar31(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'a_fix',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%a_fix(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'b_fix',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%b_fix(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'c_fix',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%c_fix(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'s_fix',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%s_fix(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'akc_active',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%akc_active(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'akn_active',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%akn_active(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'ekc_active',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%ekc_active(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'ekn_active',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%ekn_active(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'kc_nonmyc',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%kc_nonmyc(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'kn_nonmyc',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%kn_nonmyc(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'kr_resorb',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%kr_resorb(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'perecm',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%perecm(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fun_cn_flex_a',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fun_cn_flex_a(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fun_cn_flex_b',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fun_cn_flex_b(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fun_cn_flex_c',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fun_cn_flex_c(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'FUN_fracfixers',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%FUN_fracfixers(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'manunitro',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%manunitro(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fleafcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fleafcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'ffrootcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%ffrootcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fstemcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fstemcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rootprof_beta',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_2)
    this%rootprof_beta(:,:) = read_tmp_2(0:mxpft,:)

    ierr = NF90_INQ_VARID(ncid,'pconv',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pconv(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pprod10',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pprod10(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pprodharv10',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pprodharv10(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'pprod100',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%pprod100(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'graincn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%graincn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'mxtmp',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%mxtmp(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'baset',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%baset(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'declfact',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%declfact(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'bfact',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%bfact(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'aleaff',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%aleaff(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'arootf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%arootf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'astemf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%astemf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'arooti',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%arooti(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fleafi',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fleafi(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'allconsl',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%allconsl(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'allconss',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%allconss(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'crop',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%crop(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'mergetoclmpft',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%mergetoclmpft(:) = read_tmp_3(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'irrigated',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%irrigated(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'ztopmx',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%ztopmx(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'laimx',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%laimx(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'gddmin',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%gddmin(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'hybgdd',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%hybgdd(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'lfemerg',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%lfemerg(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'grnfill',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%grnfill(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'mbbopt',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%mbbopt(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'medlynslope',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%medlynslope(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'medlynintercept',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%medlynintercept(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'mxmat',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%mxmat(:) = read_tmp_3(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'cc_leaf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%cc_leaf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'cc_lstem',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%cc_lstem(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'cc_dstem',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%cc_dstem(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'cc_other',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%cc_other(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fstemcn',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fstemcn(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_leaf',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_leaf(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_lstem',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_lstem(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_dstem',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_dstem(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_other',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_other(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_root',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_root(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_lroot',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_lroot(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fm_droot',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fm_droot(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fsr_pft',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fsr_pft(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'fd_pft',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%fd_pft(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rswf_min',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rswf_min(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'rswf_max',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%rswf_max(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'min_planting_temp',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%min_planting_temp(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'min_NH_planting_date',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%min_NH_planting_date(:) = read_tmp_3(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'min_SH_planting_date',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%min_SH_planting_date(:) = read_tmp_3(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'max_NH_planting_date',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%max_NH_planting_date(:) = read_tmp_3(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'max_SH_planting_date',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_3)
    this%max_SH_planting_date(:) = read_tmp_3(0:mxpft)

    do m = 0,mxpft
       this%dwood(m) = dwood
       this%root_radius(m)  = root_radius
       this%root_density(m) = root_density
    end do

    if (use_flexibleCN) then
       ierr = NF90_INQ_VARID(ncid,'i_vcad',clm_varid)
       ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
       this%i_vcad(:) = read_tmp_1(0:mxpft)

       ierr = NF90_INQ_VARID(ncid,'s_vcad',clm_varid)
       ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
       this%s_vcad(:) = read_tmp_1(0:mxpft)

       ierr = NF90_INQ_VARID(ncid,'i_flnr',clm_varid)
       ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
       this%i_flnr(:) = read_tmp_1(0:mxpft)

       ierr = NF90_INQ_VARID(ncid,'s_flnr',clm_varid)
       ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
       this%s_flnr(:) = read_tmp_1(0:mxpft)

    end if

    if ( use_crop .and. use_dynroot )then
       ierr = NF90_INQ_VARID(ncid,'root_dmx',clm_varid)
       ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
       this%root_dmx(:) = read_tmp_1(0:mxpft)
    end if

    ierr = NF90_INQ_VARID(ncid,'nstem',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%nstem(:) = read_tmp_1(0:mxpft)

    ierr = NF90_INQ_VARID(ncid,'taper',clm_varid)
    ierr = NF90_GET_VAR(ncid, clm_varid, read_tmp_1)
    this%taper(:) = read_tmp_1(0:mxpft)

    ierr = NF90_CLOSE(ncid)
   
    ! jkolassa, Dec 2021: not using biomass heat storage module, so set the following 4 parameters to 0
    this%dbh = 0.0_r8
    this%fbw = 0.0_r8
    this%rstem_per_dbh = 0.0_r8
    this%wood_density = 0.0_r8

    ! Set vegetation family identifier (tree/shrub/grass)
    do m = 0,mxpft 
       if (m == ndllf_evr_tmp_tree .or. m == ndllf_evr_brl_tree &
            .or. m == ndllf_dcd_brl_tree .or. m == nbrdlf_evr_trp_tree &
            .or. m == nbrdlf_evr_tmp_tree .or. m == nbrdlf_dcd_trp_tree &
            .or. m == nbrdlf_dcd_tmp_tree .or. m == nbrdlf_dcd_brl_tree) then 
          this%is_tree(m) = .true.
       else 
          this%is_tree(m) = .false.
       endif
       if(m == nbrdlf_evr_shrub .or. m == nbrdlf_dcd_tmp_shrub .or. m == nbrdlf_dcd_brl_shrub) then 
          this%is_shrub(m) = .true.
       else 
          this%is_shrub(m) = .false.
       endif
       if(m == nc3_arctic_grass .or. m == nc3_nonarctic_grass .or. m == nc4_grass) then 
          this%is_grass(m) = .true.
       else 
          this%is_grass(m) = .false.
       endif

    end do

    if (use_cndv) then
       this%fcur(:) = this%fcurdv(:)
    end if

    ! jk, Dec 2021:  we are not using the crop or irrigation modules at this point, so set the flags to 0 everywhere
    
    this%irrigated(:) = 0.0_r8
    this%crop(:)      = 0.0_r8

    ! jk Dec 2021: all PFTs are known to model since we are not using the crop model, so set flag to true everywhere
    this%is_pft_known_to_model(:) = .true.

 end subroutine init_pftcon_type

end module CNCLM_pftconMod
