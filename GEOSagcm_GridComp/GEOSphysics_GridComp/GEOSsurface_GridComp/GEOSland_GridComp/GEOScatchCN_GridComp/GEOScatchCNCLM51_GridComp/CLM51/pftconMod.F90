#include "MAPL_Generic.h"

module pftconMod

  use shr_kind_mod     , only: r8 => shr_kind_r8
  use nanMod           , only : nan, bigint
  use clm_varpar       , only : mxpft, numrad,nvariants, ivis, inir
  use clm_varctl       , only : use_flexibleCN, use_cndv
  use netcdf 
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use MAPL             , only : NetCDF4_FileFormatter, pFIO_READ
  use MAPL_ExceptionHandling
  use ncdio_pio        , only : ncd_io
  use ESMF

  ! !PUBLIC TYPES:
  implicit none

  INCLUDE 'netcdf.inc'
  save
!
! !PUBLIC MEMBER FUNCTIONS:
 

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
  integer, public :: npcropmin             = 16   ! value for first crop functional type (not including the more generic C3 crop PFT)

  ! variables that do not apply here, but are needed; set to mxpft + 1 in initialization routine

  integer, public :: ntmp_corn                    ! value for temperate corn, rain fed (rf)
  integer, public :: nirrig_tmp_corn        ! value for temperate corn, irrigated (ir)
  integer, public :: nswheat                ! value for spring temperate cereal (rf)
  integer, public :: nirrig_swheat          ! value for spring temperate cereal (ir)
  integer, public :: nwwheat                ! value for winter temperate cereal (rf)
  integer, public :: nirrig_wwheat          ! value for winter temperate cereal (ir)
  integer, public :: ntmp_soybean           ! value for temperate soybean (rf)
  integer, public :: nirrig_tmp_soybean     ! value for temperate soybean (ir)
  integer, public :: nbarley                ! value for spring barley (rf)
  integer, public :: nirrig_barley          ! value for spring barley (ir)
  integer, public :: nwbarley               ! value for winter barley (rf)
  integer, public :: nirrig_wbarley         ! value for winter barley (ir)
  integer, public :: nrye                   ! value for spring rye (rf)
  integer, public :: nirrig_rye             ! value for spring rye (ir)
  integer, public :: nwrye                  ! value for winter rye (rf)
  integer, public :: nirrig_wrye            ! value for winter rye (ir)
  integer, public :: ncassava               ! ...and so on
  integer, public :: nirrig_cassava
  integer, public :: ncitrus
  integer, public :: nirrig_citrus
  integer, public :: ncocoa
  integer, public :: nirrig_cocoa
  integer, public :: ncoffee
  integer, public :: nirrig_coffee
  integer, public :: ncotton
  integer, public :: nirrig_cotton
  integer, public :: ndatepalm
  integer, public :: nirrig_datepalm
  integer, public :: nfoddergrass
  integer, public :: nirrig_foddergrass
  integer, public :: ngrapes
  integer, public :: nirrig_grapes
  integer, public :: ngroundnuts
  integer, public :: nirrig_groundnuts
  integer, public :: nmillet
  integer, public :: nirrig_millet
  integer, public :: noilpalm
  integer, public :: nirrig_oilpalm
  integer, public :: npotatoes
  integer, public :: nirrig_potatoes
  integer, public :: npulses
  integer, public :: nirrig_pulses
  integer, public :: nrapeseed
  integer, public :: nirrig_rapeseed
  integer, public :: nrice
  integer, public :: nirrig_rice
  integer, public :: nsorghum
  integer, public :: nirrig_sorghum
  integer, public :: nsugarbeet
  integer, public :: nirrig_sugarbeet
  integer, public :: nsugarcane
  integer, public :: nirrig_sugarcane
  integer, public :: nsunflower
  integer, public :: nirrig_sunflower
  integer, public :: nmiscanthus
  integer, public :: nirrig_miscanthus
  integer, public :: nswitchgrass
  integer, public :: nirrig_switchgrass
  integer, public :: ntrp_corn              !value for tropical corn (rf)
  integer, public :: nirrig_trp_corn        !value for tropical corn (ir)
  integer, public :: ntrp_soybean           !value for tropical soybean (rf)
  integer, public :: nirrig_trp_soybean     !value for tropical soybean (ir)
  integer, public :: npcropmax              ! value for last prognostic crop in list
  integer, public :: nc3irrig               ! value for irrigated generic crop (ir)

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

  contains

     procedure, public ::  init_pftcon_type

  end type pftcon_type

type(pftcon_type), public, target, save :: pftcon

  integer, public, parameter :: pftname_len = 40         ! max length of pftname       
  character(len=pftname_len), public :: pftname(0:mxpft) ! PFT description

  real(r8), public, parameter :: reinickerp = 1.6_r8     ! parameter in allometric equation
  real(r8), public, parameter :: dwood  = 2.5e5_r8       ! cn wood density (gC/m3); lpj:2.0e5
  real(r8), public, parameter :: allom1 = 100.0_r8       ! parameters in
  real(r8), public, parameter :: allom2 =  40.0_r8       ! ...allometric
  real(r8), public, parameter :: allom3 =   0.5_r8       ! ...equations
  real(r8), public, parameter :: allom1s = 250.0_r8      ! modified for shrubs by
  real(r8), public, parameter :: allom2s =   8.0_r8      ! X.D.Z
! root radius, density from Bonan, GMD, 2014
  real(r8), public, parameter :: root_density = 0.31e06_r8 !(g biomass / m3 root)
  real(r8), public, parameter :: root_radius = 0.29e-03_r8 !(m)

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

!--------------------------------
  subroutine init_pftcon_type(this,paramfile)

  ! !DESCRIPTION:
! Initialize CTSM PFT constants                                
!                                                                              
  use abortutils  , only : endrun

! !ARGUMENTS:                                                                                                           
    implicit none
    !INPUT/OUTPUT
    class(pftcon_type) :: this
    character(len=ESMF_MAXSTR), intent(in) :: paramfile


    integer            :: ierr, clm_varid,  status, m, n
    logical            :: readv ! has variable been read in or not
    type(Netcdf4_fileformatter) :: ncid

    real(r8), allocatable, dimension(:)   :: read_tmp_1
    real(r8), allocatable, dimension(:,:) :: read_tmp_2
    integer , allocatable, dimension(:)   :: read_tmp_3

    character(len=512) :: msg

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
    allocate( this%mergetoclmpft (0:mxpft) ); this%mergetoclmpft (:) = bigint !#
    allocate( this%is_pft_known_to_model  (0:mxpft) ); this%is_pft_known_to_model(:) = .false. !#
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
    allocate( this%mxmat         (0:mxpft) ); this%mxmat     (:) = bigint
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

    call ncid%open(trim(paramfile),pFIO_READ, RC=status)

    call ncd_io('pftname',pftname, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('z0mr', this%z0mr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('displar', this%displar, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('dleaf', this%dleaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('c3psn', this%c3psn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rholvis', this%rhol(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rholnir', this%rhol(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rhosvis', this%rhos(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rhosnir', this% rhos(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('taulvis', this%taul(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('taulnir', this%taul(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('tausvis', this%taus(:,ivis), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('tausnir', this%taus(:,inir), 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('xl', this%xl, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('roota_par', this%roota_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rootb_par', this%rootb_par, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('slatop', this%slatop, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('dsladlai', this%dsladlai, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('leafcn', this%leafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('biofuel_harvfrac', this%biofuel_harvfrac, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('flnr', this%flnr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('smpso', this%smpso, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('smpsc', this%smpsc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fnitr', this%fnitr, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('woody', this%woody, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lflitcn', this%lflitcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('frootcn', this%frootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('livewdcn', this%livewdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('deadwdcn', this%deadwdcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('grperc', this%grperc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('grpnow', this%grpnow, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('froot_leaf', this%froot_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('stem_leaf', this%stem_leaf, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('croot_stem', this%croot_stem, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('flivewd', this%flivewd, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fcur', this%fcur, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fcurdv', this%fcurdv, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lf_flab', this%lf_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lf_fcel', this%lf_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lf_flig', this%lf_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fr_flab', this%fr_flab, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fr_fcel', this%fr_fcel, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fr_flig', this%fr_flig, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('leaf_long', this%leaf_long, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('evergreen', this%evergreen, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('stress_decid', this%stress_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('season_decid', this%season_decid, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

!KO
    call ncd_io('season_decid_temperate', this%season_decid_temperate, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
!KO

    call ncd_io('pftpar20', this%pftpar20, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar28', this%pftpar28, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar29', this%pftpar29, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar30', this%pftpar30, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pftpar31', this%pftpar31, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('a_fix', this%a_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('b_fix', this%b_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('c_fix', this%c_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('s_fix', this%s_fix, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('akc_active', this%akc_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('akn_active', this%akn_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ekc_active', this%ekc_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ekn_active', this%ekn_active, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('kc_nonmyc', this%kc_nonmyc, 'read', ncid, readvar=readv,   posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('kn_nonmyc', this%kn_nonmyc, 'read', ncid, readvar=readv,   posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('kr_resorb', this%kr_resorb, 'read', ncid, readvar=readv,   posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('perecm', this%perecm, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fun_cn_flex_a', this%fun_cn_flex_a, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fun_cn_flex_b', this%fun_cn_flex_b, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fun_cn_flex_c', this%fun_cn_flex_c, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('FUN_fracfixers', this%FUN_fracfixers, 'read', ncid, readvar=readv,         posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('manunitro', this%manunitro, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fleafcn', this%fleafcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ffrootcn', this%ffrootcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fstemcn', this%fstemcn, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rootprof_beta', this%rootprof_beta, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pconv', this%pconv, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pprod10', this%pprod10, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pprodharv10', this%pprodharv10, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('pprod100', this%pprod100, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('graincn', this%graincn, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mxtmp', this%mxtmp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('baset', this%baset, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('declfact', this%declfact, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('bfact', this%bfact, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('aleaff', this%aleaff, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('arootf', this%arootf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('astemf', this%astemf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('arooti', this%arooti, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fleafi', this%fleafi, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('allconsl', this%allconsl, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('allconss', this%allconss, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('crop', this%crop, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mergetoclmpft', this%mergetoclmpft, 'read', ncid, readvar=readv)
    if ( .not. readv ) then
       call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    end if

    call ncd_io('irrigated', this%irrigated, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('ztopmx', this%ztopmx, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('laimx', this%laimx, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('gddmin', this%gddmin, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('hybgdd', this%hybgdd, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('lfemerg', this%lfemerg, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('grnfill', this%grnfill, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mbbopt', this%mbbopt, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('medlynslope', this%medlynslope, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('medlynintercept', this%medlynintercept, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('mxmat', this%mxmat, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_leaf', this% cc_leaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_lstem', this%cc_lstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_dstem', this%cc_dstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('cc_other', this%cc_other, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_leaf', this% fm_leaf, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_lstem', this%fm_lstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_dstem', this%fm_dstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_other', this%fm_other, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_root', this% fm_root, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_lroot', this%fm_lroot, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fm_droot', this%fm_droot, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fsr_pft', this% fsr_pft, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('fd_pft', this%  fd_pft, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rswf_min', this% rswf_min, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('rswf_max', this% rswf_max, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    call ncd_io('planting_temp', this%planttemp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('min_planting_temp', this%minplanttemp, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('min_NH_planting_date', this%mnNHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('min_SH_planting_date', this%mnSHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('max_NH_planting_date', this%mxNHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

    call ncd_io('max_SH_planting_date', this%mxSHplantdate, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    
   npcropmax          = mxpft            ! last prognostic crop in list

   ! jkolassa Jan 2023: we do not model these crops, but the below variables are needed
   ! for some checks; setting them to mxpft + 1
   ntmp_corn          = mxpft + 1        ! value for temperate corn, rain fed (rf)
   nirrig_tmp_corn    = mxpft + 1        ! value for temperate corn, irrigated (ir)
   nswheat            = mxpft + 1        ! value for spring temperate cereal (rf)
   nirrig_swheat      = mxpft + 1        ! value for spring temperate cereal (ir)
   nwwheat            = mxpft + 1        ! value for winter temperate cereal (rf)
   nirrig_wwheat      = mxpft + 1        ! value for winter temperate cereal (ir)
   ntmp_soybean       = mxpft + 1        ! value for temperate soybean (rf)
   nirrig_tmp_soybean = mxpft + 1        ! value for temperate soybean (ir)
   nbarley            = mxpft + 1        ! value for spring barley (rf)
   nirrig_barley      = mxpft + 1        ! value for spring barley (ir)
   nwbarley           = mxpft + 1        ! value for winter barley (rf)
   nirrig_wbarley     = mxpft + 1        ! value for winter barley (ir)
   nrye               = mxpft + 1        ! value for spring rye (rf)
   nirrig_rye         = mxpft + 1        ! value for spring rye (ir)
   nwrye              = mxpft + 1        ! value for winter rye (rf)
   nirrig_wrye        = mxpft + 1        ! value for winter rye (ir)
   ncassava           = mxpft + 1        ! ...and so on
   nirrig_cassava     = mxpft + 1
   ncitrus            = mxpft + 1
   nirrig_citrus      = mxpft + 1
   ncocoa             = mxpft + 1
   nirrig_cocoa       = mxpft + 1
   ncoffee            = mxpft + 1
   nirrig_coffee      = mxpft + 1
   ncotton            = mxpft + 1
   nirrig_cotton      = mxpft + 1
   ndatepalm          = mxpft + 1
   nirrig_datepalm    = mxpft + 1
   nfoddergrass       = mxpft + 1
   nirrig_foddergrass = mxpft + 1
   ngrapes            = mxpft + 1
   nirrig_grapes      = mxpft + 1
   ngroundnuts        = mxpft + 1
   nirrig_groundnuts  = mxpft + 1
   nmillet            = mxpft + 1
   nirrig_millet      = mxpft + 1
   noilpalm           = mxpft + 1
   nirrig_oilpalm     = mxpft + 1
   npotatoes          = mxpft + 1
   nirrig_potatoes    = mxpft + 1
   npulses            = mxpft + 1
   nirrig_pulses      = mxpft + 1
   nrapeseed          = mxpft + 1
   nirrig_rapeseed    = mxpft + 1
   nrice              = mxpft + 1
   nirrig_rice        = mxpft + 1
   nsorghum           = mxpft + 1
   nirrig_sorghum     = mxpft + 1
   nsugarbeet         = mxpft + 1
   nirrig_sugarbeet   = mxpft + 1
   nsugarcane         = mxpft + 1
   nirrig_sugarcane   = mxpft + 1
   nsunflower         = mxpft + 1
   nirrig_sunflower   = mxpft + 1
   nmiscanthus        = mxpft + 1
   nirrig_miscanthus  = mxpft + 1
   nswitchgrass       = mxpft + 1
   nirrig_switchgrass = mxpft + 1
   ntrp_corn          = mxpft + 1    !value for tropical corn (rf)
   nirrig_trp_corn    = mxpft + 1    !value for tropical corn (ir)
   ntrp_soybean       = mxpft + 1    !value for tropical soybean (rf)
   nirrig_trp_soybean = mxpft + 1    !value for tropical soybean (ir)
   npcropmax          = mxpft + 1    ! value for last prognostic crop in list
   nc3irrig           = mxpft + 1    ! value for irrigated generic crop (ir)

    do m = 0,mxpft
       this%dwood(m) = dwood
       this%root_radius(m)  = root_radius
       this%root_density(m) = root_density
    end do

    !
    ! clm 5 nitrogen variables
    !
    if (use_flexibleCN) then
       call ncd_io('i_vcad', this%i_vcad, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

       call ncd_io('s_vcad', this%s_vcad, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

       call ncd_io('i_flnr', this%i_flnr, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

       call ncd_io('s_flnr', this%s_flnr, 'read', ncid, readvar=readv)
       if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    end if

    call ncd_io('nstem',this%nstem, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))
    call ncd_io('taper',this%taper, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=' ERROR: error in reading in pft data'//errMsg(sourcefile, __LINE__))

   
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

end module pftconMod
