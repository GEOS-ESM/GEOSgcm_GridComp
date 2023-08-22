module clm_varpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
!
 use clm_varpar_shared, only : VAR_COL =>VAR_COL_45, VAR_PFT => VAR_PFT_45, &
                               numpft => numpft_CN, NUM_ZON => NUM_ZON_CN, &
                               NUM_VEG => NUM_VEG_CN
! !PUBLIC TYPES:
  implicit none
  save
!
! Define number of levels

  integer, parameter :: nlevsoi     =   1     ! number of hydrologically active soil layers
  integer, parameter :: nlevgrnd    =   1     ! number of ground layers (includes lower layers that are hydrologically inactive)
  integer, parameter :: nlevsno     =   0     ! maximum number of snow layers
  integer, parameter :: nlevcan     =   1     ! number of canopy layers
  integer            :: nlevurb               ! number of urban layers
  logical, public :: more_vertlayers = .false. ! true => run with more vertical soil layers. ".false." is the default setting in CLM4.5

! Define indices used in surface file read
! maxpatch_pft     = max number of plant functional types in naturally vegetated landunit

  integer            :: maxpatch_pft

! clm_varpar_init seems to do something similar; less prone to error to move
! these three lines there? (slevis)
  integer, parameter :: max_pft_per_col   = numpft+1
  
  integer            :: nlevdecomp                    ! number of biogeochemically active soil layers
  integer            :: nlevdecomp_full               ! number of biogeochemical layers (includes lower layers that are biogeochemically inactive)

! For CH4 code
  integer, parameter :: ngases = 3 ! CH4, O2, & CO2

! !PUBLIC MEMBER FUNCTIONS:
  public clm_varpar_init          ! set parameters

! !REVISION HISTORY:
! Created by Mariana Vertenstein

  ! CatchCN parameters
  ! ------------------

! CN types:
! https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/lnd/clm2/pftdata/pft-physiology.c130503.nc
!
! PFT Description 
!  0   0  Bare                                              
!  1   1  Needleleaf evergreen temperate tree                        
!  2   2  Needleleaf evergreen boreal tree                          
!  3   3  Needleleaf deciduous boreal tree                           
!  4   4  Broadleaf evergreen tropical tree                          
!  5   5  Broadleaf evergreen temperate tree                         
!  6   6  Broadleaf deciduous tropical tree                         
!  7   7  Broadleaf deciduous temperate tree                         
!  8   8  Broadleaf deciduous boreal tree                            
!  9   9  Broadleaf evergreen temperate shrub                        
! 10  10  Broadleaf deciduous temperate shrub [moisture + deciduous]                      
! 11   -  Broadleaf deciduous temperate shrub [moisture stress only] 
! 12  11  Broadleaf deciduous boreal shrub                           
! 13  12  Arctic c3 grass                                            
! 14  13  Cool c3 grass [moisture + deciduous]                                            
! 15   -  Cool c3 grass [moisture stress only]                       
! 16  14  Warm c4 grass [moisture + deciduous]                                             
! 17   -  Warm c4 grass [moisture stress only]                       
! 18  15  C3 crop [moisture + deciduous]                                                   
! 19  16  C3 crop [moisture stress only]                       
  
! Catchment types and potential PFT mapping:
! 
! 1:  BROADLEAF EVERGREEN TREES  => 4,5
! 2:  BROADLEAF DECIDUOUS TREES  => 6,7,8
! 3:  NEEDLELEAF TREES           => 1,2,3
! 4:  GROUND COVER               => 0,13-19
! 5:  BROADLEAF SHRUBS           => 9,10,11
! 6:  DWARF TREES (TUNDRA)       =>  12
! 7:  BARE SOIL                  =>  0
! 8:  DESERT                     =>  0
! 9:  ICE                        => n/a  

  real, parameter, PUBLIC, dimension(NUM_ZON) :: CN_zone_weight = (/0.10,0.45,0.45/) ! gkw: tunable; must sum to 1
  integer, parameter, PUBLIC :: map_cat(0:numpft) = (/4,3,3,3,1,1,2,2,2,5,5,5,6,4,4,4,4,4,4,4/) 
  
! -------------------------------------------------------
! Module Varaibles (initialized in clm_varpar_init)
! -------------------------------------------------------

#ifndef CENTURY_DECOMP
  ! parameters for decomposition cascade
  integer, parameter :: ndecomp_pools = 8
  integer, parameter :: ndecomp_cascade_transitions = 9
  integer, parameter :: i_met_lit = 1
  integer, parameter :: i_cel_lit = 2
  integer, parameter :: i_lig_lit = 3
  integer, parameter :: i_cwd = 4
  integer, parameter :: nsompools = 4
#else
  ! parameters for decomposition cascade
  integer, parameter :: ndecomp_pools = 7
  integer, parameter :: ndecomp_cascade_transitions = 10
  integer, parameter :: i_met_lit = 1
  integer, parameter :: i_cel_lit = 2
  integer, parameter :: i_lig_lit = 3
  integer, parameter :: i_cwd = 4
  integer, parameter :: nsompools = 3
#endif  

!EOP
!-----------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_varpar_init
!
! !INTERFACE:
  subroutine clm_varpar_init()
!
! !DESCRIPTION:
! This subroutine initializes parameters in clm_varpar
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  maxpatch_pft   = numpft + 1 ! MAXPATCH_PFT ! gkw: this was set via compiler directive
  nlevurb        = 5          ! use the value from CLM4.5 for now, change later if needed, fzeng, 28 Mar 2017
  
  ! here is a switch to set the number of soil levels for the biogeochemistry calculations.
  ! currently it works on either a single level or on nlevsoi and nlevgrnd levels
#ifdef VERTSOILC
  nlevdecomp      = nlevsoi
  nlevdecomp_full = nlevgrnd
#else
  nlevdecomp      = 1
  nlevdecomp_full = 1
#endif

  end subroutine clm_varpar_init

!------------------------------------------------------------------------------
end module clm_varpar
