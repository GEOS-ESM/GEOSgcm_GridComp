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
  use clm_varpar_shared, only :    &
       VAR_COL => VAR_COL_51,      &
       VAR_PFT => VAR_PFT_51,      &
       numpft  => NUM_PFT_CN_51,   &
       NUM_ZON => NUM_ZON_CN,      &
       NUM_VEG => NUM_VEG_CN_51
  
! !PUBLIC TYPES:
  implicit none
  save
  
  ! Define number of levels
  
  integer,         parameter :: nlevsoi           =  1      ! number of hydrologically active soil layers
  integer,         parameter :: nlevgrnd          =  1      ! number of ground layers (includes lower layers that are hydrologically inactive)
  integer,         parameter :: nlevsno           =  0      ! maximum number of snow layers
  integer, public            :: nlevurb           =  0      ! number of urban layers; jkolassa Aug 2022: changed b/c having more urban than ground layers caused an issue w/ initialization of soil layers in column type
  integer, public            :: nlevmaxurbgrnd              ! maximum of the number of ground and urban layers
  integer, public, parameter :: nlayer            =  3      ! number of VIC soil layer --Added by AWang

  integer, public            :: nlevlak                     ! number of lake layers
  integer, public            :: nlevdecomp                  ! number of biogeochemically active soil layers
  integer, public            :: nlevdecomp_full             ! number of biogeochemical layers 
                                                            ! (includes lower layers that are biogeochemically inactive)

  integer, public            :: ndecomp_pools
  integer, public            :: ndecomp_cascade_transitions
  integer, public            :: ndecomp_cascade_outtransitions

  ! for soil matrix 
  integer, public            :: ndecomp_pools_vr            ! total number of pools ndecomp_pools*vertical levels

  integer,         parameter :: mxpft             = 15
  integer, public            :: maxveg                          ! # of pfts + cfts
  integer, public            :: maxsoil_patches   = numpft + 1  ! # of pfts + cfts + bare ground; replaces maxpatch_pft, which is obsolete

  integer, public            :: natpft_lb         =  0      ! In PATCH arrays, lower bound of Patches on the natural veg landunit (i.e., bare ground index)
                                                       
  integer, public, parameter :: nvariants         =  2      ! number of variants of PFT constants
                                                       
  integer, public, parameter :: numrad            =  2      ! number of solar radiation bands: vis, nir
  integer, public, parameter :: ivis              =  1      ! index for visible band
  integer, public, parameter :: inir              =  2      ! index for near-infrared band
  integer, public, parameter :: nlevcan           =  1      ! number of leaf layers in canopy layer
  integer, public, parameter :: nvegwcs           =  4      ! number of vegetation water conductance segments

  real,    public, parameter, dimension(NUM_ZON) :: CN_zone_weight    = (/0.10,0.45,0.45/) ! gkw: tunable; must sum to 1
  
  integer, public, parameter                     :: map_cat(0:numpft) = (/4,3,3,3,1,1,2,2,2,5,5,6,4,4,4,4/)

  ! constants for decomposition cascade

  integer, public, parameter :: i_met_lit         = 1
  integer, public, parameter :: i_cel_lit         = i_met_lit + 1
  integer, public, parameter :: i_lig_lit         = i_cel_lit + 1
  integer, public            :: i_cwd

   !Matrix index (when use_matrixcn)
  integer, public, parameter :: ileaf             =  1  ! leaf pool index
  integer, public, parameter :: ileaf_st          =  2  ! leaf storage pool index
  integer, public, parameter :: ileaf_xf          =  3  ! leaf transfer pool index
  integer, public, parameter :: ifroot            =  4  ! fine root pool index
  integer, public, parameter :: ifroot_st         =  5  ! fine root storage pool index
  integer, public, parameter :: ifroot_xf         =  6  ! fine root transfer pool index
  integer, public, parameter :: ilivestem         =  7  ! live stem pool index
  integer, public, parameter :: ilivestem_st      =  8  ! live stem storage pool index
  integer, public, parameter :: ilivestem_xf      =  9  ! live stem transfer pool index
  integer, public, parameter :: ideadstem         = 10  ! dead stem pool index
  integer, public, parameter :: ideadstem_st      = 11  ! dead stem storage pool index
  integer, public, parameter :: ideadstem_xf      = 12  ! dead stem transfer pool index
  integer, public, parameter :: ilivecroot        = 13  ! live coarse root pool index
  integer, public, parameter :: ilivecroot_st     = 14  ! live coarse root storage pool index
  integer, public, parameter :: ilivecroot_xf     = 15  ! live coarse root transfer pool index
  integer, public, parameter :: ideadcroot        = 16  ! dead coarse root pool index
  integer, public, parameter :: ideadcroot_st     = 17  ! dead coarse root storage pool index
  integer, public, parameter :: ideadcroot_xf     = 18  ! dead coarse root transfer pool index
  integer, public, parameter :: igrain            = 19  ! grain pool index
  integer, public, parameter :: igrain_st         = 20  ! grain storage pool index
  integer, public, parameter :: igrain_xf         = 21  ! grain transfer pool

  integer, public            :: ncphtrans               ! maximum number of vegetation C transfers through phenology
  integer, public            :: ncphouttrans            ! maximum number of vegetation C transfers out of vegetation through phenology
  integer, public            :: ncgmtrans               ! maximum number of vegetation C transfers through gap mortality
  integer, public            :: ncgmouttrans            ! maximum number of vegetation C transfers out of vegetation through gap mortality
  integer, public            :: ncfitrans               ! maximum number of vegetation C transfers through fire
  integer, public            :: ncfiouttrans            ! maximum number of vegetation C transfers out of vegetation trhough fire
  integer, public            :: nnphtrans               ! maximum number of vegetation N transfers through phenology
  integer, public            :: nnphouttrans            ! maximum number of vegetation N transfers out of vegetation through phenology
  integer, public            :: nngmtrans               ! maximum number of vegetation N transfers through gap mortality
  integer, public            :: nngmouttrans            ! maximum number of vegetation N transfers out of vegetation through gap mortality
  integer, public            :: nnfitrans               ! maximum number of vegetation N transfers through fire
  integer, public            :: nnfiouttrans            ! maximum number of vegetation N transfers out of vegetation trhough fire
                                                        
  integer, public            :: iretransn               ! retranslocation pool index
  integer, public            :: ioutc                   ! external C pool index
  integer, public            :: ioutn                   ! external N pool index


  integer, public, parameter :: nvegpool_natveg   = 18  ! number of vegetation matrix pool without crop
  integer, public, parameter :: nvegpool_crop     =  3  ! number of vegetation matrix pool with crop
  integer, public, parameter :: nveg_retransn     =  1  ! number of vegetation retranslocation pool
  integer, public            :: nvegcpool               ! number of vegetation C pools
  integer, public            :: nvegnpool               ! number of vegetation N pools


  ! For CH4 code 
  integer,         parameter :: ngases            = 3   ! CH4, O2, & CO2

  integer, public            :: max_patch_per_col

  integer, public            :: maxpatch_glcmec   = 0   ! max number of elevation classes (set to 0 here, not specified in CLM clm_varpar.F90)

contains

  !------------------------------------
  subroutine clm_varpar_init()
    !
    ! !DESCRIPTION:
    ! This subroutine initializes parameters in clm_varpar
    !
    use clm_varctl, only : use_vertsoilc, use_extralakelayers, use_fates, &
                         use_century_decomp, use_crop
    !
    ! !ARGUMENTS:
    implicit none

    !----------------------------

    nlevmaxurbgrnd    = max0(nlevurb,nlevgrnd)
    nlevmaxurbgrnd    = nlevgrnd                ! jkolassa: set this here, since we are not modelling urban tiles for now
    max_patch_per_col = maxsoil_patches         ! since we don't have CFTs or urban patches
    maxveg            = maxsoil_patches - 1     ! # of patches without bare ground

    ! here is a switch to set the number of soil levels for the biogeochemistry calculations.
    ! currently it works on either a single level or on nlevsoi and nlevgrnd levels
    if (use_vertsoilc) then
       nlevdecomp      = nlevsoi
       nlevdecomp_full = nlevgrnd
       ! nlevdecomp_full = nlevdecomp + 1 !jkolassa Nov 2024: nlevdecomp_full needs to be larger than nlevdecomp 
       ! when use_vertsoilc is true
    else
       nlevdecomp      = 1
       nlevdecomp_full = 1
    end if
    
    if (.not. use_extralakelayers) then
       nlevlak     =  10     ! number of lake layers
    else
       nlevlak     =  25     ! number of lake layers (Yields better results for site simulations)
    end if
    
    if ( use_fates ) then
       i_cwd = 0
       if (use_century_decomp) then
          ndecomp_pools               = 6
          ndecomp_cascade_transitions = 8
       else
          ndecomp_pools               = 7
          ndecomp_cascade_transitions = 7
       end if
    else
       i_cwd = 4
       if (use_century_decomp) then
          ndecomp_pools                  =  7
          ndecomp_cascade_transitions    = 10
          ndecomp_cascade_outtransitions =  0
       else
          ndecomp_pools                  =  8
          ndecomp_cascade_transitions    =  9
          ndecomp_cascade_outtransitions =  1
       end if
    endif
    ndecomp_pools_vr = ndecomp_pools * nlevdecomp

    if (use_crop)then
       nvegcpool = nvegpool_natveg + nvegpool_crop
       ncphtrans    = 18
       nnphtrans    = 37
       ncphouttrans =  4
       nnphouttrans =  5
    else
       nvegcpool = nvegpool_natveg
       ncphtrans    = 17
       nnphtrans    = 34
       ncphouttrans =  3
       nnphouttrans =  4
    end if
    ncgmtrans    = 18
    ncgmouttrans = 18
    ncfitrans    = 20
    ncfiouttrans = 18
    nngmtrans    = 19
    nngmouttrans = 19
    nnfitrans    = 21
    nnfiouttrans = 19
    nvegnpool    = nvegcpool + 1
    iretransn    = nvegnpool
    ioutc        = nvegcpool + 1
    ioutn        = nvegnpool + 1
    
  end subroutine clm_varpar_init
  
end module clm_varpar

! =================== EOF ============================================================
