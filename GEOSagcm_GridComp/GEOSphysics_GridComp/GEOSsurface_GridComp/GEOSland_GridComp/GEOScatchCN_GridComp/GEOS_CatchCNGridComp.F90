!  $Id$

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_CatchCNGridCompMod

!BOP
! !MODULE: GEOS_CatchCN --- ESMF gridded component implementing CatchmentCN LSM

! !DESCRIPTION:
!
!   {\tt Catch} is a gridded component to compute the energy and water
!   fluxes due to land-surface processes, using the Catchment LSM
!   of Koster et al. (2014). 
!   Koster, R. D., G. Walker, G. J. Collatz, and P. E. Thornton, 2014. 
!   Hydroclimatic controls on the means and variability of vegetation 
!   phenology and carbon uptake. J. Climate, 27, 5632-5652. doi: 
!   10.1175/JCLI-D-13-00477.1. 
!   All of its calculations are done
!   in a tile space defined by the inherited location stream.
!   It has a two-stage run method. The first stage obtains
!   drag coefficients at all the subtiles and defines
!   effective tile-mean surface quantities. The second
!   stage calls the Catchment-CN LSM. {\tt CatchCN} has no children.

!
! !USES:

  use sfclayer  ! using module that contains sfc layer code
  use ESMF
  use GEOS_Mod
  use GEOS_UtilsMod
  use DragCoefficientsMod
  use CATCHMENT_CN_MODEL
  use compute_rc_mod
  use CN_DriverMod
  USE STIEGLITZSNOW,   ONLY :                 &
       snow_albedo, StieglitzSnow_calc_tpsnow, N_CONSTIT,   &
       NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD, &
       NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD, &
       NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD, &
       NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD, &
       NUM_SSDP, NUM_SSSV, NUM_SSWT, NUM_SSSD, &
       StieglitzSnow_calc_asnow

  USE CATCH_CONSTANTS, ONLY :                 &
       N_GT           => CATCH_N_GT,          &
       N_SNOW         => CATCH_N_SNOW,        &
       RHOFS          => CATCH_SNWALB_RHOFS,  &
       SNWALB_VISMAX  => CATCH_SNWALB_VISMAX, &
       SNWALB_NIRMAX  => CATCH_SNWALB_NIRMAX, &
       SLOPE          => CATCH_SNWALB_SLOPE

  USE  clm_varpar, ONLY :                     &
       NUM_ZON, NUM_VEG, VAR_COL, VAR_PFT,    &
       CN_zone_weight, map_cat, firefac
 
  USE MAPL_BaseMod
  use MAPL_ConstantsMod,only: Tzero => MAPL_TICE, pi => MAPL_PI 
  use clm_time_manager, only: get_days_per_year, get_step_size
  use pftvarcon,        only: noveg
  USE lsm_routines,     ONLY : sibalb, catch_calc_soil_moist     

implicit none
private

! !PUBLIC MEMBER FUNCTIONS:

public SetServices

!
!EOP

integer,parameter :: FSAT=1  !  Saturated subtile
integer,parameter :: FTRN=2  !  Transition subtile
integer,parameter :: FWLT=3  !  Wilting subtile
integer,parameter :: FSNW=4  !  Snowcover subtile

integer,parameter :: NUM_SUBTILES=4

! Vegetation type as follows:
!                  1:  BROADLEAF EVERGREEN TREES
!                  2:  BROADLEAF DECIDUOUS TREES
!                  3:  NEEDLELEAF TREES
!                  4:  GROUND COVER
!                  5:  BROADLEAF SHRUBS
!                  6:  DWARF TREES (TUNDRA)
!===================================================
!ALT: we currently use only 6 types (see above)
!     in the legacy code we used to have 8 
!     (or 10 with the sea and land ice) with
!     these additional entries
!                  7:  BARE SOIL
!                  8:  DESERT

integer           :: NUM_ENSEMBLE, USE_ASCATZ0, DO_CO2SC
integer,parameter :: NTYPS = MAPL_NUMVEGTYPES

real,   parameter :: HPBL           = 1000.
real,   parameter :: MIN_VEG_HEIGHT = 0.01
real,   parameter :: Z0_BY_ZVEG     = 0.13
real,   parameter :: D0_BY_ZVEG     = 0.66

! Emissivity values from Wilber et al (1999, NATA-TP-1999-209362)
! Fu-Liou bands have been combined to Chou bands (though these are broadband only)
! IGBP veg types have been mapped to Sib-Mosaic types 
! Details in ~suarez/Emiss on cerebus

real,   parameter :: EMSVEG(NTYPS) = (/ 0.99560, 0.99000, 0.99560, 0.99320, &
                                        0.99280, 0.99180 /)
real,   parameter :: EMSBARESOIL   =    0.94120
real,   parameter :: EMSSNO        =    0.99999

! moved SURFLAY from catchment.F90 to enable run-time changes for off-line system
! - reichle, 29 Oct 2010

! real,   parameter :: SURFLAY = 20.  ! moved to GetResource in RUN2  LLT:12Jul3013

! ROOTL import from GEOS_VegdynGridComp was disabled and brought the look up table 
! in order to obtain ROOTL for primary and secondary types.

!  map catchment type into PFT
!  ---------------------------
!PFT 	Description
! 0 	bare
! 1 	needleleaf evergreen temperate tree
! 2 	needleleaf evergreen boreal tree
! 3 	needleleaf deciduous boreal tree
! 4 	broadleaf evergreen tropical tree
! 5 	broadleaf evergreen temperate tree
! 6 	broadleaf deciduous tropical tree
! 7 	broadleaf deciduous temperate tree
! 8 	broadleaf deciduous boreal tree
! 9 	broadleaf evergreen temperate shrub
! 10 	broadleaf deciduous temperate shrub [moisture + deciduous]
! 11 	broadleaf deciduous temperate shrub [moisture stress only]
! 12 	broadleaf deciduous boreal shrub
! 13 	arctic c3 grass
! 14 	cool c3 grass [moisture + deciduous]
! 15 	cool c3 grass [moisture stress only]
! 16 	warm c4 grass [moisture + deciduous]
! 17 	warm c4 grass [moisture stress only]
! 18 	crop          [moisture + deciduous]
! 19 	crop          [moisture stress only]

! Catchment types and PFT mapping:
! 
! 1:  BROADLEAF EVERGREEN TREES  => 4,5
! 2:  BROADLEAF DECIDUOUS TREES  => 6,7,8
! 3:  NEEDLELEAF TREES           => 1,2,3
! 4:  GROUND COVER               => 13-19
! 5:  BROADLEAF SHRUBS           => 9,10,11
! 6:  DWARF TREES (TUNDRA)       => 12
! 7:  BARE SOIL                  => 0
! 8:  DESERT                     => 0
! 9:  ICE                        => n/a

! index map for CLM PFTs --> catchment veg types

! pchakrab: save the logical variable OFFLINE
! Internal state and its wrapper
type T_OFFLINE_MODE
   private
   logical :: CATCH_OFFLINE
end type T_OFFLINE_MODE
type OFFLINE_WRAP
   type(T_OFFLINE_MODE), pointer :: ptr=>null()
end type OFFLINE_WRAP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for component
! !INTERFACE:

subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),intent(INOUT) :: GC
    integer, optional,  intent(  OUT) :: RC

! !DESCRIPTION:
! This version uses GEOS\_GenericSetServices, overriding
! only the run method. It also relies on MAPL\_Generic to
! handle data services. 

!EOP
!
! ErrLog Variables

    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: COMP_NAME
    integer                    :: STATUS

! Local Variables

    type(MAPL_MetaComp), pointer :: MAPL=>null()
    type(T_OFFLINE_MODE), pointer :: internal=>null()
    type(OFFLINE_WRAP) :: wrap
    integer :: OFFLINE_MODE
    logical :: is_OFFLINE
    integer :: RESTART
    integer :: DO_GOSWIM

! Begin...
! --------

! Get my name and set-up traceback handle
! ------------------------------------------------------------------------------

    Iam='SetServices'
    call ESMF_GridCompGet ( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//trim(Iam)

! pchakrab: Read CATCHMENT_OFFLINE from resource file and save
! it in the private internal state of the GridComp. It is a little
! unusual to read resource file in SetServices, but we need to know
! at this stage where we are running Catch in the offline mode or not

    allocate(internal, stat=status)
    VERIFY_(status)
    wrap%ptr => internal
    call ESMF_UserCompSetInternalState(gc, 'OfflineMode', wrap, status)

    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)
    call MAPL_GetResource ( MAPL, OFFLINE_MODE, Label="CATCHMENT_OFFLINE:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    wrap%ptr%CATCH_OFFLINE = OFFLINE_MODE /= 0

    is_OFFLINE = wrap%ptr%CATCH_OFFLINE

    call MAPL_GetResource ( MAPL, NUM_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, USE_ASCATZ0, Label="USE_ASCATZ0:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, DO_CO2SC,        Label="USE_CO2SC:",DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, DO_GOSWIM, Label="N_CONST_LAND4SNWALB:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

! Set the Run entry points
! ------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN1, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN2, RC=STATUS )
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

!BOS

!  !IMPORT STATE:

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_pressure'            ,&
         UNITS              = 'Pa'                          ,&
         SHORT_NAME         = 'PS'                          ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 

    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_air_temperature'     ,&
         UNITS              = 'K'                           ,&
         SHORT_NAME         = 'TA'                          ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 

    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_air_specific_humidity',&
         UNITS              = 'kg kg-1'                     ,&
         SHORT_NAME         = 'QA'                          ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_wind_speed'          ,&
         UNITS              = 'm s-1'                       ,&
         SHORT_NAME         = 'UU'                          ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_uwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'UWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'levellm_vwind',                     &
        UNITS              = 'm s-1',                             &
        SHORT_NAME         = 'VWINDLMTILE',                       &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'liquid_water_convective_precipitation',&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'PCU'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'liquid_water_large_scale_precipitation',&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'PLS'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'snowfall'                    ,&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'SNO'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'icefall'                     ,&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'ICE'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'freezing_rain_fall'          ,&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'FRZR'                        ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RC=STATUS  )

    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_par_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRPAR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_par_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFPAR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_nir_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFNIR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_uvr_diffuse_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DFUVR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_longwave_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'LWDNSRF'                     ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'ALW'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'linearization_of_surface_upwelling_longwave_flux',&
         UNITS              = 'W_m-2 K-1'                   ,&
         SHORT_NAME         = 'BLW'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    IF (DO_CO2SC /= 0) THEN
       call MAPL_AddImportSpec(GC,                              &
            SHORT_NAME         = 'CO2SC',                             &
            LONG_NAME          = 'CO2 Surface Concentration Bin 001', &
            UNITS              = '1e-6',                              &
            DIMS               = MAPL_DimsTileOnly,                   &
            VLOCATION          = MAPL_VLocationNone,                  &
            RC=STATUS  )
       VERIFY_(STATUS)
    ENDIF

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'leaf_area_index'             ,&
         UNITS              = '1'                           ,&
         SHORT_NAME         = 'LAI'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'greeness_fraction'           ,&
         UNITS              = '1'                           ,&
         SHORT_NAME         = 'GRN'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'evaporation'                 ,&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'EVAP'                        ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'derivative_of_evaporation_wrt_QS',&
         UNITS              = 'kg m-2 s-1'                  ,&
         SHORT_NAME         = 'DEVAP'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'upward_sensible_heat_flux'   ,&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'SH'                          ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'derivative_of_sensible_heat_wrt_Ts',&
         UNITS              = 'W m-2 K-1'                   ,&
         SHORT_NAME         = 'DSH'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_layer_height'        ,&
         UNITS              = 'm'                           ,&
         SHORT_NAME         = 'DZ'                          ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC                         ,&
        LONG_NAME          = 'vegetation_root_length'      ,&
        UNITS              = 'm'                           ,&
        SHORT_NAME         = 'ROOTL'                       ,&
        DIMS               = MAPL_DimsTileOnly             ,&
        VLOCATION          = MAPL_VLocationNone            ,&
                                                 RC=STATUS  ) 
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC                         ,&
        LONG_NAME          = 'canopy_height'               ,&
        UNITS              = 'm'                           ,&
        SHORT_NAME         = 'Z2CH'                        ,&
        DIMS               = MAPL_DimsTileOnly             ,&
        VLOCATION          = MAPL_VLocationNone            ,&
                                                 RC=STATUS  ) 
   VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         SHORT_NAME         = 'THATM',                       &
         LONG_NAME          = 'effective_surface_skin_temperature',&
         UNITS              = 'K',                           &
         DIMS               = MAPL_DimsTileOnly,             &
         VLOCATION          = MAPL_VLocationNone,            &
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                         &
         SHORT_NAME         = 'QHATM',                       &
         LONG_NAME          = 'effective_surface_specific_humidity',&
         UNITS              = 'kg kg-1',                     &
         DIMS               = MAPL_DimsTileOnly,             &
         VLOCATION          = MAPL_VLocationNone,            &
                                                  RC=STATUS  )
    VERIFY_(STATUS)
 
    call MAPL_AddImportSpec(GC,                         &
         SHORT_NAME         = 'CTATM',                       &
         LONG_NAME          = 'surface_exchange_coefficient_for_heat', &
         UNITS              = 'kg m-2 s-1',                  &
         DIMS               = MAPL_DimsTileOnly,             &
         VLOCATION          = MAPL_VLocationNone,            &
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                         &
         SHORT_NAME         = 'CQATM',                       &
         LONG_NAME          = 'surface_exchange_coefficient_for_moisture', &
         UNITS              = 'kg m-2 s-1',                  &
         DIMS               = MAPL_DimsTileOnly,             &
         VLOCATION          = MAPL_VLocationNone,            &
                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                ,&
       SHORT_NAME = 'ASCATZ0'                                 ,&
       LONG_NAME  = 'ASCAT_roughness_length'		      ,&
       UNITS      = 'm'                                       ,&
       DIMS       = MAPL_DimsTileOnly                         ,&
       VLOCATION  = MAPL_VLocationNone                        ,&
       RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddImportSpec(GC                             ,&
        SHORT_NAME         = 'NDVI'                        ,& 
        LONG_NAME          = 'normalized_difference_vegetation_index' ,&
        UNITS              = '1'                           ,&        
        DIMS               = MAPL_DimsTileOnly             ,&
        VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'dust_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'dust_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'dust_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'dust_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'DUSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_DUSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
 
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'black_carbon_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'BCSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_BCSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'organic_carbon_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'OCSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_OCSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sulfate_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sulfate_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sulfate_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sulfate_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SUSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SUSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_dry_depos_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSDP',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSDP/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_wet_depos_conv_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSSV',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSSV/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_wet_depos_ls_scav_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSWT',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSWT/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                         &
         LONG_NAME          = 'sea_salt_gravity_sett_all_bins', &
         UNITS              = 'kg m-2 s-1',                  &
         SHORT_NAME         = 'SSSD',                        &
         DIMS               = MAPL_DimsTileOnly,             &
         UNGRIDDED_DIMS     = (/NUM_SSSD/),                  &
         VLOCATION          = MAPL_VLocationNone,            &
         RC=STATUS  ) 
    VERIFY_(STATUS)

!  !INTERNAL STATE:

! if is_offline, some variables ( in the last) are not required
  if (is_OFFLINE) then
     RESTART = MAPL_RestartSkip
  else
     RESTART = MAPL_RestartRequired
  endif

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'topo_baseflow_param_1'     ,&
    UNITS              = 'kg m-4'                    ,&
    SHORT_NAME         = 'BF1'                       ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'topo_baseflow_param_2'     ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'BF2'                       ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'topo_baseflow_param_3'     ,&
    UNITS              = 'log(m)'                    ,&
    SHORT_NAME         = 'BF3'                       ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'max_rootzone_water_content',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'VGWMAX'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'moisture_threshold'        ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CDCR1'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'max_water_content'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CDCR2'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'saturated_matric_potential',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'PSIS'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'clapp_hornberger_b'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BEE'                       ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_porosity'             ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'POROS'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'wetness_at_wilting_point'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WPWET'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'sfc_sat_hydraulic_conduct' ,&
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'COND'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'vertical_transmissivity'   ,&
    UNITS              = 'm-1'                       ,&
    SHORT_NAME         = 'GNU'                       ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'wetness_param_1'           ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARS1'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'wetness_param_2'           ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARS2'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'wetness_param_3'           ,&
    UNITS              = 'm+4 kg-2'                  ,&
    SHORT_NAME         = 'ARS3'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'shape_param_1'             ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARA1'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'shape_param_2'             ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ARA2'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'shape_param_3'             ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARA3'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'shape_param_4'             ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ARA4'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_1'         ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARW1'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_2'         ,&
    UNITS              = 'm+2 kg-1'                  ,&
    SHORT_NAME         = 'ARW2'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_3'         ,&
    UNITS              = 'm+4 kg-2'                  ,&
    SHORT_NAME         = 'ARW3'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'min_theta_param_4'         ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ARW4'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_1'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSA1'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_2'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSA2'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_3'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSB1'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_4'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TSB2'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_5'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ATAU'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'water_transfer_param_6'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BTAU'                      ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'vegetation_type'           ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ITY'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_VEG/)                 ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'vegetation_fraction'       ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FVG'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_VEG/)                 ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'canopy_temperature'        ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TC'                        ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'canopy_specific_humidity'  ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QC'                        ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'ground_temperature'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TG'                        ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'interception_reservoir_capac',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CAPAC'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'catchment_deficit'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CATDEF'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'root_zone_excess'          ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'RZEXC'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'surface_excess'            ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SRFEXC'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_1' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT1'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_2' ,&
    UNITS              = 'J_m-2'                     ,&
    SHORT_NAME         = 'GHTCNT2'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_3' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT3'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_4' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT4'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_5' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT5'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_heat_content_layer_6' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'GHTCNT6'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'mean_catchment_temp_incl_snw',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TSURF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'snow_mass_layer_1'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN1'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'snow_mass_layer_2'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN2'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'snow_mass_layer_3'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WESNN3'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'heat_content_snow_layer_1' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN1'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'heat_content_snow_layer_2' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN2'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'heat_content_snow_layer_3' ,&
    UNITS              = 'J m-2'                     ,&
    SHORT_NAME         = 'HTSNNN3'                   ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'snow_depth_layer_1'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN1'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'snow_depth_layer_2'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN2'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'snow_depth_layer_3'        ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZN3'                    ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'surface_heat_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CH'                        ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = RESTART                     ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'surface_momentum_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CM'                        ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = RESTART                     ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'surface_moisture_exchange_coffiecient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CQ'                        ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = RESTART                     ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'subtile_fractions'         ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FR'                        ,&
    DIMS               = MAPL_DimsTileTile           ,&
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = RESTART                     ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC,                  &
    SHORT_NAME         = 'WW',                        &
    LONG_NAME          = 'vertical_velocity_scale_squared', &
    UNITS              = 'm+2 s-2',                   &
    DIMS               = MAPL_DimsTileTile,           &
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone,          &
    RESTART            = RESTART                     ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC,                  &
    SHORT_NAME         = 'DCH',                        &
    LONG_NAME          = 'ch difference, optional in louissurface', &
    UNITS              = '1',                   &
    DIMS               = MAPL_DimsTileTile,           &
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone,          &
    RESTART            = MAPL_RestartSkip            ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC,                  &
    SHORT_NAME         = 'DCQ',                        &
    LONG_NAME          = 'cq difference, optional in louissurface', &
    UNITS              = '1',                   &
    DIMS               = MAPL_DimsTileTile,           &
    NUM_SUBTILES       = NUM_SUBTILES                ,&
    VLOCATION          = MAPL_VLocationNone,          &
    RESTART            = MAPL_RestartSkip            ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'catchment_tile_id'         ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'TILE_ID'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CLM_nitrogen_deposition'   ,&
    UNITS              = 'g m-2 s-1'                 ,&
    SHORT_NAME         = 'NDEP'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'cli_2m_T_(MERRA2)'      ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'CLI_T2M'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'MODIS soil albedo vis dir' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BGALBVR'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'MODIS soil albedo vis dif' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BGALBVF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'MODIS soil albedo nir dir' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BGALBNR'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'MODIS soil albedo nir dif' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BGALBNF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'column_rst_vars'           ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNCOL'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_ZON*VAR_COL/)         ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'PFT_rst_vars'              ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNPFT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_ZON*NUM_VEG*VAR_PFT/) ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for ground temp'    ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TGWM'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_ZON/)                 ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for soil moisture'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'RZMM'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_ZON/)                 ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for sfc soil moist' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'SFMCM'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for baseflow'       ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'BFLOWM'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for total water'    ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'TOTWATM'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for air temperature',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TAIRM'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for soil temp'      ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPM'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN summing counter'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNSUM'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for sunlit photosyn',&
    UNITS              = 'umol m-2 s-1'              ,&
    SHORT_NAME         = 'PSNSUNM'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_VEG,NUM_ZON/)         ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for shaded photosyn',&
    UNITS              = 'umol m-2 s-1'              ,&
    SHORT_NAME         = 'PSNSHAM'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    UNGRIDDED_DIMS     = (/NUM_VEG,NUM_ZON/)         ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for snow depth'     ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNDZM'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CN sum for area snow cover',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOWM'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  !---------- GOSWIM snow impurity related variables ----------

  if (DO_GOSWIM /= 0) then 
  
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_1'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU001'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_2'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU002'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_3'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU003'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_4'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU004'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'dust_mass_in_snow_bin_5'   ,&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RDU005'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'hydrophobic_black_carbon_mass_in_snow_bin_1',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RBC001'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'hydrophilic_black_carbon_mass_in_snow_bin_2',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'RBC002'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'hydrophobic_organic_carbon_mass_in_snow_bin_1',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'ROC001'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                  ,&
          LONG_NAME          = 'hydrophilic_organic_carbon_mass_in_snow_bin_2',&
          UNITS              = 'kg m-2'                    ,&
          SHORT_NAME         = 'ROC002'                    ,&
          FRIENDLYTO         = trim(COMP_NAME)             ,&
          DIMS               = MAPL_DimsTileOnly           ,&
          UNGRIDDED_DIMS     = (/N_SNOW/)                  ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartOptional        ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)
     
  endif

  !  !EXPORT STATE:
  
  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'evaporation'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'EVAPOUT'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                     &
     LONG_NAME          = 'sublimation'               ,&
     UNITS              = 'kg m-2 s-1'                ,&
     SHORT_NAME         = 'SUBLIM'                    ,&
     DIMS               = MAPL_DimsTileOnly           ,&
     VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
     VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'upward_sensible_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHOUT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'runoff_flux'               ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNOFF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'interception_loss_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPINT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baresoil_evap_energy_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'transpiration_energy_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPVEG'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_ice_evaporation_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPICE'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil moisture in Upper 10cm'     ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WAT10CM'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'totoal soil moisture'      ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'WATSOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil frozen water content' ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'ICESOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_evaporation_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSNO'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baseflow_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'BASEFLOW'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'overland_runoff_including_throughflow'  ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNSURF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowmelt_flux'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'SMELT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_outgoing_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'HLWUP'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_longwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'LWNDSRF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
    VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                    ,&
    LONG_NAME          = 'surface_net_downward_shortwave_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SWNDSRF'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
    VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'total_latent_energy_flux'  ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'HLATN'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'rainwater_infiltration_flux',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'QINFIL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'areal_fraction_saturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'AR1'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'areal_fraction_transpiration_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'AR2'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'root_zone_equilibrium_moisture',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'RZEQ'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ground_energy_flux'        ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'GHFLX'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ave_catchment_temp_incl_snw',& 
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSURF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_top_snow_layer',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSNOW'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_unsaturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPUNST'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_saturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'temperature_wilted_zone'   ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPWLT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  if (is_OFFLINE) then
     call MAPL_AddInternalSpec(GC,                       &
       LONG_NAME          = 'fractional_area_of_land_snowcover',&
       UNITS              = '1'                         ,&
       SHORT_NAME         = 'ASNOW'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartSkip            ,&
       FRIENDLYTO         = trim(COMP_NAME)             ,&
                                           RC=STATUS  )
     VERIFY_(STATUS)
  else
  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_land_snowcover',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOW'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  endif

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'downward_heat_flux_into_snow',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'SHSNOW'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'averaged_snow_temperature' ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'AVETSNOW'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_saturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_unsaturated_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRUST'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_wilting_zone',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'FRWLT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_mass'                 ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'SNOWMASS'                  ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snow_depth_in_snow_covered_area' ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNOWDP'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_soil_wetness'      ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'root_zone_soil_wetness'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'ave_prof_soil__moisture'   ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET3'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_surface_layer'       ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCSF'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_root_zone'           ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCRZ'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'water_ave_prof'            ,&
    UNITS              = 'm3 m-3'                   ,&
    SHORT_NAME         = 'WCPR'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_1' ,&
    UNITS              = 'C'                         ,&
    SHORT_NAME         = 'TP1'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_2' ,&
    UNITS              = 'C'                         ,&
    SHORT_NAME         = 'TP2'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_3' ,&
    UNITS              = 'C'                         ,&
    SHORT_NAME         = 'TP3'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_4' ,&
    UNITS              = 'C'                         ,&
    SHORT_NAME         = 'TP4'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_5' ,&
    UNITS              = 'C'                         ,&
    SHORT_NAME         = 'TP5'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperatures_layer_6' ,&
    UNITS              = 'C'                         ,&
    SHORT_NAME         = 'TP6'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  if (is_OFFLINE) then
     call MAPL_AddInternalSpec(GC,                       &
        LONG_NAME          = 'surface_emissivity'        ,&
        UNITS              = '1'                         ,&
        SHORT_NAME         = 'EMIS'                      ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RESTART            = MAPL_RestartSkip            ,&
        FRIENDLYTO         = trim(COMP_NAME)             ,&
                                           RC=STATUS  )
     VERIFY_(STATUS)
  else
  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_emissivity'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'EMIS'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  endif

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_visible_beam',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBVR'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_visible_diffuse',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBVF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_near_infrared_beam',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBNR'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_albedo_near_infrared_diffuse',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBNF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_surface_skin_temperature',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'DELTS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_surface_specific_humidity',&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'DELQS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_evaporation'        ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'DELEVAP'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'change_upward_sensible_energy_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'DELSH'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_skin_temperature'  ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TST'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'land_surface_skin_temperature'  ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'LST'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_specific_humidity' ,&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QST'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'turbulence_surface_skin_temperature',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TH'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'turbulence_surface_skin_specific_hum',&
    UNITS              = 'kg kg-1'                   ,&
    SHORT_NAME         = 'QH'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_heat_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CHT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_momentum_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CMT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_moisture_exchange_coefficient',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CQT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'neutral_drag_coefficient'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_bulk_richardson_number',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'RIT'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_roughness'         ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'Z0'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT2M',                     &
        LONG_NAME          = 'temperature 2m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ2M',                     &
        LONG_NAME          = 'humidity 2m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU2M',                    &
        LONG_NAME          = 'zonal 2m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV2M',                    &
        LONG_NAME          = 'meridional 2m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOT10M',                     &
        LONG_NAME          = 'temperature 10m wind from MO sfc', &
        UNITS              = 'K',                         &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOQ10M',                     &
        LONG_NAME          = 'humidity 10m wind from MO sfc',    &
        UNITS              = 'kg kg-1',                   &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU10M',                    &
        LONG_NAME          = 'zonal 10m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV10M',                    &
        LONG_NAME          = 'meridional 10m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOU50M',                    &
        LONG_NAME          = 'zonal 50m wind from MO sfc',&
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                    &
        SHORT_NAME         = 'MOV50M',                    &
        LONG_NAME          = 'meridional 50m wind from MO sfc', &
        UNITS              = 'm s-1',                     &
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone,          &
                                               RC=STATUS  )
     VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_roughness_for_heat',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'Z0H'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'zero_plane_displacement_height',&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'D0'                        ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GUST',                      &
    LONG_NAME          = 'gustiness',                 &
    UNITS              = 'm s-1',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'VENT',                      &
    LONG_NAME          = 'surface_ventilation_velocity',&
    UNITS              = 'm s-1',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


     call MAPL_AddExportSpec(GC,                             &
        SHORT_NAME         = 'ACCUM',                             &
        LONG_NAME          = 'net_ice_accumulation_rate',         &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'EVLAND',                    &
    LONG_NAME          = 'Evaporation_land',          &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'PRLAND',                    &
    LONG_NAME          = 'Total_precipitation_land',  &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SNOLAND',                   &
    LONG_NAME          = 'snowfall_land',             &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DRPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_par_beam_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DFPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_par_diffuse_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LHSNOW',                    &
    LONG_NAME          = 'Latent_heat_flux_snow',     &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWNETSNOW',                    &
    LONG_NAME          = 'Net_shortwave_snow',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWUPSNOW',                    &
    LONG_NAME          = 'Net_longwave_snow',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWDNSNOW',                    &
    LONG_NAME          = 'Net_longwave_snow',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TCSORIG',                   &
    LONG_NAME          = 'Input_tc_for_snow',         &
    UNITS              = 'K',                         &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TPSN1IN',                   &
    LONG_NAME          = 'Input_temp_of_top_snow_lev',&
    UNITS              = 'K',                         &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TPSN1OUT',                  &
    LONG_NAME          = 'Output_temp_of_top_snow_lev',&
    UNITS              = 'K',                         &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHSNOW',                    &
    LONG_NAME          = 'Ground_heating_snow',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LHLAND',                    &
    LONG_NAME          = 'Latent_heat_flux_land',     &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SHLAND',                    &
    LONG_NAME          = 'Sensible_heat_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWLAND',                    &
    LONG_NAME          = 'Net_shortwave_land',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWDOWNLAND',                &
    LONG_NAME          = 'Incident_shortwave_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWLAND',                    &
    LONG_NAME          = 'Net_longwave_land',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHLAND',                    &
    LONG_NAME          = 'Ground_heating_land',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHTSKIN',                   &
    LONG_NAME          = 'Ground_heating_skin_temp',  &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SMLAND',                    &
    LONG_NAME          = 'Snowmelt_flux_land',        &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TWLAND',                    &
    LONG_NAME          = 'Avail_water_storage_land',  &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TELAND',                    &
    LONG_NAME          = 'Total_energy_storage_land', &
    UNITS              = 'J m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'TSLAND',                    &
    LONG_NAME          = 'Total_snow_storage_land',   &
    UNITS              = 'kg m-2',                    &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DWLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_water',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DHLAND',                    &
    LONG_NAME          = 'rate_of_change_of_total_land_energy',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPLAND',                    &
    LONG_NAME          = 'rate_of_spurious_land_energy_source',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPWATR',                    &
    LONG_NAME          = 'rate_of_spurious_land_water_source',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPSNOW',                    &
    LONG_NAME          = 'rate_of_spurious_snow_energy',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'vegetation_type'           ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ITY'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_exposed_leaf-area_index',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNLAI'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_leaf-area_index'  ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNTLAI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_exposed_stem-area_index',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNSAI'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_carbon'           ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNTOTC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_vegetation_carbon',&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNVEGC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_root_carbon'      ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNROOT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_net_primary_production' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNPP'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_gross_primary_production',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNGPP'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_soil_respiration' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNSR'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_net_ecosystem_exchange' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNEE'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'abstract_C_pool_to_meet_excess_MR_demand' ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNXSMR'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_added_to_maintain_positive_C' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNADD'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_carbon_loss_to_fire'    ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNLOSS'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_fractional_area_burn_rate' ,&
    UNITS              = 's-1'                       ,&
    SHORT_NAME         = 'CNBURN'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'absorbed_PAR'              ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'PARABS'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'incident_PAR'              ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'PARINC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'saturated_stomatal_conductance' ,&
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SCSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'unstressed_stomatal_conductance' ,&
    UNITS              = 'm s-1'                     ,&
    SHORT_NAME         = 'SCUNS'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'transpiration coefficient' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'BTRANT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'solar induced fluorescence',&
    UNITS              = 'umol m-2 sm s-1'           ,&
    SHORT_NAME         = 'SIF'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'fire season length'        ,&
    UNITS              = 'days'                      ,&
    SHORT_NAME         = 'CNFSEL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU001'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU002'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_3',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU003'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_4',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU004'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_5',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTDU005'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTBC001'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTBC002'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_1',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTOC001'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_2',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RMELTOC002'                ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
  VERIFY_(STATUS)

!EOS

    call MAPL_TimerAdd(GC,    name="RUN1"  ,RC=STATUS)
    VERIFY_(STATUS)
    if (is_OFFLINE) then
       call MAPL_TimerAdd(GC,    name="-RUN0"  ,RC=STATUS)
       VERIFY_(status)
    end if
    call MAPL_TimerAdd(GC,    name="-SURF" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-CATCH" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-ALBEDO" ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final method
! ---------------------------------

    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!BOP
! !IROUTINE: RUN1 -- First Run stage for the catchment component
! !INTERFACE:

subroutine RUN1 ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),intent(inout) :: GC     !Gridded component
    type(ESMF_State),   intent(inout) :: IMPORT !Import state
    type(ESMF_State),   intent(inout) :: EXPORT !Export state
    type(ESMF_Clock),   intent(inout) :: CLOCK  !The clock
    integer,optional,   intent(out  ) :: RC     !Error code:

! !DESCRIPTION: Does the cds computation and roughness length
!EOP
! ErrLog Variables

    character(len=ESMF_MAXSTR) :: IAm
    integer :: STATUS

    character(len=ESMF_MAXSTR) :: COMP_NAME

! Locals

    type(MAPL_MetaComp),pointer :: MAPL
    type(ESMF_State)                :: INTERNAL
    type(ESMF_Alarm)                :: ALARM
    type(ESMF_Config)               :: CF
    type(ESMF_VM)                   :: VM

! -----------------------------------------------------
! IMPORT Pointers
! ----------------------------------------------------  -

    real, dimension(:),     pointer :: PS
    real, dimension(:),     pointer :: TA
    real, dimension(:),     pointer :: QA
    real, dimension(:),     pointer :: UU
    real, pointer, dimension(:)    :: UWINDLMTILE
    real, pointer, dimension(:)    :: VWINDLMTILE
    real, dimension(:),     pointer :: DZ
    real, dimension(:),     pointer :: LAI
    real, dimension(:),     pointer :: Z2CH
    real, dimension(:),     pointer :: PCU
    real, dimension(:),     pointer :: ASCATZ0
    real, dimension(:),     pointer :: NDVI

! -----------------------------------------------------
! INTERNAL Pointers
! -----------------------------------------------------

    real, dimension(:,:), pointer :: ITY
    real, dimension(:,:), pointer :: FVG
    real, dimension(:,:), pointer :: TC
    real, dimension(:,:), pointer :: QC
    real, dimension(:,:), pointer :: CH
    real, dimension(:,:), pointer :: CM
    real, dimension(:,:), pointer :: CQ
    real, dimension(:,:), pointer :: FR
    real, dimension(:,:), pointer :: WW
    real, dimension(:,:), pointer :: cncol
    real, dimension(:,:), pointer :: cnpft
    real, dimension(:,:), pointer :: DCH
    real, dimension(:,:), pointer :: DCQ

! -----------------------------------------------------
! EXPORT Pointers
! -----------------------------------------------------

    real, dimension(:),   pointer :: TH
    real, dimension(:),   pointer :: QH
    real, dimension(:),   pointer :: CHT
    real, dimension(:),   pointer :: CMT
    real, dimension(:),   pointer :: CQT
    real, dimension(:),   pointer :: CNT
    real, dimension(:),   pointer :: RIT
    real, dimension(:),   pointer :: Z0
    real, dimension(:),   pointer :: Z0H
    real, dimension(:),   pointer :: D0
    real, dimension(:),   pointer :: GST
    real, dimension(:),   pointer :: VNT
   real, pointer, dimension(:  )  :: MOT2M
   real, pointer, dimension(:  )  :: MOQ2M
   real, pointer, dimension(:  )  :: MOU2M
   real, pointer, dimension(:  )  :: MOV2M
   real, pointer, dimension(:  )  :: MOT10M
   real, pointer, dimension(:  )  :: MOQ10M
   real, pointer, dimension(:  )  :: MOU10M
   real, pointer, dimension(:  )  :: MOV10M
   real, pointer, dimension(:  )  :: MOU50M
   real, pointer, dimension(:  )  :: MOV50M
    real, dimension(:),   pointer :: ITYO


! From old bucket version of CDS calculation
! ------------------------------------------

    integer             :: N
    integer             :: NT
    real,   allocatable :: UCN(:)
    real,   allocatable :: TVA(:)
    real,   allocatable :: TVS(:)
    real,   allocatable :: URA(:)
    real,   allocatable :: UUU(:)
    real,   allocatable :: ZVG(:)
    real,   allocatable :: DZE(:)
    real,   allocatable :: D0T(:)
    real,   allocatable :: CHX(:)
    real,   allocatable :: CQX(:)
    real,   allocatable :: CN(:)
    real,   allocatable :: RE(:)
    real,   allocatable :: ZT(:)
    real,   allocatable :: ZQ(:)
    integer,allocatable :: VEG1(:)
    integer,allocatable :: VEG2(:)
    real,   allocatable :: FVG1(:)
    real,   allocatable :: FVG2(:)
    real,   allocatable :: Z0T(:,:)
   real, allocatable              :: U50M (:)
   real, allocatable              :: V50M (:)
   real, allocatable              :: T10M (:)
   real, allocatable              :: Q10M (:)
   real, allocatable              :: U10M (:)
   real, allocatable              :: V10M (:)
   real, allocatable              :: T2M (:)
   real, allocatable              :: Q2M (:)
   real, allocatable              :: U2M (:)
   real, allocatable              :: V2M (:)
   real, allocatable              :: RHOH(:)
   real, allocatable              :: VKH(:)
   real, allocatable              :: VKM(:)
   real, allocatable              :: USTAR(:)
   real, allocatable              :: XX(:)
   real, allocatable              :: YY(:)
   real, allocatable              :: CU(:)
   real, allocatable              :: CT(:)
   real, allocatable              :: RIB(:)
   real, allocatable              :: ZETA(:)
   real, allocatable              :: WS(:)
   integer, allocatable           :: IWATER(:)
   real, allocatable              :: PSMB(:)
   real, allocatable              :: PSL(:)
   integer                        :: niter

   integer                        :: CHOOSEMOSFC
   integer                        :: CHOOSEZ0
   real                           :: SCALE4Z0

! gkw: for CN model
! -----------------
    integer, parameter :: nveg  = num_veg ! number of vegetation types
    integer, parameter :: nzone = num_zon ! number of stress zones

    integer, allocatable :: ityp(:,:,:)
    real,    allocatable :: fveg(:,:,:), elai(:,:,:), esai(:,:,:), wtzone(:,:), lai1(:), lai2(:), wght(:)

    integer :: nv, nz, ib
    real    :: bare
    logical, save :: first = .true.
    integer*8, save :: istep = 1 ! gkw: legacy variable from offline

  ! Offline mode

   type(OFFLINE_WRAP)             :: wrap
   logical                        :: is_OFFLINE

!=============================================================================
! Begin...
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
! Get the target component's name and set-up traceback handle.
! ------------------------------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//"::RUN1"

    ! Get component's offline mode from its pvt internal state
    call ESMF_UserCompGetInternalState(gc, 'OfflineMode', wrap, status)
    VERIFY_(status)
    is_OFFLINE = wrap%ptr%CATCH_OFFLINE

    call ESMF_VMGetCurrent ( VM, RC=STATUS ) 
    ! if (MAPL_AM_I_Root(VM)) print *, trim(Iam)//'::OFFLINE mode: ', is_OFFLINE
    
! Get my internal MAPL_Generic state
! ----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN1")

! Get parameters from generic state
! ---------------------------------

    call MAPL_Get ( MAPL                          ,&
                                INTERNAL_ESMF_STATE=INTERNAL   ,&
                                                      RC=STATUS )
    VERIFY_(STATUS)

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call MAPL_GetResource ( MAPL, CHOOSEMOSFC, Label="CHOOSEMOSFC:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, CHOOSEZ0, Label="CHOOSEZ0:", DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, SCALE4Z0, Label="SCALE4Z0:", DEFAULT=0.5, RC=STATUS)
    VERIFY_(STATUS)

! Pointers to inputs
!-------------------

   call MAPL_GetPointer(IMPORT,UU     , 'UU'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,UWINDLMTILE     , 'UWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,VWINDLMTILE     , 'VWINDLMTILE'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,DZ     , 'DZ'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,TA     , 'TA'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,QA     , 'QA'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PS     , 'PS'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,LAI    , 'LAI'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,Z2CH   , 'Z2CH'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,PCU    , 'PCU'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,ASCATZ0, 'ASCATZ0',    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(IMPORT,NDVI   , 'NDVI'   ,    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to internals
!----------------------
 
   call MAPL_GetPointer(INTERNAL,ITY  , 'ITY'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,FVG  , 'FVG'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,TC   , 'TC'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,QC   , 'QC'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,FR   , 'FR'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CH   , 'CH'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CM   , 'CM'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CQ   , 'CQ'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,WW   , 'WW'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CNCOL  ,'CNCOL'  ,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CNPFT  ,'CNPFT'  ,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,DCH  , 'DCH'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,DCQ  , 'DCQ'     ,    RC=STATUS)
   VERIFY_(STATUS)

! Pointers to outputs
!--------------------

   call MAPL_GetPointer(EXPORT,QH    , 'QH'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,TH    , 'TH'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CHT   , 'CHT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CMT   , 'CMT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CQT   , 'CQT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,CNT   , 'CNT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,RIT   , 'RIT'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0    , 'Z0'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,Z0H   , 'Z0H'     ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,D0    , 'D0'      ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,GST   , 'GUST'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,VNT   , 'VENT'    ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOT2M, 'MOT2M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOQ2M, 'MOQ2M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU2M, 'MOU2M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV2M, 'MOV2M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOT10M, 'MOT10M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOQ10M, 'MOQ10M'   ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU10M, 'MOU10M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV10M, 'MOV10M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOU50M, 'MOU50M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,MOV50M, 'MOV50M'  ,    RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,ITYO  , 'ITY'     ,    RC=STATUS)
   VERIFY_(STATUS)

   NT = size(TA)

   allocate(TVA(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(TVS(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(URA(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UUU(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VEG1(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VEG2(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(FVG1(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(FVG2(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(DZE(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZVG(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Z0T(NT,NUM_SUBTILES),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(D0T(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CHX(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CQX(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RE (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CN (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZT (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZQ (NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(UCN(NT),STAT=STATUS)
   VERIFY_(STATUS)
   allocate(T2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Q2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(v2M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(T10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(Q10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(v10M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(U50M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(v50M (NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RHOH(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PSMB(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(PSL(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VKH(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(VKM(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(USTAR(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(XX(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(YY(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CU(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(CT(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(RIB(NT)  ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(ZETA(NT) ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(WS(NT)   ,STAT=STATUS)
   VERIFY_(STATUS)
   allocate(IWATER(NT),STAT=STATUS)
   VERIFY_(STATUS)

   allocate(   ityp(nt,nveg,nzone) )
   allocate(   fveg(nt,nveg,nzone) )
   allocate( wtzone(nt,nzone) )
   allocate(   elai(nt,nveg,nzone) )
   allocate(   esai(nt,nveg,nzone) )

   allocate ( lai1(nt) )
   allocate ( lai2(nt) )
   allocate ( wght(nt) )

!  Vegetation types used to index into tables
!--------------------------------------------

   where(ITY(:,1) > 0.)
     VEG1 = map_cat(nint(ITY(:,1)))  ! gkw: primary   CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
    elsewhere
     VEG1 = map_cat(nint(ITY(:,2)))  ! gkw: primary   CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   where(ITY(:,3) > 0.)
     VEG2 = map_cat(nint(ITY(:,3)))  ! gkw: secondary CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
    elsewhere
     VEG2 = map_cat(nint(ITY(:,4)))  ! gkw: secondary CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   ASSERT_((count(VEG1>NTYPS.or.VEG1<1)==0))
   ASSERT_((count(VEG2>NTYPS.or.VEG2<1)==0))

   ! At this point, bare soil is not allowed in CatchCN. FVEG in BCs 
   ! files do not have bare soil either. However, at times, tiny bare 
   ! fractions appear due to truncation. We add that tiny fraction to the
   ! largest of the 4 fractions and ensure bare is zero. (Sarith 3/3/16)

   DO N = 1, NT
      BARE = 1.      

      DO NV = 1, NVEG
         BARE = BARE - FVG(N,NV)! subtract vegetated fractions 
      END DO

      if (BARE /= 0.) THEN
         IB = MAXLOC(FVG(N,:),1)
         FVG (N,IB) = FVG(N,IB) + BARE ! This also corrects cases sum gt 1.
      ENDIF
      
   END DO
   
   FVG1 = fvg(:,1) + fvg(:,2)  ! gkw: primary   vegetation fraction
   FVG2 = fvg(:,3) + fvg(:,4)  ! gkw: secondary vegetation fraction

! set CLM CN PFT & fraction, set carbon zone weights
! --------------------------------------------------
   do nz = 1,nzone
     ityp(:,:,nz) = nint(ity(:,:))
     fveg(:,:,nz) = fvg(:,:)
     wtzone(:,nz) = CN_zone_weight(nz)
   end do

! initialize CN model and transfer restart variables on startup
! -------------------------------------------------------------
   if(first) then
     call CN_init(istep,nt,nveg,nzone,ityp,fveg,cncol,var_col,cnpft,var_pft)     
     first = .false.
   endif

   ! For the OFFLINE case, first update some diagnostic vars
   if (is_OFFLINE) then
      call MAPL_TimerOn(MAPL, "-RUN0")
      call RUN0(gc, import, export, clock, rc)
      call MAPL_TimerOff(MAPL, "-RUN0")
   end if

! obtain LAI from previous time step (from CN model)
! --------------------------------------------------
   call get_CN_LAI(nt,nveg,nzone,ityp,fveg,elai,esai=esai)

   lai1 = 0.
   wght = 0.
   do nz = 1,nzone
     do nv = 1,2
       lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
       wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
     end do
   end do
   lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

   lai2 = 0.
   wght = 0.
   do nz = 1,nzone
     do nv = 3,4
       lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
       wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
     end do
   end do
   lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

   lai = fvg1*lai1 + fvg2*lai2   ! gkw: this is a VEGDYN import

   deallocate ( ityp )
   deallocate ( fveg )
   deallocate ( elai )
   deallocate ( esai )
   deallocate ( wtzone )

!  Clear the output tile accumulators
!------------------------------------
   
   CHX = 0.0
   CQX = 0.0

   if(associated(TH )) TH  = 0.0
   if(associated(QH )) QH  = 0.0
   if(associated(CMT)) CMT = 0.0
   if(associated(CNT)) CNT = 0.0
   if(associated(RIT)) RIT = 0.0
   if(associated(Z0H)) Z0H = 0.0
   if(associated(GST)) GST = 0.0
   if(associated(VNT)) VNT = 0.0
   if(associated(MOU50M)) MOU50M = 0.0
   if(associated(MOV50M)) MOV50M = 0.0
   if(associated(MOT10M)) MOT10M = 0.0
   if(associated(MOQ10M)) MOQ10M = 0.0
   if(associated(MOU10M)) MOU10M = 0.0
   if(associated(MOV10M)) MOV10M = 0.0
   if(associated( MOT2M))  MOT2M = 0.0
   if(associated( MOQ2M))  MOQ2M = 0.0
   if(associated( MOU2M))  MOU2M = 0.0
   if(associated( MOV2M))  MOV2M = 0.0

   SUBTILES: do N=1,NUM_SUBTILES

!  Effective vegetation height. In catchment, LAI dependence 
!   includes the effect of partially vegetated areas,
!   as well as the phenology of the deciduous types. These
!   effects will be separated in future formulations.

   ZVG  = fvg1*(Z2CH - SCALE4Z0*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI1)) + &
          fvg2*(Z2CH - SCALE4Z0*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI2)) 


!  For now roughnesses and displacement heights
!   are the same for all subtiles.

   Z0T(:,N)  = Z0_BY_ZVEG*ZVG
   IF (USE_ASCATZ0 == 1) THEN
      WHERE (NDVI <= 0.2)
         Z0T(:,N)  = ASCATZ0
      END WHERE
   ENDIF
   D0T  = D0_BY_ZVEG*ZVG

   DZE  = max(DZ - D0T, 10.)

   if(associated(Z0 )) Z0  = Z0T(:,N)
   if(associated(D0 )) D0  = D0T

!  Compute the three surface exchange coefficients
!-------------------------------------------------

! Choose sfc layer: if CHOOSEMOSFC is 1, choose helfand MO,
!                   if CHOOSEMOSFC is 0 (default), choose louis

   call MAPL_TimerOn(MAPL,"-SURF")
   if(CHOOSEMOSFC.eq.0) then
   WW(:,N) = 0.
   CM(:,N) = 0.

    call louissurface(3,N,UU,WW,PS,TA,TC,QA,QC,PCU,LAI,Z0T,DZE,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE,DCH,DCQ)

   elseif (CHOOSEMOSFC.eq.1)then
  
    niter = 6   ! number of internal iterations in the helfand MO surface layer routine
    IWATER = 3
  
    PSMB = PS * 0.01            ! convert to MB
! Approximate pressure at top of surface layer: hydrostatic, eqn of state using avg temp and press
    PSL = PSMB * (1. - (DZE*MAPL_GRAV)/(MAPL_RGAS*(TA+TC(:,N)) ) ) /   &
               (1. + (DZE*MAPL_GRAV)/(MAPL_RGAS*(TA+TC(:,N)) ) )
  
    CALL helfsurface( UWINDLMTILE,VWINDLMTILE,TA,TC(:,N),QA,QC(:,N),PSL,PSMB,Z0T(:,N),lai,  &
                      IWATER,DZE,niter,nt,RHOH,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,  &
                      t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0)
  
    CM(:,N)  = VKM
    CH(:,N)  = VKH
    CQ(:,N)  = VKH
  
    CN = (MAPL_KARMAN/ALOG(DZE/Z0T(:,N) + 1.0)) * (MAPL_KARMAN/ALOG(DZE/Z0T(:,N) + 1.0))
    ZT = Z0T(:,N)
    ZQ = Z0T(:,N)
    RE = 0.
    UUU = UU  
    UCN = 0.
  
!  Aggregate to tiles for MO only diagnostics
!--------------------------------------------
      if(associated(MOU50M))MOU50M = MOU50M + U50M(:)*FR(:,N)
      if(associated(MOV50M))MOV50M = MOV50M + V50M(:)*FR(:,N)
      if(associated(MOT10M))MOT10M = MOT10M + T10M(:)*FR(:,N)
      if(associated(MOQ10M))MOQ10M = MOQ10M + Q10M(:)*FR(:,N)
      if(associated(MOU10M))MOU10M = MOU10M + U10M(:)*FR(:,N)
      if(associated(MOV10M))MOV10M = MOV10M + V10M(:)*FR(:,N)
      if(associated(MOT2M))MOT2M = MOT2M + T2M(:)*FR(:,N)
      if(associated(MOQ2M))MOQ2M = MOQ2M + Q2M(:)*FR(:,N)
      if(associated(MOU2M))MOU2M = MOU2M + U2M(:)*FR(:,N)
      if(associated(MOV2M))MOV2M = MOV2M + V2M(:)*FR(:,N)

    endif
    call MAPL_TimerOff(MAPL,"-SURF")

!  Aggregate to tile
!-------------------

      CHX     = CHX + CH(:,N)*FR(:,N)
      CQX     = CQX + CQ(:,N)*FR(:,N)

      if(associated(CMT)) CMT     = CMT + CM(:,N)        *FR(:,N)
      if(associated(CNT)) CNT     = CNT + CN(:  )        *FR(:,N)
      if(associated(RIT)) RIT     = RIT + RIB(:  )       *FR(:,N)
      if(associated( TH)) TH      = TH  + CH(:,N)*TC(:,N)*FR(:,N)
      if(associated( QH)) QH      = QH  + CQ(:,N)*QC(:,N)*FR(:,N)
      if(associated(Z0H)) Z0H     = Z0H + ZT             *FR(:,N)
      if(associated(GST)) GST     = GST + WW(:,N)        *FR(:,N)
      if(associated(VNT)) VNT     = VNT + UUU            *FR(:,N)

      WW(:,N) = max(CH(:,N)*(TC(:,N)-TA-(MAPL_GRAV/MAPL_CP)*DZE)/TA + MAPL_VIREPS*CQ(:,N)*(QC(:,N)-QA),0.0)
      WW(:,N) = (HPBL*MAPL_GRAV*WW(:,N))**(2./3.)

   end do SUBTILES

   if(associated( TH)) TH  = TH /CHX
   if(associated( QH)) QH  = QH /CQX
   if(associated(CHT)) CHT = CHX
   if(associated(CQT)) CQT = CQX
   if(associated(GST)) GST = sqrt(max(GST+UCN,0.0))
   if(associated(ITYO)) ITYO = real(VEG1)   ! gkw: primary type exported... where it is used?

   deallocate ( lai1 )
   deallocate ( lai2 )
   deallocate ( wght )

   deallocate(TVA)
   deallocate(TVS)
   deallocate(URA)
   deallocate(UUU)
   deallocate(ZVG)
   deallocate(DZE)
   deallocate(Z0T)
   deallocate(D0T)
   deallocate(CHX)
   deallocate(CQX)
   deallocate(VEG1)
   deallocate(VEG2)
   deallocate(FVG1)
   deallocate(FVG2)
   deallocate(RE )
   deallocate(CN )
   deallocate(ZT )
   deallocate(ZQ )
   deallocate(UCN)
   deallocate(U50M )
   deallocate(V50M )
   deallocate(T10M )
   deallocate(Q10M )
   deallocate(U10M )
   deallocate(V10M )
   deallocate(T2M )
   deallocate(Q2M )
   deallocate(U2M )
   deallocate(V2M )
   deallocate(RHOH)
   deallocate(VKH)
   deallocate(VKM)
   deallocate(USTAR)
   deallocate(XX)
   deallocate(YY)
   deallocate(CU)
   deallocate(CT)
   deallocate(RIB)
   deallocate(ZETA)
   deallocate(WS)
   deallocate(IWATER)
   deallocate(PSMB)
   deallocate(PSL)

!  All done
! ------------------------------------------------------------------------------

    call MAPL_TimerOff ( MAPL, "RUN1"  )
    call MAPL_TimerOff ( MAPL, "TOTAL" )

    RETURN_(ESMF_SUCCESS)

end subroutine RUN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

subroutine RUN2 ( GC, IMPORT, EXPORT, CLOCK, RC )

! ------------------------------------------------------------------------------
! !ARGUMENTS:
! ------------------------------------------------------------------------------

    type(ESMF_GridComp),intent(inout) :: GC
    type(ESMF_State),   intent(inout) :: IMPORT
    type(ESMF_State),   intent(inout) :: EXPORT
    type(ESMF_Clock),   intent(inout) :: CLOCK
    integer,optional,   intent(out  ) :: RC

! ------------------------------------------------------------------------------
! ErrLog Variables
! ------------------------------------------------------------------------------

    character(len=ESMF_MAXSTR) :: Iam="RUN2"
    integer :: STATUS
    character(len=ESMF_MAXSTR) :: COMP_NAME

! ------------------------------------------------------------------------------
! Local derived type aliases
! ------------------------------------------------------------------------------

    type(MAPL_MetaComp),pointer :: MAPL
    type(ESMF_Alarm)                :: ALARM

    integer :: IM,JM
    real    :: SURFLAY              ! Default (Ganymed-3 and earlier) SURFLAY=20.0 for Old Soil Params
                                    !         (Ganymed-4 and later  ) SURFLAY=50.0 for New Soil Params
    integer :: CHOOSEMOSFC
    integer :: incl_Louis_extra_derivs

    real    :: SCALE4Z0

! ------------------------------------------------------------------------------
! Begin: Get the target components name and
! set-up traceback handle.
! ------------------------------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL, RUNALARM=ALARM, RC=STATUS)
    VERIFY_(STATUS)

! Get parameters (0:Louis, 1:Monin-Obukhov)
! -----------------------------------------
    call MAPL_GetResource ( MAPL, CHOOSEMOSFC, Label="CHOOSEMOSFC:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, incl_Louis_extra_derivs, Label="INCL_LOUIS_EXTRA_DERIVS:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, SURFLAY, Label="SURFLAY:", DEFAULT=50.0, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, SCALE4Z0, Label="SCALE4Z0:", DEFAULT=0.5, RC=STATUS)
    VERIFY_(STATUS)

! ------------------------------------------------------------------------------
! If its time, recalculate the LSM tile routine
! ------------------------------------------------------------------------------

    call MAPL_TimerOn ( MAPL,"TOTAL" )
    call MAPL_TimerOn ( MAPL,"RUN2"  )

    if(ESMF_AlarmIsRinging(ALARM, RC=STATUS))then
       call ESMF_AlarmRingerOff(ALARM, RC=STATUS)
       VERIFY_(STATUS)
       call Driver ( RC=STATUS )
       VERIFY_(STATUS)
    endif

    call MAPL_TimerOff ( MAPL, "RUN2"  )
    call MAPL_TimerOff ( MAPL, "TOTAL" )

    RETURN_(ESMF_SUCCESS)

    contains

! ------------------------------------------------------------------------------
! ------------------------------------------------------------------------------

      subroutine Driver ( RC )
        integer,optional,intent(OUT) :: RC

        character(len=ESMF_MAXSTR) :: IAm
        integer :: STATUS

        ! --------------------------------------------------------------------------
        ! Local derived type aliases
        ! --------------------------------------------------------------------------

        type(ESMF_STATE) :: INTERNAL

        ! -----------------------------------------------------
        ! IMPORT Pointers
        ! -----------------------------------------------------

        real, dimension(:),   pointer :: PS
        real, dimension(:),   pointer :: TA
        real, dimension(:),   pointer :: QA
        real, dimension(:),   pointer :: UU
        real, dimension(:),   pointer :: DZ
        real, dimension(:),   pointer :: PCU
        real, dimension(:),   pointer :: PLS
        real, dimension(:),   pointer :: SNO
        real, dimension(:),   pointer :: ICE
        real, dimension(:),   pointer :: FRZR
        real, dimension(:),   pointer :: THATM
        real, dimension(:),   pointer :: QHATM
        real, dimension(:),   pointer :: CTATM
        real, dimension(:),   pointer :: CQATM

        real, dimension(:),   pointer :: drpar
        real, dimension(:),   pointer :: dfpar
        real, dimension(:),   pointer :: drnir
        real, dimension(:),   pointer :: dfnir
        real, dimension(:),   pointer :: druvr
        real, dimension(:),   pointer :: dfuvr
        real, dimension(:),   pointer :: lwdnsrf
        real, dimension(:),   pointer :: alw
        real, dimension(:),   pointer :: blw
        real, dimension(:),   pointer :: CO2SC

        real, dimension(:),   pointer :: evap
        real, dimension(:),   pointer :: devap
        real, dimension(:),   pointer :: sh
        real, dimension(:),   pointer :: dsh

        real, dimension(:),   pointer :: ROOTL
        real, dimension(:),   pointer :: Z2CH
        real, dimension(:),   pointer :: LAI
        real, dimension(:),   pointer :: GRN
        real, dimension(:),   pointer :: ASCATZ0
        real, dimension(:),   pointer :: NDVI

        real, dimension(:,:), pointer :: DUDP
        real, dimension(:,:), pointer :: DUSV
        real, dimension(:,:), pointer :: DUWT
        real, dimension(:,:), pointer :: DUSD
        real, dimension(:,:), pointer :: BCDP
        real, dimension(:,:), pointer :: BCSV
        real, dimension(:,:), pointer :: BCWT
        real, dimension(:,:), pointer :: BCSD
        real, dimension(:,:), pointer :: OCDP
        real, dimension(:,:), pointer :: OCSV
        real, dimension(:,:), pointer :: OCWT
        real, dimension(:,:), pointer :: OCSD
        real, dimension(:,:), pointer :: SUDP
        real, dimension(:,:), pointer :: SUSV
        real, dimension(:,:), pointer :: SUWT
        real, dimension(:,:), pointer :: SUSD
        real, dimension(:,:), pointer :: SSDP
        real, dimension(:,:), pointer :: SSSV
        real, dimension(:,:), pointer :: SSWT
        real, dimension(:,:), pointer :: SSSD

        ! -----------------------------------------------------
        ! INTERNAL Pointers
        ! -----------------------------------------------------

        real, dimension(:),   pointer :: bf1
        real, dimension(:),   pointer :: bf2
        real, dimension(:),   pointer :: bf3
        real, dimension(:),   pointer :: vgwmax
        real, dimension(:),   pointer :: cdcr1
        real, dimension(:),   pointer :: cdcr2
        real, dimension(:),   pointer :: psis
        real, dimension(:),   pointer :: bee
        real, dimension(:),   pointer :: poros
        real, dimension(:),   pointer :: wpwet
        real, dimension(:),   pointer :: cond
        real, dimension(:),   pointer :: gnu
        real, dimension(:),   pointer :: ars1
        real, dimension(:),   pointer :: ars2
        real, dimension(:),   pointer :: ars3
        real, dimension(:),   pointer :: ara1
        real, dimension(:),   pointer :: ara2
        real, dimension(:),   pointer :: ara3
        real, dimension(:),   pointer :: ara4
        real, dimension(:),   pointer :: arw1
        real, dimension(:),   pointer :: arw2
        real, dimension(:),   pointer :: arw3
        real, dimension(:),   pointer :: arw4
        real, dimension(:),   pointer :: tsa1
        real, dimension(:),   pointer :: tsa2
        real, dimension(:),   pointer :: tsb1
        real, dimension(:),   pointer :: tsb2
        real, dimension(:),   pointer :: atau
        real, dimension(:),   pointer :: btau
        real, dimension(:,:), pointer :: ity
        real, dimension(:,:), pointer :: fvg
        real, dimension(:),   pointer :: capac
        real, dimension(:),   pointer :: catdef
        real, dimension(:),   pointer :: rzexc
        real, dimension(:),   pointer :: srfexc
        real, dimension(:),   pointer :: ghtcnt1
        real, dimension(:),   pointer :: ghtcnt2
        real, dimension(:),   pointer :: ghtcnt3
        real, dimension(:),   pointer :: ghtcnt4
        real, dimension(:),   pointer :: ghtcnt5
        real, dimension(:),   pointer :: ghtcnt6
        real, dimension(:),   pointer :: tsurf
        real, dimension(:),   pointer :: wesnn1
        real, dimension(:),   pointer :: wesnn2
        real, dimension(:),   pointer :: wesnn3
        real, dimension(:),   pointer :: htsnnn1
        real, dimension(:),   pointer :: htsnnn2
        real, dimension(:),   pointer :: htsnnn3
        real, dimension(:),   pointer :: sndzn1
        real, dimension(:),   pointer :: sndzn2
        real, dimension(:),   pointer :: sndzn3
        real, dimension(:,:), pointer :: tc
        real, dimension(:,:), pointer :: tg
        real, dimension(:,:), pointer :: qc
        real, dimension(:,:), pointer :: ch
        real, dimension(:,:), pointer :: cm
        real, dimension(:,:), pointer :: cq
        real, dimension(:,:), pointer :: fr
        real, dimension(:,:), pointer :: dcq
        real, dimension(:,:), pointer :: dch
        real, dimension(:),   pointer :: tile_id
        real, dimension(:),   pointer :: ndep
        real, dimension(:),   pointer :: cli_t2m
        real, dimension(:),   pointer :: bgalbvr
        real, dimension(:),   pointer :: bgalbvf
        real, dimension(:),   pointer :: bgalbnr
        real, dimension(:),   pointer :: bgalbnf
        real, dimension(:,:), pointer :: cncol
        real, dimension(:,:), pointer :: cnpft
        real, dimension(:,:), pointer :: tgwm
        real, dimension(:,:), pointer :: rzmm
        real, dimension(:),   pointer :: sfmcm
        real, dimension(:),   pointer :: bflowm
        real, dimension(:),   pointer :: totwatm
        real, dimension(:),   pointer :: tairm
        real, dimension(:),   pointer :: tpm
        real, dimension(:),   pointer :: cnsum
        real, dimension(:,:,:), pointer :: psnsunm
        real, dimension(:,:,:), pointer :: psnsham
        real, dimension(:),   pointer :: sndzm
        real, dimension(:),   pointer :: asnowm
        real, dimension(:,:), pointer :: RDU001
        real, dimension(:,:), pointer :: RDU002
        real, dimension(:,:), pointer :: RDU003
        real, dimension(:,:), pointer :: RDU004
        real, dimension(:,:), pointer :: RDU005
        real, dimension(:,:), pointer :: RBC001
        real, dimension(:,:), pointer :: RBC002
        real, dimension(:,:), pointer :: ROC001
        real, dimension(:,:), pointer :: ROC002

        ! -----------------------------------------------------
        ! EXPORT Pointers
        ! -----------------------------------------------------

        real, dimension(:),   pointer :: evapout
        real, dimension(:),   pointer :: sublim
        real, dimension(:),   pointer :: shout
        real, dimension(:),   pointer :: runoff
        real, dimension(:),   pointer :: evpint
        real, dimension(:),   pointer :: evpsoi
        real, dimension(:),   pointer :: evpveg
        real, dimension(:),   pointer :: evpice
        real, dimension(:),   pointer :: evpsno
        real, dimension(:),   pointer :: bflow
        real, dimension(:),   pointer :: runsurf
        real, dimension(:),   pointer :: smelt
        real, dimension(:),   pointer :: accum
        real, dimension(:),   pointer :: hlwup
        real, dimension(:),   pointer :: swndsrf
        real, dimension(:),   pointer :: lwndsrf
        real, dimension(:),   pointer :: hlatn
        real, dimension(:),   pointer :: qinfil
        real, dimension(:),   pointer :: ar1
        real, dimension(:),   pointer :: ar2
        real, dimension(:),   pointer :: rzeq
        real, dimension(:),   pointer :: ghflx
        real, dimension(:),   pointer :: tpsurf
        real, dimension(:),   pointer :: tpsn1
        real, dimension(:),   pointer :: tpust
        real, dimension(:),   pointer :: tpsat
        real, dimension(:),   pointer :: tpwlt
        real, dimension(:),   pointer :: asnow
        real, dimension(:),   pointer :: frsat
        real, dimension(:),   pointer :: frust
        real, dimension(:),   pointer :: frwlt
        real, dimension(:),   pointer :: tp1
        real, dimension(:),   pointer :: tp2
        real, dimension(:),   pointer :: tp3
        real, dimension(:),   pointer :: tp4
        real, dimension(:),   pointer :: tp5
        real, dimension(:),   pointer :: tp6
        real, dimension(:),   pointer :: emis
        real, dimension(:),   pointer :: albvr
        real, dimension(:),   pointer :: albvf
        real, dimension(:),   pointer :: albnr
        real, dimension(:),   pointer :: albnf
        real, dimension(:),   pointer :: delts
        real, dimension(:),   pointer :: delqs
        real, dimension(:),   pointer :: delevap
        real, dimension(:),   pointer :: delsh
        real, dimension(:),   pointer :: tst
        real, dimension(:),   pointer :: lst
        real, dimension(:),   pointer :: qst

        real, dimension(:),   pointer :: WET1
        real, dimension(:),   pointer :: WET2
        real, dimension(:),   pointer :: WET3
        real, dimension(:),   pointer :: WCSF
        real, dimension(:),   pointer :: WCRZ
        real, dimension(:),   pointer :: WCPR
        real, dimension(:),   pointer :: SNOMAS
        real, dimension(:),   pointer :: SNOWDP

        real, dimension(:),   pointer :: EVLAND
        real, dimension(:),   pointer :: PRLAND
        real, dimension(:),   pointer :: SNOLAND
        real, dimension(:),   pointer :: DRPARLAND
        real, dimension(:),   pointer :: DFPARLAND
        real, dimension(:),   pointer :: LHSNOW
        real, dimension(:),   pointer :: SWNETSNOW1
        real, dimension(:),   pointer :: LWUPSNOW
        real, dimension(:),   pointer :: LWDNSNOW
        real, dimension(:),   pointer :: TCSORIG
        real, dimension(:),   pointer :: TPSN1IN
        real, dimension(:),   pointer :: TPSN1OUT
        real, dimension(:),   pointer :: GHSNOW
        real, dimension(:),   pointer :: LHLAND
        real, dimension(:),   pointer :: SHLAND
        real, dimension(:),   pointer :: SWLAND
        real, dimension(:),   pointer :: SWDOWNLAND
        real, dimension(:),   pointer :: LWLAND
        real, dimension(:),   pointer :: GHLAND
        real, dimension(:),   pointer :: GHTSKIN
        real, dimension(:),   pointer :: SMLAND
        real, dimension(:),   pointer :: TWLAND
        real, dimension(:),   pointer :: TELAND
        real, dimension(:),   pointer :: TSLAND
        real, dimension(:),   pointer :: DWLAND
        real, dimension(:),   pointer :: DHLAND
        real, dimension(:),   pointer :: SPLAND
        real, dimension(:),   pointer :: SPWATR
        real, dimension(:),   pointer :: SPSNOW

        real, dimension(:),   pointer :: CNLAI
        real, dimension(:),   pointer :: CNTLAI
        real, dimension(:),   pointer :: CNSAI
        real, dimension(:),   pointer :: CNTOTC
        real, dimension(:),   pointer :: CNVEGC
        real, dimension(:),   pointer :: CNROOT
        real, dimension(:),   pointer :: CNNPP
        real, dimension(:),   pointer :: CNGPP
        real, dimension(:),   pointer :: CNSR
        real, dimension(:),   pointer :: CNNEE
        real, dimension(:),   pointer :: CNXSMR
        real, dimension(:),   pointer :: CNADD
        real, dimension(:),   pointer :: CNLOSS
        real, dimension(:),   pointer :: CNBURN
        real, dimension(:),   pointer :: PARABS
        real, dimension(:),   pointer :: PARINC
        real, dimension(:),   pointer :: SCSAT
        real, dimension(:),   pointer :: SCUNS
        real, dimension(:),   pointer :: BTRANT
        real, dimension(:),   pointer :: SIF
        real, dimension(:),   pointer :: CNFSEL

        real, dimension(:),   pointer :: WAT10CM
        real, dimension(:),   pointer :: WATSOI
        real, dimension(:),   pointer :: ICESOI
        real, dimension(:),   pointer :: SHSNOW
        real, dimension(:),   pointer :: AVETSNOW
        real, pointer, dimension(:)   :: RMELTDU001
        real, pointer, dimension(:)   :: RMELTDU002
        real, pointer, dimension(:)   :: RMELTDU003
        real, pointer, dimension(:)   :: RMELTDU004
        real, pointer, dimension(:)   :: RMELTDU005
        real, pointer, dimension(:)   :: RMELTBC001
        real, pointer, dimension(:)   :: RMELTBC002
        real, pointer, dimension(:)   :: RMELTOC001
        real, pointer, dimension(:)   :: RMELTOC002

        ! --------------------------------------------------------------------------
        ! Local pointers for tile variables
        ! --------------------------------------------------------------------------

        INTEGER,pointer,dimension(:) :: CAT_ID
        real,pointer,dimension(:) :: dzsf
        real,pointer,dimension(:) :: swnetfree
        real,pointer,dimension(:) :: swnetsnow
        real,pointer,dimension(:) :: qa1
        real,pointer,dimension(:) :: qa2
        real,pointer,dimension(:) :: qa4
        real,pointer,dimension(:) :: tilezero
        real,pointer,dimension(:) :: zth
        real,pointer,dimension(:) :: lats
        real,pointer,dimension(:) :: lons
        real,pointer,dimension(:) :: slr
        real,pointer,dimension(:) :: rdc
	real,pointer,dimension(:) :: PRECU
	real,pointer,dimension(:) :: PRELS
	real,pointer,dimension(:) :: SNOW
	real,pointer,dimension(:) :: UUU, RHO
	real,pointer,dimension(:) :: LAI0,GRN0,ZVG
	real,pointer,dimension(:) :: Z0, D0
	real,pointer,dimension(:) :: sfmc, rzmc, prmc, entot, wtot
	real,pointer,dimension(:) :: ghflxsno, ghflxtskin
        real,pointer,dimension(:) :: SHSNOW1, AVETSNOW1, WAT10CM1, WATSOI1, ICESOI1
        real,pointer,dimension(:) :: LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW
        real,pointer,dimension(:) :: TCSORIG1, TPSN1IN1, TPSN1OUT1
	real,pointer,dimension(:) :: WCHANGE, ECHANGE, HSNACC, EVACC, SHACC
	real,pointer,dimension(:) :: SNOVR, SNOVF, SNONR, SNONF
	real,pointer,dimension(:) :: VSUVR, VSUVF
	real,pointer,dimension(:) :: ALWX, BLWX
        real,pointer,dimension(:) :: LHACC, SUMEV
	real,pointer,dimension(:) :: fveg1, fveg2
        real,pointer,dimension(:) :: FICE1
        real,pointer,dimension(:) :: SLDTOT

!       real*8,pointer,dimension(:) :: fsum

        real,pointer,dimension(:,:) :: ghtcnt
        real,pointer,dimension(:,:) :: wesnn
        real,pointer,dimension(:,:) :: htsnnn
        real,pointer,dimension(:,:) :: sndzn
        real,pointer,dimension(:,:) :: shsbt
        real,pointer,dimension(:,:) :: dshsbt
        real,pointer,dimension(:,:) :: evsbt
        real,pointer,dimension(:,:) :: devsbt
        real,pointer,dimension(:,:) :: DEDTC 
        real,pointer,dimension(:,:) :: DHSDQA
        real,pointer,dimension(:,:) :: CFT
        real,pointer,dimension(:,:) :: RA
        real,pointer,dimension(:,:) :: CFQ
        real,pointer,dimension(:,:) :: TCO
        real,pointer,dimension(:,:) :: QCO
        real,pointer,dimension(:,:) :: DQS
        real,pointer,dimension(:,:) :: QSAT

        integer,dimension(:),pointer :: veg1
        integer,dimension(:),pointer :: veg2

        real,pointer,dimension(:) :: RCSAT 
	real,pointer,dimension(:) :: DRCSDT
	real,pointer,dimension(:) :: DRCSDQ
	real,pointer,dimension(:) :: RCUNS 
	real,pointer,dimension(:) :: DRCUDT
	real,pointer,dimension(:) :: DRCUDQ

        real,pointer,dimension(:,:,:) :: RCONSTIT
        real,pointer,dimension(:,:)   :: TOTDEPOS
        real,pointer,dimension(:,:)   :: RMELT

        ! --------------------------------------------------------------------------
        ! Locals for parameter lookup
        ! --------------------------------------------------------------------------

        ! vegetation calculations

        real,dimension(NTYPS) :: VGRF11
        real,dimension(NTYPS) :: VGRF12
        real,dimension(NTYPS) :: VGTR11
        real,dimension(NTYPS) :: VGTR12
        real,dimension(NTYPS) :: VGROCA
        real,dimension(NTYPS) :: VGROTD
        real,dimension(NTYPS) :: VGRDRS
        real,dimension(NTYPS) :: VGDDA, VGDDB, VGDDC
        real,dimension(NTYPS) :: VGRDA, VGRDB

        real,dimension(:),allocatable :: RSL1, RSL2
        real,dimension(:),allocatable :: SQSCAT

        ! albedo calculation stuff

        type(ESMF_Config)           :: CF
        type(MAPL_SunOrbit)         :: ORBIT
        type(ESMF_Time)             :: CURRENT_TIME, StopTime, NextTime
        type(ESMF_Time)             :: BEFORE
        type(ESMF_Time)             :: NOW
        type(ESMF_Time)             :: MODELSTART
        type(ESMF_Time)             :: AFTER
        type(ESMF_TimeInterval)     :: DELT
        type(ESMF_TimeInterval)     :: TINT
        real                        :: DT_SOLAR
        type(ESMF_Alarm)            :: SOLALARM
        logical                     :: solalarmison
        logical                     :: debugzth
        real                        :: FAC

        real,parameter              :: PRECIPFRAC=1.0
        real                        :: DT
        integer                     :: NTILES
        integer                     :: I, N
        integer                     :: AEROSOL_DEPOSITION
        integer                     :: N_CONST_LAND4SNWALB
        integer                     :: DO_GOSWIM

	! dummy variables for call to get snow temp

        real    :: FICE
        logical :: DUMFLAG1,DUMFLAG2
        integer                         :: nmax
        type(ESMF_VM)                   :: VM

   ! variables from sun orbit
        integer :: years_per_cycle, days_per_cycle
        real    :: year_length
        real, pointer, dimension(:) :: zco => null()
        real, pointer, dimension(:) :: zso => null()

#ifdef DBG_CATCH_INPUTS
        ! vars for debugging purposes
        type(ESMF_Grid)                 :: TILEGRID
        type (MAPL_LocStream)           :: LOCSTREAM
        integer, pointer                :: mask(:)
        integer                         :: nt
        integer, save                   :: unit_i=0
        logical, save                   :: firsttime=.true.
        integer                         :: unit
	integer 			:: NT_GLOBAL

#endif

       ! Offline case

        type(OFFLINE_WRAP)          :: wrap
        logical                     :: is_OFFLINE
        real,dimension(:,:),allocatable :: ALWN, BLWN
        ! un-adelterated TC's and QC's
        real, pointer               :: TC1_0(:), TC2_0(:),  TC4_0(:)
        real, pointer               :: QA1_0(:), QA2_0(:),  QA4_0(:)

        ! --------------------------------------------------------------------------
        ! Lookup tables
        ! --------------------------------------------------------------------------

        data VGRF11 / 0.100, 0.100, 0.070, 0.105, 0.100, 0.100 /
        data VGRF12 / 0.160, 0.160, 0.160, 0.360, 0.160, 0.160 /
        data VGTR11 / 0.050, 0.050, 0.050, 0.070, 0.050, 0.050 /
        data VGTR12 / 0.001, 0.001, 0.001, 0.220, 0.001, 0.001 /
        data VGROTD / 1.000, 1.000, 0.500, 0.500, 0.500, 0.200 / 

        data VGROCA / 0.384E-6, 0.384E-6, 0.384E-6, 0.384E-6, 0.384E-6, 0.384E-6/
        data VGRDRS / 0.750E13, 0.750E13, 0.750E13, 0.400E13, 0.750E13, 0.750E13/

! Correction to RDC formulation -Randy Koster, 4/1/2011
!        data VGRDA / 285.9, 294.9, 652.9,  25.8,  100.7,  22.9,  23.8, 23.8/
!        data VGRDB / 5.1 ,  7.2, 10.8,  4.8,  1.8,  5.1,  .000, .000/

        data VGRDA / 285.9, 355.18, 660.24,  30.06,  100.7,  24.36/
        data VGRDB / 5.1 ,  7.2, 10.5,  4.8,  1.8,  5.1/

! gkw: following is for CN model
! ------------------------------
    integer, parameter :: nveg  = num_veg ! number of vegetation types
    integer, parameter :: nzone = num_zon ! number of stress zones

    real, allocatable, dimension(:) :: btran, wgt, bt1, bt2, bt4, wpp, fwet
    real, allocatable, dimension(:) :: sm1, sm2, sm4
    real, allocatable, dimension(:) :: btran1, btran2, btran3, tcx, qax
    real, allocatable, dimension(:) :: rcx, rcxdt, rcxdq, tx1, tx2, tx3, qx1, qx2, qx3
    real, allocatable, dimension(:,:) :: tgw, rzm, rc00, rcdt,rcdq, totcolc, wtzone
    real, allocatable, dimension(:,:,:) :: elai,esai,fveg,tlai,psnsun,psnsha,laisun,laisha
    integer, allocatable, dimension(:,:,:) :: ityp
    real, allocatable, dimension(:) :: car1, car2, car4
    real, allocatable, dimension(:) :: parzone, para
    real, allocatable, dimension(:) :: npp, gpp, sr, nee, padd, root, vegc, xsmr
    real, allocatable, dimension(:) :: burn, fsel, closs
    real, allocatable, dimension(:) :: dayl, dayl_fac, CO2V

    real, parameter :: dtc = 0.03 ! canopy temperature perturbation (K) [approx 1:10000]
    real, parameter :: dea = 0.10 ! vapor pressure perturbation (Pa) [approx 1:10000]

    real, allocatable, dimension(:) :: totwat ! total soil liquid water (kg/m2)
    real, save :: ashift = 0. ! for baseflow. gkw: this should match value in routine "base" in catchment

    integer :: nz, iv
    real :: cn1, cn2, cn3, cn12, cn23, ar, f1, f2, f3, f4, ax1, ax2, ax4

    real, allocatable, dimension(:,:) :: psnsunx,psnshax,sifsunx,sifshax,laisunx,laishax,elaz,esaz,fvez
    integer, allocatable, dimension(:,:) :: ityz

    real, allocatable, dimension(:) :: ht, tp, soilice
    real :: zbar, frice

    real, allocatable, dimension(:,:)  :: col
    real, allocatable, dimension(:,:,:) :: pft

    integer :: ntile, nv, dpy, ierr, iok, ndt
    integer, save :: year_prev = -9999
    real :: dtcn ! carbon model time step

    integer :: AGCM_YY, AGCM_MM, AGCM_DD, AGCM_S, dofyr
    logical, save :: first = .true.
    integer*8, save :: istep = 1 ! gkw: legacy variable from offline

! solar declination related
    real :: ob, declin, zs, zc, max_decl, max_dayl
    integer :: year, iday, idayp1

    real :: co2
    real, external :: getco2

! temporaries for call to SIBALB for each type
! --------------------------------------------
    real, allocatable, dimension(:) :: lai1, lai2, wght
    real, allocatable, dimension(:) :: ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp
    real, allocatable, dimension(:) :: SNOVR_tmp, SNONR_tmp, SNOVF_tmp, SNONF_tmp

        IAm=trim(COMP_NAME)//"::RUN2::Driver"

        ! Begin

        IAm=trim(COMP_NAME)//"Driver"

        ! --------------------------------------------------------------------------
        ! Get time step from configuration
        ! --------------------------------------------------------------------------

        call ESMF_GridCompGet  ( GC, CONFIG=CF, RC=STATUS )
        VERIFY_(STATUS)

        ! --------------------------------------------------------------------------
        ! Get my internal MAPL_Generic state
        ! --------------------------------------------------------------------------

        call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
        VERIFY_(STATUS)

        call MAPL_Get(MAPL, HEARTBEAT = DT, RC=STATUS)
        VERIFY_(STATUS)

        call ESMF_ConfigGetAttribute ( CF, DT                  ,&
             Label   = trim(COMP_NAME)//"_DT:"     ,&
             Default = DT                          ,&
             RC=STATUS )
        VERIFY_(STATUS)

        call MAPL_GetResource ( MAPL, DO_GOSWIM, Label="N_CONST_LAND4SNWALB:", DEFAULT=0, RC=STATUS)
        VERIFY_(STATUS)

        ! Get component's private internal state
        call ESMF_UserCompGetInternalState(gc, 'OfflineMode', wrap, status)
        VERIFY_(status)

        call ESMF_VMGetCurrent ( VM, RC=STATUS )
        ! Component's offline mode
        is_OFFLINE = wrap%ptr%CATCH_OFFLINE
        ! if (MAPL_AM_I_Root(VM)) print *, trim(Iam)//'::OFFLINE mode: ', is_OFFLINE

        ! --------------------------------------------------------------------------
        ! Get parameters from generic state.
        ! --------------------------------------------------------------------------

        call MAPL_Get ( MAPL                 ,&
             RUNALARM  = ALARM                            ,&
             ORBIT     = ORBIT                            ,&
             TILELATS  = LATS                             ,&
             TILELONS  = LONS                             ,&
             INTERNAL_ESMF_STATE = INTERNAL               ,&
             RC=STATUS )
        VERIFY_(STATUS)

        ! Get parameters to zero the deposition rate 
        ! 0: Use all GOCART aerosol values, 1: turn OFF everythying, 
	! 2: turn off dust ONLY,3: turn off Black Carbon ONLY,4: turn off Organic Carbon ONLY
        ! __________________________________________

        call MAPL_GetResource ( MAPL, AEROSOL_DEPOSITION, Label="AEROSOL_DEPOSITION:", DEFAULT=0, RC=STATUS)
        VERIFY_(STATUS)

        ! GOSWIM ANOW_ALBEDO 
        ! 0 : GOSWIM snow albedo scheme is turned off
        ! 9 : i.e. N_CONSTIT in Stieglitz to turn on GOSWIM snow albedo scheme
 
        call MAPL_GetResource ( MAPL, N_CONST_LAND4SNWALB, Label="N_CONST_LAND4SNWALB:", DEFAULT=0, RC=STATUS)
        VERIFY_(STATUS)

        ! -----------------------------------------------------
        ! IMPORT Pointers
        ! -----------------------------------------------------

        call MAPL_GetPointer(IMPORT,PS     ,'PS'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,TA     ,'TA'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,QA     ,'QA'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,UU     ,'UU'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DZ     ,'DZ'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,PCU    ,'PCU'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,PLS    ,'PLS'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SNO    ,'SNO'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,ICE    ,'ICE'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,FRZR   ,'FRZR'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DRPAR  ,'DRPAR'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DFPAR  ,'DFPAR'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DRNIR  ,'DRNIR'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DFNIR  ,'DFNIR'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DRUVR  ,'DRUVR'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DFUVR  ,'DFUVR'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,LWDNSRF,'LWDNSRF',RC=STATUS); VERIFY_(STATUS)

        call MAPL_GetPointer(IMPORT,ALW    ,'ALW'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,BLW    ,'BLW'    ,RC=STATUS); VERIFY_(STATUS)

        call MAPL_GetPointer(IMPORT,EVAP   ,'EVAP'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DEVAP  ,'DEVAP'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SH     ,'SH'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DSH    ,'DSH'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,THATM  ,'THATM'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,QHATM  ,'QHATM'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,CTATM  ,'CTATM'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,CQATM  ,'CQATM'  ,RC=STATUS); VERIFY_(STATUS)        
        IF (DO_CO2SC /= 0) call MAPL_GetPointer(IMPORT,CO2SC  ,'CO2SC'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,LAI    ,'LAI'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,GRN    ,'GRN'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,ROOTL  ,'ROOTL'  ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,Z2CH   ,'Z2CH'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,ASCATZ0,'ASCATZ0',RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,NDVI   ,'NDVI'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DUDP   ,'DUDP'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DUSV   ,'DUSV'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DUWT   ,'DUWT'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,DUSD   ,'DUSD'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,BCDP   ,'BCDP'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,BCSV   ,'BCSV'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,BCWT   ,'BCWT'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,BCSD   ,'BCSD'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,OCDP   ,'OCDP'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,OCSV   ,'OCSV'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,OCWT   ,'OCWT'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,OCSD   ,'OCSD'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SUDP   ,'SUDP'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SUSV   ,'SUSV'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SUWT   ,'SUWT'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SUSD   ,'SUSD'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SSDP   ,'SSDP'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SSSV   ,'SSSV'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SSWT   ,'SSWT'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT,SSSD   ,'SSSD'   ,RC=STATUS); VERIFY_(STATUS)

        ! -----------------------------------------------------
        ! INTERNAL Pointers
        ! -----------------------------------------------------

        call MAPL_GetPointer(INTERNAL,BF1        ,'BF1'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BF2        ,'BF2'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BF3        ,'BF3'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,VGWMAX     ,'VGWMAX'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CDCR1      ,'CDCR1'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CDCR2      ,'CDCR2'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,PSIS       ,'PSIS'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BEE        ,'BEE'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,POROS      ,'POROS'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,WPWET      ,'WPWET'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,COND       ,'COND'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GNU        ,'GNU'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARS1       ,'ARS1'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARS2       ,'ARS2'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARS3       ,'ARS3'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARA1       ,'ARA1'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARA2       ,'ARA2'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARA3       ,'ARA3'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARA4       ,'ARA4'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARW1       ,'ARW1'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARW2       ,'ARW2'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARW3       ,'ARW3'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ARW4       ,'ARW4'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TSA1       ,'TSA1'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TSA2       ,'TSA2'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TSB1       ,'TSB1'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TSB2       ,'TSB2'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ATAU       ,'ATAU'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BTAU       ,'BTAU'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ITY        ,'ITY'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,FVG        ,'FVG'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TC         ,'TC'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,QC         ,'QC'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TG         ,'TG'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CAPAC      ,'CAPAC'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CATDEF     ,'CATDEF'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RZEXC      ,'RZEXC'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SRFEXC     ,'SRFEXC'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GHTCNT1    ,'GHTCNT1'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GHTCNT2    ,'GHTCNT2'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GHTCNT3    ,'GHTCNT3'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GHTCNT4    ,'GHTCNT4'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GHTCNT5    ,'GHTCNT5'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GHTCNT6    ,'GHTCNT6'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TSURF      ,'TSURF'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,WESNN1     ,'WESNN1'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,WESNN2     ,'WESNN2'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,WESNN3     ,'WESNN3'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,HTSNNN1    ,'HTSNNN1'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,HTSNNN2    ,'HTSNNN2'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,HTSNNN3    ,'HTSNNN3'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNDZN1     ,'SNDZN1'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNDZN2     ,'SNDZN2'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNDZN3     ,'SNDZN3'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CH         ,'CH'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CM         ,'CM'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CQ         ,'CQ'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,FR         ,'FR'         ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,DCQ        ,'DCQ'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,DCH        ,'DCH'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TILE_ID    ,'TILE_ID'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,NDEP       ,'NDEP'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CLI_T2M    ,'CLI_T2M'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBVR    ,'BGALBVR'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBVF    ,'BGALBVF'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBNR    ,'BGALBNR'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBNF    ,'BGALBNF'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CNCOL      ,'CNCOL'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CNPFT      ,'CNPFT'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TGWM       ,'TGWM'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RZMM       ,'RZMM'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SFMCM      ,'SFMCM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BFLOWM     ,'BFLOWM'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TOTWATM    ,'TOTWATM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TAIRM      ,'TAIRM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TPM        ,'TPM'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CNSUM      ,'CNSUM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,PSNSUNM    ,'PSNSUNM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,PSNSHAM    ,'PSNSHAM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNDZM      ,'SNDZM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ASNOWM     ,'ASNOWM'     ,RC=STATUS); VERIFY_(STATUS)
        if (DO_GOSWIM /= 0) then
           call MAPL_GetPointer(INTERNAL,RDU001     ,'RDU001'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,RDU002     ,'RDU002'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,RDU003     ,'RDU003'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,RDU004     ,'RDU004'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,RDU005     ,'RDU005'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,RBC001     ,'RBC001'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,RBC002     ,'RBC002'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,ROC001     ,'ROC001'     , RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,ROC002     ,'ROC002'     , RC=STATUS); VERIFY_(STATUS)
        endif

        ! -----------------------------------------------------
        ! EXPORT POINTERS
        ! -----------------------------------------------------

        call MAPL_GetPointer(EXPORT,EVAPOUT,'EVAPOUT',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SUBLIM,'SUBLIM',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SHOUT,  'SHOUT'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RUNOFF, 'RUNOFF' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPINT, 'EVPINT' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPSOI, 'EVPSOI' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPVEG, 'EVPVEG' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPICE, 'EVPICE' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WAT10CM,'WAT10CM',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WATSOI, 'WATSOI' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ICESOI, 'ICESOI' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPSNO, 'EVPSNO'              ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,BFLOW,  'BASEFLOW',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RUNSURF,'RUNSURF',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SMELT,  'SMELT'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,HLWUP,  'HLWUP'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWNDSRF,'SWNDSRF',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWNDSRF,'LWNDSRF',ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,HLATN,  'HLATN'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,QINFIL, 'QINFIL' ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,AR1,    'AR1'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,AR2,    'AR2'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RZEQ,   'RZEQ'   ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHFLX,  'GHFLX'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSURF, 'TPSURF' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSN1,  'TPSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPUST,  'TPUNST' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSAT,  'TPSAT'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPWLT,  'TPWLT'  ,             RC=STATUS); VERIFY_(STATUS)
        if (is_OFFLINE) then
           call MAPL_GetPointer(INTERNAL, ASNOW, 'ASNOW', RC=STATUS); VERIFY_(STATUS)
        else
           call MAPL_GetPointer(EXPORT,ASNOW,  'ASNOW'  ,             RC=STATUS); VERIFY_(STATUS)
        endif
        call MAPL_GetPointer(EXPORT,SHSNOW, 'SHSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,AVETSNOW,'AVETSNOW',           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FRSAT,  'FRSAT'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FRUST,  'FRUST'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FRWLT,  'FRWLT'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP1,    'TP1'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP2,    'TP2'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP3,    'TP3'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP4,    'TP4'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP5,    'TP5'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP6,    'TP6'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        if (is_OFFLINE) then
           call MAPL_GetPointer(INTERNAL, EMIS, 'EMIS', RC=STATUS);   VERIFY_(STATUS)
        else
           call MAPL_GetPointer(EXPORT,EMIS,   'EMIS'   ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        endif
        call MAPL_GetPointer(EXPORT,ALBVR,  'ALBVR'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBVF,  'ALBVF'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBNR,  'ALBNR'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBNF,  'ALBNF'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DELTS,  'DELTS'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DELQS,  'DELQS'  ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TST  ,  'TST'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,QST  ,  'QST'    ,ALLOC=.true.,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LST  ,  'LST'    ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WET1 ,  'WET1'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WET2 ,  'WET2'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WET3 ,  'WET3'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WCSF ,  'WCSF'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WCRZ ,  'WCRZ'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WCPR ,  'WCPR'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ACCUM,  'ACCUM'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SNOMAS,'SNOWMASS',             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SNOWDP, 'SNOWDP' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVLAND, 'EVLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PRLAND, 'PRLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SNOLAND, 'SNOLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DRPARLAND, 'DRPARLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DFPARLAND, 'DFPARLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LHSNOW, 'LHSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWNETSNOW1, 'SWNETSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWUPSNOW, 'LWUPSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWDNSNOW, 'LWDNSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TCSORIG, 'TCSORIG' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSN1IN, 'TPSN1IN' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSN1OUT, 'TPSN1OUT' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LHLAND, 'LHLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SHLAND, 'SHLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWLAND, 'SWLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWDOWNLAND, 'SWDOWNLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWLAND, 'LWLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHLAND, 'GHLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHSNOW, 'GHSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHTSKIN,'GHTSKIN',             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SMLAND, 'SMLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TWLAND, 'TWLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TELAND, 'TELAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TSLAND, 'TSLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DWLAND, 'DWLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DHLAND, 'DHLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPLAND, 'SPLAND' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPWATR, 'SPWATR' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPSNOW, 'SPSNOW' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNLAI,  'CNLAI'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTLAI, 'CNTLAI' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSAI,  'CNSAI'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTOTC, 'CNTOTC' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNVEGC, 'CNVEGC' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNROOT, 'CNROOT' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNPP,  'CNNPP'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNGPP,  'CNGPP'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSR,   'CNSR'   ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNEE,  'CNNEE'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNXSMR, 'CNXSMR' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNADD,  'CNADD'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNLOSS, 'CNLOSS' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNBURN, 'CNBURN' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PARABS, 'PARABS' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PARINC, 'PARINC' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SCSAT,  'SCSAT'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SCUNS,  'SCUNS'  ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,BTRANT, 'BTRANT' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SIF,    'SIF'    ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFSEL, 'CNFSEL' ,             RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU001,'RMELTDU001',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU002,'RMELTDU002',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU003,'RMELTDU003',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU004,'RMELTDU004',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU005,'RMELTDU005',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTBC001,'RMELTBC001',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTBC002,'RMELTBC002',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTOC001,'RMELTOC001',  RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTOC002,'RMELTOC002',  RC=STATUS); VERIFY_(STATUS)

        NTILES = size(PS)

    allocate(   ityp(ntiles,nveg,nzone) )
    allocate(   fveg(ntiles,nveg,nzone) )
    allocate( wtzone(ntiles,nzone) )
    allocate(   elai(ntiles,nveg,nzone) )
    allocate(   esai(ntiles,nveg,nzone) )
    allocate(   tlai(ntiles,nveg,nzone) )

! set CLM CN PFT & fraction, set carbon zone weights
! --------------------------------------------------
    do nz = 1,nzone
      ityp(:,:,nz) = nint(ity(:,:))
      fveg(:,:,nz) = fvg(:,:)
      wtzone(:,nz) = CN_zone_weight(nz)
    end do

! obtain LAI from previous time step (from CN model)
! --------------------------------------------------
    call get_CN_LAI(ntiles,nveg,nzone,ityp,fveg,elai,esai=esai)


        ! --------------------------------------------------------------------------
        ! ALLOCATE LOCAL POINTERS
        ! --------------------------------------------------------------------------

        allocate(GHTCNT (6,NTILES))
        allocate(WESNN  (3,NTILES))
        allocate(HTSNNN (3,NTILES))
        allocate(SNDZN  (3,NTILES))
        allocate(TILEZERO (NTILES))
        allocate(DZSF     (NTILES))
        allocate(SWNETFREE(NTILES))
        allocate(SWNETSNOW(NTILES))
        allocate(VEG1     (NTILES))
        allocate(VEG2     (NTILES))
        allocate(RCSAT    (NTILES))
        allocate(DRCSDT   (NTILES))
        allocate(DRCSDQ   (NTILES))
        allocate(RCUNS    (NTILES))
        allocate(DRCUDT   (NTILES))
        allocate(DRCUDQ   (NTILES))
        allocate(ZTH      (NTILES))  
        allocate(SLR      (NTILES))  
        allocate(RSL1     (NTILES)) 
        allocate(RSL2     (NTILES)) 
        allocate(SQSCAT   (NTILES))
        allocate(RDC      (NTILES))  
	allocate(UUU      (NTILES))
	allocate(RHO      (NTILES))
	allocate(ZVG      (NTILES))
	allocate(LAI0     (NTILES))
	allocate(GRN0     (NTILES))
	allocate(Z0       (NTILES))
	allocate(D0       (NTILES))
	allocate(SFMC     (NTILES))
	allocate(RZMC     (NTILES))
	allocate(PRMC     (NTILES))
	allocate(ENTOT    (NTILES))
	allocate(ghflxsno (NTILES))
	allocate(ghflxtskin(NTILES))
	allocate(WTOT     (NTILES))
	allocate(WCHANGE  (NTILES))
	allocate(ECHANGE  (NTILES))
	allocate(HSNACC   (NTILES))
	allocate(EVACC    (NTILES))
	allocate(SHACC    (NTILES))
	allocate(VSUVR    (NTILES))
	allocate(VSUVF    (NTILES))
	allocate(SNOVR    (NTILES))
	allocate(SNOVF    (NTILES))
	allocate(SNONR    (NTILES))
	allocate(SNONF    (NTILES))
	allocate(CAT_ID   (NTILES))
	allocate(ALWX     (NTILES))
	allocate(BLWX     (NTILES))
        allocate(SHSNOW1   (NTILES))
        allocate(AVETSNOW1 (NTILES))
        allocate(WAT10CM1  (NTILES))
        allocate(WATSOI1   (NTILES))
        allocate(ICESOI1   (NTILES))
        allocate(LHSNOW1   (NTILES))
        allocate(LWUPSNOW1 (NTILES))
        allocate(LWDNSNOW1 (NTILES))
        allocate(NETSWSNOW (NTILES))
        allocate(TCSORIG1  (NTILES))
        allocate(TPSN1IN1  (NTILES))
        allocate(TPSN1OUT1 (NTILES))
        allocate(LHACC     (NTILES))
        allocate(SUMEV     (NTILES))
	allocate(fveg1     (NTILES))
	allocate(fveg2     (NTILES))
        allocate(FICE1     (NTILES)) 
        allocate(SLDTOT    (NTILES))             ! total solid precip

        allocate(SHSBT    (NTILES,NUM_SUBTILES))
        allocate(DSHSBT   (NTILES,NUM_SUBTILES))
        allocate(EVSBT    (NTILES,NUM_SUBTILES))
        allocate(DEVSBT   (NTILES,NUM_SUBTILES))
        allocate(DEDTC    (NTILES,NUM_SUBTILES))
        allocate(DHSDQA   (NTILES,NUM_SUBTILES))
        allocate(CFT      (NTILES,NUM_SUBTILES))
        allocate(CFQ      (NTILES,NUM_SUBTILES))
        allocate(TCO      (NTILES,NUM_SUBTILES))
        allocate(QCO      (NTILES,NUM_SUBTILES))
        allocate(DQS      (NTILES,NUM_SUBTILES))
        allocate(QSAT     (NTILES,NUM_SUBTILES))
        allocate(RA       (NTILES,NUM_SUBTILES))
        allocate(RCONSTIT (NTILES,N_SNOW,N_constit))
        allocate(TOTDEPOS (NTILES,N_constit))
        allocate(RMELT    (NTILES,N_constit))
        allocate(ALWN     (NTILES,NUM_SUBTILES))
        allocate(BLWN     (NTILES,NUM_SUBTILES))
        allocate(TC1_0    (NTILES))
        allocate(TC2_0    (NTILES))
        allocate(TC4_0    (NTILES))
        allocate(QA1_0    (NTILES))
        allocate(QA2_0    (NTILES))
        allocate(QA4_0    (NTILES))

        call ESMF_VMGetCurrent ( VM, RC=STATUS )
        ! --------------------------------------------------------------------------
        ! Catchment Id and vegetation types used to index into tables
        ! --------------------------------------------------------------------------

        CAT_ID = nint(tile_id)

        where(ITY(:,1) > 0.)             ! gkw: account for "split" types
          VEG1 = map_cat(nint(ITY(:,1))) ! map  primary  CN PFT to catchment type
         elsewhere
          VEG1 = map_cat(nint(ITY(:,2))) ! map  primary  CN PFT to catchment type
        endwhere
        where(ITY(:,3) > 0.)
          VEG2 = map_cat(nint(ITY(:,3))) ! map secondary CN PFT to catchment type
         elsewhere
          VEG2 = map_cat(nint(ITY(:,4))) ! map secondary CN PFT to catchment type
        endwhere

        fveg1(:) = fvg(:,1) + fvg(:,2) ! sum veg fractions (primary)   gkw: NVEG specific
        fveg2(:) = fvg(:,3) + fvg(:,4) ! sum veg fractions (secondary) gkw: fveg1+fveg2=1

        allocate ( lai1(ntiles) )
        allocate ( lai2(ntiles) )
        allocate ( wght(ntiles) )

        lai1 = 0.
        wght = 0.
        do nz = 1,nzone
          do nv = 1,2
            lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
            wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
          end do
        end do
        lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

        lai2 = 0.
        wght = 0.
        do nz = 1,nzone
          do nv = 3,4
            lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
            wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
          end do
        end do
        lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

! LAI seen by the land model
! --------------------------
        lai = fveg1*lai1 + fveg2*lai2 ! gkw: prognostic LAI on catch_internal_rst (overwrite VEGDYN import)

        ! --------------------------------------------------------------------------
        ! surface layer depth for soil moisture
        ! --------------------------------------------------------------------------
        
        DZSF(    :) = SURFLAY

        ! --------------------------------------------------------------------------
        ! build arrays from internal state
        ! --------------------------------------------------------------------------

        GHTCNT(1,:) = GHTCNT1
        GHTCNT(2,:) = GHTCNT2
        GHTCNT(3,:) = GHTCNT3
        GHTCNT(4,:) = GHTCNT4
        GHTCNT(5,:) = GHTCNT5
        GHTCNT(6,:) = GHTCNT6

        WESNN (1,:) = WESNN1
        WESNN (2,:) = WESNN2
        WESNN (3,:) = WESNN3

        HTSNNN(1,:) = HTSNNN1
        HTSNNN(2,:) = HTSNNN2
        HTSNNN(3,:) = HTSNNN3

        SNDZN (1,:) = SNDZN1
        SNDZN (2,:) = SNDZN2
        SNDZN (3,:) = SNDZN3

        debugzth = .false.

        ! --------------------------------------------------------------------------
        ! Get the current time. 
        ! --------------------------------------------------------------------------

        call ESMF_ClockGet( CLOCK, currTime=CURRENT_TIME, startTime=MODELSTART, TIMESTEP=DELT,  RC=STATUS )
        VERIFY_(STATUS)
        if (MAPL_AM_I_Root(VM).and.debugzth) then
         print *,' start time of clock '
         CALL ESMF_TimePrint ( MODELSTART, OPTIONS="string", RC=STATUS )
        endif

        ! --------------------------------------------------------------------------
        ! retrieve the zenith angle
        ! --------------------------------------------------------------------------

!! The next sequence is to make sure that the albedo here and in solar are in sync
!!
! Need to know when Solar was called last, so first get the solar alarm
        call ESMF_ClockGetAlarm ( CLOCK, alarmname="SOLAR_Alarm", ALARM=SOLALARM, RC=STATUS )
!        VERIFY_(STATUS)
      if(status==0) then 
! Get the interval of the solar alarm - first get it in seconds
        call ESMF_ConfigGetAttribute ( CF, DT_SOLAR, Label="SOLAR_DT:", DEFAULT=DT, RC=STATUS )
        VERIFY_(STATUS)
! Now make an ESMF interval from the increment in seconds
        CALL ESMF_TimeIntervalSet ( TINT, S=NINT(DT_SOLAR), RC=STATUS )
        VERIFY_(STATUS)
! Now print out the solar alarm interval
        if (MAPL_AM_I_Root(VM).and.debugzth) CALL ESMF_TimeIntervalPrint ( TINT, OPTIONS="string", RC=STATUS )
! Now find out if it is ringing now: if so, set "BEFORE" to last time it rang before now
         solalarmison = ESMF_AlarmIsRinging(SOLALARM,RC=STATUS)
         VERIFY_(STATUS)
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' logical for solar alarm ',solalarmison
!     if so, set "BEFORE" to last time it rang before now
        if(solalarmison) then
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm is ringing '
         NOW = CURRENT_TIME
         BEFORE = NOW - TINT
! Now print out the last time solar alarm rang
         if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", RC=STATUS )
!     If alarm is not ringing now, find out when it rang last
        else
         if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm is not ringing '
         call ESMF_AlarmGet ( SOLALARM, prevRingTime=BEFORE, RC=STATUS )
         VERIFY_(STATUS)
! PrevRingTime can lie: if alarm never went off yet it gives next alarm time, not prev.
         if(BEFORE > CURRENT_TIME) then
          BEFORE = BEFORE-TINT
          if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm not ringing, prev time lied '
          if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", RC=STATUS )
         else
          if (MAPL_AM_I_Root(VM).and.debugzth)print *,' In catch, solar alarm not ringing, prev time okay '
          if (MAPL_AM_I_Root(VM).and.debugzth)CALL ESMF_TimePrint ( BEFORE, OPTIONS="string", RC=STATUS )
         endif
! Now print out the last time solar alarm rang
        endif
     else
        BEFORE = CURRENT_TIME
        TINT = DELT
     end if
! Get the zenith angle at the center of the time between the last solar call and the next one
        call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR, &
            INTV = TINT,     &
            currTime=BEFORE+DELT,  &
            RC=STATUS )
        VERIFY_(STATUS)

        ZTH = max(0.0,ZTH)
        ZVG  = fveg1*(Z2CH - SCALE4Z0*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI1)) + &
             fveg2*(Z2CH - SCALE4Z0*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI2)) 


        !  For now roughnesses and displacement heights
        !   are the same for all subtiles.
        !---------------------------------------------------

        Z0   = Z0_BY_ZVEG*ZVG
        IF (USE_ASCATZ0 == 1) WHERE (NDVI <= 0.2) Z0 = ASCATZ0        
        D0   = D0_BY_ZVEG*ZVG

        UUU = max(UU,MAPL_USMIN) * (log((ZVG-D0+Z0)/Z0) &
             / log((max(DZ-D0,10.)+Z0)/Z0))

        !--------------- GOSWIM IMPORTS FROM GOCART ---------------
        ! Initialization 
        RCONSTIT(:,:,:)  = 0.0
        TOTDEPOS(:,:) = 0.0
        RMELT(:,:)  = 0.0
        !------------------------------------------------------------------

        ! Zero the light-absorbing aerosol (LAA) deposition rates from  GOCART:

        select case (AEROSOL_DEPOSITION)
        case (0)
           DUDP(:,:)=0.
           DUSV(:,:)=0.
           DUWT(:,:)=0.
           DUSD(:,:)=0.
           BCDP(:,:)=0.
           BCSV(:,:)=0.
           BCWT(:,:)=0.
           BCSD(:,:)=0.
           OCDP(:,:)=0.
           OCSV(:,:)=0.
           OCWT(:,:)=0.
           OCSD(:,:)=0.
           
        case (2)
           DUDP(:,:)=0.
           DUSV(:,:)=0.
           DUWT(:,:)=0.
           DUSD(:,:)=0.
           
        case (3)
           BCDP(:,:)=0.
           BCSV(:,:)=0.
           BCWT(:,:)=0.
           BCSD(:,:)=0.
           
        case (4)
           OCDP(:,:)=0.
           OCSV(:,:)=0.
           OCWT(:,:)=0.
           OCSD(:,:)=0.
           
        end select

! Convert the dimentions for LAAs from GEOS_SurfGridComp.F90 to GEOS_LandIceGridComp.F90
! Note: Explanations of each variable
! TOTDEPOS(:,1): Combined dust deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,2): Combined dust deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,3): Combined dust deposition from size bin 3 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,4): Combined dust deposition from size bin 4 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,5): Combined dust deposition from size bin 5 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,6): Combined hydrophobic BC deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,7): Combined hydrophilic BC deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,8): Combined hydrophobic OC deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,9): Combined hydrophilic OC deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
!============================= Possible future applications ====================================
! TOTDEPOS(:,10): Combined sulfate deposition from size bin 3 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,11): Combined sea salt deposition from size bin 1 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,12): Combined sea salt deposition from size bin 2 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,13): Combined sea salt deposition from size bin 3 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,14): Combined sea salt deposition from size bin 4 (dry, conv-scav, ls-scav, sed)
! TOTDEPOS(:,15): Combined sea salt deposition from size bin 5 (dry, conv-scav, ls-scav, sed)

        TOTDEPOS(:,1) = DUDP(:,1) + DUSV(:,1) + DUWT(:,1) + DUSD(:,1)
        TOTDEPOS(:,2) = DUDP(:,2) + DUSV(:,2) + DUWT(:,2) + DUSD(:,2)
        TOTDEPOS(:,3) = DUDP(:,3) + DUSV(:,3) + DUWT(:,3) + DUSD(:,3)
        TOTDEPOS(:,4) = DUDP(:,4) + DUSV(:,4) + DUWT(:,4) + DUSD(:,4)
        TOTDEPOS(:,5) = DUDP(:,5) + DUSV(:,5) + DUWT(:,5) + DUSD(:,5)
        TOTDEPOS(:,6) = BCDP(:,1) + BCSV(:,1) + BCWT(:,1) + BCSD(:,1)
        TOTDEPOS(:,7) = BCDP(:,2) + BCSV(:,2) + BCWT(:,2) + BCSD(:,2)
        TOTDEPOS(:,8) = OCDP(:,1) + OCSV(:,1) + OCWT(:,1) + OCSD(:,1)
        TOTDEPOS(:,9) = OCDP(:,2) + OCSV(:,2) + OCWT(:,2) + OCSD(:,2)
!============================= Possible future applications ====================================
!        TOTDEPOS(:,10) = SUDP(:,1) + SUSV(:,1) + SUWT(:,1) + SUSD(:,1)
!        TOTDEPOS(:,11) = SSDP(:,1) + SSSV(:,1) + SSWT(:,1) + SSSD(:,1)
!        TOTDEPOS(:,12) = SSDP(:,2) + SSSV(:,2) + SSWT(:,2) + SSSD(:,2)
!        TOTDEPOS(:,13) = SSDP(:,3) + SSSV(:,3) + SSWT(:,3) + SSSD(:,3)
!        TOTDEPOS(:,14) = SSDP(:,4) + SSSV(:,4) + SSWT(:,4) + SSSD(:,4)
!        TOTDEPOS(:,15) = SSDP(:,5) + SSSV(:,5) + SSWT(:,5) + SSSD(:,5)

! --------------- GOSWIM PROGRNOSTICS ---------------------------

! Conversion of the masses of the snow impurities
! Note: Explanations of each variable
! Number of snow layer is 15: N = 1-15
! RCONSTIT(NTILES,N,1): Dust mass from bin 1 in layer N
! RCONSTIT(NTILES,N,2): Dust mass from bin 2 in layer N
! RCONSTIT(NTILES,N,3): Dust mass from bin 3 in layer N
! RCONSTIT(NTILES,N,4): Dust mass from bin 4 in layer N
! RCONSTIT(NTILES,N,5): Dust mass from bin 5 in layer N
! RCONSTIT(NTILES,N,6): Hydrophobic BC mass from bin 1 in layer N
! RCONSTIT(NTILES,N,7): Hydrophilic BC mass from bin 2 in layer N
! RCONSTIT(NTILES,N,8): Hydrophobic OC mass from bin 1 in layer N
! RCONSTIT(NTILES,N,9): Hydrophilic OC mass from bin 2 in layer N
!============================= Possible future applications ====================================
! RCONSTIT(NTILES,N,10): Sulfate mass from size bin 3 in layer N
! RCONSTIT(NTILES,N,11): Sea salt mass from size bin 1 in layer N
! RCONSTIT(NTILES,N,12): Sea salt mass from size bin 2 in layer N
! RCONSTIT(NTILES,N,13): Sea salt mass from size bin 3 in layer N
! RCONSTIT(NTILES,N,14): Sea salt mass from size bin 4 in layer N
! RCONSTIT(NTILES,N,15): Sea salt mass from size bin 5 in layer N

        if (DO_GOSWIM /= 0) then
           RCONSTIT(:,:,1) = RDU001(:,:)
           RCONSTIT(:,:,2) = RDU002(:,:)
           RCONSTIT(:,:,3) = RDU003(:,:)
           RCONSTIT(:,:,4) = RDU004(:,:)
           RCONSTIT(:,:,5) = RDU005(:,:)
           RCONSTIT(:,:,6) = RBC001(:,:)
           RCONSTIT(:,:,7) = RBC002(:,:)
           RCONSTIT(:,:,8) = ROC001(:,:)
           RCONSTIT(:,:,9) = ROC002(:,:)

!============================= Possible future applications ====================================
!        RCONSTIT(:,:,10) = RSU003(:,:)
!        RCONSTIT(:,:,11) = RSS001(:,:)
!        RCONSTIT(:,:,12) = RSS002(:,:)
!        RCONSTIT(:,:,13) = RSS003(:,:)
!        RCONSTIT(:,:,14) = RSS004(:,:)
!        RCONSTIT(:,:,15) = RSS005(:,:)
        endif

        ! --------------------------------------------------------------------------
        ! Update raditation exports
        ! --------------------------------------------------------------------------

        allocate ( ALBVR_tmp(ntiles) )
        allocate ( ALBNR_tmp(ntiles) )
        allocate ( ALBVF_tmp(ntiles) )
        allocate ( ALBNF_tmp(ntiles) )
        allocate ( SNOVR_tmp(ntiles) )
        allocate ( SNONR_tmp(ntiles) )
        allocate ( SNOVF_tmp(ntiles) )
        allocate ( SNONF_tmp(ntiles) )

        call    SIBALB(NTILES, VEG1,LAI1,GRN, ZTH,         & 
                       BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
                       ALBVR, ALBNR, ALBVF, ALBNF, MODIS_SCALE=.TRUE.  )         ! instantaneous snow-free albedos on tiles
        call STIEGLITZSNOW_CALC_TPSNOW(NTILES, HTSNNN(1,:), WESNN(1,:), TPSN1OUT1, FICE1)
        TPSN1OUT1 =  TPSN1OUT1 + Tzero

        call   SNOW_ALBEDO(NTILES,N_snow, N_CONST_LAND4SNWALB, VEG1, LAI1, ZTH,        &
                 RHOFS,                                              &   
                 SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                & 
                 WESNN, HTSNNN, SNDZN,                               &
                 ALBVR, ALBNR, ALBVF, ALBNF,                         & ! instantaneous snow-free albedos on tiles
                 SNOVR, SNONR, SNOVF, SNONF,                         & ! instantaneous snow albedos on tiles
                 RCONSTIT, UUU, TPSN1OUT1, DRPAR, DFPAR)           

        call    SIBALB(NTILES, VEG2,LAI2,GRN, ZTH,         & 
                       BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
                       ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp, MODIS_SCALE=.TRUE. ) ! instantaneous snow-free albedos on tiles

        call   SNOW_ALBEDO(NTILES,N_snow, N_CONST_LAND4SNWALB, VEG2, LAI2, ZTH,        &
                 RHOFS,                                              &   
                 SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                & 
                 WESNN, HTSNNN, SNDZN,                               &
                 ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp, & ! instantaneous snow-free albedos on tiles
                 SNOVR_tmp, SNONR_tmp, SNOVF_tmp, SNONF_tmp, & ! instantaneous snow albedos on tiles
                 RCONSTIT, UUU, TPSN1OUT1, DRPAR, DFPAR )  

        ALBVR(:) = ALBVR(:)*fveg1(:) + ALBVR_tmp(:)*fveg2(:)
        ALBNR(:) = ALBNR(:)*fveg1(:) + ALBNR_tmp(:)*fveg2(:)
        ALBVF(:) = ALBVF(:)*fveg1(:) + ALBVF_tmp(:)*fveg2(:)
        ALBNF(:) = ALBNF(:)*fveg1(:) + ALBNF_tmp(:)*fveg2(:)

        SNOVR(:) = SNOVR(:)*fveg1(:) + SNOVR_tmp(:)*fveg2(:)
        SNONR(:) = SNONR(:)*fveg1(:) + SNONR_tmp(:)*fveg2(:)
        SNOVF(:) = SNOVF(:)*fveg1(:) + SNOVF_tmp(:)*fveg2(:)
        SNONF(:) = SNONF(:)*fveg1(:) + SNONF_tmp(:)*fveg2(:)

        ! --------------------------------------------------------------------------
        ! albedo/swnet partitioning
        ! --------------------------------------------------------------------------

        VSUVR = DRPAR + DRUVR
        VSUVF = DFPAR + DFUVR

        if(associated(SWDOWNLAND)) SWDOWNLAND = DRPAR + DFPAR + DRUVR + DFUVR + DRNIR + DFNIR

        SWNETFREE = (1.-ALBVR)*VSUVR + (1.-ALBVF)*VSUVF + (1.-ALBNR)*DRNIR + (1.-ALBNF)*DFNIR 
        SWNETSNOW = (1.-SNOVR)*VSUVR + (1.-SNOVF)*VSUVF + (1.-SNONR)*DRNIR + (1.-SNONF)*DFNIR 

        ! --------------------------------------------------------------------------
        ! Parameters that depend on vegetation type only                             gkw: these are not used in unified
        ! --------------------------------------------------------------------------

        RSL1   = VGRDRS(VEG1)/(ROOTL*VGROTD(VEG1)) 

        RSL2   = ROOTL*VGROCA(VEG1)
        RSL2   = (RSL2 - 3.0 - 2.*alog(RSL2/(1.-RSL2)))/(8.*MAPL_PI*ROOTL*VGROTD(VEG1))

        ! --------------------------------------------------------------------------
        ! Greenness and type dependent parameters
        ! --------------------------------------------------------------------------

        SQSCAT = fveg1*((VGTR11(VEG1)+VGRF11(VEG1)) * GRN  + (VGTR12(VEG1)+VGRF12(VEG1)) * (1.-GRN)) + & 
                 fveg2*((VGTR11(VEG2)+VGRF11(VEG2)) * GRN  + (VGTR12(VEG2)+VGRF12(VEG2)) * (1.-GRN))
        SQSCAT = sqrt(1.0 - SQSCAT)

        ! --------------------------------------------------------------------------
        ! LAI and type dependent parameters; RDC formulation now uses veg fractions gkw: 2013-11-25, see note from Randy
        ! --------------------------------------------------------------------------

        RDC = max(VGRDA(VEG1),VGRDA(VEG2))*min(1.,lai/2.)
        RDC = max(RDC,0.001)

        RHO = PS/(MAPL_RGAS*(TA*(1+MAPL_VIREPS*QA)))

        DEDTC=0.0
        DHSDQA=0.0

        if(is_OFFLINE) then
           do N=1,NUM_SUBTILES
              CFT   (:,N) = 1.0
              CFQ   (:,N) = 1.0
              SHSBT (:,N) = MAPL_CP*CH(:,N)*(TC(:,N)-TA)
              EVSBT (:,N) = CQ(:,N)*(QC(:,N)-QA)
              DSHSBT(:,N) = MAPL_CP*CH(:,N)
              DEVSBT(:,N) = CQ(:,N)
              BLWN(:,N) = EMIS*MAPL_STFBOL*TC(:,N)*TC(:,N)*TC(:,N)
              ALWN(:,N) = -3.0*BLWN(:,N)*TC(:,N)
              BLWN(:,N) =  4.0*BLWN(:,N)
           end do
           if(CHOOSEMOSFC==0 .and. incl_Louis_extra_derivs ==1) then
              do N=1,NUM_SUBTILES
                 DEVSBT(:,N)=CQ(:,N)+max(0.0,-DCQ(:,N)*MAPL_VIREPS*TC(:,N)*(QC(:,N)-QA))
                 DEDTC(:,N) =max(0.0,-DCQ(:,N)*(1.+MAPL_VIREPS*QC(:,N))*(QC(:,N)-QA))
                 DSHSBT(:,N)=MAPL_CP*(CH(:,N)+max(0.0,-DCH(:,N)*(1.+MAPL_VIREPS*QC(:,N))*(TC(:,N)-TA)))
                 DHSDQA(:,N)=max(0.0,-MAPL_CP*DCH(:,N)*MAPL_VIREPS*TC(:,N)*(TC(:,N)-TA))
              enddo
           endif
        else
           do N=1,NUM_SUBTILES
              CFT   (:,N) = (CH(:,N)/CTATM)
              CFQ   (:,N) = (CQ(:,N)/CQATM)
              SHSBT (:,N) = (SH  + DSH  *(TC(:,N)-THATM))*CFT(:,N)
              EVSBT (:,N) = (EVAP+ DEVAP*(QC(:,N)-QHATM))*CFQ(:,N)
              DSHSBT(:,N) =  DSH  *CFT(:,N)
              DEVSBT(:,N) =  DEVAP*CFQ(:,N)
              ALWN(:,N)=ALW
              BLWN(:,N)=BLW
           end do
        end if

        do N=1,NUM_SUBTILES
           DQS(:,N) = GEOS_DQSAT ( TC(:,N), PS, QSAT=QSAT(:,N), PASCALS=.true., RAMP=0.0 )
           QC (:,N) = min(max(QA(:),QSAT(:,N)),QC(:,N))
           QC (:,N) = max(min(QA(:),QSAT(:,N)),QC(:,N))
           RA (:,N) = RHO/CH(:,N)
        end do


        QC(:,FSNW) = QSAT(:,FSNW)

	! --------------------------------------------------------------------------
        ! get total solid precip
        ! --------------------------------------------------------------------------

        SLDTOT = SNO+ICE+FRZR
        
	! --------------------------------------------------------------------------
	! protect the forcing from unsavory values, as per practice in offline
	! driver
	! --------------------------------------------------------------------------

        ASSERT_(count(PLS<0.)==0)
        ASSERT_(count(PCU<0.)==0)
        ASSERT_(count(SLDTOT<0.)==0)

        LAI0  = max(0.0001     , LAI)
        GRN0  = max(0.0001     , GRN)		
        ZTH   = max(0.0001     , ZTH)

        TCO   = TC
        QCO   = QC

        ! --------------------------------------------------------------------------
        ! actual CATCHMENT call
        ! --------------------------------------------------------------------------

        TILEZERO = 0.0

        call MAPL_TimerOn  ( MAPL, "-CATCH" )

#ifdef DBG_CATCH_INPUTS
        call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_LocStreamGet(LOCSTREAM, TILEGRID=TILEGRID, RC=STATUS)
        VERIFY_(STATUS)

        call MAPL_TileMaskGet(tilegrid,  mask, rc=status)
        VERIFY_(STATUS)

         if (UNIT_i == 0) then
           unit_i = GETFILE( "catch_inputs.data", form="unformatted", RC=STATUS )
           VERIFY_(STATUS)
        endif
        unit = unit_i

! Inputs
        call MAPL_VarWrite(unit, tilegrid, PCU, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, PLS, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SNO, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ICE, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, FRZR, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, UUU, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, EVSBT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, EVSBT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT(:,FTRN),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, EVSBT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, EVSBT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, TILEZERO      ,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        
        call MAPL_VarWrite(unit, tilegrid, TA, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QA, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, RA(:,FSAT),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FTRN),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FWLT),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, ZTH,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DRPAR,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DFPAR,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SWNETFREE,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SWNETSNOW,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, LWDNSRF, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, PS*.01, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, LAI0,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, GRN0,  mask=mask, rc=status); VERIFY_(STATUS)
!       call MAPL_VarWrite(unit, tilegrid, Z2CH,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SQSCAT,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RSL1,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RSL2,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RDC, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FTRN) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FWLT) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FSNW) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FSNW) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWX, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWX, mask=mask, rc=status); VERIFY_(STATUS)

! params
        if (firsttime) then
            firsttime = .false.
           unit = GETFILE( "catch_params.data", form="unformatted", RC=STATUS )
           VERIFY_(STATUS)

           NT_GLOBAL = size(mask)

           call WRITE_PARALLEL(NT_GLOBAL, UNIT)
           call WRITE_PARALLEL(DT, UNIT)
           call WRITE_PARALLEL(PRECIPFRAC, UNIT)
           call MAPL_VarWrite(unit, tilegrid, VEG1, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, VEG2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BF1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BF2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BF3,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, VGWMAX,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CDCR1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CDCR2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, PSIS, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BEE,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, POROS,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WPWET,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, COND,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GNU, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARS1, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARS2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARS3, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA1, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA3, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA4, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW3,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW4,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSA1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSA2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSB1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSB2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ATAU, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BTAU, mask=mask, rc=status); VERIFY_(STATUS)

           call FREE_FILE(unit, RC=STATUS)
           VERIFY_(STATUS)

! Updates
           unit = GETFILE( "catch_updates.data", form="unformatted", RC=STATUS )
           VERIFY_(STATUS)


           call MAPL_VarWrite(unit, tilegrid, TC(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TC(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TC(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, QC(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, QC(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, QC(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TG(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TG(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TG(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CAPAC, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CATDEF, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, RZEXC, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SRFEXC, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(1,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(2,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(3,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(4,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(5,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(6,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSURF, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WESNN(1,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WESNN(2,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WESNN(3,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, HTSNNN(1,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, HTSNNN(2,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, HTSNNN(3,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SNDZN(1,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SNDZN(2,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SNDZN(3,:), mask=mask, rc=status); VERIFY_(STATUS)
           
           call FREE_FILE(unit, RC=STATUS)
           VERIFY_(STATUS)

        end if
        deallocate(mask)
#endif

! ----------------------------------------------------------------------------------------

! gkw: start on main CN block

    allocate(   btran(ntiles) )
    allocate(     wgt(ntiles) )
    allocate(     bt1(ntiles) )
    allocate(     bt2(ntiles) )
    allocate(     bt4(ntiles) )
    allocate(     wpp(ntiles) )
    allocate(    fwet(ntiles) )
    allocate(     sm1(ntiles) )
    allocate(     sm2(ntiles) )
    allocate(     sm4(ntiles) )
    allocate(  btran1(ntiles) )
    allocate(  btran2(ntiles) )
    allocate(  btran3(ntiles) )
    allocate(     tcx(ntiles) )
    allocate(     qax(ntiles) )
    allocate(     rcx(ntiles) )
    allocate(   rcxdt(ntiles) )
    allocate(   rcxdq(ntiles) )
    allocate(     tx1(ntiles) )
    allocate(     tx2(ntiles) )
    allocate(     tx3(ntiles) )
    allocate(     qx1(ntiles) )
    allocate(     qx2(ntiles) )
    allocate(     qx3(ntiles) )
    allocate(    car1(ntiles) )
    allocate(    car2(ntiles) )
    allocate(    car4(ntiles) )
    allocate( parzone(ntiles) )
    allocate(    para(ntiles) )
    allocate(  totwat(ntiles) )
    allocate(     npp(ntiles) )
    allocate(     gpp(ntiles) )
    allocate(      sr(ntiles) )
    allocate(     nee(ntiles) )
    allocate(    padd(ntiles) )
    allocate(    root(ntiles) )
    allocate(    vegc(ntiles) )
    allocate(    xsmr(ntiles) )
    allocate(    burn(ntiles) )
    allocate(    fsel(ntiles) )
    allocate(   closs(ntiles) )
    allocate(    dayl(ntiles) )
    allocate(dayl_fac(ntiles) )
    allocate(CO2V    (ntiles) )

    allocate(     tgw(ntiles,nzone) )
    allocate(     rzm(ntiles,nzone) )
    allocate(    rc00(ntiles,nzone) )
    allocate(    rcdt(ntiles,nzone) )
    allocate(    rcdq(ntiles,nzone) )
    allocate( totcolc(ntiles,nzone) )

    allocate( psnsunx(ntiles,nveg) )
    allocate( psnshax(ntiles,nveg) )
    allocate( sifsunx(ntiles,nveg) )
    allocate( sifshax(ntiles,nveg) )
    allocate( laisunx(ntiles,nveg) )
    allocate( laishax(ntiles,nveg) )
    allocate(    elaz(ntiles,nveg) )
    allocate(    esaz(ntiles,nveg) )
    allocate(    fvez(ntiles,nveg) )
    allocate(    ityz(ntiles,nveg) )

    allocate( psnsun(ntiles,nveg,nzone) )
    allocate( psnsha(ntiles,nveg,nzone) )
    allocate( laisun(ntiles,nveg,nzone) )
    allocate( laisha(ntiles,nveg,nzone) )

    allocate(      ht(N_gt) )
    allocate(      tp(N_gt) )
    allocate( soilice(N_gt) )

! get current date & time gkw: this is used to transfer CN restart vars & set declination
! -----------------------
    call ESMF_TimeGet  ( CURRENT_TIME, YY = AGCM_YY, &
				       MM = AGCM_MM, &
				       DD = AGCM_DD, &
				       S  = AGCM_S , &  ! S is second of day (because H & M are absent)
				       dayOfYear = dofyr , &
				       rc=status )
    VERIFY_(STATUS)

! declination  gkw: this is ugly... get someone to make ZS & ZC available as optional arg for MAPL_SunGetInsolation
! -----------                                                                              or MAPL_SunOrbitQuery
    call MAPL_SunOrbitQuery(Orbit,zc=zco,zs=zso,years_per_cycle=years_per_cycle, &
           days_per_cycle=days_per_cycle,year_length=year_length,rc=status)
    VERIFY_(STATUS)
    YEAR = mod(AGCM_YY-1,YEARS_PER_CYCLE)  ! gkw: made ORBIT public in MAPL_sun_uc.P90 (temporary solution)

    IDAY = YEAR*int(YEAR_LENGTH)+dofyr
    IDAYP1 = mod(IDAY,DAYS_PER_CYCLE) + 1

    ASSERT_(IDAY   <= 1461 .AND. IDAY   > 0)
    ASSERT_(IDAYP1 <= 1461 .AND. IDAYP1 > 0)

    FAC = real(AGCM_S)/86400.

    !ZS = ORBIT%ZS(IDAYP1)*FAC + ORBIT%ZS(IDAY)*(1.-FAC) !   sine of solar declination
    !ZC = ORBIT%ZC(IDAYP1)*FAC + ORBIT%ZC(IDAY)*(1.-FAC) ! cosine of solar declination
    ZS = ZSO(IDAYP1)*FAC + ZSO(IDAY)*(1.-FAC) !   sine of solar declination
    ZC = ZCO(IDAYP1)*FAC + ZCO(IDAY)*(1.-FAC) ! cosine of solar declination

    nullify(ZSO,ZCO) 
    declin = asin(ZS)

!if( MAPL_AM_I_Root() ) write(6,444) declin,iday,idayp1,fac,zs,dofyr,AGCM_S
!444 format('orbit:',f10.6,2i5,f7.4,f9.5,i4,i6)

! get ending time; determine if this is last call before ending time
! ------------------------------------------------------------------
    call ESMF_ClockGet ( clock,  StopTime=StopTime ,rc=STATUS )
    VERIFY_(STATUS)

    NextTime = CURRENT_TIME + DELT

    if(ntiles > 0) then ! gkw: skip threads with no land tiles

! gkw: assign new vegetation types and fractions
! ----------------------------------------------
    cat_id = nint(tile_id) ! gkw: temporary for debugging

! compute daylength (and daylength factor) gkw: LATS is in radians, OB in degrees
! ----------------------------------------
    call MAPL_SunOrbitQuery(ORBIT,           & ! gkw: this is the correct way to obtain obliquity
                            OBLIQUITY=OB,    & ! gkw: ideally, ZS & ZC would be obtained this way
			    rc=status )
    VERIFY_(STATUS)

    do n = 1,ntiles
!!!   fac = -(sin(lats(n))*sin(declin))/(cos(lats(n))*cos(declin))
      fac = -(sin(lats(n))*zs)/(cos(lats(n))*zc)
      fac = min(1.,max(-1.,fac))
      dayl(n) = (86400./MAPL_PI) * acos(fac)   ! daylength (seconds)

      max_decl = ob * (MAPL_PI/180.)
      if (lats(n) < 0.) max_decl = -max_decl
      fac = -(sin(lats(n))*sin(max_decl))/(cos(lats(n))*cos(max_decl))
      fac = min(1.,max(-1.,fac))
      max_dayl = (86400./MAPL_PI) * acos(fac)  ! maximum daylength (sec; function of latitude & obliquity)

! calculate dayl_factor as the ratio of (current:max dayl)^2
! set a minimum of 0.01 (1%) for the dayl_factor [gkw: from CLM4]

      dayl_fac(n) = min(1.,max(0.01,(dayl(n)*dayl(n))/(max_dayl*max_dayl)))
    end do

!write(6,845) cat_id(1),lats(1),lons(1),dayl(1),-(sin(lats(1))*sin(declin))/(cos(lats(1))*cos(declin)), &
!                                               -(sin(lats(1))*zs)/(cos(lats(1))*zc),ob,dayl_fac(1)
!845 format('dayl:',i6,2f10.5,f9.2,2f10.6,f7.3,f8.5)

! gkw: obtain catchment area fractions and soil moisture
! ------------------------------------------------------
    call catch_calc_soil_moist( ntiles, veg1, dzsf, vgwmax, cdcr1, cdcr2, psis, bee, poros, wpwet, &
                              ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1, arw2, arw3, arw4,    &
                              srfexc, rzexc, catdef, car1, car2, car4, sfmc, rzmc, prmc )
                              
! obtain saturated canopy resistance following Farquhar, CLM4 implementation    

! compute RC & PSN in each of the CN zones
! ----------------------------------------

! "btran" in the catchment zones; map into CN zones
! -------------------------------------------------
       sm1 = 1.0
       bt1 = 1.0          ! saturated area

! gkw: bt2 is unstressed region only (subtract saturated and wilting areas)
       do n = 1,ntiles
         if(car2(n) > 0.) then
           bt2(n)=(rzmc(n)/poros(n) - car1(n) - car4(n)*wpwet(n))/car2(n)
          else
           bt2(n)= rzmc(n)/poros(n)
         endif
         bt2(n) = max(bt2(n),wpwet(n))
         bt2(n) = min(bt2(n),1.)
       end do

       sm2 = bt2
       bt2 = bt2**(-bee)

       wpp = wpwet ** (-bee)    ! soil water potential at wilting point
       bt2 = (bt2-wpp)/(1.-wpp) ! unstressed area (using soil water potential; PSIS common to all terms)

!      sm4 = wpwet - 1.e-7      ! wilting point minus epsilon (to account for truncation)

       do n = 1,ntiles
         if(car4(n) > 0.) then
           sm4(n)=(rzmc(n)/poros(n) - car1(n) - car2(n)*sm2(n))/car4(n)
          else
           sm4(n)= wpwet(n)
         endif
         sm4(n) = max(sm4(n),1.e-3)
         sm4(n) = min(sm4(n),wpwet(n)-1.e-7)
       end do

       bt4 = 0.0          ! stressed area

       do n = 1,ntiles

         ax1 = car1(n)
	 ax2 = car2(n)
         ax4 = 1. - ax1 - ax2

	 cn1 = wtzone(n,1)
	 cn2 = wtzone(n,2)
	 cn3 = wtzone(n,3)

! CN zone 1
         if(ax1 .gt. cn1) then
           f1 = cn1 ; f2 = 0. ; f4 = 0.
          else
           if((ax1+ax2) .gt. cn1) then
             f1 = ax1 ; f2 = cn1-ax1 ; f4 = 0.
            else
             f1 = ax1 ; f2 = ax2 ; f4 = cn1-ax1-ax2
           endif
         endif

         btran1(n) = (f1*bt1(n)     + f2*bt2(n)     + f4*bt4(n)    )/cn1
	 tgw(n,1) =  (f1*tg(n,fsat) + f2*tg(n,ftrn) + f4*tg(n,fwlt))/cn1
	 tx1(n)   =  (f1*tc(n,fsat) + f2*tc(n,ftrn) + f4*tc(n,fwlt))/cn1
	 qx1(n)   =  (f1*qc(n,fsat) + f2*qc(n,ftrn) + f4*qc(n,fwlt))/cn1
         rzm(n,1) =  (f1*sm1(n)     + f2*sm2(n)     + f4*sm4(n)    )/cn1

! CN zone 2
         if(ax1 .gt. cn1) then
           cn12 = cn1 + cn2
           if(car1(n) .gt. cn12) then
             f1 = cn2 ; f2 = 0. ; f4 = 0.
            else
             if((ax1+ax2) .lt. cn12) then
               f1 = ax1-cn1 ; f2 = ax2 ; f4 = cn12-ax1-ax2
              else
               f1 = ax1-cn1 ; f2 = cn12-ax1 ; f4 = 0.
             endif
           endif
          else
           cn23 = cn2 + cn3
           if(ax4 .gt. cn23) then
             f1 = 0. ; f2 = 0. ; f4 = cn2
            else
             if(ax4 .lt. cn3) then
               f1 = 0. ; f2 = cn2 ; f4 = 0.
              else
               f1 = 0. ; f2 = cn23-ax4 ; f4 = ax4-cn3
             endif
           endif
         endif

         btran2(n) = (f1*bt1(n)     + f2*bt2(n)     + f4*bt4(n)    )/cn2
	 tgw(n,2) =  (f1*tg(n,fsat) + f2*tg(n,ftrn) + f4*tg(n,fwlt))/cn2
	 tx2(n)   =  (f1*tc(n,fsat) + f2*tc(n,ftrn) + f4*tc(n,fwlt))/cn2
	 qx2(n)   =  (f1*qc(n,fsat) + f2*qc(n,ftrn) + f4*qc(n,fwlt))/cn2
         rzm(n,2) =  (f1*sm1(n)     + f2*sm2(n)     + f4*sm4(n)    )/cn2

! CN zone 3
         if(ax4 .gt. cn3) then
           f1 = 0. ; f2 = 0. ; f4 = cn3
          else
           if((ax4+ax2) .gt. cn3) then
             f1 = 0. ; f2 = cn3-ax4 ; f4 = ax4
            else
             f1 = cn3-ax4-ax2 ; f2 = ax2 ; f4 = ax4
           endif
         endif

         btran3(n) = (f1*bt1(n)     + f2*bt2(n)     + f4*bt4(n)    )/cn3
	 tgw(n,3) =  (f1*tg(n,fsat) + f2*tg(n,ftrn) + f4*tg(n,fwlt))/cn3
	 tx3(n)   =  (f1*tc(n,fsat) + f2*tc(n,ftrn) + f4*tc(n,fwlt))/cn3
	 qx3(n)   =  (f1*qc(n,fsat) + f2*qc(n,ftrn) + f4*qc(n,fwlt))/cn3
         rzm(n,3) =  (f1*sm1(n)     + f2*sm2(n)     + f4*sm4(n)    )/cn3

       end do

! soil temperature and hydrologic state
! -------------------------------------
    DO N=1,ntiles

! soil temperatures
! -----------------
      zbar = -sqrt(1.e-20+catdef(n)/bf1(n))+bf2(n)
      HT(:)=GHTCNT(:,N)
      CALL GNDTMP_CN(poros(n),zbar,ht,frice,tp,soilice)
      tp1(n) = tp(1) + Tzero
      tp2(n) = tp(2) + Tzero
      tp3(n) = tp(3) + Tzero
      tp4(n) = tp(4) + Tzero
      tp5(n) = tp(5) + Tzero
      tp6(n) = tp(6) + Tzero

! total soil liquid water
! -----------------------
      totwat(n) = cdcr2(n) - catdef(n) + rzexc(n) + srfexc(n)
      totwat(n) = totwat(n)*(1. - frice)

! baseflow
! --------
      zbar = sqrt(1.e-20+catdef(n)/bf1(n))-bf2(n)
      bflow(n) = (1.-frice)*1000.* &
    	    cond(n)*exp(-(bf3(n)-ashift)-gnu(n)*zbar)/gnu(n)
      IF(catdef(n) >= cdcr1(n)) bflow(n) = 0.
      bflow(n) = min(cond(n),bflow(n))
    end do

! get CO2
! -------
    call MAPL_GetResource( MAPL, CO2, 'CO2:', default=350.e-6, RC=STATUS)
    VERIFY_(STATUS)

    if(CO2 < 0.0) CO2 = GETCO2(AGCM_YY,dofyr)

    CO2V (:) = CO2
    
! use CO2SC from GOCART/CO2
! -------------------------
    
    IF (DO_CO2SC /= 0) THEN
       
       where ((CO2SC >= 0.) .and. (CO2SC <= 1000.))
          CO2V = CO2SC * 1e-6
       end where
       
    endif

    if(associated(BTRANT)) btrant = 0. 
    if(associated(SIF))    sif    = 0. 

! fraction of foliage that is wet  gkw 20140327
! -------------------------------
    do n = 1,ntiles
      if(lai(n) > 1.e-4) then
        fwet(n) = min(1.,max(0.,capac(n)/(0.2*lai(n))))
       else
        fwet(n) = 0.
      endif
    end do

    para(:) = 0. ! zero out absorbed PAR summing array

    do nz = 1,nzone

      if(nz == 1) then
        btran = btran1
        tcx = tx1
        qax = qx1
      endif

      if(nz == 2) then
        btran = btran2
        tcx = tx2
        qax = qx2
      endif

      if(nz == 3) then
        btran = btran3
        tcx = tx3
        qax = qx3
      endif

      do nv = 1,nveg
        elaz(:,nv) = elai(:,nv,nz)
        esaz(:,nv) = esai(:,nv,nz)
        ityz(:,nv) = ityp(:,nv,nz)
        fvez(:,nv) = fveg(:,nv,nz)
      end do

      do n = 1,ntiles
        if(tp1(n) < (Tzero-0.01)) btran(n) = 0. ! no photosynthesis if ground fully frozen
      end do

      call compute_rc(ntiles,nveg,TCx,QAx,TA,PS,ZTH,DRPAR,DFPAR,            &
                      elaz,esaz,ityz,fvez,btran,fwet,                       &
                      RCx,RCxDT,RCxDQ,psnsunx,psnshax,laisunx,laishax,      &
                      dayl_fac,co2v,dtc,dea,parzone,sifsunx,sifshax)

      rc00(:,nz) = rcx(:)
      rcdt(:,nz) = rcxdt(:)
      rcdq(:,nz) = rcxdq(:)

      psnsun(:,:,nz) = psnsunx(:,:)
      psnsha(:,:,nz) = psnshax(:,:)
      laisun(:,:,nz) = laisunx(:,:)
      laisha(:,:,nz) = laishax(:,:)

      para(:)   = para(:) + parzone(:)*wtzone(:,nz)
      if(associated(BTRANT)) btrant(:) = btrant(:) + btran(:)*wtzone(:,nz)
      if(associated(SIF)) then
        do nv = 1,nveg
          sif(:) = sif(:) + wtzone(:,nz)*fvez(:,nv)*(sifsunx(:,nv)*laisunx(:,nv) + sifshax(:,nv)*laishax(:,nv))
        end do
      endif

    end do

    deallocate (co2v)

    if(associated(PARABS)) parabs = para
    if(associated(PARINC)) parinc = drpar + dfpar

! set the number of days per year when crossing year boundary or on restart gkw: use GEOS5/MAPL value
! -------------------------------------------------------------------------
    if(AGCM_YY .ne. year_prev) then
      dpy = get_days_per_year(AGCM_YY) ! set the number of days for current year
      year_prev = AGCM_YY
    endif

! set time step for CN model
! --------------------------

    call MAPL_GetResource ( MAPL, DTCN, Label="CATCHCN_DT:", DEFAULT=5400., RC=STATUS)
    VERIFY_(STATUS)

! CN time step over 4 hours may fail; limit to 4 hours; verify that DTCN (CATCHCN_DT) is a multiple of DT
! ------------------------------------------------------------------------------------------
    dtcn = min(dtcn,14400.)
    if(mod(dtcn,dt) /= 0) stop 'dtcn'

    ndt = get_step_size( nint(dtcn) ) ! gkw: get_step_size must be called here to set CN model time step

! sum over interval for CN
! ------------------------
    tgwm    = tgwm    + tgw
    tpm     = tpm     + tp1
    sfmcm   = sfmcm   + sfmc
    rzmm    = rzmm    + rzm
    bflowm  = bflowm  + bflow
    totwatm = totwatm + totwat
    tairm   = tairm   + TA
    psnsunm = psnsunm + psnsun*laisun
    psnsham = psnsham + psnsha*laisha
    do n = 1,N_snow
      sndzm(:) = sndzm(:) + sndzn(n,:)
    end do
    asnowm = asnowm + asnow
    cnsum   = cnsum   + 1.

! call CN model every CATCHCN_DT seconds
! --------------------------------
    if(mod(AGCM_S,nint(dtcn)) == 0) then

! compute mean state over interval
! --------------------------------
    do nz = 1,nzone
      tgwm(:,nz) = tgwm(:,nz) / cnsum(:)
      rzmm(:,nz) = rzmm(:,nz) / cnsum(:)
      do nv = 1,nveg
        psnsunm(:,nv,nz) = psnsunm(:,nv,nz) / cnsum(:)
        psnsham(:,nv,nz) = psnsham(:,nv,nz) / cnsum(:)
      end do
    end do
    tpm     = tpm     / cnsum
    sfmcm   = sfmcm   / cnsum
    bflowm  = bflowm  / cnsum
    totwatm = totwatm / cnsum
    tairm   = tairm   / cnsum
    sndzm   = sndzm   / cnsum
    asnowm  = asnowm  / cnsum

    laisun = 1.
    laisha = 1.

    call CN_Driver(istep,ntiles,nveg,nzone,dayl, &
                   tgwm,tpm,tp2,tp3,tp4,tp5,tp6, &
                   sfmcm,rzmm,wpwet, &
 	           psis,bee,poros,vgwmax,bflowm,totwatm, &
 	           tairm,psnsunm,psnsham,laisun,laisha, &
                   ityp,fveg,wtzone,sndzm,asnowm,ndep,elai,esai,tlai,totcolc,cat_id,cli_t2m, &
 		   npp,gpp,sr,nee,root,padd,vegc,xsmr,burn,fsel,closs,firefac)

! save scaled CN diagnostics
! --------------------------
    if(associated(CNLAI)) then
      cnlai(:) = 0.
      do nz = 1,nzone
        do nv = 1,nveg
          cnlai(:) = cnlai(:) + elai(:,nv,nz)*fveg(:,nv,nz)*wtzone(:,nz)
        end do
      end do
      cnlai(:) = cnlai(:) * cnsum
    endif

    if(associated(CNTLAI)) then
      cntlai(:) = 0.
      do nz = 1,nzone
        do nv = 1,nveg
          cntlai(:) = cntlai(:) + tlai(:,nv,nz)*fveg(:,nv,nz)*wtzone(:,nz)
        end do
      end do
      cntlai(:) = cntlai(:) * cnsum
    endif

    if(associated(CNSAI)) then
      cnsai(:) = 0.
      do nz = 1,nzone
        do nv = 1,nveg
          cnsai(:) = cnsai(:) + esai(:,nv,nz)*fveg(:,nv,nz)*wtzone(:,nz)
        end do
      end do
      cnsai(:) = cnsai(:) * cnsum
    endif

    if(associated(CNTOTC)) then
      cntotc(:) = 0.
      do nz = 1,nzone
        cntotc(:) = cntotc(:) + 1.e-3*totcolc(:,nz)*wtzone(:,nz)
      end do
      cntotc(:) = cntotc(:) * cnsum
    endif

    if(associated(CNVEGC)) cnvegc = 1.e-3*vegc  * cnsum
    if(associated(CNROOT)) cnroot = 1.e-3*root  * cnsum
    if(associated(CNNPP )) cnnpp  = 1.e-3*npp   * cnsum
    if(associated(CNGPP )) cngpp  = 1.e-3*gpp   * cnsum
    if(associated(CNSR  )) cnsr   = 1.e-3*sr    * cnsum
    if(associated(CNNEE )) cnnee  = 1.e-3*nee   * cnsum
    if(associated(CNXSMR)) cnxsmr = 1.e-3*xsmr  * cnsum
    if(associated(CNADD )) cnadd  = 1.e-3*padd  * cnsum
    if(associated(CNLOSS)) cnloss = 1.e-3*closs * cnsum ! total fire C loss (kg/m2/s)
    if(associated(CNBURN)) cnburn = burn        * cnsum ! area fractional fire burn rate (s-1)
    if(associated(CNFSEL)) cnfsel = fsel        * cnsum ! fire season length (days)

! reset summing arrays
! --------------------
    tgwm    = 0.
    tpm     = 0.
    sfmcm   = 0.
    rzmm    = 0.
    bflowm  = 0.
    totwatm = 0.
    tairm   = 0.
    psnsunm = 0.
    psnsham = 0.
    sndzm   = 0.
    asnowm  = 0.
    cnsum   = 0.

    else ! CN diags set to zero

    if(associated(CNLAI )) cnlai  = 0.
    if(associated(CNTLAI)) cntlai = 0.
    if(associated(CNSAI )) cnsai  = 0.
    if(associated(CNTOTC)) cntotc = 0.
    if(associated(CNVEGC)) cnvegc = 0.
    if(associated(CNROOT)) cnroot = 0.
    if(associated(CNNPP )) cnnpp  = 0.
    if(associated(CNGPP )) cngpp  = 0.
    if(associated(CNSR  )) cnsr   = 0.
    if(associated(CNNEE )) cnnee  = 0.
    if(associated(CNXSMR)) cnxsmr = 0.
    if(associated(CNADD )) cnadd  = 0.
    if(associated(CNLOSS)) cnloss = 0.
    if(associated(CNBURN)) cnburn = 0.
    if(associated(CNFSEL)) cnfsel = 0.

    endif
    
! copy CN_restart vars to catch_internal_rst gkw: only do if stopping
! ------------------------------------------
  if(NextTime == StopTime) then

    call CN_exit(ntiles,nveg,nzone,ityp,fveg,cncol,var_col,cnpft,var_pft)     
       i = 1
       do iv = 1,VAR_PFT
          do nv = 1,NUM_VEG
             do nz = 1, NUM_ZON
                do  n = 1,ntiles
                   ! to ensure unused array elements don't have crazy numbers in restart files.
                   if(fveg (n,nv,nz) == 0.) cnpft (n,i) = 0.
                end do
                i = i + 1
             end do
          end do
       end do
  endif

! update LAI for primary & secondary vegetation types
! ---------------------------------------------------
    lai1 = 0.
    wght = 0.
    do nz = 1,nzone
      do nv = 1,2
        lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
        wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
      end do
    end do
    lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

    lai2 = 0.
    wght = 0.
    do nz = 1,nzone
      do nv = 3,4
        lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
        wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
      end do
    end do
    lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

    lai = fveg1*lai1 + fveg2*lai2 ! gkw: prognostic LAI on catch_internal_rst (overwrite VEGDYN import)

! have stomatal resistance in the CN zones; map as conductance into catchment zones
! ---------------------------------------------------------------------------------
    do n = 1,ntiles

      ax1 = car1(n)
      ax2 = car2(n)
      ax4 = 1. - ax1 - ax2

      cn1 = wtzone(n,1)
      cn2 = wtzone(n,2)
      cn3 = wtzone(n,3)

! catchment: saturated area

      if(ax1 .lt. cn1) then
        f1 = ax1 ; f2 = 0. ; f3 = 0.
       else
        if(ax1 .lt. (cn1+cn2)) then
          f1 = cn1 ; f2 = ax1-cn1 ; f3 = 0.
         else
          f1 = cn1 ; f2 = cn2 ; f3 = ax1-cn1-cn2
        endif
      endif

      if(ax1 .gt. 0.) then
        rcsat(n) = ax1/(f1/rc00(n,1)+f2/rc00(n,2)+f3/rc00(n,3))
        rcxdt(n) = ax1/(f1/rcdt(n,1)+f2/rcdt(n,2)+f3/rcdt(n,3))
        rcxdq(n) = ax1/(f1/rcdq(n,1)+f2/rcdq(n,2)+f3/rcdq(n,3))
       else
        rcsat(n) = 1.e3
        rcxdt(n) = 1.e3
        rcxdq(n) = 1.e3
      endif

! compute deriviatives
      drcsdt(n) = (rcxdt(n) - rcsat(n)) / dtc
      drcsdq(n) = (rcxdq(n) - rcsat(n)) / (0.622*dea/PS(n))

! catchment: unstressed area

      if(ax1 .lt. cn1) then
        ar = ax1 + ax2
        if(ar .lt. cn1) then
          f1 = ax2 ; f2 = 0. ; f3 = 0.
         else
          if(ar .lt. (cn1+cn2)) then
            f1 = cn1-ax1 ; f2 = ar-cn1 ; f3 = 0.
           else
            f1 = cn1-ax1 ; f2 = cn2 ; f3 = ar-cn1-cn2
          endif
        endif
       else
        ar = ax2 + ax4
        if(ar .lt. cn3) then
          f1 = 0. ; f2 = 0. ; f3 = ax2
         else
          if(ax4 .gt. cn3) then
            f1 = 0. ; f2 = ax2 ; f3 = 0.
           else
            f1 = 0. ; f2 = ar-cn3 ; f3 = cn3-ax4
          endif
        endif
      endif

      if(ax2 .gt. 0.) then
        rcuns(n) = ax2/(f1/rc00(n,1)+f2/rc00(n,2)+f3/rc00(n,3))
        rcxdt(n) = ax2/(f1/rcdt(n,1)+f2/rcdt(n,2)+f3/rcdt(n,3))
        rcxdq(n) = ax2/(f1/rcdq(n,1)+f2/rcdq(n,2)+f3/rcdq(n,3))
       else
        rcuns(n) = 1.e3
        rcxdt(n) = 1.e3
        rcxdq(n) = 1.e3
      endif

! compute deriviatives
      drcudt(n) = (rcxdt(n) - rcuns(n)) / dtc
      drcudq(n) = (rcxdq(n) - rcuns(n)) / (0.622*dea/PS(n))

    end do

    if(associated(SCSAT )) scsat  = 1. / rcsat
    if(associated(SCUNS )) scuns  = 1. / rcuns

    endif ! end of check for zero tiles

! gkw: end of main CN block

! Andrea Molod (Oct 21, 2016):
 
        do N=1,NUM_SUBTILES
           DQS(:,N) = GEOS_DQSAT ( TC(:,N), PS, QSAT=QSAT(:,N),PASCALS=.true., RAMP=0.0 )
           QC (:,N) = min(max(QA(:),QSAT(:,N)),QC(:,N))
           QC (:,N) = max(min(QA(:),QSAT(:,N)),QC(:,N))
           RA (:,N) = RHO/CH(:,N)
        end do

! call unified land model
! -----------------------

        if (ntiles > 0) then

           call CATCHCN ( NTILES, LONS, LATS                         ,&
                DT	, PRECIPFRAC, cat_id, VEG1,VEG2,FVEG1,FVEG2,DZSF   ,&
                PCU     , PLS       , SNO, ICE, FRZR                 ,&
                UUU                                                  ,&

                EVSBT(:,FSAT),     DEVSBT(:,FSAT),     DEDTC(:,FSAT) ,&
                SHSBT(:,FSAT),     DHSDQA(:,FSAT),     DSHSBT(:,FSAT),&
                EVSBT(:,FTRN),     DEVSBT(:,FTRN),     DEDTC(:,FTRN) ,&
                SHSBT(:,FTRN),     DHSDQA(:,FTRN),     DSHSBT(:,FTRN),&
                EVSBT(:,FWLT),     DEVSBT(:,FWLT),     DEDTC(:,FWLT) ,&
                SHSBT(:,FWLT),     DHSDQA(:,FWLT),     DSHSBT(:,FWLT),&
                EVSBT(:,FSNW),     DEVSBT(:,FSNW),     DEDTC(:,FSNW) ,&
                SHSBT(:,FSNW),     DHSDQA(:,FSNW),     DSHSBT(:,FSNW),&

                TA           ,     QA                                ,&

                RA(:,FSAT), RA(:,FTRN), RA(:,FWLT), RA(:,FSNW)       ,&

                ZTH,  SWNETFREE, SWNETSNOW, LWDNSRF                  ,&

                PS*.01                                               ,&

                LAI0, GRN0,       SQSCAT, RSL1, RSL2, RDC            ,&

                QSAT(:,FSAT) ,    DQS(:,FSAT) ,   ALWN(:,1),    BLWN(:,1) ,&
                QSAT(:,FTRN) ,    DQS(:,FTRN) ,   ALWN(:,2),    BLWN(:,2) ,&
                QSAT(:,FWLT) ,    DQS(:,FWLT) ,   ALWN(:,3),    BLWN(:,3) ,&
                QSAT(:,FSNW) ,    DQS(:,FSNW) ,   ALWN(:,4),    BLWN(:,4) ,&

           RCSAT,DRCSDT,DRCSDQ,  RCUNS,DRCUDT,DRCUDQ,              &
                BF1, BF2, BF3, VGWMAX, CDCR1, CDCR2, PSIS	          ,&
                BEE, POROS, WPWET, COND, GNU	                  ,&
                ARS1, ARS2, ARS3, ARA1, ARA2, ARA3, ARA4	          ,&
                ARW1, ARW2, ARW3, ARW4, TSA1, TSA2, TSB1, TSB2	  ,&
                ATAU, BTAU, .false.			          ,&

                TG(:,FSAT), TG(:,FTRN), TG(:,FWLT)		          ,& 
                TC(:,FSAT), TC(:,FTRN), TC(:,FWLT)		          ,& 
                QC(:,FSAT), QC(:,FTRN), QC(:,FWLT)		          ,&

                CAPAC, CATDEF, RZEXC, SRFEXC, GHTCNT                 ,&
                WESNN, HTSNNN, SNDZN                                 ,&

                EVAPOUT, SHOUT, RUNOFF, EVPINT, EVPSOI, EVPVEG       ,&
                EVPICE                                               ,&
                BFLOW                                                ,&
                RUNSURF                                              ,&
                SMELT                                                ,&
                HLWUP                                                ,&
                SWNDSRF                                              ,&
                HLATN                                                ,&
                QINFIL                                               ,&
                AR1                                                  ,&
                AR2                                                  ,&
                RZEQ                                                 ,&
                GHFLX                                                ,&
                GHFLXSNO                                             ,&
                GHFLXTSKIN                                           ,&
                TC(:,FSNW)                                           ,&
                ASNOW                                                ,&
                TP1, TP2, TP3, TP4, TP5, TP6,  SFMC, RZMC, PRMC      ,&
                ENTOT,WTOT, WCHANGE, ECHANGE, HSNACC, EVACC, SHACC   ,&
           TSURF                                                ,&
                SHSNOW1, AVETSNOW1, WAT10CM1, WATSOI1, ICESOI1       ,&
                LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW             ,&
                TCSORIG1, TPSN1IN1, TPSN1OUT1                        ,&
                TC1_0=TC1_0, TC2_0=TC2_0, TC4_0=TC4_0                ,&
                QA1_0=QA1_0, QA2_0=QA2_0, QA4_0=QA4_0                ,&
                RCONSTIT=RCONSTIT, RMELT=RMELT, TOTDEPOS=TOTDEPOS, LHACC=LHACC)

        end if

        call MAPL_TimerOff ( MAPL, "-CATCHCN" )

        if (is_OFFLINE) then
           TC(:,FSAT) = TC1_0
           TC(:,FTRN) = TC2_0
           TC(:,FWLT) = TC4_0
           QC(:,FSAT) = QA1_0
           QC(:,FTRN) = QA2_0
           QC(:,FWLT) = QA4_0
           EVACC = 0.0
           SHACC = 0.0
        endif

        QC(:,FSNW) =  GEOS_QSAT ( TC(:,FSNW), PS, PASCALS=.true., RAMP=0.0 )

        ! --------------------------------------------------------------------------
        ! update subtile fractions
        ! --------------------------------------------------------------------------

        EMIS    = fveg1*(EMSVEG(VEG1) + (EMSBARESOIL - EMSVEG(VEG1))*exp(-LAI1)) + &
                  fveg2*(EMSVEG(VEG2) + (EMSBARESOIL - EMSVEG(VEG2))*exp(-LAI2))

        EMIS    = EMIS     *(1.-ASNOW) + EMSSNO   *ASNOW

        call MAPL_SunGetInsolation(LONS, LATS,      &
            ORBIT, ZTH, SLR, &
            INTV = TINT,     &
            currTime=CURRENT_TIME+DELT,  &
            RC=STATUS )
        VERIFY_(STATUS)

        ZTH = max(0.0,ZTH)

        ! --------------------------------------------------------------------------
        ! Update raditation exports
        ! --------------------------------------------------------------------------


        call MAPL_TimerOn(MAPL,"-ALBEDO")
        call    SIBALB(NTILES, VEG1,LAI1,GRN, ZTH,         & 
                       BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
                       ALBVR, ALBNR, ALBVF, ALBNF, MODIS_SCALE=.TRUE.  )         ! instantaneous snow-free albedos on tiles
        call STIEGLITZSNOW_CALC_TPSNOW(NTILES, HTSNNN(1,:), WESNN(1,:), TPSN1OUT1, FICE1)
        TPSN1OUT1 =  TPSN1OUT1 + Tzero        

        call   SNOW_ALBEDO(NTILES,N_snow, N_CONST_LAND4SNWALB, VEG1, LAI1, ZTH,        &
                 RHOFS,                                              &   
                 SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                & 
                 WESNN, HTSNNN, SNDZN,                               &
                 ALBVR, ALBNR, ALBVF, ALBNF,                         & ! instantaneous snow-free albedos on tiles
                 SNOVR, SNONR, SNOVF, SNONF,                         & ! instantaneous snow albedos on tiles
                 RCONSTIT, UUU, TPSN1OUT1,DRPAR, DFPAR)

        call    SIBALB(NTILES, VEG2,LAI2,GRN, ZTH,         & 
                       BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
                       ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp, MODIS_SCALE=.TRUE.  ) ! instantaneous snow-free albedos on tiles

        call   SNOW_ALBEDO(NTILES,N_snow, N_CONST_LAND4SNWALB, VEG2, LAI2, ZTH,        &
                 RHOFS,                                              &   
                 SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                & 
                 WESNN, HTSNNN, SNDZN,                               &
                 ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp, & ! instantaneous snow-free albedos on tiles
                 SNOVR_tmp, SNONR_tmp, SNOVF_tmp, SNONF_tmp, & ! instantaneous snow albedos on tiles
                 RCONSTIT, UUU, TPSN1OUT1,DRPAR, DFPAR ) 

        ALBVR(:) = ALBVR(:)*fveg1(:) + ALBVR_tmp(:)*fveg2(:)
        ALBNR(:) = ALBNR(:)*fveg1(:) + ALBNR_tmp(:)*fveg2(:)
        ALBVF(:) = ALBVF(:)*fveg1(:) + ALBVF_tmp(:)*fveg2(:)
        ALBNF(:) = ALBNF(:)*fveg1(:) + ALBNF_tmp(:)*fveg2(:)

        SNOVR(:) = SNOVR(:)*fveg1(:) + SNOVR_tmp(:)*fveg2(:)
        SNONR(:) = SNONR(:)*fveg1(:) + SNONR_tmp(:)*fveg2(:)
        SNOVF(:) = SNOVF(:)*fveg1(:) + SNOVF_tmp(:)*fveg2(:)
        SNONF(:) = SNONF(:)*fveg1(:) + SNONF_tmp(:)*fveg2(:)

        ALBVR   = ALBVR    *(1.-ASNOW) + SNOVR    *ASNOW
        ALBVF   = ALBVF    *(1.-ASNOW) + SNOVF    *ASNOW
        ALBNR   = ALBNR    *(1.-ASNOW) + SNONR    *ASNOW
        ALBNF   = ALBNF    *(1.-ASNOW) + SNONF    *ASNOW
        call MAPL_TimerOff(MAPL,"-ALBEDO")

        LWNDSRF = LWDNSRF - HLWUP

        ! --------------------------------------------------------------------------
        ! update outputs
        ! --------------------------------------------------------------------------

        DELTS = 0.0
        DELQS = 0.0

        do N=1,NUM_SUBTILES
           DELTS   = DELTS + CFT(:,N)*(TC(:,N)-TCO(:,N))*FR(:,N)
           DELQS   = DELQS + CFQ(:,N)*(QC(:,N)-QCO(:,N))*FR(:,N)
        end do

        FR(:,FSAT) =           AR1  * (1-ASNOW)
        FR(:,FTRN) =           AR2  * (1-ASNOW)
        FR(:,FWLT) = (1.0-(AR1+AR2))* (1-ASNOW)
        FR(:,FSNW) =                     ASNOW

        FR = min( max( fr,0.0 ), 1.0 )

        TST   = 0.0
        QST   = 0.0
        do N=1,NUM_SUBTILES
           TST     = TST   +           TC(:,N)          *FR(:,N)
           QST     = QST   +           QC(:,N)          *FR(:,N)
        end do

        if ( .not. is_OFFLINE) then
!amm add correction term to latent heat diagnostics (HLATN is always allocated)
!    this will impact the export LHLAND
            HLATN = HLATN - LHACC
! also add some portion of the correction term to evap from soil, int, veg and snow
           SUMEV = EVPICE+EVPSOI+EVPVEG+EVPINT
           where(SUMEV>0.)
           EVPICE = EVPICE - EVACC*EVPICE/SUMEV
           EVPSOI = EVPSOI - EVACC*EVPSOI/SUMEV
           EVPINT = EVPINT - EVACC*EVPINT/SUMEV
           EVPVEG = EVPVEG - EVACC*EVPVEG/SUMEV
           endwhere
        endif

        if(associated( LST  )) LST    = TST
        if(associated( TPSURF))TPSURF = TSURF
        if(associated( WET1 )) WET1   = max(min(SFMC / POROS,1.0),0.0)
        if(associated( WET2 )) WET2   = max(min(RZMC / POROS,1.0),0.0)
        if(associated( WET3 )) WET3   = max(min(PRMC / POROS,1.0),0.0)
        if(associated( WCSF )) WCSF   = SFMC
        if(associated( WCRZ )) WCRZ   = RZMC
        if(associated( WCPR )) WCPR   = PRMC

        if(associated( ACCUM)) ACCUM  = SLDTOT - EVPICE*(1./MAPL_ALHS) - SMELT 

        if(associated(EVPSNO)) EVPSNO = EVPICE
        if(associated(SUBLIM)) SUBLIM = EVPICE*(1./MAPL_ALHS)*FR(:,FSNW)
        if(associated(EVLAND)) EVLAND = EVAPOUT-EVACC
        if(associated(PRLAND)) PRLAND = PCU+PLS+SLDTOT
        if(associated(SNOLAND)) SNOLAND = SLDTOT
        if(associated(DRPARLAND)) DRPARLAND = DRPAR
        if(associated(DFPARLAND)) DFPARLAND = DFPAR
        if(associated(LHLAND)) LHLAND = HLATN
        if(associated(SHLAND)) SHLAND = SHOUT-SHACC
        if(associated(SWLAND)) SWLAND = SWNDSRF
        if(associated(LWLAND)) LWLAND = LWNDSRF
        if(associated(GHLAND)) GHLAND = GHFLX
        if(associated(GHSNOW)) GHSNOW = GHFLXSNO
        if(associated(SHSNOW)) SHSNOW = SHSNOW1                   
        if(associated(AVETSNOW)) AVETSNOW = AVETSNOW1             
        if(associated(WAT10CM)) WAT10CM = WAT10CM1                
        if(associated(WATSOI)) WATSOI = WATSOI1                   
        if(associated(ICESOI)) ICESOI = ICESOI1                   
        if(associated(LHSNOW)) LHSNOW = LHSNOW1                   
        if(associated(LWUPSNOW)) LWUPSNOW = LWUPSNOW1             
        if(associated(LWDNSNOW)) LWDNSNOW = LWDNSNOW1             
        if(associated(SWNETSNOW1)) SWNETSNOW1 = NETSWSNOW         
        if(associated(TCSORIG)) TCSORIG = TCSORIG1                
        if(associated(TPSN1IN)) TPSN1IN = TPSN1IN1                
        if(associated(TPSN1OUT)) TPSN1OUT = TPSN1OUT1
        if(associated(GHTSKIN))GHTSKIN = GHFLXTSKIN
        if(associated(SMLAND)) SMLAND = SMELT
        if(associated(TWLAND)) TWLAND = WTOT
        if(associated(TELAND)) TELAND = ENTOT
        if(associated(TSLAND)) TSLAND = WESNN (1,:) + WESNN (2,:) + WESNN (3,:)
        if(associated(DWLAND)) DWLAND = WCHANGE
        if(associated(DHLAND)) DHLAND = ECHANGE
        if(associated(SPLAND)) SPLAND = SHACC
        if(associated(SPWATR)) SPWATR = EVACC
        if(associated(SPSNOW)) SPSNOW = HSNACC

        if(associated(FRSAT )) FRSAT  = max( min( FR(:,FSAT),1.0 ), 0.0 )
        if(associated(FRUST )) FRUST  = max( min( FR(:,FTRN),1.0 ), 0.0 )
        if(associated(FRWLT )) FRWLT  = max( min( FR(:,FWLT),1.0 ), 0.0 )

        if(associated(SNOMAS)) SNOMAS = WESNN (1,:) + WESNN (2,:) + WESNN (3,:)
        if(associated(SNOWDP)) SNOWDP = SNDZN (1,:) + SNDZN (2,:) + SNDZN (3,:)

        if(associated(RMELTDU001)) RMELTDU001 = RMELT(:,1) 
        if(associated(RMELTDU002)) RMELTDU002 = RMELT(:,2) 
        if(associated(RMELTDU003)) RMELTDU003 = RMELT(:,3) 
        if(associated(RMELTDU004)) RMELTDU004 = RMELT(:,4) 
        if(associated(RMELTDU005)) RMELTDU005 = RMELT(:,5) 
        if(associated(RMELTBC001)) RMELTBC001 = RMELT(:,6) 
        if(associated(RMELTBC002)) RMELTBC002 = RMELT(:,7) 
        if(associated(RMELTOC001)) RMELTOC001 = RMELT(:,8) 
        if(associated(RMELTOC002)) RMELTOC002 = RMELT(:,9) 

        if(associated(TPSN1OUT)) then
           where(WESNN(1,:)>0.)
            TPSN1OUT = TPSN1OUT1                
           elsewhere
            TPSN1OUT = MAPL_UNDEF
           end where              
        end if

        if(associated(TPSN1)) then
           where(WESNN(1,:)>0.)
              TPSN1  = TC(:,FSNW)
           elsewhere
              TPSN1  = MAPL_UNDEF
           end where
        end if

        if(associated(TPSAT)) then
           where(FR(:,FSAT)>0.)
              TPSAT  = TC(:,FSAT)
           elsewhere
              TPSAT  = MAPL_UNDEF
           end where
        end if

        if(associated(TPWLT)) then
           where(FR(:,FWLT)>0.)
              TPWLT  = TC(:,FWLT)
           elsewhere
              TPWLT  = MAPL_UNDEF
           end where
        end if

        if(associated(TPUST)) then
           where(FR(:,FTRN)>0.)
              TPUST  = TC(:,FTRN)
           elsewhere
              TPUST  = MAPL_UNDEF
           end where
        end if


        ! --------------------------------------------------------------------------
        ! update internal state arrays
        ! --------------------------------------------------------------------------

        GHTCNT1 = GHTCNT(1,:)
        GHTCNT2 = GHTCNT(2,:)
        GHTCNT3 = GHTCNT(3,:)
        GHTCNT4 = GHTCNT(4,:)
        GHTCNT5 = GHTCNT(5,:)
        GHTCNT6 = GHTCNT(6,:)

        WESNN1  = WESNN (1,:)
        WESNN2  = WESNN (2,:)
        WESNN3  = WESNN (3,:)

        HTSNNN1 = HTSNNN(1,:)
        HTSNNN2 = HTSNNN(2,:)
        HTSNNN3 = HTSNNN(3,:)

        SNDZN1  = SNDZN (1,:)
        SNDZN2  = SNDZN (2,:)
        SNDZN3  = SNDZN (3,:)

        if (DO_GOSWIM /= 0) then
           RDU001(:,:) = RCONSTIT(:,:,1) 
           RDU002(:,:) = RCONSTIT(:,:,2) 
           RDU003(:,:) = RCONSTIT(:,:,3) 
           RDU004(:,:) = RCONSTIT(:,:,4) 
           RDU005(:,:) = RCONSTIT(:,:,5) 
           RBC001(:,:) = RCONSTIT(:,:,6) 
           RBC002(:,:) = RCONSTIT(:,:,7) 
           ROC001(:,:) = RCONSTIT(:,:,8) 
           ROC002(:,:) = RCONSTIT(:,:,9) 
        end if

        ! --------------------------------------------------------------------------

        deallocate ( wght )
        deallocate ( lai1 )
        deallocate ( lai2 )
        deallocate ( ALBVR_tmp )
        deallocate ( ALBNR_tmp )
        deallocate ( ALBVF_tmp )
        deallocate ( ALBNF_tmp )
        deallocate ( SNOVR_tmp )
        deallocate ( SNONR_tmp )
        deallocate ( SNOVF_tmp )
        deallocate ( SNONF_tmp )

        deallocate(GHTCNT   )
        deallocate(WESNN    )
        deallocate(HTSNNN   )
        deallocate(SNDZN    )
	deallocate(TILEZERO )
        deallocate(DZSF     )
	deallocate(SWNETFREE)
	deallocate(SWNETSNOW)
	deallocate(VEG1     )
	deallocate(VEG2     )
	deallocate(RCSAT    )
	deallocate(DRCSDT   )
	deallocate(DRCSDQ   )
	deallocate(RCUNS    )
	deallocate(DRCUDT   )
	deallocate(DRCUDQ   )
	deallocate(ZTH      )
	deallocate(SLR      )
	deallocate(RSL1     )
	deallocate(RSL2     )
	deallocate(SQSCAT   )
	deallocate(RDC      )
        deallocate(UUU      )
        deallocate(RHO      )
        deallocate(ZVG      )
        deallocate(LAI0     )
        deallocate(GRN0     )
        deallocate(Z0       )
	deallocate(D0       )
	deallocate(SFMC     )
	deallocate(RZMC     )
        deallocate(PRMC     )
        deallocate(ENTOT    )
        deallocate(WTOT     )
        deallocate(GHFLXSNO )
        deallocate(SHSNOW1  )
        deallocate(AVETSNOW1)
        deallocate(WAT10CM1 )
        deallocate(WATSOI1  )
        deallocate(ICESOI1  )
        deallocate(LHSNOW1  )
        deallocate(LWUPSNOW1)
        deallocate(LWDNSNOW1)
        deallocate(NETSWSNOW)
        deallocate(TCSORIG1 )
        deallocate(LHACC )
        deallocate(SUMEV )
        deallocate(TPSN1IN1 )
        deallocate(TPSN1OUT1)
        deallocate(GHFLXTSKIN)
        deallocate(WCHANGE  )
        deallocate(ECHANGE  )
        deallocate(HSNACC   )
        deallocate(EVACC    )
        deallocate(SHACC    )
        deallocate(VSUVR    )
        deallocate(VSUVF    )
        deallocate(SNOVR    )
        deallocate(SNOVF    )
        deallocate(SNONR    )
        deallocate(SNONF    )
        deallocate(SHSBT    )
        deallocate(DSHSBT   )
        deallocate(EVSBT    )
        deallocate(DEVSBT   )
        deallocate(DEDTC    )
        deallocate(DHSDQA   )
        deallocate(CFT      )
        deallocate(CFQ      )
        deallocate(TCO      )
        deallocate(QCO      )
        deallocate(DQS      )
        deallocate(QSAT     )
        deallocate(RA       )
        deallocate(CAT_ID   )
        deallocate(ALWX     )
        deallocate(BLWX     )
        deallocate(ALWN     )
        deallocate(BLWN     )
        deallocate(TC1_0    )
        deallocate(TC2_0    )
        deallocate(TC4_0    )
        deallocate(QA1_0    )
        deallocate(QA2_0    )
        deallocate(QA4_0    )
        deallocate(fveg1    )
        deallocate(fveg2    )
        deallocate(RCONSTIT )
        deallocate(TOTDEPOS )
        deallocate(RMELT    )
        deallocate(FICE1    )
        deallocate(SLDTOT   )
        deallocate(   btran )
	deallocate(     wgt )
	deallocate(     bt1 )
	deallocate(     bt2 )
	deallocate(     bt4 )
	deallocate(     wpp )
	deallocate(    fwet )
	deallocate(     sm1 )
	deallocate(     sm2 )
	deallocate(     sm4 )
	deallocate(  btran1 )
	deallocate(  btran2 )
	deallocate(  btran3 )
	deallocate(     tcx )
	deallocate(     qax )
	deallocate(     rcx )
	deallocate(   rcxdt )
	deallocate(   rcxdq )
	deallocate(     tx1 )
	deallocate(     tx2 )
	deallocate(     tx3 )
	deallocate(     qx1 )
	deallocate(     qx2 )
	deallocate(     qx3 )
	deallocate(    car1 )
	deallocate(    car2 )
	deallocate(    car4 )
	deallocate( parzone )
	deallocate(    para )
	deallocate(  totwat )
	deallocate(     npp )
	deallocate(     gpp )
	deallocate(      sr )
	deallocate(     nee )
	deallocate(    padd )
	deallocate(    root )
	deallocate(    vegc )
	deallocate(    xsmr )
	deallocate(    burn )
	deallocate(    fsel )
	deallocate(   closs )
	deallocate(    dayl )
	deallocate(dayl_fac )

	deallocate(     tgw )
	deallocate(     rzm )
	deallocate(    rc00 )
	deallocate(    rcdt )
	deallocate(    rcdq )
	deallocate( totcolc )
	deallocate(  wtzone )

	deallocate( psnsunx )
	deallocate( psnshax )
	deallocate( sifsunx )
	deallocate( sifshax )
	deallocate( laisunx )
	deallocate( laishax )
	deallocate(    elaz )
	deallocate(    esaz )
	deallocate(    fvez )
	deallocate(    ityz )

	deallocate(   elai )
	deallocate(   esai )
	deallocate(   fveg )
	deallocate(   tlai )
	deallocate( psnsun )
	deallocate( psnsha )
	deallocate( laisun )
	deallocate( laisha )
	deallocate(   ityp )

	deallocate(      ht )
	deallocate(      tp )
	deallocate( soilice )

        RETURN_(ESMF_SUCCESS)

      end subroutine Driver

end subroutine RUN2

!BOP
! !IROUTINE: RUN0 -- Extra run method for the OFFLINE case, called by RUN1
! !INTERFACE:

subroutine RUN0(gc, import, export, clock, rc)
  
  ! !ARGUMENTS:
  type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
  type(ESMF_State),    intent(inout) :: import ! Import state
  type(ESMF_State),    intent(inout) :: export ! Export state
  type(ESMF_Clock),    intent(inout) :: clock  ! The clock
  integer, optional,   intent(  out) :: rc     ! Error code

  ! !DESCRIPTION: In the OFFLINE case, some diagnostic vars (INTERNAL states
  ! asnow and emis) are updated here.
  !EOP

  ! ErrLog variables
  integer :: status
  character(len=ESMF_MAXSTR) :: Iam
  character(len=ESMF_MAXSTR) :: comp_name

  ! Local variables

  !! ESMF/MAPL variables
  type(MAPL_MetaComp), pointer :: MAPL
  type(ESMF_State) :: INTERNAL

  !! IMPORT pointers
  real, pointer :: ps (:)=>null()

  !! INTERNAL pointers
  !! -asnow-emis-ww-fr-
  real, pointer :: ity(:,:)=>null()
  real, pointer :: fvg(:,:)=>null()
  real, pointer :: asnow(:)=>null()
  real, pointer :: emis(:)=>null()
  real, pointer :: ww(:,:)=>null()
  real, pointer :: fr(:,:)=>null()
  real, pointer :: DCQ(:,:)=>null()
  real, pointer :: DCH(:,:)=>null()
  !! -prognostic-variables-
  real, pointer :: tc(:,:)=>null()
  real, pointer :: qc(:,:)=>null()
  real, pointer :: htsnnn1(:)=>null()
  real, pointer :: wesnn1(:)=>null()
  real, pointer :: wesnn2(:)=>null()
  real, pointer :: wesnn3(:)=>null()
  real, pointer :: srfexc(:)=>null()
  real, pointer :: rzexc(:)=>null()
  real, pointer :: catdef(:)=>null()
  !! -parameters-
  real, pointer :: vgwmax(:)=>null()
  real, pointer :: cdcr1(:)=>null()
  real, pointer :: cdcr2(:)=>null()
  real, pointer :: psis(:)=>null()
  real, pointer :: bee(:)=>null()
  real, pointer :: poros(:)=>null()
  real, pointer :: wpwet(:)=>null()
  real, pointer :: ars1(:)=>null()
  real, pointer :: ars2(:)=>null()
  real, pointer :: ars3(:)=>null()
  real, pointer :: ara1(:)=>null()
  real, pointer :: ara2(:)=>null()
  real, pointer :: ara3(:)=>null()
  real, pointer :: ara4(:)=>null()
  real, pointer :: arw1(:)=>null()
  real, pointer :: arw2(:)=>null()
  real, pointer :: arw3(:)=>null()
  real, pointer :: arw4(:)=>null()

  !! Miscellaneous
  integer :: ntiles, nv, nz
  real, allocatable :: dummy(:)
  real :: SURFLAY
  real, allocatable :: dzsf(:), ar1(:), ar2(:), wesnn(:,:)
  real, allocatable :: catdefcp(:), srfexccp(:), rzexccp(:)
  real, allocatable :: VEG1(:), VEG2(:)
  integer, allocatable :: ityp(:,:,:)
  real,    allocatable :: fveg(:,:,:), elai(:,:,:), esai(:,:,:), wtzone(:,:), lai1(:), lai2(:), wght(:)
  real, allocatable,dimension(:) :: fveg1, fveg2

  ! Begin...

  ! Get component name and setup traceback handle
  call ESMF_GridCompGet(gc, name=comp_name, rc=status)
  VERIFY_(status)
  Iam = trim(comp_name)//"::RUN0"

  ! Get MAPL object
  call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
  VERIFY_(status)

  ! Get component's internal ESMF state
  call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, rc=status)
  VERIFY_(status)

  ! Pointers to IMPORTs
  call MAPL_GetPointer(import, ps, 'PS', rc=status)
  VERIFY_(status)

  ! Pointers to INTERNALs
  call MAPL_GetPointer(INTERNAL, ITY, 'ITY', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, FVG, 'FVG', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, asnow, 'ASNOW', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, emis, 'EMIS', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, fr, 'FR', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ww, 'WW', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, DCQ, 'DCQ', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, DCH, 'DCH', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, tc, 'TC', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, qc, 'QC', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, htsnnn1, 'HTSNNN1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, wesnn1, 'WESNN1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, wesnn2, 'WESNN2', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, wesnn3, 'WESNN3', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, vgwmax, 'VGWMAX', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, cdcr1, 'CDCR1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, cdcr2, 'CDCR2', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, psis, 'PSIS', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, bee, 'BEE', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, poros, 'POROS', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, wpwet, 'WPWET', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ars1, 'ARS1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ars2, 'ARS2', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ars3, 'ARS3', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ara1, 'ARA1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ara2, 'ARA2', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ara3, 'ARA3', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ara4, 'ARA4', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, arw1, 'ARW1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, arw2, 'ARW2', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, arw3, 'ARW3', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, arw4, 'ARW4', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, srfexc, 'SRFEXC', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, rzexc, 'RZEXC', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, catdef, 'CATDEF', rc=status)
  VERIFY_(status)

  ! Number of tiles and a dummy real array
  ntiles = size(HTSNNN1)
  allocate(dummy(ntiles), stat=status)
  VERIFY_(status)
  ! Reset WW
  WW = 0.

  ! get CNLAI to compute emmissivity
  allocate(fveg1     (NTILES))
  allocate(fveg2     (NTILES))
  allocate(veg1(ntiles), stat=status)
  VERIFY_(status)
  allocate(veg2(ntiles), stat=status)
  VERIFY_(status)
  allocate(   ityp(ntiles,num_veg,num_zon) )
  allocate(   fveg(ntiles,num_veg,num_zon) )
  allocate( wtzone(ntiles,num_zon) )
  allocate(   elai(ntiles,num_veg,num_zon) )
  allocate(   esai(ntiles,num_veg,num_zon) )

  allocate ( lai1(ntiles) )
  allocate ( lai2(ntiles) )
  allocate ( wght(ntiles) )

! set CLM CN PFT & fraction, set carbon zone weights
! --------------------------------------------------
   do nz = 1,num_zon
     ityp(:,:,nz) = nint(ity(:,:))
     fveg(:,:,nz) = fvg(:,:)
     wtzone(:,nz) = CN_zone_weight(nz)
   end do

! obtain LAI from previous time step (from CN model)
! --------------------------------------------------

   call get_CN_LAI(ntiles,num_veg,num_zon,ityp,fveg,elai,esai=esai)

   lai1 = 0.
   wght = 0.
   do nz = 1,num_zon
     do nv = 1,2
       lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
       wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
     end do
   end do
   lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

   lai2 = 0.
   wght = 0.
   do nz = 1,num_zon
     do nv = 3,4
       lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
       wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
     end do
   end do
   lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

   deallocate ( ityp )
   deallocate ( fveg )
   deallocate ( elai )
   deallocate ( esai )
   deallocate ( wtzone )

!  Vegetation types used to index into tables
!--------------------------------------------

   where(ITY(:,1) > 0.)
     VEG1 = map_cat(nint(ITY(:,1)))  ! gkw: primary   CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
    elsewhere
     VEG1 = map_cat(nint(ITY(:,2)))  ! gkw: primary   CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   where(ITY(:,3) > 0.)
     VEG2 = map_cat(nint(ITY(:,3)))  ! gkw: secondary CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
    elsewhere
     VEG2 = map_cat(nint(ITY(:,4)))  ! gkw: secondary CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   ASSERT_((count(VEG1>NTYPS.or.VEG1<1)==0)) 
   ASSERT_((count(VEG2>NTYPS.or.VEG2<1)==0))
   fveg1(:) = fvg(:,1) + fvg(:,2) ! sum veg fractions (primary)   gkw: NUM_VEG specific
   fveg2(:) = fvg(:,3) + fvg(:,4) ! sum veg fractions (secondary) gkw: fveg1+fveg2=1

  ! Compute ASNOW and EMIS
  allocate(wesnn(3,ntiles), stat=status)
  VERIFY_(status)
  wesnn(1,:) = wesnn1
  wesnn(2,:) = wesnn2
  wesnn(3,:) = wesnn3
  call StieglitzSnow_calc_asnow(3, ntiles, wesnn, asnow)

  EMIS    = fveg1*(EMSVEG(NINT(VEG1)) + (EMSBARESOIL - EMSVEG(NINT(VEG1)))*exp(-LAI1)) + &
       fveg2*(EMSVEG(NINT(VEG2)) + (EMSBARESOIL - EMSVEG(NINT(VEG2)))*exp(-LAI2))

  emis = emis*(1.-asnow) + EMSSNO*asnow

  ! Compute FR
  ! Step 1: set dzsf
  ! Step 2: compute ar1, ar2 via call to catch_calc_soil_moist()
  ! Step 3: compute fr

  ! -step-1-
  call MAPL_GetResource(MAPL, SURFLAY, Label="SURFLAY:", DEFAULT=50.0, rc=status)
  VERIFY_(status)
  allocate(dzsf(ntiles), stat=status)
  VERIFY_(status)
  dzsf = SURFLAY

  ! -step-2-
  allocate(ar1(ntiles), stat=status)
  VERIFY_(status)
  allocate(ar2(ntiles), stat=status)
  VERIFY_(status)
  ! -we-don't-want-to-modify-srfexc-rzexc-and-catdef-
  ! -so-we-create-local-copies-
  allocate(catdefcp(ntiles), stat=status)
  VERIFY_(status)
  allocate(srfexccp(ntiles), stat=status)
  VERIFY_(status)
  allocate(rzexccp(ntiles), stat=status)
  VERIFY_(status)
  catdefcp = catdef
  srfexccp = srfexc
  rzexccp = rzexc
  call catch_calc_soil_moist(                                                   &
       ! intent(in)
       ntiles, nint(veg1), dzsf, vgwmax, cdcr1, cdcr2,                           &
       psis, bee, poros, wpwet,                                                 &
       ars1, ars2, ars3,                                                        &
       ara1, ara2, ara3, ara4,                                                  &
       arw1, arw2, arw3, arw4,                                                  &
       ! intent(inout)
       ! from process_cat
       srfexccp, rzexccp, catdefcp,                                             &
       ! use this one can match process_cat
       ! srfexc, rzexc, catdef,                                             &
       ! intent(out)
       ar1, ar2, dummy                                                          &
       )

  ! -step-3-
  fr(:,FSAT) =           ar1  * (1-asnow)
  fr(:,FTRN) =           ar2  * (1-asnow)
  fr(:,FWLT) = (1.0-(ar1+ar2))* (1-asnow)
  fr(:,FSNW) =                     asnow
  fr = min(max(fr,0.0),1.0)

  ! Overwrite the top layer snow temperature tc(4) with its diagnosed value
  call StieglitzSnow_calc_tpsnow(ntiles, htsnnn1, wesnn1, tc(:,4), dummy)
  tc(:,FSNW) = tc(:,FSNW) + MAPL_TICE ! Convert to K

  ! Overwrite qc(4)
  !qc(:,FSNW) = GEOS_QSAT(tc(:,FSNW), PS, PASCALS=.true., RAMP=0.0)
   qc(:,FSNW) = MAPL_EQsat(tc(:,FSNW),PS,OverIce=.true.)

  ! Clean up
  if (allocated(catdefcp)) deallocate(catdefcp)
  if (allocated(srfexccp)) deallocate(srfexccp)
  if (allocated(rzexccp)) deallocate(rzexccp)
  if (allocated(dummy)) deallocate(dummy)
  if (allocated(dzsf)) deallocate(dzsf)
  if (allocated(ar1)) deallocate(ar1)
  if (allocated(ar2)) deallocate(ar2)
  if (allocated(wesnn)) deallocate(wesnn)
  if (allocated(fveg1)) deallocate (fveg1)
  if (allocated(fveg2)) deallocate (fveg2)
  if (allocated(veg1)) deallocate (veg1)
  if (allocated(veg2)) deallocate (veg2)
  if (allocated(lai1)) deallocate (lai1)
  if (allocated(lai2)) deallocate (lai2)
  if (allocated(wght)) deallocate (wght)

  ! All done
  RETURN_(ESMF_SUCCESS)

end subroutine RUN0

end module GEOS_CatchCNGridCompMod

