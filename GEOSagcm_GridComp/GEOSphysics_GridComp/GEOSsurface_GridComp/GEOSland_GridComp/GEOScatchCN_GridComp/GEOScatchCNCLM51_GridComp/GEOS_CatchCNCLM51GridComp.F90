!  $Id$ 

#include "MAPL_Generic.h"
#define DEALLOC_(A) if(associated(A))then;A=0;if(MAPL_ShmInitialized)then; call MAPL_DeAllocNodeArray(A,rc=STATUS);else; deallocate(A,stat=STATUS);endif;_VERIFY(STATUS);NULLIFY(A);endif

!=============================================================================
module GEOS_CatchCNCLM51GridCompMod

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

  use, intrinsic :: iso_fortran_env, only: INT64
  use sfclayer  ! using module that contains sfc layer code
  use ESMF
  use GEOS_Mod
  use GEOS_UtilsMod
  use DragCoefficientsMod
  use CATCHMENT_CN_MODEL
  use CNCLM_DriverMod
  use CNCLM_Photosynthesis
  use CNCLM_initMod
  USE STIEGLITZSNOW,   ONLY :                  &
       StieglitzSnow_snow_albedo,              &
       StieglitzSnow_calc_tpsnow,              &
       N_CONSTIT,                              &
       NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD, &
       NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD, &
       NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD, &
       NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD, &
       NUM_SSDP, NUM_SSSV, NUM_SSWT, NUM_SSSD, &
       StieglitzSnow_calc_asnow

  USE CATCH_CONSTANTS, ONLY :                  &
       N_SNOW         => CATCH_N_SNOW,         &
       N_GT           => CATCH_N_GT,           &
       DZGT           => CATCH_DZGT,           &
       DZTSURF        => CATCH_DZTSURF,        &
       RHOFS          => CATCH_SNOW_RHOFS,     &
       SNWALB_VISMAX  => CATCH_SNOW_VISMAX,    &
       SNWALB_NIRMAX  => CATCH_SNOW_NIRMAX,    &
       SLOPE          => CATCH_SNOW_SLOPE,     &
       PEATCLSM_POROS_THRESHOLD

  USE  clm_varpar, ONLY :                      &
       NUM_ZON, NUM_VEG, VAR_COL, VAR_PFT,     &
       CN_zone_weight, map_cat, numpft
 
  USE MAPL
  use MAPL_ConstantsMod,only: Tzero => MAPL_TICE, pi => MAPL_PI, MAPL_RHOWTR 
  use clm_time_manager, only: get_days_per_year, get_step_size, get_nstep, is_first_step
  use pftconMod,        only: noveg
  use lsm_routines,     only: sibalb, catch_calc_soil_moist, catch_calc_peatclsm_waterlevel, catch_calc_zbar, gndtmp

  use update_model_para4cn, only : upd_curr_date_time
  use WaterType
  use CNVegetationFacade
  use catch_wrap_stateMod

implicit none
private

  include "netcdf.inc"

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

integer,parameter :: NTYPS = MAPL_NUMVEGTYPES

real,   parameter :: SAI4ZVG(NTYPS) = (/ 0.60653, 0.60653, 0.60653, 1.0, 1.0, 1.0 /)
real,   parameter :: HPBL           = 1000.
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
    type(T_CATCHCN_STATE), pointer :: CATCHCN_INTERNAL
    class(T_CATCHCN_STATE), pointer :: statePtr
    type(CATCHCN_WRAP) :: wrap
    integer :: OFFLINE_MODE
    integer :: RESTART

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

    allocate(CATCHCN_INTERNAL, stat=status)
    VERIFY_(status)
    statePtr => CATCHCN_INTERNAL

    ! resource variables for offline GEOSldas; for documentation, see GEOSldas/src/Applications/LDAS_App/GEOSldas_LDAS.rc
    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)
    call MAPL_GetResource ( MAPL, CATCHCN_INTERNAL%CATCH_OFFLINE, Label="CATCHMENT_OFFLINE:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, CATCHCN_INTERNAL%CATCH_SPINUP,  Label="CATCHMENT_SPINUP:",  DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    OFFLINE_MODE = CATCHCN_INTERNAL%CATCH_OFFLINE    ! shorthand

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
         LONG_NAME          = 'surface_downwelling_PAR_beam_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'DRPAR'                       ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
         RC=STATUS  ) 
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'surface_downwelling_PAR_diffuse_flux',&
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
         LONG_NAME          = 'surface_absorbed_longwave_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'LWDNSRF'                     ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux',&
         UNITS              = 'W m-2'                       ,&
         SHORT_NAME         = 'ALW'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,&
         LONG_NAME          = 'linearization_of_surface_emitted_longwave_flux',&
         UNITS              = 'W_m-2 K-1'                   ,&
         SHORT_NAME         = 'BLW'                         ,&
         DIMS               = MAPL_DimsTileOnly             ,&
         VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
    VERIFY_(STATUS)

    IF (catchcn_internal%ATM_CO2 == 4) THEN
       call MAPL_AddImportSpec(GC,                                    &
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
         LONG_NAME          = 'vegetation_greenness_fraction'           ,&
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
  if      ( OFFLINE_MODE == 1 ) then
     RESTART = MAPL_RestartSkip
  elseif  ( OFFLINE_MODE == 2 ) then
     RESTART = MAPL_RestartOptional
  elseif  ( OFFLINE_MODE == 0 ) then
     RESTART = MAPL_RestartRequired
  else
     ASSERT_(.FALSE.)
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
    LONG_NAME          = 'maximum_soil_water_content_above_wilting_point'         ,&
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
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'POROS'                     ,&
    FRIENDLYTO         = trim(COMP_NAME)             ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                  ,&
    LONG_NAME          = 'soil_wilting_point_in_degree_of_saturation_units'  ,&
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
    LONG_NAME          = 'surface_layer_soil_temperature',&
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
    LONG_NAME          = 'vegetation_interception_water_storage',&
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
    RESTART            = RESTART                     ,&
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

  if (catchcn_internal%SNOW_ALBEDO_INFO == 1) then
    call MAPL_AddInternalSpec(GC                  ,&
       LONG_NAME          = 'effective_snow_reflectivity',&
       UNITS              = '1'                         ,&
       SHORT_NAME         = 'SNOWALB'                   ,&
       FRIENDLYTO         = trim(COMP_NAME)             ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  )
     VERIFY_(STATUS)
  endif

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
    LONG_NAME          = 'surface_moisture_exchange_coefficient',&
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

  if     (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND == 1) then
     
     ! for *analytical* extra derivatives in louissurface
     
     call MAPL_AddInternalSpec(GC,                          &
          SHORT_NAME         = 'delCH_delTVA',              &
          LONG_NAME          = 'partial_derivative_of_CH_wrt_virtual_Tair', & 
          UNITS              = '1',                         &
          DIMS               = MAPL_DimsTileTile,           &
          NUM_SUBTILES       = NUM_SUBTILES                ,&
          VLOCATION          = MAPL_VLocationNone,          &
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  )
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC,                          &
          SHORT_NAME         = 'delCQ_delTVA',              &  
          LONG_NAME          = 'partial_derivative_of_CQ_wrt_virtual_Tair', &
          UNITS              = '1',                         &
          DIMS               = MAPL_DimsTileTile,           &
          NUM_SUBTILES       = NUM_SUBTILES                ,&
          VLOCATION          = MAPL_VLocationNone,          &
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  )
     VERIFY_(STATUS)

  elseif (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND >= 2) then 
     
     ! for *numerical* extra derivatives in helfsurface and louissurface
     
     call MAPL_AddInternalSpec(GC                          ,&
          LONG_NAME          = 'partial_derivative_of_CH_wrt_canopy_temperature', &
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'delCH_delTC'               ,&
          DIMS               = MAPL_DimsTileTile           ,&
          NUM_SUBTILES       = NUM_SUBTILES                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  )
     VERIFY_(STATUS)
     
     call MAPL_AddInternalSpec(GC                          ,&
          LONG_NAME          = 'partial_derivative_of_CQ_wrt_canopy_specific_humidity', &
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'delCQ_delQC'               ,&
          DIMS               = MAPL_DimsTileTile           ,&
          NUM_SUBTILES       = NUM_SUBTILES                ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RESTART            = MAPL_RestartSkip            ,&
          RC=STATUS  )
     VERIFY_(STATUS)

  end if

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
    LONG_NAME          = 'CLM_peak_month_agricultural_fire',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ABM'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CLM_peatland_fraction'     ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'PEATF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CLM_gross_domestic_product',&
    UNITS              = 'K 1995US$/capita'          ,&
    SHORT_NAME         = 'GDP'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'CLM_human_density_2010'    ,&
    UNITS              = 'individual/km2'            ,&
    SHORT_NAME         = 'HDM'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
    RESTART            = MAPL_RestartRequired        ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
    LONG_NAME          = 'field_capacity'            ,&
    UNITS              = 'm3/m3'                     ,&
    SHORT_NAME         = 'FIELDCAP'                  ,&
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
       LONG_NAME          = 'CN sum for surface layer soil temperature',&
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
       SHORT_NAME         = 'SFMM'                      ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       UNGRIDDED_DIMS     = (/NUM_ZON/)                 ,&
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
       LONG_NAME          = 'CN sum for relative humidity',&
       UNITS              = '%'                         ,&
       SHORT_NAME         = 'RHM'                       ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
     
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for wind speed'     ,&
       UNITS              = 'm s-1'                     ,&
       SHORT_NAME         = 'WINDM'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for rainfall'       ,&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RAINFM'                    ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for snow fall'      ,&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'SNOWFM'                    ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for surface runoff' ,&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RUNSRFM'                   ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for frac saturated area',&
       UNITS              = '1'                         ,&
       SHORT_NAME         = 'AR1M'                      ,&
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
       LONG_NAME          = 'CN sum for sunlit leaf maintenance respiration',&
       UNITS              = 'umol CO2 m-2 s-1'                              ,&
       SHORT_NAME         = 'LMRSUNM'                   ,&
       DIMS               = MAPL_DimsTileOnly           ,& 
       VLOCATION          = MAPL_VLocationNone          ,&
       UNGRIDDED_DIMS     = (/NUM_VEG,NUM_ZON/)         ,& 
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for shaded leaf maintenance respiration',&
       UNITS              = 'umol CO2 m-2 s-1'                   ,&
       SHORT_NAME         = 'LMRSHAM'                  ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       UNGRIDDED_DIMS     = (/NUM_VEG,NUM_ZON/)         ,&
       RESTART            = MAPL_RestartOptional        ,&  
       RC=STATUS  )         
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for sunlit leaf area index',&
       UNITS              = '1'                              ,&
       SHORT_NAME         = 'LAISUNM'                   ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       UNGRIDDED_DIMS     = (/NUM_VEG,NUM_ZON/)         ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = 'CN sum for shaded leaf area index',&
       UNITS              = '1'                   ,&
       SHORT_NAME         = 'LAISHAM'                  ,&
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
       LONG_NAME          = '5-day running mean of CN sum for snow depth',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'SNDZM5D'                   ,&
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

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '10-day running mean of 2-m temperature',&
       UNITS              = 'K'                         ,&
       SHORT_NAME         = 'T2M10D'                    ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '10-day running mean of surface layer soil temperature',&
       UNITS              = 'K'                         ,&
       SHORT_NAME         = 'TG10D'                    ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '5-day running mean of daily minimum 2-m temperature',&
       UNITS              = 'K'                         ,&
       SHORT_NAME         = 'T2MMIN5D'                    ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '30-day running mean of surface relative humidity',&
       UNITS              = '%'                         ,&
       SHORT_NAME         = 'RH30D'                    ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '10-day running mean of total precipitation',&
       UNITS              = 'mm H2O/s'                         ,&
       SHORT_NAME         = 'TPREC10D'                  ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '60-day running mean of total precipitation',&
       UNITS              = 'mm H2O/s'                         ,&
       SHORT_NAME         = 'TPREC60D'                  ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddInternalSpec(GC                       ,&
       LONG_NAME          = '365-day running mean of total ET',&
       UNITS              = 'W m-2'                         ,&
       SHORT_NAME         = 'ET365D'                  ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)


  call MAPL_AddInternalSpec(GC,                    &
       LONG_NAME          = 'overland_runoff_including_throughflow'  ,&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'RUNSURF'                   ,&
       FRIENDLYTO         = trim(COMP_NAME)             ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartOptional        ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
    
  !---------- GOSWIM snow impurity related variables ----------

  if (catchcn_internal%N_CONST_LAND4SNWALB /= 0) then 
  
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

!EOS

  !  EXPORT STATE:
  
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
    LONG_NAME          = 'runoff_total_flux'         ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'RUNOFF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'interception_loss_latent_heat_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPINT'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baresoil_evaporation_latent_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSOI'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'transpiration_latent_heat_flux' ,&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPVEG'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'snowpack_evaporation_latent_heat_flux',&
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
    LONG_NAME          = 'total_soil_moisture'      ,&
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
    LONG_NAME          = 'snowpack_evaporation_latent_heat_flux',&
    UNITS              = 'W m-2'                     ,&
    SHORT_NAME         = 'EVPSNO'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'baseflow_flux_land'             ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'BASEFLOW'                  ,&
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

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'snow_frozen_fraction_layer_1' ,&
    UNITS              = '1'                            ,&
    SHORT_NAME         = 'FICE1'                        ,&
    DIMS               = MAPL_DimsTileOnly              ,&
    VLOCATION          = MAPL_VLocationNone             ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'snow_frozen_fraction_layer_2' ,&
    UNITS              = '1'                            ,&
    SHORT_NAME         = 'FICE2'                        ,&
    DIMS               = MAPL_DimsTileOnly              ,&
    VLOCATION          = MAPL_VLocationNone             ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'snow_frozen_fraction_layer_3' ,&
    UNITS              = '1'                            ,&
    SHORT_NAME         = 'FICE3'                        ,&
    DIMS               = MAPL_DimsTileOnly              ,&
    VLOCATION          = MAPL_VLocationNone             ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_emitted_longwave_flux',&
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
    LONG_NAME          = 'total_latent_heat_flux_consistent_with_evaporation_from_turbulence'  ,&
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
    LONG_NAME          = 'surface_temperature_of_land_incl_snow',& 
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSURF'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_snow_on_land',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSNOW'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_unsaturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPUNST'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_saturated_zone',&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPSAT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_temperature_of_wilting_zone'   ,&
    UNITS              = 'K'                         ,&
    SHORT_NAME         = 'TPWLT'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'fractional_area_of_snow_on_land',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ASNOW'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

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
    LONG_NAME          = 'snow_depth_within_snow_covered_area_fraction_on_land' ,&
    UNITS              = 'm'                         ,&
    SHORT_NAME         = 'SNOWDP'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_wetness_surface'      ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET1'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_wetness_rootzone'    ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET2'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_wetness_profile'   ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'WET3'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_moisture_surface'       ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCSF'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_moisture_rootzone'           ,&
    UNITS              = 'm3 m-3'                    ,&
    SHORT_NAME         = 'WCRZ'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_moisture_profile'            ,&
    UNITS              = 'm3 m-3'                   ,&
    SHORT_NAME         = 'WCPR'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperature_layer_1' ,&
    UNITS              = 'K'                         ,&  ! units now K, rreichle & borescan, 6 Nov 2020
    SHORT_NAME         = 'TP1'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperature_layer_2' ,&
    UNITS              = 'K'                         ,&  ! units now K, rreichle & borescan, 6 Nov 2020
    SHORT_NAME         = 'TP2'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperature_layer_3' ,&
    UNITS              = 'K'                         ,&  ! units now K, rreichle & borescan, 6 Nov 2020
    SHORT_NAME         = 'TP3'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperature_layer_4' ,&
    UNITS              = 'K'                         ,&  ! units now K, rreichle & borescan, 6 Nov 2020
    SHORT_NAME         = 'TP4'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperature_layer_5' ,&
    UNITS              = 'K'                         ,&  ! units now K, rreichle & borescan, 6 Nov 2020
    SHORT_NAME         = 'TP5'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'soil_temperature_layer_6' ,&
    UNITS              = 'K'                         ,&  ! units now K, rreichle & borescan, 6 Nov 2020
    SHORT_NAME         = 'TP6'                       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_emissivity'        ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'EMIS'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_reflectivity_visible_beam',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBVR'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_reflectivity_visible_diffuse',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBVF'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_reflectivity_near_infrared_beam',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'ALBNR'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    LONG_NAME          = 'surface_reflectivity_near_infrared_diffuse',&
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
    LONG_NAME          = 'change_upward_sensible_heat_flux',&
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

  call MAPL_AddExportSpec(GC,                         &
    SHORT_NAME         = 'ACCUM',                     &
    LONG_NAME          = 'net_ice_accumulation_rate', &
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'EVLAND',                    &
    LONG_NAME          = 'total_evapotranspiration_land',          &
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
    LONG_NAME          = 'surface_downwelling_PAR_beam_flux', &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'DFPARLAND',                 &
    LONG_NAME          = 'surface_downwelling_PAR_diffuse_flux', &
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
    LONG_NAME          = 'Net_shortwave_flux_snow',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWUPSNOW',                    &
    LONG_NAME          = 'surface_emitted_longwave_flux_snow',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWDNSNOW',                    &
    LONG_NAME          = 'surface_absorbed_longwave_flux_snow',         &
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
    LONG_NAME          = 'Net_shortwave_flux_land',        &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SWDOWNLAND',                &
    LONG_NAME          = 'Incident_shortwave_flux_land',   &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'LWLAND',                    &
    LONG_NAME          = 'Net_longwave_flux_land',         &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHLAND',                    &
    LONG_NAME          = 'Ground_heating_flux_land',       &
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'GHTSKIN',                   &
    LONG_NAME          = 'Ground_heating_flux_for_skin_temp_land',  &
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
    LONG_NAME          = 'total_water_storage_land',  &
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
    SHORT_NAME         = 'SPLAND',                    &              ! a.k.a. SPSHLAND
    LONG_NAME          = 'Spurious_sensible_heat_flux_land',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPLH',                      &
    LONG_NAME          = 'Spurious_latent_heat_flux_land',&
    UNITS              = 'W m-2',                     &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPWATR',                    &              ! a.k.a. SPEVLAND
    LONG_NAME          = 'Spurious_evapotranspiration_flux_land',&
    UNITS              = 'kg m-2 s-1',                &
    DIMS               = MAPL_DimsTileOnly,           &
    VLOCATION          = MAPL_VLocationNone,          &
                                           RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                    &
    SHORT_NAME         = 'SPSNOW',                    &
    LONG_NAME          = 'Spurious_snow_energy_flux_land',&
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
    LONG_NAME          = 'CN_fine_root_carbon'       ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNFROOTC'                  ,&
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
    LONG_NAME          = 'CN_autotrophic_respiration' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNAR'                      ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_heterotrophic_respiration' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNHR'                      ,&
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
    LONG_NAME          = 'CN_fire_count'             ,&
    UNITS              = 'count km-2 s-1'            ,&
    SHORT_NAME         = 'CNFIRE_CNT'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_peat_C_loss_to_fire'    ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNSOM_CLOSS'               ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_N_deployed_to_growth_storage',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNDEPLOY'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_denitrification_rate '  ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNDENIT'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_soil_min_N_loss_to_leaching',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNSMINN_LEACHED'           ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)             

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_soil_mineral_N'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNSMINN'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_N_loss_to_fire'         ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNFIRE_NLOSS'              ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_leaf_N'                 ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNLEAFN'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_leaf_C'                 ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNLEAFC'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_gross_N_mineralization_rate',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNGROSS_NMIN'              ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)               

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_net_N_mineralization_rate',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNET_NMIN'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_N_fixation_to_soil_min_N',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNFIX_TO_SMINN'           ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_actual_N_immobilization',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNACTUAL_IMMOB'            ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_fraction_potential_gpp' ,&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNFPG'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_fraction_potential_immobilization',&
    UNITS              = '1'                         ,&
    SHORT_NAME         = 'CNFPI'                     ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)                         

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_soil_min_N_plant_uptake',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNSMINN_TO_PLANT'          ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_deployment_soil_min_N_uptake' ,&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNSMINN_TO_NPOOL'          ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_atm_N_dep_to_soil_min_N',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNNDEP_TO_SMINN'           ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_vegetation_N'     ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNTOTVEGN'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_litter_N'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNTOTLITN'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)             

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_soil_organic_N'   ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNTOTSOMN'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_plant_retranslocated_N' ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNRETRANSN'                ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_deployment_retranslocated_N',&
    UNITS              = 'kg m-2 s-1'                ,&
    SHORT_NAME         = 'CNRETRANSN_TO_NPOOL'       ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_fuel_C'                 ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNFUELC'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_litter_C'         ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNTOTLITC'                 ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)      

  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_coarse_woody_debris_C'  ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNCWDC'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  ) 
  VERIFY_(STATUS)
  
  call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'CN_total_root_C'           ,&
    UNITS              = 'kg m-2'                    ,&
    SHORT_NAME         = 'CNROOT'                   ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
 
   call MAPL_AddExportSpec(GC                         ,&
    LONG_NAME          = 'fire season length'        ,&
    UNITS              = 'days'                      ,&
    SHORT_NAME         = 'CNFSEL'                    ,&
    DIMS               = MAPL_DimsTileOnly           ,&
    VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )  

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
    LONG_NAME          = 'CO2 Surface Concentration used'  ,&
    UNITS              = '1e-6'                      ,&
    SHORT_NAME         = 'CNCO2'                     ,&
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

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'depth_to_water_table_from_surface_in_peat',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'PEATCLSM_WATERLEVEL'       ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'change_in_free_surface_water_reservoir_on_peat',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'PEATCLSM_FSWCHANGE'        ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,& 
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL1',&
       UNITS              = 'm'                         ,&     
       SHORT_NAME         = 'DZGT1'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL2',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT2'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL3',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT3'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL4',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT4'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL5',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT5'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSOIL6',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZGT6'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_PRMC_and_GWETPROF',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZPR'                      ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_RZMC_and_GWETROOT',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZRZ'                      ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_SFMC_and_GWETTOP',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZSF'                      ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'thickness_of_soil_layer_associated_with_TSATLAND_TUNSTLAND_and_TWLTLAND',&
       UNITS              = 'm'                         ,&
       SHORT_NAME         = 'DZTS'                      ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'soil_wilting_point_in_equivalent_mass_of_total_profile_water',&
       UNITS              = 'kg m-2'                    ,&
       SHORT_NAME         = 'WPEMW'                     ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                  ,&
       LONG_NAME          = 'soil_wilting_point_in_volumetric_units',&
       UNITS              = 'm3 m-3'                    ,&
       SHORT_NAME         = 'WPMC'                      ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RC=STATUS  )
  VERIFY_(STATUS)

  
!EOS

    call MAPL_TimerAdd(GC,    name="RUN1"  ,RC=STATUS)
    VERIFY_(STATUS)
    if (OFFLINE_MODE /=0) then
       call MAPL_TimerAdd(GC,    name="-RUN0"  ,RC=STATUS)
       VERIFY_(status)
    end if
    call MAPL_TimerAdd(GC,    name="-SURF" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="-CATCHCNCLM51" ,RC=STATUS)
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

! !DESCRIPTION: Compute roughness length and exchange coefficients ("cds"), incl. derivatives
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

    ! for analytical extra derivatives (louissurface)

    real, dimension(:,:), pointer :: delCH_delTVA
    real, dimension(:,:), pointer :: delCQ_delTVA

    ! for numerical  extra derivatives (louissurface, helfsurface)

    real, dimension(:,:), pointer :: delCH_delTC
    real, dimension(:,:), pointer :: delCQ_delQC

    real, dimension(:,:), pointer :: cncol
    real, dimension(:,:), pointer :: cnpft

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

   integer                        :: CHOOSEZ0
   real                           :: SCALE4Z0
   real                           :: SCALE4ZVG
   real                           :: SCALE4Z0_u
   real                           :: MIN_VEG_HEIGHT

   ! ------------------------------------- 
   !
   ! for numerical extra derivatives (louissurface, helfsurface)
   
   real, parameter    :: MOSFC_pert_fac = 0.001   ! size of multiplicative pert for numerical derivatives
   
   ! Louis needs 2d arrays;  Helfand would work with 1d arrays but use 2d arrays to avoid "if" statements
   
   real, dimension(:,:), allocatable :: DeltaTC,  CHpert
   real, dimension(:,:), allocatable :: DeltaQC,  CQpert
   
   real, dimension(:,:), allocatable :: DummyZ0T, DummyCM
   
   ! -------------------------------------

! gkw: for CN model
! -----------------
    integer, parameter :: nveg  = num_veg ! number of vegetation types
    integer, parameter :: nzone = num_zon ! number of stress zones

    integer, allocatable :: ityp(:,:,:)
    real,    allocatable :: fveg(:,:,:), elai(:,:,:), esai(:,:,:), tlai(:,:,:), wtzone(:,:), lai1(:), lai2(:), wght(:)

    real,pointer,dimension(:) :: lats
    real,pointer,dimension(:) :: lons

    integer :: nv, nz, ib
    real    :: bare
    logical, save :: first = .true.
    integer(INT64), save :: istep_cn = 0 ! gkw: legacy variable from offline
    real :: ndt
    integer(INT64) :: nstep_cn

  ! Offline mode

   type(CATCHCN_WRAP)             :: wrap
   type(T_CATCHCN_STATE), pointer :: catchcn_internal
   integer                        :: OFFLINE_MODE

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
    call ESMF_UserCompGetInternalState(gc, 'CatchcnInternal', wrap, status)
    VERIFY_(status)
    catchcn_internal => wrap%ptr
    OFFLINE_MODE = catchcn_internal%CATCH_OFFLINE          ! shorthand

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

    call MAPL_Get ( MAPL                                 ,&
            TILELATS  = LATS                             ,&      ! [radians]
            TILELONS  = LONS                             ,&      ! [radians]
            INTERNAL_ESMF_STATE = INTERNAL               ,&
            RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetResource ( MAPL, CHOOSEZ0, Label="CHOOSEZ0:", DEFAULT=3, RC=STATUS)
    VERIFY_(STATUS)     
    call ESMF_VMGetCurrent(VM,       rc=STATUS)
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
   
   if     (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND == 1) then   
      
      call MAPL_GetPointer(INTERNAL,delCH_delTVA , 'delCH_delTVA' ,    RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(INTERNAL,delCQ_delTVA , 'delCQ_delTVA' ,    RC=STATUS)
      VERIFY_(STATUS)
      
   elseif (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND >= 2) then 
      
      call MAPL_GetPointer(INTERNAL,delCQ_delQC  , 'delCQ_delQC'  ,    RC=STATUS)
      VERIFY_(STATUS)   
      call MAPL_GetPointer(INTERNAL,delCH_delTC  , 'delCH_delTC'  ,    RC=STATUS)
      VERIFY_(STATUS)
      
   end if

   call MAPL_GetPointer(INTERNAL,CNCOL  ,'CNCOL'  ,   RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(INTERNAL,CNPFT  ,'CNPFT'  ,   RC=STATUS)
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

   if (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND>=2) then
       
      ! allocate variables for numerical extra derivatives (louissurface, helfsurface)
       
      allocate(DeltaTC( NT,NUM_SUBTILES),STAT=STATUS);  VERIFY_(STATUS)
      allocate(DeltaQC( NT,NUM_SUBTILES),STAT=STATUS);  VERIFY_(STATUS)
      
      allocate(CHpert(  NT,NUM_SUBTILES),STAT=STATUS);  VERIFY_(STATUS)
      allocate(CQpert(  NT,NUM_SUBTILES),STAT=STATUS);  VERIFY_(STATUS)
      
      allocate(DummyZ0T(NT,NUM_SUBTILES),STAT=STATUS);  VERIFY_(STATUS)
      allocate(DummyCM( NT,NUM_SUBTILES),STAT=STATUS);  VERIFY_(STATUS)
      
   end if
   
   allocate(   ityp(nt,nveg,nzone) )
   allocate(   fveg(nt,nveg,nzone) )
   allocate( wtzone(nt,nzone) )
   allocate(   elai(nt,nveg,nzone) )
   allocate(   esai(nt,nveg,nzone) )
   allocate(   tlai(nt,nveg,nzone) )

   allocate ( lai1(nt) )
   allocate ( lai2(nt) )
   allocate ( wght(nt) )

!  Vegetation types used to index into tables
!--------------------------------------------

   where(ITY(:,1) > 0.)
     VEG1 = map_cat(nint(ITY(:,1)))  ! gkw: primary   CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   where(ITY(:,2) > 0.)
     VEG2 = map_cat(nint(ITY(:,2)))  ! gkw: secondary CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   _ASSERT((count(VEG1>NTYPS.or.VEG1<1)==0),'needs informative message')
   _ASSERT((count(VEG2>NTYPS.or.VEG2<1)==0),'needs informative message')

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
   
   FVG1 = fvg(:,1) 
   FVG2 = fvg(:,2)

! set CLM CN PFT & fraction, set carbon zone weights
! --------------------------------------------------
   do nz = 1,nzone
     ityp(:,:,nz) = nint(ity(:,:))
     fveg(:,:,nz) = fvg(:,:)
     wtzone(:,nz) = CN_zone_weight(nz)
   end do

! call to set CN time step before any other CN routines are called (jkolassa May 2023)
! ------------------------------------------------------------------------------------------
  catchcn_internal%DTCN = min(catchcn_internal%DTCN,14400.)
  ndt = get_step_size( nint(catchcn_internal%DTCN) ) ! gkw: get_step_size must be called here to set CN model time step

! update CN time step number
! --------------------------
  nstep_cn = get_nstep(istep_cn) 

! initialize CN model and transfer restart variables on startup
! -------------------------------------------------------------
   if(first) then
      call CN_init(nt,ityp,fveg,cncol,cnpft,lats,lons,catchcn_internal%DTCN,water_inst,bgc_vegetation_inst,.true.) 
      call get_CN_LAI(nt,ityp,fveg,elai,esai=esai)
      first = .false.
   endif

   ! For the OFFLINE case, first update some diagnostic vars
   if (OFFLINE_MODE /=0) then
      call MAPL_TimerOn(MAPL, "-RUN0")
      call RUN0(gc, import, export, clock, rc)
      call MAPL_TimerOff(MAPL, "-RUN0")
   end if

! obtain LAI from previous time step (from CN model)
! --------------------------------------------------

   call get_CN_LAI(nt,ityp,fveg,elai,esai=esai,tlai=tlai)   
   
   lai1 = 0.
   wght = 0.
   do nz = 1,nzone
     nv = 1
     lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
     wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
   end do
   lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

   lai2 = 0.
   wght = 0.
   do nz = 1,nzone
     nv = 2
     lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
     wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
   end do
   lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

   lai = fvg1*lai1 + fvg2*lai2   ! gkw: this is a VEGDYN import

   deallocate ( ityp )
   deallocate ( fveg )
   deallocate ( elai )
   deallocate ( esai )
   deallocate ( wtzone )
   deallocate ( tlai   )
   
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

    select case (CATCHCN_INTERNAL%Z0_FORMULATION)
       case (0)  ! no scaled at all
          SCALE4ZVG   = 1
          SCALE4Z0    = 1
          SCALE4Z0_u  = 1
          MIN_VEG_HEIGHT = 0.01
       case (1) ! This case is bugged
          SCALE4ZVG   = 1
          SCALE4Z0    = 2
          SCALE4Z0_u  = 1
          MIN_VEG_HEIGHT = 0.01
       case (2)
          SCALE4ZVG   = 1
          SCALE4Z0    = 2
          SCALE4Z0_u  = 2
          MIN_VEG_HEIGHT = 0.01
       case (3)
          SCALE4ZVG   = 0.5
          SCALE4Z0    = 1
          SCALE4Z0_u  = 1
          MIN_VEG_HEIGHT = 0.01
       case (4)
          SCALE4ZVG   = 1
          SCALE4Z0    = 2
          SCALE4Z0_u  = 2
          MIN_VEG_HEIGHT = 0.1
    end select

    SUBTILES: do N=1,NUM_SUBTILES

      ! jkolassa Jul 2025: For Z0-formulation == 4 use old (6-class) veg type
      !                    and mix two veg types
      !                    Consider updating to a 15-PFT resolution in the future
      if (CATCHCN_INTERNAL%Z0_FORMULATION == 4) then
         ! make canopy height >= min veg height:
         Z2CH = max(Z2CH,MIN_VEG_HEIGHT)
         ZVG  = fvg1*(Z2CH - SAI4ZVG(VEG1)*SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI1)) + &
                fvg2*(Z2CH - SAI4ZVG(VEG2)*SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI2))
      else
         ZVG  = fvg1*(Z2CH - SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI1)) + &
                fvg2*(Z2CH - SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI2))
         !Z0T(:,N)  = Z0_BY_ZVEG*ZVG*SCALE4Z0 
      endif

!  For now roughnesses and displacement heights
!   are the same for all subtiles.

      Z0T(:,N)  = Z0_BY_ZVEG*ZVG*SCALE4Z0
   IF (catchcn_internal%USE_ASCATZ0 == 1) THEN
      WHERE (NDVI <= 0.2)
         Z0T(:,N)  = ASCATZ0
      END WHERE
   ENDIF
   D0T  = D0_BY_ZVEG*ZVG

   DZE  = max(DZ - D0T, 10.)

   if(associated(Z0 )) Z0  = Z0T(:,N)
   if(associated(D0 )) D0  = D0T


!  Compute surface exchange coefficients
!---------------------------------------
      
    call MAPL_TimerOn(MAPL,"-SURF")    ! timer for computation of MOSFC exchange coeffs and derivs (Louis or Helfand)
      
    if (CATCHCN_INTERNAL%CHOOSEMOSFC.eq.0) then

       ! Louis surface turbulence
         
       WW(:,N) = 0.
       CM(:,N) = 0.
       
       if (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND==1) then
          
          ! analytical extra derivatives (default for Louis)
          
          call louissurface(3,N,UU,WW,PS,TA,TC,QA,QC,PCU,LAI,Z0T,DZE,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE,delCH_delTVA,delCQ_delTVA)
          
       else
          
          ! none .or. numerical extra derivatives
          
          if (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND>=2) then
             
             ! Prep calculation of numerical extra derivatives.  Start with calling louissurface with perturbed inputs and
             !  save only the perturbed exchange coeffs.  The final call with nominal inputs produces the unperturbed
             !  exchange coeffs and other outputs (CN, RIB, ZT, ZQ, etc).
             ! Must use properly initialized dummmies for Z0T and CM because these are intent(inout).
             
             ! perturb TC: send in (TC+DeltaTC), get back CHpert
             
             DeltaTC  = MOSFC_pert_fac*TC
             
             DummyZ0T = Z0T
             DummyCM  = CM
             
             call louissurface(   3,N,UU,WW,PS,TA,TC+DeltaTC,QA,QC        ,PCU,LAI,DummyZ0T,DZE,DummyCM,CN,RIB,ZT,ZQ,CHpert,CQ    ,UUU,UCN,RE)
             
             if (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND==2) then
                
                ! perturb QC: send in (QC+DeltaQC), get back CQpert
                
                DeltaQC  = MOSFC_pert_fac*QC
                
                DummyZ0T = Z0T
                DummyCM  = CM
                
                call louissurface(3,N,UU,WW,PS,TA,TC        ,QA,QC+DeltaQC,PCU,LAI,DummyZ0T,DZE,DummyCM,CN,RIB,ZT,ZQ,CH    ,CQpert,UUU,UCN,RE)
                
             end if
             
          end if
          
          ! Call with nominal inputs [after calls with perturbed inputs to obtain correct outputs (CN, RIB, ZT, ZQ, etc.)]
          
          call louissurface(3,N,UU,WW,PS,TA,TC,QA,QC,PCU,LAI,Z0T,DZE,CM,CN,RIB,ZT,ZQ,CH,CQ,UUU,UCN,RE)
          
       end if  ! MOSFC_EXTRA_DERIVS_OFFL_LAND
       
    elseif (CATCHCN_INTERNAL%CHOOSEMOSFC.eq.1)then
       
       ! Helfand surface turbulence
       
       niter  = 6                  ! number of internal iterations in the Helfand MO surface layer routine
       IWATER = 3
       
       PSMB = PS * 0.01            ! convert to MB
       
       ! Approximate pressure at top of surface layer: hydrostatic, eqn of state using avg temp and press
       PSL = PSMB * (1. - (DZE*MAPL_GRAV)/(MAPL_RGAS*(TA+TC(:,N)) ) ) /   &
            (1. + (DZE*MAPL_GRAV)/(MAPL_RGAS*(TA+TC(:,N)) ) )
       
       if (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND==2) then
          
          ! Prep calculation of numerical extra derivatives.  Start with calling louissurface with perturbed inputs and
          !  save only the perturbed exchange coeffs.  The final call with nominal inputs produces the unperturbed
          !  exchange coeffs and other outputs (CN, RIB, ZT, ZQ, etc).
          ! Must use properly initialized dummmies for Z0T and CM because these are intent(inout).
          
          ! perturb TC: send in (TC+DeltaTC), get back CHpert
          
          DeltaTC( :,N) = MOSFC_pert_fac*TC(:,N)
          
          DummyZ0T(:,N) = Z0T(:,N)
          
          CALL helfsurface(UWINDLMTILE,VWINDLMTILE,TA,TC(:,N)+DeltaTC(:,N),QA,QC(:,N)             ,PSL,PSMB,DummyZ0T(:,N),lai,  &
               IWATER,DZE,niter,nt,RHOH,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,                                                  &
               t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0)
          
          CHpert(  :,N) = VKH
          
          ! perturb QC: send in (QC+DeltaQC), get back CQpert
          
          DeltaQC( :,N) = MOSFC_pert_fac*QC(:,N)
          
          DummyZ0T(:,N) = Z0T(:,N)
          
          CALL helfsurface(UWINDLMTILE,VWINDLMTILE,TA,TC(:,N)             ,QA,QC(:,N)+DeltaQC(:,N),PSL,PSMB,DummyZ0T(:,N),lai,  &
               IWATER,DZE,niter,nt,RHOH,VKH,VKM,USTAR,XX,YY,CU,CT,RIB,ZETA,WS,                                                  &
               t2m,q2m,u2m,v2m,t10m,q10m,u10m,v10m,u50m,v50m,CHOOSEZ0)
          
          CQpert(  :,N) = VKH
          
       end if  ! MOSFC_EXTRA_DERIVS_OFFL_LAND==2
       
       ! Call with nominal inputs [after calls with perturbed inputs to obtain correct outputs (Z0T, [*]2m, [*]10m, etc.)]
       
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

    endif  ! CHOOSEMOSFC

    if     (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND==2) then                        
       
       ! finalize numerical derivatives
       
       delCH_delTC(:,N) = (CHpert(:,N) - CH(:,N)) / DeltaTC(:,N)
       delCQ_delQC(:,N) = (CQpert(:,N) - CQ(:,N)) / DeltaQC(:,N)
         
    elseif (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND==3) then                        
       
       ! finalize numerical derivatives (valid for Louis only!)
       
       ! For Louis, the exchange coeffs depend only on the *virtual* temperature (not true for Helfand).
       !
       ! This lets us compute the derivatives of the exchange coefficients w.r.t. both TC and QC from 
       ! just one additional call to louissurface() with perturbed TC.
       !
       ! In the following, "del" indicates a derivative and "Delta" indicates a difference term.
       !
       ! We have:
       !
       ! (1) CH = CQ
       !
       ! (2) TVC = TC*(1 + eps*QC)    (virtual temperature; eps = MAPL_VIREPS)
       !
       !     (2a) ==> delTVC_delQC = eps*TC
       !
       !     (2b) ==> DeltaTVC = (TVCpert - TVC) = TCpert*(1 + eps*QC) - TC*(1 + eps*QC) = DeltaTC*(1 + eps*QC)
       !
       ! (3) delCH_delTC  = (CHpert - CH)/DeltaTC
       !
       ! (4) delCH_delTVC = (CHpert - CH)/DeltaTVC = (CHpert - CH)/(DeltaTC*(1 + eps*QC))
       !
       ! Using (1)-(4), we have:
       !
       ! delCQ_delQC = delCH_delQC                                                   using (1)
       !                                                                             
       !             = delCH_delTVC                         * delTVC_delQC           using chain rule
       !                                                                             
       !             = (CHpert - CH)/(DeltaTC*(1 + eps*QC)) * delTVC_delQC           using (4)      
       !                          
       !             = (CHpert - CH)/DeltaTC * 1/(1+eps*QC) * eps*TC                 using (2a)
       !
       !             = delCH_delTC           * 1/(1+eps*QC) * eps*TC                 using (3)
       
       delCH_delTC(:,N) = (CHpert(:,N) - CH(:,N)) / DeltaTC(:,N)
       
       delCQ_delQC(:,N) = delCH_delTC(:,N) * MAPL_VIREPS*TC(:,N)/(1.+MAPL_VIREPS*QC(:,N))
       
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
      if(associated(VNT)) VNT     = VNT + UUU            *FR(:,N)

      WW(:,N) = max(CH(:,N)*(TC(:,N)-TA-(MAPL_GRAV/MAPL_CP)*DZE)/TA + MAPL_VIREPS*CQ(:,N)*(QC(:,N)-QA),0.0)
      WW(:,N) = (HPBL*MAPL_GRAV*WW(:,N))**(2./3.)
      if(associated(GST)) GST     = GST + WW(:,N)        *FR(:,N)

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
   if (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND>=2) then
      deallocate(DeltaTC )
      deallocate(DeltaQC )
      deallocate(CHpert  )
      deallocate(CQpert  )
      deallocate(dummyZ0T)
      deallocate(dummyCM )
   end if
   
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

    type(MAPL_MetaComp),pointer      :: MAPL
    type(ESMF_Alarm)                 :: ALARM

    integer                          :: IM,JM

    real                             :: SCALE4ZVG
    real                             :: SCALE4Z0_u
    real                             :: MIN_VEG_HEIGHT

    type(ESMF_VM)                    :: VM
    type (T_CATCHCN_STATE), pointer  :: CATCHCN_INTERNAL
    type (CATCHCN_WRAP)              :: wrap
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

     ! Get component's private internal state
    call ESMF_UserCompGetInternalState(gc, 'CatchcnInternal', wrap, status)
    VERIFY_(status)
    CATCHCN_INTERNAL=>wrap%ptr

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL, RUNALARM=ALARM, RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_VMGetCurrent(VM,       rc=STATUS)

   select case (CATCHCN_INTERNAL%Z0_FORMULATION)
      case (0)  ! no scaled at all
         SCALE4ZVG   = 1
         SCALE4Z0_u  = 1
         MIN_VEG_HEIGHT = 0.01
      case (1) ! This case is bugged
         SCALE4ZVG   = 1
         SCALE4Z0_u  = 1
         MIN_VEG_HEIGHT = 0.01
      case (2)
         SCALE4ZVG   = 1
         SCALE4Z0_u  = 2
         MIN_VEG_HEIGHT = 0.01
      case (3)
         SCALE4ZVG   = 0.5
         SCALE4Z0_u  = 1
         MIN_VEG_HEIGHT = 0.01
      case (4)
         SCALE4ZVG   = 1
         SCALE4Z0_u  = 2
         MIN_VEG_HEIGHT = 0.1
   end select

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

        real, dimension(:),   pointer :: THATM
        real, dimension(:),   pointer :: QHATM
        real, dimension(:),   pointer :: CTATM
        real, dimension(:),   pointer :: CQATM
        real, dimension(:),   pointer :: ICE
        real, dimension(:),   pointer :: FRZR
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
        real, dimension(:),   pointer :: snowalb
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
        real, dimension(:,:), pointer :: delCQ_delTVA
        real, dimension(:,:), pointer :: delCH_delTVA
        real, dimension(:,:), pointer :: delCH_delTC
        real, dimension(:,:), pointer :: delCQ_delQC
        real, dimension(:),   pointer :: tile_id
        real, dimension(:),   pointer :: ndep
        real, dimension(:),   pointer :: abm
        real, dimension(:),   pointer :: peatf
        real, dimension(:),   pointer :: gdp
        real, dimension(:),   pointer :: hdm
        real, dimension(:),   pointer :: fieldcap        
        real, dimension(:),   pointer :: cli_t2m
        real, dimension(:),   pointer :: bgalbvr
        real, dimension(:),   pointer :: bgalbvf
        real, dimension(:),   pointer :: bgalbnr
        real, dimension(:),   pointer :: bgalbnf
        real, dimension(:,:), pointer :: cncol
        real, dimension(:,:), pointer :: cnpft
        real, dimension(:,:), pointer :: tgwm
        real, dimension(:,:), pointer :: rzmm
        real, dimension(:,:), pointer :: sfmm
        real, dimension(:),   pointer :: bflowm
        real, dimension(:),   pointer :: totwatm
        real, dimension(:),   pointer :: tairm
        real, dimension(:),   pointer :: rhm    
        real, dimension(:),   pointer :: windm  
        real, dimension(:),   pointer :: rainfm 
        real, dimension(:),   pointer :: snowfm 
        real, dimension(:),   pointer :: runsrfm
        real, dimension(:),   pointer :: ar1m   
        real, dimension(:),   pointer :: tpm
        real, dimension(:),   pointer :: cnsum
        real, dimension(:,:,:), pointer :: psnsunm
        real, dimension(:,:,:), pointer :: psnsham
        real, dimension(:,:,:), pointer :: lmrsunm
        real, dimension(:,:,:), pointer :: lmrsham
        real, dimension(:,:,:), pointer :: laisunm
        real, dimension(:,:,:), pointer :: laisham
        real, dimension(:),   pointer :: sndzm
        real, dimension(:),   pointer :: sndzm5d
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

        real, dimension(:),   pointer :: T2M10D
        real, dimension(:),   pointer :: TG10D
        real, dimension(:),   pointer :: T2MMIN5D
        real, dimension(:),   pointer :: RH30D
        real, dimension(:),   pointer :: TPREC10D
        real, dimension(:),   pointer :: TPREC60D
        real, dimension(:),   pointer :: ET365D

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
        real, dimension(:),   pointer :: fice1 
        real, dimension(:),   pointer :: fice2 
        real, dimension(:),   pointer :: fice3 
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
        real, dimension(:),   pointer :: SPLH
        real, dimension(:),   pointer :: SPWATR
        real, dimension(:),   pointer :: SPSNOW

        real, dimension(:),   pointer :: CNLAI
        real, dimension(:),   pointer :: CNTLAI
        real, dimension(:),   pointer :: CNSAI
        real, dimension(:),   pointer :: CNTOTC
        real, dimension(:),   pointer :: CNVEGC
        real, dimension(:),   pointer :: CNFROOTC
        real, dimension(:),   pointer :: CNNPP
        real, dimension(:),   pointer :: CNGPP
        real, dimension(:),   pointer :: CNSR
        real, dimension(:),   pointer :: CNAR
        real, dimension(:),   pointer :: CNHR
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
        real, dimension(:),   pointer :: CNCO2 
        real, dimension(:),   pointer :: CNFIRE_CNT         
        real, dimension(:),   pointer :: CNSOM_CLOSS        
        real, dimension(:),   pointer :: CNNDEPLOY            
        real, dimension(:),   pointer :: CNDENIT            
        real, dimension(:),   pointer :: CNSMINN_LEACHED    
        real, dimension(:),   pointer :: CNSMINN            
        real, dimension(:),   pointer :: CNFIRE_NLOSS       
        real, dimension(:),   pointer :: CNLEAFN            
        real, dimension(:),   pointer :: CNLEAFC            
        real, dimension(:),   pointer :: CNGROSS_NMIN       
        real, dimension(:),   pointer :: CNNET_NMIN         
        real, dimension(:),   pointer :: CNNFIX_TO_SMINN    
        real, dimension(:),   pointer :: CNACTUAL_IMMOB     
        real, dimension(:),   pointer :: CNFPG              
        real, dimension(:),   pointer :: CNFPI              
        real, dimension(:),   pointer :: CNSMINN_TO_PLANT   
        real, dimension(:),   pointer :: CNSMINN_TO_NPOOL   
        real, dimension(:),   pointer :: CNNDEP_TO_SMINN    
        real, dimension(:),   pointer :: CNTOTVEGN          
        real, dimension(:),   pointer :: CNTOTLITN          
        real, dimension(:),   pointer :: CNTOTSOMN          
        real, dimension(:),   pointer :: CNRETRANSN         
        real, dimension(:),   pointer :: CNRETRANSN_TO_NPOOL
        real, dimension(:),   pointer :: CNFUELC            
        real, dimension(:),   pointer :: CNTOTLITC          
        real, dimension(:),   pointer :: CNCWDC             
        real, dimension(:),   pointer :: CNROOT            
        real, dimension(:),   pointer :: CNFSEL

        real, dimension(:),   pointer :: WAT10CM
        real, dimension(:),   pointer :: WATSOI
        real, dimension(:),   pointer :: ICESOI
        real, dimension(:),   pointer :: SHSNOW
        real, dimension(:),   pointer :: AVETSNOW

        real, dimension(:),   pointer :: RMELTDU001
        real, dimension(:),   pointer :: RMELTDU002
        real, dimension(:),   pointer :: RMELTDU003
        real, dimension(:),   pointer :: RMELTDU004
        real, dimension(:),   pointer :: RMELTDU005
        real, dimension(:),   pointer :: RMELTBC001
        real, dimension(:),   pointer :: RMELTBC002
        real, dimension(:),   pointer :: RMELTOC001
        real, dimension(:),   pointer :: RMELTOC002
                                   
        real, dimension(:),   pointer :: PEATCLSM_WATERLEVEL
        real, dimension(:),   pointer :: PEATCLSM_FSWCHANGE

        real, dimension(:),   pointer :: DZGT1
        real, dimension(:),   pointer :: DZGT2
        real, dimension(:),   pointer :: DZGT3
        real, dimension(:),   pointer :: DZGT4
        real, dimension(:),   pointer :: DZGT5
        real, dimension(:),   pointer :: DZGT6
        real, dimension(:),   pointer :: DZPR
        real, dimension(:),   pointer :: DZRZ
        real, dimension(:),   pointer :: DZSF
        real, dimension(:),   pointer :: DZTS 
        real, dimension(:),   pointer :: WPEMW
        real, dimension(:),   pointer :: WPMC

        ! --------------------------------------------------------------------------
        ! Local pointers for tile variables
        ! --------------------------------------------------------------------------

        INTEGER,pointer,dimension(:) :: CAT_ID
        real,pointer,dimension(:) :: DZSF_in_mm
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
        real,pointer,dimension(:) :: TCSORIG1, TPSN1IN1, TPSN1OUT1, FSW_CHANGE
	real,pointer,dimension(:) :: WCHANGE, ECHANGE, HSNACC, LHOUT, EVACC, LHACC, SHACC
	real,pointer,dimension(:) :: SNOVR, SNOVF, SNONR, SNONF
	real,pointer,dimension(:) :: VSUVR, VSUVF
	real,pointer,dimension(:) :: ALWX, BLWX
	real,pointer,dimension(:) :: fveg1, fveg2
        real,pointer,dimension(:) :: FICE1TMP
        real,pointer,dimension(:) :: SLDTOT

!       real*8,pointer,dimension(:) :: fsum

        real,pointer,dimension(:,:) :: ghtcnt
        real,pointer,dimension(:,:) :: wesnn
        real,pointer,dimension(:,:) :: htsnnn
        real,pointer,dimension(:,:) :: sndzn
        real,pointer,dimension(:,:) :: ficesout
        real,pointer,dimension(:,:) :: shsbt
        real,pointer,dimension(:,:) :: dshsbt
        real,pointer,dimension(:,:) :: evsbt
        real,pointer,dimension(:,:) :: devsbt
        real,pointer,dimension(:,:) :: DEDTC 
        real,pointer,dimension(:,:) :: DHSDQC
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
        real,allocatable,dimension(:) :: rdc_tmp_1, rdc_tmp_2

        ! albedo calculation stuff

        type(ESMF_Config)           :: CF
        type(MAPL_SunOrbit)         :: ORBIT
        type(ESMF_Time)             :: CURRENT_TIME, StopTime, NextTime, NextRecordTime
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
        real                        :: DT
        integer                     :: NTILES
        integer                     :: I, J, K, N

	! dummy variables for call to get snow temp

        real    :: FICE
        logical :: DUMFLAG1,DUMFLAG2
        integer                         :: nmax
        type(ESMF_VM)                   :: VM

#ifdef DBG_CNLSM_INPUTS
        ! vars for debugging purposes
        type(ESMF_Grid)                 :: TILEGRID
        type (MAPL_LocStream)           :: LOCSTREAM
        integer, pointer                :: mask(:)
        integer                         :: nt
        integer, save                   :: unit_i=0
        logical, save                   :: firsttime=.true.
        integer                         :: unit
        integer                         :: NT_GLOBAL

#endif

       ! Offline case

        type(CATCHCN_WRAP)          :: wrap
        type(T_CATCHCN_STATE), pointer :: catchcn_internal
        integer                     :: OFFLINE_MODE
        real,dimension(:,:),allocatable :: ALWN, BLWN
        ! unadulterated TC's and QC's
        real, pointer               :: TC1_0(:), TC2_0(:),  TC4_0(:)
        real, pointer               :: QA1_0(:), QA2_0(:),  QA4_0(:)

        ! CATCHMENT_SPINUP
        integer                     :: CurrMonth, CurrDay, CurrHour, CurrMin, CurrSec

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

    real, allocatable, dimension(:) ::  wgt, wpp, fwet, wet_in 
    real, allocatable, dimension(:,:) :: sm   ! soil water as frac of WHC for the 3 dydrological zones at root depth
    real, allocatable, dimension(:) :: SWSRF1, SWSRF2, SWSRF4    ! soil water as frac of WHC for the 3 dydrological zones at surface soil
    real, allocatable, dimension(:,:) :: tcx, qax
    real, allocatable, dimension(:,:) :: tgw, rzm, sfm,rc00, rcdt,rcdq, totcolc, wtzone
    real, allocatable, dimension(:,:) :: btran_fire, bt
    real, allocatable, dimension(:,:,:) :: btran,elai,esai,fveg,tlai,psnsun,psnsha,laisun,laisha,lmrsun,lmrsha
    integer, allocatable, dimension(:,:,:) :: ityp
    real, allocatable, dimension(:) :: car1, car2, car4
    real, allocatable, dimension(:) :: para
    real, allocatable, dimension(:) :: rcxdt, rcxdq
    real, allocatable, dimension(:) :: dayl, dayl_fac
    real, allocatable, dimension(:), save :: nee, npp, gpp, sr, aresp, hresp, padd, frootc, vegc, xsmr,burn, closs
    real, allocatable, dimension(:) :: nfire, som_closs, fsnow
    real, allocatable, dimension(:) :: ndeploy, denit, sminn_leached, sminn, fire_nloss
    real, allocatable, dimension(:) :: leafn, leafc, gross_nmin, net_nmin, nfix_to_sminn, actual_immob
    real, allocatable, dimension(:) :: fpg, fpi, sminn_to_plant, sminn_to_npool, ndep_to_sminn
    real, allocatable, dimension(:) :: totvegn, totlitn, totsomn, retransn, retransn_to_npool  
    real, allocatable, dimension(:) :: fuelc, totlitc, cwdc, rootc 

    ! ***************************************************************************************************************************************************************
    ! Begin Carbon Tracker variables
    !
    ! use EEA global average CO2 to scale 2001-2014 CarbonTracker CO2 monthly mean diurnal cycle to obtain CO2 for 1850-2000.
    ! extended from the last cycle when carbon reaches equilibrium with the 2001-2014 CarbonTracker CO2 monthly mean diurnal 
    ! cycle * 280ppm/389.8899ppm, fzeng, Apr 2017.
    ! EEA global average CO2 is from http://www.eea.europa.eu/data-and-maps/figures/atmospheric-concentration-of-co2-ppm-1  
    ! --------------------------------------------------------------------------------------------------------------------  

    real               :: co2g                 ! global average atmospheric carbon dioxide concentration, varies after 1850
    integer, parameter :: byr_co2g  = 1851     ! year global average atmospheric CO2 concentration began to increase from 280.e-6 
    integer, parameter :: myr_co2g  = 1950     ! year global average atmospheric CO2 concentration reached 311.e-6 
    integer, parameter :: eyr_co2g  = 2012     ! year global average atmospheric CO2 concentration reached 391.e-6 
    real, parameter    :: co2g_byr  = 280.e-6  ! pre-industrial global average atmospheric carbon dioxide concentration (i.e. before byr_co2g)
    real, parameter    :: co2g_myr  = 311.e-6  ! global average atmospheric CO2 concentration in myr_co2g
    real, parameter    :: co2g_eyr  = 391.e-6  ! global average atmospheric CO2 concentration in eyr_co2g
    real, parameter    :: dco2g_1   = (co2g_myr-co2g_byr)/(myr_co2g-byr_co2g) ! yearly atmospheric CO2 concentration increment for period 1 (byr_co2g to myr_co2g)
    real, parameter    :: dco2g_2   = (co2g_eyr-co2g_myr)/(eyr_co2g-myr_co2g) ! yearly atmospheric CO2 concentration increment for period 2 (myr_co2g to eyr_co2g)
    real, parameter    :: CTco2g = 389.8899e-6     ! Spatial (tile area weighted) and temporal average of 2001-2014 CarbonTracker CO2
    real, allocatable, dimension(:) :: co2v  ! spatial varying atmospheric carbon dioxide concentration

    ! parameters for calculating CT indices for tiles
    ! -----------------------------------------------
    integer, parameter :: CT_grid_N_lon  = 120  ! lon dimension CarbonTracker CO2 data
    integer, parameter :: CT_grid_N_lat  =  90  ! lat dimension CarbonTracker CO2 data
    real, parameter    :: CT_grid_dlon = 360./real(CT_grid_N_lon), CT_grid_dlat = 180./real(CT_grid_N_lat)
    INTEGER            ::  info, comm, CTfile, Y1, M1, This3H, ThisCO2_Year, NUNQ, CO2_YEAR
    logical, allocatable, dimension (:)        :: unq_mask
    integer, allocatable, dimension (:,:)      :: CT_index
    integer, allocatable, dimension (:)        :: ct2cat, ThisIndex, loc_int
    integer, allocatable, dimension (:), save  :: ct_tid
    real, dimension (:,:,:,:),     allocatable :: CTCO2_TMP
    real, dimension (:,:,:), save, allocatable :: CT_CO2V
    logical, save      :: first_ct = .true.
    integer, save      :: FIRST_YY
        
    ! End Carbon Tracker variables
    ! *************************************************************************************************************************************************************

    ! prescribe DYNVEG parameters
    ! ---------------------------

    real, parameter :: dtc = 0.03 ! canopy temperature perturbation (K) [approx 1:10000]
    real, parameter :: dea = 0.10 ! vapor pressure perturbation (Pa) [approx 1:10000]

    real, allocatable, dimension(:) :: totwat ! total soil liquid water (kg/m2)
    real, save :: ashift = 0. ! for baseflow. gkw: this should match value in routine "base" in catchment
    real :: Qair_sat                                 ! saturated specific humidity (kg/kg)
    real, allocatable, dimension(:) :: Qair_relative          ! relative humidity (%)

    integer :: nz, iv
    real :: cn1, cn2, cn3, cn12, cn23, ar, ax1, ax2, ax4
    real, dimension(fsat:fwlt) :: f1, f2, f3, f4

    real, allocatable, dimension(:,:,:,:) ::  albdir, albdif
    integer, allocatable, dimension(:) :: ityp_tmp
    
    ! static summing arrays for CN
    ! ----------------------------
    real, allocatable, dimension(:) :: ht, tp, soilice
    real :: zbar, frice

    real, allocatable, dimension(:,:)  :: col
    real, allocatable, dimension(:,:,:) :: pft
    
    real, allocatable, dimension(:) :: lnfm
    character(len=ESMF_MAXSTR)      :: LNFMFile, CO2_CycleFile

    integer :: ntile, nv, dpy, ierr, iok, ndt
    integer, save :: year_prev = -9999
    
    integer, save :: n1d                ! number of land model steps in a 1-day period
    integer, save :: n5d                ! number of land model steps in a 5-day period
    integer, save :: n10d               ! number of land model steps in a 10-day period
    integer, save :: n30d               ! number of land model steps in a 30-day period
    integer, save :: n60d               ! number of land model steps in a 60-day period    
    integer, save :: n365d              ! number of land model steps in a 365-day period

    ! For accumulated fields
    ! NOTE: In CNPhenologyMod.F90, init_gdd20 is always set to .false. as well. For GEOS-5 runs, need to discard at least the first 2 years.
    ! This is not a problem for offline runs because we always spin up the model whenever we change meterology. fzeng, July 2017 
    ! -------------------------------------------------------------------------------------------------------------------------------------- 
    logical, parameter :: init_accum = .false.! jkolassa May 2023: needs to be set to true if no CNCLM51 restart is available
    logical, parameter :: init_accum_365 = .false.! jkolassa May 2023: needs to be set to true if no CNCLM51 restart is available
    integer, save :: istep                    ! model time step index
    integer, save :: istep_365                ! model time step index
    integer :: accper                         ! number of time steps accumulated in a period of XX days, increases from 1 to nXXd in the first XX days,
    ! and remains as nXXd thereafter
    integer, allocatable, dimension(:) :: ta_count
    real, allocatable, dimension(:)    :: TA_MIN                
   
    integer :: AGCM_YY, AGCM_MM, AGCM_DD, AGCM_MI, AGCM_S, AGCM_HH, dofyr, AGCM_S_ofday
    logical, save :: first = .true.
    integer(INT64), save :: istep_cn = 1 ! gkw: legacy variable from offline

    ! solar declination related
    real :: ob, declin, zs, zc, max_decl, max_dayl
    integer :: year, iday, idayp1

   ! real :: co2
    real, external :: getco2

    ! temporaries for call to SIBALB for each type
    ! --------------------------------------------
    real, allocatable, dimension(:) :: lai1, lai2, wght
    real, allocatable, dimension(:) :: ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp
    real, allocatable, dimension(:) :: SNOVR_tmp, SNONR_tmp, SNOVF_tmp, SNONF_tmp

    logical           :: record
    type(ESMF_Alarm)  :: RecordAlarm

    ! Variables for FPAR
    real   , allocatable, dimension (:,:,:)      :: parzone

    integer :: cn_count = 0 
    logical :: first_cn

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

        ! Get component's private internal state
        call ESMF_UserCompGetInternalState(gc, 'CatchcnInternal', wrap, status)
        VERIFY_(status)
        catchcn_internal => wrap%ptr
        OFFLINE_MODE = catchcn_internal%CATCH_OFFLINE            ! shorthand
        ! if (MAPL_AM_I_Root(VM)) print *, trim(Iam)//'::OFFLINE mode: ', is_OFFLINE

        call ESMF_VMGetCurrent ( VM, RC=STATUS )

        ! --------------------------------------------------------------------------
        ! Get parameters from generic state.
        ! --------------------------------------------------------------------------

        call MAPL_Get ( MAPL                 ,&
             RUNALARM  = ALARM                            ,&
             ORBIT     = ORBIT                            ,&
             TILELATS  = LATS                             ,&      ! [radians]
             TILELONS  = LONS                             ,&      ! [radians]
             INTERNAL_ESMF_STATE = INTERNAL               ,&
             RC=STATUS )
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
        call MAPL_GetPointer(IMPORT,FRZR   ,'FRZR'   ,RC=STATUS); VERIFY_(STATUS)
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
        IF (catchcn_internal%ATM_CO2 == 4) &
             call MAPL_GetPointer(IMPORT,CO2SC  ,'CO2SC'  ,RC=STATUS); VERIFY_(STATUS)
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

        if     (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND == 1) then
           
           call MAPL_GetPointer(INTERNAL,delCQ_delTVA ,'delCQ_delTVA'       ,RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,delCH_delTVA ,'delCH_delTVA'       ,RC=STATUS); VERIFY_(STATUS)
           
        elseif (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND >= 2) then 
           
           call MAPL_GetPointer(INTERNAL,delCH_delTC  ,'delCH_delTC'        ,RC=STATUS); VERIFY_(STATUS)
           call MAPL_GetPointer(INTERNAL,delCQ_delQC  ,'delCQ_delQC'        ,RC=STATUS); VERIFY_(STATUS)
           
        end if

        call MAPL_GetPointer(INTERNAL,TILE_ID    ,'TILE_ID'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,NDEP       ,'NDEP'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ABM        ,'ABM'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,PEATF      ,'PEATF'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,GDP        ,'GDP'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,HDM        ,'HDM'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,FIELDCAP   ,'FIELDCAP'   ,RC=STATUS); VERIFY_(STATUS)        
        call MAPL_GetPointer(INTERNAL,CLI_T2M    ,'CLI_T2M'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBVR    ,'BGALBVR'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBVF    ,'BGALBVF'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBNR    ,'BGALBNR'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BGALBNF    ,'BGALBNF'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CNCOL      ,'CNCOL'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CNPFT      ,'CNPFT'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TGWM       ,'TGWM'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RZMM       ,'RZMM'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SFMM       ,'SFMM'       ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,BFLOWM     ,'BFLOWM'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TOTWATM    ,'TOTWATM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TAIRM      ,'TAIRM'      ,RC=STATUS); VERIFY_(STATUS)  
        call MAPL_GetPointer(INTERNAL,RHM        ,'RHM'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,WINDM      ,'WINDM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RAINFM     ,'RAINFM'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNOWFM     ,'SNOWFM'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RUNSRFM    ,'RUNSRFM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,AR1M       ,'AR1M'       ,RC=STATUS); VERIFY_(STATUS)           
        call MAPL_GetPointer(INTERNAL,TPM        ,'TPM'        ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,CNSUM      ,'CNSUM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,PSNSUNM    ,'PSNSUNM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,PSNSHAM    ,'PSNSHAM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,LMRSUNM    ,'LMRSUNM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,LMRSHAM    ,'LMRSHAM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,LAISUNM    ,'LAISUNM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,LAISHAM    ,'LAISHAM'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNDZM      ,'SNDZM'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,SNDZM5D    ,'SNDZM5D'    ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ASNOWM     ,'ASNOWM'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,T2M10D     ,'T2M10D'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TG10D      ,'TG10D'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,T2MMIN5D   ,'T2MMIN5D'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RH30D      ,'RH30D'      ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TPREC10D   ,'TPREC10D'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,TPREC60D   ,'TPREC60D'   ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,ET365D     ,'ET365D'     ,RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(INTERNAL,RUNSURF    ,'RUNSURF'    ,RC=STATUS); VERIFY_(STATUS)
 
        if (catchcn_internal%N_CONST_LAND4SNWALB /= 0) then
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

        call MAPL_GetPointer(EXPORT,EVAPOUT            , 'EVAPOUT',ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SUBLIM             , 'SUBLIM' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SHOUT              , 'SHOUT'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RUNOFF             , 'RUNOFF' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPINT             , 'EVPINT' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPSOI             , 'EVPSOI' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPVEG             , 'EVPVEG' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPICE             , 'EVPICE' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WAT10CM            , 'WAT10CM',ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WATSOI             , 'WATSOI' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ICESOI             , 'ICESOI' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVPSNO             , 'EVPSNO'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,BFLOW              , 'BASEFLOW',ALLOC=.true.,          RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SMELT              , 'SMELT'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FICE1              , 'FICE1'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FICE2              , 'FICE2'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FICE3              , 'FICE3'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,HLWUP              , 'HLWUP'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWNDSRF            , 'SWNDSRF',ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWNDSRF            , 'LWNDSRF',ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,HLATN              , 'HLATN'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,QINFIL             , 'QINFIL' ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,AR1                , 'AR1'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,AR2                , 'AR2'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RZEQ               , 'RZEQ'   ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHFLX              , 'GHFLX'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSURF             , 'TPSURF'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSN1              , 'TPSNOW'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPUST              , 'TPUNST'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSAT              , 'TPSAT'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPWLT              , 'TPWLT'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ASNOW              , 'ASNOW'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SHSNOW             , 'SHSNOW'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,AVETSNOW           , 'AVETSNOW'            ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FRSAT              , 'FRSAT'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FRUST              , 'FRUST'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,FRWLT              , 'FRWLT'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP1                , 'TP1'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP2                , 'TP2'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP3                , 'TP3'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP4                , 'TP4'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP5                , 'TP5'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TP6                , 'TP6'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EMIS               , 'EMIS'   ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBVR              , 'ALBVR'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBVF              , 'ALBVF'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBNR              , 'ALBNR'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ALBNF              , 'ALBNF'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DELTS              , 'DELTS'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DELQS              , 'DELQS'  ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TST                , 'TST'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,QST                , 'QST'    ,ALLOC=.true.,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LST                , 'LST'                 ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WET1               , 'WET1'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WET2               , 'WET2'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WET3               , 'WET3'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WCSF               , 'WCSF'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WCRZ               , 'WCRZ'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WCPR               , 'WCPR'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,ACCUM              , 'ACCUM'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SNOMAS             , 'SNOWMASS'            ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SNOWDP             , 'SNOWDP'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,EVLAND             , 'EVLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PRLAND             , 'PRLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SNOLAND            , 'SNOLAND'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DRPARLAND          , 'DRPARLAND'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DFPARLAND          , 'DFPARLAND'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LHSNOW             , 'LHSNOW'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWNETSNOW1         , 'SWNETSNOW'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWUPSNOW           , 'LWUPSNOW'            ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWDNSNOW           , 'LWDNSNOW'            ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TCSORIG            , 'TCSORIG'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSN1IN            , 'TPSN1IN'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TPSN1OUT           , 'TPSN1OUT'            ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LHLAND             , 'LHLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SHLAND             , 'SHLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWLAND             , 'SWLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SWDOWNLAND         , 'SWDOWNLAND'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,LWLAND             , 'LWLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHLAND             , 'GHLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHSNOW             , 'GHSNOW'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,GHTSKIN            , 'GHTSKIN'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SMLAND             , 'SMLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TWLAND             , 'TWLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TELAND             , 'TELAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,TSLAND             , 'TSLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DWLAND             , 'DWLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DHLAND             , 'DHLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPLAND             , 'SPLAND'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPLH               , 'SPLH'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPWATR             , 'SPWATR'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SPSNOW             , 'SPSNOW'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNLAI              , 'CNLAI'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTLAI             , 'CNTLAI'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSAI              , 'CNSAI'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTOTC             , 'CNTOTC'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNVEGC             , 'CNVEGC'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFROOTC           , 'CNFROOTC'            ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNPP              , 'CNNPP'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNGPP              , 'CNGPP'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSR               , 'CNSR'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNAR               , 'CNAR'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNHR               , 'CNHR'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNEE              , 'CNNEE'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNXSMR             , 'CNXSMR'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNADD              , 'CNADD'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNLOSS             , 'CNLOSS'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNBURN             , 'CNBURN'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PARABS             , 'PARABS'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PARINC             , 'PARINC'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SCSAT              , 'SCSAT'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SCUNS              , 'SCUNS'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,BTRANT             , 'BTRANT'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,SIF                , 'SIF'                 ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNCO2              , 'CNCO2'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFIRE_CNT         , 'CNFIRE_CNT'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSOM_CLOSS        , 'CNSOM_CLOSS'         ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNDEPLOY          , 'CNNDEPLOY'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNDENIT            , 'CNDENIT'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSMINN_LEACHED    , 'CNSMINN_LEACHED'     ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSMINN            , 'CNSMINN'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFIRE_NLOSS       , 'CNFIRE_NLOSS'        ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNLEAFN            , 'CNLEAFN'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNLEAFC            , 'CNLEAFC'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNGROSS_NMIN       , 'CNGROSS_NMIN'        ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNET_NMIN         , 'CNNET_NMIN'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNFIX_TO_SMINN    , 'CNNFIX_TO_SMINN'     ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNACTUAL_IMMOB     , 'CNACTUAL_IMMOB'      ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFPG              , 'CNFPG'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFPI              , 'CNFPI'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSMINN_TO_PLANT   , 'CNSMINN_TO_PLANT'    ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNSMINN_TO_NPOOL   , 'CNSMINN_TO_NPOOL'    ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNNDEP_TO_SMINN    , 'CNNDEP_TO_SMINN'     ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTOTVEGN          , 'CNTOTVEGN'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTOTLITN          , 'CNTOTLITN'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTOTSOMN          , 'CNTOTSOMN'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNRETRANSN         , 'CNRETRANSN'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNRETRANSN_TO_NPOOL, 'CNRETRANSN_TO_NPOOL' ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNFUELC            , 'CNFUELC'             ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNTOTLITC          , 'CNTOTLITC'           ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNCWDC             , 'CNCWDC'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,CNROOT             , 'CNROOT'              ,           RC=STATUS); VERIFY_(STATUS)        
        call MAPL_GetPointer(EXPORT,CNFSEL             , 'CNFSEL'              ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU001         , 'RMELTDU001'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU002         , 'RMELTDU002'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU003         , 'RMELTDU003'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU004         , 'RMELTDU004'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTDU005         , 'RMELTDU005'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTBC001         , 'RMELTBC001'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTBC002         , 'RMELTBC002'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTOC001         , 'RMELTOC001'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,RMELTOC002         , 'RMELTOC002'          ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PEATCLSM_WATERLEVEL, 'PEATCLSM_WATERLEVEL' ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,PEATCLSM_FSWCHANGE , 'PEATCLSM_FSWCHANGE'  ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZGT1              , 'DZGT1'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZGT2              , 'DZGT2'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZGT3              , 'DZGT3'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZGT4              , 'DZGT4'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZGT5              , 'DZGT5'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZGT6              , 'DZGT6'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZPR               , 'DZPR'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZRZ               , 'DZRZ'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZSF               , 'DZSF'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,DZTS               , 'DZTS'                ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WPEMW              , 'WPEMW'               ,           RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,WPMC               , 'WPMC'                ,           RC=STATUS); VERIFY_(STATUS)

        NTILES = size(PS)
        
    allocate(   ityp(ntiles,nveg,nzone) )
    allocate(   fveg(ntiles,nveg,nzone) )
    allocate(   wtzone   (ntiles,nzone) )
    allocate(   elai(ntiles,nveg,nzone) )
    allocate(   esai(ntiles,nveg,nzone) )
    allocate(   tlai(ntiles,nveg,nzone) )

! initialize CN model and transfer restart variables on startup
! -------------------------------------------------------------
    if(first) then
    
      ! set number of time steps within a XX-day/hour period for 2m temperature XX-day/hour "running mean"
      ! --------------------------------------------------------------------------------------------------
      n1d  = 86400/dt
      n5d  = 5*86400/dt
      n10d = 10*86400/dt
      n30d = 30*86400/dt
      n60d = 60*86400/dt
      n365d = 365*86400/dt
      ! fzeng: this is done in such way to exclude istep in the restart file
      if(init_accum) then
        istep = 0                             ! set model time step index to 0 when begin to accumulate the cumulative variables, fzeng, 21 Apr 2017 
       else
        istep = maxval((/n10d,n30d,n60d/)) ! otherwise, set model time step index to the maximum of these nXX
      end if           
      
      ! jkolassa: implement this separately for 365-day running mean of ET
      if(init_accum_365) then 
        istep_365 = 0                             ! set model time step index to 0 when begin to accumulate the cumulative variables, fzeng, 21 Apr 2017 
       else
        istep_365 = maxval((/n10d,n30d,n60d,n365d/)) ! otherwise, set model time step index to the maximum of these nXX
      end if


      ! variables used for summing CN inputs over multiple land model calls; not saved on restart 
      ! fzeng: run must end on a CN call step
      ! -----------------------------------------------------------------------------------------
!      allocate( lmrsunm(ntiles,nveg,nzone) )
!      allocate( lmrsham(ntiles,nveg,nzone) )
!      allocate(             runsrf(ntiles) )      
!      
!      lmrsunm = 0.
!      lmrsham = 0.      
!      runsrf = 0.
           
      first = .false.
      
    endif

! set CLM CN PFT & fraction, set carbon zone weights
! --------------------------------------------------
    do nz = 1,nzone
      ityp(:,:,nz) = nint(ity(:,:))
      fveg(:,:,nz) = fvg(:,:)
      wtzone(:,nz) = CN_zone_weight(nz)
    end do

    ! obtain LAI from previous time step (from CN model)
    ! --------------------------------------------------
    call get_CN_LAI(ntiles,ityp,fveg,elai,esai=esai,tlai = tlai)

! OPTIONAL IMPOSE MONTHLY MEAN DIURNAL CYCLE FROM NOAA CARBON TRACKER
! -------------------------------------------------------------------

    IF ((catchcn_internal%ATM_CO2 == 1).OR.(catchcn_internal%ATM_CO2 == 2)) THEN
       READ_CT_CO2: IF(first_ct) THEN

          ! Carbon Tracker grid tiles mapping

          allocate (CT_INDEX  (1:CT_grid_N_lon, 1:CT_grid_N_lat))
          do j = 1, CT_grid_N_lat
             do i = 1, CT_grid_N_lon
                CT_INDEX (i,j) = (j - 1) * CT_grid_N_lon + i
             end do
          end do

          allocate (ct2cat (1:  NTILES))
          allocate (ct_tid (1:  NTILES))
          
          ct_tid = -9999
          ct2cat = 0

          do N = 1, NTILES             
             I =  NINT ((CEILING (lons(n)*90./MAPL_PI)*2 + 180.) / CT_grid_dlon)
             J =  NINT ((CEILING (lats(n)*90./MAPL_PI)*2 +  90.) / CT_grid_dlat)
             CT2CAT  (N) = ct_index (i,j)
          end do

          N = count(ct2cat > 0)
          
          allocate (unq_mask(1:N ))
          allocate (loc_int (1:N ))
          
          loc_int = pack(ct2cat ,mask = (ct2cat > 0))
          call MAPL_Sort (loc_int)

          unq_mask = .true.
          
          do i = 2,N 
             unq_mask(i) = .not.(loc_int(i) == loc_int(i-1)) 
          end do
          
          NUNQ = count(unq_mask)  
          
          allocate (ThisIndex (1:NUNQ))
          ThisIndex =  pack(loc_int, mask = unq_mask )

          do i = 1, NUNQ             
             where (ct2cat ==  ThisIndex(i)) ct_tid = i                
          end do
          
          ! Reading Carbon Tracker CO2_MonthlyMean_DiurnalCycle

          call ESMF_ClockGet( CLOCK, startTime=MODELSTART, RC=STATUS ); VERIFY_(STATUS)
          call ESMF_TimeGet ( MODELSTART, YY = FIRST_YY,  rc=status )  ; VERIFY_(STATUS)
          CALL ESMF_VMGet(vm, MPICOMMUNICATOR=comm, rc=status); VERIFY_(status)
          call MPI_Info_create(info, STATUS); VERIFY_(status)
          call MPI_Info_set(info, "romio_cb_read", "automatic", STATUS); VERIFY_(status)

          call MAPL_GetResource (MAPL, CO2_CycleFile, label = 'CO2_MonthlyMean_DiurnalCycle_FILE:', default = 'CO2_MonthlyMean_DiurnalCycle.nc4', RC=STATUS )
          VERIFY_(STATUS) 

          STATUS = NF_OPEN (trim(CO2_CycleFile), NF_NOWRITE, CTfile); VERIFY_(status) 

          allocate (CT_CO2V (1: NUNQ, 1:12, 1:8))
          allocate (CTCO2_TMP (1:CT_grid_N_lon, 1:CT_grid_N_lat, 1:12, 1:8))

          STATUS = NF_GET_VARA_REAL (CTfile, VarID(CTfile,'CO2'),  (/1,1,1,1/), &
               (/CT_grid_N_lon, CT_grid_N_lat, 12, 8/), CTCO2_TMP);VERIFY_(STATUS) 

          do N = 1, NUNQ
             I = MOD (ThisIndex(N), CT_grid_N_lon)
             IF(I == 0) I = CT_grid_N_lon
             J = (ThisIndex(N) -I) / CT_grid_N_lon + 1 

             CT_CO2V (N,:,:) = CTCO2_TMP (I,J,:,:)
             
          end do

          status = NF_CLOSE (CTFile); VERIFY_(status) 
          first_ct = .false.

          deallocate (CTCO2_TMP,ct2cat, unq_mask, loc_int, ct_index, ThisIndex)

       ENDIF READ_CT_CO2       
    ENDIF


        ! --------------------------------------------------------------------------
        ! ALLOCATE LOCAL POINTERS
        ! --------------------------------------------------------------------------

        allocate(GHTCNT  (N_GT,  NTILES))
        allocate(WESNN   (N_SNOW,NTILES))
        allocate(HTSNNN  (N_SNOW,NTILES))
        allocate(SNDZN   (N_SNOW,NTILES))
        allocate(FICESOUT(N_SNOW,NTILES))

        allocate(TILEZERO (NTILES))
        allocate(DZSF_in_mm(NTILES))
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
	allocate(RDC_TMP_1(NTILES))
        allocate(RDC_TMP_2(NTILES))
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
        allocate(LHOUT    (NTILES))
	allocate(EVACC    (NTILES))
        allocate(LHACC    (NTILES))
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
	allocate(fveg1     (NTILES))
	allocate(fveg2     (NTILES))
        allocate(FICE1TMP  (NTILES)) 
        allocate(SLDTOT    (NTILES))             ! total solid precip
        allocate(FSW_CHANGE(NTILES))

        allocate(SHSBT    (NTILES,NUM_SUBTILES))
        allocate(DSHSBT   (NTILES,NUM_SUBTILES))
        allocate(EVSBT    (NTILES,NUM_SUBTILES))
        allocate(DEVSBT   (NTILES,NUM_SUBTILES))
        allocate(DEDTC    (NTILES,NUM_SUBTILES))
        allocate(DHSDQC   (NTILES,NUM_SUBTILES))
        allocate(CFT      (NTILES,NUM_SUBTILES))
        allocate(CFQ      (NTILES,NUM_SUBTILES))
        allocate(TCO      (NTILES,NUM_SUBTILES))
        allocate(QCO      (NTILES,NUM_SUBTILES))
        allocate(DQS      (NTILES,NUM_SUBTILES))
        allocate(QSAT     (NTILES,NUM_SUBTILES))
        allocate(RA       (NTILES,NUM_SUBTILES))
        allocate(ALWN     (NTILES,NUM_SUBTILES))
        allocate(BLWN     (NTILES,NUM_SUBTILES))
        
        allocate(TC1_0    (NTILES))
        allocate(TC2_0    (NTILES))
        allocate(TC4_0    (NTILES))
        allocate(QA1_0    (NTILES))
        allocate(QA2_0    (NTILES))
        allocate(QA4_0    (NTILES))
        allocate(RCONSTIT (NTILES,N_SNOW,N_constit))
        allocate(TOTDEPOS (NTILES,N_constit))
        allocate(RMELT    (NTILES,N_constit))

        allocate(TA_MIN   (NTILES))
        allocate(ta_count (NTILES))

        call ESMF_VMGetCurrent ( VM, RC=STATUS )

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
        ! Offline land spin-up.
        ! --------------------------------------------------------------------------

        if (CATCHCN_INTERNAL%CATCH_SPINUP /= 0) then

           ! remove snow every Aug 1, 0z (Northern Hemisphere) or Feb 1, 0z (Southern Hemisphere)
           !
           ! assumes that CURRENT_TIME actually hits 0z on first of month (which seems safe enough)

           call ESMF_TimeGet(CURRENT_TIME, mm=CurrMonth, dd=CurrDay, h=CurrHour, m=CurrMin, s=CurrSec, rc=STATUS)
           VERIFY_(STATUS)

           if (CurrDay==1 .and. CurrHour==0 .and. CurrMin==0 .and. CurrSec==0) then

              if      (CurrMonth==8) then

                 where ( LATS >= 0. )    ! [radians]

                    WESNN1  = 0.
                    WESNN2  = 0.
                    WESNN3  = 0.
                    HTSNNN1 = 0.
                    HTSNNN2 = 0.
                    HTSNNN3 = 0.
                    SNDZN1  = 0.
                    SNDZN2  = 0.
                    SNDZN3  = 0.

                 end where

              else if (CurrMonth==2) then

                 where ( LATS <  0. )    ! [radians]

                    WESNN1  = 0.
                    WESNN2  = 0.
                    WESNN3  = 0.
                    HTSNNN1 = 0.
                    HTSNNN2 = 0.
                    HTSNNN3 = 0.
                    SNDZN1  = 0.
                    SNDZN2  = 0.
                    SNDZN3  = 0.

                 end where

              end if

           end if  ! 0z on first of month

        end if     ! if (CATCHCN_INTERNAL%CATCH_SPINUP /= 0)

        ! --------------------------------------------------------------------------
        ! Catchment Id and vegetation types used to index into tables
        ! --------------------------------------------------------------------------

        CAT_ID = nint(tile_id)

        where(ITY(:,1) > 0.)           
          VEG1 = map_cat(nint(ITY(:,1))) ! map  primary  CN PFT to catchment type
        endwhere
        where(ITY(:,2) > 0.)
          VEG2 = map_cat(nint(ITY(:,2))) ! map secondary CN PFT to catchment type
        endwhere

        fveg1(:) = fvg(:,1) 
        fveg2(:) = fvg(:,2)

        allocate ( lai1(ntiles) )
        allocate ( lai2(ntiles) )
        allocate ( wght(ntiles) )

        lai1 = 0.
        wght = 0.
        do nz = 1,nzone
          nv = 1
          lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
          wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
        end do
        lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

        lai2 = 0.
        wght = 0.
        do nz = 1,nzone
          nv = 2
          lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
          wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
        end do
        lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

! LAI seen by the land model
! --------------------------
        lai = fveg1*lai1 + fveg2*lai2 ! gkw: prognostic LAI on catch_internal_rst (overwrite VEGDYN import)

        ! --------------------------------------------------------------------------
        ! surface layer depth for soil moisture
        ! --------------------------------------------------------------------------
        
        DZSF_in_mm(:) = catchcn_internal%SURFLAY             ! same as DZSF but in units of [mm]

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

        if (CATCHCN_INTERNAL%Z0_FORMULATION == 4) then 
           ! make canopy height >= min veg height:
           Z2CH = max(Z2CH,MIN_VEG_HEIGHT)
           ZVG  = fveg1*(Z2CH - SAI4ZVG(VEG1)*SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI1)) + &
                  fveg2*(Z2CH - SAI4ZVG(VEG2)*SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI2))     
        else 
           ZVG  = fveg1*(Z2CH - SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI1)) + &
                  fveg2*(Z2CH - SCALE4ZVG*(Z2CH - MIN_VEG_HEIGHT)*exp(-LAI2))
        endif

        Z0   = Z0_BY_ZVEG*ZVG*SCALE4Z0_u
        IF (catchcn_internal%USE_ASCATZ0 == 1) WHERE (NDVI <= 0.2) Z0 = ASCATZ0        
        D0   = D0_BY_ZVEG*ZVG

        UUU = max(UU,MAPL_USMIN) * (log((ZVG-D0+Z0)/Z0) &
             / log((max(DZ-D0,10.)+Z0)/Z0))

        !--------------- GOSWIM IMPORTS FROM GOCART ---------------
        ! Initialization
        if (N_CONSTIT > 0) then 
           RCONSTIT(:,:,:)  = 0.0
           TOTDEPOS(:,:) = 0.0
           RMELT(:,:)  = 0.0
        endif
        !------------------------------------------------------------------

        ! Zero the light-absorbing aerosol (LAA) deposition rates from  GOCART:

        select case (catchcn_internal%AEROSOL_DEPOSITION)
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

        if (CATCHCN_INTERNAL%N_CONST_LAND4SNWALB /= 0) then
        
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
        
        ! old RDC formulation implemented in orginial GEOScatchCN_GridCom
        ! RDC = max(VGRDA(VEG1),VGRDA(VEG2))*min(1.,lai/2.)

        ! new RDC formulation used to reproduce Fanwei Zeng's LDASsa Catchment-CN.4.0 and Eunjee Lee's Catchment-CN.4.5 simulations
        rdc_tmp_1 = max( VGRDA(VEG1)*min( 1., LAI1/VGRDB(VEG1) ), 0.001)
        rdc_tmp_2 = max( VGRDA(VEG2)*min( 1., LAI2/VGRDB(VEG2) ), 0.001)
        RDC = max(rdc_tmp_1,rdc_tmp_2)*min(1.,lai/2.)
        RDC = max(RDC,0.001)

        RHO = PS/(MAPL_RGAS*(TA*(1+MAPL_VIREPS*QA)))

        !--------------------------------------------------------------------------------------------------------
        !                                                      
        ! MOSFC variable names and description:                
        !                                                      
        !--------------------------------------------------------------------------------------------------------
        ! GEOS_CatchGridComp.F90   | catchment.F90 | dimension  | Description
        !--------------------------------------------------------------------------------------------------------        
        !  TA                      |     TM        | NT         | surface (lowest model level) air temperature
        !  QA                      |     QM        | NT         | surface (lowest model level) air spec humidity
        !--------------------------------------------------------------------------------------------------------
        !  TC                      |     TC        | NT-by-NSBT | canopy (air) temperature
        !  QC                      |     QA (!)    | NT-by-NSBT | canopy (air) specific humidity
        !  CH                      |     -         | NT-by-NSBT | exchange coeff for heat
        !  CQ                      |     -         | NT-by-NSBT | exchange coeff for humidity
        !  EVSBT                   |     ETURB     | NT-by-NSBT | evaporation
        ! DEVSBT                   |    DEDQA      | NT-by-NSBT | deriv of evap w.r.t. canopy spec humidity
        ! DEDTC                    |    DEDTC      | NT-by-NSBT | deriv of evap w.r.t. canopy temperature
        !  SHSBT                   |     HSTURB    | NT-by-NSBT | sensible heat flux (SH)
        ! DHSDQC (formerly DHSDQA) |    DHSDQA     | NT-by-NSBT | deriv of SH   w.r.t. canopy spec humidity 
        ! DSHSBT                   |    DHSDTC     | NT-by-NSBT | deriv of SH   w.r.t. canopy temperature 
        !--------------------------------------------------------------------------------------------------------
        !   *SBT = sub-tile (?)
        !   NT   = number of tiles
        !   NSBT = number of subtiles (per tile)
        !
        ! For land, CH = CQ in Helfand and Louis.
        !
        !
        ! MOSFC equations:
        !
        !  EVSBT =      CQ * (QC - QA)
        !  SHSBT = Cp * CH * (TC - TA)           [ Cp  = MAPL_CP ]
        !
        ! Derivatives obtained via product rule.  See equations below.
        !
        ! For analytical derivatives, additionally use the following identities:
        !
        !    virtual TC: TVC = TC*(1 + eps)*QC   [ eps = MAPL_VIREPS ]
        !    virtual TA: TVA = TA*(1 + eps)*QA
        !
        !    delTVC_delQC =  TC*eps
        !    delTVC_delTC = (1 + eps)
        !
        !    delCQ_delQC = delCQ_delTVC * delTVC_delQC 
        !    delCH_delTC = delCH_delTVC * delTVC_delTC
        !
        !    CQ=CQ(Ri) where Ri is proportional to deltaTVA=TVA-TVC --> produces a minus sign
        !
        !    delCQ_delTVC = -1*delCQ_delTVA
        !    delCH_delTVC = -1*delCH_delTVA
        !
        !--------------------------------------------------------------------------------------------------
        ! reichle, 9/9/2024
        !--------------------------------------------------------------------------------------------------
        
        ! initialize derivatives that may not be filled later
        
        DEDTC =0.0
        DHSDQC=0.0

        if(OFFLINE_MODE /=0) then

           ! Catchment in offline (land-only) mode
        
           do N=1,NUM_SUBTILES

              CFT   (:,N) = 1.0
              CFQ   (:,N) = 1.0

              SHSBT (:,N) = MAPL_CP*CH(:,N)*(TC(:,N)-TA)
              EVSBT (:,N) = CQ(:,N)*(QC(:,N)-QA)

              BLWN(:,N) = EMIS*MAPL_STFBOL*TC(:,N)*TC(:,N)*TC(:,N)
              ALWN(:,N) = -3.0*BLWN(:,N)*TC(:,N)
              BLWN(:,N) =  4.0*BLWN(:,N)

           end do
           
           select case (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND)
              
           case (0)    ! ignore derivatives of exchange coeffs w.r.t. canopy temp and specific humidity
              
              do N=1,NUM_SUBTILES
                 DEVSBT(:,N) =           CQ(:,N)
                 DSHSBT(:,N) = MAPL_CP*  CH(:,N)
              end do
              
           case (1)    ! Louis only: analytical derivatives of exchange coeffs w.r.t. canopy temp and specific humidity
              
              _ASSERT( CATCHCN_INTERNAL%CHOOSEMOSFC==0, 'must use Louis scheme for MOSFC analytical derivatives' )
              
              do N=1,NUM_SUBTILES
                 DEVSBT(:,N) =           CQ(:,N) + max( 0.0,         -delCQ_delTVA(:,N)*    MAPL_VIREPS*TC(:,N) *(QC(:,N)-QA) )
                 DEDTC( :,N) =                     max( 0.0,         -delCQ_delTVA(:,N)*(1.+MAPL_VIREPS*QC(:,N))*(QC(:,N)-QA) )
                 DSHSBT(:,N) = MAPL_CP*( CH(:,N) + max( 0.0,         -delCH_delTVA(:,N)*(1.+MAPL_VIREPS*QC(:,N))*(TC(:,N)-TA) ) )
                 DHSDQC(:,N) =                     max( 0.0, -MAPL_CP*delCH_delTVA(:,N)*    MAPL_VIREPS*TC(:,N) *(TC(:,N)-TA) )
              end do
              
           case (2,3)  ! numerical derivatives of exchange coeffs w.r.t. canopy temp and specific humidity
              
              do N=1,NUM_SUBTILES
                 DEVSBT(:,N) =           CQ(:,N) + max( 0.0,          delCQ_delQC( :,N)*                         (QC(:,N)-QA) )
                 DEDTC( :,N) =                     max( 0.0,          delCH_delTC( :,N)*                         (QC(:,N)-QA) )
                 DSHSBT(:,N) = MAPL_CP*( CH(:,N) + max( 0.0,          delCH_delTC( :,N)*                         (TC(:,N)-TA) ) )   
                 DHSDQC(:,N) =                     max( 0.0,  MAPL_CP*delCQ_delQC( :,N)*                         (TC(:,N)-TA) )
              end do
              
           case default
              
              _ASSERT(.false., 'unknown MOSFC_EXTRA_DERIVS_OFFL_LAND')
              
           end select
                      
        else

           ! GCM: Catchment coupled to atmosphere
           
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

        end if   ! Catchment offline

        ! Compute DQS; make sure QC is between QA and QSAT; compute RA.
        !
        !   Some 1,000 lines below, duplicate code was present and removed in Jan 2022. 
        !   - reichle, 14 Jan 2022.
        !
        ! reichle, 9/9/2024:
        !
        ! WHY IS QC RESET *AFTER* IT WAS USED TO CALCULATE THE EXCHANGE COEFFS (AND DERIVS) ABOVE???
        !
        ! For reference, the following comment was copied from from LDASsa m3-16, specifically, from 
        ! reichle-LDASsa_m3-16_6/src/Components/GEOSlana_GridComp/process_cat.F90 (Lines 330-333):
        !
        ! ! compute surface exchange coefficients etc BEFORE possibly resetting
        ! ! profile of Qair-QAx-Qsat(surf) -- for consistency with two-stage 
        ! ! run-method in GEOS_CatchGridComp.F90
        ! ! reichle+qliu,  9 Oct 2008

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

        SLDTOT = SNO+ICE      ! do *not* add FRZR (freezing rain) to solid precip, see comment below

        ! FRZR (freezing rain) is rain that falls as super-cooled liquid water, which freezes upon
        !      impact on a sufficiently cold surface.  As such, FRZR is *not* solid precipitation
        !      and should be considered rainfall.
        !
        ! As of Jun 2025, FRZR is identical to 0 and can be ignored.  Looking ahead, make sure to
        !      account correctly for FRZR in the input precipitation variables.  Once it's filled
        !      with non-zero values, FRZR will (probably) be included in PLS+PCU.  It is (probably)
        !      better to replace PCU & PLS with RAIN and FRZR, where RAIN (probably) does *not*
        !      include FRZR and PCU+PLS=RAIN+FRZR (TO BE CONFIRMED!).
        !
        ! - reichle, 6/6/2025

	! --------------------------------------------------------------------------
	! protect the forcing from unsavory values, as per practice in offline
	! driver
	! --------------------------------------------------------------------------

        _ASSERT(count(PLS<0.)   ==0, 'encountered neg precip value (PLS)'   )
        _ASSERT(count(PCU<0.)   ==0, 'encountered neg precip value (PCU)'   )
        _ASSERT(count(SLDTOT<0.)==0, 'encountered neg precip value (SLDTOT)')

        LAI0  = max(0.0001     , LAI)
        GRN0  = max(0.0001     , GRN)		
        ZTH   = max(0.0001     , ZTH)

        TCO   = TC
        QCO   = QC

        ! --------------------------------------------------------------------------
        ! actual CATCHMENT call
        ! --------------------------------------------------------------------------

        TILEZERO = 0.0

        call MAPL_TimerOn  ( MAPL, "-CATCHCNCLM51" )


! ----------------------------------------------------------------------------------------

! gkw: start on main CN block

    allocate(   btran(ntiles,nveg,nzone) )
    allocate(   btran_fire(ntiles,nzone) )
    allocate(     wgt(ntiles) )
    allocate(     wpp(ntiles) )
    allocate(    fwet(ntiles) )
    allocate(  wet_in(ntiles) )
    allocate(     bt(ntiles,fsat:fwlt))
    allocate(     sm(ntiles,fsat:fwlt))
    allocate(  SWSRF1(ntiles) )
    allocate(  SWSRF2(ntiles) )
    allocate(  SWSRF4(ntiles) )
    allocate(     tcx(ntiles,nzone) )
    allocate(     qax(ntiles,nzone) )
    allocate(   rcxdt(ntiles) )
    allocate(   rcxdq(ntiles) )
    allocate(    car1(ntiles) )
    allocate(    car2(ntiles) )
    allocate(    car4(ntiles) )
    allocate( parzone(ntiles,nveg,nzone) )
    allocate(    para(ntiles) )
    allocate (  totwat(ntiles) )
    if(.not. allocated(npp )) allocate(     npp(ntiles) )
    if(.not. allocated(gpp )) allocate(     gpp(ntiles) )
    if(.not. allocated(sr  )) allocate(      sr(ntiles) )
    if(.not. allocated(aresp)) allocate(  aresp(ntiles) )
    if(.not. allocated(hresp)) allocate(  hresp(ntiles) )
    if(.not. allocated(nee )) allocate(     nee(ntiles) )
    if(.not. allocated(padd)) allocate(    padd(ntiles) )
    if(.not. allocated(frootc)) allocate(frootc(ntiles) )
    if(.not. allocated(vegc)) allocate(    vegc(ntiles) )
    if(.not. allocated(xsmr)) allocate(    xsmr(ntiles) )
    if(.not. allocated(burn)) allocate(    burn(ntiles) )
    if(.not. allocated(closs))allocate(   closs(ntiles) )

    allocate(             nfire(ntiles) )
    allocate(         som_closs(ntiles) )
    allocate(    dayl(ntiles) )
    allocate(dayl_fac(ntiles) )
    allocate(CO2V    (ntiles) )
    allocate(             fsnow(ntiles) )
    allocate(          ityp_tmp(ntiles) )
    allocate(     Qair_relative(ntiles) )    
    allocate(           ndeploy(ntiles) )
    allocate(             denit(ntiles) )
    allocate(     sminn_leached(ntiles) )
    allocate(             sminn(ntiles) )
    allocate(        fire_nloss(ntiles) )
    allocate(             leafn(ntiles) )
    allocate(             leafc(ntiles) )
    allocate(        gross_nmin(ntiles) )
    allocate(          net_nmin(ntiles) )
    allocate(     nfix_to_sminn(ntiles) )
    allocate(      actual_immob(ntiles) )
    allocate(               fpg(ntiles) )
    allocate(               fpi(ntiles) )
    allocate(    sminn_to_plant(ntiles) )
    allocate(    sminn_to_npool(ntiles) )
    allocate(     ndep_to_sminn(ntiles) )
    allocate(           totvegn(ntiles) )
    allocate(           totlitn(ntiles) )
    allocate(           totsomn(ntiles) )
    allocate(          retransn(ntiles) )
    allocate( retransn_to_npool(ntiles) )
    allocate(             fuelc(ntiles) )
    allocate(           totlitc(ntiles) )
    allocate(              cwdc(ntiles) )
    allocate(             rootc(ntiles) )
    allocate(              lnfm(ntiles) )

    allocate(     tgw(ntiles,nzone) )
    allocate(     rzm(ntiles,nzone) )
    allocate(    rc00(ntiles,nzone) )
    allocate(    rcdt(ntiles,nzone) )
    allocate(    rcdq(ntiles,nzone) )
    allocate( totcolc(ntiles,nzone) )
    allocate(         sfm(ntiles,nzone) )    

    allocate(       albdir(ntiles,nveg,nzone,2) )
    allocate(       albdif(ntiles,nveg,nzone,2) )

    allocate( psnsun(ntiles,nveg,nzone) )
    allocate( psnsha(ntiles,nveg,nzone) )
    allocate( laisun(ntiles,nveg,nzone) )
    allocate( laisha(ntiles,nveg,nzone) )
    allocate( lmrsun(ntiles,nveg,nzone) )
    allocate( lmrsha(ntiles,nveg,nzone) )
    allocate(      ht(N_gt) )
    allocate(      tp(N_gt) )
    allocate( soilice(N_gt) )

! get current date & time gkw: this is used to transfer CN restart vars & set declination
! -----------------------
    call ESMF_TimeGet  ( CURRENT_TIME, YY = AGCM_YY,       &
				       MM = AGCM_MM,       &
				       DD = AGCM_DD,       &
                                       H  = AGCM_HH,       &
                                       M  = AGCM_MI,       &
				       S  = AGCM_S ,       &  
				       dayOfYear = dofyr , &
				       rc=status )
    VERIFY_(STATUS)

    AGCM_S_ofday = AGCM_S + 60 * AGCM_MI + 3600 * AGCM_HH

! get ending time; determine if this is last call before ending time
! ------------------------------------------------------------------
    call ESMF_ClockGet ( clock,  StopTime=StopTime ,rc=STATUS )
    VERIFY_(STATUS)

    NextTime = CURRENT_TIME + DELT

    ! 0-land tiles processors hang in MAPL_ReadForcing
    ! Thus moved reading lnfm here
    ! ------------------------------------------------
    
    if(mod(AGCM_S_ofday,nint(catchcn_internal%DTCN)) == 0) then
        ! Get lightening frequency clim file name from configuration
       call MAPL_GetResource ( MAPL, LNFMFILE, label = 'LNFM_FILE:', default = 'lnfm.dat', RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_ReadForcing(MAPL,'LNFM',LNFMFILE,CURRENT_TIME,lnfm,ON_TILES=.true.,RC=STATUS)
       VERIFY_(STATUS)
    endif

    if(ntiles > 0) then ! gkw: skip threads with no land tiles

! gkw: assign new vegetation types and fractions
! ----------------------------------------------
    cat_id = nint(tile_id) ! gkw: temporary for debugging

! compute daylength (and daylength factor)
! ----------------------------------------

    ! current daylight duration
    call MAPL_SunGetDaylightDuration(ORBIT,lats,dayl,currTime=CURRENT_TIME,RC=STATUS)
    VERIFY_(STATUS)
    ! maximum daylight duration (at solstice)
    call MAPL_SunGetDaylightDurationMax(ORBIT,lats,dayl_fac,currTime=CURRENT_TIME,RC=STATUS)
    VERIFY_(STATUS)
    ! dayl_fac is ratio current:maximum dayl squared (min 0.01 [gkw: from CLM4])
    dayl_fac = min(1.,max(0.01,(dayl/dayl_fac)**2))

! gkw: obtain catchment area fractions and soil moisture
! ------------------------------------------------------
    call catch_calc_soil_moist( ntiles, DZSF_in_mm, vgwmax, cdcr1, cdcr2, psis, bee, poros, wpwet,&
         ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1, arw2, arw3, arw4, bf1, bf2,              &
         srfexc, rzexc, catdef, car1, car2, car4, sfmc, rzmc, prmc,                               &
         SWSRF1OUT=SWSRF1, SWSRF2OUT=SWSRF2, SWSRF4OUT=SWSRF4 )
    
                            
                              
! obtain saturated canopy resistance following Farquhar, CLM4 implementation    

! compute RC & PSN in each of the CN zones
! ----------------------------------------

! "btran" in the catchment zones; map into CN zones
! -------------------------------------------------
       sm(:,fsat) = 1.0
       
! gkw: bt2 is unstressed region only (subtract saturated and wilting areas)
       do n = 1,ntiles
         if(car2(n) > 0.) then
           sm(n,ftrn)=(rzmc(n)/poros(n) - car1(n) - car4(n)*wpwet(n))/car2(n)
          else
           sm(n,ftrn)= rzmc(n)/poros(n)
         endif
         sm(n,ftrn) = max(sm(n,ftrn),wpwet(n))
         sm(n,ftrn) = min(sm(n,ftrn),1.)

         if(car4(n) > 0.) then
           sm(n,fwlt)=(rzmc(n)/poros(n) - car1(n) - car2(n)*sm(n,ftrn))/car4(n)
          else
           sm(n,fwlt)= wpwet(n)
         endif
         sm(n,fwlt) = max(sm(n,fwlt),1.e-3)
         sm(n,fwlt) = min(sm(n,fwlt),wpwet(n)-1.e-7)
       end do
       
       bt(:,fsat) = 1.0
       bt(:,ftrn) = sm(:,ftrn)**(-bee)
       wpp = wpwet ** (-bee)
       bt(:,ftrn) = (bt(:,ftrn)-wpp)/(1.-wpp)
       bt(:,fwlt) = 0.
 
       do n = 1,ntiles

         ax1 = car1(n)
	 ax2 = car2(n)
         ax4 = 1. - ax1 - ax2

	 cn1 = wtzone(n,1)
	 cn2 = wtzone(n,2)
	 cn3 = wtzone(n,3)

! CN zone 1
         if(ax1 .gt. cn1) then
           f1(1) = cn1 ; f2(1) = 0. ; f4(1) = 0.
          else
           if((ax1+ax2) .gt. cn1) then
             f1(1) = ax1 ; f2(1) = cn1-ax1 ; f4(1) = 0.
            else
             f1(1) = ax1 ; f2(1) = ax2 ; f4(1) = cn1-ax1-ax2
           endif
         endif

! CN zone 2
         if(ax1 .gt. cn1) then
           cn12 = cn1 + cn2
           if(car1(n) .gt. cn12) then
             f1(2) = cn2 ; f2(2) = 0. ; f4(2) = 0.
            else
             if((ax1+ax2) .lt. cn12) then
               f1(2) = ax1-cn1 ; f2(2) = ax2 ; f4(2) = cn12-ax1-ax2
              else
               f1(2) = ax1-cn1 ; f2(2) = cn12-ax1 ; f4(2) = 0.
             endif
           endif
          else
           cn23 = cn2 + cn3
           if(ax4 .gt. cn23) then
             f1(2) = 0. ; f2(2) = 0. ; f4(2) = cn2
            else
             if(ax4 .lt. cn3) then
               f1(2) = 0. ; f2(2) = cn2 ; f4(2) = 0.
              else
               f1(2) = 0. ; f2(2) = cn23-ax4 ; f4(2) = ax4-cn3
             endif
           endif
         endif

! CN zone 3
         if(ax4 .gt. cn3) then
           f1(3) = 0. ; f2(3) = 0. ; f4(3) = cn3
          else
           if((ax4+ax2) .gt. cn3) then
             f1(3) = 0. ; f2(3) = cn3-ax4 ; f4(3) = ax4
            else
             f1(3) = cn3-ax4-ax2 ; f2(3) = ax2 ; f4(3) = ax4
           endif
         endif

         do nz = 1,nzone
            btran_fire(n,nz) = (f1(nz)*bt(n,fsat)     + f2(nz)*bt(n,ftrn)     + f4(nz)*bt(n,fwlt)    )/wtzone(n,nz)
	    tgw(n,nz) =  (f1(nz)*tg(n,fsat) + f2(nz)*tg(n,ftrn) + f4(nz)*tg(n,fwlt))/wtzone(n,nz)
	    tcx(n,nz) =  (f1(nz)*tc(n,fsat) + f2(nz)*tc(n,ftrn) + f4(nz)*tc(n,fwlt))/wtzone(n,nz)
	    qax(n,nz) =  (f1(nz)*qc(n,fsat) + f2(nz)*qc(n,ftrn) + f4(nz)*qc(n,fwlt))/wtzone(n,nz)
            rzm(n,nz) =  (f1(nz)*sm(n,fsat) + f2(nz)*sm(n,ftrn) + f4(nz)*sm(n,fwlt))/wtzone(n,nz)
            sfm(n,nz) =  (f1(nz)*SWSRF1(n)  + f2(nz)*SWSRF2(n)  + f4(nz)*SWSRF4(n) )/wtzone(n,nz)
         end do
         
       end do !n 

! soil temperature and hydrologic state
! -------------------------------------
    DO N=1,ntiles

! soil temperatures
! -----------------
       
      ! zbar function - reichle, 29 Jan 2022 (minus sign applied in call to GNDTMP)
      ZBAR = catch_calc_zbar( bf1(n), bf2(n), catdef(n) )  
      HT(:)=GHTCNT(:,N)
      CALL GNDTMP(poros(n),-1.*zbar,ht,frice,tp,soilice)  ! note minus sign for zbar

      ! At the CatchCNGridComp level, tp1, tp2, .., tp6 are export variables in units of Kelvin,
      ! - rreichle & borescan, 6 Nov 2020

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
      bflow(n) = (1.-frice)*1000.* &
    	    cond(n)*exp(-(bf3(n)-ashift)-gnu(n)*zbar)/gnu(n)
      IF(catdef(n) >= cdcr1(n)) bflow(n) = 0.
      bflow(n) = min(cond(n),bflow(n))
    end do

! compute relative humidity (%) used in CNFireMod
! -----------------------------------------------
    do n = 1,ntiles
       Qair_sat = MAPL_EQsat(TA(n), PS(n) )
       Qair_relative(n) = QA(n) / Qair_sat * 100. 
    end do

    Qair_relative(:) = min(max(0., Qair_relative(:)), 100.)

! compute accumulated fields, fzeng
! following the methods in accFldsMod.F90 and accumulMod.F90 in CLM4.5
! --------------------------------------------------------------------

    istep = istep + 1
    istep_365 = istep_365 + 1
    TA_MIN(:) = 1000.

    ! running mean - reset accumulation period until greater than nstep
    ! fzeng & gkw: may not be exactly 2m, but it is consistent with t_ref2m in CN model
    ! T2M10D   (T10    in CLM4.5): 10-day running mean of 2-m temperature (K)
    ! TPREC10D (PREC10 in CLM4.5): 10-day running mean of total precipitation (mm H2O/s)
    ! TPREC60D (PREC60 in CLM4.5): 60-day running mean of total precipitation (mm H2O/s)
    ! ---------------------------------------------------------------------------------   
     if(init_accum) then     
        
        ! (1) 5-day running mean of snow depth      
         accper = min(istep,n5d)
         SNDZM5D   = ((accper-1)*SNDZM5D + SNDZM) / accper 

        ! (1) 10-day running mean of 2-m temperature (K) and total precipitation (mm H2O/s)     
         accper = min(istep,n10d)
         T2M10D   = ((accper-1)*T2M10D + TA) / accper
         TPREC10D = ((accper-1)*TPREC10D + PCU + PLS + SNO) / accper      
         TG10D    = ((accper-1)*TG10D + TG(:,1)) / accper         

        ! (2) 30-day running mean of relative humidity [%]     
         accper = min(istep,n30d)
         RH30D   = ((accper-1)*RH30D + Qair_relative) / accper


        ! (2) 60-day running mean of total precipitation (mm H2O/s)
         accper = min(istep,n60d)
         TPREC60D = ((accper-1)*TPREC60D + PCU + PLS + SNO) / accper      


        ! jkolassa: for T2MMIN5D compute minimum T2M once per day, then use that value to compute new 5-day running mean of minimum T2M

         do n = 1,ntiles
            ta_count(n) = ta_count(n) + 1
            TA_MIN(n) = min(TA_MIN(n),TA(n))

            if (ta_count(n) == n1d) then 
                T2MMIN5D(n)   = ((accper-1)*T2MMIN5D(n)   + TA_MIN(n))  / accper
                TA_MIN(n) = 1000.
                ta_count(n) = 0
            end if 
         end do

      else
 
         SNDZM5D  = ((n5d-1)*SNDZM5D   + SNDZM)  / n5d
         T2M10D   = ((n10d-1)*T2M10D   + TA)  / n10d
         TG10D    = ((n10d-1)*TG10D   + TG(:,1))  / n10d
         TPREC10D = ((n10d-1)*TPREC10D + PCU + PLS + SNO) / n10d
         RH30D    = ((n30d-1)*RH30D    + Qair_relative)  / n30d
         TPREC60D = ((n60d-1)*TPREC60D + PCU + PLS + SNO) / n60d


         ! jkolassa: for T2MMIN5D compute minimum T2M once per day, then use that value to compute new 5-day running mean of minimum T2M

         do n = 1,ntiles
            ta_count(n) = ta_count(n) + 1
            TA_MIN(n) = min(TA_MIN(n),TA(n))

            if (ta_count(n) == n1d) then
                T2MMIN5D(n)   = ((accper-1)*T2MMIN5D(n)   + TA_MIN(n))  / accper
                TA_MIN(n) = 1000.
                ta_count(n) = 0
            end if 
         end do 

     endif


! get CO2
! -------

    if(catchcn_internal%ATM_CO2 == 3) catchcn_internal%CO2 = GETCO2(AGCM_YY,dofyr)

    CO2V (:) = catchcn_internal%CO2
    
! use CO2SC from GOCART/CO2
! -------------------------
    
    IF (catchcn_internal%ATM_CO2 == 4) THEN
       
       where ((CO2SC >= 0.) .and. (CO2SC <= 1000.))
          CO2V = CO2SC * 1e-6
       end where
       
    endif

    IF(catchcn_internal%ATM_CO2 == 1) co2g = 1.    ! DO NOT SCALE USE CT CLIMATOLOGY

    CALC_CTCO2_SF: IF(catchcn_internal%ATM_CO2 == 2)  THEN

       ! Compute scale factor to scale CarbonTracker CO2 monthly mean diurnal cycle (3-hourly)
       CO2_YEAR = AGCM_YY
       IF(catchcn_internal%CO2_YEAR_IN > 0) CO2_YEAR = catchcn_internal%CO2_YEAR_IN

       ! update EEA global average CO2 and co2 scalar at the beginning of each year, fz, 26 Sep 2016
       ! -------------------------------------------------------------------------------------------

       IF (AGCM_YY /= CO2_YEAR) CO2_YEAR = CO2_YEAR + AGCM_YY - FIRST_YY
          
       if (CO2_YEAR  < byr_co2g) then
          co2g = co2g_byr
       elseif ((CO2_YEAR >= byr_co2g).AND.(CO2_YEAR <= myr_co2g)) then
          co2g = co2g_byr + dco2g_1 * (CO2_YEAR - byr_co2g)
       else
          co2g = co2g_myr + dco2g_2 * (CO2_YEAR - myr_co2g)
       endif
     
       co2g = co2g / CTco2g  ! = co2g/CTco2g, is used to scale CarbonTracker CO2 monthly mean diurnal cycle (3-hourly)
     
    ENDIF CALC_CTCO2_SF

      USE_CT_CO2: IF((catchcn_internal%ATM_CO2 == 1).OR.(catchcn_internal%ATM_CO2 == 2)) THEN

         IF(AGCM_DD < 16) THEN

            ! interpolate between AGCM_MM - 1 and  AGCM_MM

            M1 = AGCM_MM -1
            Y1 = AGCM_YY
            if(M1 == 0) then ; M1 = 12 ; Y1 = AGCM_YY -1 ; endif

            call ESMF_TimeSet(BEFORE, YY = Y1, MM = M1, DD = 16, &
                 H =  0, M =  0, S = 0, rc = STATUS) ; VERIFY_(STATUS)
            call ESMF_TimeSet(AFTER , YY = AGCM_YY, MM = AGCM_MM,   DD = 15, &
                 H = 23, M = 59, S = 59, rc = STATUS); VERIFY_(STATUS)
         
            call MAPL_Interp_Fac (CURRENT_TIME,BEFORE,AFTER,FAC,RC=STATUS ) ; VERIFY_(STATUS)
            ASSERT_(FAC >= 0.0)
            ASSERT_(FAC <= 1.0)        

            DO N = 1,NTILES
               CO2V (N) = (FAC * CT_CO2V (CT_TID (N),M1, AGCM_HH/3+1)  +  (1.0-FAC) * CT_CO2V (CT_TID (N),AGCM_MM, AGCM_HH/3+1)) * &
                 CO2G * 1.e-6 ! scale by EEA global average CO2 * convert from ppm
            END DO
         ELSE
            
            ! interpolate between AGCM_MM and  AGCM_MM + 1

            M1 = AGCM_MM +1
            Y1 = AGCM_YY
            if(M1 == 13) then ; M1 = 1 ; Y1 = AGCM_YY +1 ; endif

            call ESMF_TimeSet(BEFORE , YY = AGCM_YY, MM = AGCM_MM,   DD = 16, &
                 H =  0, M =  0, S =  0, rc = STATUS) ; VERIFY_(STATUS)           
            call ESMF_TimeSet(AFTER, YY = Y1, MM = M1, DD = 15, &
                 H = 23, M = 59, S = 59, rc = STATUS) ; VERIFY_(STATUS)
         
            call MAPL_Interp_Fac (CURRENT_TIME,BEFORE,AFTER,FAC,RC=STATUS ) ; VERIFY_(STATUS)
            ASSERT_(FAC >= 0.0)
            ASSERT_(FAC <= 1.0)        
            DO N = 1,NTILES
               CO2V (N) = (FAC * CT_CO2V (CT_TID (N),AGCM_MM, AGCM_HH/3+1)  +  (1.0-FAC) * CT_CO2V (CT_TID (N),M1 , AGCM_HH/3+1)) * &
                 CO2G * 1.e-6 ! scale by EEA global average CO2 * convert from ppm
            END DO
         ENDIF
          
      ENDIF USE_CT_CO2

    if(associated(BTRANT)) btrant = 0. 

! fraction of foliage that is wet  gkw 20140327
! -------------------------------
    do n = 1,ntiles
      if(lai(n) > 1.e-4) then 
        fwet(n) = min(1.,max(0.,capac(n)/(0.2*lai(n))))
       else
        fwet(n) = 0. 
      endif
    end do

! compute snow-free albedo for each PFT in each zone  gkw: assume the snow albedo is not very important
! --------------------------------------------------
     do nz = 1,nzone
      do nv = 1,nveg
        ityp_tmp(:) = map_cat(ityp(:,nv,nz))

! fzeng: note that this is not exactly the same as calling sibalb_vis in the unified model because the 
! "if(fveg(i)>1.e-4 .and. zth(i)>0.01)" branch in subroutine sibalb_vis is absent in the current subroutine sibalb.
! -----------------------------------------------------------------------------------------------------------------

        call SIBALB(ntiles, ityp_tmp, elai(:,nv,nz), GRN, ZTH,            &
                    BGALBVR, BGALBVF, BGALBNR, BGALBNF,                   & ! gkw: MODIS soil background albedo
                    ALBVR, ALBNR, ALBVF, ALBNF, MODIS_SCALE=.TRUE.  )       ! instantaneous snow-free albedos on tiles

        ! Get TPSN1OUT1 for SNOW_ALBEDO parameterization
        
        call STIEGLITZSNOW_CALC_TPSNOW(NTILES, HTSNNN(1,:), WESNN(1,:), TPSN1OUT1, FICE1TMP )    
        TPSN1OUT1 =  TPSN1OUT1 + MAPL_TICE
        
        call StieglitzSnow_snow_albedo(ntiles, N_snow, catchcn_internal%N_CONST_LAND4SNWALB, ityp_tmp,   &
                    elai(:,nv,nz), ZTH,                                      &
                    RHOFS,                                                &   
                    SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                  & 
                    WESNN, HTSNNN, SNDZN,                                 &
                    ALBVR, ALBNR, ALBVF, ALBNF,                           & ! instantaneous snow-free albedos on tiles
                    SNOVR, SNONR, SNOVF, SNONF,                           & ! instantaneous snow albedos on tiles
                    RCONSTIT, UUU, TPSN1OUT1, DRPAR, DFPAR) 

! fsnow: pft-level; asnow: grid-level
! -----------------------------------
        where(tlai(:,nv,nz) > 0.)
          fsnow(:) = 1. - elai(:,nv,nz)/tlai(:,nv,nz)
          fsnow(:) = min(max(fsnow(:),0.),1.)
         elsewhere
          fsnow(:) = 0.
        endwhere

        ! visible
        albdir(:,nv,nz,1) = albvr(:)*(1.-fsnow(:)) + snovr(:)*fsnow(:)
        albdif(:,nv,nz,1) = albvf(:)*(1.-fsnow(:)) + snovf(:)*fsnow(:)
        
        ! NIR
        albdir(:,nv,nz,2) = albnr(:)*(1.-fsnow(:)) + snonr(:)*fsnow(:)
        albdif(:,nv,nz,2) = albnf(:)*(1.-fsnow(:)) + snonf(:)*fsnow(:)

      end do ! nv
    end do ! nz

      wet_in = max(min(PRMC / POROS,1.0),0.0)  
 
      call catchcn_calc_rc(ntiles,fveg,TCx,QAx,PS,co2v,dayl_fac, &
            T2M10D,TA,cond,psis,rzm,bee,capac,fwet,ZTH,ityp,&
            DRPAR,DFPAR,albdir,albdif,dtc,dea,water_inst,bgc_vegetation_inst,rc00,rcdq,rcdt,&
            laisun,laisha,psnsun,psnsha,lmrsun,lmrsha,parzone,&
            btran)
  
      para(:) = 0. ! zero out absorbed PAR summing array
      do nz = 1,nzone
         do nv = 1,nveg
            do n = 1,ntiles
               if (fveg(n,nv,nz)>1.e-4) then ! account for fact that parzone is undefined if fveg = 0
                 para(n)     = para(n) + parzone(n,nv,nz)*wtzone(n,nz)*fveg(n,nv,nz)
                 if(associated(BTRANT)) then
                    btrant(n) = btrant(n) + btran(n,nv,nz)*fveg(n,nv,nz)*wtzone(n,nz)
                 end if 
               end if
            end do
         end do
      end do

    if(associated(CNCO2)) CNCO2 = CO2V * 1e6
    deallocate (co2v)

    if(associated(PARABS)) parabs = para
    if(associated(PARINC)) parinc = drpar + dfpar

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

    ! Get TPSN1OUT1 for SNOW_ALBEDO parameterization
    
    call STIEGLITZSNOW_CALC_TPSNOW(NTILES, HTSNNN(1,:), WESNN(1,:), TPSN1OUT1, FICE1TMP)
    TPSN1OUT1 =  TPSN1OUT1 + Tzero

    call   StieglitzSnow_snow_albedo(NTILES,N_snow, catchcn_internal%N_CONST_LAND4SNWALB, VEG1, LAI1, ZTH,        &
         RHOFS,                                              &   
         SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                & 
         WESNN, HTSNNN, SNDZN,                               &
         ALBVR, ALBNR, ALBVF, ALBNF,                         & ! instantaneous snow-free albedos on tiles
         SNOVR, SNONR, SNOVF, SNONF,                         & ! instantaneous snow albedos on tiles
         RCONSTIT, UUU, TPSN1OUT1, DRPAR, DFPAR)           
    
    call    SIBALB(NTILES, VEG2,LAI2,GRN, ZTH,         & 
         BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
         ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp, MODIS_SCALE=.TRUE. ) ! instantaneous snow-free albedos on tiles
  
    call   StieglitzSnow_snow_albedo(NTILES,N_snow, catchcn_internal%N_CONST_LAND4SNWALB, VEG2, LAI2, ZTH,        &
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
 
    if (catchcn_internal%SNOW_ALBEDO_INFO == 1) then

        ! use MODIS-derived snow albedo from bcs (via Catch restart)
        !
        ! as a restart parameter from the bcs, snow albedo must not have no-data-values
        ! (checks for unphysical values should be in the make_bcs package)

        call MAPL_GetPointer(INTERNAL,SNOWALB,'SNOWALB',RC=STATUS); VERIFY_(STATUS)

        SNOVR = SNOWALB
        SNONR = SNOWALB
        SNOVF = SNOWALB
        SNONF = SNOWALB 

    endif

    ! --------------------------------------------------------------------------
    ! albedo/swnet partitioning
    ! --------------------------------------------------------------------------
    
    VSUVR = DRPAR + DRUVR
    VSUVF = DFPAR + DFUVR
    
    if(associated(SWDOWNLAND)) SWDOWNLAND = DRPAR + DFPAR + DRUVR + DFUVR + DRNIR + DFNIR
    
    SWNETFREE = (1.-ALBVR)*VSUVR + (1.-ALBVF)*VSUVF + (1.-ALBNR)*DRNIR + (1.-ALBNF)*DFNIR 
    SWNETSNOW = (1.-SNOVR)*VSUVR + (1.-SNOVF)*VSUVF + (1.-SNONR)*DRNIR + (1.-SNONF)*DFNIR 

! set the number of days per year when crossing year boundary or on restart gkw: use GEOS5/MAPL value
! -------------------------------------------------------------------------
    if(AGCM_YY .ne. year_prev) then
      dpy = get_days_per_year(AGCM_YY) ! set the number of days for current year
      year_prev = AGCM_YY
    endif

    ! CN time step over 4 hours may fail; limit to 4 hours; verify that DTCN is a multiple of DT
    ! ------------------------------------------------------------------------------------------
    catchcn_internal%DTCN = min(catchcn_internal%DTCN,14400.)
    if(mod(catchcn_internal%DTCN,dt) /= 0) stop 'dtcn'
    
    ! sum over interval for CN
    ! ------------------------

    tgwm    = tgwm    + tgw
    tpm     = tpm     + tp1
    sfmm    = sfmm    + sfm
    rzmm    = rzmm    + rzm
    bflowm  = bflowm  + bflow
    totwatm = totwatm + totwat
    
    tairm   = tairm   + TA
    rhm     = rhm     + Qair_relative
    windm   = windm   + UU
    rainfm  = rainfm  + (PCU + PLS)
    snowfm  = snowfm  + SNO
    runsrfm = runsrfm + RUNSURF
    ar1m    = ar1m    + car1        
    psnsunm = psnsunm + psnsun
    psnsham = psnsham + psnsha
    lmrsunm = lmrsunm + lmrsun
    lmrsham = lmrsham + lmrsha
    laisunm = laisunm + laisun
    laisham = laisham + laisha
    do n = 1,N_snow
       sndzm(:) = sndzm(:) + sndzn(n,:)
    end do
    asnowm = asnowm + asnow
    cnsum   = cnsum   + 1.
    
    ! call CN model every DTCN seconds
    ! --------------------------------

    if(mod(AGCM_S_ofday,nint(catchcn_internal%DTCN)) == 0) then
       
       cn_count = cn_count + 1

       ! check whether CN is on its first 1.5 hours; since CN_Driver is called once right at the beginning, we set this variable to true when CN_Driver is called for the second time
       if (cn_count .le. 2) then
           first_cn = is_first_step(.true.)
       else
           first_cn = is_first_step(.false.)
       end if

       ! fzeng: pass current date_time to the CN routines.
       call upd_curr_date_time( AGCM_YY, AGCM_MM, AGCM_DD, dofyr, &
            AGCM_HH, AGCM_MI, AGCM_S )
       
       ! compute mean state over interval
       ! --------------------------------
       do nz = 1,nzone
          tgwm(:,nz) = tgwm(:,nz) / cnsum(:)
          rzmm(:,nz) = rzmm(:,nz) / cnsum(:)
          sfmm(:,nz) = sfmm(:,nz) / cnsum(:)
          do nv = 1,nveg
             psnsunm(:,nv,nz) = psnsunm(:,nv,nz) / cnsum(:)
             psnsham(:,nv,nz) = psnsham(:,nv,nz) / cnsum(:)
             lmrsunm(:,nv,nz) = lmrsunm(:,nv,nz) / cnsum(:)
             lmrsham(:,nv,nz) = lmrsham(:,nv,nz) / cnsum(:)
             laisunm(:,nv,nz) = laisunm(:,nv,nz) / cnsum(:)
             laisham(:,nv,nz) = laisham(:,nv,nz) / cnsum(:)
          end do
       end do
       tpm     = tpm     / cnsum
       bflowm  = bflowm  / cnsum
       totwatm = totwatm / cnsum
       tairm   = tairm   / cnsum
       rhm     = rhm     / cnsum
       windm   = windm   / cnsum
       rainfm  = rainfm  / cnsum
       snowfm  = snowfm  / cnsum
       runsrfm = runsrfm / cnsum
       ar1m    = ar1m    / cnsum          
       sndzm   = sndzm   / cnsum
       asnowm  = asnowm  / cnsum

       call CN_Driver(istep_cn,ntiles,ityp,fveg,ndep,tpm,tairm,psis,bee,dayl,btran_fire,ar1m,&
                      rzmm,sfmm,rhm,windm,rainfm,snowfm,TPREC10D,TPREC60D,ET365D,gdp,&
                      abm,peatf,hdm,lnfm,poros,RH30D,totwatm,bflowm,runsrfm,sndzm,&
                      asnowm,TG10D,T2MMIN5D,SNDZM5D,water_inst, first_cn, &
                      psnsunm, psnsham, lmrsunm, lmrsham, laisunm, laisham, wpwet, &
                      elai,esai,tlai,totcolc,npp,gpp,sr,aresp,hresp,nee,burn,closs,nfire,&
                      som_closs,frootc,vegc,xsmr,ndeploy,denit,sminn_leached,sminn,&
                      fire_nloss,leafn,leafc,gross_nmin,net_nmin,&
                      nfix_to_sminn,actual_immob,fpg,fpi,sminn_to_plant,&
                      sminn_to_npool,ndep_to_sminn,totvegn,totlitn,totsomn,&
                      retransn,retransn_to_npool,fuelc,totlitc,cwdc,rootc)

       istep_cn = istep_cn + 1

       ! jkolassa: padd is a correction term that we may no longer need;
       !           I am setting it to zero here in order to avoid having to change
       !           the restart file for now       

       padd(:) = 0.

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
       
       if(associated(CNFIRE_CNT         )) cnfire_cnt          = nfire                   * cnsum ! fire count (s-1)
       if(associated(CNSOM_CLOSS        )) cnsom_closs         = 1.e-3*som_closs         * cnsum ! peat fire C loss (kg/m2/s)
       if(associated(CNNDEPLOY          )) cnndeploy           = 1.e-3*ndeploy           * cnsum           
       if(associated(CNDENIT            )) cndenit             = 1.e-3*denit             * cnsum 
       if(associated(CNSMINN_LEACHED    )) cnsminn_leached     = 1.e-3*sminn_leached     * cnsum 
       if(associated(CNSMINN            )) cnsminn             = 1.e-3*sminn             * cnsum           
       if(associated(CNFIRE_NLOSS       )) cnfire_nloss        = 1.e-3*fire_nloss        * cnsum           
       if(associated(CNLEAFN            )) cnleafn             = 1.e-3*leafn             * cnsum 
       if(associated(CNLEAFC            )) cnleafc             = 1.e-3*leafc             * cnsum 
       if(associated(CNGROSS_NMIN       )) cngross_nmin        = 1.e-3*gross_nmin        * cnsum 
       if(associated(CNNET_NMIN         )) cnnet_nmin          = 1.e-3*net_nmin          * cnsum 
       if(associated(CNNFIX_TO_SMINN    )) cnnfix_to_sminn     = 1.e-3*nfix_to_sminn     * cnsum 
       if(associated(CNACTUAL_IMMOB     )) cnactual_immob      = 1.e-3*actual_immob      * cnsum 
       if(associated(CNFPG              )) cnfpg               = fpg                     * cnsum 
       if(associated(CNFPI              )) cnfpi               = fpi                     * cnsum 
       if(associated(CNSMINN_TO_PLANT   )) cnsminn_to_plant    = 1.e-3*sminn_to_plant    * cnsum 
       if(associated(CNSMINN_TO_NPOOL   )) cnsminn_to_npool    = 1.e-3*sminn_to_npool    * cnsum 
       if(associated(CNNDEP_TO_SMINN    )) cnndep_to_sminn     = 1.e-3*ndep_to_sminn     * cnsum 
       if(associated(CNTOTVEGN          )) cntotvegn           = 1.e-3*totvegn           * cnsum 
       if(associated(CNTOTLITN          )) cntotlitn           = 1.e-3*totlitn           * cnsum 
       if(associated(CNTOTSOMN          )) cntotsomn           = 1.e-3*totsomn           * cnsum 
       if(associated(CNRETRANSN         )) cnretransn          = 1.e-3*retransn          * cnsum 
       if(associated(CNRETRANSN_TO_NPOOL)) cnretransn_to_npool = 1.e-3*retransn_to_npool * cnsum 
       if(associated(CNFUELC            )) cnfuelc             = 1.e-3*fuelc             * cnsum 
       if(associated(CNTOTLITC          )) cntotlitc           = 1.e-3*totlitc           * cnsum 
       if(associated(CNCWDC             )) cncwdc              = 1.e-3*cwdc              * cnsum 
       if(associated(CNROOT             )) cnroot              = 1.e-3*rootc             * cnsum          
       if(associated(CNFSEL             )) cnfsel              = 0. 
       ! reset summing arrays
       ! --------------------
       tgwm    = 0.
       tpm     = 0.
       sfmm    = 0.
       rzmm    = 0.
       bflowm  = 0.
       totwatm = 0.
       tairm   = 0.
       rhm     = 0.
       windm   = 0.
       rainfm  = 0.
       snowfm  = 0.
       runsrfm = 0.
       ar1m    = 0.           
       psnsunm = 0.
       psnsham = 0.
       lmrsunm = 0.
       lmrsham = 0.
       laisunm = 0.
       laisham = 0.
       sndzm   = 0.
       asnowm  = 0.
       cnsum   = 0.
       
    else ! CN diags set to zero
          
       if(associated(CNLAI )) cnlai  = 0.
       if(associated(CNTLAI)) cntlai = 0.
       if(associated(CNSAI )) cnsai  = 0.
       if(associated(CNTOTC)) cntotc = 0.
       if(associated(CNFIRE_CNT         )) cnfire_cnt          = 0. 
       if(associated(CNSOM_CLOSS        )) cnsom_closs         = 0. 
       if(associated(CNNDEPLOY          )) cnndeploy           = 0. 
       if(associated(CNDENIT            )) cndenit             = 0. 
       if(associated(CNSMINN_LEACHED    )) cnsminn_leached     = 0. 
       if(associated(CNSMINN            )) cnsminn             = 0. 
       if(associated(CNFIRE_NLOSS       )) cnfire_nloss        = 0. 
       if(associated(CNLEAFN            )) cnleafn             = 0. 
       if(associated(CNLEAFC            )) cnleafc             = 0. 
       if(associated(CNGROSS_NMIN       )) cngross_nmin        = 0. 
       if(associated(CNNET_NMIN         )) cnnet_nmin          = 0. 
       if(associated(CNNFIX_TO_SMINN    )) cnnfix_to_sminn     = 0. 
       if(associated(CNACTUAL_IMMOB     )) cnactual_immob      = 0. 
       if(associated(CNFPG              )) cnfpg               = 0. 
       if(associated(CNFPI              )) cnfpi               = 0. 
       if(associated(CNSMINN_TO_PLANT   )) cnsminn_to_plant    = 0. 
       if(associated(CNSMINN_TO_NPOOL   )) cnsminn_to_npool    = 0. 
       if(associated(CNNDEP_TO_SMINN    )) cnndep_to_sminn     = 0. 
       if(associated(CNTOTVEGN          )) cntotvegn           = 0. 
       if(associated(CNTOTLITN          )) cntotlitn           = 0. 
       if(associated(CNTOTSOMN          )) cntotsomn           = 0. 
       if(associated(CNRETRANSN         )) cnretransn          = 0. 
       if(associated(CNRETRANSN_TO_NPOOL)) cnretransn_to_npool = 0. 
       if(associated(CNFUELC            )) cnfuelc             = 0. 
       if(associated(CNTOTLITC          )) cntotlitc           = 0. 
       if(associated(CNCWDC             )) cncwdc              = 0. 
       if(associated(CNROOT             )) cnroot              = 0.
       
    endif

       ! CN_Driver outputs at DTCN are saved and used to populate below exports
       !                        uniformly outside DTCN.  
       ! -----------------------------------------------------------------------

       if(associated(CNVEGC)) cnvegc = 1.e-3*vegc  ! * cnsum
       if(associated(CNFROOTC)) cnfrootc = 1.e-3*frootc ! * cnsum
       if(associated(CNNPP )) cnnpp  = 1.e-3*npp   ! * cnsum
       if(associated(CNGPP )) cngpp  = 1.e-3*gpp   ! * cnsum
       if(associated(CNSR  )) cnsr   = 1.e-3*sr    ! * cnsum
       if(associated(CNAR  )) cnar   = 1.e-3*aresp    ! * cnsum
       if(associated(CNHR  )) cnhr   = 1.e-3*hresp    ! * cnsum
       if(associated(CNNEE )) cnnee  = 1.e-3*nee   ! * cnsum
       if(associated(CNXSMR)) cnxsmr = 1.e-3*xsmr  ! * cnsum
       if(associated(CNADD )) cnadd  = 1.e-3*padd  ! * cnsum
       if(associated(CNLOSS)) cnloss = 1.e-3*closs ! * cnsum ! total fire C loss (kg/m2/s)
       if(associated(CNBURN)) cnburn = burn        ! * cnsum ! area fractional fire burn rate (s-1)
       
       ! copy CN_restart vars to catch_internal_rst gkw: only do if stopping
       ! ------------------------------------------
       record = .false. 
       call ESMF_ClockGetAlarm ( CLOCK, alarmname="RecordAlarm001", ALARM=RecordAlarm, RC=STATUS )
       if (status == 0) then
           call ESMF_AlarmGet( RecordAlarm, RingTime = NextRecordTime, _RC)
           if (NextTime == NextRecordTime) record = .true.
       endif

       if(NextTime == StopTime .or. record) then
          
          call CN_exit(ntiles,ityp,fveg,cncol,cnpft)    
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
      nv = 1
      lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
      wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
    end do
    lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

    lai2 = 0.
    wght = 0.
    do nz = 1,nzone
      nv = 2
      lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
      wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
    end do
    lai2 = lai2 / max(wght,1.e-8) ! LAI for secondary vegetation type

    lai = fveg1*lai1 + fveg2*lai2 ! gkw: prognostic LAI on catch_internal_rst (overwrite VEGDYN import)
    LAI0  = max(0.0001     , LAI)

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
        f1(1) = ax1 ; f2(1) = 0. ; f3(1) = 0.
       else
        if(ax1 .lt. (cn1+cn2)) then
          f1(1) = cn1 ; f2(1) = ax1-cn1 ; f3(1) = 0.
         else
          f1(1) = cn1 ; f2(1) = cn2 ; f3(1) = ax1-cn1-cn2
        endif
      endif

      if(ax1 .gt. 0.) then
        rcsat(n) = ax1/(f1(1)/rc00(n,1)+f2(1)/rc00(n,2)+f3(1)/rc00(n,3))
        rcxdt(n) = ax1/(f1(1)/rcdt(n,1)+f2(1)/rcdt(n,2)+f3(1)/rcdt(n,3))
        rcxdq(n) = ax1/(f1(1)/rcdq(n,1)+f2(1)/rcdq(n,2)+f3(1)/rcdq(n,3))
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
          f1(2) = ax2 ; f2(2) = 0. ; f3(2) = 0.
         else
          if(ar .lt. (cn1+cn2)) then
            f1(2) = cn1-ax1 ; f2(2) = ar-cn1 ; f3(2) = 0.
           else
            f1(2) = cn1-ax1 ; f2(2) = cn2 ; f3(2) = ar-cn1-cn2
          endif
        endif
       else
        ar = ax2 + ax4
        if(ar .lt. cn3) then
          f1(2) = 0. ; f2(2) = 0. ; f3(2) = ax2
         else
          if(ax4 .gt. cn3) then
            f1(2) = 0. ; f2(2) = ax2 ; f3(2) = 0.
           else
            f1(2) = 0. ; f2(2) = ar-cn3 ; f3(2) = cn3-ax4
          endif
        endif
      endif

      if(ax2 .gt. 0.) then
        rcuns(n) = ax2/(f1(2)/rc00(n,1)+f2(2)/rc00(n,2)+f3(2)/rc00(n,3))
        rcxdt(n) = ax2/(f1(2)/rcdt(n,1)+f2(2)/rcdt(n,2)+f3(2)/rcdt(n,3))
        rcxdq(n) = ax2/(f1(2)/rcdq(n,1)+f2(2)/rcdq(n,2)+f3(2)/rcdq(n,3))
       else
        rcuns(n) = 1.e3
        rcxdt(n) = 1.e3
        rcxdq(n) = 1.e3
      endif

! compute derivatives
      drcudt(n) = (rcxdt(n) - rcuns(n)) / dtc
      drcudq(n) = (rcxdq(n) - rcuns(n)) / (0.622*dea/PS(n))

    end do

    if(associated(SCSAT )) scsat  = 1. / rcsat
    if(associated(SCUNS )) scuns  = 1. / rcuns

    endif ! end of check for zero tiles

! gkw: end of main CN block

#ifdef DBG_CNLSM_INPUTS
        call MAPL_Get(MAPL, LocStream=LOCSTREAM, RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_LocStreamGet(LOCSTREAM, NT_GLOBAL=NT_GLOBAL, TILEGRID=TILEGRID, RC=STATUS)
        VERIFY_(STATUS)

        call MAPL_TileMaskGet(tilegrid,  mask, rc=status)
        VERIFY_(STATUS)

         if (UNIT_i == 0) then
           unit_i = GETFILE( "catchcnclm51_inputs.data", form="unformatted", RC=STATUS )
           VERIFY_(STATUS)
        endif
        unit = unit_i

! Inputs

        call MAPL_VarWrite(unit, tilegrid, PCU, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, PLS, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SNO, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ICE,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, FRZR, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, UUU, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, EVSBT (:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEDTC (:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT (:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DHSDQC(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, EVSBT (:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEDTC (:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT (:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DHSDQC(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, EVSBT (:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEDTC (:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT (:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DHSDQC(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, EVSBT (:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEVSBT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DEDTC (:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SHSBT (:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DHSDQC(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DSHSBT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, TA, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QA, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FSAT),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FTRN),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FWLT),  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RA(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ZTH,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SWNETFREE,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SWNETSNOW,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, LWDNSRF, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, PS*.01, mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, LAI0,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, GRN0,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, SQSCAT,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RSL1,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RSL2,  mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RDC, mask=mask, rc=status); VERIFY_(STATUS)

        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FSAT) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWN(:,1)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWN(:,1)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FTRN) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWN(:,2)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWN(:,2)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FWLT) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWN(:,3)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWN(:,3)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, QSAT(:,FSNW), mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DQS(:,FSNW) , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, ALWN(:,4)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, BLWN(:,4)   , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RCSAT       , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DRCSDT      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DRCSDQ      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, RCUNS       , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DRCUDT      , mask=mask, rc=status); VERIFY_(STATUS)
        call MAPL_VarWrite(unit, tilegrid, DRCUDQ      , mask=mask, rc=status); VERIFY_(STATUS)

! params
        if (firsttime) then
            firsttime = .false.
           unit = GETFILE( "catchcnclm51_params.data", form="unformatted", RC=STATUS )
           VERIFY_(STATUS)

           call WRITE_PARALLEL(NT_GLOBAL, UNIT)
           call WRITE_PARALLEL(DT, UNIT)
           call WRITE_PARALLEL(catchcn_internal%USE_FWET_FOR_RUNOFF, UNIT)
           call MAPL_VarWrite(unit, tilegrid, LONS,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, LATS,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, VEG1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, VEG2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, FVEG1, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, FVEG2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, DZSF_in_mm,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BF1,   mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BF2,   mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BF3,   mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, VGWMAX,mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CDCR1, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CDCR2, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, PSIS,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BEE,   mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, POROS, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WPWET, mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, COND,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GNU,   mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARS1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARS2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARS3,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA3,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARA4,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW3,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ARW4,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSA1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSA2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSB1,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TSB2,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, ATAU,  mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, BTAU,  mask=mask, rc=status); VERIFY_(STATUS)

           call FREE_FILE(unit, RC=STATUS)
           VERIFY_(STATUS)

! Updates
           unit = GETFILE( "catchcnclm51_updates.data", form="unformatted", RC=STATUS )
           VERIFY_(STATUS)

           call MAPL_VarWrite(unit, tilegrid, TG(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TG(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TG(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TC(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TC(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, TC(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, QC(:,FSAT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, QC(:,FTRN), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, QC(:,FWLT), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CAPAC,      mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, CATDEF,     mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, RZEXC,      mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SRFEXC,     mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(1,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(2,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(3,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(4,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(5,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, GHTCNT(6,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WESNN(1,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WESNN(2,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, WESNN(3,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, HTSNNN(1,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, HTSNNN(2,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, HTSNNN(3,:),mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SNDZN(1,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SNDZN(2,:), mask=mask, rc=status); VERIFY_(STATUS)
           call MAPL_VarWrite(unit, tilegrid, SNDZN(3,:), mask=mask, rc=status); VERIFY_(STATUS)
           
           call FREE_FILE(unit, RC=STATUS)
           VERIFY_(STATUS)

        end if
        DEALLOC_(mask)
#endif

! call unified land model
! -----------------------
        if (ntiles > 0) then

           call CATCHCN ( NTILES, LONS, LATS, DT,catchcn_internal%USE_FWET_FOR_RUNOFF, &                    ! LONS, LATS are in [radians] !!!
                catchcn_internal%FWETC, catchcn_internal%FWETL, cat_id, VEG1,VEG2,FVEG1,FVEG2,DZSF_in_mm,&
                PCU      ,     PLS ,     SNO, ICE, FRZR              ,&
                UUU                                                  ,&

                EVSBT(:,FSAT),     DEVSBT(:,FSAT),     DEDTC(:,FSAT) ,&
                SHSBT(:,FSAT),     DHSDQC(:,FSAT),     DSHSBT(:,FSAT),&
                EVSBT(:,FTRN),     DEVSBT(:,FTRN),     DEDTC(:,FTRN) ,&
                SHSBT(:,FTRN),     DHSDQC(:,FTRN),     DSHSBT(:,FTRN),&
                EVSBT(:,FWLT),     DEVSBT(:,FWLT),     DEDTC(:,FWLT) ,&
                SHSBT(:,FWLT),     DHSDQC(:,FWLT),     DSHSBT(:,FWLT),&
                EVSBT(:,FSNW),     DEVSBT(:,FSNW),     DEDTC(:,FSNW) ,&
                SHSBT(:,FSNW),     DHSDQC(:,FSNW),     DSHSBT(:,FSNW),&

                TA           ,     QA                                ,&

                RA(:,FSAT), RA(:,FTRN), RA(:,FWLT), RA(:,FSNW)       ,&

                ZTH,  SWNETFREE, SWNETSNOW, LWDNSRF                  ,& ! LWDNSRF = *absorbed* longwave only (excl reflected)

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

                EVAPOUT, SHOUT, RUNOFF                               ,&  ! EVAPOUT:                        kg/m2/s
                EVPINT, EVPSOI, EVPVEG, EVPICE                       ,&  ! EVPINT, EVPSOI, EVPVEG, EVPICE: W/m2
                BFLOW                                                ,&
                RUNSURF                                              ,&
                SMELT                                                ,&
                HLWUP                                                ,&  ! *emitted* longwave only (excl reflected)
                SWNDSRF                                              ,&
                LHOUT                                                ,&  ! renamed from HLATN to avoid confusion w/ HLATN=LHFX in SurfaceGC
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
                ENTOT,WTOT, WCHANGE, ECHANGE, HSNACC                 ,&
                EVACC                                                ,&  ! kg/m2/
                LHACC, SHACC                                         ,&  ! W/m2
                TSURF                                                ,&
                SHSNOW1, AVETSNOW1, WAT10CM1, WATSOI1, ICESOI1       ,&
                LHSNOW1, LWUPSNOW1, LWDNSNOW1, NETSWSNOW             ,&
                TCSORIG1, TPSN1IN1, TPSN1OUT1, FSW_CHANGE, FICESOUT  ,&
                TC1_0=TC1_0, TC2_0=TC2_0, TC4_0=TC4_0                ,&
                QA1_0=QA1_0, QA2_0=QA2_0, QA4_0=QA4_0                ,&
                RCONSTIT=RCONSTIT, RMELT=RMELT, TOTDEPOS=TOTDEPOS)

           ! Change units of TP1, TP2, .., TP6 export variables from Celsius to Kelvin.
           ! This used to be done at the level the Surface GridComp.
           ! With this change, gridded TSOIL[n] exports from Surface and tile-space TP[n] exports
           ! from Catch are now consistently in units of Kelvin.
           ! - rreichle, borescan, 6 Nov 2020
           
           TP1 = TP1 + MAPL_TICE
           TP2 = TP2 + MAPL_TICE
           TP3 = TP3 + MAPL_TICE
           TP4 = TP4 + MAPL_TICE
           TP5 = TP5 + MAPL_TICE
           TP6 = TP6 + MAPL_TICE

        end if

        ! compute 365-day running mean of total ET (excluding sublimation from snow)
        if(init_accum_365) then     
           ! 365-day running mean of total ET (W m-2)
           accper = min(istep_365,n365d)
           ET365D = ((accper-1)*ET365D + EVPSOI + EVPINT + EVPVEG) / accper     
        else            
           ET365D = ((n365d-1)*ET365D + EVPSOI + EVPINT + EVPVEG) / n365d
        endif 

        if (OFFLINE_MODE /=0) then
           
           ! in offline mode, disregard the GCM-specific TC/QC modifications and the 
           ! accounting terms for turbulent fluxes (but not snow heat accounting term)
           
           TC(:,FSAT) = TC1_0
           TC(:,FTRN) = TC2_0
           TC(:,FWLT) = TC4_0
           QC(:,FSAT) = QA1_0
           QC(:,FTRN) = QA2_0
           QC(:,FWLT) = QA4_0
           
           EVACC = 0.0
           LHACC = 0.0
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
        if (ntiles > 0) then
        call    SIBALB(NTILES, VEG1,LAI1,GRN, ZTH,         & 
                       BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
                       ALBVR, ALBNR, ALBVF, ALBNF, MODIS_SCALE=.TRUE.  )         ! instantaneous snow-free albedos on tiles

        call STIEGLITZSNOW_CALC_TPSNOW(NTILES, HTSNNN(1,:), WESNN(1,:), TPSN1OUT1, FICE1TMP)
        TPSN1OUT1 =  TPSN1OUT1 + Tzero        

        call    StieglitzSnow_snow_albedo(NTILES,N_snow, catchcn_internal%N_CONST_LAND4SNWALB, VEG1, LAI1, ZTH,        &
                 RHOFS,                                              &   
                 SNWALB_VISMAX, SNWALB_NIRMAX, SLOPE,                & 
                 WESNN, HTSNNN, SNDZN,                               &
                 ALBVR, ALBNR, ALBVF, ALBNF,                         & ! instantaneous snow-free albedos on tiles
                 SNOVR, SNONR, SNOVF, SNONF,                         & ! instantaneous snow albedos on tiles
                 RCONSTIT, UUU, TPSN1OUT1,DRPAR, DFPAR)

        call    SIBALB(NTILES, VEG2,LAI2,GRN, ZTH,         & 
                       BGALBVR, BGALBVF, BGALBNR, BGALBNF, & ! gkw: MODIS soil background albedo
                       ALBVR_tmp, ALBNR_tmp, ALBVF_tmp, ALBNF_tmp, MODIS_SCALE=.TRUE.  ) ! instantaneous snow-free albedos on tiles


        call    StieglitzSnow_snow_albedo(NTILES,N_snow, catchcn_internal%N_CONST_LAND4SNWALB, VEG2, LAI2, ZTH,        &
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

        if (catchcn_internal%SNOW_ALBEDO_INFO == 1) then

           ! use MODIS-derived snow albedo from bcs (via Catch restart)
           !  
           ! as a restart parameter from the bcs, snow albedo must not have no-data-values 
           ! (checks for unphysical values should be in the make_bcs package)

           SNOVR = SNOWALB
           SNONR = SNOWALB
           SNOVF = SNOWALB
           SNONF = SNOWALB

        endif

        ALBVR   = ALBVR    *(1.-ASNOW) + SNOVR    *ASNOW
        ALBVF   = ALBVF    *(1.-ASNOW) + SNOVF    *ASNOW
        ALBNR   = ALBNR    *(1.-ASNOW) + SNONR    *ASNOW
        ALBNF   = ALBNF    *(1.-ASNOW) + SNONF    *ASNOW
        endif        
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
        if(associated(PRLAND)) PRLAND = PCU+PLS+SLDTOT
        if(associated(SNOLAND)) SNOLAND = SLDTOT
        if(associated(DRPARLAND)) DRPARLAND = DRPAR
        if(associated(DFPARLAND)) DFPARLAND = DFPAR

       
        ! -----------------------------------------------------------------------------------------
        !
        ! IMPORTANT: Surface turbulent fluxes in [*]LAND exports are returned as calculated by Catchment! 
        !
        ! For completeness, also return the "accounting" terms that represent the difference between
        ! the flux calculated by Catchment and the flux calculated by the atmosphere (Turbulence GC).
        !
        !   EVAP__calculated_by_atmosphere = EVLAND - EVACC     kg/m2/s
        !   EFLUX_calculated_by_atmosphere = LHLAND - LHACC     W/m2      
        !   HFLUX_calculated_by_atmosphere = SHLAND - SHACC     W/m2
        !
        ! Note: LHACC added for completeness and consistency with the mass flux term (EVACC);
        !       strictly speaking, the atmosphere only receives the evap mass flux and the sensible 
        !       heat flux, but not the latent heat flux.
        !
        ! In previous model versions, the [*]ACC accouting terms were subtracted from the
        ! Catchment-calculated fluxes, thus returning the turbulent fluxes as calculated by 
        ! the atmosphere.
        !
        ! Note: In offline mode, the [*]ACC accounting terms are zeroed out above. 
        !
        ! - reichle, 17 July 2024
         
        if(associated(EVLAND)) EVLAND = EVAPOUT   ! EVLAND is what Catchment thinks it should be 
        if(associated(LHLAND)) LHLAND = LHOUT     ! LHLAND is what Catchment thinks it should be
        if(associated(SHLAND)) SHLAND = SHOUT     ! SHLAND is what Catchment thinks it should be
        
        if(associated(SPWATR)) SPWATR = EVACC
        if(associated(SPLH  )) SPLH   = LHACC
        if(associated(SPLAND)) SPLAND = SHACC
        
        ! Compute latent heat flux that is consistent with the evap mass flux as calculated 
        ! by the atmosphere (Turbulence GC).  In the "flx" HISTORY collection, EFLUX is 
        ! the all-surface latent heat flux, with HLATN (as below) being the land contribution.

        HLATN   = LHOUT   - LHACC

        ! In previous model versions, the evap mass flux EVAPOUT and the sensible heat flux SHOUT
        ! were returned as calculated by Catchment.  These diagnostic are not written in MERRA-2.
        ! In FP and GEOSIT, they are only included in ocean ("ocn") Collections.  They appear to 
        ! be used also in the "gmichem" and "S2S" collections. 
        ! For consistency with previous model versions, keep returning EVAPOUT and SHOUT as 
        ! calculated by Catchment.

        ! Overview of surface turbulent flux variables in subroutine catchment(), the Catch, Land, and Surface Gridded Components, and the M21 "lnd" HISTORY collection
        !
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        ! Description        | Units   | catchment() | Catch[CN]GC                 | LandGC  | SurfaceGC              | HISTORY    | Notes                                                        |
        !                    |         | ArgName     | VarName | export            | export  | export  | averaged     | M21C       |                                                              |
        !                    |         |             |         | (tile)            | (tile)  | (grid)  | over         | "lnd"      |                                                              |
        !=========================================================================================================================================================================================|
        ! evap mass flux     | kg/m2/s | EVAP        | EVAPOUT | EVLAND            | EVLAND  | EVLAND  | land         | EVLAND     |                                                              |
        ! evap mass flux     | kg/m2/s | -           | -       | EVAPOUT           | EVAPOUT | EVAPOUT | all surfaces | n/a        |                                                              |
        ! EV accounting term | kg/m2/s | EVACC       | EVACC   | SPWATR            | SPWATR  | SPWATR  | land         | SPEVLAND   | EVLAND-SPEVLAND = EVAP_from_TurbGC [100% land]               |
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        ! latent heat flux   | W/m2    | LHFLUX      | LHOUT   | LHLAND            | LHLAND  | LHLAND  | land         | LHLAND     |                                                              |
        ! latent heat flux   | W/m2    | -           | -       | HLATN=LHOUT-LHACC | HLATN   | LHFX    | all surfaces | n/a        | consistent w/ EVAP_from_TurbGC [100% land]                   |
        ! LH accounting term | W/m2    | LHACC       | LHACC   | SPLH              | SPLH    | SPLH    | land         | SPLHLAND   | (LHLAND-SPLHLAND) consistent w/ EVAP_from_TurbGC [100% land] |
        ! LH component       | W/m2    | EINT        | EVPINT  | EVPINT            | EVPINT  | EVPINT  | land         | LHLANDINTR |                                                              |
        ! LH component       | W/m2    | ESOI        | EVPSOI  | EVPSOI            | EVPSOI  | EVPSOI  | land         | LHLANDSOIL |                                                              | 
        ! LH component       | W/m2    | EVEG        | EVPVEG  | EVPVEG            | EVPVEG  | EVPVEG  | land         | LHLANDTRNS |                                                              |
        ! LH component       | W/m2    | ESNO        | EVPICE  | EVPICE            | EVPICE  | EVPICE  | land         | LHLANDSBLN |                                                              |
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
        ! sensible heat flux | W/m2    | SHFLUX      | SHOUT   | SHLAND            | SHLAND  | SHLAND  | land         | SHLAND     |                                                              |
        ! sensible heat flux | W/m2    | -           | -       | SHOUT             | SHOUT   | SHOUT   | all surfaces | n/a        |                                                              |
        ! SH accounting term | W/m2    | SHACC       | SHACC   | SPLAND            | SPLAND  | SPLAND  | land         | SPSHLAND   | SHLAND-SPSHLAND = SH_from_TurbGC [for 100% land]             |
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|


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
        if(associated(SPSNOW)) SPSNOW = HSNACC

        if(associated(FRSAT )) FRSAT  = max( min( FR(:,FSAT),1.0 ), 0.0 )
        if(associated(FRUST )) FRUST  = max( min( FR(:,FTRN),1.0 ), 0.0 )
        if(associated(FRWLT )) FRWLT  = max( min( FR(:,FWLT),1.0 ), 0.0 )

        if(associated(SNOMAS)) SNOMAS = WESNN (1,:) + WESNN (2,:) + WESNN (3,:)
        if(associated(SNOWDP)) SNOWDP = SNDZN (1,:) + SNDZN (2,:) + SNDZN (3,:)

        if(associated(FICE1 )) FICE1  = max( min( FICESOUT(1,:),1.0 ), 0.0 )
        if(associated(FICE2 )) FICE2  = max( min( FICESOUT(2,:),1.0 ), 0.0 )
        if(associated(FICE3 )) FICE3  = max( min( FICESOUT(3,:),1.0 ), 0.0 )

        if (N_CONSTIT > 0) then
           if(associated(RMELTDU001)) RMELTDU001 = RMELT(:,1) 
           if(associated(RMELTDU002)) RMELTDU002 = RMELT(:,2) 
           if(associated(RMELTDU003)) RMELTDU003 = RMELT(:,3) 
           if(associated(RMELTDU004)) RMELTDU004 = RMELT(:,4) 
           if(associated(RMELTDU005)) RMELTDU005 = RMELT(:,5) 
           if(associated(RMELTBC001)) RMELTBC001 = RMELT(:,6) 
           if(associated(RMELTBC002)) RMELTBC002 = RMELT(:,7) 
           if(associated(RMELTOC001)) RMELTOC001 = RMELT(:,8) 
           if(associated(RMELTOC002)) RMELTOC002 = RMELT(:,9)
        end if

        if(associated(DZGT1 )) DZGT1  = DZGT(1)                               ! [m]
        if(associated(DZGT2 )) DZGT2  = DZGT(2)                               ! [m]
        if(associated(DZGT3 )) DZGT3  = DZGT(3)                               ! [m]
        if(associated(DZGT4 )) DZGT4  = DZGT(4)                               ! [m]
        if(associated(DZGT5 )) DZGT5  = DZGT(5)                               ! [m]
        if(associated(DZGT6 )) DZGT6  = DZGT(6)                               ! [m]
                            
        if(associated(DZPR  )) DZPR   = CDCR2/(1.-WPWET)/POROS/MAPL_RHOWTR    ! [m] 
        if(associated(DZRZ  )) DZRZ   = VGWMAX/POROS/MAPL_RHOWTR              ! [m]
        if(associated(DZSF  )) DZSF   = DZSF_in_mm/1000.                      ! [m]
        if(associated(DZTS  )) DZTS   = DZTSURF                               ! [m]
                            
        if(associated(WPEMW )) WPEMW  = WPWET*POROS*DZPR*MAPL_RHOWTR          ! [kg/m2]
        if(associated(WPMC  )) WPMC   = WPWET*POROS                           ! [m3/m3]
        
        if(associated(PEATCLSM_FSWCHANGE))  then
           where (POROS >= PEATCLSM_POROS_THRESHOLD)
              PEATCLSM_FSWCHANGE = FSW_CHANGE
           elsewhere
              PEATCLSM_FSWCHANGE = MAPL_UNDEF
           end where
        end if

        if(associated(PEATCLSM_WATERLEVEL)) then
           PEATCLSM_WATERLEVEL = catch_calc_peatclsm_waterlevel( BF1, BF2, CDCR2, POROS, WPWET, CATDEF )
        endif

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

        if (catchcn_internal%N_CONST_LAND4SNWALB /= 0 ) then
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
        if (allocated (ALBVR_tmp)) deallocate ( ALBVR_tmp )
        if (allocated (ALBNR_tmp)) deallocate ( ALBNR_tmp )
        if (allocated (ALBVF_tmp)) deallocate ( ALBVF_tmp )
        if (allocated (ALBNF_tmp)) deallocate ( ALBNF_tmp )
        if (allocated (SNOVR_tmp)) deallocate ( SNOVR_tmp )
        if (allocated (SNONR_tmp)) deallocate ( SNONR_tmp )
        if (allocated (SNOVF_tmp)) deallocate ( SNOVF_tmp )
        if (allocated (SNONF_tmp)) deallocate ( SNONF_tmp )

        deallocate(GHTCNT   )
        deallocate(WESNN    )
        deallocate(HTSNNN   )
        deallocate(SNDZN    )
        deallocate(FICESOUT )
        deallocate(TILEZERO )
        deallocate(DZSF_in_mm)
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
        deallocate(RDC_TMP_1)
        deallocate(RDC_TMP_2)
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
        deallocate(TPSN1IN1 )
        deallocate(TPSN1OUT1)
        deallocate(GHFLXTSKIN)
        deallocate(WCHANGE  )
        deallocate(ECHANGE  )
        deallocate(HSNACC   )
        deallocate(LHOUT    )
        deallocate(EVACC    )
        deallocate(LHACC    )
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
        deallocate(DHSDQC   )
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
        deallocate(FICE1TMP )
        deallocate(SLDTOT   )
        deallocate(FSW_CHANGE)
        deallocate(   btran )
        deallocate(     wgt )
        deallocate(     wpp )
        deallocate(    fwet )
        deallocate(  wet_in )
        deallocate(     sm  )
        deallocate(  SWSRF1 )
        deallocate(  SWSRF2 )
        deallocate(  SWSRF4 )
        deallocate(     tcx )
        deallocate(     qax )
        deallocate(   rcxdt )
        deallocate(   rcxdq )
        deallocate(    car1 )
        deallocate(    car2 )
        deallocate(    car4 )
        deallocate( parzone )
        deallocate(    para )
        deallocate(  totwat )
        deallocate(   nfire )
        deallocate(som_closs)
        deallocate(    dayl )
        deallocate(dayl_fac )
        deallocate(             fsnow )
        deallocate(          ityp_tmp )
        deallocate(     Qair_relative )        
        deallocate(           ndeploy )
        deallocate(             denit )
        deallocate(     sminn_leached )
        deallocate(             sminn )
        deallocate(        fire_nloss )
        deallocate(             leafn )
        deallocate(             leafc )
        deallocate(        gross_nmin )
        deallocate(          net_nmin )
        deallocate(     nfix_to_sminn )
        deallocate(      actual_immob )
        deallocate(               fpg )
        deallocate(               fpi )
        deallocate(    sminn_to_plant )
        deallocate(    sminn_to_npool )
        deallocate(     ndep_to_sminn )
        deallocate(           totvegn )
        deallocate(           totlitn )
        deallocate(           totsomn )
        deallocate(          retransn )
        deallocate( retransn_to_npool )
        deallocate(             fuelc )
        deallocate(           totlitc )
        deallocate(              cwdc )
        deallocate(             rootc )
        deallocate(              lnfm )
        
        deallocate(     tgw )
        deallocate(     rzm )
        deallocate(    rc00 )
        deallocate(    rcdt )
        deallocate(    rcdq )
        deallocate( totcolc )
        deallocate(  wtzone )
        deallocate(               sfm )
        deallocate(            bt )    
        deallocate(     btran_fire )
        deallocate(            albdir )
        deallocate(            albdif )
        deallocate(   elai )
        deallocate(   esai )
        deallocate(   fveg )
        deallocate(   tlai )
        deallocate( psnsun )
        deallocate( psnsha )
        deallocate( laisun )
        deallocate( laisha )
        deallocate(   ityp )
        deallocate(            lmrsun )
        deallocate(            lmrsha )
        deallocate(      ht )
        deallocate(      tp )
        deallocate( soilice )
        deallocate(TA_MIN)
        deallocate(ta_count)
        call MAPL_TimerOff  ( MAPL, "-CATCHCNCLM51" )
        RETURN_(ESMF_SUCCESS)

      end subroutine Driver

! Commented out functions betai(), betacf(), and gammln().
! These functions are not used and were reproduced identically in  
! GEOS_CatchCNCLM40GridComp.F90 and in GEOS_CatchCNCLM45GridComp.F90.
! Another copy was in GEOScatchCN_GridComp/utils/math_routines.F90 but
! there function betai() was missing the restriction 0.0125<x<0.9875.
! - reichle, 23 May 2022
!
!      ! ---------------------------------------------------------
!      
!      FUNCTION betai(a,b,x)
!        REAL betai,a,b,x
!        REAL bt
!        !external gammln
!        
!        if (x < 0.0125) x = 0.0125
!        if (x > 0.9875) x = 0.9875
!        
!        if(x.lt.0..or.x.gt.1.)print *, 'bad argument x in betai',x
!        if(x.lt.0..or.x.gt.1.)stop
!        if(x.eq.0..or.x.eq.1.)then
!           bt=0.
!        else 
!           bt=exp(gammln(a+b)-gammln(a)-gammln(b) &
!                +a*log(x)+b*log(1.-x))
!        endif
!        
!        if(x.lt.(a+1.)/(a+b+2.))then 
!           betai=bt*betacf(a,b,x)/a
!           return
!        else
!           betai=1.-bt*betacf(b,a,1.-x)/b 
!           return
!        endif
!        
!      END FUNCTION betai
!
!      ! -------------------------------------------------------
!
!      FUNCTION betacf(a,b,x)
!
!        INTEGER MAXIT
!        REAL betacf,a,b,x,EPS,FPMIN
!        PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
!        INTEGER m,m2
!        REAL aa,c,d,del,h,qab,qam,qap
!        
!        qab=a+b 
!        qap=a+1. 
!        qam=a-1.
!        c=1. 
!        d=1.-qab*x/qap
!        
!        if(abs(d).lt.FPMIN)d=FPMIN
!        d=1./d
!        h=d
!        do m=1,MAXIT
!           m2=2*m
!           aa=m*(b-m)*x/((qam+m2)*(a+m2))
!           d=1.+aa*d 
!           if(abs(d).lt.FPMIN)d=FPMIN
!           c=1.+aa/c
!           if(abs(c).lt.FPMIN)c=FPMIN
!           d=1./d
!           h=h*d*c
!           aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
!           d=1.+aa*d 
!           if(abs(d).lt.FPMIN)d=FPMIN
!           c=1.+aa/c
!           if(abs(c).lt.FPMIN)c=FPMIN
!           d=1./d
!           del=d*c
!           h=h*del
!           if(abs(del-1.).lt.EPS)exit 
!        enddo
!        betacf=h
!        return
!        
!      END FUNCTION betacf
!      
!      ! --------------------------------------------------------------
!      
!      FUNCTION gammln(xx)
!        
!        REAL gammln,xx
!        INTEGER j
!        DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
!        
!        SAVE cof,stp
!        DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,          &
!             24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
!             -.5395239384953d-5,2.5066282746310005d0/
!        x=xx
!        y=x
!        tmp=x+5.5d0
!        tmp=(x+0.5d0)*log(tmp)-tmp
!        ser=1.000000000190015d0
!        do  j=1,6
!           y=y+1.d0
!           ser=ser+cof(j)/y
!        enddo
!        gammln=tmp+log(stp*ser/x)
!        return
!        
!      END FUNCTION gammln

      ! --------------------------------------------------------------
      
      integer function VarID (NCFID, VNAME) 
        
        integer, intent (in)      :: NCFID
        character(*), intent (in) :: VNAME
        integer                   :: status
        
        STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID); VERIFY_(STATUS)
        
      end function VarID

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
  !! -asnow-emis-ww-fr-D[xx]
  real, pointer :: ity(:,:)=>null()
  real, pointer :: fvg(:,:)=>null()
  real, pointer :: asnow(:)=>null()
  real, pointer :: emis(:)=>null()
  real, pointer :: ww(:,:)=>null()
  real, pointer :: fr(:,:)=>null()
  real, pointer :: delCQ_delTVA(:,:)=>null()
  real, pointer :: delCH_delTVA(:,:)=>null()
  real, pointer :: delCH_delTC(:,:)=>null()
  real, pointer :: delCQ_delQC(:,:)=>null()

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
  real, pointer :: bf1(:)=>null()
  real, pointer :: bf2(:)=>null()

  !! Miscellaneous
  integer :: ntiles, nv, nz
  real, allocatable :: dummy(:)
  real, allocatable :: DZSF_in_mm(:), ar1(:), ar2(:), wesnn(:,:)
  real, allocatable :: catdefcp(:), srfexccp(:), rzexccp(:)
  real, allocatable :: VEG1(:), VEG2(:)
  integer, allocatable :: ityp(:,:,:)
  real,    allocatable :: fveg(:,:,:), elai(:,:,:), esai(:,:,:), wtzone(:,:), lai1(:), lai2(:), wght(:)
  real, allocatable,dimension(:) :: fveg1, fveg2

  type(T_CATCHCN_STATE), pointer :: catchcn_internal
  type(CATCHCN_WRAP) :: wrap

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

  call ESMF_UserCompGetInternalState(gc, 'CatchcnInternal', wrap, status)
  VERIFY_(status)
  catchcn_internal => wrap%ptr
  ! Pointers to IMPORTs
  call MAPL_GetPointer(import, ps, 'PS', rc=status)
  VERIFY_(status)

  ! Pointers to EXPORTs
  call MAPL_GetPointer(export, asnow, 'ASNOW', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(export, emis, 'EMIS', rc=status)
  VERIFY_(status)

  ! Pointers to INTERNALs
  call MAPL_GetPointer(INTERNAL, ITY, 'ITY', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, FVG, 'FVG', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, fr, 'FR', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, ww, 'WW', rc=status)
  VERIFY_(status)

  if     (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND == 1) then
 
     call MAPL_GetPointer(INTERNAL, delCQ_delTVA, 'delCQ_delTVA', rc=status)
     VERIFY_(status)
     call MAPL_GetPointer(INTERNAL, delCH_delTVA, 'delCH_delTVA', rc=status)
     VERIFY_(status)

  elseif (CATCHCN_INTERNAL%MOSFC_EXTRA_DERIVS_OFFL_LAND >= 2) then 
     
     call MAPL_GetPointer(INTERNAL, delCH_delTC, 'delCH_delTC', rc=status)
     VERIFY_(status)
     call MAPL_GetPointer(INTERNAL, delCQ_delQC, 'delCQ_delQC', rc=status)
     VERIFY_(status)
     
  end if

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
  call MAPL_GetPointer(INTERNAL, bf1, 'BF1', rc=status)
  VERIFY_(status)
  call MAPL_GetPointer(INTERNAL, bf2, 'BF2', rc=status)
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

   call get_CN_LAI(ntiles,ityp,fveg,elai,esai=esai)

   lai1 = 0.
   wght = 0.
   do nz = 1,num_zon
     nv = 1
     lai1(:) = lai1(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
     wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
   end do
   lai1 = lai1 / max(wght,1.e-8) ! LAI for primary vegetation type

   lai2 = 0.
   wght = 0.
   do nz = 1,num_zon
     nv = 2
     lai2(:) = lai2(:) + max(elai(:,nv,nz),0.)*fveg(:,nv,nz)*wtzone(:,nz)
     wght(:) = wght(:) +                       fveg(:,nv,nz)*wtzone(:,nz)
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
   endwhere
   where(ITY(:,2) > 0.)
     VEG2 = map_cat(nint(ITY(:,2)))  ! gkw: secondary CN PFT type mapped to catchment type; ITY should be > 0 even if FVEG=0
   endwhere
   _ASSERT((count(VEG1>NTYPS.or.VEG1<1)==0),'needs informative message') 
   _ASSERT((count(VEG2>NTYPS.or.VEG2<1)==0),'needs informative message')
   fveg1(:) = fvg(:,1) 
   fveg2(:) = fvg(:,2) 

  ! Compute ASNOW and EMIS
  allocate(wesnn(3,ntiles), stat=status)
  VERIFY_(status)
  wesnn(1,:) = wesnn1
  wesnn(2,:) = wesnn2
  wesnn(3,:) = wesnn3
  call StieglitzSnow_calc_asnow(N_snow, ntiles, wesnn, asnow)

  EMIS    = fveg1*(EMSVEG(NINT(VEG1)) + (EMSBARESOIL - EMSVEG(NINT(VEG1)))*exp(-LAI1)) + &
       fveg2*(EMSVEG(NINT(VEG2)) + (EMSBARESOIL - EMSVEG(NINT(VEG2)))*exp(-LAI2))

  emis = emis*(1.-asnow) + EMSSNO*asnow

  ! Compute FR
  ! Step 1: set dzsf
  ! Step 2: compute ar1, ar2 via call to catch_calc_soil_moist()
  ! Step 3: compute fr

  ! -step-1-
  allocate(DZSF_in_mm(ntiles), stat=status)
  VERIFY_(status)
  DZSF_in_mm = catchcn_internal%SURFLAY

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
       ntiles, DZSF_in_mm, vgwmax, cdcr1, cdcr2,                                &
       psis, bee, poros, wpwet,                                                 &
       ars1, ars2, ars3,                                                        &
       ara1, ara2, ara3, ara4,                                                  &
       arw1, arw2, arw3, arw4, bf1, bf2,                                        &
       srfexccp, rzexccp, catdefcp,                                             &
       ar1, ar2, dummy                                                          &
       )

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
  if (allocated(catdefcp))   deallocate(catdefcp)
  if (allocated(srfexccp))   deallocate(srfexccp)
  if (allocated(rzexccp))    deallocate(rzexccp)
  if (allocated(dummy))      deallocate(dummy)
  if (allocated(DZSF_in_mm)) deallocate(DZSF_in_mm)
  if (allocated(ar1))        deallocate(ar1)
  if (allocated(ar2))        deallocate(ar2)
  if (allocated(wesnn))      deallocate(wesnn)
  if (allocated(fveg1))      deallocate(fveg1)
  if (allocated(fveg2))      deallocate(fveg2)
  if (allocated(veg1))       deallocate(veg1)
  if (allocated(veg2))       deallocate(veg2)
  if (allocated(lai1))       deallocate(lai1)
  if (allocated(lai2))       deallocate(lai2)
  if (allocated(wght))       deallocate(wght)

  ! All done
  RETURN_(ESMF_SUCCESS)

end subroutine RUN0

end module GEOS_CatchCNCLM51GridCompMod

subroutine SetServices(gc, rc)
   use ESMF
   use GEOS_CatchCNCLM51GridCompMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine

