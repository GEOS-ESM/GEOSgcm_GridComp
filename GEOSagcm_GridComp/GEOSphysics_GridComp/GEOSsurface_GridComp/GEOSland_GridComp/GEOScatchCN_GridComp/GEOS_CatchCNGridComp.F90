! It is a proxy of clm 4.0 and clm 4.5

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_CatchCNGridCompMod

  use ESMF
  use MAPL
  USE STIEGLITZSNOW,   ONLY :                 &
       N_CONSTIT,   &
       NUM_DUDP, NUM_DUSV, NUM_DUWT, NUM_DUSD, &
       NUM_BCDP, NUM_BCSV, NUM_BCWT, NUM_BCSD, &
       NUM_OCDP, NUM_OCSV, NUM_OCWT, NUM_OCSD, &
       NUM_SUDP, NUM_SUSV, NUM_SUWT, NUM_SUSD, &
       NUM_SSDP, NUM_SSSV, NUM_SSWT, NUM_SSSD

  use  catch_wrap_stateMod

  implicit none
  private
! !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices

  integer :: CATCHCN

contains
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
    character(len=ESMF_MAXSTR)   :: CATCHCN_VERSION
    type(ESMF_GridComp), pointer :: gcs(:)
    type(T_CATCHCN_STATE), pointer :: CATCHCN_INTERNAL_STATE
    class(T_CATCH_STATE),  pointer :: statePtr
    type(CATCHCN_WRAP)             :: wrap

    character(len=ESMF_MAXSTR)              :: SURFRC
    type(ESMF_Config)                       :: SCF, CF
    integer                                 :: LSM_CHOICE
    character(len=ESMF_MAXSTR)              :: tmp
    integer                                 :: NUM_LDAS_ENSEMBLE, ens_id_width

! Begin...
! --------
    Iam='SetServices'
    call ESMF_GridCompGet ( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//trim(Iam)


    call MAPL_GetObjectFromGC(gc, MAPL, rc=status)
    VERIFY_(status)

    call MAPL_GetResource ( MAPL, NUM_LDAS_ENSEMBLE, Label="NUM_LDAS_ENSEMBLE:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, ens_id_width, Label="ENS_ID_WIDTH:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    allocate(CATCHCN_INTERNAL_STATE)
    statePtr =>CATCHCN_INTERNAL_STATE

    ! resource variables for offline GEOSldas; for documentation, see GEOSldas/src/Applications/LDAS_App/GEOSldas_LDAS.rc
    call MAPL_GetResource ( MAPL, CATCHCN_INTERNAL_STATE%CATCH_OFFLINE, Label="CATCHMENT_OFFLINE:", DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, CATCHCN_INTERNAL_STATE%CATCH_SPINUP,  Label="CATCHMENT_SPINUP:",  DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    ! resource variables from GEOS_SurfaceGridComp.rc
    call MAPL_GetResource ( MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)

    call surface_params_to_wrap_state(statePtr, SCF,  _RC)

    call ESMF_ConfigDestroy(SCF, _RC)

    call MAPL_Get (MAPL, CF=CF, _RC)
    call ESMF_ConfigSetAttribute(CF, value=CATCHCN_INTERNAL_STATE%ATM_CO2, Label='ATM_CO2:', _RC)
    call ESMF_ConfigSetAttribute(CF, value=CATCHCN_INTERNAL_STATE%N_CONST_LAND4SNWALB, Label='N_CONST_LAND4SNWALB:', _RC)
    call ESMF_ConfigSetAttribute(CF, value=CATCHCN_INTERNAL_STATE%RUN_IRRIG, Label='RUN_IRRIG:', _RC)
    call ESMF_ConfigSetAttribute(CF, value=CATCHCN_INTERNAL_STATE%PRESCRIBE_DVG, Label='PRESCRIBE_DVG:', _RC)
    call ESMF_ConfigSetAttribute(CF, value=CATCHCN_INTERNAL_STATE%SNOW_ALBEDO_INFO, Label='SNOW_ALBEDO_INFO:', _RC)
    call MAPL_Set (MAPL, CF=CF, _RC)

    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=2, RC=STATUS)
    VERIFY_(STATUS)
    tmp = ''
    if (NUM_LDAS_ENSEMBLE >1) then
        !catchcn_exxxx
        tmp(1:ens_id_width)=COMP_NAME(8:8+ens_id_width-1)
    endif
    if ( LSM_CHOICE == 2 ) then
       CATCHCN = MAPL_AddChild('CATCHCNCLM40'//trim(tmp), 'setservices_', parentGC=GC, sharedObj='libGEOScatchCNCLM40_GridComp.so', RC=STATUS)
       VERIFY_(STATUS)
    else if ( LSM_CHOICE == 4 ) then
       CATCHCN = MAPL_AddChild('CATCHCNCLM51'//trim(tmp), 'setservices_', parentGC=GC, sharedObj='libGEOScatchCNCLM51_GridComp.so', RC=STATUS)
       VERIFY_(STATUS)
    else
       _ASSERT( .false., " LSM_CHOICE should equal 2 (CLM40) or 4 (CLM51)")
    endif

    wrap%ptr =>CATCHCN_INTERNAL_STATE
    call ESMF_UserCompSetInternalState(gc, 'CatchcnInternal', wrap, status)
    VERIFY_(status)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN1, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN, RUN2, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE, Finalize, RC=status)
    VERIFY_(status)

! Set the state variable specs. ( should be the combinations of clm4.0 and clm4.5
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

    IF (CATCHCN_INTERNAL_STATE%ATM_CO2 == 4) THEN
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
         LONG_NAME          = 'vegetation_greenness_fraction',&
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

  !  EXPORT STATE:
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LST',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TST',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QST',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DELTS',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DELQS',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBVR',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBVF',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBNR',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ALBNF',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EMIS',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNOWMASS', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNOWDP',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHFLX',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPUNST',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSURF',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSNOW',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPWLT',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSAT',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ASNOW',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SHSNOW',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'AVETSNOW', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FRSAT',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FRUST',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'FRWLT',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WET1',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WET2',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WET3',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WCSF',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WCRZ',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WCPR',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP1',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP2',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP3',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP4',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP5',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TP6',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVAPOUT',  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SUBLIM' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SHOUT'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RUNOFF' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPINT' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPSOI' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPVEG' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPICE' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WAT10CM',  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WATSOI' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ICESOI' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVPSNO' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BASEFLOW', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RUNSURF',  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SMELT'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HLWUP'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWNDSRF',  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWNDSRF',  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HLATN'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QINFIL' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ACCUM'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'EVLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PRLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNOLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DRPARLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DFPARLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LHSNOW' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWNETSNOW' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWUPSNOW' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWDNSNOW' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TCSORIG' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSN1IN' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TPSN1OUT' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LHLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SHLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SWDOWNLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'LWLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHSNOW' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTSKIN' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SMLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TWLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TELAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TSLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DWLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'DHLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPLAND' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPWATR' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SPSNOW' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WESNN1' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WESNN2' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'WESNN3' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CAPAC'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC, SHORT_NAME = 'POROS', CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'COND' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'PSIS' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BEE'  , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'WPWET', CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'GNU'  , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'VGWMAX',CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF1'  , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF2'  , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BF3'  , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'CDCR1', CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'CDCR2', CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS1' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS2' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARS3' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA1' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA2' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA3' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARA4' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW1' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW2' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW3' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ARW4' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSA1' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSA2' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSB1' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'TSB2' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'ATAU' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC, SHORT_NAME = 'BTAU' , CHILD_ID = CATCHCN, RC=STATUS); VERIFY_(STATUS)

!   From catmentcn grid internal to be perturbed by land_pert grid
!   WESNN1-3 are originally exported
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TC'     , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TG'     , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QC'     , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CATDEF' , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RZEXC'  , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SRFEXC' , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HTSNNN1', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HTSNNN2', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'HTSNNN3', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNDZN1' , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNDZN2' , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SNDZN3' , CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT1', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT2', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT3', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT4', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT5', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GHTCNT6', CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)

    ! Unified CN from RUN2 of the first catchment instance

    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNLAI'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNTLAI' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNSAI'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNTOTC' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNVEGC' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNROOT' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    if (LSM_CHOICE >= 3) then ! jkolassa: needed for CNCLM51
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNFROOTC' ,  CHILD_ID = CATCHCN, RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNAR'   ,  CHILD_ID = CATCHCN, RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNHR'   ,  CHILD_ID = CATCHCN, RC=STATUS  )
       VERIFY_(STATUS)
    endif
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNNPP'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNGPP'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNSR'   ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNNEE'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNXSMR' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNADD'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNLOSS' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNBURN' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PARABS' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PARINC' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SCSAT'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SCUNS'  ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'BTRANT' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'SIF'    ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNFSEL' ,  CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'TH',       CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'QH',       CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CHT',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CQT',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CMT',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'CNT',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RIT',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'Z0',       CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'D0',       CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'Z0H',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'VENT',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'GUST',     CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOU50M',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOV50M',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOT10M',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOQ10M',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOU10M',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOV10M',   CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOT2M',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOQ2M',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOU2M',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'MOV2M',    CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ITY',      CHILD_ID = CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU001', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU002', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU003', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU004', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTDU005', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTBC001', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTBC002', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTOC001', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RMELTOC002', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PEATCLSM_WATERLEVEL',CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC, SHORT_NAME = 'PEATCLSM_FSWCHANGE' ,CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)

    if (CATCHCN_INTERNAL_STATE%N_CONST_LAND4SNWALB /= 0) then
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU001', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU002', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU003', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU004', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RDU005', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RBC001', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'RBC002', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ROC001', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ROC002', CHILD_ID = CATCHCN, RC=STATUS) ; VERIFY_(STATUS)
    endif

!EOS

    call MAPL_TimerAdd(GC,    name="INITIALIZE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN1"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN2"  ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final method
! ---------------------------------

    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

end subroutine SetServices



! !IROUTINE: Initialize -- Initialize method for the composite Surface Gridded Component

! !INTERFACE:

subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)           :: IAm
    integer                              :: STATUS
    character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp    ), pointer       :: MAPL
    type (MAPL_MetaComp    ), pointer       :: CHILD_MAPL
    type (MAPL_LocStream       )            :: LOCSTREAM
    type (ESMF_GridComp        ), pointer   :: GCS(:)
    type (CATCHCN_WRAP)                     :: wrap
    integer                                 :: I

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOn(MAPL,"INITIALIZE", RC=STATUS ); VERIFY_(STATUS)

    call MAPL_Get (MAPL, LOCSTREAM=LOCSTREAM, GCS=GCS, RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_UserCompGetInternalState(gc, 'CatchcnInternal', wrap, status)
    VERIFY_(status)

! Place the land tilegrid in the generic state of each child component
!---------------------------------------------------------------------
    do I = 1, SIZE(GCS)
       call MAPL_GetObjectFromGC( GCS(I), CHILD_MAPL, RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_Set (CHILD_MAPL, LOCSTREAM=LOCSTREAM, RC=STATUS )
       VERIFY_(STATUS)
       call ESMF_UserCompSetInternalState(gcs(I), 'CatchcnInternal', wrap, status)
       VERIFY_(status)
    end do

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)

    call MAPL_TimerOff(MAPL,"INITIALIZE", RC=STATUS ); VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"TOTAL", RC=STATUS ); VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

end subroutine Initialize

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
    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    character(len=ESMF_MAXSTR),pointer  :: GCnames(:)
    integer                             :: I

! ------------------------------------------------------------------------------
! Get the target component's name and set-up traceback handle.
! ------------------------------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//"::RUN1"

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN1")

    call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCnames=GCnames,rc=STATUS)
    VERIFY_(STATUS)

! Call the children's RUN methods
!--------------------------------

    DO I = 1, size(GCS)
       call MAPL_TimerOn(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
       call ESMF_GridCompRun(GCS(I), importState=GIM(I), exportState=GEX(I), &
                             CLOCK=CLOCK, PHASE=1, userRC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
    END DO

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

    character(len=ESMF_MAXSTR) :: Iam
    integer :: STATUS
    character(len=ESMF_MAXSTR) :: COMP_NAME

! ------------------------------------------------------------------------------
! Local derived type aliases
! ------------------------------------------------------------------------------

    type(MAPL_MetaComp),pointer     :: MAPL


    type (ESMF_GridComp),      pointer  :: GCS(:)
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    character(len=ESMF_MAXSTR),pointer  :: GCnames(:)
    integer                             :: I

! ------------------------------------------------------------------------------
! Get the target component's name and set-up traceback handle.
! ------------------------------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam=trim(COMP_NAME)//"::RUN2"

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

! Start timers
! ------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN2")

    call MAPL_Get (MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCnames=GCnames,rc=STATUS)
    VERIFY_(STATUS)

! Call the children's RUN methods
!--------------------------------

    DO I = 1, size(GCS)
       call MAPL_TimerOn(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
       call ESMF_GridCompRun(GCS(I), importState=GIM(I), exportState=GEX(I), &
                             CLOCK=CLOCK, PHASE=2, userRC=STATUS)
       VERIFY_(STATUS)
       call MAPL_TimerOff(MAPL,trim(GCnames(i)), RC=STATUS ); VERIFY_(STATUS)
    END DO

!  All done
! ------------------------------------------------------------------------------

    call MAPL_TimerOff ( MAPL, "RUN2"  )
    call MAPL_TimerOff ( MAPL, "TOTAL" )
    RETURN_(ESMF_SUCCESS)
end subroutine RUN2

subroutine Finalize(gc, import, export, clock, rc)
    type(ESMF_GridComp), intent(inout) :: gc     ! Gridded component
    type(ESMF_State),    intent(inout) :: import ! Import state
    type(ESMF_State),    intent(inout) :: export ! Export state
    type(ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,   intent(  out) :: rc     ! Error code

    ! !DESCRIPTION:
    ! Clean-up.
    !EOP
    ! ErrLog variables
    integer :: status
    character(len=ESMF_MAXSTR) :: Iam
    character(len=ESMF_MAXSTR) :: comp_name
    ! Begin...
    ! Get component's name and setup traceback handle
    call ESMF_GridCompget(gc, name=comp_name, rc=status)
    VERIFY_(status)
    Iam = trim(comp_name) // "::Finalize"
    ! Call Finalize for every child
    call MAPL_GenericFinalize(gc, import, export, clock, rc=status)
    VERIFY_(status)
    ! End
    RETURN_(ESMF_SUCCESS)

end subroutine Finalize

end module GEOS_CatchCNGridCompMod

