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
    character(len=ESMF_MAXSTR)              :: SURFRC
    type(ESMF_Config)                       :: SCF
    integer :: DO_GOSWIM, LSM_CHOICE, ATM_CO2
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

    tmp = ''
    if (NUM_LDAS_ENSEMBLE >1) then
        !catchcn_exxxx
        !ENS_ID_WIDTH is the width of the digit. Add 2 for '_e'
        tmp(1:ens_id_width+2)=COMP_NAME(8:8+ens_id_width-1 +2)
    endif


    call MAPL_GetResource ( MAPL, LSM_CHOICE, Label="LSM_CHOICE:", DEFAULT=2, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, SURFRC, label = 'SURFRC:', default = 'GEOS_SurfaceGridComp.rc', RC=STATUS) ; VERIFY_(STATUS)
    SCF = ESMF_ConfigCreate(rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(SCF,SURFRC,rc=status) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='ATM_CO2:', value=ATM_CO2, DEFAULT=2, RC=STATUS) ; VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute (SCF, label='N_CONST_LAND4SNWALB:'  , value=DO_GOSWIM  , DEFAULT=0, RC=STATUS); VERIFY_(STATUS)

    if ( LSM_CHOICE == 2 ) then
       CATCHCN = MAPL_AddChild('CATCHCNCLM40'//trim(tmp), 'setservices_', parentGC=GC, sharedObj='libGEOScatchCNCLM40_GridComp.so', RC=STATUS)
       VERIFY_(STATUS)       
    else if ( LSM_CHOICE == 3 ) then
       CATCHCN = MAPL_AddChild('CATCHCNCLM45'//trim(tmp), 'setservices_', parentGC=GC, sharedObj='libGEOScatchCNCLM45_GridComp.so', RC=STATUS)
       VERIFY_(STATUS)       
    else
       _ASSERT( .false., " LSM_CHOICE should equal 2 (CLM40) or 3 (CLM45)")
    endif


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

    IF (ATM_CO2 == 4) THEN
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

  !  EXPORT STATE:
    call MAPL_AddExportSpec ( GC, CATCHCN, RC=STATUS  )
    VERIFY_(STATUS)

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

! Place the land tilegrid in the generic state of each child component
!---------------------------------------------------------------------
    do I = 1, SIZE(GCS)
       call MAPL_GetObjectFromGC( GCS(I), CHILD_MAPL, RC=STATUS )
       VERIFY_(STATUS)
       call MAPL_Set (CHILD_MAPL, LOCSTREAM=LOCSTREAM, RC=STATUS )
       VERIFY_(STATUS)
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

