    ! !INTERNAL STATE:

    !BOS

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS',    &
         default    = 1.0e-6,                                      &
         RESTART    = MAPL_RestartRequired,                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                        &
         SHORT_NAME = 'QLLS',                                            &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                         &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                             &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLLS',                                      &
         LONG_NAME  = 'large_scale_cloud_area_fraction',           &
         UNITS      = '1',                                         &
         FRIENDLYTO = 'DYNAMICS',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLCN',                                      &
         LONG_NAME  = 'convective_cloud_area_fraction',            &
         UNITS      = '1',                                         &
         FRIENDLYTO = 'DYNAMICS',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddInternalSpec(GC,                                     &
         SHORT_NAME = 'QILS',                                         &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                      &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME = 'QICN',                                        &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                     &
         FRIENDLYTO = 'DYNAMICS:TURBULENCE',                         &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

 

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPL',                                       &
         LONG_NAME  ='particle_number_for_liquid_cloud',           &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NCPL),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,          &
         DEFAULT = 50.0e6 ,   RC=STATUS  )  
     
         VERIFY_(STATUS)                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPI',                                       &
         LONG_NAME  ='particle_number_for_ice_cloud',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NCPI),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,        &
         DEFAULT = 1.0e3,   RC=STATUS  )  
    
        VERIFY_(STATUS)                                                   


    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NRAIN',                                       &
         LONG_NAME  ='particle_number_for_rain',           &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NRAIN),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,          &
         DEFAULT = 0.0 ,   RC=STATUS  )  
     
         VERIFY_(STATUS)                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NSNOW',                                       &
         LONG_NAME  ='particle_number_for_snow',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NSNOW),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,        &
         DEFAULT = 0.0,   RC=STATUS  )  
    
        VERIFY_(STATUS)                                                   

       call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NGRAUPEL',                                       &
         LONG_NAME  ='particle_number_for_graupel',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES_NGRAUPEL),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,        &
         DEFAULT = 0.0,   RC=STATUS  )  
    
        VERIFY_(STATUS)


    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QRAIN',                                     &
         LONG_NAME  = 'mass_fraction_of_rain',                     & 
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES_QRAIN),                       &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QSNOW',                                     &
         LONG_NAME  = 'mass_fraction_of_snow',                     &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES_QSNOW),                       &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QGRAUPEL',                                  &
         LONG_NAME  = 'mass_fraction_of_graupel',                  &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES_QGRAUPEL),                       &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QW',                                        &
         LONG_NAME  = 'mass_fraction_of_wet_air',                  &
         UNITS      = 'kg kg-1',                                   &
         RESTART    = MAPL_RestartSkip,                            &
         FRIENDLYTO = 'TURBULENCE:' // trim(COMP_NAME),            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


     
    
  !if (DOSHLW /= 0) then
  !        
  !  call MAPL_AddInternalSpec(GC,                                  &
  !       SHORT_NAME ='CUSH',                                       &
  !        LONG_NAME = 'Cumulus_scale_height_from_UW_shlw_convection',&
  !       UNITS      ='m',                                           &       
  !       DIMS      = MAPL_DimsHorzOnly,                             &
  !       VLOCATION = MAPL_VLocationNone,                            &
  !       DEFAULT= 1000.0,                                           &
  !       RC=STATUS  )  
  !       VERIFY_(STATUS)                         
  !  
  !end if 
  
    ! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'geopotential_height',                       &
         UNITS      = 'm',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
                                                        RC=STATUS  )

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PLE',                                         &
         LONG_NAME  = 'air_pressure',                                &
         UNITS      = 'Pa',                                          &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)                                                                          

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PKE',                                         &
         LONG_NAME  = 'edge_p$^\kappa$',                      &
         UNITS      = 'Pa$^\kappa$',                               &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PLK',                                       &
         LONG_NAME  = 'mid-layer_p$^\kappa$',                      &
         UNITS      = 'Pa$^\kappa$',                               &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationCenter,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'PREF',                                       &
         LONG_NAME  = 'reference_air_pressure',                     &
         UNITS      = 'Pa',                                         &
         DIMS       = MAPL_DimsVertOnly,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'KH',                                         &
         LONG_NAME  = 'scalar_diffusivity',                         &
         UNITS      = 'm+2 s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                              &
         SHORT_NAME = 'TKE',                                        &
         LONG_NAME  = 'turbulent_kinetic_energy',                   &
         UNITS      = 'm+2 s-2',                                    &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TH',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'U',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'V',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'W',                                         &
         LONG_NAME  = 'vertical_velocity',                         &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TS',                                        &
         LONG_NAME  = 'surface temperature',                       &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'SNOMAS',                                    &
         LONG_NAME  = 'snow_mass',                       &
         UNITS      = 'kg/m2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'FRLANDICE',                                    &
         LONG_NAME  = 'areal_landice_fraction',                       &
         UNITS      = '1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'FRLAND',                                    &
         LONG_NAME  = 'areal_land_fraction',                       &
         UNITS      = '1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'FROCEAN',                                   &
         LONG_NAME  = 'areal_ocean_fraction',                      &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'FRACI',                                     &
         LONG_NAME  = 'ice_covered_fraction_of_tile',              &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         AVERAGING_INTERVAL = AVRGNINT,                            &
         REFRESH_INTERVAL   = RFRSHINT,                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                   &
       SHORT_NAME = 'ZPBL',                                       &
       LONG_NAME  = 'planetary_boundary_layer_height',            &
       UNITS      = 'm',                                          &
       DIMS       = MAPL_DimsHorzOnly,                            &
       VLOCATION  = MAPL_VLocationNone,                           &
       RC=STATUS  )
    VERIFY_(STATUS)

    ! These bundles should be changed when we merge w/ the head.

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'MTR',                                        &
         LONG_NAME  = 'tracers_for_moist',                          &
         UNITS      = 'X',                                          &
         DATATYPE   = MAPL_BundleItem,                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RESTART    = MAPL_RestartSkip,                            &
         RC=STATUS )
    VERIFY_(STATUS)                                                                           

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME = 'KPBL',                                       &
        LONG_NAME  = 'planetary_boundary_layer_level',             &
        UNITS      = '1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                            &
        VLOCATION  = MAPL_VLocationNone,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
                                                        RC=STATUS  )
     VERIFY_(STATUS)      

     call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME = 'KPBL_SC',                                    &
        LONG_NAME  = 'boundary_layer_level_for_UW_shlw',           &
        UNITS      = '1',                                          &
        DIMS       = MAPL_DimsHorzOnly,                            &
        VLOCATION  = MAPL_VLocationNone,                           &
        AVERAGING_INTERVAL = AVRGNINT,                             &
        REFRESH_INTERVAL   = RFRSHINT,                             &
        DEFAULT    = 72.,                                          &
                                                        RC=STATUS  )
     VERIFY_(STATUS)      

    call MAPL_AddImportSpec(GC,                                         &
         LONG_NAME  = 'aerosol_cloud_interaction',                             &
         UNITS      = '1',                                                     &
         SHORT_NAME = 'AERO_ACI',                                              &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationCenter,                                    &
         DATATYPE   = MAPL_StateItem,                                         &
         RC=STATUS  )
    VERIFY_(STATUS)      


    !new imports required for Aer-Cloud Interactions

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUGWX',                                    &
         LONG_NAME  = 'surface_eastward_gravity_wave_stress',      &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)   

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUGWY',                                    &
         LONG_NAME  = 'surface_northward_gravity_wave_stress',     &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)   


    call MAPL_AddImportSpec(GC,                             &
         LONG_NAME          = 'eastward_surface_stress_on_air',    &
         UNITS              = 'N m-2',                             &
         SHORT_NAME         = 'TAUX',                              &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         LONG_NAME          = 'northward_surface_stress_on_air',   &
         UNITS              = 'N m-2',                             &
         SHORT_NAME         = 'TAUY',                              &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUOROX',                                   &
         LONG_NAME  = 'surface_eastward_orographic_gravity_wave_stress',      &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TAUOROY',                                   &
         LONG_NAME  = 'surface_northward_orographic_gravity_wave_stress',     &
         UNITS      = 'N m-2',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddImportSpec ( gc,                                  &
         SHORT_NAME = 'OMEGA',                                     &
         LONG_NAME  = 'vertical_pressure_velocity',                &
         UNITS      = 'Pa s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
         LONG_NAME  = 'Blackadar_length_scale_for_scalars',                    &
         UNITS      = 'm',                                                     &
         SHORT_NAME = 'ALH',                                                   &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationEdge,                                      &
         RC=STATUS  )                                                   
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                   &
         SHORT_NAME = 'RADLW',                                           &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave',        &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( GC,                                   &
         SHORT_NAME = 'RADSW',                                           &
         LONG_NAME  = 'air_temperature_tendency_due_to_longwave',        &
         UNITS      = 'K s-1',                                           &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    !============================================

   call MAPL_AddImportSpec ( GC,                                          & !USe the nature run to force cirrus
         SHORT_NAME = 'WSUB_NATURE',                                 &
         LONG_NAME  =  'variance in wsub from the nature run',     &
         UNITS      = 'm2 s-2',                                    &
         RESTART    = MAPL_RestartSkip,                            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'TROPP',                                          &
         LONG_NAME  = 'tropopause_pressure_based_on_blended_estimate',  &
         UNITS      = 'Pa',                                             &
         DIMS       = MAPL_DimsHorzOnly,                                &
         VLOCATION  = MAPL_VLocationNone,                               &
         AVERAGING_INTERVAL = AVRGNINT,                                 &
         REFRESH_INTERVAL   = RFRSHINT,                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    !if ( IQVAINC /=0 ) then
    !   ! The following import is only for offline purposes
    !   ! NOTE: This is only used in offline applications so when adding new 
    !   !       fields to IMPORT state the suggestion is to add them BEFORE
    !   !       this state - unlike the usual procedure of always appending
    !   !       to the end of the state.
    !   call MAPL_AddImportSpec(GC,                                    &
    !        SHORT_NAME = 'QVAINC',                                    &
    !        LONG_NAME  = 'specific_humidity_analysis_increment',      &
    !        UNITS      = 'kg kg-1',                                   &
    !        default    = 0.0,                                         &
    !        DIMS       = MAPL_DimsHorzVert,                           &
    !        VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    !   VERIFY_(STATUS)                                                                          
    !endif

   call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'air_temperature',                           &
         UNITS      = 'K',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationCenter,                       &
                                                        RC=STATUS  )
   VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                    &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS    )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                      &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/s',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'QV_DYN_IN',                                 &
         LONG_NAME  = 'spec_humidity_at_begin_of_time_step',       &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'T_DYN_IN',                                 &
         LONG_NAME  = 'temperature_at_begin_of_time_step',       &
         UNITS      = 'K',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'U_DYN_IN',                                 &
         LONG_NAME  = 'u_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'V_DYN_IN',                                 &
         LONG_NAME  = 'v_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec ( gc,                                 &
         SHORT_NAME = 'PLE_DYN_IN',                                 &
         LONG_NAME  = 'edge_pressure_at_begin_of_time_step',       &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME         = 'AREA',                              &
        LONG_NAME          = 'agrid_cell_area',            &
        UNITS              = 'm+2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME         = 'USTAR',                             &
        LONG_NAME          = 'surface_velocity_scale',            &
        UNITS              = 'm s-1',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'TSTAR',                             &
        LONG_NAME          = 'surface_temperature_scale',         &
        UNITS              = 'K',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'QSTAR',                             &
        LONG_NAME          = 'surface_moisture_scale',            &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = '2-meter_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'T2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                         &
       LONG_NAME  = '2-meter_specific_humidity',                             &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'Q2M',                                                   &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'TA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
       LONG_NAME  = 'surface_air_specific_humidity',                         &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'QA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                                   &
        LONG_NAME          = 'sensible_heat_flux',                &
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SH',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation',                       &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'EVAP',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface geopotential height',       &
        UNITS              = 'm+2 s-2',                           &
        SHORT_NAME         = 'PHIS',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)
!-srf-gf-scheme

    ! !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                              &
       SHORT_NAME = 'WTHV2',                                                 &
       LONG_NAME  = 'Buoyancy_flux_for_SHOC',                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QCTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_total_cloud_water',         &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QLTOT',                                      &
         LONG_NAME  = 'grid_box_mass_fraction_of_cloud_liquid_water',        &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QITOT',                                      &
         LONG_NAME  = 'grid_box_mass_fraction_of_cloud_ice_water',           &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QRTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_falling_rain',              &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QSTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_falling_snow',              &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QGTOT',                                      &
         LONG_NAME  = 'mass_fraction_of_falling_graupel',              &
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'QPTOTLS',                                    &
         LONG_NAME  = 'mass_fraction_of_large_scale_falling_precip', & 
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTHDT ',                                     &
         LONG_NAME = 'pressure_weighted_potential_temperature_tendency_due_to_moist',&
         UNITS     = 'Pa K s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTDTFRIC',                                   &
         LONG_NAME = 'pressure_weighted_temperature_tendency_due_to_moist_friction',&
         UNITS     = 'Pa K s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DQDT  ',                                     &
         LONG_NAME = 'specific_humidity_tendency_due_to_moist',    &
         UNITS     = 'kg kg-1 s-1',                                &
         DIMS      = MAPL_DimsHorzVert,                           &
         VLOCATION = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DUDT  ',                                      &
         LONG_NAME = 'zonal_wind_tendency_due_to_moist',            &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &                  
         SHORT_NAME = 'DVDT  ',                                      &
         LONG_NAME = 'meridional_wind_tendency_due_to_moist',       &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DWDT  ',                                      &
         LONG_NAME = 'vertical_velocity_tendency_due_to_moist',       &
         UNITS     = 'm s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DTHDTCN',                                     &
         LONG_NAME = 'potential_temperature_tendency_due_to_convection',&
         UNITS     = 'K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQDTCN',                                      &
         LONG_NAME = 'specific_humidity_tendency_due_to_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQLDT ',                                      &
         LONG_NAME = 'total_liq_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DQIDT ',                                      &
         LONG_NAME = 'total_ice_water_tendency_due_to_moist',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQCDTCN ',                                    &
         LONG_NAME = 'condensate_tendency_due_to_convection',       &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_DQLDT ',                                  &
         LONG_NAME = 'convective_condensate_source',                &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'CNV_PRC3 ',                                   &
         LONG_NAME = 'convective_precipitation_from_RAS',           &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQRL   ',                                     &
         LONG_NAME = 'large_scale_rainwater_source',                &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DQRC   ',                                     &
         LONG_NAME = 'convective_rainwater_source',                 &
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_MF0',                                     &
         LONG_NAME = 'cloud_base_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               & 
         SHORT_NAME = 'CNV_MFD',                                     & 
         LONG_NAME = 'detraining_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &    
         DIMS      = MAPL_DimsHorzVert,                            &  
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &                  
         SHORT_NAME = 'CNV_MFC',                                     & 
         LONG_NAME = 'cumulative_mass_flux',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_FREQ',                                    & 
         LONG_NAME = 'convective_frequency',                        &
         UNITS     = 'fraction',                                    &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_BASEP',                                   & 
         LONG_NAME = 'pressure_at_convective_cloud_base',           &
         UNITS     = 'Pa',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'CNV_TOPP',                                    & 
         LONG_NAME = 'pressure_at_convective_cloud_top',            &
         UNITS     = 'Pa',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
         RC=STATUS  )
    VERIFY_(STATUS)



    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNV_UPDF',                                    &
         LONG_NAME = 'updraft_areal_fraction',                      &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNV_CVW',                                     &
         LONG_NAME = 'updraft_vertical_velocity',                   &
         UNITS     = 'hPa s-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )

    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CNV_QC',                                      &
         LONG_NAME ='grid_mean_convective_condensate',             &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'SHLW_PRC3 ',                                   &
         LONG_NAME = 'shallow_convective_rain',           &
         UNITS     = 'kg kg-1 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME= 'SHLW_SNO3 ',                                   &
         LONG_NAME = 'shallow_convective_snow',           &
         UNITS     = 'kg kg-1 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CUFRC_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_cloud_fraction',            &
         UNITS      = 'fraction',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQRDT_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_precipitating_condensate',            &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQCDT_SC',                             &
         LONG_NAME  = 'shallow_cumulus_condensate_tendency',            &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQSDT_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_precipitating_frozen_condensate',            &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DTHDT_SC',                                  &
         LONG_NAME  = 'Potential_temperature_tendency_from_shallow_convection',            &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQIDT_SC',                                  &
         LONG_NAME  = 'Ice_tendency_from_shallow_convection',      &
         UNITS      = 'kg ks-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQLDT_SC',                                      &
         LONG_NAME = 'Liquid_water_tendency_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQVDT_SC',                                      &
         LONG_NAME = 'Specific_humidity_tendency_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DUDT_SC',                                      &
         LONG_NAME = 'Zonal_wind_tendency_from_shallow_convection',&
         UNITS     = 'm s-2',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DVDT_SC',                                      &
         LONG_NAME = 'Meridional_wind_tendency_from_shallow_convection',&
         UNITS     = 'm s-2',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'WLCL_SC',                                      &
         LONG_NAME = 'Vertical_velocity_at_LCL_for_shallow_convection',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QTSRC_SC',                                      &
         LONG_NAME = 'Total_water_in_source_air_for_shallow_convection',&
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THLSRC_SC',                                      &
         LONG_NAME = 'Liquid_potential_temperature_of_source_air_for_shallow_convection',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THVLSRC_SC',                                      &
         LONG_NAME = 'Liquid_virtual_potential_temperature_source_air_shallow_convection',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'TKEAVG_SC',                                      &
         LONG_NAME = 'Average_boundary_layer_TKE_used_for_shallow_convection',&
         UNITS     = 'm2 s-2',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CLDTOP_SC',                                      &
         LONG_NAME = 'Cloud_top_height_from_shallow_convection',&
         UNITS     = 'm',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'WUP_SC',                                      &
         LONG_NAME = 'Vertical_velocity_in_shallow_convection_updraft',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QTUP_SC',                                      &
         LONG_NAME = 'Total_water_in_shallow_convection_updraft',&
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THLUP_SC',                                      &
         LONG_NAME = 'Liquid_potential_temperature_in_shallow_convection_updraft',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'THVUP_SC',                                      &
         LONG_NAME = 'Virtual_potential_temperature_in_shallow_convection_updraft',&
         UNITS     = 'K',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'UUP_SC',                                      &
         LONG_NAME = 'Zonal_wind_in_shallow_convection_updraft',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'VUP_SC',                                      &
         LONG_NAME = 'Meridional_wind_in_shallow_convection_updraft',&
         UNITS     = 'm s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ENTR_SC',                                      &
         LONG_NAME = 'Lateral_entrainment_rate_in_shallow_convection_updraft',&
         UNITS     = 'Pa-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DETR_SC',                                      &
         LONG_NAME = 'Lateral_detrainment_rate_in_shallow_convection_updraft',&
         UNITS     = 'Pa-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'XC_SC',                                      &
         LONG_NAME = 'Critical_mixing_fraction_in_shallow_updraft',&
         UNITS     = 'Pa-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLDET_SC',                                      &
         LONG_NAME = 'Detrained_liquid_water_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QIDET_SC',                                      &
         LONG_NAME = 'Detrained_ice_water_from_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLENT_SC',                                      &
         LONG_NAME = 'Sink_from_entrained_liquid_water_by_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QIENT_SC',                                      &
         LONG_NAME = 'Sink_from_entrained_ice_water_by_shallow_convection',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLSUB_SC',                                      &
         LONG_NAME = 'Shallow_convective_subsidence_liquid_water_tendency',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QISUB_SC',                                      &
         LONG_NAME = 'Shallow_convective_subsidence_ice_water_tendency',&
         UNITS     = 'kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNT_SC',                                      &
         LONG_NAME = 'Shallow_cloud_top_interface',            &
         UNITS     = 'index',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CNB_SC',                                      &
         LONG_NAME = 'Shallow_cloud_bottom_interface',            &
         UNITS     = 'index',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QCU_SC',                                      &
         LONG_NAME = 'Shallow_updraft_condensate',            &
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QLU_SC',                                      &
         LONG_NAME = 'Shallow_updraft_liquid_condensate',            &
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QIU_SC',                                      &
         LONG_NAME = 'Shallow_updraft_frozen_condensate',            &
         UNITS     = 'kg kg-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'UMF_SC',                                      &
         LONG_NAME = 'Shallow_updraft_mass_flux_at_interfaces', &
         UNITS     = 'kg m-2 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'MFD_SC',                                      &
         LONG_NAME = 'Shallow_updraft_detrained_mass_flux', &
         UNITS     = 'kg m-2 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CIN_SC',                                      &
         LONG_NAME = 'Convective_inhibition_for_shallow_convection', &
         UNITS     = 'J kg-1',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLCL_SC',                                      &
         LONG_NAME = 'Pressure_at_lift_condensation_level',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLFC_SC',                                      &
         LONG_NAME = 'Pressure_at_level_of_free_convection',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PINV_SC',                                      &
         LONG_NAME = 'Pressure_at_inversion_level',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PREL_SC',                                      &
         LONG_NAME = 'Pressure_at_release_level',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PBUP_SC',                                      &
         LONG_NAME = 'Pressure_at_level_neutral_buoyancy',   &
         UNITS     = 'Pa',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CBMF_SC',                                      &
         LONG_NAME = 'cloud_base_mass_flux_due_to_shallow_convection',&
         UNITS     = 'kg s-1 m-2',                                 &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='RL',                                          & 
         LONG_NAME ='liquid_cloud_particle_effective_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RI',                                          & 
         LONG_NAME = 'ice_phase_cloud_particle_effective_radius',   &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RR',                                          & 
         LONG_NAME = 'falling_rain_particle_effective_radius',      &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RS',                                          & 
         LONG_NAME  = 'falling_snow_particle_effective_radius',       &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

 call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RG',                                          & 
         LONG_NAME  = 'falling_graupel_particle_effective_radius',       &
         UNITS     = 'm',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
    
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CLDNCCN',                                     & 
         LONG_NAME ='number_concentration_of_cloud_particles',     &
         UNITS     ='m-3',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QSATI'  ,                                     & 
         LONG_NAME = 'saturation_spec_hum_over_ice',                &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QSATL'  ,                                     & 
         LONG_NAME = 'saturation_spec_hum_over_liquid',             &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ALPHT'  ,                                     & 
         LONG_NAME = 'pdf_spread_for_condensation_over_qsat_total', &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ALPH1'  ,                                     & 
         LONG_NAME ='pdf_spread_for_condensation_over_qsat_term1', &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ALPH2'  ,                                     & 
         LONG_NAME = 'pdf_spread_for_condensation_over_qsat_term2', &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CFPDFX'  ,                                    & 
         LONG_NAME = 'cloud_fraction_internal_in_PDF_scheme',       &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RHCLR'  ,                                     & 
         LONG_NAME = 'RH_clear_sky',                                &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CFPDF'  ,                                     & 
         LONG_NAME = 'cloud_fraction_after_PDF',                    &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FCLD'  ,                                      & 
         LONG_NAME = 'cloud_fraction_for_radiation',                &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='QV',                                          & 
         LONG_NAME ='water_vapor_for_radiation',                   &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QL',                                          & 
         LONG_NAME = 'in_cloud_cloud_liquid_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QI',                                          & 
         LONG_NAME = 'in_cloud_cloud_ice_for_radiation',                     &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QR',                                          & 
         LONG_NAME = 'Falling_rain_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QS',                                          & 
         LONG_NAME = 'Falling_snow_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'QG',                                          &
         LONG_NAME = 'Falling_graupel_for_radiation',                  &
         UNITS     = 'kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DBZ',                                          &
         LONG_NAME = 'Simulated_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DBZ_MAX',                                          &
         LONG_NAME = 'Maximum_simulated_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_RAIN',                                     &
         LONG_NAME ='falling_rain_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_SNOW',                                     &
         LONG_NAME ='falling_snow_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_ICE',                                     &
         LONG_NAME ='falling_ice_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='PRCP_GRAUPEL',                                     &
         LONG_NAME ='falling_graupel_precipitation_at_surface',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='LS_PRCP',                                     & 
         LONG_NAME ='nonanvil_large_scale_precipitation',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_PRCP',                                     & 
         LONG_NAME = 'anvil_precipitation',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CN_PRCP',                                     & 
         LONG_NAME ='convective_precipitation',                    &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='SC_PRCP',                                    & 
         LONG_NAME ='shallow_convective_precipitation',            &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ER_PRCP',                                     & 
         LONG_NAME = 'spurious_rain_from_RH_cleanup',          &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FILLNQV',                                     & 
         LONG_NAME = 'filling_of_negative_Q',          &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PGENTOT',                                     & 
         LONG_NAME = 'Total_column_production_of_precipitation',    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PREVTOT',                                     & 
         LONG_NAME = 'Total_column_re-evap/subl_of_precipitation',    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'LS_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_nonanvil_large_scale_showers',&
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_anvil_showers',             &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CN_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_convective_showers',        &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'SC_ARF',                                      & 
         LONG_NAME = 'areal_fraction_of_shallow_showers',        &
         UNITS     = '1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PTYPE',                                         &
         LONG_NAME = 'surface_precipitation_type',                 &
         UNITS     = '1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FRZR',                                         &
         LONG_NAME = 'freezing_rain_fall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'ICE',                                         &
         LONG_NAME = 'icefall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'SNO',                                         & 
         LONG_NAME = 'snowfall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PCU',                                         & 
         LONG_NAME = 'convective_rainfall',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLS',                                         & 
         LONG_NAME = 'large_scale_rainfall',                        &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TPREC',                                       & 
         LONG_NAME ='total_precipitation',                         &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='HOURNORAIN',                                  & 
         LONG_NAME ='time-during_an_hour_with_no_precipitation',   &
         UNITS     ='s',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TPW',                                         & 
         LONG_NAME ='total_precipitable_water',                    &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CCWP',                                        & 
         LONG_NAME ='grid_mean_conv_cond_water_path_diagnostic',   &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CWP',                                         & 
         LONG_NAME ='condensed_water_path',                        &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='LWP',                                         & 
         LONG_NAME ='liquid_water_path',                           &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='IWP',                                         & 
         LONG_NAME ='ice_water_path',                              &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='BYNCY',                                       & 
         LONG_NAME ='buoyancy_of surface_parcel',                  &
         UNITS     ='m s-2',                                       &
         DIMS      = MAPL_DimsHorzVert,                            & 
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CAPE',                                        & 
         LONG_NAME ='cape_for_surface_parcel',                     &
         UNITS     ='J kg-1',                                      &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='INHB',                                        & 
         LONG_NAME ='inhibition_for_surface_parcel',               &
         UNITS     ='J kg-1',                                      &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQ0',                                        & 
         LONG_NAME ='Total_Water_Substance_Before',                &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQ1',                                        & 
         LONG_NAME ='Total_Water_Substance_After',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCPTE',                                        & 
         LONG_NAME ='Total_VI_DcpT',                         &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVE0',                                        & 
         LONG_NAME ='Total_VI_MSE_Before',                         &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVE1',                                        & 
         LONG_NAME ='Total_VI_MSE_After',                          &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVEX',                                        & 
         LONG_NAME ='Total_VI_MSE_Somewhere',                      &
         UNITS     ='J m-2'  ,                                     &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZPBLCN',                                      & 
         LONG_NAME ='boundary_layer_depth',                        &
         UNITS     ='m'   ,                                        &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLCL',                                        & 
         LONG_NAME ='lifting_condensation_level',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLFC',                                        & 
         LONG_NAME ='level_of_free_convection',                    &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZCBL',                                        & 
         LONG_NAME ='height_of_cloud_base_layer',                  &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='MXDIAM',                                      & 
         LONG_NAME ='diameter_of_largest_RAS_plume',               &
         UNITS     ='m'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH600',                                      &
         LONG_NAME ='rh_at_600mb',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='Q600',                                      &
         LONG_NAME ='q_at_600mb',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QCBL',                                      &
         LONG_NAME ='q_at_cloud_base_level',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CNV_FRC',                                      &
         LONG_NAME ='convective_fraction',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='STOCH_CNV',                                     &
         LONG_NAME ='stochastic_factor_for_convection',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SIGMA_DEEP',                                     &
         LONG_NAME ='sigma_for_deep_in_convection',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SIGMA_MID',                                     &
         LONG_NAME ='sigma_for_mid_in_convection',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_TIME',                                     & 
         LONG_NAME ='timescale_for_RAS_plumes',               &
         UNITS     ='s'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_TRG',                                     &
         LONG_NAME ='rh_trigger_for_RAS_plumes',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_TOKI',                                     &
         LONG_NAME ='tokioka_factor_for_RAS_plumes',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_PBL',                                     &
         LONG_NAME ='pbl_fraction_for_RAS_plumes',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RAS_WFN',                                     &
         LONG_NAME ='RAS_work_function_before_scaling',               &
         UNITS     =''  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ENTLAM',                                      &
         LONG_NAME ='entrainment parameter',                       &
         UNITS     ='m-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

#if 0 
    ! taken out since they are now friendly to dynamics
    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='QLCN',                                       &
         LONG_NAME  ='mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='QICN',                                       &
         LONG_NAME  ='mass_fraction_of_convective_cloud_ice_water', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CLLS',                                       &
         LONG_NAME  ='large_scale_cloud_volume_fraction',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CLCN',                                       &
         LONG_NAME  ='convective_cloud_volume_fraction',             &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          
#endif


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH1',                                         & 
         LONG_NAME ='relative_humidity_before_moist',              &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RH2',                                         & 
         LONG_NAME ='relative_humidity_after_moist',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    !Outputs to give model trajectory in the moist TLM/ADJ
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TH_moist',                                    & 
         LONG_NAME ='potential_temp_before_moist',                 &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='Q_moist',                                     & 
         LONG_NAME ='specific_humidity_before_moist',              &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KCBL_moist',                                  & 
         LONG_NAME ='KCBL_before_moist',                           &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='ctop_moist',                                  & 
         LONG_NAME ='ctop_after_ras',                              &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TS_moist',                                    & 
         LONG_NAME ='surface_temp_before_moist',                   &
         UNITS     ='K',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KHu_moist',                                   &
         LONG_NAME ='upper_index_where_Kh_greater_than_2',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KHl_moist',                                   &
         LONG_NAME ='lower_index_where_Kh_greater_than_2',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    !End outputs for trajectory of moist TLM/ADJ

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='KHX',                                         &
         LONG_NAME ='scalar_diffusivity_for_pbl_clouds',               &
         UNITS     ='m+2 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='DTSX',                                   &
         LONG_NAME ='dts_for_pbl_clouds',         &
         UNITS     ='1',                                           &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RHX',                                         & 
         LONG_NAME ='relative_humidity_after_PDF',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REVSU_CN',                                    & 
         LONG_NAME ='evap_subl_of_convective_precipitation',       &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REVSU_LSAN',                                    & 
         LONG_NAME ='evap_subl_of_non_convective_precipitation',       &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_CN',                                      & 
         LONG_NAME ='evaporation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_CN_GF',                                   & 
         LONG_NAME ='evaporation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_SC',                                      & 
         LONG_NAME ='evaporation_of_shallow_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_AN',                                      & 
         LONG_NAME ='evaporation_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='REV_LS',                                      & 
         LONG_NAME ='evaporation_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_CN',                                      & 
         LONG_NAME ='sublimation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_CN_GF',                                   & 
         LONG_NAME ='sublimation_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_SC',                                      & 
         LONG_NAME ='sublimation_of_shallow_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_AN',                                      & 
         LONG_NAME ='sublimation_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RSU_LS',                                      & 
         LONG_NAME ='sublimation_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACR_TOT',                                      & 
         LONG_NAME ='total_accretion_of__precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_CN',                                      & 
         LONG_NAME ='liq_liq_accretion_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_SC',                                      & 
         LONG_NAME ='liq_liq_accretion_of_shallow_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_AN',                                      & 
         LONG_NAME ='liq_liq_accretion_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRLL_LS',                                      & 
         LONG_NAME ='liq_liq_accretion_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_CN',                                      & 
         LONG_NAME ='ice_liq_accretion_of_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_SC',                                      & 
         LONG_NAME ='ice_liq_accretion_of_shallow_convective_precipitation',     &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_AN',                                      & 
         LONG_NAME ='ice_liq_accretion_of_anvil_precipitation',          &
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ACRIL_LS',                                      & 
         LONG_NAME ='ice_liq_accretion_of_nonanvil_large_scale_precipitation',&
         UNITS     ='kg kg-1 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_CN',                                      & 
         LONG_NAME ='3D_flux_of_ice_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_CN_GF',                                   & 
         LONG_NAME ='3D_flux_of_ice_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_SC',                                      & 
         LONG_NAME ='3D_flux_of_ice_shallow_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_AN',                                      & 
         LONG_NAME ='3D_flux_of_ice_anvil_precipitation',          &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_LS',                                      & 
         LONG_NAME ='3D_flux_of_ice_nonanvil_large_scale_precipitation',&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFI_LSAN',                                      & 
         LONG_NAME ='3D_flux_of_ice_nonconvective_precipitation'  ,&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_CN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_CN_GF',                                    & 
         LONG_NAME ='3D_flux_of_liquid_convective_precipitation',   &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_SC',                                      & 
         LONG_NAME ='3D_flux_of_liquid_shallow_convective_precipitation',     &
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_AN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_anvil_precipitation',          &
         UNITS     ='kg m-2 s-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_LS',                                      & 
         LONG_NAME ='3D_flux_of_liquid_nonanvil_large_scale_precipitation',&
         UNITS     ='kg m-2 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFL_LSAN',                                      & 
         LONG_NAME ='3D_flux_of_liquid_nonconvective_precipitation'  ,&
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME='DPDTMST',                                      & 
         LONG_NAME ='layer_pressure_thickness_tendency_from_moist', &
         UNITS     ='Pa s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                             &
         VLOCATION = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCNVL',                                       & 
         LONG_NAME ='convective_source_of_cloud_liq',   &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DCNVI',                                       & 
         LONG_NAME ='convective_source_of_cloud_ice',        &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DLPDF',                                    & 
         LONG_NAME ='pdf_source_sink_of_cloud_liq',    &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIPDF',                                    & 
         LONG_NAME ='pdf_source_sink_of_cloud_ice',    &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DLFIX',                                    & 
         LONG_NAME ='fix_source_sink_of_cloud_liq',          &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIFIX',                                    & 
         LONG_NAME ='fix_source_sink_of_cloud_ice',          &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='AUT',                                      & 
         LONG_NAME ='autoconv_sink_of_cloud_liq',               &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='EVAPC',                                  & 
         LONG_NAME ='evaporation_of_cloud_liq',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SDM',                                    & 
         LONG_NAME ='sedimentation_sink_of_cloud_ice',        &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLICE_AN',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_anvil_snow',       &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLICE_LS',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_largescale_snow',  &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLWAT_AN',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_anvil_rain',     &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLWAT_LS',                              & 
         LONG_NAME ='autoconversion_fall_velocity_of_largescale_rain',&
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_AN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_anvil_rain',              &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_LS',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_largescale_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_CN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_convective_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLRN_SC',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_shallow_rain',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_AN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_anvil_snow',              &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_LS',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_largescale_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_CN',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_convective_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
         SHORT_NAME='VFALLSN_SC',                               & 
         LONG_NAME ='reevaporation_fall_velocity_of_shallow_snow',         &
         UNITS     ='m s-1',                                    &
         DIMS      = MAPL_DimsHorzVert,                         &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SUBLC',                                  & 
         LONG_NAME ='sublimation_of_cloud_ice',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='FRZ_TT',                                 & 
         LONG_NAME ='freezing_of_cloud_condensate',           &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='FRZ_PP',                                 & 
         LONG_NAME ='freezing_of_precip_condensate',          &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFRZ',                                 & 
         LONG_NAME ='Probability_of_freezing_of_aerosol_part',          &
         UNITS     ='1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='LIQANMOVE',                              &
!!$         LONG_NAME ='move2anv_source_of_anvil_liq',           &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='ICEANMOVE',                              & 
!!$         LONG_NAME ='move2anv_source_of_anvil_ice',           &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                                   &
!!$         SHORT_NAME = 'DLSCLD'  ,                                 & 
!!$         LONG_NAME = 'move2anv_change_in_large_scale_cloud_fraction',   &
!!$         UNITS     = '1',                                         &
!!$         DIMS      = MAPL_DimsHorzVert,                           &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                                   &
!!$         SHORT_NAME = 'DANCLD'  ,                                 & 
!!$         LONG_NAME = 'move2anv_change_in_anvil_cloud_fraction',   &
!!$         UNITS     = '1',                                         &
!!$         DIMS      = MAPL_DimsHorzVert,                           &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='CURAINMOVE',                             &
!!$         LONG_NAME ='movels2conv_source_of_cnv_rain',         &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)
!!$
!!$    call MAPL_AddExportSpec(GC,                               &
!!$         SHORT_NAME='CUSNOWMOVE',                             & 
!!$         LONG_NAME ='movels2conv_source_of_cnv_snow',         &
!!$         UNITS     ='kg kg-1 s-1',                            &
!!$         DIMS      = MAPL_DimsHorzVert,                       &
!!$         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
!!$    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFLCNMOVE',                              &
         LONG_NAME ='moved_source_of_cnv_rain',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PFICNMOVE',                              & 
         LONG_NAME ='moved_source_of_cnv_snow',               &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzVert,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CU2DRAINMOVE',                           &
         LONG_NAME ='moved_2d_source_of_cnv_rain',            &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzOnly,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='CU2DSNOWMOVE',                           & 
         LONG_NAME ='moved_2d_source_of_cnv_snow',            &
         UNITS     ='kg kg-1 s-1',                            &
         DIMS      = MAPL_DimsHorzOnly,                       &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! Vertically integrated water substance conversions


    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PDFLZ',                                        & 
         LONG_NAME ='statistical_source_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PDFIZ',                                        & 
         LONG_NAME ='statistical_source_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVRNZ',                                        & 
         LONG_NAME ='convective_production_of_rain_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVLZ',                                        & 
         LONG_NAME ='convective_source_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='CNVIZ',                                        & 
         LONG_NAME ='convective_source_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='EVPCZ',                                        & 
         LONG_NAME ='evaporation_loss_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SUBCZ',                                        & 
         LONG_NAME ='sublimation_loss_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='EVPPZ',                                        & 
         LONG_NAME ='evaporation_loss_of_precip_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SUBPZ',                                        & 
         LONG_NAME ='sublimation_loss_of_precip_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='AUTZ',                                        & 
         LONG_NAME ='autoconversion_loss_of_cloud_water',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SDMZ',                                        & 
         LONG_NAME ='sedimentation_loss_of_cloud_ice',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='COLLLZ',                                        & 
         LONG_NAME ='accretion_loss_of_cloud_water_to_rain',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='COLLIZ',                                        & 
         LONG_NAME ='accretion_loss_of_cloud_water_to_snow',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='FRZCZ',                                        & 
         LONG_NAME ='net_freezing_of_cloud_condensate',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='FRZPZ',                                        & 
         LONG_NAME ='net_freezing_of_precip_condensate',          &
         UNITS     ='kg m-2 s-1'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RCCODE',                                      & 
         LONG_NAME ='Convection_return_codes',                     &
         UNITS     ='codes',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TRIEDLV',                                     & 
         LONG_NAME ='Tested_for_convection_at_this_level',         &
         UNITS     ='0 or 1',                                      &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! MATMAT Exports for after-RAS inoutputs

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QVRAS',                                       & 
         LONG_NAME ='water_vapor_after_ras',                       &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THRAS',                                       & 
         LONG_NAME ='potential_temperature_after_ras',             &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='URAS',                                        & 
         LONG_NAME ='eastward_wind_after_ras',                     &
         UNITS     ='m s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='VRAS',                                        & 
         LONG_NAME ='northward_wind_after_ras',                    &
         UNITS     ='m s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    ! MATMAT Exports for before-RAS inputs for RAStest

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THOI',                                        & 
         LONG_NAME ='potential_temperature_before_ras',            &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QHOI',                                        & 
         LONG_NAME ='specific_humidity_before_ras',                &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QSSI',                                        & 
         LONG_NAME ='saturation_specific_humidity_before_ras',     &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='DQSI',                                        & 
         LONG_NAME ='deriv_sat_specific_humidity_wrt_t_before_ras',&
         UNITS     ='kg kg-1 K-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='PLEI',                                        & 
         LONG_NAME ='air_pressure_before_ras',                     &
         UNITS     ='Pa',                                          &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TPERTI',                                      & 
         LONG_NAME ='temperature_perturbation_before_ras',         &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='KCBLI',                                       & 
         LONG_NAME ='cloud_base_layer_before_ras',                 &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    ! !RECORD IMPORTS AND INTERNALS AT TOP OF MOIST:

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QX0',                                          &
         LONG_NAME  ='specific_humidity',                          &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLLSX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLLSX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLCNX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QLCNX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='CLLSX0',                                       &
         LONG_NAME  ='large_scale_cloud_area_fraction',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='CLCNX0',                                       &
         LONG_NAME  ='convective_cloud_area_fraction',             &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QILSX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QILSX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          


    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QICNX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_convective_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QICNX1',                                       &
         LONG_NAME  ='final_mass_fraction_of_convective_cloud_ice_water', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QCCNX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_convective_cloud_condensate', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='QCLSX0',                                       &
         LONG_NAME  ='initial_mass_fraction_of_large_scale_cloud_condensate', &
         UNITS      ='kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          



    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'KHX0',                                         &
         LONG_NAME  = 'scalar_diffusivity',                         &
         UNITS      = 'm+2 s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'THX0',                                        &
         LONG_NAME  = 'potential_temperature',                     &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'UX0',                                         &
         LONG_NAME  = 'eastward_wind',                             &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'VX0',                                         &
         LONG_NAME  = 'northward_wind',                            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'TSX0',                                        &
         LONG_NAME  = 'surface temperature',                       &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'FRLANDX0',                                    &
         LONG_NAME  = 'areal_land_fraction',                       &
         UNITS      = '1',                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         RC=STATUS  )
    VERIFY_(STATUS)


!!! downdraft diags
    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_MFC',                                   &
         LONG_NAME  = 'Downdraft_mass_flux',                       &
         UNITS      = 'kg m-2 s-1',                                &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_RH1',                                  &
         LONG_NAME  = 'Downdraft_in_cloud_RH_before',       &
         UNITS      = '1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_RH2',                                  &
         LONG_NAME  = 'Downdraft_in_cloud_RH_after',       &
         UNITS      = '1',                                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_TC',                                    &
         LONG_NAME  = 'Temperature_excess_in_DDF',                 &
         UNITS      = 'K',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_QVC',                                   &
         LONG_NAME  = 'Spec_hum_excess_in_DDF',                    &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_BYNC',                                  &
         LONG_NAME  = 'Buoyancy_of_DDF',                           &
         UNITS      = 'm s-2',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_MUPH',                                  &
         LONG_NAME  = 'Downdraft_moistening_from_evap_subl',       &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_DQDT',                                  &
         LONG_NAME  = 'Total_Downdraft_moistening',                      &
         UNITS      = 'kg kg-1 s-1',                               &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_DTDT',                                  &
         LONG_NAME  = 'Total_Downdraft_heating',                         &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'DDF_ZSCALE',                                &
         LONG_NAME  = 'vertical_scale_for_downdraft',              &
         UNITS      = 'm',                                         &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,                        &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEDISS',                                    &
         LONG_NAME  = 'kinetic_energy_diss_in_RAS',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEMST',                                    &
         LONG_NAME  = 'vertically_integrated_kinetic_energy_tendency_across_moist',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'KEMST2',                                    &
         LONG_NAME  = 'vertically_integrated_KE_dissipation_in_RAS',       &
         UNITS      = 'W m-2',                                    &
         DIMS       = MAPL_DimsHorzOnly,                          &
         VLOCATION  = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    ! CAR 12/5/08
    ! Aerosol Scavenging diagnostics/export states
    ! ------------------------------
    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDUDT ',                                         & 
         LONG_NAME ='dust_tendency_due_to_conv_scav',                 &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSSDT ',                                         &
         LONG_NAME ='sea_salt_tendency_due_to_conv_scav',             &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DOCDT ',                                         &
         LONG_NAME ='organic_carbon_tendency_due_to_conv_scav',       &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DBCDT ',                                         &
         LONG_NAME ='black_carbon_tendency_due_to_conv_scav',         &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSUDT ',                                         &
         LONG_NAME ='sulfate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DNIDT',                                          &
         LONG_NAME ='nitrate_tendency_due_to_conv_scav',              &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DDUDTcarma ',                                    &
         LONG_NAME ='carma_dust_tendency_due_to_conv_scav',           &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='DSSDTcarma ',                                    &
         LONG_NAME ='carma_seasalt_tendency_due_to_conv_scav',        &
         UNITS     ='kg m-2 s-1',                                     &
         DIMS      = MAPL_DimsHorzOnly,                               &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='LFR',                                            &
         LONG_NAME ='lightning_flash_rate',                           &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A1X1',                                           &
         LONG_NAME ='LFR_Term_number_1',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A2X2',                                           &
         LONG_NAME ='LFR_Term_number_2',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A3X3',                                           &
         LONG_NAME ='LFR_Term_number_3',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A4X4',                                           &
         LONG_NAME ='LFR_Term_number_4',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='A5X5',                                           &
         LONG_NAME ='LFR_Term_number_5',                              &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='LFR_GCC',                                        &
         LONG_NAME ='lightning_flash_rate_for_GEOSCHEMchem',          &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THMOIST',                                     & 
         LONG_NAME ='potential_temperature_after_all_of_moist',   &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SMOIST',                                      & 
         LONG_NAME ='dry_static_energy_after_all_of_moist',        &
         UNITS     ='m+2 s-2',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    !-------Aerosol Cloud Interactions Diagnostics  

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SMAX_LIQ',                                          & 
         LONG_NAME ='Maximum incloud supersaturation for liquid',        &
         UNITS     ='%',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='WSUB',                                          & 
         LONG_NAME ='Subgrid Scale in-cloud vertical velocity',          &
         UNITS     ='m s-1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                        &
         SHORT_NAME='CCN01',                                             & 
         LONG_NAME ='CCN conc at 0.1 % supersaturation (grid_avg)',                 &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CCN04',                                             & 
         LONG_NAME ='CCN conc at 0.4 % supersaturation (grid_avg)',                 &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CCN1',                                              & 
         LONG_NAME ='CCN conc at 1.0 % supersaturation (grid_avg)',                 &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SMAX_ICE',                                          & 
         LONG_NAME ='Maximum incloud supersaturation for ice',        &
         UNITS     ='%',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CDNC_NUC',                                          & 
         LONG_NAME ='Nucleated cloud droplet concentration (grid_avg)',        &
         UNITS     ='m-3',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='INC_NUC',                                          & 
         LONG_NAME ='Nucleated ice crystal concentration (grid_avg)',        &
         UNITS     ='m-3',                                                 &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NCPL_VOL',                                       &
         LONG_NAME  ='particle_number_for_liquid_cloud', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                                          

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NCPI_VOL',                                       &
         LONG_NAME  ='particle_number_for_ice_cloud', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SO4',                                       &
         LONG_NAME  ='Sulfate number conc.', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='ORG',                                       &
         LONG_NAME  ='Organic number conc. (hydrophilic)',        &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='BCARBON',                                       &
         LONG_NAME  ='Black carbon number conc. (hydrophilic)',        &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='DUST',                                       &
         LONG_NAME  ='Total dust number conc', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SEASALT',                                       &
         LONG_NAME  ='Total sea number conc', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NHET_NUC',                                       &
         LONG_NAME  ='Nucleated ice crystal concentration by het freezing (grid_avg)', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NLIM_NUC',                                       &
         LONG_NAME  ='Limiting IN concentration allowing hom freezing', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SAT_RAT',                                         & 
         LONG_NAME ='saturation_ratio_after_moist',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 



    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQVDT_micro',                                         & 
         LONG_NAME ='Q tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQLDT_micro',                                         & 
         LONG_NAME ='QL tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQIDT_micro',                                         & 
         LONG_NAME ='QI tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQRDT_micro',                                         &
         LONG_NAME ='QR tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQSDT_micro',                                         &
         LONG_NAME ='QS tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQGDT_micro',                                         &
         LONG_NAME ='QG tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQADT_micro',                                         &
         LONG_NAME ='QA tendency due to microphysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_micro',                                         & 
         LONG_NAME ='T tendency due to microphysics ',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DUDT_micro',                                         &
         LONG_NAME ='U tendency due to microphysics ',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DVDT_micro',                                         &
         LONG_NAME ='V tendency due to microphysics ',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQX1',                                        & 
         LONG_NAME ='Total_Water_Substance_bef_macro',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVQX2',                                        & 
         LONG_NAME ='Total_Water_Substance_bef_micro',                 &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RL_MASK',                                          & 
         LONG_NAME ='volumetric_liquid_cloud_particle_volume_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='RI_MASK',                                          & 
         LONG_NAME ='volumetric_ice_cloud_particle_volume_radius',      &
         UNITS     ='m',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='KAPPA',                                       & 
         LONG_NAME ='kappa parameter for activation',              &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CFLIQ',                                       &
         LONG_NAME  ='liquid_cloud_area_fraction (LS)',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)      

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='CFICE',                                       &
         LONG_NAME  ='ice_cloud_area_fraction (LS)',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS) 

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='RHICE',                                       &
         LONG_NAME  ='Relative humidity wrt ice',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='RHLIQ',                                       &
         LONG_NAME  ='Relative humidity wrt liquid',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME ='SC_ICE',                                       &
         LONG_NAME  ='Effective Freezing RHi',            &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NHET_IMM',                                       &
         LONG_NAME  ='Immersion_IN', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='NHET_DEP',                                       &
         LONG_NAME  ='Deposition_IN', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='DUST_IMM',                                       &
         LONG_NAME  ='Immersion IN from dust', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='DUST_DEP',                                       &
         LONG_NAME  ='deposition IN from dust', &
         UNITS      ='m-3',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SCF',                                       &
         LONG_NAME  ='Supercooled cloud fraction', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                           &
         SHORT_NAME ='SCF_ALL',                                       &
         LONG_NAME  ='Supercooled cloud fraction including snow and rain', &
         UNITS      ='1',                                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_GW',                                          & 
         LONG_NAME ='Subgrid Scale vertical velocity variance from GW',  &
         UNITS     ='m2 s-2',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_CNV',                                          & 
         LONG_NAME ='Subgrid Scale vertical velocity variance from convection',  &
         UNITS     ='m2 s-2',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_TURB',                                          & 
         LONG_NAME ='Subgrid Scale vertical velocity variance from turbulence', &
         UNITS     ='m2 s-2',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SIGW_RC',                                          & 
         LONG_NAME ='Mean subgrid Scale vertical velocity from rad cooling', &
         UNITS     ='m s-1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CNV_FICE',                                          & 
         LONG_NAME ='Ice fraction in convective tower', &
         UNITS     ='1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CNV_NDROP',                                          & 
         LONG_NAME ='Droplet number conc. in conv. detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='CNV_NICE',                                          & 
         LONG_NAME ='Ice crystal number conc. in conv. detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SC_NDROP',                                          & 
         LONG_NAME ='Droplet number conc. in shallow detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='SC_NICE',                                          & 
         LONG_NAME ='Ice crystal number conc. in shallow detrainment', &
         UNITS     ='m-3',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='RHCmicro',                                          & 
         LONG_NAME ='Corrected RHc after micro', &
         UNITS     ='1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='CCNCOLUMN',                                   & 
         LONG_NAME ='Vertically integrated CCN at 1% ssat',        &
         UNITS     ='m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NDCOLUMN',                                   & 
         LONG_NAME ='Vertically integrated NCPL',        &
         UNITS     ='m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCCOLUMN',                                   & 
         LONG_NAME ='Vertically integrated NCPI',        &
         UNITS     ='m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='BERG',                                                 &
         LONG_NAME  ='ice mixing ration tendency due to Bergeron process ',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='BERGS',                                                 &
         LONG_NAME  ='Snow mixing ration tendency due to Bergeron process ',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='MELT',                                                 &
         LONG_NAME  ='Melting of cloud ice',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='QCRES',                                                 &
         LONG_NAME  ='Residual cloud tendency in micro',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='QIRES',                                                 &
         LONG_NAME  ='Residual ice cloud tendency in micro',  &
         UNITS      ='Kg Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='AUTICE',                                      & 
         LONG_NAME ='autoconv_sink_of_cloud_ice',               &
         UNITS     ='kg kg-1 s-1',                                   &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DT_RASP',                                                 &
         LONG_NAME  ='T tendency from ras precip',  &
         UNITS      ='K s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNHET_IMM',                                                 &
         LONG_NAME  ='Ice number tendency due to immersion freezing',  &
         UNITS      ='Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNHET_CT',                                                 &
         LONG_NAME  ='Ice number tendency due to contact freezing',  &
         UNITS      ='Kg-1 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_macro',                                         & 
         LONG_NAME ='T tendency due to macrophysics ',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQVDT_macro',                                         &
         LONG_NAME ='QV tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQLDT_macro',                                         &
         LONG_NAME ='QL tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQIDT_macro',                                         &
         LONG_NAME ='QI tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQADT_macro',                                         &
         LONG_NAME ='QA tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_moist',                                         & 
         LONG_NAME ='T tendency due to moist',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDTCN',                                         & 
         LONG_NAME ='T tendency due to convection',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)  

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='FRZPP_LS',                                                 &
         LONG_NAME  ='Tendency in T from micro precip freezing',  &
         UNITS      ='K s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='SNOWMELT_LS',                                                 &
         LONG_NAME  ='Tendency in T from micro snow melting',  &
         UNITS      ='K s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCNUC',                                                 &
         LONG_NAME  ='Ice number tendency due to nucleation on aerosol',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCHMSPLIT',                                                 &
         LONG_NAME  ='Ice number tendency due to H-M splittering',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCSUBL',                                                 &
         LONG_NAME  ='Ice number tendency due to sublimation',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCAUTICE',                                                 &
         LONG_NAME  ='Ice number tendency due to autoconversion to snow',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCACRIS',                                                 &
         LONG_NAME  ='Ice number tendency due to accretion by snow',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNCCNV',                                                 &
         LONG_NAME  ='Ice crystal number tendency from convective detrainment',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDCCN',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to activation',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDACRLS',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to accretion by snow',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDEVAPC',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to evaporation',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDACRLR',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to accretion by rain',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDAUTLIQ',                                                 &
         LONG_NAME  ='Cloud droplet number tendency due to autoconversion',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME ='DNDCNV',                                                 &
         LONG_NAME  ='Cloud droplet number tendency from convective detrainment',  &
         UNITS      ='m-3 s-1',                                            &
         DIMS       = MAPL_DimsHorzVert,                                     &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='CLDREFFL_TOP',                                   & 
         LONG_NAME ='Droplet effective radius at cloud top',        &
         UNITS     ='m'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='CLDREFFI_TOP',                                   & 
         LONG_NAME ='ice crystal effective radius at cloud top',        &
         UNITS     ='m'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCPL_TOP',                                   & 
         LONG_NAME ='Grid-averaged NCPL at cloud top',        &
         UNITS     ='m-3'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCPI_TOP',                                   & 
         LONG_NAME ='Grid-averaged NCPI at cloud top',        &
         UNITS     ='m-3'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME='NCPL_CLDBASE',                                   & 
         LONG_NAME ='IN-CLOUD NCPL at cloud base',        &
         UNITS     ='m-3'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)
    


    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'LWC',                                      &
         LONG_NAME  = 'liquid water content',        &
         UNITS      = 'kg m-3',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME = 'IWC',                                      &
         LONG_NAME  = 'ice water content',        &
         UNITS      = 'kg m-3',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                             &
         SHORT_NAME='QCVAR_EXP',                                   & 
         LONG_NAME ='inverse relative variance of cloud water',        &
         UNITS      = '1',                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)
        
                
    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='LTS',                                          & 
         LONG_NAME ='Lower tropospheric stability', &
         UNITS     ='K',                                             &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)
        
      call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='EIS',                                          & 
         LONG_NAME ='Estimated Inversion Strength', &
         UNITS     ='K',                                             &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
     VERIFY_(STATUS)
    
    
     call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='RAS_ALPHA',                                          & 
         LONG_NAME ='RAS relaxation parameter', &
         UNITS     ='1',                                             &
         DIMS      = MAPL_DimsHorzVert,                                  &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                        &
         SHORT_NAME='RAS_TAU',                                          & 
         LONG_NAME ='RAS total relaxation timescale',                   &
         UNITS     ='1',                                                &
         DIMS      = MAPL_DimsHorzVert,                                 &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                                        &
         SHORT_NAME = 'DTDT_BL',                                       &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_bound_layer', &
         UNITS      = 'K s-1',                                         &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddImportSpec(GC,                                       &
         SHORT_NAME = 'DQDT_BL',                                      &
         LONG_NAME  = 'tendency_of_spec_humidity_due_to_bound_layer',  &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DNDDT ',                                      &
         LONG_NAME = 'total_liq_droplet_number_tendency_due_to_moist',       &
         UNITS     = 'm-3 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'DNCDT ',                                      &
         LONG_NAME = 'total_ice_crystal_number_tendency_due_to_moist',       &
         UNITS     = 'm-3 s-1',                                 &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'surface_air_temperature',                               &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'TA',                                                    &
       DIMS       = MAPL_DimsHorzOnly,                                       &
       VLOCATION  = MAPL_VLocationNone,                                      &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

!-srf-gf-scheme
    !IF(DEBUG_GF==1 .and. (CONVPAR_OPTION=='GF' .or. CONVPAR_OPTION=='BOTH')) THEN

    !   call MAPL_AddExportSpec(GC,                                     &
    !     SHORT_NAME = 'DTRDT_GF',                                    &
    !     LONG_NAME  = 'tendency_of_tracer_due_GF',        &
    !     UNITS      = 'kg kg-1 s-1',                                  &
    !     DIMS       = MAPL_DimsHorzVert,                              &
    !     VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    !   VERIFY_(STATUS)    
    !   call MAPL_AddExportSpec(GC,                                     &
    !     SHORT_NAME = 'DQDT_GF',                                    &
    !     LONG_NAME  = 'tendency_of_spec_humidity_due_GF',        &
    !     UNITS      = 'kg kg-1 s-1',                                  &
    !     DIMS       = MAPL_DimsHorzVert,                              &
    !     VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                     &
    !     SHORT_NAME = 'DTDT_GF',                                    &
    !     LONG_NAME  = 'tendency_of_temp_due_GF',        &
    !     UNITS      = 'K s-1',                                  &
    !     DIMS       = MAPL_DimsHorzVert,                              &
    !     VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                     &
    !     SHORT_NAME = 'MUPDP',                                    &
    !     LONG_NAME  = 'Mass_flux_deep_GF',        &
    !     UNITS      = 'kg m-2 s-1',                                  &
    !     DIMS       = MAPL_DimsHorzVert,                              &
    !     VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                     &
    !     SHORT_NAME = 'MUPSH',                                    &
    !     LONG_NAME  = 'Mass_flux_shallow_GF',        &
    !     UNITS      = 'kg m-2 s-1',                                  &
    !     DIMS       = MAPL_DimsHorzVert,                              &
    !     VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                     &
    !     SHORT_NAME = 'MUPMD',                                    &
    !     LONG_NAME  = 'Mass_flux_congestus_GF',        &
    !     UNITS      = 'kg m-2 s-1',                                  &
    !     DIMS       = MAPL_DimsHorzVert,                              &
    !     VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    !   VERIFY_(STATUS)
    !   !-2d
    !   call MAPL_AddExportSpec(GC,                                &
    !    SHORT_NAME         = 'MFDP',                              &
    !    LONG_NAME          = 'mass_flux_cloud_base_deep_GF',      &
    !    UNITS              = 'kg m-2 s-1',                        &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                &
    !    SHORT_NAME         = 'MFSH',                              &
    !    LONG_NAME          = 'mass_flux_cloud_base_shal_GF',      &
    !    UNITS              = 'kg m-2 s-1',                        &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                &
    !    SHORT_NAME         = 'MFMD',                              &
    !    LONG_NAME          = 'mass_flux_cloud_base_mid_GF',       &
    !    UNITS              = 'kg m-2 s-1',                        &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                &
    !    SHORT_NAME         = 'ERRDP',                             &
    !    LONG_NAME          = 'convection_code_deep_GF',                   & 
    !    UNITS              = '1',                                &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !   VERIFY_(STATUS)
    !   call MAPL_AddExportSpec(GC,                                &
    !    SHORT_NAME         = 'ERRSH',                             &
    !    LONG_NAME          = 'convection_code_shallow_GF',                   & 
    !    UNITS              = '1',                                &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)
    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'ERRMD',                             &
    !    LONG_NAME          = 'convection_code_mid_GF',            & 
    !    UNITS              = '1',                                &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)
    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'AA0',                               &
    !    LONG_NAME          = 'cloud work function 0',             & 
    !    UNITS              = 'J kg-1',                              &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)

    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'AA1',                               &
    !    LONG_NAME          = 'cloud work function 1',             & 
    !    UNITS              = 'J kg-1',                              &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)

    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'AA2',                               &
    !    LONG_NAME          = 'cloud work function 2',             & 
    !    UNITS              = 'J kg-1',                              &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)

    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'AA3',                               &
    !    LONG_NAME          = 'cloud work function 3',             & 
    !    UNITS              = 'J kg-1',                              &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)
    !    
    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'AA1_CIN',                           &
    !    LONG_NAME          = 'cloud work function CIN',           & 
    !    UNITS              = 'J kg-1',                              &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)

    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'AA1_BL',                            &
    !    LONG_NAME          = 'Bound layer AA1',                   & 
    !    UNITS              = 'J kg-1 s-1',                            &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)
    !    
    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'TAU_BL',                             &
    !    LONG_NAME          = 'Bound layer time scale',            & 
    !    UNITS              = 's',                                &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)
    !    
    !    call MAPL_AddExportSpec(GC,                               &
    !    SHORT_NAME         = 'TAU_EC',                             &
    !    LONG_NAME          = 'cape removal time scale',           & 
    !    UNITS              = 's',                                &
    !    DIMS               = MAPL_DimsHorzOnly,                   &
    !    VLOCATION          = MAPL_VLocationNone,                  &
    !                                                   RC=STATUS  )
    !    VERIFY_(STATUS)

    !ENDIF
!-srf-gf-scheme
!========================================================================================
!--kml--- activation for single-moment uphysics
        !call MAPL_AddExportSpec(GC,                               &
        !SHORT_NAME         = 'NACTL',                             &
        !LONG_NAME          = 'activ aero # conc liq phase for 1-mom',              & 
        !UNITS              = 'm-3',                               &
        !DIMS               = MAPL_DimsHorzVert,                   &
        !VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
        !VERIFY_(STATUS)

        !call MAPL_AddExportSpec(GC,                               &
        !SHORT_NAME         = 'NACTI',                             &
        !LONG_NAME          = 'activ aero # conv ice phase for 1-mom',               & 
        !UNITS              = 'm-3',                               &
        !DIMS               = MAPL_DimsHorzVert,                   &
        !VLOCATION          = MAPL_VLocationCenter,     RC=STATUS  )
        !VERIFY_(STATUS)
!--kml--- activation for single-moment uphysics

!========================================================================================             


    !EOS
