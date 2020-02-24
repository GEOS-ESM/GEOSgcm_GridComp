! $Id$

! VERIFY_ and RETURN_ macros for error handling.

#include "MAPL_Generic.h"

!=============================================================================
module GEOS_DatmoDynGridCompMod
!BOP

! !MODULE: GEOS_Singcol -- A Module to drive single column model with profile data.


! !USES:

  use ESMF
  use MAPL_Mod
  use PPM
  use cfmip_data_mod
  
  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt MOIST} 
!

!EOP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code
    
! !DESCRIPTION: This version uses the GEOS\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the GEOS\_GenericState.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (ESMF_Config          )            :: CF

    integer      :: MY_STEP
    integer      :: ACCUMINT
    real         :: DT


!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    
    Iam = trim(COMP_NAME) // 'SetServices'
    call write_parallel(trim(IAM))

! Register services for this component
! ------------------------------------
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   Run , __RC__)
    call MAPL_GridCompSetEntryPoint ( gc, ESMF_METHOD_RUN,   RunAddIncs, __RC__)
    
! Get the configuration from the component
!-----------------------------------------
    call ESMF_GridCompGet( GC, CONFIG = CF, __RC__ )
    

! Set the state variable specs.
! -----------------------------
    call ESMF_ConfigGetAttribute ( CF, DT,  Label="RUN_DT:", __RC__)
    call ESMF_ConfigGetAttribute ( CF, DT,  Label=trim(COMP_NAME)//"_DT:", default=DT, __RC__)

    MY_STEP = nint(DT)

    call ESMF_ConfigGetAttribute ( CF, DT, Label=trim(COMP_NAME)//'Avrg:',default=DT, __RC__)

    ACCUMINT = nint(DT)


!BOS
! !INTERNAL STATE:

     call MAPL_AddInternalSpec(GC,                                 &
         SHORT_NAME = 'PREF',                                      &
         LONG_NAME  = 'reference_air_pressure',                    &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         AVERAGING_INTERVAL = ACCUMINT,                            &
         REFRESH_INTERVAL   = MY_STEP,                             &
         ADD2EXPORT = .TRUE.,                                      &
                                                           __RC__  )

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME='PLE',                                         &
         LONG_NAME ='Pressure_at_the_edges ',                      &
         UNITS     ='Pascals',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                   __RC__  )
  
    call MAPL_AddInternalSpec(GC,                               &
         SHORT_NAME='T ',                                      &
         LONG_NAME ='air_temperature',                    &
         UNITS     ='K',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddInternalSpec(GC,                               &
         SHORT_NAME='U  ',                                      &
         LONG_NAME ='Zonal wind ',        &
         UNITS     ='m/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddInternalSpec(GC,                               &
         SHORT_NAME='V  ',                                      &
         LONG_NAME ='meridional wind',        &
         UNITS     ='m/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddInternalSpec(GC,                               &
         SHORT_NAME='OM ',                                      &
         LONG_NAME ='pressure velocity',                        &
         UNITS     ='Pa/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
                                                        __RC__  )
    

! !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME ='DTDT',                                       &
         LONG_NAME  ='T_tendency',                                 &
         UNITS      ='K s-1',                                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                        __RC__  )
                                                                              
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME ='DUDT',                                       &
         LONG_NAME  ='later',                                      &
         UNITS      ='m s-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                         __RC__  )
                                                                              
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME ='DVDT',                                       &
         LONG_NAME  ='later',                                      &
         UNITS      ='m s-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                        __RC__  )
                                                                              
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME ='DWDT',                                       &
         LONG_NAME  ='later',                                      &
         UNITS      ='m s-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
                                                        __RC__  )
                                                                              
    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME ='DPEDT',                                      &
         LONG_NAME  ='air_pressure',                               &
         UNITS      ='Pa',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
                                                        __RC__  )
                                                                              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Not sure why these are really needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQVANA',                                    &
         LONG_NAME  = 'specific_humidity_vapor_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQLANA',                                    &
         LONG_NAME  = 'specific_humidity_liquid_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQIANA',                                    &
         LONG_NAME  = 'specific_humidity_ice_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQRANA',                                    &
         LONG_NAME  = 'specific_humidity_rain_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQSANA',                                    &
         LONG_NAME  = 'specific_humidity_snow_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DQGANA',                                    &
         LONG_NAME  = 'specific_humidity_graupel_increment_from_analysis', &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'DOXANA',                                    &
         LONG_NAME  = 'ozone_increment_from_analysis',             &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
    call MAPL_AddImportSpec ( gc,                             &
         SHORT_NAME = 'PHIS',                                      &
         LONG_NAME  = 'surface_geopotential_height',               &
         UNITS      = 'm+2 sec-2',                                 &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               __RC__  )
    
!!This bundle should come from PHYSICS(MOIST)
    call MAPL_AddImportSpec( gc,                              &
        SHORT_NAME = 'TRADV',                                        &
        LONG_NAME  = 'advected_quantities',                        &
        UNITS      = 'unknown',                                    &
        DATATYPE   = MAPL_BundleItem,               &
        __RC__  )
    
! !EXPORT STATE:
 
! first copies of exportable internal variables
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PLE',                                         &
         LONG_NAME ='Pressure at the edges ',                      &
         UNITS     ='Pascals',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='T  ',                                           &
         LONG_NAME ='air_temperature',                    &
         UNITS     ='K',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='U  ',                                      &
         LONG_NAME ='Zonal wind ',        &
         UNITS     ='m/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='V',                                      &
         LONG_NAME ='meridional wind',        &
         UNITS     ='m/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='W',                                      &
         LONG_NAME ='meridional wind',        &
         UNITS     ='m/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='Q',                                      &
         LONG_NAME ='meridional wind',        &
         UNITS     ='m/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='OMEGA',                                      &
         LONG_NAME ='pressure velocity',                        &
         UNITS     ='Pa/s',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DIV',                                        &
         LONG_NAME ='divergence',        &
         UNITS     ='s-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='LHOBS',                                       &
         LONG_NAME ='Obs. latent heat flux (surface)',             &
         UNITS     ='W m-2',                                       &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SHOBS',                                       &
         LONG_NAME ='Obs. sensible heat flux (surface)',           &
         UNITS     ='W m-2',                                       &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PCPOBS',                                      &
         LONG_NAME ='Obs. precipitation rate',                     &
         UNITS     ='mm/d',                                        &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TSAIROBS',                                    &
         LONG_NAME ='Obs. sfc(2m?) air temp.',                     &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        __RC__  )
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TGSOILOBS',                                   &
         LONG_NAME ='Obs. Soil Temp.',                             &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PSFCOBS',                                     &
         LONG_NAME ='Obs. Sfc. Pressure',                          &
         UNITS     ='[Pa,hPa]',                                    &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TH ',                                      &
         LONG_NAME ='potential_temperature',                    &
         UNITS     ='K',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='ZLE',                                         &
         LONG_NAME ='Geop. height at the edges',                      &
         UNITS     ='m',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QOBS',                                        &
         LONG_NAME ='Obs. Spec. Humidity',                         &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TOBS',                                        &
         LONG_NAME ='Obs. Temperature',                            &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TVATOBS',                                     &
         LONG_NAME ='Obs. Temperature Tendency V adv',             &
         UNITS     ='K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='THATOBS',                                     &
         LONG_NAME ='Obs. Temperature Tendency H adv',             &
         UNITS     ='K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QVATOBS',                                     &
         LONG_NAME ='Obs. Moisture Tendency V adv',             &
         UNITS     ='K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QHATOBS',                                     &
         LONG_NAME ='Obs. Moisture Tendency H adv',             &
         UNITS     ='K s-1',                                       &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QTEST ',                                      &
         LONG_NAME ='test_tracer',                    &
         UNITS     ='1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='UE',                                      &
         LONG_NAME ='Diagnosed_Edge_Winds',                    &
         UNITS     ='m s-1',                                        &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='WWTG',                                      &
         LONG_NAME ='weak_T-gradient_compensating_W',                    &
         UNITS     ='m s-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
                                                        __RC__  )
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='WTOT',                                      &
         LONG_NAME ='total_vertical_velocity',                    &
         UNITS     ='m s-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationEdge,                         &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QLOBS',                                       &
         LONG_NAME ='Obs. Cloud liquid',                           &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='QIOBS',                                       &
         LONG_NAME ='Obs. Cloud ice',                              &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                __RC__  )
    

!JTB: New exports for Surface comp.
!------------------------------------------
    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='PS',                                           &
         LONG_NAME ='Surface Pressure',                             &
         UNITS     ='Pa',                                           &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='TA',                                           &
         LONG_NAME ='Surface_air_temperature',                      &
         UNITS     ='K',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='SPEED',                                        &
         LONG_NAME ='Surface_wind_speed',                           &
         UNITS     ='m s-1',                                        &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='DZ',                                           &
         LONG_NAME ='Surface_layer_height',                         &
         UNITS     ='m',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='QA',                                           &
         LONG_NAME ='Surface_air_spec_humidity',                    &
         UNITS     ='1',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='QSKINOBS',                                     &
         LONG_NAME ='obs_skin_spec_humidity_whatever_that_means',   &
         UNITS     ='1',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='TSKINOBS',                                     &
         LONG_NAME ='obs_skin_temperature',                         &
         UNITS     ='K',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='PHIS',                                     &
         LONG_NAME ='obs_skin_temperature',                         &
         UNITS     ='K',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='VARFLT',                                     &
         LONG_NAME ='obs_skin_temperature',                         &
         UNITS     ='K',                                            &
         DIMS      = MAPL_DimsHorzOnly,                             &
         VLOCATION = MAPL_VLocationNone,                            &
                                                        __RC__  )
 !-srf-gf-scheme
    call MAPL_AddExportSpec(GC,                                                    &
         SHORT_NAME = 'DYNF_Q',                                                    &
         LONG_NAME  = 'dynamics_forcing_to_convection_for_specific_humidity',      &
         UNITS      = 'kg kg-1 s-1',                                               &
         DIMS       =  MAPL_DimsHorzVert,                                          &
         VLOCATION  =  MAPL_VLocationCenter,                                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                    &
         SHORT_NAME = 'DYNF_T',                                                    &
         LONG_NAME  = 'dynamics_forcing_to_convection_for_air_temperature',        &
         UNITS      = 'K s-1',                                                     &
         DIMS       =  MAPL_DimsHorzVert,                                          &
         VLOCATION  =  MAPL_VLocationCenter,                                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                    &
         SHORT_NAME = 'DYNF_PLE',                                                  &
        LONG_NAME  = 'dynamics_forcing_to_convection_for_pressure',               &
         UNITS      = 'Pa s-1',                                                    &
         DIMS       =  MAPL_DimsHorzVert,                                          &
         VLOCATION  =  MAPL_VLocationEdge,                                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                    &
         SHORT_NAME = 'DYNF_UA',                                                    &
         LONG_NAME  = 'dynamics_forcing_to_convection_for_eastward_wind_Agrid',        &
         UNITS      = 'm s-2',                                                     &
         DIMS       =  MAPL_DimsHorzVert,                                          &
         VLOCATION  =  MAPL_VLocationCenter,                                       &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                                    &
         SHORT_NAME = 'DYNF_VA',                                                    &
         LONG_NAME  = 'dynamics_forcing_to_convection_for_norhtward_wind_Agrid',        &
         UNITS      = 'm s-2',                                                     &
         DIMS       =  MAPL_DimsHorzVert,                                          &
         VLOCATION  =  MAPL_VLocationCenter,                                       &
         RC=STATUS  )
    VERIFY_(STATUS)
    
!    call MAPL_AddExportSpec ( gc,                                  &
!         SHORT_NAME = 'T_N',                                       &
!         LONG_NAME  = 'air_temperature at begin of time step',     &
!         UNITS      = 'K',                                         &
!         DIMS       = MAPL_DimsHorzVert,                           &
!         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
!     VERIFY_(STATUS)
     call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'Q_N',                                       &
         LONG_NAME  = 'spec_humidity_at_begin_of_time_step',       &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
!     call MAPL_AddExportSpec ( gc,                                 &
!         SHORT_NAME = 'PLE_N',                                     &
!         LONG_NAME  = 'edge_pressure at begin of time step',       &
!         UNITS      = 'Pa',                                        &
!         DIMS       = MAPL_DimsHorzVert,                           &
!         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
!     VERIFY_(STATUS)
!     call MAPL_AddExportSpec ( gc,                                 &
!         SHORT_NAME = 'U_N',                                       &
!         LONG_NAME  = 'eastward_wind at begin of time step',       &
!         UNITS      = 'm s-1',                                     &
!         DIMS       = MAPL_DimsHorzVert,                           &
!         FIELD_TYPE = MAPL_VectorField,                            &
!         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
!     VERIFY_(STATUS)

!    call MAPL_AddExportSpec ( gc,                                  &
!         SHORT_NAME = 'V_N',                                       &
!         LONG_NAME  = 'northward_wind at begin of time step',      &
!         UNITS      = 'm s-1',                                     &
!         DIMS       = MAPL_DimsHorzVert,                           &
!         FIELD_TYPE = MAPL_VectorField,                            &
!         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
!     VERIFY_(STATUS)
!-srf-gf-scheme
    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'QV_DYN_IN',                                 &
         LONG_NAME  = 'spec_humidity_at_begin_of_time_step',       &
         UNITS      = 'kg kg-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'T_DYN_IN',                                 &
         LONG_NAME  = 'temperature_at_begin_of_time_step',       &
         UNITS      = 'K',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'U_DYN_IN',                                 &
         LONG_NAME  = 'u_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'V_DYN_IN',                                 &
         LONG_NAME  = 'v_wind_at_begin_of_time_step',       &
         UNITS      = 'm s-1',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'PLE_DYN_IN',                                 &
         LONG_NAME  = 'edge_pressure_at_begin_of_time_step',       &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                 &
         SHORT_NAME = 'PKE',                                 &
         LONG_NAME  = 'edge_pressure_with_k',       &
         UNITS      = 'Pa',                                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,             RC=STATUS  )
    VERIFY_(STATUS)

   

!! Exports added for consistency with Superdyn

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME         = 'PEANA',                                             &
       LONG_NAME          = 'total_potential_energy_tendency_due_to_analysis',   &
       UNITS              = 'W m-2',                                             &
       DIMS               = MAPL_DimsHorzOnly,                                   &
       VLOCATION          = MAPL_VLocationNone,                        __RC__ )
    

    call MAPL_AddExportSpec ( gc,                                                &
       SHORT_NAME   = 'PEPHY',                                                   &
       LONG_NAME    = 'total_potential_energy_tendency_due_to_physics',          &
       UNITS        = 'W m-2',                                                   &
       DIMS         = MAPL_DimsHorzOnly,                                         &
       VLOCATION    = MAPL_VLocationNone,                              __RC__ )
    

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DOXDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ozone_tendency_due_to_analysis',               &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQVDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_analysis',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQLDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_liquid_water_tendency_due_to_analysis',        &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQIDTANAINT',                                                        &
         LONG_NAME  = 'vertically_integrated_ice_water_tendency_due_to_analysis',           &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTANAINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_analysis',                 &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DTHVDTPHYINT',                                                       &
         LONG_NAME  = 'vertically_integrated_THV_tendency_due_to_physics',                  &
         UNITS      = 'K kg m-2 s-1',                                                       &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                           &
         SHORT_NAME = 'DQVDTDYNINT',                                                        &
         LONG_NAME  = 'vertically_integrated_water_vapor_tendency_due_to_dynamics',         &
         UNITS      = 'kg m-2 s-1',                                                         &
         DIMS       = MAPL_DimsHorzOnly,                                                    &
         VLOCATION  = MAPL_VLocationNone,                                        __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'PV',                                        &
         LONG_NAME  = 'ertels_isentropic_potential_vorticity',     &
         UNITS      = 'm+2 kg-1 sec-1',                            &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'EPV',                                       &
         LONG_NAME  = 'ertels_potential_vorticity',                &
         UNITS      = 'K m+2 kg-1 sec-1',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                                        &
       SHORT_NAME         = 'TROPP_BLENDED',                                             &
       LONG_NAME          = 'tropopause_pressure_based_on_blended_estimate',             &
       UNITS              = 'Pa',                                                        &
       DIMS               = MAPL_DimsHorzOnly,                                           &
       VLOCATION          = MAPL_VLocationNone,                                __RC__ )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='S',                                      &
         LONG_NAME ='static_energy',                    &
         UNITS     ='J kg-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='TV',                                      &
         LONG_NAME ='air_virtual_temperature',                    &
         UNITS     ='K',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PLK',                                      &
         LONG_NAME ='Exner_quantity',                    &
         UNITS     ='J kg-1',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='PL',                                      &
         LONG_NAME ='midlevel_pressures',                    &
         UNITS     ='Pa',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
                                                        __RC__  )
    

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQVDTDYN',                                      &
         LONG_NAME  = 'tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQLLSDTDYN',                                    &
         LONG_NAME  = 'tendency_of_large_scale_cloud_water_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQILSDTDYN',                                    &
         LONG_NAME  = 'tendency_of_large_scale_cloud_ice_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQLCNDTDYN',                                    &
         LONG_NAME  = 'tendency_of_convective_cloud_water_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DQICNDTDYN',                                    &
         LONG_NAME  = 'tendency_of_convective_cloud_ice_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DCLLSDTDYN',                                    &
         LONG_NAME  = 'tendency_of_large_scale_cloud_fraction_due_to_dynamics', &
         UNITS      = 'fraction',                                      &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'DCLCNDTDYN',                                    &
         LONG_NAME  = 'tendency_of_convective_cloud_fraction_due_to_dynamics', &
         UNITS      = 'fraction',                                      &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'DTDTDYN',                                     &
         LONG_NAME  = 'tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K sec-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             __RC__    )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'HDQDTDYN',                                      &
         LONG_NAME  = 'horiz_tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                      &
         SHORT_NAME = 'VDQDTDYN',                                      &
         LONG_NAME  = 'vertical_tendency_of_specific_humidity_due_to_dynamics', &
         UNITS      = 'kg/kg/sec',                                     &
         DIMS       = MAPL_DimsHorzVert,                               &
         VLOCATION  = MAPL_VLocationCenter,                 __RC__  )
     

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'HDTDTDYN',                                    &
         LONG_NAME  = 'horiz_tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K sec-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             __RC__    )
     

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'HDTHDTDYN',                                   &
         LONG_NAME  = 'horiz_tendency_of_air_pot_temp_due_to_dynamics', &
         UNITS      = 'K sec-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             __RC__    )
     

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'VDTDTDYN',                                    &
         LONG_NAME  = 'vertical_tendency_of_air_temperature_due_to_dynamics', &
         UNITS      = 'K sec-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             __RC__    )
     

    call MAPL_AddExportSpec ( gc,                                    &
         SHORT_NAME = 'VDTHDTDYN',                                   &
         LONG_NAME  = 'vertical_tendency_of_air_pot_temp_due_to_dynamics', &
         UNITS      = 'K sec-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,             __RC__    )
     
    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DXC',                                       &
         LONG_NAME  = 'cgrid_delta_x',                             &
         UNITS      = 'm'  ,                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DYC',                                       &
         LONG_NAME  = 'cgrid_delta_y',                             &
         UNITS      = 'm'  ,                                       &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               RC=STATUS  )
     VERIFY_(STATUS)


    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AREA',                                      &
         LONG_NAME  = 'agrid_cell_area',                           &
         UNITS      = 'm+2'  ,                                     &
         DIMS       = MAPL_DimsHorzOnly,                           &
         VLOCATION  = MAPL_VLocationNone,               __RC__  )
     
    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'AK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_a',                   &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'BK',                                        &
         LONG_NAME  = 'hybrid_sigma_pressure_b',                   &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsVertOnly,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_CGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_C-Grid',                   &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_CGRID',                                   &
         LONG_NAME  = 'northward_wind_on_C-Grid',                  &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'U_DGRID',                                   &
         LONG_NAME  = 'eastward_wind_on_native_D-Grid',            &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'V_DGRID',                                   &
         LONG_NAME  = 'northward_wind_on_native_D-Grid',           &
         UNITS      = 'm s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'PT',                                        &
         LONG_NAME  = 'scaled_potential_temperature',              &
         UNITS      = 'K Pa$^{-\kappa}$',                          &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                &
         SHORT_NAME = 'PE',                                        &
         LONG_NAME  = 'air_pressure',                              &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,               RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                                  &
         SHORT_NAME = 'DELP',                                      &
         LONG_NAME  = 'pressure_thickness',                        &
         UNITS      = 'Pa',                                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
     VERIFY_(STATUS)
!EOS



! Set the Profiling timers
! ------------------------
    call MAPL_TimerAdd(GC, name="DRIVER"    ,__RC__)
    call MAPL_TimerAdd(GC, name="MISC"      ,__RC__)
     
! Set generic init and final methods
! ----------------------------------
    call MAPL_GenericSetServices    ( GC, __RC__)

    RETURN_(ESMF_SUCCESS)
     
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: RUN -- Run method for the SINGLE COLUMN component

! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:
    
! !DESCRIPTION: This version uses the GEOS\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),  pointer  :: STATE
    !type (GEOS_GenericState), pointer   :: STATE
    type (ESMF_Config)                   :: CF
    type (ESMF_State)                    :: INTERNAL
    type (ESMF_Alarm)                    :: ALARM
    type (ESMF_FieldBundle)              :: qbundle
    type (ESMF_Field)                    :: field

    logical                              :: FRIENDLY
    logical ,allocatable                 :: FRIENDLY_arr(:)


! Local variables

    character(len=ESMF_MAXSTR)          :: DATA
    character(len=ESMF_MAXSTR)          :: QNAME
    character(len=ESMF_MAXSTR), allocatable  :: QNAMEarr(:)
    integer                :: SCMmonth, SCMday, SCMhour, SCMminute, SCMyear, SCMsecond
    integer                :: STRTmonth, STRTday, STRThour, STRTminute, STRTyear, STRTsecond
    real, dimension(6)     :: SCMtime
    type (ESMF_Time)       :: currentTime, startTime
    type (ESMF_TimeInterval) :: timeStep

    real    :: DT,Fac0,Fac1,DTXX,RELAX_TO_OBS
    real    :: OROGSGH 

    integer :: IM,JM,LM,L,K,NQ,ii,NOT1,COLDSTART,Ktrc,iip1,itr,ntracs

    real, pointer, dimension(:,:,:) :: PLE,PLEOUT
    real, pointer, dimension(:,:,:) :: ZLE
    real, pointer, dimension(:)     :: PREF,PREF_IN
    real, pointer, dimension(:,:,:) :: Q, QOBS, QHATOBS, QVATOBS
    real, pointer, dimension(:,:,:) :: T, TOBS, THATOBS, TVATOBS,TOUT
    real, pointer, dimension(:,:,:) :: QLOBS, QIOBS
    real, pointer, dimension(:,:,:) :: TH,PL
    real, pointer, dimension(:,:,:) :: U,UE,UOUT
    real, pointer, dimension(:,:,:) :: V,VOUT
    real, pointer, dimension(:,:,:) :: TTPHYS,QTEST,OM,OMOUT,div
    real, pointer, dimension(:,:,:) :: DUDT,DVDT,DTDT,PV,EPV
    real, pointer, dimension(:,:,:) :: DQVDTDYN,DTDTDYN
    real, pointer, dimension(:,:)   :: DQVDTDYNINT
    real, pointer, dimension(:,:,:) :: DQLLSDTDYN,DQILSDTDYN,DQLCNDTDYN,DQICNDTDYN,DCLLSDTDYN,DCLCNDTDYN
    real, pointer, dimension(:,:,:) :: HDQDTDYN,HDTDTDYN,VDQDTDYN,VDTDTDYN
    real, pointer, dimension(:,:,:) :: HDTHDTDYN,VDTHDTDYN

    real, pointer, dimension(:,:)   :: PSFCOBS
    real, pointer, dimension(:,:)   :: PCPOBS
    real, pointer, dimension(:,:)   :: TSAIROBS
    real, pointer, dimension(:,:)   :: TGSOILOBS
    real, pointer, dimension(:,:)   :: PS, SGH, PHIS, PHISOU
    real, pointer, dimension(:,:)   :: DZ
    real, pointer, dimension(:,:)   :: TA
    real, pointer, dimension(:,:)   :: SPEED
    real, pointer, dimension(:,:)   :: QA
    real, pointer, dimension(:,:)   :: TSKINOBS
    real, pointer, dimension(:,:)   :: QSKINOBS
    real, pointer, dimension(:,:)   :: LHOBS
    real, pointer, dimension(:,:)   :: SHOBS
    real, pointer, dimension(:,:)   :: VARFLT
    real, pointer, dimension(:,:)   :: DUMMYAREA
    real, pointer, dimension(:,:)   :: DUMMYDXC,DUMMYDYC
    real, pointer, dimension(:,:,:) :: DUMMYW,DUMMYPLK
    real, pointer, dimension(:,:,:) :: PKEOUT
    real, pointer, dimension(:)     :: AK,BK
    real, pointer, dimension(:,:,:) :: U_CGRID,V_CGRID
    real, pointer, dimension(:,:,:) :: U_DGRID,V_DGRID
    real, pointer, dimension(:,:,:)   :: PT,PE
    real, pointer, dimension(:,:,:)   :: DELTAP

    real, pointer, dimension(:,:,:) :: STATICEN

    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS

    type three_d_ptr
         real, pointer, dimension(:,:,:) :: ptr3
    end type three_d_ptr

    real, allocatable, dimension(:,:,:) :: PKE,ZLO,PLO,QTDYN,TTDYN,OMOBS,THOBS
    real, allocatable, dimension(:,:,:) :: UdQdx,VdQdy,UdTdx,VdTdy, CFMIPRLX,CFMIPRLX1
    type(three_d_ptr), allocatable :: TRCarr(:)
    real, pointer , dimension(:,:,:) :: qdum
            real,allocatable, dimension(:) :: WF, XXX
!=======================================================================

  ! temporary garbage dump for profile data

      real, ALLOCATABLE, SAVE, DIMENSION(:       ) :: time        !Calenday day
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::   yy        !Year
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  mo        !Month
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  dd        !Day
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  hh        !Hour
      real, ALLOCATABLE, SAVE, DIMENSION(:       ) ::  mm        !Minutes
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: LHF
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: SHF 
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: PSFC
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: PCP_OBS 
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: TS_AIR
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: TG_SOIL
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: TSKIN
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: QSKIN
      REAL, ALLOCATABLE, SAVE, DIMENSION(:       ) :: QSFCAIR
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: tt
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: qq , QiQi, QLQL
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: uu
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: vv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: T_H_adv 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: T_V_adv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q_H_adv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q_V_adv
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: T_t_dyn 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q_t_dyn
      
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: T_ana
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q_ana

      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q1 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: Q2 
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: PLE_DATA
      REAL, ALLOCATABLE, SAVE, DIMENSION(: , :   ) :: OMEGA
      REAL, SAVE :: ptop

      real, dimension(4) :: VC 

      INTEGER :: NT, NLEVEL,I,J,VERTADV, useana, advscheme
      real :: zrel,zrelp,qfloor       

      LOGICAL :: USE_ASCII_DATA, AT_START, CFMIP, CFMIP2,CFMIP3
      LOGICAL, SAVE :: ALREADY_HAVE_DATA
      integer, save :: I_time_step,cfcse
      real blendwgt

!=============================================================================

! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    
    Iam = trim(COMP_NAME) // 'Run'

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, __RC__)
    

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( STATE, IM=IM, JM=JM, LM=LM,   &
                               RUNALARM = ALARM,             &
                               CF       = CF,                &
                               LONS     = LONS,              &
                               LATS     = LATS,              &
                               INTERNAL_ESMF_STATE=INTERNAL, &
                                                   __RC__ )
    

#ifdef DEBUG_SCM
           print *," LONS, LATS in SCM " ,lons*180./MAPL_PI ,lats*180./MAPL_PI
#endif

    call ESMF_ConfigGetAttribute ( CF, DTXX,  Label="RUN_DT:", __RC__)
    

    call ESMF_ConfigGetAttribute ( CF, RELAX_TO_OBS,  Label="RELAX_TO_OBS:", &
                                         DEFAULT=0.00,  __RC__)

    call ESMF_ConfigGetAttribute ( CF, VERTADV ,  Label="VERTICAL_ADVECTION:", &
                                         DEFAULT=0,  __RC__)

    call ESMF_ConfigGetAttribute ( CF, USEANA ,  Label="USE_ANALYSIS_INC:", &
                                         DEFAULT=0,  __RC__)

    call ESMF_ConfigGetAttribute ( CF, ADVSCHEME ,  Label="ADVECTION_SCHEME:", &
                                         DEFAULT=3,  __RC__)

    call ESMF_ConfigGetAttribute ( CF, COLDSTART ,  Label="COLDSTART:", &
                                         DEFAULT=1,  __RC__)

    call ESMF_ConfigGetAttribute( cf, value=DATA, label ='DRIVER_DATA:', &
                                    DEFAULT='none', rc = status )

    call ESMF_ConfigGetAttribute( cf, NT, label ='SCM_TIME_LENGTH:', &
                                    DEFAULT=0, rc = status )

    call ESMF_ConfigGetAttribute( cf, NLEVEL, label ='SCM_NLEVEL:', &
                                    DEFAULT=0, rc = status )

    call ESMF_ConfigGetAttribute( cf, USE_ASCII_DATA, label ='SCM_DATA_DRIVER:', &
                                    DEFAULT=.true., rc = status )

    call ESMF_ConfigGetAttribute( cf, CFMIP, label ='SCM_CFMIP:', &
                                    DEFAULT=.false., rc = status )

    call ESMF_ConfigGetAttribute( cf, CFMIP2, label ='SCM_CFMIP2:', &
                                    DEFAULT=.false., rc = status )
                                    
    call ESMF_ConfigGetAttribute( cf, CFMIP3, label ='SCM_CFMIP3:', &
                                    DEFAULT=.false., rc = status )

    call ESMF_ConfigGetAttribute( cf, CFCSE, label ='CGILS_CASE:', &
                                    DEFAULT=0, rc = status )


    call ESMF_ConfigGetAttribute ( CF, OROGSGH,  Label="OROG_STDEV:", &
                                         DEFAULT=100.,  __RC__)

    if ( CFMIP .and. CFMIP2) then
            print *, " Error - SCM_CFMIP and SCM_CFMIP2 cannot be set at the same time  "  ! This should never happen
            RETURN_(ESMF_FAILURE)
    end if

    if ( ( CFMIP .or. CFMIP2) .and. USE_ASCII_DATA ) then
            print *, " Error - USE_ASCII_DATA cannot be set at the same time as  SCM_CFMIP or SCM_CFMIP2  "  ! This should never happen
            RETURN_(ESMF_FAILURE)
    end if

    if ( NT == 0 ) then
            print *, " Error - SCM_NT == 0 -- experiment configuration is wrong (should be set in AGCM.rc ) "  ! This should never happen
            RETURN_(ESMF_FAILURE)
    end if

    if ( NLEVEL == 0 ) then
            print *, " Error - SCM_NLEVEL == 0 -- experiment configuration is wrong (should be set in AGCM.rc )"  ! This should never happen
            RETURN_(ESMF_FAILURE)
    end if

    if ( NT > 1 .and. ( CFMIP .or. CFMIP2 ) ) then
            print *, " Error - SCM_NT > 1 not allowed with CFMIP cases -- experiment configuration is wrong (set in AGCM.rc )"  ! This should never happen
            RETURN_(ESMF_FAILURE)
    end if


    DT  = DTXX

     VC    = 1.0
     vc(3) = 0.0
     if(vertadv.eq.1) vc(2) = 0.
! afe add message here?
     if ( useana .eq. 1) then
      vc(4) = 1.
     else
      vc(4) = 0.
     endif

! Allocate arrays for driver data
!----------------------------------------
 ! Here key off of whether TIME is allocated to determine
 ! whether data has been allocated and read in or not. Data arrays
 ! have "SAVE" attribute, so are allocated once and not deallocated.
 ! Kludgey but seems to work - JTB 7/21/04

     if (.not.(ALLOCATED(TIME)) ) THEN 
        ALLOCATE(time (NT)        , __STAT__ )       
        ALLOCATE(yy (NT)          , __STAT__)        !Year
        ALLOCATE(mo (NT)          , __STAT__)        !Month
        ALLOCATE(dd (NT)          , __STAT__)        !Day
        ALLOCATE(hh (NT)          , __STAT__)        !Hour
        ALLOCATE(mm (NT)          , __STAT__)        !Minutes
        ALLOCATE(PSFC (NT)        , __STAT__)
        ALLOCATE(PCP_OBS (NT)     , __STAT__ ) 
        ALLOCATE(TS_AIR (NT)      , __STAT__ )
        ALLOCATE(TG_SOIL (NT)     , __STAT__ )
        ALLOCATE(TSKIN (NT)       , __STAT__ )
        ALLOCATE(QSKIN (NT)       , __STAT__ )
        ALLOCATE(QSFCAIR(NT)      , __STAT__ )
        ALLOCATE(LHF(NT)          , __STAT__ )
        ALLOCATE(SHF(NT)          , __STAT__ )
        ALLOCATE(tt (NT,LM)       , __STAT__ )
        ALLOCATE(qq (NT,LM)       , __STAT__ )
        ALLOCATE(uu (NT,LM)       , __STAT__ )
        ALLOCATE(vv (NT,LM)       , __STAT__ )
        ALLOCATE(T_H_adv(NT,LM)   , __STAT__ )  
        T_h_adv = 0.
        ALLOCATE(T_V_adv (NT,LM)  , __STAT__)  
        T_v_adv = 0.
        ALLOCATE(Q_H_adv (NT,LM)  , __STAT__)  
        Q_h_adv = 0.
        ALLOCATE(Q_V_adv (NT,LM)  , __STAT__)  
        Q_v_adv = 0.
        ALLOCATE(T_t_dyn (NT,LM)  , __STAT__)  
        T_t_dyn = 0.
        ALLOCATE(Q_t_dyn (NT,LM)  , __STAT__)  
        Q_t_dyn = 0.
        ALLOCATE(T_ana (NT,LM)    , __STAT__)  
        T_ana = 0. 
        ALLOCATE(Q_ana (NT,LM)    , __STAT__)  
        Q_ana = 0.
        ALLOCATE(Q1 (NT,LM)       , __STAT__) 
        ALLOCATE(Q2 (NT,LM)       , __STAT__) 
        ALLOCATE(OMEGA(NT,0:LM)   , __STAT__)
        ALLOCATE(PLE_DATA(NT,0:LM), __STAT__)
        ALLOCATE(qiqi(NT,LM)      , __STAT__)
        ALLOCATE(qlql(NT,LM)      , __STAT__)
        
        I_time_step=0
        ALREADY_HAVE_DATA=.FALSE.
        
      endif

      call MAPL_GetPointer(INTERNAL, PLE, 'PLE', __RC__)
      call MAPL_GetPointer(INTERNAL, PREF_IN,'PREF',__RC__)
      call MAPL_GetPointer(INTERNAL, U,   'U'  , __RC__)
      call MAPL_GetPointer(INTERNAL, V,   'V'  , __RC__)
      call MAPL_GetPointer(INTERNAL, T,   'T'  , __RC__)
      call MAPL_GetPointer(INTERNAL, OM,  'OM' , __RC__)
      call ESMF_StateGet(IMPORT, 'TRADV' , QBUNDLE,   __RC__)

! If its time, get the profile data
! --------------------------------------------

    if ( ESMF_AlarmIsRinging(ALARM, RC=STATUS) ) then
       call ESMF_AlarmRingerOff(ALARM, __RC__)
       call MAPL_TimerOn(STATE,"DRIVER")
       call GET_THE_DATA()
       call MAPL_TimerOff(STATE,"DRIVER")
    endif
    
! We need to get Q and any other "advected" water quantities
! from a bundle
      call ESMF_FieldBundleGet(QBUNDLE, FieldCount=NQ, __RC__)
      
      ALLOCATE( TRCarr(NQ-1), __STAT__ )
      ALLOCATE( FRIENDLY_arr(NQ-1), __STAT__ )
      ALLOCATE( QNAMEarr(NQ-1), __STAT__ )
      Ktrc=1

      do K=1,NQ

         call ESMF_FieldBundleGet    (QBUNDLE, fieldIndex=K, field=FIELD, __RC__)
         call ESMF_FieldGet (FIELD, NAME=QNAME,  __RC__)
         
         ! Get item's friendly status (default is not friendly)
         !-----------------------------------------------------

         call ESMF_AttributeGet  (FIELD, NAME="FriendlyToDYNAMICS",VALUE=FRIENDLY, __RC__)
         if(STATUS /= ESMF_SUCCESS) FRIENDLY = .false.

         if( trim(QNAME).eq.'Q') then
           call ESMF_FieldGet (FIELD, 0, Q, __RC__)
          else
           call ESMF_FieldGet (FIELD, 0, trcarr(ktrc)%ptr3 , __RC__)
            !!TRCarr(Ktrc)%ptr3 = Qdum
           FRIENDLY_arr(Ktrc) = FRIENDLY 
           QNAMEarr(Ktrc)     = QNAME 
           Ktrc=Ktrc+1
         end if
     end do

       do K=1,nq-1 
         if (.not. FRIENDLY_arr(K) ) then
            print *, trim(QNAMEarr(k)), " unfriendly var in dynamics bundle - Error "  ! This should never happen
            RETURN_(ESMF_FAILURE)
         endif
       enddo                    
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Figure out what time it is. Check to see if we 
! are at the beginning of the run
!--------------------------------------------------
      call ESMF_ClockGet(CLOCK, currTime=currentTime, startTime=startTime, timeStep=timeStep, __RC__)
      call ESMF_TimeGet(currentTime, YY=SCMyear, MM= SCMmonth, &
                     DD=SCMday, H=SCMhour, M=SCMminute, &
                     S=SCMsecond, __RC__)
      call ESMF_TimeGet(startTime, YY=STRTyear, MM= STRTmonth, &
                     DD=STRTday, H=STRThour, M=STRTminute, &
                     S=STRTsecond, __RC__)
      !!AT_START   =  ( currentTime == startTime + timeStep)
      AT_START   =  ( I_time_step == 0 )
      AT_START   =  ( AT_START .AND. (COLDSTART == 1) )

#ifdef DEBUG_SCM
      print *," This is the time yyyy mm dd:hh min sec in SCM "
      print *,  SCMyear,SCMmonth,SCMday,SCMhour,SCMminute,SCMsecond
      print *," AT_START: ",AT_START
#endif

      call MAPL_GetPointer(EXPORT, PCPOBS,   'PCPOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, LHOBS,   'LHOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, SHOBS,   'SHOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TSAIROBS,   'TSAIROBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TGSOILOBS,   'TGSOILOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, PSFCOBS,   'PSFCOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TH,  'TH'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, ZLE,  'ZLE'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QTEST,  'QTEST'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, VARFLT,'VARFLT'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, PLEOUT, 'PLE', __RC__)
      call MAPL_GetPointer(EXPORT, UOUT,   'U'  , __RC__)
      call MAPL_GetPointer(EXPORT, VOUT,   'V'  , __RC__)
      call MAPL_GetPointer(EXPORT, TOUT,   'T'  , __RC__)
      call MAPL_GetPointer(EXPORT, OMOUT,  'OMEGA' , __RC__)
      call MAPL_GetPointer(EXPORT, DIV,    'DIV' , ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, PL ,    'PL' , __RC__)
      call MAPL_GetPointer(EXPORT, DTDTDYN,   'DTDTDYN'  , __RC__)
      call MAPL_GetPointer(EXPORT, HDTDTDYN,   'HDTDTDYN'  , __RC__)
      call MAPL_GetPointer(EXPORT, VDTDTDYN,   'VDTDTDYN'  , __RC__)
      call MAPL_GetPointer(EXPORT, HDTHDTDYN,   'HDTHDTDYN'  , __RC__)
      call MAPL_GetPointer(EXPORT, VDTHDTDYN,   'VDTHDTDYN'  , __RC__)
      call MAPL_GetPointer(EXPORT, DQVDTDYN,  'DQVDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DQLLSDTDYN,  'DQLLSDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DQILSDTDYN,  'DQILSDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DQLCNDTDYN,  'DQLCNDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DQICNDTDYN,  'DQICNDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DCLLSDTDYN,  'DCLLSDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DCLCNDTDYN,  'DCLCNDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, HDQDTDYN,  'HDQDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, VDQDTDYN,  'VDQDTDYN' , __RC__)
      call MAPL_GetPointer(EXPORT, DQVDTDYNINT,  'DQVDTDYNINT' , __RC__)
      call MAPL_GetPointer(EXPORT, DUMMYAREA,  'AREA' , __RC__)
      call MAPL_GetPointer(EXPORT, DUMMYW,  'W' , __RC__)
      call MAPL_GetPointer(EXPORT, DUMMYPLK,  'PLK' , __RC__)
      call MAPL_GetPointer(EXPORT, PKEOUT,  'PKE' , __RC__)
      call MAPL_GetPointer(EXPORT, DUMMYDXC,  'DXC' , __RC__)
      call MAPL_GetPointer(EXPORT, DUMMYDYC,  'DYC' , __RC__)
      call MAPL_GetPointer(EXPORT, AK,  'AK' , __RC__)
      call MAPL_GetPointer(EXPORT, BK,  'BK' , __RC__)
      call MAPL_GetPointer(EXPORT, U_DGRID,  'U_DGRID' , __RC__)
      call MAPL_GetPointer(EXPORT, V_DGRID,  'V_DGRID' , __RC__)
      call MAPL_GetPointer(EXPORT, U_CGRID,  'U_CGRID' , __RC__)
      call MAPL_GetPointer(EXPORT, V_CGRID,  'V_CGRID' , __RC__)
      call MAPL_GetPointer(EXPORT, PT,  'PT' , __RC__)
      call MAPL_GetPointer(EXPORT, PE,  'PE' , __RC__)
      call MAPL_GetPointer(EXPORT, DELTAP, 'DELP' , ALLOC=.true., __RC__)

! whenever datmodyn starts using gocart, this will need to be a real value --
! see fvdycore for example
      if(associated(DUMMYAREA)) then
          DUMMYAREA=1.0 
      end if

      if(associated(DUMMYDXC)) then
          DUMMYDXC=1.0 
      end if

      if(associated(DUMMYW)) then
          DUMMYW=1.0 
      end if

      if(associated(DUMMYPLK)) then
          DUMMYPLK=1.0 
      end if

      if(associated(DUMMYDYC)) then
          DUMMYDYC=1.0 
      end if

! added to satisfy desires of da
      if(associated(AK)) then
          AK=0.0 
      end if

     if(associated(BK)) then
          BK=0.0 
      end if

! these could in principle be hooked to U and V (with interpolation),
! but setting to 0.0 to catch unexpectedness

      if(associated(U_CGRID)) then
          U_CGRID=0.0 
      end if

      if(associated(V_CGRID)) then
          V_CGRID=0.0 
      end if

      if(associated(U_DGRID)) then
          U_DGRID=0.0 
      end if

      if(associated(V_DGRID)) then
          V_DGRID=0.0 
      end if

      if(associated(PT)) then
          PT=0.0 
      end if

      if(associated(PE)) PE=0.0 






!JTB: Outputs for surface component
!------------------------------------
      call MAPL_GetPointer(EXPORT, PS,    'PS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TA,    'TA'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, DZ,    'DZ'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, SPEED, 'SPEED' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QA,    'QA'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QSKINOBS,    'QSKINOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TSKINOBS,    'TSKINOBS'    , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TOBS,  'TOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QOBS,  'QOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QLOBS,  'QLOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QIOBS,  'QIOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, THATOBS,  'THATOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, TVATOBS,  'TVATOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QHATOBS,  'QHATOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, QVATOBS,  'QVATOBS' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, UE,  'UE' , &
                           ALLOC=.true., __RC__)
      call MAPL_GetPointer(EXPORT, PHISOU,  'PHIS' , &
                           ALLOC=.true., __RC__)

      ALLOCATE( VdTdy(IM,JM,1:LM), __STAT__ )
      ALLOCATE( VdQdy(IM,JM,1:LM), __STAT__ )
      ALLOCATE( UdTdx(IM,JM,1:LM), __STAT__ )
      ALLOCATE( UdQdx(IM,JM,1:LM), __STAT__ )

      ALLOCATE( CFMIPRLX(IM,JM,1:LM), __STAT__ )
      ALLOCATE( CFMIPRLX1(IM,JM,1:LM), __STAT__ )
     
      ALLOCATE( THOBS(IM,JM,1:LM), __STAT__ )
      ALLOCATE( QTDYN(IM,JM,1:LM), __STAT__ )
      ALLOCATE( TTDYN(IM,JM,1:LM), __STAT__ )
      ALLOCATE( OMOBS(IM,JM,0:LM), __STAT__ )
      ALLOCATE( SGH(IM,JM), __STAT__ )

      SGH  = OROGSGH           
      PHISOU = MAPL_GRAV * SGH

      I_time_step = I_time_step+1

      ASSERT_( .NOT. ((USE_ASCII_DATA .eqv. .FALSE.) .AND. (RELAX_TO_OBS > 0.) ) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If USE_ASCII_DATA=.T. then data exists to drive run, i.e. NT>1
! If USE_ASCII_DATA=.F. then NT=1 and only ICs are provided 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

     if (USE_ASCII_DATA) then
        SCMtime = 1.0* (/SCMyear,SCMmonth,SCMday,SCMhour,SCMminute,SCMsecond /)
        ii   = WhereInTime ( yy,mo,dd,hh,mm, nt, DT, SCMtime, Fac0, Fac1 , __RC__)
        iip1 = ii+1
        
#ifdef DEBUG_SCM
         print *, " Data times and weights"
         print *, yy(ii),mo(ii),dd(ii),hh(ii),mm(ii), Fac0
         print *, SCMtime
         print *, yy(iip1),mo(iip1),dd(iip1),hh(iip1),mm(iip1), Fac1
#endif
     else
        ii   = 1
        iip1 = 1
        Fac0 = 0.5
        Fac1 = 0.5
        SCMtime = 1.0* (/SCMyear,SCMmonth,SCMday,SCMhour,SCMminute,SCMsecond /)
#ifdef DEBUG_SCM
         print *, " NO multi-time Driver Data "
         print *, SCMtime
#endif
     endif

     ALLOCATE( PKE(IM,JM,0:LM), __STAT__ )
     ALLOCATE( PLO(IM,JM,1:LM), __STAT__ )
   ALLOCATE( ZLO(IM,JM,1:LM), __STAT__ )    

 
     do l=0,lm
        PLE(:,:,l)   = ( Fac0*PLE_DATA(ii,l) + Fac1*PLE_DATA(iip1,l) )
     end do

      PS(:,:)     =   PLE(:,:,LM)
      PLO(:,:,1:LM) = ( PLE(:,:,0:LM-1)+PLE(:,:,1:LM) )/2.0
      PKE =  ( PLE / MAPL_P00 )**MAPL_KAPPA 
      DELTAP = ( PLE(:,:,1:LM) - PLE(:,:,0:LM-1) )
      
      do l=1,lm
        QOBS(:,:,l)=( Fac0*QQ(ii,l) + Fac1*QQ(iip1,l) )/1000.
        TOBS(:,:,l)=( Fac0*TT(ii,l) + Fac1*TT(iip1,l) )
        U(:,:,l)   =( Fac0*UU(ii,l) + Fac1*UU(iip1,l) )
        V(:,:,l)   =( Fac0*VV(ii,l) + Fac1*VV(iip1,l) )
        QLOBS(:,:,l)=( Fac0*QLQL(ii,l) + Fac1*QLQL(iip1,l) ) 
        QIOBS(:,:,l)=( Fac0*QiQi(ii,l) + Fac1*QiQi(iip1,l) )
      end do
      THOBS  = TOBS * ( ( MAPL_P00 / PLO )**MAPL_KAPPA )

      PSFCOBS(:,:)   = ( Fac0*PSFC(ii)    + Fac1*PSFC(iip1) )
      PCPOBS(:,:)    = ( Fac0*PCP_OBS(ii) + Fac1*PCP_OBS(iip1) )
      TSAIROBS(:,:)  = ( Fac0*TS_AIR(ii)  + Fac1*TS_AIR(iip1) )
      TGSOILOBS(:,:) = ( Fac0*TG_SOIL(ii) + Fac1*TG_SOIL(iip1) )
      LHOBS(:,:)     = ( Fac0*LHF(ii)   + Fac1*LHF(iip1) ) 
      SHOBS(:,:)     = ( Fac0*SHF(ii)   + Fac1*SHF(iip1) ) 
      QSKINOBS(:,:)  = ( Fac0*QSKIN(ii)   + Fac1*QSKIN(iip1) )/1000.
      TSKINOBS(:,:)  = ( Fac0*TSKIN(ii)   + Fac1*TSKIN(iip1) )
      SPEED(:,:)  =   SQRT( U(:,:,LM)**2  + V(:,:,LM)**2 ) 

      do l=0,lm
        OMOBS(:,:,L) = ( Fac0*OMEGA(ii,l)    + Fac1*OMEGA(iip1,l) )
      end do

       do l=1,lm
        QTDYN(:,:,l) = &
         Fac0*( vc(1)*Q_H_ADV(ii,l)  + vc(2)*Q_V_ADV(ii,l)   + vc(3)*Q_t_DYN(ii,l) + vc(4)*Q_ana(ii,l) ) & 
      +  Fac1*( vc(1)*Q_H_ADV(iip1,l)+ vc(2)*Q_V_ADV(iip1,l) + vc(3)*Q_t_DYN(iip1,l) + vc(4)*Q_ana(iip1,l) ) 
        TTDYN(:,:,l) = &
         Fac0*( vc(1)*T_H_ADV(ii,l)  + vc(2)*T_V_ADV(ii,l)   + vc(3)*T_t_DYN(ii,l) + vc(4)*T_ana(ii,l) ) & 
     +   Fac1*( vc(1)*T_H_ADV(iip1,l)+ vc(2)*T_V_ADV(iip1,l) + vc(3)*T_t_DYN(iip1,l) + vc(4)*T_ana(iip1,l) ) 
      enddo!

      TTDYN = -1.0 * TTDYN 
      QTDYN = -1.0 * QTDYN / 1000.

! diagnostics for tracers that are friendly to dynamics (and advected)
   if ( SIZE(qnamearr) > 0 ) then
    ntracs = SIZE(qnamearr)
    do itr = 1,ntracs
     if(trim(QNAMEarr(itr)).eq.'QLLS'.and.associated(DQLLSDTDYN))DQLLSDTDYN = TRCarr(itr)%ptr3
     if(trim(QNAMEarr(itr)).eq.'QLCN'.and.associated(DQLCNDTDYN))DQLCNDTDYN = TRCarr(itr)%ptr3
     if(trim(QNAMEarr(itr)).eq.'QILS'.and.associated(DQILSDTDYN))DQILSDTDYN = TRCarr(itr)%ptr3
     if(trim(QNAMEarr(itr)).eq.'QICN'.and.associated(DQICNDTDYN))DQICNDTDYN = TRCarr(itr)%ptr3
     if(trim(QNAMEarr(itr)).eq.'CLLS'.and.associated(DCLLSDTDYN))DCLLSDTDYN = TRCarr(itr)%ptr3
     if(trim(QNAMEarr(itr)).eq.'CLCN'.and.associated(DCLCNDTDYN))DCLCNDTDYN = TRCarr(itr)%ptr3
    end do
   end if
     if ( AT_START .or. RELAX_TO_OBS.gt.0) then
! Above top pressures in driving dataset, use model T (blend gradually)
      do L = 1,LM
      do J = 1,JM
      do I = 1,IM
       if(plo(i,j,L).lt.ptop) then
        blendwgt = max( 0., 1. - 0.001*(ptop-plo(i,j,L) )  )
       else
        blendwgt = 1.
       endif
       TOBS(i,j,L) = TOBS(i,j,L) * blendwgt + T(i,j,L) * (1.-blendwgt)
       QOBS(i,j,L) = QOBS(i,j,L) * blendwgt + Q(i,j,L) * (1.-blendwgt)
      enddo
      enddo
      enddo
     endif
      if ( AT_START ) then
        Q  = QOBS
        T  = TOBS
        OM = OMOBS
        QTEST        =0.
        QTEST(:,:,55:60)=1.
      endif     
      TH = T * ( ( MAPL_P00 / PLO )**MAPL_KAPPA )
      OM = OMOBS

      call MAPL_GetPointer(IMPORT, PHIS,  'PHIS' , __RC__)
      ZLE(:,:,LM) = PHIS / MAPL_GRAV
      do l=LM,1,-1
         ZLE(:,:,L-1) = ZLE(:,:,L) + MAPL_CP * TH(:,:,L)*(PKE(:,:,L)-PKE(:,:,L-1))/MAPL_GRAV
      end do
      DZ(:,:)     =  0.5 * ( ZLE(:,:,LM-1) -  ZLE(:,:,LM) )

     ZLO = 0.5*(ZLE(:,:,0:LM-1)+ZLE(:,:,1:LM))

!       print *,'ZLE0',ZLE
!       print *,'ZLO0',ZLO
!       print *,'PKE',PKE
!       print *,'PLE',PLE
!       print *,'DELTAP',DELTAP 
     

        ! These quantities are special. Used in SFC to calculate CDs 
        ! In 3D GEOS5 they are just set == values_at_LM
      QA(:,:)     = (1.0-RELAX_TO_OBS)*Q(:,:,LM) + RELAX_TO_OBS*     &
                                    ( Fac0*QSFCAIR(ii) + Fac1*QSFCAIR(iip1) )/1000.
      TA(:,:)     = (1.0-RELAX_TO_OBS)*T(:,:,LM) + RELAX_TO_OBS*TSAIROBS(:,:)

! Load tendency diagnostics with 'before' values

      if (associated(DTDTDYN))   DTDTDYN  = T
      if (associated(DQVDTDYN))  DQVDTDYN = Q
      if (associated(DQVDTDYNINT))  then
       do k=1,lm
        DQVDTDYNINT = DQVDTDYNINT - Q(:,:,k)*deltap(:,:,k)
       enddo
      endif

      if (associated(HDTDTDYN).and.VERTADV==1)   HDTDTDYN = T
      if (associated(HDTHDTDYN).and.VERTADV==1)   HDTHDTDYN = TH
      if (associated(HDQDTDYN).and.VERTADV==1)   HDQDTDYN = Q

! Increment T and Q with  tendencies from dataset that we are using
      if ( .NOT. AT_START ) then
        Q = Q +  (1.0-RELAX_TO_OBS) * DT * QTDYN +  RELAX_TO_OBS*( QOBS - Q) 
        T = T +  (1.0-RELAX_TO_OBS) * DT * TTDYN +  RELAX_TO_OBS*( TOBS - T)                             
      endif     

      call MAPL_GetPointer(EXPORT  , STATICEN  , 'S'  , __RC__)
      STATICEN  = MAPL_GRAV * ( ZLE(:,:,1:LM) + ZLE(:,:,0:LM-1) ) * 0.5 + MAPL_CP * T

      TH     = T    * ( ( MAPL_P00 / PLO )**MAPL_KAPPA )

      if (associated(HDTHDTDYN).and.VERTADV==1)  HDTHDTDYN = ( TH - HDTHDTDYN ) / DT
      if (associated(HDTDTDYN).and.VERTADV==1)  HDTDTDYN = ( T - HDTDTDYN ) / DT
      if (associated(HDQDTDYN).and.VERTADV==1)  HDQDTDYN = ( Q - HDQDTDYN ) / DT

! Load vertical tendency diagnostics with 'before' values (values after horiz advection)
      if (associated(VDTDTDYN).and.VERTADV==1)   VDTDTDYN = T
      if (associated(VDTHDTDYN).and.VERTADV==1)   VDTHDTDYN = TH
      if (associated(VDQDTDYN).and.VERTADV==1)   VDQDTDYN = Q


!! Do vertical advection by using omega - this routine updates TH and Q

      if (VERTADV==1.and. .not.AT_START) THEN
        call ADVECT_VERTICAL(IM,JM,LM,DIV)
      endif
      T  = TH * ( ( PLO / MAPL_P00 )**MAPL_KAPPA )

      if (associated(DTDTDYN))   DTDTDYN  = ( T - DTDTDYN  ) / DT
      if (associated(DQVDTDYN))  DQVDTDYN = ( Q - DQVDTDYN ) / DT
      if (associated(DQVDTDYNINT))  then
       do k=1,lm
        DQVDTDYNINT = DQVDTDYNINT + Q(:,:,k)*deltap(:,:,k)
       enddo
       DQVDTDYNINT = DQVDTDYNINT / (MAPL_GRAV*DT)
      endif

      if (associated(VDTDTDYN).and.VERTADV==1)  VDTDTDYN = ( T - VDTDTDYN ) / DT
      if (associated(VDTHDTDYN).and.VERTADV==1)  VDTHDTDYN = ( TH - VDTHDTDYN ) / DT
      if (associated(VDQDTDYN).and.VERTADV==1)  VDQDTDYN = ( Q - VDQDTDYN ) / DT
  
! diagnostics for tracers that are friendly to dynamics (and advected)
   if ( SIZE(qnamearr) > 0 ) then
    ntracs = SIZE(qnamearr)
    do itr = 1,ntracs
     if(trim(QNAMEarr(itr)).eq.'QLLS'.and.associated(DQLLSDTDYN))DQLLSDTDYN = (TRCarr(itr)%ptr3-DQLLSDTDYN)/DT
     if(trim(QNAMEarr(itr)).eq.'QLCN'.and.associated(DQLCNDTDYN))DQLCNDTDYN = (TRCarr(itr)%ptr3-DQLCNDTDYN)/DT
     if(trim(QNAMEarr(itr)).eq.'QILS'.and.associated(DQILSDTDYN))DQILSDTDYN = (TRCarr(itr)%ptr3-DQILSDTDYN)/DT
     if(trim(QNAMEarr(itr)).eq.'QICN'.and.associated(DQICNDTDYN))DQICNDTDYN = (TRCarr(itr)%ptr3-DQICNDTDYN)/DT
     if(trim(QNAMEarr(itr)).eq.'CLLS'.and.associated(DCLLSDTDYN))DCLLSDTDYN = (TRCarr(itr)%ptr3-DCLLSDTDYN)/DT
     if(trim(QNAMEarr(itr)).eq.'CLCN'.and.associated(DCLCNDTDYN))DCLCNDTDYN = (TRCarr(itr)%ptr3-DCLCNDTDYN)/DT
    end do
   end if

      CFMIPRLX = 0.00

!    print *,'CFMIP3',CFMIP3
!    Forcing based on Phase 2 of CGILS intercomparison. See Blossey et al. (2016)

      if ( CFMIP3 ) then  

         ZLO = 0.5*(ZLE(:,:,0:LM-1)+ZLE(:,:,1:LM))

!        print *,'cfmip3'         
!        print *,'cfcse',cfcse
         if (CFCSE .eq. 12) then
           zrel=1200.
           zrelp=1500.
!           qfloor=0.003561
           qfloor=0.003581839
         elseif (CFCSE .eq. 11) then
           zrel=2500.
           zrelp=3000.
           qfloor=3.55e-3
         elseif (CFCSE .eq. 6) then
           zrel=4000.
           zrelp=4800.
           qfloor=0.
         else
           print *,'error - define the right case'
           RETURN_(ESMF_FAILURE)
         endif
         

!        print *,'CFCSE',CFCSE
!       print *,'ZLO',ZLO      
!       print *,'ZLE',ZLE
          where((ZLO>zrel).and.(ZLO<zrelp))
              CFMIPRLX=1./(2.*10800.)*(1.-cos(3.14*(ZLO-zrel)/(zrelp-zrel)))
          endwhere

          where(ZLO>=zrelp) 
              CFMIPRLX=1./10800.
          endwhere
          
          
           CFMIPRLX1=0. 
           where((Q<qfloor) .and. (ZLO<1300.))
           CFMIPRLX1=1./3600.
           endwhere


!        print *,'Q',Q
!        print *,'CFMIPRLX1',CFMIPRLX1
!        print *,'CFMIPRLX',CFMIPRLX

         T = T - CFMIPRLX * ( T - TOBS ) * DT 
         Q = Q - CFMIPRLX * ( Q - QOBS ) * DT-CFMIPRLX1*(Q-qfloor)*DT  
         if (associated(DTDTDYN))   DTDTDYN  = DTDTDYN  - CFMIPRLX * ( T - TOBS )
         if (associated(DQVDTDYN))  DQVDTDYN = DQVDTDYN - CFMIPRLX * ( Q - QOBS )-CFMIPRLX1*(Q-qfloor)*DT
!      print *,'modified T and Q'   
   end if  ! CFMIP3 switch





!       if ( NT ==  1 .and. .not. CFMIP2 ) then  ! cgils case
       if ( CFMIP ) then  ! cgils case

         where(( PLO < 40000. ).and.( PLO >= 20000. ))
             CFMIPRLX = 86400. - (86400.-14400.)*(40000. - PLO)/20000.
         endwhere
         where( PLO < 20000. )
             CFMIPRLX = 14400. 
         endwhere
         where( CFMIPRLX > 0.000 )
             CFMIPRLX = 1./ CFMIPRLX 
         endwhere
         T = T - CFMIPRLX * ( T - TOBS ) * DT 
         Q = Q - CFMIPRLX * ( Q - QOBS ) * DT  
         if (associated(DTDTDYN))   DTDTDYN  = DTDTDYN  - CFMIPRLX * ( T - TOBS )
         if (associated(DQVDTDYN))  DQVDTDYN = DQVDTDYN - CFMIPRLX * ( Q - QOBS )
      end if

       if ( CFMIP2 ) then
         where(( PLO < 60000. ).and.( PLO >= 40000. ))
             CFMIPRLX = 86400. - (86400.-14400.)*(60000. - PLO)/20000.
         endwhere
         where( PLO < 40000. )
             CFMIPRLX = 14400. 
         endwhere
         where( CFMIPRLX > 0.000 )
             CFMIPRLX = 1./ CFMIPRLX 
         endwhere
         T = T - CFMIPRLX * ( T - TOBS ) * DT 
         Q = Q - CFMIPRLX * ( Q - QOBS ) * DT  
         if (associated(DTDTDYN))   DTDTDYN  = DTDTDYN  - CFMIPRLX * ( T - TOBS )
         if (associated(DQVDTDYN))  DQVDTDYN = DQVDTDYN - CFMIPRLX * ( Q - QOBS )
      end if

      call MAPL_GetPointer(EXPORT, PREF,   'PREF'    , &
                           ALLOC=.true., __RC__)
      

      PREF = PREF_IN

      VARFLT = 0.
      if (associated(PL    )) PL     = PLO
      if (associated(PLEOUT)) PLEOUT = PLE
      if (associated(TOUT))   TOUT   = T
      if (associated(VOUT))   VOUT   = V
      if (associated(UOUT))   UOUT   = U
      if (associated(OMOUT))  OMOUT  = OM
      if (associated(PKEOUT))  PKEOUT  = PKE

      DEALLOCATE(THOBS)
      DEALLOCATE(PKE)
      DEALLOCATE(PLO)
      DEALLOCATE(TTDYN)
      DEALLOCATE(QTDYN)
      DEALLOCATE(CFMIPRLX)

      DEALLOCATE(OMOBS)
      DEALLOCATE(SGH)
 
      DEALLOCATE( VdTdy )
      DEALLOCATE( VdQdy )
      DEALLOCATE( UdTdx )
      DEALLOCATE( UdQdx )

      DEALLOCATE( TRCarr )
      DEALLOCATE( FRIENDLY_arr )
      DEALLOCATE( QNAMEarr )

   RETURN_(ESMF_SUCCESS)

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ADVECT_VERTICAL(IM,JM,LM , DIV )
      implicit none
      integer, intent(in):: IM,JM,LM
!
! advscheme is a flag to indicate which vertical advection scheme to use
!            1 - third order upwind
!            2 - first order upwind
!            3 - first order upwind - new routine
! (future)   4 - use remapping scheme like fv core
!
      real, dimension(IM,JM,1:LM), intent(out)   :: DIV

      real, dimension(IM,JM,1:LM+1) :: WF
      real, dimension(LM) :: XXX,DP,PP
      integer :: ntracs,itr

      do l=0,LM
       do j=1,JM
        do i=1,IM
         WF(I,J,L+1)  = OM(I,J,L)
        end do
       end do
      end do
 
      DIV(:,:,1:LM)  = ( OM(:,:,1:LM) - OM(:,:,0:LM-1) ) / ( PLE(:,:,1:LM) - PLE(:,:,0:LM-1) )

      do j=1,JM
       do i=1,IM
          DP    =  PLE(i,j,1:LM)   - PLE(i,j,0:LM-1) 
          PP    =( PLE(i,j,1:LM)   + PLE(i,j,0:LM-1) ) / 2.0 
       end do
      end do

      do j=1,JM
       do i=1,IM
          XXX   = QTEST(i,j,1:LM)
         if(advscheme.eq.1) then
          call upstad3( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.2) then
          call upstad1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.3) then
          call upvert1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         else
          print *,' INVALID VERTICAL ADVECTION SCHEME - DID NOT DO ADVECTION '
          RETURN_(ESMF_FAILURE)
         endif
          QTEST(i,j,1:LM) = XXX
       end do
      end do

      do j=1,JM
       do i=1,IM
          XXX   = Q(i,j,1:LM)
         if(advscheme.eq.1) then
          call upstad3( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.2) then
          call upstad1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.3) then
          call upvert1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         else
          print *,' INVALID VERTICAL ADVECTION SCHEME - DID NOT DO ADVECTION '
          RETURN_(ESMF_FAILURE)
         endif
          Q(i,j,1:LM) = XXX
       end do
      end do

      do j=1,JM
       do i=1,IM
          XXX   = TH(i,j,1:LM)
         if(advscheme.eq.1) then
          call upstad3( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.2) then
          call upstad1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.3) then
          call upvert1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         else
          print *,' INVALID VERTICAL ADVECTION SCHEME - DID NOT DO ADVECTION '
          RETURN_(ESMF_FAILURE)
         endif
          TH(i,j,1:LM) = XXX
       end do
      end do

      if ( SIZE(qnamearr) > 0 ) then 
        ntracs = SIZE(qnamearr)
        do itr = 1,ntracs
#ifdef DEBUG_SCM
            print *, trim(QNAMEarr(itr)), " being advected !!!!! "
#endif
         do j=1,JM
         do i=1,IM
                 XXX   = TRCarr(itr)%ptr3(i,j,1:LM)
         if(advscheme.eq.1) then
          call upstad3( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.2) then
          call upstad1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         elseif (advscheme.eq.3) then
          call upvert1( WF(i,j,1:LM+1) , XXX, PP, DT, 1 , LM  )
         else
            print *,' INVALID VERTICAL ADVECTION SCHEME - DID NOT DO ADVECTION '
            RETURN_(ESMF_FAILURE)
         endif
         TRCarr(itr)%ptr3(i,j,1:LM) = XXX
         end do
         end do
        end do
      end if

  end subroutine ADVECT_VERTICAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine GET_THE_DATA

     
     if (.not.(ALREADY_HAVE_DATA)) THEN
  
      QLQL = 0.
      QIQI = 0.     

       if ( CFMIP .or. CFMIP2 .or. CFMIP3 ) then
           CALL CFMIP_IC(trim(DATA)//'.bin', NT, NLEVEL,     &
                     LM ,                                &
                     PREF_IN,                            &
                     time,                               &
                     yy,                                 &
                     mo,                                 &
                     dd,                                 &
                     hh,                                 &
                     mm,                                 &
                     PCP_OBS,                            &
                     TS_AIR,                             &
                     TG_SOIL,                            &
                     TSKIN,                              &
                     QSFCAIR,                            &
                     QSKIN,                              &
                     PSFC,                               &
		     LHF,                                &
		     SHF,                                &
                     ptop,                               &
                     tt,                                 &
                     qq,                                 &
                     uu,                                 &
                     vv,                                 &
                     omega,                              &
                     T_H_adv,                            &
                     T_V_adv,                            &
                     Q_H_adv,                            &
                     Q_V_adv,                            &
                     PLE_DATA                             ) 


                     LHF = 0.
                     SHF = 0.

         else if (USE_ASCII_DATA) then      
                CALL ARM2(trim(DATA)//'.dat', NT, NLEVEL,     &
                     LM ,                                &
                     PREF_IN,                            &
                     time,                               &
                     yy,                                 &
                     mo,                                 &
                     dd,                                 &
                     hh,                                 &
                     mm,                                 &
                     PCP_OBS,                            &
                     TS_AIR,                             &
                     TG_SOIL,                            &
                     TSKIN,                              &
                     QSFCAIR,                            &
                     QSKIN,                              &
                     PSFC,                               &
		     LHF,                                &
		     SHF,                                &
		     ptop,                               &
                     tt,                                 &
                     qq,                                 &
                     uu,                                 &
                     vv,                                 &
                     omega,                              &
                     T_H_adv,                            &
                     T_V_adv,                            &
                     Q_H_adv,                            &
                     Q_V_adv,                            &
                     useana,                             &
                     T_ana,                              &
                     Q_ana,                              &
                     Q1,                                 &
                     Q2,                                 & 
                     PLE_DATA                             ) 

     
         else 
         write(*,*) " No good data "
         RETURN_(ESMF_FAILURE)
         
        end if
      ALREADY_HAVE_DATA=.TRUE.

      else 
#ifdef DEBUG_SCM
        write(*,*) " ALREADY HAVE DATA "
#endif
      endif

    end subroutine GET_THE_DATA

  end subroutine RUN

! !IROUTINE: RunAddIncs Run method 2 for Data Dynamics. Adds in physics tendencies

! !INTERFACE:

   subroutine RunAddIncs ( GC, IMPORT, EXPORT, CLOCK, RC )
   !!subroutine RunDummy ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code

! Local derived type aliases

    type (MAPL_MetaComp),  pointer       :: STATE
    type (ESMF_Config)                   :: CF
    type (ESMF_State)                    :: INTERNAL

! !DESCRIPTION: The run method to add increments
!

!EOP

! ErrLog Variables

   character(len=ESMF_MAXSTR)          :: IAm
   integer                             :: STATUS
   character(len=ESMF_MAXSTR)          :: COMP_NAME


  real, pointer, dimension(:,:,:) :: T, TTPHYS,PLE,ZLE,S
  integer :: IM,JM,LM

  real    :: RELAX_TO_OBS,DTXX,DT
  real, allocatable, dimension(:,:,:) :: DELTAP,ZLO

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "RunAddIncs"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, __RC__ )
    
    Iam = trim(COMP_NAME) // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, __RC__)
 
! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get( STATE, IM=IM, JM=JM, LM=LM,   &
                   CF       = CF,                &
                   INTERNAL_ESMF_STATE=INTERNAL, &
                    __RC__ )
    
    call ESMF_ConfigGetAttribute ( CF, DTXX,  Label="RUN_DT:", __RC__)
    call ESMF_ConfigGetAttribute ( CF, RELAX_TO_OBS,  Label="RELAX_TO_OBS:", &
                                   DEFAULT=0.00,  __RC__)
    DT=DTXX

    call MAPL_GetPointer(INTERNAL, PLE, 'PLE', __RC__)
    call MAPL_GetPointer(INTERNAL, T,   'T'  , __RC__)
    call MAPL_GetPointer(IMPORT  , TTPHYS, 'DTDT', __RC__)
    call MAPL_GetPointer(EXPORT  , ZLE, 'ZLE', __RC__)
    call MAPL_GetPointer(EXPORT  , S  , 'S'  , __RC__)

    ALLOCATE( DELTAP(IM,JM,1:LM), __STAT__ )
    ALLOCATE( ZLO(IM,JM,1:LM), __STAT__ )

    DELTAP = ( PLE(:,:,1:LM) - PLE(:,:,0:LM-1) )
    ZLO    = ( ZLE(:,:,1:LM) + ZLE(:,:,0:LM-1) ) * 0.5

    ! This needs to be protected against having
    ! RELAX_TO_OBS=1 when USE_ASCII_DATA=.F.

    T = T +  (1.0-RELAX_TO_OBS) * DT * TTPHYS / DELTAP 

    S  = MAPL_GRAV * ZLO + MAPL_CP * T

    RETURN_(ESMF_SUCCESS)

  end subroutine RunAddIncs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function WhereInTime ( yy,mo,dd,hh,mm, nt, DT, SCMtime,Fac0,Fac1, RC ) result(ii)
integer,                intent(in) :: nt
real,    dimension(nt), intent(in) :: yy,mo,dd,hh,mm
real,    dimension(6),  intent(in) :: SCMtime
real,                   intent(in) :: DT
real,                   intent(out):: Fac0,Fac1
integer, optional,      intent(out):: RC
integer                            :: ii

integer :: i

real    ::  SCM_SecOfDay,Da1_SecOfDay,Da0_SecOfDay,hh_i,hh_i1
real    ::  SCM_SecN,Da1_SecN,Da0_SecN

ii=-1
Fac0=0.
Fac1=0.
RC=0

SCM_SecN = SecOfYear( SCMtime(1), SCMtime(2), SCMtime(3), SCMtime(4), SCMtime(5), SCMtime(6) )

if ( nt == 1) then
   i=1
   Da0_SecN = SecOfYear( yy(i),   mo(i),   dd(i),   hh(i)  , mm(i)  , 0. ) 
   if (ABS(Da0_SecN-SCM_SecN) < 60.) then
      ii=1
   else
      ii=-99
   endif
   RETURN
end if

do i=1,nt-1

   if ( yy(i) .ne. SCMtime(1) ) cycle

   Da0_SecN = SecOfYear( yy(i),   mo(i),   dd(i),   hh(i)  , mm(i)  , 0. ) 
   Da1_SecN = SecOfYear( yy(i+1), mo(i+1), dd(i+1), hh(i+1), mm(i+1), 0. ) 

   if ( yy(i+1) .gt. yy(i) ) then
         Da1_SecN = Da1_SecN + SecOfYear( yy(i), 12., 31., 23., 59., 59. )   
   endif

   if ( ( Da0_SecN .gt. SCM_SecN ) .or. ( Da1_SecN .le. SCM_SecN ) ) cycle

         !! now we should be in the correct Data interval
         !! Lets stop here
       ii=i
  
       Fac0 = ( Da1_SecN - SCM_SecN ) / ( Da1_SecN - Da0_SecN ) 
       Fac1 = 1.0-Fac0
  
end do

IF (II.lt.1) THEN
  RC=-1
  print *, ' BAD Time in SCM ',SCMtime
  print *, ' ii = ',ii
  print *," This data starts on "
  print *, yy(1),mo(1),dd(1),hh(1),mm(1)
  print *," This data ends on "
  print *, yy(nt),mo(nt),dd(nt),hh(nt),mm(nt)
  print *, Da0_SecN, SCM_SecN, Da1_SecN
  print *," this year, that year "
  print *, yy(i),yy(i+1)
  print *, Fac0,Fac1,nt

  !do i=1,nt
  !print *, yy(i),mo(i),dd(i),hh(i),mm(i)
  !enddo 

ENDIF

end function WhereInTime 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function SecOfYear(yr,mo,dy,hh,mm,ss) result(SecN)
real, intent(in) :: yr,mo,dy,hh,mm,ss
real             :: SecN 
integer, dimension(12)      :: monlen
integer, dimension(0:12)    :: DaysSoFar
integer :: i, DayN, iyr, imo

imo=INT(mo)
iyr=INT(yr)

if (MOD(iyr,4).ne.0) then ! yes, this happens to work for us, don't overthink it
   !!        J   F   M   A   M   J   J   A   S   O   N   D
   monlen=(/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
ELSE
   monlen=(/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
endif

DaysSoFar(0)=0
do i=1,12
   DaysSoFar(i)=DaysSoFar(i-1)+monlen(i)
end do

DayN=DaysSoFar(imo-1)+dy-1

SecN = DayN*24.*3600. + hh*3600. + mm*60. + ss

end function SecOfYear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GEOS_DatmoDynGridCompMod

