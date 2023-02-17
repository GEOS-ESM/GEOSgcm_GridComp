
! VERIFY_ and RETURN_ macros for error handling.

!#define UWDIAG 1
!#define PDFDIAG 1

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_Moist -- A Module to compute moist processes, including convection,
!   large-scale condensation and precipitation and cloud parameters.

! !INTERFACE:

module GEOS_MoistGridCompMod

  ! !USES:

  use ESMF
  use MAPL
  use GEOS_GFDL_1M_InterfaceMod
  use GEOS_BACM_1M_InterfaceMod
  use GEOS_MGB2_2M_InterfaceMod
  use GEOS_RAS_InterfaceMod
  use GEOS_GF_InterfaceMod
  use GEOS_UW_InterfaceMod

  use aer_cloud
  use Aer_Actv_Single_Moment
  use Lightning_mod, only: HEMCO_FlashRate
  use GEOSmoist_Process_Library
  use GEOS_UtilsMod

  implicit none

  character(LEN=ESMF_MAXSTR):: CONVPAR_OPTION  ! GF, RAS, NONE
  character(LEN=ESMF_MAXSTR):: SHALLOW_OPTION  ! UW, NONE
  character(LEN=ESMF_MAXSTR):: CLDMICR_OPTION  ! BACM_1M, GFDL_1M, MGB2_2M

  private

  logical :: DEBUG = .false.
  logical :: LDIAGNOSE_PRECIP_TYPE
  logical :: LUPDATE_PRECIP_TYPE
  logical :: USE_AERO_BUFFER
  real    :: CCN_OCN
  real    :: CCN_LND
  
  ! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  ! !DESCRIPTION:
  ! 
  !   {\tt GEOS\_MoistGridCompMod} implements moist processes in GEOS-5. These
  !   include all processes that involve phase changes in the atmosphere, such
  !   as large-scale condensation, convective clouds, and all rain and cloud
  !   formation. Its state consists of water vapor, various types of condensate,
  !   and fractions of various cloud types. 
  !   two moment cloud microphysics (Barahona et al., GMD, 2014.) can be run by setting CLDMACRO==MGB2_2M. 
  !   When using 2-moment microphysics the number concentration of ice crystals and cloud droplets 
  !   are part of the state.
  !

  !EOP

contains

  !BOP

  ! !IROUTINE: SetServices -- Sets ESMF services for this component

  ! !INTERFACE:

  subroutine SetServices ( GC, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

    ! !DESCRIPTION:  {\tt GEOS\_MoistGridCompMod} uses the default Initialize and Finalize 
    !                services, but registers its own Run method.

    !EOP

    !=============================================================================
    !
    ! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    ! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: STATE 
    type (ESMF_Config      )            :: CF

    integer      :: RFRSHINT
    integer      :: AVRGNINT
    real         :: DT
  
    logical :: LCONVPAR
    logical :: LSHALLOW
    logical :: LCLDMICR
 
    !=============================================================================

    ! Begin...

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! Set the Init entry point
    !-----------------------------------------
    call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_INITIALIZE,  Initialize, RC=status )
    VERIFY_(STATUS)

    ! Set the Run entry point
    ! -----------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=status )
    VERIFY_(STATUS)

    ! Get the configuration from the component
    !-----------------------------------------
    call ESMF_GridCompGet( GC, CONFIG = CF, RC=STATUS )
    VERIFY_(STATUS)

    ! Set the state variable specs.
    ! -----------------------------
    call ESMF_ConfigGetAttribute( CF, DT,       Label="RUN_DT:",                              RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute( CF, RFRSHINT, Label="REFRESH_INTERVAL:",  default=nint(DT), RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute( CF, AVRGNINT, Label='AVERAGING_INTERVAL:',default=RFRSHINT, RC=STATUS)
    VERIFY_(STATUS)

    ! Inititialize deep convective parameterizations (Options: RAS, GF or NONE)
    !----------------------------------------------------------------------
    call ESMF_ConfigGetAttribute( CF, CONVPAR_OPTION, Label='CONVPAR_OPTION:', default="GF", RC=STATUS)
    VERIFY_(STATUS)
    LCONVPAR = adjustl(CONVPAR_OPTION)=="RAS" .or. &
               adjustl(CONVPAR_OPTION)=="GF" .or. &
               adjustl(CONVPAR_OPTION)=="NONE"
    _ASSERT( LCONVPAR, 'Unsupported Deep Convection Option' )

    ! Inititialize shallow convective parameterizations (Options: UW or NONE)
    !----------------------------------------------------------------------
    call ESMF_ConfigGetAttribute( CF, SHALLOW_OPTION, Label="SHALLOW_OPTION:",  default="UW", RC=STATUS)
    VERIFY_(STATUS)
    LSHALLOW = adjustl(SHALLOW_OPTION)=="UW" .or. &
               adjustl(SHALLOW_OPTION)=="NONE"
    _ASSERT( LSHALLOW, 'Unsupported Shallow Convection Option' )

    ! Inititialize cloud microphysics (Options: BACM_1M, MGB2_2M or GFDL_1M)
    !--------------------------------------------------------------
    call ESMF_ConfigGetAttribute( CF, CLDMICR_OPTION, Label="CLDMICR_OPTION:",  default="GFDL_1M", RC=STATUS)
    VERIFY_(STATUS)
    LCLDMICR = adjustl(CLDMICR_OPTION)=="BACM_1M" .or. &
               adjustl(CLDMICR_OPTION)=="MGB2_2M" .or. &
               adjustl(CLDMICR_OPTION)=="GFDL_1M"
    _ASSERT( LCLDMICR, 'Unsupported Cloud Microphysics Option' )


    if (adjustl(CONVPAR_OPTION)=="RAS"    ) call     RAS_Setup(GC, CF, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CONVPAR_OPTION)=="GF"     ) call      GF_Setup(GC, CF, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(SHALLOW_OPTION)=="UW"     ) call      UW_Setup(GC, CF, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CLDMICR_OPTION)=="BACM_1M") call BACM_1M_Setup(GC, CF, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CLDMICR_OPTION)=="MGB2_2M") call MGB2_2M_Setup(GC, CF, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CLDMICR_OPTION)=="GFDL_1M") call GFDL_1M_Setup(GC, CF, RC=STATUS) ; VERIFY_(STATUS)

    ! !IMPORT STATE:

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
         SHORT_NAME = 'PREF',                                       &
         LONG_NAME  = 'reference_air_pressure',                     &
         UNITS      = 'Pa',                                         &
         DIMS       = MAPL_DimsVertOnly,                            &
         VLOCATION  = MAPL_VLocationEdge,                           &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME = 'ZLE',                                       &
         LONG_NAME  = 'geopotential_height',                       &
         UNITS      = 'm',                                         &
         DIMS       =  MAPL_DimsHorzVert,                          &
         VLOCATION  =  MAPL_VLocationEdge,                         &
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

    call MAPL_AddImportSpec(GC,                                              &
         SHORT_NAME = 'WQT',                                                   &
         LONG_NAME  = 'Total_water_flux',                                      &
         UNITS      = 'kg kg-1 m s-1',                                               &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationCenter,                                    &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                              &
         SHORT_NAME = 'WHL',                                                   &
         LONG_NAME  = 'Liquid_water_static_energy_flux',                       &
         UNITS      = 'K m s-1',                                               &
         DIMS       = MAPL_DimsHorzVert,                                       &
         VLOCATION  = MAPL_VLocationCenter,                                    &
         AVERAGING_INTERVAL = AVRGNINT,                             &
         REFRESH_INTERVAL   = RFRSHINT,                             &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'W2',                                       &
            LONG_NAME  = 'variance_of_vertical_velocity',             &
            UNITS      = 'm2 s-2',                                       &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'W3',                                       &
            LONG_NAME  = 'third_moment_of_vertical_velocity',         &
            UNITS      = 'm3 s-3',                                     &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'HL3',                                       &
            LONG_NAME  = 'third_moment_of_liquid_water_static_energy',    &
            UNITS      = 'K+3',                                       &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'EDMF_FRC',                                       &
            LONG_NAME  = 'Mass_Flux_Fractional_Area',                 &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'HL2',                                       &
            LONG_NAME  = 'variance_of_liquid_water_static_energy',    &
            UNITS      = 'K+2',                                       &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'QT2',                                       &
            LONG_NAME  = 'variance_of_total_water_specific_humidity', &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'QT3',                                       &
            LONG_NAME  = 'third_moment_of_total_water_specific_humidity', &
            UNITS      = '1',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

       call MAPL_AddImportSpec(GC,                                  &
            SHORT_NAME = 'HLQT',                                      &
            LONG_NAME  = 'covariance_of_liquid_water_static_energy_and_total_water_specific_humidity', &
            UNITS      = 'K',                                         &
            DIMS       = MAPL_DimsHorzVert,                           &
            VLOCATION  = MAPL_VLocationCenter,                          &
            AVERAGING_INTERVAL = AVRGNINT,                            &
            REFRESH_INTERVAL   = RFRSHINT,                            &
            RC=STATUS )
       VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                             &
         SHORT_NAME = 'T',                                         &
         LONG_NAME  = 'temperature',                               &
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

    call MAPL_AddImportSpec(GC,                                     &
         LONG_NAME  = 'aerosols',                                   &
         UNITS      = '1',                                          &
         SHORT_NAME = 'AERO',                                   &
         DIMS       = MAPL_DimsHorzVert,                            &
         VLOCATION  = MAPL_VLocationCenter,                         &
         DATATYPE   = MAPL_StateItem,                               &
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
       LONG_NAME  = 'HYSTPDF_iterations',                                    &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDFITERS',                                              &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_relative_area_fraction',                       &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_A',                                                 &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

#ifdef PDFDIAG
    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_vertical_velocity_standard_deviation_first_plume', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_SIGW1',                                             &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_vertical_velocity_standard_deviation_second_plume', &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_SIGW2',                                             &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_vertical_velocity_of_first_plume',         &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'PDF_W1',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_vertical_velocity_of_second_plume',        &
       UNITS      = 'm s-1',                                                 &
       SHORT_NAME = 'PDF_W2',                                                &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_liq_wat_pot_temp_of_first_plume',       &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_SIGTH1',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_liq_wat_pot_temp_of_second_plume',      &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_SIGTH2',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_liq_wat_pot_temp_of_first_plume',          &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_TH1',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_liq_wat_pot_temp_of_second_plume',         &
       UNITS      = 'K',                                                     &
       SHORT_NAME = 'PDF_TH2',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_total_water_of_first_plume',            &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_SIGQT1',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_stddev_total_water_of_second_plume',           &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_SIGQT2',                                            &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_total_water_of_first_plume',               &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_QT1',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_avg_total_water_of_second_plume',              &
       UNITS      = 'kg kg-1',                                               &
       SHORT_NAME = 'PDF_QT2',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_corr_total_water_liq_wat_pot_temp',            &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_RQTTH',                                               &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_corr_vertical_velocity_liq_wat_pot_temp',      &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_RWTH',                                              &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                              &
       LONG_NAME  = 'SHOC_PDF_corr_vertical_velocity_total_water',           &
       UNITS      = '1',                                                     &
       SHORT_NAME = 'PDF_RWQT',                                              &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)
#endif

    call MAPL_AddExportSpec(GC,                                              &
       SHORT_NAME = 'WTHV2',                                                 &
       LONG_NAME  = 'Buoyancy_flux_for_SHOC',                                &
       UNITS      = '1',                                                     &
       DIMS       = MAPL_DimsHorzVert,                                       &
       VLOCATION  = MAPL_VLocationCenter,                                    &
                                                                  RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                              &
       SHORT_NAME = 'WQL',                                                   &
       LONG_NAME  = 'Liquid_water_flux',                                     &
       UNITS      = 'kg kg-1 m s-1',                                         &
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
         SHORT_NAME = 'QPTOTLS',                                    &
         LONG_NAME  = 'mass_fraction_of_large_scale_falling_precip', & 
         UNITS      = 'kg kg-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                              &
         SHORT_NAME = 'DTDT',                                     &
         LONG_NAME = 'pressure_weighted_temperature_tendency_due_to_moist',&
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
         SHORT_NAME = 'DQADT ',                                      &
         LONG_NAME = 'total_cloud_tendency_due_to_moist',       &
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
         SHORT_NAME='DQRDT',                                         &
         LONG_NAME ='QRAIN tendency due to moist ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQSDT',                                         &
         LONG_NAME ='QSNOW tendency due to moist ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQGDT',                                         &
         LONG_NAME ='QGRAUPEL tendency due to moist ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'SHL_DQCDT',                                  &
         LONG_NAME = 'shallow_cu_condensate_source',                &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,                         &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME= 'CNV_DQCDT',                                  &
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

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CLDBASEHGT',                                &
         LONG_NAME = 'Height_of_cloud_base',                       &
         UNITS     = 'm',                                          &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                           &
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

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'QTFLX_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_total_water_flux',          &
         UNITS      = 'kg kg-1 m s-1',                             &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'SLFLX_SC',                                  &
         LONG_NAME  = 'shallow_cumulus_liquid_static_energy_flux', &
         UNITS      = 'K m s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'UFLX_SC',                                   &
         LONG_NAME  = 'shallow_cumulus_u_wind_flux',               &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME = 'VFLX_SC',                                   &
         LONG_NAME  = 'shallow_cumulus_v_wind_flux',               &
         UNITS      = '1',                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationEdge,                          &
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
         SHORT_NAME = 'DQADT_SC',                             &
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
         SHORT_NAME = 'DTDT_SC',                                  &
         LONG_NAME  = 'Temperature_tendency_from_shallow_convection',            &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQIDT_SC',                                  &
         LONG_NAME  = 'Ice_tendency_from_shallow_convection',      &
         UNITS      = 'kg kg-1 s-1',                               &
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
         SHORT_NAME = 'DCM_SC',                                      &
         LONG_NAME = 'Shallow_convection_detrained_cloud_mass', &
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
         SHORT_NAME = 'CUSH_SC',                                      &
         LONG_NAME = 'cumulus_scale_height_for_shallow_convection',&
         UNITS     = 'm',                                          &
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
         SHORT_NAME ='NCCN_LIQ',                                     &
         LONG_NAME ='number_concentration_of_cloud_liquid_particles',     &
         UNITS     ='cm-3',                                         &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='NCCN_ICE',                                     &
         LONG_NAME ='number_concentration_of_ice_cloud_particles',     &
         UNITS     ='cm-3',                                         &
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
         LONG_NAME ='deep_convective_precipitation',               &
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
         SHORT_NAME ='CNPCPRATE',                                  &
         LONG_NAME ='convective_precipitation_from_GF',            &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='LS_SNR',                                     &
         LONG_NAME ='large_scale_snow',                            &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'AN_SNR',                                     &
         LONG_NAME = 'anvil_snow',                         &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='CN_SNR',                                     &
         LONG_NAME ='deep_convective_snow',                        &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='SC_SNR',                                    &
         LONG_NAME ='shallow_convective_snow',            &
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
         SHORT_NAME = 'DTDT_ER',                                  &
         LONG_NAME  = 'temperature_tendency_from_RH_cleanup',            &
         UNITS      = 'K s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'DQVDT_ER',                                  &
         LONG_NAME  = 'specific_humidty_tendency_from_RH_cleanup',            &
         UNITS      = 'kg m-2 s-1',                                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='SC_MSE',                                    & 
         LONG_NAME ='shallow_convective_column_MSE_tendency',      &
         UNITS     ='W m-2',                                       &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME ='SC_QT',                                       & 
         LONG_NAME ='shallow_convective_column_QT_tendency',      &
         UNITS     ='kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'FILLNQV_IN',                                     & 
         LONG_NAME = 'filling_of_negative_Q_on_entry_to_moist',          &
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
         SHORT_NAME = 'PREC_STRAT',                                &
         LONG_NAME = 'all_stratiform_precipitation',               &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PREC_CONV',                                 &
         LONG_NAME = 'all_convective_precipitation',               &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'RAIN',                                         &
         LONG_NAME = 'rainfall',                                    &
         UNITS     = 'kg m-2 s-1',                                  &
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
         SHORT_NAME = 'SNOWTOTAL',                                &
         LONG_NAME = 'snowfall_total',                       &
         UNITS     = 'mm',                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PRECTOTAL',                                &
         LONG_NAME = 'precipitation_total',                        &     
         UNITS     = 'mm',                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'KUCHERA_RATIO',                       &
         LONG_NAME = 'kuchera_snow_to_liquid_ratio',     &
         UNITS     = 'unitless',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PCU',                                         & 
         LONG_NAME = 'liquid_convective_precipitation',            &
         UNITS     = 'kg m-2 s-1',                                  &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'PLS',                                         & 
         LONG_NAME = 'liquid_large_scale_precipitation',           &
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
         SHORT_NAME='CLWP',                                        & 
         LONG_NAME ='cloud_liquid_water_path',                     &
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
         SHORT_NAME='QCBL',                                      &
         LONG_NAME ='q_at_cloud_base_level',               &
         UNITS     ='kg/kg'  ,                                         &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='SRF_TYPE',                                      &
         LONG_NAME ='surface type for ice_fraction',               &
         UNITS     =''  ,                                         &
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
         SHORT_NAME='RHCRIT',                                         &
         LONG_NAME ='critical_relative_humidity_for_PDF',               &
         UNITS     ='1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
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
         SHORT_NAME = 'UMST0',                                    &
         LONG_NAME  = 'zonal_wind_before_moist_processes',       &
         UNITS      = 'm s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( gc,                            &
         SHORT_NAME = 'VMST0',                                    &
         LONG_NAME  = 'meridonal_wind_before_moist_processes',       &
         UNITS      = 'm s-1',                                    &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,              RC=STATUS  )
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

    call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME='LFR_GCC',                                        &
         LONG_NAME ='lightning_flash_rate_for_GEOSCHEMchem',          &
         UNITS     ='km-2 s-1',                                       &
         DIMS      = MAPL_DimsHorzOnly,                               &
         VLOCATION = MAPL_VLocationNone,                              &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='UAFMOIST',                                     &
         LONG_NAME ='zonal_wind_after_all_of_moist',   &
         UNITS     ='m s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='VAFMOIST',                                     &
         LONG_NAME ='meridional_wind_after_all_of_moist',   &
         UNITS     ='m s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='THAFMOIST',                                     & 
         LONG_NAME ='potential_temperature_after_all_of_moist',   &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='TAFMOIST',                                     &
         LONG_NAME ='temperature_after_all_of_moist',   &
         UNITS     ='K',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='QAFMOIST',                                     &
         LONG_NAME ='specific_humidity_after_all_of_moist',   &
         UNITS     ='kg kg-1',                                     &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME='SAFMOIST',                                      & 
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
    
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='NWFA',                                      &
         LONG_NAME ='Number concentration of water-friendly aerosol',               &
         UNITS     ='Kg-1'  ,                                         &
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
         SHORT_NAME='TPERT_SC',                                          & 
         LONG_NAME ='Shallow_convection_source_air_temperature_perturbation', &
         UNITS     ='K',                                             &
         DIMS      = MAPL_DimsHorzOnly,                                  &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                          &
         SHORT_NAME='QPERT_SC',                                          & 
         LONG_NAME ='Shallow_convection_source_air_humidity_perturbation', &
         UNITS     ='kg kg-1',                                             &
         DIMS      = MAPL_DimsHorzOnly,                                  &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
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
         SHORT_NAME='DUDT_macro',                                         &
         LONG_NAME ='U tendency due to macrophysics ',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DVDT_macro',                                         &
         LONG_NAME ='V tendency due to macrophysics ',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
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
         SHORT_NAME='DQRDT_macro',                                         &
         LONG_NAME ='QRAIN tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQSDT_macro',                                         &
         LONG_NAME ='QSNOW tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQGDT_macro',                                         &
         LONG_NAME ='QGRAUPEL tendency due to macrophysics ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DUDT_DC',                                         &
         LONG_NAME ='U-wind tendency due to deep convection',               &
         UNITS     ='m s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DVDT_DC',                                         &
         LONG_NAME ='V-wind tendency due to deep convection',               &
         UNITS     ='K s-2',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DTDT_DC',                                         &
         LONG_NAME ='T tendency due to deep convection',               &
         UNITS     ='K s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQVDT_DC',                                         &
         LONG_NAME ='QV tendency due to deep convection ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQLDT_DC',                                         &
         LONG_NAME ='QL tendency due to deep convection ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQIDT_DC',                                         &
         LONG_NAME ='QI tendency due to deep convection ',               &
         UNITS     ='kg kg-1 s-1',                                           &
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME='DQADT_DC',                                         &
         LONG_NAME ='QA tendency due to deep convection ',               &
         UNITS     ='kg kg-1 s-1',                                           &
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
    
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DQDT_GF',                                    &
         LONG_NAME  = 'tendency_of_spec_humidity_due_GF',        &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DTDT_GF',                                    &
         LONG_NAME  = 'tendency_of_temp_due_GF',        &
         UNITS      = 'K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DUDT_GF',                                    &
         LONG_NAME  = 'tendency_of_zonal_wind_due_GF',        &
         UNITS      = 'kg kg-1 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'DVDT_GF',                                    &
         LONG_NAME  = 'tendency_of_meridional_wind_due_GF',        &
         UNITS      = 'K s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MUPDP',                                    &
         LONG_NAME  = 'Mass_flux_deep_GF_updraft',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MDNDP',                                    &
         LONG_NAME  = 'Mass_flux_deep_GF_downdraft',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MUPSH',                                    &
         LONG_NAME  = 'Mass_flux_shallow_GF',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'MUPMD',                                    &
         LONG_NAME  = 'Mass_flux_congestus_GF',        &
         UNITS      = 'kg m-2 s-1',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'VAR3d_a',                                    &
         LONG_NAME  = 'dummy_array_for_output_GF',        &
         UNITS      = 'unknown',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'VAR3d_b',                                  &
         LONG_NAME  = 'dummy_array_for_output_GF',                &
         UNITS      = 'unknown',                                  &
         DIMS       = MAPL_DimsHorzVert,                          &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'VAR3d_c',                                    &
         LONG_NAME  = 'dummy_array_for_output_GF',        &
         UNITS      = 'unknown',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                     &
         SHORT_NAME = 'VAR3d_d',                                    &
         LONG_NAME  = 'dummy_array_for_output_GF',        &
         UNITS      = 'unknown',                                  &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
       VERIFY_(STATUS)
       !-2d
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'MFDP',                              &
        LONG_NAME          = 'mass_flux_cloud_base_deep_GF',      &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'MFSH',                              &
        LONG_NAME          = 'mass_flux_cloud_base_shal_GF',      &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'MFMD',                              &
        LONG_NAME          = 'mass_flux_cloud_base_mid_GF',       &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'ERRDP',                             &
        LONG_NAME          = 'convection_code_deep_GF',                   & 
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'ERRSH',                             &
        LONG_NAME          = 'convection_code_shallow_GF',                   & 
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'ERRMD',                             &
        LONG_NAME          = 'convection_code_mid_GF',            & 
        UNITS              = '1',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA0',                               &
        LONG_NAME          = 'cloud work function 0',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA1',                               &
        LONG_NAME          = 'cloud work function 1',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA2',                               &
        LONG_NAME          = 'cloud work function 2',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA3',                               &
        LONG_NAME          = 'cloud work function 3',             & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA1_CIN',                           &
        LONG_NAME          = 'cloud work function CIN',           & 
        UNITS              = 'J kg-1',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'AA1_BL',                            &
        LONG_NAME          = 'Bound layer AA1',                   & 
        UNITS              = 'J kg-1 s-1',                            &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'TAU_BL',                             &
        LONG_NAME          = 'Bound layer time scale',            & 
        UNITS              = 's',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)
        
        call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'TAU_EC',                             &
        LONG_NAME          = 'cape removal time scale',           & 
        UNITS              = 's',                                &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
        VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='TPWI',                                        & 
         LONG_NAME ='initial_total_precipitable_water',            &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
         VERIFY_(STATUS)
        
        call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='TPWI_star',                                   & 
         LONG_NAME ='saturation_initial_total_precipitable_water', &
         UNITS     ='kg m-2'  ,                                    &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
         VERIFY_(STATUS)

        call MAPL_AddExportSpec(GC,                                &
         SHORT_NAME='LFR_GF',                                 & 
         LONG_NAME ='lightning_flash_density ',                    &
         UNITS     ='km-2 day-1'  ,                                &
         DIMS      = MAPL_DimsHorzOnly,                            & 
         VLOCATION = MAPL_VLocationNone,                RC=STATUS  )
        VERIFY_(STATUS)

    !EOS

    ! Set the Profiling timers
    ! ------------------------

    call MAPL_TimerAdd(GC,name="---CONV_TRACERS"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---AERO_ACTIVATE"   ,RC=STATUS)
    VERIFY_(STATUS)    
    call MAPL_TimerAdd(GC,name="---CLDMACRO"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="---CLDMACRO"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,name="--FLASH"    ,RC=STATUS)
    VERIFY_(STATUS)


    ! Set generic init and final methods
    ! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!-!-!-!!!!!!!!!!!!!!!

  !EOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !BOP

  ! !IROUTINE: Initialize -- Initialize method for the composite Moist Gridded Component

  ! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    ! !DESCRIPTION: The Initialize method of the Moist Physics Gridded Component first 
    !   calls the Initialize method of the child Dynamics.  The Dynamics Initialize method will
    !   create the ESMF GRID, which will then be used to set the GRID associated with the
    !   SuperDyn Composite Component itself.  It should be noted that the 
    !   SuperDyn Initialize method also invokes the GEOS Topo Utility which creates all
    !   topography related quantities.

    !EOP


    ! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

    ! Local derived type aliases

    type (MAPL_MetaComp),      pointer  :: MAPL
    type (ESMF_Grid )                   :: GRID
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL

    type (ESMF_Config)                  :: CF

    !=============================================================================

    ! Begin... 

    ! Get the target components name and set-up traceback handle.
    ! -----------------------------------------------------------
    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, config=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // trim(Iam)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS )
    VERIFY_(STATUS)

    ! Call Generic Initialize for MOIST GC
    !----------------------------------------

    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK, RC=STATUS )
    VERIFY_(STATUS)

    ! Get parameters from generic state.
    !-----------------------------------
    call MAPL_GetResource( MAPL, LDIAGNOSE_PRECIP_TYPE, Label="DIAGNOSE_PRECIP_TYPE:",  default=.FALSE.,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LUPDATE_PRECIP_TYPE,   Label="UPDATE_PRECIP_TYPE:",    default=.FALSE., RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, USE_AEROSOL_NN  , 'USE_AEROSOL_NN:'  , DEFAULT=.TRUE.        , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, USE_BERGERON    , 'USE_BERGERON:'    , DEFAULT=USE_AEROSOL_NN, RC=STATUS); VERIFY_(STATUS)
    if (USE_AEROSOL_NN) then
      call MAPL_GetResource( MAPL, USE_AERO_BUFFER , 'USE_AERO_BUFFER:' , DEFAULT=.TRUE. , RC=STATUS); VERIFY_(STATUS)
      call aer_cloud_init()
      call WRITE_PARALLEL ("INITIALIZED aer_cloud_init")
    else
      call MAPL_GetResource( MAPL, CCN_OCN, 'NCCN_OCN:', DEFAULT= 300., RC=STATUS); VERIFY_(STATUS) ! #/cm^3
      call MAPL_GetResource( MAPL, CCN_LND, 'NCCN_LND:', DEFAULT= 100., RC=STATUS); VERIFY_(STATUS) ! #/cm^3
    endif

    if (adjustl(CONVPAR_OPTION)=="RAS"    ) call     RAS_Initialize(MAPL, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CONVPAR_OPTION)=="GF"     ) call      GF_Initialize(MAPL, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(SHALLOW_OPTION)=="UW"     ) call      UW_Initialize(MAPL, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CLDMICR_OPTION)=="BACM_1M") call BACM_1M_Initialize(MAPL, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CLDMICR_OPTION)=="GFDL_1M") call GFDL_1M_Initialize(MAPL, RC=STATUS) ; VERIFY_(STATUS)
    if (adjustl(CLDMICR_OPTION)=="MGB2_2M") call MGB2_2M_Initialize(MAPL, RC=STATUS) ; VERIFY_(STATUS)

    ! All done
    !---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


  !===================================================================================

  !BOP

  ! !IROUTINE: RUN -- Run method for the CONVECT component

  ! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

    ! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
    !                the Initialize and Finalize services, as well as allocating

    !EOP


    ! ErrLog Variables

    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type( ESMF_VM )                 :: VMG

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: DT_MOIST

    ! Local variables
    real                                :: Tmax
    real, allocatable, dimension(:,:,:) :: PLEmb, PKE, ZLE0, PK, MASS
    real, allocatable, dimension(:,:,:) :: PLmb,  ZL0, DZET
    real, allocatable, dimension(:,:,:) :: QST3, DQST3, MWFA
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D
    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W
    real, pointer, dimension(:,:)   :: FRLAND, FRLANDICE, FRACI, SNOMAS
    real, pointer, dimension(:,:)   :: SH, EVAP, KPBL
    real, pointer, dimension(:,:,:) :: OMEGA
    type(ESMF_State)                :: AERO
    type(ESMF_FieldBundle)          :: TR
    ! Exports
    real, pointer, dimension(:,:,:) :: DQDT, DQADT, DQIDT, DQLDT, DQRDT, DQSDT, DQGDT 
    real, pointer, dimension(:,:,:) :: DTDT, DUDT,  DVDT,  DWDT
    real, pointer, dimension(:,:,:) :: DPDTMST, PFL_LSAN, PFI_LSAN
    real, pointer, dimension(:,:,:) :: DTDT_ER, DQVDT_ER
    real, pointer, dimension(:,:  ) :: PTYPE, TPREC, CN_PRCP, LS_PRCP, AN_PRCP, SC_PRCP, PLS, PCU
    real, pointer, dimension(:,:  ) :: RAIN, SNOW, ICE, FRZR, PREC_STRAT, PREC_CONV
    real, pointer, dimension(:,:,:) :: BYNCY
    real, pointer, dimension(:,:  ) :: CAPE, INHB
    real, pointer, dimension(:,:  ) :: CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:) :: CFICE, CFLIQ
    real, pointer, dimension(:,:,:  ) :: NWFA
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D
    

    
    real :: CNV_NUMLIQ_SC, CNV_NUMICE_SC
    

    integer :: IM,JM,LM
    integer :: I, J, L

    !=============================================================================

    ! Begin... 

    ! Get my name and set-up traceback handle
    ! ---------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, VM=VMG, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS) ; VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"TOTAL")

    ! If its time, call run methods
    ! --------------------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    if ( ESMF_AlarmIsRinging( ALARM, RC=STATUS) ) then

       call ESMF_AlarmRingerOff(ALARM, RC=STATUS) ; VERIFY_(STATUS)

       ! Internal State
       call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, NACTL,   'NACTL'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(INTERNAL, NACTI,   'NACTI'    , RC=STATUS); VERIFY_(STATUS)

       ! Import State
       call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, W,       'W'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, KPBL,    'KPBL'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, SH,      'SH'      , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, EVAP,    'EVAP'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   , RC=STATUS); VERIFY_(STATUS)
       call   ESMF_StateGet(IMPORT,'AERO',    AERO     , RC=STATUS); VERIFY_(STATUS)
       call   ESMF_StateGet(IMPORT,'MTR',     TR       , RC=STATUS); VERIFY_(STATUS)

       ! Update SRF_TYPE for ice_fraction
       call MAPL_GetPointer(IMPORT, FRLAND,    'FRLAND'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FRLANDICE, 'FRLANDICE' , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FRACI,     'FRACI'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, SNOMAS,    'SNOMAS'    , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SRF_TYPE,  'SRF_TYPE'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       SRF_TYPE = 0.0 ! Ocean
       where (FRLAND > 0.1)
         SRF_TYPE = 1.0 ! Land
       end where
       where ( (SNOMAS > 0.1) .OR. (FRLANDICE > 0.5) .OR. (FRACI > 0.5) )
         SRF_TYPE = 2.0 ! Ice/Snow
       end where

       ! Allocatables
        ! Edge variables 
       ALLOCATE ( ZLE0 (IM,JM,0:LM) )
       ALLOCATE ( PLEmb(IM,JM,0:LM) )
       ALLOCATE ( PKE  (IM,JM,0:LM) )
        ! Layer variables
       ALLOCATE ( ZL0  (IM,JM,LM  ) )
       ALLOCATE ( DZET (IM,JM,LM  ) )
       ALLOCATE ( PLmb (IM,JM,LM  ) )
       ALLOCATE ( PK   (IM,JM,LM  ) )
       ALLOCATE ( DQST3(IM,JM,LM  ) )
       ALLOCATE (  QST3(IM,JM,LM  ) )
       ALLOCATE ( MASS (IM,JM,LM  ) )
       ALLOCATE ( TMP3D(IM,JM,LM  ) )  
       ALLOCATE ( TMP2D(IM,JM     ) )

       ! Save input winds
       call MAPL_GetPointer(EXPORT, PTR3D, 'UMST0', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       PTR3D = U
       call MAPL_GetPointer(EXPORT, PTR3D, 'VMST0', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       PTR3D = V

       ! Derived States
       MASS     = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )/MAPL_GRAV
       call FILLQ2ZERO(Q, MASS, TMP2D)
       call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV_IN', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = TMP2D
       PLEmb    =  PLE*.01
       PKE      = (PLE/MAPL_P00)**(MAPL_KAPPA)
       PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
       PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
       DO L=0,LM
          ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
       END DO
       ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
       DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
       DQST3    = GEOS_DQSAT(T, PLmb, QSAT=QST3)

       ! These may be used by children
       call MAPL_GetPointer(EXPORT, NWFA, 'NWFA', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)        
       call MAPL_GetPointer(EXPORT, CNV_FRC, 'CNV_FRC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, BYNCY,   'BYNCY'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, CAPE,    'CAPE'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, INHB,    'INHB'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call BUOYANCY( T, Q, QST3, DQST3, DZET, ZL0, BYNCY, CAPE, INHB)
       CNV_FRC = 0.0
       if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
         WHERE (CAPE .ne. MAPL_UNDEF)
            CNV_FRC =(MAX(1.e-6,MIN(1.0,(CAPE-CNV_FRACTION_MIN)/(CNV_FRACTION_MAX-CNV_FRACTION_MIN))))
         END WHERE
       endif
       if (CNV_FRACTION_EXP /= 1.0) then
          CNV_FRC = CNV_FRC**CNV_FRACTION_EXP
       endif

       ! Extract convective tracers from the TR bundle
       call MAPL_TimerOn (MAPL,"---CONV_TRACERS")
       call CNV_Tracers_Init(TR, RC)
       call MAPL_TimerOff(MAPL,"---CONV_TRACERS")

       ! Get aerosol activation properties
       call MAPL_TimerOn (MAPL,"---AERO_ACTIVATE")
       if (USE_AEROSOL_NN) then
         allocate ( AeroProps(IM,JM,LM) )
         ! Pressures in Pa
         call Aer_Activation(IM,JM,LM, Q, T, PLmb*100.0, PLE, ZL0, ZLE0, QLCN, QICN, QLLS, QILS, &
                             SH, EVAP, KPBL, OMEGA, FRLAND, USE_AERO_BUFFER, &
                             AeroProps, AERO, NACTL, NACTI, NWFA)
       else
         do L=1,LM
           NACTL(:,:,L) = (CCN_LND*FRLAND + CCN_OCN*(1.0-FRLAND))*1.e6 ! #/m^3
           NACTI(:,:,L) = (CCN_LND*FRLAND + CCN_OCN*(1.0-FRLAND))*1.e6 ! #/m^3
         end do
       endif
       call MAPL_GetPointer(EXPORT, PTR3D, 'NCCN_LIQ', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = NACTL*1.e-6
       call MAPL_GetPointer(EXPORT, PTR3D, 'NCCN_ICE', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = NACTI*1.e-6

       call MAPL_TimerOff(MAPL,"---AERO_ACTIVATE")

       if (adjustl(CONVPAR_OPTION)=="RAS"    ) call     RAS_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CONVPAR_OPTION)=="GF"     ) call      GF_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(SHALLOW_OPTION)=="UW"     ) call      UW_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CLDMICR_OPTION)=="BACM_1M") call BACM_1M_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CLDMICR_OPTION)=="GFDL_1M") call GFDL_1M_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)
       if (adjustl(CLDMICR_OPTION)=="MGB2_2M") call MGB2_2M_Run(GC, IMPORT, EXPORT, CLOCK, RC=STATUS) ; VERIFY_(STATUS)

       ! Exports
         ! Cloud fraction exports
         call MAPL_GetPointer(EXPORT, CFICE, 'CFICE', ALLOC=.true., RC=STATUS); VERIFY_(STATUS)
         if (associated(CFICE)) then
           CFICE=0.0
           WHERE (QILS+QICN .gt. 1.0e-12)
              CFICE=(CLLS+CLCN)*(QILS+QICN)/(QLLS+QLCN+QILS+QICN)
           END WHERE
           CFICE=MAX(MIN(CFICE, 1.0), 0.0)
         endif
         call MAPL_GetPointer(EXPORT, CFLIQ, 'CFLIQ', RC=STATUS); VERIFY_(STATUS)
         if (associated(CFLIQ)) then
           CFLIQ=0.0
           WHERE (QLLS+QLCN .gt. 1.0e-12)
              CFLIQ=(CLLS+CLCN)*(QLLS+QLCN)/(QLLS+QLCN+QILS+QICN)
           END WHERE
           CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         endif
         ! Rain-out and Relative Humidity where RH > 110%
         call MAPL_GetPointer(EXPORT,  DTDT_ER,  'DTDT_ER', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DQVDT_ER, 'DQVDT_ER', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
          DTDT_ER = T
         DQVDT_ER = Q
         ! some diagnostics to export
         QST3=GEOS_QsatICE (T, PLmb*100.0)
         call MAPL_GetPointer(EXPORT, PTR3D, 'SAT_RAT', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           where (CFICE .lt. 0.99 .and. QST3 .gt. 1.0e-20)
            TMP3D = max((Q - QST3*CFICE), 0.0)/(1.0-CFICE)
            PTR3D = min(TMP3D/QST3, 2.0)
           elsewhere
            PTR3D = 1.0
           end where
         endif
         call MAPL_GetPointer(EXPORT, PTR3D, 'RHICE', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           PTR3D = Q/QST3
           where (T>MAPL_TICE)
             PTR3D=0.0
           end where
         endif
     !!! QST3  = GEOS_QsatLQU (T, PLmb*100.0, DQ=DQST3) !clean up only with respect to liquid water
         DQST3 = GEOS_DQSAT   (T, PLmb, QSAT=QST3)      ! this qsat function expects hPa...
         call MAPL_GetPointer(EXPORT, PTR3D, 'RHLIQ', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) PTR3D = Q/QST3
         where ( Q > 1.1*QST3 )
            TMP3D = (Q - 1.1*QST3)/( 1.0 + 1.1*DQST3*MAPL_ALHL/MAPL_CP )
         elsewhere
            TMP3D = 0.0
         endwhere
         call MAPL_GetPointer(EXPORT, LS_PRCP, 'LS_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, PTR2D,   'ER_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         PTR2D = SUM(TMP3D*MASS,3)/DT_MOIST
         LS_PRCP = LS_PRCP + PTR2D
         Q = Q - TMP3D
         T = T + (MAPL_ALHL/MAPL_CP)*TMP3D
          DTDT_ER = (T -  DTDT_ER)/DT_MOIST
         DQVDT_ER = (Q - DQVDT_ER)/DT_MOIST
         ! cleanup any negative QV/QC/CF
         call FILLQ2ZERO(Q, MASS, TMP2D)
         call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR2D)) PTR2D = TMP2D/DT_MOIST

       if (USE_AEROSOL_NN) then
         deallocate ( AeroProps )
       endif
       

       ! Export Total Moist Tendencies

       call MAPL_GetPointer(EXPORT, DUDT, 'DUDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DUDT)) then
          DUDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DUDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DUDT = DUDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DVDT, 'DVDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DVDT)) then
          DVDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DVDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DVDT = DVDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DTDT, 'DTDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DTDT)) then
          DTDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DTDT_ER'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DTDT = DTDT + PTR3D
          DTDT = DTDT*(PLE(:,:,1:LM)-PLE(:,:,0:LM-1)) ! Pressure weighted tendency
       endif

       call MAPL_GetPointer(EXPORT, DQDT, 'DQDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQDT)) then
          DQDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQVDT_ER'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQDT = DQDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQLDT, 'DQLDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQLDT)) then
          DQLDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQLDT = DQLDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQIDT, 'DQIDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQIDT)) then
          DQIDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQIDT = DQIDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQRDT, 'DQRDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQRDT)) then
          DQRDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQRDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQRDT = DQRDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQSDT, 'DQSDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQSDT)) then
          DQSDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQSDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQSDT = DQSDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQGDT, 'DQGDT', RC=STATUS); VERIFY_(STATUS)
       if (associated(DQGDT)) then
          DQGDT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQGDT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQGDT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQGDT = DQGDT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DQADT, 'DQADT'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(DQADT)) then
          DQADT = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_DC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_macro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'DQADT_micro', RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DQADT = DQADT + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, DPDTMST, 'DPDTMST'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(DPDTMST)) then
          DPDTMST = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_CN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_CN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_SC'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) DPDTMST = DPDTMST + PTR3D(:,:,0:LM-1)-PTR3D(:,:,1:LM)
          DPDTMST = MAPL_GRAV * DPDTMST
        endif

       call MAPL_GetPointer(EXPORT, PFL_LSAN, 'PFL_LSAN'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(PFL_LSAN)) then
          PFL_LSAN = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFL_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFL_LSAN = PFL_LSAN + PTR3D
       endif

       call MAPL_GetPointer(EXPORT, PFI_LSAN, 'PFI_LSAN'  , RC=STATUS); VERIFY_(STATUS)
       if (associated(PFI_LSAN)) then
          PFI_LSAN = 0.0
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_AN'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
          call MAPL_GetPointer(EXPORT, PTR3D, 'PFI_LS'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR3D)) PFI_LSAN = PFI_LSAN + PTR3D
       endif

       ! Combine Precip Exports

       ! liquid convective precip
       call MAPL_GetPointer(EXPORT, PCU, 'PCU', RC=STATUS); VERIFY_(STATUS)
       if (associated(PCU)) then
          PCU = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PCU = PCU + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PCU = PCU + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE' , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PCU = PCU + PTR2D
          PCU = MAX(PCU, 0.0)
       endif

       ! liquid large-scale precip
       call MAPL_GetPointer(EXPORT, PLS, 'PLS', RC=STATUS); VERIFY_(STATUS)
       if (associated(PLS)) then
          PLS = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PLS = PLS + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) PLS = PLS + PTR2D
          PLS = MAX(PLS, 0.0)
       endif

       ! all liquid precip
       call MAPL_GetPointer(EXPORT, RAIN, 'RAIN', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(RAIN)) then
          RAIN = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE' , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) RAIN = RAIN + PTR2D
          RAIN = MAX(RAIN, 0.0)
       endif

       ! all frozen precip (snow at this point)
       call MAPL_GetPointer(EXPORT, SNOW, 'SNO', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(SNOW)) then
          SNOW = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_SNR'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SNOW = SNOW + PTR2D
          SNOW = MAX(SNOW, 0.0)
       endif

       ! all deep convective precip (rain+snow)
       call MAPL_GetPointer(EXPORT, CN_PRCP, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(CN_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'CNPCPRATE' , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) CN_PRCP = CN_PRCP + PTR2D
       endif

       ! all large-scale precip (rain+snow)
       call MAPL_GetPointer(EXPORT, LS_PRCP, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(LS_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) LS_PRCP = LS_PRCP + PTR2D
       endif

       ! all anvil precip (rain+snow)
       call MAPL_GetPointer(EXPORT, AN_PRCP, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(AN_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) AN_PRCP = AN_PRCP + PTR2D
       endif

       ! all shallow precip (rain+snow)
       call MAPL_GetPointer(EXPORT, SC_PRCP, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
       if (associated(SC_PRCP)) then
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_SNR'    , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) SC_PRCP = SC_PRCP + PTR2D
       endif

       ! Total - all precip (rain+snow)
       call MAPL_GetPointer(EXPORT, TPREC, 'TPREC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(TPREC)) then
          TPREC = 0.0
          call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
          if (associated(PTR2D)) TPREC = TPREC + PTR2D
          TPREC = MAX(TPREC, 0.0)
       endif

       ! diagnosed stratiform precip (rain+snow)
       call MAPL_GetPointer(EXPORT, PREC_STRAT, 'PREC_STRAT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PREC_STRAT)) then
          if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
             PREC_STRAT = (1.0-CNV_FRC)*TPREC
          else
             PREC_STRAT = 0.0
             call MAPL_GetPointer(EXPORT, PTR2D, 'LS_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_STRAT = PREC_STRAT + PTR2D
             call MAPL_GetPointer(EXPORT, PTR2D, 'AN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_STRAT = PREC_STRAT + PTR2D
             PREC_STRAT = MAX(PREC_STRAT, 0.0)
          endif
       endif

       ! diagnosed convective precip (rain+snow)
       call MAPL_GetPointer(EXPORT, PREC_CONV, 'PREC_CONV', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PREC_CONV)) then
          if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
             PREC_CONV = CNV_FRC*TPREC
          else
             PREC_CONV = 0.0
             call MAPL_GetPointer(EXPORT, PTR2D, 'CN_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_CONV = PREC_CONV + PTR2D
             call MAPL_GetPointer(EXPORT, PTR2D, 'SC_PRCP'   , RC=STATUS); VERIFY_(STATUS)
             if (associated(PTR2D)) PREC_CONV = PREC_CONV + PTR2D
             PREC_CONV = MAX(PREC_CONV, 0.0)
          endif
       endif

     ! Diagnostic precip types: 
       call MAPL_GetPointer(EXPORT, ICE,   'ICE',   ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, FRZR,  'FRZR',  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (LUPDATE_PRECIP_TYPE .OR. LDIAGNOSE_PRECIP_TYPE) then
          call MAPL_GetPointer(EXPORT, PTYPE, 'PTYPE', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
          call DIAGNOSE_PRECIP_TYPE(IM, JM, LM, TPREC, PLS, PCU, RAIN, SNOW, ICE, FRZR, &
                                    PTYPE, PLE, T/PK, PK, PKE, ZL0, LUPDATE_PRECIP_TYPE)
       endif 
     ! Get Kuchera snow:rain ratios
       do I = 1,IM
          do J = 1,JM
              Tmax = 0.0
              do L =  LM, 1, -1
                 if (PLmb(I,J,L).gt.500.) then
                    Tmax = MAX(Tmax,T(I,J,L))
                 end if
              end do
              if (Tmax <= 271.16) then
                 TMP2D(I,J) = 12.0 + (271.16 - Tmax)
              else
                 TMP2D(I,J) = 12.0 + 2*(271.16 - Tmax)
              end if
              TMP2D(I,J) = max(0.0,TMP2D(I,J))
          end do
       end do
       call MAPL_GetPointer(EXPORT, PTR2D,'KUCHERA_RATIO', RC=STATUS); VERIFY_(STATUS)
       if(associated(PTR2D)) PTR2D=TMP2D
     ! Accumulated precip totals (mm), apply KUCHERA_RATIO for SNOW
       call MAPL_GetPointer(EXPORT, PTR2D,'SNOWTOTAL', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = TMP2D*DT_MOIST*(SNOW+ICE)
       call MAPL_GetPointer(EXPORT, PTR2D,'PRECTOTAL', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = DT_MOIST*TPREC

       call MAPL_GetPointer(EXPORT, PTR3D, 'QLTOT', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QLLS+QLCN

       call MAPL_GetPointer(EXPORT, PTR3D, 'QITOT', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QILS+QICN

       call MAPL_GetPointer(EXPORT, PTR3D, 'QCTOT', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = MIN(CLLS+CLCN,1.0)

       ! Cloud condensate exports
       call MAPL_GetPointer(EXPORT, PTR3D, 'QLLSX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QLLS

       call MAPL_GetPointer(EXPORT, PTR3D, 'QILSX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QILS

       call MAPL_GetPointer(EXPORT, PTR3D, 'QLCNX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QLCN

       call MAPL_GetPointer(EXPORT, PTR3D, 'QICNX1', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = QICN

       ! Fill wind, temperature & RH exports needed for SYNCTQ

       call MAPL_GetPointer(EXPORT, PTR3D, 'UAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = U

       call MAPL_GetPointer(EXPORT, PTR3D, 'VAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = V

       call MAPL_GetPointer(EXPORT, PTR3D, 'TAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T

       call MAPL_GetPointer(EXPORT, PTR3D, 'QAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = Q

       call MAPL_GetPointer(EXPORT, PTR3D, 'THAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T/PK

       call MAPL_GetPointer(EXPORT, PTR3D, 'SAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) then
          do L=1,LM
            PTR3D(:,:,L) = MAPL_CP*T(:,:,L) + MAPL_GRAV*(ZL0(:,:,L)+ZLE(:,:,LM))
          enddo
       endif

       call MAPL_GetPointer(EXPORT, PTR3D, 'RH2', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = MAX(MIN( Q/GEOS_QSAT (T, PLmb) , 1.02 ),0.0)
       call MAPL_GetPointer(EXPORT, PTR2D, 'CWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS+QICN+QILS )*MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'CLWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS ) *MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'LWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS) *MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'IWP', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( QICN+QILS ) *MASS , 3 )
       call MAPL_GetPointer(EXPORT, PTR2D, 'TPW', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = SUM( ( Q         ) *MASS , 3 )

       ! Lightning Exports
       call MAPL_GetPointer(EXPORT, PTR2D, 'LFR_GCC', NotFoundOk=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = 0.0

    else

       ! Internal State
       call MAPL_GetPointer(INTERNAL, Q,        'Q'    , RC=STATUS); VERIFY_(STATUS)
       ! Import State
       call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
       ! Allocatables
        ! Edge variables 
       ALLOCATE ( ZLE0 (IM,JM,0:LM) )
       ALLOCATE ( PLEmb(IM,JM,0:LM) )
        ! Layer variables
       ALLOCATE ( ZL0  (IM,JM,LM  ) )
       ALLOCATE ( PLmb (IM,JM,LM  ) )
       ALLOCATE ( PK   (IM,JM,LM  ) )
       ALLOCATE ( MASS (IM,JM,LM  ) )
       ALLOCATE ( TMP2D(IM,JM     ) )
       ! dervied states
       MASS     = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )/MAPL_GRAV
       call FILLQ2ZERO(Q, MASS, TMP2D)
       call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV_IN', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR2D)) PTR2D = TMP2D
       PLEmb    = PLE*.01
       PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
       PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
       DO L=0,LM
          ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
       END DO
       ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface

       ! Fill Wind, Temperature & RH exports needed for SYNCTQ

       call MAPL_GetPointer(EXPORT, PTR3D, 'UAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = U

       call MAPL_GetPointer(EXPORT, PTR3D, 'VAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = V

       call MAPL_GetPointer(EXPORT, PTR3D, 'TAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T
      
       call MAPL_GetPointer(EXPORT, PTR3D, 'QAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = Q
 
       call MAPL_GetPointer(EXPORT, PTR3D, 'THAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = T/PK
       
       call MAPL_GetPointer(EXPORT, PTR3D, 'SAFMOIST', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) then
          do L=1,LM
            PTR3D(:,:,L) = MAPL_CP*T(:,:,L) + MAPL_GRAV*(ZL0(:,:,L)+ZLE(:,:,LM))
          enddo
       endif
       
       call MAPL_GetPointer(EXPORT, PTR3D, 'RH2', RC=STATUS); VERIFY_(STATUS)
       if (associated(PTR3D)) PTR3D = MAX(MIN( Q/GEOS_QSAT (T, PLmb) , 1.02 ),0.0)

    endif

    call MAPL_TimerOff(MAPL,"TOTAL")

    RETURN_(ESMF_SUCCESS)

  end subroutine RUN

end module GEOS_MoistGridCompMod

