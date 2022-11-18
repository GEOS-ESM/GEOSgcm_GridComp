!$Id$

#include "MAPL_Generic.h"

module GEOS_OceanGridCompMod

!BOP
! !MODULE: GEOS_OceanGridCompMod -- Implements ESMF wrapper to invoke the DATASEA/MIT/MOM ocean models.

! !USES:
  use ESMF
  use MAPL
#ifdef BUILD_MIT_OCEAN
  use MIT_GEOS5PlugMod, only: MITSetServices => SetServices  ! this sets IRF
#endif
  use GEOS_DataSeaGridCompMod, only: DataSeaSetServices  => SetServices

  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices

  character(len=ESMF_MAXSTR)          :: OCEAN_NAME
  integer            :: DO_DATASEA
  real :: OrphanDepth

! !DESCRIPTION:
!
!   {\tt GuestOcean\_GridComp} is a light-weight gridded component that serves an
!   interface to ocean/data\_ocean components.
!
!EOP

  type :: T_PrivateState
     type(ESMF_Clock)  :: CLOCK
  end type T_PrivateState

  type :: T_PrivateState_Wrap
     type(T_PrivateState), pointer :: ptr
  end type T_PrivateState_Wrap

  integer ::          OCN 
  integer ::          OCNd 
  logical ::      DUAL_OCEAN


contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for GuestOcean

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices,
!       which sets the Run, Initialize, and Finalize services,
!       as well as allocating our instance of a generic state and putting it in the
!   gridded component (GC). Here we override all three methods and declare
!       the specs for the Imports and Export States (no MAPL controlled Internal State).
!       GuestOcean state variables (the bulletin board and the time) are kept
!       in the GuestOcean's Private Internal state.
!
!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Local vars
    type  (MAPL_MetaComp), pointer     :: MAPL
    type  (ESMF_Config)                :: CF
    integer ::      iDUAL_OCEAN
    character(len=ESMF_MAXSTR)         :: charbuf_
    character(len=ESMF_MAXSTR)         :: sharedObj

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL,       DO_DATASEA,     Label="USE_DATASEA:" ,       DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL,       OrphanDepth,    Label="OGCM_TOP_LAYER:" ,    DEFAULT=10.0, RC=STATUS)
    VERIFY_(STATUS)

    if(DO_DATASEA/=0) then
       OCEAN_NAME="DATASEA"
       OCN = MAPL_AddChild(GC, NAME=OCEAN_NAME, SS=DataSeaSetServices, RC=STATUS)
       VERIFY_(STATUS)
    else
       call MAPL_GetResource ( MAPL, OCEAN_NAME, Label="OCEAN_NAME:", DEFAULT="MOM", __RC__ )
       select case (trim(OCEAN_NAME))
          case ("MOM")
             call MAPL_GetResource ( MAPL, sharedObj,  Label="MOM_GEOS5PLUGMOD:", DEFAULT="libMOM_GEOS5PlugMod.so", __RC__ )
             OCN = MAPL_AddChild(OCEAN_NAME,'setservices_', parentGC=GC, sharedObj=sharedObj,  __RC__)
          case ("MOM6")
             call MAPL_GetResource ( MAPL, sharedObj,  Label="MOM6_GEOSPLUG:", DEFAULT="libMOM6_GEOSPlug.so", __RC__ )
             OCN = MAPL_AddChild(OCEAN_NAME,'setservices_', parentGC=GC, sharedObj=sharedObj,  __RC__)
#ifdef BUILD_MIT_OCEAN
          case ("MIT")
             OCN = MAPL_AddChild(GC, NAME=OCEAN_NAME, SS=MITSetServices, __RC__)
#endif
          case default
             charbuf_ = "OCEAN_NAME: " // trim(OCEAN_NAME) // " is not implemented, ABORT!"
             _ASSERT(.false., charbuf_)
       end select
    endif

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, __RC__ )
    DUAL_OCEAN = iDUAL_OCEAN /= 0

    OCNd = 0
    if (dual_ocean) then
       OCNd = MAPL_AddChild(GC, NAME="DATASEA", SS=DataSeaSetServices, RC=STATUS)
       VERIFY_(STATUS)
    endif

! Set the state variable specs.
! -----------------------------

!BOS

!  !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'FROCEAN',                           &
         LONG_NAME          = 'fraction_of_gridbox_covered_by_ocean',&
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'TAUX',                              &
         LONG_NAME          = 'Agrid_eastward_stress_on_ocean',     &
         UNITS              = 'N m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'TAUY',                              &
         LONG_NAME          = 'Agrid_northward_stress_on_ocean',    &
         UNITS              = 'N m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
        SHORT_NAME         = 'PENUVR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                            &
        SHORT_NAME         = 'PENPAR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                                &
        SHORT_NAME         = 'PENUVF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_UV_flux',  &
        UNITS              = 'W m-2',                             &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                                &
        SHORT_NAME         = 'PENPAF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_PAR_flux', &
        UNITS              = 'W m-2',                             &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
          LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DRNIR'                       ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
          LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DFNIR'                       ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
          SHORT_NAME         = 'SWHEAT',                            &
          LONG_NAME          = 'solar_heating_rate',                &
          UNITS              = 'W m-2',                             &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
          LONG_NAME          = 'river_discharge_at_ocean_points',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'DISCHARGE'                   ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    if (trim(OCEAN_NAME) == "MOM") then
       call MAPL_AddImportSpec(GC,                             &
           SHORT_NAME         = 'TR',                                &
           LONG_NAME          = 'tracer_mixing_ratios',              &
           UNITS              = '1',                                 &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
           DATATYPE           = MAPL_BundleItem,                     &
                                                          RC=STATUS  )
       VERIFY_(STATUS)
   
       call MAPL_AddImportSpec(GC,                             &
           SHORT_NAME         = 'TRFLUX',                            &
           LONG_NAME          = 'surface_fluxes_of_tracers',         &
           UNITS              = 'X',                                 &
           DIMS               = MAPL_DimsHorzOnly,                   &
           VLOCATION          = MAPL_VLocationNone,                  &
           DATATYPE           = MAPL_BundleItem,                     &
                                                          RC=STATUS  )
       VERIFY_(STATUS)
    endif

    call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'LWFLX'                   ,&
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SHFLX'                     ,&
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'QFLUX'                   ,&
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'ocean_snowfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SNOW'                   ,&
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                               &
        LONG_NAME          = 'ocean_rainfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'RAIN'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                    &
        SHORT_NAME         = 'FRESH',                         &
        LONG_NAME          = 'fresh_water_flux_due_to_ice_dynamics', &
          UNITS              = 'kg m-2 s-1'                ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FSALT',                         &
        LONG_NAME          = 'salt_flux_due_to_ice_dynamics', &
        UNITS              = 'kg m-2 s-1',                            &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'FHOCN',                         &
        LONG_NAME          = 'heat_flux_due_to_ice_dynamics', &
        UNITS              = 'W m-2',                            &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PEN_OCN',                           &
        LONG_NAME          = 'penetrated_shortwave_flux_at_the_bottom_of_first_ocean_model_layer',     &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RC=STATUS  )
     VERIFY_(STATUS)

    if (dual_ocean) then
       call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'FRACICEd',                           &
         LONG_NAME          = 'fractional_cover_of_seaice',        &
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
       VERIFY_(STATUS)
    endif

!  ! Need to have this internal state to fill in orphan points:

    call MAPL_AddInternalSpec(GC,                                 &
         SHORT_NAME         = 'TS_FOUND',                          &
         LONG_NAME          = 'foundation_temperature_for_interface_layer',&
         UNITS              = 'K',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         FRIENDLYTO         = trim(COMP_NAME), &
         DEFAULT            = 280.0, &
         RC=STATUS  )
    VERIFY_(STATUS)

!ALT Note the FRACICE from datasea is inhereted (in AMIP or dual_ocean)

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'MASKO',                             &
         LONG_NAME          = 'ocean_mask',                        &
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'SS_FOUND',                          &
         LONG_NAME          = 'foundation_salinity_for_interface_layer',&
         UNITS              = 'PSU',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'FRZMLT',                            &
         LONG_NAME          = 'freeze_melt_potential',             &
         UNITS              = 'W m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

! Diagnostics exports

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'TAUX',                              &
         LONG_NAME          = 'Agrid_eastward_stress_on_ocean',     &
         UNITS              = 'N m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'TAUY',                              &
         LONG_NAME          = 'Agrid_northward_stress_on_ocean',    &
         UNITS              = 'N m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'SWHEAT',                            &
         LONG_NAME          = 'solar_heating_rate',                &
         UNITS              = 'W m-2',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'RFLUX',                             &
         LONG_NAME          = 'downward_radiative_heat_flux_at_ocean_bottom',&
         UNITS              = 'W m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         LONG_NAME          = 'river_discharge_at_ocean_points',&
         UNITS              = 'kg m-2 s-1'                ,&
         SHORT_NAME         = 'DISCHARGE'                 ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'FROCEAN',                           &
         LONG_NAME          = 'fraction_of_gridbox_covered_by_ocean',&
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWFLX'                   ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
        LONG_NAME          = 'surface_net_downward_shortwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SWFLX'                   ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
         UNITS              = 'W m-2',                             &
        SHORT_NAME         = 'SHFLX'                     ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'QFLUX'                   ,&
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
        SHORT_NAME         = 'SFLX',                         &
        LONG_NAME          = 'salt_flux_due_to_ice_dynamics', &
        UNITS              = 'kg m-2 s-1',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
          SHORT_NAME         = 'RAIN',                              &
          LONG_NAME          = 'ocean_rainfall',&
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                 &
          SHORT_NAME         = 'SNOW',                              &
          LONG_NAME          = 'ocean_snowfall',&
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                     &    
          SHORT_NAME         = 'PEN_OCN',                           &    
          LONG_NAME          = 'penetrated_shortwave_flux_at_the_bottom_of_first_ocean_model_layer',&
          UNITS              = 'W m-2',                             &    
          DIMS               = MAPL_DimsHorzOnly,                   &    
          VLOCATION          = MAPL_VLocationNone,                  &    
          RC=STATUS  )
     VERIFY_(STATUS)

! Exports of child

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'TW',                                        &
         CHILD_ID   = OCN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'SW',                                        &
         CHILD_ID   = OCN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'UW',                                        &
         CHILD_ID   = OCN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'VW',                                        &
         CHILD_ID   = OCN,                                         &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    if(DO_DATASEA==0) then
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'DH',                                        &
            CHILD_ID   = OCN,                                         &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'UWB',                                        &
            CHILD_ID   = OCN,                                         &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'VWB',                                        &
            CHILD_ID   = OCN,                                         &
            RC=STATUS  )
       VERIFY_(STATUS)
       if (trim(OCEAN_NAME) == "MOM") then
          call MAPL_AddExportSpec ( GC   ,                          &
               SHORT_NAME = 'SSH',                                       &
               CHILD_ID   = OCN,                                         &
               RC=STATUS  )
          VERIFY_(STATUS)
       endif
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'SLV',                                       &
            CHILD_ID   = OCN,                                         &
            RC=STATUS  )
       VERIFY_(STATUS)
       if (trim(OCEAN_NAME) == "MOM") then
          call MAPL_AddExportSpec ( GC   ,                               &
               SHORT_NAME = 'PBO',                                       &
               CHILD_ID   = OCN,                                         &
               RC=STATUS  )
          VERIFY_(STATUS)
       endif

       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'T',                                         &
            CHILD_ID   = OCN,                                         &
            RC=STATUS  )
       VERIFY_(STATUS)
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'S',                                         &
            CHILD_ID   = OCN,                                         &
            RC=STATUS  )
       VERIFY_(STATUS)
    end if

!EOS

    if(DO_DATASEA==0) then
       call MAPL_TerminateImport    ( GC, SHORT_NAME=               &
          [character(len=9) :: 'TAUX  ','TAUY  ',                   &
            'PENUVR','PENPAR','PENUVF','PENPAF', 'DRNIR', 'DFNIR',  &
            'DISCHARGE', 'LWFLX', 'SHFLX', 'QFLUX', 'RAIN', 'SNOW', &
            'SFLX','SWHEAT'],                                       &  ! do not terminate import of PEN_OCN since it is not used in the `plug'
            CHILD=OCN,                          RC=STATUS  )
       VERIFY_(STATUS)
    end if

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=status)
    VERIFY_(STATUS)
! phase 1
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	  Run,        RC=status)
    VERIFY_(STATUS)
! phase 2 - this is only used in the predictor part of the replay for dual ocean
    if (DUAL_OCEAN) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	  Run,        RC=status)
       VERIFY_(STATUS)
    end if

! terminate child's import for a temperature correction - we will fill it here, if we run in dual_ocean mode, otherwise nobody needs this variable
    if(DUAL_OCEAN) then
       call MAPL_TerminateImport    ( GC,   &
            SHORT_NAME = (/'DEL_TEMP'/),                                    &
            CHILD      = OCN,                                        &
            RC=STATUS  )
       VERIFY_(STATUS)
    endif


!=============================================================================
! Generic SetServices--This creates the generic state and calls SetServices for children
!---------------------------------------------------------------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
  VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="INITIALIZE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUN"        ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="--ModRun"   ,RC=STATUS)
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices

! -----------------------------------------------------------------

!BOP

! !IROUTINE: INITIALIZE -- Initialize method for ExternalOcean wrapper

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),      intent(INOUT) :: GC     ! Gridded component
    type(ESMF_State),         intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),         intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),         intent(INOUT) :: CLOCK  ! The clock
    integer, optional,        intent(  OUT) :: RC     ! Error code:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)      :: IAm
    integer             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),     pointer   :: State
    type (ESMF_Grid)                    :: Grid
    type (T_PrivateState),    pointer   :: PrivateSTATE
    type (T_PrivateState_Wrap)          :: WRAP
    integer                             :: IM, JM, LM
    real                                :: DT

    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)
    type (ESMF_TimeInterval)            :: timeStep
    type (ESMF_Time)                    :: currTime

    real, pointer ::   MASK(:,:)
    real, pointer ::  MASKO(:,:)
    real, pointer :: MASK3D(:,:,:)
    real, pointer :: DH(:,:,:)

!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, grid=GRID, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, State, RC=STATUS)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(STATE,"INITIALIZE")
    call MAPL_TimerOn(STATE,"TOTAL"     )

! Get info from the Generic state
!--------------------------------

    call MAPL_Get(STATE,             &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
                                       RC=STATUS )
    VERIFY_(STATUS)


! Allocate the private state...
!------------------------------

    allocate( PrivateSTATE , stat=STATUS )
    VERIFY_(STATUS)

    wrap%ptr => PrivateState

! And put it in the GC
!---------------------

    CALL ESMF_UserCompSetInternalState( GC, TRIM(OCEAN_NAME)//'_internal_state', WRAP, STATUS )
    VERIFY_(status)

! Initialize the PrivateState. First the time...
!-----------------------------------------------
    call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    RC=STATUS)             ! Get AGCM Heartbeat
    VERIFY_(status)
    call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, RC=STATUS) ! set Default OCEAN_DT to AGCM Heartbeat
    VERIFY_(status)

    CALL ESMF_TimeIntervalSet(timeStep, S=NINT(DT), RC=status)
    VERIFY_(status)

    call ESMF_ClockGet(CLOCK, currTIME=currTime, RC=STATUS)
    VERIFY_(STATUS)

!ALT: check with Max about moving the clock 1 step forward
    PrivateState%clock = ESMF_ClockCreate(NAME = TRIM(OCEAN_NAME)//"Clock", &
         timeStep=timeStep, startTime=currTime, rc=status)
    VERIFY_(status)


! Initialize the Ocean Model.
!
!  This verifies the Grid and the Time against the private restarts.
!
!  Verifying that the GuestOcean grid and decomposition match those
!  inherited by the component. This is simply asserting that the
!  local im, jm, and lm are the same and making sure its internal
!  communication is consistent with the VM in the host's grid, probably
!  by initializing the guests's internal communication with the
!  communicator that comes from the VM in the Grid's Layout.
!
!  Note thet tha bulletin board states , GIM(:) and GEX(:), have been created and
!  populated with nodata fields. The ESMF arrays will be filled later.

!-----------------------------------------------------------------------

! Get sizes from my internal state
!---------------------------------
    call MAPL_Get(STATE, &
         IM=IM, &
         JM=JM, &
         LM=LM, &
         RC=STATUS)
    VERIFY_(STATUS)

! Once we know we have a valid  ESMF grid, we can call MAPL_GenericInitialize.
! This will allow us to use the built-in checkpoint/restart for our states.
!----------------------------------------------------------------------------


    call MAPL_TimerOff(STATE,"TOTAL"     )
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)
    call MAPL_TimerOn (STATE,"TOTAL"     )

    if(DO_DATASEA==0) then
       call MAPL_GetPointer(EXPORT, MASKO, 'MASKO'  , alloc=.true.,__RC__)

       select case (trim(OCEAN_NAME))
          case ("MOM")
             call MAPL_GetPointer(GEX(OCN), MASK3D, 'MOM_3D_MASK', __RC__)
             MASK => MASK3D(:,:,1)
          case ("MOM6")
             call MAPL_GetPointer(GEX(OCN), MASK, 'MOM_2D_MASK', __RC__)
          case ("MIT")
             call MAPL_GetPointer(GEX(OCN), MASK3D, 'MIT_3D_MASK', __RC__)
             MASK => MASK3D(:,:,1)
       end select
       if(associated(MASKO)) MASKO = MASK
    end if
 
    call MAPL_TimerOff(STATE,"TOTAL"     )
    call MAPL_TimerOff(STATE,"INITIALIZE")

! All Done
!---------
    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

! ========================================================

!BOP

! !IROUTINE: Run        -- Run method for ExternalModel wrapper

! !INTERFACE:

  subroutine Run ( gc, import, export, clock, rc )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component
    type(ESMF_State),    intent(INOUT) :: import ! Import state
    type(ESMF_State),    intent(INOUT) :: export ! Export state
    type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
    integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),     pointer   :: STATE
    type (ESMF_Time)                    :: EndTime
    type (ESMF_Time)                    :: MyTime,ct
    type (T_PrivateState),    pointer   :: PrivateSTATE
    type (T_PrivateState_Wrap)          :: WRAP
    type (ESMF_GridComp    ), pointer   :: GCS(:)
    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)

! Pointers to Imports

    real, pointer :: FROCEAN(:,:)
    real, pointer :: TAUXi(:,:)
    real, pointer :: TAUYi(:,:)
    real, pointer :: PENUVRi(:,:)
    real, pointer :: PENPARi(:,:)
    real, pointer :: PENUVFi(:,:)
    real, pointer :: PENPAFi(:,:)
    real, pointer :: DRNIRi(:,:)
    real, pointer :: DFNIRi(:,:)
    real, pointer :: HEATi(:,:,:)
    real, pointer :: DISCHARGEi(:,:)
    real, pointer :: LWFLXi(:,:)
    real, pointer :: SHFLXi(:,:)
    real, pointer :: QFLUXi(:,:)
    real, pointer :: SNOWi(:,:)
    real, pointer :: RAINi(:,:)
    real, pointer :: FHOCN(:,:)
    real, pointer :: FRESH(:,:)
    real, pointer :: FSALT(:,:)
    real, pointer :: PEN_OCN(:,:)

! Pointers to Exports

    real, pointer :: TS_FOUND (:,:)
    real, pointer :: SS_FOUND (:,:)
    real, pointer :: FRZMLTe(:,:)

! Diagnostics exports

    real, pointer :: RFLUX (:,:)
    real, pointer :: TAUXe (:,:)
    real, pointer :: TAUYe (:,:)
    real, pointer :: HEATe (:,:,:)
    real, pointer :: FROCEANe (:,:)
    real, pointer :: DISCHARGEe(:,:)
    real, pointer :: LWFLXe(:,:)
    real, pointer :: SWFLXe(:,:)
    real, pointer :: SHFLXe(:,:)
    real, pointer :: QFLUXe(:,:)
    real, pointer :: RAINe(:,:)
    real, pointer :: SNOWe(:,:)
    real, pointer :: SFLXe(:,:)
    real, pointer :: PEN_OCNe(:,:)


! Pointers to imports of child

    real, pointer :: TAUX(:,:)
    real, pointer :: TAUY(:,:)
    real, pointer :: PENUVR(:,:)
    real, pointer :: PENPAR(:,:)
    real, pointer :: PENUVF(:,:)
    real, pointer :: PENPAF(:,:)
    real, pointer :: DRNIR(:,:)
    real, pointer :: DFNIR(:,:)
    real, pointer :: HEAT(:,:,:)
    real, pointer :: DISCHARGE(:,:)
    real, pointer :: LWFLX(:,:)
    real, pointer :: SHFLX(:,:)
    real, pointer :: QFLUX(:,:)
    real, pointer :: RAIN(:,:)
    real, pointer :: SNOW(:,:)
    real, pointer :: SFLX(:,:)
    real, pointer :: FI(:,:)
    real, pointer :: FId(:,:)

! Pointers to exports of child

    real, pointer :: TW  (:,:)
    real, pointer :: SW  (:,:)
    real, pointer ::   MASK(:,:)
    real, pointer :: MASK3D(:,:,:)
    real, pointer :: FRZMLT(:,:)

    real, pointer :: TWd  (:,:)
    real, pointer :: DEL_TEMP (:,:)

! Locals

    integer           :: I,J,L
    integer           :: IM
    integer           :: JM
    integer           :: LM
    integer           :: NUM
    real, allocatable :: WGHT(:,:)
    real, pointer     :: WGHTi(:,:)
    real              :: DT, TAU_SST
    real              :: TAU_SST_UNDER_ICE
    real, pointer     :: LONS  (:,:)
    real, pointer     :: LATS  (:,:)
    real, parameter   :: OrphanSalinity=34.0  ! SA: ought to revisit this in context of OGCM rewrite. In reality one should update S as well.
    real              :: Tfreeze

    integer :: ID
    integer :: PHASE
    integer, allocatable :: PREDICTOR_CHLD(:)
    real, pointer     :: FR(:,:) => null()
    real, pointer     :: FRI(:,:,:) => null()
    character(len=ESMF_MAXSTR) :: replayMode

! Get the component's name and set-up traceback handle.
! -----------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name, currentPhase=PHASE, RC=status )
    VERIFY_(status)
    if (PHASE >= 10) PHASE = PHASE - 10 ! to be replaced by MAPL get_phase 
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, RC=status)
    VERIFY_(status)

! Profilers
!----------

    call MAPL_TimerOn (STATE,"RUN"  )
    call MAPL_TimerOn (STATE,"TOTAL")

! Get child's import ad export to use as a bulletin board
!--------------------------------------------------------
    call MAPL_Get(STATE,             &
         GCS       = GCS,                        &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
         LONS      = LONS,                       &
         LATS      = LATS,                       &
         IM        = IM,                         &
         JM        = JM,                         &
         LM        = LM,                         &
                                       RC=STATUS )
    VERIFY_(STATUS)


! Check the clocks to set set-up the "run-to" time
!-------------------------------------------------

    call ESMF_ClockGet( CLOCK, currTime=endTime, RC=STATUS)
    VERIFY_(status)

! Get GuestModel's private internal state
!---------------------------------

    CALL ESMF_UserCompGetInternalState( GC, TRIM(OCEAN_NAME)//'_internal_state', WRAP, STATUS )
    VERIFY_(STATUS)

    PrivateSTATE => WRAP%PTR

    call ESMF_ClockGet( PrivateState%CLOCK, currTime=myTime, RC=STATUS)
    VERIFY_(status)

    if (myTime > EndTime) then
       call ESMF_ClockSet(PrivateState%Clock,direction=ESMF_DIRECTION_REVERSE,rc=status)
       VERIFY_(status)
       do
         call ESMF_ClockAdvance(PrivateState%Clock,rc=status)
         VERIFY_(status)
         call ESMF_ClockGet(PrivateState%Clock,currTime=ct,rc=status)
         VERIFY_(status)
         if (ct==endTime) exit
       enddo
       call ESMF_ClockSet(PrivateState%Clock,direction=ESMF_DIRECTION_FORWARD,rc=status)
       VERIFY_(status)
       call ESMF_ClockGet( PrivateState%CLOCK, currTime=myTime, RC=STATUS)
       VERIFY_(status)
    end if

    if( MyTime <= EndTime ) then ! Time to run

! We get the ocean-land mask (now computed in Initialize of Plug)
! ---------------------------------------------------------------
       if(DO_DATASEA==0) then
          select case(trim(OCEAN_NAME))
             case ("MOM")
                call MAPL_GetPointer(GEX(OCN), MASK3D, 'MOM_3D_MASK', __RC__)
                MASK => MASK3D(:,:,1)
             case ("MIT")
                call MAPL_GetPointer(GEX(OCN), MASK3D, 'MIT_3D_MASK', __RC__)
                MASK => MASK3D(:,:,1)
             case ("MOM6")
                call MAPL_GetPointer(GEX(OCN), MASK, 'MOM_2D_MASK', __RC__)
             end select
       else
          allocate(MASK3D(IM,JM,LM), STAT=STATUS); VERIFY_(STATUS)
          MASK3D=1.0
          allocate(MASK(IM,JM), STAT=STATUS); VERIFY_(STATUS)
          MASK=1.0
       end if

! Get ocean time step and misc. parameters
!-----------------------------------------

       call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    RC=STATUS)             ! Get AGCM Heartbeat
       VERIFY_(status)
       call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, RC=STATUS) ! set Default OCEAN_DT to AGCM Heartbeat
       VERIFY_(status)

! Get pointers to imports
!--------------------------------------------------------------------------------
       call MAPL_GetPointer(IMPORT, FROCEAN, 'FROCEAN', RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, TAUXi, 'TAUX'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, TAUYi, 'TAUY'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, PENUVRi, 'PENUVR'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, PENPARi, 'PENPAR'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, PENUVFi, 'PENUVF'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, PENPAFi, 'PENPAF'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, DRNIRi, 'DRNIR'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, DFNIRi, 'DFNIR'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, HEATi, 'SWHEAT' , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, DISCHARGEi, 'DISCHARGE'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, LWFLXi, 'LWFLX'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, SHFLXi, 'SHFLX'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, QFLUXi, 'QFLUX'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, SNOWi, 'SNOW'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, RAINi, 'RAIN'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FHOCN, 'FHOCN'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FRESH, 'FRESH'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, FSALT, 'FSALT'   , RC=STATUS); VERIFY_(STATUS)

! Get pointers from ImExState
!----------------------------
       if(DO_DATASEA==0) then
          call MAPL_GetPointer(GIM(OCN), TAUX, 'TAUX'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), TAUY, 'TAUY'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), PENUVR, 'PENUVR'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), PENPAR, 'PENPAR'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), PENUVF, 'PENUVF'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), PENPAF, 'PENPAF'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), DRNIR, 'DRNIR'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), DFNIR, 'DFNIR'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), HEAT, 'SWHEAT', RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), DISCHARGE, 'DISCHARGE'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), LWFLX, 'LWFLX'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), SHFLX, 'SHFLX'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), QFLUX, 'QFLUX'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), RAIN, 'RAIN'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), SNOW, 'SNOW'  , RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(OCN), SFLX, 'SFLX'  , RC=STATUS); VERIFY_(STATUS) ! and do not add import of PEN_OCN here since it is not used in the `plug'
       end if

       call MAPL_GetPointer(IMPORT, PEN_OCN, 'PEN_OCN',RC=STATUS); VERIFY_(STATUS)

       call MAPL_GetPointer(GEX(OCN), TW,   'TW'  , alloc=.true., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(GEX(OCN), SW,   'SW'  , alloc=.true., RC=STATUS); VERIFY_(STATUS)

       if (dual_ocean) then
          call MAPL_GetPointer(GEX(OCNd), TWd,   'TW'  , alloc=.true., RC=STATUS); VERIFY_(STATUS)
          call MAPL_GetPointer(IMPORT, FId, 'FRACICEd'   , RC=STATUS); VERIFY_(STATUS)
       end if
       
       if(DO_DATASEA==0) then
          call MAPL_GetPointer(GEX(OCN), FRZMLT,   'FRZMLT'  , alloc=.true., RC=STATUS); VERIFY_(STATUS)
       end if

! Get pointers to exports
!--------------------------------------------------------

       call MAPL_GetPointer(EXPORT, TS_FOUND,'TS_FOUND', RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SS_FOUND,'SS_FOUND', RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, FRZMLTe,'FRZMLT', RC=STATUS); VERIFY_(STATUS)

! Diagnostics exports
!---------------------------------------------------------
       call MAPL_GetPointer(EXPORT, RFLUX,  'RFLUX' , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, FROCEANe,'FROCEAN', RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, TAUXe, 'TAUX'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, TAUYe, 'TAUY'   , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, HEATe, 'SWHEAT' , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, DISCHARGEe, 'DISCHARGE'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, LWFLXe, 'LWFLX'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SWFLXe, 'SWFLX'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SHFLXe, 'SHFLX'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, QFLUXe, 'QFLUX'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, RAINe, 'RAIN'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SNOWe, 'SNOW'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, SFLXe, 'SFLX'  , RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetPointer(EXPORT, PEN_OCNe,'PEN_OCN', RC=STATUS); VERIFY_(STATUS)

       if(associated(FROCEANe)) FROCEANe = FROCEAN

! Allocate space for temporary arrays
!------------------------------------

       allocate(WGHT(IM,JM), STAT=STATUS); VERIFY_(STATUS)

! Weight for ocean grid
!----------------------
       wght=0.0
       where(MASK>0 .and. FROCEAN /= MAPL_UNDEF)
          WGHT = FROCEAN/MASK
       elsewhere
          WGHT = 0.0
       end where

#ifdef BUILD_MIT_OCEAN
       call MAPL_GetPointer(GIM(OCN), WGHTi, 'WGHT', __RC__)
       WGHTi = WGHT
#endif

       if(DO_DATASEA==0) then
! Copy imports into ImEx variables
!---------------------------------
          PENUVR = PENUVRi * WGHT
          PENPAR = PENPARi * WGHT
          PENUVF = PENUVFi * WGHT
          PENPAF = PENPAFi * WGHT
          DRNIR = DRNIRi * WGHT
          DFNIR = DFNIRi * WGHT
          DISCHARGE = DISCHARGEi * WGHT
          LWFLX = LWFLXi * WGHT
          QFLUX = QFLUXi * WGHT
          SHFLX = (SHFLXi- FHOCN) * WGHT
          RAIN = (RAINi+FRESH) * WGHT
          SNOW = SNOWi * WGHT
          SFLX = FSALT * WGHT

! This stress forces the ocean, combined with sea ice bottom stress later
!------------------------------------------------------------------------
          TAUX = TAUXi * WGHT
          TAUY = TAUYi * WGHT


! Prepare radiative heating for ocean
!------------------------------------

          if(associated(RFLUX )) RFLUX  = 0.0
          select case (trim(OCEAN_NAME))
             case ("MOM", "DATASEA")
                do L=1,LM
                   HEAT(:,:,L) = HEATi(:,:,L)*WGHT
                   if(associated(RFLUX)) then
                      RFLUX = RFLUX + (1.0-MASK3D(:,:,L))*HEAT(:,:,L)
                   end if
                end do
             case ("MOM6")
                ! No 3D Mask from MOM6. Do nothing for now!
          end select

          if (associated(HEATe)) HEATe = HEATi
          if (associated(TAUXe)) TAUXe = TAUXi
          if (associated(TAUYe)) TAUYe = TAUYi
          if (associated(DISCHARGEe)) DISCHARGEe = DISCHARGEi
          if (associated(LWFLXe)) LWFLXe = LWFLXi
          if (associated(SWFLXe)) SWFLXe = PENUVRi+PENPARi+PENUVFi+PENPAFi+DRNIRi+DFNIRi-PEN_OCN
          if (associated(SHFLXe)) SHFLXe = SHFLXi
          if (associated(QFLUXe)) QFLUXe = QFLUXi
          if (associated(RAINe)) RAINe = RAINi
          if (associated(SNOWe)) SNOWe = SNOWi
          if (associated(SFLXe)) SFLXe = SFLX
          if (associated(PEN_OCNe)) PEN_OCNe = PEN_OCN
       end if !DO_DATASEA

! Loop the ocean model
!---------------------

       NUM = 0
       do while ( MyTime <= endTime )

! Run ocean for one time step (DT)
!---------------------------------

          call MAPL_TimerOff(STATE,"TOTAL")
          call MAPL_TimerOn (STATE,"--ModRun")

          if (.not. DUAL_OCEAN) then
             call MAPL_GenericRunChildren(GC, IMPORT, EXPORT, PrivateState%CLOCK, RC=STATUS)
             VERIFY_(STATUS)
          else
             if (PHASE == 1) then
                ! corrector
                call ESMF_GridCompRun( GCS(OCNd), importState=GIM(OCNd), &
                     exportState=GEX(OCNd), clock=CLOCK, phase=1, userRC=STATUS)
                VERIFY_(STATUS)
                call MAPL_GenericRunCouplers( STATE, CHILD=OCNd, CLOCK=CLOCK, RC=STATUS )
                VERIFY_(STATUS)
                call ESMF_GridCompRun( GCS(OCN), importState=GIM(OCN), &
                     exportState=GEX(OCN), clock=CLOCK, phase=1, userRC=STATUS)
                VERIFY_(STATUS)
                call MAPL_GenericRunCouplers( STATE, CHILD=OCN, CLOCK=CLOCK, RC=STATUS )
                VERIFY_(STATUS)
             else
                ! predictor
                call ESMF_GridCompRun( GCS(OCNd), importState=GIM(OCNd), &
                     exportState=GEX(OCNd), clock=CLOCK, phase=1, userRC=STATUS)
                VERIFY_(STATUS)
                call MAPL_GenericRunCouplers( STATE, CHILD=OCNd, CLOCK=CLOCK, RC=STATUS )
                VERIFY_(STATUS)
             end if
          end if

          if (DUAL_OCEAN .and. PHASE == 1) then
             ! calculate temperature correction to send back to MOM
             call MAPL_GetPointer(GIM(OCN), DEL_TEMP, 'DEL_TEMP', RC=STATUS)
             VERIFY_(STATUS)
             call MAPL_GetPointer(GIM(OCNd), FI , 'FRACICE'  , RC=STATUS)
             VERIFY_(STATUS)
             
             call MAPL_GetResource(STATE,TAU_SST, Label="TAU_SST:", default=432000.0 ,RC=STATUS)
             VERIFY_(status)

             call MAPL_GetResource(STATE,TAU_SST_UNDER_ICE, Label="TAU_SST_UNDER_ICE:", default=86400.0 ,RC=STATUS)
             VERIFY_(status)

             ! we should have valid pointers to TW and TWd by now
             
             DEL_TEMP = 0.0 ! we do not want uninitiazed variables
             where(MASK > 0.0 .and. FI < 0.05)

                ! what about relaxation
                DEL_TEMP = (TWd - TW)*DT/(DT+TAU_SST)

             end where

             where(MASK > 0.0 .and. FI >= 0.05 .and. FId > FI)
                ! 0.054 (C/psu) is the ratio between the freezing temperature and salinity of brine.
                ! -0.054*SW gives salinity dependent freezing temperature 
                ! ideally this const should be from the ocean model, but doing so is difficult here
                DEL_TEMP = ((-0.054*SW+MAPL_TICE) - TW)*DT/(DT+TAU_SST_UNDER_ICE)

             end where

             ! put it back to MOM
             call ESMF_GridCompRun( GCS(OCN), importState=GIM(OCN), &
               exportState=GEX(OCN), clock=CLOCK, phase=2, userRC=STATUS )
             VERIFY_(STATUS)
          end if
          
          call MAPL_TimerOff(STATE,"--ModRun")
          call MAPL_TimerOn (STATE,"TOTAL")

! Bump the time in the internal state
!------------------------------------

          call ESMF_ClockAdvance( PrivateState%clock,                    rc=status)
          VERIFY_(status)
          call ESMF_ClockGet    ( PrivateState%clock, currTime= myTime , rc=status)
          VERIFY_(status)

          NUM = NUM + 1

       end do

       if(associated(SS_FOUND)) then
          SS_FOUND = OrphanSalinity
          where(WGHT > 0.0)
             SS_FOUND = SW
          end where
       end if

       if(associated(FRZMLTe)) then
          if(DO_DATASEA == 0) then
            if(trim(OCEAN_NAME) == "MIT") then
! for now, when MIT fill frzmlt here (not filled inside)
! and also for now we set the depth of the top ocean level to 50. m in MIT plug
             where(WGHT > 0.0 )
               FRZMLTe = (MAPL_TICE-0.054*SW - TS_FOUND) * (MAPL_RHO_SEAWATER*MAPL_CAPWTR*50.)/DT
             end where
            else
! when NOT MIT assume frzmlt filled inside and get from ocean export
             where(WGHT > 0.0 )
                FRZMLTe = FRZMLT
             end where
            endif
          else
             FRZMLTe = 0.0
          end if          
       end if

       if (DUAL_OCEAN) then
          !ALT we might not have FI yet, so let get it again
          call MAPL_GetPointer(GIM(OCNd), FI , 'FRACICE'  , RC=STATUS)
          VERIFY_(STATUS)
          where(WGHT > 0.0)
             where(FI < 0.05)
                TS_FOUND = TWd
             elsewhere
                TS_FOUND = TW
             end where
          end where
       else
          where(WGHT > 0.0)
             TS_FOUND = TW
          end where
       end if

! Update orphan points
       if(DO_DATASEA == 0) then
          WGHT=FROCEAN*(1-MASK)
          Tfreeze=MAPL_TICE-0.054*OrphanSalinity

          where(wght>0.0)
             TS_FOUND=TS_FOUND+ &
                      DT*(LWFLXi+(PENUVRi+PENPARi+PENUVFi+PENPAFi+DRNIRi+DFNIRi - PEN_OCN)-SHFLXi-QFLUXi*MAPL_ALHL-MAPL_ALHF*SNOWi+FHOCN)/(OrphanDepth*MAPL_RHO_SEAWATER*MAPL_CAPWTR) ! explicit update in time
             FRZMLTe = (Tfreeze - TS_FOUND) * (MAPL_RHO_SEAWATER*MAPL_CAPWTR*OrphanDepth)/DT
             TS_FOUND=max(TS_FOUND, Tfreeze)
          end where

       end if

       deallocate(WGHT, STAT=STATUS); VERIFY_(STATUS)

       if(DO_DATASEA/=0) then
          deallocate(MASK3D, STAT=STATUS); VERIFY_(STATUS)
          deallocate(MASK, STAT=STATUS); VERIFY_(STATUS)
       end if

    end if ! Time to run

! Profilers
!----------

    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOff(STATE,"RUN"  )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

end module GEOS_OceanGridCompMod
