!  $Id$
#include "MAPL_Generic.h"
!=============================================================================
!BOP

! !MODULE: GEOS_Ogcm   -- A composite component for the ogcm components.

! !INTERFACE:

module GEOS_OgcmGridCompMod

! !USES:

  use ESMF
  use MAPL

  use GEOS_OceanBioGeoChemGridCompMod,   only : ObioSetServices    => SetServices
  use GEOS_OradBioGridCompMod,           only : OradBioSetServices => SetServices
  use GEOS_OradGridCompMod,              only : OradSetServices    => SetServices

  use GEOS_OceanGridCompMod,             only : OceanSetServices => SetServices
  use GEOS_SeaIceGridCompMod,            only : SeaIceSetServices => SetServices

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

! !DESCRIPTION:
! 
!   {\tt GEOS\_Ogcm} is a light-weight gridded component that implements the
!      interface to the ogcm components. The ogcm computational components
!      (MOMx/Ocean, CICEx/GEOS_Seaice, OceanRadiation, OceanBioGeochemistry, etc)
!      are its children.
!      This component currently serves as an interface between the exchange
!      grid and the ocean's grid. Its ``natural'' grid is the ocean part of the 
!      exchange grid, and all its imports and exports are on this grid. The natural
!      grid of all of its children is currently the ocean's rectangular grid.
!      The ESMF grid that is in the gridded component is created by the parent
!      and it is the ocean's rectangular grid. At present the exchange grid information
!      is kept in the generic state.
!
!      The fact that some of these are friendlies---all the ``skin'' 
!      components---means that it cannot be a ``no-work-no-change'' component.
!      The interpolation of these to the ocean grid leave an ocean grid imprint
!      on them. No such happens on the atmospheric side. So we should think of these
!      exchange grid friendlies as ocean variables.  
!
!EOP

  integer, parameter :: NUM_SNOW_LAYERS=1
  integer            :: NUM_ICE_CATEGORIES
  integer            :: NUM_ICE_LAYERS

  integer            :: DO_CICE_THERMO
  integer            :: DO_DATASEAONLY
  integer            :: DO_DATAICE
  integer            :: DO_OBIO
  logical            :: DO_DATA_ATM4OCN

  logical          :: ocean_extData
  logical          :: ocean_sssData
  logical          :: seaIceT_extData

  integer, parameter :: NUM_DUDP = 5
  integer, parameter :: NUM_DUWT = 5
  integer, parameter :: NUM_DUSD = 5
  integer, parameter :: NUM_BCDP = 2
  integer, parameter :: NUM_BCWT = 2
  integer, parameter :: NUM_OCDP = 2
  integer, parameter :: NUM_OCWT = 2
  integer, parameter :: NB_OBIO  = 33 !total number of bands for OradBio
  integer, parameter :: NB_CHOU_UV   = 5 ! Number of UV bands
  integer, parameter :: NB_CHOU_NIR  = 3 ! Number of near-IR bands
  integer, parameter :: NB_CHOU      = NB_CHOU_UV + NB_CHOU_NIR ! Total number of bands
!--------

  character(len=ESMF_MAXSTR)          :: OCEAN_NAME
!=============================================================================

  integer ::        OBIO
  integer ::        ORAD
  integer ::      SEAICE
  integer ::       OCEAN

  logical ::      DUAL_OCEAN

  type T_OGCM_STATE
     private
     logical :: useInterp = .false.
  end type T_OGCM_STATE

! Wrapper for extracting internal state
! -------------------------------------
  type OGCM_WRAP
     type (T_OGCM_STATE), pointer :: PTR => null()
  end type OGCM_WRAP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION: This version uses the MAPL\_GenericSetServices, which in addition
!                to setting default IRF methods, also allocates
!   our instance of a generic state and puts it in the 
!   gridded component (GC). Here we override the Initialize and Run methods.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    type (MAPL_MetaComp),  pointer          :: MAPL
    type (ESMF_Config)                      :: CF

    type (T_OGCM_STATE), pointer            :: ogcm_internal_state => null() 
    type (OGCM_wrap)                        :: wrap

    integer ::      iDUAL_OCEAN
!=============================================================================

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

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, RC=STATUS )
    VERIFY_(STATUS)
    DUAL_OCEAN = iDUAL_OCEAN /= 0

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL,       DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,       DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    if (DO_CICE_THERMO /= 0) then
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
       VERIFY_(STATUS)
       if (DO_CICE_THERMO == 1) then
          call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS,     Label="CICE_N_ICE_LAYERS:" ,     RC=STATUS)
          VERIFY_(STATUS)
       endif 
    else
       NUM_ICE_CATEGORIES = 1
       NUM_ICE_LAYERS     = 1
    endif

! this get resource is repeated in Ocean - change both together!
    call MAPL_GetResource ( MAPL, DO_DATASEAONLY, Label="USE_DATASEA:" ,        DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DO_DATAICE,     Label="USE_DATASEAICE:" ,     DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DO_OBIO,        Label="USE_OCEANOBIOGEOCHEM:",DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DO_DATA_ATM4OCN, Label="USE_DATA_ATM4OCN:" ,  DEFAULT=.FALSE.,     __RC__ )
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, OCEAN_NAME,     Label="OCEAN_NAME:",          DEFAULT="MOM", __RC__ )
    
! following logic is to make sure: configuration of cetain components (CICE, OBIO, etc) has associated models "active."
    if (DO_DATA_ATM4OCN) then
       _ASSERT(DO_DATASEAONLY==0,'needs informative message')
    end if
    if (DO_DATASEAONLY/=0) then
       _ASSERT(DO_CICE_THERMO==0,'needs informative message')
       _ASSERT(DO_DATAICE    /=0,'needs informative message')
       _ASSERT(DO_OBIO       ==0,'needs informative message')
    end if

    call MAPL_GetResource (MAPL,   ocean_extData, Label="OCEAN_EXT_DATA:",   DEFAULT=.FALSE., __RC__ ) ! .TRUE. or .FALSE.
    if (DO_DATASEAONLY==1) then ! Fake-ocean (i.e., data ocean). 
      call MAPL_GetResource (MAPL, ocean_sssData,  Label="OCEAN_SSS_DATA:",  DEFAULT=.FALSE., __RC__ ) ! .TRUE. or .FALSE.
      call MAPL_GetResource (MAPL, seaIceT_extData,Label="SEAICE_THICKNESS_EXT_DATA:", DEFAULT=.FALSE., _RC ) ! .TRUE. or .FALSE.
    endif

! Set the Run and initialize entry points
!----------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run       , RC=STATUS )
    VERIFY_(STATUS)
    if (DUAL_OCEAN) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run       , RC=STATUS )
       VERIFY_(STATUS)
    end if

! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------

    if (DO_OBIO/=0) then
       OBIO = MAPL_AddChild(GC, NAME='OBIO', SS=ObioSetServices, RC=STATUS)
       VERIFY_(STATUS)
       ORAD = MAPL_AddChild(GC, NAME='ORAD', SS=OradBioSetServices, RC=STATUS)
       VERIFY_(STATUS)
    else
       OBIO = 0
       ORAD = MAPL_AddChild(GC, NAME='ORAD', SS=OradSetServices, RC=STATUS)
       VERIFY_(STATUS)
    end if
    
    SEAICE = MAPL_AddChild(GC, NAME='SEAICE', SS=SeaIceSetServices, RC=STATUS)
    VERIFY_(STATUS)

    OCEAN = MAPL_AddChild(GC, NAME='OCEAN', SS=OceanSetServices, RC=STATUS) 
    VERIFY_(STATUS)
   
! Set the state variable specs.
! -----------------------------

!BOS

!  !IMPORT STATE:

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'eastward_stress_on_ocean'          ,&
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUXW'                             ,&
    DIMS               = MAPL_DimsTileOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'northward_stress_on_ocean',         &
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUYW'                             ,&
    DIMS               = MAPL_DimsTileOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'eastward_stress_on_ice'            ,&
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUXI'                             ,&
    DIMS               = MAPL_DimsTileOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'northward_stress_on_ice',           &
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUYI'                             ,&
    DIMS               = MAPL_DimsTileOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'ocean_ustar_cubed',                 &
    UNITS              = 'm+3 s-3'                           ,&
    SHORT_NAME         = 'OUSTAR3'                           ,&
    DIMS               = MAPL_DimsTileOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          ) 
  VERIFY_(STATUS)


  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME           = 'surface_air_pressure',             &
    UNITS               = 'Pa',                               &
    SHORT_NAME          = 'PS',                               &
    DIMS                = MAPL_DimsTileOnly,                  &
    VLOCATION           = MAPL_VLocationNone,                 &
    DEFAULT             = 100000.,                            &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                             &
       SHORT_NAME         = 'PENUVR',                            &
       LONG_NAME          = 'net_downward_penetrating_direct_UV_flux',  &
       UNITS              = 'W m-2',                             &
       DIMS               = MAPL_DimsTileOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)
  
  call MAPL_AddImportSpec(GC,                             &
       SHORT_NAME         = 'PENPAR',                            &
       LONG_NAME          = 'net_downward_penetrating_direct_PAR_flux', &
       UNITS              = 'W m-2',                             &
       DIMS               = MAPL_DimsTileOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)
  
  call MAPL_AddImportSpec(GC,                             &
       SHORT_NAME         = 'PENUVF',                            &
       LONG_NAME          = 'net_downward_penetrating_diffuse_UV_flux',  &
       UNITS              = 'W m-2',                             &
       DIMS               = MAPL_DimsTileOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                             &
       SHORT_NAME         = 'PENPAF',                            &
       LONG_NAME          = 'net_downward_penetrating_diffuse_PAR_flux', &
       UNITS              = 'W m-2',                             &
       DIMS               = MAPL_DimsTileOnly,                   &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                         ,&
       LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
       UNITS              = 'W m-2'                       ,&
       SHORT_NAME         = 'DRNIR'                       ,&
       DIMS               = MAPL_DimsTileOnly             ,&
       VLOCATION          = MAPL_VLocationNone            ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC                         ,&
       LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
       UNITS              = 'W m-2'                       ,&
       SHORT_NAME         = 'DFNIR'                       ,&
       DIMS               = MAPL_DimsTileOnly             ,&
       VLOCATION          = MAPL_VLocationNone            ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                    &
       LONG_NAME          = 'river_discharge_at_ocean_points',&
       UNITS              = 'kg m-2 s-1'                ,&
       SHORT_NAME         = 'DISCHRG'                   ,&
       DIMS               = MAPL_DimsTileOnly           ,&
       VLOCATION          = MAPL_VLocationNone          ,&
       RESTART            = MAPL_RestartSkip            ,&
       RC=STATUS  ) 
  VERIFY_(STATUS)

  if (DO_OBIO/=0) then
    call OBIO_SetServices(DO_DATA_ATM4OCN, RC)
  end if
  
! These are supposed to be friendly to us
!----------------------------------------------

  if (.not. seaIceT_extData) then
    if (DO_CICE_THERMO <= 1) then  
         call MAPL_AddImportSpec(GC,                               &
              SHORT_NAME         = 'HI',                           &
              LONG_NAME          = 'seaice_skin_layer_mass',       &
              UNITS              = 'kg',                           &
              DIMS               = MAPL_DimsTileOnly,              &
              VLOCATION          = MAPL_VLocationNone,             &
              DEFAULT            = 0.0,                            &
              _RC)
    endif

    call MAPL_AddImportSpec(GC,                            &
      SHORT_NAME         = 'SI',                           &
      LONG_NAME          = 'seaice_skin_salinity',         &
      UNITS              = 'psu',                          &
      DIMS               = MAPL_DimsTileOnly,              &
      VLOCATION          = MAPL_VLocationNone,             &
      DEFAULT            = 0.0,                            &
      _RC)
  endif

  if (DO_CICE_THERMO /= 0) then  
     call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME         = 'TI',                                &
          LONG_NAME          = 'seaice_skin_temperature',           &
          UNITS              = 'K',                                 &
          DIMS               = MAPL_DimsTileOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          DEFAULT            = MAPL_TICE,                           &
          _RC)
  else
     if (.not. seaIceT_extData) then
        call MAPL_AddImportSpec(GC,                                 &
             SHORT_NAME         = 'TI',                             &
             LONG_NAME          = 'seaice_skin_temperature',        &
             UNITS              = 'K',                              &
             DIMS               = MAPL_DimsTileOnly,                &
             VLOCATION          = MAPL_VLocationNone,               &
             DEFAULT            = MAPL_TICE,                        &
             _RC)
     endif
  endif

  if (DO_CICE_THERMO == 1) then  
     call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME         = 'FRACICE',                         &
          LONG_NAME          = 'fractional_cover_of_seaice',        &
          UNITS              = '1',                                 &
          DIMS               = MAPL_DimsTileOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          _RC)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'VOLICE',                            &
          LONG_NAME          = 'ice_category_volume_per_unit_area_of_grid_cell',&
          UNITS              = 'm',                                 &
          DIMS               = MAPL_DimsTileOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          DEFAULT            = 0.0,                                 &
          _RC)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'VOLSNO',                            &
          LONG_NAME          = 'sno_category_volume_per_unit_area_of_grid_cell',&
          UNITS              = 'm',                                 &
          DIMS               = MAPL_DimsTileOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          DEFAULT            = 0.0,                                 &
          _RC)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'ERGICE',                            &
          LONG_NAME          = 'ice_category_layer_internal_energy',&
          UNITS              = 'J m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS,NUM_ICE_CATEGORIES/),&
          DEFAULT            = 0.0,                                 &
          _RC)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'ERGSNO',                            &
          LONG_NAME          = 'snow_category_layer_internal_energy',&
          UNITS              = 'J m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS,NUM_ICE_CATEGORIES/),&
          DEFAULT            = 0.0,                                 &
          _RC)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'TAUAGE',                            &
          LONG_NAME          = 'volume_weighted_mean_ice_age',      &
          UNITS              = 's',                                 &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          DIMS               = MAPL_DimsTileOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          DEFAULT            = 0.0,                                 &
          _RC)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'MPOND',                            &
          LONG_NAME          = 'pond_volume',                       &
          UNITS              = 'm',                                 &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          DIMS               = MAPL_DimsTileOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          DEFAULT            = 0.0,                                 &
          _RC)
  endif

  call MAPL_AddImportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWFLX'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                     &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHFLX'                     ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                     &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'QFLUX'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                     &
        LONG_NAME          = 'ocean_snowfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'SNOW'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
                                               RC=STATUS  ) 
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                     &
        LONG_NAME          = 'ocean_rainfall'            ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'RAIN'                   ,&
        DIMS               = MAPL_DimsTileOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
  VERIFY_(STATUS)

  if (DO_CICE_THERMO <= 1) then  
    call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'FRESH',                           &
         LONG_NAME          = 'fresh_water_flux_due_to_thermodynamics', &
         UNITS              = 'kg m-2 s-1',                      &
         DIMS               = MAPL_DimsTileOnly,                     &
         VLOCATION          = MAPL_VLocationNone,                &
         _RC)

  call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'FSALT',                           &
         LONG_NAME          = 'salt_flux_due_to_thermodynamics', &
         UNITS              = 'kg m-2 s-1',                      &
         DIMS               = MAPL_DimsTileOnly,                     &
         VLOCATION          = MAPL_VLocationNone,                &
         _RC)

  call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'FHOCN',                           &
         LONG_NAME          = 'heat_flux_due_to_thermodynamics', &
         UNITS              = 'W m-2',                           &
         DIMS               = MAPL_DimsTileOnly,                     &
         VLOCATION          = MAPL_VLocationNone,                &
         _RC)
  endif

  call MAPL_AddImportSpec(GC,                                  &
         SHORT_NAME         = 'PEN_OCN',                         &
         LONG_NAME          = 'penetrated_shortwave_flux_at_the_bottom_of_first_ocean_model_layer', &
         UNITS              = 'W m-2',                           &
         DIMS               = MAPL_DimsTileOnly,                 &
         VLOCATION          = MAPL_VLocationNone,                &
         RC=STATUS  )
  VERIFY_(STATUS)

!  !EXPORT STATE:

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'UW',                                &
    LONG_NAME          = 'zonal_velocity_of_surface_water',   &
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VW',                                &
    LONG_NAME          = 'meridional_velocity_of_surface_water',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'UI',                                &
    LONG_NAME          = 'zonal_velocity_of_surface_seaice',  &
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VI',                                &
    LONG_NAME          = 'meridional_velocity_of_surface_seaice',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TILELONS',                          &
    LONG_NAME          = 'longitude',                         &
    UNITS              = 'degrees',                           &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TILELATS',                          &
    LONG_NAME          = 'latitude',                          &
    UNITS              = 'degrees',                           &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'KPAR',                              &
    LONG_NAME          = 'PAR_extinction_coefficient',        &
    UNITS              = 'm-1 ',                              &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'TS_FOUND',                          &
    LONG_NAME          = 'foundation_temperature_for_interface_layer',&
    UNITS              = 'K',                                 &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'SS_FOUND',                          &
    LONG_NAME          = 'foundation_salinity_for_interface_layer',&
    UNITS              = 'PSU',                               &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                 &
    SHORT_NAME         = 'FRZMLT',                            &
    LONG_NAME          = 'freeze_melt_potential',             &
    UNITS              = 'W m-2',                             &
    DIMS               = MAPL_DimsTileOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  if (DO_CICE_THERMO == 0) then  
     call MAPL_AddExportSpec(GC,                            &
          SHORT_NAME         = 'FRACICE',                   &
          LONG_NAME          = 'fractional_cover_of_seaice',&
          UNITS              = '1',                         &
          DIMS               = MAPL_DimsTileOnly,           &
          VLOCATION          = MAPL_VLocationNone,          &
          _RC)

     if (seaIceT_extData) then
       call MAPL_AddExportSpec(GC,                            &
            SHORT_NAME         = 'SEAICETHICKNESS',           &
            LONG_NAME          = 'seaice_thickness',          &
            UNITS              = 'm',                         &
            DIMS               = MAPL_DimsTileOnly,           &
            VLOCATION          = MAPL_VLocationNone,          &
            _RC)
     endif

  elseif (DO_CICE_THERMO == 1) then  

     call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = 'TAUXIBOT',                           &
          LONG_NAME          = 'eastward_stress_at_base_of_ice',    &
          UNITS              = 'N m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          _RC)
     
     call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = 'TAUYIBOT',                           &
          LONG_NAME          = 'northward_stress_at_base_of_ice',   &
          UNITS              = 'N m-2',                             &
          DIMS               = MAPL_DimsTileOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          _RC)
  else  
     call MAPL_AddExportSpec(GC,                            &
          SHORT_NAME         = 'FRACICE',                           &
          LONG_NAME          = 'fractional_cover_of_seaice',        &
          UNITS              = '1',                                 &
          DIMS               = MAPL_DimsTileOnly,                   &
          UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
          VLOCATION          = MAPL_VLocationNone,                  &
          _RC)
  
  end if

  if (DO_CICE_THERMO == 2) then
     call MAPL_AddExportSpec ( GC   ,                          &
          SHORT_NAME = 'SURFSTATE',                            &
          CHILD_ID   = SEAICE ,                                &
                                                           _RC )
  endif
  
!EOS

! Connections between the children
!---------------------------------

#ifdef BUILD_MIT_OCEAN
  call MAPL_AddConnectivity ( GC,   &
       SHORT_NAME  = (/'ICESTATES'/), &
       DST_ID = OCEAN,               &
       SRC_ID = SEAICE,             &
       RC=STATUS  )
  VERIFY_(STATUS)
#endif

  if(DO_DATASEAONLY==0) then
!   if (trim(OCEAN_NAME) == "MOM") then  ! MOM5 only
       ! Radiation to Ocean
       call MAPL_AddConnectivity ( GC,  &
            SHORT_NAME  = (/'SWHEAT'/), &
            DST_ID = OCEAN,             &
            SRC_ID = ORAD,              &
            RC=STATUS  )
       VERIFY_(STATUS)
     
       ! Ocean to Radiation
       call MAPL_AddConnectivity ( GC,  &
            SHORT_NAME  = (/'DH   ',    &
                            'MASKO'/),  &
            DST_ID = ORAD,              &
            SRC_ID = OCEAN,             &
            RC=STATUS  )
       VERIFY_(STATUS)
!   end if
  end if

  call MAPL_AddConnectivity ( GC,  &
          SHORT_NAME  = (/'FRACICE'/), & 
          DST_ID = OCEAN,             &
          SRC_ID = SEAICE,            &
          _RC)

  if (seaIceT_extData) then
    call MAPL_AddConnectivity ( GC,  &
         SHORT_NAME  = (/'SEAICETHICKNESS'/), & 
         DST_ID = OCEAN,             &
         SRC_ID = SEAICE,            &
         _RC)
  endif

  if(DUAL_OCEAN) then 
     call MAPL_AddConnectivity ( GC,  &
          SRC_NAME  = (/'FRACICEd'/), & 
          DST_NAME  = (/'FRACICEd'/), & 
          DST_ID = OCEAN,             &
          SRC_ID = SEAICE,            &
          RC=STATUS  )
     VERIFY_(STATUS)
     ! ice nudging needs it 
     call MAPL_AddConnectivity ( GC,  &
          SHORT_NAME  = (/'SS_FOUND'/), & 
          DST_ID = SEAICE,            &
          SRC_ID = OCEAN,             &
          RC=STATUS  )
     VERIFY_(STATUS)
  endif
  
  if(DO_DATASEAONLY==0) then
     call MAPL_AddConnectivity ( GC,  &
          SHORT_NAME  = (/'UWB','VWB','UW ','VW ','SLV'/), &
          SRC_ID = OCEAN,             &
          DST_ID = SEAICE,             &
          RC=STATUS  )
     VERIFY_(STATUS)
     
     if (trim(OCEAN_NAME) == "MOM") then  ! MOM5
       call MAPL_AddConnectivity ( GC,  &
            SHORT_NAME  = (/'TAUXBOT ','TAUYBOT ', 'HICE    ', 'HSNO    ', &
                            'STROCNXB','STROCNYB', 'AICEU   ', 'FRESH   ', &
                            'FSALT   ','FHOCN   '/), &
            DST_ID = OCEAN,             &
            SRC_ID = SEAICE,            &
            RC=STATUS  )
       VERIFY_(STATUS)
     else ! MOM6
       call MAPL_AddConnectivity ( GC,  &
            SHORT_NAME  = (/'TAUXBOT ','TAUYBOT ', 'HICE    ', 'HSNO    ', &
                            'FRESH   ','FSALT   ', 'FHOCN   ', 'AICE    '/), &
            DST_ID = OCEAN,             &
            SRC_ID = SEAICE,            &
            RC=STATUS  )
       VERIFY_(STATUS)
       if (trim(OCEAN_NAME) /= "MIT") then  !
          call MAPL_AddConnectivity ( GC,   &
               SHORT_NAME  = (/'UWC','VWC'/), &
               SRC_ID = OCEAN,                &
               DST_ID = SEAICE,               &
               _RC)
       endif
    end if
  end if

  if (DO_CICE_THERMO > 1) then
     call MAPL_AddConnectivity ( GC,                          &
          SRC_NAME  = (/'TS_FOUND', 'SS_FOUND', 'FRZMLT  '/), & 
          DST_NAME  = (/'SST     ', 'SSS     ', 'FRZMLT  '/), & 
          DST_ID    = SEAICE,                                 &
          SRC_ID    = OCEAN,                                  &
          _RC )
  endif

! Children's imports are in the ocean grid and are all satisfied
! by OGCM from exchange grid quantities.
  
  if (ocean_extData) then

    if (DO_DATASEAONLY==1) then ! fake-ocean (i.e., data ocean)
      if (ocean_sssData .or. seaIceT_extData) then
        if (ocean_sssData .and. (.not. seaIceT_extData)) then ! only data_sss only
         call MAPL_TerminateImport   ( GC, ["DATA_SST  ","DATA_SSS  ", "DATA_ICEC ","DATA_KPAR "], [ocean,ocean,seaice,orad], RC=STATUS  )
        endif
        if ((.not. ocean_sssData) .and. seaIceT_extData) then ! only data_sit only
         call MAPL_TerminateImport   ( GC, ["DATA_SST  ","DATA_ICEC ", "DATA_SIT  ","DATA_KPAR "], [ocean,seaice,seaice,orad], RC=STATUS  )
        endif
        ! both data_sss and data_sit
        call MAPL_TerminateImport   ( GC, ["DATA_SST  ","DATA_SSS  ", "DATA_ICEC ", "DATA_SIT  ", "DATA_KPAR "], [ocean,ocean,seaice,seaice,orad], RC=STATUS  )
      else
        ! no data_sss and data_sit
        call MAPL_TerminateImport    ( GC, ["DATA_SST  ",             "DATA_ICEC ","DATA_KPAR "], [ocean,seaice,orad], RC=STATUS  )
      endif

    else ! we got real ocean and sea ice in case of coupled model, and only data KPAR is used.
      call MAPL_TerminateImport    ( GC, ["DATA_KPAR "], [orad], RC=STATUS  ) ! need to terminate others as well: cosz, discharge, frocean, pice, taux, tauy
    endif
  else
    call MAPL_TerminateImport(GC, ['DATA_UW', 'DATA_VW'], [OCEAN, OCEAN], _RC)
  endif

! Set the Profiling timers
! ------------------------

  call MAPL_TimerAdd(GC,    name="INITIALIZE"   ,RC=STATUS)
  VERIFY_(STATUS)
  call MAPL_TimerAdd(GC,    name="RUN"          ,RC=STATUS)
  VERIFY_(STATUS)

  call MAPL_TimerAdd(GC, name="InitChild"    ,RC=STATUS)
  VERIFY_(STATUS)

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( ogcm_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => ogcm_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'OGCM_state',wrap,status )
    VERIFY_(STATUS)

! Call SetServices 
!------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)
 
    RETURN_(ESMF_SUCCESS)
 
    contains

    subroutine OBIO_SetServices(DO_DATA_ATM4OCN, RC)
    
      logical,                intent(IN   ) ::  DO_DATA_ATM4OCN
      integer, optional,      intent(  OUT) ::  RC
      integer                               :: STATUS

      call MAPL_AddImportSpec(GC,                            &
         LONG_NAME          = 'surface_wind_speed'        ,&
         UNITS              = 'm s-1'                     ,&
         SHORT_NAME         = 'UU'                        ,&
         DIMS               = MAPL_DimsTileOnly           ,&
         VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS          ) 
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                             &
         LONG_NAME          = 'CO2 Surface Concentration Bin 001', &
         UNITS              = '1e-6'                       ,&
         SHORT_NAME         = 'CO2SC'                      ,&
         DIMS               = MAPL_DimsTileOnly            ,&
         VLOCATION          = MAPL_VLocationNone           ,&
         RESTART            = MAPL_RestartSkip             ,&
         RC=STATUS  ) 
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                             &
           LONG_NAME          = 'Dust Dry Deposition'        ,&
           UNITS              = 'kg m-2 s-1'                 ,&
           SHORT_NAME         = 'DUDP'                       ,&
           DIMS               = MAPL_DimsTileOnly            ,&
           UNGRIDDED_DIMS     = (/NUM_DUDP/)                 ,&
           VLOCATION          = MAPL_VLocationNone           ,&
           RESTART            = MAPL_RestartSkip             ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                             &
           LONG_NAME          = 'Dust Wet Deposition'        ,&
           UNITS              = 'kg m-2 s-1'                 ,&
           SHORT_NAME         = 'DUWT'                       ,&
           DIMS               = MAPL_DimsTileOnly            ,&
           UNGRIDDED_DIMS     = (/NUM_DUWT/)                 ,&
           VLOCATION          = MAPL_VLocationNone           ,&
           RESTART            = MAPL_RestartSkip             ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)
     
      call MAPL_AddImportSpec(GC,                             &
           LONG_NAME          = 'Dust Sedimentation'         ,&
           UNITS              = 'kg m-2 s-1'                 ,&
           SHORT_NAME         = 'DUSD'                       ,&
           DIMS               = MAPL_DimsTileOnly            ,&
           UNGRIDDED_DIMS     = (/NUM_DUSD/)                 ,&
           VLOCATION          = MAPL_VLocationNone           ,&
           RESTART            = MAPL_RestartSkip             ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                             &
           LONG_NAME          = 'Black Carbon Dry Deposition',&
           UNITS              = 'kg m-2 s-1'                 ,&
           SHORT_NAME         = 'BCDP'                       ,&
           DIMS               = MAPL_DimsTileOnly            ,&
           UNGRIDDED_DIMS     = (/NUM_BCDP/)                 ,&
           VLOCATION          = MAPL_VLocationNone           ,&
           RESTART            = MAPL_RestartSkip             ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)
     
      call MAPL_AddImportSpec(GC,                             &
           LONG_NAME          = 'Black Carbon Wet Deposition',&
           UNITS              = 'kg m-2 s-1'                 ,&
           SHORT_NAME         = 'BCWT'                       ,&
           DIMS               = MAPL_DimsTileOnly            ,&
           UNGRIDDED_DIMS     = (/NUM_BCWT/)                 ,&
           VLOCATION          = MAPL_VLocationNone           ,&
           RESTART            = MAPL_RestartSkip             ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)
     
      call MAPL_AddImportSpec(GC,                               &
           LONG_NAME          = 'Organic Carbon Dry Deposition',&
           UNITS              = 'kg m-2 s-1'                   ,&
           SHORT_NAME         = 'OCDP'                         ,&
           DIMS               = MAPL_DimsTileOnly              ,&
           UNGRIDDED_DIMS     = (/NUM_OCDP/)                   ,&
           VLOCATION          = MAPL_VLocationNone             ,&
           RESTART            = MAPL_RestartSkip               ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)
     
      call MAPL_AddImportSpec(GC,                               &
           LONG_NAME          = 'Organic Carbon Wet Deposition',&
           UNITS              = 'kg m-2 s-1'                   ,&
           SHORT_NAME         = 'OCWT'                         ,&
           DIMS               = MAPL_DimsTileOnly              ,&
           UNGRIDDED_DIMS     = (/NUM_OCWT/)                   ,&
           VLOCATION          = MAPL_VLocationNone             ,&
           RESTART            = MAPL_RestartSkip               ,&
           RC=STATUS  ) 
      VERIFY_(STATUS)
   
      call MAPL_AddImportSpec(GC,                               &
           SHORT_NAME         = 'DRBAND'                       ,&
           LONG_NAME          = 'surface_downwelling_shortwave_beam_flux_per_OBIO_band', &
           UNITS              = 'W m-2'                        ,&
           DIMS               = MAPL_DimsTileOnly              ,&
           UNGRIDDED_DIMS     = (/NB_OBIO/)                    ,&
           VLOCATION          = MAPL_VLocationNone             ,&
           RESTART            = MAPL_RestartSkip               ,&
           RC=STATUS  )
      VERIFY_(STATUS)

      call MAPL_AddImportSpec(GC,                               &
           SHORT_NAME         = 'DFBAND'                       ,&
           LONG_NAME          = 'surface_downwelling_shortwave_diffuse_flux_per_OBIO_band' ,&
           UNITS              = 'W m-2'                        ,&
           DIMS               = MAPL_DimsTileOnly              ,&
           UNGRIDDED_DIMS     = (/NB_OBIO/)                    ,&
           VLOCATION          = MAPL_VLocationNone             ,&
           RESTART            = MAPL_RestartSkip               ,&
           RC=STATUS  )
      VERIFY_(STATUS)

!     if (trim(OCEAN_NAME) == "MOM") then  ! MOM5 only
        ! Ocean to OceanBio
        call MAPL_AddConnectivity ( GC,   &
             SHORT_NAME  = (/'DH   ', 'T    ',       & 
                             'S    ', 'MASKO'/),     &
             DST_ID = OBIO,                          &
             SRC_ID = OCEAN,                         &
             RC=STATUS  )
        VERIFY_(STATUS)
!     end if
     
      ! OceanRad to OceanBio
      call MAPL_AddConnectivity ( GC,    &
           SHORT_NAME  = (/'TIRRQ   ',   &
                           'CDOMABSQ'/), &
           DST_ID = OBIO,                &
           SRC_ID = ORAD,                &
           RC=STATUS  )
      VERIFY_(STATUS)

      ! OceanBio to OceanRad
      call MAPL_AddConnectivity ( GC,   &
           SHORT_NAME  = (/'DIATOM','CHLORO','CYANO ','DINO  ',&
                           'PHAEO ','COCCO ','CDET  ','PIC   ',&
                           'CDC   ','AVGQ  '/), &
           DST_ID = ORAD,               &
           SRC_ID = OBIO,               &
           RC=STATUS  )
      VERIFY_(STATUS)
     
      ! Seaice to OceanBio
      call MAPL_AddConnectivity ( GC,   &
           SHORT_NAME  = (/'FRACICE'/), &
           DST_ID = OBIO,               &
           SRC_ID = SEAICE,             &
           RC=STATUS  )
      VERIFY_(STATUS)

    end subroutine OBIO_SetServices
 
  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!BOP

! !IROUTINE: Initialize -- Initialize method for the GEOS Ogcm component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the Ogcm Composite Gridded Component.
!   It reads the tiling file that defines the exchange grid 
!   It then does a Generic\_Initialize

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm 
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME
    
! Local derived type aliases

    type (MAPL_MetaComp    ), pointer   :: MAPL => null()
    type (MAPL_LocStream       )            :: EXCH
    type (ESMF_State           ), pointer   :: GIM(:) => null()
    type (ESMF_GridComp        ), pointer   :: GCS(:) => null()
    type(ESMF_FIELDBUNDLE      )            :: BUNDLE
    type(ESMF_FIELD            )            :: FIELD
    type(ESMF_Grid             )            :: grid
    integer, pointer, dimension(:)          :: TYPE => null()
    integer                                 :: I
    integer                                 :: J
    integer                                 :: N_CHILDREN
    integer                                 :: COUNTS(3)

    real, pointer                           :: FROCEAN  (:,:) => null()
    real, pointer                           :: LONS     (:  ) => null()
    real, pointer                           :: LATS     (:  ) => null()
    real, pointer                           :: TLONS    (:  ) => null()
    real, pointer                           :: TLATS    (:  ) => null()

    real, pointer                           :: PTR2d   (:,:) => null()
    logical                                 :: found
    logical, allocatable                    :: hasThisImport(:)

    integer                     :: iUseInterp
    logical                     :: UseInterp
    logical :: ACUBE, OCUBE
    integer :: NGRIDS
    integer :: A_IDX, O_IDX
    integer :: ARES, ORES
    integer :: iInterp
    integer, pointer :: GRIDIM(:)=> null()
    integer, pointer :: GRIDJM(:)=> null()
    character(len=MAPL_TileNameLength)          :: GRIDNAME
    character(len=MAPL_TileNameLength), pointer :: GNAMES(:)=> null()

    type (T_OGCM_STATE), pointer        :: ogcm_internal_state => null() 
    type (OGCM_wrap)                    :: wrap
    type(ESMF_State)                    :: SURFST

    type (ESMF_StateItem_Flag) :: itemType

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE")
    call MAPL_TimerOn(MAPL,"TOTAL"     )

    call ESMF_UserCompGetInternalState(gc, 'OGCM_state', wrap, status)
    VERIFY_(STATUS)
    ogcm_internal_state => wrap%ptr

! Get the Ocean part of the Xchg grid
!------------------------------------

    call MAPL_GenericMakeXchgNatural(MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get(MAPL, GCS = GCS, RC=STATUS )
    VERIFY_(STATUS)

! Call Initialize for every Child
!--------------------------------

    call MAPL_TimerOff(MAPL,"TOTAL"     )
    call MAPL_TimerOn(MAPL,"InitChild") 
    call MAPL_GenericInitialize ( GC, IMPORT, EXPORT, CLOCK,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerOff(MAPL,"InitChild")
    call MAPL_TimerOn (MAPL,"TOTAL"     )

! Get info from the Generic state
!--------------------------------

    call MAPL_Get(MAPL,             &
         TILETYPES = TYPE,                       &
         TILELONS  = LONS,                       &
         TILELATS  = LATS,                       &
         GIM       = GIM,                        &
                                       RC=STATUS )
    VERIFY_(STATUS)


! These are static exports
!-------------------------

    call MAPL_GetPointer(EXPORT, TLONS, 'TILELONS',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TLATS, 'TILELATS',  RC=STATUS)
    VERIFY_(STATUS)

    if(associated(TLONS)) TLONS = LONS
    if(associated(TLATS)) TLATS = LATS

! Manipulate friendly imports
!    get grid and sizes
    call ESMF_GridCompGet(GC, grid=grid, RC=status) 
    VERIFY_(STATUS)
    call MAPL_GridGet(grid, localCellCountPerDim=COUNTS, RC=STATUS)
    VERIFY_(STATUS)


! Fill the ocean fraction exposed to atmosphere (skin area) 
!   in the childrens import, if they want it
!----------------------------------------------------------

    call MAPL_Get (MAPL, ExchangeGrid=EXCH,   RC=STATUS )
    VERIFY_(STATUS)

    do I = 1, size(GIM)
       call ESMF_StateGet(GIM(I), 'FROCEAN', itemType=itemType, RC=STATUS)
       VERIFY_(STATUS)
       if (itemType == ESMF_STATEITEM_FIELD) then
          call ESMF_StateGet(GIM(I), 'FROCEAN', FIELD, RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_GetPointer(GIM(I), FROCEAN, 'FROCEAN',   RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_LocStreamFracArea( EXCH, MAPL_OCEAN, FROCEAN, RC=STATUS) 
          VERIFY_(STATUS)
       end if
    end do

    if (DO_CICE_THERMO > 1) then
        call ESMF_StateGet(EXPORT, 'SURFSTATE', SURFST, __RC__)
        call MAPL_GetPointer(SURFST, FROCEAN, 'FROCEAN', __RC__)
        call MAPL_LocStreamFracArea( EXCH, MAPL_OCEAN, FROCEAN, RC=STATUS) 
        VERIFY_(STATUS)
    endif 

! Put OBIO tracers into the OCEAN's tracer bundle.
!-------------------------------------------------

    if (DO_DATASEAONLY==0) then
       if (trim(OCEAN_NAME) == "MOM") then
         call ESMF_StateGet(GIM(OCEAN), 'TR', BUNDLE, RC=STATUS)
         VERIFY_(STATUS)
         if (DO_OBIO/=0) then
            call MAPL_GridCompGetFriendlies(GCS(OBIO),"OCEAN", BUNDLE, RC=STATUS )
            VERIFY_(STATUS)
         end if
      endif
    end if

!   The section below attempts to make an intellegent guess of the default
!   for INTERPOLATE_SST 

!   Get the name of the ocean grid
    call ESMF_GridGet(GRID, name=gridname, rc=status)
    VERIFY_(STATUS)
!   first query the exchange grid for the names of the 2 grids (ATM and OCN)

    call MAPL_LocStreamGet(Exch, GRIDNAMES = GNAMES, &
         GRIDIM=GRIDIM, GRIDJM=GRIDJM, RC=STATUS)
    VERIFY_(STATUS)
!   query exchange grid for ngrids
    ngrids = size(gnames)
    _ASSERT(ngrids==2,'needs informative message')

    
!   validate that gridname is there
    found = .false.
    DO I = 1, NGRIDS
       IF (GNAMES(I) == GRIDNAME) THEN
          FOUND = .TRUE.
          exit
       ENDIF
    ENDDO
    _ASSERT(FOUND,'needs informative message')

    O_IDX = I
    A_IDX = 3-I
! we pick the "other" gridname (i.e. ATM). 
! this logic works only when ngrids==2; 3-1=2;3-2=1

! Check if any of the grids is on a cubed-sphere
    OCUBE = (  GRIDIM(O_IDX)==6*GRIDJM(O_IDX)) ! MIT Ocean uses "different" cubed-sphere (it is fatter in the IM direction)
    ACUBE = (6*GRIDIM(A_IDX)==  GRIDJM(A_IDX))
    ARES = GRIDIM(A_IDX)
    if( ACUBE ) ARES = 4*ARES
    ORES = GRIDIM(O_IDX)
    if( OCUBE ) then
       ORES = 4*(ORES/6)
    end if

    if (ARES > ORES) then
       ! Ocean grid is coarser; we should interpolate
       iInterp = 1
    else
       ! No interpolation needed
       iInterp = 0
    end if

    call MAPL_GetResource(MAPL, iUseInterp, 'INTERPOLATE_SST:', &
         default=iInterp, RC=STATUS )
    VERIFY_(STATUS)
    useInterp = (iUseInterp /= 0)
    ogcm_internal_state%useInterp = useInterp

! All Done
!---------

    call MAPL_TimerOff(MAPL,"TOTAL"     )
    call MAPL_TimerOff(MAPL,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!BOP

! !IROUTINE: RUN -- Run method for the Ogcm component

! !INTERFACE:

  subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )

  use ice_constants,   only: puny, Tocnfrz

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Periodically refreshes the ozone mixing ratios.

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

    type (MAPL_MetaComp),      pointer  :: MAPL   => null()
    type (ESMF_GridComp),      pointer  :: GCS(:) => null()
    type (ESMF_State),         pointer  :: GIM(:) => null()
    type (ESMF_State),         pointer  :: GEX(:) => null()
    type (MAPL_LocStream)               :: EXCHGrid
    logical                             :: FRIENDLY
    type(ESMF_FIELD)                    :: FIELD

! Pointers to imports (All tiled)

    real, pointer, dimension(:) :: TAUXW => null()
    real, pointer, dimension(:) :: TAUYW => null()
    real, pointer, dimension(:) :: TAUXI => null()
    real, pointer, dimension(:) :: TAUYI => null()
    real, pointer, dimension(:) :: USTR3 => null()
    real, pointer, dimension(:) :: UU => null()
    real, pointer, dimension(:) :: PS => null()
    real, pointer, dimension(:) :: PENUVR => null()
    real, pointer, dimension(:) :: PENUVF => null()
    real, pointer, dimension(:) :: PENPAR => null()
    real, pointer, dimension(:) :: PENPAF => null()
    real, pointer, dimension(:) :: DRNIR => null()
    real, pointer, dimension(:) :: DFNIR => null()
    real, pointer, dimension(:) :: HI => null()
    real, pointer, dimension(:) :: SI => null()
    real, pointer, dimension(:) :: DISCHARGE => null() 
    real, pointer, dimension(:) :: CO2SC => null()

    real, pointer, dimension(:,:) :: DUDP => null()
    real, pointer, dimension(:,:) :: DUWT => null()
    real, pointer, dimension(:,:) :: DUSD => null()
    real, pointer, dimension(:,:) :: BCDP => null()
    real, pointer, dimension(:,:) :: BCWT => null()
    real, pointer, dimension(:,:) :: OCDP => null()
    real, pointer, dimension(:,:) :: OCWT => null()
    real, pointer, dimension(:,:) :: DRBAND => null()          
    real, pointer, dimension(:,:) :: DFBAND => null()
    real, pointer, dimension(:)   :: TI => null()
    real, pointer, dimension(:)   :: FR => null()
    real, pointer, dimension(:,:) :: TI8 => null()
    real, pointer, dimension(:,:) :: FR8 => null()
    real, pointer, dimension(:,:) :: VOLICE => null()
    real, pointer, dimension(:,:) :: VOLSNO => null()
    real, pointer, dimension(:,:) :: TAUAGE => null()
    real, pointer, dimension(:,:) :: MPOND => null()

    real, pointer, dimension(:,:,:) :: ERGICE => null()
    real, pointer, dimension(:,:,:) :: ERGSNO => null()

    real, pointer, dimension(:) :: TAUXIBOT => null() 
    real, pointer, dimension(:) :: TAUYIBOT => null()

    real, pointer, dimension(:) :: LWFLX => null()
    real, pointer, dimension(:) :: SHFLX => null()
    real, pointer, dimension(:) :: QFLUX => null()
    real, pointer, dimension(:) :: SNOW => null()
    real, pointer, dimension(:) :: RAIN => null()
    real, pointer, dimension(:) :: FRESH => null()
    real, pointer, dimension(:) :: FSALT => null()
    real, pointer, dimension(:) :: FHOCN => null()
    real, pointer, dimension(:) :: PEN_OCN => null()

! Pointers to ocn grid versions

    real, pointer, dimension(:,:) :: TAUXIO => null()
    real, pointer, dimension(:,:) :: TAUYIO => null()
    real, pointer, dimension(:,:) :: TAUXWO => null()
    real, pointer, dimension(:,:) :: TAUYWO => null()
    real, pointer, dimension(:,:) :: USTR3O => null()
    real, pointer, dimension(:,:) :: PSO    => null()
    real, pointer, dimension(:,:) :: USTR3B => null()
    real, pointer, dimension(:,:) :: UUB    => null()
    real, pointer, dimension(:,:) :: UUR    => null()
    real, pointer, dimension(:,:) :: PSB    => null()
    real, pointer, dimension(:,:) :: PSR    => null()
    real, pointer, dimension(:,:) :: CO2SCB => null()

    real, pointer, dimension(:,:,:) :: DUDPB => null()
    real, pointer, dimension(:,:,:) :: DUWTB => null()
    real, pointer, dimension(:,:,:) :: DUSDB => null()
    real, pointer, dimension(:,:,:) :: BCDPB => null()
    real, pointer, dimension(:,:,:) :: BCWTB => null()
    real, pointer, dimension(:,:,:) :: OCDPB => null()
    real, pointer, dimension(:,:,:) :: OCWTB => null()
    real, pointer, dimension(:,:,:) :: DRBANDR => null()
    real, pointer, dimension(:,:,:) :: DFBANDR => null()
    real, pointer, dimension(:,:) :: PENUVRO => null()
    real, pointer, dimension(:,:) :: PENUVFO => null()
    real, pointer, dimension(:,:) :: PENPARO => null()
    real, pointer, dimension(:,:) :: PENPAFO => null()
    real, pointer, dimension(:,:) :: DRNIRO  => null()
    real, pointer, dimension(:,:) :: DFNIRO  => null()
    real, pointer, dimension(:,:) :: DISCHARGEOB => null()

    real, pointer, dimension(:,:) :: PENUVRM    => null()
    real, pointer, dimension(:,:) :: PENUVFM    => null()
    real, pointer, dimension(:,:) :: PENPARM    => null()
    real, pointer, dimension(:,:) :: PENPAFM    => null()
    real, pointer, dimension(:,:) :: DRNIRM    => null()
    real, pointer, dimension(:,:) :: DFNIRM    => null()    
    real, pointer, dimension(:,:) :: DISCHARGEO => null()
    real, pointer, dimension(:,:) :: UWO => null()
    real, pointer, dimension(:,:) :: VWO => null()

    real, pointer, dimension(:,:) :: HIO => null()
    real, pointer, dimension(:,:) :: SIO => null()
    real, pointer, dimension(:,:) :: UIO => null()
    real, pointer, dimension(:,:) :: SEAICETHICKNESSO  => null()
    real, pointer, dimension(:,:) :: VIO => null()
    real, pointer, dimension(:,:) :: KPARO => null()
    real, pointer, dimension(:,:) :: TS_FOUNDO => null()
    real, pointer, dimension(:,:) :: SS_FOUNDO => null()
    real, pointer, dimension(:,:) :: FRZMLTO   => null()
    real, pointer, dimension(:,:) :: TAUXIBOTO => null() 
    real, pointer, dimension(:,:) :: TAUYIBOTO => null()
    real, pointer, dimension(:,:) :: STROCNXB => null() 
    real, pointer, dimension(:,:) :: STROCNYB => null()
    real, pointer, dimension(:,:) :: AICEU => null()
    real, pointer, dimension(:,:) :: AICE  => null()
    real, pointer, dimension(:,:) :: UWBO => null()
    real, pointer, dimension(:,:) :: VWBO => null()

    real, pointer, dimension(:,:,:) :: TIO8 => null()
    real, pointer, dimension(:,:)   :: FRI => null()
    real, pointer, dimension(:,:,:) :: FRO8 => null()
    real, pointer, dimension(:,:)   :: FRO => null()
    real, pointer, dimension(:,:)   :: TIO => null()
    real, pointer, dimension(:,:,:) :: VOLICEO => null()
    real, pointer, dimension(:,:,:) :: VOLSNOO => null()
    real, pointer, dimension(:,:,:) :: TAUAGEO => null()
    real, pointer, dimension(:,:,:) :: MPONDO => null()
    real, pointer, dimension(:,:,:) :: ERGICEO => null()
    real, pointer, dimension(:,:,:) :: ERGSNOO => null()

    real, pointer, dimension(:,:) :: LWFLXO => null()
    real, pointer, dimension(:,:) :: SHFLXO => null()
    real, pointer, dimension(:,:) :: QFLUXO => null()
    real, pointer, dimension(:,:) :: SNOWO => null()
    real, pointer, dimension(:,:) :: RAINO => null()
    real, pointer, dimension(:,:) :: FRESHO   => null()
    real, pointer, dimension(:,:) :: FSALTO   => null()
    real, pointer, dimension(:,:) :: FHOCNO   => null()
    real, pointer, dimension(:,:) :: PEN_OCNO => null()

! Pointers to exports

    real, pointer, dimension(:) :: UW       => null()
    real, pointer, dimension(:) :: VW       => null()
    real, pointer, dimension(:) :: UI       => null()
    real, pointer, dimension(:) :: VI       => null()
    real, pointer, dimension(:) :: KPAR     => null()
    real, pointer, dimension(:) :: TS_FOUND => null()
    real, pointer, dimension(:) :: SS_FOUND => null()
    real, pointer, dimension(:) :: FRZMLT   => null()
    real, pointer, dimension(:) :: SEAICETHICKNESS   => null()

    real, allocatable, dimension(:) :: VARTILE

    integer                     :: N, K
    integer                     :: iUseInterp
    logical                     :: UseInterp

    type (T_OGCM_STATE), pointer        :: ogcm_internal_state => null() 
    type (OGCM_wrap)                    :: wrap

    integer :: ID
    integer :: PHASE
    integer :: PHASE_
    integer, allocatable :: CHLD(:)
    type(ESMF_Info) :: infoh

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = 'Run'
    call ESMF_GridCompGet( GC, name=COMP_NAME, currentPhase=PHASE, RC=status)
    VERIFY_(STATUS)
    if (PHASE >= 10) PHASE = PHASE - 10 ! to be replaced with MAPL get_phase   
    Iam = trim(COMP_NAME) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"RUN"  )
    call MAPL_TimerOn(MAPL,"TOTAL")

! Get parameters from generic state.
!-----------------------------------

    call MAPL_Get(MAPL,             &
         ExchangeGrid  = ExchGrid,  &
         GIM       = GIM,           & 
         GEX       = GEX,           &
         GCS       = GCS,           &
    RC=STATUS )
    VERIFY_(STATUS)

! Pointers to imports
!--------------------

    call MAPL_GetPointer(IMPORT, TAUXW   ,  'TAUXW'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUYW   ,  'TAUYW'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUXI   ,  'TAUXI'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUYI   ,  'TAUYI'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, USTR3   ,  'OUSTAR3', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PS      ,  'PS'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENUVR  ,  'PENUVR' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENUVF  ,  'PENUVF' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENPAR  ,  'PENPAR' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENPAF  ,  'PENPAF' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DRNIR, 'DRNIR', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DFNIR, 'DFNIR', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DISCHARGE, 'DISCHRG', RC=STATUS)
    VERIFY_(STATUS)
    if (.not. seaIceT_extData) then
      if (DO_CICE_THERMO <= 1) then  
          call MAPL_GetPointer(IMPORT, HI      ,  'HI'     , _RC)
      endif
      call MAPL_GetPointer(IMPORT, SI      ,  'SI'     , _RC)
    endif

    if (DO_CICE_THERMO /= 0) then  
       call MAPL_GetPointer(IMPORT, TI8     ,  'TI'     , _RC)
       VERIFY_(STATUS)
    else
       if (.not. seaIceT_extData) then
         call MAPL_GetPointer(IMPORT, TI      ,  'TI'     , _RC)
       endif 
    endif
    
    if (DO_CICE_THERMO == 1) then  
       call MAPL_GetPointer(IMPORT, FR8     , 'FRACICE' , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, VOLICE  , 'VOLICE'  , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, VOLSNO  , 'VOLSNO'  , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, ERGICE  , 'ERGICE'  , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, ERGSNO  , 'ERGSNO'  , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, TAUAGE  , 'TAUAGE'  , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(IMPORT, MPOND   , 'MPOND'   , RC=STATUS)
       VERIFY_(STATUS)
    endif 
    
    if (DO_OBIO/=0) then
      call OBIO_RunTransforms(DO_DATA_ATM4OCN, RC)
    endif

    call MAPL_GetPointer(IMPORT, LWFLX, 'LWFLX', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SHFLX, 'SHFLX', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QFLUX, 'QFLUX', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SNOW,  'SNOW' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RAIN,  'RAIN' , RC=STATUS)
    VERIFY_(STATUS)
    if (DO_CICE_THERMO <= 1) then  
       call MAPL_GetPointer(IMPORT, FRESH, 'FRESH', _RC)
       call MAPL_GetPointer(IMPORT, FSALT, 'FSALT', _RC)
       call MAPL_GetPointer(IMPORT, FHOCN, 'FHOCN', _RC)
    endif 
    call MAPL_GetPointer(IMPORT, PEN_OCN,'PEN_OCN',RC=STATUS)
    VERIFY_(STATUS)

! Verify that the saltwater ice variables are friendly to seaice
!---------------------------------------------------------------

    if (.not. seaIceT_extData) then
       call ESMF_StateGet (IMPORT, 'TI', FIELD, _RC)
       call ESMF_InfoGetFromHost(FIELD,infoh,_RC)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,_RC)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'SI', FIELD, _RC)
       call ESMF_InfoGetFromHost(FIELD,infoh,_RC)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,_RC)
       _ASSERT(FRIENDLY,'needs informative message')

      if(DO_CICE_THERMO <= 1) then
         call ESMF_StateGet (IMPORT, 'HI', FIELD, _RC)
         call ESMF_InfoGetFromHost(FIELD,infoh,_RC)
         call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,_RC)
         _ASSERT(FRIENDLY,'needs informative message')
      endif
    endif 

    if(DO_CICE_THERMO==1) then
       call ESMF_StateGet (IMPORT, 'FRACICE', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'VOLICE', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'VOLSNO', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'ERGICE', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'ERGSNO', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'TAUAGE', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')

       call ESMF_StateGet (IMPORT, 'MPOND', FIELD, RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_InfoGetFromHost(FIELD,infoh,RC=STATUS)
       call ESMF_InfoGet(infoh,'FriendlyToSEAICE',FRIENDLY,RC=STATUS)
       VERIFY_(STATUS)
       _ASSERT(FRIENDLY,'needs informative message')
    end if
    
! Children's Imports
!-------------------

    call MAPL_GetPointer(GIM(SEAICE), TAUXIO  ,  'TAUX', notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(SEAICE), TAUYIO  ,  'TAUY', notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(GIM(OCEAN ), TAUXWO  ,  'TAUX',  notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(OCEAN ), TAUYWO  ,  'TAUY',  notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(OCEAN ), USTR3O  ,  'OUSTAR3', notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(OCEAN ), PSO     ,  'PS'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)

    if(DO_DATASEAONLY==0) then
       call MAPL_GetPointer(GIM(OCEAN  ), PENUVRM ,  'PENUVR',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(OCEAN  ), PENUVFM ,  'PENUVF',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(OCEAN  ), PENPARM ,  'PENPAR',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(OCEAN  ), PENPAFM ,  'PENPAF',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(OCEAN ), DRNIRM   ,  'DRNIR' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(OCEAN ), DFNIRM   ,  'DFNIR' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(OCEAN ), DISCHARGEO, 'DISCHARGE',  RC=STATUS)
       VERIFY_(STATUS)
    end if
    
    call MAPL_GetPointer(GIM(ORAD  ), PENUVRO ,  'PENUVR',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(ORAD  ), PENUVFO ,  'PENUVF',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(ORAD  ), PENPARO ,  'PENPAR',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(ORAD  ), PENPAFO ,  'PENPAF',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(ORAD  ), DRNIRO ,  'DRNIR',  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(GIM(ORAD  ), DFNIRO ,  'DFNIR',  RC=STATUS)
    VERIFY_(STATUS)

    if (.not. seaIceT_extData) then
      if (DO_CICE_THERMO <= 1) then  
         call MAPL_GetPointer(GIM(SEAICE), HIO     ,  'HI'    ,  _RC)
      endif
      call MAPL_GetPointer(GIM(SEAICE), SIO     ,  'SI'    ,  _RC)
    endif
    if (DO_CICE_THERMO == 0) then  
       if (.not. seaIceT_extData) then
         call MAPL_GetPointer(GIM(SEAICE), TIO   ,  'TI'    ,  _RC)
       endif
    elseif(DO_CICE_THERMO == 1) then
       call MAPL_GetPointer(GIM(SEAICE), TIO8    ,  'TI'    ,  _RC)
    else
       call MAPL_GetPointer(GEX(SEAICE), TIO8    ,  'TI'    ,  _RC)
    endif
    if (DO_CICE_THERMO == 1) then  
       call MAPL_GetPointer(GIM(SEAICE), FRO8    , 'FRACICE',  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), VOLICEO , 'VOLICE' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), VOLSNOO , 'VOLSNO' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), ERGICEO , 'ERGICE' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), ERGSNOO , 'ERGSNO' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), TAUAGEO , 'TAUAGE' ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), MPONDO  , 'MPOND'  ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), FRESHO  , 'FRESH'  ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), FSALTO  , 'FSALT'  ,  RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_GetPointer(GIM(SEAICE), FHOCNO  , 'FHOCN'  ,  RC=STATUS)
       VERIFY_(STATUS)
    endif

   call MAPL_GetPointer(GIM(OCEAN), LWFLXO, 'LWFLX',  RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(GIM(OCEAN), SHFLXO, 'SHFLX',  RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(GIM(OCEAN), QFLUXO, 'QFLUX',  RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(GIM(OCEAN), SNOWO, 'SNOW',  RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(GIM(OCEAN), RAINO, 'RAIN',  RC=STATUS)
   VERIFY_(STATUS)
   call MAPL_GetPointer(GIM(OCEAN), PEN_OCNO,'PEN_OCN',RC=STATUS)
   VERIFY_(STATUS)

! Transform imports to the ocean grid
!------------------------------------

    if(associated(TAUXIO)) then
       call MAPL_LocStreamTransform( ExchGrid, TAUXIO  ,  TAUXI  , RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated(TAUYIO)) then
       call MAPL_LocStreamTransform( ExchGrid, TAUYIO  ,  TAUYI  , RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated(TAUXWO)) then
       call MAPL_LocStreamTransform( ExchGrid, TAUXWO  ,  TAUXW  , RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated(TAUYWO)) then
       call MAPL_LocStreamTransform( ExchGrid, TAUYWO  ,  TAUYW  , RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated(USTR3O)) then
       call MAPL_LocStreamTransform( ExchGrid, USTR3O  ,  USTR3  , RC=STATUS) 
       VERIFY_(STATUS)
    endif
    if(associated(PSO)) then
       call MAPL_LocStreamTransform( ExchGrid, PSO     ,  PS     , RC=STATUS) 
       VERIFY_(STATUS)
    endif
    
    call MAPL_LocStreamTransform( ExchGrid, PENUVRO,  PENUVR, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, PENUVFO,  PENUVF, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, PENPARO,  PENPAR, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, PENPAFO,  PENPAF, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, DRNIRO,   DRNIR,  RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, DFNIRO,   DFNIR,  RC=STATUS) 
    VERIFY_(STATUS)

    if(DO_DATASEAONLY==0) then
       call MAPL_LocStreamTransform( ExchGrid, DISCHARGEO, DISCHARGE, RC=STATUS) 
       VERIFY_(STATUS)
     
       PENUVRM= PENUVRO
       PENUVFM= PENUVFO
       PENPARM= PENPARO 
       PENPAFM= PENPAFO
       DRNIRM= DRNIRO
       DFNIRM= DFNIRO
    end if

    if ( associated(DISCHARGEOB) ) then
       call MAPL_LocStreamTransform( ExchGrid, DISCHARGEOB, DISCHARGE, RC=STATUS)
       VERIFY_(STATUS)
    end if
    
    if (.not. seaIceT_extData) then
      call MAPL_LocStreamTransform( ExchGrid, SIO    ,  SI    , _RC)
      if (DO_CICE_THERMO <= 1) then  
         call MAPL_LocStreamTransform( ExchGrid, HIO    ,  HI    , _RC)
      endif
    endif
    
    if (DO_CICE_THERMO == 0) then  
      if (.not. seaIceT_extData) then
        call MAPL_LocStreamTransform( ExchGrid, TIO  ,  TI    , _RC)
      endif
    elseif (DO_CICE_THERMO == 1) then
       allocate(VARTILE(size(TI8,dim=1)), STAT=STATUS) 
       VERIFY_(STATUS)
       do n=1,NUM_ICE_CATEGORIES 
! When running dual ocean, sea ice data needs the import also
          VARTILE = TI8(:,N) * FR8(:,N)
          call MAPL_LocStreamTransform( ExchGrid, TIO8(:,:,N),  VARTILE, RC=STATUS) 
          VERIFY_(STATUS)
          call MAPL_LocStreamTransform( ExchGrid, FRO8(:,:,N),  FR8(:,N), RC=STATUS) 
          VERIFY_(STATUS)
          where(FRO8(:,:,N) /= MAPL_UNDEF .and. FRO8(:,:,N) < puny)
            TIO8(:,:,N) = MAPL_TICE+Tocnfrz 
          endwhere 
          where(FRO8(:,:,N) /= MAPL_UNDEF .and. FRO8(:,:,N) >= puny)
            TIO8(:,:,N) = TIO8(:,:,N) / FRO8(:,:,N) 
          endwhere 
          call MAPL_LocStreamTransform( ExchGrid, VOLICEO(:,:,N),  &
               VOLICE (:,  N),  & 
               RC=STATUS) 
          VERIFY_(STATUS)
          call MAPL_LocStreamTransform( ExchGrid, VOLSNOO(:,:,N),  &
               VOLSNO (:,  N), RC=STATUS) 
          VERIFY_(STATUS)
          VARTILE = TAUAGE(:,N) * FR8(:,N) 
          call MAPL_LocStreamTransform( ExchGrid, TAUAGEO(:,:,N),  &
               VARTILE, RC=STATUS) 
          VERIFY_(STATUS)
          where(FRO8(:,:,N) /= MAPL_UNDEF .and. FRO8(:,:,N) < puny)
            TAUAGEO(:,:,N) = 0.0 
          endwhere 
          where(FRO8(:,:,N) /= MAPL_UNDEF .and. FRO8(:,:,N) >= puny)
            TAUAGEO(:,:,N) = TAUAGEO(:,:,N) / FRO8(:,:,N) 
          endwhere 
          call MAPL_LocStreamTransform( ExchGrid, MPONDO(:,:,N),  &
               MPOND (:,  N), RC=STATUS) 
          VERIFY_(STATUS)
          do k=1,NUM_ICE_LAYERS 
             call MAPL_LocStreamTransform( ExchGrid,                             &  
                  ERGICEO(:,:,NUM_ICE_LAYERS*(N-1)+K),  &
                  ERGICE (:,  K,N), RC=STATUS) 
             VERIFY_(STATUS)
          enddo
          do k=1,NUM_SNOW_LAYERS 
             call MAPL_LocStreamTransform( ExchGrid,                              &
                  ERGSNOO(:,:,NUM_SNOW_LAYERS*(N-1)+K),  &
                  ERGSNO (:,  K,N), RC=STATUS) 
             VERIFY_(STATUS)
          enddo
       enddo
       deallocate(VARTILE, STAT=STATUS)
       VERIFY_(STATUS)
       if (.not. (dual_ocean .and. phase==2) ) then
! for efficiency, dont need these if running predictor in dual ocean
        call MAPL_LocStreamTransform( ExchGrid, FRESHO,    FRESH,   RC=STATUS) 
        VERIFY_(STATUS)
        call MAPL_LocStreamTransform( ExchGrid, FSALTO,    FSALT,   RC=STATUS) 
        VERIFY_(STATUS)
        call MAPL_LocStreamTransform( ExchGrid, FHOCNO,    FHOCN,   RC=STATUS) 
        VERIFY_(STATUS)
       endif
    endif

    call MAPL_LocStreamTransform( ExchGrid, LWFLXO,  LWFLX, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, SHFLXO,  SHFLX, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, QFLUXO,  QFLUX, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, SNOWO,  SNOW, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, RAINO,  RAIN, RC=STATUS) 
    VERIFY_(STATUS)
    call MAPL_LocStreamTransform( ExchGrid, PEN_OCNO,  PEN_OCN, RC=STATUS) 
    VERIFY_(STATUS)

! Pointers to tile outputs
!-------------------------
    if (DO_CICE_THERMO == 0) then  
       call MAPL_GetPointer(EXPORT, FR      ,  'FRACICE', _RC)
       if (seaIceT_extData) then
         call MAPL_GetPointer(EXPORT, SEAICETHICKNESS,  'SEAICETHICKNESS', _RC)
       endif
    elseif (DO_CICE_THERMO == 1) then  
       call MAPL_GetPointer(EXPORT, TAUXIBOT,  'TAUXIBOT', _RC)
       call MAPL_GetPointer(EXPORT, TAUYIBOT,  'TAUYIBOT', _RC)
    else  
       call MAPL_GetPointer(EXPORT, FR8     ,  'FRACICE',  _RC)
    end if

    call MAPL_GetPointer(EXPORT, UW      ,  'UW'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VW      ,  'VW'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UI      ,  'UI'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VI      ,  'VI'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, KPAR    ,  'KPAR'   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TS_FOUND,  'TS_FOUND' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SS_FOUND,  'SS_FOUND' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZMLT  ,  'FRZMLT'   , RC=STATUS)
    VERIFY_(STATUS)

! Mark as needed in the children by allocating
!---------------------------------------------

    if(associated(UW)) then
       call MAPL_GetPointer(GEX(OCEAN ), UWO ,  'UW'    , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(VW)) then
       call MAPL_GetPointer(GEX(OCEAN ), VWO ,  'VW'    , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(UI)) then
       call MAPL_GetPointer(GEX(SEAICE), UIO ,  'UI'    , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(VI)) then
       call MAPL_GetPointer(GEX(SEAICE), VIO ,  'VI'    , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(KPAR)) then
       call MAPL_GetPointer(GEX(ORAD ), KPARO ,  'KPAR' , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(TS_FOUND)) then
       call MAPL_GetPointer(GEX(OCEAN ), TS_FOUNDO ,  'TS_FOUND' , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(SS_FOUND)) then
       call MAPL_GetPointer(GEX(OCEAN ), SS_FOUNDO ,  'SS_FOUND' , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(FRZMLT)) then
       call MAPL_GetPointer(GEX(OCEAN ), FRZMLTO ,  'FRZMLT' , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(DO_DATASEAONLY==0) then
       call MAPL_GetPointer(GEX(OCEAN ), UWBO ,  'UWB'    , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
       
       call MAPL_GetPointer(GEX(OCEAN ), VWBO ,  'VWB'    , alloc=.true., RC=STATUS)
       VERIFY_(STATUS)
    end if

    if (DO_CICE_THERMO == 0) then  
       call MAPL_GetPointer(GEX(SEAICE), FRO  ,  'FRACICE', alloc=.true., _RC)
       if (seaIceT_extData) then
         call MAPL_GetPointer(GEX(SEAICE), SEAICETHICKNESSO, 'SEAICETHICKNESS', alloc=.true., _RC)
       endif
    elseif (DO_CICE_THERMO == 1) then  
       call MAPL_GetPointer(GEX(SEAICE), FRI  ,  'FRACICE', alloc=.true., _RC)

       if(associated(TAUXIBOT)) then
          call MAPL_GetPointer(GEX(SEAICE), TAUXIBOTO , 'TAUXBOT' , alloc=.true., _RC)
       end if
       
       if(associated(TAUYIBOT)) then
          call MAPL_GetPointer(GEX(SEAICE), TAUYIBOTO , 'TAUYBOT' , alloc=.true., _RC)
       end if
    else
       call MAPL_GetPointer(GEX(SEAICE), FRO8  ,  'FRSEAICE', alloc=.true.,    _RC)
    endif

    call MAPL_TimerOff(MAPL,"TOTAL"     )

    if (.not. DUAL_OCEAN) then
       call MAPL_GenericRunChildren(GC, IMPORT, EXPORT, CLOCK, RC=STATUS)
       VERIFY_(STATUS)
    else
       if (PHASE == 1) then
          ! corrector
          ! run explicitly phase 1 of all the children
          allocate(CHLD(4), stat=status)
          VERIFY_(STATUS)
          CHLD = (/OBIO,ORAD,SEAICE,OCEAN/)
          DO N=1, size(CHLD)
             ID = CHLD(N)
             if (ID <= 0) cycle
             call ESMF_GridCompRun( GCS(ID), importState=GIM(ID), &
                  exportState=GEX(ID), clock=CLOCK, phase=1, userRC=STATUS )
             VERIFY_(STATUS)
             call MAPL_GenericRunCouplers( MAPL, CHILD=ID, CLOCK=CLOCK, RC=STATUS )
             VERIFY_(STATUS)
          END DO
          deallocate(CHLD)

       else
       ! run explicitly the children excluding "real" seaice (ocean has the data part inside ocean)
          allocate(CHLD(4), stat=status)
          VERIFY_(STATUS)
          CHLD = (/OBIO,ORAD,SEAICE,OCEAN/)
          DO N=1, size(CHLD)
             ID = CHLD(N)
             if (ID <= 0) cycle
             if (ID /= OCEAN .and. ID /= SEAICE) then
                phase_ = 1
             else
                phase_ = phase
             end if
             call ESMF_GridCompRun( GCS(ID), importState=GIM(ID), &
                  exportState=GEX(ID), clock=CLOCK, phase=phase_, userRC=STATUS )
             VERIFY_(STATUS)
             call MAPL_GenericRunCouplers( MAPL, CHILD=ID, CLOCK=CLOCK, RC=STATUS )
             VERIFY_(STATUS)
          END DO
          deallocate(CHLD)

       end if
    end if

    call MAPL_TimerOn (MAPL,"TOTAL"     )

    call ESMF_UserCompGetInternalState(gc, 'OGCM_state', wrap, status)
    VERIFY_(STATUS)
    ogcm_internal_state => wrap%ptr

    useInterp = ogcm_internal_state%useInterp

    if (.not. seaIceT_extData) then
      call MAPL_LocStreamTransform( ExchGrid, SI     ,  SIO   , _RC)
      if (DO_CICE_THERMO <= 1) then  
         call MAPL_LocStreamTransform( ExchGrid, HI     ,  HIO   , _RC)
      endif
    endif

    ! call Run2 of SEAICE to do ice nudging 
    if (dual_ocean) then
        if(PHASE==1) then
             ! phase 3 is corrector stage same as phase 1
             call ESMF_GridCompRun( GCS(SEAICE), importState=GIM(SEAICE), &
                  exportState=GEX(SEAICE), clock=CLOCK, phase=3, userRC=STATUS )
             VERIFY_(STATUS)
             call MAPL_GenericRunCouplers( MAPL, CHILD=SEAICE, CLOCK=CLOCK, RC=STATUS )
             VERIFY_(STATUS)
        else
             ! phase 4 predictor stage same as phase 2
             call ESMF_GridCompRun( GCS(SEAICE), importState=GIM(SEAICE), &
                  exportState=GEX(SEAICE), clock=CLOCK, phase=4, userRC=STATUS )
             VERIFY_(STATUS)
             call MAPL_GenericRunCouplers( MAPL, CHILD=SEAICE, CLOCK=CLOCK, RC=STATUS )
             VERIFY_(STATUS)
        endif
    endif

    if (DO_CICE_THERMO == 0) then  
       call MAPL_LocStreamTransform( ExchGrid, FR     ,  FRO   , INTERP=useInterp, _RC)

       if (.not. seaIceT_extData) then
         call MAPL_LocStreamTransform( ExchGrid, TI     ,  TIO   , _RC)
       else
         call MAPL_LocStreamTransform( ExchGrid, SEAICETHICKNESS,  SEAICETHICKNESSO, INTERP=useInterp, _RC)
       endif
    elseif (DO_CICE_THERMO == 1) then
       do n=1,NUM_ICE_CATEGORIES
           call MAPL_LocStreamTransform( ExchGrid, TI8(:,N),  TIO8(:,:,N), RC=STATUS) 
           VERIFY_(STATUS)
           if (DO_CICE_THERMO == 1) then
           call MAPL_LocStreamTransform( ExchGrid, FR8(:,N),  FRO8(:,:,N),  & 
                INTERP=useInterp, RC=STATUS) 
           VERIFY_(STATUS)
           call MAPL_LocStreamTransform( ExchGrid, VOLICE (:,  N), &
               VOLICEO(:,:,N), RC=STATUS) 
           VERIFY_(STATUS)
           call MAPL_LocStreamTransform( ExchGrid, VOLSNO (:,  N), &
               VOLSNOO(:,:,N), RC=STATUS) 
          VERIFY_(STATUS)
          call MAPL_LocStreamTransform( ExchGrid, TAUAGE (:,  N), &
               TAUAGEO(:,:,N), RC=STATUS) 
          VERIFY_(STATUS)
          call MAPL_LocStreamTransform( ExchGrid, MPOND (:,  N), &
               MPONDO(:,:,N), RC=STATUS) 
           VERIFY_(STATUS)
          do k=1,NUM_ICE_LAYERS 
             call MAPL_LocStreamTransform( ExchGrid, ERGICE (:,  K,N),  &
                  ERGICEO(:,:,NUM_ICE_LAYERS*(N-1)+K),RC=STATUS) 
             VERIFY_(STATUS)
          enddo
          do k=1,NUM_SNOW_LAYERS 
             call MAPL_LocStreamTransform( ExchGrid, ERGSNO (:,  K,N),  &
                  ERGSNOO(:,:,NUM_SNOW_LAYERS*(N-1)+K),RC=STATUS) 
             VERIFY_(STATUS)
          enddo
          endif
       enddo
       
       if(associated(TAUXIBOT)) then
          call MAPL_LocStreamTransform( ExchGrid, TAUXIBOT, TAUXIBOTO, RC=STATUS) 
          VERIFY_(STATUS)
       end if
       
       if(associated(TAUYIBOT)) then
          call MAPL_LocStreamTransform( ExchGrid, TAUYIBOT, TAUYIBOTO, RC=STATUS) 
          VERIFY_(STATUS)
       end if
    else
       do n=1,NUM_ICE_CATEGORIES
          call MAPL_LocStreamTransform( ExchGrid, TI8(:,N),  TIO8(:,:,N),   _RC)
          call MAPL_LocStreamTransform( ExchGrid, FR8(:,N),  FRO8(:,:,N),   _RC)
       enddo
    endif

    call MAPL_GetResource(MAPL, iUseInterp, 'INTERPOLATE_OCEAN_ICE_CURRENTS:', &
         default=0, RC=STATUS )
    VERIFY_(STATUS)
    useInterp = (iUseInterp /= 0)

    if(associated(UW)) then
       call MAPL_LocStreamTransform( ExchGrid, UW     ,  UWO   , &
         INTERP=useInterp, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(VW)) then
       call MAPL_LocStreamTransform( ExchGrid, VW     ,  VWO   , &
         INTERP=useInterp, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(UI)) then
       call MAPL_LocStreamTransform( ExchGrid, UI     ,  UIO   , &
         INTERP=useInterp, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(VI)) then
       call MAPL_LocStreamTransform( ExchGrid, VI     ,  VIO   , &
         INTERP=useInterp, RC=STATUS)
       VERIFY_(STATUS)
    end if

    if(associated(KPAR)) then
       call MAPL_LocStreamTransform( ExchGrid, KPAR   ,  KPARO , RC=STATUS) 
       VERIFY_(STATUS)
    end if

    if(associated(TS_FOUND)) then
       call MAPL_LocStreamTransform( ExchGrid, TS_FOUND   ,  TS_FOUNDO , &
         INTERP=useInterp, RC=STATUS) 
       VERIFY_(STATUS)
    end if

    if(associated(SS_FOUND)) then
       call MAPL_LocStreamTransform( ExchGrid, SS_FOUND   ,  SS_FOUNDO , &
         INTERP=useInterp, RC=STATUS) 
       VERIFY_(STATUS)
    end if

    if(associated(FRZMLT)) then
       call MAPL_LocStreamTransform( ExchGrid, FRZMLT ,  FRZMLTO , RC=STATUS) 
       VERIFY_(STATUS)
    end if

!  All done
!-----------

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"RUN" )

    RETURN_(ESMF_SUCCESS)
   
    contains

    subroutine OBIO_RunTransforms(DO_DATA_ATM4OCN, RC)

      logical,                    intent(IN   ) ::  DO_DATA_ATM4OCN
      integer, optional,          intent(  OUT) ::  RC
      
      character(len=ESMF_MAXSTR), parameter :: IAm="OBIO_RunTransforms"
      integer                               :: STATUS

      call MAPL_GetPointer(IMPORT, UU      ,  'UU',      RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, CO2SC   ,  'CO2SC'  , RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetPointer(IMPORT, DUDP    ,  'DUDP'   , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, DUWT    ,  'DUWT'   , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, DUSD    ,  'DUSD'   , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, DRBAND    , 'DRBAND'    , RC=STATUS)
      VERIFY_(STATUS)
      call MAPL_GetPointer(IMPORT, DFBAND    , 'DFBAND'    , RC=STATUS)
      VERIFY_(STATUS)

      call MAPL_GetPointer(GIM(OBIO ), USTR3B  ,  'OUSTAR3'  , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(OBIO ), UUB     ,  'UU'       , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(OBIO ), PSB     ,  'PS'       , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(ORAD ), UUR     ,  'UU'       , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(ORAD ), PSR     ,  'PS'       , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(GIM(OBIO ), CO2SCB  ,  'CO2SC'    , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetPointer(GIM(OBIO ), DUDPB   ,  'DUDP'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(OBIO ), DUWTB   ,  'DUWT'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(OBIO ), DUSDB   ,  'DUSD'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(ORAD ), DRBANDR  , 'DRBAND'   , notfoundOK=.true.,  RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(GIM(ORAD ), DFBANDR  , 'DFBAND'   , notfoundOK=.true.,  RC=STATUS); VERIFY_(STATUS)          
      call MAPL_GetPointer(GIM(OBIO ), DISCHARGEOB   ,  'DISCHARGE'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)

      if(.not. DO_DATA_ATM4OCN) then
        call MAPL_GetPointer(IMPORT, BCDP      , 'BCDP'      , RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT, BCWT      , 'BCWT'      , RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT, OCDP      , 'OCDP'      , RC=STATUS)
        VERIFY_(STATUS)
        call MAPL_GetPointer(IMPORT, OCWT      , 'OCWT'      , RC=STATUS)
        VERIFY_(STATUS)
      end if

      if(associated(UUB)) then
         call MAPL_LocStreamTransform( ExchGrid, UUB     ,  UU     , RC=STATUS) 
         VERIFY_(STATUS)
      endif
      if(associated(UUR)) then
         call MAPL_LocStreamTransform( ExchGrid, UUR     ,  UU     , RC=STATUS)
         VERIFY_(STATUS)
      endif
      if(associated(PSB)) then
         call MAPL_LocStreamTransform( ExchGrid, PSB     ,  PS     , RC=STATUS)    
         VERIFY_(STATUS)
      endif
      if(associated(PSR)) then
         call MAPL_LocStreamTransform( ExchGrid, PSR     ,  PS     , RC=STATUS)    
         VERIFY_(STATUS)
      endif

      if(associated(USTR3B)) then
         call MAPL_LocStreamTransform( ExchGrid, USTR3B  ,  USTR3  , RC=STATUS)
         VERIFY_(STATUS)
      endif
      if(associated(CO2SCB)) then
         call MAPL_LocStreamTransform( ExchGrid, CO2SCB  ,  CO2SC  , RC=STATUS) 
         VERIFY_(STATUS)
      endif

      if(associated(DUDPB)) then
       do N = 1, NUM_DUDP
          call MAPL_LocStreamTransform( ExchGrid, DUDPB(:,:,N), DUDP(:,N), RC=STATUS )
          VERIFY_(STATUS)
       end do
      endif
      if(associated(DUWTB)) then
       do N = 1, NUM_DUWT
          call MAPL_LocStreamTransform( ExchGrid, DUWTB(:,:,N), DUWT(:,N), RC=STATUS )
          VERIFY_(STATUS)
       end do
      endif
      if(associated(DUSDB)) then
       do N = 1, NUM_DUSD
          call MAPL_LocStreamTransform( ExchGrid, DUSDB(:,:,N), DUSD(:,N), RC=STATUS )
          VERIFY_(STATUS)
       end do
      endif
      if(associated(DRBANDR)) then
        do N = 1, NB_OBIO
           call MAPL_LocStreamTransform( ExchGrid, DRBANDR(:,:,N), DRBAND(:,N),   RC=STATUS )
           VERIFY_(STATUS)
        end do
      endif
      if(associated(DFBANDR)) then
        do N = 1, NB_OBIO
           call MAPL_LocStreamTransform( ExchGrid, DFBANDR(:,:,N), DFBAND(:,N),   RC=STATUS )
           VERIFY_(STATUS)
        end do
      endif

      if(.not. DO_DATA_ATM4OCN) then
        call MAPL_GetPointer(GIM(OBIO ), BCDPB   ,  'BCDP'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(GIM(OBIO ), BCWTB   ,  'BCWT'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(GIM(OBIO ), OCDPB   ,  'OCDP'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(GIM(OBIO ), OCWTB   ,  'OCWT'     , notfoundOK=.true., RC=STATUS); VERIFY_(STATUS)
      end if

      if(.not. DO_DATA_ATM4OCN) then
       if(associated(BCDPB)) then
          do N = 1, NUM_BCDP
             call MAPL_LocStreamTransform( ExchGrid, BCDPB(:,:,N), BCDP(:,N), RC=STATUS )
             VERIFY_(STATUS)
          end do
       endif
       if(associated(BCWTB)) then
          do N = 1, NUM_BCWT
             call MAPL_LocStreamTransform( ExchGrid, BCWTB(:,:,N), BCWT(:,N), RC=STATUS )
             VERIFY_(STATUS)
          end do
       endif
       if(associated(OCDPB)) then
          do N = 1, NUM_OCDP
             call MAPL_LocStreamTransform( ExchGrid, OCDPB(:,:,N), OCDP(:,N), RC=STATUS )
             VERIFY_(STATUS)
          end do
       endif
       if(associated(OCWTB)) then
          do N = 1, NUM_OCWT
             call MAPL_LocStreamTransform( ExchGrid, OCWTB(:,:,N), OCWT(:,N), RC=STATUS )
             VERIFY_(STATUS)
          end do
       endif
      end if

      RETURN_(ESMF_SUCCESS)
    end subroutine OBIO_RunTransforms

  end subroutine RUN

end module GEOS_OgcmGridCompMod

