!  $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_MITDynaGridCompMod

!BOP

! !MODULE: GEOS_MITDyna -- The MIT seaice Dynamics model

! !DESCRIPTION:
! 
!   {\tt GEOS\_DataSeaIce is a gridded component that reads the 
!   ocean\_bcs file 
!   This module interpolates the SST and sea ice data from 
!   either daily or monthly values to the correct time of the simulation.
!   Data are read only if the simulation time is not in the save interval.
!   Surface Albedo and Surface roughness calculations are also takencare of in
!   this module.
!

! !USES: 

  use ESMF
  use MAPL

  use ice_domain_size,    only: init_domain_size
  use ice_init,           only: alloc_dyna_arrays, dealloc_dyna_arrays
  use ice_work,           only: init_work

  implicit none
  private


! !PUBLIC MEMBER FUNCTIONS:

  public SetServices


  integer, parameter :: NUM_3D_ICE_TRACERS=3
  integer            :: NUM_ICE_CATEGORIES
  integer            :: NUM_ICE_LAYERS
  integer, parameter :: NUM_SNOW_LAYERS=1
  integer            :: NUM_ICE_LAYERS_ALL
  integer            :: NUM_SNOW_LAYERS_ALL


!=============================================================================


   contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! ! INTERFACE:

  subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

!  !DESCRIPTION: This version uses the MAPL\_GenericSetServices. This function sets
!                the Initialize and Finalize services, as well as allocating
!   our instance of a generic state and putting it in the 
!   gridded component (GC). Here we only need to set the run method and
!   add the state variable specifications (also generic) to our instance
!   of the generic state. This is the way our true state variables get into
!   the ESMF\_State INTERNAL, which is in the MAPL\_MetaComp.
!
!EOP

!=============================================================================
!
! ErrLog Variables


    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Local derived type aliases

    type (MAPL_MetaComp),  pointer          :: MAPL
    type (ESMF_Config)                      :: CF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = "SetServices"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get constants from CF
! ---------------------

    call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS,     Label="CICE_N_ICE_LAYERS:" ,     RC=STATUS)
    VERIFY_(STATUS)

    NUM_ICE_LAYERS_ALL  = NUM_ICE_LAYERS  * NUM_ICE_CATEGORIES
    NUM_SNOW_LAYERS_ALL = NUM_SNOW_LAYERS * NUM_ICE_CATEGORIES

! Set the Run entry point
! -----------------------
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,        Run,        RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,   Finalize,   RC=status)
    VERIFY_(STATUS)


! Set the state variable specs.
! -----------------------------

!BOS

! !Import state:

  call MAPL_AddImportSpec(GC,                            &
    SHORT_NAME         = 'FRACICE',                           &
    LONG_NAME          = 'fractional_cover_of_seaice',        &
    UNITS              = '1',                                 &
   ! PRECISION          = ESMF_KIND_R8,                        &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    RESTART            = MAPL_RestartOptional,                &
    DEFAULT            = 0.0,                                 &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    SHORT_NAME         = 'TI',                                &
    LONG_NAME          = 'seaice_skin_temperature',           &
    UNITS              = 'K',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    RESTART            = MAPL_RestartOptional,                &
    DEFAULT            = MAPL_TICE,                           &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    SHORT_NAME         = 'SI',                                &
    LONG_NAME          = 'seaice_skin_salinity',              &
    UNITS              = 'psu',                               &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    RESTART            = MAPL_RestartOptional,                &
    DEFAULT            = 4.,                                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                                &
    SHORT_NAME         = 'VOLICE',                            &
    LONG_NAME          = 'ice_category_volume_per_unit_area_of_grid_cell',&
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    RESTART            = MAPL_RestartOptional,                &
    DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                                &
    SHORT_NAME         = 'VOLSNO',                            &
    LONG_NAME          = 'sno_category_volume_per_unit_area_of_grid_cell',&
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    RESTART            = MAPL_RestartOptional,                &
    DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                                &
        SHORT_NAME         = 'ERGICE',                            &
        LONG_NAME          = 'ice_category_layer_internal_energy',&
        UNITS              = 'J m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS_ALL/),              &
        !VLOCATION          = MAPL_VLocationCenter,                 &
    !    DEFAULT            = 0.0,                                 &
        RESTART            = MAPL_RestartSkip,                    &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                &
        SHORT_NAME         = 'ERGSNO',                            &
        LONG_NAME          = 'snow_category_layer_internal_energy',&
        UNITS              = 'J m-2',                             &
        !DIMS               = MAPL_DimsHorzVert,                   &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS_ALL/),             &
        !VLOCATION          = MAPL_VLocationCenter,                 &
        RESTART            = MAPL_RestartSkip,                    &
        !DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                                     ,&
        LONG_NAME          = 'melt_pond_volume'                     ,&
        UNITS              = 'm'                                ,&
        SHORT_NAME         = 'MPOND'                                 ,&
        !DIMS               = MAPL_DimsHorzVert,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        !VLOCATION          = MAPL_VLocationCenter,                 &
        !DEFAULT            = 0.0                                    ,&
        RESTART            = MAPL_RestartSkip,                    &
        RC=STATUS                                                 )

     VERIFY_(STATUS)

   call MAPL_AddImportSpec(GC,                                &
        SHORT_NAME         = 'TAUAGE',                            &
        LONG_NAME          = 'volume_weighted_mean_ice_age',      &
        UNITS              = 's',                                 &
        !DIMS               = MAPL_DimsHorzVert,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        !VLOCATION          = MAPL_VLocationCenter,                 &
        RESTART            = MAPL_RestartSkip,                    &
        !DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'eastward_stress_on_ice'            ,&
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUX'                              ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
    RESTART            = MAPL_RestartSkip,                    &
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                            &
    LONG_NAME          = 'northward_stress_on_ice',           &
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUY'                              ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
    RESTART            = MAPL_RestartSkip,                    &
                                           RC=STATUS          )
  VERIFY_(STATUS)


    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'UW',                                &
         LONG_NAME          = 'water_skin_eastward_velocity',      &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         DEFAULT            = 0.0,                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'VW',                                &
         LONG_NAME          = 'water_skin_northward_velocity',&
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         DEFAULT            = 0.0,                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'UWB',                                &
         LONG_NAME          = 'water_skin_eastward_velocity', &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         DEFAULT            = 0.0,                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'VWB',                                &
         LONG_NAME          = 'water_skin_northward_velocity',&
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         DEFAULT            = 0.0,                                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'SLV',                           &
         LONG_NAME          = 'sea_level_with_ice_loading',     &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         DEFAULT            = 0.0,                                 &
         RC=status  )
    VERIFY_(status)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'FROCEAN',                           &
         LONG_NAME          = 'fraction_of_gridbox_covered_by_skin',&
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RESTART            = MAPL_RestartSkip,                    &
         RC=status  )
    VERIFY_(status)

    call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'HI',                                &
         LONG_NAME          = 'seaice_skin_layer_depth',            &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'FRESH',                         &
         LONG_NAME          = 'fresh_water_flux_due_to_thermodynamics', &
         UNITS              = 'kg m-2 s-1',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RESTART            = MAPL_RestartSkip,                    &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'FSALT',                         &
         LONG_NAME          = 'salt_flux_due_to_thermodynamics', &
         UNITS              = 'kg m-2 s-1',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RESTART            = MAPL_RestartSkip,                    &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                            &
         SHORT_NAME         = 'FHOCN',                         &
         LONG_NAME          = 'heat_flux_due_to_thermodynamics', &
         UNITS              = 'W m-2',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RESTART            = MAPL_RestartSkip,                    &
                                                   RC=STATUS  )
    VERIFY_(STATUS)

! !Internal state:

   call MAPL_AddInternalSpec(GC                                     ,&
        LONG_NAME          = 'x_component_of_velocity'              ,&
        UNITS              = 'm s-1'                                ,&
        SHORT_NAME         = 'UVEL'                                 ,&
        PRECISION          = ESMF_KIND_R8,                        &
        DIMS               = MAPL_DimsHorzOnly                      ,&
        VLOCATION          = MAPL_VLocationNone                     ,&
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                     )

     VERIFY_(STATUS)

   call MAPL_AddInternalSpec(GC                                        ,&
        LONG_NAME          = 'y_component_of_velocity'              ,&
        UNITS              = 'm s-1'                                   ,&
        SHORT_NAME         = 'VVEL'                                     ,&
        PRECISION          = ESMF_KIND_R8,                        &
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)


   call MAPL_AddInternalSpec(GC                                        ,&
        LONG_NAME          = 'strain_rate_I_component_velocity_divergence' ,&
        UNITS              = 's-1'                                            ,&
        SHORT_NAME         = 'DIVU'                                     ,&
        PRECISION          = ESMF_KIND_R8,                        &
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

   call MAPL_AddInternalSpec(GC                                        ,&
        LONG_NAME          = 'strain_rate_II_component'                ,&
        UNITS              = 's-1'                                     ,&
        PRECISION          = ESMF_KIND_R8,                        &
        SHORT_NAME         = 'SHEAR'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                            )

   VERIFY_(STATUS)

   call MAPL_AddInternalSpec(GC                                        ,&
        LONG_NAME          = 'ice_strength'                ,&
        UNITS              = 'N m-1'                                     ,&
        SHORT_NAME         = 'STRENGTH'                                ,&
        PRECISION          = ESMF_KIND_R8,                        &
        DIMS               = MAPL_DimsHorzOnly                         ,&
        VLOCATION          = MAPL_VLocationNone                        ,&
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                            )

    VERIFY_(STATUS)

   call MAPL_AddInternalSpec(GC                                        ,&
        LONG_NAME          = 'ice_extent_mask_at_U_cell'               ,&
        UNITS              = '1'                                       ,&
        SHORT_NAME         = 'ICEUMASK'                                ,&
        DIMS               = MAPL_DimsHorzOnly                         ,&
        VLOCATION          = MAPL_VLocationNone                        ,&
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                            )

    VERIFY_(STATUS)


   ! this internal state has been reconfigured to have a dimension 
   ! of 2D + UNGRIDDED or (IM, JM, 50).
   ! It having UNGRIDDED_DIMS=50 makes the current restarts compatible.
   ! the first 12 slices are for restart of ice_stress_tensor_component,
   ! the rest will have some other uses or be left unused 
   call MAPL_AddInternalSpec(GC                                        ,&
        LONG_NAME          = 'ice_stress_tensor_component'             ,&
        UNITS              = 's-1'                                     ,&
        PRECISION          = ESMF_KIND_R8,                        &
        SHORT_NAME         = 'STRESSCOMP'                                ,&
        DIMS               = MAPL_DimsHorzVert                              ,&
        VLOCATION          = MAPL_VLocationCenter                       ,&
        !UNGRIDDED_DIMS     = (/50/),                              &
        !DIMS               = MAPL_DimsHorzOnly,                   &
        !VLOCATION          = MAPL_VLocationNone,                  &
        DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                            )

   VERIFY_(STATUS)


!  !Export state:

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'UI',                                &
    LONG_NAME          = 'zonal_velocity_of_surface_seaice',   &
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VI',                                &
    LONG_NAME          = 'meridional_velocity_of_surface_seaice',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VEL',                                &
    LONG_NAME          = 'ice_drift_speed',&
    UNITS              = 'm s-1 ',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'TAUXBOT',                           &
        LONG_NAME          = 'eastward_stress_at_base_of_ice',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'TAUYBOT',                           &
        LONG_NAME          = 'northward_stress_at_base_of_ice',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'TAUXOCNB',                           &
        LONG_NAME          = 'x_stress_at_base_of_ice',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'TAUYOCNB',                           &
        LONG_NAME          = 'y_stress_at_base_of_ice',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'STROCNXB',                           &
        LONG_NAME          = 'x_stress_at_base_of_ice_weighted_by_aiu',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                                  &
        SHORT_NAME         = 'STROCNYB',                           &
        LONG_NAME          = 'y_stress_at_base_of_ice_weighted_by_aiu',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC,                            &
        SHORT_NAME         = 'FRACICE',                           &
        LONG_NAME          = 'fractional_cover_of_seaice',        &
        UNITS              = '1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'ice_strength'                            ,&
        UNITS              = 'N m-1'                                   ,&
        SHORT_NAME         = 'STRENGTH'                                ,&
        DIMS               = MAPL_DimsHorzOnly                         ,&
        VLOCATION          = MAPL_VLocationNone                        ,&
        RC=STATUS                                                       )

     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'strain_rate_I_component_velocity_divergence' ,&
        UNITS              = 's-1'                                            ,&
        SHORT_NAME         = 'DIVU'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'strain_rate_I_component_ocean_velocity_divergence' ,&
        UNITS              = 's-1'                                            ,&
        SHORT_NAME         = 'DIVUOCN'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'faked_velocity_divergence' ,&
        UNITS              = 's-1'                                            ,&
        SHORT_NAME         = 'FAKEDIV'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'air_ice_stress_divergence' ,&
        UNITS              = 'N m-3'                                            ,&
        SHORT_NAME         = 'TAUADIV'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'air_ice_stress_no_aice_scaling_divergence' ,&
        UNITS              = 'N m-3'                                            ,&
        SHORT_NAME         = 'TAU1DIV'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

     call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'ocean_ice_stress_divergence' ,&
        UNITS              = 'N m-3'                                            ,&
        SHORT_NAME         = 'TAUODIV'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

     VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC                                        ,&
        LONG_NAME          = 'strain_rate_II_component'                ,&
        UNITS              = 's-1'                                     ,&
        SHORT_NAME         = 'SHEAR'                                     ,&
        DIMS               = MAPL_DimsHorzOnly                              ,&
        VLOCATION          = MAPL_VLocationNone                             ,&
        RC=STATUS                                                            )

   VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'HICE',                            &
    LONG_NAME          = 'mean_ice_thickness_of_grid_cell',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'HICE0',                            &
    LONG_NAME          = 'mean_ice_thickness_of_grid_cell_covered_by_ice',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'HSNO',                            &
    LONG_NAME          = 'mean_snow_thickness_of_grid_cell',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'HSNO0',                            &
    LONG_NAME          = 'mean_snow_thickness_of_grid_cell_covered_by_ice',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'VICE0',                            &
    LONG_NAME          = 'ice_volume_of_grid_cell_before_cleaup',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'VICE1',                            &
    LONG_NAME          = 'ice_volume_of_grid_cell_after_cleaup',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'DRAFT',                            &
    LONG_NAME          = 'mean_ice_draft_of_grid_cell',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'DRAFT0',                            &
    LONG_NAME          = 'mean_ice_draft_of_grid_cell_covered_by_ice',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'AICE',                            &
    LONG_NAME          = 'ice_concentration_of_grid_cell',   &
    UNITS              = '1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'AICEU',                            &
    LONG_NAME          = 'ice_concentration_of_grid_cell_Bgrid',   &
    UNITS              = '1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)


  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'DAIDTD',                                &
        LONG_NAME          = 'ice_area_tendency_dueto_dynamics', &
        UNITS              = '% day-1',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'DVIDTD',                                &
        LONG_NAME          = 'ice_volume_tendency_dueto_dynamics', &
        UNITS              = 'cm day-1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'DVIRDGDT',                                &
        LONG_NAME          = 'rate_of_ice_volume_ridged',               &
        UNITS              = 'cm day-1',                                 &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'STRCORX',                                &
        LONG_NAME          = 'stress_due_to_coriolis_effect_x_direction', &
        UNITS              = 'N m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'STRCORY',                                &
        LONG_NAME          = 'stress_due_to_coriolis_effect_y_direction', &
        UNITS              = 'N m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'STRTLTX',                                &
        LONG_NAME          = 'stress_due_to_sea_surface_slope_x_direction', &
        UNITS              = 'N m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'STRTLTY',                                &
        LONG_NAME          = 'stress_due_to_sea_surface_slope_y_direction', &
        UNITS              = 'N m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'STRINTX',                                &
        LONG_NAME          = 'divergence_of_internal_ice_stress_x_direction', &
        UNITS              = 'N m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                           &
        SHORT_NAME         = 'STRINTY',                                &
        LONG_NAME          = 'divergence_of_internal_ice_stress_y_direction', &
        UNITS              = 'N m-2',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'eastward_stress_on_ice'            ,&
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUXI'                             ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'northward_stress_on_ice',           &
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUYI'                             ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'smoothed_eastward_stress_on_ice'   ,&
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUXISM'                           ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'smoothed_northward_stress_on_ice',  &
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUYISM'                           ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'x_stress_on_ice_Bgrid'            ,&
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUXIB'                             ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    LONG_NAME          = 'y_stress_on_ice_Bgrid',           &
    UNITS              = 'N m-2'                             ,&
    SHORT_NAME         = 'TAUYIB'                             ,&
    DIMS               = MAPL_DimsHorzOnly                   ,&
    VLOCATION          = MAPL_VLocationNone                  ,&
                                           RC=STATUS          )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                               &
         LONG_NAME          = 'water_skin_eastward_velocity', &
         UNITS              = 'm s-1 ',                            &
         SHORT_NAME         = 'UOCN',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                               &
         LONG_NAME          = 'water_skin_northward_velocity',&
         UNITS              = 'm s-1 ',                            &
         SHORT_NAME         = 'VOCN',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                               &
         LONG_NAME          = 'sea_surface_height',                &
         UNITS              = 'm',                                 &
         SHORT_NAME         = 'SSH',                           &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                               &
         LONG_NAME          = 'sea_level_with_ice_loading',                &
         UNITS              = 'm',                                 &
         SHORT_NAME         = 'SLV',                           &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    LONG_NAME          = 'open_water_fraction_of_grid_cell',   &
    UNITS              = '1',                                 &
    SHORT_NAME         = 'FROCEAN',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'AREA',                            &
    LONG_NAME          = 'area_of_grid_cell',   &
    UNITS              = 'm2',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'TMASK',                            &
    LONG_NAME          = 'ocean_mask_for_sea_ice',   &
    UNITS              = '1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                                        ,&
    SHORT_NAME         = 'ICEUMASK0'                                ,&
    LONG_NAME          = 'ice_extent_mask_at_U_cell_at_the_beginning' ,&
    UNITS              = '1'                                       ,&
    DIMS               = MAPL_DimsHorzOnly                         ,&
    VLOCATION          = MAPL_VLocationNone                        ,&
    RC=STATUS                                                            )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                                        ,&
    SHORT_NAME         = 'ICEUMASK1'                                ,&
    LONG_NAME          = 'ice_extent_mask_at_U_cell_at_the_end' ,&
    UNITS              = '1'                                       ,&
    DIMS               = MAPL_DimsHorzOnly                         ,&
    VLOCATION          = MAPL_VLocationNone                        ,&
    RC=STATUS                                                            )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                                        ,&
    SHORT_NAME         = 'UVEL0'                                ,&
    LONG_NAME          = 'x_velocity_at_U_cell_at_the_beginning' ,&
    UNITS              = 'm s-1'                                       ,&
    DIMS               = MAPL_DimsHorzOnly                         ,&
    VLOCATION          = MAPL_VLocationNone                        ,&
    RC=STATUS                                                            )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                                        ,&
    SHORT_NAME         = 'UVEL1'                                ,&
    LONG_NAME          = 'x_velocity_at_U_cell_at_the_end' ,&
    UNITS              = 'm s-1'                                       ,&
    DIMS               = MAPL_DimsHorzOnly                         ,&
    VLOCATION          = MAPL_VLocationNone                        ,&
    RC=STATUS                                                            )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                                        ,&
    SHORT_NAME         = 'VVEL0'                                ,&
    LONG_NAME          = 'y_velocity_at_U_cell_at_the_beginning' ,&
    UNITS              = 'm s-1'                                       ,&
    DIMS               = MAPL_DimsHorzOnly                         ,&
    VLOCATION          = MAPL_VLocationNone                        ,&
    RC=STATUS                                                            )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC                                        ,&
    SHORT_NAME         = 'VVEL1'                                ,&
    LONG_NAME          = 'y_velocity_at_U_cell_at_the_end' ,&
    UNITS              = 'm s-1'                                       ,&
    DIMS               = MAPL_DimsHorzOnly                         ,&
    VLOCATION          = MAPL_VLocationNone                        ,&
    RC=STATUS                                                            )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'ANGLET',                            &
    LONG_NAME          = 'rotation_angle_on_T_cell_center',   &
    UNITS              = 'radians',                           &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'ANGLEU',                            &
    LONG_NAME          = 'rotation_angle_on_U_cell_center',   &
    UNITS              = 'radians',                           &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TRANSIX',                         &
    LONG_NAME          = 'x_component_of_sea_ice_volume_transport', &
    UNITS              = 'm3 s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TRANSIY',                         &
    LONG_NAME          = 'y_component_of_sea_ice_volume_transport', &
    UNITS              = 'm3 s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC                    ,&
      SHORT_NAME         = 'HICEUNT',                     &
      LONG_NAME          = 'grid_cell_mean_ice_thickness_untouched_by_run',&
      UNITS              = 'm'                         ,&
      DIMS               = MAPL_DimsHorzOnly           ,&
      VLOCATION          = MAPL_VLocationNone          ,&
                                           RC=STATUS  )
   VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TRANSIMX',                         &
    LONG_NAME          = 'x_component_of_sea_ice_mass_transport_tri', &
    UNITS              = 'kg s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TRANSIMY',                         &
    LONG_NAME          = 'y_component_of_sea_ice_mass_transport_tri', &
    UNITS              = 'kg s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VOLICE0',                         &
    LONG_NAME          = 'total_ice_volume_before_dynamics', &
    UNITS              = 'm',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VOLSNO0',                         &
    LONG_NAME          = 'total_snow_volume_before_dynamics', &
    UNITS              = 'm',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VOLICE1',                         &
    LONG_NAME          = 'total_ice_volume_after_dynamics', &
    UNITS              = 'm',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VOLSNO1',                         &
    LONG_NAME          = 'total_snow_volume_after_dynamics', &
    UNITS              = 'm',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'FRESH',                         &
    LONG_NAME          = 'fresh_water_flux_into_ocean', &
    UNITS              = 'kg m-2 s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'FSALT',                         &
    LONG_NAME          = 'salt_flux_into_ocean', &
    UNITS              = 'kg m-2 s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'FHOCN',                         &
    LONG_NAME          = 'heat_flux_into_ocean', &
    UNITS              = 'W m-2',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'SIG1',                         &
    LONG_NAME          = 'normalized_principal_stress_1', &
    UNITS              = '1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'SIG2',                         &
    LONG_NAME          = 'normalized_principal_stress_2', &
    UNITS              = '1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'DAIDTNUDG',                         &
    LONG_NAME          = 'ice_area_tendency_due_to_nudging',  &
    UNITS              = '% day-1',                           &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'DVIDTNUDG',                         &
    LONG_NAME          = 'ice_volume_tendency_due_to_nudging',&
    UNITS              = 'cm day-1',                          &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

! category dimensional exports
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'AICEN',                           &
    LONG_NAME          = 'seaice_area_for_each_category',     &
    UNITS              = '1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VICEN',                           &
    LONG_NAME          = 'seaice_volume_for_each_category',   &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'HIFLXE',                           &
    LONG_NAME          = 'ice_volume_transports_across_E_cell_edges',     &
    UNITS              = 'm3 s-1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'HIFLXN',                           &
    LONG_NAME          = 'ice_volume_transports_across_N_cell_edges',     &
    UNITS              = 'm3 s-1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'STRESSES0',                           &
    LONG_NAME          = 'ice_stress_components_at_the_beginning',     &
    UNITS              = 'N m-1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/12/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'STRESSES1',                           &
    LONG_NAME          = 'ice_stress_components_at_the_end',     &
    UNITS              = 'N m-1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/12/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  !
  ! the following exports are specifically for CMIP5 
  !
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'AICE_CMIP5',                        &
    LONG_NAME          = 'sea_ice_area_fraction',             &
    UNITS              = '%',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'HICE_CMIP5',                        &
    LONG_NAME          = 'sea_ice_thickness',             &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'HSNO_CMIP5',                        &
    LONG_NAME          = 'surface_snow_thickness',             &
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'FWM_CMIP5',                         &
    LONG_NAME          = 'frozen_water_mass',                 &
    UNITS              = 'kg m-2',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'WIST_CMIP5',                         &
    LONG_NAME          = 'weighted_ice_surface_temperature',   &
    UNITS              = 'K',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'WTSNINT_CMIP5',                         &
    LONG_NAME          = 'weighted_temperature_at_interface_between_sea_ice_and_snow',   &
    UNITS              = 'K',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TRANSIX_CMIP5',                         &
    LONG_NAME          = 'x_component_of_sea_ice_mass_transport', &
    UNITS              = 'kg s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TRANSIY_CMIP5',                         &
    LONG_NAME          = 'y_component_of_sea_ice_mass_transport', &
    UNITS              = 'kg s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'UI_CMIP5',                         &
    LONG_NAME          = 'x_component_of_sea_ice_velocity', &
    UNITS              = 'm s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)
  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'VI_CMIP5',                         &
    LONG_NAME          = 'y_component_of_sea_ice_velocity', &
    UNITS              = 'm s-1',                            &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'ICESTATES',                      &
    LONG_NAME          = 'container_for_seaice_variables_for_MITgcm', &
    UNITS              = 'N/A',                            &
    DIMS               = MAPL_DimsHorzOnly,                &
    VLOCATION          = MAPL_VLocationNone,               &
    DATATYPE           = MAPL_StateItem,                   &
!    RESTART            = MAPL_RestartSkip,                 &
                                                   RC=STATUS  )
  VERIFY_(STATUS)


! !Export state - increments for seaice:

  call MAPL_AddExportSpec(GC,                                    &
       SHORT_NAME         = 'DEL_FRACICE',                        &
       LONG_NAME          = 'delta_fractional_cover_of_seaice',  &
       UNITS              = '1',                                 &
       DIMS               = MAPL_DimsHorzOnly,                   &
       UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
       VLOCATION          = MAPL_VLocationNone,                  &
       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'DEL_TI',                                &
    LONG_NAME          = 'delta_seaice_skin_temperature',           &
    UNITS              = 'K',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    DEFAULT            = MAPL_TICE,                           &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'DEL_SI',                                &
    LONG_NAME          = 'delta_seaice_skin_salinity',              &
    UNITS              = 'psu',                               &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
    DEFAULT            = 4.,                                  &
                                                   RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'DEL_VOLICE',                            &
    LONG_NAME          = 'delta_ice_category_volume_per_unit_area_of_grid_cell',&
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'DEL_VOLSNO',                            &
    LONG_NAME          = 'delta_sno_category_volume_per_unit_area_of_grid_cell',&
    UNITS              = 'm',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
   DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'DEL_ERGICE',                            &
        LONG_NAME          = 'delta_ice_category_layer_internal_energy',&
        UNITS              = 'J m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        UNGRIDDED_DIMS     = (/NUM_ICE_LAYERS_ALL/),              &
        !VLOCATION          = MAPL_VLocationCenter,                 &
    ! DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'DEL_ERGSNO',                            &
        LONG_NAME          = 'delta_snow_category_layer_internal_energy',&
        UNITS              = 'J m-2',                             &
        !DIMS               = MAPL_DimsHorzVert,                   &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        UNGRIDDED_DIMS     = (/NUM_SNOW_LAYERS_ALL/),             &
        !VLOCATION          = MAPL_VLocationCenter,                 &
       !DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC                                     ,&
        LONG_NAME          = 'delta_melt_pond_volume'                     ,&
        UNITS              = 'm'                                ,&
        SHORT_NAME         = 'DEL_MPOND'                                 ,&
        !DIMS               = MAPL_DimsHorzVert,                   &
       UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        !VLOCATION          = MAPL_VLocationCenter,                 &
        !DEFAULT            = 0.0                                    ,&
        RC=STATUS                                                 )

     VERIFY_(STATUS)

   call MAPL_AddExportSpec(GC,                                &
        SHORT_NAME         = 'DEL_TAUAGE',                            &
        LONG_NAME          = 'delta_volume_weighted_mean_ice_age',      &
        UNITS              = 's',                                 &
        !DIMS               = MAPL_DimsHorzVert,                   &
        UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        !VLOCATION          = MAPL_VLocationCenter,                 &
        !DEFAULT            = 0.0,                                 &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                            &
         SHORT_NAME         = 'DEL_HI',                                &
         LONG_NAME          = 'delta_seaice_skin_layer_depth',            &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
                                                   RC=STATUS  )
    VERIFY_(STATUS)


!EOS

    call MAPL_TimerAdd(GC,    name="INITIALIZE",RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="RUN"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,    name="FINALIZE"  ,RC=STATUS)
    VERIFY_(STATUS)

! Set generic init and final methods
! ----------------------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices

! !IROUTINE: Initialize -- Initialize method for the GEOS CICE Dynamics

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: The Initialize method of the CICE dynamics Gridded Component.
!   It does some initializing work on CICE data structures 
!   It then does a Generic_Initialize

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

    integer                                 :: Comm
    integer                                 :: NDTE   ! number of subcycles
    integer                                 :: EVP_DAMPING   
    integer                                 :: NDYN_DT  
    integer                                 :: STRENGTH, RDG_PARTIC, RDG_REDIST
    real                                    :: MU_RDG, PSTAR, IODRAG
    real                                    :: DTI 
    real(kind=ESMF_KIND_R8)                 :: DTR
    character(len=ESMF_MAXSTR)              :: GRIDTYPE
    character(len=ESMF_MAXSTR)              :: GRIDFILE
    character(len=ESMF_MAXSTR)              :: GRIDFILEFORM
    character(len=ESMF_MAXSTR)              :: KMTFILE
    character(len=ESMF_MAXSTR)              :: KMTFILEFORM
    character(len=ESMF_MAXSTR)              :: PROCESSOR_SHAPE
    character(len=ESMF_MAXSTR)              :: DISTRIBUTION_TYPE
    character(len=ESMF_MAXSTR)              :: DISTRIBUTION_WGHT 
    character(len=ESMF_MAXSTR)              :: EW_BOUNDARY_TYPE
    character(len=ESMF_MAXSTR)              :: NS_BOUNDARY_TYPE
    character(len=ESMF_MAXSTR)              :: ADVECTION

! Local derived type aliases

    type (ESMF_VM)                         :: VM
    type (MAPL_MetaComp    ), pointer      :: MAPL
    type (ESMF_State)                      :: INTERNAL

    type(ESMF_Alarm)                       :: ALARM_OCNEXCH   
    type(ESMF_Alarm)                       :: ALARM   
    type(ESMF_TimeInterval)                :: RING_INTERVAL  
    type(ESMF_Time)                        :: RING_TIME  
    type(ESMF_Time)                        :: currTime  
    integer                                :: CPLS 


#ifdef MODIFY_TOPOGRAPHY
    type (MAPL_LocStream       )            :: EXCH
    real, allocatable                       :: FROCEAN(:,:)
#endif

! pointers to internal
   real                   , pointer, dimension(:,:  ) :: ICEUMASK
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:,:) :: STRESSCOMP
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:  ) :: UVEL
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:  ) :: VVEL

   logical,  allocatable          :: ICEUM(:,:) 

    integer                       :: IM, JM
    integer                       :: NXG, NYG
    integer                       :: NPES
    integer                       :: OGCM_IM, OGCM_JM
    integer                       :: OGCM_NX, OGCM_NY
    type (ESMF_State) :: state
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

! Generic initialize
! ------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL"     )

    call MAPL_Get(MAPL,             &
                  INTERNAL_ESMF_STATE = INTERNAL,   &
                  RUNALARM = ALARM,                 &
                  IM=IM, &
                  JM=JM, & 
                  NX=NXG, & 
                  NY=NYG, & 
                                RC=STATUS )
    VERIFY_(STATUS)


! CICE grid initialization using the communicator from the VM
!------------------------------------------------------
    call ESMF_VMGetCurrent(VM, rc=STATUS)
    VERIFY_(STATUS)
    
    call ESMF_VMGet(VM, mpiCommunicator=Comm, petCount=NPES, rc=STATUS)
    VERIFY_(STATUS)


    ! Get the ocean layout from the VM
    !---------------------------------

    call MAPL_GetResource( MAPL, OGCM_IM, Label="OGCM.IM_WORLD:", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, OGCM_JM, Label="OGCM.JM_WORLD:", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, OGCM_NX, Label="OGCM.NX:", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, OGCM_NY, Label="OGCM.NY:", RC=STATUS)
    VERIFY_(STATUS)

    ! CICE grid initialization
    !-------------------------

    ASSERT_(mod(OGCM_IM,OGCM_NX)==0)
    ASSERT_(mod(OGCM_JM,OGCM_NY)==0)
 
    if(MAPL_AM_I_ROOT()) then
       print*, 'Initializing CICE domain size with: '
       print*, 'OGCM_NX = ', OGCM_NX, 'OGCM_NY = ', OGCM_NY
       print*, 'OGCM_IM = ', OGCM_IM, 'OGCM_JM = ', OGCM_JM
    endif

    ! Do some necessary CICE initialization
    !--------------------------------------

    call init_domain_size(OGCM_IM, OGCM_JM, OGCM_NX, OGCM_NY, NPES)

    call alloc_dyna_arrays( MAPL_AM_I_Root(), Iam )
    call init_work            ! work arrays

    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE work array initialized in '//trim(comp_name)
    endif
 
    call ESMF_StateGet(EXPORT, 'ICESTATES', state, __RC__)
    call ESMF_StateAdd(state, [import,export], __RC__)
    ! ALT: the statement above looks funny (at best): 
    ! the state item ICESTATES 
    ! itself is contained in the export state of this component, and yet
    ! it contains the export (and import). It seems cyclical, 
    ! but no harm is done, and it is a convenient way to encapsulate
    ! the seaice variables and "ship" them to MITgcm
    call ESMF_AttributeSet(state, name='ICECOMPNAME', value=trim(comp_name),__RC__)

! All Done
!---------

    call MAPL_TimerOff(MAPL,"TOTAL"     )
    call MAPL_TimerOff(MAPL,"INITIALIZE")

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!BOP

! ! IROUTINE: RUN -- Run stage for the DataSeaIce component

! !INTERFACE:

subroutine RUN ( GC, IMPORT, EXPORT, CLOCK, RC )


! !ARGUMENTS:

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
  integer, optional,   intent(  out) :: RC     ! Error code:

! ! DESCRIPTION: Periodically refreshes the SST and Ice information.

!EOP

! ErrLog Variables

  character(len=ESMF_MAXSTR)          :: IAm
  integer                             :: STATUS
  character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals

  type (MAPL_MetaComp), pointer       :: MAPL
         
! pointers to export

   real, pointer, dimension(:,:)  :: FRESH
   real, pointer, dimension(:,:)  :: FSALT
   real, pointer, dimension(:,:)  :: FHOCN

! pointers to import

   real, pointer, dimension(:,:)  :: FROCEAN
   real, pointer, dimension(:,:)  :: FRESHTHM
   real, pointer, dimension(:,:)  :: FSALTTHM
   real, pointer, dimension(:,:)  :: FHOCNTHM


!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, _RC )
    Iam = trim(COMP_NAME) // Iam


! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, MAPL, _RC)

! Start Total timer
!------------------

    call MAPL_TimerOn(MAPL,"TOTAL")
    call MAPL_TimerOn(MAPL,"RUN" )
   
    call MAPL_GetPointer(IMPORT, FROCEAN, 'FROCEAN', _RC)
    call MAPL_GetPointer(IMPORT, FRESHTHM ,  'FRESH'     , _RC)
    call MAPL_GetPointer(IMPORT, FSALTTHM ,  'FSALT'     , _RC)
    call MAPL_GetPointer(IMPORT, FHOCNTHM ,  'FHOCN'     , _RC)
   
    call MAPL_GetPointer(EXPORT, FSALT   ,  'FSALT'  , _RC)
    call MAPL_GetPointer(EXPORT, FHOCN   ,  'FHOCN'  , _RC)

     if(associated(FRESH)) then
         FRESH = 0.0
!@         call get_ocn_fluxes(fre=FRESH) 
         where(FROCEAN /= MAPL_UNDEF)
            FRESH = FRESH + FRESHTHM 
         end where
     endif  
     if(associated(FSALT)) then
         FSALT = 0.0
!@         call get_ocn_fluxes(sal=FSALT) 
         where(FROCEAN /= MAPL_UNDEF)
            FSALT = FSALT + FSALTTHM
         end where
     endif  
     if(associated(FHOCN)) then
         FHOCN = 0.0
 !@        call get_ocn_fluxes(hea=FHOCN) 
         where(FROCEAN /= MAPL_UNDEF)
            FHOCN = FHOCN + FHOCNTHM
         end where
     endif  

!  All done
!-----------

   call MAPL_TimerOff(MAPL,"RUN"  )
   call MAPL_TimerOff(MAPL,"TOTAL")

   RETURN_(ESMF_SUCCESS)
end subroutine RUN

! !IROUTINE: Finalize        -- Finalize method for CICEDyna wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp),    pointer                   :: MAPL 


! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Finalize"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"   )
    call MAPL_TimerOn(MAPL,"FINALIZE")

    call dealloc_dyna_arrays( MAPL_AM_I_ROOT(), Iam )

    call MAPL_TimerOff(MAPL,"FINALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"   )

! Generic Finalize
! ------------------
    
    call MAPL_GenericFinalize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Finalize

end module GEOS_MITDynaGridCompMod
