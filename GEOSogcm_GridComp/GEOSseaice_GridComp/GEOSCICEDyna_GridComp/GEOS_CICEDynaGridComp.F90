!  $Id$

#include "MAPL_Generic.h"

!=============================================================================

module GEOS_CICEDynaGridCompMod

!BOP

! !MODULE: GEOS_CICEDyna -- The CICE Dynamics model

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
  use ice_constants,      only: rhoi, rhos, rhow, puny, Tffresh
  use ice_state,          only: nt_Tsfc, nt_iage, nt_volpn
  use ice_communicate,    only: init_communicate
  use ice_domain,         only: init_domain_blocks    
  use ice_domain_size,    only: init_domain_size
  use ice_grid,           only: init_grid1, &
                                init_grid2, &
                                get_tarea,  &               
                                get_angle,  &               
                                get_angleT, &               
                                get_tmask  
  use ice_dyn_evp,        only: init_evp, &
                                get_ice_principal_stresses
  use ice_mechred,        only: init_ice_mechred
  use ice_transport_driver, only: init_transport
  use ice_timers,         only: init_ice_timers
  use ice_fileunits,      only: init_fileunits
  use ice_calendar,       only: init_calendar
  use ice_flux,           only: init_coupler_flux, &
                                set_atm_fluxes, &
                                get_ocn_fluxes, &
                                get_bottom_stresses, &
                                get_strtilt_term,   &
                                get_strint_term,    &
                                get_strair_term,    &
                                get_momentum_terms, &
                                get_ice_flux_vars,  &
                                get_smooth_windstress, &
                                get_ice_tendencies, &
                                get_ice_transport_flux,   &
                                get_ice_ridge_tendencies, &
                                init_stress_tensor, &
                                gather_scatter_stress, &
                                finalize_stress_tensor
  use ice_state,          only: init_trcr_depend
  use ice_step_mod,       only: step_dynamics,                & 
                                cleanup_itd_before_dynamics
  use ice_state,          only: get_ice_state, set_ice_state, &
                                set_dyn_ice_state,            &
                                set_therm_ice_state,          &
                                get_aggregate_ice_state,      &
                                get_ice_fr_thickness_tgrid,   &
                                get_snow_thickness_tgrid,     &
                                get_ice_vol_transport,        &
                                get_ice_mass_transport,       &
                                get_ice_vel_tgrid_lon_lat
  use ice_itd,            only: aggregate_ice_states
  use ice_therm_vertical, only: get_is_interface_temp
  use ice_init,           only: alloc_dyna_arrays, dealloc_dyna_arrays
  use ice_work,           only: init_work

  implicit none
  private

  integer, parameter :: FILE_HEADER_SIZE=14

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices


  integer, parameter :: NUM_3D_ICE_TRACERS=3
  integer            :: NUM_ICE_CATEGORIES
  integer            :: NUM_ICE_LAYERS
  integer, parameter :: NUM_SNOW_LAYERS=1
  integer            :: NUM_ICE_LAYERS_ALL
  integer            :: NUM_SNOW_LAYERS_ALL
  integer, parameter :: IDEB=15
  integer, parameter :: JDEB=25
  integer, parameter :: TARPE=63


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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

    call MAPL_GenericMakeXchgNatural(MAPL, RC=STATUS)
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

    call MAPL_GetPointer(INTERNAL, UVEL,        'UVEL'     , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(INTERNAL, VVEL,        'VVEL'     , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(INTERNAL, ICEUMASK  ,  'ICEUMASK' , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(INTERNAL, STRESSCOMP,  'STRESSCOMP', RC=STATUS)
    VERIFY_(STATUS) 

    call MAPL_GetResource (MAPL, CPLS, Label="CICE_COUPLE_OCEAN:", DEFAULT=1, RC=STATUS)
    VERIFY_(STATUS)
    if(MAPL_AM_I_ROOT()) then
       print*,'CPLS = ', CPLS
    endif

    call ESMF_AlarmGet(ALARM, ringTime=RING_TIME,  &
                    ringInterval=RING_INTERVAL,    &
                    RC=STATUS) 
    VERIFY_(STATUS) 

    call ESMF_ClockGet(CLOCK, currTime=currTime, RC=STATUS)
    VERIFY_(STATUS)

    ALARM_OCNEXCH = ESMF_AlarmCreate(NAME = trim(COMP_NAME)//"OCN_COUP_Alarm", &
                   CLOCK        = CLOCK,                 &
                   !RingTime     = RING_TIME,             &
                   RingTime     = currTime,             &
                   ringInterval = CPLS*RING_INTERVAL,    & 
                   RC=STATUS)
    VERIFY_(STATUS) 


!!$    if(MAPL_AM_I_ROOT()) then
!!$       print*, 'currTime = '    
!!$       call ESMF_TimePrint(currTime, "string", RC=STATUS)
!!$       print*, 'RING_TIME = '    
!!$       call ESMF_TimePrint(RING_TIME, "string", RC=STATUS)
!!$    endif

    if(RING_TIME == currTime) then
      call ESMF_AlarmRingerOn(ALARM_OCNEXCH, RC=STATUS); VERIFY_(STATUS)
    end if

    call MAPL_StateAlarmAdd(MAPL, ALARM_OCNEXCH, RC=STATUS)
    VERIFY_(STATUS)


    !call MAPL_Get (MAPL, ExchangeGrid=EXCH,   RC=STATUS )
    !VERIFY_(STATUS)

    !allocate(FROCEAN(IM,JM), STAT=STATUS)
    !VERIFY_(STATUS)
    !FROCEAN = 0.0

    !call MAPL_LocStreamFracArea( EXCH, MAPL_OCEAN, FROCEAN, RC=STATUS)
    !VERIFY_(STATUS)

    call MAPL_Get(MAPL, HEARTBEAT = DTI, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DTI, Label="CICE_DT:", DEFAULT=DTI, RC=STATUS)
    VERIFY_(STATUS)
    DTR = DTI
    if(MAPL_AM_I_ROOT()) then
       print*,'DTI = ', DTI, ' DTR = ', DTR
    endif
! Get the layout from the grid
!-----------------------------

    call ESMF_VMGetCurrent(VM, rc=STATUS)
    VERIFY_(STATUS)

! CICE grid initialization using the communicator from the VM
!------------------------------------------------------
    
    call ESMF_VMGet(VM, mpiCommunicator=Comm, petCount=NPES, rc=STATUS)
    VERIFY_(STATUS)

    ! call CICE initializing MPI communicator routine
    call init_communicate(Comm)
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI comm initialized'
    endif

    call init_fileunits       ! unit numbers
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE file units initialized'
    endif

    call MAPL_GetResource (MAPL, PROCESSOR_SHAPE,   Label="CICE_PROCESSOR_SHAPE:", DEFAULT="square-pop", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, DISTRIBUTION_TYPE, Label="CICE_DISTRIBUTION_TYPE", DEFAULT="cartesian", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, DISTRIBUTION_WGHT, Label="CICE_DISTRIBUTION_WGHT:", DEFAULT="block", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, EW_BOUNDARY_TYPE,   Label="CICE_EW_BOUNDARY_TYPE:", DEFAULT="cyclic", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, NS_BOUNDARY_TYPE,  Label="CICE_NS_BOUNDARY_TYPE:", DEFAULT="tripole", RC=STATUS)
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
       print*, 'CICE work array initialized'
    endif
 
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE DISTRIBUTION_TYPE = ', DISTRIBUTION_TYPE
       print*, 'CICE PROCESSOR_SHAPE = ', PROCESSOR_SHAPE
    endif
    call init_domain_blocks(NPES, PROCESSOR_SHAPE, DISTRIBUTION_TYPE, &
                        DISTRIBUTION_WGHT, EW_BOUNDARY_TYPE, NS_BOUNDARY_TYPE)

    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE domain blocks initialized'
    endif

    call init_ice_timers      ! initialize all timers
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE timer initialized'
    endif

    call MAPL_GetResource (MAPL, GRIDTYPE, Label="CICE_GRIDTYPE:",  DEFAULT="tripole",  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, GRIDFILE, Label="CICE_GRID:", DEFAULT="unknown_grid", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL,  KMTFILE, Label="CICE_KMT:",  DEFAULT="unknown_kmt",  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, GRIDFILEFORM, Label="CICE_GRIDFORM:", DEFAULT="bin", RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL,  KMTFILEFORM, Label="CICE_KMTFORM:",  DEFAULT="bin",  RC=STATUS)
    VERIFY_(STATUS)

    if(MAPL_AM_I_ROOT()) then
       print*, 'gridtype =', GRIDTYPE
       print*, 'gridfile =', GRIDFILE
       print*, 'kmtfile =',  KMTFILE
       print*, 'gridfileform =', GRIDFILEFORM
       print*, 'kmtfileform =', KMTFILEFORM
    endif

    call init_grid1(GRIDFILE, GRIDFILEFORM, GRIDTYPE, KMTFILE, KMTFILEFORM)
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI grid1 initialized'
    endif

#ifdef MODIFY_TOPOGRAPHY
    call MAPL_Get (MAPL, ExchangeGrid=EXCH,   RC=STATUS )
    VERIFY_(STATUS)
    allocate(FROCEAN(IM,JM), STAT=STATUS)
    VERIFY_(STATUS)
    FROCEAN = 0.0
    call MAPL_LocStreamFracArea( EXCH, MAPL_OCEAN, FROCEAN, RC=STATUS)
    VERIFY_(STATUS)
    !*** we pass frocean to init_grid2 such that some ocean grid points
    !***  as decided by tmask, will be set to land 
    call init_grid2 (FROCEAN)
    deallocate(FROCEAN)
#else
    ! note: we do NOT need to pass arguments again
    call init_grid2
#endif
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI grid2 initialized'
    endif

    !deallocate(FROCEAN)

    call MAPL_GetResource (MAPL, STRENGTH,   Label="CICE_STRENGTH:",    DEFAULT=1,       RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetResource (MAPL, RDG_PARTIC, Label="CICE_RDG_PARTIC:",  DEFAULT=1,       RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, RDG_REDIST, Label="CICE_RDG_REDIST:",  DEFAULT=1,       RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetResource (MAPL, MU_RDG,     Label="CICE_MU_RDG:",      DEFAULT=4.0,     RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (MAPL, PSTAR,      Label="CICE_PSTAR:",       DEFAULT=27500.0, RC=STATUS)
    VERIFY_(STATUS)

    call init_ice_mechred(STRENGTH, RDG_PARTIC, RDG_REDIST, MU_RDG, PSTAR)
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI ice_mechred initialized'
       print*, 'mu_rdg = ', MU_RDG
       print*, 'Pstar(Hibler79) = ', PSTAR
    endif

    call MAPL_GetResource (MAPL,  NDTE, Label="CICE_NDTE:",  DEFAULT=300,  RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetResource (MAPL,  EVP_DAMPING, Label="CICE_EVP_DAMPING:",  DEFAULT=0,  RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetResource (MAPL,  IODRAG, Label="CICE_IO_DRAG:",  DEFAULT=0.00536,  RC=STATUS)
    VERIFY_(STATUS) 
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE NDTE = ', NDTE
       print*, 'CICE EVP DAMPING ', EVP_DAMPING
       print*, 'CICE ICE-OCEAN DRAG ', IODRAG
       if(EVP_DAMPING > 0) then 
         print*, 'CICE evp damping is ENABLED'
       else
         print*, 'CICE evp damping is DISABLED'
       endif 
    endif
    ! define evp dynamics parameters, variables 
    call init_evp (DTR, NDTE, EVP_DAMPING, IODRAG)  
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI evp initialized'
    endif

    call init_coupler_flux 
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI coupler flux initialized'
    endif
    

    !call init_trcr_depend(.true., .true.)

    call MAPL_GetResource (MAPL,  ADVECTION, Label="CICE_ADVECTION:",  DEFAULT="remap",  RC=STATUS)
    VERIFY_(STATUS)
    ! initialize horizontal transport
    call init_transport(ADVECTION)
    if(MAPL_AM_I_ROOT()) then
       print*, 'CICE MPI transport initialized'
    endif

    call MAPL_GetResource (MAPL,  NDYN_DT, Label="CICE_NDYN_DT:",  DEFAULT=1,  RC=STATUS)
    VERIFY_(STATUS)
    call init_calendar(NDYN_DT)      ! initialize some calendar stuff
#if 0
    allocate(ICEUM(IM,JM),STAT=STATUS)
    VERIFY_(STATUS)
    where(ICEUMASK > 0.5)
        ICEUM = .true.
    elsewhere
        ICEUM = .false.
    endwhere 
#endif
    call init_stress_tensor(ICEUMASK, STRESSCOMP)
    ! up to now, cice stress tensors contain valid values in the interior domain;
    ! unlike other state variables, stress tensors are distributed in a special routine
    ! which also fills in the ghost values correctly
    ! we have to do gather and then a scatter to properly initialzie stress tensors
    call gather_scatter_stress
    !call set_dyn_ice_state(UVEL, VVEL)
#if 0
    deallocate(ICEUM)
#endif

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

  type (MAPL_MetaComp), pointer       :: STATE
  type (ESMF_State       )            :: INTERNAL
  type(ESMF_logical)                  :: FRIENDLY
  type(ESMF_FIELD)                    :: FIELD
  type (ESMF_Time)                    :: CurrentTime
  real                                :: DT
  real                                :: RUN_DT

  integer                             :: k, N, NDYN_DT
  integer                             :: NX, NY, I, J
  real(kind=ESMF_KIND_R8)             :: DYN_DT

  type(ESMF_Alarm)                    :: OCN_COUP 
  type(ESMF_Time)                     :: currTime 
  logical                             :: TO_CPL
  integer                             :: CPLS
         
! pointers to export

   real, pointer, dimension(:,:)  :: UI
   real, pointer, dimension(:,:)  :: VI
   real, pointer, dimension(:,:)  :: TAUXBOT
   real, pointer, dimension(:,:)  :: TAUYBOT
   real, pointer, dimension(:,:)  :: TAUXOCNB
   real, pointer, dimension(:,:)  :: TAUYOCNB
   real, pointer, dimension(:,:)  :: STROCNXB
   real, pointer, dimension(:,:)  :: STROCNYB
   real, pointer, dimension(:,:)  :: FRACICEOUT 
   real, pointer, dimension(:,:)  :: DIVUO
   real, pointer, dimension(:,:)  :: DIVUOCN
   real, pointer, dimension(:,:)  :: FAKEDIV
   real, pointer, dimension(:,:)  :: TAUADIV
   real, pointer, dimension(:,:)  :: TAU1DIV
   real, pointer, dimension(:,:)  :: TAUODIV
   real, pointer, dimension(:,:)  :: STRENGTHO
   real, pointer, dimension(:,:)  :: SHEARO
   real, pointer, dimension(:,:)  :: HICE
   real, pointer, dimension(:,:)  :: HICEUNT
   real, pointer, dimension(:,:)  :: HICE0
   real, pointer, dimension(:,:)  :: VICE0
   real, pointer, dimension(:,:)  :: VICE1
   real, pointer, dimension(:,:)  :: HSNO
   real, pointer, dimension(:,:)  :: HSNO0
   real, pointer, dimension(:,:)  :: DRAFT
   real, pointer, dimension(:,:)  :: DRAFT0
   real, pointer, dimension(:,:)  :: AICE
   real, pointer, dimension(:,:)  :: AICEU
   real, pointer, dimension(:,:)  :: DAIDTD
   real, pointer, dimension(:,:)  :: DVIDTD
   real, pointer, dimension(:,:)  :: DVIRDGDT
   real, pointer, dimension(:,:)  :: STRCORX
   real, pointer, dimension(:,:)  :: STRCORY
   real, pointer, dimension(:,:)  :: STRTLTX
   real, pointer, dimension(:,:)  :: STRTLTY
   real, pointer, dimension(:,:)  :: STRINTX
   real, pointer, dimension(:,:)  :: STRINTY
   real, pointer, dimension(:,:)  :: TAUXIOUT
   real, pointer, dimension(:,:)  :: TAUYIOUT
   real, pointer, dimension(:,:)  :: TAUXISM
   real, pointer, dimension(:,:)  :: TAUYISM
   real, pointer, dimension(:,:)  :: TAUXIB
   real, pointer, dimension(:,:)  :: TAUYIB
   real, pointer, dimension(:,:)  :: UOCN
   real, pointer, dimension(:,:)  :: VOCN
   real, pointer, dimension(:,:)  :: SSHOUT
   real, pointer, dimension(:,:)  :: SLVOUT
   real, pointer, dimension(:,:)  :: FROCEANOUT
   real, pointer, dimension(:,:)  :: AREA
   real, pointer, dimension(:,:)  :: TMSK
   real, pointer, dimension(:,:)  :: UMSK0
   real, pointer, dimension(:,:)  :: UMSK1
   real, pointer, dimension(:,:)  :: UVEL0
   real, pointer, dimension(:,:)  :: UVEL1
   real, pointer, dimension(:,:)  :: VVEL0
   real, pointer, dimension(:,:)  :: VVEL1
   real, pointer, dimension(:,:)  :: VEL
   real, pointer, dimension(:,:)  :: ANGLET
   real, pointer, dimension(:,:)  :: ANGLEU
   real, pointer, dimension(:,:)  :: TRANSIX
   real, pointer, dimension(:,:)  :: TRANSIY
   real, pointer, dimension(:,:)  :: TRANSIMX
   real, pointer, dimension(:,:)  :: TRANSIMY
   real, pointer, dimension(:,:)  :: FRESH
   real, pointer, dimension(:,:)  :: FSALT
   real, pointer, dimension(:,:)  :: FHOCN
   real, pointer, dimension(:,:)  :: VOLICE0
   real, pointer, dimension(:,:)  :: VOLSNO0
   real, pointer, dimension(:,:)  :: VOLICE1
   real, pointer, dimension(:,:)  :: VOLSNO1
   real, pointer, dimension(:,:)  :: SIG1
   real, pointer, dimension(:,:)  :: SIG2

!  pointers to category dimensional exports
   real, pointer, dimension(:,:,:)  :: AICEN
   real, pointer, dimension(:,:,:)  :: VICEN
   real, pointer, dimension(:,:,:)  :: HIFLXE
   real, pointer, dimension(:,:,:)  :: HIFLXN
   real, pointer, dimension(:,:,:)  :: STRESS0
   real, pointer, dimension(:,:,:)  :: STRESS1

   ! pointers to CMIP5 exports
   real, pointer, dimension(:,:)  :: AICE_C5
   real, pointer, dimension(:,:)  :: HICE_C5
   real, pointer, dimension(:,:)  :: HSNO_C5
   real, pointer, dimension(:,:)  :: FWM_C5
   real, pointer, dimension(:,:)  :: WIST_C5
   real, pointer, dimension(:,:)  :: WTSNINT_C5
   real, pointer, dimension(:,:)  :: TRANSIX_C5
   real, pointer, dimension(:,:)  :: TRANSIY_C5
   real, pointer, dimension(:,:)  :: UI_C5
   real, pointer, dimension(:,:)  :: VI_C5

! pointers to internal

   real(kind=ESMF_KIND_R8), pointer, dimension(:,:)  :: UVEL
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:)  :: VVEL
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:)  :: DIVU
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:)  :: SHEAR
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:)  :: STRENGTH
   real                   , pointer, dimension(:,:)  :: ICEUMASK
   real(kind=ESMF_KIND_R8), pointer, dimension(:,:,:):: STRESSCOMP

! pointers to import

   real, pointer, dimension(:,:,:)  :: TI
   real, pointer, dimension(:,:,:)  :: FR
   real, pointer, dimension(:,:,:)  :: VOLICE
   real, pointer, dimension(:,:,:)  :: VOLSNO
   real, pointer, dimension(:,:,:)  :: TAUAGE
   real, pointer, dimension(:,:,:)  :: MPOND

   real, pointer, dimension(:,:,:)  :: ERGICE
   real, pointer, dimension(:,:,:)  :: ERGSNO

   real, pointer, dimension(:,:)  :: TAUXI
   real, pointer, dimension(:,:)  :: TAUYI
   real, pointer, dimension(:,:)  :: HI
   real, pointer, dimension(:,:)  :: SI
   real, pointer, dimension(:,:)  :: UW
   real, pointer, dimension(:,:)  :: VW
   real, pointer, dimension(:,:)  :: UWB
   real, pointer, dimension(:,:)  :: VWB
   real, pointer, dimension(:,:)  :: SSH
   real, pointer, dimension(:,:)  :: SLV
   real, pointer, dimension(:,:)  :: FROCEAN
   real, pointer, dimension(:,:)  :: FRESHTHM
   real, pointer, dimension(:,:)  :: FSALTTHM
   real, pointer, dimension(:,:)  :: FHOCNTHM

   real,     allocatable          :: TRCRN(:,:,:,:) 
   real,     allocatable          :: FRALL(:,:) 
   real,     allocatable          :: TAUXBOT0(:,:) 
   real,     allocatable          :: TAUYBOT0(:,:) 
   logical,  allocatable          :: TMASK(:,:) 

   integer                        :: myPE, nDEs, steady_state_ocean = 0
   integer                        :: taui_smooth = 0
   type( ESMF_VM )                :: VMG
   real, pointer, dimension(:,:)  :: LATS
   real, pointer, dimension(:,:)  :: LONS

   real                           :: HOLD, scaling = 1.0
   real                           :: scaling_currents = 1.0

!  Begin...
!----------

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam


! Get my internal MAPL_Generic state
!----------------------------------

    call MAPL_GetObjectFromGC(GC, STATE, STATUS)
    VERIFY_(STATUS)


    call MAPL_Get(STATE,                     &
         INTERNAL_ESMF_STATE = INTERNAL,     &
         !RUNALARM = OCN_COUP,                &
         LATS  = LATS ,                      &
         LONS  = LONS ,                      &
                                RC=STATUS )
    VERIFY_(STATUS)


! Start Total timer
!------------------

   call MAPL_TimerOn(STATE,"TOTAL")
   call MAPL_TimerOn(STATE,"RUN" )

    call MAPL_GetPointer(IMPORT, TI      ,  'TI'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FR      ,  'FRACICE', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VOLICE  ,  'VOLICE'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VOLSNO  ,  'VOLSNO'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ERGICE  ,  'ERGICE'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ERGSNO  ,  'ERGSNO'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUAGE  ,  'TAUAGE'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, MPOND   ,  'MPOND'      , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(IMPORT, TAUXI   ,  'TAUX'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUYI   ,  'TAUY'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, UW      ,  'UW'         , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VW      ,  'VW'         , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, UWB     ,  'UWB'        , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, VWB     ,  'VWB'        , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SLV     ,  'SLV'        , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FROCEAN ,  'FROCEAN'    , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FRESHTHM ,  'FRESH'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FSALTTHM ,  'FSALT'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FHOCNTHM ,  'FHOCN'     , RC=STATUS)
    VERIFY_(STATUS)
   
    call MAPL_GetPointer(INTERNAL, UVEL      ,  'UVEL'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, VVEL      ,  'VVEL'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, DIVU      ,  'DIVU'     , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, SHEAR     ,  'SHEAR'    , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(INTERNAL, STRENGTH  ,  'STRENGTH' , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(INTERNAL, STRESSCOMP, 'STRESSCOMP', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(INTERNAL, ICEUMASK  ,  'ICEUMASK' , RC=STATUS)
    VERIFY_(STATUS) 

    call MAPL_GetPointer(EXPORT,      UI  , 'UI'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,      VI  , 'VI'       , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUXBOT  , 'TAUXBOT'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUYBOT  , 'TAUYBOT'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUXOCNB , 'TAUXOCNB' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUYOCNB , 'TAUYOCNB' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, STROCNXB , 'STROCNXB' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, STROCNYB , 'STROCNYB' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRACICEOUT,'FRACICE'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DIVUO    ,  'DIVU'    , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DIVUOCN  ,  'DIVUOCN' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FAKEDIV  ,  'FAKEDIV' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUADIV  ,  'TAUADIV' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAU1DIV  ,  'TAU1DIV' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUODIV  ,  'TAUODIV' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, STRENGTHO, 'STRENGTH' , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, SHEARO,    'SHEAR' , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, HICE,   'HICE'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, HICEUNT,'HICEUNT', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, HICE0,  'HICE0' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, HSNO,   'HSNO'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, HSNO0,  'HSNO0' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VICE0,  'VICE0' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VICE1,  'VICE1' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DRAFT,  'DRAFT' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DRAFT0, 'DRAFT0', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AICE,   'AICE'  , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, AICEU,  'AICEU' , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, DAIDTD, 'DAIDTD', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, DVIDTD, 'DVIDTD', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, DVIRDGDT,'DVIRDGDT',RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, STRCORX, 'STRCORX', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, STRCORY, 'STRCORY', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, STRTLTX, 'STRTLTX', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, STRTLTY, 'STRTLTY', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, STRINTX, 'STRINTX', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, STRINTY, 'STRINTY', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, TAUXIOUT,'TAUXI'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUYIOUT,'TAUYI'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUXISM,'TAUXISM' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUYISM,'TAUYISM' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUXIB, 'TAUXIB'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAUYIB, 'TAUYIB'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UOCN    , 'UOCN'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VOCN    , 'VOCN'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FROCEANOUT,'FROCEAN', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA    , 'AREA'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TMSK    , 'TMASK'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UMSK0   , 'ICEUMASK0', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UMSK1   , 'ICEUMASK1', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UVEL0   , 'UVEL0', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UVEL1   , 'UVEL1', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VVEL0   , 'VVEL0', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VVEL1   , 'VVEL1', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VEL     ,  'VEL' , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ANGLET  , 'ANGLET', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ANGLEU  , 'ANGLEU', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SSHOUT  , 'SSH'   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SLVOUT  , 'SLV'   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TRANSIX ,'TRANSIX', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, TRANSIY ,'TRANSIY', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, TRANSIMX,'TRANSIMX', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, TRANSIMY,'TRANSIMY', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, VOLICE0 ,  'VOLICE0', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VOLSNO0 ,  'VOLSNO0', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VOLICE1 ,  'VOLICE1', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VOLSNO1 ,  'VOLSNO1', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRESH   ,  'FRESH'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FSALT   ,  'FSALT'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FHOCN   ,  'FHOCN'  , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SIG1    ,  'SIG1'   , RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SIG2    ,  'SIG2'   , RC=STATUS)
    VERIFY_(STATUS)

! category dimensional exports
    call MAPL_GetPointer(EXPORT, AICEN   ,'AICEN',    RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, VICEN   ,'VICEN',    RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, HIFLXE  ,'HIFLXE',    RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, HIFLXN  ,'HIFLXN',    RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, STRESS0,'STRESSES0',    RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, STRESS1,'STRESSES1',    RC=STATUS)
    VERIFY_(STATUS) 

    ! CMIP5 exports
    call MAPL_GetPointer(EXPORT, AICE_C5,   'AICE_CMIP5'  , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, HICE_C5,   'HICE_CMIP5'  , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, HSNO_C5,   'HSNO_CMIP5'  , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, FWM_C5,    'FWM_CMIP5'   , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, WIST_C5,   'WIST_CMIP5'  , RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, WTSNINT_C5,'WTSNINT_CMIP5', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, TRANSIX_C5,'TRANSIX_CMIP5', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, TRANSIY_C5,'TRANSIY_CMIP5', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, UI_C5,     'UI_CMIP5', RC=STATUS)
    VERIFY_(STATUS) 
    call MAPL_GetPointer(EXPORT, VI_C5,     'VI_CMIP5', RC=STATUS)
    VERIFY_(STATUS) 
    
    if(associated(HICEUNT )) HICEUNT = sum(VOLICE(:,:,1:NUM_ICE_CATEGORIES), dim=3)

    if(associated(UMSK0   )) UMSK0   = ICEUMASK
    if(associated(UVEL0   )) UVEL0   = UVEL
    if(associated(VVEL0   )) VVEL0   = VVEL
    if(associated(STRESS0 )) STRESS0 = STRESSCOMP(:,:,1:12)

     NX = size(ERGICE,dim=1)
     NY = size(ERGICE,dim=2)

     call ESMF_GridCompGet( GC, VM=VMG, RC=STATUS )
     VERIFY_(STATUS)
     call ESMF_VMGet       (VMG, localpet=MYPE, petcount=nDEs,  RC=STATUS)
     VERIFY_(STATUS)


     allocate(TRCRN(NX,NY,NUM_3D_ICE_TRACERS,NUM_ICE_CATEGORIES),STAT=STATUS)
     VERIFY_(STATUS)
     allocate(FRALL(NX,NY),STAT=STATUS)
     VERIFY_(STATUS)
     allocate(TMASK(NX,NY),STAT=STATUS)
     VERIFY_(STATUS)
     TMASK = .false.  
     call  get_tmask(TMASK)

    call MAPL_GetResource (STATE,  NDYN_DT, Label="CICE_NDYN_DT:",  DEFAULT=1,  RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_Get(STATE, HEARTBEAT = DT, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource (STATE, DT, Label="CICE_DT:", DEFAULT=DT, RC=STATUS)
    VERIFY_(STATUS)
    DYN_DT = DT / real(max(1, NDYN_DT),kind=8) ! dynamics et al timestep

     if(associated(VOLICE0)) VOLICE0 = sum(VOLICE,dim=3)
     if(associated(VOLSNO0)) VOLSNO0 = sum(VOLSNO,dim=3)

     do I=1,NX
       do J=1,NY
         do N=1,NUM_ICE_CATEGORIES 
         TRCRN(I,J,nt_Tsfc,N)   = TI    (I,J,N) - MAPL_TICE
         TRCRN(I,J,nt_iage,N)   = TAUAGE(I,J,N)
         TRCRN(I,J,nt_volpn,N)  = MPOND (I,J,N)
         enddo
       enddo
     enddo

     call set_ice_state (FR(:,:,1:NUM_ICE_CATEGORIES),     &     
                         TRCRN,                            &
                         VOLICE(:,:,1:NUM_ICE_CATEGORIES), & 
                         VOLSNO(:,:,1:NUM_ICE_CATEGORIES), &
                         ERGICE(:,:,1:NUM_ICE_LAYERS_ALL), &
                         ERGSNO(:,:,1:NUM_SNOW_LAYERS_ALL),&
                         UVEL,    VVEL,   &
                         FROCEAN       )
     !call set_therm_ice_state (FR(:,:,1:NUM_ICE_CATEGORIES),     &     
     !                    TRCRN,                            &
     !                    VOLICE(:,:,1:NUM_ICE_CATEGORIES), & 
     !                    VOLSNO(:,:,1:NUM_ICE_CATEGORIES), &
     !                    ERGICE(:,:,1:NUM_ICE_LAYERS_ALL), &
     !                    ERGSNO(:,:,1:NUM_SNOW_LAYERS_ALL) &
     !                        )

     call cleanup_itd_before_dynamics
     call aggregate_ice_states

     !call ESMF_ClockGet(CLOCK, currTime=currTime, RC=STATUS)
     !VERIFY_(STATUS)
     !if(MAPL_AM_I_ROOT()) then
     !   print*, 'in RUN currTime = '    
     !   call ESMF_TimePrint(currTime, "string", RC=STATUS)
     !endif

     call MAPL_GetResource (STATE, CPLS, Label="CICE_COUPLE_OCEAN:", DEFAULT=1, RC=STATUS)
     VERIFY_(STATUS)
     !where (FROCEAN > 0.0)
     !    STRESSCOMP(:,:,14) = STRESSCOMP(:,:,14) + UWB
     !    STRESSCOMP(:,:,16) = STRESSCOMP(:,:,16) + VWB
     !    STRESSCOMP(:,:,18) = STRESSCOMP(:,:,18) + SSH
     !endwhere

     call  MAPL_StateAlarmGet(STATE, OCN_COUP, & 
                             trim(COMP_NAME)//"OCN_COUP_Alarm", &
                             RC=STATUS)
     VERIFY_(STATUS)

     TO_CPL = ESMF_AlarmIsRinging(OCN_COUP, RC=STATUS)
     VERIFY_(STATUS)

     if ( TO_CPL ) then
         call ESMF_AlarmRingerOff(OCN_COUP, RC=STATUS); VERIFY_(STATUS)
         !if(MAPL_AM_I_ROOT()) then
         !    print*, 'OCN COUP Alarm is ringing!'
         !endif

         !STRESSCOMP(:,:,13) = STRESSCOMP(:,:,14) / CPLS
         !STRESSCOMP(:,:,14) = 0.0
         !STRESSCOMP(:,:,15) = STRESSCOMP(:,:,16) / CPLS
         !STRESSCOMP(:,:,16) = 0.0
         !STRESSCOMP(:,:,17) = STRESSCOMP(:,:,18) / CPLS
         !STRESSCOMP(:,:,18) = 0.0

     end if

    call mapl_getresource(state, scaling, Label = "cice_sea_level_scaling:", default = 1.0, rc = status)
    VERIFY_(status)
    call mapl_getresource(state, scaling_currents, Label = "cice_ocean_currents_scaling:", default = 1.0, rc = status)
    VERIFY_(status)
    call mapl_getresource(state, taui_smooth, Label = "cice_wind_stress_smooth:", default = 0, rc = status)
    VERIFY_(status)
    !call set_atm_fluxes(TAUXI, TAUYI, UWB, VWB, scaling*SLV, taui_smooth > 0)
    call set_atm_fluxes(TAUXI, TAUYI, scaling_currents*UWB, scaling_currents*VWB, scaling*SLV, MAPL_UNDEF)
    !call set_atm_fluxes(TAUXI, TAUYI,         &
    !                    STRESSCOMP(:,:,13),   & 
    !                    STRESSCOMP(:,:,15),   & 
    !                    STRESSCOMP(:,:,17),   & 
    !                    FROCEAN)

    do K=1, NDYN_DT
      call mapl_getresource(state, steady_state_ocean, Label = "steady_state_ocean:", default = 0, rc = status)
      VERIFY_(status)
      !if(steady_state_ocean == 0) call step_dynamics(DYN_DT, FROCEAN, WORK1, WORK2, WORK3)
      if(steady_state_ocean == 0) call step_dynamics(DYN_DT)
    enddo

     call get_ice_state (FR(:,:,1:NUM_ICE_CATEGORIES),     &     
                         TRCRN,                            &
                         VOLICE(:,:,1:NUM_ICE_CATEGORIES), & 
                         VOLSNO(:,:,1:NUM_ICE_CATEGORIES), &
                         ERGICE(:,:,1:NUM_ICE_LAYERS_ALL), &
                         ERGSNO(:,:,1:NUM_SNOW_LAYERS_ALL),&
                         UVEL,    VVEL,   &
                         DIVU,    SHEAR,  &
                         STRENGTH, FROCEAN)

     FRALL = sum(FR, dim=3)  

     if(associated(VOLICE1)) VOLICE1 = sum(VOLICE,dim=3)
     if(associated(VOLSNO1)) VOLSNO1 = sum(VOLSNO,dim=3)

     
     if(associated(UMSK1   )) UMSK1 = ICEUMASK
     if(associated(UVEL1   )) UVEL1   = UVEL
     if(associated(VVEL1   )) VVEL1   = VVEL
     if(associated(STRESS1 )) STRESS1 = STRESSCOMP(:,:,1:12)

     
     if(associated(AICEN)) then
         AICEN = FR
         do N=1,NUM_ICE_CATEGORIES  
            where(.not. TMASK)
               AICEN(:,:,N) = MAPL_UNDEF
            endwhere 
         enddo   
     endif 

     if(associated(VICEN)) then
         VICEN = VOLICE
         do N=1,NUM_ICE_CATEGORIES  
            where(.not. TMASK)
               VICEN(:,:,N) = MAPL_UNDEF
            endwhere 
         enddo   
     endif 
      
     if(associated(FRESH)) then
         FRESH = 0.0
         call get_ocn_fluxes(fre=FRESH) 
         FRESH = FRESH + FRESHTHM 
     endif  
     if(associated(FSALT)) then
         FSALT = 0.0
         call get_ocn_fluxes(sal=FSALT) 
         FSALT = FSALT + FSALTTHM
     endif  
     if(associated(FHOCN)) then
         FHOCN = 0.0
         call get_ocn_fluxes(hea=FHOCN) 
         FHOCN = FHOCN + FHOCNTHM
     endif  

     if(associated(HICE)) then
        HICE = sum(VOLICE, dim=3)
        where(FROCEAN==0.0)
            HICE = MAPL_UNDEF
        endwhere
     endif  
     if(associated(HSNO)) then
        HSNO = sum(VOLSNO, dim=3)
        where(FROCEAN==0.0)
            HSNO = MAPL_UNDEF
        endwhere
     endif  
     if(associated(AICE)) then
        AICE = sum(FR, dim=3) 
        where(FROCEAN==0.0)
            AICE = MAPL_UNDEF
        endwhere
     endif  

     if(associated(VICE0)) then
        VICE0 = MAPL_UNDEF
        call get_ice_flux_vars (v0 = VICE0) 
     endif
     if(associated(VICE1)) then
        VICE1 = MAPL_UNDEF
        call get_ice_flux_vars (v1 = VICE1) 
     endif

     if(associated(HICE0)) then
       HICE0 = 0.0
       if(associated(HICE) .and. associated(AICE)) then
         where(AICE > 0.0)
            HICE0 = HICE / AICE     
         elsewhere
            HICE0 = MAPL_UNDEF
         endwhere 
       endif
       where(.not. TMASK)
            HICE0 = MAPL_UNDEF
       endwhere 
     endif

     if(associated(HSNO0)) then
       HSNO0 = 0.0
       if(associated(HSNO) .and. associated(AICE)) then
         where(AICE > 0.0)
            HSNO0 = HSNO / AICE     
         elsewhere
            HSNO0 = MAPL_UNDEF
         endwhere 
       endif
       where(.not. TMASK)
            HSNO0 = MAPL_UNDEF
       endwhere 
     endif

     if(associated(DRAFT)) then
        if(associated(HICE) .and. associated(HSNO)) then
           DRAFT = (rhoi * HICE + rhos * HSNO) / rhow 
        else
           DRAFT = 0. 
        endif   
        where(.not. TMASK)
           DRAFT = MAPL_UNDEF 
        endwhere 
     endif

     if(associated(DRAFT0)) then
        if(associated(HICE0) .and. associated(HSNO0)) then
           DRAFT0 = (rhoi * HICE0 + rhos * HSNO0) / rhow 
        else
           DRAFT0 = 0. 
        endif   
        where(.not. TMASK)
           DRAFT0 = MAPL_UNDEF 
        endwhere 
     endif

     if(associated(DAIDTD) .and. associated(DVIDTD)) then
        call get_ice_tendencies(DAIDTD, DVIDTD) 
        DAIDTD = DAIDTD * 8640000  ! 1 s-1 -> %  day-1
        DVIDTD = DVIDTD * 8640000  ! m s-1 -> cm day-1
     endif

     if(associated(HIFLXE) .and. associated(HIFLXN)) then
        HIFLXE = 0.0
        HIFLXN = 0.0
        call get_ice_transport_flux(HIFLXE, HIFLXN) 
        HIFLXE = HIFLXE / DT  
        HIFLXN = HIFLXN / DT  
     endif

     if(associated(DVIRDGDT)) then
        call get_ice_ridge_tendencies(DVIRDGDT, 3) 
        DVIRDGDT = DVIRDGDT * 8640000  ! m s-1 -> cm day-1
     endif

     do I=1,NX
       do J=1,NY
         do N=1,NUM_ICE_CATEGORIES 
            TI    (I,J,N)   = TRCRN(I,J,nt_Tsfc,N) + MAPL_TICE   
            TAUAGE(I,J,N)   = TRCRN(I,J,nt_iage,N)   
            MPOND (I,J,N)   = TRCRN(I,J,nt_volpn,N) 
         enddo
       enddo
     enddo

     if(associated(UI) .and. associated(VI)) then
        call get_ice_vel_tgrid_lon_lat(UI, VI, MAPL_UNDEF, .false.) 
     endif

     if(associated(VEL     )) then
        if(associated(UI) .and. associated(VI)) then
           where (TMASK)
             VEL   = sqrt(UI*UI+VI*VI)
           elsewhere
             VEL   = MAPL_Undef
           endwhere
        else
           VEL   = MAPL_Undef
        endif
     endif

     if(associated(DIVUO    )) DIVUO     = DIVU
     if(associated(STRENGTHO)) STRENGTHO = STRENGTH
     if(associated(SHEARO   )) SHEARO    = SHEAR

     if(associated(DIVUOCN    )) then
         DIVUOCN = MAPL_UNDEF
         call get_momentum_terms(divuo = DIVUOCN)
     endif

     if(associated(FAKEDIV    )) then
         FAKEDIV = MAPL_UNDEF
         call get_momentum_terms(fakedivo = FAKEDIV)
     endif

     if(associated(TAUADIV    )) then
         TAUADIV = MAPL_UNDEF
         call get_momentum_terms(tauadivo = TAUADIV)
     endif

     if(associated(TAU1DIV    )) then
         TAU1DIV = MAPL_UNDEF
         call get_momentum_terms(tau1divo = TAU1DIV)
     endif

     if(associated(TAUODIV    )) then
         TAUODIV = MAPL_UNDEF
         call get_momentum_terms(tauodivo = TAUODIV)
     endif

     if(associated(DIVUO    )) then
        where(FROCEAN==0.0)
          DIVUO = MAPL_UNDEF
        endwhere
     endif
     if(associated(STRENGTHO    )) then
        where(FROCEAN==0.0)
          STRENGTHO = MAPL_UNDEF
        endwhere
     endif
     if(associated(SHEARO    )) then
        where(FROCEAN==0.0)
          SHEARO = MAPL_UNDEF
        endwhere
     endif

     if(associated(FRACICEOUT)) then
        FRACICEOUT = sum(FR, dim=3)
        where(FROCEAN==0.0)
           FRACICEOUT = MAPL_UNDEF
        endwhere
     endif 

     if(associated(TAUXBOT) .and. associated(TAUYBOT)) then
        call get_bottom_stresses(TAUXBOT, TAUYBOT) 
     endif

     if(associated(STRTLTX) .and. associated(STRTLTY)) then
        STRTLTX  = MAPL_UNDEF
        STRTLTY  = MAPL_UNDEF
        call get_strtilt_term(STRTLTX, STRTLTY, MAPL_UNDEF) 
     endif

     if(associated(STRCORX)) then
         STRCORX  = MAPL_UNDEF
         call get_momentum_terms(strcorxo = STRCORX)
     endif

     if(associated(STRCORY)) then
         STRCORY  = MAPL_UNDEF
         call get_momentum_terms(strcoryo = STRCORY)
     endif

     if(associated(STRINTX) .and. associated(STRINTY)) then
        STRINTX = MAPL_UNDEF
        STRINTY = MAPL_UNDEF
        call get_strint_term(STRINTX, STRINTY, MAPL_UNDEF) 
     endif

     if(associated(TAUXOCNB)) then
         TAUXOCNB  = MAPL_UNDEF
         call get_momentum_terms(strocnxo = TAUXOCNB)
     endif

     if(associated(TAUYOCNB)) then
         TAUYOCNB  = MAPL_UNDEF
         call get_momentum_terms(strocnyo = TAUYOCNB)
     endif

     if(associated(STROCNXB)) then
         STROCNXB  = 0.0
         call get_momentum_terms(strocnxo = STROCNXB)
     endif

     if(associated(STROCNYB)) then
         STROCNYB  = 0.0
         call get_momentum_terms(strocnyo = STROCNYB)
     endif

     if(associated(AICEU)) then
         AICEU  = 0.0
         call get_momentum_terms(aiuo = AICEU)
     endif

     if(associated(SIG1) .and. associated(SIG2)) then
        SIG1 = MAPL_UNDEF
        SIG2 = MAPL_UNDEF
        call get_ice_principal_stresses(SIG1, SIG2)  
        where(SIG1 > MAPL_UNDEF)
           SIG1 = MAPL_UNDEF
        endwhere
        where(SIG2 > MAPL_UNDEF)
           SIG2 = MAPL_UNDEF
        endwhere
     endif

     if(associated(TAUXIOUT    )) then
        TAUXIOUT = TAUXI
     endif

     if(associated(TAUYIOUT    )) then
        TAUYIOUT = TAUYI
     endif

     if(associated(TAUXISM) .and. associated(TAUYISM)) then
        TAUXISM = MAPL_UNDEF
        TAUYISM = MAPL_UNDEF
        call get_smooth_windstress(TAUXISM, TAUYISM)
     endif

     if(associated(TAUXIB) .and. associated(TAUYIB)) then
        TAUXIB = MAPL_UNDEF
        TAUYIB = MAPL_UNDEF
        call get_strair_term(TAUXIB, TAUYIB)
     endif

     if(associated(SSHOUT    )) then
        SSHOUT = SLV
     endif

     if(associated(SLVOUT    )) then
        SLVOUT = SLV
     endif

     if(associated(UOCN    )) then
        UOCN = UWB
     endif

     if(associated(VOCN    )) then
        VOCN = VWB
     endif

     if(associated(FROCEANOUT    )) then
        FROCEANOUT = FROCEAN
     endif

     if(associated(AREA    )) then
        call get_tarea(AREA)
     endif

     if(associated(TMSK    )) then
          TMSK = 0.0
          where(TMASK)
             TMSK = 1.0
          endwhere
     endif

     if(associated(ANGLET    )) then
        ANGLET = 0.0
        call get_angleT(ANGLET)
     endif

     if(associated(ANGLEU    )) then
        ANGLEU = 0.0
        call get_angle(ANGLEU)
     endif

     if(associated(TRANSIX) .and. associated(TRANSIY)) then
        call get_ice_vol_transport(TRANSIX, TRANSIY, .false.)  
        where(TMASK .and. FRALL == 0.0)
           TRANSIX = MAPL_UNDEF 
           TRANSIY = MAPL_UNDEF
        endwhere 
        where(.not. TMASK)
           TRANSIX = MAPL_UNDEF 
           TRANSIY = MAPL_UNDEF
        endwhere 
     endif

     if(associated(TRANSIMX) .and. associated(TRANSIMY)) then
        call get_ice_mass_transport(TRANSIMX, TRANSIMY, .false.)  
        where(TMASK .and. FRALL == 0.0)
           TRANSIMX = MAPL_UNDEF
           TRANSIMY = MAPL_UNDEF
        endwhere 
        where(.not. TMASK)
           TRANSIMX = MAPL_UNDEF 
           TRANSIMY = MAPL_UNDEF
        endwhere 
     endif

     ! CMIP5 exports
     if(associated(AICE_C5)) then
        AICE_C5  = 100.0 * sum(FR, dim=3)
        where(.not. TMASK)
            AICE_C5 = MAPL_UNDEF 
        endwhere 
     endif
     if(associated(HICE_C5)) then
        HICE_C5 = sum(VOLICE, dim=3)
        where(.not. TMASK)
            HICE_C5 = MAPL_UNDEF
        endwhere 
     endif
     if(associated(HSNO_C5)) then
        HSNO_C5  = sum(VOLSNO, dim=3)
        where(.not. TMASK)
            HSNO_C5 = MAPL_UNDEF 
        endwhere 
     endif
     if(associated(FWM_C5) .and. associated(HICE) .and. associated(HSNO)) then
        FWM_C5  = rhoi * HICE + rhos * HSNO
        where(.not. TMASK)
            FWM_C5 = MAPL_UNDEF 
        endwhere 
     endif

     if(associated(WIST_C5)) then
        WIST_C5  = sum(FR(:,:,1:NUM_ICE_CATEGORIES)*TI(:,:,1:NUM_ICE_CATEGORIES), dim=3) 
        where(TMASK .and. FRALL > puny)
            WIST_C5 = WIST_C5 / FRALL
        elsewhere
            WIST_C5 = MAPL_UNDEF
        endwhere 
     endif

     if(associated(WTSNINT_C5)) then
        call get_is_interface_temp(WTSNINT_C5, MAPL_UNDEF)
        !call get_is_interface_temp(TSNINT, MAPL_UNDEF)
        !WTSNINT_C5  = sum(FR(:,:,1:NUM_ICE_CATEGORIES)*TSNINT, dim=3) 
        !WTSNINT_C5  = WTSNINT_C5 + Tffresh  ! deg C -> Kelvin 
        where(.not. TMASK)
           WTSNINT_C5 = MAPL_UNDEF
        endwhere 
     endif

     if(associated(TRANSIX_C5) .and. associated(TRANSIY_C5)) then
        call get_ice_mass_transport(TRANSIX_C5, TRANSIY_C5, .true.)  
        where(TMASK .and. FRALL == 0.0)
           TRANSIX_C5 = 0.0 
           TRANSIY_C5 = 0.0 
        endwhere 
        where(.not. TMASK)
           TRANSIX_C5 = MAPL_UNDEF 
           TRANSIY_C5 = MAPL_UNDEF
        endwhere 
     endif

     if(associated(UI_C5) .and. associated(VI_C5)) then
        call get_ice_vel_tgrid_lon_lat(UI_C5, VI_C5, MAPL_UNDEF, .true.) 
        where(TMASK .and. FRALL == 0.0)
           UI_C5 = MAPL_UNDEF
           VI_C5 = MAPL_UNDEF
        endwhere 
        where(.not. TMASK) 
           UI_C5 = MAPL_UNDEF
           VI_C5 = MAPL_UNDEF
        endwhere 
     endif


     deallocate(TRCRN)
     deallocate(FRALL)  
     deallocate(TMASK)

!  All done
!-----------

   call MAPL_TimerOff(STATE,"RUN"  )
   call MAPL_TimerOff(STATE,"TOTAL")

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
    type (ESMF_State   )                               :: INTERNAL

    real                   , pointer, dimension(:,:)   :: ICEUMASK
    real(kind=ESMF_KIND_R8), pointer, dimension(:,:,:) :: STRESSCOMP


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

     call MAPL_Get(MAPL,                               &
                  INTERNAL_ESMF_STATE = INTERNAL,     &
                  RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, STRESSCOMP, 'STRESSCOMP', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, ICEUMASK  ,  'ICEUMASK' , RC=STATUS)
    VERIFY_(STATUS)

    call gather_scatter_stress
    call finalize_stress_tensor(ICEUMASK, STRESSCOMP)

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



end module GEOS_CICEDynaGridCompMod
