!  $Id$

#include "MAPL_Generic.h"

! GEOS-5 default real kind

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

module MOM_GEOS5PlugMod

!BOP
! !MODULE: MOM5_GEOS5PlugMod -- wrapper for MOM5.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a wrapper for MOM.
! It uses ESMF AND MAPL. It has heavy dependencies on FMS and MOM.
!
! This should be built like MOM, so that its default reals
! are the same as MOM's. It may also be an adequate plug for HIM.
!
! It does not use the configuration.
! Its time step is the clocks time step.
! Each run invocation runs one time step.
!

!YV: Notes: MLD here is not MLD, but boundary layer depth (hblt) from KPP
!           Barotropic stream function (PSI) is not filled. These need to be fixed

!USES:
  use ESMF
  use MAPL
  use MAPL_ConstantsMod,        only: MAPL_TICE

! These MOM dependencies are all we are currently using.

  use constants_mod,            only: constants_init
  use diag_manager_mod,         only: diag_manager_init, diag_manager_end
  use field_manager_mod,        only: field_manager_init, field_manager_end

  use fms_mod,                  only: fms_init, fms_end
  use fms_io_mod,               only: fms_io_exit

  use mpp_domains_mod,          only: domain2d, mpp_update_domains

  use time_manager_mod,         only: set_calendar_type, time_type
  use time_manager_mod,         only: set_time, set_date
  use time_manager_mod,         only: JULIAN
  use time_manager_mod,         only: operator( + )

  use ocean_model_mod,          only: ocean_model_init, update_ocean_model, ocean_model_end, ocean_model_restart
  use ocean_types_mod,          only: ocean_public_type, ice_ocean_boundary_type
  
! MAT ocean_state_type renamed due to GNU build issue with simultaneous MOM5/MOM6 model
  use ocean_model_mod,          only: get_ocean_domain, mom5_ocean_state_type 

! mjs added these two

  use ocean_model_mod,          only: mom4_get_dimensions
  use ocean_model_mod,          only: ocean_model_data_get
  use ocean_model_mod,          only: mom4_get_latlon_UVsurf, mom4_get_UVsurfB
  use ocean_model_mod,          only: mom4_get_thickness, mom4_get_tsurf, mom4_get_ssurf
  use ocean_model_mod,          only: mom4_get_pointers_to_variables, mom4_get_streamfunction,  mom4_get_mld
  use ocean_model_mod,          only: mom4_get_prog_tracer_index, mom4_put_prog_tracer, mom4_get_prog_tracer 
  use ocean_model_mod,          only: mom4_get_diag_tracer_index, mom4_get_diag_tracer 
  use ocean_model_mod,          only: mom4_get_temperature_index, mom4_get_salinity_index, &
       mom4_get_uv, mom4_get_latlon_uv, mom4_get_density
  use ocean_model_mod,          only: mom4_get_3D_tmask, mom4_set_swheat, mom4_set_swheat_fr

! This was added for a to b; Balaji was reluctant to expose ice_grid_mod.

  use mpp_parameter_mod,          only: AGRID, SCALAR_PAIR
  use mpp_io_mod,                 only: MPP_RDONLY, MPP_NETCDF 
  use mpp_io_mod,                 only: mpp_open, mpp_close
  use fms_mod,                    only: read_data

! Nothing on the MOM side is visible through this module.

  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices
!EOP

! These are the MOM-side bulletin boards, where things are in
! MOM's precision and the B grid

  type MOM_MAPL_Type
     type(ocean_public_type)           :: Ocean
     type(ice_ocean_boundary_type)   :: Ice_ocean_boundary
  end type MOM_MAPL_Type

  type MOM_MAPLWrap_Type
     type(MOM_MAPL_Type), pointer :: Ptr
  end type MOM_MAPLWrap_Type

  logical :: DUAL_OCEAN

contains


!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the PhysicsGcm GC needs to register its
!   Initialize and Run.  It uses the MAPL_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (AGCM and OGCM) and runs their
!   respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals
    type (MAPL_MetaComp),  pointer     :: MAPL  
    integer                            :: iDUAL_OCEAN

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Get the MAPL object
! -------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, RC=STATUS )
    DUAL_OCEAN = iDUAL_OCEAN /= 0


!BOS

!  !IMPORT STATE:

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
         SHORT_NAME         = 'PS',                                &
         LONG_NAME          = 'Surface Atmospheric Pressure',      &
         UNITS              = 'Pa',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'PICE',                              &
         LONG_NAME          = 'pressure due to ice weight',        &
         UNITS              = 'Pa',                                &
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

     call MAPL_AddImportSpec(GC                     ,&
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'LWFLX'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                     &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHFLX'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                     &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'QFLUX'                   ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  ) 
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'RAIN',                              &
         LONG_NAME          = 'ocean_rainfall',&
         UNITS              = 'kg m-2 s-1',                        &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'SNOW',                              &
         LONG_NAME          = 'ocean_snowfall',&
         UNITS              = 'kg m-2 s-1',                        &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'SFLX',                              &
         LONG_NAME          = 'salt_flux_from_sea_ice_to_ocean',      &
         UNITS              = 'kg m-2 s-1',                        &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

        call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENUVR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENPAR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENUVF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                             &
        SHORT_NAME         = 'PENPAF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DRNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,&
          LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DFNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  ) 
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                               &
          SHORT_NAME          = 'TR',                                &
          LONG_NAME           = 'tracer_mixing_ratios',              &
          UNITS               = '1',                                 &
          DIMS                = MAPL_DimsHorzVert,                   &
          VLOCATION           = MAPL_VLocationCenter,                &
          DATATYPE            = MAPL_BundleItem,                     &
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                    &
          LONG_NAME          = 'river_discharge_at_ocean_points',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'DISCHARGE'                 ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  ) 
     VERIFY_(STATUS)


    call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'STROCNXB',                           &
        LONG_NAME          = 'x_stress_at_base_of_ice_weighted_by_aiu',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'STROCNYB',                           &
        LONG_NAME          = 'y_stress_at_base_of_ice_weighted_by_aiu',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )    

     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                &
          SHORT_NAME         = 'AICEU',                            &
          LONG_NAME          = 'ice_concentration_of_grid_cell_Bgrid',   &
          UNITS              = '1',                                 &
          DIMS               = MAPL_DimsHorzOnly,                   &
          VLOCATION          = MAPL_VLocationNone,                  &
          RC=STATUS  )
     VERIFY_(STATUS)

    if (dual_ocean) then
       call MAPL_AddImportSpec(GC,                                &
            SHORT_NAME         = 'DEL_TEMP',                          &
            LONG_NAME          = 'temperature correction to top level MOM (Tsst-Tmom',   &
            UNITS              = 'K',                                 &
            DIMS               = MAPL_DimsHorzOnly,                   &
            VLOCATION          = MAPL_VLocationNone,                  &
            RC=STATUS  )
       VERIFY_(STATUS)
    end if

!  !EXPORT STATE:

! Run1 exports

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'UW',                                &
         LONG_NAME          = 'surface_Agrid_eastward_velocity',   &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'VW',                                &
         LONG_NAME          = 'surface_Agrid_northward_velocity',  &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'UWB',                                &
         LONG_NAME          = 'surface_Bgrid_X_velocity',    &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'VWB',                                &
         LONG_NAME          = 'surface_Bgrid_Y_velocity',  &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'TW',                                &
         LONG_NAME          = 'surface_temperature',               &
         UNITS              = 'K',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'SW',                                &
         LONG_NAME          = 'surface_salinity',                  &
         UNITS              = 'psu',                               &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'MOM_3D_MASK',                       &
         LONG_NAME          = 'Mom4_ocean_mask_at_t-points',       &
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'AREA',                              &
         LONG_NAME          = 'Mom4_ocean_area_at_t-points',       &
         UNITS              = 'm+2',                               &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)


    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'SSH',                              &
         LONG_NAME          = 'sea_level_height',                 &
         UNITS              = 'm',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'SLV',                              &
         LONG_NAME          = 'sea_level_with_ice_loading',       &
         UNITS              = 'm',                                &
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

! Diagnostic exports


    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'DH',                                &
         LONG_NAME          = 'layer_thickness',                   &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'RHO',                               &
         LONG_NAME          = 'density',                           &
         UNITS              = 'kg m-3',                            &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'MASSCELLO',                         &
         LONG_NAME          = 'mass_per_unit_area',                &
         UNITS              = 'kg m-2',                            &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'HC',                                &
         LONG_NAME          = 'heat_content',                      &
         UNITS              = 'J m-2',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'U',                                 &
         LONG_NAME          = 'eastward_current',                  &
         UNITS              = 'm s-1',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'V',                                 &
         LONG_NAME          = 'northward_current',                 &
         UNITS              = 'm s-1',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'UX',                                &
         LONG_NAME          = 'x_current',                         &
         UNITS              = 'm s-1',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'VX',                                &
         LONG_NAME          = 'y_current',                         &
         UNITS              = 'm s-1',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'T',                                 &
         LONG_NAME          = 'potential_temperature',          &
         UNITS              = 'K',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'TCON',                            &
         LONG_NAME          = 'conservative_temperature',             &
         UNITS              = 'K',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'S',                                 &
         LONG_NAME          = 'salinity',                          &
         UNITS              = 'psu',                               &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'WMO',                               &
         LONG_NAME          = 'upward_mass_transport',             &
         UNITS              = 'tonne s-1',                         &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'WMOSQ',                           &
         LONG_NAME          = 'upward_mass_transport_squared',     &
         UNITS              = 'tonne2 s-2',                        &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'TOSSQ',                              &
         LONG_NAME          = 'surface_temperature_squared',       &
         UNITS              = 'K2',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'PBO',                          &
         LONG_NAME          = 'pressure_at_sea_floor',        &
         UNITS              = 'dbar',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'OMLDAMAX',                          &
         LONG_NAME          = 'maximum_mixed_layer_thickness',     &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME         = 'DEPTH',                        &
         LONG_NAME          = 'layer_depth',                  &
         UNITS              = 'm',                            &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'MLD',                               &
         LONG_NAME          = 'mixed_layer_depth',                 &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'PSI',                               &
         LONG_NAME          = 'barotropic_streamfunction',         &
         UNITS              = 'kg s-1',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'SWFRAC',                            &
         LONG_NAME          = 'shortwave_fractional_decay',        &
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)
    
!EOS

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	    Run,        RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_WRITERESTART, Record,     RC=status)
    VERIFY_(STATUS)
    if (dual_ocean) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	    Run2,        RC=status)
       VERIFY_(STATUS)
    end if

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="INITIALIZE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUN"        ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUN2"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="FINALIZE"   ,RC=STATUS)
    VERIFY_(STATUS)

! Generic SetServices
! -------------------

    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

! All done
! --------

    RETURN_(ESMF_SUCCESS)
  
  end subroutine SetServices


!BOP

! !IROUTINE: INITIALIZE -- Initialize method for ExternalOcean wrapper

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp),     intent(INOUT) :: GC     ! Gridded component 
    type(ESMF_State),        intent(INOUT) :: IMPORT ! Import state
    type(ESMF_State),        intent(INOUT) :: EXPORT ! Export state
    type(ESMF_Clock),        intent(INOUT) :: CLOCK  ! The clock
    integer, optional,       intent(  OUT) :: RC     ! Error code:

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    integer                                :: FMSlayout(2)
    integer                                :: counts(7)
    integer                                :: Comm
    integer                                :: isc,iec,jsc,jec
    integer                                :: IM, JM, LM
    integer                                :: IMW, JMW
    integer                                :: YEAR,MONTH,DAY,HR,MN,SC

! Locals with MOM types

    type(time_type)                        :: Time        
    type(time_type)                        :: DT 

! Locals with ESMF and MAPL types

    type(ESMF_VM)                          :: VM
    type (MAPL_MetaComp), pointer          :: MAPL 
    type(ESMF_Grid)                        :: Grid
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: boundary
    type(ocean_public_type),       pointer :: Ocean
    type(mom5_ocean_state_type),   pointer :: Ocean_State
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state 
    type(MOM_MAPLWrap_Type)                :: wrap

    integer                                :: DT_OCEAN

    REAL_, pointer                         :: TW  (:,:)
    REAL_, pointer                         :: SW  (:,:)
    REAL_, pointer                         :: OMLDAMAX  (:,:)
    REAL_, pointer                         :: DH  (:,:,:)
    REAL_, pointer                         :: AREA(:,:)
    REAL_, pointer                         :: MASK(:,:,:)

    real, allocatable                      :: Tmp3(:,:,:), Tmp2(:,:)

    REAL_, pointer, dimension(:, :)        ::  mld, psi, sea_lev, ssh, pbo
    REAL_, pointer, dimension(:, :, :)     :: TL, SL
    integer                                :: i,j

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // trim(Iam)


! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL"     )
    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Get the grid, configuration
!----------------------------

    call ESMF_GridCompGet( GC, grid=Grid,  RC=status )
    VERIFY_(STATUS)

! Get the layout from the grid
!-----------------------------

    call ESMF_VMGetCurrent(VM, rc=STATUS)
    VERIFY_(STATUS)

! Get the dimensions of the DElayout
!-----------------------------------

    call MAPL_Get(MAPL, NX=FMSlayout(1), NY=FMSlayout(2), RC=STATUS )
    VERIFY_(STATUS)

! Set the time for MOM
!---------------------

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =MN,    S =SC,  &
                                        RC=STATUS )
    VERIFY_(STATUS)

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, RC=status)
    VERIFY_(status)

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( MOM_MAPL_internal_state, stat=status )
    VERIFY_(STATUS)

    wrap%ptr => MOM_MAPL_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'MOM_MAPL_state',wrap,status )
    VERIFY_(STATUS)

    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean    => MOM_MAPL_internal_state%Ocean

! FMS initialization using the communicator from the VM
!------------------------------------------------------

    call ESMF_VMGet(VM, mpiCommunicator=Comm, rc=STATUS)
    VERIFY_(STATUS)

    call fms_init(Comm)

! Do MOM stuff
!-------------

    call constants_init
    call field_manager_init
    call diag_manager_init
    call set_calendar_type (JULIAN                )
    DT   = set_time (DT_OCEAN, 0)
    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)
    call ocean_model_init  (Ocean, Ocean_state, Time, Time)

! Check local sizes of two horizontal dimensions
!-----------------------------------------------

    call mom4_get_dimensions(isc, iec, jsc, jec, nk_out=LM)
    call MAPL_GridGet(GRID, localCellCountPerDim=counts, RC=status)
    VERIFY_(STATUS)

    IM=iec-isc+1
    JM=jec-jsc+1
    
    ASSERT_(counts(1)==IM)
    ASSERT_(counts(2)==JM)

! Allocate MOM's flux bulletin board.
!------------------------------------

    allocate ( Boundary% u_flux          (isc:iec,jsc:jec), &
               Boundary% v_flux          (isc:iec,jsc:jec), &
               Boundary% t_flux          (isc:iec,jsc:jec), &
               Boundary% q_flux          (isc:iec,jsc:jec), &
               Boundary% salt_flux       (isc:iec,jsc:jec), &
               Boundary% lw_flux         (isc:iec,jsc:jec), &
               Boundary% sw_flux_vis_dir (isc:iec,jsc:jec), &
               Boundary% sw_flux_vis_dif (isc:iec,jsc:jec), &
               Boundary% sw_flux_nir_dir (isc:iec,jsc:jec), &
               Boundary% sw_flux_nir_dif (isc:iec,jsc:jec), &
               Boundary% lprec           (isc:iec,jsc:jec), &
               Boundary% fprec           (isc:iec,jsc:jec), &
               Boundary% runoff          (isc:iec,jsc:jec), &
               Boundary% calving         (isc:iec,jsc:jec), &
               Boundary% runoff_hflx     (isc:iec,jsc:jec), &
               Boundary% calving_hflx    (isc:iec,jsc:jec), &
               Boundary% p               (isc:iec,jsc:jec), &
                                                stat=STATUS )
    VERIFY_(STATUS)

! Clear the fluxes we will not be using
!--------------------------------------

    Boundary%u_flux          = 0.0
    Boundary%v_flux          = 0.0
    Boundary%t_flux          = 0.0
    Boundary%q_flux          = 0.0
    Boundary%salt_flux       = 0.0
    Boundary%lw_flux         = 0.0
    Boundary%sw_flux_vis_dir = 0.0
    Boundary%sw_flux_vis_dif = 0.0
    Boundary%sw_flux_nir_dir = 0.0
    Boundary%sw_flux_nir_dif = 0.0
    Boundary%lprec           = 0.0
    Boundary%fprec           = 0.0
    Boundary%runoff          = 0.0
    Boundary%calving         = 0.0
    Boundary%runoff_hflx     = 0.0
    Boundary%calving_hflx    = 0.0
    Boundary%p               = 0.0

! Profilers
! ---------

    call MAPL_TimerOff(MAPL,"INITIALIZE")
    call MAPL_TimerOff(MAPL,"TOTAL"     )

! Generic initialize
! ------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

! Make sure exports neede by the parent prior to our run call are initialized
!----------------------------------------------------------------------------

    call MAPL_GetPointer(EXPORT, TW,   'TW'  ,  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SW,   'SW'  ,  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DH,   'DH'  ,  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MASK, 'MOM_3D_MASK',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA, 'AREA',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, OMLDAMAX,   'OMLDAMAX'  ,  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, mld, 'MLD', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, psi, 'PSI', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, sea_lev, 'SLV', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ssh, 'SSH', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, pbo, 'PBO', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TL, 'T',  alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SL, 'S', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)

! Get the 3-D MOM data
!---------------------

    allocate(Tmp3(IM,JM,LM),stat=status); VERIFY_(STATUS)

    call mom4_get_3D_tmask(Tmp3)
    MASK = real(Tmp3,kind=G5KIND)

    call mom4_get_thickness(Tmp3)
    where(MASK > 0.0)
       DH = real(Tmp3,kind=G5KIND)
    elsewhere
       DH = MAPL_UNDEF
    end where

    call mom4_get_salinity_index(i)
    call mom4_get_prog_tracer(i,fld=Tmp3) 
    where(MASK > 0.0)
       SL = real(Tmp3,kind=G5KIND)
    elsewhere
       SL = MAPL_UNDEF
    end where

! Convert conservative temp to potential if necessary
    i=-1
    call mom4_get_diag_tracer_index(i,'pot_temp')
    if(i .eq. -1) then
       call mom4_get_temperature_index(i)
       call mom4_get_prog_tracer(i,fld=Tmp3)
       where(MASK > 0.0)
          TL = real(Tmp3+MAPL_TICE,kind=G5KIND)
       elsewhere
          TL = MAPL_UNDEF
       end where
    else
       call mom4_get_diag_tracer(i,fld=Tmp3)
       where(MASK > 0.0)
          TL = real(Tmp3+MAPL_TICE,kind=G5KIND)
       elsewhere
          TL = MAPL_UNDEF
       end where
    end if

    deallocate(Tmp3)

! Get the 2-D MOM data
!---------------------
    allocate(Tmp2(IM,JM),stat=status); VERIFY_(STATUS)

    call mom4_get_Tsurf(Ocean,Tmp2)
    where(MASK(:,:,1) > 0.0)
       TW = real(Tmp2,kind=G5KIND)
    elsewhere
       TW = MAPL_UNDEF
    end where
    
    call mom4_get_Ssurf(Ocean,Tmp2)
    where(MASK(:,:,1) > 0.0)
       SW = real(Tmp2,kind=G5KIND)
    elsewhere
       SW = MAPL_UNDEF
    end where
 
    if(associated(area)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'area', Tmp2, isc, jsc)
       AREA = real(Tmp2,kind=G5KIND)
    end if 

    if(associated(sea_lev)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', Tmp2, isc, jsc)
       sea_lev = real(merge(tsource = tmp2, fsource = real(MAPL_UNDEF), mask = (mask(:, :, 1) > 0.0)), kind=G5KIND)
    end if 

    if(associated(ssh)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'eta_t', Tmp2, isc, jsc)
       ssh = real(merge(tsource = tmp2, fsource = real(MAPL_UNDEF), mask = (mask(:, :, 1) > 0.0)), kind=G5KIND)
    end if 

    if(associated(pbo)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'pbot_t', tmp2, isc, jsc)
       pbo = real(merge(tsource = 1.0e-04*tmp2, fsource = real(MAPL_UNDEF), mask = (mask(:, :, 1) > 0.0)),kind=G5KIND)
    end if 

    tmp2=mom4_get_mld()
    OMLDAMAX   = real(merge(tsource = tmp2, fsource = real(MAPL_UNDEF), mask = (mask(:, :, 1) > 0.0)), kind=G5KIND)

    deallocate(Tmp2)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize


!=================================================================================

!BOP

! !IROUTINE: Run  -- Run method for External Model Plug

! !INTERFACE:

  subroutine Run  ( gc, import, export, clock, rc )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: import ! Import state
    type(ESMF_State),    intent(INOUT) :: export ! Export state
    type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
    integer, optional,   intent(  OUT) :: rc     ! Error code:
    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    integer                                :: IM, JM, LM
    integer                                :: IMw, JMw
    integer                                :: I
    integer                                :: J

    integer                                :: steady_state_ocean = 0
    logical                                :: ocean_seg_start = .true.
    logical                                :: ocean_seg_end   = .true.

! Required exports

    REAL_, pointer                         :: TW  (:,:)
    REAL_, pointer                         :: SW  (:,:)
    REAL_, pointer                         :: UW  (:,:)
    REAL_, pointer                         :: VW  (:,:)
    REAL_, pointer                         :: UWB (:,:)
    REAL_, pointer                         :: VWB (:,:)
    REAL_, pointer                         :: DH  (:,:,:)
    REAL_, pointer                         :: SSH  (:,:)
    REAL_, pointer                         :: SLV  (:,:)
    REAL_, pointer                         :: FRZMLT  (:,:)
    REAL_, pointer                         :: MASK(:,:,:)
    REAL_, pointer                         :: AREA(:,:)
    REAL_, pointer                         :: DEPTH(:,:,:)

! Optional Exports

    REAL_, pointer                         :: TL  (:,:,:)
    REAL_, pointer                         :: SL  (:,:,:)
    REAL_, pointer                         :: RL  (:,:,:)
    REAL_, pointer                         :: UL  (:,:,:)
    REAL_, pointer                         :: VL  (:,:,:)
    REAL_, pointer                         :: UX  (:,:,:)
    REAL_, pointer                         :: VX  (:,:,:)
    REAL_, pointer                         :: WRHOT(:,:,:)
    REAL_, pointer                         :: WRHOTSQ(:,:,:)
    REAL_, pointer                         :: TSSQ(:,:)
    REAL_, pointer                         :: TCON(:,:,:)
    REAL_, pointer                         :: PBO(:,:)
    REAL_, pointer                         :: MASS(:,:,:)
    REAL_, pointer                         :: HC(:,:,:)
    REAL_, pointer                         :: OMLDAMAX(:,:)
    REAL_, pointer                         :: SWFRAC(:,:,:)

! Imports
    REAL_, pointer                         :: TAUX(:,:)
    REAL_, pointer                         :: TAUY(:,:)
    REAL_, pointer                         :: PS  (:,:)
    REAL_, pointer                         :: PICE(:,:)
    REAL_, pointer                         :: HEAT(:,:,:)
    REAL_, pointer                         :: TRACER(:,:,:)
    REAL_, pointer                         :: LWFLX(:,:)
    REAL_, pointer                         :: SHFLX(:,:)
    REAL_, pointer                         :: QFLUX(:,:)
    REAL_, pointer                         :: RAIN(:,:)
    REAL_, pointer                         :: SNOW(:,:)
    REAL_, pointer                         :: SFLX(:,:)
    REAL_, pointer                         :: PENUVR(:,:)
    REAL_, pointer                         :: PENPAR(:,:)
    REAL_, pointer                         :: PENUVF(:,:)
    REAL_, pointer                         :: PENPAF(:,:)
    REAL_, pointer                         :: DRNIR(:,:)
    REAL_, pointer                         :: DFNIR(:,:)
    REAL_, pointer                         :: DISCHARGE(:,:)
    REAL_, pointer                         :: STROCNXB(:,:)
    REAL_, pointer                         :: STROCNYB(:,:)
    REAL_, pointer                         :: AICEU(:,:)

! Temporaries

    real, allocatable                      :: U(:,:)
    real, allocatable                      :: V(:,:)
    real, allocatable                      :: H(:,:,:)
    real, allocatable                      :: G(:,:,:)
    real, allocatable                      :: cos_rot(:,:)
    real, allocatable                      :: sin_rot(:,:)
    real                                   :: EPSLN

    type(MAPL_MetaComp),           pointer :: MAPL 
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state 
    type(MOM_MAPLWrap_Type)                :: wrap
    type(ice_ocean_boundary_type), pointer :: boundary
    type(ocean_public_type),       pointer :: Ocean
    type(mom5_ocean_state_type),   pointer :: Ocean_State
    type(domain2d)                         :: OceanDomain
    integer                                :: isc,iec,jsc,jec
    integer                                :: isd,ied,jsd,jed

    integer                                :: YEAR,MONTH,DAY,HR,MN,SC
    type(time_type)                        :: Time        
    type(time_type)                        :: DT

    integer                                :: ii, jj
    integer                                :: cnt,l
    integer                                :: tracer_index
    real                                   :: sum, pice_scaling = 1.0
    integer :: DT_OCEAN
    real, parameter   :: CW = 3992.10322329649

    integer                                :: na
    type(ESMF_FieldBundle)                 :: TR
    type(ESMF_Field)                       :: Field
    type(ESMF_Array)                       :: Array
    character(len=13)                      :: TRNAME
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

    REAL_, pointer, dimension(:, :)        :: mld, psi

    REAL_, pointer, dimension(:,:)         :: LATS
    REAL_, pointer, dimension(:,:)         :: LONS
! Begin
!------


! Get the component's name and set-up traceback handle.
! -----------------------------------------------------
    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(status)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)


    call MAPL_Get(MAPL,                     &
         INTERNAL_ESMF_STATE = INTERNAL,     &
         LATS  = LATS ,                      &
         LONS  = LONS ,                      &
                                RC=STATUS )
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN"  )

! Get the Plug's private internal state
!--------------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR 

! Aliases to MOM types
!---------------------

    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean    => MOM_MAPL_internal_state%Ocean


    call get_ocean_domain(OceanDomain)
    call mom4_get_dimensions(isc, iec, jsc, jec, isd, ied, jsd, jed, LM)

    IM=iec-isc+1
    JM=jec-jsc+1

! Temporaries with MOM default reals
!-----------------------------------

    allocate(U(IM,JM   ), stat=STATUS); VERIFY_(STATUS)
    allocate(V(IM,JM   ), stat=STATUS); VERIFY_(STATUS)
    allocate(H(IM,JM,LM), stat=STATUS); VERIFY_(STATUS)
    allocate(G(IM,JM,LM), stat=STATUS); VERIFY_(STATUS)
    allocate(cos_rot(IM,JM   ), stat=STATUS); VERIFY_(STATUS)
    allocate(sin_rot(IM,JM   ), stat=STATUS); VERIFY_(STATUS)

! Get IMPORT pointers
!--------------------

    call MAPL_GetPointer(IMPORT, TAUX, 'TAUX'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUY, 'TAUY'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PS,   'PS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PICE, 'PICE'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, HEAT, 'SWHEAT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, LWFLX, 'LWFLX'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SHFLX, 'SHFLX'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QFLUX, 'QFLUX'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RAIN, 'RAIN'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SNOW, 'SNOW'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SFLX, 'SFLX'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENUVR, 'PENUVR'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENPAR, 'PENPAR'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENUVF, 'PENUVF'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENPAF, 'PENPAF'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DRNIR, 'DRNIR'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DFNIR, 'DFNIR'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DISCHARGE, 'DISCHARGE', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, STROCNXB, 'STROCNXB', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, STROCNYB, 'STROCNYB', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, AICEU, 'AICEU', RC=STATUS); VERIFY_(STATUS)

! Get EXPORT pointers
!--------------------

    call MAPL_GetPointer(EXPORT, UW,   'UW'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VW,   'VW'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UWB,  'UWB' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VWB,  'VWB' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TW,   'TW'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SW,   'SW'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SSH,   'SSH', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SLV,   'SLV', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZMLT,   'FRZMLT', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DEPTH, 'DEPTH', RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, DH,   'DH'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MASK, 'MOM_3D_MASK', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA, 'AREA', RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, TL,   'T'   ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RL,   'RHO' ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MASS, 'MASSCELLO' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, HC,   'HC' ,        RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SL,   'S'   ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UL,   'U'   ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VL,   'V'   ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UX,   'UX'  ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VX,   'VX'  ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WRHOT,'WMO'   ,     RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WRHOTSQ,'WMOSQ',    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TSSQ,   'TOSSQ'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PBO,    'PBO'   ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TCON, 'TCON'   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, OMLDAMAX,   'OMLDAMAX', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, mld, 'MLD', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, psi, 'PSI', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SWFRAC, 'SWFRAC', RC=STATUS); VERIFY_(STATUS)

    call mapl_getresource(mapl, pice_scaling, Label = "MOM_PICE_SCALING:", default = 1.0, rc = status); VERIFY_(status)
    boundary%P         = pice_scaling*real(PICE, kind = KIND(Boundary%P))
    boundary%lw_flux   = real(LWFLX, kind=KIND(Boundary%P)) ! Long wave flux: both positive down
    boundary%t_flux    = real(SHFLX, kind=KIND(Boundary%P)) ! Sencible heat flux: both positive up
    boundary%q_flux    = real(QFLUX, kind=KIND(Boundary%P)) ! Evaporation: both positive up
    boundary%lprec     = real(RAIN, kind=KIND(Boundary%P))  ! Liquid precipitation: both positive down
    boundary%fprec     = real(SNOW, kind=KIND(Boundary%P))  ! Frozen precipitation: both positive down
    boundary%salt_flux =-real(SFLX, kind=KIND(Boundary%P))  ! Salt flux: MOM positive up, GEOS positive down

    boundary%runoff    = real(DISCHARGE,kind=KIND(Boundary%P))

! All shortwave components are positive down  in MOM and in GEOS
!---------------------------------------------------------------
    boundary%sw_flux_vis_dir = real(PENUVR+PENPAR, kind=KIND(Boundary%P))
    boundary%sw_flux_vis_dif = real(PENUVF+PENPAF, kind=KIND(Boundary%P))
    boundary%sw_flux_nir_dir = real(DRNIR, kind=KIND(Boundary%P))
    boundary%sw_flux_nir_dif = real(DFNIR, kind=KIND(Boundary%P))

! Precision conversion
!---------------------

    H = real(HEAT, kind=KIND(H))

! KPP requires fractional shortwave decay (i.e. penetrated shortwave at T levels normalized by surface 
! shortwave flux). We compute it from surface shortwave flux and shotwave heating here.
!-----------------------------------------------------------------------------------------------------
    U=real(PENUVR+PENPAR+PENUVF+PENPAF+DRNIR+DFNIR, kind=KIND(U))
    
    where(U>0.0)
       V=1.0
       U=1.0/U
    elsewhere
       V=0.0 ! short wave fraction should be 0 when surface flux is 0
    end where
    
    do l=1,LM
       G(:,:,l)=V-0.5*H(:,:,l)*U
       V=V-H(:,:,l)*U
    end do
    G=max(0.0,G) ! this protects from tiny negatives at depth where sw heating is very small
    call mom4_set_swheat_fr(G)
    
    if(associated(SWFRAC)) then
       where(MASK > 0.0)
          SWFRAC = real(G,kind=G5KIND)
       elsewhere
          SWFRAC = MAPL_UNDEF
       end where
    end if

! Subtract surface flux from the top level heating rate, because MOM ocean_sbc adds surface sw flux to 
! Tprog(index_temp)%stf (surface temp tracer flux).
!--------------------------------------------------
    U=real(PENUVR+PENPAR+PENUVF+PENPAF+DRNIR+DFNIR, kind=KIND(U))
    H(:,:,1) = H(:,:,1)-U
    call mom4_set_swheat(H)

! Convert input stresses over water to B grid
!--------------------------------------------
    U = 0.0
    V = 0.0
    call transformA2B(real(TAUX,kind=kind(U)), real(TAUY,kind=kind(V)), U, V)

! Rotate input stress over water along i,j of tripolar grid, and combine with stress under ice 
!---------------------------------------------------------------------------------------------
    call ocean_model_data_get(Ocean_State, Ocean, 'cos_rot', cos_rot, isc, jsc)
    call ocean_model_data_get(Ocean_State, Ocean, 'sin_rot', sin_rot, isc, jsc)

    Boundary%U_flux = (U*cos_rot + V*sin_rot)*(1-AICEU) - STROCNXB
    Boundary%V_flux = (-U*sin_rot + V*cos_rot)*(1-AICEU) - STROCNYB

! Set the time for MOM
!---------------------

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =MN,    S =SC,  &
                                        RC=STATUS )
    VERIFY_(STATUS)

    CALL ESMF_TimeIntervalGet(TINT, S=DT_OCEAN, RC=status)
    VERIFY_(status)

    DT   = set_time (DT_OCEAN, 0)
    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)

! Copy tracers from IMPORT bundle to MOM internal state
!------------------------------------------------------

    call ESMF_StateGet(IMPORT, 'TR', TR, RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_FieldBundleGet(TR,FieldCount=NA, RC=STATUS)
    VERIFY_(STATUS)

    TRNAME = 'GENtracer_XXX'

    do I=1,NA
       call ESMF_FieldBundleGet(TR,   I,   FIELD,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet(field=field, array=Array,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ArrayGet(Array, farrayPtr=Tracer,  RC=STATUS)
       VERIFY_(STATUS)

       H = real(TRACER, kind=KIND(U))

       write(TRNAME(11:13),'(I3.3)') I
       call mom4_get_prog_tracer_index(tracer_index,TRNAME)
       if(tracer_index > 0) call mom4_put_prog_tracer(tracer_index,H)
    end do

! Run MOM for one time step
!--------------------------

    call mapl_getresource(mapl, steady_state_ocean, Label = "steady_state_ocean:", default = 0, rc = status); VERIFY_(status)
    if(steady_state_ocean == 0) call update_ocean_model(Boundary, Ocean_State, Ocean, Time, DT)
     
! Copy tracers from MOM internal state to IMPORT bundle
!------------------------------------------------------

    do I=1,NA
       call ESMF_FieldBundleGet(TR,   I,   FIELD,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_FieldGet(field=field, array=Array,  RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ArrayGet(Array, farrayPtr=Tracer,  RC=STATUS)
       VERIFY_(STATUS)

       write(TRNAME(11:13),'(I3.3)') I
       call mom4_get_prog_tracer_index(tracer_index,TRNAME)
       if(tracer_index > 0) call mom4_get_prog_tracer(tracer_index,H)

       where(MASK > 0.0)
          TRACER = real(H, kind=KIND(TRACER))
       elsewhere
          TRACER = MAPL_UNDEF
       end where
       
    end do


! Get export fields

! Required Exports at G5 precision
!---------------------------------

    call mom4_get_Tsurf(Ocean,U)
    where(MASK(:,:,1) > 0.0)
       TW = real(U, kind=G5KIND)
    elsewhere
       TW = MAPL_UNDEF
    end where

    call mom4_get_Ssurf(Ocean,U)
    where(MASK(:,:,1) > 0.0)
       SW = real(U, kind=G5KIND)
    elsewhere
       SW = MAPL_UNDEF
    end where

    call mom4_get_salinity_index(i)
    call mom4_get_prog_tracer(i,fld=H) 
    where(MASK > 0.0)
       SL = real(H,kind=G5KIND)
    elsewhere
       SL = MAPL_UNDEF
    end where
 
! Convert conservative temp to potential if necessary
    i=-1
    call mom4_get_diag_tracer_index(i,'pot_temp')
    if(i .eq. -1) then
       call mom4_get_temperature_index(i)
       call mom4_get_prog_tracer(i,fld=H)
       where(MASK > 0.0)
          TL = real(H+MAPL_TICE,kind=G5KIND)
       elsewhere
          TL = MAPL_UNDEF
       end where
    else
       call mom4_get_diag_tracer(i,fld=H)
       where(MASK > 0.0)
          TL = real(H+MAPL_TICE,kind=G5KIND)
       elsewhere
          TL = MAPL_UNDEF
       end where
    end if

    call mom4_get_thickness(H)
    where(MASK > 0.0)
       DH = real(H, kind=G5KIND)
    end where

! Optional Exports at G5 precision
!---------------------------------
! Get the A grid currents at MOM precision
!-----------------------------------------
    if(associated(UW) .or. associated(VW)) then
       call mom4_get_latlon_UVsurf(OCEAN, U, V, STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(UW  )) then
       where(MASK(:,:,1) > 0.0)
          UW = real(U, kind=G5KIND)
       elsewhere
          UW=0.0
       end where
     endif

    if(associated(VW  )) then
       where(MASK(:,:,1) > 0.0)
          VW = real(V, kind=G5KIND)
       elsewhere
          VW=0.0
       end where       
    end if

! Get the B grid currents at MOM precision
!-----------------------------------------
    if(associated(UWB) .or. associated(VWB)) then
       call mom4_get_UVsurfB(OCEAN, U, V, STATUS)
       VERIFY_(STATUS)
    endif

    if(associated(UWB  )) then
       where(MASK(:,:,1)>0.0)
          UWB = real(U, kind=G5KIND)
       elsewhere
          UWB=0.0
       end where
     endif

    if(associated(VWB  )) then
       where(MASK(:,:,1)>0.0)
          VWB = real(V, kind=G5KIND)
       elsewhere
          VWB=0.0
       end where       
    end if

    if(associated(RL  )) then
       call mom4_get_density  (H)
       where(MASK > 0.0)
          RL = real(H,kind=G5KIND)
       elsewhere
          RL = MAPL_UNDEF
       end where
    end if

    if(associated(MASS  )) then
       call mom4_get_density  (H)
       where(MASK > 0.0)
          MASS = DH*real(H,kind=G5KIND)
       elsewhere
          MASS = MAPL_UNDEF
       end where
    end if

    if(associated(HC  )) then
       call mom4_get_density  (H)
       call mom4_get_temperature_index(i)
       call mom4_get_prog_tracer(i,fld=G)
       where(MASK > 0.0)
          HC = CW*DH*real(H*(G+MAPL_TICE),kind=G5KIND)
       elsewhere
          HC = MAPL_UNDEF
       end where
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   if(associated(TCON  )) then
! Convert potential temp to conservativ if necessary
       i=-1
       call mom4_get_diag_tracer_index(i,'con_temp')
       if(i .eq. -1) then
          call mom4_get_temperature_index(i)
          call mom4_get_prog_tracer(i,fld=H)
          where(MASK > 0.0)
             TCON = real(H+MAPL_TICE,kind=G5KIND)
          elsewhere
             TCON = MAPL_UNDEF
          end where
       else
          call mom4_get_diag_tracer(i,fld=H)
          where(MASK > 0.0)
             TCON = real(H+MAPL_TICE,kind=G5KIND)
          elsewhere
             TCON = MAPL_UNDEF
          end where
       end if
    end if

    if(associated(UL) .or. associated(VL)) then
       call mom4_get_latlon_UV(G, H, STATUS)
       VERIFY_(STATUS)

       if(associated(UL  )) then
          where(MASK > 0.0)
             UL = real(G,kind=G5KIND)
          elsewhere
             UL = MAPL_UNDEF
          end where
       end if

       if(associated(VL  )) then
          where(MASK > 0.0)
             VL = real(H,kind=G5KIND)
          elsewhere
             VL = MAPL_UNDEF
          end where
       end if
    endif

   if(associated(UX) .or. associated(VX)) then
       call mom4_get_UV(G, H, STATUS)
       VERIFY_(STATUS)

       if(associated(UX  )) then
          where(MASK > 0.0)
             UX = real(G,kind=G5KIND)
          elsewhere
             UX = MAPL_UNDEF
          end where
       end if

       if(associated(VX  )) then
          where(MASK > 0.0)
             VX = real(H,kind=G5KIND)
          elsewhere
             VX = MAPL_UNDEF
          end where
       end if
    endif

    if(associated(WRHOT)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'wrhot', H, isc, jsc) 
       do l=1,LM
          where(MASK(:,:,l) > 0.0)
             WRHOT(:,:,l) = real(H(:,:,l), kind = G5KIND)
             WRHOT(:,:,l) = WRHOT(:,:,l)*AREA*1.e-3
          elsewhere
             WRHOT(:,:,l) = MAPL_UNDEF
          end where
       end do
    end if

    if(associated(WRHOTSQ)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'wrhot', H, isc, jsc) 
       do l=1,LM
          where(MASK(:,:,l) > 0.0)
             WRHOTSQ(:,:,l) = real(H(:,:,l), kind = G5KIND)
             WRHOTSQ(:,:,l) = WRHOTSQ(:,:,l)*AREA*1e-3
             WRHOTSQ(:,:,l) = WRHOTSQ(:,:,l)*WRHOTSQ(:,:,l)
          elsewhere
             WRHOTSQ(:,:,l) = MAPL_UNDEF
          end where
       end do
    end if

    if(associated(SSH)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'eta_t', U, isc, jsc) 
       where(MASK(:,:,1) > 0.0)
          SSH = real(U, kind = G5KIND)
       elsewhere
          SSH = MAPL_UNDEF
       end where
    end if

    if(associated(SLV)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', U, isc, jsc) 
       where(MASK(:,:,1)>0.0)
          SLV = real(U, kind = G5KIND)
       elsewhere
          SLV=0.0
       end where       
    end if

    if(associated(FRZMLT)) then
       ! frazil in mom5 already contains melt potential
       call ocean_model_data_get(Ocean_State, Ocean, 'frazil', U, isc, jsc) 
       where(MASK(:,:,1)>0.0)
          FRZMLT = real(U, kind = G5KIND)
       elsewhere
          FRZMLT = 0.0
       end where       
    end if
    
    if(associated(PBO)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'pbot_t', U, isc, jsc) 
       where(MASK(:,:,1) > 0.0)
          PBO = real(U*1e-4, kind = G5KIND)
       elsewhere
          PBO = MAPL_UNDEF
       end where
    end if

    if(associated(DEPTH)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'geodepth_zt', H, isc, jsc) 
          DEPTH = real(H, kind = G5KIND)
    end if

    if(associated(TSSQ)) then
       where(MASK(:,:,1) > 0.0)
          TSSQ = TW*TW
       elsewhere
          TSSQ = MAPL_UNDEF
       end where
    end if

    U=mom4_get_mld()
    if(HR==0 .and.  MN==0 .and. SC==0) then
       OMLDAMAX = real(U, kind = G5KIND)
    else
       OMLDAMAX  = max(OMLDAMAX , real(U, kind = G5KIND))
    endif

    if(associated(mld)) then
      mld = mask(:, :, 1)*u
    endif
    if(associated(psi)) then
      psi = mask(:, :, 1)*mom4_get_streamfunction()
    endif

    deallocate(H)
    deallocate(G)
    deallocate(U,V)
    deallocate(cos_rot,sin_rot)

    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

! All Done
!---------
    RETURN_(ESMF_SUCCESS)
  contains

    subroutine transformA2B(U, V, uvx, uvy)
      real, dimension(:,:) :: U, V
      real                 :: uvx(isc:,jsc:)
      real                 :: uvy(isc:,jsc:)

      integer              :: i, j
      real, allocatable    :: TX(:,:), TY(:,:)
      
      
      allocate(tx(isd:ied,jsd:jed), stat=STATUS); VERIFY_(STATUS)
      allocate(ty(isd:ied,jsd:jed), stat=STATUS); VERIFY_(STATUS)

      tx(isc:iec, jsc:jec) = U
      ty(isc:iec, jsc:jec) = V

      call mpp_update_domains(tx, ty, OceanDomain, gridtype=AGRID, flags=SCALAR_PAIR)

      do j = jsc, jec
         do i = isc, iec
            sum = 0.0
            cnt = 0
            do ii = 0,1
               do jj = 0,1
                  if (tx(i+ii,j+jj) /= MAPL_Undef) then
                     cnt = cnt+1
                     sum = sum + tx(i+ii,j+jj)
                  end if
               end do
            end do
            if (cnt /= 0) then
               uvx(i,j) = sum/real(cnt)
            else
               uvx(i,j) = 0.0
            end if

            sum = 0.0
            cnt = 0
            do ii = 0,1
               do jj = 0,1
                  if (ty(i+ii,j+jj) /= MAPL_Undef) then
                     cnt = cnt+1
                     sum = sum + ty(i+ii,j+jj)
                  end if
               end do
            end do
            if (cnt /= 0) then
               uvy(i,j) = sum/real(cnt)
            else
               uvy(i,j) = 0.0
            end if
         enddo
      enddo
      deallocate(ty, tx)
    end subroutine transformA2B
  end subroutine Run

!=================================================================================

!BOP

! !IROUTINE: Run2  -- Run2 method, needed only when in dual_ocean mode. Apply correction to top-level MOM temperature, based on DEL_TEMP

! !INTERFACE:

  subroutine Run2  ( gc, import, export, clock, rc )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
    type(ESMF_State),    intent(INOUT) :: import ! Import state
    type(ESMF_State),    intent(INOUT) :: export ! Export state
    type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
    integer, optional,   intent(  OUT) :: rc     ! Error code:
    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				   :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    integer                                :: IM, JM, LM
    integer                                :: tracer_index


! Imports
    REAL_, pointer                         :: DEL_TEMP(:,:)

! Temporaries

    real, allocatable                      :: T(:,:,:)

! Pointers to export    
    REAL_, pointer                         :: MASK(:,:,:)

    type(MAPL_MetaComp),           pointer :: MAPL 
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state 
    type(MOM_MAPLWrap_Type)                :: wrap
!    type(ice_ocean_boundary_type), pointer :: boundary
!    type(ocean_public_type),       pointer :: Ocean
!    type(mom5_ocean_state_type),   pointer :: Ocean_State
!    type(domain2d)                         :: OceanDomain
    integer                                :: isc,iec,jsc,jec
    integer                                :: isd,ied,jsd,jed

! Begin
!------


! Get the component's name and set-up traceback handle.
! -----------------------------------------------------
    Iam = "Run2"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(status)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)


! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN2"  )

! Get the Plug's private internal state
!--------------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR 

! Aliases to MOM types
!---------------------

!    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
!    Ocean    => MOM_MAPL_internal_state%Ocean


!    call get_ocean_domain(OceanDomain)
    call mom4_get_dimensions(isc, iec, jsc, jec, isd, ied, jsd, jed, LM)

    IM=iec-isc+1
    JM=jec-jsc+1

! Temporaries with MOM default reals
!-----------------------------------

    allocate(T(IM,JM,LM), stat=STATUS); VERIFY_(STATUS)

! Get IMPORT pointers
!--------------------

    call MAPL_GetPointer(IMPORT, DEL_TEMP, 'DEL_TEMP', RC=STATUS); VERIFY_(STATUS)

! Get EXPORT pointers
!--------------------
    ! by now this should be allocated, so 'alloc=.true.' is needed
    call MAPL_GetPointer(EXPORT, MASK, 'MOM_3D_MASK', RC=STATUS)
    VERIFY_(STATUS)


    call mom4_get_temperature_index(tracer_index)
    ASSERT_(tracer_index > 0) ! temperature index is valid
    call mom4_get_prog_tracer(tracer_index,fld=T)

    where(MASK(:,:,1) > 0.0) ! correct only ocean points
       !ALT: Note that we modify only top level of T
       !     we do not need to worry about temperature units
       !     since we are applying difference

       !     some relaxation ??? here or in guest ???

       T(:,:,1) = T(:,:,1) + DEL_TEMP

    end where
    call mom4_put_prog_tracer(tracer_index,T)

    deallocate(T)

    call MAPL_TimerOff(MAPL,"RUN2"  )
    call MAPL_TimerOff(MAPL,"TOTAL")

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Run2

!BOP
    
! !IROUTINE: Finalize        -- Finalize method for GuestOcean wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer:: MAPL 
    type(ESMF_Time)                  :: MyTime
    type(MOM_MAPL_Type),     pointer :: MOM_MAPL_internal_state 
    type(MOM_MAPLWrap_Type)          :: wrap
    type(ocean_public_type),     pointer :: Ocean
    type(mom5_ocean_state_type), pointer :: Ocean_State

! ErrLog Variables

    character(len=ESMF_MAXSTR)	     :: IAm
    integer			     :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals with MOM types

    type(time_type)                  :: Time        
    integer                          :: YEAR,MONTH,DAY,HR,MN,SC

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

! Get the Plug's private internal state
!--------------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR 

    Ocean => MOM_MAPL_internal_state%Ocean

! Set the times  for MOM
!----------------------

    call ESMF_ClockGet( CLOCK, currTime=MyTime, RC=STATUS)
    VERIFY_(status)
    
    call ESMF_TimeGet (MyTime,                    &
         YY=YEAR, MM=MONTH, DD=DAY, &
         H=HR,    M =MN,    S =SC,  &
         RC=STATUS )
    VERIFY_(STATUS)

    Time = set_date(YEAR,MONTH,DAY,HR,MN,SC)

    call ocean_model_end (Ocean, Ocean_State, Time)
    call diag_manager_end(Time )
    call field_manager_end

    call fms_io_exit
!    call fms_end

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


!====================================================================

! !IROUTINE: Record        -- Record method for GuestOcean wrapper (write intermediate restarts)

! !INTERFACE:

  subroutine Record ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer    :: MAPL 
    type(MOM_MAPL_Type),     pointer :: MOM_MAPL_internal_state 
    type(MOM_MAPLWrap_Type)          :: wrap
    type(mom5_ocean_state_type),  pointer :: Ocean_State

! ErrLog Variables

    character(len=ESMF_MAXSTR)	     :: IAm
    integer			     :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals
    character(len=14)                :: timeStamp
    logical                          :: doRecord

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Record"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL")

    doRecord = MAPL_RecordAlarmIsRinging(MAPL, MODE=MAPL_Write2Disk, RC=status)
    VERIFY_(STATUS)

    if (doRecord) then

! Get the Plug's private internal state
!--------------------------------------

       CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
       VERIFY_(STATUS)

       MOM_MAPL_internal_state => WRAP%PTR 

       call MAPL_DateStampGet(clock, timeStamp, rc=status)
       VERIFY_(STATUS)

       call ocean_model_restart (Ocean_State, timeStamp)
        VERIFY_(STATUS)

    end if

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

!====================================================================

end module MOM_GEOS5PlugMod

subroutine SetServices(gc, rc)
   use ESMF
   use MOM_GEOS5PlugMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc,rc=rc)
end subroutine
