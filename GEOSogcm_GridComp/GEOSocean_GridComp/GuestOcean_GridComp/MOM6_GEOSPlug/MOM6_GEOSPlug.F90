!  $Id$

#include "MAPL_Generic.h"

! GEOS   default real kind

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

module MOM6_GEOSPlugMod

!BOP
! !MODULE: MOM6_GEOSPlugMod -- to couple with MOM6.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a coupler for MOM.
! It uses ESMF AND MAPL. It has heavy dependencies on FMS and MOM.
!
! This should be built like MOM, so that its default reals
! are the same as for MOM. It may also be an adequate plug for it.
!
! It does not use the configuration.
! Its time step is the clocks time step.
! Each run invocation runs one time step.
!

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

  use mpp_domains_mod,          only: domain2d, mpp_update_domains, &
                                      mpp_get_compute_domain,       &
                                      mpp_get_data_domain

  use mpp_parameter_mod,        only: AGRID, BGRID_NE, CGRID_NE, SCALAR_PAIR

  use time_manager_mod,         only: set_calendar_type, time_type
  use time_manager_mod,         only: set_time, set_date
  use time_manager_mod,         only: JULIAN

  use ocean_model_mod,          only: ocean_model_init,     &
                                      ocean_model_init_sfc, &
                                      update_ocean_model,   &
                                      ocean_model_end,      &
                                      ocean_model_restart

  use ocean_model_mod,          only: ocean_model_data_get, &
                                      ocean_public_type,    &
                                      ocean_state_type,     &
                                      ocean_model_get_UV_surf

  use MOM_surface_forcing_gfdl, only: ice_ocean_boundary_type

  use ocean_model_mod,          only: get_ocean_grid
  use MOM_grid,                 only: ocean_grid_type
  use MOM_domains,              only: pass_vector

! Nothing on the MOM side is visible through this module.

  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices
!EOP

! These are the MOM-side bulletin boards, where things are in
! MOM precision and on its grid

  type MOM_MAPL_Type
   type(ocean_public_type),       pointer   :: Ocean
   type(ice_ocean_boundary_type), pointer   :: Ice_ocean_boundary
   type(ocean_state_type),        pointer   :: Ocean_state
  end type MOM_MAPL_Type

! A wrapper-derived data type to connect our internal state with MOM
  type MOM_MAPLWrap_Type
     type(MOM_MAPL_Type), pointer :: Ptr
  end type MOM_MAPLWrap_Type

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the PhysicsGCM GC needs to register its
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

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

!BOS

!  !IMPORT STATE:

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'TAUX',                              &
         LONG_NAME          = 'Agrid_eastward_stress_on_ocean',    &
         UNITS              = 'N m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'TAUY',                              &
         LONG_NAME          = 'Agrid_northward_stress_on_ocean',   &
         UNITS              = 'N m-2',                             &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'PS',                                &
         LONG_NAME          = 'Surface Atmospheric Pressure',      &
         UNITS              = 'Pa',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'PICE',                              &
         LONG_NAME          = 'pressure due to ice weight',        &
         UNITS              = 'Pa',                                &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'SWHEAT',                            &
         LONG_NAME          = 'solar_heating_rate',                &
         UNITS              = 'W m-2',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                     ,             &
        LONG_NAME          = 'surface_net_downward_longwave_flux',&
        UNITS              = 'W m-2'                     ,        &
        SHORT_NAME         = 'LWFLX'                   ,          &
        DIMS               = MAPL_DimsHorzOnly           ,        &
        VLOCATION          = MAPL_VLocationNone          ,        &
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                          &
        LONG_NAME          = 'upward_sensible_heat_flux' ,&
        UNITS              = 'W m-2'                     ,&
        SHORT_NAME         = 'SHFLX'                     ,&
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                          &
        LONG_NAME          = 'evaporation'               ,&
        UNITS              = 'kg m-2 s-1'                ,&
        SHORT_NAME         = 'QFLUX'                   ,  &
        DIMS               = MAPL_DimsHorzOnly           ,&
        VLOCATION          = MAPL_VLocationNone          ,&
        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'RAIN',                              &
         LONG_NAME          = 'ocean_rainfall',                    &
         UNITS              = 'kg m-2 s-1',                        &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'SNOW',                              &
         LONG_NAME          = 'ocean_snowfall',                    &
         UNITS              = 'kg m-2 s-1',                        &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
         SHORT_NAME         = 'SFLX',                              &
         LONG_NAME          = 'salt_flux_from_sea_ice_to_ocean',   &
         UNITS              = 'kg m-2 s-1',                        &
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

    call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PENPAR',                            &
        LONG_NAME          = 'net_downward_penetrating_direct_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PENUVF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_UV_flux',  &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'PENPAF',                            &
        LONG_NAME          = 'net_downward_penetrating_diffuse_PAR_flux', &
        UNITS              = 'W m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC                         ,     &
          LONG_NAME          = 'net_surface_downwelling_nir_beam_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DRNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
          RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC                         ,     &
          LONG_NAME          = 'net_surface_downwelling_nir_diffuse_flux',&
          UNITS              = 'W m-2'                       ,&
          SHORT_NAME         = 'DFNIR'                       ,&
          DIMS               = MAPL_DimsHorzOnly             ,&
          VLOCATION          = MAPL_VLocationNone            ,&
                                                  RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                            &
          LONG_NAME          = 'river_discharge_at_ocean_points',&
          UNITS              = 'kg m-2 s-1'                ,&
          SHORT_NAME         = 'DISCHARGE'                 ,&
          DIMS               = MAPL_DimsHorzOnly           ,&
          VLOCATION          = MAPL_VLocationNone          ,&
          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
          SHORT_NAME         = 'AICE',                             &
          LONG_NAME          = 'ice_concentration_of_grid_cell',   &
          UNITS              = '1',                                &
          DIMS               = MAPL_DimsHorzOnly,                  &
          VLOCATION          = MAPL_VLocationNone,                 &
          RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
        SHORT_NAME         = 'TAUXBOT',                            &
        LONG_NAME          = 'eastward_stress_at_base_of_ice_Agrid',    &
        UNITS              = 'N m-2',                              &
        DIMS               = MAPL_DimsHorzOnly,                    &
        VLOCATION          = MAPL_VLocationNone,                   &
        RC=STATUS  )
     VERIFY_(STATUS)

    call MAPL_AddImportSpec(GC,                                    &
        SHORT_NAME         = 'TAUYBOT',                            &
        LONG_NAME          = 'northward_stress_at_base_of_ice_Agrid',    &
        UNITS              = 'N m-2',                              &
        DIMS               = MAPL_DimsHorzOnly,                    &
        VLOCATION          = MAPL_VLocationNone,                   &
        RC=STATUS  )
     VERIFY_(STATUS)

!  !EXPORT STATE:

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'UW',                                &
         LONG_NAME          = 'surface_Agrid_eastward_velocity',   &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'VW',                                &
         LONG_NAME          = 'surface_Agrid_northward_velocity',  &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'UWB',                               &
         LONG_NAME          = 'surface_Bgrid_X_velocity',          &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'VWB',                               &
         LONG_NAME          = 'surface_Bgrid_Y_velocity',          &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'TW',                                &
         LONG_NAME          = 'surface_temperature',               &
         UNITS              = 'K',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'SW',                                &
         LONG_NAME          = 'surface_salinity',                  &
         UNITS              = 'psu',                               &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'MOM_2D_MASK',                       &
         LONG_NAME          = 'MOM_ocean_mask_at_t-points',        &
         UNITS              = '1',                                 &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'AREA',                              &
         LONG_NAME          = 'MOM_ocean_area_at_t-points',        &
         UNITS              = 'm+2',                               &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'SLV',                              &
         LONG_NAME          = 'sea_level_with_ice_loading_and_invBaro',       &
         UNITS              = 'm',                                &
         DIMS               = MAPL_DimsHorzOnly,                  &
         VLOCATION          = MAPL_VLocationNone,                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'FRAZIL',                           &
         LONG_NAME          = 'heating_from_frazil_formation',    &
         UNITS              = 'W m-2',                            &
         DIMS               = MAPL_DimsHorzOnly,                  &
         VLOCATION          = MAPL_VLocationNone,                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'MELT_POT',                         &
         LONG_NAME          = 'heat_that_can_be_used_to_melt_sea_ice',    &
         UNITS              = 'W m-2',                            &
         DIMS               = MAPL_DimsHorzOnly,                  &
         VLOCATION          = MAPL_VLocationNone,                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                   &
         SHORT_NAME         = 'FRZMLT',                           &
         LONG_NAME          = 'freeze_melt_potential',            &
         UNITS              = 'W m-2',                            &
         DIMS               = MAPL_DimsHorzOnly,                  &
         VLOCATION          = MAPL_VLocationNone,                 &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'DUM1',                              &
         LONG_NAME          = 'dummy_export1',        &
         UNITS              = '1',                               &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)
    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'DUM2',                              &
         LONG_NAME          = 'dummy_export2',        &
         UNITS              = '1',                               &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         RC=STATUS  )
    VERIFY_(STATUS)

!  !Diagnostic exports
!Get rid of following 3D exports

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'U',                                 &
         LONG_NAME          = 'eastward_current',                  &
         UNITS              = 'm s-1',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'V',                                 &
         LONG_NAME          = 'northward_current',                 &
         UNITS              = 'm s-1',                             &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'DH',                                &
         LONG_NAME          = 'layer_thickness',                   &
         UNITS              = 'm OR kg m-2',                       &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    & 
         SHORT_NAME         = 'DEPTH',                             &
         LONG_NAME          = 'layer_depth',                       &
         UNITS              = 'm',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'T',                                 &
         LONG_NAME          = 'potential_temperature',             &
         UNITS              = 'K',                                 &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'S',                                 &
         LONG_NAME          = 'salinity',                          &
         UNITS              = 'psu',                               &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
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
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,        RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   RC=status)
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_WRITERESTART, Record,     RC=status)
    VERIFY_(STATUS)

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="INITIALIZE" ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC,   name="RUN"        ,RC=STATUS)
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

!=============================================================================

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

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals

    integer                                :: counts(7)
    integer                                :: Comm
    integer                                :: isc,iec,jsc,jec
    integer                                :: isd,ied,jsd,jed
    integer                                :: IM, JM
    integer                                :: g_isc,g_iec,g_jsc,g_jec
    integer                                :: g_isd,g_ied,g_jsd,g_jed

    integer                                :: YEAR,MONTH,DAY,HR,MN,SC

! Locals with MOM types

    type(time_type)                        :: Time
    type(time_type)                        :: DT

! Locals with ESMF and MAPL types

    type(ESMF_VM)                          :: VM
    type(MAPL_MetaComp), pointer           :: MAPL
    type(ESMF_Grid)                        :: Grid
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: Boundary                => null()
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap

    type(ocean_grid_type),         pointer :: Ocean_grid              => null()

    REAL_, pointer                         :: TW  (:,:)        => null()
    REAL_, pointer                         :: SW  (:,:)        => null()
    REAL_, pointer                         :: AREA(:,:)        => null()
    REAL_, pointer                         :: MASK(:,:)        => null()

    real, allocatable                      :: Tmp2(:,:)

    REAL_, pointer, dimension(:, :)        :: sea_lev          => null()

    integer                                :: DT_OCEAN
    character(len=7)                       :: wind_stagger     ! 'AGRID' or 'BGRID' or 'CGRID'
    integer                                ::iwind_stagger     !  AGRID  or  BGRID  or  CGRID : integer values

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

! Allocate this instance of the internal state and wrap
! -----------------------------------------------------

    allocate ( MOM_MAPL_internal_state, stat=status )
    VERIFY_(STATUS)

    wrap%ptr => MOM_MAPL_internal_state

! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    allocate ( Boundary, stat=STATUS); VERIFY_(STATUS)
    allocate ( Ocean,    stat=STATUS); VERIFY_(STATUS)

    MOM_MAPL_internal_state%Ice_ocean_boundary => Boundary
    MOM_MAPL_internal_state%Ocean              => Ocean

! FMS initialization using the communicator from the VM
!------------------------------------------------------

    call ESMF_VMGet(VM, mpiCommunicator=Comm, rc=STATUS)
    VERIFY_(STATUS)

    call fms_init(Comm)

! Init MOM stuff
!---------------

    call constants_init
    call field_manager_init
    call set_calendar_type ( JULIAN)
    call diag_manager_init                                   !SA: could pass time_init, not available before (MOM5)

    DT   = set_time (DT_OCEAN, 0)
    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)

! Check run time wind stagger option set in AGCM.rc
! to make sure it matches what is expected here
!----------------------------------------------------

    call MAPL_GetResource( MAPL, wind_stagger, Label="ocean_wind_stagger:", DEFAULT="AGRID", RC=STATUS)
    VERIFY_(STATUS)

    if ( trim(wind_stagger) == "AGRID") then
      iwind_stagger = AGRID
      if (MAPL_AM_I_Root()) print *, ' Surface stress stagger for MOM6: AGRID. Its value= ', AGRID
    elseif ( ( trim(wind_stagger) == "BGRID") .or. ( trim(wind_stagger) == "CGRID")) then
      print *, ' Surface stress stagger for MOM6: BGRID_NE or CGRID_NE. These options are not supported. Exiting!'
      ASSERT_(.false.)
    else
      print *, ' Surface stress stagger for MOM6 is invalid, stopping.'
      ASSERT_(.false.)
    endif

! Initialize ocean model
!-----------------------

    Ocean%is_ocean_pe = .true.
    call ocean_model_init  (Ocean, Ocean_state, Time, Time, iwind_stagger)
 
    MOM_MAPL_internal_state%Ocean_State => Ocean_State

    call ocean_model_init_sfc(Ocean_state, Ocean)

! Get the ocean grid and sizes of global and computational domains
!-----------------------------------------------------------------

    call get_ocean_grid (Ocean_state, Ocean_grid)
    isc  = Ocean_grid%isc; iec  = Ocean_grid%iec
    isd  = Ocean_grid%isd; ied  = Ocean_grid%ied

    jsc  = Ocean_grid%jsc; jec  = Ocean_grid%jec
    jsd  = Ocean_grid%jsd; jed  = Ocean_grid%jed

! -------
! instead of:
!   call mpp_get_compute_domain(Ocean%Domain, g_isc, g_iec, g_jsc, g_jec)
!   g_isd  = Ocean_grid%isd_global;      g_jsd  = Ocean_grid%jsd_global
!!  g_ied = ied +g_isd -1;               g_jed = jed +g_jsd -1               ! local + global -1
!   g_ied = ied+(Ocean_grid%idg_offset); g_jed = jed+(Ocean_grid%jdg_offset)
! do this:
    call mpp_get_compute_domain(Ocean_grid%Domain%mpp_domain, g_isc, g_iec, g_jsc, g_jec)
    call mpp_get_data_domain   (Ocean_grid%Domain%mpp_domain, g_isd, g_ied, g_jsd, g_jed)
! -------

! Check local sizes of horizontal dimensions
!--------------------------------------------
    call MAPL_GridGet(GRID, localCellCountPerDim=counts, RC=status)
    VERIFY_(STATUS)

    IM=iec-isc+1
    JM=jec-jsc+1

    ASSERT_(counts(1)==IM)
    ASSERT_(counts(2)==JM)

! Check run time surface current stagger option set in MOM_input 
! to make sure it matches what is expected here
!---------------------------------------------------------------

    if (MAPL_AM_I_Root()) then
     if ( (Ocean%stagger == AGRID) .or. (Ocean%stagger == BGRID_NE)) then
       print *, ' Surface velocity stagger set in ocean model: (MOM6) AGRID or BGRID_NE. These are not supported, try CGRID_NE. Exiting!'
       ASSERT_(.false.)
     elseif (Ocean%stagger == CGRID_NE) then
       print *, ' Surface velocity stagger set in ocean model: (MOM6) CGRID_NE.'
     else
       print *, ' Surface velocity stagger set in ocean model: (MOM6) is invalid, stopping.'
       ASSERT_(.false.)
     endif
    endif

! Allocate MOM flux bulletin board.
!------------------------------------

    allocate ( Boundary% u_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% v_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% t_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% q_flux          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% salt_flux       (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% lw_flux         (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_vis_dir (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_vis_dif (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_nir_dir (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% sw_flux_nir_dif (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% lprec           (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% fprec           (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% runoff          (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% calving         (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% stress_mag      (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% ustar_berg      (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% area_berg       (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% mass_berg       (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% runoff_hflx     (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% calving_hflx    (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% p               (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% mi              (g_isd:g_ied,g_jsd:g_jed), &
               Boundary% ice_rigidity    (g_isd:g_ied,g_jsd:g_jed), &
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
    Boundary%stress_mag      = 0.0
    Boundary%ustar_berg      = 0.0
    Boundary%area_berg       = 0.0
    Boundary%mass_berg       = 0.0
    Boundary%runoff_hflx     = 0.0
    Boundary%calving_hflx    = 0.0
    Boundary%p               = 0.0
    Boundary%mi              = 0.0
    Boundary% ice_rigidity   = 0.0

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

    call MAPL_GetPointer(EXPORT, MASK,     'MOM_2D_MASK', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TW,       'TW'  ,        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SW,       'SW'  ,        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA,     'AREA',        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, sea_lev,  'SLV',         alloc=.true., RC=STATUS)
    VERIFY_(STATUS)

! Get the 2-D MOM data
!---------------------
    allocate(Tmp2(IM,JM), stat=status); VERIFY_(STATUS)

    call ocean_model_data_get(Ocean_State, Ocean, 'mask', Tmp2, isc, jsc)
    MASK = real(Tmp2, kind=G5KIND)

    call ocean_model_data_get(Ocean_State, Ocean, 't_surf', Tmp2, g_isc, g_jsc) ! this comes to us in deg C
    where(MASK(:,:) > 0.0)
       TW = real(Tmp2, kind=G5KIND) + MAPL_TICE                                 ! because C to K was subtracted in MOM
    elsewhere
       TW = MAPL_UNDEF
    end where

    call ocean_model_data_get(Ocean_State, Ocean, 's_surf', Tmp2, g_isc, g_jsc) ! comes to us in PSU
    where(MASK(:,:) > 0.0)
       SW = real(Tmp2, kind=G5KIND)
    elsewhere
       SW = MAPL_UNDEF
    end where

    if(associated(area)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'area', Tmp2, isc, jsc)
       AREA = real(Tmp2, kind=G5KIND)
    end if

    deallocate(Tmp2)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine Initialize

!=================================================================================

!BOP

! !IROUTINE: Run  -- Run method for External Model Plug

! !INTERFACE:

  subroutine Run  ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC       ! Gridded component
    type(ESMF_State),    intent(INOUT) :: IMPORT   ! Import state
    type(ESMF_State),    intent(INOUT) :: EXPORT   ! Export state
    type(ESMF_Clock),    intent(INOUT) :: CLOCK    ! The supervisor clock
    integer, optional,   intent(  OUT) :: RC       ! Error code:
    type(ESMF_State)                   :: INTERNAL ! Internal state

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals with ESMF and MAPL types

    type(MAPL_MetaComp),       pointer :: MAPL               => null()
    type(ESMF_Time)                    :: MyTime
    type(ESMF_TimeInterval)            :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: Boundary                 => null()
    type(ocean_public_type),       pointer :: Ocean                    => null()
    type(ocean_state_type),        pointer :: Ocean_State              => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state  => null()
    type(MOM_MAPLWrap_Type)                :: wrap

!   type(ocean_grid_type),         pointer :: Ocean_grid               => null()

! Required exports

    REAL_, pointer                     :: TW    (:,:)        => null()
    REAL_, pointer                     :: SW    (:,:)        => null()
    REAL_, pointer                     :: UW    (:,:)        => null()
    REAL_, pointer                     :: VW    (:,:)        => null()
    REAL_, pointer                     :: UWB   (:,:)        => null()
    REAL_, pointer                     :: VWB   (:,:)        => null()
    REAL_, pointer                     :: SLV   (:,:)        => null()
    REAL_, pointer                     :: FRAZIL(:,:)        => null()
    REAL_, pointer                     :: MELT_POT(:,:)      => null()
    REAL_, pointer                     :: FRZMLT(:,:)        => null()
    REAL_, pointer                     :: MASK  (:,:)        => null()
    REAL_, pointer                     :: AREA  (:,:)        => null()

! Optional Exports
! none

! Imports
    REAL_, pointer                     :: TAUX(:,:)          => null()
    REAL_, pointer                     :: TAUY(:,:)          => null()
    REAL_, pointer                     :: PS  (:,:)          => null()
    REAL_, pointer                     :: PICE(:,:)          => null()
    REAL_, pointer                     :: LWFLX(:,:)         => null()
    REAL_, pointer                     :: SHFLX(:,:)         => null()
    REAL_, pointer                     :: QFLUX(:,:)         => null()
    REAL_, pointer                     :: RAIN(:,:)          => null()
    REAL_, pointer                     :: SNOW(:,:)          => null()
    REAL_, pointer                     :: SFLX(:,:)          => null()
    REAL_, pointer                     :: PENUVR(:,:)        => null()
    REAL_, pointer                     :: PENPAR(:,:)        => null()
    REAL_, pointer                     :: PENUVF(:,:)        => null()
    REAL_, pointer                     :: PENPAF(:,:)        => null()
    REAL_, pointer                     :: DRNIR(:,:)         => null()
    REAL_, pointer                     :: DFNIR(:,:)         => null()
    REAL_, pointer                     :: DISCHARGE(:,:)     => null()
    REAL_, pointer                     :: AICE(:,:)          => null()
    REAL_, pointer                     :: TAUXBOT(:,:)       => null()
    REAL_, pointer                     :: TAUYBOT(:,:)       => null()

! Temporaries

    real, allocatable                  :: U (:,:),  V(:,:)
    real, allocatable                  :: cos_rot(:,:)
    real, allocatable                  :: sin_rot(:,:)

    integer                            :: IM, JM

    integer                            :: steady_state_ocean = 0       ! SA: Per Atanas T, "name" of this var is misleading
                                                                       ! We run ocean model only when it = 0

    character(len=7)                   :: pres_loading                 ! yes or no

    integer                            :: isc,iec,jsc,jec

    integer                            :: YEAR,MONTH,DAY,HR,MN,SC
    type(time_type)                    :: Time
    type(time_type)                    :: DT

    real                               :: pice_scaling = 1.0
    integer                            :: DT_OCEAN


    REAL_, pointer, dimension(:,:)     :: LATS  => null()
    REAL_, pointer, dimension(:,:)     :: LONS  => null()

! Begin
!------

! Get the component name and set-up traceback handle.
! -----------------------------------------------------
    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(status)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)


    call MAPL_Get(MAPL,                      &
         INTERNAL_ESMF_STATE = INTERNAL,     &
         LATS  = LATS ,                      &
         LONS  = LONS ,                      &
                                RC=STATUS )
    VERIFY_(STATUS)

! Profilers
!----------

    call MAPL_TimerOn (MAPL,"TOTAL")
    call MAPL_TimerOn (MAPL,"RUN"  )

! Get the Plug private internal state
!--------------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR

! Aliases to MOM types
!---------------------

    Boundary    => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean       => MOM_MAPL_internal_state%Ocean
    Ocean_State => MOM_MAPL_internal_state%Ocean_State

! Get domain size
!----------------

! -------
! do this:
    call mpp_get_compute_domain(Ocean%Domain, isc, iec, jsc, jec)
! instead of:
!   call get_ocean_grid (Ocean_state, Ocean_grid)
!   isc  = Ocean_grid%isc; iec  = Ocean_grid%iec
!   jsc  = Ocean_grid%jsc; jec  = Ocean_grid%jec
! -------

    IM=iec-isc+1
    JM=jec-jsc+1

! Temporaries with MOM default reals
!-----------------------------------

    allocate(U(IM,JM   ),    stat=STATUS); VERIFY_(STATUS)
    allocate(V(IM,JM   ),    stat=STATUS); VERIFY_(STATUS)
    allocate(cos_rot(IM,JM), stat=STATUS); VERIFY_(STATUS)
    allocate(sin_rot(IM,JM), stat=STATUS); VERIFY_(STATUS)

! Get IMPORT pointers
!--------------------

    call MAPL_GetPointer(IMPORT, TAUX,     'TAUX'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUY,     'TAUY'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PS,       'PS'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PICE,     'PICE'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, LWFLX,    'LWFLX'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SHFLX,    'SHFLX'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QFLUX,    'QFLUX'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RAIN,     'RAIN'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SNOW,     'SNOW'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SFLX,     'SFLX'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENUVR,   'PENUVR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENPAR,   'PENPAR'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENUVF,   'PENUVF'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PENPAF,   'PENPAF'  ,  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DRNIR,    'DRNIR'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DFNIR,    'DFNIR'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DISCHARGE,'DISCHARGE', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, AICE,     'AICE',      RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUXBOT,  'TAUXBOT',   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUYBOT,  'TAUYBOT',   RC=STATUS); VERIFY_(STATUS)

! Get EXPORT pointers
!--------------------

    call MAPL_GetPointer(EXPORT, UW,    'UW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VW,    'VW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UWB,   'UWB' ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VWB,   'VWB' ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TW,    'TW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SW,    'SW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SLV,   'SLV',    RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, FRAZIL,  'FRAZIL',   alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MELT_POT,'MELT_POT', alloc=.true., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZMLT,  'FRZMLT',                 RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, MASK, 'MOM_2D_MASK', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA, 'AREA',        RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, pice_scaling, Label = "MOM_PICE_SCALING:", default = 1.0, rc = status); VERIFY_(status)

! Fill in ocean boundary fluxes/forces
!-------------------------------------

    call MAPL_GetResource( MAPL, pres_loading, Label="pres_loading:", DEFAULT="NO", RC=STATUS)
    VERIFY_(STATUS)

    ! NOTE: PICE that is available here is all = 0. This should be made realistic, for now it is from MOM5 legacy
    !       Need to study with zero pressure loading (CTL: as now), exp1 ( with PS only), exp2 (with PS and PICE), exp3 (PICE only).
    if ( pres_loading == "YES") then
      Boundary%P              (isc:iec,jsc:jec)= real(PS,          kind=KIND(Boundary%p)) + & ! Pressure of overlying atmospheric
                                  pice_scaling * real(PICE,        kind=KIND(Boundary%p))     ! and ice
    else
      Boundary%P              (isc:iec,jsc:jec)= &
                                   pice_scaling* real(PICE,        kind=KIND(Boundary%p)) ! Pressure of overlying ice only
    end if

    Boundary%lw_flux        (isc:iec,jsc:jec)= real(LWFLX,         kind=KIND(Boundary%p)) ! Long wave flux: both positive down
    Boundary%t_flux         (isc:iec,jsc:jec)= real(SHFLX,         kind=KIND(Boundary%p)) ! Sensible heat flux: both positive up
    Boundary%q_flux         (isc:iec,jsc:jec)= real(QFLUX,         kind=KIND(Boundary%p)) ! specific humidity flux [kg m-2 s-1] ( OR evaporation flux ?)
    Boundary%lprec          (isc:iec,jsc:jec)= real(RAIN,          kind=KIND(Boundary%p)) ! Liquid precipitation: both positive down
    Boundary%fprec          (isc:iec,jsc:jec)= real(SNOW,          kind=KIND(Boundary%p)) ! Frozen precipitation: both positive down
    Boundary%salt_flux      (isc:iec,jsc:jec)=-real(SFLX,          kind=KIND(Boundary%p)) ! Salt flux: MOM positive up, GEOS positive down
    Boundary%runoff         (isc:iec,jsc:jec)= real(DISCHARGE,     kind=KIND(Boundary%p)) ! mass flux of liquid runoff [kg m-2 s-1]

! All shortwave components are positive down in MOM and in GEOS
!--------------------------------------------------------------
    Boundary%sw_flux_vis_dir(isc:iec,jsc:jec)= real(PENUVR+PENPAR, kind=KIND(Boundary%p)) ! direct visible sw radiation        [W m-2]
    Boundary%sw_flux_vis_dif(isc:iec,jsc:jec)= real(PENUVF+PENPAF, kind=KIND(Boundary%p)) ! diffuse visible sw radiation       [W m-2]
    Boundary%sw_flux_nir_dir(isc:iec,jsc:jec)= real(DRNIR,         kind=KIND(Boundary%p)) ! direct  Near InfraRed sw radiation [W m-2]
    Boundary%sw_flux_nir_dif(isc:iec,jsc:jec)= real(DFNIR,         kind=KIND(Boundary%p)) ! diffuse Near InfraRed sw radiation [W m-2]

! Convert input stresses over water to MOM wind stagger
!------------------------------------------------------
    U = 0.0; V = 0.0
!   Using A-grid ice stress, note ice (CICE) stress has opposite sign to atmosphere
    U = real( TAUX*(1.-AICE) - TAUXBOT*AICE, kind=kind(U))
    V = real( TAUY*(1.-AICE) - TAUYBOT*AICE, kind=kind(V))

! Grid rotation angles - these could be saved in the first instance, rather doing every time
!   cos_rot = 1. ! A-grid
    call ocean_model_data_get(Ocean_State, Ocean, 'cos_rot', cos_rot, isc, jsc)
!   sin_rot = 0. ! A-grid
    call ocean_model_data_get(Ocean_State, Ocean, 'sin_rot', sin_rot, isc, jsc)

! Rotate input stress over water along i,j of tripolar grid, and combine with stress under ice
!---------------------------------------------------------------------------------------------
    Boundary%U_flux = 0.0;  Boundary%V_flux = 0.0 ! Initialize stress

    Boundary%U_flux  (isc:iec,jsc:jec)= real( (U*cos_rot - V*sin_rot), kind=KIND(Boundary%p))
    Boundary%V_flux  (isc:iec,jsc:jec)= real( (U*sin_rot + V*cos_rot), kind=KIND(Boundary%p))

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

! Run MOM for one time step
!--------------------------

    ! set following to non-zero only if the coupled model becomes unstable 
    ! (inconsistent atmosphere or bad restart! or some instabilities) - per Atanas T
    call MAPL_GetResource(MAPL, steady_state_ocean, Label = "steady_state_ocean:", default = 0, rc = status); VERIFY_(status)

    if(steady_state_ocean == 0) then
      call update_ocean_model(Boundary, Ocean_State, Ocean, Time, DT)
    else
      ! SA: steady_state_ocean: Call update_ocean_model with additional args now available with MOM6
      !     though this will cause MOM_FATAL error
      call update_ocean_model(Boundary, Ocean_State, Ocean, Time, DT, .false., .false.)
    endif

! Get required Exports at GEOS precision
!---------------------------------------

!   mask
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 'mask', U, isc, jsc)
    MASK = real(U, kind=G5KIND)

!   surface (potential) temperature (K)
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 't_surf', U, isc, jsc) ! this comes to us in deg C
    where(MASK(:,:) > 0.0)
       TW = real(U, kind=G5KIND) + MAPL_TICE                             ! because C to K was subtracted in MOM
    elsewhere
       TW = MAPL_UNDEF
    end where

!   surface salinity (PSU)
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 's_surf', U, isc, jsc) ! this comes to us in PSU
    where(MASK(:,:) > 0.0)
       SW = real(U, kind=G5KIND)
    elsewhere
       SW = MAPL_UNDEF
    end where

!   sea level (m)
    U = 0.0
    if(associated(SLV)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', U, isc, jsc) ! this comes to us in m
       where(MASK(:,:)>0.0)
          SLV = real(U, kind = G5KIND)
       elsewhere
          SLV = 0.0
       end where
    end if

!   frazil (W/m^2)
    U = 0.0
    if(associated(FRAZIL)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'frazil',   U, isc, jsc)  ! this comes to us in J/m2
       where(MASK(:,:)>0.0)
          FRAZIL =  real( (U)/dt_ocean, kind = G5KIND) ! relying on fortran to promote the int (dt_ocean) to real
       elsewhere
          FRAZIL =  0.0
       end where
    endif

!   melt potential (W/m^2)
    U = 0.0
    if(associated(MELT_POT)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'melt_pot', U, isc, jsc)  ! this comes to us in J/m2
       where(MASK(:,:)>0.0)
          MELT_POT = -real( (U)/dt_ocean, kind = G5KIND) ! relying on fortran to promote the int (dt_ocean) to real
       elsewhere
          MELT_POT =  0.0
       end where
       MELT_POT = MIN ( MELT_POT, 0.0) ! make sure melt potential is <= 0
    endif

!   freezing melt potential (W/m^2)
    if(associated(FRZMLT)) then
       if ( (.not.associated(FRAZIL)) .or. (.not.associated(MELT_POT))) then
         print *, 'You are asking for freeze melt potential, without asking for frazil and melt potential. You must ask for all. Exiting!'
         ASSERT_(.false.)
       endif

       where(MASK(:,:)>0.0)
          FRZMLT = FRAZIL + MELT_POT
       elsewhere
          FRZMLT = 0.0
       end where
    end if

! currents (m/s)
!---------------
!   A-grid currents (for the atmospheric model)
    U = 0.0; V = 0.0
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'ua', U, isc, jsc) ! this comes to us in m/s
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'va', V, isc, jsc) ! this comes to us in m/s

    if(associated(UW )) then
      where(MASK(:,:) > 0.0)
        UW = real(U, kind=G5KIND)
      elsewhere
        UW=0.0
      end where
    endif

    if(associated(VW )) then
      where(MASK(:,:) > 0.0)
        VW = real(V, kind=G5KIND)
      elsewhere
        VW=0.0
      end where
    end if

!   B-grid currents (for CICE dynamics)
    U = 0.0; V = 0.0
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'ub', U, isc, jsc) ! this comes to us in m/s
    call ocean_model_get_UV_surf(Ocean_State, Ocean, 'vb', V, isc, jsc) ! this comes to us in m/s

    if(associated(UWB )) then
      where(MASK(:,:) > 0.0)
        UWB = real(U, kind=G5KIND)
      elsewhere
        UWB =0.0
      end where
    endif

    if(associated(VWB )) then
      where(MASK(:,:) > 0.0)
        VWB = real(V, kind=G5KIND)
      elsewhere
        VWB =0.0
      end where
    end if

! Optional Exports at GEOS precision
!-----------------------------------
! none
!   3d exports with MOM6, such as depths, T, S, U, V, etc
!   will not be exported. If needed, write them on tri-polar grid directly from MOM6

    deallocate(U, V)
    deallocate(cos_rot,sin_rot)

    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

! All Done
!---------
    RETURN_(ESMF_SUCCESS)
  end subroutine Run

!BOP

!====================================================================

! !IROUTINE: Finalize        -- Finalize method for GuestOcean wrapper

! !INTERFACE:

  subroutine Finalize ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
  type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
  integer, optional,   intent(  OUT) :: RC     ! Error code

!EOP

    type(MAPL_MetaComp),           pointer :: MAPL
    type(ESMF_Time)                        :: MyTime
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(ice_ocean_boundary_type), pointer :: Boundary                => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
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

! Get the Plug private internal state
!--------------------------------------

    CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    MOM_MAPL_internal_state => WRAP%PTR

    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean    => MOM_MAPL_internal_state%Ocean
    Ocean_State => MOM_MAPL_internal_state%Ocean_State

! Set the times for MOM
!----------------------

    call ESMF_ClockGet( CLOCK, currTime=MyTime, RC=STATUS)
    VERIFY_(status)

    call ESMF_TimeGet (MyTime,      &
         YY=YEAR, MM=MONTH, DD=DAY, &
         H=HR,    M =MN,    S =SC,  &
         RC=STATUS )
    VERIFY_(STATUS)

    Time = set_date(YEAR,MONTH,DAY,HR,MN,SC)

    call ocean_model_end (Ocean, Ocean_State, Time) ! SA: this also calls ocean_model_save_restart(...)

    call diag_manager_end(Time )
    call field_manager_end
    call fms_io_exit

    deallocate ( Boundary% u_flux          , &
                 Boundary% v_flux          , &
                 Boundary% t_flux          , &
                 Boundary% q_flux          , &
                 Boundary% salt_flux       , &
                 Boundary% lw_flux         , &
                 Boundary% sw_flux_vis_dir , &
                 Boundary% sw_flux_vis_dif , &
                 Boundary% sw_flux_nir_dir , &
                 Boundary% sw_flux_nir_dif , &
                 Boundary% lprec           , &
                 Boundary% fprec           , &
                 Boundary% runoff          , &
                 Boundary% calving         , &
                 Boundary% stress_mag      , &
                 Boundary% ustar_berg      , &
                 Boundary% area_berg       , &
                 Boundary% mass_berg       , &
                 Boundary% runoff_hflx     , &
                 Boundary% calving_hflx    , &
                 Boundary% p               , &
                 Boundary% mi              , &
                 Boundary% ice_rigidity    , &
                                stat=STATUS )
    VERIFY_(STATUS)

    deallocate ( Ocean,                   STAT=STATUS); VERIFY_(STATUS)
    deallocate ( Boundary,                STAT=STATUS); VERIFY_(STATUS)
    deallocate ( MOM_MAPL_internal_state, STAT=STATUS); VERIFY_(STATUS)
!

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

! !IROUTINE: Record -- Record method for GuestOcean wrapper (write intermediate restarts)

! !INTERFACE:

  subroutine Record ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
  type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
  integer, optional,   intent(  OUT) :: RC     ! Error code

!EOP

    type(MAPL_MetaComp),     pointer :: MAPL
    type(MOM_MAPL_Type),     pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)          :: wrap
    type(ocean_state_type),  pointer :: Ocean_State             => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
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

    doRecord = MAPL_RecordAlarmIsRinging(MAPL, RC=status)
    VERIFY_(STATUS)

    if (doRecord) then

! Get the Plug private internal state
!--------------------------------------

       CALL ESMF_UserCompGetInternalState( GC, 'MOM_MAPL_state', WRAP, STATUS )
       VERIFY_(STATUS)

       MOM_MAPL_internal_state => WRAP%PTR
       Ocean_State             => MOM_MAPL_internal_state%Ocean_State

       call MAPL_DateStampGet(clock, timeStamp, rc=status)
       VERIFY_(STATUS)

! Write a restart
!-----------------

       call ocean_model_restart (Ocean_State, timeStamp)
       VERIFY_(STATUS)

    end if

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

!====================================================================

end module MOM6_GEOSPlugMod

subroutine SetServices(gc, rc)
   use ESMF
   use MOM6_GEOSPlugMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine

