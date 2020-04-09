!  $Id$

#include "MAPL_Generic.h"

! GEOS   default real kind

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

module MOM6_GEOSPlugMod

!BOP
! !MODULE: MOM6_GEOSPlugMod -- to couple with MOM6.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a couple for MOM.
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
  use MAPL_Mod
  use MAPL_ConstantsMod,        only: MAPL_TICE

! These MOM dependencies are all we are currently using.

  use constants_mod,            only: constants_init
  use diag_manager_mod,         only: diag_manager_init, diag_manager_end
  use field_manager_mod,        only: field_manager_init, field_manager_end

  use fms_mod,                  only: fms_init, fms_end
  use fms_io_mod,               only: fms_io_exit

  use mpp_domains_mod,          only: domain2d, mpp_update_domains, mpp_get_compute_domain
  use mpp_parameter_mod,        only: AGRID, SCALAR_PAIR

  use time_manager_mod,         only: set_calendar_type, time_type
  use time_manager_mod,         only: set_time, set_date
  use time_manager_mod,         only: JULIAN

  use ocean_model_mod,          only: ocean_model_init,   ocean_model_init_sfc, &
                                      update_ocean_model, ocean_model_end, ocean_model_restart
  use ocean_model_mod,          only: ocean_model_data_get
  use ocean_model_mod,          only: ocean_public_type, ocean_state_type
  use MOM_surface_forcing,      only: ice_ocean_boundary_type

! Nothing on the MOM side is visible through this module.

  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices
!EOP

! These are the MOM-side bulletin boards, where things are in
! MOM precision and the B grid

  type MOM_MAPL_Type
     type(ocean_public_type)         :: Ocean               ! mom6/config_src/coupled_driver/ocean_model_MOM.F90     [MOM6 GEOS]
                                                            ! mom/src/mom5/ocean_core/ocean_types.F90                [MOM5 GEOS]

     type(ice_ocean_boundary_type)   :: Ice_ocean_boundary  ! mom6/config_src/coupled_driver/MOM_surface_forcing.F90 [MOM6 GEOS]
                                                            ! mom/src/mom5/ocean_core/ocean_types.F90                [MOM5 GEOS]
  end type MOM_MAPL_Type

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


    call MAPL_AddImportSpec(GC,                                   &
        SHORT_NAME         = 'STROCNXB',                          &
        LONG_NAME          = 'x_stress_at_base_of_ice_weighted_by_aiu',    &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                  &
        SHORT_NAME         = 'STROCNYB',                          &
        LONG_NAME          = 'y_stress_at_base_of_ice_weighted_by_aiu',   &
        UNITS              = 'N m-2',                             &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )

     VERIFY_(STATUS)

     call MAPL_AddImportSpec(GC,                                   &
          SHORT_NAME         = 'AICEU',                            &
          LONG_NAME          = 'ice_concentration_of_grid_cell_Bgrid',   &
          UNITS              = '1',                                &
          DIMS               = MAPL_DimsHorzOnly,                  &
          VLOCATION          = MAPL_VLocationNone,                 &
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
!SA: following horz dims is wrong.
!        DIMS               = MAPL_DimsHorzVert,                   &
!
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
         UNITS              = 'J m-2',                            &
         DIMS               = MAPL_DimsHorzOnly,                  &
         VLOCATION          = MAPL_VLocationNone,                 &
         RC=STATUS  )
    VERIFY_(STATUS)

!  !Diagnostic exports

    call MAPL_AddExportSpec(GC,                                    &
         SHORT_NAME         = 'DH',                                &
         LONG_NAME          = 'layer_thickness',                   &
         UNITS              = 'm OR kg m-2',                       &
         DIMS               = MAPL_DimsHorzVert,                   &
         VLOCATION          = MAPL_VLocationCenter,                &
         RC=STATUS  )
    VERIFY_(STATUS)

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

    call MAPL_AddExportSpec(GC,                               &   ! SA: already have DH. Get rid of it. SOON!
         SHORT_NAME         = 'DEPTH',                        &
         LONG_NAME          = 'layer_depth',                  &
         UNITS              = 'm',                            &
         DIMS               = MAPL_DimsHorzVert,              &
         VLOCATION          = MAPL_VLocationCenter,           &
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

    character(len=ESMF_MAXSTR)	       :: IAm
    integer			       :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals

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
    type(MAPL_MetaComp), pointer           :: MAPL
    type(ESMF_Grid)                        :: Grid
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

! Locals

    type(ice_ocean_boundary_type), pointer :: Boundary                => null()
    type(ocean_public_type),       pointer :: Ocean                   => null()
    type(ocean_state_type),        pointer :: Ocean_State             => null()
    type(domain2d),                pointer :: OceanDomain             => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)                :: wrap

    integer                                :: DT_OCEAN

    REAL_, pointer                         :: TW  (:,:)        => null()
    REAL_, pointer                         :: SW  (:,:)        => null()
    REAL_, pointer                         :: DH  (:,:,:)      => null()
    REAL_, pointer                         :: AREA(:,:)        => null()
    REAL_, pointer                         :: MASK(:,:)        => null()

    real, allocatable                      :: Tmp2(:,:)

    REAL_, pointer, dimension(:, :)        :: sea_lev => null()
    REAL_, pointer, dimension(:, :, :)     :: TL      => null()
    REAL_, pointer, dimension(:, :, :)     :: SL      => null()

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

    call ESMF_UserCompSetInternalState ( GC, 'MOM_MAPL_state', WRAP, STATUS )
    VERIFY_(STATUS)

    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean    => MOM_MAPL_internal_state%Ocean
!SA: what to be done for:
!   Ocean_State ??

! FMS initialization using the communicator from the VM
!------------------------------------------------------

    call ESMF_VMGet(VM, mpiCommunicator=Comm, rc=STATUS)
    VERIFY_(STATUS)

    call fms_init(Comm)                                      ! FMS/fms/fms.F90                               [MOM6 GEOS]
                                                             ! mom/src/shared/fms/fms.F90                    [MOM5 GEOS]
! Init MOM stuff
!---------------

    call constants_init                                      ! FMS/constants/constants.F90                   [MOM6 GEOS] SA: sync with MAPL_Constants
                                                             ! mom/src/shared/constants/constants.F90        [MOM5 GEOS]

    call field_manager_init                                  ! FMS/field_manager/field_manager.F90           [MOM6 GEOS]
                                                             ! mom/src/shared/field_manager/field_manager.F90[MOM5 GEOS]

    call set_calendar_type ( JULIAN)                         ! FMS/time_manager/time_manager.F90             [MOM6 GEOS]
                                                             ! mom/src/shared/time_manager/time_manager.F90  [MOM5 GEOS]

    call diag_manager_init                                   ! FMS/diag_manager/diag_manager.F90             [MOM6 GEOS] SA: could pass time_init, not available before (MOM5)
                                                             ! mom/src/shared/diag_manager/diag_manager.F90  [MOM5 GEOS]

    DT   = set_time (DT_OCEAN, 0)                            ! FMS/time_manager/time_manager.F90             [MOM6 GEOS]
                                                             ! mom/src/shared/time_manager/time_manager.F90  [MOM5 GEOS]

    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)                ! FMS/time_manager/time_manager.F90             [MOM6 GEOS]
                                                             ! mom/src/shared/time_manager/time_manager.F90  [MOM5 GEOS]

    Ocean%is_ocean_pe = .true.
    call ocean_model_init  (Ocean, Ocean_state, Time, Time)  ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
                                                             ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]

! you probably need this as well? Check with Alistair, Bob.
    call ocean_model_init_sfc(Ocean_state, Ocean)

! Check local sizes of horizontal dimensions
!--------------------------------------------

!   call mom4_get_dimensions(isc, iec, jsc, jec, nk_out=LM)          ! mom/src/mom5/ocean_core/ocean_model.F90    [MOM5 GEOS]

    OceanDomain => Ocean%Domain
    call mpp_get_compute_domain(OceanDomain, isc, iec, jsc, jec)            ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]

!   print *, '[isc, iec], [jsc, jec] = ', '[', isc, iec, ']', '[', jsc, jec, ']'

!   following is incorrect
!   isc = Ocean_state%grid%isc; iec = Ocean_state%grid%iec           ! SA: Seems to be no mom6_get_dimensions in MOM6.
!   jsc = Ocean_state%grid%jsc; jec = Ocean_state%grid%jec

!   print *, 'Ocean_state%grid%isc = ', Ocean_state%grid%isc
!   print *, 'Ocean_state%grid%jsc = ', Ocean_state%grid%jsc
!   print *, 'Ocean_state%grid%iec = ', Ocean_state%grid%iec
!   print *, 'Ocean_state%grid%jec = ', Ocean_state%grid%jec

    call MAPL_GridGet(GRID, localCellCountPerDim=counts, RC=status)
    VERIFY_(STATUS)

    IM=iec-isc+1
    JM=jec-jsc+1
    LM=Ocean_state%GV%ke

    ASSERT_(counts(1)==IM)
    ASSERT_(counts(2)==JM)

! Allocate MOM flux bulletin board.
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
               Boundary% stress_mag      (isc:iec,jsc:jec), &        ! SA: additions in MOM6
               Boundary% ustar_berg      (isc:iec,jsc:jec), &
               Boundary% area_berg       (isc:iec,jsc:jec), &
               Boundary% mass_berg       (isc:iec,jsc:jec), &
               Boundary% runoff_hflx     (isc:iec,jsc:jec), &
               Boundary% calving_hflx    (isc:iec,jsc:jec), &
               Boundary% p               (isc:iec,jsc:jec), &
               Boundary% mi              (isc:iec,jsc:jec), &
               Boundary% ice_rigidity    (isc:iec,jsc:jec), &
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
    Boundary%stress_mag      = 0.0   ! SA: additions in MOM6
    Boundary%ustar_berg      = 0.0
    Boundary%area_berg       = 0.0
    Boundary%mass_berg       = 0.0
    Boundary%runoff_hflx     = 0.0
    Boundary%calving_hflx    = 0.0
    Boundary%p               = 0.0
    Boundary% mi             = 0.0
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

    call MAPL_GetPointer(EXPORT, TW,       'TW'  ,        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SW,       'SW'  ,        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DH,       'DH'  ,        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MASK,     'MOM_2D_MASK', alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA,     'AREA',        alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, sea_lev,  'SLV',         alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TL,       'T',           alloc=.true., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SL,       'S',           alloc=.true., RC=STATUS)
    VERIFY_(STATUS)

! Get the 3-D MOM data
!---------------------
                                                           ! SA: cut and loose with MASK! Anyway, it is 2D now.
    DH = real(Ocean_state%MOM_CSp%h, kind=G5KIND)          ! mom6/src/core/MOM.F90                   [MOM6 GEOS] layer thickness [H ~> m or kg m-2]
!   call mom4_get_thickness(Tmp3)                          ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS] thickness (in meters) of each layer

    SL = real(Ocean_state%MOM_CSp%S, kind=G5KIND)          ! mom6/src/core/MOM.F90                   [MOM6 GEOS] in ppt
!   call mom4_get_salinity_index(i)                        ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS]
!   call mom4_get_prog_tracer(i,fld=Tmp3)

    TL = real(Ocean_state%MOM_CSp%T + MAPL_TICE, kind=G5KIND) ! mom6/src/core/MOM.F90                   [MOM6 GEOS] potential temperature [degC]
!   call mom4_get_prog_tracer(i,fld=Tmp3)                     ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS]

! Get the 2-D MOM data
!---------------------
    allocate(Tmp2(IM,JM), stat=status); VERIFY_(STATUS)

    call ocean_model_data_get(Ocean_State, Ocean, 'mask', Tmp2, isc, jsc)       ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
!   call mom4_get_3D_tmask(Tmp3)                                                ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS] ! SA: No need for 3D mask, 2D is good enough!
    MASK = real(Tmp2, kind=G5KIND)

    call ocean_model_data_get(Ocean_State, Ocean, 't_surf', Tmp2, isc, jsc)     ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS] ! SA: in K
!   call mom4_get_Tsurf(Ocean,Tmp2)                                             ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
    where(MASK(:,:) > 0.0)
       TW = real(Tmp2, kind=G5KIND)
    elsewhere
       TW = MAPL_UNDEF
    end where

    call ocean_model_data_get(Ocean_State, Ocean, 's_surf', Tmp2, isc, jsc)     ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS] ! SA: in PSU
!   call mom4_get_Ssurf(Ocean,Tmp2)                                             ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
    where(MASK(:,:) > 0.0)
       SW = real(Tmp2, kind=G5KIND)
    elsewhere
       SW = MAPL_UNDEF
    end where

    if(associated(area)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'area', Tmp2, isc, jsc)     ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
                                                                                 ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
       AREA = real(Tmp2, kind=G5KIND)
    end if

    if(associated(sea_lev)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', Tmp2, isc, jsc)  ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS] ! SA: includes Inv Baro in M
                                                                                 ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
       sea_lev = real(merge(tsource = Tmp2, fsource = real(MAPL_UNDEF), mask = (MASK(:, :) > 0.0)), kind=G5KIND)
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

    character(len=ESMF_MAXSTR)		   :: IAm
    integer				                  :: STATUS
    character(len=ESMF_MAXSTR)         :: COMP_NAME

! Locals

    integer                            :: IM, JM, LM
    integer                            :: IMw, JMw
    integer                            :: I
    integer                            :: J

    integer                            :: steady_state_ocean = 0     ! SA: Per Atanas T, "name" of this var is misleading! We run ocean model only when it = 0 !
    logical                            :: ocean_seg_start    = .true.
    logical                            :: ocean_seg_end      = .true.

! Required exports

    REAL_, pointer                     :: TW    (:,:)        => null()
    REAL_, pointer                     :: SW    (:,:)        => null()
    REAL_, pointer                     :: UW    (:,:)        => null()
    REAL_, pointer                     :: VW    (:,:)        => null()
    REAL_, pointer                     :: UWB   (:,:)        => null()
    REAL_, pointer                     :: VWB   (:,:)        => null()
    REAL_, pointer                     :: DH    (:,:,:)      => null()
    REAL_, pointer                     :: SLV   (:,:)        => null()
    REAL_, pointer                     :: FRAZIL(:,:)        => null()
    REAL_, pointer                     :: MASK  (:,:)        => null()
    REAL_, pointer                     :: AREA  (:,:)        => null()
    REAL_, pointer                     :: DEPTH(:,:,:)       => null()

! Optional Exports

    REAL_, pointer                     :: TL  (:,:,:)       => null()
    REAL_, pointer                     :: SL  (:,:,:)       => null()
    REAL_, pointer                     :: SWFRAC(:,:,:)     => null()

! Imports
    REAL_, pointer                     :: TAUX(:,:)         => null()
    REAL_, pointer                     :: TAUY(:,:)         => null()
    REAL_, pointer                     :: PS  (:,:)         => null()
    REAL_, pointer                     :: PICE(:,:)         => null()
    REAL_, pointer                     :: HEAT(:,:,:)       => null()
    REAL_, pointer                     :: LWFLX(:,:)        => null()
    REAL_, pointer                     :: SHFLX(:,:)        => null()
    REAL_, pointer                     :: QFLUX(:,:)        => null()
    REAL_, pointer                     :: RAIN(:,:)         => null()
    REAL_, pointer                     :: SNOW(:,:)         => null()
    REAL_, pointer                     :: SFLX(:,:)         => null()
    REAL_, pointer                     :: PENUVR(:,:)       => null()
    REAL_, pointer                     :: PENPAR(:,:)       => null()
    REAL_, pointer                     :: PENUVF(:,:)       => null()
    REAL_, pointer                     :: PENPAF(:,:)       => null()
    REAL_, pointer                     :: DRNIR(:,:)        => null()
    REAL_, pointer                     :: DFNIR(:,:)        => null()
    REAL_, pointer                     :: DISCHARGE(:,:)    => null()
    REAL_, pointer                     :: STROCNXB(:,:)     => null()
    REAL_, pointer                     :: STROCNYB(:,:)     => null()
    REAL_, pointer                     :: AICEU(:,:)        => null()

! Temporaries

    real, allocatable                  :: U(:,:)
    real, allocatable                  :: V(:,:)
!   real, allocatable                  :: H(:,:,:) ! SA: not used
!   real, allocatable                  :: G(:,:,:) ! SA: not used
    real, allocatable                  :: cos_rot(:,:)
    real, allocatable                  :: sin_rot(:,:)
    real                               :: EPSLN

    type(MAPL_MetaComp),           pointer :: MAPL                     => null()
    type(MOM_MAPL_Type),           pointer :: MOM_MAPL_internal_state  => null()
    type(MOM_MAPLWrap_Type)                :: wrap

    type(ice_ocean_boundary_type), pointer :: Boundary     => null()
    type(ocean_public_type),       pointer :: Ocean        => null()
    type(ocean_state_type),        pointer :: Ocean_State  => null()
    type(domain2d),                pointer :: OceanDomain  => null()

    integer                                :: isc,iec,jsc,jec
    integer                                :: isd,ied,jsd,jed

    integer                                :: YEAR,MONTH,DAY,HR,MN,SC
    type(time_type)                        :: Time
    type(time_type)                        :: DT

    integer                                :: l
    real                                   :: pice_scaling = 1.0
    integer           :: DT_OCEAN
    real, parameter   :: CW = 3992.10322329649                          ! SA: should use MAPL_Constants or add to it!

    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

    REAL_, pointer, dimension(:,:)         :: LATS  => null()
    REAL_, pointer, dimension(:,:)         :: LONS  => null()

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

    Boundary => MOM_MAPL_internal_state%Ice_ocean_boundary
    Ocean    => MOM_MAPL_internal_state%Ocean

! Get domain size
!----------------

    OceanDomain => Ocean%Domain
    call mpp_get_compute_domain(OceanDomain, isc, iec, jsc, jec)            ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
!   call get_ocean_domain(OceanDomain)                                      ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
!   call mom4_get_dimensions(isc, iec, jsc, jec, isd, ied, jsd, jed, LM)    ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
    print *, '[isc, iec], [jsc, jec]:', '[', isc, iec, ']', '[', jsc, jec, ']'

! these [ (isd, ied), (jsd, jed)] need to filled but NOT following way.
! since ocean_state isn't yet specified.  
!   isd = Ocean_state%grid%isd; ied = Ocean_state%grid%ied
!   jsd = Ocean_state%grid%jsd; jed = Ocean_state%grid%jed

!   print *, 'Ocean_state%grid%isd = ', Ocean_state%grid%isd
!   print *, 'Ocean_state%grid%jsd = ', Ocean_state%grid%jsd
!   print *, 'Ocean_state%grid%ied = ', Ocean_state%grid%ied
!   print *, 'Ocean_state%grid%jed = ', Ocean_state%grid%jed

    IM=iec-isc+1
    JM=jec-jsc+1
!   LM=Ocean_state%GV%ke
    LM=50
    print *, 'IM, JM, LM=', IM, JM, LM

! Temporaries with MOM default reals
!-----------------------------------

    allocate(U(IM,JM   ),    stat=STATUS); VERIFY_(STATUS)
    allocate(V(IM,JM   ),    stat=STATUS); VERIFY_(STATUS)
!   allocate(H(IM,JM,LM),    stat=STATUS); VERIFY_(STATUS)  ! SA: this is not needed ??
!   allocate(G(IM,JM,LM),    stat=STATUS); VERIFY_(STATUS)  ! SA: this is not needed ??
    allocate(cos_rot(IM,JM), stat=STATUS); VERIFY_(STATUS)
    allocate(sin_rot(IM,JM), stat=STATUS); VERIFY_(STATUS)

! Get IMPORT pointers
!--------------------

    call MAPL_GetPointer(IMPORT, TAUX,     'TAUX'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUY,     'TAUY'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PS,       'PS'    ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PICE,     'PICE'  ,    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, HEAT,     'SWHEAT',    RC=STATUS); VERIFY_(STATUS)
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
    call MAPL_GetPointer(IMPORT, STROCNXB, 'STROCNXB',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, STROCNYB, 'STROCNYB',  RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, AICEU,    'AICEU',     RC=STATUS); VERIFY_(STATUS)

! Get EXPORT pointers
!--------------------

    call MAPL_GetPointer(EXPORT, UW,    'UW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VW,    'VW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UWB,   'UWB' ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VWB,   'VWB' ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TW,    'TW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SW,    'SW'  ,   RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SLV,   'SLV',    RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRAZIL,'FRAZIL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DEPTH, 'DEPTH',  RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, DH,   'DH'  ,        RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MASK, 'MOM_2D_MASK', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AREA, 'AREA',        RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, TL,   'T'   ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SL,   'S'   ,       RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SWFRAC, 'SWFRAC',   RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, pice_scaling, Label = "MOM_PICE_SCALING:", default = 1.0, rc = status); VERIFY_(status)

    Boundary%P         = pice_scaling* &
                         real(PICE,      kind=KIND(Boundary%p)) ! Pressure of overlying ice and atmosphere
    Boundary%lw_flux   = real(LWFLX,     kind=KIND(Boundary%p)) ! Long wave flux: both positive down
    Boundary%t_flux    = real(SHFLX,     kind=KIND(Boundary%p)) ! Sensible heat flux: both positive up
    Boundary%q_flux    = real(QFLUX,     kind=KIND(Boundary%p)) ! specific humidity flux [kg m-2 s-1] ( OR evaporation flux ?)
    Boundary%lprec     = real(RAIN,      kind=KIND(Boundary%p)) ! Liquid precipitation: both positive down
    Boundary%fprec     = real(SNOW,      kind=KIND(Boundary%p)) ! Frozen precipitation: both positive down
    Boundary%salt_flux =-real(SFLX,      kind=KIND(Boundary%p)) ! Salt flux: MOM positive up, GEOS positive down
    Boundary%runoff    = real(DISCHARGE, kind=KIND(Boundary%p)) ! mass flux of liquid runoff [kg m-2 s-1]

! All shortwave components are positive down  in MOM and in GEOS
!---------------------------------------------------------------
    Boundary%sw_flux_vis_dir = real(PENUVR+PENPAR, kind=KIND(Boundary%p)) ! direct visible sw radiation        [W m-2]
    Boundary%sw_flux_vis_dif = real(PENUVF+PENPAF, kind=KIND(Boundary%p)) ! diffuse visible sw radiation       [W m-2]
    Boundary%sw_flux_nir_dir = real(DRNIR,         kind=KIND(Boundary%p)) ! direct Near InfraRed sw radiation  [W m-2]
    Boundary%sw_flux_nir_dif = real(DFNIR,         kind=KIND(Boundary%p)) ! diffuse Near InfraRed sw radiation [W m-2]

! Convert input stresses over water to B grid
!--------------------------------------------
    U = 0.0
    V = 0.0
! comment for now
!   call transformA2B( real(TAUX,kind=kind(U)), real(TAUY,kind=kind(V)), U, V)

! Rotate input stress over water along i,j of tripolar grid, and combine with stress under ice
!---------------------------------------------------------------------------------------------
    call ocean_model_data_get(Ocean_State, Ocean, 'cos_rot', cos_rot, isc, jsc)    ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
    call ocean_model_data_get(Ocean_State, Ocean, 'sin_rot', sin_rot, isc, jsc)    ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]

! comment for now
!   Boundary%U_flux =  (U*cos_rot + V*sin_rot)*(1.-AICEU) - STROCNXB
!   Boundary%V_flux = (-U*sin_rot + V*cos_rot)*(1.-AICEU) - STROCNYB

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

    DT   = set_time (DT_OCEAN, 0)                       ! FMS/time_manager/time_manager.F90             [MOM6 GEOS]
                                                        ! mom/src/shared/time_manager/time_manager.F90  [MOM5 GEOS]

    Time = set_date (YEAR,MONTH,DAY,HR,MN,SC)           ! FMS/time_manager/time_manager.F90             [MOM6 GEOS]
                                                        ! mom/src/shared/time_manager/time_manager.F90  [MOM5 GEOS]

! Run MOM for one time step
!--------------------------

    ! steady_state_ocean is always (by default) = 0 when we have an ocean model in GEOS.
    ! set it to non-zero only if the coupled model becomes unstable (inconsistent atmosphere - bad restart! or some instabilities) - per Atanas T
    call MAPL_GetResource(MAPL, steady_state_ocean, Label = "steady_state_ocean:", default = 0, rc = status); VERIFY_(status)

    ! SA: steady_state_ocean is not needed. Call update_ocean_model with additional args now available with MOM6
    if(steady_state_ocean == 0) then

      print *, 'plug: Run: call update_ocean_model'
      call update_ocean_model(Boundary, Ocean_State, Ocean, Time, DT)  ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
      print *, 'plug: Run: called update_ocean_model'
      ASSERT_(.FALSE.)
                                                                       ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
    endif

! Get export fields

! Required Exports at GEOS precision
!-----------------------------------

!   surface (potential) temperature (K)
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 't_surf', U, isc, jsc) ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS] ! SA: in K
!   call mom4_get_Tsurf(Ocean,U)                                         ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
    where(MASK(:,:) > 0.0)
       TW = real(U, kind=G5KIND)
    elsewhere
       TW = MAPL_UNDEF
    end where

!   surface salinity (PSU)
    U = 0.0
    call ocean_model_data_get(Ocean_State, Ocean, 's_surf', U, isc, jsc) ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS] ! SA: in PSU
!   call mom4_get_Ssurf(Ocean,U)                                         ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
    where(MASK(:,:) > 0.0)
       SW = real(U, kind=G5KIND)
    elsewhere
       SW = MAPL_UNDEF
    end where

!   3D T, S and depth, they are not 'masked'
!   ocean salinity (3D; PSU)
    SL = real( Ocean_state%MOM_CSp%S, kind=G5KIND) ! mom6/src/core/MOM.F90                   [MOM6 GEOS]
!   call mom4_get_salinity_index(i)                ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS]
!   call mom4_get_prog_tracer(i,fld=H)

!   ocean potential temperature (3D; K)
    TL = real(Ocean_state%MOM_CSp%T + MAPL_TICE, kind=G5KIND) ! mom6/src/core/MOM.F90                   [MOM6 GEOS] potential temperature [degC]
!   call mom4_get_prog_tracer(i,fld=H)                        ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS]

    DH = real(Ocean_state%MOM_CSp%h, kind=G5KIND)  ! mom6/src/core/MOM.F90                   [MOM6 GEOS] layer thickness [H ~> m or kg m-2]
!   call mom4_get_thickness(H)                     ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS] thickness (in meters) of each layer

! Optional Exports at GEOS precision
!-----------------------------------

! Get the A grid currents at MOM precision
!-----------------------------------------
! SA: like MOM5, MOM6 also runs with Bgrid staggering (default) for currents. So currents are on Bgrid. Need B to A grid.
!     something like... mct_driver/ocn_cap_methods.F90 "rotate ssh gradients from local coordinates..." which is what happens in mom4_get_latlon_UVsurf
    if(associated(UW) .or. associated(VW)) then
       U = 0.0; V = 0.0 ! SA: for now!
!      call mom4_get_latlon_UVsurf(OCEAN, U, V, STATUS)    ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
!      VERIFY_(STATUS)
    endif

    if(associated(UW  )) then
       where(MASK(:,:) > 0.0)
          UW = real(U, kind=G5KIND)
       elsewhere
          UW=0.0
       end where
     endif

    if(associated(VW  )) then
       where(MASK(:,:) > 0.0)
          VW = real(V, kind=G5KIND)
       elsewhere
          VW=0.0
       end where
    end if

! Get the B grid currents at MOM precision
!-----------------------------------------

    if(associated(UWB) .or. associated(VWB)) then
       UWB = Ocean%u_surf  ! MOM6 run with Bgrid staggering (its default) ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
       VWB = Ocean%v_surf  ! MOM6 run with Bgrid staggering (its default)
!      call mom4_get_UVsurfB(OCEAN, U, V, STATUS)                         ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
!      VERIFY_(STATUS)
    endif

    if(associated(UWB  )) then
       where(MASK(:,:)>0.0)
          UWB = real(U, kind=G5KIND)
       elsewhere
          UWB=0.0
       end where
     endif

    if(associated(VWB  )) then
       where(MASK(:,:)>0.0)
          VWB = real(V, kind=G5KIND)
       elsewhere
          VWB=0.0
       end where
    end if

    if(associated(SLV)) then
       call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', U, isc, jsc) ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS] ! SA: in PSU
!      call ocean_model_data_get(Ocean_State, Ocean, 'sea_lev', U, isc, jsc) ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
       where(MASK(:,:)>0.0)
          SLV = real(U, kind = G5KIND)
       elsewhere
          SLV=0.0
       end where
    end if

    if(associated(FRAZIL)) then
       FRAZIL = Ocean%frazil                                                ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
!      call ocean_model_data_get(Ocean_State, Ocean, 'frazil', U, isc, jsc) ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
       where(MASK(:,:)>0.0)
          FRAZIL = real(U, kind = G5KIND)
       elsewhere
          FRAZIL=0.0
       end where
    end if

    if(associated(DEPTH)) then
       DEPTH = real(Ocean_state%MOM_CSp%h, kind=G5KIND)                          ! mom6/src/core/MOM.F90                   [MOM6 GEOS] layer thickness [H ~> m or kg m-2]
!      call ocean_model_data_get(Ocean_State, Ocean, 'geodepth_zt', H, isc, jsc) ! mom/src/mom5/ocean_core/ocean_model.F90 [MOM5 GEOS]
    end if

!   deallocate(H)
!   deallocate(G)
    deallocate(U,V)
    deallocate(cos_rot,sin_rot)

    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

! All Done
!---------
    RETURN_(ESMF_SUCCESS)
  contains

    subroutine transformA2B(U, V, uvx, uvy) ! SA: Per Atanas T, GEOS "likes" winds/currents on A-grid, so we need stuff like this...

      real, intent(IN)     :: U (:,:)
      real, intent(IN)     :: V (:,:)

      real, INTENT(INOUT)  :: uvx(isc:,jsc:)
      real, INTENT(INOUT)  :: uvy(isc:,jsc:)

      integer              :: i, j, ii, jj, cnt
      real, allocatable    :: tx(:,:), ty(:,:)
      real                 :: sum


      allocate(tx(isd:ied,jsd:jed), stat=STATUS); VERIFY_(STATUS)
      allocate(ty(isd:ied,jsd:jed), stat=STATUS); VERIFY_(STATUS)

      tx(isc:iec, jsc:jec) = U
      ty(isc:iec, jsc:jec) = V

      call mpp_update_domains(tx, ty, OceanDomain, gridtype=AGRID, flags=SCALAR_PAIR)   ! FMS/mpp/mpp_domains.F90            [MOM6 GEOS]
                                                                                        ! mom/src/shared/mpp/mpp_domains.F90 [MOM5 GEOS]

      do j = jsc, jec
         do i = isc, iec
            sum = 0.0
            cnt = 0
            do ii = 0,1
               do jj = 0,1
                  if (tx(i+ii,j+jj) /= MAPL_UNDEF) then
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
                  if (ty(i+ii,j+jj) /= MAPL_UNDEF) then
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

    type(MAPL_MetaComp),       pointer :: MAPL
    type(ESMF_Time)                    :: MyTime
    type(MOM_MAPL_Type),       pointer :: MOM_MAPL_internal_state => null()
    type(MOM_MAPLWrap_Type)            :: wrap
    type(ocean_public_type),   pointer :: Ocean                   => null()
    type(ocean_state_type),    pointer :: Ocean_State             => null()
    type(ice_ocean_boundary_type), pointer :: Boundary            => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)	    :: IAm
    integer			                   :: STATUS
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

    Ocean => MOM_MAPL_internal_state%Ocean

! Set the times for MOM
!----------------------

    call ESMF_ClockGet( CLOCK, currTime=MyTime, RC=STATUS)
    VERIFY_(status)

    call ESMF_TimeGet (MyTime,      &
         YY=YEAR, MM=MONTH, DD=DAY, &
         H=HR,    M =MN,    S =SC,  &
         RC=STATUS )
    VERIFY_(STATUS)

    Time = set_date(YEAR,MONTH,DAY,HR,MN,SC)        ! FMS/time_manager/time_manager.F90                  [MOM6 GEOS]
                                                    ! mom/src/shared/time_manager/time_manager.F90       [MOM5 GEOS]

    call ocean_model_end (Ocean, Ocean_State, Time) ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
                                                    ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]
                                                    ! SA: unlike in mom5, this also calls ocean_model_save_restart(...)
                                                    !     possible duplication of ocean_model_restart(...) called by below Record

    call diag_manager_end(Time )                    ! FMS/diag_manager/diag_manager.F90                  [MOM6 GEOS]
                                                    ! mom/src/shared/diag_manager/diag_manager.F90       [MOM5 GEOS]

    call field_manager_end                          ! FMS/field_manager/field_manager.F90                [MOM6 GEOS]
                                                    ! mom/src/shared/field_manager/field_manager.F90     [MOM5 GEOS]

    call fms_io_exit                                ! FMS/fms/fms_io.F90                                 [MOM6 GEOS]
                                                    ! mom/src/shared/fms/fms_io.F90                      [MOM5 GEOS]

! deallocate

    deallocate ( Boundary% u_flux        , &
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
               Boundary% stress_mag      , &        ! SA: additions in MOM6
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

    deallocate( MOM_MAPL_internal_state, STAT=STATUS); VERIFY_(STATUS)
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

! !IROUTINE: Record        -- Record method for GuestOcean wrapper (write intermediate restarts)

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

    character(len=ESMF_MAXSTR)	    :: IAm
    integer			                   :: STATUS
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

       call MAPL_DateStampGet(clock, timeStamp, rc=status)
       VERIFY_(STATUS)

! Write a restart
!-----------------

       call ocean_model_restart (Ocean_State, timeStamp) ! mom6/config_src/coupled_driver/ocean_model_MOM.F90 [MOM6 GEOS]
       VERIFY_(STATUS)                                   ! mom/src/mom5/ocean_core/ocean_model.F90            [MOM5 GEOS]

    end if

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

!====================================================================

end module MOM6_GEOSPlugMod
