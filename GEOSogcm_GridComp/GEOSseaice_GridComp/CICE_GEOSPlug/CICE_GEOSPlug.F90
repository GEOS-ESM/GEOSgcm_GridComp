!  $Id$

#include "MAPL_Generic.h"

! GEOS   default real kind

#define G5KIND      4
#define REAL_       real(kind=G5KIND)

module CICE_GEOSPlugMod

!BOP
! !MODULE: CICE_GEOSPlugMod -- to couple with CICE6 and later.

!DESCRIPTION:
! A  MAPL/ESMF Gridded Component that acts as a coupler for CICE6.
!

!USES:
  use ESMF
  use MAPL
  use CICE_InitMod                 
  use CICE_FinalMod                 
  use CICE_RunMod                 
  use ice_import_export


  implicit none
  private

  !PUBLIC MEMBER FUNCTIONS:
  public :: SetServices

  integer            :: NUM_ICE_CATEGORIES
  logical, private   :: ice_grid_init2 

contains

  
!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the CICE_GEOSplug GC needs to register its
!   Initialize, Run and Finalize.  It uses the MAPL_Generic construct for defining
!   state specs and couplings among its children.  

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
    character(len=ESMF_MAXSTR)          :: COMP_NAME

! Locals
    type (ESMF_Config)                  :: CF

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME,  CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ice_grid_init2 = .false.

    call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
    VERIFY_(STATUS)

! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,   Initialize, _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,          Run,        _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,     Finalize,   _RC)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_WRITERESTART, Record,     _RC)

! Set the state variable specs.
! -----------------------------

!BOS
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
         SHORT_NAME         = 'UWC',                                &
         LONG_NAME          = 'water_skin_eastward_cgrid_velocity', &
         UNITS              = 'm s-1 ',                            &
         DIMS               = MAPL_DimsHorzOnly,                   &
         VLOCATION          = MAPL_VLocationNone,                  &
         DEFAULT            = 0.0,                                 &
         RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddImportSpec(GC,                               &
         SHORT_NAME         = 'VWC',                                &
         LONG_NAME          = 'water_skin_northward_cgrid_velocity',&
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


  call MAPL_AddImportSpec(GC                         ,&
          SHORT_NAME         = 'FRZMLT'                    ,&
          LONG_NAME          = 'freeze_melt_potential',     &
          UNITS              = 'W m-2'                     ,&
          DIMS               = MAPL_DimsHorzOnly,           &
          VLOCATION          = MAPL_VLocationNone          ,&
          DEFAULT            = 0.0,                         &
          _RC  )

  call MAPL_AddImportSpec(GC                         ,&
          SHORT_NAME         = 'SST'                       ,&
          LONG_NAME          = 'sea_surface_temperature'   ,&
          UNITS              = 'K'                         ,&
          DIMS               = MAPL_DimsHorzOnly,           &
          VLOCATION          = MAPL_VLocationNone          ,&
          !DEFAULT            = MAPL_TICE - 1.8,             &
          _RC  )

  call MAPL_AddImportSpec(GC                         ,&
          SHORT_NAME         = 'SSS'                       ,&
          LONG_NAME          = 'sea_surface_salinity'      ,&
          UNITS              = 'psu'                       ,&
          DIMS               = MAPL_DimsHorzOnly,           &
          VLOCATION          = MAPL_VLocationNone          ,&
          DEFAULT            = 33.0,                        &
          _RC  )


  !call MAPL_AddImportSpec(GC,                               &
  !       SHORT_NAME         = 'FROCEAN',                           &
  !       LONG_NAME          = 'fraction_of_gridbox_covered_by_skin',&
  !       UNITS              = '1',                                 &
  !       DIMS               = MAPL_DimsHorzOnly,                   &
  !       VLOCATION          = MAPL_VLocationNone,                  &
  !       RESTART            = MAPL_RestartSkip,                    &
  !       RC=status  )
  !VERIFY_(status)

  call MAPL_AddImportSpec(GC,                                 &
        SHORT_NAME         = 'SI',                                &
        LONG_NAME          = 'seaice_skin_salinity',              &
        UNITS              = 'psu',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
        RESTART            = MAPL_RestartSkip,                    &
        DEFAULT            = 4.,                                  &
        _RC  )


  ! === Exports

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'TI',                                &
    LONG_NAME          = 'seaice_skin_temperature',           &
    UNITS              = 'K',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    _RC )

  call MAPL_AddExportSpec(GC,                            &
    SHORT_NAME         = 'FRSEAICE',                           &
    LONG_NAME          = 'fractional_cover_of_seaice',        &
    UNITS              = '1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    UNGRIDDED_DIMS     = (/NUM_ICE_CATEGORIES/),              &
    VLOCATION          = MAPL_VLocationNone,                  &
    _RC )

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

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'AICE',                            &
    LONG_NAME          = 'ice_concentration_of_grid_cell',   &
    UNITS              = '1',                                 &
    DIMS               = MAPL_DimsHorzOnly,                   &
    VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'FRACICE',                       &
    LONG_NAME          = 'fractional_cover_of_seaice',    &
    UNITS              = '1',                             &
    DIMS               = MAPL_DimsHorzOnly,               &
    VLOCATION          = MAPL_VLocationNone,              &
                                                     _RC  )

  call MAPL_AddExportSpec(GC,                             &
    SHORT_NAME         = 'HICE',                          &
    LONG_NAME          = 'mean_ice_thickness_of_grid_cell', &
    UNITS              = 'm',                             &
    DIMS               = MAPL_DimsHorzOnly,               &
    VLOCATION          = MAPL_VLocationNone,              &
                                                      _RC )
  VERIFY_(STATUS)

  call MAPL_AddExportSpec(GC,                                &
    SHORT_NAME         = 'HSNO',                            &
    LONG_NAME          = 'mean_snow_thickness_of_grid_cell',   &
    UNITS              = 'm',                                 &
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


  !*CALLBACK*
  !=================================================================================
  ! an ESMF state to pass information b.w. GCs using callback
  ! to be connected to the import in SeaiceInterface
  !
  call MAPL_AddExportSpec(GC                                ,&
          SHORT_NAME         = 'SURFSTATE'                  ,&
          LONG_NAME          = 'surface_state_for_seaice_thermo_coupling',  &
          UNITS              = 'W m-2'                      ,&
          DATATYPE           = MAPL_StateItem               ,&
                                                       __RC__)

 !=================================================================================

!EOS


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

! !IROUTINE: INITIALIZE -- Initialize method for CICE wrapper

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

    character(len=ESMF_MAXSTR)             :: IAm
    integer                                :: STATUS
    character(len=ESMF_MAXSTR)             :: COMP_NAME

! Locals

    integer                                :: YEAR,MONTH,DAY,HR,MN,SC

! Locals 
    integer                                :: Comm 

    integer                                :: NPES
    integer                                :: OGCM_IM, OGCM_JM
    integer                                :: OGCM_NX, OGCM_NY
    integer                                :: BLK_NX,  BLK_NY
    integer                                :: counts(7)

! Locals with ESMF and MAPL types

    type(ESMF_VM)                          :: VM
    type(MAPL_MetaComp), pointer           :: MAPL
    type(ESMF_Grid)                        :: Grid
    type(ESMF_Time)                        :: MyTime
    type(ESMF_TimeInterval)                :: TINT

! Locals

    type(ESMF_State)                       :: SURFST


    REAL_, pointer                         :: FROCEAN(:,:)        => null()

    integer                                :: DT_SEAICE

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

    call MAPL_TimerOn(MAPL,"INITIALIZE")

! Generic initialize
! ------------------

    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL"     )

! Get the grid, configuration
!----------------------------
    !call MAPL_GetPointer(IMPORT, FROCEAN,     'FROCEAN'  ,    __RC__)

    call ESMF_GridCompGet( GC, grid=Grid,  RC=status )
    VERIFY_(STATUS)

! Get the layout from the grid
!-----------------------------

    call ESMF_VMGetCurrent(VM, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL,DT_SEAICE,  Label="RUN_DT:",    _RC)             ! Get AGCM Heartbeat
    call MAPL_GetResource(MAPL,DT_SEAICE,  Label="OCEAN_DT:",  DEFAULT=DT_SEAICE, _RC) ! set Default OCEAN_DT to AGCM Heartbeat

! Set the time for CICE
!---------------------

    call ESMF_ClockGet(CLOCK, currTIME=MyTime, TimeStep=TINT,  RC=STATUS)
    VERIFY_(STATUS)

    call ESMF_TimeGet (MyTime,                    &
                       YY=YEAR, MM=MONTH, DD=DAY, &
                       H=HR,    M =MN,    S =SC,  &
                                               _RC)


    ! Get the ocean layout from the VM
    !---------------------------------

    call MAPL_GetResource( MAPL, OGCM_IM, Label="OGCM.IM_WORLD:", __RC__)
    call MAPL_GetResource( MAPL, OGCM_JM, Label="OGCM.JM_WORLD:", __RC__)
    call MAPL_GetResource( MAPL, OGCM_NX, Label="OGCM.NX:",       __RC__)
    call MAPL_GetResource( MAPL, OGCM_NY, Label="OGCM.NY:",       __RC__)

    !ASSERT_(mod(OGCM_IM,OGCM_NX)==0)
    !ASSERT_(mod(OGCM_JM,OGCM_NY)==0)
    if (mod(OGCM_IM,OGCM_NX) == 0) then
        BLK_NX = OGCM_IM/OGCM_NX
    else
        BLK_NX = OGCM_IM/OGCM_NX + 1
    endif 
    if (mod(OGCM_JM,OGCM_NY) == 0) then
        BLK_NY = OGCM_JM/OGCM_NY
    else
        BLK_NY = OGCM_JM/OGCM_NY + 1
    endif 

    call ESMF_VMGet(VM, mpiCommunicator=Comm,  petCount=NPES, rc=STATUS)
    VERIFY_(STATUS)

    call MAPL_GridGet(GRID, localCellCountPerDim=counts, RC=status)
    VERIFY_(STATUS)


! Init CICE 
!---------------
    call cice_init1(Comm, NPES, BLK_NX, BLK_NY, &
                    DT_SEAICE, MAPL_TICE, MAPL_ALHL, MAPL_ALHS)

    !call ice_import_grid(FROCEAN, __RC__)
    !if (counts(2) /= BLK_NY) &
    !     print*,counts(2),BLK_NY
    ! there should be an ASSERT here to make sure the block size from CICE and MAPL agree
    ! block sizes from MAPL are contained in counts
    ! there should be a call to CICE which returns the block size for current PE
    ! the two should match
    !ASSERT_(counts(1) == BLK_NX)  
    !ASSERT_(counts(2) == BLK_NY)  
     
    !call cice_init2(YEAR, MONTH, DAY, HR, MN, SC) ! init cice calendar here
    call cice_cal_init(YEAR, MONTH, DAY, HR, MN, SC) ! init cice calendar here

    !*CALLBACK*
    !=====================================================================================
    call ESMF_StateGet(EXPORT, 'SURFSTATE', SURFST, _RC)

    !!attach the thermo coupling method
    !
    call ESMF_MethodAdd(SURFST, label='thermo_coupling', userRoutine=thermo_coupling, _RC)
    call ESMF_MethodAdd(SURFST, label='prep_albedo',     userRoutine=prep_albedo,     _RC)

    call LoadSurfaceStates(SURFST, _RC)

    ! the value of len below must be the maximum string length in VARLIST
    call LoadOcnVars(SURFST, IMPORT, VARLIST=[character(len=6)::'SST','SSS','FRZMLT'], _RC)
    !=====================================================================================



    call MAPL_TimerOff(MAPL,"TOTAL"     )


! Profilers
! ---------
    call MAPL_TimerOff(MAPL,"INITIALIZE")

! Make sure exports neede by the parent prior to our run call are initialized
!----------------------------------------------------------------------------


! All Done
!---------

    RETURN_(ESMF_SUCCESS)

contains

  subroutine LoadSurfaceStates(SURFST, RC)

! !ARGUMENTS:
    type(ESMF_State), intent(INOUT)    :: SURFST
    integer, optional, intent(OUT)     :: RC


! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS

    type(ESMF_Grid)                    :: grid

! Get the component name and set-up traceback handle.
! -----------------------------------------------------
    Iam = trim(comp_name) // "LoadSurfaceStates"

    ! create and add some fields to the callback state (i.e. SURFSTATE)
    call ESMF_GridCompGet(gc, grid=grid, __RC__)
    !call AddSurfField('TSKINICE', SURFST, GRID,     &    
    !                 UGRID=NUM_ICE_CATEGORIES, __RC__)
    call AddSurfField('EVAP', SURFST, GRID,         &
                     UGRID=NUM_ICE_CATEGORIES, __RC__)
    call AddSurfField('LHF', SURFST, GRID,          &
                     UGRID=NUM_ICE_CATEGORIES, __RC__)
    call AddSurfField('FSURF', SURFST, GRID,        &
                     UGRID=NUM_ICE_CATEGORIES, __RC__)
    call AddSurfField('DFSURFDTS', SURFST, GRID,    &
                     UGRID=NUM_ICE_CATEGORIES, __RC__)
    call AddSurfField('DLHFDTS', SURFST, GRID,      &
                     UGRID=NUM_ICE_CATEGORIES, __RC__)
    call AddSurfField('SNOW',  SURFST, GRID,  __RC__)
    call AddSurfField('RAIN',  SURFST, GRID,  __RC__)
    call AddSurfField('DRPAR', SURFST, GRID,  __RC__)
    call AddSurfField('DFPAR', SURFST, GRID,  __RC__)
    call AddSurfField('DRNIR', SURFST, GRID,  __RC__)
    call AddSurfField('DFNIR', SURFST, GRID,  __RC__)
    call AddSurfField('DRUVR', SURFST, GRID,  __RC__)
    call AddSurfField('DFUVR', SURFST, GRID,  __RC__)
    call AddSurfField('COSZ',  SURFST, GRID,  __RC__)

    call AddSurfField('FROCEAN',  SURFST, GRID,  __RC__)

    ! callback return fields

    ! fields with categories
    call AddSurfField('DTS',   SURFST, GRID,        & 
                     UGRID=NUM_ICE_CATEGORIES, __RC__)
    ! aggregated fields
    call AddSurfField('ALBVR',  SURFST, GRID,  __RC__)
    call AddSurfField('ALBVF',  SURFST, GRID,  __RC__)
    call AddSurfField('ALBNR',  SURFST, GRID,  __RC__)
    call AddSurfField('ALBNF',  SURFST, GRID,  __RC__)
    call AddSurfField('PENUVR', SURFST, GRID,  __RC__)
    call AddSurfField('PENUVF', SURFST, GRID,  __RC__)
    call AddSurfField('PENPAR', SURFST, GRID,  __RC__)
    call AddSurfField('PENPAF', SURFST, GRID,  __RC__)
    call AddSurfField('GHTSKIN', SURFST, GRID,  __RC__)

    RETURN_(ESMF_SUCCESS)
    
  end subroutine LoadSurfaceStates 

  subroutine AddSurfField(FLD_NAME, SURFST, GRID, UGRID, RC)

! !ARGUMENTS:
    character(len=*),    intent(IN)    :: FLD_NAME
    type(ESMF_State), intent(INOUT)    :: SURFST
    type(ESMF_Grid),  intent(INOUT)    :: GRID
    integer, optional,   intent(IN)    :: UGRID
    integer, optional,  intent(OUT)    :: RC

! ErrLog Variables

    character(len=ESMF_MAXSTR)         :: IAm
    integer                            :: STATUS

    type(ESMF_Field)                   :: fld

    Iam = trim(comp_name) // "AddSurfField"

    fld = MAPL_FieldCreateEmpty(FLD_NAME, GRID, __RC__)
    if (present(UGRID)) then 
       call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzOnly,           &
           location=MAPL_VLocationNone, typekind=MAPL_R4, hw=0,          &
           ungrid=[UGRID], __RC__)
    else
       call MAPL_FieldAllocCommit(fld, dims=MAPL_DimsHorzOnly,           &
           location=MAPL_VLocationNone, typekind=MAPL_R4, hw=0,          &
           __RC__)
    endif
    call MAPL_StateAdd(SURFST, fld, __RC__)

    RETURN_(ESMF_SUCCESS)
  
  end subroutine AddSurfField  

  subroutine LoadOcnVars(SURFST, IMPORT, VARLIST, RC)
     type(ESMF_State), intent(in) :: import
     type(ESMF_State), intent(inout) :: surfst
     character(len=*), intent(in) :: varlist(:)
     integer, optional, intent(out) :: RC
 
     type(ESMF_Field) :: field
     integer :: n, status
 
     do n=1, size(varlist)
         call ESMF_StateGet(import, varlist(n), field, _RC)
         call ESMF_StateAdd(surfst, [field], _RC)
     enddo
      _RETURN(ESMF_SUCCESS)
  end subroutine LoadOcnVars

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


!   type(ocean_grid_type),         pointer :: Ocean_grid               => null()

! Required exports

    REAL_, pointer                     :: MASK  (:,:)        => null()
    REAL_, pointer                     :: AREA  (:,:)        => null()
    REAL_, pointer                     :: FRI   (:,:)        => null()
    REAL_, pointer                     :: AICE  (:,:)        => null()
    REAL_, pointer                     :: FRESH (:,:)        => null()
    REAL_, pointer                     :: FHOCN (:,:)        => null()
    REAL_, pointer                     :: FSALT (:,:)        => null()
    REAL_, pointer                     :: TI    (:,:,:)      => null()
    REAL_, pointer                     :: FI    (:,:,:)      => null()
    REAL_, pointer                     :: TAUXBOT(:,:)       => null()
    REAL_, pointer                     :: TAUYBOT(:,:)       => null()
    REAL_, pointer                     :: UI    (:,:)        => null()
    REAL_, pointer                     :: VI    (:,:)        => null()

! Optional Exports
! none

! Imports
    REAL_, pointer                     :: TAUX(:,:)          => null()
    REAL_, pointer                     :: TAUY(:,:)          => null()
    REAL_, pointer                     :: UWB(:,:)           => null()
    REAL_, pointer                     :: VWB(:,:)           => null()
    REAL_, pointer                     :: UWC(:,:)           => null()
    REAL_, pointer                     :: VWC(:,:)           => null()
    REAL_, pointer                     :: SLV(:,:)           => null()

! Temporaries

    integer                            :: IM, JM

    integer                            :: YEAR,MONTH,DAY,HR,MN,SC



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


! 
!---------------------

    call MAPL_GetPointer(IMPORT, TAUX,     'TAUX'        ,                 _RC)
    call MAPL_GetPointer(IMPORT, TAUY,     'TAUY'        ,                 _RC)
    call MAPL_GetPointer(IMPORT,  SLV,     'SLV'         ,                 _RC)
    call MAPL_GetPointer(IMPORT,  UWB,     'UWB'         ,                 _RC)
    call MAPL_GetPointer(IMPORT,  VWB,     'VWB'         ,                 _RC)
    call MAPL_GetPointer(IMPORT,  UWC,     'UWC'         ,                 _RC)
    call MAPL_GetPointer(IMPORT,  VWC,     'VWC'         ,                 _RC)

    call MAPL_GetPointer(EXPORT,   TI,     'TI'          ,  alloc=.true.,  _RC)
    call MAPL_GetPointer(EXPORT,   FI,     'FRSEAICE'    ,  alloc=.true.,  _RC)
    call MAPL_GetPointer(EXPORT,   UI,     'UI'          ,  alloc=.true.,  _RC)
    call MAPL_GetPointer(EXPORT,   VI,     'VI'          ,  alloc=.true.,  _RC)
    call MAPL_GetPointer(EXPORT,   TAUXBOT,'TAUXBOT'     ,  alloc=.true.,  _RC)
    call MAPL_GetPointer(EXPORT,   TAUYBOT,'TAUYBOT'     ,  alloc=.true.,  _RC)

    call MAPL_GetPointer(EXPORT,   FRI,    'FRACICE'     ,                 _RC)
    call MAPL_GetPointer(EXPORT,   AICE,   'AICE'        ,                 _RC)
    call MAPL_GetPointer(EXPORT,   FHOCN,  'FHOCN'       ,                 _RC)
    call MAPL_GetPointer(EXPORT,   FRESH,  'FRESH'       ,                 _RC)
    call MAPL_GetPointer(EXPORT,   FSALT,  'FSALT'       ,                 _RC)


    !call ice_import_thermo2()

    call ice_import_dyna(TAUX, TAUY, SLV, UWB, VWB, UWC, VWC, _RC)


    call CICE_Run

! Get exports needed by GEOS
!---------------------

    !call ice_export_thermo2

    call ice_export_dyna(TAUXBOT, TAUYBOT, UI, VI, _RC)


    if(associated(TI)) then
      call ice_export_field('TI',       TI, _RC)
    endif 

    if(associated(FI)) then
      call ice_export_field('FRSEAICE', FI, _RC)
    endif 

    if(associated(FRI)) then
      call ice_export_field('FRACICE', FRI, _RC)
    endif 

    if(associated(AICE)) then
      call ice_export_field('FRACICE', AICE, _RC)
    endif 

    if(associated(FHOCN)) then
      call ice_export_field('FHOCN', FHOCN, _RC)
    endif 

    if(associated(FRESH)) then
      call ice_export_field('FRESH', FRESH, _RC)
    endif 

    if(associated(FSALT)) then
      call ice_export_field('FSALT', FSALT, _RC)
    endif 


    call MAPL_TimerOff(MAPL,"RUN"   )
    call MAPL_TimerOff(MAPL,"TOTAL" )

! All Done
!---------
    RETURN_(ESMF_SUCCESS)
  end subroutine Run

!BOP

!====================================================================

! !IROUTINE: Record -- Record method (write intermediate restarts)

! !INTERFACE:

  subroutine Record ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: GC     ! Gridded component
  type(ESMF_State),    intent(INOUT) :: IMPORT ! Import state
  type(ESMF_State),    intent(INOUT) :: EXPORT ! Export state
  type(ESMF_Clock),    intent(INOUT) :: CLOCK  ! The supervisor clock
  integer, optional,   intent(  OUT) :: RC     ! Error code

!EOP

    type(MAPL_MetaComp),     pointer :: MAPL

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals
    character(len=14)                :: timeStamp
    logical                          :: doRecord

    __Iam__('Record')

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, _RC)
    Iam = trim(COMP_NAME)//'::'//Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, _RC)

! Profilers
!----------

    call MAPL_TimerOn(MAPL,"TOTAL")

    doRecord = MAPL_RecordAlarmIsRinging(MAPL, MODE=MAPL_Write2Disk, _RC)

    if (doRecord) then

! Get the private internal state
!--------------------------------


       call MAPL_DateStampGet(clock, timeStamp, _RC)

! Write a restart
!-----------------

       call ice_checkpoint(timeStamp)

    end if

    call MAPL_TimerOff(MAPL,"TOTAL")
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

!====================================================================

! !IROUTINE: Finalize        -- Finalize method for CICE wrapper

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

! ErrLog Variables

    character(len=ESMF_MAXSTR)       :: IAm
    integer                          :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals 

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


! Set the times for CICE
!----------------------

    call ESMF_ClockGet( CLOCK, currTime=MyTime, RC=STATUS)
    VERIFY_(status)

    call ESMF_TimeGet (MyTime,      &
         YY=YEAR, MM=MONTH, DD=DAY, &
         H=HR,    M =MN,    S =SC,  &
         RC=STATUS )
    VERIFY_(STATUS)

    call CICE_Finalize !BZ: note save restarts is in ice_step??


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


  !*CALLBACK*
  !=====================================================================================
  subroutine thermo_coupling(state, rc)

  !! Arguments
  !! ---------
     type(ESMF_State)                      :: state
     integer,          intent(OUT)         :: rc

!EOP

     integer                               :: status
     REAL_, pointer                        :: FRO(:,:)          => null()


! ErrLog Variables

     character(len=ESMF_MAXSTR), parameter   :: IAm=' thermo_coupling'

     if (.not. ice_grid_init2) then
        call MAPL_GetPointer(state, FRO, 'FROCEAN', __RC__)
        call ice_import_grid(FRO, rc=STATUS)
        VERIFY_(STATUS)
        call cice_init2
        ice_grid_init2 = .TRUE.
     endif

     ! unpack fields and send them to cice
     call ice_import_thermo1(state, rc=STATUS)
     
     ! let cice update surface temperature and fluxes 
     call ice_fast_physics     

     ! export the relevant fields from cice
     call ice_export_thermo1(state, rc=STATUS)

     ! pack them back into state 

     RETURN_(ESMF_SUCCESS)

  end subroutine thermo_coupling

  subroutine prep_albedo(state, rc)

  !! Arguments
  !! ---------
     type(ESMF_State)                        :: state
     integer,           intent(OUT)          :: rc

!EOP

     integer                               :: status


! ErrLog Variables

     character(len=ESMF_MAXSTR), parameter   :: IAm=' prep_albedo'

     ! unpack fields and send them to cice
     call ice_import_radiation(state, rc=STATUS)
     
     ! let cice update surface temperature and fluxes 
     call ice_radiation     

     ! export the relevant fields from cice
     call ice_export_radiation(state, rc=STATUS)

     ! pack them back into state 

     RETURN_(ESMF_SUCCESS)

  end subroutine prep_albedo 

  !=====================================================================================


!====================================================================


end module CICE_GEOSPlugMod

subroutine SetServices(gc, rc)
   use ESMF
   use CICE_GEOSPlugMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine
