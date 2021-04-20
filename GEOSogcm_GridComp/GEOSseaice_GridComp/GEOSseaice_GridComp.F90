

#include "MAPL_Generic.h"

module GEOSseaice_GridCompMod

!BOP
! !MODULE:  GEOSseaice_GridCompMod -- Implements ESMF wrapper to invoke the DATASEAICE/CICE4/CICE6 seaice models.

! !USES:

  use ESMF
  use MAPL
  use GEOS_CICEDynaGridCompMod,          only : CICE4SeaIceSetServices  => SetServices
  use GEOS_DataSeaIceGridCompMod,        only : DataSeaIceSetServices   => SetServices
  use ice_prescribed_mod,                only : ice_nudging

   

  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices

  character(len=ESMF_MAXSTR)          :: SEAICE_NAME
  integer                             :: DO_DATASEAICE

! !DESCRIPTION:
!
!   {\tt GEOSseaice\_GridComp} is a light-weight gridded component that serves an
!   interface to seaice/data\_seaice components.
!
!EOP

  type :: T_PrivateState
     type(ESMF_Clock)  :: CLOCK
  end type T_PrivateState

  type :: T_PrivateState_Wrap
     type(T_PrivateState), pointer :: ptr
  end type T_PrivateState_Wrap




  integer ::          ICE
  integer ::          ICEd
  logical ::      DUAL_OCEAN


contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for GEOSseaice

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
    !character(len=ESMF_MAXSTR)         :: sharedObj

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )
    Iam = trim(COMP_NAME) // Iam


! Set the state variable specs.
! -----------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__ )

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL,  DO_DATASEAICE,  Label="USE_DATASEAICE:" ,  DEFAULT=1, __RC__ )

! Initialize these IDs (0 means not used)
! ---------------------------------------
    ICE = 0
    ICEd = 0

    if(DO_DATASEAICE/=0) then
       SEAICE_NAME="DATASEAICE"
       ICE = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=DataSeaiceSetServices, __RC__)
    else
       call MAPL_GetResource ( MAPL, SEAICE_NAME, Label="SEAICE_NAME:", DEFAULT="CICE4", __RC__ )
       select case (trim(SEAICE_NAME))
          case ("CICE4")
             ICE = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=CICE4SeaIceSetServices, __RC__)
          case default
             charbuf_ = "SEAICE_NAME: " // trim(SEAICE_NAME) // " is not implemented, ABORT!"
             call WRITE_PARALLEL(charbuf_)
             VERIFY_(999)
       end select
    endif

    call MAPL_GetResource( MAPL, iDUAL_OCEAN, 'DUAL_OCEAN:', default=0, __RC__ )
    DUAL_OCEAN = iDUAL_OCEAN /= 0

    if (dual_ocean) then
       SEAICE_NAME="DATASEAICE"
       ICEd = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=DataSeaiceSetServices, __RC__) 
    endif

! Set the state variable specs.
! -----------------------------
!BOS

!  !IMPORT STATE:

   ! import states from child components will be promoted here 



!  !EXPORT STATE:


! Exports of child

    if(DO_DATASEAICE==0) then

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'UI',                                   &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'VI',                                   &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'VEL',                                  &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUXBOT',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUYBOT',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUXOCNB',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUYOCNB',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STROCNXB',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STROCNYB',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'FRACICE',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRENGTH',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DIVU',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DIVUOCN',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUADIV',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAU1DIV',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUODIV',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'SHEAR',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HICE',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HICE0',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HSNO',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HSNO0',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DRAFT',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DRAFT0',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'AICE',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'AICEU',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DAIDTD',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DVIDTD',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DVIRDGDT',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRCORX',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRCORY',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRTLTX',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRTLTY',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRINTX',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'STRINTY',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUXI',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUYI',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUXIB',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TAUYIB',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'UOCN',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'VOCN',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'SSH',                                  &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'SLV',                                  &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'FROCEAN',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'AREA',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TMASK',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'ICEUMASK0',                            &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'ICEUMASK1',                            &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'ANGLET',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'ANGLEU',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TRANSIX',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TRANSIY',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TRANSIMX',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'TRANSIMY',                             &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'FRESH',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'FSALT',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'FHOCN',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'SIG1',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'SIG2',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DAIDTNUDG',                            &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'DVIDTNUDG',                            &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'AICEN',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'VICEN',                                &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HIFLXE',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HIFLXN',                               &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )
    endif

!EOS

    
    if(DUAL_OCEAN) then
       ! in dual ocean mode, both Dataseaice and CICEDyna are running,
       ! but HI, SI and TI have dimensions consistent with CICE4ColumnPhys,
       ! not SimpleSeaice; the other CICE vars are copies of those from
       ! CICEDyna and they will be filled  
       call MAPL_TerminateImport    ( GC, SHORT_NAME=              &
          [character(len=9) :: 'HI'     ,     'SI'    ,            &     
                               'FRACICE',     'TI'    ,            &
                               'VOLICE' ,     'VOLSNO',            &    
                               'ERGICE ',     'ERGSNO',            & 
                               'MPOND'  ,     'TAUAGE'],           &
            CHILD=ICEd,          __RC__  )
    end if


! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__ )
! phase 1
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	  Run,        __RC__ )
! phase 2 - this is only used in the predictor part of the replay for dual ocean
    if (DUAL_OCEAN) then
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,	  Run,     __RC__ )
    end if


!=============================================================================
! Generic SetServices--This creates the generic state and calls SetServices for children
!---------------------------------------------------------------------------------------

    call MAPL_GenericSetServices    ( GC, __RC__ )

! Set the Profiling timers
! ------------------------

    call MAPL_TimerAdd(GC,   name="INITIALIZE" ,__RC__)
    call MAPL_TimerAdd(GC,   name="RUN"        ,__RC__)
    call MAPL_TimerAdd(GC,   name="--ModRun"   ,__RC__)

! All Done
!---------

    RETURN_(ESMF_SUCCESS)
  end subroutine SetServices

! -----------------------------------------------------------------

!BOP

! !IROUTINE: INITIALIZE -- Initialize method for ExternalSeaice wrapper

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

    character(len=ESMF_MAXSTR)          :: IAm
    integer                             :: STATUS
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

    CALL ESMF_UserCompSetInternalState( GC, TRIM(SEAICE_NAME)//'_internal_state', WRAP, STATUS )
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
         RC=STATUS)
    VERIFY_(STATUS)

! Once we know we have a valid  ESMF grid, we can call MAPL_GenericInitialize.
! This will allow us to use the built-in checkpoint/restart for our states.
!----------------------------------------------------------------------------


    call MAPL_TimerOff(STATE,"TOTAL"     )
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, RC=status )
    VERIFY_(STATUS)

 
    call MAPL_TimerOff(STATE,"INITIALIZE")

! All Done
!---------
    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

! ========================================================
!BOP

! !IROUTINE: Run        -- Run method for ExternalSeaice wrapper

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



! Locals

    integer           :: I,J,L
    integer           :: IM
    integer           :: JM
    integer           :: LM
    integer           :: NUM
    real, allocatable :: WGHT(:,:)
    real              :: DT, TAU_SST
    real              :: TAU_SST_UNDER_ICE
    real, pointer     :: LONS  (:,:)
    real, pointer     :: LATS  (:,:)

    integer           :: ID
    integer           :: PHASE
    real, pointer     :: FR(:,:) => null()
    real, pointer     :: FRI(:,:,:) => null()
    integer, allocatable :: PREDICTOR_CHLD(:)
    character(len=ESMF_MAXSTR) :: replayMode

! Get the component's name and set-up traceback handle.
! -----------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name, currentPhase=PHASE, RC=status )
    VERIFY_(status)
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
       where(MASK>0)
          WGHT = FROCEAN/MASK
       elsewhere
          WGHT = 0.0
       end where

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

          if (associated(HEATe)) HEATe = HEAT
          if (associated(TAUXe)) TAUXe = TAUX
          if (associated(TAUYe)) TAUYe = TAUY
          if (associated(DISCHARGEe)) DISCHARGEe = DISCHARGE
          if (associated(LWFLXe)) LWFLXe = LWFLX
          if (associated(SWFLXe)) SWFLXe = PENUVR+PENPAR+PENUVF+PENPAF+DRNIR+DFNIR
          if (associated(SHFLXe)) SHFLXe = SHFLX
          if (associated(QFLUXe)) QFLUXe = QFLUX
          if (associated(RAINe)) RAINe = RAIN
          if (associated(SNOWe)) SNOWe = SNOW
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
             where(WGHT > 0.0 )
                FRZMLTe = FRZMLT
             end where
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
                      DT*(LWFLXi+(PENUVR+PENPAR+PENUVF+PENPAF+DRNIR+DFNIR - PEN_OCN)-SHFLXi-QFLUXi*MAPL_ALHL-MAPL_ALHF*SNOWi+FHOCN)/(OrphanDepth*MAPL_RHO_SEAWATER*MAPL_CAPWTR) ! explicit update in time
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

end module GEOSseaice_GridCompMod
