

#include "MAPL_Generic.h"

module GEOS_SeaIceGridCompMod

!BOP
! !MODULE:  GEOSseaice_GridCompMod -- Implements ESMF wrapper to invoke the DATASEAICE/CICE4/CICE6 seaice models.

! !USES:

  use ESMF
  use MAPL
#ifdef BUILD_MIT_OCEAN
  use GEOS_MITDynaGridCompMod,           only : GEOSMITSeaIceSetServices  => SetServices
#endif
  use GEOS_DataSeaIceGridCompMod,        only : DataSeaIceSetServices     => SetServices
  use ice_prescribed_mod,                only : ice_nudging

   

  implicit none
  private

! !PUBLIC ROUTINES:

  public SetServices

  character(len=ESMF_MAXSTR)          :: SEAICE_NAME
  integer                             :: DO_DATASEAICE
  logical                             :: seaIceT_extData

! !DESCRIPTION:
!
!   {\tt GEOSseaice\_GridComp} is a light-weight gridded component that serves an
!   interface to seaice/data\_seaice components.
!
!EOP

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
    integer                            :: iDUAL_OCEAN
    character(len=ESMF_MAXSTR)         :: charbuf_
    character(len=ESMF_MAXSTR)         :: sharedObj

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

    call MAPL_GetResource ( MAPL,  seaIceT_extData, Label="SEAICE_THICKNESS_EXT_DATA:",  DEFAULT=.FALSE., _RC ) ! .TRUE. or .FALSE.

! Initialize these IDs (0 means not used)
! ---------------------------------------
    ICE = 0
    ICEd = 0

    if(DO_DATASEAICE/=0) then
       SEAICE_NAME="DATASEAICE"
       ICE = MAPL_AddChild(GC, NAME=SEAICE_NAME, SS=DataSeaiceSetServices, __RC__)
    else
#ifdef BUILD_MIT_OCEAN
       ICE = MAPL_AddChild(GC, NAME="MITSEAICEDYNA", SS=GEOSMITSeaIceSetServices, __RC__)
       call MAPL_AddExportSpec ( GC, SHORT_NAME = 'ICESTATES',    &
                                 CHILD_ID = ICE, __RC__  )

#else             
       call MAPL_GetResource ( MAPL, SEAICE_NAME, Label="SEAICE_NAME:", DEFAULT="CICE4", __RC__ )
       select case (trim(SEAICE_NAME))
          case ("CICE4")
             call MAPL_GetResource ( MAPL, sharedObj,  Label="GEOSCICEDyna_GridComp:", DEFAULT="libGEOSCICEDyna_GridComp.so", __RC__ )
             ICE = MAPL_AddChild(SEAICE_NAME,'setservices_', parentGC=GC, sharedObj=sharedObj,  __RC__)
          case ("CICE6")
             call MAPL_GetResource ( MAPL, sharedObj,  Label="CICE_GEOSPLUG:", DEFAULT="libCICE_GEOSPlug.so", __RC__ )
             ICE = MAPL_AddChild(SEAICE_NAME,'setservices_', parentGC=GC, sharedObj=sharedObj,  __RC__)

          case default
             charbuf_ = "SEAICE_NAME: " // trim(SEAICE_NAME) // " is not implemented, ABORT!"
             call WRITE_PARALLEL(charbuf_)
             VERIFY_(999)
       end select
#endif             
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

   if (DUAL_OCEAN) then
        
      call MAPL_AddImportSpec(GC,                                       &
         SHORT_NAME         = 'SS_FOUND',                               &
         LONG_NAME          = 'foundation_salinity_for_interface_layer',&
         UNITS              = 'PSU',                                    &
         DIMS               = MAPL_DimsHorzOnly,                        &
         VLOCATION          = MAPL_VLocationNone,                       &
         RESTART            = MAPL_RestartSkip,                         &
                                                                 __RC__ )
   endif 

!  !EXPORT STATE:
   if (DUAL_OCEAN) then

      call MAPL_AddExportSpec(GC,                                       &
         SHORT_NAME         = 'FRACICEd',                               &
         LONG_NAME          = 'fractional_cover_of_seaice',             &
         UNITS              = '1',                                      &
         DIMS               = MAPL_DimsHorzOnly,                        &
         VLOCATION          = MAPL_VLocationNone,                       &
                                                                 __RC__ )

   endif 


! Exports of child
    call MAPL_AddExportSpec ( GC   ,                              &
             SHORT_NAME = 'FRACICE',                              &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )
    call MAPL_AddExportSpec ( GC   ,                              &
             SHORT_NAME = 'UI',                                   &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

    call MAPL_AddExportSpec ( GC   ,                              &
             SHORT_NAME = 'VI',                                   &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )
    if (seaIceT_extData) then
      call MAPL_AddExportSpec ( GC   ,                          &
           SHORT_NAME = 'SEAICETHICKNESS',                      &
           CHILD_ID   = ICE ,                                   &
                                                           __RC__ )
    endif

    if(DO_DATASEAICE==0) then

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'AICE',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HICE',                                 &
             CHILD_ID   = ICE ,                                   &
                                                           __RC__ )

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'HSNO',                                 &
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
    endif

    if(trim(SEAICE_NAME) == 'CICE4') then

        call MAPL_AddExportSpec ( GC   ,                          &
             SHORT_NAME = 'VEL',                                  &
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
             SHORT_NAME = 'HICE0',                                &
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

    if (trim(SEAICE_NAME) == 'CICE6') then
       call MAPL_AddExportSpec ( GC   ,                           &
            SHORT_NAME = 'SURFSTATE',                             &
            CHILD_ID   =  ICE ,                                   &
                                                              _RC )

       call MAPL_AddExportSpec ( GC   ,                           &
            SHORT_NAME = 'TI',                                    &
            CHILD_ID   =  ICE ,                                   &
                                                              _RC )

       call MAPL_AddExportSpec ( GC   ,                           &
            SHORT_NAME = 'FRSEAICE',                              &
            CHILD_ID   =  ICE ,                                   &
                                                              _RC )
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
                               CHILD=ICEd,                 __RC__  )
    end if


! Set the Initialize, Run, Finalize entry points
! ----------------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__ )
! phase 1
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,        __RC__ )
    if (DUAL_OCEAN) then
! phase 2 - this is only used in the predictor part of the replay for dual ocean
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run,     __RC__ )
! phase 3 - this is only used in the corrector part of the replay for dual ocean
!           ice nudging only
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2,    __RC__ )
! phase 4 - this is only used in the predictor part of the replay for dual ocean
!           ice nudging only
       call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run2,    __RC__ )
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
    if (DUAL_OCEAN) then
       call MAPL_TimerAdd(GC,   name="--IceNudging"   ,__RC__)
    endif

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
    real                                :: DT

    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)


!=============================================================================

! Begin...

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Initialize"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME,  __RC__ )
    Iam = trim(comp_name) // trim(Iam)

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, State, __RC__ )

! Profilers
!----------

    call MAPL_TimerOn(STATE,"INITIALIZE")
    call MAPL_TimerOn(STATE,"TOTAL"     )

! Get info from the Generic state
!--------------------------------

    call MAPL_Get(STATE,             &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
                                          __RC__ )


! Initialize the PrivateState. First the time...
!-----------------------------------------------
    call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    __RC__)             ! Get AGCM Heartbeat
    call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, __RC__) ! set Default OCEAN_DT to AGCM Heartbeat



! Once we know we have a valid  ESMF grid, we can call MAPL_GenericInitialize.
! This will allow us to use the built-in checkpoint/restart for our states.
!----------------------------------------------------------------------------


    call MAPL_TimerOff(STATE,"TOTAL"     )
    call MAPL_GenericInitialize( GC, IMPORT, EXPORT, CLOCK, __RC__ )
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
    type (ESMF_Config)                  :: CF
    type (ESMF_GridComp    ), pointer   :: GCS(:)
    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)


    integer, parameter                  :: NUM_SNOW_LAYERS=1
    integer                             :: NUM_ICE_CATEGORIES
    integer                             :: NUM_ICE_LAYERS

    integer                             :: DO_CICE_THERMO


! Pointers to Imports

! Pointers to Exports
    real, pointer, dimension(:,:)   :: FRd       => null()


! Diagnostics exports



! Pointers to imports of child

    real, pointer, dimension(:,:,:) :: TIO8    => null()
    real, pointer, dimension(:,:,:) :: FRO8    => null()
    real, pointer, dimension(:,:,:) :: VOLICEO => null()
    real, pointer, dimension(:,:,:) :: VOLSNOO => null()
    real, pointer, dimension(:,:,:) :: TAUAGEO => null()
    real, pointer, dimension(:,:,:) :: MPONDO  => null()
    real, pointer, dimension(:,:,:) :: ERGICEO => null()
    real, pointer, dimension(:,:,:) :: ERGSNOO => null()

    real, pointer, dimension(:,:,:) :: TIO8d   => null()
    real, pointer, dimension(:,:,:) :: FRO8d   => null()
    real, pointer, dimension(:,:,:) :: VOLICEOd=> null()
    real, pointer, dimension(:,:,:) :: VOLSNOOd=> null()
    real, pointer, dimension(:,:,:) :: TAUAGEOd=> null()
    real, pointer, dimension(:,:,:) :: MPONDOd => null()
    real, pointer, dimension(:,:,:) :: ERGICEOd=> null()
    real, pointer, dimension(:,:,:) :: ERGSNOOd=> null()

! Pointers to exports of child
    real, pointer, dimension(:,:)   :: FId       => null()


! Locals

    integer           :: IM
    integer           :: JM
    real              :: DT 


    integer           :: PHASE

! Get the component's name and set-up traceback handle.
! -----------------------------------------------------

    Iam = "Run"
    call ESMF_GridCompGet( gc, NAME=comp_name,  CONFIG=CF, currentPhase=PHASE, __RC__ )
    if (PHASE >= 10) PHASE = PHASE - 10 ! to be replaced with MAPL get_phase   
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, __RC__ )

! Profilers
!----------

    call MAPL_TimerOn (STATE,"RUN"  )
    call MAPL_TimerOn (STATE,"TOTAL")

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( STATE,       DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,       DEFAULT=0, RC=STATUS)
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



! Get child's import ad export to use as a bulletin board
!--------------------------------------------------------
    call MAPL_Get(STATE,                         &
         GCS       = GCS,                        &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
         IM        = IM,                         &
         JM        = JM,                         &
                                       __RC__    )


    if (dual_ocean) then

       call MAPL_GetPointer(EXPORT   , FRd ,       'FRACICEd', __RC__)

       call MAPL_GetPointer(GIM(ICE) , FRO8     , 'FRACICE',  __RC__)
       call MAPL_GetPointer(GIM(ICE) , TIO8     , 'TI'     ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , VOLICEO  , 'VOLICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , VOLSNOO  , 'VOLSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , ERGICEO  , 'ERGICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , ERGSNOO  , 'ERGSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , TAUAGEO  , 'TAUAGE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , MPONDO   , 'MPOND'  ,  __RC__)

       call MAPL_GetPointer(GIM(ICEd), FRO8d    , 'FRACICE',  __RC__)
       call MAPL_GetPointer(GIM(ICEd), TIO8d    , 'TI'     ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), VOLICEOd , 'VOLICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), VOLSNOOd , 'VOLSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), ERGICEOd , 'ERGICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), ERGSNOOd , 'ERGSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), TAUAGEOd , 'TAUAGE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), MPONDOd  , 'MPOND'  ,  __RC__)

       call MAPL_GetPointer(GEX(ICEd), FId       , 'FRACICE'   , alloc=.TRUE., __RC__)

    end if

    call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    __RC__)             ! Get AGCM Heartbeat
    call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, __RC__) ! set Default OCEAN_DT to AGCM Heartbeat

    ! Loop the sea ice model
    !---------------------

    
    ! Run sea ice for one time step (DT)
    !---------------------------------
    
    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOn (STATE,"--ModRun")

    if (.not. DUAL_OCEAN) then
       call MAPL_GenericRunChildren(GC, IMPORT, EXPORT, CLOCK, __RC__)
    else
       if (PHASE == 1) then
          ! corrector
          call ESMF_GridCompRun( GCS(ICEd), importState=GIM(ICEd), &
               exportState=GEX(ICEd), clock=CLOCK, phase=1, userRC=STATUS)
          VERIFY_(STATUS)
          call MAPL_GenericRunCouplers( STATE, CHILD=ICEd, CLOCK=CLOCK, __RC__ )
          call ESMF_GridCompRun( GCS(ICE), importState=GIM(ICE), &
               exportState=GEX(ICE), clock=CLOCK, phase=1, userRC=STATUS)
          VERIFY_(STATUS)
          call MAPL_GenericRunCouplers( STATE, CHILD=ICE, CLOCK=CLOCK, __RC__ )

          if(associated(FRd)) then
             FRd = FId 
          endif 

       else
          ! predictor
          call ESMF_GridCompRun( GCS(ICEd), importState=GIM(ICEd), &
               exportState=GEX(ICEd), clock=CLOCK, phase=1, userRC=STATUS)
          VERIFY_(STATUS)
          call MAPL_GenericRunCouplers( STATE, CHILD=ICEd, CLOCK=CLOCK, __RC__ )

       end if

    end if

    call MAPL_TimerOff(STATE,"--ModRun")

! Profilers
!----------
    call MAPL_TimerOff(STATE,"RUN"  )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

! !IROUTINE: Run2        -- Run method for nudging sea ice in dual-ocean corrector stage

! !INTERFACE:

  subroutine Run2 ( gc, import, export, clock, rc )

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
    type (ESMF_Config)                  :: CF
    type (ESMF_GridComp    ), pointer   :: GCS(:)
    type (ESMF_State       ), pointer   :: GIM(:)
    type (ESMF_State       ), pointer   :: GEX(:)


    integer, parameter                  :: NUM_SNOW_LAYERS=1
    integer                             :: NUM_ICE_CATEGORIES
    integer                             :: NUM_ICE_LAYERS

    integer                             :: DO_CICE_THERMO


! Pointers to Imports
    real, pointer, dimension(:,:)   :: SS_FOUNDi => null()

! Pointers to Exports

! Diagnostics exports



! Pointers to imports of child

    real, pointer, dimension(:,:,:) :: TIO8    => null()
    real, pointer, dimension(:,:,:) :: FRO8    => null()
    real, pointer, dimension(:,:,:) :: VOLICEO => null()
    real, pointer, dimension(:,:,:) :: VOLSNOO => null()
    real, pointer, dimension(:,:,:) :: TAUAGEO => null()
    real, pointer, dimension(:,:,:) :: MPONDO  => null()
    real, pointer, dimension(:,:,:) :: ERGICEO => null()
    real, pointer, dimension(:,:,:) :: ERGSNOO => null()
    real, pointer, dimension(:,:)   :: AICEDO  => null()
    real, pointer, dimension(:,:)   :: HICEDO  => null()
    real, pointer, dimension(:,:)   :: HSNODO  => null()

    real, pointer, dimension(:,:,:) :: TIO8d   => null()
    real, pointer, dimension(:,:,:) :: FRO8d   => null()
    real, pointer, dimension(:,:,:) :: VOLICEOd=> null()
    real, pointer, dimension(:,:,:) :: VOLSNOOd=> null()
    real, pointer, dimension(:,:,:) :: TAUAGEOd=> null()
    real, pointer, dimension(:,:,:) :: MPONDOd => null()
    real, pointer, dimension(:,:,:) :: ERGICEOd=> null()
    real, pointer, dimension(:,:,:) :: ERGSNOOd=> null()

! Pointers to exports of child
    real, pointer, dimension(:,:)   :: FId       => null()
    real, pointer, dimension(:,:)   :: DAIDTNUDG => null()
    real, pointer, dimension(:,:)   :: DVIDTNUDG => null()



! Locals

    integer           :: IM
    integer           :: JM
    real              :: DT 
    integer           :: CAT_DIST               ! parameters for sea ice nudging
    real              :: HIN, RN,  TAU_SIT      ! parameters for sea ice nudging

    integer           :: PHASE

! Get the component's name and set-up traceback handle.
! -----------------------------------------------------

    Iam = "Run2"
    call ESMF_GridCompGet( gc, NAME=comp_name,  CONFIG=CF, currentPhase=PHASE, __RC__ )
    Iam = trim(comp_name) // Iam
    if (PHASE >= 10) PHASE = PHASE - 10 ! to be replaced with MAPL get_phase   

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, STATE, __RC__ )

! Profilers
!----------

    call MAPL_TimerOn (STATE,"RUN"  )
    call MAPL_TimerOn (STATE,"TOTAL")

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( STATE,       DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,       DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    if (DO_CICE_THERMO /= 0) then
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS,     Label="CICE_N_ICE_LAYERS:" ,     RC=STATUS)
       VERIFY_(STATUS)
    else
       NUM_ICE_CATEGORIES = 1
       NUM_ICE_LAYERS     = 1
    endif



! Get child's import ad export to use as a bulletin board
!--------------------------------------------------------
    call MAPL_Get(STATE,                         &
         GCS       = GCS,                        &
         GIM       = GIM,                        &
         GEX       = GEX,                        &
         IM        = IM,                         &
         JM        = JM,                         &
                                       __RC__    )


    if (dual_ocean) then


       call MAPL_GetPointer(GIM(ICE) , FRO8     , 'FRACICE',  __RC__)
       call MAPL_GetPointer(GIM(ICE) , TIO8     , 'TI'     ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , VOLICEO  , 'VOLICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , VOLSNOO  , 'VOLSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , ERGICEO  , 'ERGICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , ERGSNOO  , 'ERGSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , TAUAGEO  , 'TAUAGE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICE) , MPONDO   , 'MPOND'  ,  __RC__)

       call MAPL_GetPointer(GIM(ICEd), FRO8d    , 'FRACICE',  __RC__)
       call MAPL_GetPointer(GIM(ICEd), TIO8d    , 'TI'     ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), VOLICEOd , 'VOLICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), VOLSNOOd , 'VOLSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), ERGICEOd , 'ERGICE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), ERGSNOOd , 'ERGSNO' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), TAUAGEOd , 'TAUAGE' ,  __RC__)
       call MAPL_GetPointer(GIM(ICEd), MPONDOd  , 'MPOND'  ,  __RC__)


       call MAPL_GetPointer(GEX(ICE) , AICEDO    , 'AICE'    , __RC__)  ! SA: AICE needs to be looked into
       call MAPL_GetPointer(GEX(ICE) , HICEDO    , 'HICE'    , __RC__)  ! BZ: These diags need to be updated
       call MAPL_GetPointer(GEX(ICE) , HSNODO    , 'HSNO'    , __RC__)  ! after ice nudging
       call MAPL_GetPointer(GEX(ICE) , DAIDTNUDG , 'DAIDTNUDG' , alloc=.TRUE., __RC__)
       call MAPL_GetPointer(GEX(ICE) , DVIDTNUDG , 'DVIDTNUDG' , alloc=.TRUE., __RC__)
       call MAPL_GetPointer(GEX(ICEd), FId       , 'FRACICE'   , alloc=.TRUE., __RC__)

       ! copy to dataseaice imports as the ice nudging is done via dataseaice states 
       TIO8d    = TIO8
       FRO8d    = FRO8
       VOLICEOd = VOLICEO
       VOLSNOOd = VOLSNOO
       TAUAGEOd = TAUAGEO
       MPONDOd  = MPONDO
       ERGICEOd = ERGICEO
       ERGSNOOd = ERGSNOO

    end if

    call MAPL_GetResource(STATE,DT,  Label="RUN_DT:",    __RC__)             ! Get AGCM Heartbeat
    call MAPL_GetResource(STATE,DT,  Label="OCEAN_DT:",  DEFAULT=DT, __RC__) ! set Default OCEAN_DT to AGCM Heartbeat

    ! Loop the sea ice model
    !---------------------

    
    ! Run ocean for one time step (DT)
    !---------------------------------
    
    call MAPL_TimerOff(STATE,"TOTAL")
    call MAPL_TimerOn (STATE,"--IceNudging")

    
    _ASSERT(DUAL_OCEAN, 'This method should only run in dual ocean mode')

    call MAPL_GetResource( STATE, HIN     , Label="SEA_ICE_NUDGING_HINEW:"    , DEFAULT=0.5, __RC__ )
    call MAPL_GetResource( STATE, CAT_DIST, Label="SEA_ICE_NUDGING_CAT_DIST:" , DEFAULT=1  , __RC__ )

    if (PHASE == 3) then
 
        call MAPL_GetResource( STATE, TAU_SIT , LABEL="SEA_ICE_NUDGING_RELAX:"    , default=86400.0, __RC__)
        call MAPL_GetResource( STATE, RN      , Label="SEA_ICE_NUDGING_R:"        , DEFAULT=0.1,     __RC__)
        call MAPL_GetPointer(IMPORT   , SS_FOUNDi , 'SS_FOUND', __RC__)

        call ice_nudging(  FRO8d,         TIO8d,          &
                           VOLICEOd,      VOLSNOOd,       &
                           ERGICEOd,      ERGSNOOd,       &
                           TAUAGEOd,      MPONDOd,        &
                           FId,           HIN,            &
                           NUM_ICE_CATEGORIES,            &
                           TAU_SIT,       RN,             &
                           NUM_ICE_LAYERS,                &
                           NUM_SNOW_LAYERS,               &
                           CAT_DIST,      DT,             &
                           salinity = SS_FOUNDi,          &
                           ai_tend = DAIDTNUDG,           &
                           vi_tend = DVIDTNUDG )

          if(associated(AICEDO)) then
             where(AICEDO/=MAPL_UNDEF)
                AICEDO = sum(FRO8d, dim=3)
             endwhere
          endif
          if(associated(HICEDO)) then
             where(HICEDO/=MAPL_UNDEF)
                HICEDO = sum(VOLICEOd, dim=3)
             endwhere
          endif
          if(associated(HSNODO)) then
             where(HSNODO/=MAPL_UNDEF)
                HSNODO = sum(VOLSNOOd, dim=3)
             endwhere
          endif
    else  
          call ice_nudging(FRO8d,         TIO8d,          &
                           VOLICEOd,      VOLSNOOd,       &
                           ERGICEOd,      ERGSNOOd,       &
                           TAUAGEOd,      MPONDOd,        &
                           FId,           HIN,            &
                           NUM_ICE_CATEGORIES,            &
                           DT,            0.0,            &
                           NUM_ICE_LAYERS,                &
                           NUM_SNOW_LAYERS,               &
                           CAT_DIST,      DT)
    end if

    TIO8    = TIO8d
    FRO8    = FRO8d
    VOLICEO = VOLICEOd
    VOLSNOO = VOLSNOOd
    TAUAGEO = TAUAGEOd
    MPONDO  = MPONDOd
    ERGICEO = ERGICEOd
    ERGSNOO = ERGSNOOd

    call MAPL_TimerOff (STATE,"--IceNudging")

! Profilers
!----------
    call MAPL_TimerOff(STATE,"RUN"  )

! All Done
!---------

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

end module GEOS_SeaIceGridCompMod
