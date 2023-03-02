! $Id: GEOS_GcmGridComp.F90,v 1.30.2.1.4.4.2.2.2.1.2.11.6.1.4.1.2.1.2.1.14.2 2019/12/18 20:23:25 ltakacs Exp $

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GcmGridCompMod -- A Module to combine Agcm Gridded Components

! !INTERFACE:

module GEOS_GcmGridCompMod

! !USES:

   use ESMF
   use MAPL
   use dist_grid_mod

   use ModelE_MAPL,     only:  AGCM_SetServices => SetServices
   use GEOS_mkiauGridCompMod,    only:  AIAU_SetServices => SetServices
   use MAPL_HistoryGridCompMod,  only:  Hist_SetServices => SetServices
   use MAPL_HistoryGridCompMod,  only:  HISTORY_ExchangeListWrap
   use iso_fortran_env

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines the Agcm GC, 

 
!EOP

integer ::       AGCM
integer ::       AIAU
integer ::       hist

integer ::       k
character(len = 2) :: suffix

type T_GCM_STATE
   private
   type(ESMF_Alarm)           :: replayStartAlarm
   type(ESMF_Alarm)           :: replayStopAlarm
   type(ESMF_Alarm)           :: replayCycleAlarm
   type(ESMF_Alarm)           :: Regular_Replay09Alarm
   type(ESMF_Alarm)           :: replayMKIAUAlarm
   type(ESMF_Alarm)           :: replayShutoffAlarm
   integer                    :: cordur
   integer                    :: predur
   integer                    :: rplfreq
   integer                    :: rplreft
   logical                    :: rplRegular
   logical                    :: checkpointRequested = .false.
   character(len=ESMF_MAXSTR) :: checkpointFilename = ''
   character(len=ESMF_MAXSTR) :: checkpointFileType = ''
   type(ESMF_GridComp)        :: history_parent
   logical                    :: run_history = .false.
end type T_GCM_STATE

! Wrapper for extracting internal state
! -------------------------------------
type GCM_WRAP
   type (T_GCM_STATE), pointer :: PTR => null()
end type GCM_WRAP

contains

!BOP

! !IROUTINE: SetServices -- Sets ESMF services for this component

! !INTERFACE:

    subroutine SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer, intent(out)               :: RC  ! return code

! !DESCRIPTION:  The SetServices for the PhysicsGcm GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (AGCM) and runs their
!   respective SetServices.

!EOP

!=============================================================================
!
! ErrLog Variables

    character(len=ESMF_MAXSTR)              :: IAm
    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: COMP_NAME

! Locals

    type (T_GCM_STATE), pointer         :: gcm_internal_state => null() 
    type (GCM_wrap)                     :: wrap
    type (T_EXTDATA_STATE), pointer     :: extdata_internal_state 
    type (EXTDATA_wrap)                 :: ExtDataWrap
    type (MAPL_MetaComp),  pointer      :: MAPL
    character(len=ESMF_MAXSTR)          :: ReplayMode
    logical                             :: rplRegular

    type (ESMF_Config)                  :: CF
    integer                             :: heartbeat

    integer                             :: ASSIMILATION_CYCLE
    integer                             :: CORRECTOR_DURATION
    integer                             :: MKIAU_FREQUENCY
    integer                             :: MKIAU_REFERENCE_TIME
    integer                             :: REPLAY_FILE_FREQUENCY
    integer                             :: REPLAY_FILE_REFERENCE_TIME

!=============================================================================

! Begin...

! Get my name and set-up traceback handle
! ---------------------------------------

    Iam = 'SetServices'
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

! Register services for this component
! ------------------------------------

    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE, Initialize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_RUN,  Run, RC=STATUS )
    VERIFY_(STATUS)

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------

    AGCM = MAPL_AddChild(GC, NAME='AGCM', SS=Agcm_SetServices, RC=STATUS)
    VERIFY_(STATUS)
    AIAU = MAPL_AddChild(GC, NAME='AIAU', SS=AIAU_SetServices, RC=STATUS)
    VERIFY_(STATUS)

! Get RUN Parameters (MERRA-2 Defaults) and Initialize for use in other Components (e.g., AGCM_GridComp and MKIAU_GridComp)
!--------------------------------------------------------------------------------------------------------------------------

  ! Get HEARTBEAT
  ! -------------
    call MAPL_GetResource( MAPL, heartbeat, Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

  ! Get ASSIMILATION_CYCLE Duration.  If NOT_PRESENT, Initialize to MERRA-2 Default and Set Attribute.
  !---------------------------------------------------------------------------------------------------
    call MAPL_GetResource( MAPL, ASSIMILATION_CYCLE, Label="ASSIMILATION_CYCLE:", default=-999, RC=STATUS)
    VERIFY_(STATUS)
    if( ASSIMILATION_CYCLE .eq. -999 ) then
        ASSIMILATION_CYCLE = 21600
        call MAPL_ConfigSetAttribute(  CF, ASSIMILATION_CYCLE, Label="ASSIMILATION_CYCLE:", RC=STATUS)
        VERIFY_(STATUS)
        IF(MAPL_AM_I_ROOT()) PRINT *,'Setting ASSIMILATION_CYCLE to ',ASSIMILATION_CYCLE
    endif

  ! Get CORRECTOR_DURATION.  If NOT_PRESENT, Initialize to MERRA-2 Default and Set Attribute.
  ! -----------------------------------------------------------------------------------------
    call MAPL_GetResource( MAPL, CORRECTOR_DURATION, Label="CORRECTOR_DURATION:", default=-999, RC=STATUS)
    VERIFY_(STATUS)
    if( CORRECTOR_DURATION .eq. -999 ) then
        CORRECTOR_DURATION = ASSIMILATION_CYCLE
        call MAPL_ConfigSetAttribute(  CF, CORRECTOR_DURATION, Label="CORRECTOR_DURATION:", RC=STATUS)
        VERIFY_(STATUS)
        IF(MAPL_AM_I_ROOT()) PRINT *,'Setting CORRECTOR_DURATION, to ',CORRECTOR_DURATION
    endif




  ! Get REPLAY Frequency and Reference Time Associated with Forcing Files
  ! ---------------------------------------------------------------------
    call MAPL_GetResource(MAPL,REPLAY_FILE_FREQUENCY,      Label="REPLAY_FILE_FREQUENCY:",      default=-999, rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,REPLAY_FILE_REFERENCE_TIME, Label="REPLAY_FILE_REFERENCE_TIME:", default=-999, rc=STATUS )
    VERIFY_(STATUS)

  ! Get MKIAU Frequency and Reference Time Associated with MKIAU Updates
  ! --------------------------------------------------------------------
    call MAPL_GetResource(MAPL,MKIAU_FREQUENCY,      Label="MKIAU_FREQUENCY:",      default=-999, rc=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,MKIAU_REFERENCE_TIME, Label="MKIAU_REFERENCE_TIME:", default=-999, rc=STATUS )
    VERIFY_(STATUS)



  ! Set MKIAU_FREQUENCY and REPLAY_FILE_FREQUENCY Defaults
  ! ------------------------------------------------------
    if( (MKIAU_FREQUENCY .eq. -999) .and. (REPLAY_FILE_FREQUENCY .eq. -999) ) then
            REPLAY_FILE_FREQUENCY = CORRECTOR_DURATION
                  MKIAU_FREQUENCY = CORRECTOR_DURATION
            call MAPL_ConfigSetAttribute(  CF, MKIAU_FREQUENCY, Label="MKIAU_FREQUENCY:", RC=STATUS)
            VERIFY_(STATUS)
            call MAPL_ConfigSetAttribute(  CF, REPLAY_FILE_FREQUENCY, Label="REPLAY_FILE_FREQUENCY:", RC=STATUS)
            VERIFY_(STATUS)
            IF(MAPL_AM_I_ROOT()) PRINT *,'Setting MKIAU_FREQUENCY, to ',MKIAU_FREQUENCY
            IF(MAPL_AM_I_ROOT()) PRINT *,'Setting REPLAY_FILE_FREQUENCY, to ',REPLAY_FILE_FREQUENCY

    else if( MKIAU_FREQUENCY .eq. -999 ) then
             MKIAU_FREQUENCY = REPLAY_FILE_FREQUENCY
             call MAPL_ConfigSetAttribute(  CF, MKIAU_FREQUENCY, Label="MKIAU_FREQUENCY:", RC=STATUS)
             VERIFY_(STATUS)
             IF(MAPL_AM_I_ROOT()) PRINT *,'Setting MKIAU_FREQUENCY, to ',MKIAU_FREQUENCY

    else if( REPLAY_FILE_FREQUENCY .eq. -999 ) then
             REPLAY_FILE_FREQUENCY = MKIAU_FREQUENCY
             call MAPL_ConfigSetAttribute(  CF, REPLAY_FILE_FREQUENCY, Label="REPLAY_FILE_FREQUENCY:", RC=STATUS)
             VERIFY_(STATUS)
             IF(MAPL_AM_I_ROOT()) PRINT *,'Setting REPLAY_FILE_FREQUENCY, to ',REPLAY_FILE_FREQUENCY
    endif


  ! Set MKIAU_REFERENCE_TIME and REPLAY_FILE_REFERENCE_TIME Defaults
  ! ----------------------------------------------------------------
    if( (MKIAU_REFERENCE_TIME .eq. -999) .and. (REPLAY_FILE_REFERENCE_TIME .eq. -999) ) then
            REPLAY_FILE_REFERENCE_TIME = 000000
                  MKIAU_REFERENCE_TIME = 000000
            call MAPL_ConfigSetAttribute(  CF, MKIAU_REFERENCE_TIME, Label="MKIAU_REFERENCE_TIME:", RC=STATUS)
            VERIFY_(STATUS)
            call MAPL_ConfigSetAttribute(  CF, REPLAY_FILE_REFERENCE_TIME, Label="REPLAY_FILE_REFERENCE_TIME:", RC=STATUS)
            VERIFY_(STATUS)
            IF(MAPL_AM_I_ROOT()) PRINT *,'Setting MKIAU_REFERENCE_TIME, to ',MKIAU_REFERENCE_TIME
            IF(MAPL_AM_I_ROOT()) PRINT *,'Setting REPLAY_FILE_REFERENCE_TIME, to ',REPLAY_FILE_REFERENCE_TIME

    else if( MKIAU_REFERENCE_TIME .eq. -999 ) then
             MKIAU_REFERENCE_TIME = REPLAY_FILE_REFERENCE_TIME
             call MAPL_ConfigSetAttribute(  CF, MKIAU_REFERENCE_TIME, Label="MKIAU_REFERENCE_TIME:", RC=STATUS)
             VERIFY_(STATUS)
             IF(MAPL_AM_I_ROOT()) PRINT *,'Setting MKIAU_REFERENCE_TIME, to ',MKIAU_REFERENCE_TIME

    else if( REPLAY_FILE_REFERENCE_TIME .eq. -999 ) then
             REPLAY_FILE_REFERENCE_TIME = MKIAU_REFERENCE_TIME
             call MAPL_ConfigSetAttribute(  CF, REPLAY_FILE_REFERENCE_TIME, Label="REPLAY_FILE_REFERENCE_TIME:", RC=STATUS)
             VERIFY_(STATUS)
             IF(MAPL_AM_I_ROOT()) PRINT *,'Setting REPLAY_FILE_REFERENCE_TIME, to ',REPLAY_FILE_REFERENCE_TIME
    endif



! Get REPLAY Parameters
! ---------------------
    call MAPL_GetResource(MAPL, ReplayMode,  'REPLAY_MODE:',  default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

! -------------------------------------------
    ! We need to know if we are doing "regular" replay
    ! if yes, then we are modifing MKIAU import specs
    ! to perform averaging of all  MKIAU imports

    if(adjustl(ReplayMode)=="Regular") then
       rplRegular = .true.
    else    
       rplRegular = .false.
    endif

! Borrow exports from children

! Export for IAU and/or Analysis purposes
! ---------------------------------------
       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'AK',       &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'BK',       &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'PS',       &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'DELP',     &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TV',       &
            CHILD_ID = AGCM,         &
            RC = STATUS              )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'U',        &
            CHILD_ID = AGCM,         &
            RC = STATUS              )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'V',        &
            CHILD_ID = AGCM,         &
            RC = STATUS              )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'PHIS',     &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'Q',        &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'QCTOT',    &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'QITOT',    &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'QLTOT',    &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'O3PPMV',   &
            CHILD_ID   = AGCM,       &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'U10N',     &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'V10N',     &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'SNOMAS',   &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'WET1',     &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TSOIL1',   &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'LWI',      &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'TS',       &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'Z0',       &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRLAND',   &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRLANDICE',&
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRLAKE',   &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FROCEAN',  &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddExportSpec ( GC, &
            SHORT_NAME = 'FRACI',    &
            CHILD_ID = AGCM,         &
            RC=STATUS                )
       VERIFY_(STATUS)

       call MAPL_AddConnectivity ( GC,                              &
            SHORT_NAME  = (/'PHIS  ','AK    ','BK    ','U     ','V     ','TV    ','PS    ','DELP  ','O3PPMV', 'TS    ','AREA  '/),       &
            DST_ID = AIAU,                                          &
            SRC_ID = AGCM,                                          &
            RC=STATUS  )
       VERIFY_(STATUS)

       call MAPL_AddConnectivity ( GC,                              &
            SHORT_NAME  = (/'DUDT ','DVDT ','DTDT ','DQVDT','DPEDT','DO3DT','DTSDT'/), &
            DST_ID = AGCM,                                          &
            SRC_ID = AIAU,                                          &
            RC=STATUS  )
       VERIFY_(STATUS)

       call MAPL_AddConnectivity ( GC,                              &
            SRC_NAME  = (/'Q            ','TROPP_BLENDED'/),        &
            DST_NAME  = (/'QV           ','TROPP_BLENDED'/),        &
            DST_ID = AIAU,                                          &
            SRC_ID = AGCM,                                          &
            RC=STATUS  )
       VERIFY_(STATUS)

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( gcm_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => gcm_internal_state

    gcm_internal_state%rplRegular = rplRegular
    gcm_internal_state%cordur     = CORRECTOR_DURATION
    gcm_internal_state%rplfreq    = MKIAU_FREQUENCY
    gcm_internal_state%rplreft    = MKIAU_REFERENCE_TIME

! If doing regular replay make state to "borrow" gc and import to ExtData
! -----------------------------------------------------------------------
    if (rplRegular) then
       allocate( extdata_internal_state, stat=status)
       VERIFY_(STATUS)
       ExtDataWrap%ptr => extdata_internal_state
       call ESMF_UserCompSetInternalState( GC, 'ExtData_state',ExtDataWrap,status)
       VERIFY_(STATUS)
      
       call history_setservice(_RC)
    end if

 
! Save pointer to the wrapped internal state in the GC
! ----------------------------------------------------

    call ESMF_UserCompSetInternalState ( GC, 'GCM_state',wrap,status )
    VERIFY_(STATUS)


    call MAPL_GenericSetServices    ( GC, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC, name="INITIALIZE"    ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="RUN"           ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--OCEAN"       ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--ATMOSPHERE"  ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--REPLAY"      ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--A2O"         ,RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_TimerAdd(GC, name="--O2A"         ,RC=STATUS)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)

    contains

      subroutine history_setservice(rc)
         integer, intent(out), optional :: rc

         integer :: status
         type(ESMF_Config) :: hist_cf, gcm_cf
         type(MAPL_MetaComp), pointer :: history_metaobj
         type(StubComponent) :: stub_component
         integer :: run_dt
         character(len=ESMF_MAXSTR) :: replay_history
         logical :: is_present

         call ESMF_GridCompGet(gc,config=gcm_cf,_RC)
         call ESMF_ConfigFindLabel(gcm_cf,"REPLAY_HISTORY_RC:",isPresent=is_present,_RC)

         if (is_present) then
            gcm_internal_state%run_history = .true. 
            call MAPL_GetResource(MAPL,replay_history,"REPLAY_HISTORY_RC:",_RC)
            hist_cf = ESMF_ConfigCreate(_RC)
            call ESMF_ConfigLoadFile(hist_cf,trim(replay_history),_RC)
            call MAPL_GetResource(MAPL,run_dt,"RUN_DT:",_RC)
            call MAPL_ConfigSetAttribute(hist_cf,value=run_dt,label="RUN_DT:",_RC)
            call MAPL_ConfigSetAttribute(hist_cf,value=replay_history,label="HIST_CF:",_RC)
            gcm_internal_state%history_parent = ESMF_GridCompCreate(name="History_GCM_parent",config=hist_cf,_RC)
            history_metaobj => null()
            call MAPL_InternalStateCreate(gcm_internal_state%history_parent,history_metaobj,_RC)
            call MAPL_Set(history_metaobj,cf=hist_cf,name="History_GCM_parent",component=stub_component,_RC)
            hist = MAPL_AddChild(history_metaobj,name="History_GCM",ss=hist_setservices,_RC) 
         end if
         _RETURN(_SUCCESS)
      end subroutine history_setservice

  end subroutine SetServices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Initialize -- Initialize method for the composite Gcm Gridded Component

! !INTERFACE:

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)           :: IAm 
    integer                              :: STATUS
    character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases

    type(ESMF_Grid)     :: agrid

    type(ESMF_DELayout) :: layout
    type (ESMF_Config)  :: CF

    type (ESMF_VM)                      :: VM
    type (MAPL_MetaComp),      pointer  :: MAPL => null()
    type (MAPL_MetaComp),      pointer  :: CMAPL => null()
    type (ESMF_GridComp),      pointer  :: GCS(:) => null()
    type (ESMF_State),         pointer  :: GIM(:) => null()
    type (ESMF_State),         pointer  :: GEX(:) => null()
    character(len=ESMF_MAXSTR), pointer :: GCNAMES(:) => null()
    character(len=ESMF_MAXSTR)          :: tmpStr
    character(len=ESMF_MAXSTR)          :: ReplayMode
    type (T_GCM_STATE), pointer         :: gcm_internal_state => null() 
    type (GCM_wrap)                     :: wrap
    
    type(ESMF_Calendar)                 :: cal
    type(ESMF_Time)                     :: rep_StartTime
    type(ESMF_Time)                     :: rep_RefTime
    type(ESMF_Time)                     :: currTime
    type(ESMF_Time)                     :: ringTime
    type(ESMF_TimeInterval)             :: CORRECTOR_DURATION
    type(ESMF_TimeInterval)             :: PREDICTOR_DURATION
    type(ESMF_TimeInterval)             :: MKIAU_FREQUENCY
    type(ESMF_TimeInterval)             :: Shutoff
    type(ESMF_Alarm)                    :: replayStartAlarm
    type(ESMF_Alarm)                    :: replayStopAlarm
    type(ESMF_Alarm)                    :: replayCycleAlarm
    type(ESMF_Alarm)                    ::   Exact_Replay09Alarm
    type(ESMF_Alarm)                    :: Regular_Replay09Alarm
    type(ESMF_Alarm)                    :: replayMKIAUAlarm
    type(ESMF_Alarm)                    :: replayShutoffAlarm
    TYPE(ESMF_Alarm)                    :: PredictorIsActive
    type(ESMF_Alarm)                    :: alarms(1)
    type(ESMF_DistGrid)                 :: distGRID
    integer                             :: i_CORRECTOR_DURATION
    integer                             :: i_PREDICTOR_DURATION
    integer                             :: i_MKIAU_FREQUENCY
    integer                             ::   MKIAU_REFERENCE_TIME
    integer                             :: rplshut
    integer                             :: rep_startdate(2)
    integer                             :: rep_YY
    integer                             :: rep_MM
    integer                             :: rep_DD
    integer                             :: rep_H
    integer                             :: rep_M
    integer                             :: rep_S
    integer                             :: MKIAU_RingDate
    integer                             :: MKIAU_RingTime
    logical                             :: FileExists

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Initialize"

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"INITIALIZE")
    call MAPL_TimerOn(MAPL,"TOTAL")

! Get my internal private state. This contains the transforms
!  between the exchange grid and the atmos grid.
!-------------------------------------------------------------
    call ESMF_UserCompGetInternalState(gc, 'GCM_state', wrap, status)
    VERIFY_(STATUS)
    gcm_internal_state => wrap%ptr


! Get children and their im/ex states from my generic state.
!----------------------------------------------------------
    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCNAMES=GCNAMES, RC=STATUS )
    VERIFY_(STATUS)


! Create Atmospheric grid
!------------------------
!   this is very klugy, but if we are running modele it has its way of making grid distro
!   and we have ours and ours is easier to make match modele's than the other way around
    call write_modele_distribution_to_file(cf,_RC)
    call MAPL_GridCreate(GCS(AGCM), rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompGet(GCS(AGCM),  grid=agrid, rc=status)
    VERIFY_(STATUS)

    call ESMF_GridCompSet(GC, grid=agrid, rc=status)
    VERIFY_(STATUS)

    GCM_INTERNAL_STATE%checkpointRequested = .false.

    call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)

! ALT replay: create alarms
! -------------------------
       i_CORRECTOR_DURATION         = gcm_internal_state%corDur

       i_MKIAU_FREQUENCY      = gcm_internal_state%rplfreq
         MKIAU_REFERENCE_TIME = gcm_internal_state%rplreft

       call MAPL_GetResource(MAPL, rplshut, Label="REPLAY_SHUTOFF:", default=-3600, rc=status)
       VERIFY_(STATUS)
       call ESMF_ClockGet(clock, currTime=currTime, calendar=cal, rc=status)
       VERIFY_(STATUS)

       call ESMF_TimeIntervalSet(CORRECTOR_DURATION,    S=i_CORRECTOR_DURATION, rc=status)
       VERIFY_(STATUS)
       call ESMF_TimeIntervalSet(MKIAU_FREQUENCY, S=i_MKIAU_FREQUENCY, rc=status)
       VERIFY_(STATUS)

       call ESMF_TimeIntervalSet(Shutoff, S=abs(rplshut), rc=status)
       VERIFY_(STATUS)

       rep_StartTime = currTime

       ! UNPACK REPLAY_STARTTIME
       ! -----------------------
       call ESMF_TimeGet  ( rep_StartTime, YY = rep_YY, &
                                           MM = rep_MM, &
                                           DD = rep_DD, &
                                           H  = rep_H , &
                                           M  = rep_M , &
                                           S  = rep_S, rc=status )
       VERIFY_(STATUS)

       ! Initialize REPLAY_STARTDATE in YYYYMMDD & HHMMSS format
       ! ------------------------------------------------------- 
       rep_startdate(1) = 10000*rep_YY + 100*rep_MM + rep_DD
       rep_startdate(2) = 10000*rep_H  + 100*rep_M  + rep_S 

       ! Update REPLAY_STARTDATE with User-Supplied values from CONFIG
       ! ------------------------------------------------------------- 
       call MAPL_GetResource( MAPL, rep_startdate(1), label='REPLAY_STARTDATE:', default=rep_startdate(1), rc=STATUS )
       call MAPL_GetResource( MAPL, rep_startdate(2), label='REPLAY_STARTTIME:', default=rep_startdate(2), rc=STATUS )

       ! REPACK REPLAY_STARTTIME
       ! -----------------------
       rep_YY =     rep_startdate(1)/10000
       rep_MM = mod(rep_startdate(1),10000)/100
       rep_DD = mod(rep_startdate(1),100)
       rep_H  =     rep_startdate(2)/10000
       rep_M  = mod(rep_startdate(2),10000)/100
       rep_S  = mod(rep_startdate(2),100)

       call ESMF_TimeSet( rep_StartTime, YY = rep_YY, &
                                         MM = rep_MM, &
                                         DD = rep_DD, &
                                          H = rep_H , &
                                          M = rep_M , &
                                          S = rep_S , &
                          calendar=cal,  rc = STATUS  )
       VERIFY_(STATUS)

       rep_H =     MKIAU_REFERENCE_TIME/10000
       rep_M = mod(MKIAU_REFERENCE_TIME,10000)/100
       rep_S = mod(MKIAU_REFERENCE_TIME,100)

     ! Set ESMF Replay Reference Time
     ! ------------------------------
       call ESMF_TimeSet( rep_RefTime,  YY = rep_YY, &
                                        MM = rep_MM, &
                                        DD = rep_DD, &
                                         H = rep_H , &
                                         M = rep_M , &
                                         S = rep_S , &
                          calendar=cal, rc = STATUS  )
       VERIFY_(STATUS)


     ! Compute Time of First Call to MKIAU
     ! -----------------------------------
       if (rep_RefTime < currTime ) then
           rep_RefTime = rep_RefTime + ( INT( (currTime-rep_RefTime)/MKIAU_FREQUENCY ) )*MKIAU_FREQUENCY
       if (rep_RefTime < currTime ) then
           rep_RefTime = rep_RefTime + MKIAU_FREQUENCY
       endif
       endif
       replayMKIAUAlarm = ESMF_AlarmCreate( name="ReplayMKIAU", clock=clock,                                &
                                            RingTime     = rep_RefTime,                                     &
                                            RingInterval = MKIAU_FREQUENCY, sticky=.false., rc=status )
       VERIFY_(STATUS)
       if(rep_RefTime == currTime) then
          call ESMF_AlarmRingerOn( replayMKIAUAlarm, rc=status )
          VERIFY_(STATUS)
       end if


     ! Create Alarms for REPLAYs Requiring 09 files
     ! --------------------------------------------
              ! EXACT REPLAY may require Correct_Step: agcm09_import
              ! ----------------------------------------------------
                Exact_Replay09Alarm = ESMF_AlarmCreate( name="ExactReplay09", clock=clock,                                            &
                                                        RingTime     = rep_StartTime + CORRECTOR_DURATION - MKIAU_FREQUENCY/2 , &
                                                        RingInterval = CORRECTOR_DURATION, sticky=.false., rc=status )
                VERIFY_(STATUS)

              ! REGULAR REPLAY may require Predictor_Step: ana09.eta 
              ! ----------------------------------------------------
                Regular_Replay09Alarm = ESMF_AlarmCreate( name="RegularReplay09", clock=clock,               &
                                                          RingTime     = rep_StartTime + CORRECTOR_DURATION, &
                                                          RingInterval = CORRECTOR_DURATION, sticky=.false., rc=status )
                VERIFY_(STATUS)


     ! Compute PREDICTOR Duration
     ! --------------------------
       PREDICTOR_DURATION = rep_RefTime - currTime
                 RingTime = rep_RefTime
       do while( RingTime + MKIAU_FREQUENCY <= rep_StartTime + CORRECTOR_DURATION )
                 RingTime  = RingTime + MKIAU_FREQUENCY
       PREDICTOR_DURATION  = PREDICTOR_DURATION + MKIAU_FREQUENCY
       end do


     ! Check for User-Defined PREDICTOR_DURATION, and Check for Consistency with Computed Value
     ! ----------------------------------------------------------------------------------------
       call MAPL_GetResource( MAPL, i_PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", default=-999, RC=STATUS)
       VERIFY_(STATUS)
       if( i_PREDICTOR_DURATION .eq. -999 ) then
           call ESMF_TimeIntervalGet( PREDICTOR_DURATION, S=i_PREDICTOR_DURATION, rc=status )
           call MAPL_ConfigSetAttribute(  CF, i_PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", RC=STATUS)
           VERIFY_(STATUS)
           IF(MAPL_AM_I_ROOT()) PRINT *,'Setting PREDICTOR_DURATION based on IAU Frequency: ',i_PREDICTOR_DURATION
       else
           call ESMF_TimeIntervalGet( PREDICTOR_DURATION, S=rep_S, rc=status )
           if( rep_S .gt. i_PREDICTOR_DURATION ) then
               IF(MAPL_AM_I_ROOT()) then
                  PRINT * 
                  PRINT *,'ERROR!             User-Defined PREDICTOR_DURATION: ',i_PREDICTOR_DURATION
                  PRINT *,'is SMALLER THAN required Calculated PREDICTOR_DURATION: ',rep_S
                  PRINT * 
               endif
               VERIFY_(ESMF_FAILURE)
           else
               if( rep_S .ne. i_PREDICTOR_DURATION ) then
                   call ESMF_TimeIntervalSet( PREDICTOR_DURATION, S=i_PREDICTOR_DURATION, rc=status )
                   VERIFY_(STATUS)
                   IF(MAPL_AM_I_ROOT()) PRINT *,' (Note: PREDICTOR_DURATION based on USER Input) '
               endif
           endif
       endif
       gcm_internal_state%predur = i_PREDICTOR_DURATION

       ! NOTE:  When using PREDICTOR_DURATION=0 or calling MKIAU at Beginning of Run,
       !        you must include AIAU_IMPORT_RESTARTS to begin the cycle
       ! ----------------------------------------------------------------------------
       if( (adjustl(ReplayMode)=="Regular") .and.                         &
           ( (i_PREDICTOR_DURATION == 0) .or. (rep_RefTime == currTime) ) ) then
           call MAPL_GetResource( MAPL, tmpStr, LABEL="AIAU_IMPORT_RESTART_FILE:", RC=STATUS)
           if (STATUS /= ESMF_SUCCESS) then
               IF(MAPL_AM_I_ROOT()) then
                                                  PRINT *,' '
                                                  PRINT *,'ERROR! ...'
                  if( i_PREDICTOR_DURATION == 0 ) PRINT *,'Using PREDICTOR_DURATION = 0 requires AIAU_IMPORT_RESTART'
                  if( rep_RefTime == currTime   ) PRINT *,'Calling MKIAU at Beginning of Run requires AIAU_IMPORT_RESTART'
                                                  PRINT *,' '
               endif
               VERIFY_(STATUS)
           else
               call ESMF_VMGetCurrent ( VM=vm, RC=STATUS )
               VERIFY_(STATUS)
               IF(MAPL_AM_I_ROOT()) then
                  Inquire(FILE = tmpStr, EXIST=FileExists)
               ENDIF
               call MAPL_CommsBcast(vm, FileExists, n=1, ROOT=MAPL_Root, rc=status)
               VERIFY_(status)
               if( .not.FileExists ) then
                   IF(MAPL_AM_I_ROOT()) then
                                                      PRINT *,' '
                                                      PRINT *,'ERROR! ...'
                      if( i_PREDICTOR_DURATION == 0 ) PRINT *,'Using PREDICTOR_DURATION = 0 requires AIAU_IMPORT_RESTART'
                      if( rep_RefTime == currTime   ) PRINT *,'Calling MKIAU at Beginning of Run requires AIAU_IMPORT_RESTART'
                                                      PRINT *,' '
                   endif
                   RETURN_(ESMF_FAILURE)
               endif
           endif
       endif


       ! Creating REPLAY Alarms
       ! ----------------------
       RingTime = rep_StartTime

       replayStartAlarm = ESMF_AlarmCreate( name="startReplay", clock=clock, &
                                            RingTime     = ringTime,         & 
                                            RingInterval = CORRECTOR_DURATION, sticky=.true., rc=status )
       VERIFY_(STATUS)

       if(ringTime == currTime) then
          call ESMF_AlarmRingerOn( replayStartAlarm, rc=status )
          VERIFY_(STATUS)
       end if

       replayCycleAlarm = ESMF_AlarmCreate( name="replayCycle", clock=clock, &
                                            RingTime     = ringTime,                                     &
                                            RingInterval = CORRECTOR_DURATION, sticky=.false., rc=status )
       VERIFY_(STATUS)

       replayStopAlarm = ESMF_AlarmCreate( name="stopReplay", clock=clock, &
                                           RingTime     = ringTime + PREDICTOR_DURATION,                &
                                           RingInterval = CORRECTOR_DURATION, sticky=.false., rc=status )
       VERIFY_(STATUS)

       replayShutoffAlarm = ESMF_AlarmCreate( name='ReplayShutOff', clock=clock, RingInterval = Shutoff, sticky=.true., RC=STATUS )
       VERIFY_(STATUS)


     ! Put MKIAU_RingDate and MKIAU_RingTime into MAPL Config
     ! ------------------------------------------------------
       call MAPL_GetResource(MAPL, MKIAU_RingDate, Label="MKIAU_RingDate:", default=-999, RC=STATUS)
       VERIFY_(STATUS)
       if( MKIAU_RingDate .eq. -999 ) then
           call ESMF_TimeGet  ( rep_RefTime, YY = rep_YY, &
                                             MM = rep_MM, &
                                             DD = rep_DD, &
                                             H  = rep_H , &
                                             M  = rep_M , &
                                             S  = rep_S, rc=status )
           VERIFY_(STATUS)

           MKIAU_RingDate = 10000*rep_YY + 100*rep_MM + rep_DD
           MKIAU_RingTime = 10000*rep_H  + 100*rep_M  + rep_S 

           call MAPL_ConfigSetAttribute(  CF, MKIAU_RingDate, Label="MKIAU_RingDate:", RC=STATUS)
           VERIFY_(STATUS)
           call MAPL_ConfigSetAttribute(  CF, MKIAU_RingTime, Label="MKIAU_RingTime:", RC=STATUS)
           VERIFY_(STATUS)

           IF(MAPL_AM_I_ROOT()) PRINT *,'Setting MKIAU_RingDate/Time to ',MKIAU_RingDate,MKIAU_RingTime
       else
           IF(MAPL_AM_I_ROOT()) PRINT *,'Setting MKIAU_RingDate and MKIAU_RingTime not allowed!'
           VERIFY_(ESMF_FAILURE)
       endif

       if (rplshut <= 0) then ! this is a "flag" to never use Shutoff alarm
          call ESMF_AlarmDisable(replayShutoffAlarm, RC=STATUS)
          VERIFY_(STATUS)
       end if

       ! Create Active PREDICTOR_STEP Alarm (Note: Sticky, so OFF & ON are set Explicitly)
       ! ---------------------------------------------------------------------------------
       PredictorIsActive = ESMF_AlarmCreate(NAME="PredictorActive", CLOCK=clock, RingTime=currTime, STICKY=.TRUE., RC=STATUS)
       VERIFY_(STATUS)

       ! Set Active PREDICTOR_STEP Alarm to OFF
       ! --------------------------------------
       CALL ESMF_AlarmRingerOff(PredictorIsActive, RC=STATUS)
       VERIFY_(STATUS)
     ! IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Predictor Alarm is ringing? ",ESMF_AlarmIsRinging(PredictorIsActive)
       

    if(gcm_internal_state%rplRegular) then

       call MAPL_GetObjectFromGC ( GCS(AGCM), CMAPL, RC=STATUS)
       VERIFY_(STATUS)
       Alarms(1) = replayStartAlarm
       call MAPL_AddRecord(CMAPL, ALARMS, (/MAPL_Write2Ram/), rc=status)
       VERIFY_(STATUS)

       GCM_INTERNAL_STATE%replayStartAlarm      = replayStartAlarm
       GCM_INTERNAL_STATE%replayStopAlarm       = replayStopAlarm
       GCM_INTERNAL_STATE%replayCycleAlarm      = replayCycleAlarm
       GCM_INTERNAL_STATE%replayMKIAUAlarm      = replayMKIAUAlarm
       GCM_INTERNAL_STATE%replayShutoffAlarm    = replayShutoffAlarm
       GCM_INTERNAL_STATE%Regular_Replay09Alarm = Regular_Replay09Alarm

       call MAPL_GetResource(MAPL, tmpStr, "MKIAU_CHECKPOINT_FILE:", default="NULL", rc=status)
       VERIFY_(STATUS)

       if (tmpStr /= "NULL" .or. (PREDICTOR_DURATION.gt.MKIAU_FREQUENCY/2) ) then
          GCM_INTERNAL_STATE%checkpointRequested = .true.
          if (tmpStr /= "NULL" ) then
              GCM_INTERNAL_STATE%checkpointFilename = tmpStr
              call MAPL_GetResource(MAPL, tmpStr, "MKIAU_CHECKPOINT_TYPE:", default="pnc4", rc=status)
              VERIFY_(STATUS)
              GCM_INTERNAL_STATE%checkpointFileType = tmpStr
          else
              GCM_INTERNAL_STATE%checkpointFilename = "mkiau_checkpoint"
              GCM_INTERNAL_STATE%checkpointFileType = "pnc4"
              call MAPL_ConfigSetAttribute(  CF, trim(GCM_INTERNAL_STATE%checkpointFilename), Label="MKIAU_CHECKPOINT_FILE:", RC=STATUS)
              VERIFY_(STATUS)
              call MAPL_ConfigSetAttribute(  CF, trim(GCM_INTERNAL_STATE%checkpointFileType), Label="MKIAU_CHECKPOINT_TYPE:", RC=STATUS)
              VERIFY_(STATUS)
              IF(MAPL_AM_I_ROOT()) then
                 PRINT *,'Setting MKIAU_CHECKPOINT_FILE to ',trim(GCM_INTERNAL_STATE%checkpointFilename)
                 PRINT *,'Setting MKIAU_CHECKPOINT_TYPE to ',trim(GCM_INTERNAL_STATE%checkpointFileType)
              ENDIF
          endif
       endif

    end if

! **********************************************************************
! ****                 Initialize Gridded Components                ****
! **********************************************************************

    call MAPL_TimerOff(MAPL,"TOTAL")

    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)

    call MAPL_TimerOn(MAPL,"TOTAL")

    if(gcm_internal_state%rplRegular .and. gcm_internal_state%run_history) then
       call initialize_history(_RC)
    end if

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"INITIALIZE")

   RETURN_(ESMF_SUCCESS)

 contains

   subroutine initialize_history(rc)
     integer, optional, intent(Out) :: rc

     integer :: status,user_status
     type(ESMF_State), allocatable :: gcm_exports(:),hist_imports(:),hist_exports(:)
     type(ESMF_GridComp), allocatable :: hist_gcs(:)
     type(ESMF_GridComp), allocatable :: gcm_gcs(:)
     type(MAPL_MetaComp), pointer :: history_metaobj
     type(HISTORY_ExchangeListWrap) :: lswrap
     integer(kind=INT64), pointer   :: LSADDR(:) => null()

     call MAPL_GetObjectFromGC ( gcm_internal_state%history_parent, history_metaobj, _RC)

     call MAPL_Get(mapl,childrens_export_states = gcm_exports, childrens_gridcomps = gcm_gcs,  _RC)

     call MAPL_Get(history_metaobj, &
          childrens_export_states = hist_exports, &
          childrens_import_states = hist_imports, &
          childrens_gridcomps = hist_gcs, _RC)

     allocate(lswrap%ptr, stat = status)
     _VERIFY(STATUS)
     call ESMF_UserCompSetInternalState(hist_gcs(hist), 'MAPL_LocStreamList', &
          lswrap, _RC)
     call MAPL_GetAllExchangeGrids(gcm_gcs(agcm), LSADDR, _RC)
     lswrap%ptr%LSADDR_PTR => LSADDR

     call ESMF_StateAdd(hist_imports(hist),[gcm_exports(agcm)],_RC)
     call ESMF_GridCompInitialize(hist_gcs(hist),importState=hist_imports(hist),&
                                       exportState=hist_exports(hist),&
                                       clock=clock,userRC=user_status,_RC)

     _RETURN(_SUCCESS)
   end subroutine initialize_history

   subroutine AllocateExports(STATE, NAMES, RC)
     type(ESMF_State)          , intent(INOUT) ::  STATE
     character(len=*)          , intent(IN   ) ::  NAMES(:)
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="AllocateExports"
     integer                               :: STATUS

     real,    pointer :: ptr(:) => null()
     integer          :: I

     DO I = 1, SIZE(NAMES)
        call MAPL_GetPointer(STATE, ptr, NAMES(I), alloc=.true.,  RC=STATUS)
        VERIFY_(STATUS)
     END DO

     RETURN_(ESMF_SUCCESS)

   end subroutine AllocateExports
   subroutine AllocateExports_UGD(STATE, NAMES, RC)
     type(ESMF_State)          , intent(INOUT) ::  STATE
     character(len=*)          , intent(IN   ) ::  NAMES(:)
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="AllocateExports_UGD"
     integer                               :: STATUS

     real,    pointer :: ptr(:,:) => null()
     integer          :: I

     DO I = 1, SIZE(NAMES)
        call MAPL_GetPointer(STATE, ptr, NAMES(I), alloc=.true.,  RC=STATUS)
        VERIFY_(STATUS)
     END DO

     RETURN_(ESMF_SUCCESS)

   end subroutine AllocateExports_UGD

 end subroutine Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !IROUTINE: Run -- Run method for the composite Gcm Gridded Component

! !INTERFACE:

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

! !DESCRIPTION: 
 

!EOP

! ErrLog Variables

    character(len=ESMF_MAXSTR)           :: IAm 
    integer                              :: STATUS
    character(len=ESMF_MAXSTR)           :: COMP_NAME

! Local derived type aliases
    type(ESMF_VM)       :: VM

    type (MAPL_MetaComp),  pointer  :: MAPL => null()
    type (MAPL_MetaComp),  pointer  :: MAPL_AGCM
    type (ESMF_GridComp),      pointer  :: GCS(:) => null()
    type (ESMF_State),         pointer  :: GIM(:) => null()
    type (ESMF_State),         pointer  :: GEX(:) => null()

    type (T_GCM_STATE), pointer         :: gcm_internal_state => null() 
    type (GCM_wrap)                     :: wrap
    type (T_ExtData_STATE), pointer     :: ExtData_internal_state  => null()
    type (ExtData_wrap)                 :: ExtDatawrap
    type (ESMF_State)                   :: Dummy

    type(ESMF_Alarm)              :: alarm
    type(ESMF_Alarm), allocatable :: AlarmList(:)
    type(ESMF_Time),  allocatable :: AlarmRingTime(:)
    type(ESMF_Time)               :: ct, replayTime
    type (ESMF_Time)              :: REPLAY_TIME

    logical,          allocatable :: ringingState(:)
    integer                       :: i, nalarms
    logical                       :: done
    logical                       :: shutoffRpl
    logical                       :: timeForRpl
    logical                       :: timeForMKIAU
    character(len=ESMF_MAXSTR)    :: FileName
    character(len=ESMF_MAXSTR)    :: FileType
    character(len=ESMF_MAXSTR)    :: ReplayMode
    character(len=ESMF_MAXSTR)    :: record_fname
    character(len=14)             :: DATESTAMP !YYYYMMDD_HHMMz

    integer                                :: N

    TYPE(ESMF_Alarm)              :: PredictorIsActive

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // "Run"

    call ESMF_VMGetCurrent ( VM=vm, RC=STATUS )
    VERIFY_(STATUS)
    
! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)


    call MAPL_TimerON(MAPL,"RUN"  )
    call MAPL_TimerON(MAPL,"TOTAL")


! Get my internal private state. This contains the transforms
!  between the exchange grid and the atmos grid.
!-------------------------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'GCM_state', wrap, status)
    VERIFY_(STATUS)
    gcm_internal_state => wrap%ptr

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, RC=STATUS )
    VERIFY_(STATUS)

    ! Check for Default DO_DATAATM=0 (FALSE) mode
    ! -------------------------------------------

       call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
       VERIFY_(STATUS)

 REPLAY: if(adjustl(ReplayMode)=="Regular") then

          call MAPL_TimerON(MAPL,"--REPLAY"  )

          call ESMF_UserCompGetInternalState(gc, 'ExtData_state', ExtDatawrap, status)
          VERIFY_(STATUS)
          ExtData_internal_state => ExtDatawrap%ptr

          shutoffRpl = ESMF_AlarmIsRinging(GCM_INTERNAL_STATE%replayShutoffAlarm, RC=STATUS)
          VERIFY_(STATUS)
          timeForRpl = ESMF_AlarmIsRinging(GCM_INTERNAL_STATE%replayStartAlarm, RC=STATUS)
          VERIFY_(STATUS)
          if (timeForRpl .and. .not. shutoffRpl ) then

             ! clear Atm IAU tendencies first
             ! ------------------------------
             if( MAPL_AM_I_Root() ) then
                 PRINT *
                 PRINT *,TRIM(Iam)//":  Zeroing IAU forcing ..."
             endif
             call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=2, userRC=status )
             VERIFY_(STATUS)

             ! we should save import and internal states
             ! this is already done by record in MAPL_CAP

             call ESMF_ClockGet(clock, currTime=ct, rc=status)
             VERIFY_(STATUS)

             ! Set Active PREDICTOR_STEP Alarm to ON
             ! -------------------------------------
             CALL ESMF_ClockGetAlarm(CLOCK, "PredictorActive", PredictorIsActive, RC=STATUS)
             VERIFY_(STATUS)
             CALL ESMF_AlarmRingerOn(PredictorIsActive, RC=STATUS)
             VERIFY_(STATUS)

             ! save alarms' states
             ! -------------------
             call ESMF_ClockGet(clock, alarmCount = nalarms, rc=status)
             VERIFY_(STATUS)
             allocate (alarmList(nalarms), alarmRingTime(nalarms), ringingState(nalarms), stat = status)
             VERIFY_(STATUS)
             call ESMF_ClockGetAlarmList(clock, alarmListFlag=ESMF_ALARMLIST_ALL, alarmList=alarmList, alarmCount = nalarms, rc=status)
             VERIFY_(STATUS)
             DO I = 1, nalarms
                call ESMF_AlarmGet(alarmList(I), ringTime=alarmRingTime(I), ringing=ringingState(I), rc=status)
                VERIFY_(STATUS)
             END DO

             ! fix the predictor alarm (huh?)
             ! ------------------------------
             call ESMF_ClockGetAlarm(clock, 'PredictorAlarm', alarm, rc=status)
             VERIFY_(STATUS)
             call ESMF_AlarmRingerOn(alarm, RC=STATUS)
             VERIFY_(STATUS)


             replayTime = ct

           ! Check Need for MKIAU at beginning of REPLAY Time Loop
           ! -----------------------------------------------------
                timeForMKIAU = ESMF_AlarmIsRinging( GCM_INTERNAL_STATE%replayMKIAUAlarm, RC=STATUS )
                VERIFY_(STATUS)

                if( timeForMKIAU ) then

                    ! Ensure Regular_Replay09Alarm is OFF at Beginning of PREDICTOR Loop
                    ! ------------------------------------------------------------------
                    call ESMF_AlarmRingerOff(GCM_INTERNAL_STATE%Regular_Replay09Alarm, RC=STATUS)
                    VERIFY_(STATUS)

                    if( MAPL_AM_I_Root() ) then
                        PRINT *,TRIM(Iam)//":  Running MKIAU at Beginning of Predictor Loop, Creating IAU forcing ..."
                        PRINT *
                    endif
                    call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=1, userRC=status )
                    VERIFY_(STATUS)

                    if( gcm_internal_state%checkpointRequested ) then
                        call MAPL_GetObjectFromGC ( GCS(AGCM), MAPL_AGCM, RC=STATUS)
                        VERIFY_(STATUS)

                        call MAPL_GetResource( MAPL_AGCM, FileName, "MKIAU_CHECKPOINT_FILE:", rc=status)
                        VERIFY_(STATUS)
                        call MAPL_GetResource( MAPL_AGCM, FileType, "MKIAU_CHECKPOINT_TYPE:", rc=status)
                        VERIFY_(STATUS)

                        call MAPL_DateStampGet(clock, datestamp, rc=status)
                        VERIFY_(STATUS)
                        if( FileType == 'binary' ) RECORD_FNAME = trim(FileName) // '.' // DATESTAMP // '.bin'
                        if( FileType == 'pnc4'   ) RECORD_FNAME = trim(Filename) // '.' // DATESTAMP // '.nc4'

                        if( MAPL_AM_I_Root() ) PRINT *,TRIM(Iam)//":  Creating MKIAU Checkpoint ..."
                        call MAPL_CheckpointState( GIM(AGCM), CLOCK, RECORD_FNAME, Filetype, MAPL_AGCM, .FALSE., RC=STATUS )
                        if( MAPL_AM_I_Root() ) PRINT *
                        VERIFY_(STATUS)
                    endif
                endif

             if( gcm_internal_state%predur .ne. 0 ) then
             PREDICTOR_TIME_LOOP: do
                ! Run the ExtData Gridded Component
                ! ---------------------------------
                call ESMF_GridCompRun ( ExtData_internal_state%gc, importState=dummy, exportState=ExtData_internal_state%ExpState, clock=clock, userRC=status )
                VERIFY_(STATUS)

                ! Run the AGCM Gridded Component
                ! ------------------------------
                !    Ensure IAU tendencies are zero when running Predictor
                !    if( MAPL_AM_I_Root() ) print *, 'Zeroing   IAU forcing ...'
                     call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=2, userRC=status )
                     VERIFY_(STATUS)
                     call ESMF_GridCompRun ( GCS(AGCM), importState=GIM(AGCM), exportState=GEX(AGCM), clock=clock, userRC=status )
                     VERIFY_(STATUS)

                
                ! Advance the Clock
                ! -----------------

                call ESMF_ClockAdvance ( clock, rc=status )
                VERIFY_(STATUS)
                call MAPL_DateStampGet(clock, datestamp, rc=status)
                VERIFY_(STATUS)
                if( MAPL_AM_I_Root() ) PRINT *,TRIM(Iam)//":  Advancing AGCM  1-Model TimeStep ... ",datestamp

                call MAPL_GenericRunCouplers( MAPL, CHILD=AGCM, CLOCK=clock, RC=status )
                VERIFY_(STATUS)

                ! check for MKIAU
                ! ---------------
                timeForMKIAU = ESMF_AlarmIsRinging( GCM_INTERNAL_STATE%replayMKIAUAlarm, RC=STATUS )
                VERIFY_(STATUS)

                if( timeForMKIAU ) then
                    if( MAPL_AM_I_Root() ) then
                        PRINT *
                        PRINT *,TRIM(Iam)//":  Running MKIAU Within Predictor Loop, Creating IAU forcing ..."
                    endif
                    call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=1, userRC=status )
                    VERIFY_(STATUS)

                    if( gcm_internal_state%checkpointRequested ) then
                        call MAPL_GetObjectFromGC ( GCS(AGCM), MAPL_AGCM, RC=STATUS)
                        VERIFY_(STATUS)

                        call MAPL_GetResource( MAPL_AGCM, FileName, "MKIAU_CHECKPOINT_FILE:", rc=status)
                        VERIFY_(STATUS)
                        call MAPL_GetResource( MAPL_AGCM, FileType, "MKIAU_CHECKPOINT_TYPE:", rc=status)
                        VERIFY_(STATUS)

                        call MAPL_DateStampGet(clock, datestamp, rc=status)
                        VERIFY_(STATUS)
                        if( FileType == 'binary' ) RECORD_FNAME = trim(FileName) // '.' // DATESTAMP // '.bin'
                        if( FileType == 'pnc4'   ) RECORD_FNAME = trim(Filename) // '.' // DATESTAMP // '.nc4'

                        if( MAPL_AM_I_Root() ) PRINT *,TRIM(Iam)//":  Creating MKIAU Checkpoint ..."
                        call MAPL_CheckpointState( GIM(AGCM), CLOCK, RECORD_FNAME, Filetype, MAPL_AGCM, .FALSE., RC=STATUS )
                        VERIFY_(STATUS)
                        if( MAPL_AM_I_Root() ) PRINT *
                    endif
                endif

                call ESMF_VMBarrier(VM, rc=status)
                VERIFY_(STATUS)

                if (gcm_internal_state%run_history) then
                   call run_history(_RC)
                end if

                DONE = ESMF_AlarmIsRinging(GCM_INTERNAL_STATE%replayStopAlarm, RC=STATUS)
                VERIFY_(STATUS)
                if ( DONE ) exit

             enddo PREDICTOR_TIME_LOOP
             endif


             ! rewind the clock
             ! --------------------------------------------------------------
             call ESMF_ClockSet(clock, direction=ESMF_DIRECTION_REVERSE, rc=status)
             VERIFY_(STATUS)

             if( gcm_internal_state%predur .ne. 0 ) then
                 if( MAPL_AM_I_Root() ) print *
                 do
                   call ESMF_ClockAdvance ( clock, rc=status )
                   VERIFY_(STATUS)
                   call MAPL_DateStampGet(clock, datestamp, rc=status)
                   VERIFY_(STATUS)
                   if( MAPL_AM_I_Root() ) PRINT *,TRIM(Iam)//":  Rewinding Clock 1-Model TimeStep ... ",datestamp
                   call ESMF_ClockGet(clock, currTime=ct, rc=status)
                   VERIFY_(STATUS)
                   if (ct ==replayTime) exit
                 enddo
             endif
             if( MAPL_AM_I_Root() ) then
                print *
                print *, 'Continue  AGCM Replay ...'
                print *
             endif

             ! restore the state of the alarms
             ! -------------------------------
             call ESMF_ClockSet(clock, direction=ESMF_DIRECTION_FORWARD, rc=status)
             VERIFY_(STATUS)
             DO I = 1, nalarms
                call ESMF_AlarmSet(alarmList(I), ringTime=alarmRingTime(I), ringing=ringingState(I), rc=status)
                VERIFY_(STATUS)
             END DO

             deallocate(alarmList, alarmRingTime, ringingState)

             ! restore import and internal states
             ! ----------------------------------
             call MAPL_GenericRefresh( GCS(AGCM), GIM(AGCM), GEX(AGCM), clock, rc=status )
             VERIFY_(STATUS)

             call ESMF_GridCompRun ( ExtData_internal_state%gc, importState=dummy, &
                  exportState=ExtData_internal_state%ExpState, clock=clock, userRC=status )
             VERIFY_(STATUS)

          end if

          if (shutoffRpl ) then
             ! clear IAU tendencies
             ! --------------------
             ! if( MAPL_AM_I_Root() ) PRINT *,TRIM(Iam)//":  Zeroing IAU forcing ..."
             call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=2, userRC=status )
             VERIFY_(STATUS)
             call MAPL_GetObjectFromGC ( GCS(AGCM), MAPL_AGCM, RC=STATUS)
             VERIFY_(STATUS)
             call MAPL_DisableRecord(MAPL_AGCM,"startReplay",rc=status)
             VERIFY_(STATUS)
          end if

          ! Set Active PREDICTOR_STEP Alarm to OFF
          ! --------------------------------------
          call ESMF_ClockGetAlarm(clock, 'PredictorAlarm', alarm, rc=status)
          VERIFY_(STATUS)
          call ESMF_AlarmRingerOff(alarm, RC=STATUS)
          VERIFY_(STATUS)

          call ESMF_AlarmRingerOff(GCM_INTERNAL_STATE%replayStartAlarm, rc=status)
          VERIFY_(STATUS)

          call MAPL_TimerOff(MAPL,"--REPLAY"  )
       end if REPLAY

    ! Ensure Active PREDICTOR_STEP Alarm of OFF
    ! -----------------------------------------
    if(ESMF_AlarmIsCreated(PredictorIsActive)) then
       IF(ESMF_AlarmIsRinging(PredictorIsActive)) CALL ESMF_AlarmRingerOff(PredictorIsActive, RC=STATUS)
    end if
    VERIFY_(STATUS)

    ! the usual time step
    !--------------------

    call MAPL_TimerOn(MAPL,"--ATMOSPHERE"  )
    call MAPL_TimerOn(MAPL,"AGCM"        )
   
    call ESMF_GridCompRun ( GCS(AGCM), importState=GIM(AGCM), exportState=GEX(AGCM), clock=clock, userRC=status )
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"AGCM"        )
    call MAPL_TimerOff(MAPL,"--ATMOSPHERE"  )

     call MAPL_TimerOff(MAPL,"TOTAL")
     call MAPL_TimerOff(MAPL,"RUN"  )


     RETURN_(ESMF_SUCCESS)
   contains

     subroutine run_history(rc)
       integer, optional, intent(out) :: rc

       integer :: user_status,status
       type(ESMF_State), allocatable :: gcm_exports(:),hist_imports(:),hist_exports(:)
       type(ESMF_GridComp), allocatable :: hist_gcs(:)
       type(MAPL_MetaComp), pointer :: history_metaobj

       call MAPL_GetObjectFromGC ( gcm_internal_state%history_parent, history_metaobj, _RC)

       call MAPL_Get(history_metaobj, &
            childrens_export_states = hist_exports, &
            childrens_import_states = hist_imports, &
            childrens_gridcomps = hist_gcs, _RC)

       call ESMF_GridCompRun(hist_gcs(hist),importState=hist_imports(hist),&
                                       exportState=hist_exports(hist),&
                                       clock=clock,userRC=user_status,_RC)
       _RETURN(_SUCCESS)
 
     end subroutine run_history
 end subroutine Run

 subroutine write_modele_distribution_to_file(cf,rc)
    type(ESMF_Config), intent(inout) :: cf
    integer, optional, intent(out) :: rc

    integer :: status,npet,mypet,jm_world
    integer, allocatable :: jms(:)
    type(ESMF_VM) :: vm
    character(len=ESMF_MAXSTR) :: grid_type
    integer :: localAdjustment, excess, p

    call ESMF_VMGetCurrent(vm,_RC)
    call ESMF_VMGet(vm,localPet=mypet,peCount=npet,_RC)

    call ESMF_ConfigGetAttribute(cf,value=jm_world,label='AGCM.JM_WORLD:',_RC)
    call ESMF_ConfigGetAttribute(cf,value=grid_type,label='AGCM.GRID_TYPE:',_RC)
    _ASSERT(trim(grid_type) == 'LatLon',"model and GEOS Grid defs do not match, must be LatLon")
    allocate(jms(0:npet-1))
    select case (npet)
    case (1)
       jms(1) = JM_WORLD
       return
    case (2)
       jms(0) = JM_WORLD/2
       jms(1) = JM_WORLD - (JM_WORLD/2)
       return
    case (3:)
       ! 1st cut - round down
       jms(0:npet-1) = JM_WORLD/npet

       ! Fix at poles
       ! redistrute excess
       excess = JM_WORLD - sum(jms(0:npet-1))
       if (excess > 0) then
         jms(0) = max(2, jms(0))
       end if
       if (excess > 1) then
         jms(npet-1) = max(2, jms(npet-1))
       end if
       ! re-calculate excess again, since jms may have changed
       excess = JM_WORLD - sum(jms(0:npet-1))
       ! redistribute any remaining excess among interior processors
       do p = 1, npet - 2
          localAdjustment = (p+1)*excess/(npet-2) - (p*excess)/(npet-2)
          jms(p) = jms(p) + localAdjustment
       end do
    end select

    call MAPL_ConfigSetAttribute(cf,value=jms,label='AGCM.JMS:',_RC)

 end subroutine

end module GEOS_GcmGridCompMod

