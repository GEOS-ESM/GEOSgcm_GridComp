! $Id: GEOS_GcmGridComp.F90,v 1.34.30.1 2018/02/26 17:44:50 wputman Exp $

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GcmGridCompMod -- A Module to combine Agcm and Ogcm Gridded Components

! !INTERFACE:

module GEOS_GcmGridCompMod

! !USES:

   use ESMF
   use MAPL_Mod

   use GEOS_dataatmGridCompMod,  only:  DATAATM_SetServices => SetServices
   use ModelE_MAPL, only: AGCM_SetServices => SetServices
   use GEOS_mkiauGridCompMod,    only:  AIAU_SetServices => SetServices
   use DFI_GridCompMod,          only:  ADFI_SetServices => SetServices
   use GEOS_OgcmGridCompMod,     only:  OGCM_SetServices => SetServices
   use m_chars,                  only:  uppercase

   use modele_grid, only: MODELE_ATMGridCreate
   use dist_grid_mod, only: npes_world

  implicit none
  private

#ifdef USE_ODAS
    integer :: chdir, flag = 0
    external chdir 
#endif

! !PUBLIC MEMBER FUNCTIONS:

  public SetServices

  integer            :: NUM_ICE_CATEGORIES
  integer            :: NUM_ICE_LAYERS
  integer, parameter :: NUM_SNOW_LAYERS=1
  integer            :: DO_CICE_THERMO  
  integer            :: DO_DATAATM
  integer            :: DO_OBIO 

!=============================================================================

! !DESCRIPTION: This gridded component (GC) combines the Agcm GC, 
!   and Ogcm GC into a new composite Gcm GC.

 
!EOP

integer ::       AGCM
integer ::       OGCM
integer ::       AIAU
integer ::       ADFI

integer :: bypass_ogcm

type T_GCM_STATE
   private
   type (MAPL_LocStreamXFORM) :: XFORM_A2O
   type (MAPL_LocStreamXFORM) :: XFORM_O2A
   type(ESMF_State)           :: impSKIN ! Ocean thin layer
   type(ESMF_State)           :: expSKIN
   type(ESMF_Alarm)           :: alarmOcn
   type(ESMF_Alarm)           :: replayStartAlarm
   type(ESMF_Alarm)           :: replayStopAlarm
   type(ESMF_Alarm)           :: replayShutoffAlarm
   type(ESMF_TimeInterval)    :: replayDuration
   integer                    :: rpldur
   integer                    :: rplfreq
   integer                    :: rplbkgave
   logical                    :: rplRegular
   logical                    :: checkpointRequested = .false.
   character(len=ESMF_MAXSTR) :: checkpointFilename = ''
   character(len=ESMF_MAXSTR) :: checkpointFileType = ''
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
    integer, optional                  :: RC  ! return code

! !DESCRIPTION:  The SetServices for the PhysicsGcm GC needs to register its
!   Initialize and Run.  It uses the MAPL\_Generic construct for defining 
!   state specs and couplings among its children.  In addition, it creates the   
!   children GCs (AGCM and OGCM) and runs their
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
    character(len=ESMF_MAXSTR)          :: ReplayNudge
    logical                             :: rplRegular

    type (ESMF_Config)                  :: CF
    integer                             :: heartbeat

    integer                             :: ASSIMILATION_CYCLE
    integer                             :: PREDICTOR_DURATION
    integer                             :: CORRECTOR_DURATION
    integer                             :: REPLAY_FILE_FREQUENCY
    integer                             :: REPLAY_BKGAVE

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
#ifdef USE_ODAS
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_FINALIZE,  Finalize, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_WRITERESTART, Record,   RC=status)
    VERIFY_(STATUS)
#endif

! Get my MAPL_Generic state
!--------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

! Get constants from CF
! ---------------------

    call MAPL_GetResource ( MAPL,       DO_CICE_THERMO,     Label="USE_CICE_Thermo:" ,       DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

    if (DO_CICE_THERMO /= 0) then     ! Before merging CICEthermo with SaltWater, following were set as integer parameters
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_CATEGORIES, Label="CICE_N_ICE_CATEGORIES:" , RC=STATUS)
       VERIFY_(STATUS)
       call ESMF_ConfigGetAttribute(CF, NUM_ICE_LAYERS,     Label="CICE_N_ICE_LAYERS:" ,     RC=STATUS)
       VERIFY_(STATUS)
    else 
       NUM_ICE_CATEGORIES = 1
       NUM_ICE_LAYERS     = 1  
    endif

    call MAPL_GetResource ( MAPL, DO_OBIO,        Label="USE_OCEANOBIOGEOCHEM:",DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource ( MAPL, DO_DATAATM,     Label="USE_DATAATM:" ,     DEFAULT=0, RC=STATUS)
    VERIFY_(STATUS)

! Create childrens gridded components and invoke their SetServices
! ----------------------------------------------------------------

    if(DO_DATAATM/=0) then
       AGCM = MAPL_AddChild(GC, NAME='DATAATM', SS=DATAATM_SetServices, RC=STATUS)
       VERIFY_(STATUS)
    else
       AGCM = MAPL_AddChild(GC, NAME='AGCM', SS=Agcm_SetServices, RC=STATUS)
       VERIFY_(STATUS)
       AIAU = MAPL_AddChild(GC, NAME='AIAU', SS=AIAU_SetServices, RC=STATUS)
       VERIFY_(STATUS)
       ADFI = MAPL_AddChild(GC, NAME='ADFI', SS=ADFI_SetServices, RC=STATUS)
       VERIFY_(STATUS)
    endif
    ! OGCM = MAPL_AddChild(GC, NAME='OGCM', SS=Ogcm_SetServices, RC=STATUS)
    ! VERIFY_(STATUS)

! Get RUN Parameters (MERRA-2 Defaults) and Initialize for use in other Components (e.g., AGCM_GridComp and MKIAU_GridComp)
!--------------------------------------------------------------------------------------------------------------------------

  ! Get HEARTBEAT
  ! -------------
    call MAPL_GetResource( MAPL, heartbeat, Label="RUN_DT:", RC=STATUS)
    VERIFY_(STATUS)

  ! Get ASSIMILATION_CYCLE Duration.  If NOT_PRESENT, Initialize and Set Attribute.
  !--------------------------------------------------------------------------------
    call MAPL_GetResource( MAPL, ASSIMILATION_CYCLE, Label="ASSIMILATION_CYCLE:", default=-999, RC=STATUS)
    VERIFY_(STATUS)
    if( ASSIMILATION_CYCLE .eq. -999 ) then
        ASSIMILATION_CYCLE = 21600
        call MAPL_ConfigSetAttribute(  CF, ASSIMILATION_CYCLE, Label="ASSIMILATION_CYCLE:", RC=STATUS)
        VERIFY_(STATUS)
        IF(MAPL_AM_I_ROOT()) PRINT *,'Setting ASSIMILATION_CYCLE to ',ASSIMILATION_CYCLE
    endif

  ! Get CORRECTOR_DURATION.  If NOT_PRESENT, Initialize and Set Attribute.
  ! ----------------------------------------------------------------------
    call MAPL_GetResource( MAPL, CORRECTOR_DURATION, Label="CORRECTOR_DURATION:", default=-999, RC=STATUS)
    VERIFY_(STATUS)
    if( CORRECTOR_DURATION .eq. -999 ) then
        CORRECTOR_DURATION = ASSIMILATION_CYCLE
        call MAPL_ConfigSetAttribute(  CF, CORRECTOR_DURATION, Label="CORRECTOR_DURATION:", RC=STATUS)
        VERIFY_(STATUS)
        IF(MAPL_AM_I_ROOT()) PRINT *,'Setting CORRECTOR_DURATION, to ',CORRECTOR_DURATION
    endif

  ! Get PREDICTOR_DURATION.  If NOT_PRESENT, Initialize and Set Attribute.
  ! ----------------------------------------------------------------------
    call MAPL_GetResource( MAPL, PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", default=-999, RC=STATUS)
    VERIFY_(STATUS)
    VERIFY_(STATUS)
    if( PREDICTOR_DURATION .eq. -999 ) then
        PREDICTOR_DURATION = CORRECTOR_DURATION/2
        call MAPL_ConfigSetAttribute(  CF, PREDICTOR_DURATION, Label="PREDICTOR_DURATION:", RC=STATUS)
        VERIFY_(STATUS)
        IF(MAPL_AM_I_ROOT()) PRINT *,'Setting PREDICTOR_DURATION, to ',PREDICTOR_DURATION
    endif

  ! Get REPLAY_FILE_FREQUENCY.  If NOT_PRESENT, Initialize and Set Attribute.
  ! -------------------------------------------------------------------------
    call MAPL_GetResource(MAPL, REPLAY_FILE_FREQUENCY, Label="REPLAY_FILE_FREQUENCY:", default=-999, RC=STATUS)
    VERIFY_(STATUS)
    if( REPLAY_FILE_FREQUENCY .eq. -999 ) then
        REPLAY_FILE_FREQUENCY = CORRECTOR_DURATION
        call MAPL_ConfigSetAttribute(  CF, REPLAY_FILE_FREQUENCY, Label="REPLAY_FILE_FREQUENCY:", RC=STATUS)
        VERIFY_(STATUS)
        IF(MAPL_AM_I_ROOT()) PRINT *,'Setting REPLAY_FILE_FREQUENCY, to ',REPLAY_FILE_FREQUENCY
    endif


! Get REPLAY Parameters
! ---------------------
    call MAPL_GetResource(MAPL, ReplayMode,  'REPLAY_MODE:',  default="NoReplay", RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, ReplayNudge, 'REPLAY_NUDGE:', default="NO",       RC=STATUS )
    VERIFY_(STATUS)

    ReplayNudge = uppercase(ReplayNudge)
    ASSERT_( adjustl(ReplayNudge)=='YES' .or.  adjustl(ReplayNudge)=='NO' )

! -------------------------------------------
    ! We need to know if we are doing "regular" replay
    ! if yes, then we are modifing MKIAU import specs
    ! to perform averaging of all  MKIAU imports

    if(adjustl(ReplayMode)=="Regular") then
       rplRegular = .true.

       call MAPL_GetResource(MAPL, REPLAY_BKGAVE, Label="REPLAY_BKGAVE:", default=0, rc=status)
       VERIFY_(STATUS)

       if( adjustl(ReplayNudge) == "YES" ) then
           CORRECTOR_DURATION = heartbeat
           PREDICTOR_DURATION = 0
                REPLAY_BKGAVE = 0
       endif

       ! we need to add/change all of MKIAU import specs
       if (REPLAY_BKGAVE /= 0) then
          call MAPL_ChildAddAttribToImportSpec ( MAPL,   CHILD_ID = AIAU ,                           &
                                                 REFRESH_INTERVAL = CORRECTOR_DURATION ,             &
                                               AVERAGING_INTERVAL =      REPLAY_BKGAVE + heartbeat , &
                                                           OFFSET = PREDICTOR_DURATION +(REPLAY_BKGAVE/2), RC=status )
          VERIFY_(STATUS)
       end if
    else    
       rplRegular = .false.
          REPLAY_BKGAVE = 0
    endif

! Borrow exports from children
#ifdef USE_ODAS
    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'MOM_3D_MASK',                          &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'skin_T',                               &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'skin_S',                               &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'T',                                    &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'S',                                    &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'U',                                    &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'V',                                    &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'Z',                                    &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'RHO',                                  &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'SLV',                                  &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'SSH',                                  &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'PBO',                                  &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'TX',                                   &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'TY',                                   &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'MLD',                                  &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec ( GC   ,                          &
         SHORT_NAME = 'PSI',                                  &
         CHILD_ID   = OGCM,                                   &
                                                        RC=STATUS  )
    VERIFY_(STATUS)
 
    if (DO_CICE_THERMO /= 0) then
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'AICE',                                 &
            CHILD_ID   = OGCM,                                   &
            RC=STATUS  )
       VERIFY_(STATUS)
       
       call MAPL_AddExportSpec ( GC   ,                          &
            SHORT_NAME = 'HICE',                                 &
            CHILD_ID   = OGCM,                                   &
            RC=STATUS  )
       VERIFY_(STATUS)
    end if
 
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'CHLOROPHYLL',                          &
         CHILD_ID   = OGCM,                                   &

                                                        RC=STATUS  )
    VERIFY_(STATUS)
#else

! Export for IAU and/or Analysis purposes
! ---------------------------------------
    if(DO_DATAATM==0) then
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
    endif
#endif

    if(DO_DATAATM==0) then
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
            SRC_NAME  = 'Q', &
            DST_NAME  = 'QV', &
            DST_ID = AIAU,                                          &
            SRC_ID = AGCM,                                          &
            RC=STATUS  )
       VERIFY_(STATUS)

       call MAPL_AddConnectivity ( GC,                              &
            SHORT_NAME  = (/'U_DGRID','V_DGRID','PT     ',          &
            'PE     ','Q      ','OX     ' /),       &
            DST_ID = ADFI,                                          &
            SRC_ID = AGCM,                                          &
            RC=STATUS  )
       VERIFY_(STATUS)
    endif

! Next 17 vars are explicitly connected through exchange grid transforms Run
!---------------------------------------------------------------------------


     ! call MAPL_TerminateImport    ( GC,   &
     !      SHORT_NAME = (/'TAUXW  ','TAUYW  ','TAUXI  ','TAUYI  ',   &
     !                     'OUSTAR3','PS     ',                       &
     !                     'HW     ','TW     ','SW     ',             &
     !                     'HI     ','TI     ','SI     ' ,            &
     !                     'PENUVR ','PENUVF ','PENPAR ','PENPAF ',   &
     !                     'DISCHRG'/), &
     !      CHILD      = OGCM,           &
     !      RC=STATUS  )
     ! VERIFY_(STATUS)

     ! call MAPL_TerminateImport    ( GC,   &
     !      SHORT_NAME = (/'UU'/),                                    &
     !      CHILD      = OGCM,                                        &
     !      RC=STATUS  )
     ! VERIFY_(STATUS)

     ! if(DO_DATAATM==0) then
     !    call MAPL_TerminateImport    ( GC,                             &
     !         SHORT_NAME = (/'CO2SC'/),                                 &
     !         CHILD      = OGCM,                                        &
     !         RC=STATUS  )
     !    VERIFY_(STATUS)

     !    call MAPL_TerminateImport    ( GC,                             &
     !         SHORT_NAME = (/'DUDP', 'DUWT', 'DUSD',                    &
     !         'BCDP', 'BCWT',                            &
     !         'OCDP', 'OCWT' /),                         &
     !         CHILD      = OGCM,                                        &
     !         RC=STATUS  )
     !    VERIFY_(STATUS)

     !    call MAPL_TerminateImport    ( GC,                             &
     !         SHORT_NAME = (/'FSWBAND  ', 'FSWBANDNA'/),                &
     !         CHILD      = OGCM,                                        &
     !         RC=STATUS  )
     !    VERIFY_(STATUS)
     ! else
     !    call MAPL_TerminateImport    ( GC,   & 
     !         SHORT_NAME = (/'KPAR   ','UW     ','VW     ','UI     ', &
     !         'VI     ','TAUXBOT','TAUYBOT'/),         &
     !         CHILD      = AGCM,           &
     !         RC=STATUS  )
     !    VERIFY_(STATUS)
     ! end if

   !  if (DO_CICE_THERMO /= 0) then  
   !     call MAPL_TerminateImport    ( GC,   &
   !        SHORT_NAME = (/ &
   !                       'FRACICE', 'VOLICE ', 'VOLSNO ',              &
   !                       'ERGICE ', 'ERGSNO ', 'TAUAGE ', 'MPOND  '/),   &
   !        CHILD      = OGCM,           &
   !        RC=STATUS  )
   !     VERIFY_(STATUS)
   ! endif 

! Allocate this instance of the internal state and put it in wrapper.
! -------------------------------------------------------------------

    allocate( gcm_internal_state, stat=status )
    VERIFY_(STATUS)
    wrap%ptr => gcm_internal_state

    gcm_internal_state%rplRegular = rplRegular
    gcm_internal_state%rpldur     = PREDICTOR_DURATION
    gcm_internal_state%rplfreq    = CORRECTOR_DURATION
    gcm_internal_state%rplbkgave  = REPLAY_BKGAVE

! If doing regular replay make state to "borrow" gc and import to ExtData
! -----------------------------------------------------------------------
    if (rplRegular) then
       allocate( extdata_internal_state, stat=status)
       VERIFY_(STATUS)
       ExtDataWrap%ptr => extdata_internal_state
       call ESMF_UserCompSetInternalState( GC, 'ExtData_state',ExtDataWrap,status)
       VERIFY_(STATUS) 
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
    type(ESMF_Grid)     :: ogrid
    type(MAPL_LocStream):: exchA
    type(MAPL_LocStream):: exchO
    type(MAPL_LocStream):: locstA

    type(ESMF_DELayout) :: layout

    type (MAPL_MetaComp),  pointer  :: MAPL => null()
    type (MAPL_MetaComp),  pointer  :: CMAPL => null()
    type (ESMF_GridComp),      pointer  :: GCS(:) => null()
    type (ESMF_State),         pointer  :: GIM(:) => null()
    type (ESMF_State),         pointer  :: GEX(:) => null()
    character(len=ESMF_MAXSTR), pointer :: GCNAMES(:) => null()
    character(len=ESMF_MAXSTR)          :: skinname
    character(len=ESMF_MAXSTR)          :: tilingfile
    character(len=ESMF_MAXSTR)          :: tmpStr
    character(len=ESMF_MAXSTR)          :: ReplayMode
    type (T_GCM_STATE), pointer         :: gcm_internal_state => null() 
    type (GCM_wrap)                     :: wrap
    
    type(ESMF_Calendar)                 :: cal
    type(ESMF_Time)                     :: rep_StartTime
    type(ESMF_Time)                     :: currTime
    type(ESMF_Time)                     :: ringTime
    type(ESMF_TimeInterval)             :: Frequency
    type(ESMF_TimeInterval)             :: Duration
    type(ESMF_TimeInterval)             :: Shutoff
    type(ESMF_TimeInterval)             :: halfRplBkgAve
    type(ESMF_Alarm)                    :: replayStartAlarm
    type(ESMF_Alarm)                    :: replayStopAlarm
    type(ESMF_Alarm)                    :: replayShutoffAlarm
! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
    TYPE(ESMF_Alarm)                    :: PredictorIsActive
! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
    type(ESMF_Alarm)                    :: alarms(1)
    type(ESMF_DistGrid)                 :: distGRID
    integer                             :: PREDICTOR_DURATION
    integer                             :: CORRECTOR_DURATION
    integer                             :: REPLAY_BKGAVE
    integer                             :: rplshut
    integer                             :: rep_startdate(2)
    integer                             :: rep_YY
    integer                             :: rep_MM
    integer                             :: rep_DD
    integer                             :: rep_H
    integer                             :: rep_M
    integer                             :: rep_S

    integer :: nx, ny, nz

!=============================================================================

! Begin... 

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
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

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX, GCNAMES=GCNAMES,&
                                RC=STATUS )
    VERIFY_(STATUS)


! Create Atmospheric grid
!------------------------

    call MAPL_GetResource(mapl, nx, label = "AGCM_IM:", rc = rc)
    VERIFY_(rc)
    
    call MAPL_GetResource(mapl, ny, label = "AGCM_JM:", rc = rc)
    VERIFY_(rc)

    call MAPL_GetResource(mapl, nz, label = "AGCM_LM:", rc = rc)
    VERIFY_(rc)
    
    agrid = MODELE_ATMGridCreate(-180._ESMF_KIND_R8, & 
         -90._ESMF_KIND_R8, & 
         180._ESMF_KIND_R8, 90._ESMF_KIND_R8, nx, ny, npes_world, & 
         half_polar_cell=.true., rc=rc)
    VERIFY_(rc)
    
    call ESMF_AttributeSet(agrid, name='GRID_LM', value=nz, rc=rc)
    VERIFY_(rc)

! Create Ocean grid
!------------------
    ! call MAPL_GridCreate(GCS(OGCM), rc=status)
    ! VERIFY_(STATUS)
    ! call ESMF_GridCompGet(GCS(OGCM),  grid=ogrid, rc=status)
    ! VERIFY_(STATUS)

    ! call ESMF_GridCompGet(GCS(AGCM),  grid=agrid, rc=status)
    ! VERIFY_(STATUS)
    ! call ESMF_GridGet(agrid, DistGrid=distgrid, rc=status)
    ! VERIFY_(STATUS)
    ! call ESMF_DistGridGet(distGRID, deLayout=layout, RC=STATUS)
    ! VERIFY_(STATUS)

! Create exchange grids from tile file
!-------------------------------------

    ! call MAPL_GetResource(MAPL, TILINGFILE, 'TILING_FILE:', &
    !      default="tile.data", RC=STATUS)
    ! VERIFY_(STATUS)

    ! call MAPL_LocStreamCreate(exchA, LAYOUT=layout, FILENAME=TILINGFILE, &
    !                           NAME='MAIN_Atm',                           &
    !                           grid=agrid, RC=STATUS)
    ! VERIFY_(STATUS)

    ! call MAPL_LocStreamCreate(exchO, LAYOUT=layout, FILENAME=TILINGFILE, &
    !                           NAME='MAIN_Ocn',  mask=(/ MAPL_OCEAN /),   &
    !                           grid=ogrid, RC=STATUS)
    ! VERIFY_(STATUS)

! ------------------------------------------------------
! Add default exchange grid to both Atm and Ocn
! ------------------------------------------------------

    ! call MAPL_ExchangeGridSet(GCS(AGCM), exchA, rc=status)
    ! VERIFY_(STATUS)
    ! call MAPL_ExchangeGridSet(GCS(OGCM), exchO, rc=status)
    ! VERIFY_(STATUS)

! Recursive setup of grids (should be disabled)
    call ESMF_GridCompSet(GCS(AGCM),  grid=agrid, rc=status)
    VERIFY_(STATUS)
    
    ! call ESMF_GridCompSet(GCS(OGCM),  grid=ogrid, rc=status)
    ! VERIFY_(STATUS)
    if(DO_DATAATM==0) then
       call ESMF_GridCompSet(GCS(AIAU),  grid=agrid, rc=status)
       VERIFY_(STATUS)
       call ESMF_GridCompSet(GCS(ADFI),  grid=agrid, rc=status)
       VERIFY_(STATUS)
    endif

!ALT we need a grid for GCM - we put either the Agrid or the Ogrid
! depending of what exports are we propagating up
#ifdef USE_ODAS
    call ESMF_GridCompSet(GC, grid=ogrid, rc=status)
    VERIFY_(STATUS)
#else
    call ESMF_GridCompSet(GC, grid=agrid, rc=status)
    VERIFY_(STATUS)
#endif

    GCM_INTERNAL_STATE%checkpointRequested = .false.

    if(gcm_internal_state%rplRegular) then

       call MAPL_GetResource(MAPL, ReplayMode, 'REPLAY_MODE:', default="NoReplay", RC=STATUS )
       VERIFY_(STATUS)

! ALT replay: create alarms
       PREDICTOR_DURATION = gcm_internal_state%rplDur
       CORRECTOR_DURATION = gcm_internal_state%rplFreq
          REPLAY_BKGAVE = gcm_internal_state%rplbkgave
       call MAPL_GetResource(MAPL, rplshut, Label="REPLAY_SHUTOFF:", default=-3600, rc=status)
       VERIFY_(STATUS)
       call ESMF_ClockGet(clock, currTime=currTime, calendar=cal, rc=status)
       VERIFY_(STATUS)

       call ESMF_TimeIntervalSet(Duration,  S=PREDICTOR_DURATION, rc=status)
       VERIFY_(STATUS)
       call ESMF_TimeIntervalSet(Frequency, S=CORRECTOR_DURATION, rc=status)
       VERIFY_(STATUS)
       call ESMF_TimeIntervalSet(Shutoff, S=abs(rplshut), rc=status)
       VERIFY_(STATUS)
       call ESMF_TimeIntervalSet(halfRplBkgAve, S=REPLAY_BKGAVE/2, rc=status)
       VERIFY_(STATUS)

       ! Offset REPLAY_STARTTIME by REPLAY_FREQUENCY when running PREDICTOR_DURATION=0 with AGCM_IMPORT
       ! --------------------------------------------------------------------------------------------
       rep_StartTime = currTime
       if( PREDICTOR_DURATION == 0 ) then
           call MAPL_GetResource( MAPL, tmpStr, LABEL="AGCM_IMPORT_RESTART_FILE:", RC=STATUS)
           if (STATUS == ESMF_SUCCESS) rep_StartTime = currTime + Frequency
       endif

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

       RingTime = rep_StartTime !LLT: this depends on type of REPLAY being used

       replayStartAlarm = ESMF_AlarmCreate( name="startReplay", clock=clock, &
            RingTime=ringTime, RingInterval=Frequency, sticky=.true., rc=status )
       VERIFY_(STATUS)
       if(ringTime == currTime) then
          call ESMF_AlarmRingerOn(replayStartAlarm, rc=status)
          VERIFY_(STATUS)
       end if

       replayStopAlarm = ESMF_AlarmCreate(name="stopReplay", clock=clock, &
            RingTime=ringTime+Duration+halfRplBkgAve, RingInterval=Frequency, &
            sticky=.false., rc=status )
       VERIFY_(STATUS)

       replayShutoffAlarm = ESMF_AlarmCreate(name='ReplayShutOff', &
            clock=clock, ringInterval=Shutoff, sticky=.false., RC=STATUS )
       VERIFY_(STATUS)
       if (rplshut <= 0) then ! this is a "flag" to never use Shutoff alarm
          call ESMF_AlarmDisable(replayShutoffAlarm, RC=STATUS)
          VERIFY_(STATUS)
       end if

! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
       PredictorIsActive = ESMF_AlarmCreate(NAME="PredictorActive", CLOCK=clock, RingTime=currTime, STICKY=.TRUE., RC=STATUS)
       VERIFY_(STATUS)
       CALL ESMF_AlarmRingerOff(PredictorIsActive, RC=STATUS)
       VERIFY_(STATUS)
!      IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Predictor Alarm is ringing? ",ESMF_AlarmIsRinging(PredictorIsActive)
! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
       
       call MAPL_GetObjectFromGC ( GCS(AGCM), CMAPL, RC=STATUS)
       VERIFY_(STATUS)
       Alarms(1) = replayStartAlarm
       call MAPL_AddRecord(CMAPL, ALARMS, (/MAPL_Write2Ram/), rc=status)
       VERIFY_(STATUS)

       GCM_INTERNAL_STATE%replayStartAlarm = replayStartAlarm
       GCM_INTERNAL_STATE%replayStopAlarm = replayStopAlarm
       GCM_INTERNAL_STATE%replayShutoffAlarm = replayShutoffAlarm
       GCM_INTERNAL_STATE%replayDuration = Duration

       call MAPL_GetResource(MAPL, tmpStr, "REPLAY_CHECKPOINT_FILE:", &
            default="", rc=status)
       VERIFY_(STATUS)
       if (tmpStr /= "") then
          GCM_INTERNAL_STATE%checkpointRequested = .true.
          GCM_INTERNAL_STATE%checkpointFilename = tmpStr
          call MAPL_GetResource(MAPL, tmpStr, "REPLAY_CHECKPOINT_TYPE:", &
               default="binary", rc=status)
          VERIFY_(STATUS)
          GCM_INTERNAL_STATE%checkpointFileType = tmpStr
       endif

    end if

! **********************************************************************
! ****                 Initialize Gridded Components                ****
! **********************************************************************

    call MAPL_TimerOff(MAPL,"TOTAL")
#ifdef USE_ODAS
   if(COMP_NAME(1:3) == "ens") then
      status = chdir("restarts." // comp_name(5:7));
      VERIFY_(STATUS) 
   end if
#endif
    call MAPL_GenericInitialize ( GC, import, export, clock, rc=status )
    VERIFY_(STATUS)
#ifdef USE_ODAS
   if(COMP_NAME(1:3) == "ens") then 
      status = chdir("..");
      VERIFY_(STATUS)
   end if
#endif

    call MAPL_TimerOn(MAPL,"TOTAL")

! Create XFORMs
!--------------

    if(DO_DATAATM/=0) then
       skinname = 'DATAATM'
    else
       skinname = 'SALTWATER'
    endif

   ! call MAPL_GetResource(MAPL, bypass_ogcm, "BYPASS_OGCM:", &
   !      default=0, rc=status)
   ! VERIFY_(STATUS)

   ! if (bypass_ogcm == 0) then
   ! call MAPL_GetChildLocstream(GCS(AGCM), locstA, skinname, rc=STATUS)
   ! VERIFY_(STATUS)

   ! call MAPL_LocStreamCreateXform ( XFORM=GCM_INTERNAL_STATE%XFORM_A2O, &
   !                                  LocStreamOut=exchO, &
   !                                  LocStreamIn=locstA, &
   !                                  NAME='XFORM_A2O', &
   !                                  RC=STATUS )
   ! VERIFY_(STATUS)

   ! call MAPL_LocStreamCreateXform ( XFORM=GCM_INTERNAL_STATE%XFORM_O2A, &
   !                                  LocStreamOut=locstA, &
   !                                  LocStreamIn=exchO, &
   !                                  NAME='XFORM_O2A', &
   !                                  RC=STATUS )
   ! VERIFY_(STATUS)
   ! end if

! This part has some explicit hierarchy build in...
!--------------------------------------------------------------

   ! call MAPL_ExportStateGet ( (/ GEX(AGCM) /), skinname, &
   !                             GCM_INTERNAL_STATE%expSKIN, rc=status )
   ! VERIFY_(STATUS)

   ! call MAPL_ImportStateGet ( GCS(AGCM) , GIM(AGCM), skinname,  &
   !                            GCM_INTERNAL_STATE%impSKIN, rc=status )
   ! VERIFY_(STATUS)

   ! call AllocateExports(GCM_INTERNAL_STATE%expSKIN, &
   !      (/'TAUXO    ', 'TAUYO    ','TAUXI    ', 'TAUYI    ', &
   !      'PENPAR   ', 'PENPAF   ','PENUVR   ', 'PENUVF   ', &
   !      'OUSTAR3  ', 'PS       '/), &
   !      RC=STATUS)
   ! VERIFY_(STATUS)

   ! if (DO_CICE_THERMO /= 0) then  
   !    call AllocateExports(GCM_INTERNAL_STATE%expSKIN, &
   !         (/'TAUXW    ', 'TAUYW    ', &
   !         'DISCHARGE'/), &
   !         RC=STATUS)
   !    VERIFY_(STATUS)
   ! end if

   ! if(DO_OBIO/=0) then
   !    call AllocateExports( GCM_INTERNAL_STATE%expSKIN,                   &
   !         (/'UU'/),                                     &
   !         RC=STATUS )
   !    VERIFY_(STATUS)
   !    if(DO_DATAATM==0) then
   !       call AllocateExports( GCM_INTERNAL_STATE%expSKIN,                   &
   !            (/'CO2SC'/),                                  &
   !            RC=STATUS )
   !       VERIFY_(STATUS)
         
   !       call AllocateExports_UGD( GCM_INTERNAL_STATE%expSKIN,               &
   !            (/'DUDP', 'DUWT', 'DUSD',                 &
   !            'BCDP', 'BCWT',                         &
   !            'OCDP', 'OCWT' /),                      &
   !            RC=STATUS )
   !       VERIFY_(STATUS)
         
   !       call AllocateExports_UGD( GCM_INTERNAL_STATE%expSKIN,               &
   !            (/'FSWBAND  ', 'FSWBANDNA'/),               &
   !            RC=STATUS )
   !       VERIFY_(STATUS)
   !    endif
   ! endif

   ! call AllocateExports(GEX(OGCM), (/'UW      ', 'VW      ', &
   !                                   'UI      ', 'VI      ', & 
   !                                   'KPAR    ', 'TS_FOUND' /), RC=STATUS)
   ! VERIFY_(STATUS)
   ! if (DO_CICE_THERMO == 0) then  
   !    call AllocateExports(GEX(OGCM), (/'FRACICE '/), RC=STATUS)
   !    VERIFY_(STATUS)
   ! else
   !    call AllocateExports(GEX(OGCM), (/'TAUXIBOT', 'TAUYIBOT'/), RC=STATUS)
   ! end if
    
   ! call ESMF_ClockGetAlarm(clock, alarmname=trim(GCNAMES(OGCM)) // '_Alarm', &
   !                         alarm=GCM_INTERNAL_STATE%alarmOcn, rc=status)
   ! VERIFY_(STATUS)

#ifdef PRINT_STATES
    call WRITE_PARALLEL ( trim(Iam)//": IMPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( IMPORT, rc=STATUS )
    call WRITE_PARALLEL ( trim(Iam)//": EXPORT State" )
    if ( MAPL_am_I_root() ) call ESMF_StatePrint ( EXPORT, rc=STATUS )
#endif

    call MAPL_TimerOff(MAPL,"TOTAL")
    call MAPL_TimerOff(MAPL,"INITIALIZE")

   RETURN_(ESMF_SUCCESS)

 contains
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

    type (MAPL_LocStreamXFORM)           :: XFORM_A2O
    type (MAPL_LocStreamXFORM)           :: XFORM_O2A
#ifdef USE_ODAS
    type(MAPL_LocStream)                 :: LOCSTREAM_O, LOCSTREAM_A
#endif
    type(ESMF_State)           :: impSKIN
    type(ESMF_State)           :: expSKIN

    type (T_GCM_STATE), pointer         :: gcm_internal_state => null() 
    type (GCM_wrap)                     :: wrap
    type (T_ExtData_STATE), pointer     :: ExtData_internal_state  => null()
    type (ExtData_wrap)                 :: ExtDatawrap
    type (ESMF_State)                   :: Dummy

    type(ESMF_Alarm)              :: alarm
    type(ESMF_Alarm), allocatable :: AlarmList(:)
    type(ESMF_Time),  allocatable :: AlarmRingTime(:)
    type(ESMF_Time)               :: ct, replayTime
    logical,          allocatable :: ringingState(:)
    integer                       :: i, nalarms
    logical                       :: done
    logical                       :: shutoffRpl
    logical                       :: timeForRpl
    character(len=ESMF_MAXSTR)    :: ReplayMode
    character(len=ESMF_MAXSTR)    :: record_fname
    character(len=14)             :: DATESTAMP !YYYYMMDD_HHMMz

    integer, dimension(NUM_ICE_CATEGORIES) :: SUBINDEXO, SUBINDEXA
    integer                                :: N

! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
    TYPE(ESMF_Alarm)                    :: PredictorIsActive
! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
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


    call MAPL_TimerON(MAPL,"TOTAL")
    call MAPL_TimerON(MAPL,"RUN"  )


! Get my internal private state. This contains the transforms
!  between the exchange grid and the atmos grid.
!-------------------------------------------------------------

    call ESMF_UserCompGetInternalState(gc, 'GCM_state', wrap, status)
    VERIFY_(STATUS)
    gcm_internal_state => wrap%ptr

    XFORM_A2O = GCM_INTERNAL_STATE%XFORM_A2O
    XFORM_O2A = GCM_INTERNAL_STATE%XFORM_O2A
    impSKIN   = GCM_INTERNAL_STATE%impSKIN
    expSKIN   = GCM_INTERNAL_STATE%expSKIN

    do N=1,NUM_ICE_CATEGORIES
       SUBINDEXO(N) = N
       SUBINDEXA(N) = N+1
    enddo

! Get children and their im/ex states from my generic state.
!----------------------------------------------------------

    call MAPL_Get ( MAPL, GCS=GCS, GIM=GIM, GEX=GEX,&
                                RC=STATUS )
    VERIFY_(STATUS)
#ifdef USE_ODAS
    call MAPL_ExchangeGridGet(GCS(OGCM), EXCH=LOCSTREAM_O, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_ExchangeGridGet(GCS(AGCM), EXCH=LOCSTREAM_A, RC=STATUS)
    VERIFY_(STATUS)
#endif
    if(DO_DATAATM==0) then
       ! replay-mode stuff
       ! -----------------
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
             if( MAPL_AM_I_Root() ) print *, 'Zeroing   IAU forcing ...'
             call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=2, userRC=status )
             VERIFY_(STATUS)

             ! we should save import and internal states
             ! this is already done by record in MAPL_CAP

             call ESMF_ClockGet(clock, currTime=ct, rc=status)
             VERIFY_(STATUS)

! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
          CALL ESMF_ClockGetAlarm(CLOCK, "PredictorActive", PredictorIsActive, RC=STATUS)
          VERIFY_(STATUS)
          CALL ESMF_AlarmRingerOn(PredictorIsActive, RC=STATUS)
          VERIFY_(STATUS)
!         IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Start of Predictor ringing is ",ESMF_AlarmIsRinging(PredictorIsActive)
! ---------------- Alarm for detecting active PREDICTOR step ---------------- !

             replayTime = ct

             ! save alarms' states
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

             ! fix the predictor alarm
             call ESMF_ClockGetAlarm(clock, 'PredictorAlarm', alarm, rc=status)
             VERIFY_(STATUS)
             call ESMF_AlarmRingerOn(alarm, RC=STATUS)
             VERIFY_(STATUS)


           ! time loop
             if( gcm_internal_state%rplDur.ne.0 ) then
             TIME_LOOP: do

                ! Run the ExtData Gridded Component
                ! ---------------------------------
                call ESMF_GridCompRun ( ExtData_internal_state%gc, importState=dummy, exportState=ExtData_internal_state%ExpState, clock=clock, userRC=status )
                VERIFY_(STATUS)

                ! Run the AGCM Gridded Component
                ! ------------------------------
                call ESMF_GridCompRun ( GCS(AGCM), importState=GIM(AGCM), exportState=GEX(AGCM), clock=clock, userRC=status )
                VERIFY_(STATUS)

                ! call RUN_OCEAN(RC=STATUS)
                ! VERIFY_(STATUS)

                ! Advance the Clock
                ! -----------------

                call ESMF_ClockAdvance ( clock, rc=status )
                VERIFY_(STATUS)
                call MAPL_DateStampGet(clock, datestamp, rc=status)
                VERIFY_(STATUS)
                if( MAPL_AM_I_Root() ) print *, 'Advancing AGCM  1-Model TimeStep ... ',datestamp

                call MAPL_GenericRunCouplers( MAPL, CHILD=AGCM, CLOCK=clock, RC=status )
                VERIFY_(STATUS)

                call ESMF_VMBarrier(VM, rc=status)
                VERIFY_(STATUS)

                DONE = ESMF_AlarmIsRinging(GCM_INTERNAL_STATE%replayStopAlarm, RC=STATUS)
                VERIFY_(STATUS)
                if ( DONE ) exit

             enddo TIME_LOOP ! end of time loop
             endif

             ! set clock to mid-point on bkg-ave
             call ESMF_ClockSet(clock, currTime=replayTime+GCM_INTERNAL_STATE%replayDuration, rc=status)
             VERIFY_(STATUS)

             ! call make Atm IAU
             if( MAPL_AM_I_Root() ) print *, 'Creating  IAU forcing ...'
             call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=1, userRC=status )
             VERIFY_(STATUS)

             ! rewind the clock
             ! --------------------------------------------------------------
             call ESMF_ClockSet(clock, direction=ESMF_DIRECTION_REVERSE, rc=status)
             VERIFY_(STATUS)

             if( gcm_internal_state%rplDur.ne.0 ) then
                 do
                   call ESMF_ClockAdvance ( clock, rc=status )
                   VERIFY_(STATUS)
                   call MAPL_DateStampGet(clock, datestamp, rc=status)
                   VERIFY_(STATUS)
                   if( MAPL_AM_I_Root() ) print *, 'Rewinding Clock 1-Model TimeStep ... ',datestamp
                   call ESMF_ClockGet(clock, currTime=ct, rc=status)
                   VERIFY_(STATUS)
                   if (ct ==replayTime) exit
                 enddo
             endif
             if( MAPL_AM_I_Root() ) then
                print *, 'Continue  AGCM Replay ...'
                print *
             endif

             ! restore the state of the alarms
             call ESMF_ClockSet(clock, direction=ESMF_DIRECTION_FORWARD, rc=status)
             VERIFY_(STATUS)
             DO I = 1, nalarms
                call ESMF_AlarmSet(alarmList(I), ringTime=alarmRingTime(I), ringing=ringingState(I), rc=status)
                VERIFY_(STATUS)
             END DO

             deallocate(alarmList, alarmRingTime, ringingState)

             ! restore import and internal states
             call MAPL_GenericRefresh( GCS(AGCM), GIM(AGCM), GEX(AGCM), clock, rc=status )
             VERIFY_(STATUS)
             call MAPL_GenericRefresh( GCS(OGCM), GIM(OGCM), GEX(OGCM), clock, rc=status )
             VERIFY_(STATUS)

             call ESMF_GridCompRun ( ExtData_internal_state%gc, importState=dummy, &
                  exportState=ExtData_internal_state%ExpState, clock=clock, userRC=status )
             VERIFY_(STATUS)

          end if
          if (shutoffRpl ) then

             ! clear IAU tendencies
             if( MAPL_AM_I_Root() ) print *, 'Zeroing   IAU forcing ...'
             call ESMF_GridCompRun ( GCS(AIAU), importState=GIM(AIAU), exportState=GEX(AIAU), clock=clock, phase=2, userRC=status )
             VERIFY_(STATUS)
             call MAPL_GetObjectFromGC ( GCS(AGCM), MAPL_AGCM, RC=STATUS)
             VERIFY_(STATUS)
             call MAPL_DisableRecord(MAPL_AGCM,"startReplay",rc=status)
             VERIFY_(STATUS)

          end if

          if (timeForRpl .and. gcm_internal_state%checkpointRequested) then
             !ALT: this is done here, 
             !     so that the tendensies are already either computed or cleared,
             !     and this guaranties identical results after "EXACT" replay(s).
             !     Also the clock had been rewound, and time stamping is correct.

             ! add timestamp to filename
             call MAPL_DateStampGet(clock, datestamp, rc=status)
             VERIFY_(STATUS)

             RECORD_FNAME = trim(gcm_internal_state%checkpointFilename) // '.' // DATESTAMP // '.bin'

             call MAPL_GetObjectFromGC ( GCS(AGCM), MAPL_AGCM, RC=STATUS)
             VERIFY_(STATUS)

             call MAPL_CheckpointState(GIM(AGCM), CLOCK, RECORD_FNAME, &
                  gcm_internal_state%checkpointFiletype, &
                  MAPL_AGCM, .FALSE., RC=STATUS)
             VERIFY_(STATUS)
          end if

          ! fix the predictor alarm
          call ESMF_ClockGetAlarm(clock, 'PredictorAlarm', alarm, rc=status)
          VERIFY_(STATUS)
          call ESMF_AlarmRingerOff(alarm, RC=STATUS)
          VERIFY_(STATUS)

          call ESMF_AlarmRingerOff(GCM_INTERNAL_STATE%replayStartAlarm, rc=status)
          VERIFY_(STATUS)

          call MAPL_TimerOff(MAPL,"--REPLAY"  )
       end if REPLAY
    endif

! ---------------- Alarm for detecting active PREDICTOR step ---------------- !
          IF(ESMF_AlarmIsRinging(PredictorIsActive)) CALL ESMF_AlarmRingerOff(PredictorIsActive, RC=STATUS)
          VERIFY_(STATUS)
!         IF(MAPL_AM_I_ROOT()) PRINT *,TRIM(Iam)//": Start of Predictor ringing is ",ESMF_AlarmIsRinging(PredictorIsActive)
! ---------------- Alarm for detecting active PREDICTOR step ---------------- !

    ! the usual time step
    !--------------------


    call MAPL_TimerOn(MAPL,"--ATMOSPHERE"  )
    call MAPL_TimerOn(MAPL,"AGCM"          )
   
    call ESMF_GridCompRun ( GCS(AGCM), importState=GIM(AGCM), exportState=GEX(AGCM), clock=clock, userRC=status )
    VERIFY_(STATUS)

    call MAPL_TimerOff(MAPL,"AGCM"          )
    call MAPL_TimerOff(MAPL,"--ATMOSPHERE"  )

    ! if(DO_DATAATM==0) then
    !    ! Accumulate for digital filter
    !    ! -----------------------------
    !    call ESMF_GridCompRun ( GCS(ADFI), importState=GIM(ADFI), exportState=GEX(ADFI), clock=clock, userRC=status )
    !    VERIFY_(STATUS)
    ! endif

    ! call RUN_OCEAN(RC=STATUS)
    ! VERIFY_(STATUS)
    

     call MAPL_TimerOff(MAPL,"RUN"  )
     call MAPL_TimerOff(MAPL,"TOTAL")


     RETURN_(ESMF_SUCCESS)
   contains
     subroutine RUN_OCEAN(rc)
       integer, optional, intent(OUT) :: rc
       integer :: status
       character(len=ESMF_MAXSTR) :: Iam='Run_Ocean'

    if ( bypass_ogcm /= 0) then
       RETURN_(ESMF_SUCCESS)
    end if

    if ( ESMF_AlarmIsRinging(GCM_INTERNAL_STATE%alarmOcn, RC=STATUS) ) then
       VERIFY_(STATUS)

!       time to couple and run ocean

! Synchronize for Next TimeStep
! -----------------------------

       call ESMF_VMBarrier(VM, rc=status)
       VERIFY_(STATUS)


       call MAPL_TimerOn(MAPL,"--A2O"  )

! get tilevars and transform them
! SURFACE exports to OGCM imports
! Example how to do TW, we need to do all of the OCGMimports

! Copy attributes to deal with friendliness
!------------------------------------------

       call MAPL_CopyFriendliness(GIM(OGCM),'TW',expSKIN,'TSKINW', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_CopyFriendliness(GIM(OGCM),'HW',expSKIN,'HSKINW', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_CopyFriendliness(GIM(OGCM),'SW',expSKIN,'SSKINW', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_CopyFriendliness(GIM(OGCM),'TI',expSKIN,'TSKINI' , RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_CopyFriendliness(GIM(OGCM),'HI',expSKIN,'HSKINI', RC=STATUS)
       VERIFY_(STATUS)
       call MAPL_CopyFriendliness(GIM(OGCM),'SI',expSKIN,'SSKINI', RC=STATUS)
       VERIFY_(STATUS)
       if (DO_CICE_THERMO /= 0) then  
          call MAPL_CopyFriendliness(GIM(OGCM),'FRACICE',expSKIN,'FR', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_CopyFriendliness(GIM(OGCM),'VOLICE',expSKIN,'VOLICE', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_CopyFriendliness(GIM(OGCM),'VOLSNO',expSKIN,'VOLSNO', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_CopyFriendliness(GIM(OGCM),'ERGICE',expSKIN,'ERGICE', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_CopyFriendliness(GIM(OGCM),'ERGSNO',expSKIN,'ERGSNO', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_CopyFriendliness(GIM(OGCM),'TAUAGE',expSKIN,'TAUAGE', RC=STATUS)
          VERIFY_(STATUS)
          call MAPL_CopyFriendliness(GIM(OGCM),'MPOND',expSKIN,'VOLPOND', RC=STATUS)
          VERIFY_(STATUS)
       endif 
       
! Do the routing between the atm and ocean's decompositions of the exchage grid
!------------------------------------------------------------------------------

       call DO_A2O(GIM(OGCM),'TW'     ,expSKIN,'TSKINW' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'HW'     ,expSKIN,'HSKINW' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'SW'     ,expSKIN,'SSKINW' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'HI'     ,expSKIN,'HSKINI' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'SI'     ,expSKIN,'SSKINI' , RC=STATUS)
       VERIFY_(STATUS)
       if (DO_CICE_THERMO == 0) then
          call DO_A2O(GIM(OGCM),'TI'     ,expSKIN,'TSKINI' , RC=STATUS)
          VERIFY_(STATUS)
       endif

       if (DO_CICE_THERMO /= 0) then  
          call DO_A2O_SUBTILES_R4R4(GIM(OGCM),'TI'     , SUBINDEXO, &
               expSKIN  ,'TSKINI' , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES_R4R8(GIM(OGCM),'FRACICE', SUBINDEXO, &
               expSKIN  ,'FR'     , SUBINDEXA, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES_R4R8(GIM(OGCM),'VOLICE' , SUBINDEXO, &
               expSKIN  ,'VOLICE' , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES_R4R8(GIM(OGCM),'VOLSNO' , SUBINDEXO, &
               expSKIN  ,'VOLSNO' , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES2D_R4R8(GIM(OGCM),'ERGICE' , SUBINDEXO, &
               expSKIN  ,'ERGICE' , SUBINDEXO, &
               NUM_ICE_LAYERS, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES2D_R4R8(GIM(OGCM),'ERGSNO' , SUBINDEXO, &
               expSKIN  ,'ERGSNO' , SUBINDEXO, &
               NUM_SNOW_LAYERS, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES_R4R4(GIM(OGCM),'TAUAGE' , SUBINDEXO, &
               expSKIN  ,'TAUAGE' , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_A2O_SUBTILES_R4R8(GIM(OGCM),'MPOND'   , SUBINDEXO, &
               expSKIN  ,'VOLPOND' , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
       endif

       if (DO_CICE_THERMO /= 0) then  
          call DO_A2O(GIM(OGCM),'TAUXW'  ,expSKIN,'TAUXW'  , RC=STATUS); VERIFY_(STATUS)
          call DO_A2O(GIM(OGCM),'TAUYW'  ,expSKIN,'TAUYW'  , RC=STATUS); VERIFY_(STATUS)
       else
          call DO_A2O(GIM(OGCM),'TAUXW'  ,expSKIN,'TAUXO'  , RC=STATUS); VERIFY_(STATUS)
          call DO_A2O(GIM(OGCM),'TAUYW'  ,expSKIN,'TAUYO'  , RC=STATUS); VERIFY_(STATUS)
       endif
       
       call DO_A2O(GIM(OGCM),'TAUXI'  ,expSKIN,'TAUXI'  , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'TAUYI'  ,expSKIN,'TAUYI'  , RC=STATUS)
       VERIFY_(STATUS)
       
       call DO_A2O(GIM(OGCM),'PENPAR' ,expSKIN,'PENPAR' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'PENPAF' ,expSKIN,'PENPAF' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'PENUVR' ,expSKIN,'PENUVR' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'PENUVF' ,expSKIN,'PENUVF' , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'OUSTAR3',expSKIN,'OUSTAR3', RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'PS'     ,expSKIN,'PS'     , RC=STATUS)
       VERIFY_(STATUS)
       call DO_A2O(GIM(OGCM),'DISCHRG',expSKIN,'DISCHARGE', RC=STATUS)
       VERIFY_(STATUS)
       if(DO_OBIO/=0) then
          call DO_A2O(GIM(OGCM),'UU'     ,expSKIN,'UU'     , RC=STATUS)
          VERIFY_(STATUS)
          if(DO_DATAATM==0) then
             call DO_A2O(GIM(OGCM),'CO2SC'  ,expSKIN,'CO2SC'  , RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'DUDP', expSKIN, 'DUDP', RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'DUWT', expSKIN, 'DUWT', RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'DUSD', expSKIN, 'DUSD', RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'BCDP', expSKIN, 'BCDP', RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'BCWT', expSKIN, 'BCWT', RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'OCDP', expSKIN, 'OCDP', RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'OCWT', expSKIN, 'OCWT', RC=STATUS)
             VERIFY_(STATUS)

             call DO_A2O_UGD(GIM(OGCM), 'FSWBAND',   expSKIN, 'FSWBAND',   RC=STATUS)
             VERIFY_(STATUS)
             call DO_A2O_UGD(GIM(OGCM), 'FSWBANDNA', expSKIN, 'FSWBANDNA', RC=STATUS)
             VERIFY_(STATUS)
          endif
       endif
       call MAPL_TimerOff(MAPL,"--A2O"  )
       
!--
! OGCM exports to SURFACE imports

! Example how to do UW, we need to do all of the SURFACE 'friendly' tilevars


       call MAPL_TimerOn(MAPL,"--OCEAN"  )
       call MAPL_TimerOn(MAPL,"OGCM"     )

       call ESMF_GridCompRun ( GCS(OGCM), importState=gim(OGCM), exportState=gex(OGCM), clock=clock, userRC=status )
       VERIFY_(STATUS)

       call MAPL_TimerOff(MAPL,"OGCM"     )
       call MAPL_TimerOff(MAPL,"--OCEAN"  )


       call MAPL_TimerOn (MAPL,"--O2A"  )
#ifdef USE_ODAS
       call MAPL_GetResource(MAPL, flag, Label = "surface_winds_to_iOdas:", default = 0, rc = status)
       VERIFY_(STATUS)
       if(flag > 0) then
           call DO_A2OG(GEX(OGCM), 'TX', expSKIN, 'U10M', RC=STATUS)
           VERIFY_(STATUS)
           call DO_A2OG(GEX(OGCM), 'TY', expSKIN, 'V10M', RC=STATUS)
           VERIFY_(STATUS)
       endif  
#endif
       call DO_O2A(expSKIN, 'TSKINW'   , GIM(OGCM), 'TW'    , RC=STATUS)
       VERIFY_(STATUS)
       if (DO_CICE_THERMO == 0) then
         call DO_O2A(expSKIN, 'TSKINI'   , GIM(OGCM), 'TI'    , RC=STATUS)
         VERIFY_(STATUS)
       endif
       call DO_O2A(expSKIN, 'HSKINW'   , GIM(OGCM), 'HW'    , RC=STATUS)
       VERIFY_(STATUS)
       call DO_O2A(expSKIN, 'HSKINI'   , GIM(OGCM), 'HI'    , RC=STATUS)
       VERIFY_(STATUS)
       call DO_O2A(expSKIN, 'SSKINW'   , GIM(OGCM), 'SW'    , RC=STATUS)
       VERIFY_(STATUS)
       call DO_O2A(expSKIN, 'SSKINI'   , GIM(OGCM), 'SI'    , RC=STATUS)
       VERIFY_(STATUS)

       if (DO_CICE_THERMO /= 0) then  
          call DO_O2A_SUBTILES_R4R4(expSKIN  , 'TSKINI'     , SUBINDEXO,  &
               GIM(OGCM), 'TI'         , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES_R8R4(expSKIN  , 'FR'         , SUBINDEXA,  &
               GIM(OGCM), 'FRACICE'    , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES_R8R4(expSKIN  , 'VOLICE'     , SUBINDEXO,  &
               GIM(OGCM), 'VOLICE'     , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES_R8R4(expSKIN  , 'VOLSNO'     , SUBINDEXO,  &
               GIM(OGCM), 'VOLSNO'     , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES2D_R8R4(expSKIN  , 'ERGICE'     , SUBINDEXO,  &
               GIM(OGCM), 'ERGICE'     , SUBINDEXO,  &
               NUM_ICE_LAYERS, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES2D_R8R4(expSKIN  , 'ERGSNO'     , SUBINDEXO,  &
               GIM(OGCM), 'ERGSNO'     , SUBINDEXO,  &
               NUM_SNOW_LAYERS, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES_R4R4(expSKIN  , 'TAUAGE'      , SUBINDEXO,  &
               GIM(OGCM), 'TAUAGE'      , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A_SUBTILES_R8R4(expSKIN  , 'VOLPOND'     , SUBINDEXO,  &
               GIM(OGCM), 'MPOND'       , SUBINDEXO, RC=STATUS)
          VERIFY_(STATUS)
       endif  

       call DO_O2A(impSKIN, 'UW'       , GEX(OGCM), 'UW'    , RC=STATUS)
       VERIFY_(STATUS)
       call DO_O2A(impSKIN, 'VW'       , GEX(OGCM), 'VW'    , RC=STATUS)
       VERIFY_(STATUS)
       call DO_O2A(impSKIN, 'KPAR'     , GEX(OGCM), 'KPAR'  , RC=STATUS)
       VERIFY_(STATUS)

       if (DO_CICE_THERMO == 0) then
          call DO_O2A(impSKIN, 'FRACICE'  , GEX(OGCM), 'FRACICE', RC=STATUS)
          VERIFY_(STATUS)
       else
          call DO_O2A(impSKIN, 'TAUXBOT'  , GEX(OGCM), 'TAUXIBOT', RC=STATUS)
          VERIFY_(STATUS)
          call DO_O2A(impSKIN, 'TAUYBOT'  , GEX(OGCM), 'TAUYIBOT', RC=STATUS)
          VERIFY_(STATUS)
       end if

       call DO_O2A(impSKIN, 'UI'       , GEX(OGCM), 'UI'    , RC=STATUS)
       VERIFY_(STATUS)
       call DO_O2A(impSKIN, 'VI'       , GEX(OGCM), 'VI'    , RC=STATUS)
       VERIFY_(STATUS)

! OGCM export of TS_FOUND to SKIN
!---------------------------------
        call DO_O2A(impSKIN, 'TS_FOUND' , GEX(OGCM), 'TS_FOUND' , RC=STATUS)
        VERIFY_(STATUS)

        call ESMF_AlarmRingerOff(GCM_INTERNAL_STATE%alarmOcn, RC=STATUS)
        VERIFY_(STATUS)

       call MAPL_TimerOff(MAPL,"--O2A"  )

     endif
     RETURN_(ESMF_SUCCESS)

   end subroutine RUN_OCEAN

   subroutine DO_A2O(STATEO,NAMEO,STATEA,NAMEA,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_A2O"
     integer                               :: STATUS

     real,    pointer :: ptrA(:) => null()
     real,    pointer :: ptrO(:) => null()

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        call MAPL_LocStreamTransform( ptrO, XFORM_A2O, ptrA, RC=STATUS ) 
        VERIFY_(STATUS)
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_A2O

   subroutine DO_A2O_UGD(STATEO,NAMEO,STATEA,NAMEA,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_A2O_UGD"
     integer                               :: STATUS
     integer                               :: N

     real,    pointer :: ptrA(:,:) => null()
     real,    pointer :: ptrO(:,:) => null()

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N = 1, size(ptrA,2)
           call MAPL_LocStreamTransform( ptrO(:,N), XFORM_A2O, ptrA(:,N), RC=STATUS ) 
           VERIFY_(STATUS)
        end do
     end if

     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_A2O_UGD

   subroutine DO_O2A(STATEA,NAMEA,STATEO,NAMEO,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_O2A"
     integer                               :: STATUS

     real,    pointer :: ptrA(:) => null()
     real,    pointer :: ptrO(:) => null()

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        call MAPL_LocStreamTransform( ptrA, XFORM_O2A, ptrO, RC=STATUS ) 
        VERIFY_(STATUS)
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_O2A

   subroutine DO_A2O_SUBTILES_R4R4(STATEO,NAMEO,SUBINDEXO,STATEA,NAMEA,SUBINDEXA,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     integer                   , intent(IN   ) ::  SUBINDEXO(:)
     integer                   , intent(IN   ) ::  SUBINDEXA(:)
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_A2O_SUBTILES_R4R4"
     integer                               :: STATUS

     real,                       pointer :: ptrA(:,:) => null()
     real,                       pointer :: ptrO(:,:) => null()
     integer                             :: N, DIMSO, DIMSA  

     DIMSO = size(SUBINDEXO)
     DIMSA = size(SUBINDEXA)
     ASSERT_(DIMSO == DIMSA)

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N=1,DIMSO
         call MAPL_LocStreamTransform( ptrO(:,SUBINDEXO(N)), XFORM_A2O, &
                                       ptrA(:,SUBINDEXA(N)), RC=STATUS ) 
         VERIFY_(STATUS)
        enddo
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_A2O_SUBTILES_R4R4

   subroutine DO_A2O_SUBTILES_R4R8(STATEO,NAMEO,SUBINDEXO,STATEA,NAMEA,SUBINDEXA,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     integer                   , intent(IN   ) ::  SUBINDEXO(:)
     integer                   , intent(IN   ) ::  SUBINDEXA(:)
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_A2O_SUBTILES_R4R8"
     integer                               :: STATUS

     real(kind=ESMF_KIND_R8),    pointer :: ptrA(:,:) => null()
     real,                       pointer :: ptrO(:,:) => null()
     integer                             :: N, DIMSO, DIMSA  

     DIMSO = size(SUBINDEXO)
     DIMSA = size(SUBINDEXA)
     ASSERT_(DIMSO == DIMSA)

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N=1,DIMSO
         call MAPL_LocStreamTransform( ptrO(:,SUBINDEXO(N)), XFORM_A2O, &
                                       ptrA(:,SUBINDEXA(N)), RC=STATUS ) 
         VERIFY_(STATUS)
        enddo
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_A2O_SUBTILES_R4R8

   subroutine DO_A2O_SUBTILES2D_R4R8(STATEO,NAMEO,SUBINDEXO,STATEA,NAMEA,SUBINDEXA, &
                                     DIMS, RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     integer                   , intent(IN   ) ::  SUBINDEXO(:)
     integer                   , intent(IN   ) ::  SUBINDEXA(:)
     integer                   , intent(IN   ) ::  DIMS
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_A2O_SUBTILES2D_R4R8"
     integer                               :: STATUS

     real(kind=ESMF_KIND_R8),    pointer :: ptrA(:,:,:) => null()
     real,                       pointer :: ptrO(:,:,:) => null()
     integer                             :: N, K, DIMSO, DIMSA  
    

     DIMSO = size(SUBINDEXO)
     DIMSA = size(SUBINDEXA)
     ASSERT_(DIMSO == DIMSA)

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N=1,DIMSO
           do K=1,DIMS
             call MAPL_LocStreamTransform( ptrO(:,K,SUBINDEXO(N)), XFORM_A2O, &
                                           ptrA(:,K,SUBINDEXA(N)), RC=STATUS ) 
             VERIFY_(STATUS)
           enddo
        enddo
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_A2O_SUBTILES2D_R4R8

   subroutine DO_O2A_SUBTILES_R4R4(STATEA,NAMEA,SUBINDEXA,STATEO,NAMEO,SUBINDEXO,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     integer                   , intent(IN   ) ::  SUBINDEXO(:)
     integer                   , intent(IN   ) ::  SUBINDEXA(:)
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_O2A_SUBTILES_R4R4"
     integer                               :: STATUS

     real,                       pointer :: ptrA(:,:) => null()
     real,                       pointer :: ptrO(:,:) => null()
     integer                             :: N, DIMSO, DIMSA  

     DIMSO = size(SUBINDEXO)
     DIMSA = size(SUBINDEXA)
     ASSERT_(DIMSO == DIMSA)

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N=1,DIMSO
           call MAPL_LocStreamTransform( ptrA(:, SUBINDEXA(N)), XFORM_O2A, &
                                         ptrO(:, SUBINDEXO(N)), RC=STATUS ) 
           VERIFY_(STATUS)
        enddo
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_O2A_SUBTILES_R4R4

   subroutine DO_O2A_SUBTILES_R8R4(STATEA,NAMEA,SUBINDEXA,STATEO,NAMEO,SUBINDEXO,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     integer                   , intent(IN   ) ::  SUBINDEXO(:)
     integer                   , intent(IN   ) ::  SUBINDEXA(:)
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_O2A_SUBTILES_R8R4"
     integer                               :: STATUS

     real(kind=ESMF_KIND_R8),    pointer :: ptrA(:,:) => null()
     real,                       pointer :: ptrO(:,:) => null()
     integer                             :: N, DIMSO, DIMSA  

     DIMSO = size(SUBINDEXO)
     DIMSA = size(SUBINDEXA)
     ASSERT_(DIMSO == DIMSA)

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N=1,DIMSO
           call MAPL_LocStreamTransform( ptrA(:, SUBINDEXA(N)), XFORM_O2A, &
                                         ptrO(:, SUBINDEXO(N)), RC=STATUS ) 
           VERIFY_(STATUS)
        enddo
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_O2A_SUBTILES_R8R4

   subroutine DO_O2A_SUBTILES2D_R8R4(STATEA,NAMEA,SUBINDEXA,STATEO,NAMEO,SUBINDEXO, &
                                     DIMS, RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     integer                   , intent(IN   ) ::  SUBINDEXO(:)
     integer                   , intent(IN   ) ::  SUBINDEXA(:)
     integer                   , intent(IN   ) ::  DIMS
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_O2A_SUBTILES2D_R8R4"
     integer                               :: STATUS

     real(kind=ESMF_KIND_R8),    pointer :: ptrA(:,:,:) => null()
     real,                       pointer :: ptrO(:,:,:) => null()
     integer                             :: N, K, DIMSO, DIMSA  

     DIMSO = size(SUBINDEXO)
     DIMSA = size(SUBINDEXA)
     ASSERT_(DIMSO == DIMSA)

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     if (associated(ptrO) .and. associated(ptrA)) then
        do N=1,DIMSO
          do K=1,DIMS
            call MAPL_LocStreamTransform( ptrA(:,K,SUBINDEXA(N)), XFORM_O2A, &
                                          ptrO(:,K,SUBINDEXO(N)), RC=STATUS ) 
            VERIFY_(STATUS)
          enddo
        enddo
     end if
     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_O2A_SUBTILES2D_R8R4

#ifdef USE_ODAS
   subroutine DO_A2OG(STATEO,NAMEO,STATEA,NAMEA,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_A2O"
     integer                               :: STATUS

     real,    pointer  :: ptrA(:) => null()
     real,    pointer  :: ptrO_2d(:,:) => null()
     real, allocatable :: varO(:)
     integer           :: nt

     call MAPL_GetPointer(STATEO, ptrO_2d, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     if (.not. associated(ptrO_2d)) then
        RETURN_(ESMF_SUCCESS)
     end if

     call MAPL_GetPointer(STATEA, ptrA, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     if (.not. associated(ptrA)) then
        RETURN_(ESMF_SUCCESS)
     end if

     call MAPL_LocStreamGet( LOCSTREAM_O, nt_local=nt, RC=STATUS ) 
     VERIFY_(STATUS)
     allocate(varO(NT), stat=status)
     VERIFY_(STATUS)

     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
!    T2T (atmosphere-to-ocean tile)
     call MAPL_LocStreamTransform( varO, XFORM_A2O, ptrA, RC=STATUS ) 
     VERIFY_(STATUS)
!    T2G (ocean-tile-to-ocean-grid)
     call MAPL_LocStreamTransform( LOCSTREAM_O, ptrO_2d, varO, RC=STATUS ) 
     VERIFY_(STATUS)

     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     deallocate(varO)
     

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_A2OG



   subroutine DO_O2AG(STATEA,NAMEA,STATEO,NAMEO,RC)
     type(ESMF_State)          , intent(INOUT) ::  STATEA
     type(ESMF_State)          , intent(INOUT) ::  STATEO
     character(len=*)          , intent(IN   ) ::  NAMEA
     character(len=*)          , intent(IN   ) ::  NAMEO
     integer, optional,          intent(  OUT) ::  RC

     character(len=ESMF_MAXSTR), parameter :: IAm="DO_O2A"
     integer                               :: STATUS

     real,    pointer  :: ptrO(:) => null()
     real,    pointer  :: ptrA_2d(:,:) => null()
     real, allocatable :: varA(:)
     integer           :: nt

     call MAPL_GetPointer(STATEA, ptra_2d, NAMEA, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     if (.not. associated(ptra_2d)) then
        RETURN_(ESMF_SUCCESS)
     end if

     call MAPL_GetPointer(STATEO, ptrO, NAMEO, notFoundOK=.true., RC=STATUS)
     VERIFY_(STATUS)
     if (.not. associated(ptrO)) then
        RETURN_(ESMF_SUCCESS)
     end if

     call MAPL_LocStreamGet( LOCSTREAM_A, nt_local=nt, RC=STATUS ) 
     VERIFY_(STATUS)
     allocate(varA(NT), stat=status)
     VERIFY_(STATUS)

     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
!    T2T (ocean-to-atmosphere tile)
     call MAPL_LocStreamTransform( varA, XFORM_O2A, ptrO, RC=STATUS ) 
     VERIFY_(STATUS)
!    T2G (atmosphere-tile-to-atmosphere-grid)
     call MAPL_LocStreamTransform( LOCSTREAM_A, ptrA_2d, varA, RC=STATUS ) 
     VERIFY_(STATUS)

     call ESMF_VMBarrier(VM, rc=status)
     VERIFY_(STATUS)
     deallocate(varA)
     

     RETURN_(ESMF_SUCCESS)
   end subroutine DO_O2AG
#endif

 end subroutine Run
#ifdef USE_ODAS
   subroutine Finalize ( GC, IMPORT, EXPORT, CLOCK, RC )

! !ARGUMENTS:

     type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
     type(ESMF_State),    intent(inout) :: IMPORT ! Import state
     type(ESMF_State),    intent(inout) :: EXPORT ! Export state
     type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
     integer, optional,   intent(  out) :: RC     ! Error code
     
! !DESCRIPTION: The Finalize method of GCM.

!EOP

! ErrLog Variables

     character(len=ESMF_MAXSTR)          :: IAm 
     integer                             :: STATUS
     character(len=ESMF_MAXSTR)          :: COMP_NAME
     
      Iam = "Finalize"
      call ESMF_GridCompGet ( GC, name=COMP_NAME, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // Iam


      if(COMP_NAME(1:3) == "ens") then 
         status = chdir("restarts." // comp_name(5:7));
         VERIFY_(STATUS)
      end if
      
      call MAPL_GenericFinalize ( GC, IMPORT, EXPORT, CLOCK, rc = STATUS)
      VERIFY_(STATUS)

      if(COMP_NAME(1:3) == "ens") then 
         status = chdir("..");
         VERIFY_(STATUS) 
      end if
      
      RETURN_(ESMF_SUCCESS)

   end subroutine Finalize

  subroutine Record ( gc, import, export, clock, rc )

! !ARGUMENTS:

  type(ESMF_GridComp), intent(INOUT) :: gc     ! Gridded component 
  type(ESMF_State),    intent(INOUT) :: import ! Import state
  type(ESMF_State),    intent(INOUT) :: export ! Export state
  type(ESMF_Clock),    intent(INOUT) :: clock  ! The supervisor clock
  integer, optional,   intent(  OUT) :: rc     ! Error code:

!EOP

    type (MAPL_MetaComp), pointer    :: MAPL => null()

! ErrLog Variables

    character(len=ESMF_MAXSTR)	     :: IAm
    integer			     :: STATUS
    character(len=ESMF_MAXSTR)       :: COMP_NAME

! Locals
    logical                          :: doRecord

! Get the target components name and set-up traceback handle.
! -----------------------------------------------------------

    Iam = "Record"
    call ESMF_GridCompGet( gc, NAME=comp_name, RC=status )
    VERIFY_(STATUS)
    Iam = trim(comp_name) // Iam

! Get my internal MAPL_Generic state

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=status)
    VERIFY_(STATUS)

    doRecord = MAPL_RecordAlarmIsRinging(MAPL, RC=status)
    VERIFY_(STATUS)

    if (doRecord) then
      status = chdir("restarts." // comp_name(5:7));
      VERIFY_(STATUS) 
    end if

    call MAPL_GenericRecord ( GC, IMPORT, EXPORT, CLOCK, rc = STATUS)
    VERIFY_(STATUS)

    if (doRecord) then
      status = chdir("..");
      VERIFY_(STATUS) 
    end if
    RETURN_(ESMF_SUCCESS)

  end subroutine Record

#endif
end module GEOS_GcmGridCompMod

