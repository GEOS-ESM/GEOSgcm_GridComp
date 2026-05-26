! $Id$

#include "MAPL_Generic.h"

!#define UWDIAG 1

!=============================================================================
!BOP

! !MODULE: GEOS_UW_InterfaceMod -- A Module to interface with the
!   UW convection

module GEOS_UW_InterfaceMod

  use ESMF
  use MAPL
  use UWSHCU   ! using module that contains uwshcu code
  use GEOSmoist_Process_Library

  implicit none

  integer USE_TRACER_TRANSP_UW      ! transport tracers in UW
  real    :: SCLM_SHALLOW
  logical :: JASON_UW, JASON_MFD_SC
  logical :: REPORT_UW_NEGATIVES
  logical :: USE_EIS

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  public :: UW_Setup, UW_Initialize, UW_Run
   
contains

subroutine UW_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_UW_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME ='CUSH',                                         &
         LONG_NAME  = 'Cumulus_scale_height_from_UW_shlw_convection',&
         UNITS      ='m',                                            &
         DIMS       = MAPL_DimsHorzOnly,                             &
         VLOCATION  = MAPL_VLocationNone,                            &
         DEFAULT= 1000.0,                                            &
         RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC, name="--UW", RC=STATUS)
    VERIFY_(STATUS)

end subroutine UW_Setup

subroutine UW_Initialize (MAPL, CLOCK, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional                   :: RC  ! return code
    integer :: LM

    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: MOIST_DT
    real                            :: UW_DT

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Alarm)        :: UW_RunAlarm
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

    call MAPL_Get(MAPL, RUNALARM=ALARM, LM=LM, RC=STATUS );VERIFY_(STATUS)
    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    MOIST_DT = DT_R8
    call MAPL_GetResource(MAPL, UW_DT, 'UW_DT:', default=MOIST_DT, RC=STATUS); VERIFY_(STATUS)

    call ESMF_ClockGet(CLOCK, calendar=calendar, currTime=currentTime, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalSet(ringInterval, S=nint(UW_DT), calendar=calendar, RC=STATUS); VERIFY_(STATUS)

    UW_RunAlarm = ESMF_AlarmCreate(Clock        = CLOCK,        &
                                   Name         = 'UW_RunAlarm',&
                                   RingInterval = ringInterval, &
                                   RingTime     = currentTime,  &
                                   Enabled      = .true.   ,    &
                                   Sticky       = .false.  , RC=STATUS); VERIFY_(STATUS)

    call MAPL_Get ( MAPL, LM=LM, RC=STATUS )
    VERIFY_(STATUS)

                JASON_UW = .FALSE.
    if (LM==72) JASON_UW = .TRUE.
    call MAPL_GetResource(MAPL, JASON_UW, 'JASON_UW:', default=JASON_UW, RC=STATUS) ; VERIFY_(STATUS)

                JASON_MFD_SC = .FALSE.
    if (LM==72) JASON_MFD_SC = .TRUE.
    call MAPL_GetResource(MAPL, JASON_MFD_SC, 'JASON_MFD_SC:', default=JASON_MFD_SC, RC=STATUS) ; VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, REPORT_UW_NEGATIVES, 'REPORT_UW_NEGATIVES:', default=.FALSE., RC=STATUS) ; VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, USE_TRACER_TRANSP_UW,        'USE_TRACER_TRANSP_UW:',default= 1      , RC=STATUS) ; VERIFY_(STATUS)
    if (JASON_UW) then
      call MAPL_GetResource(MAPL, SHLWPARAMS%WINDSRCAVG,       'WINDSRCAVG:'      ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%MIXSCALE,         'MIXSCALE:'        ,DEFAULT=0.0,    RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%CRIQC,            'CRIQC:'           ,DEFAULT=1.0e-3, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%THLSRC_FAC,       'THLSRC_FAC:'      ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRC_FAC,        'QTSRC_FAC:'       ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRCHGT,         'QTSRCHGT:'        ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKFRE,            'RKFRE:'           ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKM,              'RKM:'             ,DEFAULT= 12.0,  RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RMAXFRAC,         'RMAXFRAC:'        ,DEFAULT= 0.1,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%FRC_RASN,         'FRC_RASN:'        ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RPEN,             'RPEN:'            ,DEFAULT= 3.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SCLM_SHALLOW,                'SCLM_SHALLOW:'    ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%NITER_XC,         'NITER_XC:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, USE_EIS,                     'UW_USE_EIS:'      ,DEFAULT=.FALSE.,RC=STATUS) ; VERIFY_(STATUS)
    else
      call MAPL_GetResource(MAPL, SHLWPARAMS%WINDSRCAVG,       'WINDSRCAVG:'      ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%MIXSCALE,         'MIXSCALE:'        ,DEFAULT=3000.0, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%MIXSCALE_HR,      'MIXSCALE_HR:'     ,DEFAULT=3000.0, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%CRIQC,            'CRIQC:'           ,DEFAULT=0.9e-3, RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%THLSRC_FAC,       'THLSRC_FAC:'      ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRC_FAC,        'QTSRC_FAC:'       ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%QTSRCHGT,         'QTSRCHGT:'        ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKFRE,            'RKFRE:'           ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKFRE_HR,         'RKFRE_HR:'        ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKM,              'RKM:'             ,DEFAULT= 12.0,  RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RKM_HR,           'RKM_HR:'          ,DEFAULT= 12.0,  RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RMAXFRAC,         'RMAXFRAC:'        ,DEFAULT= 0.1,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RMAXFRAC_HR,      'RMAXFRAC_HR:'     ,DEFAULT= 0.1,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%FRC_RASN,         'FRC_RASN:'        ,DEFAULT= 0.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%RPEN,             'RPEN:'            ,DEFAULT= 3.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SCLM_SHALLOW,                'SCLM_SHALLOW:'    ,DEFAULT= 1.0,   RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, SHLWPARAMS%NITER_XC,         'NITER_XC:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, USE_EIS,                     'UW_USE_EIS:'      ,DEFAULT=.FALSE.,RC=STATUS) ; VERIFY_(STATUS)
    endif
    call MAPL_GetResource(MAPL, SHLWPARAMS%ITER_CIN,         'ITER_CIN:'        ,DEFAULT=2,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_CINCIN,       'USE_CINCIN:'      ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%CRIDIST_OPT,      'CRIDIST_OPT:'     ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_SELF_DETRAIN, 'USE_SELF_DETRAIN:',DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_MOMENFLX,     'USE_MOMENFLX:'    ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%USE_CUMPENENT,    'USE_CUMPENENT:'   ,DEFAULT=1,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%SCVERBOSE,        'SCVERBOSE:'       ,DEFAULT=0,      RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RLE,              'RLE:'             ,DEFAULT=0.1,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%MUMIN1,           'MUMIN1:'          ,DEFAULT=0.906,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RBUOY,            'RBUOY:'           ,DEFAULT=1.0,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RDRAG,            'RDRAG:'           ,DEFAULT=1.0,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%EPSVARW,          'EPSVARW:'         ,DEFAULT=5.e-4,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%PGFC,             'PGFC:'            ,DEFAULT=0.7,    RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%KEVP,             'KEVP:'            ,DEFAULT=2.e-6,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%RDROP,            'SHLW_RDROP:'      ,DEFAULT=8.e-6,  RC=STATUS) ; VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, SHLWPARAMS%DETRHGT,          'DETRHGT:'         ,DEFAULT=1800.0, RC=STATUS) ; VERIFY_(STATUS)

end subroutine UW_Initialize

subroutine UW_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:)   :: CUSH

    ! Imports
    real, pointer, dimension(:,:)   :: FRLAND, SH, EVAP, KPBL_SC
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, TKE

    ! allocatable derived quantities
    real,    allocatable, dimension(:,:,:) :: ZLE0, ZL0
    real,    allocatable, dimension(:,:,:) :: PL, PK, PKE, DP
    real,    allocatable, dimension(:,:,:) :: MASS
    real,    allocatable, dimension(:,:)   :: RKM2D, RKFRE, MIX2D, RMAXFRAC2D
    real,    allocatable, dimension(:,:,:) :: TMP3D

    ! Required Exports (connectivities to moist siblings)
    real, pointer, dimension(:,:)   :: CNPCPRATE
    real, pointer, dimension(:,:)   :: CNV_FRC, SRF_TYPE
    ! Exports
    real, pointer, dimension(:,:,:) :: CUFRC_SC
    real, pointer, dimension(:,:,:) :: UMF_SC, MFD_SC, DCM_SC
    real, pointer, dimension(:,:,:) :: QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC
    real, pointer, dimension(:,:,:) :: DTDT_SC, DQVDT_SC, DQRDT_SC, DQSDT_SC, &
                                       DUDT_SC,  DVDT_SC, DQIDT_SC, DQLDT_SC, DQADT_SC
    real, pointer, dimension(:,:,:) :: ENTR_SC, DETR_SC, QLDET_SC, &
                                       QIDET_SC, QLENT_SC, QIENT_SC, &
                                       QLSUB_SC, QISUB_SC, SC_NDROP, SC_NICE
    real, pointer, dimension(:,:)   :: TPERT_SC, QPERT_SC, LTS, EIS
    real, pointer, dimension(:,:)   :: CBMF_SC, PLCL_SC, PLFC_SC,     &
                                       PINV_SC, PREL_SC, PBUP_SC,    &
                                       CLDTOP_SC, SC_QT, SC_MSE
#ifdef UWDIAG
    real, pointer, dimension(:,:)   :: CIN_SC, CNT_SC, CNB_SC,      &
                                       WLCL_SC, QTSRC_SC, THLSRC_SC, &
                                       THVLSRC_SC, TKEAVG_SC
    real, pointer, dimension(:,:,:) :: SHL_DQCDT, WUP_SC, QTUP_SC,   &
                                       THLUP_SC, THVUP_SC, UUP_SC,   & 
                                       VUP_SC, XC_SC, QCU_SC,        &
                                       QLU_SC, QIU_SC
#endif
    real, pointer, dimension(:,:,:) :: QLTOT, QITOT
    real, pointer, dimension(:,:,:) ::   DQVDT_FILL
    real, pointer, dimension(:,:,:) :: DQLLSDT_FILL
    real, pointer, dimension(:,:,:) :: DQLCNDT_FILL
    real, pointer, dimension(:,:,:) :: DQILSDT_FILL
    real, pointer, dimension(:,:,:) :: DQICNDT_FILL
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:)   :: PTR2D

    real, pointer, dimension(:,:)   :: LONS, LATS

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: UW_DT, MOIST_DT
    real                            :: SIG
    type(ESMF_Alarm)                :: alarm
    logical                         :: alarm_is_ringing
    type( ESMF_VM )                 :: VMG

    ! Local variables

    real :: fac_eis                    ! Estimated enversion strength 0:1 factor
    real :: rkfre_base                 ! Base fractional entrainment rate before EIS modification
    real :: rkm_base                   ! Base momentum entrainment rate before EIS modification  
    real :: mix2d_base                 ! Base mixing length scale before EIS modification
    real :: rmaxfrac_base              ! Base maximum updraft area fraction before EIS modification
    real :: eis_rkfre_factor           ! EIS modification factor for RKFRE [0-1]
    real :: eis_rkm_factor             ! EIS modification factor for RKM
    real :: eis_mix2d_factor           ! EIS modification factor for MIX2D [0-1]
    real :: eis_rmaxfrac_factor        ! EIS modification factor for RMAXFRAC [1.0-1.1]

    integer                         :: I, J, L, K
    integer                         :: IM,JM,LM

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS); VERIFY_(STATUS)
    call MAPL_Get( MAPL, RUNALARM=ALARM, &
         INTERNAL_ESMF_STATE=INTERNAL, IM=IM, JM=JM, LM=LM, &
         LONS     = LONS,              &
         LATS     = LATS,              &
         RC=STATUS )
    VERIFY_(STATUS)
    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    MOIST_DT = DT_R8

    call ESMF_GridCompGet ( GC, VM=VMG, RC=STATUS )
    VERIFY_(STATUS)

    ! Internals
    call MAPL_GetPointer(INTERNAL, Q,      'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,   'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,   'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,   'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,   'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,   'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,   'QICN'    , RC=STATUS); VERIFY_(STATUS)
    ! Imports
    call MAPL_GetPointer(IMPORT, T         ,'T'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U         ,'U'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V         ,'V'         ,RC=STATUS); VERIFY_(STATUS)
    ! Tendency Export
    call MAPL_GetPointer(EXPORT, DUDT_SC,    'DUDT_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DVDT_SC,    'DVDT_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTDT_SC,    'DTDT_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQVDT_SC,   'DQVDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDT_SC,   'DQIDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDT_SC,   'DQLDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQRDT_SC,   'DQRDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQSDT_SC,   'DQSDT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQADT_SC,   'DQADT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Lower tropospheric stability and estimated inversion strength from MoistGC
    call MAPL_GetPointer(EXPORT, LTS,   'LTS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EIS,   'EIS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(EXPORT, PLCL_SC,   'PLCL_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PLFC_SC,   'PLFC_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PINV_SC,   'PINV_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PREL_SC,   'PREL_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PBUP_SC,   'PBUP_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CBMF_SC,   'CBMF_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CLDTOP_SC, 'CLDTOP_SC' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

#ifdef UWDIAG
      call MAPL_GetPointer(EXPORT, SHL_DQCDT, 'SHL_DQCDT' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNT_SC,    'CNT_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNB_SC,    'CNB_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, WLCL_SC,   'WLCL_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QTSRC_SC,  'QTSRC_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THLSRC_SC, 'THLSRC_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THVLSRC_SC,'THVLSRC_SC', RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, TKEAVG_SC, 'TKEAVG_SC' , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, WUP_SC,    'WUP_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QTUP_SC,   'QTUP_SC'   , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THLUP_SC,  'THLUP_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, THVUP_SC,  'THVUP_SC'  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, UUP_SC,    'UUP_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, VUP_SC,    'VUP_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, XC_SC,     'XC_SC'     , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QCU_SC,    'QCU_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QLU_SC,    'QLU_SC'    , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QIU_SC,    'QIU_SC'    , RC=STATUS); VERIFY_(STATUS)
#endif
    
    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--UW")

    ! Get parameters from generic state.
    !-----------------------------------

    ! Internals
    call MAPL_GetPointer(INTERNAL, CUSH,   'CUSH'    , RC=STATUS); VERIFY_(STATUS)

    ! Imports
    call MAPL_GetPointer(IMPORT, FRLAND    ,'FRLAND'    ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ZLE       ,'ZLE'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE       ,'PLE'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SH        ,'SH'        ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP      ,'EVAP'      ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBL_SC   ,'KPBL_SC'   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TKE       ,'TKE'       ,RC=STATUS); VERIFY_(STATUS)

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PKE  (IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PL   (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( DP   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( RKFRE  (IM,JM) )
    ALLOCATE ( RKM2D  (IM,JM) )
    ALLOCATE ( MIX2D  (IM,JM) )
    ALLOCATE ( RMAXFRAC2D (IM,JM) )

    ! Derived States
    !--------------------------------------------------------------
      
    ! 1. Compute the top-of-atmosphere edge (k = 0)
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, PKE, PLE, ZLE0, ZLE) &
    !$OMP PRIVATE(i, j)
    do j = 1, JM
       !DIR$ IVDEP
       do i = 1, IM
          PKE(i,j,0)  = (PLE(i,j,0) / MAPL_P00)**(MAPL_KAPPA)
          ZLE0(i,j,0) = ZLE(i,j,0) - ZLE(i,j,LM)
       end do
    end do
    !$OMP END PARALLEL DO

    ! 2. Compute the remaining edges and all layer variables (k = 1 to LM)
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, PKE, PLE, PL, PK, &
    !$OMP        ZLE0, ZLE, ZL0, DP, MASS) &
    !$OMP PRIVATE(i, j, k)
    do k = 1, LM
       do j = 1, JM
          !DIR$ IVDEP
          do i = 1, IM
               
             ! Edge variables
             PKE(i,j,k)  = (PLE(i,j,k) / MAPL_P00)**(MAPL_KAPPA)
             ZLE0(i,j,k) = ZLE(i,j,k) - ZLE(i,j,LM)
               
             ! Layer variables
             PL(i,j,k)   = 0.5 * (PLE(i,j,k-1) + PLE(i,j,k))
             PK(i,j,k)   = (PL(i,j,k) / MAPL_P00)**(MAPL_KAPPA)
               
             ZL0(i,j,k)  = 0.5 * (ZLE0(i,j,k-1) + ZLE0(i,j,k))
               
             DP(i,j,k)   = PLE(i,j,k) - PLE(i,j,k-1)
             MASS(i,j,k) = DP(i,j,k) / MAPL_GRAV
               
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call ESMF_ClockGetAlarm(clock, 'UW_RunAlarm', alarm, RC=STATUS); VERIFY_(STATUS)
    alarm_is_ringing = ESMF_AlarmIsRinging(alarm, RC=STATUS); VERIFY_(STATUS)

    if (alarm_is_ringing) then

!!! call WRITE_PARALLEL('UW is Running')
    call ESMF_AlarmRingerOff(alarm, RC=STATUS); VERIFY_(STATUS)
    call ESMF_AlarmGet(alarm, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    UW_DT = DT_R8

    ! Required Exports (connectivities to moist siblings)
    call MAPL_GetPointer(EXPORT, MFD_SC,     'MFD_SC'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QLDET_SC,   'QLDET_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QIDET_SC,   'QIDET_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CUFRC_SC,   'CUFRC_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNPCPRATE,  'CNPCPRATE' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_FRC ,   'CNV_FRC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SRF_TYPE,   'SRF_TYPE'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Exports
    call MAPL_GetPointer(EXPORT, UMF_SC,     'UMF_SC'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DCM_SC,     'DCM_SC'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ENTR_SC,    'ENTR_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DETR_SC,    'DETR_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QLSUB_SC,   'QLSUB_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QISUB_SC,   'QISUB_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_NDROP,   'SC_NDROP'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_NICE,    'SC_NICE'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TPERT_SC,   'TPERT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QPERT_SC,   'QPERT_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QTFLX_SC,   'QTFLX_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SLFLX_SC,   'SLFLX_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, UFLX_SC,    'UFLX_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFLX_SC,    'VFLX_SC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

    ! 1. Fetch all pointers first
    !--------------------------------------------------------------
    call MAPL_GetPointer(IMPORT, PTR2D, 'AREA', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QLTOT, 'QLTOT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QITOT, 'QITOT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

    ! 2. 2D parameters for UW
    !--------------------------------------------------------------
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, JASON_UW, SHLWPARAMS, RKFRE, RKM2D, MIX2D, RMAXFRAC2D, &
    !$OMP        USE_EIS, EIS, srf_type, PTR2D) &
    !$OMP PRIVATE(i, j, fac_eis, SIG, rkfre_base, rkm_base, mix2d_base, rmaxfrac_base, &
    !$OMP        eis_rkfre_factor, eis_rkm_factor, eis_mix2d_factor, eis_rmaxfrac_factor)
    do j = 1, JM
       !DIR$ IVDEP
       do i = 1, IM
          if (JASON_UW) then
             RKFRE(i,j)      = SHLWPARAMS%RKFRE
             RKM2D(i,j)      = SHLWPARAMS%RKM
             MIX2D(i,j)      = SHLWPARAMS%MIXSCALE
             RMAXFRAC2D(i,j) = SHLWPARAMS%RMAXFRAC
          else
             fac_eis = 0.0
             if (USE_EIS) fac_eis = get_fac_eis(EIS(i,j), srf_type(i,j))
             SIG = SIGMA(SQRT(PTR2D(i,j)))

             ! Base resolution-dependent parameters
             rkfre_base    = SHLWPARAMS%RKFRE    * SIG + SHLWPARAMS%RKFRE_HR    * (1.0 - SIG)
             rkm_base      = SHLWPARAMS%RKM      * SIG + SHLWPARAMS%RKM_HR      * (1.0 - SIG)
             mix2d_base    = SHLWPARAMS%MIXSCALE * SIG + SHLWPARAMS%MIXSCALE_HR * (1.0 - SIG)
             rmaxfrac_base = SHLWPARAMS%RMAXFRAC * SIG + SHLWPARAMS%RMAXFRAC_HR * (1.0 - SIG)

             ! EIS-based regime modifications
             eis_rkfre_factor    = 1.0 - 0.8 * fac_eis
             eis_rkm_factor      = 1.0 + 0.4 * fac_eis
             eis_mix2d_factor    = 1.0 - 0.3 * fac_eis
             eis_rmaxfrac_factor = 1.0 + 0.1 * fac_eis

             ! Apply EIS modifications
             RKFRE(i,j)      = rkfre_base * eis_rkfre_factor
             RKM2D(i,j)      = rkm_base   * eis_rkm_factor
             MIX2D(i,j)      = mix2d_base * eis_mix2d_factor
             RMAXFRAC2D(i,j) = rmaxfrac_base * eis_rmaxfrac_factor

             ! Optional: Add minimum limits
             RKFRE(i,j)      = max(RKFRE(i,j), 0.1)
             RKM2D(i,j)      = min(RKM2D(i,j), 14.0)
             MIX2D(i,j)      = max(MIX2D(i,j), 1500.0)
             RMAXFRAC2D(i,j) = max(min(RMAXFRAC2D(i,j), 0.8), 0.05)
          end if
       end do
    end do
    !$OMP END PARALLEL DO

    ! 3. Combine condensates for input (not updated within UW) 
    !--------------------------------------------------------------
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, QLTOT, QLLS, QLCN, QITOT, QILS, QICN) &
    !$OMP PRIVATE(i, j, k)
    do k = 1, LM
       do j = 1, JM
          !DIR$ IVDEP
          do i = 1, IM
             QLTOT(i,j,k) = QLLS(i,j,k) + QLCN(i,j,k)
             QITOT(i,j,k) = QILS(i,j,k) + QICN(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

      !  Call UW shallow convection
      !----------------------------------------------------------------
      call compute_uwshcu_inv(IM*JM, LM, UW_DT,           & ! IN
            PL, ZL0, PK, PLE, ZLE0, PKE, DP,              &
            U, V, Q, QLTOT, QITOT, T, TKE, RKFRE, KPBL_SC,&
            SH, EVAP, CNPCPRATE, FRLAND, RKM2D, MIX2D, RMAXFRAC2D, &
            CUSH,                                         & ! INOUT
            UMF_SC, DCM_SC, DQVDT_SC, DQLDT_SC, DQIDT_SC, & ! OUT
            DTDT_SC, DUDT_SC, DVDT_SC, DQRDT_SC,          &
            DQSDT_SC, CUFRC_SC, ENTR_SC, DETR_SC,         &
            QLDET_SC, QIDET_SC, QLSUB_SC, QISUB_SC,       &
            SC_NDROP, SC_NICE, TPERT_SC, QPERT_SC,        &
            QTFLX_SC, SLFLX_SC, UFLX_SC, VFLX_SC,         &
            CBMF_SC, PLCL_SC, PLFC_SC, PINV_SC,           & ! DIAG ONLY
            PREL_SC, PBUP_SC, CLDTOP_SC,                  &
#ifdef UWDIAG 
            QCU_SC, QLU_SC,                               & 
            QIU_SC, SHL_DQCDT, CNT_SC, CNB_SC,            &
            CIN_SC, WLCL_SC, QTSRC_SC, THLSRC_SC,         &
            THVLSRC_SC, TKEAVG_SC, CLDTOP_SC, WUP_SC,     &
            QTUP_SC, THLUP_SC, THVUP_SC, UUP_SC, VUP_SC,  &
            XC_SC,                                        &
#endif 
            USE_TRACER_TRANSP_UW)

      ! 1. Fetch ALL pointers at the top
      !--------------------------------------------------------------
      call MAPL_GetPointer(EXPORT, QLENT_SC, 'QLENT_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, QIENT_SC, 'QIENT_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_QT,   'SC_QT',     RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, SC_MSE,  'SC_MSE',    RC=STATUS); VERIFY_(STATUS)

      ! 2. Fused 3D Loop for Detrainment and Conversions
      !--------------------------------------------------------------
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(IM, JM, LM, JASON_MFD_SC, DETR_SC, UMF_SC, DP, MFD_SC, &
      !$OMP        DCM_SC, DQADT_SC, SCLM_SHALLOW, MASS, QLENT_SC, QLDET_SC, QIENT_SC, QIDET_SC) &
      !$OMP PRIVATE(i, j, k)
      do k = 1, LM
         do j = 1, JM
            !DIR$ IVDEP 
            do i = 1, IM
               ! Calculate detrained mass flux
               if (JASON_MFD_SC) then
                  if (DETR_SC(i,j,k) /= MAPL_UNDEF) then
                     MFD_SC(i,j,k) = 0.5 * (UMF_SC(i,j,k) + UMF_SC(i,j,k-1)) * DETR_SC(i,j,k) * DP(i,j,k)
                  else
                     MFD_SC(i,j,k) = 0.0
                  end if
               else
                  MFD_SC(i,j,k) = DCM_SC(i,j,k)
               end if
               
               DQADT_SC(i,j,k) = MFD_SC(i,j,k) * SCLM_SHALLOW / MASS(i,j,k)

               ! Convert detrained water units before passing to cloud
               QLENT_SC(i,j,k) = 0.0
               QIENT_SC(i,j,k) = 0.0
               
               if (QLDET_SC(i,j,k) < 0.0) then
                  QLENT_SC(i,j,k) = QLDET_SC(i,j,k)
                  QLDET_SC(i,j,k) = 0.0
               end if
               
               if (QIDET_SC(i,j,k) < 0.0) then
                  QIENT_SC(i,j,k) = QIDET_SC(i,j,k)
                  QIDET_SC(i,j,k) = 0.0
               end if

               ! Scale the detrained fluxes before exporting
               QLDET_SC(i,j,k) = QLDET_SC(i,j,k) * MASS(i,j,k)
               QIDET_SC(i,j,k) = QIDET_SC(i,j,k) * MASS(i,j,k)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      ! 3. Whole-array copies for direct 3D variables
      !--------------------------------------------------------------
      call MAPL_GetPointer(EXPORT, PTR3D, 'SHLW_PRC3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) PTR3D = DQRDT_SC    ! [kg/kg/s]
      call MAPL_GetPointer(EXPORT, PTR3D, 'SHLW_SNO3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) PTR3D = DQSDT_SC    ! [kg/kg/s]
      call MAPL_GetPointer(EXPORT, PTR2D,  'CUSH_SC', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D = CUSH

      ! 4. Fused 2D Loop for Column Integrals
      !--------------------------------------------------------------
      ! We parallelize over j,i and accumulate over k internally
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP SHARED(IM, JM, LM, SC_QT, SC_MSE, DQSDT_SC, DQRDT_SC, DQVDT_SC, &
      !$OMP        QLENT_SC, QLSUB_SC, QIENT_SC, QISUB_SC, MASS, QLDET_SC, QIDET_SC, &
      !$OMP        DTDT_SC, DQIDT_SC) &
      !$OMP PRIVATE(i, j, k)
      do j = 1, JM
         do i = 1, IM
            if (associated(SC_QT)) SC_QT(i,j) = 0.0
            if (associated(SC_MSE)) SC_MSE(i,j) = 0.0
            
            do k = 1, LM
               if (associated(SC_QT)) then
                  SC_QT(i,j) = SC_QT(i,j) + &
                     ( DQSDT_SC(i,j,k) + DQRDT_SC(i,j,k) + DQVDT_SC(i,j,k) + &
                       QLENT_SC(i,j,k) + QLSUB_SC(i,j,k) + QIENT_SC(i,j,k) + &
                       QISUB_SC(i,j,k) ) * MASS(i,j,k) + &
                       QLDET_SC(i,j,k) + QIDET_SC(i,j,k)
               end if
               
               if (associated(SC_MSE)) then
                  SC_MSE(i,j) = SC_MSE(i,j) + &
                     ( MAPL_CP   * DTDT_SC(i,j,k) + &
                       MAPL_ALHL * DQVDT_SC(i,j,k) - &
                       MAPL_ALHF * DQIDT_SC(i,j,k) ) * MASS(i,j,k)
               end if
            end do
         end do
      end do
      !$OMP END PARALLEL DO

  endif


  ! 1. Fetch all pointers FIRST before doing any math
  !--------------------------------------------------------------
  call MAPL_GetPointer(EXPORT, QLDET_SC, 'QLDET_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, QIDET_SC, 'QIDET_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, QLSUB_SC, 'QLSUB_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, QLENT_SC, 'QLENT_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, QISUB_SC, 'QISUB_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, QIENT_SC, 'QIENT_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

  ! 2. Apply tendencies in a single fused loop with OpenMP
  !--------------------------------------------------------------
  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(IM, JM, LM, Q, DQVDT_SC, MOIST_DT, T, DTDT_SC, U, DUDT_SC, V, DVDT_SC, &
  !$OMP        CLCN, DQADT_SC, QLCN, QLDET_SC, MASS, QICN, QIDET_SC, &
  !$OMP        QLLS, QLSUB_SC, QLENT_SC, QILS, QISUB_SC, QIENT_SC) &
  !$OMP PRIVATE(i, j, k)
  do k = 1, LM
     do j = 1, JM
        !DIR$ IVDEP 
        !DIR$ VECTOR ALWAYS
        do i = 1, IM
           ! Apply tendencies
           Q(i,j,k) = Q(i,j,k) + DQVDT_SC(i,j,k) * MOIST_DT
           T(i,j,k) = T(i,j,k) + DTDT_SC(i,j,k)  * MOIST_DT
           U(i,j,k) = U(i,j,k) + DUDT_SC(i,j,k)  * MOIST_DT
           V(i,j,k) = V(i,j,k) + DVDT_SC(i,j,k)  * MOIST_DT
           
           ! Tiedtke-style cloud fraction 
           CLCN(i,j,k) = MAX(0.0, MIN(CLCN(i,j,k) + DQADT_SC(i,j,k)*MOIST_DT, 1.0))
           
           ! Add detrained shallow convective ice/liquid source
           QLCN(i,j,k) = MAX(0.0, QLCN(i,j,k) + QLDET_SC(i,j,k)*MOIST_DT/MASS(i,j,k))
           QICN(i,j,k) = MAX(0.0, QICN(i,j,k) + QIDET_SC(i,j,k)*MOIST_DT/MASS(i,j,k))
           
           ! Apply condensate tendency from subsidence, and sink from
           ! condensate entrained into shallow updraft. 
           QLLS(i,j,k) = MAX(0.0, QLLS(i,j,k) + (QLSUB_SC(i,j,k)+QLENT_SC(i,j,k))*MOIST_DT)
           QILS(i,j,k) = MAX(0.0, QILS(i,j,k) + (QISUB_SC(i,j,k)+QIENT_SC(i,j,k))*MOIST_DT)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

! Cleanup negative water species
! ------------------------------
  call MAPL_GetPointer(EXPORT,   DQVDT_FILL,   'DQVDT_FILL_SC', RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_SC', RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_SC', RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_SC', RC=STATUS); VERIFY_(STATUS)
  call MAPL_GetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_SC', RC=STATUS); VERIFY_(STATUS)
  if (REPORT_UW_NEGATIVES) then
    call FILLQ2ZERO( Q       , MASS, DT=MOIST_DT, DQDT=  DQVDT_FILL, WARNING_LABEL="QV   After UW ShallowCu", VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QLLS    , MASS, DT=MOIST_DT, DQDT=DQLLSDT_FILL, WARNING_LABEL="QLLS After UW ShallowCu", VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QLCN    , MASS, DT=MOIST_DT, DQDT=DQLCNDT_FILL, WARNING_LABEL="QLCN After UW ShallowCu", VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QILS    , MASS, DT=MOIST_DT, DQDT=DQILSDT_FILL, WARNING_LABEL="QILS After UW ShallowCu", VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QICN    , MASS, DT=MOIST_DT, DQDT=DQICNDT_FILL, WARNING_LABEL="QICN After UW ShallowCu", VM=VMG, RC=STATUS); VERIFY_(STATUS)
  else
    call FILLQ2ZERO( Q       , MASS, DT=MOIST_DT, DQDT=  DQVDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QLLS    , MASS, DT=MOIST_DT, DQDT=DQLLSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QLCN    , MASS, DT=MOIST_DT, DQDT=DQLCNDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QILS    , MASS, DT=MOIST_DT, DQDT=DQILSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QICN    , MASS, DT=MOIST_DT, DQDT=DQICNDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
  endif

  if (DEBUG_TQ_ERRORS) then
        do L=1,LM                
          do J=1,JM              
           do I=1,IM             
             if (T(I,J,L) > 333.0) then
                 print *, "Temperature spike detected : ", T(I,J,L)
                 print *, "         UW Temp Increment : ", DTDT_SC(I,J,L)*MOIST_DT
                 print *, "    AFTER UW ShallowCu"
                 print *, "  Latitude       =", LATS(I,J)*180.0/MAPL_PI
                 print *, "  Longitude      =", LONS(I,J)*180.0/MAPL_PI
                 print *, "  Pressure (mb)  =", 0.5*(PLE(I,J,L)+PLE(I,J,L-1))/100.0                  
                 print *, "                       CLLS=",  CLLS(I,J,L),             "CLCN=", CLCN(I,J,L)
                 print *, "  QV=",     Q(I,J,L), "  QL=",  QLLS(I,J,L)+QLCN(I,J,L), "  QI=", QILS(I,J,L)+QICN(I,J,L)
             endif               
           end do ! IM loop      
         end do ! JM loop        
       end do ! LM loop          
  endif

  call MAPL_TimerOff (MAPL,"--UW")

end subroutine UW_Run

end module GEOS_UW_InterfaceMod
