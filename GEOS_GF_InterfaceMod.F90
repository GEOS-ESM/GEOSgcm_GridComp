! $Id$

! #include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GF_InterfaceMod -- A Module to interface with the
!   GF convection

module GEOS_GF_InterfaceMod

!   use ESMF
!   use MAPL
  use GEOS_UtilsMod
  use GEOSmoist_Process_Library
  use Aer_Actv_Single_Moment
  use ConvPar_GF_SharedParams
!   use ConvPar_GF_GEOS5
  use ConvPar_GF2020

  implicit none

  private

!   character(len=ESMF_MAXSTR)              :: IAm
!   integer                                 :: STATUS

!   ! specify how to handle friendlies with DYN:TRB:CHM:ANA
!   type FRIENDLIES_TYPE
!          character(len=ESMF_MAXSTR) :: CNV_TR
!   end type FRIENDLIES_TYPE
!   type (FRIENDLIES_TYPE) FRIENDLIES

  integer :: USE_GF2020
  logical :: STOCHASTIC_CNV
  real    :: STOCH_TOP, STOCH_BOT
  real    :: SCLM_DEEP
  real    :: GF_MIN_AREA
  logical :: FIX_CNV_CLOUD

!   public :: GF_Setup, GF_Initialize, GF_Run
  public :: GF_Run

contains

! subroutine GF_Setup (GC, CF, RC)
!     type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
!     type(ESMF_Config),   intent(inout) :: CF
!     integer, optional                  :: RC  ! return code
!     character(len=ESMF_MAXSTR)         :: COMP_NAME

!     IAm = "GEOS_GF_InterfaceMod"
!     call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
!     VERIFY_(STATUS)
!     Iam = trim(COMP_NAME) // Iam

!     call ESMF_ConfigGetAttribute( CF ,CONVECTION_TRACER, Label='CONVECTION_TRACER:', default= 0, RC=STATUS );VERIFY_(STATUS)
!     if (CONVECTION_TRACER == 1) then
!        FRIENDLIES%CNV_TR   = trim(COMP_NAME)//"DYNAMICS"
!     else
!        FRIENDLIES%CNV_TR   = trim(COMP_NAME)
!     endif
!     call MAPL_AddInternalSpec(GC,                                  &
!          SHORT_NAME ='CNV_TR',                                     &
!          LONG_NAME  ='convection_tracer',                          &
!          UNITS      ='1',                                          &
!          FRIENDLYTO = trim(FRIENDLIES%CNV_TR),                     &
!          DIMS       = MAPL_DimsHorzVert,                           &
!          VLOCATION  = MAPL_VLocationCenter,                        &
!          RESTART    = MAPL_RestartSkip,                            &
!          DEFAULT    = 0.0,   RC=STATUS  )
!          VERIFY_(STATUS)

!     call MAPL_TimerAdd(GC, name="--GF", RC=STATUS)
!     VERIFY_(STATUS)

! end subroutine GF_Setup

! subroutine GF_Initialize (MAPL, CLOCK, RC)
!     type (MAPL_MetaComp), intent(inout) :: MAPL
!     type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!     integer, optional                   :: RC  ! return code
!     integer :: LM
!     type (ESMF_Alarm   )            :: ALARM
!     type (ESMF_TimeInterval)        :: TINT
!     real(ESMF_KIND_R8)              :: DT_R8
!     real                            :: MOIST_DT
!     real                            :: GF_DT

!     type(ESMF_Calendar)     :: calendar
!     type(ESMF_Time)         :: currentTime
!     type(ESMF_Alarm)        :: GF_RunAlarm
!     type(ESMF_Time)         :: ringTime
!     type(ESMF_TimeInterval) :: ringInterval
!     integer                 :: year, month, day, hh, mm, ss

!     call MAPL_Get(MAPL, RUNALARM=ALARM, LM=LM, RC=STATUS );VERIFY_(STATUS)
!     call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
!     call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)           
!     MOIST_DT = DT_R8
!     call MAPL_GetResource(MAPL, GF_DT, 'GF_DT:', default=MOIST_DT, RC=STATUS); VERIFY_(STATUS)

!     call ESMF_ClockGet(CLOCK, calendar=calendar, currTime=currentTime, RC=STATUS); VERIFY_(STATUS)
!     call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, RC=STATUS); VERIFY_(STATUS)
!     call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, RC=STATUS); VERIFY_(STATUS)
!     call ESMF_TimeIntervalSet(ringInterval, S=nint(GF_DT), calendar=calendar, RC=STATUS); VERIFY_(STATUS)
    
!     GF_RunAlarm = ESMF_AlarmCreate(Clock        = CLOCK,        &
!                                    Name         = 'GF_RunAlarm',&
!                                    RingInterval = ringInterval, &
!                                    RingTime     = currentTime,  &
!                                    Enabled      = .true.   ,    &
!                                    Sticky       = .false.  , RC=STATUS); VERIFY_(STATUS)

!     if (LM .eq. 72) then
!       call MAPL_GetResource(MAPL, USE_GF2020                , 'USE_GF2020:'            ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!     else
!       call MAPL_GetResource(MAPL, USE_GF2020                , 'USE_GF2020:'            ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!     endif
!     IF (USE_GF2020==1) THEN
!       call MAPL_GetResource(MAPL, ICUMULUS_GF(DEEP)         , 'DEEP:'                  ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, ICUMULUS_GF(SHAL)         , 'SHALLOW:'               ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, ICUMULUS_GF(MID)          , 'CONGESTUS:'             ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CLOSURE_CHOICE(DEEP)      , 'CLOSURE_DEEP:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CLOSURE_CHOICE(SHAL)      , 'CLOSURE_SHALLOW:'       ,default= 7,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CLOSURE_CHOICE(MID)       , 'CLOSURE_CONGESTUS:'     ,default= 3,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_ENTR_RATE(DEEP)       , 'ENTR_DP:'               ,default= 1.e-4,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_ENTR_RATE(SHAL)       , 'ENTR_SH:'               ,default= 2.e-3,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_ENTR_RATE(MID)        , 'ENTR_MD:'               ,default= 9.e-4,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(DEEP)    , 'FADJ_MASSFLX_DP:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(SHAL)    , 'FADJ_MASSFLX_SH:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(MID)     , 'FADJ_MASSFLX_MD:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_TRACER_TRANSP         , 'USE_TRACER_TRANSP:'     ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_TRACER_SCAVEN         , 'USE_TRACER_SCAVEN:'     ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_SCALE_DEP             , 'USE_SCALE_DEP:'         ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_MOMENTUM_TRANSP       , 'USE_MOMENTUM_TRANSP:'   ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, DICYCLE                   , 'DICYCLE:'               ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, OUTPUT_SOUND              , 'OUTPUT_SOUND:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_MEMORY                , 'USE_MEMORY:'            ,default=-1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, TAU_OCEA_CP               , 'TAU_OCEA_CP:'           ,default= 21600., RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, TAU_LAND_CP               , 'TAU_LAND_CP:'           ,default= 21600., RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, DOWNDRAFT                 , 'DOWNDRAFT:'             ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_FLUX_FORM             , 'USE_FLUX_FORM:'         ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_TRACER_EVAP           , 'USE_TRACER_EVAP:'       ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, APPLY_SUB_MP              , 'APPLY_SUB_MP:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, ALP1                      , 'ALP1:'                  ,default= 1.,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, LIGHTNING_DIAG            , 'LIGHTNING_DIAG:'        ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, OVERSHOOT                 , 'OVERSHOOT:'             ,default= 0.,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_WETBULB               , 'USE_WETBULB:'           ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, LAMBAU_SHDN               , 'LAMBAU_SHDN:'           ,default= 2.,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, MAX_TQ_TEND               , 'MAX_TQ_TEND:'           ,default= 100., RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_CLOUD_DISSIPATION     , 'USE_CLOUD_DISSIPATION:' ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_SMOOTH_TEND           , 'USE_SMOOTH_TEND:'       ,default= 1,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, BETA_SH                   , 'BETA_SH:'               ,default= 2.2,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_LINEAR_SUBCL_MF       , 'USE_LINEAR_SUBCL_MF:'   ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CAP_MAXS                  , 'CAP_MAXS:'              ,default= 50.,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, ZERO_DIFF                 , 'ZERO_DIFF:'             ,default= 0,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, GF_ENV_SETTING            , 'GF_ENV_SETTING:'        ,default= 'DYNAMICS', RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, STOCH_TOP                 , 'STOCH_TOP:'             ,default= 2.50,  RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, STOCH_BOT                 , 'STOCH_BOT:'             ,default= 0.75,  RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, STOCHASTIC_CNV            , 'STOCHASTIC_CNV:'        ,default= .FALSE.,RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, GF_MIN_AREA               , 'GF_MIN_AREA:'           ,default= -1.e6, RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, TAU_MID                   , 'TAU_MID:'               ,default= 5400., RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, TAU_DEEP                  , 'TAU_DEEP:'              ,default= 10800.,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CLEV_GRID                 , 'CLEV_GRID:'             ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, VERT_DISCR                , 'VERT_DISCR:'            ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_FCT                   , 'USE_FCT:'               ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, SATUR_CALC                , 'SATUR_CALC:'            ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, BC_METH                   , 'BC_METH:'               ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_REBCB                 , 'USE_REBCB:'             ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, AUTOCONV                  , 'AUTOCONV:'              ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, LAMBAU_DEEP               , 'LAMBAU_DEEP:'           ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL ,MOIST_TRIGGER             , 'MOIST_TRIGGER:'         ,default= 0,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, SGS_W_TIMESCALE           , 'SGS_W_TIMESCALE:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL ,FRAC_MODIS                , 'FRAC_MODIS:'            ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL ,USE_SMOOTH_PROF           , 'USE_SMOOTH_PROF:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL ,EVAP_FIX                  , 'EVAP_FIX:'              ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, ADV_TRIGGER               , 'ADV_TRIGGER:'           ,default= 0,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(DEEP)       , 'AVE_LAYER_DP:'          ,default= 40.,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(SHAL)       , 'AVE_LAYER_SH:'          ,default= 30.,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(MID)        , 'AVE_LAYER_MD:'          ,default= 40.,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(DEEP)      , 'USE_EXCESS_DP:'         ,default= 2,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(SHAL)      , 'USE_EXCESS_SH:'         ,default= 1,     RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(MID)       , 'USE_EXCESS_MD:'         ,default= 2,     RC=STATUS );VERIFY_(STATUS)
!       if (AUTOCONV == 1) then 
!          call MAPL_GetResource(MAPL, C0_DEEP                , 'C0_DEEP:'               ,default= 2.0e-3,RC=STATUS );VERIFY_(STATUS)
!          call MAPL_GetResource(MAPL, C0_MID                 , 'C0_MID:'                ,default= 2.0e-3,RC=STATUS );VERIFY_(STATUS)
!          call MAPL_GetResource(MAPL, C0_SHAL                , 'C0_SHAL:'               ,default= 0.    ,RC=STATUS );VERIFY_(STATUS)
!          call MAPL_GetResource(MAPL, QRC_CRIT               , 'QRC_CRIT:'              ,default= 2.0e-4,RC=STATUS );VERIFY_(STATUS)
!       else        
!          call MAPL_GetResource(MAPL, C0_DEEP                , 'C0_DEEP:'               ,default= 1.5e-3,RC=STATUS );VERIFY_(STATUS)
!          call MAPL_GetResource(MAPL, C0_MID                 , 'C0_MID:'                ,default= 1.0e-3,RC=STATUS );VERIFY_(STATUS)
!          call MAPL_GetResource(MAPL, C0_SHAL                , 'C0_SHAL:'               ,default= 1.0e-3,RC=STATUS );VERIFY_(STATUS)
!          call MAPL_GetResource(MAPL, QRC_CRIT               , 'QRC_CRIT:'              ,default= 1.0e-4,RC=STATUS );VERIFY_(STATUS)
!       endif
!       call MAPL_GetResource(MAPL, C1                        , 'C1:'                    ,default= 0.,    RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(DEEP)   , 'HEI_DOWN_LAND_DP:'      ,default= 0.3,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(SHAL)   , 'HEI_DOWN_LAND_SH:'      ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(MID)    , 'HEI_DOWN_LAND_MD:'      ,default= 0.3,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(DEEP)  , 'HEI_DOWN_OCEAN_DP:'     ,default= 0.6,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(SHAL)  , 'HEI_DOWN_OCEAN_SH:'     ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(MID)   , 'HEI_DOWN_OCEAN_MD:'     ,default= 0.6,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(DEEP)   , 'HEI_UPDF_LAND_DP:'      ,default= 0.4,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(SHAL)   , 'HEI_UPDF_LAND_SH:'      ,default= 0.2,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(MID)    , 'HEI_UPDF_LAND_MD:'      ,default= 0.4,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(DEEP)  , 'HEI_UPDF_OCEAN_DP:'     ,default= 0.4,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(SHAL)  , 'HEI_UPDF_OCEAN_SH:'     ,default= 0.2,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(MID)   , 'HEI_UPDF_OCEAN_MD:'     ,default= 0.4,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(DEEP)    , 'MAX_EDT_LAND_DP:'       ,default= 0.4,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(SHAL)    , 'MAX_EDT_LAND_SH:'       ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(MID)     , 'MAX_EDT_LAND_MD:'       ,default= 0.4,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(DEEP)   , 'MAX_EDT_OCEAN_DP:'      ,default= 0.3,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(SHAL)   , 'MAX_EDT_OCEAN_SH:'      ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(MID)    , 'MAX_EDT_OCEAN_MD:'      ,default= 0.3,   RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, SCLM_DEEP                 , 'SCLM_DEEP:'             ,default= 1.0,   RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, FIX_CNV_CLOUD             , 'FIX_CNV_CLOUD:'         ,default= .FALSE., RC=STATUS); VERIFY_(STATUS)
!     ELSE
!       ! Inititialize parameters of convection scheme GF (circa 2019)
!       call MAPL_GetResource(MAPL, icumulus_gf(deep)   ,'DEEP:'             ,default= 1 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, icumulus_gf(shal)   ,'SHALLOW:'          ,default= 0 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, icumulus_gf(mid)    ,'CONGESTUS:'        ,default= 1 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, closure_choice(deep),'CLOSURE_DEEP:'     ,default= 0 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, closure_choice(shal),'CLOSURE_SHALLOW:'  ,default= 7 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, closure_choice(mid) ,'CLOSURE_CONGESTUS:',default= 3 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_ENTR_RATE(DEEP) ,'ENTR_DP:'          ,default= 1.0e-4,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_ENTR_RATE(SHAL) ,'ENTR_SH:'          ,default= 1.4e-3,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, CUM_ENTR_RATE(MID)  ,'ENTR_MD:'          ,default= 9.0e-4,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, C1                  ,'C1:'               ,default= 1.0e-3,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_TRACER_TRANSP   ,'USE_TRACER_TRANSP:',default= 1 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_TRACER_SCAVEN   ,'USE_TRACER_SCAVEN:',default= 1 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_SCALE_DEP       ,'USE_SCALE_DEP:'    ,default= 1 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, DICYCLE             ,'DICYCLE:'          ,default= 1 , RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, STOCHASTIC_CNV      ,'STOCHASTIC_CNV:'   ,default= .FALSE., RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, GF_MIN_AREA         ,'GF_MIN_AREA:'      ,default= 1.e6,  RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, TAU_DEEP            ,'TAU_DEEP:'         ,default= 5400.0, RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, TAU_MID             ,'TAU_MID:'          ,default= 3600.0, RC=STATUS ); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_FLUX_FORM       ,'USE_FLUX_FORM:'    ,default= 1 , RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_FCT             ,'USE_FCT:'          ,default= 0 , RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, USE_TRACER_EVAP     ,'USE_TRACER_EVAP:'  ,default= 1 , RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, ALP1                ,'ALP1:'             ,default= 1., RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, BC_METH             ,'BC_METH:'          ,default= 0 , RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, C0_DEEP             ,'C0_DEEP:'          ,default= 2.0e-3,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, C0_MID              ,'C0_MID:'           ,default= 2.0e-3,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, C0_SHAL             ,'C0_SHAL:'          ,default= 1.0e-3,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, QRC_CRIT            ,'QRC_CRIT:'         ,default= 2.0e-4,RC=STATUS );VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, GF_ENV_SETTING      ,'GF_ENV_SETTING:'   ,default= 'DYNAMICS', RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, SCLM_DEEP           ,'SCLM_DEEP:'        ,default= 1.0    , RC=STATUS); VERIFY_(STATUS)
!       call MAPL_GetResource(MAPL, FIX_CNV_CLOUD       ,'FIX_CNV_CLOUD:'    ,default= .FALSE., RC=STATUS); VERIFY_(STATUS)
!     ENDIF
 
! end subroutine GF_Initialize


subroutine GF_Run (IM, JM, LM, dirName, rank_str)!(GC, IMPORT, EXPORT, CLOCK, RC)
   !  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
   !  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
   !  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
   !  type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
   !  integer, optional,   intent(  out) :: RC     ! Error code:

    ! Local derived type aliases

   !  type (MAPL_MetaComp), pointer   :: MAPL
   !  type (ESMF_Config  )            :: CF
   !  type (ESMF_State   )            :: INTERNAL
   !  type (ESMF_TimeInterval)        :: TINT
   !  real(ESMF_KIND_R8)              :: DT_R8
    real                            :: GF_DT
   !  type(ESMF_Alarm)                :: alarm
   !  logical                         :: alarm_is_ringing

   character(len=100), intent(in) :: dirName
   character(len=20), intent(in) :: rank_str

    ! Local variables
    integer                         :: I, J, L, fileID
    integer, intent(in)                         :: IM,JM,LM
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS
    real                            :: minrhx

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: CNV_TR ! tracer memory

    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W, KH, OMEGA
    real, pointer, dimension(:,:  ) :: FRLAND, AREA
    real, pointer, dimension(:,:,:) :: DQVDTDYN
    real, pointer, dimension(:,:,:) :: DTDTDYN
    real, pointer, dimension(:,:,:) :: QV_DYN_IN
    real, pointer, dimension(:,:,:) :: U_DYN_IN
    real, pointer, dimension(:,:,:) :: V_DYN_IN
    real, pointer, dimension(:,:,:) :: T_DYN_IN
    real, pointer, dimension(:,:,:) :: PLE_DYN_IN
    real, pointer, dimension(:,:,:) :: DTDT_BL
    real, pointer, dimension(:,:,:) :: DQDT_BL
    real, pointer, dimension(:,:  ) :: KPBL
    real, pointer, dimension(:,:,:) :: RADSW, RADLW

    ! allocatable derived quantities
    real,    allocatable, dimension(:,:,:) :: ZLE0
    real,    allocatable, dimension(:,:,:) :: PL, PK, ZL0
    real,    allocatable, dimension(:,:,:) :: MASS, fQi, QST3
    real,    allocatable, dimension(:,:,:) :: TH
    real,    allocatable, dimension(:,:,:) :: REVSU, PRFIL
    integer, allocatable, dimension(:,:)   :: SEEDINI   
    real,    allocatable, dimension(:,:)   :: SEEDCNV 
    real,    allocatable, dimension(:,:,:) :: TMP3D
    real,    allocatable, dimension(:,:)   :: TMP2D

    ! Required Exports (connectivities to moist siblings)
    real, pointer, dimension(:,:,:) :: CNV_MFD, CNV_MFC, CNV_CVW, CNV_QC, CNV_DQCDT, CNV_PRC3, CNV_UPDF
    real, pointer, dimension(:,:,:) :: DUDT_DC, DVDT_DC, DTDT_DC, DQVDT_DC, DQIDT_DC, DQLDT_DC, DQADT_DC
    real, pointer, dimension(:,:  ) :: CNV_FRC, SRF_TYPE
    ! Exports
    real, pointer, dimension(:,:,:) :: CNV_MF0, ENTLAM
    real, pointer, dimension(:,:,:) :: MUPDP,MDNDP,MUPSH,MUPMD
    real, pointer, dimension(:,:,:) :: VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d
    real, pointer, dimension(:,:  ) :: T2M,Q2M,TA,QA,SH,EVAP,PHIS
    real, pointer, dimension(:,:  ) :: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
    real, pointer, dimension(:,:  ) :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC
    real, pointer, dimension(:,:  ) :: TPWI,TPWI_star,LFR_GF,CNPCPRATE
    real, pointer, dimension(:,:,:) :: RSU_CN,REV_CN,PFL_CN,PFI_CN
    real, pointer, dimension(:,:  ) :: SIGMA_DEEP, SIGMA_MID
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D

    ! CK : Extra variables
   integer :: CNV_tracer_size
   character(2) :: n_str

   print*,'In GF_RUN'

!     call ESMF_ClockGetAlarm(clock, 'GF_RunAlarm', alarm, RC=STATUS); VERIFY_(STATUS)
!     alarm_is_ringing = ESMF_AlarmIsRinging(alarm, RC=STATUS); VERIFY_(STATUS)

!     if (alarm_is_ringing) then

! !!! call WRITE_PARALLEL('GF is Running')
!     call ESMF_AlarmRingerOff(alarm, RC=STATUS); VERIFY_(STATUS)
!     call ESMF_AlarmGet(alarm, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
!     call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
!       GF_DT = DT_R8

!     call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ); VERIFY_(STATUS)

!     ! Get my internal MAPL_Generic state
!     !-----------------------------------

!     call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS); VERIFY_(STATUS)

!     call MAPL_TimerOn (MAPL,"--GF")

!     ! Get parameters from generic state.
!     !-----------------------------------

!     call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
!          CF       = CF,                &
!          LONS     = LONS,              &
!          LATS     = LATS,              &
!          INTERNAL_ESMF_STATE=INTERNAL, &
!          RC=STATUS )
!     VERIFY_(STATUS)

   !  ! Internals
   !  call MAPL_GetPointer(INTERNAL, Q,      'Q'       , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, QLLS,   'QLLS'    , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, QLCN,   'QLCN'    , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, CLCN,   'CLCN'    , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, CLLS,   'CLLS'    , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, QILS,   'QILS'    , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, QICN,   'QICN'    , RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(INTERNAL, CNV_TR, 'CNV_TR'  , RC=STATUS); VERIFY_(STATUS)

   !  ! Imports
   !  call MAPL_GetPointer(IMPORT, FRLAND    ,'FRLAND'    ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, AREA      ,'AREA'      ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, ZLE       ,'ZLE'       ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, PLE       ,'PLE'       ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, T         ,'T'         ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, U         ,'U'         ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, V         ,'V'         ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, W         ,'W'         ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, OMEGA     ,'OMEGA'     ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, KH        ,'KH'        ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, SH        ,'SH'        ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, EVAP      ,'EVAP'      ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, T2M       ,'T2M'       ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, Q2M       ,'Q2M'       ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, TA        ,'TA'        ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, QA        ,'QA'        ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, PHIS      ,'PHIS'      ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, DQVDTDYN  ,'DQVDTDYN'  ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, DTDTDYN   ,'DTDTDYN'   ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, QV_DYN_IN ,'QV_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, U_DYN_IN  ,'U_DYN_IN'  ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, V_DYN_IN  ,'V_DYN_IN'  ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, T_DYN_IN  ,'T_DYN_IN'  ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, PLE_DYN_IN,'PLE_DYN_IN',RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, DTDT_BL   ,'DTDT_BL'   ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, DQDT_BL   ,'DQDT_BL'   ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, RADSW     ,'RADSW'     ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, RADLW     ,'RADLW'     ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(IMPORT, KPBL      ,'KPBL'      ,RC=STATUS); VERIFY_(STATUS)

   open(newunit=fileID, file=trim(dirName) // '/USE_GF2020_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_GF2020
   close(fileID)

   open(newunit=fileID, file=trim(dirName) // '/GF_DT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) GF_DT
   close(fileID)

   allocate(LONS(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/LONS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) LONS
   close(fileID)

   allocate(LATS(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/LATS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) LATS
   close(fileID)

   open(newunit=fileID, file=trim(dirName) // '/ICUMULUS_GF_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ICUMULUS_GF(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/ICUMULUS_GF_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ICUMULUS_GF(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/ICUMULUS_GF_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ICUMULUS_GF(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CLOSURE_CHOICE_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CLOSURE_CHOICE(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CLOSURE_CHOICE_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CLOSURE_CHOICE(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CLOSURE_CHOICE_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CLOSURE_CHOICE(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_ENTR_RATE_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_ENTR_RATE(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_ENTR_RATE_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_ENTR_RATE(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_ENTR_RATE_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_ENTR_RATE(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_FADJ_MASSFLX_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_FADJ_MASSFLX(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_FADJ_MASSFLX_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_FADJ_MASSFLX(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_FADJ_MASSFLX_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_FADJ_MASSFLX(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_TRACER_TRANSP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_TRACER_TRANSP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_TRACER_SCAVEN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_TRACER_SCAVEN
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_SCALE_DEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_SCALE_DEP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_MOMENTUM_TRANSP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_MOMENTUM_TRANSP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/DICYCLE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) DICYCLE
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/OUTPUT_SOUND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) OUTPUT_SOUND
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_MEMORY_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_MEMORY
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/TAU_OCEA_CP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) TAU_OCEA_CP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/TAU_LAND_CP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) TAU_LAND_CP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/DOWNDRAFT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) DOWNDRAFT
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_FLUX_FORM_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_FLUX_FORM
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_TRACER_EVAP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_TRACER_EVAP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/APPLY_SUB_MP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) APPLY_SUB_MP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/ALP1_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ALP1
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/LIGHTNING_DIAG_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) LIGHTNING_DIAG
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/OVERSHOOT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) OVERSHOOT
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_WETBULB_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_WETBULB
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/LAMBAU_SHDN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) LAMBAU_SHDN
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/MAX_TQ_TEND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) MAX_TQ_TEND
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_CLOUD_DISSIPATION_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_CLOUD_DISSIPATION
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_SMOOTH_TEND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_SMOOTH_TEND
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/BETA_SH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) BETA_SH
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_LINEAR_SUBCL_MF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_LINEAR_SUBCL_MF
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CAP_MAXS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CAP_MAXS
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/ZERO_DIFF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ZERO_DIFF
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/GF_ENV_SETTING_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) GF_ENV_SETTING
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/STOCH_TOP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) STOCH_TOP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/STOCH_BOT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) STOCH_BOT
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/STOCHASTIC_CNV_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) STOCHASTIC_CNV
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/GF_MIN_AREA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) GF_MIN_AREA
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/TAU_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) TAU_MID
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/TAU_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) TAU_DEEP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CLEV_GRID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CLEV_GRID
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/VERT_DISCR_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) VERT_DISCR
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_FCT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_FCT
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/SATUR_CALC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) SATUR_CALC
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/BC_METH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) BC_METH
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_REBCB_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_REBCB
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/AUTOCONV_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) AUTOCONV
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/LAMBAU_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) LAMBAU_DEEP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/MOIST_TRIGGER_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) MOIST_TRIGGER
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/SGS_W_TIMESCALE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) SGS_W_TIMESCALE
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/FRAC_MODIS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) FRAC_MODIS
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/USE_SMOOTH_PROF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) USE_SMOOTH_PROF
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/EVAP_FIX_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) EVAP_FIX
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/ADV_TRIGGER_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ADV_TRIGGER
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_AVE_LAYER_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_AVE_LAYER(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_AVE_LAYER_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_AVE_LAYER(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_AVE_LAYER_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_AVE_LAYER(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_USE_EXCESS_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_USE_EXCESS(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_USE_EXCESS_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_USE_EXCESS(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_USE_EXCESS_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_USE_EXCESS(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/C0_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) C0_DEEP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/C0_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) C0_MID
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/C0_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) C0_SHAL
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/QRC_CRIT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QRC_CRIT
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/C1_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) C1
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_DOWN_LAND_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_DOWN_LAND(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_DOWN_LAND_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_DOWN_LAND(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_DOWN_LAND_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_DOWN_LAND(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_DOWN_OCEAN_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_DOWN_OCEAN(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_DOWN_OCEAN_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_DOWN_OCEAN(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_DOWN_OCEAN_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_DOWN_OCEAN(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_UPDF_LAND_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_UPDF_LAND(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_UPDF_LAND_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_UPDF_LAND(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_UPDF_LAND_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_UPDF_LAND(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_UPDF_OCEAN_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_UPDF_OCEAN(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_UPDF_OCEAN_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_UPDF_OCEAN(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_HEI_UPDF_OCEAN_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_HEI_UPDF_OCEAN(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_MAX_EDT_LAND_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_MAX_EDT_LAND(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_MAX_EDT_LAND_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_MAX_EDT_LAND(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_MAX_EDT_LAND_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_MAX_EDT_LAND(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_MAX_EDT_OCEAN_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_MAX_EDT_OCEAN(DEEP)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_MAX_EDT_OCEAN_SHAL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_MAX_EDT_OCEAN(SHAL)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/CUM_MAX_EDT_OCEAN_MID_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CUM_MAX_EDT_OCEAN(MID)
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/SCLM_DEEP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) SCLM_DEEP
   close(fileID)
   
   open(newunit=fileID, file=trim(dirName) // '/FIX_CNV_CLOUD_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) FIX_CNV_CLOUD
   close(fileID)

   allocate(Q(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) Q
   close(fileID)
   
   allocate(QLLS(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QLLS
   close(fileID)
   
   allocate(QLCN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QLCN
   close(fileID)
   
   allocate(CLLS(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CLLS
   close(fileID)
   
   allocate(CLCN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CLCN
   close(fileID)
   
   allocate(QILS(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QILS
   close(fileID)

   allocate(QICN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QICN
   close(fileID)
   
   allocate(CNV_TR(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/CNV_TR_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CNV_TR
   close(fileID)
   
   allocate(ZLE(IM, JM, 0:LM))
   open(newunit=fileID, file=trim(dirName) // '/ZLE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) ZLE
   close(fileID)
   
   allocate(PLE(IM, JM, 0:LM))
   open(newunit=fileID, file=trim(dirName) // '/PLE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) PLE
   close(fileID)
   
   allocate(T(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) T
   close(fileID)
   
   allocate(U(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/U_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) U
   close(fileID)
   
   allocate(V(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/V_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) V
   close(fileID)
   
   allocate(W(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/W_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) W
   close(fileID)
   
   allocate(KH(IM, JM, 0:LM))
   open(newunit=fileID, file=trim(dirName) // '/KH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) KH
   close(fileID)

   allocate(SH(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/SH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) SH
   close(fileID)

   allocate(EVAP(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/EVAP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) EVAP
   close(fileID)
   
   allocate(T2M(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/T2M_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) T2M
   close(fileID)

   allocate(Q2M(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/Q2M_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) Q2M
   close(fileID)

   allocate(TA(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/TA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) TA
   close(fileID)

   allocate(QA(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/QA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QA
   close(fileID)

   allocate(PHIS(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/PHIS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) PHIS
   close(fileID)

   allocate(OMEGA(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/OMEGA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) OMEGA
   close(fileID)
   
   allocate(FRLAND(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/FRLAND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) FRLAND
   close(fileID)
   
   allocate(AREA(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/AREA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) AREA
   close(fileID)
   
   allocate(DQVDTDYN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/DQVDTDYN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) DQVDTDYN
   close(fileID)

   allocate(DTDTDYN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/DTDTDYN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) DTDTDYN
   close(fileID)
   
   allocate(QV_DYN_IN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/QV_DYN_IN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) QV_DYN_IN
   close(fileID)
   
   allocate(U_DYN_IN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/U_DYN_IN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) U_DYN_IN
   close(fileID)
   
   allocate(V_DYN_IN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/V_DYN_IN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) V_DYN_IN
   close(fileID)
   
   allocate(T_DYN_IN(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/T_DYN_IN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) T_DYN_IN
   close(fileID)
   
   allocate(PLE_DYN_IN(IM, JM, 0:LM))
   open(newunit=fileID, file=trim(dirName) // '/PLE_DYN_IN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) PLE_DYN_IN
   close(fileID)
   
   allocate(DTDT_BL(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/DTDT_BL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) DTDT_BL
   close(fileID)
   
   allocate(DQDT_BL(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/DQDT_BL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) DQDT_BL
   close(fileID)
   
   allocate(KPBL(IM, JM))
   open(newunit=fileID, file=trim(dirName) // '/KPBL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) KPBL
   close(fileID)
   
   allocate(RADSW(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/RADSW_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) RADSW
   close(fileID)
   
   allocate(RADLW(IM, JM, LM))
   open(newunit=fileID, file=trim(dirName) // '/RADLW_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) RADLW
   close(fileID)

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PRFIL(IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PL   (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( TH   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( fQi  (IM,JM,LM  ) )
    ALLOCATE ( QST3 (IM,JM,LM  ) )
    ALLOCATE ( REVSU(IM,JM,LM  ) )
    ALLOCATE ( TMP3D(IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( SEEDINI(IM,JM) )
    ALLOCATE ( SEEDCNV(IM,JM) )
    ALLOCATE ( TMP2D  (IM,JM) )

    ! derived quantaties
    ! Derived States
    PL       = 0.5*(PLE(:,:,0:LM-1) + PLE(:,:,1:LM))
    PK       = (PL/MAPL_P00)**(MAPL_KAPPA)
    DO L=0,LM
       ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
    END DO
    ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
    TH       = T/PK
    MASS     = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )/MAPL_GRAV

   !  ! Required Exports (connectivities to moist siblings)
   !  call MAPL_GetPointer(EXPORT, CNV_MFD,    'CNV_MFD'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_MFC,    'CNV_MFC'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_CVW,    'CNV_CVW'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_QC,     'CNV_QC'    ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_DQCDT,  'CNV_DQCDT' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_PRC3,   'CNV_PRC3'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_UPDF,   'CNV_UPDF'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNV_FRC,    'CNV_FRC'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, SRF_TYPE,   'SRF_TYPE'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  ! Exports used below
   !  call MAPL_GetPointer(EXPORT, CNV_MF0,    'CNV_MF0'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, CNPCPRATE,  'CNPCPRATE' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, ENTLAM,     'ENTLAM'    ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, SIGMA_DEEP, 'SIGMA_DEEP',  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, SIGMA_MID,  'SIGMA_MID' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, TPWI,       'TPWI'      ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, TPWI_star,  'TPWI_star' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)

   !  call MAPL_GetPointer(EXPORT, LFR_GF   ,'LFR_GF'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MUPDP    ,'MUPDP'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MDNDP    ,'MDNDP'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MUPSH    ,'MUPSH'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MUPMD    ,'MUPMD'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, VAR3d_a  ,'VAR3d_a'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, VAR3d_b  ,'VAR3d_b'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, VAR3d_c  ,'VAR3d_c'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, VAR3d_d  ,'VAR3d_d'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MFDP     ,'MFDP'      ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MFSH     ,'MFSH'      ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, MFMD     ,'MFMD'      ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, ERRDP    ,'ERRDP'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, ERRSH    ,'ERRSH'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, ERRMD    ,'ERRMD'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, RSU_CN   ,'RSU_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, REV_CN   ,'REV_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, PFI_CN   ,'PFI_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, PFL_CN   ,'PFL_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, AA0      ,'AA0'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, AA1      ,'AA1'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, AA2      ,'AA2'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, AA3      ,'AA3'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, AA1_BL   ,'AA1_BL'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, AA1_CIN  ,'AA1_CIN'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, TAU_BL   ,'TAU_BL'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, TAU_EC   ,'TAU_EC'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)

   !  ! Initialize tendencies
   !  call MAPL_GetPointer(EXPORT,  DUDT_DC,    'DUDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT,  DVDT_DC,    'DVDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT,  DTDT_DC,    'DTDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, DQVDT_DC,   'DQVDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, DQLDT_DC,   'DQLDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, DQIDT_DC,   'DQIDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
   !  call MAPL_GetPointer(EXPORT, DQADT_DC,   'DQADT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)

   allocate(CNV_MFD(IM, JM, LM))
   allocate(CNV_MFC(IM, JM, 0:LM))
   allocate(CNV_CVW(IM, JM, LM))
   allocate(CNV_QC(IM, JM, LM))
   allocate(CNV_DQCDT(IM, JM, LM))
   allocate(CNV_PRC3(IM, JM, LM))
   allocate(CNV_UPDF(IM, JM, LM))
   allocate(DUDT_DC(IM, JM, LM))
   allocate(DVDT_DC(IM, JM, LM))
   allocate(DTDT_DC(IM, JM, LM))
   allocate(DQVDT_DC(IM, JM, LM))
   allocate(DQIDT_DC(IM, JM, LM))
   allocate(DQLDT_DC(IM, JM, LM))
   allocate(DQADT_DC(IM, JM, LM))
   allocate(CNV_FRC(IM, JM))
   allocate(SRF_TYPE(IM, JM))
   allocate(CNV_MF0(IM, JM, LM))
   allocate(ENTLAM(IM, JM, LM))
   allocate(MUPDP(IM, JM, LM))
   allocate(MDNDP(IM, JM, LM))
   allocate(MUPSH(IM, JM, LM))
   allocate(MUPMD(IM, JM, LM))
   allocate(VAR3d_a(IM, JM, LM))
   allocate(VAR3d_b(IM, JM, LM))
   allocate(VAR3d_c(IM, JM, LM))
   allocate(VAR3d_d(IM, JM, LM))
   allocate(MFDP(IM, JM))
   allocate(MFSH(IM, JM))
   allocate(MFMD(IM, JM))
   allocate(ERRDP(IM, JM))
   allocate(ERRSH(IM, JM))
   allocate(ERRMD(IM, JM))
   allocate(AA0(IM, JM))
   allocate(AA1(IM, JM))
   allocate(AA2(IM, JM))
   allocate(AA3(IM, JM))
   allocate(AA1_BL(IM, JM))
   allocate(AA1_CIN(IM, JM))
   allocate(TAU_BL(IM, JM))
   allocate(TAU_EC(IM, JM))
   allocate(TPWI(IM, JM))
   allocate(TPWI_star(IM, JM))
   allocate(LFR_GF(IM, JM))
   allocate(CNPCPRATE(IM, JM))
   allocate(RSU_CN(IM, JM, LM))
   allocate(REV_CN(IM, JM, LM))
   allocate(PFL_CN(IM, JM, 0:LM))
   allocate(PFI_CN(IM, JM, 0:LM))
   allocate(SIGMA_DEEP(IM, JM))
   allocate(SIGMA_MID(IM, JM))
   allocate(PTR3D(IM, JM, LM))
   allocate(PTR2D(IM, JM))

   ! Note : CNV_FRC and SRF_TYPE, even though these are designated as "exports", they are used at inputs to calcuations below
   open(newunit=fileID, file=trim(dirName) // '/CNV_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CNV_FRC
   close(fileID)

   open(newunit=fileID, file=trim(dirName) // '/SRF_TYPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) SRF_TYPE
   close(fileID)

   ! Note: CNV_Tracers is typically set in GEOS_moistGridComp.  This is to enable testing of GF_Run individually
   ! Read in Tracer data
   open(newunit=fileID, file=trim(dirName) // '/CNV_Tracers_Size_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
   read(fileID) CNV_tracer_size
   close(fileID)

   allocate(CNV_Tracers(CNV_tracer_size))
   
   do I = 1,CNV_tracer_size
      allocate(CNV_Tracers(I)%Q(IM, JM, LM))

      if(I>10) then
         write(n_str,'(i2)') I
      else
         write(n_str,'(i1)') I
      endif

      open(newunit=fileID, file=trim(dirName) // '/CNV_Tracers_' // trim(n_str) // '_Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
      read(fileID) CNV_tracers(I)%Q
      close(fileID)

      open(newunit=fileID, file=trim(dirName) // '/CNV_Tracers_' // trim(n_str) // '_Vect_Hcts_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
      read(fileID) CNV_tracers(I)%Vect_Hcts
      close(fileID)

      open(newunit=fileID, file=trim(dirName) // '/CNV_Tracers_' // trim(n_str) // '_fscav_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
      read(fileID) CNV_tracers(I)%fscav
      close(fileID)
   enddo

   TPWI       = SUM( Q*MASS, 3 )
   TPWI_star  = SUM( GEOS_QSAT(T,PL,PASCALS=.true.)*MASS, 3 )

    if (STOCHASTIC_CNV) then
       ! Create bit-processor-reproducible random white noise for convection [0:1]
       SEEDINI = 1000000 * ( 100*T(:,:,LM)   - INT( 100*T(:,:,LM) ) )
       SEEDCNV = SQRT(MAX(MIN(SEEDINI/1000000.0,1.0),0.0))
       ! Create stochastic variability to GF sigma
       SEEDCNV = SQRT(1.0-(1.0-SEEDCNV))*(STOCH_TOP-STOCH_BOT)+STOCH_BOT
    else
       SEEDCNV = 1.0
    endif
   !  CALL MAPL_GetPointer(EXPORT, PTR2D,  'STOCH_CNV', RC=STATUS); VERIFY_(STATUS)
    if(associated(PTR2D)) PTR2D = SEEDCNV

    ! Modify AREA (m^2) here so GF scale dependence has a CNV_FRC dependence
    if (GF_MIN_AREA > 0) then
       where (AREA > GF_MIN_AREA)
          TMP2D = GF_MIN_AREA*CNV_FRC + AREA*(1.0-CNV_FRC)
       elsewhere
          TMP2D =  AREA
       endwhere
    else if (GF_MIN_AREA < 0) then
       where (AREA > ABS(GF_MIN_AREA)) 
          TMP2D = AREA*CNV_FRC + ABS(GF_MIN_AREA)*(1.0-CNV_FRC)
       elsewhere
          TMP2D =  AREA
       endwhere
    else
       TMP2D = AREA
    endif

    call read_convpar_gf2020_module_var(dirName, rank_str)

    IF (USE_GF2020==1) THEN
         !- call GF2020 interface routine
         ! PLE and PL are passed in Pq
         call GF2020_Interface(   IM,JM,LM,LONS,LATS,GF_DT                       &
                                 ,PLE, PL, ZLE0, ZL0, PK, MASS, OMEGA, KH           &
                                 ,T, TH, Q, U, V, QLCN, QICN, QLLS, QILS, CNPCPRATE &
                                 ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQCDT, ENTLAM     &
                                 ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC, CLCN, CLLS    &
                                 ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                                 ,RADSW   ,RADLW  ,DQDT_BL  ,DTDT_BL                &
                                 ,FRLAND, TMP2D, T2M           &
                                 ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS                       &
                                 ,KPBL ,CNV_FRC, SRF_TYPE                           &
                                 ,SEEDCNV, SIGMA_DEEP, SIGMA_MID                    &
                                 ,DQVDT_DC,DTDT_DC,DUDT_DC,DVDT_DC                  &
                                 ,MUPDP,MUPSH,MUPMD,MDNDP                           &
                                 ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
                                 ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
                                 ,DTDTDYN,DQVDTDYN                                  &
                                 ,REVSU, PRFIL                                      &
                                 ,TPWI, TPWI_star, LFR_GF                           &
                                 ,VAR3d_a, VAR3d_b, VAR3d_c, VAR3d_d, CNV_TR)
    ELSE
         ! !- call GF/GEOS5 interface routine
         ! ! PLE and PL are passed in Pq
         ! call GF_GEOS5_Interface( IM,JM,LM,LONS,LATS,GF_DT                       &
         !                         ,PLE, PL, ZLE0, ZL0, PK, MASS, OMEGA               &
         !                         ,T, TH, Q, U, V, QLCN, QICN, QLLS, QILS, CNPCPRATE &
         !                         ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQCDT,ENTLAM      &
         !                         ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC , CLCN         &
         !                         ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
         !                         ,RADSW   ,RADLW  ,DQDT_BL  ,DTDT_BL                &
         !                         ,FRLAND, TMP2D, T2M, Q2M      &
         !                         ,TA ,QA ,SH ,EVAP ,PHIS                            &  
         !                         ,KPBL ,CNV_FRC, SRF_TYPE                           &
         !                         ,SEEDCNV, SIGMA_DEEP, SIGMA_MID                    &
         !                         ,DQVDT_DC,DTDT_DC,DUDT_DC,DVDT_DC                  &
         !                         ,MUPDP,MUPSH,MUPMD                                 &
         !                         ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
         !                         ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
         !                         ,DTDTDYN,DQVDTDYN                                  &
         !                         ,REVSU, PRFIL)
    ENDIF

    ! add tendencies to the moist import state
      U  = U  +  DUDT_DC*GF_DT
      V  = V  +  DVDT_DC*GF_DT
      Q  = Q  + DQVDT_DC*GF_DT
      T  = T  +  DTDT_DC*GF_DT
      TH = T/PK
    ! update DeepCu QL/QI/CF tendencies
      fQi = ice_fraction( T, CNV_FRC, SRF_TYPE )
      TMP3D    = CNV_DQCDT/MASS
      DQLDT_DC = (1.0-fQi)*TMP3D
      DQIDT_DC =      fQi *TMP3D
      DQADT_DC = CNV_MFD*SCLM_DEEP/MASS
    ! evap/subl and precip fluxes
      do L=1,LM
         !--- sublimation/evaporation tendencies (kg/kg/s)
           RSU_CN (:,:,L) = REVSU(:,:,L)*     fQi(:,:,L)
           REV_CN (:,:,L) = REVSU(:,:,L)*(1.0-fQi(:,:,L))
         !--- preciptation fluxes (kg/kg/s)
           PFI_CN (:,:,L) = PRFIL(:,:,L)*     fQi(:,:,L)
           PFL_CN (:,:,L) = PRFIL(:,:,L)*(1.0-fQi(:,:,L))
      enddo
    ! add QI/QL/CL tendencies
      QLCN =         QLCN + DQLDT_DC*GF_DT
      QICN =         QICN + DQIDT_DC*GF_DT
      CLCN = MAX(MIN(CLCN + DQADT_DC*GF_DT, 1.0), 0.0)
    ! Export
      ! call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_FICE', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) PTR3D = fQi
    ! fix 'convective' cloud fraction 
      if (FIX_CNV_CLOUD) then
      ! fix convective cloud
      TMP3D = GEOS_DQSAT(T, PL, PASCALS=.true., QSAT=QST3)
      TMP3D = QST3
      WHERE (CLCN < 1.0)
         TMP3D = ( Q - QST3 * CLCN )/(1.-CLCN)
      END WHERE
      minrhx = 0.001
      WHERE ( (( TMP3D - minrhx*QST3 ) < 0.0 ) .AND. (CLCN > 0.0) )
         CLCN = (Q  - minrhx*QST3 )/( QST3*(1.0-minrhx) )
      END WHERE
      ! If still cant make suitable env RH then destroy anvil
      WHERE ( CLCN < 0.0 )
         CLCN = 0.
         DQLDT_DC = DQLDT_DC - (QLCN       )/GF_DT
         DQIDT_DC = DQIDT_DC - (       QICN)/GF_DT
         DQVDT_DC = DQVDT_DC + (QLCN + QICN)/GF_DT
          Q       =  Q       + (QLCN + QICN)
         TMP3D    = (MAPL_ALHL*QLCN + MAPL_ALHS*QICN)/MAPL_CP
         DTDT_DC  = DTDT_DC  - TMP3D/GF_DT
          T       =  T       - TMP3D
          TH      =  T/PK
         QLCN = 0.
         QICN = 0.
      END WHERE
      endif

      ! call MAPL_GetPointer(EXPORT, PTR3D, 'DQRC', RC=STATUS); VERIFY_(STATUS)
      if(associated(PTR3D)) PTR3D = CNV_PRC3 / GF_DT

   !  call MAPL_TimerOff (MAPL,"--GF")

   !  endif

end subroutine GF_Run

end module GEOS_GF_InterfaceMod

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.