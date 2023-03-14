! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GF_InterfaceMod -- A Module to interface with the
!   GF convection

module GEOS_GF_InterfaceMod

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use GEOSmoist_Process_Library
  use Aer_Actv_Single_Moment
  use ConvPar_GF_GEOS5

  implicit none

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  ! specify how to handle friendlies with DYN:TRB:CHM:ANA
  type FRIENDLIES_TYPE
         character(len=ESMF_MAXSTR) :: CNV_TR
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  integer :: USE_GF2020
  logical :: STOCHASTIC_CNV
  real    :: STOCH_TOP, STOCH_BOT
  real    :: SCLM_DEEP
  real    :: GF_MIN_AREA
  logical :: FIX_CNV_CLOUD

  public :: GF_Setup, GF_Initialize, GF_Run

contains

subroutine GF_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_GF_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call ESMF_ConfigGetAttribute( CF ,CONVECTION_TRACER, Label='CONVECTION_TRACER:', default= 0, RC=STATUS )
    VERIFY_(STATUS)
    if (CONVECTION_TRACER == 1) then
       FRIENDLIES%CNV_TR   = trim(COMP_NAME)//"DYNAMICS"
    else
       FRIENDLIES%CNV_TR   = trim(COMP_NAME)
    endif

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='CNV_TR',                                     &
         LONG_NAME  ='convection_tracer',                          &
         UNITS      ='1',                                          &
         FRIENDLYTO = trim(FRIENDLIES%CNV_TR),                     &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         RESTART    = MAPL_RestartSkip,                            &
         DEFAULT    = 0.0,   RC=STATUS  )
         VERIFY_(STATUS)

    call MAPL_TimerAdd(GC, name="--GF", RC=STATUS)
    VERIFY_(STATUS)

end subroutine GF_Setup

subroutine GF_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    call MAPL_GetResource(MAPL, ICUMULUS_GF(DEEP)           ,'DEEP:'                  ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, ICUMULUS_GF(SHAL)           ,'SHALLOW:'               ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, ICUMULUS_GF(MID)            ,'CONGESTUS:'             ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CLOSURE_CHOICE(DEEP)        ,'CLOSURE_DEEP:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CLOSURE_CHOICE(SHAL)        ,'CLOSURE_SHALLOW:'       ,default= 7,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CLOSURE_CHOICE(MID)         ,'CLOSURE_CONGESTUS:'     ,default= 3,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_ENTR_RATE(DEEP)         ,'ENTR_DP:'               ,default= 1.e-4,RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_ENTR_RATE(SHAL)         ,'ENTR_SH:'               ,default= 2.e-3,RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_ENTR_RATE(MID)          ,'ENTR_MD:'               ,default= 5.e-4,RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(DEEP)      ,'FADJ_MASSFLX_DP:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(SHAL)      ,'FADJ_MASSFLX_SH:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(MID)       ,'FADJ_MASSFLX_MD:'       ,default= 0.5,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_TRACER_TRANSP           ,'USE_TRACER_TRANSP:'     ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_TRACER_SCAVEN           ,'USE_TRACER_SCAVEN:'     ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_SCALE_DEP               ,'USE_SCALE_DEP:'         ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_MOMENTUM_TRANSP         ,'USE_MOMENTUM_TRANSP:'   ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, DICYCLE                     ,'DICYCLE:'               ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, OUTPUT_SOUND                ,'OUTPUT_SOUND:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_MEMORY                  ,'USE_MEMORY:'            ,default=-1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, TAU_OCEA_CP                 ,'TAU_OCEA_CP:'           ,default= 21600., RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, TAU_LAND_CP                 ,'TAU_LAND_CP:'           ,default= 21600., RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, DOWNDRAFT                   , 'DOWNDRAFT:'            ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_FLUX_FORM               , 'USE_FLUX_FORM:'        ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_TRACER_EVAP             , 'USE_TRACER_EVAP:'      ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, APPLY_SUB_MP                , 'APPLY_SUB_MP:'         ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, ALP1                        , 'ALP1:'                 ,default= 1.,   RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, LIGHTNING_DIAG              , 'LIGHTNING_DIAG:'       ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, OVERSHOOT                   , 'OVERSHOOT:'            ,default= 0.,   RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_WETBULB                 , 'USE_WETBULB:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, LAMBAU_SHDN                 , 'LAMBAU_SHDN:'          ,default= 2.,   RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, MAX_TQ_TEND                 , 'MAX_TQ_TEND:'          ,default= 100., RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_CLOUD_DISSIPATION       , 'USE_CLOUD_DISSIPATION:',default= 1.0,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_SMOOTH_TEND             , 'USE_SMOOTH_TEND:'      ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, BETA_SH                     , 'BETA_SH:'              ,default= 2.2,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_LINEAR_SUBCL_MF         , 'USE_LINEAR_SUBCL_MF:'  ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CAP_MAXS                    , 'CAP_MAXS:'             ,default= 50.,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_GF2020                  , 'USE_GF2020:'           ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    IF(USE_GF2020==1) THEN
       call MAPL_GetResource(MAPL, ZERO_DIFF                 , 'ZERO_DIFF:'           ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, GF_ENV_SETTING            , 'GF_ENV_SETTING:'      ,default= 'CURRENT', RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, STOCH_TOP                 , 'STOCH_TOP:'           ,default= 1.25,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, STOCH_BOT                 , 'STOCH_BOT:'           ,default= 0.75,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, STOCHASTIC_CNV            , 'STOCHASTIC_CNV:'      ,default= .TRUE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, GF_MIN_AREA               , 'GF_MIN_AREA:'         ,default= 36.e6, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, TAU_MID                   , 'TAU_MID:'             ,default= 5400., RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, TAU_DEEP                  , 'TAU_DEEP:'            ,default= 10800.,RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CLEV_GRID                 , 'CLEV_GRID:'           ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, VERT_DISCR                , 'VERT_DISCR:'          ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, USE_FCT                   , 'USE_FCT:'             ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, SATUR_CALC                , 'SATUR_CALC:'          ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, BC_METH                   , 'BC_METH:'             ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, USE_REBCB                 , 'USE_REBCB:'           ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, AUTOCONV                  , 'AUTOCONV:'            ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, LAMBAU_DEEP               , 'LAMBAU_DEEP:'         ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,MOIST_TRIGGER             , 'MOIST_TRIGGER:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, SGS_W_TIMESCALE           , 'SGS_W_TIMESCALE:'     ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,FRAC_MODIS                , 'FRAC_MODIS:'          ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,USE_SMOOTH_PROF           , 'USE_SMOOTH_PROF:'     ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,EVAP_FIX                  , 'EVAP_FIX:'            ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, ADV_TRIGGER               , 'ADV_TRIGGER:'         ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(DEEP)       , 'AVE_LAYER_DP:'        ,default= 40.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(SHAL)       , 'AVE_LAYER_SH:'        ,default= 30.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(MID)        , 'AVE_LAYER_MD:'        ,default= 40.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(DEEP)      , 'USE_EXCESS_DP:'       ,default= 2,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(SHAL)      , 'USE_EXCESS_SH:'       ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(MID)       , 'USE_EXCESS_MD:'       ,default= 2,     RC=STATUS );VERIFY_(STATUS)
       if(AUTOCONV == 1) then 
          call MAPL_GetResource(MAPL, C0_DEEP                , 'C0_DEEP:'             ,default= 2.0e-3,RC=STATUS );VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, C0_MID                 , 'C0_MID:'              ,default= 2.0e-3,RC=STATUS );VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, C0_SHAL                , 'C0_SHAL:'             ,default= 0.    ,RC=STATUS );VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, QRC_CRIT               , 'QRC_CRIT:'            ,default= 2.e-4, RC=STATUS );VERIFY_(STATUS)
       else        
          call MAPL_GetResource(MAPL, C0_DEEP                , 'C0_DEEP:'             ,default= 1.5e-3,RC=STATUS );VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, C0_MID                 , 'C0_MID:'              ,default= 1.0e-3,RC=STATUS );VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, C0_SHAL                , 'C0_SHAL:'             ,default= 1.0e-3,RC=STATUS );VERIFY_(STATUS)
          call MAPL_GetResource(MAPL, QRC_CRIT               , 'QRC_CRIT:'            ,default= 1.0e-4,RC=STATUS );VERIFY_(STATUS)
       endif
       call MAPL_GetResource(MAPL, C1                        , 'C1:'                  ,default= 0.,    RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(DEEP)   , 'HEI_DOWN_LAND_DP:'    ,default= 0.65,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(SHAL)   , 'HEI_DOWN_LAND_SH:'    ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(MID)    , 'HEI_DOWN_LAND_MD:'    ,default= 0.65,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(DEEP)  , 'HEI_DOWN_OCEAN_DP:'   ,default= 0.65,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(SHAL)  , 'HEI_DOWN_OCEAN_SH:'   ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(MID)   , 'HEI_DOWN_OCEAN_MD:'   ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(DEEP)   , 'HEI_UPDF_LAND_DP:'    ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(SHAL)   , 'HEI_UPDF_LAND_SH:'    ,default= 0.10,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(MID)    , 'HEI_UPDF_LAND_MD:'    ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(DEEP)  , 'HEI_UPDF_OCEAN_DP:'   ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(SHAL)  , 'HEI_UPDF_OCEAN_SH:'   ,default= 0.10,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(MID)   , 'HEI_UPDF_OCEAN_MD:'   ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(DEEP)    , 'MAX_EDT_LAND_DP:'     ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(SHAL)    , 'MAX_EDT_LAND_SH:'     ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(MID)     , 'MAX_EDT_LAND_MD:'     ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(DEEP)   , 'MAX_EDT_OCEAN_DP:'    ,default= 0.2,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(SHAL)   , 'MAX_EDT_OCEAN_SH:'    ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(MID)    , 'MAX_EDT_OCEAN_MD:'    ,default= 0.2,   RC=STATUS );VERIFY_(STATUS)
   ELSE
       call MAPL_GetResource(MAPL, ZERO_DIFF                 , 'ZERO_DIFF:'           ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, GF_ENV_SETTING            , 'GF_ENV_SETTING:'      ,default= 'DYNAMICS', RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, STOCH_TOP                 , 'STOCH_TOP:'           ,default= 1.25,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, STOCH_BOT                 , 'STOCH_BOT:'           ,default= 0.75,  RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, STOCHASTIC_CNV            , 'STOCHASTIC_CNV:'      ,default= .FALSE., RC=STATUS); VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, GF_MIN_AREA               , 'GF_MIN_AREA:'         ,default= 1.e6,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, TAU_MID                   , 'TAU_MID:'             ,default= 3600., RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, TAU_DEEP                  , 'TAU_DEEP:'            ,default= 5400., RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CLEV_GRID                 , 'CLEV_GRID:'           ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, VERT_DISCR                , 'VERT_DISCR:'          ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, USE_FCT                   , 'USE_FCT:'             ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, SATUR_CALC                , 'SATUR_CALC:'          ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, BC_METH                   , 'BC_METH:'             ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, USE_REBCB                 , 'USE_REBCB:'           ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, AUTOCONV                  , 'AUTOCONV:'            ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, LAMBAU_DEEP               , 'LAMBAU_DEEP:'         ,default= 2.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,MOIST_TRIGGER             , 'MOIST_TRIGGER:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, SGS_W_TIMESCALE           , 'SGS_W_TIMESCALE:'     ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,FRAC_MODIS                , 'FRAC_MODIS:'          ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,USE_SMOOTH_PROF           , 'USE_SMOOTH_PROF:'     ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL ,EVAP_FIX                  , 'EVAP_FIX:'            ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, ADV_TRIGGER               , 'ADV_TRIGGER:'         ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C0_DEEP                   , 'C0_DEEP:'             ,default= 2.e-3, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C0_MID                    , 'C0_MID:'              ,default= 2.e-3, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C0_SHAL                   , 'C0_SHAL:'             ,default= 0.    ,RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, QRC_CRIT                  , 'QRC_CRIT:'            ,default= 2.e-4, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C1                        , 'C1:'                  ,default= 1.e-3, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(DEEP)       , 'AVE_LAYER_DP:'        ,default= 50.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(SHAL)       , 'AVE_LAYER_SH:'        ,default= 30.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(MID)        , 'AVE_LAYER_MD:'        ,default= 50.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(DEEP)      , 'USE_EXCESS_DP:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(SHAL)      , 'USE_EXCESS_SH:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_USE_EXCESS(MID)       , 'USE_EXCESS_MD:'       ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(DEEP)   , 'HEI_DOWN_LAND_DP:'    ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(SHAL)   , 'HEI_DOWN_LAND_SH:'    ,default= 0.2,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_LAND(MID)    , 'HEI_DOWN_LAND_MD:'    ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(DEEP)  , 'HEI_DOWN_OCEAN_DP:'   ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(SHAL)  , 'HEI_DOWN_OCEAN_SH:'   ,default= 0.2,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_DOWN_OCEAN(MID)   , 'HEI_DOWN_OCEAN_MD:'   ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(DEEP)   , 'HEI_UPDF_LAND_DP:'    ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(SHAL)   , 'HEI_UPDF_LAND_SH:'    ,default= 0.10,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_LAND(MID)    , 'HEI_UPDF_LAND_MD:'    ,default= 0.5,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(DEEP)  , 'HEI_UPDF_OCEAN_DP:'   ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(SHAL)  , 'HEI_UPDF_OCEAN_SH:'   ,default= 0.10,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_HEI_UPDF_OCEAN(MID)   , 'HEI_UPDF_OCEAN_MD:'   ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(DEEP)    , 'MAX_EDT_LAND_DP:'     ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(SHAL)    , 'MAX_EDT_LAND_SH:'     ,default= 0.00,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_LAND(MID)     , 'MAX_EDT_LAND_MD:'     ,default= 0.35,  RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(DEEP)   , 'MAX_EDT_OCEAN_DP:'    ,default= 0.9,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(SHAL)   , 'MAX_EDT_OCEAN_SH:'    ,default= 0.0,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_MAX_EDT_OCEAN(MID)    , 'MAX_EDT_OCEAN_MD:'    ,default= 0.9,   RC=STATUS );VERIFY_(STATUS)
    ENDIF
    call MAPL_GetResource(MAPL, SCLM_DEEP       , 'SCLM_DEEP:'       , DEFAULT= 1.0       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CNV_2MOM        , 'CNV_2MOM:'        , DEFAULT= .FALSE.   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, FIX_CNV_CLOUD   , 'FIX_CNV_CLOUD:'   , DEFAULT= .TRUE.    , RC=STATUS); VERIFY_(STATUS)

end subroutine GF_Initialize


subroutine GF_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Local derived type aliases

    type (MAPL_MetaComp), pointer   :: MAPL
    type (ESMF_Config  )            :: CF
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_Alarm   )            :: ALARM
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: DT_MOIST

    ! Local variables
    integer                         :: I, J, L
    integer                         :: IM,JM,LM
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS
    real                            :: minrhx

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
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
    real,    allocatable, dimension(:,:,:) :: TH, KE
    integer, allocatable, dimension(:,:)   :: SEEDINI   
    real,    allocatable, dimension(:,:)   :: SEEDCNV 
    real,    allocatable, dimension(:,:,:) :: TMP3D
    real,    allocatable, dimension(:,:)   :: TMP2D

    ! Required Exports (connectivities to moist siblings)
    real, pointer, dimension(:,:,:) :: CNV_MFD, CNV_MFC, CNV_CVW, CNV_QC, CNV_DQCDT, CNV_PRC3, CNV_UPDF
    real, pointer, dimension(:,:,:) :: DUDT_DC, DVDT_DC, DTDT_DC, DTHDT_DC, DQVDT_DC, DQIDT_DC, DQLDT_DC, DQADT_DC
    real, pointer, dimension(:,:  ) :: CNV_FRC
    ! Exports
    real, pointer, dimension(:,:,:) :: CNV_MF0, ENTLAM
    real, pointer, dimension(:,:,:) :: MUPDP,MDNDP,MUPSH,MUPMD
    real, pointer, dimension(:,:,:) :: VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d
    real, pointer, dimension(:,:  ) :: USTAR,TSTAR,QSTAR,T2M,Q2M,TA,QA,SH,EVAP,PHIS
    real, pointer, dimension(:,:  ) :: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
    real, pointer, dimension(:,:  ) :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC
    real, pointer, dimension(:,:  ) :: TPWI,TPWI_star,LFR_GF,CNPCPRATE
    real, pointer, dimension(:,:,:) :: RSU_CN,REV_CN,PFL_CN,PFI_CN
    real, pointer, dimension(:,:  ) :: SIGMA_DEEP, SIGMA_MID
    real, pointer, dimension(:,:,:) :: CNV_NICE, CNV_NDROP, CNV_FICE
    real, pointer, dimension(:,:,:) :: NCPL, NCPI
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D

    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ) 
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--GF")

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         CF       = CF,                &
         LONS     = LONS,              &
         LATS     = LATS,              &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    ! Internals
    call MAPL_GetPointer(INTERNAL, Q,      'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,   'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,   'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,   'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,   'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,   'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,   'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NACTL,  'NACTL'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NACTI,  'NACTI'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CNV_TR, 'CNV_TR'  , RC=STATUS); VERIFY_(STATUS)

    ! Imports
    call MAPL_GetPointer(IMPORT, FRLAND    ,'FRLAND'    ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, AREA      ,'AREA'      ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ZLE       ,'ZLE'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE       ,'PLE'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T         ,'T'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U         ,'U'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V         ,'V'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, W         ,'W'         ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OMEGA     ,'OMEGA'     ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KH        ,'KH'        ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SH        ,'SH'        ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP      ,'EVAP'      ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, USTAR     ,'USTAR'     ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TSTAR     ,'TSTAR'     ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QSTAR     ,'QSTAR'     ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T2M       ,'T2M'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, Q2M       ,'Q2M'       ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TA        ,'TA'        ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QA        ,'QA'        ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PHIS      ,'PHIS'      ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DQVDTDYN  ,'DQVDTDYN'  ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DTDTDYN   ,'DTDTDYN'   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QV_DYN_IN ,'QV_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U_DYN_IN  ,'U_DYN_IN'  ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V_DYN_IN  ,'V_DYN_IN'  ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T_DYN_IN  ,'T_DYN_IN'  ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE_DYN_IN,'PLE_DYN_IN',RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DTDT_BL   ,'DTDT_BL'   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, DQDT_BL   ,'DQDT_BL'   ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RADSW     ,'RADSW'     ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RADLW     ,'RADLW'     ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBL      ,'KPBL'      ,RC=STATUS); VERIFY_(STATUS)

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PL   (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( TH   (IM,JM,LM  ) )
    ALLOCATE ( KE   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( fQi  (IM,JM,LM  ) )
    ALLOCATE ( QST3 (IM,JM,LM  ) )
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
    KE       = (V**2+U**2)

    ! Required Exports (connectivities to moist siblings)
    call MAPL_GetPointer(EXPORT, CNV_MFD,    'CNV_MFD'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_MFC,    'CNV_MFC'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_CVW,    'CNV_CVW'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_QC,     'CNV_QC'    ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_DQCDT,  'CNV_DQCDT' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_PRC3,   'CNV_PRC3'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_UPDF,   'CNV_UPDF'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_FRC,    'CNV_FRC'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Exports used below
    call MAPL_GetPointer(EXPORT, CNV_MF0,    'CNV_MF0'   ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNPCPRATE,  'CNPCPRATE' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ENTLAM,     'ENTLAM'    ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SIGMA_DEEP, 'SIGMA_DEEP',  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SIGMA_MID,  'SIGMA_MID' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TPWI,       'TPWI'      ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TPWI_star,  'TPWI_star' ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    TPWI       = SUM( Q*MASS, 3 )
    TPWI_star  = SUM( GEOS_QSAT(T,PL,PASCALS=.true.)*MASS, 3 )
    call MAPL_GetPointer(EXPORT, LFR_GF   ,'LFR_GF'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MUPDP    ,'MUPDP'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MDNDP    ,'MDNDP'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MUPSH    ,'MUPSH'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MUPMD    ,'MUPMD'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VAR3d_a  ,'VAR3d_a'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VAR3d_b  ,'VAR3d_b'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VAR3d_c  ,'VAR3d_c'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VAR3d_d  ,'VAR3d_d'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MFDP     ,'MFDP'      ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MFSH     ,'MFSH'      ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, MFMD     ,'MFMD'      ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ERRDP    ,'ERRDP'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ERRSH    ,'ERRSH'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ERRMD    ,'ERRMD'     ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_CN   ,'RSU_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_CN   ,'REV_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_CN   ,'PFI_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_CN   ,'PFL_CN'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AA0      ,'AA0'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AA1      ,'AA1'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AA2      ,'AA2'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AA3      ,'AA3'       ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AA1_BL   ,'AA1_BL'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AA1_CIN  ,'AA1_CIN'   ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAU_BL   ,'TAU_BL'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, TAU_EC   ,'TAU_EC'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    ! 2-moment stuff
    call MAPL_GetPointer(EXPORT, NCPL     ,'NCPL_VOL'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, NCPI     ,'NCPL_VOL'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    TMP3D = PL*R_AIR/T
    NCPL = NACTL/TMP3D ! kg-1
    NCPI = NACTI/TMP3D ! kg-1
    call MAPL_GetPointer(EXPORT, CNV_FICE , 'CNV_FICE' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_NDROP, 'CNV_NDROP',ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_NICE , 'CNV_NICE' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS)

    ! Initialize tendencies
    call MAPL_GetPointer(EXPORT,  DUDT_DC,    'DUDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DVDT_DC,    'DVDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DTDT_DC,    'DTDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTHDT_DC,   'DTHDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQVDT_DC,   'DQVDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDT_DC,   'DQLDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDT_DC,   'DQIDT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQADT_DC,   'DQADT_DC'  ,  ALLOC = .TRUE., RC=STATUS); VERIFY_(STATUS)

    if (STOCHASTIC_CNV) then
       ! Create bit-processor-reproducible random white noise for convection [0:1]
       SEEDINI = 1000000 * ( 100*T(:,:,LM)   - INT( 100*T(:,:,LM) ) )
       SEEDCNV = SQRT(MAX(MIN(SEEDINI/1000000.0,1.0),0.0))
       SEEDCNV = (1.1*CNV_FRC - 0.1) * ((1.0-(1.0-SEEDCNV))*(STOCH_TOP-STOCH_BOT)+STOCH_BOT)
    else
       SEEDCNV = 1.0
    endif
    CALL MAPL_GetPointer(EXPORT, PTR2D,  'STOCH_CNV', RC=STATUS); VERIFY_(STATUS)
    if(associated(PTR2D)) PTR2D = SEEDCNV

! Modify AREA (m^2) here so GF scale dependence has a CNV_FRC dependence
    if (GF_MIN_AREA > 0) then
       where (AREA > GF_MIN_AREA)
          TMP2D = GF_MIN_AREA*CNV_FRC + AREA*(1.0-CNV_FRC)
       elsewhere
          TMP2D = GF_MIN_AREA
       endwhere
    else if (GF_MIN_AREA < 0) then
       TMP2D = ABS(GF_MIN_AREA)
    else
       TMP2D = AREA
    endif

         !- call GF/GEOS5 interface routine
         ! PLE and PL are passed in Pq
         call GF_GEOS5_Interface( IM,JM,LM,LONS,LATS,DT_MOIST                       &
                                 ,PLE, PL, ZLE0, ZL0, PK, MASS, OMEGA, KH           &
                                 ,T, TH, Q, U, V, QLCN, QICN, QLLS, QILS, CNPCPRATE &
                                 ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQCDT, ENTLAM     &
                                 ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC, CLCN, CLLS    &
                                 ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                                 ,RADSW   ,RADLW  ,DQDT_BL  ,DTDT_BL                &
                                 ,FRLAND, TMP2D, USTAR, TSTAR, QSTAR, T2M           &
                                 ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS                       &
                                 ,KPBL ,CNV_FRC                                     &
                                 ,SEEDCNV, SIGMA_DEEP, SIGMA_MID                    &
                                 ,DQVDT_DC,DTDT_DC,DUDT_DC,DVDT_DC                  &
                                 ,MUPDP,MUPSH,MUPMD,MDNDP                           &
                                 ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
                                 ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
                                 ,DTDTDYN,DQVDTDYN                                  &
                                 ,NCPL, NCPI, CNV_NICE, CNV_NDROP, CNV_FICE         &
                                 ,RSU_CN, REV_CN, PFI_CN, PFL_CN                    &
                                 ,TPWI, TPWI_star, LFR_GF                           &
                                 ,VAR3d_a, VAR3d_b, VAR3d_c, VAR3d_d, CNV_TR)

    ! Fill the TH tendency
      DTHDT_DC = DTDT_DC/PK
    ! add tendencies to the moist import state
      U  = U  +  DUDT_DC*DT_MOIST
      V  = V  +  DVDT_DC*DT_MOIST
      Q  = Q  + DQVDT_DC*DT_MOIST
      T  = T  +  DTDT_DC*DT_MOIST
      TH = TH + DTHDT_DC*DT_MOIST
    ! update DeepCu QL/QI/CF tendencies
      TMP3D= CNV_DQCDT/MASS
      fQi  = ice_fraction( T, CNV_FRC )
      DQLDT_DC = (1.0-fQi)*TMP3D
      DQIDT_DC =      fQi *TMP3D
      DQADT_DC = CNV_MFD*SCLM_DEEP/MASS
    ! 2Momoent
     !dNi = make_IceNumber (dQi, TE)
     !dNl = make_DropletNumber (dQl, 0.0, FRLAND)
    ! add QI/QL/CL tendencies
      QLCN =         QLCN + DQLDT_DC*DT_MOIST
      QICN =         QICN + DQIDT_DC*DT_MOIST
      CLCN = MAX(MIN(CLCN + DQADT_DC*DT_MOIST, 1.0), 0.0)

      ! fix 'convective' cloud fraction 
      if (FIX_CNV_CLOUD) then
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
         Q    = Q + QLCN + QICN
         T    = T - (MAPL_ALHL*QLCN + MAPL_ALHS*QICN)/MAPL_CP
         QLCN = 0.
         QICN = 0.
      END WHERE
      endif

      ! Heating from cumulus friction
      call MAPL_GetPointer(EXPORT, PTR3D, 'DTDTFRIC', RC=STATUS); VERIFY_(STATUS)
      if(associated(PTR3D)) then
          KE = (0.5/DT_MOIST)*(KE - (V**2+U**2))*MASS
          TMP3D = 1.e-4 ! KEX
          TMP2D = SUM(KE,3)/MAX(SUM(TMP3D*MASS,3), 1.0e-6) ! IKEX/IKEX2 
          do L=1,LM
             PTR3D(:,:,L) = -(1./MAPL_CP) * TMP2D * TMP3D(:,:,L) * (PLE(:,:,L)-PLE(:,:,L-1))
          end do
      end if

      call MAPL_GetPointer(EXPORT, PTR3D, 'DQRC', RC=STATUS); VERIFY_(STATUS)
      if(associated(PTR3D)) PTR3D = CNV_PRC3 / DT_MOIST

    call MAPL_TimerOff (MAPL,"--GF")

end subroutine GF_Run

end module GEOS_GF_InterfaceMod
