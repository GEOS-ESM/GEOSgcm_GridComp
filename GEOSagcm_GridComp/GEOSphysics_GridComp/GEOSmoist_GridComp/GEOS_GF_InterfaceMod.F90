! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GF_InterfaceMod -- A Module to interface with the
!   GF convection

module GEOS_GF_InterfaceMod

  use ESMF
  use MAPL

  USE ConvPar_GF_GEOS5, only: GF_GEOS5_INTERFACE, MAXIENS, ICUMULUS_GF, CLOSURE_CHOICE, DEEP, SHAL, MID &
                             ,USE_SCALE_DEP,DICYCLE,TAU_DEEP,TAU_MID,HCTS                               &
                             ,USE_TRACER_TRANSP, USE_TRACER_SCAVEN,USE_MEMORY,CONVECTION_TRACER         &
                             ,USE_FLUX_FORM,USE_TRACER_EVAP,DOWNDRAFT,USE_FCT                           &
                             ,USE_REBCB, VERT_DISCR, SATUR_CALC, CLEV_GRID, APPLY_SUB_MP, ALP1          &
                             ,SGS_W_TIMESCALE, LIGHTNING_DIAG,AUTOCONV, BC_METH,OVERSHOOT,USE_WETBULB   &
                             ,C1,C0_DEEP, QRC_CRIT,LAMBAU_DEEP,LAMBAU_SHDN,C0_MID                       &
                             ,CUM_MAX_EDT_LAND,CUM_MAX_EDT_OCEAN, CUM_HEI_DOWN_LAND, CUM_HEI_DOWN_OCEAN &
                             ,CUM_HEI_UPDF_LAND,CUM_HEI_UPDF_OCEAN,USE_MOMENTUM_TRANSP,CUM_ENTR_RATE    &
                             ,ZERO_DIFF,MOIST_TRIGGER,FRAC_MODIS,MAX_TQ_TEND,CUM_FADJ_MASSFLX           &
                             ,CUM_USE_EXCESS,CUM_AVE_LAYER,ADV_TRIGGER,EVAP_FIX,USE_SMOOTH_PROF         &
                             ,OUTPUT_SOUND,TAU_OCEA_CP,TAU_LAND_CP,USE_CLOUD_DISSIPATION,USE_SMOOTH_TEND&
                             ,BETA_SH,C0_SHAL,USE_LINEAR_SUBCL_MF,CAP_MAXS

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
         RESTART = MAPL_RestartSkip,                               &
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
    call MAPL_GetResource(MAPL, CLOSURE_CHOICE(SHAL)        ,'CLOSURE_SHALLOW:'       ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CLOSURE_CHOICE(MID)         ,'CLOSURE_CONGESTUS:'     ,default= 3,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_ENTR_RATE(DEEP)         ,'ENTR_DP:'               ,default= 1.e-4,RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_ENTR_RATE(SHAL)         ,'ENTR_SH:'               ,default= 2.e-3,RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_ENTR_RATE(MID)          ,'ENTR_MD:'               ,default= 5.e-4,RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(DEEP)      ,'FADJ_MASSFLX_DP:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(SHAL)      ,'FADJ_MASSFLX_SH:'       ,default= 1.0,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_FADJ_MASSFLX(MID)       ,'FADJ_MASSFLX_MD:'       ,default= 0.5,  RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_USE_EXCESS(DEEP)        ,'USE_EXCESS_DP:'         ,default= 2,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_USE_EXCESS(SHAL)        ,'USE_EXCESS_SH:'         ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CUM_USE_EXCESS(MID)         ,'USE_EXCESS_MD:'         ,default= 2,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_TRACER_TRANSP           ,'USE_TRACER_TRANSP:'     ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_TRACER_SCAVEN           ,'USE_TRACER_SCAVEN:'     ,default= 2,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_SCALE_DEP               ,'USE_SCALE_DEP:'         ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_MOMENTUM_TRANSP         ,'USE_MOMENTUM_TRANSP:'   ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, DICYCLE                     ,'DICYCLE:'               ,default= 1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, OUTPUT_SOUND                ,'OUTPUT_SOUND:'          ,default= 0,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, USE_MEMORY                  ,'USE_MEMORY:'            ,default=-1,    RC=STATUS );VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, CONVECTION_TRACER           ,'CONVECTION_TRACER:'     ,default= 0,    RC=STATUS );VERIFY_(STATUS)
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
       call MAPL_GetResource(MAPL, TAU_MID                   , 'TAU_MID:'             ,default= 3600., RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, TAU_DEEP                  , 'TAU_DEEP:'            ,default= 5400., RC=STATUS );VERIFY_(STATUS)
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
       call MAPL_GetResource(MAPL, ZERO_DIFF                 , 'ZERO_DIFF:'           ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, ADV_TRIGGER               , 'ADV_TRIGGER:'         ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(DEEP)       , 'AVE_LAYER_DP:'        ,default= 40.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(SHAL)       , 'AVE_LAYER_SH:'        ,default= 30.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(MID)        , 'AVE_LAYER_MD:'        ,default= 40.,   RC=STATUS );VERIFY_(STATUS)
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
       call MAPL_GetResource(MAPL, TAU_MID                   , 'TAU_MID:'             ,default= 3600., RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, TAU_DEEP                  , 'TAU_DEEP:'            ,default= 10800.,RC=STATUS );VERIFY_(STATUS)
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
       call MAPL_GetResource(MAPL, ZERO_DIFF                 , 'ZERO_DIFF:'           ,default= 1,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, ADV_TRIGGER               , 'ADV_TRIGGER:'         ,default= 0,     RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C0_DEEP                   , 'C0_DEEP:'             ,default= 2.e-3, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C0_MID                    , 'C0_MID:'              ,default= 2.e-3, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C0_SHAL                   , 'C0_SHAL:'             ,default= 0.    ,RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, QRC_CRIT                  , 'QRC_CRIT:'            ,default= 2.e-4, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, C1                        , 'C1:'                  ,default= 1.e-3, RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(DEEP)       ,'AVE_LAYER_DP:'         ,default= 50.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(SHAL)       ,'AVE_LAYER_SH:'         ,default= 30.,   RC=STATUS );VERIFY_(STATUS)
       call MAPL_GetResource(MAPL, CUM_AVE_LAYER(MID)        ,'AVE_LAYER_MD:'         ,default= 50.,   RC=STATUS );VERIFY_(STATUS)
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

    integer                         :: IM,JM,LM
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS

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

#ifdef NODISABLE

    call MAPL_GetPointer(INTERNAL, CNV_TR, 'CNV_TR', RC=STATUS); VERIFY_(STATUS)

         if (STOCHASTIC_CNV /= 0) then
           SEEDINI(:,:,1) = 1000000 * ( 100*TEMP(:,:,LM)   - INT( 100*TEMP(:,:,LM) ) )
          ! Create bit-processor-reproducible random white noise for convection [0:1]
           SEEDCNV(:,:)   = SEEDINI(:,:,1)/1000000.0
           where (SEEDCNV > 1.0)
              SEEDCNV = 1.0
           end where
           where (SEEDCNV < 0.0)
              SEEDCNV = 0.0
           end where
          !SEEDCNV = SEEDCNV*(1.875-1.0)+1.0
           SEEDCNV = SEEDCNV*(1.75-0.5)+0.5
         else
           SEEDCNV(:,:) = 1.0
         endif
         CALL MAPL_GetPointer(EXPORT, STOCH_CNV,  'STOCH_CNV', RC=STATUS)
         VERIFY_(STATUS)
         if (associated(STOCH_CNV)) STOCH_CNV = SEEDCNV
         if(associated(DTDT_DC)) DTDT_DC=TH1*PK
         if(associated(DQVDT_DC)) DQVDT_DC=Q1
         if(associated(DQLDT_DC)) DQLDT_DC=QLLS+QLCN
         if(associated(DQIDT_DC)) DQIDT_DC=QILS+QICN
         if(associated(DQADT_DC)) DQADT_DC=CLLS+CLCN
         !-initialize/reset output arrays 
         CNV_MF0  =0.0 ! 'cloud_base_mass_flux'              - 'kg m-2 s-1'
         CNV_MFD  =0.0 ! 'detraining_mass_flux',             - 'kg m-2 s-1'    
         CNV_MFC  =0.0 ! 'cumulative_mass_flux',             - 'kg m-2 s-1'   
         CNV_CVW  =0.0 ! 'updraft_vertical_velocity',        - 'hPa s-1', 
         CNV_QC   =0.0 ! 'grid_mean_convective_condensate',  - 'kg kg-1',   
         CNV_UPDF =0.0 ! 'updraft_areal_fraction',           - '1', 
         ENTLAM   =0.0 ! 'entrainment parameter',            - 'm-1',  
         CNV_PRC3 =0.0 ! 'convective_precipitation           - 'kg m-2 s-1'
         call MAPL_GetPointer(IMPORT, USTAR    ,'USTAR'   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, TSTAR    ,'TSTAR'   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, QSTAR    ,'QSTAR'   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, T2M      ,'T2M  '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, Q2M      ,'Q2M  '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, TA       ,'TA   '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, QA       ,'QA   '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, PHIS     ,'PHIS '   ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DQVDTDYN ,'DQVDTDYN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DTDTDYN  ,'DTDTDYN'  ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, QV_DYN_IN,'QV_DYN_IN',RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, U_DYN_IN ,'U_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, V_DYN_IN ,'V_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, T_DYN_IN ,'T_DYN_IN' ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, PLE_DYN_IN,'PLE_DYN_IN',RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DTDT_BL  ,'DTDT_BL'  ,RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(IMPORT, DQDT_BL  ,'DQDT_BL'  ,RC=STATUS); VERIFY_(STATUS)

         call MAPL_GetPointer(EXPORT, TPWI,     'TPWI'      ,ALLOC = .TRUE. , RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, TPWI_star,'TPWI_star' ,ALLOC = .TRUE. , RC=STATUS); VERIFY_(STATUS)
         TPWI       = SUM( Q1*MASS, 3 )
         TPWI_star  = SUM( GEOS_QSAT(TH1*PK, PLO)*MASS, 3 )
         call MAPL_GetPointer(EXPORT, LFR_GF   ,'LFR_GF'    ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);LFR_GF=0.0
         call MAPL_GetPointer(EXPORT, DTRDT_GF,'DTRDT_GF' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);DTRDT_GF=0.0
         call MAPL_GetPointer(EXPORT, DQDT_GF,'DQDT_GF' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);DQDT_GF=0.0
         call MAPL_GetPointer(EXPORT, DTDT_GF,'DTDT_GF' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);DTDT_GF=0.0
         call MAPL_GetPointer(EXPORT, MUPDP  ,'MUPDP'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MUPDP=0.0
         call MAPL_GetPointer(EXPORT, MDNDP  ,'MDNDP'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MDNDP=0.0
         call MAPL_GetPointer(EXPORT, MUPSH  ,'MUPSH'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MUPSH=0.0
         call MAPL_GetPointer(EXPORT, MUPMD  ,'MUPMD'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MUPMD=0.0
         call MAPL_GetPointer(EXPORT, VAR3d_a,'VAR3d_a'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);VAR3d_a=0.0
         call MAPL_GetPointer(EXPORT, VAR3d_b,'VAR3d_b'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);VAR3d_b=0.0
         call MAPL_GetPointer(EXPORT, VAR3d_c,'VAR3d_c'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);VAR3d_c=0.0
         call MAPL_GetPointer(EXPORT, VAR3d_d,'VAR3d_d'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);VAR3d_d=0.0
         call MAPL_GetPointer(EXPORT, MFDP   ,'MFDP'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MFDP=0.0
         call MAPL_GetPointer(EXPORT, MFSH   ,'MFSH'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MFSH=0.0
         call MAPL_GetPointer(EXPORT, MFMD   ,'MFMD'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);MFMD=0.0
         call MAPL_GetPointer(EXPORT, ERRDP  ,'ERRDP'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);ERRDP=0.0
         call MAPL_GetPointer(EXPORT, ERRSH  ,'ERRSH'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);ERRSH=0.0
         call MAPL_GetPointer(EXPORT, ERRMD  ,'ERRMD'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);ERRMD=0.0
         call MAPL_GetPointer(EXPORT, RSU_CN_GF  ,'RSU_CN_GF'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);RSU_CN_GF=0.0
         call MAPL_GetPointer(EXPORT, REV_CN_GF  ,'REV_CN_GF'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);REV_CN_GF=0.0
         call MAPL_GetPointer(EXPORT, PFI_CN_GF  ,'PFI_CN_GF'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);PFI_CN_GF=0.0
         call MAPL_GetPointer(EXPORT, PFL_CN_GF  ,'PFL_CN_GF'        ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);PFL_CN_GF=0.0
         call MAPL_GetPointer(EXPORT, AA0      ,'AA0'          ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA0=0.0
         call MAPL_GetPointer(EXPORT, AA1      ,'AA1'          ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA1=0.0
         call MAPL_GetPointer(EXPORT, AA2      ,'AA2'          ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA2=0.0
         call MAPL_GetPointer(EXPORT, AA3      ,'AA3'          ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA3=0.0
         call MAPL_GetPointer(EXPORT, AA1_BL   ,'AA1_BL'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA1_BL =0.0
         call MAPL_GetPointer(EXPORT, AA1_CIN  ,'AA1_CIN' ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);AA1_CIN=0.0
         call MAPL_GetPointer(EXPORT, TAU_BL   ,'TAU_BL'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);TAU_BL =0.0
         call MAPL_GetPointer(EXPORT, TAU_EC   ,'TAU_EC'  ,ALLOC = .TRUE. ,RC=STATUS); VERIFY_(STATUS);TAU_EC =0.0

! Modify AREA (m^2) here so GF scale dependence has a CNV_FRACTION dependence
         if (GF_MIN_AREA > 0) then
           GF_AREA = AREA
           WHERE (AREA > GF_MIN_AREA)
             GF_AREA = GF_MIN_AREA*CNV_FRACTION + AREA*(1.0-CNV_FRACTION)
           END WHERE
         else if (GF_MIN_AREA < 0) then
           GF_AREA = ABS(GF_MIN_AREA)
         else
          !GF_AREA = AREA
           GF_AREA = (111000.00*(360.0/imsize))**2 ! KM^2
         endif
         !- call GF/GEOS5 interface routine
         call GF_GEOS5_Interface( IM,JM,LM,NTR,ITRCR,LONS,LATS,DT_MOIST              &
                                 ,T, PLE, PLO, ZLE, ZLO, PK, U, V, OMEGA, KH        &
                                 ,TH1, Q1, U1, V1, QLCN, QICN,QLLS,QILS, CNVPRCP    &
                                 ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQLDT,ENTLAM      &
                                 ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC , CLCN,CLLS    &
                                 ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                                 ,RADSW   ,RADLW  ,DQDT_BL  ,DTDT_BL                &
                                 ,FRLAND, GF_AREA,USTAR,TSTAR,QSTAR,T2M             &
                                 ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS                       &
                                 ,KPBL                                              &
                                 ,SEEDCNV, SIGMA_DEEP, SIGMA_MID                    &
                                 ,DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD,MDNDP           &
                                 ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
                                 ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
                                 ,DTDTDYN,DQVDTDYN                                  &
                                 ,NCPL, NCPI, CNV_NICE, CNV_NDROP, CNV_FICE, CLDMICR_OPTION &
                                 ,XHO,FSCAV,CNAMES,QNAMES,DTRDT_GF                  &
                                 ,RSU_CN_GF,REV_CN_GF, PFI_CN_GF, PFL_CN_GF         &
                                 ,TPWI,TPWI_star,LFR_GF                             &
                                 ,VAR3d_a,VAR3d_b,VAR3d_c,VAR3d_d,CNV_TR)
#endif

    call MAPL_TimerOff (MAPL,"--GF")

end subroutine GF_Run

end module GEOS_GF_InterfaceMod
