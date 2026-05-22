! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GFDL_1M_InterfaceMod -- A Module to interface with the
!   GFDL_1M cloud microphysics

module GEOS_GFDL_1M_InterfaceMod

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use GEOS_RadarMod
  use GEOSmoist_Process_Library
  use Aer_Actv_Single_Moment
  use gfdl2_cloud_microphys_mod, only : gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver, ICE_LSC_VFALL_PARAM, ICE_CNV_VFALL_PARAM
  use gfdl_mp_mod, only : gfdl_mp_init, gfdl_mp_driver, do_ref, do_hail, do_sedi_heat, do_sedi_melt_qi, do_sedi_melt_qs, do_sedi_melt_qg, ifflag

  implicit none

  private

  character(len=ESMF_MAXSTR)              :: IAm
  integer                                 :: STATUS

  ! specify how to handle friendlies with DYN:TRB:CHM:ANA
  type FRIENDLIES_TYPE
         character(len=ESMF_MAXSTR) :: QV
         character(len=ESMF_MAXSTR) :: CLLS
         character(len=ESMF_MAXSTR) :: CLCN
         character(len=ESMF_MAXSTR) :: QLLS
         character(len=ESMF_MAXSTR) :: QLCN
         character(len=ESMF_MAXSTR) :: QILS
         character(len=ESMF_MAXSTR) :: QICN
         character(len=ESMF_MAXSTR) :: QRAIN
         character(len=ESMF_MAXSTR) :: QSNOW
         character(len=ESMF_MAXSTR) :: QGRAUPEL
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  character(len=ESMF_MAXSTR)        :: COMP_NAME

  ! Local resource variables
  real    :: TURNRHCRIT_PARAM
  real    :: MIN_RH_CRIT, MAX_RH_CRIT, MIN_RH_UNSTABLE, MIN_RH_STABLE
  real    :: TAU_EVAP, CCW_EVAP_EFF
  real    :: TAU_SUBL, CCI_EVAP_EFF
  integer :: PDFSHAPE
  real    :: ANV_ICEFALL
  real    :: LS_ICEFALL
  real    :: FAC_RL
  real    :: MIN_RL
  real    :: MAX_RL
  real    :: FAC_RI
  real    :: MIN_RI
  real    :: MAX_RI
  logical :: LHYDROSTATIC
  logical :: LPHYS_HYDROSTATIC
  logical :: LMELTFRZ_CLDMACRO
  logical :: LMELTFRZ_CLDMICRO
  real    :: GFDL_MP_KLID

  logical :: LIQUID_SKIN_SNOW
  logical :: LIQUID_SKIN_GRAUPEL
  logical :: LIQUID_SKIN_HAIL

  real, PARAMETER :: W_START = 6.0
  real, PARAMETER :: W_FULL  = 12.0
  real :: fraction_hail

  logical :: REPORT_GFDL_1M_NEGATIVES

  logical :: GFDL_MP3

  public :: GFDL_1M_Setup, GFDL_1M_Initialize, GFDL_1M_Run

contains

subroutine GFDL_1M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_GFDL_1M_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    ! !INTERNAL STATE:

      FRIENDLIES%QV       = "DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS"
      FRIENDLIES%CLLS     = "DYNAMICS"
      FRIENDLIES%CLCN     = "DYNAMICS"
      FRIENDLIES%QLLS     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QLCN     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QILS     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QICN     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QRAIN    = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QSNOW    = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QGRAUPEL = "DYNAMICS:TURBULENCE"

    !BOS

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'Q',                                         &
         LONG_NAME  = 'specific_humidity',                         &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QV),                         &
         default    = 1.0e-6,                                      &
         RESTART    = MAPL_RestartRequired,                        &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                        &
         SHORT_NAME = 'QLLS',                                            &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%QLLS),                             &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = trim(FRIENDLIES%QLCN),                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLLS',                                      &
         LONG_NAME  = 'large_scale_cloud_area_fraction',           &
         UNITS      = '1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%CLLS),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLCN',                                      &
         LONG_NAME  = 'convective_cloud_area_fraction',            &
         UNITS      = '1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%CLCN),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                     &
         SHORT_NAME = 'QILS',                                         &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                      &
         FRIENDLYTO = trim(FRIENDLIES%QILS),                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME = 'QICN',                                        &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                     &
         FRIENDLYTO = trim(FRIENDLIES%QICN),                         &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,               RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QRAIN',                                     &
         LONG_NAME  = 'mass_fraction_of_rain',                     &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QRAIN),                      &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QSNOW',                                     &
         LONG_NAME  = 'mass_fraction_of_snow',                     &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QSNOW),                      &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QGRAUPEL',                                  &
         LONG_NAME  = 'mass_fraction_of_graupel',                  &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QGRAUPEL),                   &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                               &
         SHORT_NAME = 'NACTL',                                  &
         LONG_NAME  = 'activ aero # conc liq phase for 1-mom',  &
         UNITS      = 'm-3',                                    &
         RESTART    = MAPL_RestartSkip,                         &
         DIMS       = MAPL_DimsHorzVert,                        &
         VLOCATION  = MAPL_VLocationCenter,     RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddInternalSpec(GC,                               &
         SHORT_NAME = 'NACTI',                                  &
         LONG_NAME  = 'activ aero # conc ice phase for 1-mom',  &
         UNITS      = 'm-3',                                    &
         RESTART    = MAPL_RestartSkip,                         &
         DIMS       = MAPL_DimsHorzVert,                        &
         VLOCATION  = MAPL_VLocationCenter,     RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'REF_DBZ',                                          &
         LONG_NAME = 'Simulated_gfdl_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &  
         DIMS      = MAPL_DimsHorzVert,                            &
         VLOCATION = MAPL_VLocationCenter,              RC=STATUS  )
    VERIFY_(STATUS) 
    
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'REF_DBZ_MAX',                                          &
         LONG_NAME = 'Maximum_composite_gfdl_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &    
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )   
    VERIFY_(STATUS)
         
    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'REF_DBZ_1KM',                                          &
         LONG_NAME = 'Base_1KM_AGL_gfdl_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'REF_DBZ_TOP',                                          &
         LONG_NAME = 'Echo_top_gfdl_radar_reflectivity',                  &
         UNITS     = 'm',                                     &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_AddExportSpec(GC,                               &
         SHORT_NAME = 'REF_DBZ_M10C',                                          &
         LONG_NAME = 'Minus_10C_gfdl_radar_reflectivity',                  &
         UNITS     = 'dBZ',                                     &
         DIMS      = MAPL_DimsHorzOnly,                            &
         VLOCATION = MAPL_VLocationNone,              RC=STATUS  )
    VERIFY_(STATUS)

    call MAPL_TimerAdd(GC, name="--GFDL_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine GFDL_1M_Setup

subroutine GFDL_1M_Initialize (MAPL, CLOCK, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    type (ESMF_Alarm   )                :: ALARM
    type (ESMF_TimeInterval)            :: TINT
    real(ESMF_KIND_R8)                  :: DT_R8
    real                                :: DT_MOIST

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL

    CHARACTER(len=ESMF_MAXSTR) :: errmsg

    real                     :: DBZ_DT
    type(ESMF_Calendar)      :: calendar
    type(ESMF_Time)          :: currTime
    type(ESMF_Alarm)         :: DBZ_RunAlarm
    type(ESMF_TimeInterval)  :: ringInterval
    integer                  :: LM, year, month, day, hh, mm, ss

    call MAPL_Get(MAPL, LM=LM, RUNALARM=ALARM, RC=STATUS );VERIFY_(STATUS)
    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    DBZ_DT = max(DT_MOIST,900.0)
    call MAPL_GetResource(MAPL, DBZ_DT, 'DBZ_DT:', default=DBZ_DT, RC=STATUS); VERIFY_(STATUS)
    ! Get the current time in addition to the calendar
    call ESMF_ClockGet(CLOCK, currTime=currTime, calendar=calendar, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalSet(ringInterval, S=nint(DBZ_DT), calendar=calendar, RC=STATUS); VERIFY_(STATUS)
    ! Add RingTime = currTime to anchor the alarm
    DBZ_RunAlarm = ESMF_AlarmCreate(Clock       = CLOCK,          &
                                   Name         = 'DBZ_RunAlarm', &
                                   RingTime     = currTime-TINT,  &
                                   RingInterval = ringInterval,   &
                                   Sticky       = .false.  , RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, LPHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    LHYDROSTATIC = LPHYS_HYDROSTATIC
    call MAPL_GetResource( MAPL, LMELTFRZ_CLDMACRO, Label="MELTFRZ_CLDMACRO:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LMELTFRZ_CLDMICRO, Label="MELTFRZ_CLDMICRO:",  default=.FALSE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, REPORT_GFDL_1M_NEGATIVES, 'REPORT_GFDL_1M_NEGATIVES:', default=.FALSE., RC=STATUS) ; VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, GFDL_MP3, Label="GFDL_MP3:",  default=.TRUE., RC=STATUS); VERIFY_(STATUS)
    if (DT_R8 <= 150.0) do_hail = .true. 
    if (DT_R8 <= 150.0) do_sedi_heat = .true.
    if (DT_R8 <= 150.0) do_sedi_melt_qi = .true.
    if (DT_R8 <= 150.0) do_sedi_melt_qs = .true.
    if (DT_R8 <= 150.0) do_sedi_melt_qg = .true.
    if (DT_R8 <= 150.0) ifflag = 1

    if (GFDL_MP3) then
      call gfdl_mp_init(LHYDROSTATIC,DT_MOIST)
      call WRITE_PARALLEL ("INITIALIZED GFDL_1M gfdl_mp v3 in non-generic GC INIT")
      call MAPL_GetResource( MAPL, do_ref, Label="DO_GFDL_REFLECTIVITY:",  default=.TRUE., RC=STATUS); VERIFY_(STATUS)
    else  
      call gfdl_cloud_microphys_init()
      call WRITE_PARALLEL ("INITIALIZED GFDL_1M gfdl_cloud_microphys in non-generic GC INIT")
      do_ref = .false.  ! Force to false so MAPL DBZ Calc triggers, as older driver has no DBZ3D
    endif 

    call MAPL_GetResource( MAPL, SH_MD_DP        , 'SH_MD_DP:'        , DEFAULT= .TRUE., RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, DBZ_VAR_INTERCP , 'DBZ_VAR_INTERCP:' , DEFAULT= DBZ_VAR_INTERCP, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, LIQUID_SKIN_SNOW    , 'LIQUID_SKIN_SNOW:'    , DEFAULT= .FALSE. , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LIQUID_SKIN_GRAUPEL , 'LIQUID_SKIN_GRAUPEL:' , DEFAULT= .FALSE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LIQUID_SKIN_HAIL    , 'LIQUID_SKIN_HAIL:'    , DEFAULT= .TRUE. , RC=STATUS); VERIFY_(STATUS)

                                 refl10cm_allow_wet_graupel = .false.
    call MAPL_GetResource( MAPL, refl10cm_allow_wet_graupel , 'refl10cm_allow_wet_graupel:' , &
                        DEFAULT= refl10cm_allow_wet_graupel, RC=STATUS); VERIFY_(STATUS)
                                 refl10cm_allow_wet_snow    = .false.
    call MAPL_GetResource( MAPL, refl10cm_allow_wet_snow    , 'refl10cm_allow_wet_snow:'    , &
                        DEFAULT= refl10cm_allow_wet_snow, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, constrain_modis_ice, 'constrain_modis_ice:', DEFAULT= .FALSE., RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, TURNRHCRIT_PARAM, 'TURNRHCRIT:'      , DEFAULT= -9999., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RH_CRIT     , 'MAX_RH_CRIT:'     , DEFAULT= 1.0000, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RH_UNSTABLE , 'MIN_RH_UNSTABLE:' , DEFAULT= 0.9750, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RH_STABLE   , 'MIN_RH_STABLE:'   , DEFAULT= 0.8750, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, PDFSHAPE        , 'PDFSHAPE:'        , DEFAULT= 1     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ICE_LSC_VFALL_PARAM, 'ICE_LSC_VFALL_PARAM:',DEFAULT= 1, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ICE_CNV_VFALL_PARAM, 'ICE_CNV_VFALL_PARAM:',DEFAULT= 1, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ANV_ICEFALL     , 'ANV_ICEFALL:'     , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LS_ICEFALL      , 'LS_ICEFALL:'      , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LIQ_RADII_PARAM , 'LIQ_RADII_PARAM:' , DEFAULT= 1     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ICE_RADII_PARAM , 'ICE_RADII_PARAM:' , DEFAULT= 2     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RI          , 'FAC_RI:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RI          , 'MIN_RI:'          , DEFAULT=  5.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RI          , 'MAX_RI:'          , DEFAULT=100.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RL          , 'FAC_RL:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RL          , 'MIN_RL:'          , DEFAULT= 2.5e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RL          , 'MAX_RL:'          , DEFAULT=60.0e-6, RC=STATUS); VERIFY_(STATUS)

    ! USE_BERGERON should be .TRUE. only when USE_AEROSOL_NN is also .TRUE.
    call MAPL_GetResource( MAPL, USE_BERGERON    , 'USE_BERGERON:'    , DEFAULT=USE_AEROSOL_NN, RC=STATUS); VERIFY_(STATUS)

                                 CCW_EVAP_EFF = 4.e-3
    call MAPL_GetResource( MAPL, CCW_EVAP_EFF, 'CCW_EVAP_EFF:', DEFAULT= CCW_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

                                 CCI_EVAP_EFF = 4.e-3
    call MAPL_GetResource( MAPL, CCI_EVAP_EFF, 'CCI_EVAP_EFF:', DEFAULT= CCI_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=  500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, GFDL_MP_KLID    , 'GFDL_MP_KLID:'    , DEFAULT= -999.0, RC=STATUS); VERIFY_(STATUS)

    call init_refl10cm()

end subroutine GFDL_1M_Initialize

subroutine GFDL_1M_Run (GC, IMPORT, EXPORT, CLOCK, RC)
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
    logical                         :: alarm_is_ringing

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W, KH
    real, pointer, dimension(:,:)   :: AREA, PHIS, FRLAND, TS, DTSX, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:) :: SL2, SL3, QT2, QT3, W2, W3, SLQT, WQT, WQL, WSL, PDF_A
    real, pointer, dimension(:,:,:) :: WTHV2
    real, pointer, dimension(:,:,:) :: OMEGA
    ! Local
    real, allocatable, dimension(:,:,:) :: U0, V0
    real, allocatable, dimension(:,:,:) :: PLEmb, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  ZL0
    real, allocatable, dimension(:,:,:) :: DZ, DZET, DP, MASS, iMASS
    real, allocatable, dimension(:,:,:) :: DQST3, QST3
    real, allocatable, dimension(:,:,:) :: DBZ3D, TMP_NACTR
    real, allocatable, dimension(:,:,:) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                                           DQSDTmic, DQGDTmic, DQADTmic, &
                                            DUDTmic,  DVDTmic,  DTDTmic, DWDTmic
    real, allocatable, dimension(:,:,:) :: QHAIL_, QGRAUPEL_
    integer, allocatable, dimension(:,:):: KLCL
    integer                             :: KLID
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D
    real, allocatable, dimension(:)     :: TMP1D
    ! Exports
    real, pointer, dimension(:,:  ) :: LONS, LATS
    real, pointer, dimension(:,:,:) :: NACTR
    real, pointer, dimension(:,:  ) :: PRCP_WATER, PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ) :: LS_PRCP, LS_SNR, ICE, FRZR, CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:) ::  DUDT_macro,  DVDT_macro,  DTDT_macro
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro
    real, pointer, dimension(:,:,:) ::  DUDT_micro,  DVDT_micro,  DTDT_micro
    real, pointer, dimension(:,:,:) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, pointer, dimension(:,:,:) :: CLDREFFL, CLDREFFI
    real, pointer, dimension(:,:,:) :: EVAPC, SUBLC
    real, pointer, dimension(:,:,:) :: RHX, REV_LS, RSU_LS
    real, pointer, dimension(:,:,:) :: VFALL_ICE, VFALL_SNOW, VFALL_GRAUPEL, VFALL_RAIN
    real, pointer, dimension(:,:,:) :: PFL_LS, PFL_AN
    real, pointer, dimension(:,:,:) :: PFI_LS, PFI_AN
    real, pointer, dimension(:,:,:) :: PFR_LS, PFS_LS, PFG_LS
    real, pointer, dimension(:,:,:) :: PDFITERS
    real, pointer, dimension(:,:,:) :: RHCRIT3D
    real, pointer, dimension(:,:,:) :: CNV_PRC3 
    real, pointer, dimension(:,:)   :: EIS, LTS
    real, pointer, dimension(:,:)   :: DBZ_MAX, DBZ_1KM, DBZ_TOP, DBZ_M10C
    real, pointer, dimension(:,:,:) :: DBZ
    real, pointer, dimension(:,:,:) ::   DQVDT_FILL
    real, pointer, dimension(:,:,:) :: DQLLSDT_FILL
    real, pointer, dimension(:,:,:) :: DQLCNDT_FILL
    real, pointer, dimension(:,:,:) :: DQILSDT_FILL
    real, pointer, dimension(:,:,:) :: DQICNDT_FILL
    real, pointer, dimension(:,:,:) ::   DQRDT_FILL
    real, pointer, dimension(:,:,:) ::   DQSDT_FILL
    real, pointer, dimension(:,:,:) ::   DQGDT_FILL
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D

    ! Local variables
    real    :: tmp_val, rand1
    real    :: x_norm, safe_max_rh_crit
    real    :: ALPHA, RHCRIT
    real    :: one_minus_sigma
    real    :: fraction_hail

    real, allocatable :: facEIS_2d(:,:), minrhcrit_2d(:,:), turnrhcrit_2d(:,:)
    real, allocatable :: qg_col(:), qh_col(:), prs_col(:), dbz_col(:)

    integer :: IM,JM,LM
    integer :: I, J, L
    type( ESMF_VM ) :: VMG

    call ESMF_GridCompGet( GC, VM=VMG, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--GFDL_1M",RC=STATUS)

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         LATS     = LATS,              & ! These are in radians
         LONS     = LONS,              & ! These are in radians
         RUNALARM = ALARM,             &
         CF       = CF,                &
         INTERNAL_ESMF_STATE=INTERNAL, &
         RC=STATUS )
    VERIFY_(STATUS)

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NACTL,    'NACTL'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NACTI,    'NACTI'   , RC=STATUS); VERIFY_(STATUS)

    ! Import State
    call MAPL_GetPointer(IMPORT, AREA,    'AREA'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, W,       'W'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, FRLAND,  'FRLAND'  , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KH,      'KH'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PDF_A,   'PDF_A'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, W2,      'W2'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, W3,      'W3'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WQT,     'WQT'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WSL,     'WSL'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SL2,     'SL2'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SL3,     'SL3'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QT2,     'QT2'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QT3,     'QT3'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SLQT,    'SLQT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TS,      'TS'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBLSC,  'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SH,      'SH'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP,    'EVAP'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PHIS,    'PHIS'    , RC=STATUS); VERIFY_(STATUS)

    ! Allocatables
     ! Edge variables
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PLEmb(IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( U0   (IM,JM,LM  ) )
    ALLOCATE ( V0   (IM,JM,LM  ) )
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PLmb (IM,JM,LM  ) )
    ALLOCATE ( DZET (IM,JM,LM  ) )
    ALLOCATE ( DZ   (IM,JM,LM  ) )
    ALLOCATE ( DP   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( iMASS(IM,JM,LM  ) )
    ALLOCATE ( DQST3(IM,JM,LM  ) )
    ALLOCATE (  QST3(IM,JM,LM  ) )
    ALLOCATE ( DBZ3D(IM,JM,LM  ) )
    ALLOCATE ( TMP3D(IM,JM,LM  ) )
     ! Local tendencies
    ALLOCATE ( DQVDTmic(IM,JM,LM  ) )
    ALLOCATE ( DQLDTmic(IM,JM,LM  ) )
    ALLOCATE ( DQIDTmic(IM,JM,LM  ) )
    ALLOCATE ( DQRDTmic(IM,JM,LM  ) )
    ALLOCATE ( DQSDTmic(IM,JM,LM  ) )
    ALLOCATE ( DQGDTmic(IM,JM,LM  ) )
    ALLOCATE ( DQADTmic(IM,JM,LM  ) )
    ALLOCATE (  DUDTmic(IM,JM,LM  ) )
    ALLOCATE (  DVDTmic(IM,JM,LM  ) )
    ALLOCATE (  DTDTmic(IM,JM,LM  ) )
    ALLOCATE (  DWDTmic(IM,JM,LM  ) )      
     ! 2D Variables
    ALLOCATE ( TMP2D        (IM,JM) )
     ! 1D Variables
    ALLOCATE ( TMP1D   (      LM  ) )   

    ! Initialize to clear DBZ
    DBZ3D = -30.0

    ! Derived States
    PLEmb    =  PLE*.01
    PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
    DO L=0,LM
       ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
    END DO
    ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
    DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
    DQST3    = GEOS_DQSAT(T, PLmb, QSAT=QST3)
    DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
    MASS     = DP/MAPL_GRAV
    iMASS    = 1.0/MASS
    U0       = U
    V0       = V
    KLCL     = FIND_KLCL( T, Q, PLmb, IM, JM, LM )
    if (GFDL_MP_KLID > 0.0) then
        KLID = GFDL_MP_KLID
    else
        KLID = 1
    endif

    ! Export and/or scratch Variable
    call MAPL_GetPointer(EXPORT, RAD_CF,   'FCLD', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QV,   'QV'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QL,   'QL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QI,   'QI'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QR,   'QR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QS,   'QS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QG,   'QG'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CLDREFFL, 'RL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CLDREFFI, 'RI'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! This export MUST have been filled in the GridComp
    call MAPL_GetPointer(EXPORT, CNV_FRC,      'CNV_FRC'      , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SRF_TYPE,     'SRF_TYPE'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Exports  required below
    call MAPL_GetPointer(EXPORT, EVAPC,        'EVAPC'        , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SUBLC,        'SUBLC'        , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PRCP_WATER,   'PRCP_WATER'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PRCP_RAIN,    'PRCP_RAIN'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PRCP_SNOW,    'PRCP_SNOW'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PRCP_ICE,     'PRCP_ICE'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PRCP_GRAUPEL, 'PRCP_GRAUPEL' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Exports to be filled
    call MAPL_GetPointer(EXPORT, LS_PRCP,  'LS_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, LS_SNR,   'LS_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ICE,      'ICE'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZR,     'FRZR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RHX   ,   'RHX'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_LS,   'REV_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_LS,   'RSU_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_AN,   'PFL_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_LS,   'PFL_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_AN,   'PFI_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_LS,   'PFI_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFR_LS,   'PFR_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFS_LS,   'PFS_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFG_LS,   'PFG_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WTHV2,    'WTHV2'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WQL,      'WQL'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PDFITERS, 'PDFITERS', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALL_ICE,     'VFALL_ICE'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALL_SNOW,    'VFALL_SNOW'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALL_GRAUPEL, 'VFALL_GRAUPEL', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALL_RAIN,    'VFALL_RAIN'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Unused Exports (forced to 0.0)
    call MAPL_GetPointer(EXPORT, PTR2D,  'CN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'CN_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    ! Lower tropospheric stability and estimated inversion strength from MoistGC
    call MAPL_GetPointer(EXPORT, LTS,   'LTS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EIS,   'EIS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)


    if (DEBUG_TQ_ERRORS) then
         do L = 1, LM
      do J=1,JM
         do I=1,IM
             if (  (       T(I,J,L) > 333.0) .OR. (       T(I,J,L) /=        T(I,J,L)) .OR. &
                   (       Q(I,J,L) < 0.0  ) .OR. (       Q(I,J,L) /=        Q(I,J,L)) .OR. &
                   (    QLLS(I,J,L) < 0.0  ) .OR. (    QLLS(I,J,L) /=     QLLS(I,J,L)) .OR. &
                   (    QLCN(I,J,L) < 0.0  ) .OR. (    QLCN(I,J,L) /=     QLCN(I,J,L)) .OR. &
                   (    QILS(I,J,L) < 0.0  ) .OR. (    QILS(I,J,L) /=     QILS(I,J,L)) .OR. &
                   (    QICN(I,J,L) < 0.0  ) .OR. (    QICN(I,J,L) /=     QICN(I,J,L)) .OR. &
                   (   QRAIN(I,J,L) < 0.0  ) .OR. (   QRAIN(I,J,L) /=    QRAIN(I,J,L)) .OR. &
                   (   QSNOW(I,J,L) < 0.0  ) .OR. (   QSNOW(I,J,L) /=    QSNOW(I,J,L)) .OR. &
                   (QGRAUPEL(I,J,L) < 0.0  ) .OR. (QGRAUPEL(I,J,L) /= QGRAUPEL(I,J,L)) ) then 
                 print *, "T or Q  spike detected : ", T(I,J,L)
                 print *, "    On Entry to GFDL   : "
                 print *, "  Latitude       =", LATS(I,J)*180.0/MAPL_PI
                 print *, "  Longitude      =", LONS(I,J)*180.0/MAPL_PI
                 print *, "  Pressure (mb)  =", PLmb(I,J,L)
                 print *, "                        CLLS=",   CLLS(I,J,L), "   CLCN=",    CLCN(I,J,L)
                 print *, "    QV=",   Q(I,J,L), " QLLS=",   QLLS(I,J,L), "   QLCN=",    QLCN(I,J,L)
                 print *, "                        QILS=",   QILS(I,J,L), "   QICN=",    QICN(I,J,L)
                 print *, "    QR=", QRAIN(I,J,L), "   QS=",  QSNOW(I,J,L), "     QG=", QGRAUPEL(I,J,L)
               endif
             enddo
         end do
      end do
    endif

    call MAPL_TimerOn(MAPL,"---CLDMACRO")
      call MAPL_GetPointer(EXPORT, DQVDT_macro, 'DQVDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQIDT_macro, 'DQIDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQLDT_macro, 'DQLDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQADT_macro, 'DQADT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQRDT_macro, 'DQRDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQSDT_macro, 'DQSDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, DQGDT_macro, 'DQGDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  DUDT_macro,  'DUDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  DVDT_macro,  'DVDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT,  DTDT_macro,  'DTDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

      !$OMP parallel do default(none) &
      !$OMP shared(LM, JM, IM, DUDT_macro, DVDT_macro, DTDT_macro, DQVDT_macro, &
      !$OMP        DQLDT_macro, DQIDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, &
      !$OMP        DQGDT_macro, U, V, T, Q, QLCN, QLLS, QICN, QILS, CLCN, CLLS, &
      !$OMP        QRAIN, QSNOW, QGRAUPEL, EVAPC, SUBLC, PDFITERS, RHX) &
      !$OMP private(I, J, L)
      do L = 1, LM
         do J = 1, JM
            do I = 1, IM
               DUDT_macro(I,J,L) = U(I,J,L)
               DVDT_macro(I,J,L) = V(I,J,L)
               DTDT_macro(I,J,L) = T(I,J,L)
               DQVDT_macro(I,J,L) = Q(I,J,L)
               DQLDT_macro(I,J,L) = QLCN(I,J,L) + QLLS(I,J,L)
               DQIDT_macro(I,J,L) = QICN(I,J,L) + QILS(I,J,L)
               DQADT_macro(I,J,L) = CLCN(I,J,L) + CLLS(I,J,L)
               DQRDT_macro(I,J,L) = QRAIN(I,J,L)
               DQSDT_macro(I,J,L) = QSNOW(I,J,L)
               DQGDT_macro(I,J,L) = QGRAUPEL(I,J,L)
               ! Clear exports
               EVAPC(I,J,L) = 0.0
               SUBLC(I,J,L) = 0.0
               PDFITERS(I,J,L) = 0.0
               RHX(I,J,L) = 0.0
            enddo
         enddo
      enddo

      ! Include shallow precip condensates if present
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          !$OMP parallel do default(none) &
          !$OMP shared(LM, JM, IM, QRAIN, PTR3D, DT_MOIST) &
          !$OMP private(I, J, L)
          do L = 1, LM
             do J = 1, JM
                do I = 1, IM
                   QRAIN(I,J,L) = QRAIN(I,J,L) + PTR3D(I,J,L) * DT_MOIST
                enddo
             enddo
          enddo
        endif
        
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then 
          !$OMP parallel do default(none) &
          !$OMP shared(LM, JM, IM, QSNOW, PTR3D, DT_MOIST) &
          !$OMP private(I, J, L)
          do L = 1, LM
             do J = 1, JM
                do I = 1, IM
                   QSNOW(I,J,L) = QSNOW(I,J,L) + PTR3D(I,J,L) * DT_MOIST
                enddo
             enddo
          enddo
        endif

      ! -----------------------------------------------------------------
      ! Precalculate 2D arrays (Hoisted out of the L loop)
      ! -----------------------------------------------------------------
        allocate(facEIS_2d(IM,JM), minrhcrit_2d(IM,JM), turnrhcrit_2d(IM,JM))

        !$OMP parallel do default(none) &
        !$OMP shared(IM, JM, EIS, SRF_TYPE, MIN_RH_UNSTABLE, MIN_RH_STABLE, &
        !$OMP        TURNRHCRIT_PARAM, PLmb, KPBLSC, facEIS_2d, minrhcrit_2d, &
        !$OMP        turnrhcrit_2d) &
        !$OMP private(I, J)
        do J=1,JM
          do I=1,IM
             facEIS_2d(I,J) = get_fac_eis(EIS(I,J),SRF_TYPE(I,J))
             minrhcrit_2d(I,J) = MIN_RH_UNSTABLE*(1.0-facEIS_2d(I,J)) + MIN_RH_STABLE*facEIS_2d(I,J)
             minrhcrit_2d(I,J) = max(0.7, minrhcrit_2d(I,J))

             if (TURNRHCRIT_PARAM <= 0.0) then
                turnrhcrit_2d(I,J) = PLmb(I, J, NINT(KPBLSC(I,J))) - 50.
             else
                turnrhcrit_2d(I,J) = TURNRHCRIT_PARAM
             endif
          enddo
        enddo

       ! evap/subl/pdf
        call MAPL_GetPointer(EXPORT, RHCRIT3D,  'RHCRIT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        
        !$OMP parallel do default(none) &
        !$OMP shared(LM, JM, IM, Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, KLID, &
        !$OMP        facEIS_2d, minrhcrit_2d, turnrhcrit_2d, MAX_RH_CRIT, PLmb, PLEmb, &
        !$OMP        AREA, RHCRIT3D, DT_MOIST, PDFSHAPE, CNV_FRC, SRF_TYPE, ZL0, NACTL, &
        !$OMP        NACTI, WSL, WQT, SL2, QT2, SLQT, W3, W2, QT3, SL3, PDF_A, PDFITERS, &
        !$OMP        WTHV2, WQL, USE_BERGERON, RHX, LMELTFRZ_CLDMACRO, CCW_EVAP_EFF, &
        !$OMP        EVAPC, CCI_EVAP_EFF, SUBLC, QST3) &
        !$OMP private(L, J, I, safe_max_rh_crit, x_norm, MIN_RH_CRIT, RHCRIT, ALPHA)
        do L=1,LM
          do J=1,JM
           do I=1,IM
           ! cleanup clouds before cldmacro
             call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), &
                                                     QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L), &
                                                     REMOVE_CLOUDS=(L < KLID) )
           
           ! Use Slingo-Ritter (1985) formulation for critical relative humidity
             ! Ensure the max is never lower than the min
             safe_max_rh_crit = MAX(MAX_RH_CRIT, minrhcrit_2d(I,J))
             if (PLmb(i,j,l) .le. turnrhcrit_2d(I,J)) then 
                MIN_RH_CRIT = minrhcrit_2d(I,J)
             else if (L .eq. LM) then
                MIN_RH_CRIT = safe_max_rh_crit
             else             
                x_norm = (PLmb(i,j,l) - turnrhcrit_2d(I,J)) / (PLEmb(i,j,LM) - turnrhcrit_2d(I,J))
                ! Cubic smoothstep S-curve: x^2 * (3 - 2x)
                MIN_RH_CRIT = minrhcrit_2d(I,J) + (safe_max_rh_crit - minrhcrit_2d(I,J)) * &
                         (x_norm * x_norm * (3.0 - 2.0 * x_norm))
             endif
           ! -----------------------------------------------------------------
           ! Scale-Aware Blending for RHCRIT
           ! -----------------------------------------------------------------
             RHCRIT = MAX_RH_CRIT + (MIN_RH_CRIT-MAX_RH_CRIT)*SQRT(SQRT(AREA(I,J)/1.e10)) 
           ! limit ALPHA to < 30%
             ALPHA = max(0.0,min(0.30, (1.0-RHCRIT)))
           ! fill RHCRIT export
             if (associated(RHCRIT3D)) RHCRIT3D(I,J,L) = 1.0-ALPHA
           ! Do CLOUD MACRO below the pressure lid
             if (L >= KLID) then
           ! Put condensates in touch with the PDF
             call hystpdf( &
                      DT_MOIST       , &
                      ALPHA          , &
                      PDFSHAPE       , &
                      CNV_FRC(I,J)   , &
                      SRF_TYPE(I,J)  , &
                      PLmb(I,J,L)    , &
                      ZL0(I,J,L)     , &
                      Q(I,J,L)       , &
                      QLLS(I,J,L)    , &
                      QLCN(I,J,L)    , &
                      QILS(I,J,L)    , &
                      QICN(I,J,L)    , &
                      T(I,J,L)       , &
                      CLLS(I,J,L)    , &
                      CLCN(I,J,L)    , &
                      NACTL(I,J,L)   , &
                      NACTI(I,J,L)   , &
                      WSL(I,J,L)     , &
                      WQT(I,J,L)     , &
                      SL2(I,J,L)     , &
                      QT2(I,J,L)     , &
                      SLQT(I,J,L)    , &
                      W3(I,J,L)      , &
                      W2(I,J,L)      , &
                      QT3(I,J,L)     , &
                      SL3(I,J,L)     , &
                      PDF_A(I,J,L)   , &
                      PDFITERS(I,J,L), &
                      WTHV2(I,J,L)   , &
                      WQL(I,J,L)     , &
                      .false.        , &
                      USE_BERGERON)
             RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
             if (LMELTFRZ_CLDMACRO) then
           ! meltfrz new condensates
             call MELTFRZ ( DT_MOIST     , &
                            1.0          , & ! since we are explicitly operating on CN types pass this always as 1.0
                            SRF_TYPE(I,J), &
                            T(I,J,L)     , &
                            QLCN(I,J,L)  , &
                            QICN(I,J,L) )
             endif
           ! evaporation for CN
             if (CCW_EVAP_EFF > 0.0) then ! else evap done inside GFDL
             EVAPC(I,J,L) = Q(I,J,L)
             call EVAP3 (         &
                  DT_MOIST      , &
                  CCW_EVAP_EFF  , &
                  RHCRIT        , &
                   PLmb(I,J,L)  , &
                      T(I,J,L)  , &
                      Q(I,J,L)  , &
                   QLCN(I,J,L)  , &
                   QICN(I,J,L)  , &
                   CLCN(I,J,L)  , &
                  NACTL(I,J,L)  , &
                  NACTI(I,J,L)  , &
                   QST3(I,J,L)  )
             EVAPC(I,J,L) = ( Q(I,J,L) - EVAPC(I,J,L) ) / DT_MOIST
             endif
           ! sublimation for CN
             if (CCI_EVAP_EFF > 0.0) then ! else subl done inside GFDL
             SUBLC(I,J,L) = Q(I,J,L)
             call SUBL3 (        &
                  DT_MOIST      , &
                  CCI_EVAP_EFF  , &
                  RHCRIT        , &
                   PLmb(I,J,L)  , &
                      T(I,J,L)  , &
                      Q(I,J,L)  , &
                   QLCN(I,J,L)  , &
                   QICN(I,J,L)  , &
                   CLCN(I,J,L)  , &
                  NACTL(I,J,L)  , &
                  NACTI(I,J,L)  , &
                   QST3(I,J,L)  )
             SUBLC(I,J,L) = ( Q(I,J,L) - SUBLC(I,J,L) ) / DT_MOIST
             endif
             endif
           ! cleanup clouds after cldmacro
             call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), &
                                                     QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L), & 
                                                     REMOVE_CLOUDS=(L < KLID) )
           end do ! IM loop
         end do ! JM loop
       end do ! LM loop

        ! Clean up memory allocations
        deallocate(facEIS_2d, minrhcrit_2d, turnrhcrit_2d)

! Get fill negative export pointers if requested
! ----------------------------------------------                      
    call MAPL_GetPointer(EXPORT,   DQVDT_FILL,   'DQVDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   DQRDT_FILL,   'DQRDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   DQSDT_FILL,   'DQSDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,   DQGDT_FILL,   'DQGDT_FILL_CLDMACRO', RC=STATUS); VERIFY_(STATUS)
! Cleanup negative water species
! ------------------------------
    call FILLQ2ZERO( Q       , MASS, DT=DT_MOIST, DQDT=  DQVDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QLLS    , MASS, DT=DT_MOIST, DQDT=DQLLSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QLCN    , MASS, DT=DT_MOIST, DQDT=DQLCNDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QILS    , MASS, DT=DT_MOIST, DQDT=DQILSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QICN    , MASS, DT=DT_MOIST, DQDT=DQICNDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
    call FILLQ2ZERO( QRAIN   , MASS, DT=DT_MOIST, DQDT=  DQRDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)  
    call FILLQ2ZERO( QSNOW   , MASS, DT=DT_MOIST, DQDT=  DQSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)  
    call FILLQ2ZERO( QGRAUPEL, MASS, DT=DT_MOIST, DQDT=  DQGDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)  

    ! Update macrophysics tendencies
    !$OMP parallel do default(none) &
    !$OMP shared(LM, JM, IM, DUDT_macro, DVDT_macro, DTDT_macro, DQVDT_macro, &
    !$OMP        DQLDT_macro, DQIDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, &
    !$OMP        DQGDT_macro, U, V, T, Q, QLCN, QLLS, QICN, QILS, CLCN, CLLS, &
    !$OMP        QRAIN, QSNOW, QGRAUPEL, DT_MOIST) &
    !$OMP private(I, J, L)
    do L = 1, LM
       do J = 1, JM
          do I = 1, IM
             DUDT_macro(I,J,L)  = (U(I,J,L)         - DUDT_macro(I,J,L))  / DT_MOIST
             DVDT_macro(I,J,L)  = (V(I,J,L)         - DVDT_macro(I,J,L))  / DT_MOIST
             DTDT_macro(I,J,L)  = (T(I,J,L)         - DTDT_macro(I,J,L))  / DT_MOIST
             DQVDT_macro(I,J,L) = (Q(I,J,L)         - DQVDT_macro(I,J,L)) / DT_MOIST
             DQLDT_macro(I,J,L) = ((QLCN(I,J,L) + QLLS(I,J,L)) - DQLDT_macro(I,J,L)) / DT_MOIST
             DQIDT_macro(I,J,L) = ((QICN(I,J,L) + QILS(I,J,L)) - DQIDT_macro(I,J,L)) / DT_MOIST
             DQADT_macro(I,J,L) = ((CLCN(I,J,L) + CLLS(I,J,L)) - DQADT_macro(I,J,L)) / DT_MOIST
             DQRDT_macro(I,J,L) = (QRAIN(I,J,L)     - DQRDT_macro(I,J,L)) / DT_MOIST
             DQSDT_macro(I,J,L) = (QSNOW(I,J,L)     - DQSDT_macro(I,J,L)) / DT_MOIST
             DQGDT_macro(I,J,L) = (QGRAUPEL(I,J,L)  - DQGDT_macro(I,J,L)) / DT_MOIST
          enddo
       enddo
    enddo
    call MAPL_TimerOff(MAPL,"---CLDMACRO")

    if (DEBUG_TQ_ERRORS) then
         do L = 1, LM
           do J = 1, JM
             do I = 1, IM
             if (  (       T(I,J,L) > 333.0) .OR. (       T(I,J,L) /=        T(I,J,L)) .OR. &
                   (       Q(I,J,L) < 0.0  ) .OR. (       Q(I,J,L) /=        Q(I,J,L)) .OR. &
                   (    QLLS(I,J,L) < 0.0  ) .OR. (    QLLS(I,J,L) /=     QLLS(I,J,L)) .OR. &
                   (    QLCN(I,J,L) < 0.0  ) .OR. (    QLCN(I,J,L) /=     QLCN(I,J,L)) .OR. &
                   (    QILS(I,J,L) < 0.0  ) .OR. (    QILS(I,J,L) /=     QILS(I,J,L)) .OR. &
                   (    QICN(I,J,L) < 0.0  ) .OR. (    QICN(I,J,L) /=     QICN(I,J,L)) .OR. &
                   (   QRAIN(I,J,L) < 0.0  ) .OR. (   QRAIN(I,J,L) /=    QRAIN(I,J,L)) .OR. &
                   (   QSNOW(I,J,L) < 0.0  ) .OR. (   QSNOW(I,J,L) /=    QSNOW(I,J,L)) .OR. &
                   (QGRAUPEL(I,J,L) < 0.0  ) .OR. (QGRAUPEL(I,J,L) /= QGRAUPEL(I,J,L)) ) then
                 print *, "T or Q  spike detected : ", T(I,J,L)
                 print *, "    GFDL DTDT_macro Temp Increment : ", DTDT_macro(I,J,L) * DT_MOIST
                 print *, "  Latitude       =", LATS(I,J)*180.0/MAPL_PI
                 print *, "  Longitude      =", LONS(I,J)*180.0/MAPL_PI
                 print *, "  Pressure (mb)  =", PLmb(I,J,L)
                 print *, "                        CLLS=",   CLLS(I,J,L), "   CLCN=",    CLCN(I,J,L)
                 print *, "    QV=",   Q(I,J,L), " QLLS=",   QLLS(I,J,L), "   QLCN=",    QLCN(I,J,L)
                 print *, "                        QILS=",   QILS(I,J,L), "   QICN=",    QICN(I,J,L)
                 print *, "    QR=", QRAIN(I,J,L), "   QS=",  QSNOW(I,J,L), "     QG=", QGRAUPEL(I,J,L)               
               endif
             enddo
           enddo
         enddo
    endif

    call MAPL_TimerOn(MAPL,"---CLDMICRO")
    ! Zero-out microphysics tendencies
    call MAPL_GetPointer(EXPORT, DQVDT_micro, 'DQVDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDT_micro, 'DQIDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDT_micro, 'DQLDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQADT_micro, 'DQADT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQRDT_micro, 'DQRDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQSDT_micro, 'DQSDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQGDT_micro, 'DQGDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DUDT_micro,  'DUDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DVDT_micro,  'DVDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DTDT_micro,  'DTDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

    !$OMP parallel do default(none) &
    !$OMP shared(LM, JM, IM, Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, &
    !$OMP        CLLS, CLCN, U, V, T, DZET, &
    !$OMP        DQVDT_micro, DQLDT_micro, DQIDT_micro, DQRDT_micro, DQSDT_micro, &
    !$OMP        DQGDT_micro, DQADT_micro, DUDT_micro, DVDT_micro, DTDT_micro, &
    !$OMP        DZ, DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, DQSDTmic, DQGDTmic, &
    !$OMP        DQADTmic, DUDTmic, DVDTmic, DTDTmic, PFI_LS, PFL_LS, &
    !$OMP        RAD_CF, RAD_QL, RAD_QI, RAD_QV, RAD_QR, RAD_QS, RAD_QG) &
    !$OMP private(I, J, L)
    do L = 1, LM
       do J = 1, JM
          do I = 1, IM
             DQVDT_micro(I,J,L) = Q(I,J,L)
             DQLDT_micro(I,J,L) = QLLS(I,J,L) + QLCN(I,J,L)
             DQIDT_micro(I,J,L) = QILS(I,J,L) + QICN(I,J,L)
             DQRDT_micro(I,J,L) = QRAIN(I,J,L)
             DQSDT_micro(I,J,L) = QSNOW(I,J,L)
             DQGDT_micro(I,J,L) = QGRAUPEL(I,J,L)
             DQADT_micro(I,J,L) = CLLS(I,J,L) + CLCN(I,J,L)
             DUDT_micro(I,J,L)  = U(I,J,L)
             DVDT_micro(I,J,L)  = V(I,J,L)
             DTDT_micro(I,J,L)  = T(I,J,L)

             ! Delta-Z layer thickness (gfdl expects this to be negative)
             DZ(I,J,L) = -1.0 * DZET(I,J,L)

             ! Zero-out local microphysics tendencies
             DQVDTmic(I,J,L) = 0.0
             DQLDTmic(I,J,L) = 0.0
             DQRDTmic(I,J,L) = 0.0
             DQIDTmic(I,J,L) = 0.0
             DQSDTmic(I,J,L) = 0.0
             DQGDTmic(I,J,L) = 0.0
             DQADTmic(I,J,L) = 0.0
             DUDTmic(I,J,L)  = 0.0
             DVDTmic(I,J,L)  = 0.0
             DTDTmic(I,J,L)  = 0.0

             ! Zero-out 3D Precipitation Fluxes
             PFI_LS(I,J,L) = 0.0
             PFL_LS(I,J,L) = 0.0

             ! Cloud fractions and condensates for radiation
             RAD_CF(I,J,L) = MIN(CLCN(I,J,L) + CLLS(I,J,L), 1.0)
             RAD_QL(I,J,L) = QLCN(I,J,L) + QLLS(I,J,L)
             RAD_QI(I,J,L) = QICN(I,J,L) + QILS(I,J,L)
             RAD_QV(I,J,L) = Q(I,J,L)
             RAD_QR(I,J,L) = QRAIN(I,J,L)
             RAD_QS(I,J,L) = QSNOW(I,J,L)
             RAD_QG(I,J,L) = QGRAUPEL(I,J,L)
          enddo
       enddo
    enddo

        ! Run the driver
         if (GFDL_MP3) then
         call gfdl_mp_driver( &
                             ! Input water/cloud species and liquid+ice CCN NACTL & NACTI (#/m^3)
                               RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, DBZ3D, NACTL, NACTI, &
                             ! Input fields
                               T, W, U, V, DZ, DP, &
                             ! Other inputs
                               DT_MOIST, RHCRIT3D, PHIS, CNV_FRC, EIS, AREA, SRF_TYPE, &
                             ! Output precipitates
                               PRCP_WATER, PRCP_RAIN, PRCP_ICE, PRCP_SNOW, PRCP_GRAUPEL, &
                             ! constant grid/time information
                               LHYDROSTATIC, 1, IM*JM, 1,LM, KLID, &
                             ! Output tendencies
                               DQADTmic, &
                             ! Output rain re-evaporation and sublimation
                               REV_LS, RSU_LS, &
                             ! Output mass flux during sedimentation (Pa kg/kg)
                               PFL_LS(:,:,1:LM), PFR_LS(:,:,1:LM), PFI_LS(:,:,1:LM), PFS_LS(:,:,1:LM), PFG_LS(:,:,1:LM) )
           ! Convert evap/subl/cloud/precipitation flux exports from (mm/day) to (kg m-2 s-1)
           !$OMP parallel do default(none) &
           !$OMP shared(LM, JM, IM, REV_LS, RSU_LS, PFL_LS, PFI_LS, PFR_LS, PFS_LS, PFG_LS) &
           !$OMP private(I, J, L)
           do L = 1, LM
              do J = 1, JM
                 do I = 1, IM
                    REV_LS(I,J,L) = REV_LS(I,J,L) * (1.0 / 86400.0)
                    RSU_LS(I,J,L) = RSU_LS(I,J,L) * (1.0 / 86400.0)
                    PFL_LS(I,J,L) = PFL_LS(I,J,L) * (1.0 / 86400.0)
                    PFI_LS(I,J,L) = PFI_LS(I,J,L) * (1.0 / 86400.0)
                    PFR_LS(I,J,L) = PFR_LS(I,J,L) * (1.0 / 86400.0)
                    PFS_LS(I,J,L) = PFS_LS(I,J,L) * (1.0 / 86400.0)
                    PFG_LS(I,J,L) = PFG_LS(I,J,L) * (1.0 / 86400.0)
                 enddo
              enddo
           enddo
           if (do_ref) then
               call MAPL_TimerOn(MAPL,"---CLD_REF_DBZ")
               call MAPL_GetPointer(EXPORT, DBZ     , 'REF_DBZ'     , RC=STATUS); VERIFY_(STATUS)
               call MAPL_GetPointer(EXPORT, DBZ_MAX , 'REF_DBZ_MAX' , RC=STATUS); VERIFY_(STATUS)
               call MAPL_GetPointer(EXPORT, DBZ_1KM , 'REF_DBZ_1KM' , RC=STATUS); VERIFY_(STATUS)
               call MAPL_GetPointer(EXPORT, DBZ_TOP , 'REF_DBZ_TOP' , RC=STATUS); VERIFY_(STATUS)
               call MAPL_GetPointer(EXPORT, DBZ_M10C, 'REF_DBZ_M10C', RC=STATUS); VERIFY_(STATUS)
               if (associated(DBZ)) DBZ = DBZ3D
               if (associated(DBZ_MAX)) then
                   DBZ_MAX=-9999.0
                   DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                      DBZ_MAX(I,J) = MAX(DBZ_MAX(I,J),DBZ3D(I,J,L))
                   END DO ; END DO ; END DO
               endif
               if (associated(DBZ_1KM)) then
                   call cs_interpolator(1, IM, 1, JM, LM, DBZ3D, 1000., ZLE0, DBZ_1KM, -20.)
               endif    
               if (associated(DBZ_TOP)) then
                   DBZ_TOP=MAPL_UNDEF
                   DO J=1,JM ; DO I=1,IM
                      DO L=LM,1,-1   
                         if (ZLE0(i,j,l) >= 25000.) continue
                         if (DBZ3D(i,j,l) >= 18.5 ) then
                             DBZ_TOP(I,J) = ZLE0(I,J,L)
                             exit
                         endif       
                      END DO
                   END DO ; END DO
               endif  
               if (associated(DBZ_M10C)) then
                   DBZ_M10C=MAPL_UNDEF
                   DO J=1,JM ; DO I=1,IM
                      DO L=LM,1,-1
                         if (ZLE0(i,j,l) >= 25000.) continue
                         if (T(i,j,l) <= MAPL_TICE-10.0) then
                             DBZ_M10C(I,J) = DBZ3D(I,J,L)
                             exit
                         endif
                      END DO
                   END DO ; END DO
               endif
               call MAPL_TimerOff(MAPL,"---CLD_REF_DBZ")
           endif
         else
         call gfdl_cloud_microphys_driver( &
                             ! Input water/cloud species and liquid+ice CCN NACTL & NACTI (#/m^3)
                               RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTL, NACTI, &
                             ! Output tendencies
                               DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                               DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
                             ! Input fields
                               T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
                             ! constant inputs
                               AREA, DT_MOIST, FRLAND, CNV_FRC, SRF_TYPE, EIS, &
                               RHCRIT3D, ANV_ICEFALL, LS_ICEFALL, &
                             ! Output rain re-evaporation and sublimation
                               REV_LS, RSU_LS, &
                             ! Output fall speeds
                               VFALL_ICE, VFALL_SNOW, VFALL_GRAUPEL, VFALL_RAIN, &
                             ! Output precipitates
                               PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
                             ! Output mass flux during sedimentation (Pa kg/kg)
                               PFL_LS(:,:,1:LM), PFI_LS(:,:,1:LM), &
                             ! constant grid/time information
                               LHYDROSTATIC, LPHYS_HYDROSTATIC, &
                               1,IM, 1,JM, 1,LM, KLID, LM)
       ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
           PFL_LS = PFL_LS/(MAPL_GRAV*DT_MOIST)
           PFI_LS = PFI_LS/(MAPL_GRAV*DT_MOIST)
           PFR_LS = 0.0
           PFS_LS = 0.0
           PFG_LS = 0.0
         T = T + DTDTmic * DT_MOIST
         U = U + DUDTmic * DT_MOIST
         V = V + DVDTmic * DT_MOIST
     ! Apply moist/cloud species tendencies
         RAD_QV = RAD_QV + DQVDTmic * DT_MOIST
         RAD_QL = RAD_QL + DQLDTmic * DT_MOIST
         RAD_QR = RAD_QR + DQRDTmic * DT_MOIST
         RAD_QI = RAD_QI + DQIDTmic * DT_MOIST
         RAD_QS = RAD_QS + DQSDTmic * DT_MOIST
         RAD_QG = RAD_QG + DQGDTmic * DT_MOIST
         endif
     ! CleanUp Negative Water Species
         if (REPORT_GFDL_1M_NEGATIVES) then
           call neg_adj_external(IM, JM, LM, RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, &
                                 WARNING_LABEL="After GFDL_1M Driver", VM=VMG, RC=STATUS); VERIFY_(STATUS)
         else
           call neg_adj_external(IM, JM, LM, RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, &
                                 VM=VMG, RC=STATUS); VERIFY_(STATUS)
         endif

     ! Update cloud fraction, redistribute clouds, and fill vapor/precip states
     !$OMP parallel do default(none) &
     !$OMP shared(LM, JM, IM, RAD_QL, RAD_QI, RAD_CF, DQADTmic, DT_MOIST, &
     !$OMP        CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, T, &
     !$OMP        Q, QRAIN, RAD_QR, QSNOW, RAD_QS, QGRAUPEL, RAD_QG) &
     !$OMP private(I, J, L)
     do L = 1, LM
        do J = 1, JM
           do I = 1, IM
              ! 1. Update cloud fraction (Replaces where/elsewhere)
              if (RAD_QL(I,J,L) + RAD_QI(I,J,L) > 0.0) then
                 RAD_CF(I,J,L) = MIN(1.0, MAX(0.0, RAD_CF(I,J,L) + DQADTmic(I,J,L) * DT_MOIST))
              else
                 RAD_CF(I,J,L) = 0.0
              endif

              ! 2. Redistribute clouds (Calls the new scalar subroutine)
              call REDISTRIBUTE_CLOUDS_SCALAR(RAD_CF(I,J,L), RAD_QL(I,J,L), RAD_QI(I,J,L), &
                                       CLCN(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QLLS(I,J,L), &
                                       QICN(I,J,L), QILS(I,J,L), RAD_QV(I,J,L), T(I,J,L))

              ! 3. Fill vapor/rain/snow/graupel state (Replaces whole-array copies)
              Q(I,J,L)        = RAD_QV(I,J,L)
              QRAIN(I,J,L)    = RAD_QR(I,J,L)
              QSNOW(I,J,L)    = RAD_QS(I,J,L)
              QGRAUPEL(I,J,L) = RAD_QG(I,J,L)
           enddo
        enddo
     enddo

     ! Get fill negative export pointers if requested
     ! ----------------------------------------------                      
         call MAPL_GetPointer(EXPORT,   DQVDT_FILL,   'DQVDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT,   DQRDT_FILL,   'DQRDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT,   DQSDT_FILL,   'DQSDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT,   DQGDT_FILL,   'DQGDT_FILL_CLDMICRO', RC=STATUS); VERIFY_(STATUS)
     ! Cleanup negative water species
     ! ------------------------------
         if (REPORT_GFDL_1M_NEGATIVES) then
           call FILLQ2ZERO( Q       , MASS, DT=DT_MOIST, DQDT=  DQVDT_FILL, WARNING_LABEL="QV   After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QLLS    , MASS, DT=DT_MOIST, DQDT=DQLLSDT_FILL, WARNING_LABEL="QLLS After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QLCN    , MASS, DT=DT_MOIST, DQDT=DQLCNDT_FILL, WARNING_LABEL="QLCN After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QILS    , MASS, DT=DT_MOIST, DQDT=DQILSDT_FILL, WARNING_LABEL="QILS After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QICN    , MASS, DT=DT_MOIST, DQDT=DQICNDT_FILL, WARNING_LABEL="QICN After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QRAIN   , MASS, DT=DT_MOIST, DQDT=  DQRDT_FILL, WARNING_LABEL="QR   After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QSNOW   , MASS, DT=DT_MOIST, DQDT=  DQSDT_FILL, WARNING_LABEL="QS   After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QGRAUPEL, MASS, DT=DT_MOIST, DQDT=  DQGDT_FILL, WARNING_LABEL="QG   After GFDL_1M REDISTRIBUTE_CLOUDS", VM=VMG, RC=STATUS); VERIFY_(STATUS)
         else
           call FILLQ2ZERO( Q       , MASS, DT=DT_MOIST, DQDT=  DQVDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QLLS    , MASS, DT=DT_MOIST, DQDT=DQLLSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QLCN    , MASS, DT=DT_MOIST, DQDT=DQLCNDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QILS    , MASS, DT=DT_MOIST, DQDT=DQILSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QICN    , MASS, DT=DT_MOIST, DQDT=DQICNDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QRAIN   , MASS, DT=DT_MOIST, DQDT=  DQRDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QSNOW   , MASS, DT=DT_MOIST, DQDT=  DQSDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
           call FILLQ2ZERO( QGRAUPEL, MASS, DT=DT_MOIST, DQDT=  DQGDT_FILL, VM=VMG, RC=STATUS); VERIFY_(STATUS)
         endif

     ! -----------------------------------------------------------------
     ! 1. 2D Precip Diagnostics (Threaded + Div to Mult)
     ! -----------------------------------------------------------------
     !$OMP parallel do default(none) &
     !$OMP shared(IM, JM, PRCP_WATER, PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
     !$OMP        LS_PRCP, LS_SNR, ICE, FRZR) &
     !$OMP private(I, J)
     do J = 1, JM
        do I = 1, IM
           PRCP_WATER(I,J)   = MAX(PRCP_WATER(I,J)   * (1.0 / 86400.0), 0.0)
           PRCP_RAIN(I,J)    = MAX(PRCP_RAIN(I,J)    * (1.0 / 86400.0), 0.0)
           PRCP_SNOW(I,J)    = MAX(PRCP_SNOW(I,J)    * (1.0 / 86400.0), 0.0)
           PRCP_ICE(I,J)     = MAX(PRCP_ICE(I,J)     * (1.0 / 86400.0), 0.0)
           PRCP_GRAUPEL(I,J) = MAX(PRCP_GRAUPEL(I,J) * (1.0 / 86400.0), 0.0)

           LS_PRCP(I,J) = PRCP_RAIN(I,J)
           LS_SNR(I,J)  = PRCP_SNOW(I,J)
           ICE(I,J)     = PRCP_ICE(I,J) + PRCP_GRAUPEL(I,J)
           FRZR(I,J)    = 0.0
        enddo
     enddo

     ! -----------------------------------------------------------------
     ! 2. Fused Flux Redistribution, MeltFreeze, and RadCouple
     ! -----------------------------------------------------------------
     !$OMP parallel do default(none) &
     !$OMP shared(LM, JM, IM, QLCN, RAD_QL, PFL_AN, PFL_LS, PFR_LS, QICN, RAD_QI, &
     !$OMP        PFI_AN, PFI_LS, PFS_LS, PFG_LS, LMELTFRZ_CLDMICRO, DT_MOIST, &
     !$OMP        CNV_FRC, SRF_TYPE, T, QLLS, QILS, Q, CLLS, CLCN, KLID, PLmb, &
     !$OMP        QRAIN, QSNOW, QGRAUPEL, NACTL, NACTI, RAD_QV, RAD_QR, RAD_QS, &
     !$OMP        RAD_QG, RAD_CF, CLDREFFL, CLDREFFI, FAC_RL, MIN_RL, MAX_RL, &
     !$OMP        FAC_RI, MIN_RI, MAX_RI, AREA) &
     !$OMP private(I, J, L, tmp_val, one_minus_sigma)
     do L = 1, LM
       do J = 1, JM
         do I = 1, IM
           ! Redistribute precipitation fluxes for chemistry (TMP3D eliminated)
           tmp_val = MIN(1.0, MAX(QLCN(I,J,L) / MAX(RAD_QL(I,J,L), 1.E-8), 0.0))
           PFL_AN(I,J,L) = (PFL_LS(I,J,L) + PFR_LS(I,J,L)) * tmp_val
           PFL_LS(I,J,L) = (PFL_LS(I,J,L) + PFR_LS(I,J,L)) - PFL_AN(I,J,L)
           
           tmp_val = MIN(1.0, MAX(QICN(I,J,L) / MAX(RAD_QI(I,J,L), 1.E-8), 0.0))
           PFI_AN(I,J,L) = (PFI_LS(I,J,L) + PFS_LS(I,J,L) + PFG_LS(I,J,L)) * tmp_val
           PFI_LS(I,J,L) = (PFI_LS(I,J,L) + PFS_LS(I,J,L) + PFG_LS(I,J,L)) - PFI_AN(I,J,L)

           ! MeltFreeze and FixUp
           if (LMELTFRZ_CLDMICRO) then
             call MELTFRZ(DT_MOIST, CNV_FRC(I,J), SRF_TYPE(I,J), T(I,J,L), QLCN(I,J,L), QICN(I,J,L))    
             call MELTFRZ(DT_MOIST, CNV_FRC(I,J), SRF_TYPE(I,J), T(I,J,L), QLLS(I,J,L), QILS(I,J,L))
             call FIX_UP_CLOUDS(Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), &
                                QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L), REMOVE_CLOUDS=(L < KLID))
           endif

           ! Get radiative properties (Scale-Aware)
           one_minus_sigma = 1.0 - SIGMA(sqrt(AREA(I,J)))
           call RADCOUPLE_SCALE_AWARE(T(I,J,L), PLmb(I,J,L), CLLS(I,J,L), CLCN(I,J,L), &
                 Q(I,J,L), QLLS(I,J,L), QILS(I,J,L), QLCN(I,J,L), QICN(I,J,L), &
                 QRAIN(I,J,L), QSNOW(I,J,L), QGRAUPEL(I,J,L), NACTL(I,J,L), NACTI(I,J,L), &
                 one_minus_sigma, &
                 RAD_QV(I,J,L), RAD_QL(I,J,L), RAD_QI(I,J,L), RAD_QR(I,J,L), RAD_QS(I,J,L), &
                 RAD_QG(I,J,L), RAD_CF(I,J,L), CLDREFFL(I,J,L), CLDREFFI(I,J,L), &
                 FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
         enddo
       enddo
     enddo

     ! -----------------------------------------------------------------
     ! 3. Cleanup negative water species (Sequential, DO NOT THREAD)
     ! -----------------------------------------------------------------
     if (REPORT_GFDL_1M_NEGATIVES) then
       call FILLQ2ZERO( Q       , MASS, WARNING_LABEL="QV   After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QLLS    , MASS, WARNING_LABEL="QLLS After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QLCN    , MASS, WARNING_LABEL="QLCN After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QILS    , MASS, WARNING_LABEL="QILS After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QICN    , MASS, WARNING_LABEL="QICN After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QRAIN   , MASS, WARNING_LABEL="QR   After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QSNOW   , MASS, WARNING_LABEL="QS   After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QGRAUPEL, MASS, WARNING_LABEL="QG   After GFDL_1M RADCOUPLE", VM=VMG, RC=STATUS); VERIFY_(STATUS)
     else
       call FILLQ2ZERO( Q       , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QLLS    , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QLCN    , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QILS    , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QICN    , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QRAIN   , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QSNOW   , MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
       call FILLQ2ZERO( QGRAUPEL, MASS, VM=VMG, RC=STATUS); VERIFY_(STATUS)
     endif

     ! -----------------------------------------------------------------
     ! 4. Update microphysics tendencies (Threaded array assignments)
     ! -----------------------------------------------------------------
     !$OMP parallel do default(none) &
     !$OMP shared(LM, JM, IM, DQVDT_micro, DQLDT_micro, DQIDT_micro, DQADT_micro, &
     !$OMP        DQRDT_micro, DQSDT_micro, DQGDT_micro, DUDT_micro, DVDT_micro, &
     !$OMP        DTDT_micro, Q, QLLS, QLCN, QILS, QICN, CLLS, CLCN, QRAIN, &
     !$OMP        QSNOW, QGRAUPEL, U, V, T, DT_MOIST) &
     !$OMP private(I, J, L)
     do L = 1, LM
        do J = 1, JM
           do I = 1, IM
              DQVDT_micro(I,J,L) = (Q(I,J,L)          - DQVDT_micro(I,J,L)) / DT_MOIST
              DQLDT_micro(I,J,L) = ((QLLS(I,J,L) + QLCN(I,J,L)) - DQLDT_micro(I,J,L)) / DT_MOIST
              DQIDT_micro(I,J,L) = ((QILS(I,J,L) + QICN(I,J,L)) - DQIDT_micro(I,J,L)) / DT_MOIST
              DQADT_micro(I,J,L) = ((CLLS(I,J,L) + CLCN(I,J,L)) - DQADT_micro(I,J,L)) / DT_MOIST
              DQRDT_micro(I,J,L) = (QRAIN(I,J,L)      - DQRDT_micro(I,J,L)) / DT_MOIST
              DQSDT_micro(I,J,L) = (QSNOW(I,J,L)      - DQSDT_micro(I,J,L)) / DT_MOIST
              DQGDT_micro(I,J,L) = (QGRAUPEL(I,J,L)   - DQGDT_micro(I,J,L)) / DT_MOIST
              DUDT_micro(I,J,L)  = (U(I,J,L)          - DUDT_micro(I,J,L))  / DT_MOIST
              DVDT_micro(I,J,L)  = (V(I,J,L)          - DVDT_micro(I,J,L))  / DT_MOIST
              DTDT_micro(I,J,L)  = (T(I,J,L)          - DTDT_micro(I,J,L))  / DT_MOIST
           enddo
        enddo
     enddo
        call MAPL_TimerOff(MAPL,"---CLDMICRO")

    if (DEBUG_TQ_ERRORS) then
         do L = 1, LM
           do J = 1, JM
             do I = 1, IM
             if (  (       T(I,J,L) > 333.0) .OR. (       T(I,J,L) /=        T(I,J,L)) .OR. &
                   (       Q(I,J,L) < 0.0  ) .OR. (       Q(I,J,L) /=        Q(I,J,L)) .OR. &
                   (    QLLS(I,J,L) < 0.0  ) .OR. (    QLLS(I,J,L) /=     QLLS(I,J,L)) .OR. &
                   (    QLCN(I,J,L) < 0.0  ) .OR. (    QLCN(I,J,L) /=     QLCN(I,J,L)) .OR. &
                   (    QILS(I,J,L) < 0.0  ) .OR. (    QILS(I,J,L) /=     QILS(I,J,L)) .OR. &
                   (    QICN(I,J,L) < 0.0  ) .OR. (    QICN(I,J,L) /=     QICN(I,J,L)) .OR. &
                   (   QRAIN(I,J,L) < 0.0  ) .OR. (   QRAIN(I,J,L) /=    QRAIN(I,J,L)) .OR. &
                   (   QSNOW(I,J,L) < 0.0  ) .OR. (   QSNOW(I,J,L) /=    QSNOW(I,J,L)) .OR. &
                   (QGRAUPEL(I,J,L) < 0.0  ) .OR. (QGRAUPEL(I,J,L) /= QGRAUPEL(I,J,L)) ) then
                 print *, "T or Q  spike detected : ", T(I,J,L)
                 print *, "    GFDL DTDT_micro Temp Increment : ", DTDT_micro(I,J,L) * DT_MOIST
                 print *, "  Latitude       =", LATS(I,J)*180.0/MAPL_PI
                 print *, "  Longitude      =", LONS(I,J)*180.0/MAPL_PI
                 print *, "  Pressure (mb)  =", PLmb(I,J,L)
                 print *, "                        CLLS=",   CLLS(I,J,L), "   CLCN=",    CLCN(I,J,L)
                 print *, "                      RAD_QL=", RAD_QL(I,J,L), " RAD_QI=",  RAD_QI(I,J,L)
                 print *, "    QV=",   Q(I,J,L), " QLLS=",   QLLS(I,J,L), "   QLCN=",    QLCN(I,J,L)
                 print *, "                        QILS=",   QILS(I,J,L), "   QICN=",    QICN(I,J,L)
                 print *, "    QR=",QRAIN(I,J,L),"   QS=",  QSNOW(I,J,L), "     QG=", QGRAUPEL(I,J,L)
               endif
             enddo
           enddo
         enddo
    endif

        call MAPL_TimerOn(MAPL,"---CLDDIAGS")

        call MAPL_GetPointer(EXPORT, PTR3D, 'DQRL', RC=STATUS); VERIFY_(STATUS)
        if(associated(PTR3D)) PTR3D = DQRDT_macro + DQRDT_micro

        ! dissipative heating tendency from KE across the macro/micro physics
        call MAPL_GetPointer(EXPORT, PTR3D, 'DTDTFRIC', RC=STATUS); VERIFY_(STATUS)
        if(associated(PTR3D)) then
          call dissipative_ke_heating(IM,JM,LM, MASS,U0,V0, &
                                      DUDT_macro+DUDT_micro,&
                                      DVDT_macro+DVDT_micro,PTR3D)
        endif

        ! Compute DBZ radar reflectivity
        call ESMF_ClockGetAlarm(clock, 'DBZ_RunAlarm', alarm, RC=STATUS); VERIFY_(STATUS)
        alarm_is_ringing = ESMF_AlarmIsRinging(alarm, RC=STATUS); VERIFY_(STATUS)
            
        call MAPL_GetPointer(EXPORT,  NACTR,  'NACTR',        RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT,  PTR2D,  'REFL10CM_MAX', RC=STATUS); VERIFY_(STATUS)
    
        ! 1. If the user explicitly requested NACTR export, fill it every time (or whenever needed)
        if (associated(NACTR)) then
            NACTR = 1.e8 * QRAIN**0.8 
        endif

        ! 2. Handle the reflectivity alarm
        if (alarm_is_ringing) then
           call ESMF_AlarmRingerOff(alarm, RC=STATUS); VERIFY_(STATUS)
           
           ! Only compute if the user actually requested the reflectivity output
           if (associated(PTR2D)) then 
               call MAPL_TimerOn(MAPL,"---CLD_REFL10CM")    
               rand1 = 0.0
               TMP3D = 0.0
               
               ! If NACTR wasn't associated, we still need it for calc_refl10cm!
               ! We can use TMP3D to temporarily hold NACTR if needed, or if calc_refl10cm 
               ! requires it as a distinct array, use a locally allocated TMP_NACTR array.
               ! Assuming TMP_NACTR is an allocatable 3D array defined at the top:
               
               if (.not. associated(NACTR)) then
                   ! Fill a local temporary array to pass into the subroutine
                   ALLOCATE ( TMP_NACTR(IM,JM,LM) )
                   TMP_NACTR = 1.e8 * QRAIN**0.8
               endif
               
               DO J=1,JM ; DO I=1,IM
                 ! Pass either the Export pointer (if associated) or the local temporary array
                 if (associated(NACTR)) then
                     call calc_refl10cm(Q(I,J,:), QRAIN(I,J,:), NACTR(I,J,:), QSNOW(I,J,:), QGRAUPEL(I,J,:), &
                        T(I,J,:), 100*PLmb(I,J,:), TMP3D(I,J,:), rand1, 1, LM, I, J) 
                 else
                     call calc_refl10cm(Q(I,J,:), QRAIN(I,J,:), TMP_NACTR(I,J,:), QSNOW(I,J,:), QGRAUPEL(I,J,:), &
                        T(I,J,:), 100*PLmb(I,J,:), TMP3D(I,J,:), rand1, 1, LM, I, J) 
                 endif
               END DO ; END DO
              
               if (.not. associated(NACTR)) then
                   DEALLOCATE ( TMP_NACTR )
               endif
 
               PTR2D = -9999.0
               DO L=1,LM ; DO J=1,JM ; DO I=1,IM 
                  PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
               END DO ; END DO ; END DO
               
               call MAPL_TimerOff(MAPL,"---CLD_REFL10CM")
           endif
        endif

        call MAPL_TimerOn(MAPL,"---CLD_CALCDBZ")
        call MAPL_GetPointer(EXPORT, DBZ     , 'DBZ'     , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_MAX , 'DBZ_MAX' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_1KM , 'DBZ_1KM' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_TOP , 'DBZ_TOP' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_M10C, 'DBZ_M10C', RC=STATUS); VERIFY_(STATUS)
        if ( (associated(DBZ) .OR. &
              associated(DBZ_MAX) .OR. associated(DBZ_1KM) .OR. associated(DBZ_TOP) .OR. associated(DBZ_M10C)) ) then
            allocate ( qg_col(LM) )
            allocate ( qh_col(LM) ) 
            allocate ( prs_col(LM) ) 
            allocate ( dbz_col(LM) ) 
            !$OMP parallel do default(none) &
            !$OMP shared(IM, JM, LM, W, QGRAUPEL, PLmb, T, Q, QRAIN, QSNOW, &
            !$OMP        DBZ_VAR_INTERCP, LIQUID_SKIN_SNOW, LIQUID_SKIN_GRAUPEL, LIQUID_SKIN_HAIL, DBZ3D) &
            !$OMP private(I, J, L, fraction_hail, qg_col, qh_col, prs_col, dbz_col)
            DO J = 1, JM
                DO I = 1, IM
                    ! 1. Prepare the 1D column data for this specific (I,J) location
                    DO L = 1, LM
                        ! Calculate a fraction between 0.0 and 1.0 based on updraft W
                        fraction_hail = MAX(0.0, MIN(1.0, (W(I,J,L) - W_START) / (W_FULL - W_START)))
                        ! Partition the mass into 1D thread-private columns
                        qh_col(L)  = QGRAUPEL(I,J,L) * fraction_hail
                        qg_col(L)  = QGRAUPEL(I,J,L) * (1.0 - fraction_hail)
                        ! Pre-multiply pressure for the function
                        prs_col(L) = 100.0 * PLmb(I,J,L)
                    END DO
                    ! 2. Call the newly refactored 1D column function
                    !    Note: We pass 1D array slices like T(I,J,:) directly.
                    dbz_col = compute_radar_reflectivity( &
                              PRS = prs_col, &
                              TMK = T(I,J,:), &
                              QVP = Q(I,J,:), &
                              QRAIN = QRAIN(I,J,:), &
                              QSNOW = QSNOW(I,J,:), &
                              QGRAUPEL = qg_col, &
                              QHAIL = qh_col, &
                              disable_variable_intercept_params = (DBZ_VAR_INTERCP == 0), &
                              liqskin_snow = LIQUID_SKIN_SNOW, &
                              liqskin_graupel = LIQUID_SKIN_GRAUPEL, &
                              liqskin_hail = LIQUID_SKIN_HAIL)
                    ! 3. Store the returned column back into the 3D state
                    DO L = 1, LM
                        DBZ3D(I,J,L) = dbz_col(L)
                    END DO
                END DO
            END DO
        end if
        if (associated(DBZ)) DBZ = DBZ3D
        if (associated(DBZ_MAX)) then
           DBZ_MAX=-9999.0
           DO L=1,LM ; DO J=1,JM ; DO I=1,IM
              DBZ_MAX(I,J) = MAX(DBZ_MAX(I,J),DBZ3D(I,J,L))
           END DO ; END DO ; END DO
        endif
        if (associated(DBZ_1KM)) then
           call cs_interpolator(1, IM, 1, JM, LM, DBZ3D, 1000., ZLE0, DBZ_1KM, -20.)
        endif
        if (associated(DBZ_TOP)) then
           DBZ_TOP=MAPL_UNDEF
           DO J=1,JM ; DO I=1,IM
              DO L=LM,1,-1
                 if (ZLE0(i,j,l) >= 25000.) continue
                 if (DBZ3D(i,j,l) >= 18.5 ) then
                     DBZ_TOP(I,J) = ZLE0(I,J,L)
                     exit
                 endif
              END DO
           END DO ; END DO
        endif
        if (associated(DBZ_M10C)) then
           DBZ_M10C=MAPL_UNDEF
           DO J=1,JM ; DO I=1,IM
              DO L=LM,1,-1
                 if (ZLE0(i,j,l) >= 25000.) continue
                 if (T(i,j,l) <= MAPL_TICE-10.0) then
                     DBZ_M10C(I,J) = DBZ3D(I,J,L)
                     exit
                 endif
              END DO
           END DO ; END DO
        endif
        call MAPL_TimerOff(MAPL,"---CLD_CALCDBZ")

        call MAPL_GetPointer(EXPORT, PTR3D, 'QRTOT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = QRAIN

        call MAPL_GetPointer(EXPORT, PTR3D, 'QSTOT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = QSNOW

        call MAPL_GetPointer(EXPORT, PTR3D, 'QGTOT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = QGRAUPEL

        call MAPL_GetPointer(EXPORT, PTR2D, 'LWP', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS+QRAIN ) *MASS , 3 )

        call MAPL_GetPointer(EXPORT, PTR2D, 'IWP', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM( ( QICN+QILS+QSNOW+QGRAUPEL ) *MASS , 3 )

        call MAPL_TimerOff(MAPL,"---CLDDIAGS")

     call MAPL_TimerOff(MAPL,"--GFDL_1M",RC=STATUS)

end subroutine GFDL_1M_Run

end module GEOS_GFDL_1M_InterfaceMod
