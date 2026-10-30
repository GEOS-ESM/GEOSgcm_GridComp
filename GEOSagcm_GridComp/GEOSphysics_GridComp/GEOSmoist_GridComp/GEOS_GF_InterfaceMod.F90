! $Id$

#include "MAPL.h"

!=============================================================================
!BOP

! !MODULE: GEOS_GF_InterfaceMod -- A Module to interface with the
!   GF convection

module GEOS_GF_InterfaceMod

  use ESMF
  use MAPL
  use gFTL2_StringVector, only: StringVector
  use GEOS_UtilsMod
  use GEOSmoist_Process_Library
  use GF2020_2M_MicrophysicsMod
  use Aer_Actv_Single_Moment
  use aer_cloud, only: AeroPropsNew
  use ConvPar_GF_SharedParams
  use ConvPar_GF_GEOS5
  use ConvPar_GF2020
  use moist_dsl_workarounds

  implicit none

  private

  integer                                 :: STATUS

  integer :: USE_GF2020
  logical :: STOCHASTIC_CNV
  real    :: STOCH_TOP, STOCH_BOT
  real    :: SCLM_DEEP
  real    :: GF_MIN_AREA
  logical :: FIX_CNV_CLOUD
  logical :: REPORT_GF_NEGATIVES
  integer :: ZERO_DIFF_TAU
  integer :: ZERO_DIFF_VGRID
  integer :: ZERO_DIFF_HEI_EDT
  integer :: ZERO_DIFF_OTHER
  logical :: USE_PYMOIST_GF2020

  public :: GF_Setup, GF_Initialize, GF_Run, GF_Finalize

contains

subroutine GF_Setup (GC, ADVECTED, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC        ! gridded component
    type(StringVector),  intent(inout) :: ADVECTED  ! internals to subscribe to advection (TRADV)
    integer, optional                  :: RC        ! return code

    call MAPL_GridCompGetResource( GC, 'CONVECTION_TRACER', CONVECTION_TRACER, default= 0, _RC )
    ! MAPL2 meant CNV_TR to be FRIENDLYTO=MOIST:DYNAMICS when CONVECTION_TRACER == 1, but
    ! wrote "MOISTDYNAMICS" (friendly to no component). MAPL3 advects it as intended.
    ! Convective transport by Moist itself (MAPL2 FRIENDLYTO=MOIST) is done in Run.
    if (CONVECTION_TRACER == 1) call ADVECTED%push_back('CNV_TR')

    call MAPL_GridCompAddSpec(GC,                                          &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                     &
         SHORT_NAME       = 'CNV_TR',                                      &
         STANDARD_NAME    = 'convection_tracer',                           &
         UNITS            = '1',                                           &
         DIMS             = 'xyz',                                         &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,                  &
         RESTART_MODE     = MAPL_RESTART_SKIP,                             &
         FILL_VALUE       = 0.0,                                           &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                          &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                     &
         SHORT_NAME       = 'DSL__GF2020_LONS',                            &
         STANDARD_NAME    = 'DSL_longitude',                               &
         UNITS            = 'radians',                                     &
         DIMS             = 'xy',                                          &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_NONE,                    &
         RESTART_MODE     = MAPL_RESTART_SKIP,                             &
         FILL_VALUE       = 0.0,                                           &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                          &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                     &
         SHORT_NAME       = 'DSL__GF2020_LATS',                            &
         STANDARD_NAME    = 'DSL_longitude',                               &
         UNITS            = 'radians',                                     &
         DIMS             = 'xy',                                          &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_NONE,                    &
         RESTART_MODE     = MAPL_RESTART_SKIP,                             &
         FILL_VALUE       = 0.0,                                           &
         _RC)

end subroutine GF_Setup

subroutine GF_Initialize (GC, CLOCK, IMPORT, EXPORT, RC)
    type (ESMF_GridComp), intent(inout) :: GC     ! gridded component
    type (ESMF_State),    intent(inout) :: IMPORT
    type (ESMF_State),    intent(inout) :: EXPORT
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional                   :: RC  ! return code
    integer :: LM
    real                            :: MOIST_DT
    real                            :: GF_DT

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Alarm)        :: GF_RunAlarm
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss
    character(len=:), allocatable :: env_setting  ! GF_ENV_SETTING is a fixed-length ConvPar variable

    call MAPL_GridCompGet(GC, num_levels=LM, _RC)
    ! Moist's timestep (MAPL2: the RUNALARM ring interval)
    call MAPL_ClockGet(CLOCK, dt=MOIST_DT, _RC)
    call MAPL_GridCompGetResource( GC, 'GF_DT', GF_DT, default=MOIST_DT, _RC)

    call ESMF_ClockGet(CLOCK, calendar=calendar, currTime=currentTime, _RC)
    call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, _RC)
    call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, _RC)
    call ESMF_TimeIntervalSet(ringInterval, S=nint(GF_DT), calendar=calendar, _RC)

    GF_RunAlarm = ESMF_AlarmCreate(Clock        = CLOCK,        &
                                   Name         = 'GF_RunAlarm', &
                                   RingInterval = ringInterval, &
                                   RingTime     = currentTime,  &
                                   Enabled      = .true.   ,    &
                                   Sticky       = .false.  , _RC)

    call MAPL_GridCompGetResource( GC, 'USE_PYMOIST_GF2020', USE_PYMOIST_GF2020, default=.FALSE., _RC)

    if (USE_PYMOIST_GF2020) then
      call MAPL_ConfigSetAttribute(CF, MOIST_DT, 'DSL__GF2020_DT', _RC)
      call MAPL_GridCompGetResource( GC, 'GF_ENV_SETTING', env_setting                    ,default= 'DYNAMICS', _RC)
      GF_ENV_SETTING = env_setting
      if (trim(GF_ENV_SETTING)=='CURRENT') then
         call MAPL_ConfigSetAttribute(CF, 0, 'DSL__GF_ENV_SETTING:', _RC)
      elseif (trim(GF_ENV_SETTING)=='DYNAMICS') then
         call MAPL_ConfigSetAttribute(CF, 1, 'DSL__GF_ENV_SETTING:', _RC)
      endif
      call MAPL_pybridge_gcinit( "pyMoist.fortran.param_interfaces.convection.GF2020_interface", MAPL, IMPORT, EXPORT )
    else
    if (LM .eq. 72) then
      call MAPL_GridCompGetResource( GC, 'USE_GF2020', USE_GF2020                            ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_LAND', ZERO_DIFF_LAND                    ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_VVEL', ZERO_DIFF_VVEL                    ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_ENTR', ZERO_DIFF_ENTR                    ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_TAU', ZERO_DIFF_TAU                      ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_VGRID', ZERO_DIFF_VGRID                  ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_HEI_EDT', ZERO_DIFF_HEI_EDT              ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_OTHER', ZERO_DIFF_OTHER                  ,default= 1,    _RC )
    else
      call MAPL_GridCompGetResource( GC, 'USE_GF2020', USE_GF2020                            ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_LAND', ZERO_DIFF_LAND                    ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_VVEL', ZERO_DIFF_VVEL                    ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_ENTR', ZERO_DIFF_ENTR                    ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_TAU', ZERO_DIFF_TAU                      ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_VGRID', ZERO_DIFF_VGRID                  ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_HEI_EDT', ZERO_DIFF_HEI_EDT              ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ZERO_DIFF_OTHER', ZERO_DIFF_OTHER                  ,default= 0,    _RC )
    endif
    call MAPL_GridCompGetResource( GC, 'REPORT_GF_NEGATIVES', REPORT_GF_NEGATIVES, default=.FALSE., _RC)
    IF (USE_GF2020==1) THEN
      call MAPL_GridCompGetResource( GC, 'DEEP', ICUMULUS_GF(DEEP)                           ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'SHALLOW', ICUMULUS_GF(SHAL)                        ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'CONGESTUS', ICUMULUS_GF(MID)                       ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'CLOSURE_DEEP', CLOSURE_CHOICE(DEEP)                ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'CLOSURE_SHALLOW', CLOSURE_CHOICE(SHAL)             ,default= 7,    _RC )
      call MAPL_GridCompGetResource( GC, 'CLOSURE_CONGESTUS', CLOSURE_CHOICE(MID)            ,default= 3,    _RC )
      call MAPL_GridCompGetResource( GC, 'USE_TRACER_TRANSP', USE_TRACER_TRANSP              ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'USE_TRACER_SCAVEN', USE_TRACER_SCAVEN              ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'USE_SCALE_DEP', USE_SCALE_DEP                      ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'USE_MOMENTUM_TRANSP', USE_MOMENTUM_TRANSP          ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'DICYCLE', DICYCLE                                  ,default= 1,    _RC )

      ! -----------------------------------------------------------------------------------------
      ! GRELL-FREITAS ENTRAINMENT, MICROPHYSICS & DETRAINMENT CONTROLS
      ! -----------------------------------------------------------------------------------------


      ! ENTRAINMENT / PLUME DILUTION RATES [m^-1]
      ! Controls the fractional rate at which dry environmental air is mixed into the updraft.
      ! Entrainment acts as the primary thermodynamic brake on convective growth. Higher values 
      ! dilute the plume's latent heat, reduce buoyancy, and force the cloud to terminate at 
      ! lower altitudes. 
      !   - MIN_ENTR_RATE : Absolute physical floor to prevent unphysical adiabatic ascent.
      !   - ENTR_DP : Low rate; allows deep, wide plumes to punch through to the tropopause.
      !   - ENTR_MD : Moderate rate; forces congestus clouds to terminate in the mid-troposphere.
      !   - ENTR_SH : High rate; ensures shallow boundary-layer clouds dilute and mix out rapidly.
      if (INT(ZERO_DIFF_ENTR) == 0) then
        call MAPL_GridCompGetResource( GC, 'MIN_ENTR_RATE', MIN_ENTR_RATE                      ,default= 0.1e-4,_RC )
        call MAPL_GridCompGetResource( GC, 'ENTR_DP', CUM_ENTR_RATE(DEEP)                      ,default= 0.6e-4,_RC )
        call MAPL_GridCompGetResource( GC, 'ENTR_MD', CUM_ENTR_RATE(MID)                       ,default= 2.0e-4,_RC )
        call MAPL_GridCompGetResource( GC, 'ENTR_SH', CUM_ENTR_RATE(SHAL)                      ,default= 6.0e-4,_RC )
      else
        ! Legacy GF settings
        call MAPL_GridCompGetResource( GC, 'MIN_ENTR_RATE', MIN_ENTR_RATE                      ,default= 0.1e-4,_RC )
        call MAPL_GridCompGetResource( GC, 'ENTR_DP', CUM_ENTR_RATE(DEEP)                      ,default= 1.0e-4,_RC )
        call MAPL_GridCompGetResource( GC, 'ENTR_MD', CUM_ENTR_RATE(MID)                       ,default= 9.0e-4,_RC )
        call MAPL_GridCompGetResource( GC, 'ENTR_SH', CUM_ENTR_RATE(SHAL)                      ,default= 1.4e-3,_RC )
      endif

      ! AUTOCONV: Toggles the convective autoconversion scheme.
      !   1 = Classic phase-agnostic Kessler scheme (applies a constant autoconversion rate 
      !       to all total condensate regardless of temperature or phase).
      !   2 = Kessler scheme with temperature dependence (uses Hu et al. phase partitioning
      !       to scale efficiency and protect ice in glaciated sub-freezing updrafts).
      call MAPL_GridCompGetResource( GC, 'AUTOCONV', AUTOCONV                                ,default= 2,     _RC )

      ! C0_*: Internal Core Precipitation Efficiency [m^-1].
      !   Pure microphysics knob. Controls how aggressively the updraft core converts its 
      !   suspended cloud condensate (liquid/ice) into falling precipitation (rain/snow).
      !   Increasing C0 wrings out the plume internally, resulting in thinner, drier anvils aloft.
      !   Decreasing C0 allows the plume to transport more mass to the upper troposphere.
      call MAPL_GridCompGetResource( GC, 'C0_DEEP', C0_DEEP                                  ,default= 2.0e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'C0_MID', C0_MID                                    ,default= 1.0e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'C0_SHAL', C0_SHAL                                  ,default= 0.0   ,_RC )

      ! C0_ICE_EFF: Ice Precipitation Efficiency Multiplier [fraction].
      !   Scales the base C0 autoconversion rate for the ice phase in AUTOCONV=2. 
      !   Decreasing this value protects upper-level ice from precipitating internally, 
      !   forcing it to detrain into the grid and build thicker cirrus anvils (increases TQI).
      call MAPL_GridCompGetResource( GC, 'C0_ICE_EFF', C0_ICE_EFF                            ,default= 0.10  ,_RC )

      ! QRC_CRIT_*: Critical Cloud Liquid Water Threshold [kg/kg].
      !   The updraft must hold this much liquid before Kessler autoconversion is allowed to begin.
      !   Separated by convection type to allow independent tuning of deep vs. mid-level moisture transport.
      call MAPL_GridCompGetResource( GC, 'QRC_CRIT_OCN_DP', QRC_CRIT_OCN_DP                  ,default= 3.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'QRC_CRIT_OCN_MD', QRC_CRIT_OCN_MD                  ,default= 3.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'QRC_CRIT_LND_DP', QRC_CRIT_LND_DP                  ,default= 3.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'QRC_CRIT_LND_MD', QRC_CRIT_LND_MD                  ,default= 3.0e-4,_RC )

      ! C1_*: Lateral Detrainment / Plume Shape Parameter [m^-1].
      !   Macro-physics knob. Completely decoupled from C0 internal microphysics.
      !   Controls the continuous fractional mass shed by the plume into the environment during ascent.
      !   A value > 0.0 activates explicit lateral shedding. A value of 0.0 treats the plume 
      !   as laterally closed, deferring entirely to the bulk detrainment logic.
      call MAPL_GridCompGetResource( GC, 'C1_DEEP', C1_DEEP                                  ,default= 3.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'C1_MID', C1_MID                                    ,default= 1.5e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'C1_SHAL', C1_SHAL                                  ,default= 0.0   ,_RC )

      if (INT(ZERO_DIFF_TAU) == 0) then
         call MAPL_GridCompGetResource( GC, 'GF_MIN_AREA', GF_MIN_AREA                          ,default= 0.0,   _RC )
                                     SGS_W_TIMESCALE = 1.0 ! factor for adjusting GF2020 timescales
         call MAPL_GridCompGetResource( GC, 'SGS_W_TIMESCALE', SGS_W_TIMESCALE                  ,default= (SGS_W_TIMESCALE), _RC )
         ! These are UPPER bounds for new GF timescales
         call MAPL_GridCompGetResource( GC, 'TAU_MID', TAU_MID                                  ,default= 3600., _RC )
         call MAPL_GridCompGetResource( GC, 'TAU_DEEP', TAU_DEEP                                ,default= 5400., _RC )
        ! FADJ_MASSFLX is a fractional mass flux tuning factor (1.0 is no reduction) in low CAPE environments
        call MAPL_GridCompGetResource( GC, 'FADJ_MASSFLX_DP', CUM_FADJ_MASSFLX(DEEP)           ,default= 1.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'FADJ_MASSFLX_SH', CUM_FADJ_MASSFLX(SHAL)           ,default= 0.5,   _RC )
        call MAPL_GridCompGetResource( GC, 'FADJ_MASSFLX_MD', CUM_FADJ_MASSFLX(MID)            ,default= 1.0,   _RC )
      else
        ! Legacy GF settings
        call MAPL_GridCompGetResource( GC, 'GF_MIN_AREA', GF_MIN_AREA                          ,default= 1.e6,  _RC )
                                    SGS_W_TIMESCALE = 0.0 ! do not use new scheme for timescales
        call MAPL_GridCompGetResource( GC, 'SGS_W_TIMESCALE', SGS_W_TIMESCALE                  ,default= (SGS_W_TIMESCALE),_RC )
        ! These are LOWER bounds for classic GF timescales
        call MAPL_GridCompGetResource( GC, 'TAU_MID', TAU_MID                                  ,default= 3600., _RC )
        call MAPL_GridCompGetResource( GC, 'TAU_DEEP', TAU_DEEP                                ,default= 5400., _RC )
        ! FADJ_MASSFLX is a fractional mass flux tuning factor (1.0 is no reduction) in low CAPE environments
        call MAPL_GridCompGetResource( GC, 'FADJ_MASSFLX_DP', CUM_FADJ_MASSFLX(DEEP)           ,default= 1.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'FADJ_MASSFLX_SH', CUM_FADJ_MASSFLX(SHAL)           ,default= 1.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'FADJ_MASSFLX_MD', CUM_FADJ_MASSFLX(MID)            ,default= 1.0,   _RC )
      endif

      if (INT(ZERO_DIFF_VGRID) == 0) then
        call MAPL_GridCompGetResource( GC, 'CLEV_GRID', CLEV_GRID                              ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'VERT_DISCR', VERT_DISCR                            ,default= 1,     _RC )
      else
        ! Legacy GF settings
        call MAPL_GridCompGetResource( GC, 'CLEV_GRID', CLEV_GRID                              ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'VERT_DISCR', VERT_DISCR                            ,default= 0,     _RC )
      endif

      if (INT(ZERO_DIFF_OTHER) == 0) then
        call MAPL_GridCompGetResource( GC, 'USE_INV_LAYERS', USE_INV_LAYERS                    ,default= .FALSE., _RC )
        call MAPL_GridCompGetResource( GC, 'USE_SMOOTH_TEND', USE_SMOOTH_TEND                  ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'SATUR_CALC', SATUR_CALC                            ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'BC_METH', BC_METH                                  ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'USE_REBCB', USE_REBCB                              ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'LAMBAU_DEEP', LAMBAU_DEEP                          ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'LAMBAU_MID', LAMBAU_MID                            ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'MOIST_TRIGGER', MOIST_TRIGGER                      ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'FRAC_MODIS', FRAC_MODIS                            ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'EVAP_FIX', EVAP_FIX                                ,default= 1,     _RC )
        call MAPL_GridCompGetResource( GC, 'AVE_LAYER_DP', CUM_AVE_LAYER(DEEP)                 ,default= 80.,   _RC )
        call MAPL_GridCompGetResource( GC, 'AVE_LAYER_SH', CUM_AVE_LAYER(SHAL)                 ,default= 20.,   _RC )
        call MAPL_GridCompGetResource( GC, 'AVE_LAYER_MD', CUM_AVE_LAYER(MID)                  ,default= 50.,   _RC )
      else
        ! Legacy GF settings
        call MAPL_GridCompGetResource( GC, 'USE_INV_LAYERS', USE_INV_LAYERS                    ,default= .TRUE., _RC )
        call MAPL_GridCompGetResource( GC, 'USE_SMOOTH_TEND', USE_SMOOTH_TEND                  ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'SATUR_CALC', SATUR_CALC                            ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'BC_METH', BC_METH                                  ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'USE_REBCB', USE_REBCB                              ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'LAMBAU_DEEP', LAMBAU_DEEP                          ,default= 2.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'LAMBAU_MID', LAMBAU_MID                            ,default= 2.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'MOIST_TRIGGER', MOIST_TRIGGER                      ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'FRAC_MODIS', FRAC_MODIS                            ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'EVAP_FIX', EVAP_FIX                                ,default= 0,     _RC )
        call MAPL_GridCompGetResource( GC, 'AVE_LAYER_DP', CUM_AVE_LAYER(DEEP)                 ,default= 50.,   _RC )
        call MAPL_GridCompGetResource( GC, 'AVE_LAYER_SH', CUM_AVE_LAYER(SHAL)                 ,default= 30.,   _RC )
        call MAPL_GridCompGetResource( GC, 'AVE_LAYER_MD', CUM_AVE_LAYER(MID)                  ,default= 50.,   _RC )
      endif

      call MAPL_GridCompGetResource( GC, 'OUTPUT_SOUND', OUTPUT_SOUND                        ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'USE_MEMORY', USE_MEMORY                            ,default=-1,    _RC )
      call MAPL_GridCompGetResource( GC, 'TAU_OCEA_CP', TAU_OCEA_CP                          ,default= 21600., _RC )
      call MAPL_GridCompGetResource( GC, 'TAU_LAND_CP', TAU_LAND_CP                          ,default= 21600., _RC )
      call MAPL_GridCompGetResource( GC, 'DOWNDRAFT', DOWNDRAFT                              ,default= .TRUE., _RC )
      call MAPL_GridCompGetResource( GC, 'USE_FLUX_FORM', USE_FLUX_FORM                      ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'USE_TRACER_EVAP', USE_TRACER_EVAP                  ,default= 1,    _RC )
      call MAPL_GridCompGetResource( GC, 'APPLY_SUB_MP', APPLY_SUB_MP                        ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'ALP1', ALP1                                        ,default= 1.,   _RC )
      call MAPL_GridCompGetResource( GC, 'LIGHTNING_DIAG', LIGHTNING_DIAG                    ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'OVERSHOOT', OVERSHOOT                              ,default= 0.,   _RC )
      call MAPL_GridCompGetResource( GC, 'USE_WETBULB', USE_WETBULB                          ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'LAMBAU_SHDN', LAMBAU_SHDN                          ,default= 2.,   _RC )
      call MAPL_GridCompGetResource( GC, 'MAX_TQ_TEND', MAX_TQ_TEND                          ,default= 100., _RC )
      call MAPL_GridCompGetResource( GC, 'BETA_SH', BETA_SH                                  ,default= 2.2,  _RC )
      call MAPL_GridCompGetResource( GC, 'USE_LINEAR_SUBCL_MF', USE_LINEAR_SUBCL_MF          ,default= 0,    _RC )
      call MAPL_GridCompGetResource( GC, 'CAP_MAXS_DP', CUM_CAP_MAXS(DEEP)                   ,default= 150.0, _RC )
      call MAPL_GridCompGetResource( GC, 'CAP_MAXS_SH', CUM_CAP_MAXS(SHAL)                   ,default=  50.0, _RC )
      call MAPL_GridCompGetResource( GC, 'CAP_MAXS_MD', CUM_CAP_MAXS(MID)                    ,default= 150.0, _RC )
      call MAPL_GridCompGetResource( GC, 'GF_ENV_SETTING', env_setting                    ,default= 'DYNAMICS', _RC)
      GF_ENV_SETTING = env_setting
      call MAPL_GridCompGetResource( GC, 'STOCH_TOP', STOCH_TOP                              ,default= 2.50,  _RC)
      call MAPL_GridCompGetResource( GC, 'STOCH_BOT', STOCH_BOT                              ,default= 0.75,  _RC)
      call MAPL_GridCompGetResource( GC, 'STOCHASTIC_CNV', STOCHASTIC_CNV                    ,default= .FALSE.,_RC)
      call MAPL_GridCompGetResource( GC, 'USE_FCT', USE_FCT                                  ,default= 0,     _RC )
      call MAPL_GridCompGetResource( GC, 'ADV_TRIGGER', ADV_TRIGGER                          ,default= 0,     _RC )
      if (INT(ZERO_DIFF_HEI_EDT) == 0) then
        ! ============================================================================
        ! Grell-Freitas Convection Scheme Tuning Parameters
        ! ============================================================================
        ! These parameters control convective updraft/downdraft structure and have
        ! been optimized based on observational constraints and physical reasoning:
        !
        ! HEI_DOWN (Downdraft Initiation Height, fraction from cloud base):
        !   - Ocean: Lower values (0.35-0.45) downdrafts start higher in cloud
        !     where moist environment allows evaporative cooling aloft
        !   - Land: Higher values (0.55-0.65) downdrafts start lower in cloud
        !     where dry boundary layer concentrates evaporation near cloud base
        !   - Mid-level: Higher than deep convection (drier mid-troposphere)
        !
        ! HEI_UPDF (Updraft Detrainment Height, fraction from cloud base):
        !   - Ocean: Higher values (0.50-0.55) deeper penetration in weak-shear,
        !     moist environments produces more anvil cloud
        !   - Land: Lower values (0.40-0.45) earlier detrainment in sheared,
        !     dry environments produces less stratiform precipitation
        !   - Mid-level: Lower than deep convection (shallower cloud depth)
        !
        ! EDT (Evaporative Downdraft Efficiency, 1 - precipitation efficiency):
        !   - Ocean: Narrower range (0.10-0.50/0.60) limited evaporation in
        !     moist environment yields 50-90% precipitation efficiency
        !   - Land: Wider range (0.20-0.70/0.80) stronger evaporation in dry
        !     environment yields 30-80% precipitation efficiency
        !   - Mid-level: Higher values than deep (drier mid-levels, more evaporation)
        !
        ! ============================================================================
        call MAPL_GridCompGetResource( GC, 'USE_EXCESS_DP', CUM_USE_EXCESS(DEEP)               ,default= 3,     _RC )
        call MAPL_GridCompGetResource( GC, 'USE_EXCESS_SH', CUM_USE_EXCESS(SHAL)               ,default= 3,     _RC )
        call MAPL_GridCompGetResource( GC, 'USE_EXCESS_MD', CUM_USE_EXCESS(MID)                ,default= 3,     _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_LAND_DP', CUM_HEI_DOWN_LAND(DEEP)         ,default= 0.30,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_LAND_SH', CUM_HEI_DOWN_LAND(SHAL)         ,default= 0.00,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_LAND_MD', CUM_HEI_DOWN_LAND(MID)          ,default= 0.30,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_OCEAN_DP', CUM_HEI_DOWN_OCEAN(DEEP)       ,default= 0.60,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_OCEAN_SH', CUM_HEI_DOWN_OCEAN(SHAL)       ,default= 0.00,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_OCEAN_MD', CUM_HEI_DOWN_OCEAN(MID)        ,default= 0.60,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_LAND_DP', CUM_HEI_UPDF_LAND(DEEP)         ,default= 0.40,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_LAND_SH', CUM_HEI_UPDF_LAND(SHAL)         ,default= 0.20,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_LAND_MD', CUM_HEI_UPDF_LAND(MID)          ,default= 0.40,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_OCEAN_DP', CUM_HEI_UPDF_OCEAN(DEEP)       ,default= 0.40,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_OCEAN_SH', CUM_HEI_UPDF_OCEAN(SHAL)       ,default= 0.20,  _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_OCEAN_MD', CUM_HEI_UPDF_OCEAN(MID)        ,default= 0.40,  _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_LAND_DP', CUM_MIN_EDT_LAND(DEEP)           ,default= 0.10,  _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_LAND_SH', CUM_MIN_EDT_LAND(SHAL)           ,default= 0.00,  _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_LAND_MD', CUM_MIN_EDT_LAND(MID)            ,default= 0.10,  _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_OCEAN_DP', CUM_MIN_EDT_OCEAN(DEEP)         ,default= 0.10,  _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_OCEAN_SH', CUM_MIN_EDT_OCEAN(SHAL)         ,default= 0.00,  _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_OCEAN_MD', CUM_MIN_EDT_OCEAN(MID)          ,default= 0.10,  _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_LAND_DP', CUM_MAX_EDT_LAND(DEEP)           ,default= 0.90,  _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_LAND_SH', CUM_MAX_EDT_LAND(SHAL)           ,default= 0.00,  _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_LAND_MD', CUM_MAX_EDT_LAND(MID)            ,default= 0.90,  _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_OCEAN_DP', CUM_MAX_EDT_OCEAN(DEEP)         ,default= 0.90,  _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_OCEAN_SH', CUM_MAX_EDT_OCEAN(SHAL)         ,default= 0.00,  _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_OCEAN_MD', CUM_MAX_EDT_OCEAN(MID)          ,default= 0.90,  _RC )
      else
        ! support for legacy GF settings
        call MAPL_GridCompGetResource( GC, 'USE_EXCESS_DP', CUM_USE_EXCESS(DEEP)               ,default= 2,     _RC )
        call MAPL_GridCompGetResource( GC, 'USE_EXCESS_SH', CUM_USE_EXCESS(SHAL)               ,default= 3,     _RC )
        call MAPL_GridCompGetResource( GC, 'USE_EXCESS_MD', CUM_USE_EXCESS(MID)                ,default= 2,     _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_LAND_DP', CUM_HEI_DOWN_LAND(DEEP)         ,default= 0.3,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_LAND_SH', CUM_HEI_DOWN_LAND(SHAL)         ,default= 0.0,   _RC )    
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_LAND_MD', CUM_HEI_DOWN_LAND(MID)          ,default= 0.3,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_OCEAN_DP', CUM_HEI_DOWN_OCEAN(DEEP)       ,default= 0.4,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_OCEAN_SH', CUM_HEI_DOWN_OCEAN(SHAL)       ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_DOWN_OCEAN_MD', CUM_HEI_DOWN_OCEAN(MID)        ,default= 0.4,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_LAND_DP', CUM_HEI_UPDF_LAND(DEEP)         ,default= 0.4,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_LAND_SH', CUM_HEI_UPDF_LAND(SHAL)         ,default= 0.2,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_LAND_MD', CUM_HEI_UPDF_LAND(MID)          ,default= 0.4,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_OCEAN_DP', CUM_HEI_UPDF_OCEAN(DEEP)       ,default= 0.4,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_OCEAN_SH', CUM_HEI_UPDF_OCEAN(SHAL)       ,default= 0.2,   _RC )
        call MAPL_GridCompGetResource( GC, 'HEI_UPDF_OCEAN_MD', CUM_HEI_UPDF_OCEAN(MID)        ,default= 0.6,   _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_LAND_DP', CUM_MIN_EDT_LAND(DEEP)           ,default= 0.1,   _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_LAND_SH', CUM_MIN_EDT_LAND(SHAL)           ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_LAND_MD', CUM_MIN_EDT_LAND(MID)            ,default= 0.1,   _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_OCEAN_DP', CUM_MIN_EDT_OCEAN(DEEP)         ,default= 0.1,   _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_OCEAN_SH', CUM_MIN_EDT_OCEAN(SHAL)         ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'MIN_EDT_OCEAN_MD', CUM_MIN_EDT_OCEAN(MID)          ,default= 0.1,   _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_LAND_DP', CUM_MAX_EDT_LAND(DEEP)           ,default= 0.9,   _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_LAND_SH', CUM_MAX_EDT_LAND(SHAL)           ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_LAND_MD', CUM_MAX_EDT_LAND(MID)            ,default= 0.9,   _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_OCEAN_DP', CUM_MAX_EDT_OCEAN(DEEP)         ,default= 0.9,   _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_OCEAN_SH', CUM_MAX_EDT_OCEAN(SHAL)         ,default= 0.0,   _RC )
        call MAPL_GridCompGetResource( GC, 'MAX_EDT_OCEAN_MD', CUM_MAX_EDT_OCEAN(MID)          ,default= 0.9,   _RC )
      endif
      call MAPL_GridCompGetResource( GC, 'SCLM_DEEP', SCLM_DEEP                              ,default= 1.0,   _RC)
      call MAPL_GridCompGetResource( GC, 'FIX_CNV_CLOUD', FIX_CNV_CLOUD                      ,default= .FALSE., _RC)
    ELSE
      ! Inititialize parameters of convection scheme GF (circa 2019)
      call MAPL_GridCompGetResource( GC, 'DEEP', icumulus_gf(deep)                ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'SHALLOW', icumulus_gf(shal)             ,default= 0 , _RC )
      call MAPL_GridCompGetResource( GC, 'CONGESTUS', icumulus_gf(mid)            ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'CLOSURE_DEEP', closure_choice(deep)     ,default= 0 , _RC )
      call MAPL_GridCompGetResource( GC, 'CLOSURE_SHALLOW', closure_choice(shal)  ,default= 7 , _RC )
      call MAPL_GridCompGetResource( GC, 'CLOSURE_CONGESTUS', closure_choice(mid) ,default= 3 , _RC )
      call MAPL_GridCompGetResource( GC, 'ENTR_DP', CUM_ENTR_RATE(DEEP)           ,default= 1.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'ENTR_SH', CUM_ENTR_RATE(SHAL)           ,default= 1.4e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'ENTR_MD', CUM_ENTR_RATE(MID)            ,default= 9.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'C1', C1                                 ,default= 1.0e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'USE_TRACER_TRANSP', USE_TRACER_TRANSP   ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'USE_TRACER_SCAVEN', USE_TRACER_SCAVEN   ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'USE_SCALE_DEP', USE_SCALE_DEP           ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'DICYCLE', DICYCLE                       ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'STOCHASTIC_CNV', STOCHASTIC_CNV         ,default= .FALSE., _RC)
      call MAPL_GridCompGetResource( GC, 'GF_MIN_AREA', GF_MIN_AREA               ,default= 1.e6,  _RC )
      call MAPL_GridCompGetResource( GC, 'TAU_DEEP', TAU_DEEP                     ,default= 5400.0, _RC )
      call MAPL_GridCompGetResource( GC, 'TAU_MID', TAU_MID                       ,default= 3600.0, _RC )
      call MAPL_GridCompGetResource( GC, 'USE_FLUX_FORM', USE_FLUX_FORM           ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'USE_FCT', USE_FCT                       ,default= 0 , _RC )
      call MAPL_GridCompGetResource( GC, 'USE_TRACER_EVAP', USE_TRACER_EVAP       ,default= 1 , _RC )
      call MAPL_GridCompGetResource( GC, 'ALP1', ALP1                             ,default= 1., _RC )
      call MAPL_GridCompGetResource( GC, 'BC_METH', BC_METH                       ,default= 0 , _RC )
      call MAPL_GridCompGetResource( GC, 'C0_DEEP', C0_DEEP                       ,default= 2.0e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'C0_MID', C0_MID                         ,default= 2.0e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'C0_SHAL', C0_SHAL                       ,default= 1.0e-3,_RC )
      call MAPL_GridCompGetResource( GC, 'QRC_CRIT', QRC_CRIT                     ,default= 2.0e-4,_RC )
      call MAPL_GridCompGetResource( GC, 'GF_ENV_SETTING', env_setting         ,default= 'DYNAMICS', _RC)
      GF_ENV_SETTING = env_setting
      call MAPL_GridCompGetResource( GC, 'SCLM_DEEP', SCLM_DEEP                   ,default= 1.0    , _RC)
      call MAPL_GridCompGetResource( GC, 'FIX_CNV_CLOUD', FIX_CNV_CLOUD           ,default= .FALSE., _RC)
    ENDIF

    endif ! USE_PYMOIST_GF2020

end subroutine GF_Initialize


subroutine GF_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Local derived type aliases

    type (ESMF_Grid    )            :: GRID
    type (ESMF_State   )            :: INTERNAL
    type (ESMF_TimeInterval)        :: TINT
    real(ESMF_KIND_R8)              :: DT_R8
    real                            :: GF_DT, MOIST_DT
    type(ESMF_Alarm)                :: alarm
    logical                         :: alarm_is_ringing

    ! Local variables
    integer                         :: I, J, L
    integer                         :: IM,JM,LM
    real, allocatable, dimension(:,:) :: LONS
    real, allocatable, dimension(:,:) :: LATS
    real                            :: minrhx, fqi_local, tmp_local
    logical                         :: ptr_is_assoc

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: NACTL
    real, pointer, dimension(:,:,:) :: CNV_TR ! tracer memory

    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W, KH, OMEGA
    real, pointer, contiguous, dimension(:,:,:) :: p3d ! scratch for the edge-field bounds remap
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

    ! Re-exported internals (were ALLOC'd exports shared with the driver and the other schemes)
    real, pointer, dimension(:,:,:) :: MFD_DC, UMF_DC, BYNCY, PFL_CN, PFI_CN
    real, pointer, dimension(:,:,:) :: DUDT_DC, DVDT_DC, DTDT_DC, DQVDT_DC, DQIDT_DC, DQLDT_DC, DQADT_DC
    real, pointer, dimension(:,:  ) :: CNV_FRC, SRF_TYPE, CNPCPRATE
    ! Imports
    real, pointer, dimension(:,:  ) :: T2M,Q2M,TA,QA,SH,EVAP,PHIS
    ! Exports computed in local arrays (were ALLOC'd), copied to the export if it is connected
    real, allocatable, dimension(:,:,:) :: CNV_CVW, CNV_QC, CNV_DQCDT, CNV_PRC3, CNV_UPDF
    real, allocatable, dimension(:,:,:) :: CNV_MF0, ENTLAM
    real, allocatable, dimension(:,:,:) :: MUPDP,MDNDP,MUPSH,MUPMD,WQT_DC
    real, allocatable, dimension(:,:  ) :: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
    real, allocatable, dimension(:,:  ) :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_DP,TAU_MD
    real, allocatable, dimension(:,:  ) :: TPWI,TPWI_star,LFR_GF
    real, allocatable, dimension(:,:,:) :: RSU_CN,REV_CN
    real, allocatable, dimension(:,:  ) :: SIGMA_DEEP, SIGMA_MID
    real, allocatable, dimension(:,:,:) :: ENTR_DP, ENTR_MD, ENTR_SH
    real, allocatable, dimension(:,:,:) :: SGS_VVEL_DP, SGS_VVEL_MD, SGS_VVEL_SH
    real, allocatable, dimension(:,:  ) :: CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH
    real, allocatable, dimension(:,:,:) :: DNDCNV, DNICNV
    ! Their export pointers
    real, pointer, dimension(:,:,:) :: CNV_CVW_ptr, CNV_QC_ptr, CNV_DQCDT_ptr, CNV_PRC3_ptr, CNV_UPDF_ptr
    real, pointer, dimension(:,:,:) :: CNV_MF0_ptr, ENTLAM_ptr
    real, pointer, dimension(:,:,:) :: MUPDP_ptr,MDNDP_ptr,MUPSH_ptr,MUPMD_ptr,WQT_DC_ptr
    real, pointer, dimension(:,:  ) :: MFDP_ptr,MFSH_ptr,MFMD_ptr,ERRDP_ptr,ERRSH_ptr,ERRMD_ptr
    real, pointer, dimension(:,:  ) :: AA0_ptr,AA1_ptr,AA2_ptr,AA3_ptr,AA1_BL_ptr,AA1_CIN_ptr,TAU_BL_ptr,TAU_DP_ptr,TAU_MD_ptr
    real, pointer, dimension(:,:  ) :: TPWI_ptr,TPWI_star_ptr,LFR_GF_ptr
    real, pointer, dimension(:,:,:) :: RSU_CN_ptr,REV_CN_ptr
    real, pointer, dimension(:,:  ) :: SIGMA_DEEP_ptr, SIGMA_MID_ptr
    real, pointer, dimension(:,:,:) :: ENTR_DP_ptr, ENTR_MD_ptr, ENTR_SH_ptr
    real, pointer, dimension(:,:,:) :: SGS_VVEL_DP_ptr, SGS_VVEL_MD_ptr, SGS_VVEL_SH_ptr
    real, pointer, dimension(:,:  ) :: CNV_TOPP_DP_ptr, CNV_TOPP_MD_ptr, CNV_TOPP_SH_ptr
    real, pointer, dimension(:,:,:) :: DNDCNV_ptr, DNICNV_ptr
    ! Exports
    real, pointer, dimension(:,:,:) :: ENTR
    real, pointer, dimension(:,:,:) ::   DQVDT_FILL
    real, pointer, dimension(:,:,:) :: DQLLSDT_FILL
    real, pointer, dimension(:,:,:) :: DQLCNDT_FILL
    real, pointer, dimension(:,:,:) :: DQILSDT_FILL
    real, pointer, dimension(:,:,:) :: DQICNDT_FILL

    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D
    
    character(len=256) :: DIAG_MSG

    type( ESMF_VM )                 :: VMG

    ! DSL fields
    real, pointer, dimension(:,:) :: DSL__GF2020_LONS, DSL__GF2020_LATS

    call MAPL_GridCompGet(GC, num_levels=LM, grid=GRID, _RC)
    call MAPL_GridGet(GRID, IM=IM, JM=JM, _RC)
    call MAPL_GridGetCoordinates(GRID, longitudes=LONS, latitudes=LATS, _RC)
    call MAPL_GridCompGetInternalState(GC, INTERNAL, _RC)
    ! Moist's timestep (MAPL2: the RUNALARM ring interval)
    call MAPL_ClockGet(CLOCK, dt=MOIST_DT, _RC)

    call ESMF_GridCompGet ( GC, VM=VMG, _RC )

    ! Internals
    call MAPL_StateGetPointer(INTERNAL, Q,      'Q'       , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLLS,   'QLLS'    , _RC)  
    call MAPL_StateGetPointer(INTERNAL, QLCN,   'QLCN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, CLCN,   'CLCN'    , _RC)  
    call MAPL_StateGetPointer(INTERNAL, CLLS,   'CLLS'    , _RC)  
    call MAPL_StateGetPointer(INTERNAL, QILS,   'QILS'    , _RC)  
    call MAPL_StateGetPointer(INTERNAL, QICN,   'QICN'    , _RC)  
    call MAPL_StateGetPointer(INTERNAL, NACTL,  'NACTL'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNV_TR, 'CNV_TR'  , _RC)
    ! Imports
    call MAPL_StateGetPointer(IMPORT, T,'T'         , _RC)
    call MAPL_StateGetPointer(IMPORT, U,'U'         , _RC)
    call MAPL_StateGetPointer(IMPORT, V,'V'         , _RC)
    call MAPL_StateGetPointer(IMPORT, PLE,'PLE'       , _RC)
    ! Edge field: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
    p3d => PLE; PLE(1:IM,1:JM,0:LM) => p3d
    ! Initialize tendencies (internals: they keep their values between GF_RunAlarm rings)
    call MAPL_StateGetPointer(INTERNAL,  DUDT_DC,    'DUDT_DC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL,  DVDT_DC,    'DVDT_DC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL,  DTDT_DC,    'DTDT_DC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQVDT_DC,   'DQVDT_DC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQLDT_DC,   'DQLDT_DC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQIDT_DC,   'DQIDT_DC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQADT_DC,   'DQADT_DC'  , _RC)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GridCompTimerStart(GC, "--GF", _RC)

    ! Get parameters from generic state.
    !-----------------------------------

    if (USE_PYMOIST_GF2020) then
      call MAPL_StateGetPointer(INTERNAL, DSL__GF2020_LONS, 'DSL__GF2020_LONS', _RC)
      call MAPL_StateGetPointer(INTERNAL, DSL__GF2020_LATS, 'DSL__GF2020_LATS', _RC)
      DSL__GF2020_LONS = LONS
      DSL__GF2020_LATS = LATS
      call CNV_Tracers_To_SOA()
      call MAPL_pybridge_gcrun_with_internal( "pyMoist.fortran.param_interfaces.convection.GF2020_interface", MAPL, IMPORT, EXPORT, INTERNAL )
      call CNV_Tracers_To_AOS()
    else

    ! Imports
    call MAPL_StateGetPointer(IMPORT, FRLAND,'FRLAND'    , _RC)
    call MAPL_StateGetPointer(IMPORT, AREA,'AREA'      , _RC)
    call MAPL_StateGetPointer(IMPORT, ZLE,'ZLE'       , _RC)
    call MAPL_StateGetPointer(IMPORT, T,'T'         , _RC)
    call MAPL_StateGetPointer(IMPORT, U,'U'         , _RC)
    call MAPL_StateGetPointer(IMPORT, V,'V'         , _RC)
    call MAPL_StateGetPointer(IMPORT, W,'W'         , _RC)
    call MAPL_StateGetPointer(IMPORT, OMEGA,'OMEGA'     , _RC)
    call MAPL_StateGetPointer(IMPORT, KH,'KH'        , _RC)
    call MAPL_StateGetPointer(IMPORT, SH,'SH'        , _RC)
    call MAPL_StateGetPointer(IMPORT, EVAP,'EVAP'      , _RC)
    call MAPL_StateGetPointer(IMPORT, T2M,'T2M'       , _RC)
    call MAPL_StateGetPointer(IMPORT, Q2M,'Q2M'       , _RC)
    call MAPL_StateGetPointer(IMPORT, TA,'TA'        , _RC)
    call MAPL_StateGetPointer(IMPORT, QA,'QA'        , _RC)
    call MAPL_StateGetPointer(IMPORT, PHIS,'PHIS'      , _RC)
    call MAPL_StateGetPointer(IMPORT, DQVDTDYN,'DQVDTDYN'  , _RC)
    call MAPL_StateGetPointer(IMPORT, DTDTDYN,'DTDTDYN'   , _RC)
    call MAPL_StateGetPointer(IMPORT, QV_DYN_IN,'QV_DYN_IN' , _RC)
    call MAPL_StateGetPointer(IMPORT, U_DYN_IN,'U_DYN_IN'  , _RC)
    call MAPL_StateGetPointer(IMPORT, V_DYN_IN,'V_DYN_IN'  , _RC)
    call MAPL_StateGetPointer(IMPORT, T_DYN_IN,'T_DYN_IN'  , _RC)
    call MAPL_StateGetPointer(IMPORT, PLE_DYN_IN,'PLE_DYN_IN', _RC)
    ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
    p3d => ZLE;        ZLE(1:IM,1:JM,0:LM)        => p3d
    p3d => KH;         KH(1:IM,1:JM,0:LM)         => p3d
    p3d => PLE_DYN_IN; PLE_DYN_IN(1:IM,1:JM,0:LM) => p3d
    call MAPL_StateGetPointer(IMPORT, DTDT_BL,'DTDT_BL'   , _RC)
    call MAPL_StateGetPointer(IMPORT, DQDT_BL,'DQDT_BL'   , _RC)
    call MAPL_StateGetPointer(IMPORT, RADSW,'RADSW'     , _RC)
    call MAPL_StateGetPointer(IMPORT, RADLW,'RADLW'     , _RC)
    call MAPL_StateGetPointer(IMPORT, KPBL,'KPBL_SC'   , _RC)

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

    ! Derived States
    !--------------------------------------------------------------
    
    ! 1. Top-of-atmosphere edge (L = 0)
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, ZLE0, ZLE) &
    !$OMP PRIVATE(I, J)
    do J = 1, JM
       !DIR$ IVDEP
       do I = 1, IM
          ZLE0(I,J,0) = ZLE(I,J,0) - ZLE(I,J,LM)
       end do
    end do
    !$OMP END PARALLEL DO

    ! 2. Remaining edges and all layer variables (L = 1 to LM)
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, ZLE0, ZLE, PL, PLE, PK, &
    !$OMP        ZL0, TH, T, MASS) &
    !$OMP PRIVATE(I, J, L)
    do L = 1, LM
       do J = 1, JM
          !DIR$ IVDEP
          do I = 1, IM
             
             ! Edge variables
             ZLE0(I,J,L) = ZLE(I,J,L) - ZLE(I,J,LM)
             
             ! Layer variables
             PL(I,J,L)   = 0.5 * (PLE(I,J,L-1) + PLE(I,J,L))
             PK(I,J,L)   = (PL(I,J,L) / MAPL_P00)**(MAPL_KAPPA)
             
             ZL0(I,J,L)  = 0.5 * (ZLE0(I,J,L-1) + ZLE0(I,J,L))
             
             TH(I,J,L)   = T(I,J,L) / PK(I,J,L)
             
             MASS(I,J,L) = (PLE(I,J,L) - PLE(I,J,L-1)) / MAPL_GRAV
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call ESMF_ClockGetAlarm(clock, 'GF_RunAlarm', alarm, _RC)
    alarm_is_ringing = ESMF_AlarmIsRinging(alarm, _RC)

    if (alarm_is_ringing) then

!!! call WRITE_PARALLEL('GF is Running')
    call ESMF_AlarmRingerOff(alarm, _RC)
    call ESMF_AlarmGet(alarm, RingInterval=TINT, _RC)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,_RC)
    GF_DT = DT_R8

    ! Required Exports (connectivities to moist siblings): re-exported internals
    call MAPL_StateGetPointer(INTERNAL, MFD_DC,     'MFD_DC'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, UMF_DC,     'UMF_DC'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNV_FRC,    'CNV_FRC'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, SRF_TYPE,   'SRF_TYPE'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, BYNCY,      'BYNCY'     , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNPCPRATE,  'CNPCPRATE' , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFI_CN   ,'PFI_CN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFL_CN   ,'PFL_CN'    , _RC)
    ! Edge fields: MAPL3 pointers are 1:LM+1; MAPL2's were 0:LM, which this code assumes
    p3d => UMF_DC; UMF_DC(1:IM,1:JM,0:LM) => p3d
    p3d => PFI_CN; PFI_CN(1:IM,1:JM,0:LM) => p3d
    p3d => PFL_CN; PFL_CN(1:IM,1:JM,0:LM) => p3d
    ! Exports computed locally, zero-filled as MAPL2's ALLOC did (some are not written on every path)
    allocate(CNV_CVW(IM,JM,LM), CNV_QC(IM,JM,LM), CNV_DQCDT(IM,JM,LM), CNV_PRC3(IM,JM,LM), CNV_UPDF(IM,JM,LM), source=0.0)
    allocate(CNV_MF0(IM,JM,LM), ENTLAM(IM,JM,LM), source=0.0)
    allocate(SIGMA_DEEP(IM,JM), SIGMA_MID(IM,JM), TPWI(IM,JM), TPWI_star(IM,JM), LFR_GF(IM,JM), source=0.0)
    TPWI       = SUM( Q*MASS, 3 )
    TPWI_star  = SUM( GEOS_QSAT(T,PL,PASCALS=.true.)*MASS, 3 )
    allocate(MUPDP(IM,JM,LM), MDNDP(IM,JM,LM), MUPSH(IM,JM,LM), MUPMD(IM,JM,LM), source=0.0)
    allocate(MFDP(IM,JM), MFSH(IM,JM), MFMD(IM,JM), ERRDP(IM,JM), ERRSH(IM,JM), ERRMD(IM,JM), source=0.0)
    allocate(RSU_CN(IM,JM,LM), REV_CN(IM,JM,LM), source=0.0)
    allocate(AA0(IM,JM), AA1(IM,JM), AA2(IM,JM), AA3(IM,JM), AA1_BL(IM,JM), AA1_CIN(IM,JM), source=0.0)
    allocate(TAU_BL(IM,JM), TAU_DP(IM,JM), TAU_MD(IM,JM), source=0.0)
    allocate(WQT_DC(IM,JM,0:LM), source=0.0) ! edge
    allocate(ENTR_DP(IM,JM,LM), ENTR_MD(IM,JM,LM), ENTR_SH(IM,JM,LM), source=0.0)
    allocate(SGS_VVEL_DP(IM,JM,LM), SGS_VVEL_MD(IM,JM,LM), SGS_VVEL_SH(IM,JM,LM), source=0.0)
    allocate(CNV_TOPP_DP(IM,JM), CNV_TOPP_MD(IM,JM), CNV_TOPP_SH(IM,JM), source=0.0)
    allocate(DNDCNV(IM,JM,LM), DNICNV(IM,JM,LM), source=0.0)



    if (STOCHASTIC_CNV) then
       ! Create bit-processor-reproducible random white noise for convection [0:1]
       SEEDINI = 1000000 * ( 100*T(:,:,LM)   - INT( 100*T(:,:,LM) ) )
       SEEDCNV = SQRT(MAX(MIN(SEEDINI/1000000.0,1.0),0.0))
       ! Create stochastic variability to GF sigma
       SEEDCNV = SQRT(1.0-(1.0-SEEDCNV))*(STOCH_TOP-STOCH_BOT)+STOCH_BOT
    else
       SEEDCNV = 1.0
    endif
    call MAPL_StateGetPointer(EXPORT, PTR2D,  'STOCH_CNV', _RC)
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

    DNDCNV   = 0.0
    DNICNV   = 0.0
    DQLDT_DC = 0.0
    DQIDT_DC = 0.0
        
    IF (USE_GF2020==1) THEN
         ! Convert OMEGA (Pa/s) to W (m/s)
         TMP3D = -1*OMEGA/(MAPL_GRAV*PL/(MAPL_RDRY*T*(1.0+MAPL_VIREPS*Q)))
         ! Compute subgrid vertical velocities
         call compute_sgs_vvel(IM,JM,LM,ZLE0,TMP3D,BYNCY, &
                               SGS_VVEL_DP,SGS_VVEL_MD,SGS_VVEL_SH)
         !- call GF2020 interface routine: PLE and PL are passed in Pa
         call GF2020_Interface( &
              !--- Grid and Timing
              IM, JM, LM, LONS, LATS, GF_DT, &
              !--- 3D Host Environmental State
              PLE, PL, ZLE0, ZL0, PK, MASS, KH, &
              T, TH, Q, U, V, TMP3D, OMEGA, BYNCY, (QLCN+QLLS), (QICN+QILS), (CLCN+CLLS), NACTL, &
              !--- 3D Forcings (Dynamics, Advection, Radiation, PBL)
              QV_DYN_IN, PLE_DYN_IN, U_DYN_IN, V_DYN_IN, T_DYN_IN, &
              RADSW, RADLW, DQDT_BL, DTDT_BL, DTDTDYN, DQVDTDYN, &
              !--- 2D Surface and Environmental Inputs
              FRLAND, TMP2D, T2M, Q2M, TA, QA, SH, EVAP, PHIS, &
              KPBL, CNV_FRC, SRF_TYPE, SEEDCNV, TPWI, TPWI_star, &
              !--- 3D Tracers and Sub-grid Velocity (InOut)
              CNV_TR, SGS_VVEL_DP, SGS_VVEL_MD, SGS_VVEL_SH, &
              !--- Output Tendencies
              DQVDT_DC, DTDT_DC, DUDT_DC, DVDT_DC, &
              !--- 2D Output Diagnostics
              CNPCPRATE, LFR_GF, SIGMA_DEEP, SIGMA_MID, &
              CNV_TOPP_DP, CNV_TOPP_MD, CNV_TOPP_SH, &
              MFDP, MFSH, MFMD, ERRDP, ERRSH, ERRMD, &
              AA0, AA1, AA2, AA3, AA1_BL, AA1_CIN, TAU_BL, TAU_DP, TAU_MD, &
              !--- 3D Output Diagnostics and Profiles
              CNV_MF0, CNV_PRC3, MFD_DC, CNV_DQCDT, ENTLAM, &
              UMF_DC, CNV_UPDF, CNV_CVW, CNV_QC, WQT_DC, &
              REVSU, PRFIL, ENTR_DP, ENTR_MD, ENTR_SH, &
              MUPDP, MUPSH, MUPMD, MDNDP, &
              AeroPropsNew, DNDCNV, DNICNV, DQLDT_DC, DQIDT_DC)
    ELSE
         !- call GF/GEOS5 interface routine
         ! PLE and PL are passed in Pa
         call GF_GEOS5_Interface( IM,JM,LM,LONS,LATS,GF_DT                       &
                                 ,PLE, PL, ZLE0, ZL0, PK, MASS, OMEGA               &
                                 ,T, TH, Q, U, V, QLCN, QICN, QLLS, QILS, CNPCPRATE &
                                 ,CNV_MF0, CNV_PRC3, MFD_DC, CNV_DQCDT,ENTLAM       &
                                 ,UMF_DC, CNV_UPDF, CNV_CVW, CNV_QC , CLCN          &
                                 ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                                 ,RADSW   ,RADLW  ,DQDT_BL  ,DTDT_BL                &
                                 ,FRLAND, TMP2D, T2M, Q2M      &
                                 ,TA ,QA ,SH ,EVAP ,PHIS                            &
                                 ,KPBL ,CNV_FRC, SRF_TYPE                           &
                                 ,SEEDCNV, SIGMA_DEEP, SIGMA_MID                    &
                                 ,DQVDT_DC,DTDT_DC,DUDT_DC,DVDT_DC                  &
                                 ,MUPDP,MUPSH,MUPMD                                 &
                                 ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
                                 ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_DP,TAU_MD &
                                 ,DTDTDYN,DQVDTDYN                                  &
                                 ,REVSU, PRFIL)
    ENDIF
        
    
    ! update DeepCu QL/QI/CF tendencies
    
    IF (USE_GF2020 == 1 .AND. USE_CUP_2M_MOISTURE) THEN
    
       call WRITE_PARALLEL ('Using CUP_2M_MOISTURE') 
      
       TMP3D = DQLDT_DC + DQIDT_DC
      
       WHERE (TMP3D > 1.0e-30)
          fQi = DQIDT_DC / TMP3D
       ELSEWHERE
          fQi = ice_fraction(T + DTDT_DC*GF_DT, CNV_FRC, SRF_TYPE)
       END WHERE

       CNV_DQCDT = TMP3D * MASS   ! keep diagnostic total consistent
       
          ! Export
      call MAPL_StateGetPointer(EXPORT, PTR3D, 'CNV_FICE', _RC)
      if (associated(PTR3D)) PTR3D = fQi
      
       DQADT_DC = MFD_DC*SCLM_DEEP/MASS
    ! evap/subl and precip fluxes
      do L=1,LM
         !--- sublimation/evaporation tendencies (kg/kg/s)
           RSU_CN (:,:,L) = REVSU(:,:,L)*     fQi(:,:,L)
           REV_CN (:,:,L) = REVSU(:,:,L)*(1.0-fQi(:,:,L))
         !--- preciptation fluxes (kg/kg/s)
           PFI_CN (:,:,L) = PRFIL(:,:,L)*     fQi(:,:,L)
           PFL_CN (:,:,L) = PRFIL(:,:,L)*(1.0-fQi(:,:,L))
      enddo

        
    
    ELSE !USE GF2M
    

        call MAPL_StateGetPointer(EXPORT, PTR3D, 'CNV_FICE', _RC)
        ptr_is_assoc = associated(PTR3D)

        ! Update DeepCu QL/QI/CF tendencies, evap/subl and precip fluxes
        !--------------------------------------------------------------
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(IM, JM, LM, T, DTDT_DC, GF_DT, CNV_FRC, SRF_TYPE, &
        !$OMP        CNV_DQCDT, MASS, DQLDT_DC, DQIDT_DC, DQADT_DC, &
        !$OMP        MFD_DC, SCLM_DEEP, RSU_CN, REVSU, REV_CN, &
        !$OMP        PFI_CN, PRFIL, PFL_CN, ptr_is_assoc, PTR3D) &
        !$OMP PRIVATE(I, J, L, fQi_local, tmp_local)
        do L = 1, LM
           do J = 1, JM
              !DIR$ IVDEP
              do I = 1, IM
                 ! 1. Calculate local ice fraction and tmp scalar
                 fQi_local = ice_fraction(T(I,J,L) + DTDT_DC(I,J,L) * GF_DT, CNV_FRC(I,J), SRF_TYPE(I,J))
                 tmp_local = CNV_DQCDT(I,J,L) / MASS(I,J,L)

                 ! Fill the exported 3D pointer if associated
                 if (ptr_is_assoc) PTR3D(I,J,L) = fQi_local

                 ! 2. Update DeepCu QL/QI/CF tendencies
                 DQLDT_DC(I,J,L) = (1.0 - fQi_local) * tmp_local
                 DQIDT_DC(I,J,L) = fQi_local * tmp_local
                 DQADT_DC(I,J,L) = MFD_DC(I,J,L) * SCLM_DEEP / MASS(I,J,L)

                 ! 3. Evap/subl and precip fluxes (kg/kg/s)
                 RSU_CN(I,J,L) = REVSU(I,J,L) * fQi_local
                 REV_CN(I,J,L) = REVSU(I,J,L) * (1.0 - fQi_local)

                 PFI_CN(I,J,L) = PRFIL(I,J,L) * fQi_local
                 PFL_CN(I,J,L) = PRFIL(I,J,L) * (1.0 - fQi_local)
              end do
           end do
        end do
        !$OMP END PARALLEL DO

    END IF  !USE GF2M
    
    ! Other Exports
    call MAPL_StateGetPointer(EXPORT, PTR3D, 'DQRC', _RC)
    if(associated(PTR3D)) PTR3D = CNV_PRC3 / GF_DT
    call MAPL_StateGetPointer(EXPORT, PTR2D, 'CCWP', _RC)
    if (associated(PTR2D)) PTR2D = SUM( CNV_QC*MASS , 3 )

    ! Copy the local arrays to their exports, if connected
    call MAPL_StateGetPointer(EXPORT, CNV_CVW_ptr,     'CNV_CVW'    , _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_QC_ptr,      'CNV_QC'     , _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_DQCDT_ptr,   'CNV_DQCDT'  , _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_PRC3_ptr,    'CNV_PRC3'   , _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_UPDF_ptr,    'CNV_UPDF'   , _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_MF0_ptr,     'CNV_MF0'    , _RC)
    call MAPL_StateGetPointer(EXPORT, ENTLAM_ptr,      'ENTLAM'     , _RC)
    call MAPL_StateGetPointer(EXPORT, SIGMA_DEEP_ptr,  'SIGMA_DEEP' , _RC)
    call MAPL_StateGetPointer(EXPORT, SIGMA_MID_ptr,   'SIGMA_MID'  , _RC)
    call MAPL_StateGetPointer(EXPORT, TPWI_ptr,        'TPWI'       , _RC)
    call MAPL_StateGetPointer(EXPORT, TPWI_star_ptr,   'TPWI_star'  , _RC)
    call MAPL_StateGetPointer(EXPORT, LFR_GF_ptr,      'LFR_GF'     , _RC)
    call MAPL_StateGetPointer(EXPORT, MUPDP_ptr,       'MUPDP'      , _RC)
    call MAPL_StateGetPointer(EXPORT, MDNDP_ptr,       'MDNDP'      , _RC)
    call MAPL_StateGetPointer(EXPORT, MUPSH_ptr,       'MUPSH'      , _RC)
    call MAPL_StateGetPointer(EXPORT, MUPMD_ptr,       'MUPMD'      , _RC)
    call MAPL_StateGetPointer(EXPORT, MFDP_ptr,        'MFDP'       , _RC)
    call MAPL_StateGetPointer(EXPORT, MFSH_ptr,        'MFSH'       , _RC)
    call MAPL_StateGetPointer(EXPORT, MFMD_ptr,        'MFMD'       , _RC)
    call MAPL_StateGetPointer(EXPORT, ERRDP_ptr,       'ERRDP'      , _RC)
    call MAPL_StateGetPointer(EXPORT, ERRSH_ptr,       'ERRSH'      , _RC)
    call MAPL_StateGetPointer(EXPORT, ERRMD_ptr,       'ERRMD'      , _RC)
    call MAPL_StateGetPointer(EXPORT, RSU_CN_ptr,      'RSU_CN'     , _RC)
    call MAPL_StateGetPointer(EXPORT, REV_CN_ptr,      'REV_CN'     , _RC)
    call MAPL_StateGetPointer(EXPORT, AA0_ptr,         'AA0'        , _RC)
    call MAPL_StateGetPointer(EXPORT, AA1_ptr,         'AA1'        , _RC)
    call MAPL_StateGetPointer(EXPORT, AA2_ptr,         'AA2'        , _RC)
    call MAPL_StateGetPointer(EXPORT, AA3_ptr,         'AA3'        , _RC)
    call MAPL_StateGetPointer(EXPORT, AA1_BL_ptr,      'AA1_BL'     , _RC)
    call MAPL_StateGetPointer(EXPORT, AA1_CIN_ptr,     'AA1_CIN'    , _RC)
    call MAPL_StateGetPointer(EXPORT, TAU_BL_ptr,      'TAU_BL'     , _RC)
    call MAPL_StateGetPointer(EXPORT, TAU_DP_ptr,      'TAU_DP'     , _RC)
    call MAPL_StateGetPointer(EXPORT, TAU_MD_ptr,      'TAU_MD'     , _RC)
    call MAPL_StateGetPointer(EXPORT, WQT_DC_ptr,      'WQT_DC'     , _RC)
    call MAPL_StateGetPointer(EXPORT, ENTR_DP_ptr,     'ENTR_DP'    , _RC)
    call MAPL_StateGetPointer(EXPORT, ENTR_MD_ptr,     'ENTR_MD'    , _RC)
    call MAPL_StateGetPointer(EXPORT, ENTR_SH_ptr,     'ENTR_SH'    , _RC)
    call MAPL_StateGetPointer(EXPORT, SGS_VVEL_DP_ptr, 'SGS_VVEL_DP', _RC)
    call MAPL_StateGetPointer(EXPORT, SGS_VVEL_MD_ptr, 'SGS_VVEL_MD', _RC)
    call MAPL_StateGetPointer(EXPORT, SGS_VVEL_SH_ptr, 'SGS_VVEL_SH', _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_TOPP_DP_ptr, 'CNV_TOPP_DP', _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_TOPP_MD_ptr, 'CNV_TOPP_MD', _RC)
    call MAPL_StateGetPointer(EXPORT, CNV_TOPP_SH_ptr, 'CNV_TOPP_SH', _RC)
    call MAPL_StateGetPointer(EXPORT, DNDCNV_ptr,      'DNDCNV'     , _RC)
    call MAPL_StateGetPointer(EXPORT, DNICNV_ptr,      'DNICNV'     , _RC)
    if (associated(CNV_CVW_ptr))     CNV_CVW_ptr     = CNV_CVW
    if (associated(CNV_QC_ptr))      CNV_QC_ptr      = CNV_QC
    if (associated(CNV_DQCDT_ptr))   CNV_DQCDT_ptr   = CNV_DQCDT
    if (associated(CNV_PRC3_ptr))    CNV_PRC3_ptr    = CNV_PRC3
    if (associated(CNV_UPDF_ptr))    CNV_UPDF_ptr    = CNV_UPDF
    if (associated(CNV_MF0_ptr))     CNV_MF0_ptr     = CNV_MF0
    if (associated(ENTLAM_ptr))      ENTLAM_ptr      = ENTLAM
    if (associated(SIGMA_DEEP_ptr))  SIGMA_DEEP_ptr  = SIGMA_DEEP
    if (associated(SIGMA_MID_ptr))   SIGMA_MID_ptr   = SIGMA_MID
    if (associated(TPWI_ptr))        TPWI_ptr        = TPWI
    if (associated(TPWI_star_ptr))   TPWI_star_ptr   = TPWI_star
    if (associated(LFR_GF_ptr))      LFR_GF_ptr      = LFR_GF
    if (associated(MUPDP_ptr))       MUPDP_ptr       = MUPDP
    if (associated(MDNDP_ptr))       MDNDP_ptr       = MDNDP
    if (associated(MUPSH_ptr))       MUPSH_ptr       = MUPSH
    if (associated(MUPMD_ptr))       MUPMD_ptr       = MUPMD
    if (associated(MFDP_ptr))        MFDP_ptr        = MFDP
    if (associated(MFSH_ptr))        MFSH_ptr        = MFSH
    if (associated(MFMD_ptr))        MFMD_ptr        = MFMD
    if (associated(ERRDP_ptr))       ERRDP_ptr       = ERRDP
    if (associated(ERRSH_ptr))       ERRSH_ptr       = ERRSH
    if (associated(ERRMD_ptr))       ERRMD_ptr       = ERRMD
    if (associated(RSU_CN_ptr))      RSU_CN_ptr      = RSU_CN
    if (associated(REV_CN_ptr))      REV_CN_ptr      = REV_CN
    if (associated(AA0_ptr))         AA0_ptr         = AA0
    if (associated(AA1_ptr))         AA1_ptr         = AA1
    if (associated(AA2_ptr))         AA2_ptr         = AA2
    if (associated(AA3_ptr))         AA3_ptr         = AA3
    if (associated(AA1_BL_ptr))      AA1_BL_ptr      = AA1_BL
    if (associated(AA1_CIN_ptr))     AA1_CIN_ptr     = AA1_CIN
    if (associated(TAU_BL_ptr))      TAU_BL_ptr      = TAU_BL
    if (associated(TAU_DP_ptr))      TAU_DP_ptr      = TAU_DP
    if (associated(TAU_MD_ptr))      TAU_MD_ptr      = TAU_MD
    if (associated(WQT_DC_ptr))      WQT_DC_ptr      = WQT_DC
    if (associated(ENTR_DP_ptr))     ENTR_DP_ptr     = ENTR_DP
    if (associated(ENTR_MD_ptr))     ENTR_MD_ptr     = ENTR_MD
    if (associated(ENTR_SH_ptr))     ENTR_SH_ptr     = ENTR_SH
    if (associated(SGS_VVEL_DP_ptr)) SGS_VVEL_DP_ptr = SGS_VVEL_DP
    if (associated(SGS_VVEL_MD_ptr)) SGS_VVEL_MD_ptr = SGS_VVEL_MD
    if (associated(SGS_VVEL_SH_ptr)) SGS_VVEL_SH_ptr = SGS_VVEL_SH
    if (associated(CNV_TOPP_DP_ptr)) CNV_TOPP_DP_ptr = CNV_TOPP_DP
    if (associated(CNV_TOPP_MD_ptr)) CNV_TOPP_MD_ptr = CNV_TOPP_MD
    if (associated(CNV_TOPP_SH_ptr)) CNV_TOPP_SH_ptr = CNV_TOPP_SH
    if (associated(DNDCNV_ptr))      DNDCNV_ptr      = DNDCNV
    if (associated(DNICNV_ptr))      DNICNV_ptr      = DNICNV

    endif ! Alarm ringing

    endif ! USE_PYMOIST_GF2020

    ! Add tendencies to the moist import state
    !--------------------------------------------------------------
    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP SHARED(IM, JM, LM, U, DUDT_DC, MOIST_DT, V, DVDT_DC, Q, DQVDT_DC, &
    !$OMP        T, DTDT_DC, QLCN, DQLDT_DC, QICN, DQIDT_DC, CLCN, DQADT_DC) &
    !$OMP PRIVATE(I, J, L)
    do L = 1, LM
       do J = 1, JM
          !DIR$ IVDEP
          do I = 1, IM
             U(I,J,L) = U(I,J,L) + DUDT_DC(I,J,L) * MOIST_DT
             V(I,J,L) = V(I,J,L) + DVDT_DC(I,J,L) * MOIST_DT
             Q(I,J,L) = Q(I,J,L) + DQVDT_DC(I,J,L) * MOIST_DT
             T(I,J,L) = T(I,J,L) + DTDT_DC(I,J,L) * MOIST_DT
             
             ! Add QI/QL/CL tendencies
             QLCN(I,J,L) = QLCN(I,J,L) + DQLDT_DC(I,J,L) * MOIST_DT
             QICN(I,J,L) = QICN(I,J,L) + DQIDT_DC(I,J,L) * MOIST_DT
             CLCN(I,J,L) = MAX(0.0, MIN(CLCN(I,J,L) + DQADT_DC(I,J,L) * MOIST_DT, 1.0))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

! Cleanup negative water species
! ------------------------------
    call MAPL_StateGetPointer(EXPORT,   DQVDT_FILL,   'DQVDT_FILL_DC', _RC)
    call MAPL_StateGetPointer(EXPORT, DQLLSDT_FILL, 'DQLLSDT_FILL_DC', _RC)
    call MAPL_StateGetPointer(EXPORT, DQLCNDT_FILL, 'DQLCNDT_FILL_DC', _RC)
    call MAPL_StateGetPointer(EXPORT, DQILSDT_FILL, 'DQILSDT_FILL_DC', _RC)
    call MAPL_StateGetPointer(EXPORT, DQICNDT_FILL, 'DQICNDT_FILL_DC', _RC)
    if (REPORT_GF_NEGATIVES) then
      call FILLQ2ZERO( Q       , MASS, DT=MOIST_DT, DQDT=  DQVDT_FILL, WARNING_LABEL="QV   After GF DeepCu", VM=VMG, _RC)
      call FILLQ2ZERO( QLLS    , MASS, DT=MOIST_DT, DQDT=DQLLSDT_FILL, WARNING_LABEL="QLLS After GF DeepCu", VM=VMG, _RC)
      call FILLQ2ZERO( QLCN    , MASS, DT=MOIST_DT, DQDT=DQLCNDT_FILL, WARNING_LABEL="QLCN After GF DeepCu", VM=VMG, _RC)
      call FILLQ2ZERO( QILS    , MASS, DT=MOIST_DT, DQDT=DQILSDT_FILL, WARNING_LABEL="QILS After GF DeepCu", VM=VMG, _RC)
      call FILLQ2ZERO( QICN    , MASS, DT=MOIST_DT, DQDT=DQICNDT_FILL, WARNING_LABEL="QICN After GF DeepCu", VM=VMG, _RC)
    else
      call FILLQ2ZERO( Q       , MASS, DT=MOIST_DT, DQDT=  DQVDT_FILL, VM=VMG, _RC)
      call FILLQ2ZERO( QLLS    , MASS, DT=MOIST_DT, DQDT=DQLLSDT_FILL, VM=VMG, _RC)
      call FILLQ2ZERO( QLCN    , MASS, DT=MOIST_DT, DQDT=DQLCNDT_FILL, VM=VMG, _RC)
      call FILLQ2ZERO( QILS    , MASS, DT=MOIST_DT, DQDT=DQILSDT_FILL, VM=VMG, _RC)
      call FILLQ2ZERO( QICN    , MASS, DT=MOIST_DT, DQDT=DQICNDT_FILL, VM=VMG, _RC)
    endif

    if (DEBUG_TQ_ERRORS) then
        do L=1,LM
          do J=1,JM
           do I=1,IM
             if (T(I,J,L) > 333.0) then
                 print *, "Temperature spike detected : ", T(I,J,L)
                 print *, "         GF Temp Increment : ", DTDT_DC(I,J,L)*MOIST_DT
                 print *, "    AFTER Grell-Freitas DeepCu"
                 print *, "  Latitude       =", LATS(I,J)*180.0/MAPL_PI
                 print *, "  Longitude      =", LONS(I,J)*180.0/MAPL_PI
                 print *, "  Pressure (mb)  =", 0.5*(PLE(I,J,L-1) + PLE(I,J,L))/100.0
                 print *, "                       CLLS=",  CLLS(I,J,L),             "CLCN=", CLCN(I,J,L)
                 print *, "  QV=",     Q(I,J,L), "  QL=",  QLLS(I,J,L)+QLCN(I,J,L), "  QI=", QILS(I,J,L)+QICN(I,J,L)
             endif
           end do ! IM loop
         end do ! JM loop
       end do ! LM loop
    endif

    call MAPL_GridCompTimerStop(GC, "--GF", _RC)

end subroutine GF_Run

subroutine GF_Finalize(gc, import, export, rc)

  type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
  type(ESMF_State),    intent(inout) :: IMPORT ! Import state
  type(ESMF_State),    intent(inout) :: EXPORT ! Export state
  integer, optional,   intent(  out) :: RC     ! Error code

  type (MAPL_MetaComp), pointer   :: MAPL

  ! Get my internal MAPL_Generic state
  !-----------------------------------
  call MAPL_GetObjectFromGC ( GC, MAPL, _RC)


  if (USE_PYMOIST_GF2020) then
    call MAPL_pybridge_gcfinalize( "pyMoist.fortran.param_interfaces.convection.GF2020_interface", MAPL, IMPORT, EXPORT )
  endif

end subroutine GF_Finalize

end module GEOS_GF_InterfaceMod
