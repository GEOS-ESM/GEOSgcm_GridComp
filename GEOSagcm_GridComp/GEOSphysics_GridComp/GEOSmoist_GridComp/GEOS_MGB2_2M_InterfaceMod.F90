! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_MGB2_2M_InterfaceMod -- A Module to interface with the
!   MGB2_2M cloud microphysics

module GEOS_MGB2_2M_InterfaceMod

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use GEOSmoist_Process_Library
  use cldwat2m_micro
  use aer_cloud
  use micro_mg3_0

  implicit none

  integer :: MGVERSION

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
         character(len=ESMF_MAXSTR) :: NCPL
         character(len=ESMF_MAXSTR) :: NCPI
         character(len=ESMF_MAXSTR) :: NRAIN
         character(len=ESMF_MAXSTR) :: NSNOW
         character(len=ESMF_MAXSTR) :: NGRAUPEL
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  character(len=ESMF_MAXSTR)        :: COMP_NAME

  ! Local resource variables
  real    :: TURNRHCRIT_PARAM
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
  logical :: LMELTFRZ
  logical :: USE_AV_V 
  logical :: PREEXISITING_ICE
  logical :: USE_BERGERON

  integer :: CCN_PARAM, IN_PARAM, Immersion_PARAM, WSUB_OPTION
  real  :: DCS, WBFFACTOR, NC_CST, NI_CST, NG_CST, MUI_CST, &
           LCCIRRUS, UISCALE, LIU_MU, NPRE_FRAC, QCVAR_CST, &
           AUT_SCALE, TS_AUTO_ICE, &
           FDROP_DUST, FDROP_SOOT, &
           DUST_INFAC, ORG_INFAC, BC_INFAC, SS_INFAC, &
           MTIME,MINCDNC, ACC_ENH, ACC_ENH_ICE, DT_MICRO, URSCALE

  public :: MGB2_2M_Setup, MGB2_2M_Initialize, MGB2_2M_Run
  public :: MGVERSION

contains

subroutine MGB2_2M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "MGB2_2M_Setup"
    
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam
    
    call ESMF_ConfigGetAttribute( CF, MGVERSION, Label="MGVERSION:",  default=3, __RC__)


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
      FRIENDLIES%NCPI     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%NCPL     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%NRAIN    = "DYNAMICS:TURBULENCE"
      FRIENDLIES%NSNOW    = "DYNAMICS:TURBULENCE"
      FRIENDLIES%NGRAUPEL = "DYNAMICS:TURBULENCE"

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

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPL',                                       &
         LONG_NAME  ='particle_number_for_liquid_cloud',           &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NCPL),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 50.0e6 ,                             __RC__  )  
                                                              

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPI',                                       &
         LONG_NAME  ='particle_number_for_ice_cloud',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NCPI),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 1.0e3,                               __RC__  )  
                                                       

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NRAIN',                                      &
         LONG_NAME  ='particle_number_for_rain',                   &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NRAIN),                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 0.0 ,                                __RC__  )  
                                                              

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NSNOW',                                      &
         LONG_NAME  ='particle_number_for_snow',                   &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NSNOW),                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 0.0,                                 __RC__  )  
                                                       

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NGRAUPEL',                                   &
         LONG_NAME  ='particle_number_for_graupel',                &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NGRAUPEL),                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 0.0,                                 __RC__  )  
    
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
    
    
    
    call MAPL_TimerAdd(GC, name="--MGB2_2M", __RC__)
    VERIFY_(STATUS)

end subroutine  MGB2_2M_Setup



subroutine MGB2_2M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_State)                   :: INTERNAL

    type (ESMF_Alarm   )                :: ALARM
    type (ESMF_TimeInterval)            :: TINT
    real(ESMF_KIND_R8)                  :: DT_R8
    real                                :: DT_MOIST

    real, pointer, dimension(:,:,:)     :: Q
    real, pointer, dimension(:,:,:)     :: QLLS, QLCN, QILS, QICN, CLLS, CLCN
    real, pointer, dimension(:,:,:)     :: QRAIN, QSNOW, QGRAUPEL
    real, pointer, dimension(:,:,:)     :: NCPL, NCPI, NRAIN, NSNOW, NGRAUPEL

    logical  :: nccons, nicons, ngcons, do_graupel
    real(ESMF_KIND_R8)  Dcsr8, micro_mg_berg_eff_factor_in, ncnstr8, ninstr8, ngnstr8, mui_cnstr8

    call MAPL_GetResource( MAPL, LMELTFRZ,          Label="MELTFRZ:",          default=.TRUE.,  __RC__ )
    call MAPL_GetResource( MAPL, PREEXISITING_ICE,  Label='PREEXISITING_ICE:', default=.FALSE., __RC__ )
    call MAPL_GetResource( MAPL, USE_BERGERON,      Label='USE_BERGERON:',     default=.TRUE.,  __RC__ )

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )
    call MAPL_Get( MAPL, &
         RUNALARM = ALARM,             &
         INTERNAL_ESMF_STATE=INTERNAL, &
         __RC__ )

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, RC=STATUS); VERIFY_(STATUS)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,RC=STATUS); VERIFY_(STATUS)
    DT_MOIST = DT_R8

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL, 'NGRAUPEL', RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, SH_MD_DP        , 'SH_MD_DP:'        , DEFAULT= .TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, DBZ_LIQUID_SKIN , 'DBZ_LIQUID_SKIN:' , DEFAULT= 0     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, TURNRHCRIT_PARAM, 'TURNRHCRIT:'      , DEFAULT= -9999., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, PDFSHAPE        , 'PDFSHAPE:'        , DEFAULT= 1     , RC=STATUS); VERIFY_(STATUS)
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

    !2M==tuning and options======
    call MAPL_GetResource(MAPL, UISCALE,        'UISCALE:',        DEFAULT= 1.0,     __RC__) !Scaling factor for sed vel of ice      
    call MAPL_GetResource(MAPL, LIU_MU,         'LIU_MU:',         DEFAULT= 2.0,     __RC__) !Liu autoconversion parameter
    call MAPL_GetResource(MAPL, NPRE_FRAC,      'NPRE_FRAC:',      DEFAULT= -1.0,    __RC__) !Fraction of preexisting ice affecting ice nucleationn            
    call MAPL_GetResource(MAPL, USE_AV_V,       'USE_AV_V:',       DEFAULT= .TRUE.,  __RC__) !Set to > 0 to use an average velocity for activation
    call MAPL_GetResource(MAPL, AUT_SCALE,      'AUT_SCALE:',      DEFAULT= 0.5,     __RC__) !scale factor for critical size for drizzle
    call MAPL_GetResource(MAPL, TS_AUTO_ICE,    'TS_AUTO_ICE:',    DEFAULT= 360.,    __RC__) !Ice autoconversion time scale
    call MAPL_GetResource(MAPL, CCN_PARAM,      'CCNPARAM:',       DEFAULT= 2,       __RC__) !CCN activation param
    call MAPL_GetResource(MAPL, IN_PARAM,       'INPARAM:',        DEFAULT= 6,       __RC__) !IN param
    call MAPL_GetResource(MAPL, Immersion_PARAM,'Immersion_PARAM:',DEFAULT= 6,       __RC__) !Immersion param
    call MAPL_GetResource(MAPL, ACC_ENH,        'ACC_ENH:',        DEFAULT= 1.0,     __RC__) !accretion rain-liquid scaling for MG2
    call MAPL_GetResource(MAPL, ACC_ENH_ICE,    'ACC_ENH_ICE:',    DEFAULT= 1.0,     __RC__) !accretion snow-ice scaling for MG2
    call MAPL_GetResource(MAPL, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.5,     __RC__) !Fraction of dust within droplets for immersion freezing
    call MAPL_GetResource(MAPL, FDROP_SOOT,     'FDROP_SOOT:',     DEFAULT= 0.05,    __RC__) !Fraction of soot within droplets for immersion freezing        
    call MAPL_GetResource(MAPL, MINCDNC,        'MINCDNC:',        DEFAULT= 25.0,    __RC__) !min nucleated droplet conc. cm-3
    call MAPL_GetResource(MAPL, MTIME,          'MTIME:',          DEFAULT= -1.0,    __RC__) !Mixing time scale for aerosol within the cloud. Default is time step
    call MAPL_GetResource(MAPL, LCCIRRUS,       'LCCIRRUS:',       DEFAULT= 500.0,   __RC__) !Characteristic Length (m) of high freq gravity waves    
    call MAPL_GetResource(MAPL, QCVAR_CST,      'QCVAR_CST:',      DEFAULT= -1.,     __RC__) !Characteristic Length (m) of high freq gravity waves 
    call MAPL_GetResource(MAPL, DUST_INFAC,     'DUST_INFAC:',     DEFAULT= 1.0,     __RC__) !scalings for the INP concentrations
    call MAPL_GetResource(MAPL, BC_INFAC,       'BC_INFAC:',       DEFAULT= 0.1,     __RC__)
    call MAPL_GetResource(MAPL, ORG_INFAC,      'ORG_INFAC:',      DEFAULT= 1.0,     __RC__)
    call MAPL_GetResource(MAPL, SS_INFAC,       'SS_INFAC:',       DEFAULT= 1.0,     __RC__)    
    call MAPL_GetResource(MAPL, DT_MICRO,       'DT_MICRO:',       DEFAULT= 300.0,   __RC__) !time step of the microphysics substepping (s) (MG2) (5 min)
    call MAPL_GetResource(MAPL, URSCALE,        'URSCALE:',        DEFAULT= 1.0,     __RC__) !Scaling factor for sed vel of rain    
    call MAPL_GetResource(MAPL, DCS,            'DCS:'    ,        DEFAULT=200.0e-6, __RC__) !ice/snow separation diameter   
    Dcsr8 = DCS
    call MAPL_GetResource(MAPL, WBFFACTOR,      'WBFFACTOR:',      DEFAULT= 0.1 ,    __RC__) !scaling for the Bergeron-Findeinsen process rate    
    micro_mg_berg_eff_factor_in = WBFFACTOR
    call MAPL_GetResource(MAPL, NC_CST ,        'NC_CST:' ,        DEFAULT= 0.0 ,    __RC__) !constant nd (set if greather than zero)      
    call MAPL_GetResource(MAPL, NI_CST ,        'NI_CST:' ,        DEFAULT= 0.0 ,    __RC__) !constant nd (set if greather than zero)     
    call MAPL_GetResource(MAPL, NG_CST ,        'NG_CST:' ,        DEFAULT= 0.0 ,    __RC__) !constant ng (set if greather than zero)     
    call MAPL_GetResource(MAPL, MUI_CST,        'MUI_CST:',        DEFAULT= -1.0 ,   __RC__) !constant ng (set if greather than zero) 
    call MAPL_GetResource(MAPL, WSUB_OPTION,    'WSUB_OPTION:',    DEFAULT= 1,       __RC__) !0- param 1- Use Wsub climatology 2-Wnet
    mui_cnstr8 =  MUI_CST
    ncnstr8 = NC_CST
    if  (NC_CST .gt. 0.0)  nccons =.true.
    ninstr8 = NI_CST
    if  (NI_CST .gt. 0.0)  nicons =.true.
    ngnstr8 = NG_CST
    if  (NG_CST .gt. 0.0)  ngcons =.true.
    !============
   
    if  (MGVERSION .gt. 1) then
        do_graupel = .false.
        if (MGVERSION .gt. 2) do_graupel = .true.
        call micro_mg_init(Dcsr8, do_graupel,  micro_mg_berg_eff_factor_in, &
                           nccons, nicons, ncnstr8, ninstr8, ngcons, ngnstr8, mui_cnstr8)
    else
        call ini_micro(Dcsr8, micro_mg_berg_eff_factor_in, &
                       nccons, nicons, ncnstr8, ninstr8, 2.0)
    end if
   
    call aer_cloud_init()

    call WRITE_PARALLEL ("INITIALIZED MGB2_2M microphysics in non-generic GC INIT")

                                 CCW_EVAP_EFF = 4.e-3
    call MAPL_GetResource( MAPL, CCW_EVAP_EFF, 'CCW_EVAP_EFF:', DEFAULT= CCW_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

                                 CCI_EVAP_EFF = 4.e-3
    call MAPL_GetResource( MAPL, CCI_EVAP_EFF, 'CCI_EVAP_EFF:', DEFAULT= CCI_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=  500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_EXP, 'CNV_FRACTION_EXP:', DEFAULT=    1.0, RC=STATUS); VERIFY_(STATUS)

end subroutine MGB2_2M_Initialize

subroutine MGB2_2M_Run  (GC, IMPORT, EXPORT, CLOCK, RC)
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

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
    real, pointer, dimension(:,:,:) :: NCPL, NCPI, NRAIN, NSNOW, NGRAUPEL
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W, KH
    real, pointer, dimension(:,:)   :: AREA, FRLAND, TS, SH, EVAP, KPBL_SC
    real, pointer, dimension(:,:,:) :: SL2, SL3, QT2, QT3, W2, W3, SLQT, WQT, WQL, WSL
    real, pointer, dimension(:,:,:) :: WTHV2
    real, pointer, dimension(:,:,:) :: OMEGA, WSUB_CLIM
    real, pointer, dimension(:,:,:) :: RADLW, RADSW
    
    ! Local
    real, allocatable, dimension(:,:,:) :: U0, V0
    real, allocatable, dimension(:,:,:) :: PLEmb, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  ZL0, GZLO
    real, allocatable, dimension(:,:,:) :: DZET, DP, MASS, iMASS
    real, allocatable, dimension(:,:,:) :: DQST3, QST3
    real, allocatable, dimension(:,:,:) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                                           DQSDTmic, DQGDTmic, DQADTmic, &
                                           DUDTmic,  DVDTmic,  DTDTmic
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D
 
    ! Exports
    real, pointer, dimension(:,:  ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ) :: LS_PRCP, LS_SNR, CNV_FRC, SRF_TYPE, ICE, FRZR 
    real, pointer, dimension(:,:,:) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:) ::  DUDT_macro,  DVDT_macro,  DTDT_macro
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro
    real, pointer, dimension(:,:,:) ::  DUDT_micro,  DVDT_micro,  DTDT_micro
    real, pointer, dimension(:,:,:) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, pointer, dimension(:,:,:) :: CLDREFFL, CLDREFFI, CLDREFFR, CLDREFFS, CLDREFFG
    real, pointer, dimension(:,:,:) :: EVAPC, SUBLC
    real, pointer, dimension(:,:,:) :: RHX, REV_LS, RSU_LS
    real, pointer, dimension(:,:,:) :: PFL_LS, PFL_AN
    real, pointer, dimension(:,:,:) :: PFI_LS, PFI_AN
    real, pointer, dimension(:,:,:) :: PDF_A, PDFITERS
    real, pointer, dimension(:,:,:) :: RHCRIT3D
    real, pointer, dimension(:,:)   :: EIS, LTS
    real, pointer, dimension(:,:)   :: DBZ_MAX, DBZ_1KM, DBZ_TOP, DBZ_M10C
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D
#ifdef PDFDIAG
    real, pointer, dimension(:,:,:) :: PDF_W1, PDF_W2, PDF_SIGW1, PDF_SIGW2,     &
                                       PDF_QT1, PDF_QT2, PDF_SIGQT1, PDF_SIGQT2, &
                                       PDF_TH1, PDF_TH2, PDF_SIGTH1, PDF_SIGTH2, &
                                       PDF_RQTTH, PDF_RWTH, PDF_RWQT
#endif    
    
    !2m
    real, pointer, dimension(:,:,:) :: SC_ICE, CDNC_NUC, INC_NUC, &
       CFICE, CFLIQ, WSUB, CCN01, CCN04, CCN1, &
       NHET_DEP, SIGW_RC, SIGW_GW, SIGW_CNV, SIGW_TURB, &
       BERG, BERGS, MELT, DNHET_CT, QCRES, QIRES, AUTICE, FRZPP_LS, &
       SNOWMELT_LS, DNCNUC, DNCSUBL, DNCHMSPLIT, DNCAUTICE, DNCACRIS, DNDCCN, &
       DNDACRLS, DNDACRLR, DNDEVAPC, DNDAUTLIQ, &
       DNHET_IMM, RHCmicro, &
       NWFA
    real, pointer, dimension(:,:)   :: QCVAR
    real, allocatable, dimension(:,:,:) :: QCNTOT, CFX,  QTOT, &
       QL_TOT, QI_TOT, ACIL_LS_X, ACIL_AN_X, ACLL_LS_X, ACLL_AN_X, DLPDF_X, DIPDF_X, DLFIX_X, DIFIX_X, &
       AUT_X, SDM_X, FRZ_TT_X, FRZ_PP_X,  FQA !check how much of these we are actually using
    real, allocatable, dimension(:, :)  :: CLDREFFI_TOP_X, CLDREFFL_TOP_X,  NCPL_TOP_X, NCPI_TOP_X, NCPL_CLDBASEX, uwind_gw
    integer, allocatable, dimension(:, :)   ::   KLCL

    ! Local variables
    real    :: facEIS
    real    :: minrhcrit, turnrhcrit, ALPHA, RHCRIT
    integer :: IM,JM,LM
    integer :: I, J, L, K
  

    integer :: num_steps_micro,  pcnst, n_modes, kbmin, kcldtop, kcldbot , &
                NAUX, kcldtopcvn, nbincontactdust, index, K0, KCBLMIN, i_src_mode, i_dst_mode
                  
    real, parameter :: pmin_trop = 10.0 !mbar minimum pressure to do cloud microphysics
    integer, parameter :: KMIN_TROP =  25

    REAL, allocatable, dimension(:,:) :: SCICE_tmp, FQA_tmp, cfaux 
                        
    real (ESMF_KIND_R8), dimension(3)       :: ccn_diag
    real (ESMF_KIND_R8), allocatable, dimension(:,:,:) :: rndstr8,naconr8  !Assume maximum 5 dust bins
    real (ESMF_KIND_R8), dimension(1)       :: prectr8, precir8            
    real (ESMF_KIND_R8)  :: disp_liu, ui_scale, dcrit, tfreez,  &
                            ts_autice, mtimesc, ur_scale
    real (ESMF_KIND_R8), allocatable, dimension(:,:)  :: ttendr8, qtendr8, cwtendr8, &
           cldor8,  rpdelr8, zmr8, omegr8, rhdfdar8, rhu00r8, ficer8 , &
             qilsr8,  &
            pintr8, kkvhr8, rflxr8,  sflxr8, lflxr8, iflxr8, gflxr8,  &
            ter8,qvr8, qcr8,qir8, ncr8,nir8, qrr8,qsr8, nrr8,nsr8, &
            qgr8,ngr8, relvarr8,accre_enhanr8, plevr8, pdelr8, cldfr8,liqcldfr8, &
            icecldfr8,qsatfacr8, qcsinksum_rate1ordr8, naair8, npccninr8, &
            tlatr8, qvlatr8, qctendr8, qitendr8, nctendr8, nitendr8, qrtendr8, &
            qstendr8, qgtendr8, nrtendr8, nstendr8, ngtendr8, effcr8,effc_fnr8, &
            effir8, sadicer8, sadsnowr8, nevaprr8, evapsnowr8, am_evp_str8, &
            prainr8,prodsnowr8, cmeoutr8, deffir8, pgamradr8,lamcradr8, qsoutr8,&
            dsoutr8, qgoutr8, ngoutr8,dgoutr8, qroutr8, reff_rainr8,reff_snowr8, &
            reff_graur8, qcsevapr8, qisevapr8, qvresr8, cmeioutr8, vtrmcr8, &
            vtrmir8, umrr8,umsr8, umgr8,qgsedtendr8, qcsedtenr8, qisedtenr8, &
            qrsedtenr8, qssedtenr8, praor8, prcor8, mnucccor8, mnucctor8, &
            msacwior8, psacwsor8, bergsor8,bergor8, meltor8,homoor8, qcresor8,&
            prcior8, praior8, qirestotr8,mnuccrtotr8, mnuccritotr8, pracstotr8, &
            meltsdtr8,frzrdtr8, mnuccdor8, pracgtotr8,psacwgtotr8, pgsacwtotr8, &
            pgracstotr8, prdgtotr8, qmultgtotr8, qmultrgtotr8,psacrtotr8, &
            npracgtotr8, nscngtotr8,ngracstotr8, nmultgtotr8, nmultrgtotr8,&
            npsacwgtotr8, nroutr8, nsoutr8, reflr8,areflr8, areflzr8, freflr8, &
            csrflr8, acsrflr8, fcsrflr8, rercldr8, ncair8, ncalr8, qrout2r8, &
            qsout2r8, nrout2r8, nsout2r8, drout2r8, dsout2r8, qgout2r8, ngout2r8,&
            dgout2r8,freqgr8, freqsr8,freqrr8, nficer8,qcratr8, & 
            !errstring, 
            ! Below arguments are "optional" (pass null pointers to omit).  
            tnd_qsnow, tnd_nsnow, re_ice, &
            prer_evap, frzimmr8,frzcntr8, frzdepr8, & ! contact is not passed since it depends on the droplet size dist 
            nsootr8, rnsootr8, & ! soot for contact IN 
            npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, & ! Number tendencies for liquid 
            npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, & 
            nsacwior8, mnuccror8,pracsor8, qiresor8, rate1ord_cw2pr, accre_enhan_icer8
      real   :: tausurf_gw, aux1,aux2,aux3, npre, dpre, nact, xscale
      real (ESMF_KIND_R8)   :: autscx   
      integer,  parameter :: ncolmicro = 1

      type  (AerProps) :: AeroAux, AeroAux_b
                                    
      call init_Aer(AeroAux)
      call init_Aer(AeroAux_b)

    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ) 
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--MGB2_2M",RC=STATUS)

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
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
    call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL, 'NGRAUPEL', RC=STATUS); VERIFY_(STATUS)

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
    call MAPL_GetPointer(IMPORT, SH,      'SH'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP,    'EVAP'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBL_SC, 'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WSUB_CLIM,'WSUB_CLIM', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RADLW,  'RADLW'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RADSW,  'RADSW'    , RC=STATUS); VERIFY_(STATUS)

    ! Exports that require memory for calculations
    call MAPL_GetPointer(EXPORT, SIGW_RC,    'SIGW_RC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, WSUB,       'WSUB'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, INC_NUC,    'INC_NUC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNHET_IMM,  'DNHET_IMM'   , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NHET_DEP,   'NHET_DEP'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SC_ICE,     'SC_ICE'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CDNC_NUC,   'CDNC_NUC'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CFICE,      'CFICE'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CFLIQ,      'CFLIQ'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RHCmicro,   'RHCmicro'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QCVAR,      'QCVAR'       , ALLOC=.TRUE., __RC__)

    ! Diagnostic exports
    call MAPL_GetPointer(EXPORT, CCN01,      'CCN01'       , __RC__)
    call MAPL_GetPointer(EXPORT, CCN04,      'CCN04'       , __RC__)
    call MAPL_GetPointer(EXPORT, CCN1,       'CCN1'        , __RC__)
    call MAPL_GetPointer(EXPORT, BERG,       'BERG'        , __RC__)
    call MAPL_GetPointer(EXPORT, BERGS,      'BERGS'       , __RC__) 
    call MAPL_GetPointer(EXPORT, MELT,       'MELT'        , __RC__)         
    call MAPL_GetPointer(EXPORT, CLDREFFR,   'RR'          , __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFS,   'RS'          , __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFG,   'RG'          , __RC__)
    call MAPL_GetPointer(EXPORT, DNHET_CT,   'DNHET_CT'    , __RC__)
    call MAPL_GetPointer(EXPORT, QCRES,      'QCRES'       , __RC__)
    call MAPL_GetPointer(EXPORT, QIRES,      'QIRES'       , __RC__)
    call MAPL_GetPointer(EXPORT, AUTICE,     'AUTICE'      , __RC__)
    call MAPL_GetPointer(EXPORT, FRZPP_LS ,  'FRZPP_LS'    , __RC__)
    call MAPL_GetPointer(EXPORT, SNOWMELT_LS,'SNOWMELT_LS' , __RC__)
    call MAPL_GetPointer(EXPORT, DNCNUC,     'DNCNUC'      , __RC__)
    call MAPL_GetPointer(EXPORT, DNCSUBL,    'DNCSUBL'     , __RC__)
    call MAPL_GetPointer(EXPORT, DNCHMSPLIT, 'DNCHMSPLIT'  , __RC__)
    call MAPL_GetPointer(EXPORT, DNCAUTICE,  'DNCAUTICE'   , __RC__)
    call MAPL_GetPointer(EXPORT, DNCACRIS,   'DNCACRIS'    , __RC__)
    call MAPL_GetPointer(EXPORT, DNDCCN,     'DNDCCN'      , __RC__)
    call MAPL_GetPointer(EXPORT, DNDACRLS,   'DNDACRLS'    , __RC__)
    call MAPL_GetPointer(EXPORT, DNDACRLR,   'DNDACRLR'    , __RC__)
    call MAPL_GetPointer(EXPORT, DNDEVAPC,   'DNDEVAPC'    , __RC__)
    call MAPL_GetPointer(EXPORT, DNDAUTLIQ,  'DNDAUTLIQ'   , __RC__)

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
    ALLOCATE ( DP   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( iMASS(IM,JM,LM  ) )
    ALLOCATE ( DQST3(IM,JM,LM  ) )
    ALLOCATE (  QST3(IM,JM,LM  ) )
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
     ! 2D Variables
    ALLOCATE ( KLCL    (IM,JM) )
    ALLOCATE ( TMP2D   (IM,JM) )

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
    call MAPL_GetPointer(EXPORT, WTHV2,    'WTHV2'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WQL,      'WQL'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PDFITERS, 'PDFITERS', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Unused Exports (forced to 0.0)
    call MAPL_GetPointer(EXPORT, PTR2D,  'CN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'CN_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    ! Lowe tropospheric stability and estimated inversion strength
    call MAPL_GetPointer(EXPORT, LTS,   'LTS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EIS,   'EIS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    KLCL = FIND_KLCL( T, Q, PLmb, IM, JM, LM )
    call MAPL_GetPointer(EXPORT, PTR2D, 'ZLCL', RC=STATUS); VERIFY_(STATUS)
    if (associated(PTR2D)) then
      do J=1,JM
         do I=1,IM
           PTR2D(I,J) = ZL0(I,J,KLCL(I,J))
         end do
      end do
    endif
    TMP3D = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
    call FIND_EIS(T/TMP3D, QST3, T, ZL0, PLEmb, KLCL, IM, JM, LM, LTS, EIS)


   ! ==========================================================================================
   ! ========================Activate the aerosols ============================================
    call MAPL_TimerOn(MAPL,"---ACTIV_2MOM", RC=STATUS); VERIFY_(STATUS)
  !!=============== vertical velocity variance
    !- Determine which W is proper import
    if (all(W == 0.0)) then
       SIGW_RC = -1*OMEGA/(MAPL_GRAV*100.*PLmb/(MAPL_RDRY*T*(1.0+MAPL_VIREPS*Q)))
    else
       SIGW_RC = W
    endif
    SIGW_RC = SIGW_RC + (RADLW + RADSW)*MAPL_CP/MAPL_GRAV
    SC_ICE = 1.0
    if (WSUB_OPTION /= 1) then ! use parameterization from Barahona et al. GMD. 2014 (Appendix)
#ifdef SKIP
          call MAPL_GetPointer(EXPORT, SIGW_GW,    'SIGW_GW'     , ALLOC=.TRUE., __RC__)
          call MAPL_GetPointer(EXPORT, SIGW_CNV,   'SIGW_CNV'    , ALLOC=.TRUE., __RC__)
          call MAPL_GetPointer(EXPORT, SIGW_TURB,  'SIGW_TURB'   , ALLOC=.TRUE., __RC__)
          do J=1,JM
             do I=1,IM
                uwind_gw(1,1:LM)           = min(0.5*SQRT( U(I,J,1:LM)**2+  V(I,J,1:LM)**2), 50.0)
                tausurf_gw   = min(0.5*SQRT(TAUOROX(I , J)**2+TAUOROY(I , J)**2), 10.0) !limit to a very high value
                call vertical_vel_variance(T(I,J,1:LM), TKE(I,J,1:LM), 100.0*PLmb(I,J,1:LM), PLE(I,J,0:LM), uwind_gw(1,1:LM), &
                                            tausurf_gw, AIRDEN(I,J,1:LM), LM, LCCIRRUS, -SH (i,j), -EVAP(i,j), ZL0(I, J, NINT(KPBL_SC(I,J))), &
                                            SIGW_GW (I, J, 1:LM), SIGW_TURB (I, J, 1:LM), SIGW_CNV (I, J, 1:LM), WSUB (I, J, 1:LM), &
                                            SIGW_RC(I, J, 1:LM))
             end do        
          end do 
#endif 
    else  !WSUB climatology 
       WSUB = WSUB_CLIM
    endif
    !- Activation
    do J=1,JM
       do I=1,IM
          kbmin= min(NINT(KPBL_SC(I,J)), LM-1)-2         
          dpre = 1.0e-9
                               npre = NPRE_FRAC      
          if (NPRE_FRAC < 0.0) npre = CNV_FRC(I,J)*ABS(NPRE_FRAC) + (1-CNV_FRC(I,J))*0.05

          do K = KMIN_TROP, LM-1 !limit to troposphere and no activation at the surface
             npre = npre*NCPI(I,J,K)
             if ((npre > 0.0) .and. (QICN(I,J,K)+QILS(I,J,K) > 0.)) &
                            dpre = ((QICN(I,J,K)+QILS(I,J,K))/(5400.0*npre*MAPL_PI))**0.33 !Assume exponential distribution

             !!Subroutine aerosol_activate contains the CCN activation and ice nucleation parameterizations. Lives in aer_cloud.F90.
             call aerosol_activate(T(I,J,K), 100.*PLmb(I,J,K), WSUB(I,J,K), SIGW_RC(I,J,K), AeroProps(I,J,K), &
                                   npre, dpre, USE_AV_V, CCN_PARAM, IN_PARAM, FDROP_DUST, FDROP_SOOT,  &
                                   DUST_INFAC, BC_INFAC, ORG_INFAC, SS_INFAC, Immersion_PARAM, & 
                                   ccn_diag, nact, &
                                   INC_NUC(I,J,K), DNHET_IMM(I,J,K), NHET_DEP(I,J,K), SC_ICE(I,J,K))
             if (T(I,J,K) > 238.0) then
                SC_ICE(I,J,K) = 1.0
             endif
             SC_ICE(I,J,K) = MIN(MAX(SC_ICE(I,J,K), 1.0), 1.8)
    
             ! diagnostics
             if (associated(CCN01)) CCN01(I,J,K) = max(ccn_diag(1), 0.0)    
             if (associated(CCN04)) CCN04(I,J,K) = max(ccn_diag(2), 0.0)    
             if (associated(CCN1 )) CCN1 (I,J,K) = max(ccn_diag(3), 0.0)    
    
             if (K .ge. kbmin-4) nact = max(nact, (1.0-CNV_FRC(I,J))*MINCDNC*1.e6)
             CDNC_NUC(I, J, K) = nact   
    
          end do 
       enddo
    enddo
    ! fill WSUB export diagnostic with W + 0.8*WSUB 
    WSUB = SIGW_RC + 0.8*WSUB
    call MAPL_TimerOff(MAPL,"---ACTIV_2MOM", RC=STATUS); VERIFY_(STATUS)

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
    DUDT_macro=U
    DVDT_macro=V
    DTDT_macro=T
    DQVDT_macro=Q
    DQLDT_macro=QLCN+QLLS
    DQIDT_macro=QICN+QILS
    DQADT_macro=CLCN+CLLS
    DQRDT_macro=QRAIN
    DQSDT_macro=QSNOW
    DQGDT_macro=QGRAUPEL

      ! Include shallow precip condensates if present
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          QRAIN = QRAIN + PTR3D*DT_MOIST
        endif
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then 
          QSNOW = QSNOW + PTR3D*DT_MOIST
        endif
       ! evap/subl/pdf
        call MAPL_GetPointer(EXPORT, RHCRIT3D,  'RHCRIT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        do L=1,LM
          do J=1,JM
           do I=1,IM
           ! Send the condensates through the pdf after convection
             facEIS = MAX(0.0,MIN(1.0,EIS(I,J)/10.0))**2
           ! determine combined minrhcrit in stable/unstable regimes
             minrhcrit  = (0.9)*(1.0-facEIS) + (0.95)*facEIS
             if (TURNRHCRIT_PARAM <= 0.0) then
              ! determine the turn pressure using the LCL
                turnrhcrit  = PLmb(I, J, KLCL(I,J)) - 250.0 ! 250mb above the LCL
             else
                turnrhcrit = TURNRHCRIT_PARAM
             endif
           ! Use Slingo-Ritter (1985) formulation for critical relative humidity
             RHCRIT = 1.0
           ! lower turn from maxrhcrit=1.0
             if (PLmb(i,j,l) .le. turnrhcrit) then
                RHCRIT = minrhcrit
             else
                if (L.eq.LM) then
                   RHCRIT = 1.0
                else
                   RHCRIT = minrhcrit + (1.0-minrhcrit)/(19.) * &
                           ((atan( (2.*(PLmb(i,j,l)-turnrhcrit)/(PLEmb(i,j,LM)-turnrhcrit)-1.) * &
                           tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                endif
             endif
           ! include grid cell area scaling and limit RHcrit to > 70%
             ALPHA = max(0.0,min(0.30, (1.0-RHCRIT)*SQRT(SQRT(AREA(I,J)/1.e10)) ) )
           ! fill RHCRIT export
             if (associated(RHCRIT3D)) RHCRIT3D(I,J,L) = 1.0-ALPHA
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
                      NCPL(I,J,L)    , &
                      NCPI(I,J,L)    , &
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
                      PREEXISITING_ICE, &
                      USE_BERGERON   , &
                      SC_ICE(I,J,L))
             RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
             if (LMELTFRZ) then
           ! meltfrz new condensates
             call MELTFRZ ( DT_MOIST     , &
                            CNV_FRC(I,J) , &
                            SRF_TYPE(I,J), &
                            T(I,J,L)     , &
                            QLCN(I,J,L)  , &
                            QICN(I,J,L) )
             call MELTFRZ ( DT_MOIST     , &
                            CNV_FRC(I,J) , &
                            SRF_TYPE(I,J), &
                            T(I,J,L)     , &
                            QLLS(I,J,L)  , &
                            QILS(I,J,L) )
             endif
           ! evaporation for CN
             if (CCW_EVAP_EFF > 0.0) then 
             RHCRIT = 1.0
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
                   NCPL(I,J,L)  , &
                   NCPI(I,J,L)  , &
                   QST3(I,J,L)  )
             EVAPC(I,J,L) = ( Q(I,J,L) - EVAPC(I,J,L) ) / DT_MOIST
             endif
           ! sublimation for CN
             if (CCI_EVAP_EFF > 0.0) then
             RHCRIT = 1.0
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
                   NCPL(I,J,L)  , &
                   NCPI(I,J,L)  , &
                   QST3(I,J,L)  )
             SUBLC(I,J,L) = ( Q(I,J,L) - SUBLC(I,J,L) ) / DT_MOIST
             endif
           ! cleanup clouds
             call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
           end do ! IM loop
         end do ! JM loop
       end do ! LM loop

    ! Update macrophysics tendencies
     DUDT_macro=( U         - DUDT_macro)/DT_MOIST
     DVDT_macro=( V         - DVDT_macro)/DT_MOIST
     DTDT_macro=( T         - DTDT_macro)/DT_MOIST
    DQVDT_macro=( Q         -DQVDT_macro)/DT_MOIST
    DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
    DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
    DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
    DQRDT_macro=( QRAIN     -DQRDT_macro)/DT_MOIST
    DQSDT_macro=( QSNOW     -DQSDT_macro)/DT_MOIST
    DQGDT_macro=( QGRAUPEL  -DQGDT_macro)/DT_MOIST
    call MAPL_TimerOff(MAPL,"---CLDMACRO")

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
    DQVDT_micro = Q
    DQLDT_micro = QLLS + QLCN
    DQIDT_micro = QILS + QICN
    DQRDT_micro = QRAIN
    DQSDT_micro = QSNOW
    DQGDT_micro = QGRAUPEL
    DQADT_micro = CLLS + CLCN
    DUDT_micro = U
    DVDT_micro = V
    DTDT_micro = T

! 1D
    allocate(ttendr8(1,LM), __STAT__)
    allocate(qtendr8(1,LM), __STAT__)
    allocate(cwtendr8(1,LM), __STAT__)
    allocate(cldor8(1,LM), __STAT__)
    allocate(rpdelr8(1,LM), __STAT__)
    allocate(zmr8(1,LM), __STAT__)
    allocate(omegr8(1,LM), __STAT__)
    allocate(rhdfdar8(1,LM), __STAT__)
    allocate(rhu00r8(1,LM), __STAT__)
    allocate(ficer8(1,LM), __STAT__)
    allocate(qilsr8(1,LM), __STAT__)
    allocate(uwind_gw(1,LM), __STAT__)
    allocate(SCICE_tmp(1,LM), __STAT__)
    allocate(FQA_tmp(1,LM), __STAT__)

    allocate(ter8(1,LM), __STAT__)
    allocate(qvr8(1,LM), __STAT__)
    allocate(qcr8(1,LM), __STAT__)
    allocate(qir8(1,LM), __STAT__)
    allocate(ncr8(1,LM), __STAT__)
    allocate(nir8(1,LM), __STAT__)
    allocate(qrr8(1,LM), __STAT__)
    allocate(qsr8(1,LM), __STAT__)
    allocate(nrr8(1,LM), __STAT__)
    allocate(nsr8(1,LM), __STAT__)
    allocate(qgr8(1,LM), __STAT__)
    allocate(ngr8(1,LM), __STAT__)
    allocate(relvarr8(1,LM), __STAT__)
    allocate(accre_enhanr8(1,LM), __STAT__)
    allocate(plevr8(1,LM), __STAT__)
    allocate(pdelr8(1,LM), __STAT__)
    allocate(cldfr8(1,LM), __STAT__)
    allocate(liqcldfr8(1,LM), __STAT__)
    allocate(icecldfr8(1,LM), __STAT__)
    allocate(qsatfacr8(1,LM), __STAT__)
    allocate(qcsinksum_rate1ordr8(1,LM), __STAT__)
    allocate(naair8(1,LM), __STAT__)
    allocate(npccninr8(1,LM), __STAT__)
    allocate(tlatr8(1,LM), __STAT__)
    allocate(qvlatr8(1,LM), __STAT__)
    allocate(qctendr8(1,LM), __STAT__)
    allocate(qitendr8(1,LM), __STAT__)
    allocate(nctendr8(1,LM), __STAT__)
    allocate(nitendr8(1,LM), __STAT__)
    allocate(qrtendr8(1,LM), __STAT__)
    allocate(qstendr8(1,LM), __STAT__)
    allocate(qgtendr8(1,LM), __STAT__)
    allocate(nrtendr8(1,LM), __STAT__)
    allocate(nstendr8(1,LM), __STAT__)
    allocate(ngtendr8(1,LM), __STAT__)
    allocate(effcr8(1,LM), __STAT__)
    allocate(effc_fnr8(1,LM), __STAT__)
    allocate(effir8(1,LM), __STAT__)
    allocate(sadicer8(1,LM), __STAT__)
    allocate(sadsnowr8(1,LM), __STAT__)
    allocate(nevaprr8(1,LM), __STAT__)
    allocate(evapsnowr8(1,LM), __STAT__)
    allocate(am_evp_str8(1,LM), __STAT__)
    allocate(prainr8(1,LM), __STAT__)
    allocate(prodsnowr8(1,LM), __STAT__)
    allocate(cmeoutr8(1,LM), __STAT__)
    allocate(deffir8(1,LM), __STAT__)
    allocate(pgamradr8(1,LM), __STAT__)
    allocate(lamcradr8(1,LM), __STAT__)
    allocate(qsoutr8(1,LM), __STAT__)
    allocate(dsoutr8(1,LM), __STAT__)
    allocate(qgoutr8(1,LM), __STAT__)
    allocate(ngoutr8(1,LM), __STAT__)
    allocate(dgoutr8(1,LM), __STAT__)
    allocate(qroutr8(1,LM), __STAT__)
    allocate(reff_rainr8(1,LM), __STAT__)
    allocate(reff_snowr8(1,LM), __STAT__)
    allocate(reff_graur8(1,LM), __STAT__)
    allocate(qcsevapr8(1,LM), __STAT__)
    allocate(qisevapr8(1,LM), __STAT__)
    allocate(qvresr8(1,LM), __STAT__)
    allocate(cmeioutr8(1,LM), __STAT__)
    allocate(vtrmcr8(1,LM), __STAT__)
    allocate(vtrmir8(1,LM), __STAT__)
    allocate(umrr8(1,LM), __STAT__)
    allocate(umsr8(1,LM), __STAT__)
    allocate(umgr8(1,LM), __STAT__)
    allocate(qgsedtendr8(1,LM), __STAT__)
    allocate(qcsedtenr8(1,LM), __STAT__)
    allocate(qisedtenr8(1,LM), __STAT__)
    allocate(qrsedtenr8(1,LM), __STAT__)
    allocate(qssedtenr8(1,LM), __STAT__)
    allocate(praor8(1,LM), __STAT__)
    allocate(prcor8(1,LM), __STAT__)
    allocate(mnucccor8(1,LM), __STAT__)
    allocate(mnucctor8(1,LM), __STAT__)
    allocate(msacwior8(1,LM), __STAT__)
    allocate(psacwsor8(1,LM), __STAT__)
    allocate(bergsor8(1,LM), __STAT__)
    allocate(bergor8(1,LM), __STAT__)
    allocate(meltor8(1,LM), __STAT__)
    allocate(homoor8(1,LM), __STAT__)
    allocate(qcresor8(1,LM), __STAT__)
    allocate(prcior8(1,LM), __STAT__)
    allocate(praior8(1,LM), __STAT__)
    allocate(qirestotr8(1,LM), __STAT__)
    allocate(mnuccrtotr8(1,LM), __STAT__)
    allocate(mnuccritotr8(1,LM), __STAT__)
    allocate(pracstotr8(1,LM), __STAT__)
    allocate(meltsdtr8(1,LM), __STAT__)
    allocate(frzrdtr8(1,LM), __STAT__)
    allocate(mnuccdor8(1,LM), __STAT__)
    allocate(pracgtotr8(1,LM), __STAT__)
    allocate(psacwgtotr8(1,LM), __STAT__)
    allocate(pgsacwtotr8(1,LM), __STAT__)
    allocate(pgracstotr8(1,LM), __STAT__)
    allocate(prdgtotr8(1,LM), __STAT__)
    allocate(qmultgtotr8(1,LM), __STAT__)
    allocate(qmultrgtotr8(1,LM), __STAT__)
    allocate(psacrtotr8(1,LM), __STAT__)
    allocate(npracgtotr8(1,LM), __STAT__)
    allocate(nscngtotr8(1,LM), __STAT__)
    allocate(ngracstotr8(1,LM), __STAT__)
    allocate(nmultgtotr8(1,LM), __STAT__)
    allocate(nmultrgtotr8(1,LM), __STAT__)
    allocate(npsacwgtotr8(1,LM), __STAT__)
    allocate(nroutr8(1,LM), __STAT__)
    allocate(nsoutr8(1,LM), __STAT__)
    allocate(reflr8(1,LM), __STAT__)
    allocate(areflr8(1,LM), __STAT__)
    allocate(areflzr8(1,LM), __STAT__)
    allocate(freflr8(1,LM), __STAT__)
    allocate(csrflr8(1,LM), __STAT__)
    allocate(acsrflr8(1,LM), __STAT__)
    allocate(fcsrflr8(1,LM), __STAT__)
    allocate(rercldr8(1,LM), __STAT__)
    allocate(ncair8(1,LM), __STAT__)
    allocate(ncalr8(1,LM), __STAT__)
    allocate(qrout2r8(1,LM), __STAT__)
    allocate(qsout2r8(1,LM), __STAT__)
    allocate(nrout2r8(1,LM), __STAT__)
    allocate(nsout2r8(1,LM), __STAT__)
    allocate(drout2r8(1,LM), __STAT__)
    allocate(dsout2r8(1,LM), __STAT__)
    allocate(qgout2r8(1,LM), __STAT__)
    allocate(ngout2r8(1,LM), __STAT__)
    allocate(dgout2r8(1,LM), __STAT__)
    allocate(freqgr8(1,LM), __STAT__)
    allocate(freqsr8(1,LM), __STAT__)
    allocate(freqrr8(1,LM), __STAT__)
    allocate(nficer8(1,LM), __STAT__)
    allocate(qcratr8(1,LM), __STAT__)
    allocate(tnd_qsnow(1,LM), __STAT__)
    allocate(tnd_nsnow(1,LM), __STAT__)
    allocate(re_ice(1,LM), __STAT__)
    allocate(prer_evap(1,LM), __STAT__)
    allocate(frzimmr8(1,LM), __STAT__)
    allocate(frzcntr8(1,LM), __STAT__)
    allocate(frzdepr8(1,LM), __STAT__)
    allocate(nsootr8(1,LM), __STAT__)
    allocate(rnsootr8(1,LM), __STAT__)
    allocate(npccnor8(1,LM), __STAT__)
    allocate(npsacwsor8(1,LM), __STAT__)
    allocate(npraor8(1,LM), __STAT__)
    allocate(nsubcor8(1,LM), __STAT__)
    allocate(nprc1or8(1,LM), __STAT__)
    allocate(npraior8(1,LM), __STAT__)
    allocate(nnucctor8(1,LM), __STAT__)
    allocate(nnucccor8(1,LM), __STAT__)
    allocate(nnuccdor8(1,LM), __STAT__)
    allocate(nsubior8(1,LM), __STAT__)
    allocate(nprcior8(1,LM), __STAT__)
    allocate(nsacwior8(1,LM), __STAT__)
    allocate(mnuccror8(1,LM), __STAT__)
    allocate(pracsor8(1,LM), __STAT__)
    allocate(qiresor8(1,LM), __STAT__)
    allocate(rate1ord_cw2pr(1,LM), __STAT__)
    allocate(accre_enhan_icer8(1,LM), __STAT__)
    allocate(pintr8(1,LM+1), __STAT__)
    allocate(kkvhr8(1,LM+1), __STAT__)
    allocate(rflxr8(1,LM+1), __STAT__)
    allocate(sflxr8(1,LM+1), __STAT__)
    allocate(lflxr8(1,LM+1), __STAT__)
    allocate(iflxr8(1,LM+1), __STAT__)
    allocate(gflxr8(1,LM+1), __STAT__)
    allocate(rndstr8(1,LM,10), __STAT__)
    allocate(naconr8(1,LM,10), __STAT__)   
    allocate(cfaux(1,LM), __STAT__)
        
    allocate(FQA(IM,JM,LM ), __STAT__)
    allocate(GZLO(IM,JM,LM ), __STAT__)
    allocate(QCNTOT(IM,JM,LM), __STAT__)
    allocate(CFX(IM,JM,LM), __STAT__)

    allocate(QTOT(IM,JM,LM ), __STAT__)
    allocate(QL_TOT(IM,JM,LM ), __STAT__)
    allocate(QI_TOT(IM,JM,LM ), __STAT__)
   !allocate(ACIL_AN_X(IM,JM,LM ), __STAT__)
   !allocate(ACIL_LS_X(IM,JM,LM ), __STAT__)
   !allocate(ACLL_AN_X(IM,JM,LM ), __STAT__)
   !allocate(ACLL_LS_X(IM,JM,LM ), __STAT__)
    allocate(DLPDF_X(IM,JM,LM ), __STAT__)
    allocate(DIPDF_X(IM,JM,LM ), __STAT__)
    allocate(DLFIX_X(IM,JM,LM ), __STAT__)
    allocate(DIFIX_X(IM,JM,LM ), __STAT__)
   !allocate(AUT_X(IM,JM,LM ), __STAT__)
   !allocate(SDM_X(IM,JM,LM ), __STAT__)
   !allocate(FRZ_TT_X(IM,JM,LM ), __STAT__)
   !allocate(FRZ_PP_X(IM,JM,LM ), __STAT__)
    allocate(CLDREFFI_TOP_X(IM,JM ), __STAT__)
    allocate(CLDREFFL_TOP_X(IM,JM ), __STAT__)
    allocate(NCPL_TOP_X(IM,JM ), __STAT__)
    allocate(NCPI_TOP_X(IM,JM ), __STAT__)
    allocate(NCPL_CLDBASEX(IM,JM ), __STAT__)

    GZLO     = MAPL_GRAV*ZL0

    PFL_LS = 0.0
    PFL_AN = 0.0
    PFI_LS = 0.0
    PFI_AN = 0.0
    FQA    = 0.0 
    QCNTOT = QLCN+QICN
    QL_TOT = QLCN+QLLS 
    QI_TOT = QICN+QILS 
    QTOT   = QL_TOT+QI_TOT

    where (QTOT .gt. 0.0)
    FQA= min(max(QCNTOT/QTOT, 0.0), 1.0)    
    end where

    CFLIQ=0.0
    CFICE=0.0
    RAD_CF=min(CLLS+CLCN, 1.0)
    WHERE (QTOT .gt. 0.0) 
    CFLIQ=RAD_CF*QL_TOT/QTOT
    CFICE=RAD_CF*QI_TOT/QTOT
    END WHERE
    
    rhdfdar8   = 1.e-8
    rhu00r8    = 0.95
    ttendr8=0.
    qtendr8=0.
    cwtendr8=0.
    naair8=0.
    rndstr8 = 2.0e-7
    npccninr8 = 0.
    naconr8   = 0.

    !initialize MG variables
     cldfr8 = 0.0
     prectr8 = 0.0
     precir8 = 0.0
     qctendr8 = 0.0
     qitendr8 = 0.0
     qvlatr8 = 0.0
     tlatr8 = 0.0
     nctendr8 = 0.0
     nitendr8 = 0.0
     effcr8 = 0.0
     effir8 = 0.0
     drout2r8 =0.0
     dsout2r8 = 0.0
     dgout2r8 = 0.0
     qrout2r8 = 0.0
     qsout2r8 =0.0
     qgout2r8 =0.0
     nrout2r8 = 0.0
     nsout2r8 =0.0
     ngout2r8 =0.0
     evapsnowr8 =0.0
     nevaprr8 =0.0
     cmeioutr8 =0.0
     bergsor8 =0.0
     mnucccor8 =0.0
     mnucctor8 =0.0
     homoor8 = 0.0
     mnuccror8 = 0.0
     pracsor8 = 0.0
     meltor8 =0.0
     qisedtenr8 =0.0
     bergor8 =0.0
     psacwsor8 = 0.0
     qcresor8 =0.0
     qiresor8 = 0.0
     praor8 =0.0
     prcor8 = 0.0
     prcior8 =0.0
     praior8 = 0.0
     msacwior8 =0.0
     frzrdtr8 =0.0
     meltsdtr8 = 0.0
     nnucctor8 =0.0
     nnucccor8 = 0.0
     nnuccdor8 =0.0
     nsacwior8 =0.0
     nsubior8 = 0.0
     npraior8 =0.0
     nprcior8 =0.0
     npccnor8 = 0.0
     npsacwsor8 =0.0
     npraor8 =0.0
     nsubcor8 =0.0
     nprc1or8 =0.0
     rndstr8 = 2.0e-7
     naconr8   = 0.
     lflxr8 = 0.0
     iflxr8 = 0.0
     rflxr8 = 0.0
     sflxr8 = 0.0
     gflxr8 = 0.0
     frzcntr8 =  0.0
     qrtendr8 =  0.0
     nrtendr8 =  0.0
     qstendr8 =  0.0
     nstendr8 =  0.0
     qgtendr8 =  0.0
     ngtendr8 =  0.0

    !Tuning factors
    accre_enhanr8= ACC_ENH
    accre_enhan_icer8= ACC_ENH_ICE
    autscx = 1.0
    disp_liu = LIU_MU
    ui_scale = UISCALE
    ur_scale  = URSCALE
    ts_autice = DT_MOIST
    mtimesc  = DT_MOIST     
    if (TS_AUTO_ICE .gt. 0.) ts_autice= TS_AUTO_ICE    
    if (MTIME .gt. 0.0) mtimesc=MTIME    
    if (QCVAR_CST .gt. 0.) then 
       QCVAR =  QCVAR_CST
    else
   !!  xscale = (9000.0/real(imsize))**(-0.666)
       TMP2D = (2.517514*SQRT(AREA)/1.e3)**(-0.666)
       call estimate_qcvar(QCVAR, IM, JM, LM, PLmb, T, GZLO, Q, QST3, TMP2D)
    end if    
  
    do J=1,JM
       do I=1,IM       
               kbmin =1                          
               rndstr8 = 2.0e-7
               naconr8   = 0.
               cldfr8(1,1:LM)     = RAD_CF(I,J,1:LM) !Assume minimum overlap 
               liqcldfr8(1,1:LM)  = CFLIQ(I,J,1:LM) 
               icecldfr8(1,1:LM)  = CFICE(I,J,1:LM)                  
               cldor8             = cldfr8  
               ter8(1,1:LM)       = T(I,J,1:LM)
               qvr8(1,1:LM)       = Q(I,J,1:LM)
               qcr8(1,1:LM)       = QL_TOT(I,J,1:LM)
               qir8(1,1:LM)       = QI_TOT(I,J,1:LM)
               ncr8(1,1:LM)       = MAX(  NCPL(I,J,1:LM), 0.0) 
               nir8(1,1:LM)       = MAX(  NCPI(I,J,1:LM), 0.0) 

               ! Nucleation  tendencies 
               naair8(1,1:LM)     = max(( INC_NUC(I, J, 1:LM)*cldfr8(1,1:LM) - nir8(1,1:LM))/DT_MOIST, 0.0) 
               npccninr8(1,1:LM)  = max((CDNC_NUC(I, J, 1:LM)*cldfr8(1,1:LM) - ncr8(1,1:LM))/DT_MOIST, 0.0)
               
               where  ((naair8 .gt. 1.0e3)) ! add cloud fraction if nucleation is happening 2018
                   icecldfr8 = max(0.05,  icecldfr8)
               end where 

               where (cldfr8(1,:) .ge. 0.001) 
                  frzimmr8 (1,1:LM) = MIN(DNHET_IMM(I, J, 1:LM), ncr8(1,1:LM)/cldfr8(1,1:LM)/DT_MOIST) !tendency    
               elsewhere 
                  frzimmr8 (1,1:LM) = 0.0 
               end where
               
               nbincontactdust = 1

               DO K=kbmin, LM
                  call getINsubset(1,  AeroProps(I, J, K), AeroAux) 
                  nbincontactdust = AeroAux_b%nmods 
                  naconr8(1, K, 1:nbincontactdust) =  AeroAux_b%num(1:nbincontactdust)
                  rndstr8( 1, K, 1:nbincontactdust)=AeroAux_b%dpg(1:nbincontactdust)/2.0
                  ! Get black carbon properties for contact ice nucleation      
                  call getINsubset(2, AeroProps(I, J, K), AeroAux)          
                  nsootr8  (1, K)  = sum(AeroAux_b%num) !
                  naux = AeroAux_b%nmods 
                  rnsootr8 (1, K)  = sum(AeroAux_b%dpg(1:naux))/naux
               END DO

               pdelr8(1,1:LM)  = PLE(I,J,1:LM) - PLE(I,J,0:LM-1)  
               rpdelr8      = 1./pdelr8 
               pintr8(1,1:LM+1) = PLE(I,J,0:LM)  
               plevr8(1,1:LM)      = 100.*PLmb(I,J,1:LM)
               zmr8(1,1:LM)        = ZL0(I,J,1:LM)     
               kkvhr8(1,1:LM+1) = KH(I,J,0:LM)  
               ficer8 = qir8 /( qcr8+qir8 + 1.e-10 )  
                          
               if (AUT_SCALE .gt. 0.0) then 
                  autscx = AUT_SCALE
               else
                  autscx =  min(max(0., (300.0 - T(I,J,LM))/ABS(AUT_SCALE)), 1.0)
                  autscx  =  1.0 - 0.995*autscx
               end if
               
               relvarr8 = QCVAR(I, J)
                  
               ! for MG23 (initial values)     
               frzdepr8(1,1:LM)  = NHET_DEP(I, J, 1:LM)/DT_MOIST
               qrr8(1,1:LM)      = QRAIN(I, J,1:LM)
               qsr8(1,1:LM)      = QSNOW(I, J,1:LM)
               qgr8(1,1:LM)      = QGRAUPEL(I, J,1:LM)                        
               nrr8(1,1:LM)      = NRAIN(I, J,1:LM)
               nsr8(1,1:LM)      = NSNOW(I, J,1:LM)
               ngr8(1,1:LM)      = NGRAUPEL(I, J,1:LM)                         
               qsatfacr8         = 1.0                        
               SCICE_tmp(1,1:LM) = SC_ICE(I, J, 1:LM)
               FQA_tmp(1,1:LM)   = FQA(I, J, 1:LM) 
                  
                   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  !CALLS to MG versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               
   if (MGVERSION .lt. 2)  then          
               
               call set_qcvar (QCVAR(I, J))
               
               call mmicro_pcond (                           &
                    ncolmicro, ncolmicro, dt_r8, DT_MICRO, ter8, ttendr8,                   &
                    ncolmicro, LM ,                                      &
                    qvr8, qtendr8, cwtendr8, qcr8, qir8,          &
                    ncr8, nir8, plevr8, pdelr8, cldfr8,           &
                    liqcldfr8, icecldfr8,                         &
                    cldor8, pintr8, rpdelr8, zmr8, omegr8,        &
                    rate1ord_cw2pr,                               &  ! <= actually an output
                    naair8, npccninr8, rndstr8,naconr8,           &
                    rhdfdar8, rhu00r8, ficer8,                    & 
                    !                          ** outputs **
                    tlatr8, qvlatr8,        &
                    qctendr8, qitendr8, nctendr8, nitendr8, effcr8,    &
                    effc_fnr8, effir8, prectr8, precir8,             &  
                    nevaprr8, evapsnowr8,      &
                    prainr8, prodsnowr8,       &
                    cmeoutr8, deffir8, pgamradr8, &
                    lamcradr8,qsout2r8,dsout2r8, qrout2r8,drout2r8, &
                    qcsevapr8,qisevapr8,   &
                    qvresr8,cmeioutr8, &
                    vtrmcr8,vtrmir8,   &
                    qcsedtenr8,qisedtenr8, &
                    praor8,prcor8,mnucccor8, &
                    mnucctor8,msacwior8,psacwsor8,&
                    bergsor8,bergor8,meltor8, &
                    homoor8,qcresor8,prcior8, &
                    praior8,qiresor8,  &
                    mnuccror8,pracsor8, &
                    meltsdtr8,frzrdtr8, ncalr8, ncair8, mnuccdor8, nnucctor8, &
                    nsoutr8, nroutr8, frzimmr8, disp_liu, &
                    nsootr8, rnsootr8, ui_scale, autscx, mtimesc, &
                    nnuccdor8, nnucccor8, nsacwior8, nsubior8, nprcior8, &
                    npraior8, npccnor8, npsacwsor8, nsubcor8, npraor8, nprc1or8,  nbincontactdust, &
                    ts_autice, rflxr8, sflxr8, accre_enhanr8, accre_enhan_icer8)

    else ! MG2/3
        
         call  micro_mg_tend_interface ( DT_MICRO, INT(PDFSHAPE), 1.-RHCRIT3D(I, J, 1:LM), SCICE_tmp, FQA_tmp, &
                             ncolmicro,             LM,               dt_r8,       & 
                             CNV_FRC(I,J), SRF_TYPE(I,J), &
                             ter8,                          qvr8,                          &
                             qcr8,                          qir8,                          &
                             ncr8,                          nir8,                          &
                             qrr8,                          qsr8,                          &
                             nrr8,                          nsr8,                          &
                             qgr8,                          ngr8,                          &
                             relvarr8,                      accre_enhanr8,   accre_enhan_icer8,  &
                             plevr8,                        pdelr8,                        &
                             cldfr8, liqcldfr8, icecldfr8,  qsatfacr8,                     &
                             qcsinksum_rate1ordr8,                                         &
                             naair8,                         npccninr8,                      &
                             rndstr8,                        naconr8,                        &
                             tlatr8,                         qvlatr8,                        &
                             qctendr8,                       qitendr8,                       &
                             nctendr8,                       nitendr8,                       &
                             qrtendr8,                       qstendr8,                       &
                             nrtendr8,                       nstendr8,                       & 
                             qgtendr8,                       ngtendr8,                       &
                             effcr8,               effc_fnr8,            effir8,              &
                             sadicer8,                       sadsnowr8,                      &
                             prectr8,                        precir8,                        &
                             nevaprr8,                       evapsnowr8,                     &
                             am_evp_str8,                                                  &
                             prainr8,                        prodsnowr8,                     &
                             cmeoutr8,                       deffir8,                        &
                             pgamradr8,                      lamcradr8,                      &
                             qsoutr8,                        dsoutr8,                        &
                             qgoutr8,     ngoutr8,           dgoutr8,                        &
                             lflxr8,               iflxr8,   & 
                             gflxr8,                           &
                             rflxr8,           sflxr8,    qroutr8,          &
                             reff_rainr8,      reff_snowr8, reff_graur8,        &
                             qcsevapr8,            qisevapr8,            qvresr8,              &
                             cmeioutr8,            vtrmcr8,              vtrmir8,              &
                             umrr8,                          umsr8,                          &
                             umgr8,                          qgsedtendr8,                    &    
                             qcsedtenr8,                     qisedtenr8,                     &
                             qrsedtenr8,                     qssedtenr8,                     &
                             praor8,                       prcor8,                       &
                             mnucccor8,          mnucctor8,          msacwior8,          &
                             psacwsor8,          bergsor8,           bergor8,            &
                             meltor8,                      homoor8,                      &
                             qcresor8,           prcior8,            praior8,            &
                             qirestotr8,           mnuccrtotr8,          mnuccritotr8, pracstotr8,           &                           
                             meltsdtr8,         frzrdtr8,          mnuccdor8,          &
                             pracgtotr8,           psacwgtotr8,          pgsacwtotr8,          &
                             pgracstotr8,          prdgtotr8,           &
                             qmultgtotr8,          qmultrgtotr8,         psacrtotr8,           &
                             npracgtotr8,          nscngtotr8,           ngracstotr8,          &
                             nmultgtotr8,          nmultrgtotr8,         npsacwgtotr8,         & 
                             nroutr8,                            nsoutr8,                        &
                             reflr8,               areflr8,              areflzr8,             &
                             freflr8,              csrflr8,              acsrflr8,             &
                             fcsrflr8,                       rercldr8,                       &
                             ncair8,                         ncalr8,                         &
                             qrout2r8,                       qsout2r8,                       &
                             nrout2r8,                       nsout2r8,                       &
                             drout2r8,                       dsout2r8,                       &
                             qgout2r8,     ngout2r8,         dgout2r8,   freqgr8,                     &
                             freqsr8,                        freqrr8,                        &
                             nficer8,                        qcratr8,                        &
!                             errstring, & ! Below arguments are "optional" (pass null pointers to omit).
                      !       tnd_qsnow,          tnd_nsnow,          re_ice,    &
                             prer_evap, &
                             frzimmr8,             frzcntr8,              frzdepr8,  & ! contact is not passed since it depends on the droplet size dist
                             nsootr8, rnsootr8,  & ! soot for contact IN
                             npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, &  ! Number tendencies for liquid
                             npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, nsacwior8,  &  ! Number tendencies for ice
                             ts_autice, ui_scale, autscx , disp_liu, nbincontactdust, ur_scale)


    end if 

             ! Update state after microphysics
              IF (MGVERSION > 1) then 
                   QRAIN(I,J,1:LM)    = max(QRAIN(I,J,1:LM) + qrtendr8(1, 1:LM)*DT_R8, 0.0) ! grid average 
                   QSNOW(I,J,1:LM)    = max(QSNOW(I,J,1:LM) + qstendr8(1, 1:LM)*DT_R8, 0.0) ! grid average                     
                   QGRAUPEL(I,J,1:LM) = max(QGRAUPEL(I,J,1:LM) + qgtendr8(1, 1:LM)*DT_R8, 0.0) ! grid average 
                   NRAIN(I,J,1:LM)    = max(NRAIN(I,J,1:LM) + nrtendr8(1, 1:LM)*DT_R8, 0.0)
                   NSNOW(I,J,1:LM)    = max(NSNOW(I,J,1:LM) + nstendr8(1, 1:LM)*DT_R8, 0.0)
                   NGRAUPEL(I,J,1:LM) = max(NGRAUPEL(I,J,1:LM) + ngtendr8(1, 1:LM)*DT_R8, 0.0)
                   if (associated(CLDREFFR)) CLDREFFR(I,J,1:LM) = reff_rainr8(1,1:LM) 
                   if (associated(CLDREFFS)) CLDREFFS(I,J,1:LM) = reff_snowr8(1,1:LM) 
                   if (associated(CLDREFFG)) CLDREFFG(I,J,1:LM) = reff_graur8(1,1:LM) 
              else                    
                   QRAIN(I,J,1:LM)    = max(qrout2r8(1,1:LM), 0.d0) ! grid average 
                   QSNOW(I,J,1:LM)    = max(qsout2r8(1,1:LM), 0.d0)                      
                   QGRAUPEL(I,J,1:LM) = 0.0 ! grid average                    
                   NRAIN(I,J,1:LM)    = max(nrout2r8(1,1:LM), 0.d0)
                   NSNOW(I,J,1:LM)    = max(nsout2r8(1,1:LM), 0.d0)
                   NGRAUPEL(I,J,1:LM) = 0.0 ! grid average                 
                   if (associated(CLDREFFR)) CLDREFFR(I,J,1:LM) = drout2r8(1,1:LM)/2.d0       
                   if (associated(CLDREFFS)) CLDREFFS(I,J,1:LM) = dsout2r8(1,1:LM)/2.d0
               end if          
               QL_TOT(I,J,1:LM) = max(QL_TOT(I,J,1:LM)   + qctendr8(1,1:LM)*DT_R8, 0.0)
               QI_TOT(I,J,1:LM) = max(QI_TOT(I,J,1:LM)   + qitendr8(1,1:LM)*DT_R8, 0.0)
               Q(I,J,1:LM)      = max(Q(I,J,1:LM)        +  qvlatr8(1,1:LM)*DT_R8, 0.0)
               T(I,J,1:LM)      =     T(I,J,1:LM)        +   tlatr8(1,1:LM)*DT_R8/MAPL_CP
               NCPL(I,J,1:LM)   = max(NCPL(I,J,1:LM)     + nctendr8(1,1:LM)*DT_R8, 0.0)
               NCPI(I,J,1:LM)   = max(NCPI(I,J,1:LM)     + nitendr8(1,1:LM)*DT_R8, 0.0)
               CLDREFFL(I,J,1:LM) = max(effcr8(1,1:LM)*1.0e-6, 1.0e-6)
               CLDREFFI(I,J,1:LM) = max(effir8(1,1:LM)*1.0e-6, 1.0e-6)
 
              ! precipitation 
               LS_PRCP(I,J)     = max(1000.*(prectr8(1)-precir8(1)), 0.0)
               LS_SNR(I,J)      = max(1000.* precir8(1)            , 0.0)          
              ! precip fluxes
               PFL_LS(I, J, 1:LM) = rflxr8(1,2:LM+1) !+ lflxr8(1,1:LM)
               PFI_LS(I, J, 1:LM) = sflxr8(1,2:LM+1) + gflxr8(1,2:LM+1) !+  iflxr8(1,1:LM)
               if (MGVERSION < 2) then   
                  !normalize precip flux
                  if (PFL_LS(I, J, LM) .gt. 1.0e-7) PFL_LS(I, J, 1:LM) = PFL_LS(I, J, 1:LM)*LS_PRCP(I,J)/PFL_LS(I, J, LM)                    
                  if (PFI_LS(I, J, LM) .gt. 1.0e-7) PFI_LS(I, J, 1:LM) = PFI_LS(I, J, 1:LM)*LS_SNR(I,J)/PFI_LS(I, J, LM)                    
               end if 
               
               ! diagnostics from the microphysics********************

               if (associated(BERG))        BERG (I,J,1:LM)    = REAL(bergor8(1,1:LM)) 
               if (associated(BERGS))       BERGS(I,J,1:LM)    = REAL(bergsor8(1,1:LM))

               if (associated(RSU_LS))      RSU_LS(I,J,1:LM)   = REAL(evapsnowr8(1,1:LM))                    ! Snow evap   
               if (associated(REV_LS))      REV_LS(I,J,1:LM)   = REAL(nevaprr8(1,1:LM) )                     ! rain evap
               if (associated(SUBLC))       SUBLC(I,J,1:LM)    = REAL(cmeioutr8(1,1:LM)) + SUBLC(I,J,1:LM)   ! Ice subl already grid -ave
               if (associated(EVAPC))       EVAPC(I,J,1:LM)    = REAL(qcsevapr8(1,1:LM) ) + EVAPC(I,J,1:LM)  ! cloud evap
              !if (associated(FRZ_TT_X))    FRZ_TT_X(I,J,1:LM) = REAL(mnucccor8(1,1:LM)+ mnucctor8(1,1:LM) + homoor8(1,1:LM))!ice mixing ratio from nucleation (hom+het)
              !if (associated(FRZ_PP_X))    FRZ_PP_X(I,J,1:LM) = REAL(mnuccror8(1,1:LM) + pracsor8(1,1:LM)) !freezing of rain (hom+ het freezing and accretion by snow) 
               if (associated(MELT))        MELT(I,J,1:LM)     = REAL(meltor8(1,1:LM))         ! melting of cloud ice  and snow
              !if (associated(SDM_X))       SDM_X(I,J,1:LM)    = REAL(qisedtenr8(1,1:LM))      ! ice sed   
               if (associated(QCRES))       QCRES(I,J,1:LM)    = REAL(qcresor8(1,1:LM) )       ! residual liquid condensation
               if (associated(QIRES))       QIRES(I,J,1:LM)    = REAL(qiresor8(1,1:LM))        ! residual ice condensation   

              !if (associated(AUT_X))       AUT_X(I,J,1:LM)     = REAL(prcor8(1,1:LM))         ! Aut liquid 
               if (associated(AUTICE))      AUTICE(I,J,1:LM)    = REAL(prcior8(1,1:LM))        ! Aut ice    
              !if (associated(ACIL_LS_X))   ACIL_LS_X(I,J,1:LM) = REAL(psacwsor8(1,1:LM) + msacwior8(1,1:LM)) ! Acc + collection of cloud  by snow
              !if (associated(ACLL_LS_X))   ACLL_LS_X(I,J,1:LM) = REAL(praor8(1,1:LM) )        ! Acc cloud by rain    
              !if (associated(ACIL_AN_X))   ACIL_AN_X(I,J,1:LM) = REAL(praior8(1,1:LM))        ! Acc ice  by snow
              !if (associated(ACLL_AN_X))   ACLL_AN_X(I,J,1:LM) = REAL(msacwior8(1,1:LM))      ! HM process

               if (associated(FRZPP_LS))    FRZPP_LS(I,J,1:LM)    = REAL(frzrdtr8(1,1:LM)) / MAPL_CP !precip freezing latent heat rate
               if (associated(SNOWMELT_LS)) SNOWMELT_LS(I,J,1:LM) = REAL(meltsdtr8(1,1:LM))/ MAPL_CP !melting of snow latent heat rate

               !diagnostics for number concentration budget (all grid-average, Kg-1 s-1)
              
               if (associated(DNCACRIS))   DNCACRIS (I,J,1:LM)    = REAL(npraior8(1,1:LM))  !ice number tendency from accretion by snow
               if (associated(DNHET_CT))   DNHET_CT(I,J,1:LM)     = REAL(nnucctor8(1,1:LM)) !ice number tendency from contact IN  
               if (associated(DNHET_IMM))  DNHET_IMM(I,J,1:LM)    = REAL(nnucccor8(1,1:LM)) !ice number tendency from immersion IN    
               if (associated(DNCNUC))     DNCNUC(I,J,1:LM)       = REAL(nnuccdor8(1,1:LM)) !ice number tendency from nucleation on aerosol             
               if (associated(DNCSUBL))    DNCSUBL (I,J,1:LM)     = REAL(nsubior8(1,1:LM))  !ice number tendency from sublimation               
               if (associated(DNCAUTICE))  DNCAUTICE (I,J,1:LM)   = REAL(nprcior8(1,1:LM))  !ice number tendency from autoconversion
               if (associated(DNCHMSPLIT)) DNCHMSPLIT(I,J,1:LM)   = REAL(nsacwior8(1,1:LM)) !ice number tendency from H-M process

               if (associated(DNDCCN))    DNDCCN(I,J,1:LM)        = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               if (associated(DNDACRLS))  DNDACRLS(I,J,1:LM)      = REAL(npsacwsor8(1,1:LM) )!droplet number tendency from accretion by snow
               if (associated(DNDACRLR))  DNDACRLR(I,J,1:LM)      = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               if (associated(DNDEVAPC))  DNDEVAPC(I,J,1:LM)      = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation 
               if (associated(DNDAUTLIQ)) DNDAUTLIQ(I,J,1:LM)     = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion

            enddo !I
         enddo !J
         !============================================Finish 2-moment micro implementation===========================

         !update water tracers
         QLCN=QL_TOT*FQA
         QLLS=QL_TOT-QLCN
         QICN=QI_TOT*FQA
         QILS=QI_TOT-QICN

         !============ Put cloud fraction back in contact with the PDF and create new condensate if neccesary (Barahona et al., GMD, 2014)============

          do K= 1, LM
             do J=1,JM
                do I= 1, IM
                  
                  call update_cld( &
                         DT_MOIST                , &
                         1.- RHCRIT3D(I, J, K)   , &
                         PDFSHAPE                , &
                         CNV_FRC(I, J)           , &
                         SRF_TYPE(I, J)          , &
                         PLmb(I, J, K)           , &
                         Q (I, J, K)             , &
                         QLLS(I, J, K)           , &
                         QLCN(I, J, K)           , &
                         QILS(I, J, K)           , &
                         QICN(I, J, K)           , &
                         T(I, J, K)              , &
                         CLLS(I, J, K)           , &
                         CLCN(I, J, K)           , &
                         SC_ICE(I, J, K)         , &
                         NCPI(I, J, K)           , &
                         NCPL(I, J, K)           , &
                         RHCmicro(I, J, K))

               end do 
           end do
        end do 

        ! Make sure ice and liquid stay within T limits  
        call meltfrz_inst2M ( &
              IM,JM,LM    , &
              T           , &
              QLLS        , &
              QLCN        , &
              QILS        , &
              QICN        , &               
              NCPL        , &
              NCPI          )
                          
         QTOT   = QICN+QILS+QLCN+QLLS
         QL_TOT = QLCN+QLLS
         QI_TOT = QICN+QILS
         RAD_CF = min(CLLS+CLCN, 1.0)
         !=============================================End Stratiform cloud processes==========================================
         !======================================================================================================================

         !Calculate CFICE and CFLIQ 
         CFLIQ=0.0
         CFICE=0.0
         WHERE (QTOT .gt. 1.0e-12) 
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE
         CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         CFICE=MAX(MIN(CFICE, 1.0), 0.0)

         !======================================================================================================================
         !===========================Clean stuff and send it to radiation ======================================================
         !======================================================================================================================

         where (QI_TOT .le. 0.0)
            CFICE =0.0
            NCPI=0.0
            CLDREFFI =  MAPL_UNDEF
         end where

         where (QL_TOT .le. 0.0)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL =  MAPL_UNDEF
         end where

      ! use RAD_ variables for holding water species to process below
         RAD_QV = Q
         RAD_QL = QLLS+QLCN
         RAD_QI = QILS+QICN
         RAD_QR = QRAIN
         RAD_QS = QSNOW
         RAD_QG = QGRAUPEL
      !These will become in-cloud for radiation later============== 

      ! Fill GEOS precip diagnostics
         PRCP_RAIN    = LS_PRCP
         PRCP_SNOW    = LS_SNR
         PRCP_ICE     = 0.0
         PRCP_GRAUPEL = 0.0
         ICE          = PRCP_ICE + PRCP_GRAUPEL
         FRZR         = 0.0
     ! Redistribute precipitation fluxes for chemistry
         TMP3D =  MIN(1.0,MAX(QLCN/MAX(RAD_QL,1.E-8),0.0))
         PFL_AN(:,:,1:LM) = PFL_LS(:,:,1:LM) * TMP3D
         PFL_LS(:,:,1:LM) = PFL_LS(:,:,1:LM) - PFL_AN(:,:,1:LM)
         TMP3D =  MIN(1.0,MAX(QICN/MAX(RAD_QI,1.E-8),0.0))
         PFI_AN(:,:,1:LM) = PFI_LS(:,:,1:LM) * TMP3D
         PFI_LS(:,:,1:LM) = PFI_LS(:,:,1:LM) - PFI_AN(:,:,1:LM)
     ! cleanup suspended precipitation condensates
         call FIX_NEGATIVE_PRECIP(RAD_QR, RAD_QS, RAD_QG)
     ! Fill vapor/rain/snow/graupel state
         Q        = RAD_QV
         QRAIN    = RAD_QR
         QSNOW    = RAD_QS
         QGRAUPEL = RAD_QG
     ! Radiation Coupling
         do L = 1, LM
           do J = 1, JM
             do I = 1, IM
              ! cleanup clouds
               call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
              ! get radiative properties
               call RADCOUPLE ( T(I,J,L), PLmb(I,J,L), CLLS(I,J,L), CLCN(I,J,L), &
                     Q(I,J,L), QLLS(I,J,L), QILS(I,J,L), QLCN(I,J,L), QICN(I,J,L), QRAIN(I,J,L), QSNOW(I,J,L), QGRAUPEL(I,J,L), NCPL(I,J,L), NCPI(I,J,L), &
                     RAD_QV(I,J,L), RAD_QL(I,J,L), RAD_QI(I,J,L), RAD_QR(I,J,L), RAD_QS(I,J,L), RAD_QG(I,J,L), RAD_CF(I,J,L), &
                     CLDREFFL(I,J,L), CLDREFFI(I,J,L), &
                     FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
            enddo
          enddo
         enddo
         call FILLQ2ZERO(RAD_QV, MASS, TMP2D)
         call FILLQ2ZERO(RAD_QL, MASS, TMP2D)
         call FILLQ2ZERO(RAD_QI, MASS, TMP2D)
         call FILLQ2ZERO(RAD_QR, MASS, TMP2D)
         call FILLQ2ZERO(RAD_QS, MASS, TMP2D)
         call FILLQ2ZERO(RAD_QG, MASS, TMP2D)
         call FILLQ2ZERO(RAD_CF, MASS, TMP2D)
         RAD_QL = MIN( RAD_QL , 0.001 )  ! Still a ridiculously large
         RAD_QI = MIN( RAD_QI , 0.001 )  ! value.
         RAD_QR = MIN( RAD_QR , 0.01  )  ! value.
         RAD_QS = MIN( RAD_QS , 0.01  )  ! value.
         RAD_QG = MIN( RAD_QG , 0.01  )  ! value.
         where (QILS+QICN .le. 0.0)
            CLDREFFI = 36.0e-6
         end where
         where (QLLS+QLCN .le. 0.0)
            CLDREFFL = 14.0e-6
         end where

         ! Update microphysics tendencies
         DQVDT_micro = ( Q          - DQVDT_micro) / DT_MOIST
         DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
         DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
         DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
         DQRDT_micro = ( QRAIN      - DQRDT_micro) / DT_MOIST
         DQSDT_micro = ( QSNOW      - DQSDT_micro) / DT_MOIST
         DQGDT_micro = ( QGRAUPEL   - DQGDT_micro) / DT_MOIST
          DUDT_micro = ( U          -  DUDT_micro) / DT_MOIST
          DVDT_micro = ( V          -  DVDT_micro) / DT_MOIST
          DTDT_micro = ( T          -  DTDT_micro) / DT_MOIST
        call MAPL_TimerOff(MAPL,"---CLDMICRO")

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
        call MAPL_GetPointer(EXPORT, PTR3D   , 'DBZ'     , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_MAX , 'DBZ_MAX' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_1KM , 'DBZ_1KM' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_TOP , 'DBZ_TOP' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_M10C, 'DBZ_M10C', RC=STATUS); VERIFY_(STATUS)

        if (associated(PTR3D) .OR. &
            associated(DBZ_MAX) .OR. associated(DBZ_1KM) .OR. associated(DBZ_TOP) .OR. associated(DBZ_M10C)) then

            call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,DBZ_LIQUID_SKIN)
            if (associated(PTR3D)) PTR3D = TMP3D

            if (associated(DBZ_MAX)) then
               DBZ_MAX=-9999.0
               DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                  DBZ_MAX(I,J) = MAX(DBZ_MAX(I,J),TMP3D(I,J,L))
               END DO ; END DO ; END DO
            endif

            if (associated(DBZ_1KM)) then
               call cs_interpolator(1, IM, 1, JM, LM, TMP3D, 1000., ZLE0, DBZ_1KM, -20.)
            endif

            if (associated(DBZ_TOP)) then
               DBZ_TOP=MAPL_UNDEF
               DO J=1,JM ; DO I=1,IM
                  DO L=LM,1,-1
                     if (ZLE0(i,j,l) >= 25000.) continue
                     if (TMP3D(i,j,l) >= 18.5 ) then
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
                         DBZ_M10C(I,J) = TMP3D(I,J,L)
                         exit
                     endif
                  END DO
               END DO ; END DO
            endif

        endif

        call MAPL_GetPointer(EXPORT, PTR2D , 'DBZ_MAX_R' , RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) then
            call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,0.0*QSNOW,0.0*QGRAUPEL,IM,JM,LM,1,0,DBZ_LIQUID_SKIN)
             PTR2D=-9999.0
             DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
             END DO ; END DO ; END DO
        endif
        call MAPL_GetPointer(EXPORT, PTR2D , 'DBZ_MAX_S' , RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) then
            call CALCDBZ(TMP3D,100*PLmb,T,Q,0.0*QRAIN,QSNOW,0.0*QGRAUPEL,IM,JM,LM,1,0,DBZ_LIQUID_SKIN)
             PTR2D=-9999.0
             DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
             END DO ; END DO ; END DO 
        endif
        call MAPL_GetPointer(EXPORT, PTR2D , 'DBZ_MAX_G' , RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) then
            call CALCDBZ(TMP3D,100*PLmb,T,Q,0.0*QRAIN,0.0*QSNOW,QGRAUPEL,IM,JM,LM,1,0,DBZ_LIQUID_SKIN)
             PTR2D=-9999.0
             DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
             END DO ; END DO ; END DO  
        endif

        call MAPL_GetPointer(EXPORT, PTR3D, 'QRTOT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = QRAIN
        call MAPL_GetPointer(EXPORT, PTR3D, 'QSTOT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = QSNOW
        call MAPL_GetPointer(EXPORT, PTR3D, 'QGTOT', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) PTR3D = QGRAUPEL

     call MAPL_TimerOff(MAPL,"--MGB2_2M",RC=STATUS)

end subroutine MGB2_2M_Run

end module GEOS_MGB2_2M_InterfaceMod
