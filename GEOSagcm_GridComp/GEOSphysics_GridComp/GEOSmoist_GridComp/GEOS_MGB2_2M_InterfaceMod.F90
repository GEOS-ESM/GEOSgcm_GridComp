! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_MGB2_2M_InterfaceMod -- A Module to interface with the
!   MGB2_2M cloud microphysics

module GEOS_MGB2_2M_InterfaceMod

  use ESMF
  use MAPL, r8 => MAPL_R8
  use GEOS_UtilsMod
  use cldwat2m_micro
  use cldmacro
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
         character(len=ESMF_MAXSTR) :: QW
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

  public :: MGB2_2M_Setup, MGB2_2M_Initialize, MGB2_2M_Run
  public :: MGVERSION

contains

subroutine MGB2_2M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    
    IAm = "GEOS_MGB2_2M_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = trim(COMP_NAME) // Iam

    call ESMF_ConfigGetAttribute( CF, MGVERSION, Label="MGVERSION:",  default=2, RC=STATUS)

    ! !INTERNAL STATE:

      FRIENDLIES%QV       = trim(COMP_NAME)//":DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS"
      FRIENDLIES%QW       = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%CLLS     = trim(COMP_NAME)//":DYNAMICS"
      FRIENDLIES%CLCN     = trim(COMP_NAME)//":DYNAMICS"
      FRIENDLIES%QLLS     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QLCN     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QILS     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QICN     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QRAIN    = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QSNOW    = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QGRAUPEL = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%NCPI     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%NCPL     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%NRAIN    = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%NSNOW    = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%NGRAUPEL = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"

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
         SHORT_NAME ='NCPL',                                       &
         LONG_NAME  ='particle_number_for_liquid_cloud',           &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NCPL),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 50.0e6 ,                             RC=STATUS  )  
    VERIFY_(STATUS)                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NCPI',                                       &
         LONG_NAME  ='particle_number_for_ice_cloud',              &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NCPI),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 1.0e3,                               RC=STATUS  )  
    VERIFY_(STATUS)                                                   

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NRAIN',                                      &
         LONG_NAME  ='particle_number_for_rain',                   &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NRAIN),                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 0.0 ,                                RC=STATUS  )  
    VERIFY_(STATUS)                                                          

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NSNOW',                                      &
         LONG_NAME  ='particle_number_for_snow',                   &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NSNOW),                      &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 0.0,                                 RC=STATUS  )  
    VERIFY_(STATUS)                                                   

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME ='NGRAUPEL',                                   &
         LONG_NAME  ='particle_number_for_graupel',                &
         UNITS      ='kg-1',                                       &
         FRIENDLYTO = trim(FRIENDLIES%NGRAUPEL),                   &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,                        &
         DEFAULT = 0.0,                                 RC=STATUS  )  
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
         SHORT_NAME = 'QW',                                        &
         LONG_NAME  = 'mass_fraction_of_wet_air',                  &
         UNITS      = 'kg kg-1',                                   &
         RESTART    = MAPL_RestartSkip,                            &
         FRIENDLYTO = trim(FRIENDLIES%QW),                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                      

    call MAPL_TimerAdd(GC, name="--MGB2_2M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine MGB2_2M_Setup

subroutine MGB2_2M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, QW

    real DCS, QCVAR_, WBFFACTOR, NC_CST, NI_CST, NG_CST, MUI_CST
    logical  :: nccons, nicons, ngcons, do_graupel
    real(ESMF_KIND_R8)  Dcsr8, qcvarr8,  micro_mg_berg_eff_factor_in, ncnstr8, ninstr8, ngnstr8, mui_cnstr8

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QW,       'QW'      , RC=STATUS); VERIFY_(STATUS)
    QW = Q+QLLS+QLCN+QILS+QICN+QRAIN+QSNOW+QGRAUPEL

#ifdef NODISABLE

    call MAPL_GetResource(MAPL, DCS,      'DCS:'    , DEFAULT=350.0e-6, RC=STATUS )
    VERIFY_(STATUS)
    Dcsr8 = DCS
    call MAPL_GetResource(MAPL, QCVAR_,   'QCVAR:'  , DEFAULT= 2.0 ,RC=STATUS) !variance of the QL distribution     
    VERIFY_(STATUS)
    qcvarr8=QCVAR_
    call MAPL_GetResource(MAPL, WBFFACTOR,   'WBFFACTOR:', DEFAULT= 1.0 ,RC=STATUS) !variance of the QL distribution     
    VERIFY_(STATUS)
    micro_mg_berg_eff_factor_in = WBFFACTOR
    call MAPL_GetResource(MAPL, NC_CST ,  'NC_CST:' , DEFAULT=  0.0 ,RC=STATUS) !constant nd (set if greather than zero)     
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, NI_CST ,  'NI_CST:' , DEFAULT=  0.0 ,RC=STATUS) !constant nd (set if greather than zero) 
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, NG_CST ,  'NG_CST:' , DEFAULT=  0.0 ,RC=STATUS) !constant ng (set if greather than zero) 
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL, MUI_CST,  'MUI_CST:', DEFAULT= -1.0 ,RC=STATUS) !constant ng (set if greather than zero) 
    VERIFY_(STATUS)
    mui_cnstr8 =  MUI_CST
    ncnstr8 = NC_CST
    if  (NC_CST .gt. 0.0)  nccons =.true.
    ninstr8 = NC_CST
    if  (NI_CST .gt. 0.0)  nicons =.true.
    ngnstr8 = NC_CST
    if  (NG_CST .gt. 0.0)  ngcons =.true.
    if  (MGVERSION .gt. 1) then
        do_graupel = .false.
        if (MGVERSION .gt. 2) do_graupel = .true.
        call micro_mg_init(Dcsr8, do_graupel,  micro_mg_berg_eff_factor_in, &
                           nccons, nicons, ncnstr8, ninstr8, ngcons, ngnstr8, mui_cnstr8)
    else
        call ini_micro(Dcsr8, micro_mg_berg_eff_factor_in, &
                       nccons, nicons, ncnstr8, ninstr8, qcvarr8)
    end if
    call aer_cloud_init()
    call WRITE_PARALLEL ("INITIALIZED MG in non-generic GC INIT")

      call MAPL_GetResource(MAPL, LCCIRRUS,       'LCCIRRUS:',       DEFAULT= 500.0,  RC=STATUS) !Characteristic Length (m) of high freq gravity waves
      call MAPL_GetResource(MAPL, UISCALE,        'UISCALE:',        DEFAULT= 1.0,    RC=STATUS) !Scaling factor for sed vel of ice      
      call MAPL_GetResource(MAPL, LIU_MU,         'LIU_MU:',         DEFAULT= 2.0,    RC=STATUS) !Liu autoconversion parameter
      call MAPL_GetResource(MAPL, NPRE_FRAC,      'NPRE_FRAC:',      DEFAULT= -1.0,   RC=STATUS) !Fraction of preexisting ice affecting ice nucleationn            
      call MAPL_GetResource(MAPL, CLDPARAMS%MIN_LTS,'LTS_LOW:',      DEFAULT= 20.0,   RC=STATUS) !lower LTS for morphology correction
      call MAPL_GetResource(MAPL, LTS_UP,         'LTS_UP:',         DEFAULT= 22.0,   RC=STATUS) !Upper LTS for morphology correction
      call MAPL_GetResource(MAPL, MIN_EXP,        'MIN_EXP:',        DEFAULT= 0.5,   RC=STATUS) !Exponent of the relation CFA=CFV^n
      call MAPL_GetResource(MAPL, MAX_EXP,        'MAX_EXP:',        DEFAULT= 1.0,   RC=STATUS) !Exponent of the relation CFA=CFV^n
      call MAPL_GetResource(MAPL, USE_AV_V,       'USE_AV_V:',       DEFAULT= 1.0,    RC=STATUS) !Set to > 0 to use an average velocity for activation
      call MAPL_GetResource(MAPL, AUTSC,          'AUT_SCALE:',      DEFAULT= 0.5,    RC=STATUS) !scale factor for critical size for drizzle
      call MAPL_GetResource(MAPL, TS_AUTO_ICE,    'TS_AUTO_ICE:',    DEFAULT= 4.0, RC=STATUS) !Ice autoconversion time scale
      call MAPL_GetResource(MAPL, TMAXLL,         'TMAXLL:',         DEFAULT= 250.0,  RC=STATUS) !Liquid clouds min T
      call MAPL_GetResource(MAPL, CCN_PARAM,      'CCNPARAM:',       DEFAULT= 2.0,    RC=STATUS) !CCN activation param
      call MAPL_GetResource(MAPL, IN_PARAM,       'INPARAM:',        DEFAULT= 6.0,    RC=STATUS) !IN param
      call MAPL_GetResource(MAPL, Immersion_param,'ImmersionPARAM:', DEFAULT= 6.0,    RC=STATUS) !Immersion param
      call MAPL_GetResource(MAPL, ACC_ENH,        'ACC_ENH:',        DEFAULT= 1.0,    RC=STATUS) !accretion rain-liquid scaling for MG2
      call MAPL_GetResource(MAPL, ACC_ENH_ICE,    'ACC_ENH_ICE:',    DEFAULT= 1.0,    RC=STATUS) !accretion snow-ice scaling for MG2
      call MAPL_GetResource(MAPL, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.5,    RC=STATUS) !Fraction of dust within droplets for immersion freezing
      call MAPL_GetResource(MAPL, FDROP_SOOT,     'FDROP_SOOT:',     DEFAULT= 0.05,   RC=STATUS) !Fraction of soot within droplets for immersion freezing        
      call MAPL_GetResource(MAPL, SIGMA_NUC,      'SIGMA_NUC:',      DEFAULT= 1.0,   RC=STATUS) !Widht of the in-cloud distribution of relative humidity in cirrus
      call MAPL_GetResource(MAPL, MIN_ALH,        'MIN_ALH:',        DEFAULT= 5.0,  RC=STATUS) !scale factor for vertical velocity in sttratocumulus
      call MAPL_GetResource(MAPL, SCWST,          'SCWST:',          DEFAULT= 3.0,  RC=STATUS) !scale factor for vertical velocity in sttratocumulus
      call MAPL_GetResource(MAPL, MINCDNC,          'MINCDNC:',      DEFAULT= 0.0,  RC=STATUS) !min nucleated droplet conc. cm-3
      call MAPL_GetResource(MAPL, TMAXCFCORR,         'TMAXCFCORR:',     DEFAULT= 285.0,  RC=STATUS) !Minimum T for CF correction
      call MAPL_GetResource(MAPL, MTIME,         'MTIME:',  DEFAULT= -1.0,    RC=STATUS) !Mixing time scale for aerosol within the cloud. Default is time step
      call MAPL_GetResource(MAPL, SWCIRRUS, 'SWCIRRUS:', DEFAULT= 3.0, RC=STATUS) !Tunes vertical velocity in cirrus
      call MAPL_GetResource(MAPL, DUST_INFAC,    'DUST_INFAC:',        DEFAULT= 1.0,   RC=STATUS)  !work on this
      call MAPL_GetResource(MAPL, BC_INFAC,        'BC_INFAC:',        DEFAULT= 0.1,   RC=STATUS)
      call MAPL_GetResource(MAPL, ORG_INFAC,      'ORG_INFAC:',        DEFAULT= 1.0,   RC=STATUS)
      call MAPL_GetResource(MAPL, SS_INFAC,        'SS_INFAC:',        DEFAULT= 1.0,   RC=STATUS)
      call MAPL_GetResource(MAPL, DT_MICRO,       'DT_MICRO:',        DEFAULT= 300.0,   RC=STATUS)    ! time step of the microphysics substepping (s) (MG2) (5 min)
      call MAPL_GetResource(MAPL, UR_SCALE,        'URSCALE:',        DEFAULT= 1.0,    RC=STATUS) !Scaling factor for sed vel of rain    
      call MAPL_GetResource(MAPL, USE_NATURE_WSUB,     'USE_NAT_WSUB:',     DEFAULT= 1.0  ,RC=STATUS) !greater than zero reads wsub from nature run                     
      call MAPL_GetResource(MAPL, DCS, 'DCS:', default=350.0e-6, RC=STATUS )
      call MAPL_GetResource(MAPL, CLDPARAMS%DISP_FACTOR_LIQ,         'DISP_FACTOR_LIQ:',     DEFAULT= 1.0,   RC=STATUS) ! Scales the droplet/ice crystal number in convective detrainment 
      call MAPL_GetResource(MAPL, CLDPARAMS%DISP_FACTOR_ICE,         'DISP_FACTOR_ICE:',     DEFAULT= 1.0,   RC=STATUS) ! Scales the droplet/ice crystal number in convective detrainment 

    call MAPL_GetResource (MAPL, JASON_TUNING, trim(COMP_NAME)//"_JASON_TUNING:", default=0, RC=STATUS); VERIFY_(STATUS)
#endif

end subroutine MGB2_2M_Initialize

subroutine MGB2_2M_Run (GC, IMPORT, EXPORT, CLOCK, RC)
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

    call MAPL_TimerOn (MAPL,"--MGB2_2M",RC=STATUS); VERIFY_(STATUS)

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
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL, 'NGRAUPEL'    , RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetPointer(IMPORT, TAUGWX, 'TAUGWX'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUGWY, 'TAUGWY'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUX,   'TAUX'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUY,   'TAUY'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OMEGA,  'OMEGA'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, ALH,    'ALH'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RADLW,  'RADLW'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, RADSW,  'RADSW'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUOROX, 'TAUOROX'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TAUOROY, 'TAUOROY'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WSUB_NATURE,  'WSUB_NATURE'     , RC=STATUS); VERIFY_(STATUS)

         !============================================= Start Stratiform cloud processes==========================================
         !set up initial values

         ! Find Convective Cloud Top
         call find_l(KCT, CNV_DQLDT, 1.0e-9, IM, JM, LM, 20, LM-2)
         ! find the minimun level for cloud micro calculations
         call find_l(KMIN_TROP, PLO, pmin_trop, IM, JM, LM, 10, LM-2)

         TEMP    = TH1*PK   
         DTDT_macro=TEMP
         DQVDT_macro=Q1
         DQLDT_macro=QLCN+QLLS
         DQIDT_macro=QICN+QILS
         DQADT_macro=CLCN+CLLS
         DQRDT_macro=QRAIN
         DQSDT_macro=QSNOW
 
         SC_ICE=1.0
         NCPL=MAX( NCPL , 0. )
         NCPI=MAX( NCPI , 0. )
         CLDREFFR = 10.0e-6 
         CLDREFFS = 90.0e-6
         CLDREFFG = 90.0e-6
         CLDREFFI = 25.0e-6
         CLDREFFL = 10.0e-6
         RAD_CF   = min(CLLS+CLCN, 1.0)
         RAD_QL   = 0.0
         RAD_QI   = 0.0
         RAD_QR   = 0.0
         RAD_QS   = 0.0
         RAD_QV   = Q1
         CDNC_NUC = 0.0
         INC_NUC  = 0.0
         T_ICE_MAX = MAPL_TICE     ! -7.0+TICE DONIFF
         T_ICE_ALL = CLDPARAMS%ICE_RAMP + MAPL_TICE 
         QSNOW_CN = 0.0
         QRAIN_CN = 0.0
         PFRZ= 0.0

      ! Find estimated inversion strength (DONIF)

       call    FIND_EIS(TH1, QSS, TEMP, ZLE, PLO, KLCL, IM, JM, LM, LTS, EIS)

         ! Clean up any negative specific humidity before the microphysics scheme
         !-----------------------------------------
         call FILLQ2ZERO2( Q1, MASS, FILLQ)  

          PFL_AN_X   = 0.0 
          PFL_LS_X  = 0.0

         !=======================================================================================================================
         !=======================================================================================================================
         !===================================Nucleation of cloud droplets and ice crystals ======================================
         ! Aerosol cloud interactions. Calculate maxCCN tendency using Fountoukis and nenes (2005) or Abdul Razzak and Ghan (2002)
         ! liquid Activation Parameterization
         ! Ice activation follows the Barahona & Nenes ice activation scheme, ACP, (2008, 2009). 
         ! Written by Donifan Barahona and described in Barahona et al. (2013)
         !=======================================================================================================================
         !=======================================================================================================================
         !=======================================================================================================================

         call MAPL_TimerOn(MAPL,"---ACTIV") !Activation timer

       if (NPRE_FRAC > 0.0) then
         NPRE_FRAC_2d(:,:) = NPRE_FRAC
       else
         ! include CNV_FRACTION dependence
         DO J=1, JM
            DO I=1, IM
            NPRE_FRAC_2d(I,J) = CNV_FRACTION(I,J)*ABS(NPRE_FRAC) + (1-CNV_FRACTION(I,J))*0.05
            END DO
         END DO
       endif

         use_average_v = .false.  
         if (USE_AV_V .gt. 0.0) then   
            use_average_v = .true.
         end if
        fdust_drop =  FDROP_DUST
        fsoot_drop =  FDROP_SOOT
        sigma_nuc_r8  = SIGMA_NUC
     frachet_org =         ORG_INFAC
     frachet_dust =         DUST_INFAC
     frachet_bc =         BC_INFAC
     frachet_ss =         SS_INFAC

        CFX=0.0
          where (QSS > 0.0) 
            CFX =Q1/(QSS)
          end where 
 !recalculate bkgtau: scaling of W variance with respect to Nature run
         if (USE_NATURE_WSUB .gt. 0.) then 
            xscale = (72000.0/imsize)            
            BKGTAU=  1.472/sqrt(1.0+ (xscale/6.0)) 
            BKGTAU = max((1.71 - BKGTAU), 0.0)*SWCIRRUS
             
      end if 
        
         do J=1,JM
            do I=1,IM

                            ccn_diag  = 0.0
                     smaxliq   = 0.0
                     smaxicer8 = 0.0
                     nheticer8 = 0.0
                     sc_icer8  = 1.0
                     naair8    = 0.0
                     npccninr8 = 0.0
                     nlimicer8 = 0.0    
                     nhet_immr8 = 0.0
                     dnhet_immr8 = 0.0
                     nhet_depr8 = 0.0
                     dust_immr8 = 0.0
                     dust_depr8 = 0.0
                      so4x = 0.0
                             dustx = 0.0
                      bcx= 0.0
                          orgx=0.0
                          seasaltx=0.0
                     
                     
                     
                     uwind_gw(1,1:LM)           = min(0.5*SQRT( U1(I,J,1:LM)**2+  V1(I,J,1:LM)**2), 50.0)
                     tausurf_gw   = min(0.5*SQRT(TAUOROX(I , J)**2+TAUOROY(I , J)**2), 10.0) !limit to a very high value     
                     if (USE_NATURE_WSUB .le. 0.) then 
                     tausurf_gw   =tausurf_gw  + min(0.5*SQRT(TAUX(I , J)**2+TAUY(I , J)**2), 5.0)*BKGTAU !adds a minimum value from unresolved sources (rewritten 04/01/15)
                     end if 
                     
                     
                     
                      aux1=PLE(i,j,LM)/(287.04*(T(i,j,LM)*(1.+0.608*Q1(i,j,LM)))) ! air_dens (kg m^-3)
                      hfs = -SH  (i,j) ! W m^-2
                      hfl = -EVAP(i,j) ! kg m^-2 s^-1
                      aux2= (hfs/MAPL_CP + 0.608*T(i,j,LM)*hfl)/aux1 ! buoyancy flux (h+le)
                       aux3= ZLE(I, J, KPBL(I,J))           ! pbl height (m)
                       !-convective velocity scale W* (m/s)
                       ZWS(i,j) = max(0.,0.001-1.5*0.41*MAPL_GRAV*aux2*aux3/T(i,j,LM))
                       ZWS(i,j) = 1.2*ZWS(i,j)**0.3333 ! m/s      
             
             
                     
                     pi_gw(1, 0:LM)        = 100.0*CNV_PLE(I,J,0:LM)                     
                     theta_tr(1, 1:LM)     = TH1(I,J,1:LM)
                             rhoi_gw = 0.0  
                     ni_gw   = 0.0 
                     ti_gw   = 0.0                                         
                     ter8(1,1:LM)      = TEMP(I,J,1:LM)   
                             pi_gw(1, 0:LM)   = 100.0*CNV_PLE(I,J,0:LM) 
                     plevr8(1,1:LM)    = 100.*PLO(I,J,:)
                     ndropr8(1, 1:LM) = NCPL(I, J, 1:LM)
                            qir8(1, 1:LM)     =  QILS(I, J,1:LM)+QICN(I, J,1:LM)
                    qcr8(1, 1:LM)     =  QLLS(I, J,1:LM)+QLCN(I, J,1:LM)
                    
                  ! where (RAD_CF(I, J, 1:LM) .gt. 0.01)
                  !  npre8(1,1:LM)     = NPRE_FRAC_2d(I,J)*NCPI(I,J,1:LM)/RAD_CF(I, J, 1:LM)
                  ! elsewhere 
                    npre8(1,1:LM)     = NPRE_FRAC_2d(I,J)*NCPI(I,J,1:LM)
                  ! end where 
                    
                    
                    omegr8(1,1:LM)    = OMEGA(I,J,1:LM)  
                    
                    rad_cooling(1,1:LM) = RADLW(I,J,1:LM)+RADSW(I,J,1:LM)
                    wparc_ls = 0.0
                    wparc_gw = 0.0
                    wparc_cgw= 0.0
                    wparc_turb = 0.0
                    swparc=0.0
                    tm_gw =ter8
                    pm_gw =plevr8
                    pfrz_inc_r8 = 0.0
                    Ksa1= 1.0
                    
                    if (FRLAND(I, J) .lt. 0.1) then 
                         lc_turb(1,1:LM)   =  max(ALH(I,J,1:LM), MIN_ALH) 
                     else
                        lc_turb(1,1:LM)   =  max(ALH(I,J,1:LM), 50.0)
                    end if 
                         
                    where ((npre8 .gt. 0.0)   .and. (qir8 .gt. 0.0))
                         dpre8    = ( qir8/(5400.0*npre8*MAPL_PI))**(0.33) !Assume exponential distribution
                    elsewhere
                        dpre8=1.0e-9
                   end where

                   call   gw_prof (1, LM, 1, tm_gw, pm_gw, pi_gw, &
                                 rhoi_gw, ni_gw, ti_gw, nm_gw) !get Brunt_Vaisala Frequency and midpoint densities 

                    kcldtopcvn=KCT(I, J)
                            Nct =nm_gw(1, kcldtopcvn)      !BV frequency ar cloud top
                    Wct = max(CNV_CVW(I, J, kcldtopcvn), 0.0)
                    fcn = maxval(CNV_UPDF(I, J, kcldtopcvn:LM))    

                            kbmin= min(KCBL(I, J), LM -1)-2
                    maxkhpbl=maxval(KH(I, J, kbmin:LM-1))    
                    
                    !ksa1= max(real(kbmin - kcldtopcvn), 1.0)

               ! ==========================================================================================    
               ! ========================Activate the aerosols ============================================ 
           
               
               
                do K = KMIN_TROP(I, J), LM-1 !limit to troposphere and no activation at the surface
    
                ! find vertical velocity variance 
!                       call zeit_ci("MOIST::aero_vvar")

                       call vertical_vel_variance(omegr8(1, K), lc_turb(1, K), ter8(1, K), plevr8(1, K), rad_cooling(1,K),  uwind_gw(1,K), &
                                                         tausurf_gw, nm_gw(1, K), LCCIRRUS, Nct, Wct, &
                                                         ksa1, fcn(1, K), KH(I, J, K), FRLAND(I, J), ZPBL(I, J), ZLE(I, J, k), maxkhpbl, &
                                                            wparc_ls(1, K), wparc_gw(1, K), wparc_cgw(1, K), wparc_turb(1, K), EIS(I, J), TKE(I, J, K))
                     
!                               call zeit_co("MOIST::aero_vvar")
                                        
                   if (FRLAND(I, J) .lt. 0.1) then 
                    if (LTS(I, J) .gt. CLDPARAMS%MIN_LTS) then                     
                           if (K .ge. kbmin-2) wparc_ls(1, K)=max(wparc_ls(1,K)+ zws(i, j), 0.00)*SCWST ! add convective velocity within the PBL
                    end if 
                   else
                      if (K .ge. kbmin-2) wparc_ls(1, K)=max(wparc_ls(1,K)+ zws(i, j), 0.00) 
                   end if 

                     if (K .ge. kbmin-2) wparc_turb(1, K)=max(wparc_turb(1,K), 0.04)    !minimum velocity within the PBL (not resolved by RAS)
                                                       
                     if (K .ge.  kcldtopcvn) wparc_cgw(1, K) = 0.0                     
                    

                        
                         if (USE_NATURE_WSUB .gt. 0.) then !use climatology from the Nature run (only for cirrus)
                                 
                                  !wparc_cgw(1, k)= max(WSUB_NATURE(I, J, K)+BKGTAU*BKGTAU, 0.0)!BKG accounts for unresolved vertical velocity at 7 km                            
      
                                   wparc_cgw(1, k)= max(WSUB_NATURE(I, J, K)*BKGTAU*BKGTAU, 0.0)!BKG accounts for unresolved vertical velocity at 7 km                            
      
                                   wparc_gw(1, k) = 0.0

                        end if 

                             swparc(1, K)=sqrt(wparc_gw(1, K)+wparc_turb(1, K)+ wparc_cgw(1, K))
 
                        
                                
                       !Supersaturations to calculate CCN diagnostics
                        ccn_diag(1)=0.001
                        ccn_diag(2)=0.004
                        ccn_diag(3)=0.01

                         AeroAux%nmods = 0
                         AeroAux%num   = 0.0
                         do i_src_mode = 1, AeroProps(I,J,K)%nmods
                            if (AeroProps(I,J,K)%num(i_src_mode) > 0.1) then
                               AeroAux%nmods = AeroAux%nmods + 1
                               i_dst_mode = AeroAux%nmods

                               AeroAux%num(i_dst_mode)   = AeroProps(I,J,K)%num(i_src_mode)
                               AeroAux%dpg(i_dst_mode)   = AeroProps(I,J,K)%dpg(i_src_mode)
                               AeroAux%sig(i_dst_mode)   = AeroProps(I,J,K)%sig(i_src_mode)
                               AeroAux%den(i_dst_mode)   = AeroProps(I,J,K)%den(i_src_mode)
                               AeroAux%kap(i_dst_mode)   = AeroProps(I,J,K)%kap(i_src_mode)
                               AeroAux%fdust(i_dst_mode) = AeroProps(I,J,K)%fdust(i_src_mode)
                               AeroAux%fsoot(i_dst_mode) = AeroProps(I,J,K)%fsoot(i_src_mode)
                               AeroAux%forg(i_dst_mode)  = AeroProps(I,J,K)%forg(i_src_mode)
                             end if
                        end do

                       rh1_r8=CFX(I, J, K)
                       
                       !if ((K.gt. 2) .and. (K .lt. LM-1)) then 
                       !    tauxr8 =  (ter8(1, K+1) +  ter8(1, K) + ter8(1, K-1))/3.0
                           !else 
                             tauxr8 = ter8(1, K)
             ! end if 
                             
!                        call zeit_co("MOIST::aero_unpack")
                                   
                     !!Subroutine aerosol_activate contains the CCN activation and ice nucleation parameterizations. Lives in aer_cloud.F90.

                     call   aerosol_activate(tauxr8, plevr8(1, K), swparc(1, K), wparc_ls(1, K),  AeroAux, &
                          npre8(1, k), dpre8(1, k), ccn_diag, ndropr8(1, k), qcr8(1, K), &
                          npccninr8(1, K), smaxliq(1, K), naair8(1, K), smaxicer8(1, K), nheticer8(1, K), &
                          nhet_immr8(1, K), dnhet_immr8(1, K), nhet_depr8(1, k), sc_icer8(1, k), &
                          dust_immr8(1, K), dust_depr8(1, k), nlimicer8(1, k), use_average_v, int(CCN_PARAM), int(IN_PARAM),  &
                          so4x(1, k), seasaltx(1, k), dustx(1, k), orgx(1, K), bcx(1, k), &                          
                                      fdust_drop, fsoot_drop, pfrz_inc_r8(1, K), rh1_r8, frachet_dust, frachet_bc, frachet_org, frachet_ss, int(Immersion_PARAM))

                      CCN01(I, J, K) = max(ccn_diag(1), 0.0)
                      CCN04(I, J, K) = max(ccn_diag(2), 0.0)
                      CCN1 (I, J, K) = max(ccn_diag(3), 0.0)
                      
                      if (K .ge. kbmin-6) npccninr8(1, K) = max(npccninr8(1, K), (1.0-CNV_FRACTION(I, J))*MINCDNC*1.e6)
                       
               end do

               SMAXL(I, J, 1:LM) = real(smaxliq(1, 1:LM)*100.0)         
               SMAXI(I, J, 1:LM) = real(smaxicer8(1, 1:LM)*100.0)
               NHET_NUC(I, J, 1:LM)  = real(nheticer8(1, 1:LM))
               NLIM_NUC(I, J, 1:LM) =  real(nlimicer8(1, 1:LM))            
               SC_ICE(I, J, 1:LM) = real(sc_icer8(1, 1:LM))                  
               CDNC_NUC(I,J,1:LM)    = real(npccninr8(1, 1:LM))
               INC_NUC (I,J,1:LM)    = real(naair8(1, 1:LM)  )       
               NHET_IMM(I, J, 1:LM)  = real(max(nhet_immr8(1, 1:LM), 0.0))
               DNHET_IMM(I, J, 1:LM)  = real(max(dnhet_immr8(1, 1:LM), 0.0))
               NHET_DEP(I, J, 1:LM)  = real(nhet_depr8(1, 1:LM))
               DUST_IMM(I, J, 1:LM)  = real(max(dust_immr8(1, 1:LM), 0.0))
               DUST_DEP(I, J, 1:LM)  = real(max(dust_depr8(1, 1:LM), 0.0))
               WSUB (I, J, 1:LM) =  real(wparc_ls(1, 1:LM)+swparc(1, 1:LM)*0.8)        
               SIGW_GW (I, J, 1:LM)   = real( wparc_gw(1, 1:LM))
               SIGW_CNV (I, J, 1:LM)   =  real(wparc_cgw(1, 1:LM))
               SIGW_TURB (I, J, 1:LM) = real(wparc_turb(1, 1:LM))
               SIGW_RC (I, J, 1:LM)   =  real(wparc_ls(1, 1:LM))
                   PFRZ (I, J, 1:LM)   =  real(pfrz_inc_r8(1, 1:LM))
                
               SO4(I, J, 1:LM)=real(so4x(1, 1:LM))
                      DUST(I, J, 1:LM)=real(dustx(1, 1:LM))
               BCARBON(I, J, 1:LM)=real(bcx(1, 1:LM))
                   ORG(I, J, 1:LM)=real(orgx(1, 1:LM))
                   SEASALT(I, J, 1:LM)=real(seasaltx(1, 1:LM))

               
            enddo
         enddo

    
        
         !=============================================End cloud particle nucleation=====================================
         !===============================================================================================================

         call MAPL_TimerOff(MAPL,"---ACTIV", RC=STATUS)
         call MAPL_TimerOn(MAPL,"---CLDMACRO")

         !==========================================================================================================
         !===================================Cloud Macrophysics ====================================================
         !==========================================================================================================

         REV_CN_X  = 0.0
         REV_AN_X  = 0.0 
         REV_LS_X  = 0.0
         REV_SC_X  = 0.0
         RSU_CN_X  = 0.0
         RSU_AN_X  = 0.0
         RSU_LS_X  = 0.0     
         RSU_SC_X  = 0.0     

         CFX=INC_NUC + NHET_IMM 
      
         CN_PRC2 = 0.0 
         LS_PRC2= 0.0
         AN_PRC2 = 0.0 
         SC_PRC2 = 0.0 
         CN_SNR  =0.0 
         LS_SNR = 0.0  
         AN_SNR = 0.0
         SC_SNR = 0.0
         DTDT_macro=TEMP     
         DQVDT_macro=Q1
         DQLDT_macro=QLCN+QLLS
         DQIDT_macro=QICN+QILS
         DQADT_macro=CLCN+CLLS
         DQRDT_macro=QRAIN
         DQSDT_macro=QSNOW
         PFI_CN_X= 0.0
         PFI_AN_X=0.0
         PFI_LS_X=0.0
         PFI_SC_X=0.0
         PFL_CN_X= 0.0
         PFL_AN_X=0.0
         PFL_LS_X=0.0    
         PFL_SC_X=0.0    
      

    ! if(associated(TVQX1))  TVQX1     = SUM( (  Q1 +  QLLS + QLCN + QILS + QICN + CNV_PRC3)*DM + CNV_DQLDT*DT_MOIST , 3 )
    
         
              CNV_MFD_X      =  CNV_MFD     ! needed for cloud fraction
              CNV_DQLDT_X    =  CNV_DQLDT
              CNV_PRC3_X     = CNV_PRC3
              CNV_UPDF_X     = CNV_UPDF 
      
     IF(ADJUSTL(CONVPAR_OPTION) == 'GF') THEN    ! GF updates the state internally so we don't need to do that here
              CNV_PRC3_X     = 0.0
              CNV_MFD_X    =  0.0     
              CNV_DQLDT_X    =  0.0
              CNV_NICE_X =  0.0
              CNV_NICE_X  =  0.0
     END IF 
        
      
              if(associated(TVQX1))  TVQX1     =  SUM( (  Q1 +  QLLS + QLCN + QILS + QICN + QRAIN + QSNOW + QGRAUPEL + SHLW_PRC3 + SHLW_SNO3)*MASS &    
                 + (CNV_DQLDT)*DT_MOIST &            
                 + (QLDET_SC  + QIDET_SC)*DT_MOIST &                  
                 , 3 ) + CNVPRCP*DT_MOIST - TVQ0 ! up to here water is conserved Donif 01/2020
                 
           
 
                       
  call  macro_cloud (                    &
              IM*JM, LM         , &
              DT_MOIST          , &
              PLO               , &
              CNV_PLE           , &
              PK                , &
              SNOMAS            , &   ! <- surf
              FRLANDICE         , &   ! <- surf
              FRLAND            , &   ! <- surf
              KH                , &   ! <- turb
              CNV_MFD_X           , &   ! <- ras
              CNV_DQLDT_X         , &   ! <- ras              
              CNV_PRC3_X          , &   ! <- ras   
              CNV_UPDF_X          , &   ! <- ras
              MFD_SC            , &   ! <- shcu   
              QLDET_SC          , &   ! <- shcu   
              QIDET_SC          , &   ! <- shcu   
              SHLW_PRC3         , &   ! <- shcu   
              SHLW_SNO3         , &   ! <- shcu   
              UFRC_SC           , &   ! <- shcu 
              U1                , &
              V1                , & 
              TH1               , &              
              Q1                , &
              QLLS              , &
              QLCN              , &
              QILS              , &
              QICN              , &
              CLCN              , &
              CLLS              , &           
              CN_PRC2           , &            
              CN_ARFX           , &
              CN_SNR            , &
              CLDPARAMS         , &
              QST3              , &
              DZET              , &
              CNV_FRACTION      ,  &
              QDDF3             , &
                                ! Diagnostics
              RHX_X             , &
              REV_CN_X          , &
              RSU_CN_X          , &
              ACLL_CN_X,ACIL_CN_X   , &
              PFL_CN_X,PFI_CN_X     , &
              DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
              DCNVL_X, DCNVi_X,       &
              ALPHT_X, CFPDF_X, &
              DQRL_X,               &
              VFALLSN_CN_X,  &
              VFALLRN_CN_X,  &
              EVAPC_X , SUBLC_X ,  &
                                ! End diagnostics
            !!====2-Moment============
              CNV_FICE,  &
              CNV_NDROP_X, &
              CNV_NICE_X,  &
              SC_NDROP,  &
              SC_NICE,   &
              SC_ICE,    &
              NCPL, &
              NCPI, &
              PFRZ, &
              DNDCNV_X, &
              DNCCNV_X, &
              DT_RASP , &
              QRAIN_CN, & !grid av
              QSNOW_CN, &
              KCBL, LTS,  CONVPAR_OPTION)

      

       TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
               CN_SNR  + LS_SNR  + AN_SNR + SC_SNR
      
      if(associated(TVQX2)) TVQX2    = SUM( ( Q1 +  QLLS + QLCN + QILS + QICN +  QRAIN +  QSNOW + QGRAUPEL)*MASS , 3 )  + TPREC*DT_MOIST -TVQ0

         TEMP    = TH1*PK

         DTDT_macro=  (TEMP-DTDT_macro)/DT_MOIST
         DQVDT_macro=(Q1-DQVDT_macro)/DT_MOIST
         DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
         DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
         DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
         DQRDT_macro=(QRAIN-DQRDT_macro)/DT_MOIST
         DQSDT_macro=(QSNOW-DQSDT_macro)/DT_MOIST
      
         !make sure QI , NI stay within T limits 
         call meltfrz_inst  (     &
              IM,JM,LM    , &
              TEMP              , &
              QLLS          , &
              QLCN         , &
              QILS           , &
              QICN          , &               
              NCPL         , &
              NCPI          )

         
        call fix_up_clouds_2M( &
         Q1, &
         TEMP, &
         QLLS,&
         QILS,&
         CLLS, &
         QLCN,&
         QICN,&
         CLCN, &
         NCPL, &
         NCPI, & 
         QRAIN, &
         QSNOW, &
         QGRAUPEL, &
         NRAIN, &
         NSNOW, &
         NGRAUPEL)
         
         
         ! need to clean up small negative values. MG does can't handle them
          call FILLQ2ZERO2( Q1, MASS, FILLQ) 
          call FILLQ2ZERO2( QGRAUPEL, MASS, FILLQ) 
          call FILLQ2ZERO2( QRAIN, MASS, FILLQ) 
          call FILLQ2ZERO2( QSNOW, MASS, FILLQ) 
          call FILLQ2ZERO2( QLLS, MASS, FILLQ)
          call FILLQ2ZERO2( QLCN, MASS, FILLQ)  
          call FILLQ2ZERO2( QILS, MASS, FILLQ)
          call FILLQ2ZERO2( QICN, MASS, FILLQ)
         
         !=============================================End cloud macrophysics=====================================
         !======================================================================================================================
         !

         FQAI = 0.0
         FQAL = 0.0
         FQA  = 0.0 
         QCNTOT = QLCN+QICN
         QTOT   =  QCNTOT+ QLLS+QILS
         QL_TOT  = QLCN+QLLS 
         QI_TOT  = QICN+QILS   

         where (QTOT .gt. 0.0)
            FQA= min(max(QCNTOT/QTOT, 0.0), 1.0)    
         end where
         
         CFLIQ=0.0
         CFICE=0.0
       
         RAD_CF   = min(CLLS+CLCN, 1.0)

         WHERE (QTOT .gt. 0.0) 
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
         END WHERE

      
      
            INC_NUC = INC_NUC*PFRZ!!Accounts for ice crystal dilution after nucleation. 
            NHET_NUC = NHET_NUC*PFRZ!
      

         !==================================================================================================================
         !===============================================Two-moment stratiform microphysics ================================
         !================This is the implementation of the Morrison and Gettelman (2008) microphysics =====================
         !==================================================================================================================

         rhdfdar8   = 1.e-8_r8
         rhu00r8    = 0.95_r8
         ttendr8=0._r8
         qtendr8=0._r8
         cwtendr8=0._r8
         naair8=0.
         rndstr8 = 2.0e-7
         npccninr8 = 0.
         naconr8   = 0
         scale_ri =  1.3 ! scaling factor to account for the different definition of Ri in Chao and Suarez

        if ((RRTMG_SORAD .gt. 0.0) .or. (RRTMG_IRRAD .gt. 0.0)) then 
        scale_ri =  1.0
        end if 

         call MAPL_TimerOff(MAPL,"---CLDMACRO", RC=STATUS)
         call MAPL_TimerOn (MAPL,"---CLDMICRO", RC=STATUS)
 
             !initialize MG variables
         nimmr8 = 0.0_r8
         cldfr8 = 0.0_r8 
         prectr8 = 0.0_r8 
             precir8 = 0.0_r8
             qctendr8 = 0.0_r8
             qitendr8 = 0.0_r8
             qvlatr8 = 0.0_r8
             tlatr8 = 0.0_r8
             nctendr8 = 0.0_r8
             nitendr8 = 0.0_r8
             effcr8 = 0.0_r8
             effir8 = 0.0_r8
             drout2r8 =0.0_r8
             dsout2r8 = 0.0_r8
         dgout2r8 = 0.0_r8
             qrout2r8 = 0.0_r8
             qsout2r8 =0.0_r8
         qgout2r8 =0.0_r8
             nrout2r8 = 0.0_r8
             nsout2r8 =0.0_r8
         ngout2r8 =0.0_r8         
             evapsnowr8 =0.0_r8
             nevaprr8 =0.0_r8
             cmeioutr8 =0.0_r8
             bergsor8 =0.0_r8
             mnucccor8 =0.0_r8
             mnucctor8 =0.0_r8
             homoor8 = 0.0_r8
             mnuccror8 = 0.0_r8
             pracsor8 = 0.0_r8
             meltor8 =0.0_r8
             qisedtenr8 =0.0_r8
             bergor8 =0.0_r8
             psacwsor8 = 0.0_r8
             qcresor8 =0.0_r8
             qiresor8 = 0.0_r8
             praor8 =0.0_r8
             prcor8 = 0.0_r8
             prcior8 =0.0_r8
             praior8 = 0.0_r8
             msacwior8 =0.0_r8
             frzrdtr8 =0.0_r8
             meltsdtr8 = 0.0_r8
             nnucctor8 =0.0_r8
             nnucccor8 = 0.0_r8
             nnuccdor8 =0.0_r8
             nsacwior8 =0.0_r8
             nsubior8 = 0.0_r8
             npraior8 =0.0_r8
             nprcior8 =0.0_r8
             npccnor8 = 0.0_r8
             npsacwsor8 =0.0_r8
             npraor8 =0.0_r8
             nsubcor8 =0.0_r8
             nprc1or8 =0.0_r8
         rndstr8 = 2.0e-7
         naconr8   = 0.
     
         lflxr8 = 0.0_r8             
         iflxr8 = 0.0_r8
         rflxr8 = 0.0_r8
         sflxr8 = 0.0_r8
         gflxr8 = 0.0_r8    

         frzcntr8 =0.0_r8 
         qrtendr8 =  0.0_r8
         nrtendr8 =  0.0_r8
         qstendr8 =  0.0_r8
         nstendr8 =  0.0_r8
     
         qgtendr8 =  0.0_r8
         ngtendr8 =  0.0_r8
     
         !accre_enhanr8= 1.0_r8 
         
         accre_enhanr8= ACC_ENH
         accre_enhan_icer8= ACC_ENH_ICE
         AN_PRC2     = 0. !prectr8(1)
         AN_SNR      = 0. !precir8(1)
         AN_ARFX     = 0. !maxval( cldfr8(1,1:LM) )    
         PFL_LS_X = 0.0
         PFI_LS_X= 0.0
         QCVAR_EXP = 2.0
         autscx = 1.0
         
         do J=1,JM
            do I=1,IM

              
                   kbmin =1            
               npccninr8  = 0.0
               naair8     = 0.0
               omegr8     = 0.0
               rndstr8 = 2.0e-7
               naconr8   = 0.

               cldfr8(1,1:LM)      = RAD_CF(I,J,1:LM) !Assume minimum overlap 
             
               liqcldfr8(1, 1:LM)  = CFLIQ(I, J,1:LM) 
               icecldfr8(1, 1:LM)  = CFICE(I, J,1:LM) 
                 
               cldor8           = cldfr8  
               ter8(1,1:LM)        = TEMP(I,J,1:LM)
               qvr8(1,1:LM)        = Q1(I,J,1:LM)

               qcr8(1,1:LM)        = QL_TOT(I,J,1:LM)
               qir8(1,1:LM)        = QI_TOT(I,J,1:LM)
               ncr8(1,1:LM)        = MAX(NCPL(I,J,1:LM), 0.0) 
               nir8(1,1:LM)        = MAX(NCPI(I,J,1:LM), 0.0) 

               ! Nucleation  tendencies 
               naair8(1, 1:LM)     = max((INC_NUC(I, J, 1:LM)*cldfr8(1,1:LM) - nir8(1,1:LM))/DT_MOIST, 0.0) 
               npccninr8(1, 1:LM)  = max((CDNC_NUC(I, J, 1:LM)*cldfr8(1,1:LM) - ncr8(1,1:LM))/DT_MOIST, 0.0)
               
               where  ((naair8 .gt. 1.0e3)) ! add cloud fraction if nucleation is happening 2018
                   icecldfr8 = max(0.05,  icecldfr8)
               end where 
             

               where (cldfr8(1,:) .ge. 0.001) 
                  nimmr8(1, 1:LM)     = MIN(DNHET_IMM(I, J, 1:LM), ncr8(1,1:LM)/cldfr8(1,1:LM)/DT_MOIST) !tendency    
               elsewhere 
                  nimmr8(1, 1:LM)   = 0.0 
               end where
               
               nhet_depr8(1, 1:LM) = NHET_DEP(I, J, 1:LM)/DT_MOIST !becomes a tendency (could be done a bit better)
               nbincontactdust = 1

               DO K=kbmin, LM
                  AeroAux = AeroProps(I, J, K)
                  ! Get dust properties for contact ice nucleation 
                  call getINsubset(1,  AeroAux, AeroAux_b) 
                  naux = AeroAux_b%nmods 
                  if (nbincontactdust  .lt. naux) then 
                     nbincontactdust =  naux 
                  end if
                  naconr8(1, K, 1:naux) =  AeroAux_b%num(1:naux)
                  rndstr8( 1, K, 1:naux)=AeroAux_b%dpg(1:naux)/2.0

                  ! Get black carbon properties for contact ice nucleation      
                  call init_Aer(AeroAux_b)   
                  call getINsubset(2,  AeroAux, AeroAux_b)          
                  nsootr8  (1, K)  = sum(AeroAux_b%num) !
                  naux = AeroAux_b%nmods 
                  rnsootr8 (1, K)  = sum(AeroAux_b%dpg(1:naux))/naux
               END DO

               pdelr8(1,1:LM)  = PLE(I,J,1:LM) - PLE(I,J,0:LM-1)  
               rpdelr8      = 1./pdelr8 
               pintr8(1,1:LM+1) = PLE(I,J,0:LM)  
               plevr8(1,1:LM)      = 100.*PLO(I,J,1:LM)
               zmr8(1,1:LM)        = ZLO(I,J,1:LM)     
               kkvhr8(1,1:LM+1) = KH(I,J,0:LM)  
               ficer8 = qir8 /( qcr8+qir8 + 1.e-10 )  
               omegr8(1,1:LM)=WSUB(I, J, 1:LM)
               
               
               !Tuning factors
               disp_liu = LIU_MU
               ui_scale = UISCALE
               urscale  = URSCALE
               ts_autice = DT_R8*TS_AUTO_ICE 
               
               
               if (AUTSC .gt. 0.0) then 
                  autscx = AUTSC
               else
                autscx =  min(max(0., (300.0 - TEMP(I,J,LM))/ABS(AUTSC)), 1.0)
                autscx  =  1.0 - 0.995*autscx
               end if
               
               if (MTIME .le. 0.0) then 
                   mtimesc  = DT_MOIST
               else               
                  mtimesc=MTIME
               end if 
  
  !!!!================Estimate qcvar following Xie and Zhang, JGR, 2015
                HMOIST_950 = 0.0
                HSMOIST_500 = 0.0
                              
                 IF (PLO(I, J, LM) .le. 500.0) then                                        
                      qcvarr8  = 2.0
                 ELSEIF (PLO(I, J, LM) .lt. 950.0) then 
                   
                    DO K=LM, 1, -1       
                         if (PLO(I,J,K) .lt. 500.0) exit  
                         HSMOIST_500 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO 
                            
                      HMOIST_950 = MAPL_CP*TEMP(I, J, LM) + GZLO(I, J, LM) + Q1(I, J, LM)*MAPL_ALHL               
                      SINST = (HMOIST_950 -  HSMOIST_500)/(PLO(I,J,LM)*100.0- 50000.0)                   
                  else
                     DO K=LM, 1, -1       
                         if (PLO(I,J,K) .lt. 500.0) exit  
                         HSMOIST_500 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO 

                    DO K=LM, 1, -1       
                     if (PLO(I,J,K) .lt. 950.0) exit  
                     HMOIST_950 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + Q1(I, J, K)*MAPL_ALHL
                    END DO                                          
                     SINST = (HMOIST_950 -  HSMOIST_500)/45000.0                  
               
                   end if  
               
                  xscale = (36000.0/imsize)**(-0.666)
                  qcvarr8 =  0.67 -0.38*SINST +  4.96*xscale - 8.32*SINST*xscale  
                  qcvarr8 = min(max(qcvarr8, 0.5), 50.0)
                  if (associated(QCVAR_EXP)) QCVAR_EXP(I, J) = real(qcvarr8)
                  relvarr8 = qcvarr8
                  
                
               ! for MG23 (initial values)     
                        frzimmr8 =  nimmr8
                        frzcntr8 = nimmr8*0.0  
                        frzdepr8 = nhet_depr8
                        qrr8(1, 1:LM)     =  QRAIN(I, J,1:LM)
                        qsr8(1, 1:LM)     =  QSNOW(I, J,1:LM)
                        qgr8(1, 1:LM)     =  QGRAUPEL(I, J,1:LM)                        
                        nrr8(1, 1:LM)     =  NRAIN(I, J,1:LM)
                        nsr8(1, 1:LM)     =  NSNOW(I, J,1:LM)
                        ngr8(1, 1:LM)     =  NGRAUPEL(I, J,1:LM)                         
                        qsatfacr8 = 1.0                        
                        SCICE_tmp(1, 1:LM)  =  SC_ICE(I, J, 1:LM)
                        FQA_tmp(1, 1:LM)  = FQA(I, J, 1:LM) 
                        ALPH_tmp(1, 1:LM)  = ALPHT_X(I, J, 1:LM)
                   
                   

 
     
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
  !CALLS to MG versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
!!!Call to MG microphysics. Lives in cldwat2m_micro.F90

               
   if (MGVERSION .lt. 2)  then          
               
               call set_qcvar (qcvarr8)
               
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
                    nsoutr8, nroutr8, nimmr8, disp_liu, &
                    nsootr8, rnsootr8, ui_scale, autscx, mtimesc, &
                    nnuccdor8, nnucccor8, nsacwior8, nsubior8, nprcior8, &
                    npraior8, npccnor8, npsacwsor8, nsubcor8, npraor8, nprc1or8,  nbincontactdust, &
                    ts_autice, rflxr8, sflxr8, accre_enhanr8, accre_enhan_icer8)

    else ! MG2/3
        
         call  micro_mg_tend_interface ( DT_MICRO, INT(CLDPARAMS%PDFSHAPE), ALPH_tmp, SCICE_tmp, FQA_tmp, &
                                        CNV_FRACTION(I, J), SNOMAS(I, J), FRLANDICE(I, J), FRLAND(I, J), & 
                             ncolmicro,             LM,               dt_r8,       & 
                             ter8,                            qvr8,                              &
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
!        !                     errstring, & ! Below arguments are "optional" (pass null pointers to omit).
                      !       tnd_qsnow,          tnd_nsnow,          re_ice,    &
                             prer_evap, &
                             frzimmr8,             frzcntr8,              frzdepr8,  & ! contact is not passed since it depends on the droplet size dist
                             nsootr8, rnsootr8,  & ! soot for contact IN
                             npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, &  ! Number tendencies for liquid
                             npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, nsacwior8,  &  ! Number tendencies for ice
                             ts_autice, ui_scale, autscx , disp_liu, nbincontactdust, urscale)

    end if 

        IF (MGVERSION > 1) then 

                   QRAIN   (I,J,1:LM)  = max(REAL(qroutr8(1, 1:LM)), 0.0)
                   QSNOW   (I,J,1:LM)  = max(REAL(qsoutr8(1, 1:LM)), 0.0)
                   QGRAUPEL(I,J,1:LM)  = max(REAL(qgoutr8(1, 1:LM)), 0.0)
                   NRAIN   (I,J,1:LM)  = max(REAL(nroutr8(1, 1:LM)), 0.0)
                   NSNOW   (I,J,1:LM)  = max(REAL(nsoutr8(1, 1:LM)), 0.0)
                   NGRAUPEL(I,J,1:LM)  = max(REAL(ngoutr8(1, 1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM)  = REAL(reff_rainr8(1, 1:LM))
                   CLDREFFS(I,J,1:LM)  = REAL(reff_snowr8(1, 1:LM))/scale_ri
                   CLDREFFG(I,J,1:LM)  = REAL(reff_graur8(1, 1:LM))/scale_ri
                   DQRL_X(I,J,1:LM)    = REAL(qroutr8(1, 1:LM)/DT_R8) !rain mixing ratio tendency from micro
            
        else
                    
                   QRAIN(I,J,1:LM)  = max(REAL(qrout2r8(1, 1:LM)), 0.0) ! grid average 
                   QSNOW(I,J,1:LM)  = max(REAL(qsout2r8(1, 1:LM)), 0.0)                      
                   NRAIN(I,J,1:LM)  = max(REAL(nrout2r8(1, 1:LM)), 0.0)
                   NSNOW(I,J,1:LM)  = max(REAL(nsout2r8(1, 1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM) = REAL(drout2r8(1, 1:LM))/2.0        
                   CLDREFFS(I,J,1:LM) = REAL(dsout2r8(1, 1:LM))/2.0/scale_ri
                   DQRL_X(I,J,1:LM)   = REAL(qrout2r8(1, 1:LM)/DT_R8) !rain mixing ratio tendency from micro
                 
         end if          
         
         
         
  
               PFL_LS_X(I, J, 1:LM) = rflxr8(1, 1:LM) !+ lflxr8(1, 1:LM)
               PFI_LS_X(I, J, 1:LM) = sflxr8(1, 1:LM) !+ gflxr8(1, 1:LM) +  iflxr8(1, 1:LM)
              
               !Update state after microphysisc
               LS_PRC2(I,J)     = max(1000.*REAL((prectr8(1)-precir8(1))), 0.0)
               LS_SNR(I,J)      = max(1000.*REAL(precir8(1)), 0.0)          
               QL_TOT(I,J,1:LM) = max(QL_TOT(I,J,1:LM)   + REAL(qctendr8(1,1:LM)) * DT_R8, 0.0)
               QI_TOT(I,J,1:LM) = max(QI_TOT(I,J,1:LM)   + REAL(qitendr8(1,1:LM)) * DT_R8, 0.0)    
               Q1(I,J,1:LM)   = MAX(Q1(I,J,1:LM)     + REAL(qvlatr8(1,1:LM)) * DT_R8, 0.0)
               TEMP(I,J,1:LM) = TEMP(I,J,1:LM)   + REAL(tlatr8(1,1:LM)) * DT_R8 / (MAPL_CP)  
               NCPL(I,J,1:LM) = MAX(NCPL(I,J,1:LM)   + REAL(nctendr8(1,1:LM)) * DT_R8, 0.0) 
               NCPI(I,J,1:LM) = MAX(NCPI(I,J,1:LM)   + REAL(nitendr8(1,1:LM)) * DT_R8, 0.0)  
                           
                       
                       

               LS_ARFX(I,J)     = maxval( REAL(cldfr8(1,1:LM)) )
                            
               CLDREFFL(I,J,1:LM) = max(REAL(effcr8(1,1:LM))*1.0e-6, 1.0e-6)             
               CLDREFFI(I,J,1:LM) = max(REAL(effir8(1,1:LM))*1.0e-6, 1.0e-6)/scale_ri !scale to match the Dge definition of Fu 1996                    

               ! diagnostics from the microphysics********************

               RSU_LS_X(I,J,1:LM) = REAL(evapsnowr8(1, 1:LM))                    !Snow evap   
               REV_LS_X(I,J,1:LM) = REAL(nevaprr8(1, 1:LM) )                        !rain evap
               SUBLC_X(I,J,1:LM) = REAL(cmeioutr8(1, 1:LM)) + SUBLC_X(I,J,1:LM)     ! Ice subl already grid -ave
               BERGS(I,J,1:LM) = REAL(bergsor8(1,1:LM))                               ! Snow Bergeron
               FRZ_TT_X(I,J,1:LM) = REAL(mnucccor8(1,1:LM)+ mnucctor8(1,1:LM) + homoor8(1,1:LM))!ice mixing ratio from nucleation (hom+het)
               FRZ_PP_X(I,J,1:LM) = REAL(mnuccror8(1, 1:LM) + pracsor8(1, 1:LM)) !freezing of rain (hom+ het freezing and accretion by snow) 
               MELT(I,J,1:LM) = REAL(meltor8(1,1:LM))                            !melting of cloud ice  and snow
               SDM_X(I,J,1:LM) = REAL(qisedtenr8(1, 1:LM))                    ! ice sed   
               EVAPC_X(I,J,1:LM) = REAL(qcsevapr8(1,1:LM) )  +  EVAPC_X(I,J,1:LM)   ! cloud evap
               BERG(I,J,1:LM)  =  REAL(bergor8(1,1:LM))                          ! Bergeron process  
               ACIL_LS_X(I,J,1:LM) =REAL(psacwsor8(1, 1:LM) + msacwior8(1, 1:LM))   !Acc + collection of cloud  by snow
               QCRES(I,J,1:LM) =REAL(qcresor8(1, 1:LM) )             !residual liquid condensation
               QIRES(I,J,1:LM) =REAL(qiresor8(1, 1:LM))                  !residual ice condensation   

               ACLL_LS_X(I,J,1:LM) =REAL(praor8(1, 1:LM) )          ! Acc cloud by rain    
               AUT_X(I,J,1:LM) = REAL(prcor8(1, 1:LM))                  ! Aut liquid 
               AUTICE(I,J,1:LM) = REAL(prcior8(1, 1:LM))                !Aut  ice    
               ACIL_AN_X(I,J,1:LM) = REAL(praior8(1, 1:LM))           !Acc ice  by snow
               ACLL_AN_X(I,J,1:LM) = REAL(msacwior8(1, 1:LM))   !HM process

               FRZPP_LS(I,J,1:LM) = REAL(frzrdtr8(1, 1:LM)) / MAPL_CP !precip freezing latent heat rate
               SNOWMELT_LS(I,J,1:LM) =REAL(meltsdtr8(1, 1:LM))/ MAPL_CP !melting of snow latent heat rate

               !diagnostics for number concentration budget (all grid-average, Kg-1 s-1)
              

               DNDCCN(I,J,1:LM)       = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               DNDACRLS(I,J,1:LM)     = REAL(npsacwsor8(1,1:LM)) !droplet number tendency from accretion by snow
               DNDACRLR(I,J,1:LM)     = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               DNDEVAPC(I,J,1:LM)     = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation 
               DNDAUTLIQ(I,J,1:LM)    = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion
               
               DNCACRIS (I,J,1:LM)     = REAL(npraior8(1,1:LM))  !ice number tendency from accretion by snow
               DNHET_CT(I,J,1:LM)      = REAL(nnucctor8(1,1:LM)) !ice number tendency from contact IN  
               DNHET_IMM(I,J,1:LM)     = REAL(nnucccor8(1,1:LM)) !ice number tendency from immersion IN    
               DNCNUC(I,J,1:LM)        = REAL(nnuccdor8(1,1:LM)) !ice number tendency from nucleation on aerosol             
               DNCSUBL (I,J,1:LM)      = REAL(nsubior8(1,1:LM))  !ice number tendency from sublimation               
               DNCAUTICE (I,J,1:LM)    = REAL(nprcior8(1,1:LM))  !ice number tendency from autoconversion
               DNCHMSPLIT(I,J,1:LM)    = REAL(nsacwior8(1,1:LM)) !ice number tendency from H-M process

               ! Total tendencies

               
               DQVDT_micro(I,J,1:LM)   = REAL(qvlatr8(1,1:LM))  
               DQIDT_micro(I,J,1:LM)   = REAL(qitendr8(1,1:LM))   
               DQLDT_micro(I,J,1:LM)   = REAL(qctendr8(1,1:LM) )    
               DTDT_micro(I,J,1:LM)    = REAL((tlatr8(1,1:LM) / MAPL_CP))
               DNDCCN(I,J,1:LM)       = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               DNDACRLS(I,J,1:LM)     = REAL(npsacwsor8(1,1:LM) )!droplet number tendency from accretion by snow
               DNDACRLR(I,J,1:LM)     = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               DNDEVAPC(I,J,1:LM)     = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation 
               DNDAUTLIQ(I,J,1:LM)    = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion

            enddo !I
         enddo !J
         !============================================Finish 2-moment micro implementation===========================

    !update water tracers
         QLCN=QL_TOT*FQA
         QLLS=QL_TOT-QLCN
         QICN=QI_TOT*FQA
         QILS=QI_TOT-QICN
         PFL_AN_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) * FQA
         PFL_LS_X(:,:,1:LM) = PFL_LS_X(:,:,1:LM) - PFL_AN_X(:,:,1:LM)
         PFI_AN_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) * FQA
         PFI_LS_X(:,:,1:LM) = PFI_LS_X(:,:,1:LM) - PFI_AN_X(:,:,1:LM)
         QTOT= QICN+QILS+QLCN+QLLS

       TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
               CN_SNR  + LS_SNR  + AN_SNR + SC_SNR

   
  
         !============ Put cloud fraction back in contact with the PDF and create new condensate if neccesary (Barahona et al., GMD, 2014)============

 DLPDF_X=  QLLS +QLCN
 DIPDF_X =  QILS +QICN

do K= 1, LM
   do J=1,JM
            do I=1,IM

                  call update_cld( &
                         DT_MOIST                , &
                         ALPHT_X(I, J, K)        , &
                         INT(CLDPARAMS%PDFSHAPE) , &
                         PLO(I, J, K)            , &
                         Q1 (I, J, K)            , &
                         QLLS(I, J, K)           , &
                         QLCN(I, J, K)           , &
                         QILS(I, J, K)           , &
                         QICN(I, J, K)           , &
                         TEMP(I, J, K)           , &
                         CLLS(I, J, K)           , &
                         CLCN(I, J, K)           , &
                         SC_ICE(I, J, K)         , &
                         NCPI(I, J, K)           , &
                         NCPL(I, J, K)           , &
                         INC_NUC(I, J, K)        , &
                         RHCmicro(I, J, K), &
             CNV_FRACTION(I, J), SNOMAS(I, J), FRLANDICE(I, J), FRLAND(I, J), .TRUE.)
              
           end do 
        end do
  end do 
  
  
   DLPDF_X=((QLLS+QLCN)  - DLPDF_X)/DT_MOIST
   DIPDF_X=((QILS+QICN)  - DIPDF_X)/DT_MOIST

         ! Make sure ice and liquid stay within T limits  

  call meltfrz_inst  (     &
              IM,JM,LM    , &
              TEMP              , &
              QLLS          , &
              QLCN         , &
              QILS           , &
              QICN          , &               
              NCPL         , &
              NCPI          )
                          
         RAD_CF =min(CLLS+CLCN, 1.0)

         TPREC = CN_PRC2 + LS_PRC2 + AN_PRC2 + SC_PRC2 + &
                 CN_SNR  + LS_SNR  + AN_SNR + SC_SNR + CNVPRCP

         !=============================================End Stratiform cloud processes==========================================
         !======================================================================================================================

         !Calculate CFICE and CFLIQ 

         CFLIQ=0.0
         CFICE=0.0
         QTOT= QICN+QILS+QLCN+QLLS
         QL_TOT = QLCN+QLLS
         QI_TOT = QICN+QILS

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
            CLDREFFI = MAPL_UNDEF
         end where

         where (QL_TOT .le. 0.0)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL = MAPL_UNDEF
         end where

         WHERE  (RAD_CF > 1e-4)
            RAD_QL = min((QLLS+QLCN)/RAD_CF, 1.0e-3)
            RAD_QI = min((QILS+QICN)/RAD_CF, 1.0e-3) !
            RAD_QR =  QRAIN/RAD_CF  
            RAD_QS =  QSNOW/RAD_CF
            RAD_QG =  QGRAUPEL/RAD_CF
         ELSEWHERE 
            RAD_QL = 0.0         
            RAD_QI = 0.0
            RAD_QR = 0.0
            RAD_QS = 0.0
            RAD_QG = 0.0
         end where

         
          
       
         !Everything in-cloud for radiation============== 

         RAD_QV = MAX( Q1 , 0. )
         RAD_QL = MAX(MIN( RAD_QL , 0.001 ), 0.0)  ! Still a ridiculously large
         RAD_QI = MAX(MIN( RAD_QI , 0.001 ), 0.0)  ! value.
         RAD_QR = MAX(MIN( RAD_QR , 0.01 ), 0.0)  ! value.
         RAD_QS = MAX(MIN( RAD_QS , 0.01 ), 0.0)  ! value
         RAD_QG = MAX(MIN( RAD_QG , 0.01 ), 0.0)  ! value
         
         !=================================================================================
         !    Units conversion for diagnostics

         CFX =100.*PLO*r_air/TEMP !density times conversion factor
 
  
         !to m-3
         NCPL_VOL=NCPL*CFX !
         NCPI_VOL=NCPI*CFX
         CDNC_NUC=CDNC_NUC*CFX 
         INC_NUC =INC_NUC*CFX             
         CNV_NICE=CNV_NICE*CFX 
         CNV_NDROP=CNV_NDROP*CFX         
         
        !to m-3 s-1

         DNHET_CT    = DNHET_CT*CFX 
         DNHET_IMM   = DNHET_IMM*CFX 
         DNCNUC      = DNCNUC*CFX 
         DNCHMSPLIT  = DNCHMSPLIT*CFX
         DNCSUBL     = DNCSUBL*CFX 
         DNCACRIS    = DNCACRIS*CFX 
         DNCAUTICE   = DNCAUTICE*CFX 
         DNCCNV      = DNCCNV_X*CFX

         DNDCCN       = DNDCCN*CFX   
         DNDACRLS     = DNDACRLS*CFX
         DNDACRLR     = DNDACRLR*CFX    
         DNDEVAPC     = DNDEVAPC*CFX   
         DNDAUTLIQ    = DNDAUTLIQ*CFX 
         DNDCNV       = DNDCNV_X*CFX

         !make grid averages

         !NHET_IMM=NHET_IMM*CFICE
         !INC_NUC = INC_NUC*CFICE
         !CDNC_NUC = CDNC_NUC*CFLIQ
         !NHET_NUC = NHET_NUC*CFICE
         !DUST_DEP = DUST_DEP*CFICE
         !DUST_IMM = DUST_IMM*CFICE

         !Grid average  volumetric  radius for comparison against field data
         WHERE  ((CFICE > 0.001) .and. (NCPI .gt. 1.0))
            RI_MASK =  CFICE*((QILS+QICN)/(800.0*1.333*MAPL_PI*NCPI))**0.333
         elsewhere
            RI_MASK=0.0
         end  where

         WHERE  ((CFLIQ > 0.001) .and. (NCPL .gt. 1.0)) 
            RL_MASK =  CFLIQ*((QLLS+QLCN)/(900.0*1.333*MAPL_PI*NCPL))**0.333
         ELSEWHERE
            RL_MASK = 0.0 
         END WHERE

         TH1 = TEMP / PK

  
         !Set rain water for radiation to 0 if preciprad flag is off (set to 0)
         if(CLDPARAMS%PRECIPRAD .eq. 0.) then
            RAD_QR = 0.
            RAD_QS = 0.
            RAD_QG = 0.      
         endif

         if (associated(QRTOT)) QRTOT = QRAIN
         if (associated(QSTOT)) QSTOT = QSNOW
         if (associated(QGTOT)) QGTOT = QGRAUPEL

         CLDREFFL = MAX(4.1e-6, CLDREFFL) !DONIF Limits according to MG2008-I 
         CLDREFFL = MIN(29.e-6, CLDREFFL)
         CLDREFFI = MAX(6.e-6, CLDREFFI)   
         CLDREFFI = MIN(89.e-6, CLDREFFI)  !maximum number for the correlation and modis sim 
  
         CLDREFFR = MAX(4.1e-6, CLDREFFR) 
         CLDREFFR = MIN(29.e-6, CLDREFFR)
         CLDREFFS = MAX(6.e-6, CLDREFFS)   
         CLDREFFS = MIN(89.e-6, CLDREFFS)  !maximum number for the correlation and modis sim   
         CLDREFFG = MAX(6.e-6, CLDREFFG)   
         CLDREFFG = MIN(89.e-6, CLDREFFG)  !maximum number for the correlation and modis sim 
          

         !===========================

         ! Diagnostic cloud top/base properties
         
         CLDREFFL_TOP_X = MAPL_UNDEF
         CLDREFFI_TOP_X =MAPL_UNDEF
         NCPL_TOP_X = MAPL_UNDEF
         NCPL_CLDBASEX = MAPL_UNDEF
         NCPI_TOP_X = MAPL_UNDEF
         kbmin = LM
         LWC_X = 0.0
         IWC_X = 0.0

         CFX =100.*PLO*r_air/TEMP
         DO I=1, IM
            DO J= 1 , JM

               cfaux(1, 1:LM) =  CFLIQ(I, J, 1:LM)
               call find_cldtop(1, LM, cfaux, kbmin)

               if (kbmin .ge. LM-1) then 
                  CLDREFFL_TOP_X (I, J)  = 8.0e-6 
                  NCPL_TOP_X (I, J)  = 0.0
               else

                  CLDREFFL_TOP_X (I, J)  = CLDREFFL(I, J,  kbmin) 
                  NCPL_TOP_X (I, J)  = NCPL_VOL(I, J,  kbmin) 
               end if

                call find_cldbase(1, LM, cfaux, kbmin)
                if (kbmin .gt. 10) then 
                  NCPL_CLDBASEX (I, J)  = NCPL_VOL(I, J,  kbmin)/max(cfaux(1, kbmin), 0.01) 
                end if
                
               cfaux(1, 1:LM) =CFICE(I, J, 1:LM)   
               call find_cldtop(1, LM, cfaux, kbmin)

               if (kbmin .ge. LM-1) then 
                  CLDREFFI_TOP_X (I, J)=20.0E-6
                  NCPI_TOP_X (I, J)  = 0.0
               else       
                  CLDREFFI_TOP_X (I, J)  = CLDREFFI(I, J,  kbmin) 
                  NCPI_TOP_X (I, J)  = NCPI_VOL(I, J,  kbmin)
               end if

            END DO
         END DO

    
      !=====Tune area cloud fraction of extended PBL clouds. Area cloud fraction may be different from Volume cloud fraction 
         FQA= 0.0
         where (RAD_CF .gt. 0.0)
              FQA =  CLCN/RAD_CF
         end where
         
       
 
         DO I=1, IM
            DO J = 1, JM    

               if (FRLAND(I, J) .lt. 0.1) then 

                  
                  cfc_aux =(TEMP(I, J, LM) - TMAXCFCORR)/2.0
                  cfc_aux =  min(max(cfc_aux,-20.0), 20.0)
                  cfc_aux=   1.0/(1.0+exp(-cfc_aux))
                  
                  DO K=LM-1, 2, -1
                     if ((RAD_CF(I, J, K) .gt. 0.01) .and. (RAD_CF(I, J, K) .lt. 0.99)) then  

                           USURF=1.0                        
                           USURF= (LTS_UP-LTS(I, J))/(LTS_UP -  CLDPARAMS%MIN_LTS) 
                           USURF=min(max(USURF, MIN_EXP), MAX_EXP)                            
                           
                           fracover=min(max((TEMP(I, J, K) -TMAXLL)/2.0, -20.0), 20.0)
                           fracover = 1.0/(1.0+exp(-fracover))
                          
                           USURF = USURF*fracover + 1.0-fracover   !only near the surface                                  
                           USURF = USURF*cfc_aux + 1.0-cfc_aux !only for the subtropics                          
                           USURF =  usurf+ (1.0-USURF)*CNV_FRACTION(I, J) !only non-convective
                                            
                           RAD_CF(I, J, K)=RAD_CF(I, J, K)**USURF                                                                
                     END IF

                  END DO
               END IF

            end do
         end do
         
    
           CLCN   =  FQA*RAD_CF
           CLLS =  (1.0-FQA)*RAD_CF   

       WHERE (QTOT .gt. 1.0e-12) 
            CFLIQ=RAD_CF*QL_TOT/QTOT
            CFICE=RAD_CF*QI_TOT/QTOT
        END WHERE
         
      if(associated(CNV_FRC )) CNV_FRC  = CNV_FRACTION

      call MAPL_TimerOff (MAPL,"---CLDMICRO", RC=STATUS)

         !=======================================================================
#endif

   call MAPL_TimerOff(MAPL,"--MGB2_2M",RC=STATUS); VERIFY_(STATUS)

end subroutine MGB2_2M_Run

end module GEOS_MGB2_2M_InterfaceMod
