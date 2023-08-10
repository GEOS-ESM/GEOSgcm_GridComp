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
  integer :: imsize
  real    :: TURNRHCRIT
  real    :: MINRHCRITLND
  real    :: MINRHCRITOCN
  real    :: MAXRHCRITLND
  real    :: MAXRHCRITOCN
  real    :: CCW_EVAP_EFF
  real    :: CCI_EVAP_EFF
  integer :: PDFSHAPE
  real    :: FAC_RL
  real    :: MIN_RL
  real    :: MAX_RL
  real    :: FAC_RI
  real    :: MIN_RI
  real    :: MAX_RI
  logical :: LHYDROSTATIC
  logical :: LPHYS_HYDROSTATIC


  real  :: DCS, QCVAR_, WBFFACTOR, NC_CST, NI_CST, NG_CST, MUI_CST, PMIN_CBL
  real  :: LCCIRRUS, UISCALE, SS_SCALE, REEVAP_MICRO, LIU_MU, TFRZ, &
           NPRE_FRAC, QCVAR, ZPBLMAXLL, TMAXLL, LTS_LOW, LTS_UP, MIN_EXP,     &
           BKGTAU, DCRIT_, USE_AV_V, AUTSC, TS_AUTO_ICE, CCN_PARAM, IN_PARAM, &
           FDROP_DUST, FDROP_SOOT, USE_WSUB_CLIM, SIGMA_NUC, MIN_ALH, &
           HMOIST_950, HSMOIST_500, SINST, MAX_EXP, MAX_CAPE, MIN_CAPE,       &
           DUST_INFAC, ORG_INFAC, BC_INFAC, SS_INFAC, RRTMG_IRRAD, RRTMG_SORAD,&
           SCWST, MTIME, SWCIRRUS, MINCDNC, TMAXCFCORR,    &
           Immersion_param, ACC_ENH, ACC_ENH_ICE, DT_MICRO, DT_AUX, UR_SCALE, &
           CNV_NUMLIQ_SC, CNV_NUMICE_SC


  public :: MGB2_2M_Setup, MGB2_2M_Initialize, MGB2_2M_Run
  public :: MGVERSION
  character(LEN=ESMF_MAXSTR):: CONVPAR_OPTION

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

    call ESMF_ConfigGetAttribute( CF, MGVERSION, Label="MGVERSION:",  default=1, __RC__)
    call ESMF_ConfigGetAttribute( CF, CONVPAR_OPTION, Label='CONVPAR_OPTION:', __RC__) ! Note: Default set in GEOS_GcmGridComp.F90


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

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
    real, pointer, dimension(:,:,:)     :: NCPL, NCPI, NRAIN, NSNOW, NGRAUPEL

    logical  :: nccons, nicons, ngcons, do_graupel
    real(ESMF_KIND_R8)  Dcsr8, qcvarr8,  micro_mg_berg_eff_factor_in, ncnstr8, ninstr8, ngnstr8, mui_cnstr8



    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: nn
    real                       :: tmprhL, tmprhO


    IAm = "MGB2_2M_Initialize"

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM.GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LPHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM.GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4

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
    call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'      , __RC__)
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'      , __RC__)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL,  'NGRAUPEL'      , __RC__)

    call WRITE_PARALLEL ("INITIALIZED MGB2_2M microphysics in non-generic GC INIT")

    call MAPL_GetResource( MAPL, PDFSHAPE        , 'PDFSHAPE:'        , DEFAULT= 2     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RI          , 'FAC_RI:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RI          , 'MIN_RI:'          , DEFAULT=  5.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RI          , 'MAX_RI:'          , DEFAULT=140.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RL          , 'FAC_RL:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RL          , 'MIN_RL:'          , DEFAULT= 2.5e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RL          , 'MAX_RL:'          , DEFAULT=60.0e-6, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, CCW_EVAP_EFF, 'CCW_EVAP_EFF:', DEFAULT= 4.e-3, RC=STATUS); VERIFY_(STATUS)
	call MAPL_GetResource( MAPL, CCI_EVAP_EFF, 'CCI_EVAP_EFF:', DEFAULT= 4.e-3, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=    0.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_EXP, 'CNV_FRACTION_EXP:', DEFAULT=    0.5, RC=STATUS); VERIFY_(STATUS)


    !2M========


    call MAPL_GetResource(MAPL, LCCIRRUS,       'LCCIRRUS:',       DEFAULT= 500.0,  __RC__) !Characteristic Length (m) of high freq gravity waves
    call MAPL_GetResource(MAPL, UISCALE,        'UISCALE:',        DEFAULT= 1.0,    __RC__) !Scaling factor for sed vel of ice
    call MAPL_GetResource(MAPL, LIU_MU,         'LIU_MU:',         DEFAULT= 2.0,    __RC__) !Liu autoconversion parameter
    call MAPL_GetResource(MAPL, NPRE_FRAC,      'NPRE_FRAC:',      DEFAULT= -1.0,   __RC__) !Fraction of preexisting ice affecting ice nucleationn
    call MAPL_GetResource(MAPL, LTS_LOW,        'LTS_LOW:',        DEFAULT= 20.0,   __RC__) !lower LTS for morphology correction
    call MAPL_GetResource(MAPL, LTS_UP,         'LTS_UP:',         DEFAULT= 22.0,   __RC__) !Upper LTS for morphology correction
    call MAPL_GetResource(MAPL, MIN_EXP,        'MIN_EXP:',        DEFAULT= 0.5,    __RC__) !Exponent of the relation CFA=CFV^n
    call MAPL_GetResource(MAPL, MAX_EXP,        'MAX_EXP:',        DEFAULT= 1.0,    __RC__) !Exponent of the relation CFA=CFV^n
    call MAPL_GetResource(MAPL, USE_AV_V,       'USE_AV_V:',       DEFAULT= 1.0,    __RC__) !Set to > 0 to use an average velocity for activation
    call MAPL_GetResource(MAPL, AUTSC,          'AUT_SCALE:',      DEFAULT= 1.0,    __RC__) !scale factor for critical size for drizzle
    call MAPL_GetResource(MAPL, TS_AUTO_ICE,    'TS_AUTO_ICE:',    DEFAULT= 360.,    __RC__) !Ice autoconversion time scale
    call MAPL_GetResource(MAPL, TMAXLL,         'TMAXLL:',         DEFAULT= 250.0,  __RC__) !Liquid clouds min T
    call MAPL_GetResource(MAPL, CCN_PARAM,      'CCNPARAM:',       DEFAULT= 2.0,    __RC__) !CCN activation param
    call MAPL_GetResource(MAPL, IN_PARAM,       'INPARAM:',        DEFAULT= 6.0,    __RC__) !IN param
    call MAPL_GetResource(MAPL, Immersion_param,'ImmersionPARAM:', DEFAULT= 6.0,    __RC__) !Immersion param
    call MAPL_GetResource(MAPL, ACC_ENH,        'ACC_ENH:',        DEFAULT= 1.0,    __RC__) !accretion rain-liquid scaling for MG2
    call MAPL_GetResource(MAPL, ACC_ENH_ICE,    'ACC_ENH_ICE:',    DEFAULT= 0.05,    __RC__) !accretion snow-ice scaling for MG2
    call MAPL_GetResource(MAPL, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.5,    __RC__) !Fraction of dust within droplets for immersion freezing
    call MAPL_GetResource(MAPL, FDROP_SOOT,     'FDROP_SOOT:',     DEFAULT= 0.05,   __RC__) !Fraction of soot within droplets for immersion freezing
    call MAPL_GetResource(MAPL, SIGMA_NUC,      'SIGMA_NUC:',      DEFAULT= 1.0,    __RC__) !Widht of the in-cloud distribution of relative humidity in cirrus
    call MAPL_GetResource(MAPL, MIN_ALH,        'MIN_ALH:',        DEFAULT= 5.0,    __RC__) !scale factor for vertical velocity in sttratocumulus
    call MAPL_GetResource(MAPL, SCWST,          'SCWST:',          DEFAULT= 3.0,    __RC__) !scale factor for vertical velocity in sttratocumulus
    call MAPL_GetResource(MAPL, MINCDNC,        'MINCDNC:',        DEFAULT= 25.0,    __RC__) !min nucleated droplet conc. cm-3
    call MAPL_GetResource(MAPL, TMAXCFCORR,     'TMAXCFCORR:',     DEFAULT= 285.0,  __RC__) !Minimum T for CF correction
    call MAPL_GetResource(MAPL, MTIME,          'MTIME:',          DEFAULT= -1.0,   __RC__) !Mixing time scale for aerosol within the cloud. Default is time step
    call MAPL_GetResource(MAPL, SWCIRRUS,       'SWCIRRUS:',       DEFAULT= 3.0,    __RC__) !Tunes vertical velocity in cirrus
    call MAPL_GetResource(MAPL, DUST_INFAC,     'DUST_INFAC:',     DEFAULT= 1.0,    __RC__)  !work on this
    call MAPL_GetResource(MAPL, BC_INFAC,       'BC_INFAC:',       DEFAULT= 0.1,    __RC__)
    call MAPL_GetResource(MAPL, ORG_INFAC,      'ORG_INFAC:',      DEFAULT= 1.0,    __RC__)
    call MAPL_GetResource(MAPL, SS_INFAC,       'SS_INFAC:',       DEFAULT= 1.0,    __RC__)
    call MAPL_GetResource(MAPL, DT_MICRO,       'DT_MICRO:',       DEFAULT= 300.0,  __RC__)    ! time step of the microphysics substepping (s) (MG2) (5 min)
    call MAPL_GetResource(MAPL, UR_SCALE,       'URSCALE:',        DEFAULT= 1.0,    __RC__) !Scaling factor for sed vel of rain
    call MAPL_GetResource(MAPL, USE_WSUB_CLIM,  'USE_WSUB_CLIM:',   DEFAULT= 1.0,    __RC__) !Use Wsub climatology
    call MAPL_GetResource( MAPL, RRTMG_IRRAD ,  'USE_RRTMG_IRRAD:',DEFAULT=1.0,     __RC__)
    call MAPL_GetResource( MAPL, RRTMG_SORAD ,  'USE_RRTMG_SORAD:',DEFAULT=1.0,     __RC__)
    call MAPL_GetResource(MAPL, CNV_NUMLIQ_SC,   'CNV_NUMLIQ_SC:', DEFAULT= 0.02 ,RC=STATUS) !scaling for conv number
    call MAPL_GetResource(MAPL, CNV_NUMICE_SC,   'CNV_NUMICE_SC:', DEFAULT= 1.0 ,RC=STATUS)
    call MAPL_GetResource(MAPL, DCS,      'DCS:'    , DEFAULT=250.0e-6, __RC__ )
    Dcsr8 = DCS
    call MAPL_GetResource(MAPL, QCVAR_,   'QCVAR:'  , DEFAULT= 2.0 ,__RC__) !variance of the QL distribution

    qcvarr8=QCVAR_
    call MAPL_GetResource(MAPL, WBFFACTOR,   'WBFFACTOR:', DEFAULT= 0.05 ,__RC__) !variance of the QL distribution

    micro_mg_berg_eff_factor_in = WBFFACTOR
    call MAPL_GetResource(MAPL, NC_CST ,  'NC_CST:' , DEFAULT=  0.0 ,__RC__) !constant nd (set if greather than zero)

    call MAPL_GetResource(MAPL, NI_CST ,  'NI_CST:' , DEFAULT=  0.0 ,__RC__) !constant nd (set if greather than zero)

    call MAPL_GetResource(MAPL, NG_CST ,  'NG_CST:' , DEFAULT=  0.0 ,__RC__) !constant ng (set if greather than zero)

    call MAPL_GetResource(MAPL, MUI_CST,  'MUI_CST:', DEFAULT= -1.0 ,__RC__) !constant ng (set if greather than zero)


    mui_cnstr8 =  MUI_CST
    ncnstr8 = NC_CST
    if  (NC_CST .gt. 0.0)  nccons =.true.
    ninstr8 = NI_CST
    if  (NI_CST .gt. 0.0)  nicons =.true.
    ngnstr8 = NG_CST
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
    real, pointer, dimension(:,:,:) :: ZLE, PLE, PK, T, U, V, W, KH, TKE
    real, pointer, dimension(:,:)   :: AREA, FRLAND, TS, DTSX, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:) :: HL2, HL3, QT2, QT3, W2, W3, HLQT, WQT, WQL, WHL, EDMF_FRC
    real, pointer, dimension(:,:,:) :: WTHV2
    real, pointer, dimension(:,:,:) :: OMEGA

    real, pointer, dimension(:,:)   :: TAUOROX, TAUOROY
    real, pointer, dimension(:,:,:) :: ALH, RADLW, RADSW, WSUB_CLIM

    ! Local
    real, allocatable, dimension(:,:,:) :: U0, V0
    real, allocatable, dimension(:,:,:) :: PLEmb, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  ZL0, GZLO, PKmb
    real, allocatable, dimension(:,:,:) :: DZ, DZET, DP, MASS, iMASS
    real, allocatable, dimension(:,:,:) :: DQST3, QST3
    real, allocatable, dimension(:,:,:) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                                           DQSDTmic, DQGDTmic, DQADTmic, &
                                           DUDTmic,  DVDTmic,  DTDTmic
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: IKEX, IKEX2
    real, allocatable, dimension(:,:)   :: frland2D
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
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D


    !2m
    real, pointer, dimension(:,:,:) :: SC_ICE, CDNC_NUC, INC_NUC, PFRZ, &
       CFICE, CFLIQ, DT_RASP, SMAXL, SMAXI, WSUB, CCN01, CCN04, CCN1, &
       NHET_NUC, NLIM_NUC, SO4, ORG, BCARBON, DUST, SEASALT, NCPL_VOL, NCPI_VOL, &
       SAT_RAT, RHICE, RL_MASK, RI_MASK, &
       NHET_IMM, NHET_DEP, DUST_IMM, DUST_DEP, SIGW_GW, SIGW_CNV, SIGW_TURB, &
       SIGW_RC, BERG, BERGS, MELT, DNHET_CT, QCRES, QIRES, AUTICE, FRZPP_LS, &
       SNOWMELT_LS, DNCNUC, DNCSUBL, DNCHMSPLIT, DNCAUTICE, DNCACRIS, DNDCCN, &
       DNDACRLS, DNDACRLR, DNDEVAPC, DNDAUTLIQ, DNDCNV, DNICNV, &
       CNV_UPDF, CNV_CVW, DNHET_IMM, CNV_MFD, CNV_DQCDT, KAPPA, RHCmicro, RHLIQ, &
       CNV_NICE, CNV_NDROP, NWFA, CNV_FICE

     real, pointer, dimension(:,:)   :: EIS, LTS, QCVAR_EXP, &
       CCNCOLUMN, NDCOLUMN, NCCOLUMN



    real, allocatable, dimension(:,:,:) :: dNI, dNL, QCNTOT, CFX,  QTOT, &
       QL_TOT, QI_TOT, ACIL_LS_X, ACIL_AN_X, ACLL_LS_X, ACLL_AN_X, DLPDF_X, DIPDF_X, DLFIX_X, DIFIX_X, &
       AUT_X, SDM_X, FRZ_TT_X, FRZ_PP_X, DCNVL_X, DCNVI_X, AIRDEN, TH1, FQA, ALPH3D  !check how much of these we are actually using

    integer, allocatable, dimension(:, :)   ::  KMIN_TROP, KLCL
    real, allocatable, dimension(:, :)  :: NPRE_FRAC_2d, CLDREFFI_TOP_X, CLDREFFL_TOP_X,  NCPL_TOP_X, NCPI_TOP_X, NCPL_CLDBASEX, ZWS, ZPBL

    ! Local variables
    real    :: ALPHA, RHCRIT
    integer :: IM,JM,LM
    integer :: I, J, L, K
    real :: dw_land = 0.20 !< base value for subgrid deviation / variability over land
    real :: dw_ocean = 0.10 !< base value for ocean


    integer :: num_steps_micro,  pcnst, n_modes, kbmin, kcldtop, kcldbot , &
                NAUX, kcldtopcvn, nbincontactdust, index, K0, KCBLMIN, i_src_mode, i_dst_mode

    real, parameter :: pmin_trop = 10.0 !mbar minimum pressure to do cloud microphysics
    logical                   :: use_average_v
    REAL, allocatable, dimension(:,:) :: SCICE_tmp, FQA_tmp,   tm_gw, pm_gw, nm_gw, theta_tr,  &
            fcn, cfaux, pi_gw, rhoi_gw, ni_gw, ti_gw, h_gw, Wbreak

    real (ESMF_KIND_R8), dimension(3)       :: ccn_diag
    real(ESMF_KIND_R8), allocatable, dimension(:,:,:) :: rndstr8,naconr8  !Assume maximum 5 dust bins
    real(ESMF_KIND_R8), dimension(1)       :: prectr8, precir8
    real (ESMF_KIND_R8)  :: tauxr8, fsoot_drop, fdust_drop, sigma_nuc_r8, rh1_r8, &
                            frachet_dust, frachet_bc, frachet_org, frachet_ss, &
                            disp_liu, ui_scale, dcrit, tfreez, qcvar8, &
                            ts_autice, dcsr8, qcvarr8, scale_ri, mtimesc, urscale


    real(ESMF_KIND_R8), allocatable, dimension(:,:)  :: ttendr8, qtendr8, cwtendr8, &
           cldor8,  rpdelr8, zmr8, omegr8, rhdfdar8, rhu00r8, ficer8 , &
           ndropr8, nimmr8, wparc, smaxliq, atot,  smaxicer8, nheticer8, incr8, swparc, &
           nhetr8, nlimicer8, qilsr8, wparc_gw, wparc_ls, wparc_turb, wparc_cnv, lc_turb, rad_cooling, wparc_rc, &
            uwind_gw, wparc_cgw, pfrz_inc_r8, pintr8, kkvhr8, rflxr8,  sflxr8, lflxr8, iflxr8, gflxr8,  &
            so4x, seasaltx, dustx, &
            orgx, bcx, ter8,qvr8, qcr8,qir8, ncr8,nir8, qrr8,qsr8, nrr8,nsr8, &
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
            nsacwior8, mnuccror8,pracsor8, qiresor8, rate1ord_cw2pr, & !only MG1
            sc_icer8, nhet_immr8, dnhet_immr8, nhet_depr8, & ! activation
            dust_immr8, dust_depr8,dpre8, npre8, accre_enhan_icer8

      real   :: maxkhpbl, tausurf_gw, fracover, cfc_aux, aux1,aux2,aux3,hfs,hfl, Nct, Wct, ksa1, Xscale


      real(ESMF_KIND_R8)   :: autscx
      real, parameter :: r_air = 3.47d-3 !m3 Pa kg-1K-1
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
    allocate(ndropr8(1,LM), __STAT__)
    allocate(nimmr8(1,LM), __STAT__)
    allocate(wparc(1,LM), __STAT__)
    allocate(smaxliq(1,LM), __STAT__)
    allocate(atot(1,LM), __STAT__)
    allocate(smaxicer8(1,LM), __STAT__)
    allocate(nheticer8(1,LM), __STAT__)
    allocate(incr8(1,LM), __STAT__)
    allocate(swparc(1,LM), __STAT__)
    allocate(nhetr8(1,LM), __STAT__)
    allocate(nlimicer8(1,LM), __STAT__)
    allocate(qilsr8(1,LM), __STAT__)
    allocate(wparc_gw(1,LM), __STAT__)
    allocate(wparc_ls(1,LM), __STAT__)
    allocate(wparc_turb(1,LM), __STAT__)
    allocate(wparc_cnv(1,LM), __STAT__)
    allocate(lc_turb(1,LM), __STAT__)
    allocate(rad_cooling(1,LM), __STAT__)
    allocate(wparc_rc(1,LM), __STAT__)
    allocate(uwind_gw(1,LM), __STAT__)
    allocate(wparc_cgw(1,LM), __STAT__)
    allocate(pfrz_inc_r8(1,LM), __STAT__)
    allocate(SCICE_tmp(1,LM), __STAT__)
    allocate(FQA_tmp(1,LM), __STAT__)

    allocate(so4x(1,LM), __STAT__)
    allocate(seasaltx(1,LM), __STAT__)
    allocate(dustx(1,LM), __STAT__)
    allocate(orgx(1,LM), __STAT__)
    allocate(bcx(1,LM), __STAT__)
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
    allocate(sc_icer8(1,LM), __STAT__)
    allocate(nhet_immr8(1,LM), __STAT__)
    allocate(dnhet_immr8(1,LM), __STAT__)
    allocate(nhet_depr8(1,LM), __STAT__)
    allocate(dust_immr8(1,LM), __STAT__)
    allocate(dust_depr8(1,LM), __STAT__)
    allocate(accre_enhan_icer8(1,LM), __STAT__)
    allocate(dpre8(1,LM), __STAT__)
    allocate(npre8(1,LM), __STAT__)
    allocate(pintr8(1,LM+1), __STAT__)
    allocate(kkvhr8(1,LM+1), __STAT__)
    allocate(rflxr8(1,LM+1), __STAT__)
    allocate(sflxr8(1,LM+1), __STAT__)
    allocate(lflxr8(1,LM+1), __STAT__)
    allocate(iflxr8(1,LM+1), __STAT__)
    allocate(gflxr8(1,LM+1), __STAT__)
    allocate(rndstr8(1,LM,10), __STAT__)
    allocate(naconr8(1,LM,10), __STAT__)
    allocate(tm_gw(1,LM), __STAT__)
    allocate(pm_gw(1,LM), __STAT__)
    allocate(nm_gw(1,LM), __STAT__)
    allocate(theta_tr(1,LM), __STAT__)
    allocate(fcn(1,LM), __STAT__)
    allocate(cfaux(1,LM), __STAT__)
    allocate(pi_gw(1,0:LM), __STAT__)
    allocate(rhoi_gw(1,0:LM), __STAT__)
    allocate(ni_gw(1,0:LM), __STAT__)
    allocate(ti_gw(1,0:LM), __STAT__)
    allocate(h_gw(1,0:LM), __STAT__)

    allocate(KMIN_TROP(IM,JM), __STAT__)
    allocate(NPRE_FRAC_2d(IM,JM), __STAT__)
    allocate(ZWS(IM,JM), __STAT__)
    allocate(ZPBL(IM,JM), __STAT__)

    allocate(FQA(IM,JM,LM ), __STAT__)
    allocate(ALPH3D(IM,JM,LM ), __STAT__)
    allocate(GZLO(IM,JM,LM ), __STAT__)
    allocate(TH1(IM,JM,LM ), __STAT__)
    allocate(PK(IM,JM,LM ), __STAT__)
    allocate(QCNTOT(IM,JM,LM), __STAT__)
    allocate(CFX(IM,JM,LM), __STAT__)
    allocate(AIRDEN(IM,JM,LM), __STAT__)

    allocate(QTOT(IM,JM,LM ), __STAT__)
    allocate(QL_TOT(IM,JM,LM ), __STAT__)
    allocate(QI_TOT(IM,JM,LM ), __STAT__)
    allocate(ACIL_AN_X(IM,JM,LM ), __STAT__)
    allocate(ACIL_LS_X(IM,JM,LM ), __STAT__)
    allocate(ACLL_AN_X(IM,JM,LM ), __STAT__)
    allocate(ACLL_LS_X(IM,JM,LM ), __STAT__)
    allocate(DLPDF_X(IM,JM,LM ), __STAT__)
    allocate(DIPDF_X(IM,JM,LM ), __STAT__)
    allocate(DLFIX_X(IM,JM,LM ), __STAT__)
    allocate(DIFIX_X(IM,JM,LM ), __STAT__)
    allocate(AUT_X(IM,JM,LM ), __STAT__)
    allocate(SDM_X(IM,JM,LM ), __STAT__)
    allocate(FRZ_TT_X(IM,JM,LM ), __STAT__)
    allocate(FRZ_PP_X(IM,JM,LM ), __STAT__)
    allocate(DCNVL_X(IM,JM,LM ), __STAT__)
    allocate(DCNVI_X(IM,JM,LM ), __STAT__)
    allocate(CLDREFFI_TOP_X(IM,JM ), __STAT__)
    allocate(CLDREFFL_TOP_X(IM,JM ), __STAT__)
    allocate(NCPL_TOP_X(IM,JM ), __STAT__)
    allocate(NCPI_TOP_X(IM,JM ), __STAT__)
    allocate(NCPL_CLDBASEX(IM,JM ), __STAT__)
    !allocate(TH(IM,JM,LM ), __STAT__)


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
    call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL, 'NGRAUPEL'    , __RC__)


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
    call MAPL_GetPointer(IMPORT, EDMF_FRC,'EDMF_FRC', RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, W2,      'W2'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, W3,      'W3'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WQT,     'WQT'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, WHL,     'WHL'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, HL2,     'HL2'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, HL3,     'HL3'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QT2,     'QT2'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, QT3,     'QT3'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, HLQT,    'HLQT'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TS,      'TS'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBLSC,  'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SH,      'SH'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP,    'EVAP'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   , RC=STATUS); VERIFY_(STATUS)

    !call MAPL_GetPointer(IMPORT, KPBLIN, 'KPBL'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUOROX, 'TAUOROX'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUOROY, 'TAUOROY'     , __RC__)
    call MAPL_GetPointer(IMPORT, ALH,    'ALH'     , __RC__)
    call MAPL_GetPointer(IMPORT, RADLW,  'RADLW'     , __RC__)
    call MAPL_GetPointer(IMPORT, RADSW,  'RADSW'     , __RC__)
    call MAPL_GetPointer(IMPORT, WSUB_CLIM,  'WSUB_CLIM'     , __RC__)
    call MAPL_GetPointer(IMPORT, TKE,    'TKE'       , __RC__)

    call MAPL_GetPointer(EXPORT, CFICE,   'CFICE'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CFLIQ,   'CFLIQ'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_FICE,   'CNV_FICE'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CCNCOLUMN,   'CCNCOLUMN'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NDCOLUMN,   'NDCOLUMN'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NCCOLUMN,   'NCCOLUMN'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RHLIQ,   'RHLIQ'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RHCmicro,   'RHCmicro'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QCVAR_EXP,   'QCVAR_EXP'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SC_ICE,      'SC_ICE'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFR,    'RR'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFS,    'RS'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFG,    'RG'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CDNC_NUC,    'CDNC_NUC'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, INC_NUC,     'INC_NUC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, PFRZ,        'PFRZ'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, LTS,         'LTS'         , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, EIS,         'EIS'         , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SMAXL,       'SMAX_LIQ'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SMAXI,       'SMAX_ICE'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, WSUB,        'WSUB'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CCN01,       'CCN01'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CCN04,       'CCN04'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CCN1,        'CCN1'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NHET_NUC,    'NHET_NUC'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NLIM_NUC,    'NLIM_NUC'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SO4,         'SO4'         , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, ORG,         'ORG'         , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, BCARBON,     'BCARBON'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DUST,        'DUST'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SEASALT,     'SEASALT'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NCPL_VOL,    'NCPL_VOL'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NCPI_VOL,    'NCPI_VOL'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SAT_RAT,     'SAT_RAT'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RHICE,       'RHICE'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RL_MASK,     'RL_MASK'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RI_MASK,     'RI_MASK'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NHET_IMM,    'NHET_IMM'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NHET_DEP,    'NHET_DEP'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DUST_IMM,    'DUST_IMM'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DUST_DEP,    'DUST_DEP'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SIGW_GW,     'SIGW_GW'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SIGW_CNV,    'SIGW_CNV'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SIGW_TURB,   'SIGW_TURB'   , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SIGW_RC,     'SIGW_RC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, BERG,        'BERG'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, BERGS,       'BERGS'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, MELT,        'MELT'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNHET_CT,    'DNHET_CT'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QCRES,       'QCRES'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QIRES,       'QIRES'       , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, AUTICE,      'AUTICE'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, FRZPP_LS ,   'FRZPP_LS'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SNOWMELT_LS, 'SNOWMELT_LS' , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNCNUC,      'DNCNUC'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNCSUBL,     'DNCSUBL'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNCHMSPLIT,  'DNCHMSPLIT'  , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNCAUTICE,   'DNCAUTICE'   , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNCACRIS,    'DNCACRIS'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNDCCN,      'DNDCCN'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNDACRLS,    'DNDACRLS'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNDACRLR,    'DNDACRLR'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNDEVAPC,    'DNDEVAPC'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNDAUTLIQ,   'DNDAUTLIQ'   , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNDCNV,      'DNDCNV'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNICNV,      'DNICNV'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNHET_IMM,   'DNHET_IMM'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, KAPPA,   'KAPPA'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

  ! This export MUST have been filled in the GridComp
    call MAPL_GetPointer(EXPORT, CNV_FRC,      'CNV_FRC'      , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SRF_TYPE,     'SRF_TYPE'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)


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
    ALLOCATE ( TMP3D(IM,JM,LM  ) )

     ! 2D Variables
    ALLOCATE ( IKEX         (IM,JM) )
    ALLOCATE ( IKEX2        (IM,JM) )
    ALLOCATE ( frland2D     (IM,JM) )
    ALLOCATE ( KLCL         (IM,JM) )
    ALLOCATE ( TMP2D        (IM,JM) )

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
    PKmb       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
    TH1       = T/PKmb
    AIRDEN = 100.*PLmb/T/MAPL_RGAS
    GZLO = MAPL_GRAV*ZL0

    ! Lowe tropospheric stability and estimated inversion strength
    call MAPL_GetPointer(EXPORT, LTS,   'LTS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EIS,   'EIS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    KLCL = FIND_KLCL( T, Q, PLmb, IM, JM, LM )
    TMP3D = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
    call FIND_EIS(TH1, QST3, T, ZL0, PLEmb, KLCL, IM, JM, LM, LTS, EIS)
    call find_l(KMIN_TROP, PLmb, pmin_trop, IM, JM, LM, 10, LM-2)

!=======================================================================================================================
!=======================================================================================================================
!===================================Nucleation of cloud droplets and ice crystals ======================================
! Aerosol cloud interactions. Calculate maxCCN tendency using Fountoukis and nenes (2005) or Abdul Razzak and Ghan (2002)
! liquid Activation Parameterization
! Ice activation follows the Barahona & Nenes ice activation scheme, ACP, (2008, 2009).
! Written by Donifan Barahona and described in Barahona et al. (2013, 2017, 2023)
!=======================================================================================================================
!=======================================================================================================================
!=======================================================================================================================

      call MAPL_TimerOn(MAPL,"---ACTIV") !Activation timer
      !!!! Include Deep Cnv tendencies for number concentrations
      call MAPL_GetPointer(EXPORT, CNV_NICE,  'CNV_NICE',  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_NDROP, 'CNV_NDROP', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, NWFA,      'NWFA',      ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

      allocate( dNI(IM,JM,LM) )
      allocate( dNL(IM,JM,LM) )
      dNI = 0.0
      dNL = 0.0
      DNDCNV = 0.0
      DNICNV = 0.0

      DO I=  1, IM
          DO J =  1, JM

             KCBLMIN = MAX(NINT(KPBLSC(I, J)), NINT(LM*0.9))
             DO K =  1, LM
                          CFX(I, J, K) = NWFA(I, J, K)*PLmb(I, J, K)/PLmb(I, J, KCBLMIN)
             end DO
          end do
        end do

      DNDCNV =  dNL
      DNICNV =  dNI


       ! CNV_MFD includes Deep+Shallow mass flux
      call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFD', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) then
       dNl =  CNV_NUMLIQ_SC*CFX*PTR3D/MASS
       NCPL = NCPL + dNL*DT_MOIST
      endif


      call MAPL_GetPointer(EXPORT, PTR3D, 'DQIDT_DC', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) then
         dNI  = make_IceNumber (PTR3D, T)*CNV_NUMICE_SC
         NCPI = NCPI + dNI*DT_MOIST
      endif


    !  call MAPL_GetPointer(EXPORT, PTR3D, 'DQLDT_DC', RC=STATUS); VERIFY_(STATUS)
    !  if (associated(PTR3D)) then

    !         dNL  = make_DropletNumber (PTR3D, NWFA)*CNV_NUMLIQ_SC
    !          NCPL = NCPL + dNL*DT_MOIST
    ! endif


    ! Include Shallow Cnv tendencies for number concentrations
      call MAPL_GetPointer(EXPORT, PTR3D, 'QIDET_SC', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) then
         dNI  = make_IceNumber (PTR3D, T)*CNV_NUMICE_SC
         NCPI = NCPI + dNI*DT_MOIST
      endif
    !  call MAPL_GetPointer(EXPORT, PTR3D, 'QLDET_SC', RC=STATUS); VERIFY_(STATUS)
      ! if (associated(PTR3D)) then
      !    dNL  = make_DropletNumber (PTR3D, NWFA)*CNV_NUMLIQ_SC
       !    NCPL = NCPL + dNL*DT_MOIST
      ! endif

      DNDCNV =  dNL + DNDCNV
      DNICNV =  dNI + DNICNV

       ! CNV_MFD includes Deep+Shallow mass flux
      call MAPL_GetPointer(EXPORT, PTR3D, 'CNV_MFD', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR3D)) then
        where (PTR3D .gt. 0.)
         CNV_NICE  =  CNV_NICE + dNI*MASS/PTR3D
         CNV_NDROP =  CNV_NDROP + dNL*MASS/PTR3D
        else where
         CNV_NICE  = 0.0
         CNV_NDROP = 0.0
        end where
      endif

      deallocate( dNI )
      deallocate( dNL )

    !================  Stratiform activation ===========================================

     if (NPRE_FRAC > 0.0) then
         NPRE_FRAC_2d = NPRE_FRAC
     else
         ! include CNV_FRC dependence
         DO J=1, JM
            DO I=1, IM
            NPRE_FRAC_2d(I,J) = CNV_FRC(I,J)*ABS(NPRE_FRAC) + (1-CNV_FRC(I,J))*0.05
            END DO
         END DO
     endif

       use_average_v = .false.
       if (USE_AV_V .gt. 0.0) then
         use_average_v = .true.
       end if
        fdust_drop   =  FDROP_DUST
        fsoot_drop   =  FDROP_SOOT
        sigma_nuc_r8 =  SIGMA_NUC
        frachet_org  =  ORG_INFAC
        frachet_dust =  DUST_INFAC
        frachet_bc   =  BC_INFAC
        frachet_ss   =  SS_INFAC

         if (USE_WSUB_CLIM .gt. 0.) then
            xscale = 8.7475*(imsize**(-0.328)) ! scale for resolutions =! 50 km
         end if
        !Supersaturations to calculate CCN diagnostics
        ccn_diag(1)=0.001
        ccn_diag(2)=0.004
        ccn_diag(3)=0.01



         do J=1,JM
            do I=1,IM

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
                     wparc_ls = 0.0
                     wparc_gw = 0.0
                     wparc_cgw= 0.0
                     wparc_turb = 0.0
                     swparc=0.0
                     pfrz_inc_r8 = 0.0
                     omegr8(1,1:LM) = OMEGA(I,J,1:LM)
                     kbmin= min(NINT(KPBLSC(I, J)), LM-1)-2
                     rad_cooling(1,1:LM) = RADLW(I,J,1:LM)+RADSW(I,J,1:LM)
                     wparc_ls(1,1:LM) =-OMEGA(I,J,1:LM)/AIRDEN(I,J,1:LM)/MAPL_GRAV + MAPL_CP*rad_cooling(1,1:LM)/MAPL_GRAV

                     !!=============== find vertical velocity variance

                     if (USE_WSUB_CLIM .le. 0.) then

                         uwind_gw(1,1:LM)           = min(0.5*SQRT( U0(I,J,1:LM)**2+  V0(I,J,1:LM)**2), 50.0)
                         tausurf_gw   = min(0.5*SQRT(TAUOROX(I , J)**2+TAUOROY(I , J)**2), 10.0) !limit to a very high value
                         aux1=PLE(i,j,LM)/(287.04*(T(i,j,LM)*(1.+0.608*Q(i,j,LM)))) ! air_dens (kg m^-3)
                         hfs = -SH  (i,j) ! W m^-2
                         hfl = -EVAP(i,j) ! kg m^-2 s^-1
                         aux2= (hfs/MAPL_CP + 0.608*T(i,j,LM)*hfl)/aux1 ! buoyancy flux (h+le)
                         aux3= ZLE(I, J, NINT(KPBLSC(I,J)))           ! pbl height (m)
                         !-convective velocity scale W* (m/s)
                         ZWS(i,j) = max(0.,0.001-1.5*0.41*MAPL_GRAV*aux2*aux3/T(i,j,LM))
                         ZWS(i,j) = 1.2*ZWS(i,j)**0.3333 ! m/s
                    	 pi_gw(1, 0:LM) = PLE(I,J,0:LM)
                         theta_tr(1,1:LM) = TH1(I,J,1:LM)
                         rhoi_gw = 0.0
                         pi_gw(1, 0:LM) = 100.0*PLE(I,J,0:LM)
	                     ni_gw = 0.0
                         ti_gw = 0.0
                         tm_gw =ter8
                         pm_gw =plevr8
                         h_gw = 0.0
                         if (FRLAND(I, J) .lt. 0.1) then
        	                 lc_turb(1,1:LM)   =  max(ALH(I,J,1:LM), MIN_ALH)
	                     else
           		             lc_turb(1,1:LM)   =  max(ALH(I,J,1:LM), 50.0)
    	                 end if

                         call   gw_prof (1, LM, 1, tm_gw, pm_gw, pi_gw, &
                                  rhoi_gw, ni_gw, ti_gw, nm_gw) !get Brunt_Vaisala Frequency and midpoint densities


        	             h_gw(1,1:LM)= (2d0*MAPL_PI/LCCIRRUS)*AIRDEN(I, J,1:LM)*uwind_gw(1,1:LM)*nm_gw(1,1:LM)

                  		 where (h_gw .gt. 0.0)
                     		h_gw=sqrt(2.0*tausurf_gw/h_gw)
                  		 end where
                         Wbreak = 0.133*(2d0*MAPL_PI/LCCIRRUS)*uwind_gw/nm_gw !Vertical velocity variance at saturation

	        		     wparc_gw=(2d0*MAPL_PI/LCCIRRUS)*uwind_gw*h_gw*0.133  	        !account for gravity wave breaking

               	         wparc_gw = min(wparc_gw, Wbreak)
                         wparc_gw=wparc_gw*wparc_gw

                         wparc_turb(1,1:LM)  =TKE(I, J, 1:LM)
                         do K = KMIN_TROP(I, J), LM-1
                             if (FRLAND(I, J) .lt. 0.1) then
                       	        if (LTS(I, J) .gt. LTS_LOW) then
                                 if (K .ge. kbmin-2) wparc_ls(1, K) = max(wparc_ls(1,K)+ zws(i, j), 0.00)*SCWST ! add convective velocity within the PBL
                               end if
                             end if
                             if (K .ge. kbmin-2) wparc_ls(1, K)=max(wparc_ls(1,K)+ zws(i, j), 0.00)
                             if (K .ge. kbmin-2) wparc_turb(1, K)=max(wparc_turb(1,K), 0.04)    !minimum velocity within the PBL (not resolved by RAS)

               		         swparc(1, K)=sqrt(wparc_gw(1, K)+wparc_turb(1, K)+ wparc_cgw(1, K))
                    	 end do

                      else
                     	swparc(1,1:LM)  = WSUB_CLIM(I, j, 1:LM)
                      end if


                         ter8(1,1:LM) = T(I,J,1:LM)
                         plevr8(1,1:LM) = PLE(I,J,1:)
                         ndropr8(1,1:LM) = NCPL(I, J, 1:LM)
                         qir8(1,1:LM) =  QILS(I, J,1:LM)+QICN(I, J,1:LM)
                         qcr8(1,1:LM) =  QLLS(I, J,1:LM)+QLCN(I, J,1:LM)
                         npre8(1,1:LM) = NPRE_FRAC_2d(I,J)*NCPI(I,J,1:LM)
                         where ((npre8 .gt. 0.0)   .and. (qir8 .gt. 0.0))
                             dpre8    = ( qir8/(5400.0*npre8*MAPL_PI))**(0.33) !Assume exponential distribution
                         elsewhere
                            dpre8=1.0e-9
                         end where

               ! ==========================================================================================
               ! ========================Activate the aerosols ============================================



                do K = KMIN_TROP(I, J), LM-1 !limit to troposphere and no activation at the surface

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

                     !!Subroutine aerosol_activate contains the CCN activation and ice nucleation parameterizations. Lives in aer_cloud.F90.

                     call   aerosol_activate(ter8(1, k), plevr8(1, K), swparc(1, K), wparc_ls(1, K),  AeroAux, &
                          npre8(1, k), dpre8(1, k), ccn_diag, ndropr8(1, k), qcr8(1, K), &
                          npccninr8(1, K), smaxliq(1, K), naair8(1, K), smaxicer8(1, K), nheticer8(1, K), &
                          nhet_immr8(1, K), dnhet_immr8(1, K), nhet_depr8(1, k), sc_icer8(1, k), &
                          dust_immr8(1, K), dust_depr8(1, k), nlimicer8(1, k), use_average_v, int(CCN_PARAM), int(IN_PARAM),  &
                          so4x(1, k), seasaltx(1, k), dustx(1, k), orgx(1, K), bcx(1, k), &
                                      fdust_drop, fsoot_drop, pfrz_inc_r8(1, K), rh1_r8, frachet_dust, frachet_bc, frachet_org, frachet_ss, int(Immersion_PARAM))

                      CCN01(I, J, K) = max(ccn_diag(1), 0.0)
                      CCN04(I, J, K) = max(ccn_diag(2), 0.0)
                      CCN1 (I, J, K) = max(ccn_diag(3), 0.0)

                      if (K .ge. kbmin-6) npccninr8(1, K) = max(npccninr8(1, K), (1.0-CNV_FRC(I, J))*MINCDNC*1.e6)

               end do

               SMAXL(I, J, 1:LM) = real(smaxliq(1,1:LM)*100.0)
               SMAXI(I, J, 1:LM) = real(smaxicer8(1,1:LM)*100.0)
               NHET_NUC(I, J, 1:LM)  = real(nheticer8(1,1:LM))
               NLIM_NUC(I, J, 1:LM) =  real(nlimicer8(1,1:LM))
               SC_ICE(I, J, 1:LM) = real(sc_icer8(1,1:LM))
               CDNC_NUC(I,J,1:LM)    = real(npccninr8(1,1:LM))
               INC_NUC (I,J,1:LM)    = real(naair8(1,1:LM)  )
               NHET_IMM(I, J, 1:LM)  = real(max(nhet_immr8(1,1:LM), 0.0))
               DNHET_IMM(I, J, 1:LM)  = real(max(dnhet_immr8(1,1:LM), 0.0))
               NHET_DEP(I, J, 1:LM)  = real(nhet_depr8(1,1:LM))
               DUST_IMM(I, J, 1:LM)  = real(max(dust_immr8(1,1:LM), 0.0))
               DUST_DEP(I, J, 1:LM)  = real(max(dust_depr8(1,1:LM), 0.0))
               WSUB (I, J, 1:LM) =  real(wparc_ls(1,1:LM)+swparc(1,1:LM)*0.8)
               SIGW_GW (I, J, 1:LM)   = real( wparc_gw(1,1:LM))
               SIGW_CNV (I, J, 1:LM)   =  real(wparc_cgw(1,1:LM))
               SIGW_TURB (I, J, 1:LM) = real(wparc_turb(1,1:LM))
               SIGW_RC (I, J, 1:LM)   =  real(wparc_ls(1,1:LM))
               PFRZ (I, J, 1:LM)   =  real(pfrz_inc_r8(1,1:LM))

               SO4(I, J, 1:LM)=real(so4x(1,1:LM))
               DUST(I, J, 1:LM)=real(dustx(1,1:LM))
               BCARBON(I, J, 1:LM)=real(bcx(1,1:LM))
               ORG(I, J, 1:LM)=real(orgx(1,1:LM))
               SEASALT(I, J, 1:LM)=real(seasaltx(1,1:LM))

            enddo
         enddo

             call MAPL_TimerOff(MAPL,"---ACTIV", __RC__)

         !=============================================End cloud particle nucleation=====================================
         !===============================================================================================================



    !==========================================================================================================
    !===================================Cloud Macrophysics ====================================================
    !==========================================================================================================

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
    call MAPL_GetPointer(EXPORT, PDF_A,     'PDF_A'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WTHV2,     'WTHV2'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WQL,       'WQL'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PDFITERS, 'PDFITERS', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Unused Exports (foreced to 0.0)
    call MAPL_GetPointer(EXPORT, PTR2D,  'CN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'CN_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_SNR'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0


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
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          QRAIN = QRAIN + PTR3D*DT_MOIST
        endif
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          QSNOW = QSNOW + PTR3D*DT_MOIST
        endif


       ! evap/subl/pdf
        call MAPL_GetPointer(EXPORT, RHCRIT3D,  'RHCRIT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        do L=1,LM
          do J=1,JM
           do I=1,IM
       ! Send the condensates through the pdf after convection
             ! based on Quass 2012 https://doi.org/10.1029/2012JD017495
             if (EIS(I,J) > 5.0) then ! Stable
                ALPHA = 1.0 - ((1.0-dw_land ) + (0.99 - (1.0-dw_land ))*exp(1.0-(PLEmb(i,j,LM)/PLEmb(i,j,l))**2))
             else ! Unstable
                ALPHA = 1.0 - ((1.0-dw_ocean) + (0.99 - (1.0-dw_ocean))*exp(1.0-(PLEmb(i,j,LM)/PLEmb(i,j,l))**4))
             endif
             ! include area scaling and limit RHcrit to > 70%
             ALPHA = min( 0.30, ALPHA*SQRT(SQRT(AREA(I,J)/1.e10)) )
           ! fill RHCRIT export
           if (associated(RHCRIT3D)) RHCRIT3D(I,J,L) = 1.0-ALPHA

           ALPH3D(I, J, L) =  ALPHA

                 end do ! IM loop
         end do ! JM loop
       end do ! LM loop


       ! Put condensates in touch with the PDF


       if (.true.) then
        call MAPL_TimerOn(MAPL,"----hystpdf")

         do L=1,LM
          do J=1,JM
           do I=1,IM

            DLPDF_X(I, J, L)=  QLLS(I, J, L) +QLCN(I, J, L)
            DIPDF_X(I, J, L)=  QILS(I, J, L) +QICN(I, J, L)

             call hystpdf( &
                      DT_MOIST       , &
                      ALPH3D(I, J, L)          , &
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
                      NCPL(I,J,L)   , &
                      NCPI(I,J,L)   , &
                      WHL(I,J,L)     , &
                      WQT(I,J,L)     , &
                      HL2(I,J,L)     , &
                      QT2(I,J,L)     , &
                      HLQT(I,J,L)    , &
                      W3(I,J,L)      , &
                      W2(I,J,L)      , &
                      QT3(I,J,L)     , &
                      HL3(I,J,L)     , &
                      EDMF_FRC(I,J,L), &
                      PDF_A(I,J,L)   , &
                      PDFITERS(I,J,L), &
                      WTHV2(I,J,L)   , &
                      WQL(I,J,L)     , &
                      .false.        , &
                      .true.)

         DLPDF_X(I, J, L)=((QLLS(I, J, L)+QLCN(I, J, L)) - DLPDF_X(I, J, L))/DT_MOIST
         DIPDF_X(I, J, L)=((QILS(I, J, L)+QICN(I, J, L)) - DIPDF_X(I, J, L))/DT_MOIST

           end do ! IM loop
         end do ! JM loop
       end do ! LM loop

       call MAPL_TimerOff(MAPL,"----hystpdf")
      end if

       do L=1,LM
          do J=1,JM
           do I=1,IM


       ! evaporation for CN/LS
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
       ! sublimation for CN/LS
             RHCRIT = 1.0 - ALPH3D(I,J,L)
             SUBLC(I,J,L) =   Q(I,J,L)
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
       ! cleanup clouds
             call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
             RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )

           end do ! IM loop
         end do ! JM loop
       end do ! LM loop


	! Clean up any negative specific humidity before the microphysics scheme
      !-----------------------------------------
         !make sure QI , NI stay within T limits
         call meltfrz_inst2M  (     &
              IM,JM,LM    , &
              T              , &
              QLLS          , &
              QLCN         , &
              QILS           , &
              QICN          , &
              NCPL         , &
              NCPI          )

        call fix_up_clouds_2M( &
         Q, &
         T, &
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
          call FILLQ2ZERO( Q, MASS, TMP2D)
          call FILLQ2ZERO( QGRAUPEL, MASS, TMP2D)
          call FILLQ2ZERO( QRAIN, MASS, TMP2D)
          call FILLQ2ZERO( QSNOW, MASS, TMP2D)
          call FILLQ2ZERO( QLLS, MASS, TMP2D)
          call FILLQ2ZERO( QLCN, MASS, TMP2D)
          call FILLQ2ZERO( QILS, MASS, TMP2D)
          call FILLQ2ZERO( QICN, MASS, TMP2D)



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


 !=============================================End cloud macrophysics=====================================
 !=========================================================================================================



 !==================================================================================================================
 !===============================================Two-moment stratiform cloud microphysics ==========================
 !==================================================================================================================


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

    FQA  = 0.0
    QCNTOT = QLCN+QICN
    QL_TOT = QLCN+QLLS
    QI_TOT = QICN+QILS
    QTOT   = QL_TOT+QI_TOT

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

    ! Update TH
    TH1 = T/PKmb

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

    !Tuning factors
    accre_enhanr8= ACC_ENH
    accre_enhan_icer8= ACC_ENH_ICE
    QCVAR_EXP = 2.0
    autscx = 1.0

    disp_liu = LIU_MU
    ui_scale = UISCALE
    urscale  = URSCALE
    ts_autice = TS_AUTO_ICE
    if (MTIME .le. 0.0) then
    	mtimesc  = DT_MOIST
    else
    	mtimesc=MTIME
    end if

      do J=1,JM
            do I=1,IM

               kbmin =1
               npccninr8  = 0.0
               naair8     = 0.0
               rndstr8 = 2.0e-7
               naconr8   = 0.

               cldfr8(1,1:LM)  = RAD_CF(I,J,1:LM) !Assume minimum overlap
               liqcldfr8(1,1:LM)  =  CFLIQ(I,J,1:LM)
               icecldfr8(1,1:LM)  =  CFICE(I,J,1:LM)

               cldor8          = cldfr8
               ter8(1,1:LM)       = T(I,J,1:LM)
               qvr8(1,1:LM)       = Q(I,J,1:LM)

               qcr8(1,1:LM)        = QL_TOT(I,J,1:LM)
               qir8(1,1:LM)        = QI_TOT(I,J,1:LM)
               ncr8(1,1:LM)        = MAX(  NCPL(I,J,1:LM), 0.0)
               nir8(1,1:LM)        = MAX(  NCPI(I,J,1:LM), 0.0)

               ! Nucleation  tendencies
               naair8(1,1:LM)     = max(( INC_NUC(I, J, 1:LM)*cldfr8(1,1:LM) - nir8(1,1:LM))/DT_MOIST, 0.0)
               npccninr8(1,1:LM)  = max((CDNC_NUC(I, J, 1:LM)*cldfr8(1,1:LM) - ncr8(1,1:LM))/DT_MOIST, 0.0)

               where  ((naair8 .gt. 1.0e3)) ! add cloud fraction if nucleation is happening 2018
                   icecldfr8 = max(0.05,  icecldfr8)
               end where

               where (cldfr8(1,:) .ge. 0.001)
                  nimmr8(1,1:LM) = MIN(DNHET_IMM(I, J, 1:LM), ncr8(1,1:LM)/cldfr8(1,1:LM)/DT_MOIST) !tendency
               elsewhere
                  nimmr8(1,1:LM) = 0.0
               end where

               nhet_depr8(1,1:LM) = NHET_DEP(I, J, 1:LM)/DT_MOIST !becomes a tendency (could be done a bit better)
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
                  call getINsubset(2,  AeroAux, AeroAux_b)
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



               if (AUTSC .gt. 0.0) then
                  autscx = AUTSC
               else
                  autscx =  min(max(0., (300.0 - T(I,J,LM))/ABS(AUTSC)), 1.0)
                  autscx  =  1.0 - 0.995*autscx
               end if



  !!!!================Estimate qcvar following Xie and Zhang, JGR, 2015
                 HMOIST_950 = 0.0
                 HSMOIST_500 = 0.0
                 IF (PLmb(I, J, LM) .le. 500.0) then
                    qcvarr8  = 2.0
                 ELSEIF (PLmb(I, J, LM) .lt. 950.0) then
                    DO K=LM, 1, -1
                         if (PLmb(I,J,K) .lt. 500.0) exit
                         HSMOIST_500 = MAPL_CP*T(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO
                    HMOIST_950 = MAPL_CP*T(I, J, LM) + GZLO(I, J, LM) + Q(I, J, LM)*MAPL_ALHL
                    SINST = (HMOIST_950 -  HSMOIST_500)/(PLmb(I,J,LM)*100.0- 50000.0)
                 ELSE
                    DO K=LM, 1, -1
                         if (PLmb(I,J,K) .lt. 500.0) exit
                         HSMOIST_500 = MAPL_CP*T(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO
                    DO K=LM, 1, -1
                         if (PLmb(I,J,K) .lt. 950.0) exit
                         HMOIST_950 = MAPL_CP*T(I, J, K) + GZLO(I, J, K) + Q(I, J, K)*MAPL_ALHL
                    END DO
                    SINST = (HMOIST_950 -  HSMOIST_500)/45000.0
                  ENDIF

                  xscale = (9000.0/imsize)**(-0.666)
                  qcvarr8 =  0.67 -0.38*SINST +  4.96*xscale - 8.32*SINST*xscale
                  qcvarr8 = min(max(qcvarr8, 0.5), 50.0)
                  if (associated(QCVAR_EXP)) QCVAR_EXP(I, J) = real(qcvarr8)
                  relvarr8 = qcvarr8

               ! for MG23 (initial values)
                  frzimmr8 =  nimmr8
                  frzcntr8 = nimmr8*0.0
                  frzdepr8 = nhet_depr8
                  qrr8(1,1:LM)     =  QRAIN(I, J,1:LM)
                  qsr8(1,1:LM)     =  QSNOW(I, J,1:LM)
                  qgr8(1,1:LM)     =  QGRAUPEL(I, J,1:LM)
                  nrr8(1,1:LM)     =  NRAIN(I, J,1:LM)
                  nsr8(1,1:LM)     =  NSNOW(I, J,1:LM)
                  ngr8(1,1:LM)     =  NGRAUPEL(I, J,1:LM)
                  qsatfacr8 = 1.0
                  SCICE_tmp(1,1:LM)  =  SC_ICE(I, J, 1:LM)
                  FQA_tmp(1,1:LM)  = FQA(I, J, 1:LM)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !CALLS to MG versions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

         call  micro_mg_tend_interface ( DT_MICRO, INT(PDFSHAPE), ALPH3D(I, J, 1:LM), SCICE_tmp, FQA_tmp, &
                             ncolmicro,             LM,               dt_r8,       &
                             CNV_FRC(I,J), SRF_TYPE(I,J), &
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
!                             errstring, & ! Below arguments are "optional" (pass null pointers to omit).
                      !       tnd_qsnow,          tnd_nsnow,          re_ice,    &
                             prer_evap, &
                             frzimmr8,             frzcntr8,              frzdepr8,  & ! contact is not passed since it depends on the droplet size dist
                             nsootr8, rnsootr8,  & ! soot for contact IN
                             npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, &  ! Number tendencies for liquid
                             npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, nsacwior8,  &  ! Number tendencies for ice
                             ts_autice, ui_scale, autscx , disp_liu, nbincontactdust, urscale)


    end if

        IF (MGVERSION > 1) then

                   QRAIN(I,J,1:LM)  = max(QRAIN(I,J,1:LM) + REAL(qrtendr8(1, 1:LM)*DT_R8), 0.0) ! grid average
                   QSNOW(I,J,1:LM)  = max(QSNOW(I,J,1:LM) + REAL(qstendr8(1, 1:LM)*DT_R8), 0.0) ! grid average
                   NRAIN(I,J,1:LM)  = max(NRAIN(I,J,1:LM) + REAL(nrtendr8(1, 1:LM)*DT_R8), 0.0)
                   NSNOW(I,J,1:LM)  = max(NSNOW(I,J,1:LM) + REAL(nstendr8(1, 1:LM)*DT_R8), 0.0)
                   CLDREFFR(I,J,1:LM)  = REAL(reff_rainr8(1,1:LM))
                   CLDREFFS(I,J,1:LM)  = REAL(reff_snowr8(1,1:LM))/scale_ri
                   CLDREFFG(I,J,1:LM)  = REAL(reff_graur8(1,1:LM))/scale_ri

                  IF (MGVERSION .gt. 2) then
                      QGRAUPEL(I,J,1:LM)  = max(QGRAUPEL(I,J,1:LM) + REAL(qgtendr8(1, 1:LM)*DT_R8), 0.0) ! grid average
                      NGRAUPEL(I,J,1:LM)  = max(NGRAUPEL(I,J,1:LM) + REAL(ngtendr8(1, 1:LM)*DT_R8), 0.0)
                  else
                      QGRAUPEL(I,J,1:LM)  = 0.0 ! grid average
                      NGRAUPEL(I,J,1:LM)  = 0.0 ! grid average
                  end if


        else
                   QRAIN(I,J,1:LM)  = max(REAL(qrout2r8(1,1:LM)), 0.0) ! grid average
                   QSNOW(I,J,1:LM)  = max(REAL(qsout2r8(1,1:LM)), 0.0)
                   NRAIN(I,J,1:LM)  = max(REAL(nrout2r8(1,1:LM)), 0.0)
                   NSNOW(I,J,1:LM)  = max(REAL(nsout2r8(1,1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM) = REAL(drout2r8(1,1:LM))/2.0
                   CLDREFFS(I,J,1:LM) = REAL(dsout2r8(1,1:LM))/2.0/scale_ri
                   QGRAUPEL(I,J,1:LM)  = 0.0 ! grid average
                   NGRAUPEL(I,J,1:LM)  = 0.0 ! grid average
         end if

               PFL_LS(I, J, 1:LM) = rflxr8(1,1:LM) !+ lflxr8(1,1:LM)
               PFI_LS(I, J, 1:LM) = sflxr8(1,1:LM) !+ gflxr8(1,1:LM) +  iflxr8(1,1:LM)

               !Update state after microphysisc
               LS_PRCP(I,J)     = max(1000.*REAL((prectr8(1)-precir8(1))), 0.0)
               LS_SNR(I,J)      = max(1000.*REAL(precir8(1)), 0.0)
               QL_TOT(I,J,1:LM) = max(QL_TOT(I,J,1:LM)   + REAL(qctendr8(1,1:LM)) * DT_R8, 0.0)
               QI_TOT(I,J,1:LM) = max(QI_TOT(I,J,1:LM)   + REAL(qitendr8(1,1:LM)) * DT_R8, 0.0)
               Q(I,J,1:LM)   = MAX(Q(I,J,1:LM)     + REAL(qvlatr8(1,1:LM)) * DT_R8, 0.0)
               T(I,J,1:LM) = T(I,J,1:LM)   + REAL(tlatr8(1,1:LM)) * DT_R8 / (MAPL_CP)
               NCPL(I,J,1:LM) = MAX(NCPL(I,J,1:LM)   + REAL(nctendr8(1,1:LM)) * DT_R8, 0.0)
               NCPI(I,J,1:LM) = MAX(NCPI(I,J,1:LM)   + REAL(nitendr8(1,1:LM)) * DT_R8, 0.0)


               CLDREFFL(I,J,1:LM) = max(REAL(effcr8(1,1:LM))*1.0e-6, 1.0e-6)
               CLDREFFI(I,J,1:LM) = max(REAL(effir8(1,1:LM))*1.0e-6, 1.0e-6)/scale_ri !scale to match the Dge definition of Fu 1996

               ! diagnostics from the microphysics********************

               RSU_LS(I,J,1:LM) = REAL(evapsnowr8(1,1:LM))                    !Snow evap
               REV_LS(I,J,1:LM) = REAL(nevaprr8(1,1:LM) )                        !rain evap
               SUBLC(I,J,1:LM) = REAL(cmeioutr8(1,1:LM)) + SUBLC(I,J,1:LM)     ! Ice subl already grid -ave
               BERGS(I,J,1:LM) = REAL(bergsor8(1,1:LM))                               ! Snow Bergeron
               FRZ_TT_X(I,J,1:LM) = REAL(mnucccor8(1,1:LM)+ mnucctor8(1,1:LM) + homoor8(1,1:LM))!ice mixing ratio from nucleation (hom+het)
               FRZ_PP_X(I,J,1:LM) = REAL(mnuccror8(1,1:LM) + pracsor8(1,1:LM)) !freezing of rain (hom+ het freezing and accretion by snow)
               MELT(I,J,1:LM) = REAL(meltor8(1,1:LM))                            !melting of cloud ice  and snow
               SDM_X(I,J,1:LM) = REAL(qisedtenr8(1,1:LM))                    ! ice sed
               EVAPC(I,J,1:LM) = REAL(qcsevapr8(1,1:LM) )  +  EVAPC(I,J,1:LM)   ! cloud evap
               BERG(I,J,1:LM)  =  REAL(bergor8(1,1:LM))                          ! Bergeron process
               ACIL_LS_X(I,J,1:LM) =REAL(psacwsor8(1,1:LM) + msacwior8(1,1:LM))   !Acc + collection of cloud  by snow
               QCRES(I,J,1:LM) =REAL(qcresor8(1,1:LM) )             !residual liquid condensation
               QIRES(I,J,1:LM) =REAL(qiresor8(1,1:LM))                  !residual ice condensation

               ACLL_LS_X(I,J,1:LM) =REAL(praor8(1,1:LM) )          ! Acc cloud by rain
               AUT_X(I,J,1:LM) = REAL(prcor8(1,1:LM))                  ! Aut liquid
               AUTICE(I,J,1:LM) = REAL(prcior8(1,1:LM))                !Aut  ice
               ACIL_AN_X(I,J,1:LM) = REAL(praior8(1,1:LM))           !Acc ice  by snow
               ACLL_AN_X(I,J,1:LM) = REAL(msacwior8(1,1:LM))   !HM process

               FRZPP_LS(I,J,1:LM) = REAL(frzrdtr8(1,1:LM)) / MAPL_CP !precip freezing latent heat rate
               SNOWMELT_LS(I,J,1:LM) =REAL(meltsdtr8(1,1:LM))/ MAPL_CP !melting of snow latent heat rate

               !diagnostics for number concentration budget (all grid-average, Kg-1 s-1)

               DNDCCN(I,J,1:LM)        = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               DNDACRLS(I,J,1:LM)      = REAL(npsacwsor8(1,1:LM)) !droplet number tendency from accretion by snow
               DNDACRLR(I,J,1:LM)      = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               DNDEVAPC(I,J,1:LM)      = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation
               DNDAUTLIQ(I,J,1:LM)     = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion

               DNCACRIS (I,J,1:LM)     = REAL(npraior8(1,1:LM))  !ice number tendency from accretion by snow
               DNHET_CT(I,J,1:LM)      = REAL(nnucctor8(1,1:LM)) !ice number tendency from contact IN
               DNHET_IMM(I,J,1:LM)     = REAL(nnucccor8(1,1:LM)) !ice number tendency from immersion IN
               DNCNUC(I,J,1:LM)        = REAL(nnuccdor8(1,1:LM)) !ice number tendency from nucleation on aerosol
               DNCSUBL (I,J,1:LM)      = REAL(nsubior8(1,1:LM))  !ice number tendency from sublimation
               DNCAUTICE (I,J,1:LM)    = REAL(nprcior8(1,1:LM))  !ice number tendency from autoconversion
               DNCHMSPLIT(I,J,1:LM)    = REAL(nsacwior8(1,1:LM)) !ice number tendency from H-M process

               DNDCCN(I,J,1:LM)        = REAL(npccnor8(1,1:LM))   !droplet number tendency from CCN activation
               DNDACRLS(I,J,1:LM)      = REAL(npsacwsor8(1,1:LM) )!droplet number tendency from accretion by snow
               DNDACRLR(I,J,1:LM)      = REAL(npraor8(1,1:LM))    !droplet number tendency from accretion by rain
               DNDEVAPC(I,J,1:LM)      = REAL(nsubcor8(1,1:LM))   !droplet number tendency from evaporation
               DNDAUTLIQ(I,J,1:LM)     = REAL(nprc1or8(1,1:LM))   !droplet number tendency from autoconversion

            enddo !I
         enddo !J
         !============================================Finish 2-moment micro implementation===========================

         !update water tracers
         QLCN=QL_TOT*FQA
         QLLS=QL_TOT-QLCN
         QICN=QI_TOT*FQA
         QILS=QI_TOT-QICN
         QTOT= QICN+QILS+QLCN+QLLS

         !============ Put cloud fraction back in contact with the PDF and create new condensate if neccesary (Barahona et al., GMD, 2014)============

    do K= 1, LM
       do J=1,JM
            do I=1,IM
                  call update_cld( &
                         DT_MOIST                , &
                         ALPH3D(I, J, K)        , &
                         PDFSHAPE , &
                         CNV_FRC(I, J)           , &
                         SRF_TYPE(I, J)          , &
                         PLmb(I, J, K)            , &
                         Q (I, J, K)            , &
                         QLLS(I, J, K)           , &
                         QLCN(I, J, K)           , &
                         QILS(I, J, K)           , &
                         QICN(I, J, K)           , &
                         T(I, J, K)           , &
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

  		call meltfrz_inst2M  (     &
              IM,JM,LM    , &
              T              , &
              QLLS          , &
              QLCN         , &
              QILS           , &
              QICN          , &
              NCPL         , &
              NCPI          )

         RAD_CF =min(CLLS+CLCN, 1.0)

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

         RAD_QV = MAX( Q , 0. )
         RAD_QL = MAX(MIN( RAD_QL , 0.001 ), 0.0)  ! Still a ridiculously large
         RAD_QI = MAX(MIN( RAD_QI , 0.001 ), 0.0)  ! value.
         RAD_QR = MAX(MIN( RAD_QR , 0.01 ), 0.0)  ! value.
         RAD_QS = MAX(MIN( RAD_QS , 0.01 ), 0.0)  ! value
         RAD_QG = MAX(MIN( RAD_QG , 0.01 ), 0.0)  ! value


      ! Fill GEOS precip diagnostics
         PRCP_RAIN =  LS_PRCP
         PRCP_SNOW = LS_SNR
         ICE     = PRCP_ICE + PRCP_GRAUPEL
         FRZR    = 0.0
      ! Redistribute precipitation fluxes for chemistry
         TMP3D =  MIN(1.0,MAX(QLCN/MAX(RAD_QL,1.E-8),0.0))
         PFL_AN(:,:,1:LM) = PFL_LS(:,:,1:LM) * TMP3D
         PFL_LS(:,:,1:LM) = PFL_LS(:,:,1:LM) - PFL_AN(:,:,1:LM)
         TMP3D =  MIN(1.0,MAX(QICN/MAX(RAD_QI,1.E-8),0.0))
         PFI_AN(:,:,1:LM) = PFI_LS(:,:,1:LM) * TMP3D
         PFI_LS(:,:,1:LM) = PFI_LS(:,:,1:LM) - PFI_AN(:,:,1:LM)
      ! cleanup suspended precipitation condensates
         call FIX_NEGATIVE_PRECIP(RAD_QR, RAD_QS, RAD_QG)

         !=================================================================================
         !    Units conversion for diagnostics


         !to m-3
         NCPL_VOL=NCPL*AIRDEN !
         NCPI_VOL=NCPI*AIRDEN
         CDNC_NUC=CDNC_NUC*AIRDEN
         INC_NUC =INC_NUC*AIRDEN
        !to m-3 s-1
         DNHET_CT    = DNHET_CT*AIRDEN
         DNHET_IMM   = DNHET_IMM*AIRDEN
         DNCNUC      = DNCNUC*AIRDEN
         DNCHMSPLIT  = DNCHMSPLIT*AIRDEN
         DNCSUBL     = DNCSUBL*AIRDEN
         DNCACRIS    = DNCACRIS*AIRDEN
         DNCAUTICE   = DNCAUTICE*AIRDEN
         DNICNV      = DNICNV*AIRDEN

         DNDCCN       = DNDCCN*AIRDEN
         DNDACRLS     = DNDACRLS*AIRDEN
         DNDACRLR     = DNDACRLR*AIRDEN
         DNDEVAPC     = DNDEVAPC*AIRDEN
         DNDAUTLIQ    = DNDAUTLIQ*AIRDEN
         DNDCNV       = DNDCNV*AIRDEN

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

         TH1 = T / PKmb

       ! !Set rain water for radiation to 0 if preciprad flag is off (set to 0)
        !if(CLDPARAMS%PRECIPRAD .eq. 0.) then
        !   RAD_QR = 0.
       !    RAD_QS = 0.
       !    RAD_QG = 0.
       ! endif

         CLDREFFL = MAX(MIN_RL, CLDREFFL) !DONIF Limits according to MG2008-I
         CLDREFFL = MIN(MAX_RL, CLDREFFL)
         CLDREFFI = MAX(MIN_RI, CLDREFFI)
         CLDREFFI = MIN(MAX_RI, CLDREFFI)  !maximum number for the correlation and modis sim

         CLDREFFR = MAX(MIN_RL, CLDREFFR)
         CLDREFFR = MIN(MAX_RL, CLDREFFR)
         CLDREFFS = MAX(MIN_RI*2., CLDREFFS)
         CLDREFFS = MIN(MAX_RI*2., CLDREFFS)  !maximum number for the correlation and modis sim
         CLDREFFG = MAX(MIN_RI*2., CLDREFFG)
         CLDREFFG = MIN(MAX_RI*2., CLDREFFG)  !maximum number for the correlation and modis sim

         !===========================

         ! Diagnostic cloud top/base properties

         CLDREFFL_TOP_X = MAPL_UNDEF
         CLDREFFI_TOP_X =MAPL_UNDEF
         NCPL_TOP_X = MAPL_UNDEF
         NCPL_CLDBASEX = MAPL_UNDEF
         NCPI_TOP_X = MAPL_UNDEF
         kbmin = LM

         DO I=1, IM
            DO J= 1 , JM

               cfaux(1,1:LM) =  CFLIQ(I, J, 1:LM)
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

               cfaux(1,1:LM) =CFICE(I, J, 1:LM)
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

      call MAPL_GetPointer(EXPORT, PTR2D, 'CLDREFFI_TOP', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D =  CLDREFFI_TOP_X

      call MAPL_GetPointer(EXPORT, PTR2D, 'CLDREFFL_TOP', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D =  CLDREFFI_TOP_X

      call MAPL_GetPointer(EXPORT, PTR2D, 'NCPL_CLDBASE', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D=  NCPL_CLDBASEX

      call MAPL_GetPointer(EXPORT, PTR2D, 'NCPL_TOP', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D=  NCPL_TOP_X

      call MAPL_GetPointer(EXPORT, PTR2D, 'NCPI_TOP', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D=  NCPI_TOP_X




      ! Clean up Relative Humidity where RH > 110%
      !---------------------------------------------
      ! moved to Moist GridComp

      if (associated(CCNCOLUMN))   CCNCOLUMN  = SUM(    CCN1*MASS/AIRDEN , 3)
      if (associated(NDCOLUMN ))    NDCOLUMN  = SUM(NCPL_VOL*MASS/AIRDEN , 3)
      if (associated(NCCOLUMN ))    NCCOLUMN  = SUM(NCPI_VOL*MASS/AIRDEN , 3)

      ! Update microphysics tendencies
      if (associated(DQVDT_micro)) DQVDT_micro = ( Q         - DQVDT_micro) / DT_MOIST
      if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
      if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
      if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
      if (associated(DQRDT_micro)) DQRDT_micro = ( QRAIN      - DQRDT_micro) / DT_MOIST
      if (associated(DQSDT_micro)) DQSDT_micro = ( QSNOW      - DQSDT_micro) / DT_MOIST
      if (associated(DQGDT_micro)) DQGDT_micro = ( QGRAUPEL   - DQGDT_micro) / DT_MOIST
      if (associated( DUDT_micro))  DUDT_micro = ( U0         -  DUDT_micro) / DT_MOIST
      if (associated( DVDT_micro))  DVDT_micro = ( V0         -  DVDT_micro) / DT_MOIST
      if (associated( DTDT_micro))  DTDT_micro = ( T       -  DTDT_micro) / DT_MOIST



      call MAPL_TimerOff (MAPL,"---CLDMICRO", __RC__)

      ! Exports


       call MAPL_GetPointer(EXPORT, PTR3D, 'SCF', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           WHERE ((QLLS+QLCN+QILS+QICN) .gt. 1.0e-12)
             PTR3D = (QLLS+QLCN)/(QLLS+QLCN+QILS+QICN)
           ELSEWHERE
             PTR3D= MAPL_UNDEF
           END WHERE
         endif

         call MAPL_GetPointer(EXPORT, PTR3D, 'SCF_ALL', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           WHERE ((QLLS+QLCN+QILS+QICN + QRAIN + QSNOW + QGRAUPEL) .gt. 1.0e-12)
             PTR3D= (QLLS+QLCN+QRAIN)/(QLLS+QLCN+QILS+QICN + QSNOW + QGRAUPEL + QRAIN)
           ELSEWHERE
             PTR3D = MAPL_UNDEF
           END WHERE
         endif


        call MAPL_GetPointer(EXPORT, PTR3D, 'DQRL', RC=STATUS); VERIFY_(STATUS)
        if(associated(PTR3D)) PTR3D = DQRDT_macro + DQRDT_micro

        ! Compute DBZ radar reflectivity
        call MAPL_GetPointer(EXPORT, PTR3D, 'DBZ'    , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, PTR2D, 'DBZ_MAX', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D) .OR. associated(PTR2D)) then
           call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,0)
           if (associated(PTR3D)) PTR3D = TMP3D
           if (associated(PTR2D)) then
              PTR2D=-9999.0
              DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                 PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
              END DO ; END DO ; END DO
           endif
        endif

   call MAPL_TimerOff(MAPL,"--MGB2_2M",__RC__)

end subroutine MGB2_2M_Run



end module GEOS_MGB2_2M_InterfaceMod
