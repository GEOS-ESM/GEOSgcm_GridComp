!, autscx $Id$

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
  use cldmacro
  use aer_cloud
  use micro_mg3_0

  implicit none

  integer :: MGVERSION

  private

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
  real  :: DCS, QCVAR_, WBFFACTOR, NC_CST, NI_CST, NG_CST, MUI_CST, PMIN_CBL
  real  :: LCCIRRUS, UISCALE, SS_SCALE, REEVAP_MICRO, LIU_MU, TFRZ, &
           NPRE_FRAC, QCVAR, ZPBLMAXLL, TMAXLL, LTS_LOW, LTS_UP, MIN_EXP,     &
           BKGTAU, DCRIT_, USE_AV_V, AUTSC, TS_AUTO_ICE, CCN_PARAM, IN_PARAM, &
           FDROP_DUST, FDROP_SOOT, USE_NATURE_WSUB, SIGMA_NUC, MIN_ALH, &
           HMOIST_950, HSMOIST_500, SINST, MAX_EXP, MAX_CAPE, MIN_CAPE,       &
           DUST_INFAC, ORG_INFAC, BC_INFAC, SS_INFAC, RRTMG_IRRAD, RRTMG_SORAD,&
           SCWST, MTIME, SWCIRRUS, MINCDNC, TMAXCFCORR,    &
           Immersion_param, ACC_ENH, ACC_ENH_ICE, DT_MICRO, DT_AUX, UR_SCALE
  integer :: KSTRAP,CBL_METHOD,CLEANUP_RH
  !character(len=ESMF_MAXSTR) :: JASON_TUNING
  character(len=ESMF_MAXSTR) :: GRIDNAME
  character(len=4)           :: imchar
  character(len=2)           :: dateline
  integer                    :: imsize,nn
  character(LEN=ESMF_MAXSTR):: CONVPAR_OPTION

  public :: MGB2_2M_Setup, MGB2_2M_Initialize, MGB2_2M_Run
  public :: MGVERSION

contains

subroutine MGB2_2M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    
    character(len=ESMF_MAXSTR)              :: IAm

    IAm = "MGB2_2M_Setup"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    
    Iam = trim(COMP_NAME) // Iam

    call ESMF_ConfigGetAttribute( CF, MGVERSION, Label="MGVERSION:",  default=2, __RC__)
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
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )  
                                                                              

    call MAPL_AddInternalSpec(GC,                                        &
         SHORT_NAME = 'QLLS',                                            &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%QLLS),                             &
         DIMS       = MAPL_DimsHorzVert,                                 &
         VLOCATION  = MAPL_VLocationCenter,                   __RC__  )  
                                                                              

    call MAPL_AddInternalSpec(GC,                                       &
         SHORT_NAME = 'QLCN',                                           &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS      = 'kg kg-1',                                        &
         FRIENDLYTO = trim(FRIENDLIES%QLCN),                            &
         DIMS       = MAPL_DimsHorzVert,                                &
         VLOCATION  = MAPL_VLocationCenter,                  __RC__  )  
                                                                              

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLLS',                                      &
         LONG_NAME  = 'large_scale_cloud_area_fraction',           &
         UNITS      = '1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%CLLS),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )  
                                                                              

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'CLCN',                                      &
         LONG_NAME  = 'convective_cloud_area_fraction',            &
         UNITS      = '1',                                         &
         FRIENDLYTO = trim(FRIENDLIES%CLCN),                       &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )  
                                                                              

    call MAPL_AddInternalSpec(GC,                                     &
         SHORT_NAME = 'QILS',                                         &
         LONG_NAME  = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                      &
         FRIENDLYTO = trim(FRIENDLIES%QILS),                          &
         DIMS       = MAPL_DimsHorzVert,                              &
         VLOCATION  = MAPL_VLocationCenter,                __RC__  )  
                                                                              

    call MAPL_AddInternalSpec(GC,                                    &
         SHORT_NAME = 'QICN',                                        &
         LONG_NAME  = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS      = 'kg kg-1',                                     &
         FRIENDLYTO = trim(FRIENDLIES%QICN),                         &
         DIMS       = MAPL_DimsHorzVert,                             &
         VLOCATION  = MAPL_VLocationCenter,               __RC__  )  
                                                                              

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
    

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QRAIN',                                     &
         LONG_NAME  = 'mass_fraction_of_rain',                     & 
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QRAIN),                      &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QSNOW',                                     &
         LONG_NAME  = 'mass_fraction_of_snow',                     &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QSNOW),                      &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    

    call MAPL_AddInternalSpec(GC,                                  &
         SHORT_NAME = 'QGRAUPEL',                                  &
         LONG_NAME  = 'mass_fraction_of_graupel',                  &
         UNITS      = 'kg kg-1',                                   &
         FRIENDLYTO = trim(FRIENDLIES%QGRAUPEL),                   &
         default    = 0.0,                                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             __RC__  )
    
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
    

end subroutine MGB2_2M_Setup

subroutine MGB2_2M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL

    logical  :: nccons, nicons, ngcons, do_graupel
    real(ESMF_KIND_R8)  Dcsr8, qcvarr8,  micro_mg_berg_eff_factor_in, ncnstr8, ninstr8, ngnstr8, mui_cnstr8

    character(len=ESMF_MAXSTR)              :: IAm

    IAm = "MGB2_2M_Initialize"
    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )
   
    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4
 
    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , __RC__)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , __RC__)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , __RC__)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', __RC__)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , __RC__)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , __RC__)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , __RC__)

!#ifdef NODISABLE

    call MAPL_GetResource(MAPL, DCS,      'DCS:'    , DEFAULT=350.0e-6, __RC__ )
    
    Dcsr8 = DCS
    call MAPL_GetResource(MAPL, QCVAR_,   'QCVAR:'  , DEFAULT= 2.0 ,__RC__) !variance of the QL distribution     
    
    qcvarr8=QCVAR_
    call MAPL_GetResource(MAPL, WBFFACTOR,   'WBFFACTOR:', DEFAULT= 1.0 ,__RC__) !variance of the QL distribution     
    
    micro_mg_berg_eff_factor_in = WBFFACTOR
    call MAPL_GetResource(MAPL, NC_CST ,  'NC_CST:' , DEFAULT=  0.0 ,__RC__) !constant nd (set if greather than zero)     
    
    call MAPL_GetResource(MAPL, NI_CST ,  'NI_CST:' , DEFAULT=  0.0 ,__RC__) !constant nd (set if greather than zero) 
    
    call MAPL_GetResource(MAPL, NG_CST ,  'NG_CST:' , DEFAULT=  0.0 ,__RC__) !constant ng (set if greather than zero) 
    
    call MAPL_GetResource(MAPL, MUI_CST,  'MUI_CST:', DEFAULT= -1.0 ,__RC__) !constant ng (set if greather than zero) 
    
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
    call WRITE_PARALLEL ("INITIALIZED MG in non-generic GC INIT")

      call MAPL_GetResource(MAPL, CLDPARAMS%RH00,           'RH_CRIT:',        DEFAULT= 1.0     ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%C_ACC,          'ACCRETION:',      DEFAULT= 2.0     ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%C_EV_R,         'RAIN_REVAP_FAC:', DEFAULT= 1.00    ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%C_EV_S,         'SNOW_REVAP_FAC:', DEFAULT= 0.5     ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%REVAP_OFF_P,    'REVAP_OFF_P:',    DEFAULT= 2000.   ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%CNVENVFC,       'CNV_ENVF:',       DEFAULT= 1.0     ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%CNVDDRFC,       'CNV_DDRF:',       DEFAULT= 0.0     ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%CNVICEPARAM,    'CNV_ICEPARAM:',   DEFAULT= 1.0     ,RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%PDFSHAPE,       'PDFSHAPE:',       DEFAULT= 1       ,RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetResource(MAPL, CLDPARAMS%T_ICE_ALL,      'T_ICE_ALL:',      DEFAULT= MAPL_TICE-27.0 ,RC=STATUS); VERIFY_(STATUS)

      call MAPL_GetResource(MAPL, CLDPARAMS%SCLM_SHALLOW    , 'SCLM_SHALLOW:' , DEFAULT= 2.0, RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%SCLM_DEEP       , 'SCLM_DEEP:'    , DEFAULT= 1.0, RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%PDFSHAPE        , 'PDFSHAPE:'     , DEFAULT= 1  , RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%SLOPERHCRIT     , 'SLOPERHCRIT:'  , DEFAULT= 20.0, RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetResource(MAPL, CLDPARAMS%TURNRHCRIT_UPPER, 'TURNRHCRIT_UP:', DEFAULT= 300.0, RC=STATUS); VERIFY_(STATUS)

      ! Horizontal resolution dependant defaults for minimum RH crit
      if( imsize.le.200       ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.80, RC=STATUS)
      if( imsize.gt.200 .and. &
          imsize.le.400       ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.90, RC=STATUS)
      if( imsize.gt.400 .and. &
          imsize.le.800       ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.93, RC=STATUS)
      if( imsize.gt.800 .and. &
          imsize.le.1600      ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.95, RC=STATUS)
      if( imsize.gt.1600 .and. &
          imsize.le.3200      ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.97 ,RC=STATUS)
      if( imsize.gt.3200 .and. &
          imsize.le.6400      ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.98 ,RC=STATUS)
      if( imsize.gt.6400 .and. &
          imsize.le.12800     ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.99 ,RC=STATUS)
      if( imsize.gt.12800     ) call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT, 'MINRHCRIT:', DEFAULT=0.99 ,RC=STATUS)

      call MAPL_GetResource( MAPL, CLDPARAMS%CNV_BETA,       'CNV_BETA:',       DEFAULT= 10.0    )
      call MAPL_GetResource( MAPL, CLDPARAMS%CCW_EVAP_EFF,   'CCW_EVAP_EFF:',   DEFAULT= 5.0e-4  )
      call MAPL_GetResource( MAPL, CLDPARAMS%CCI_EVAP_EFF,   'CCI_EVAP_EFF:',   DEFAULT= 4.0e-3  )
      call MAPL_GetResource( MAPL, CLDPARAMS%TURNRHCRIT,     'TURNRHCRIT:',     DEFAULT= 884.0   )
      call MAPL_GetResource(MAPL, LCCIRRUS,       'LCCIRRUS:',       DEFAULT= 500.0,  __RC__) !Characteristic Length (m) of high freq gravity waves
      call MAPL_GetResource(MAPL, UISCALE,        'UISCALE:',        DEFAULT= 1.0,    __RC__) !Scaling factor for sed vel of ice      
      call MAPL_GetResource(MAPL, LIU_MU,         'LIU_MU:',         DEFAULT= 2.0,    __RC__) !Liu autoconversion parameter
      call MAPL_GetResource(MAPL, NPRE_FRAC,      'NPRE_FRAC:',      DEFAULT= -1.0,   __RC__) !Fraction of preexisting ice affecting ice nucleationn            
      call MAPL_GetResource(MAPL, LTS_LOW,        'LTS_LOW:',        DEFAULT= 20.0,   __RC__) !lower LTS for morphology correction
      call MAPL_GetResource(MAPL, LTS_UP,         'LTS_UP:',         DEFAULT= 22.0,   __RC__) !Upper LTS for morphology correction
      call MAPL_GetResource(MAPL, MIN_EXP,        'MIN_EXP:',        DEFAULT= 0.5,    __RC__) !Exponent of the relation CFA=CFV^n
      call MAPL_GetResource(MAPL, MAX_EXP,        'MAX_EXP:',        DEFAULT= 1.0,    __RC__) !Exponent of the relation CFA=CFV^n
      call MAPL_GetResource(MAPL, USE_AV_V,       'USE_AV_V:',       DEFAULT= 1.0,    __RC__) !Set to > 0 to use an average velocity for activation
      call MAPL_GetResource(MAPL, AUTSC,          'AUT_SCALE:',      DEFAULT= 0.5,    __RC__) !scale factor for critical size for drizzle
      call MAPL_GetResource(MAPL, TS_AUTO_ICE,    'TS_AUTO_ICE:',    DEFAULT= 4.0,    __RC__) !Ice autoconversion time scale
      call MAPL_GetResource(MAPL, TMAXLL,         'TMAXLL:',         DEFAULT= 250.0,  __RC__) !Liquid clouds min T
      call MAPL_GetResource(MAPL, CCN_PARAM,      'CCNPARAM:',       DEFAULT= 2.0,    __RC__) !CCN activation param
      call MAPL_GetResource(MAPL, IN_PARAM,       'INPARAM:',        DEFAULT= 6.0,    __RC__) !IN param
      call MAPL_GetResource(MAPL, Immersion_param,'ImmersionPARAM:', DEFAULT= 6.0,    __RC__) !Immersion param
      call MAPL_GetResource(MAPL, ACC_ENH,        'ACC_ENH:',        DEFAULT= 1.0,    __RC__) !accretion rain-liquid scaling for MG2
      call MAPL_GetResource(MAPL, ACC_ENH_ICE,    'ACC_ENH_ICE:',    DEFAULT= 1.0,    __RC__) !accretion snow-ice scaling for MG2
      call MAPL_GetResource(MAPL, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.5,    __RC__) !Fraction of dust within droplets for immersion freezing
      call MAPL_GetResource(MAPL, FDROP_SOOT,     'FDROP_SOOT:',     DEFAULT= 0.05,   __RC__) !Fraction of soot within droplets for immersion freezing        
      call MAPL_GetResource(MAPL, SIGMA_NUC,      'SIGMA_NUC:',      DEFAULT= 1.0,    __RC__) !Widht of the in-cloud distribution of relative humidity in cirrus
      call MAPL_GetResource(MAPL, MIN_ALH,        'MIN_ALH:',        DEFAULT= 5.0,    __RC__) !scale factor for vertical velocity in sttratocumulus
      call MAPL_GetResource(MAPL, SCWST,          'SCWST:',          DEFAULT= 3.0,    __RC__) !scale factor for vertical velocity in sttratocumulus
      call MAPL_GetResource(MAPL, MINCDNC,        'MINCDNC:',        DEFAULT= 0.0,    __RC__) !min nucleated droplet conc. cm-3
      call MAPL_GetResource(MAPL, TMAXCFCORR,     'TMAXCFCORR:',     DEFAULT= 285.0,  __RC__) !Minimum T for CF correction
      call MAPL_GetResource(MAPL, MTIME,          'MTIME:',          DEFAULT= -1.0,   __RC__) !Mixing time scale for aerosol within the cloud. Default is time step
      call MAPL_GetResource(MAPL, SWCIRRUS,       'SWCIRRUS:',       DEFAULT= 3.0,    __RC__) !Tunes vertical velocity in cirrus
      call MAPL_GetResource(MAPL, DUST_INFAC,     'DUST_INFAC:',     DEFAULT= 1.0,    __RC__)  !work on this
      call MAPL_GetResource(MAPL, BC_INFAC,       'BC_INFAC:',       DEFAULT= 0.1,    __RC__)
      call MAPL_GetResource(MAPL, ORG_INFAC,      'ORG_INFAC:',      DEFAULT= 1.0,    __RC__)
      call MAPL_GetResource(MAPL, SS_INFAC,       'SS_INFAC:',       DEFAULT= 1.0,    __RC__)
      call MAPL_GetResource(MAPL, DT_MICRO,       'DT_MICRO:',       DEFAULT= 300.0,  __RC__)    ! time step of the microphysics substepping (s) (MG2) (5 min)
      call MAPL_GetResource(MAPL, UR_SCALE,       'URSCALE:',        DEFAULT= 1.0,    __RC__) !Scaling factor for sed vel of rain    
      call MAPL_GetResource(MAPL, USE_NATURE_WSUB,'USE_NAT_WSUB:',   DEFAULT= 1.0,    __RC__) !greater than zero reads wsub from nature run                     
      call MAPL_GetResource( MAPL, RRTMG_IRRAD ,  'USE_RRTMG_IRRAD:',DEFAULT=0.0,     __RC__)
      call MAPL_GetResource( MAPL, RRTMG_SORAD ,  'USE_RRTMG_SORAD:',DEFAULT=0.0,     __RC__)
      call MAPL_GetResource(MAPL, CLDPARAMS%DISP_FACTOR_LIQ, 'DISP_FACTOR_LIQ:',  DEFAULT= 1.0,   __RC__) ! Scales the droplet/ice crystal number in convective detrainment 
      call MAPL_GetResource(MAPL, CLDPARAMS%DISP_FACTOR_ICE, 'DISP_FACTOR_ICE:',  DEFAULT= 1.0,   __RC__) ! Scales the droplet/ice crystal number in convective detrainment 
      call MAPL_GetResource(MAPL, CLEANUP_RH,                'CLEANUP_RH:',       DEFAULT= 0,     __RC__)
      call MAPL_GetResource(MAPL,GRIDNAME,'AGCM_GRIDNAME:', __RC__)
      call MAPL_GetResource(MAPL, PMIN_CBL,   'PMIN_CBL',   DEFAULT= 50000.0, __RC__)
      call MAPL_GetResource(MAPL,CBL_METHOD,  'CBL_METHOD:', DEFAULT= 6     , __RC__)
      call MAPL_GetResource(MAPL, KSTRAP,  'STRAPPING:',     DEFAULT=-1, __RC__)

      GRIDNAME =  AdjustL(GRIDNAME)
      nn = len_trim(GRIDNAME)
      dateline = GRIDNAME(nn-1:nn)
      imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
      read(imchar,*) imsize
      if(dateline.eq.'CF') imsize = imsize*4

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
    real, pointer, dimension(    :) :: PREF
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS
    real, pointer, dimension(:,:)   :: CNV_FRC, FRLAND, KPBLIN, SH, EVAP
    real, pointer, dimension(:,:)   :: SNOMAS, FRLANDICE
    real, pointer, dimension(:,:,:) :: Q, QRAIN, QSNOW, QGRAUPEL, QLLS, QLCN,     &
                                       CLCN, CLLS, QILS, QICN, NCPL, NCPI, NRAIN, &
                                       NSNOW, NGRAUPEL
    real, pointer, dimension(:,:)   :: TAUGWX, TAUGWY, TAUX, TAUY, TAUOROX, TAUOROY
    real, pointer, dimension(:,:,:) :: OMEGA, ALH, RADLW, RADSW, WSUB_NATURE
    real, pointer, dimension(:,:,:) :: PLE, ZLE, U, V, T, KH, TKE
    real, pointer, dimension(:,:,:) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, &
                                       DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro, &
                                       DTDT_macro,  DUDT_macro,  DVDT_macro, &
                                       DQVDT_micro, DQIDT_micro, DQLDT_micro, &
                                       DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro, &
                                       DTDT_micro,  DUDT_micro,  DVDT_micro, &
                                       SC_ICE, CLDREFFR, CLDREFFS, & 
                                       CLDREFFG, CLDREFFL, CLDREFFI, RAD_CF  , &
                                       RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QV, &
                                       CDNC_NUC, INC_NUC, PFRZ
    real, pointer, dimension(:,:,:) ::  &
       CNV_FICE, CFICE, CFLIQ, DT_RASP, SMAXL, SMAXI, WSUB, CCN01, CCN04, CCN1, &
       NHET_NUC, NLIM_NUC, SO4, ORG, BCARBON, DUST, SEASALT, NCPL_VOL, NCPI_VOL, &
       SAT_RAT, RHICE, RL_MASK, RI_MASK, &
       NHET_IMM, NHET_DEP, DUST_IMM, DUST_DEP, SIGW_GW, SIGW_CNV, SIGW_TURB, &
       SIGW_RC, BERG, BERGS, MELT, DNHET_CT, QCRES, QIRES, AUTICE, FRZPP_LS, &
       SNOWMELT_LS, DNCNUC, DNCSUBL, DNCHMSPLIT, DNCAUTICE, DNCACRIS, DNDCCN, &
       DNDACRLS, DNDACRLR, DNDEVAPC, DNDAUTLIQ, DNDCNV, DNCCNV, &
       PFLCNMOVE, PFICNMOVE, CNV_UPDF, CNV_CVW, DNHET_IMM, CNV_MFD, CNV_DQCDT, &
       CNV_PRC3, SHLW_PRC3, SHLW_SNO3, QLDET_SC, QIDET_SC, MFD_SC, EVAPC, SUBLC
    real, pointer, dimension(:,:,:) :: SC_NDROP, SC_NICE, CUFRC_SC, & 
       RAD_QG, RHCmicro, RHLIQ
    real, pointer, dimension(:,:)   :: EIS, LTS, QCVAR_EXP, &
       CCNCOLUMN, NDCOLUMN, NCCOLUMN, CU2DRAINMOVE, CU2DSNOWMOVE

    real, pointer, dimension(:,:,:) :: RHX, REV_AN, RSU_AN, REV_LS, RSU_LS, &
                                            PFL_AN, PFI_AN, PFL_LS, PFI_LS
    real, pointer, dimension(:,:,:) :: PTR3D

    real, pointer, dimension(:,:)   :: CN_PRCP, CN_SNR, ER_PRCP, LS_PRCP, LS_SNR, CN_ARF, LS_ARF
    real, pointer, dimension(:,:)   :: PTR2D

    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D

    integer, allocatable, dimension(:, :)   :: KCT, KMIN_TROP, KLCL, KPBL, KCBL
    real, allocatable, dimension(:, :)  :: NPRE_FRAC_2d
    real, allocatable, dimension(:,:,:) :: CNV_PLE, PLO, TH1, Q1, U1, V1, TEMP, PK
    real, allocatable, dimension(:,:,:) :: QSNOW_AN, QRAIN_AN, DP, MASS, DQS, QSS
    real, allocatable, dimension(:,:,:) :: QCNTOT, CFX
    real (ESMF_KIND_R8)  :: tauxr8, fsoot_drop, fdust_drop, sigma_nuc_r8, rh1_r8, &
                            frachet_dust, frachet_bc, frachet_org, frachet_ss

    real, allocatable, dimension(:,:)   :: ZWS, ZPBL

  real, allocatable, dimension(:,:,:)  ::  QTOT, QL_TOT, QI_TOT
  real, allocatable, dimension(:,:,:) :: DQST3, QST3, DZET, QDDF3, FQAI, &
                     FQAL, FQA, ZLO, GZLO, DQSDT, ZLE0, ZL0

  ! Manage diagnostic outputs for accretion
  !---------------------------------------------------
  real, allocatable, dimension(:,:,:) :: ACIL_LS_X, ACIL_AN_X, ACLL_LS_X, ACLL_AN_X

  ! Manage diagnostic outputs for 3D precip fluxes
  !---------------------------------------------------
  real, allocatable, dimension(:,:,:) :: DLPDF_X, DIPDF_X, DLFIX_X, DIFIX_X, &
       AUT_X, SDM_X, FRZ_TT_X, FRZ_PP_X, DCNVL_X, DCNVI_X

  real, allocatable, dimension(:,:) ::  CLDREFFI_TOP_X, CLDREFFL_TOP_X,  NCPL_TOP_X, NCPI_TOP_X, NCPL_CLDBASEX


  real, allocatable, dimension(:,:,:) :: ALPHT_X

  real, allocatable, dimension(:,:,:) :: TH

  real, allocatable, dimension(:,:,:) :: VFALLRN_AN_X, VFALLSN_AN_X

    integer :: I, J, L, K

    real, parameter :: pmin_trop = 10.0 !mbar minimum pressure to do cloud microphysics
    logical                   :: use_average_v
    real (ESMF_KIND_R8), dimension(3)       :: ccn_diag
    real(ESMF_KIND_R8), allocatable, dimension(:,:)  :: ttendr8, qtendr8, cwtendr8, &
           cldor8,  rpdelr8, zmr8, omegr8, rhdfdar8, rhu00r8, ficer8 , &
           ndropr8, nimmr8
    real (ESMF_KIND_R8), allocatable, dimension(:,:) :: wparc, smaxliq, atot, &
           smaxicer8, nheticer8, incr8, swparc, &
           nhetr8, nlimicer8, qilsr8, wparc_gw, wparc_ls, &
           wparc_turb, wparc_cnv, lc_turb, rad_cooling, wparc_rc, &
           uwind_gw, wparc_cgw, pfrz_inc_r8
    real(ESMF_KIND_R8), dimension(1)       :: prectr8, precir8
    real(ESMF_KIND_R8)  :: disp_liu, ui_scale, dcrit, tfreez, qcvar8, &
                           ts_autice, dcsr8, qcvarr8, scale_ri, mtimesc, urscale
    integer :: num_steps_micro,  pcnst, n_modes, kbmin, kcldtop, kcldbot , &
                NAUX, kcldtopcvn, nbincontactdust, index
    real(ESMF_KIND_R8), allocatable, dimension(:,:)  :: pintr8, kkvhr8, rflxr8, &
                                                    sflxr8, lflxr8, iflxr8, gflxr8
    real(ESMF_KIND_R8), allocatable, dimension(:,:,:) :: rndstr8,naconr8  !Assume maximum 5 dust bins
    REAL, allocatable, dimension(:,:) :: SCICE_tmp, FQA_tmp, ALPH_tmp

    real(ESMF_KIND_R8), allocatable, dimension(:,:) ::so4x, seasaltx, dustx, &
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

      real, allocatable, dimension (:, :) :: tm_gw, pm_gw, nm_gw, theta_tr,  &
            fcn, cfaux
      real, allocatable, dimension (:, :) :: pi_gw, rhoi_gw, ni_gw, ti_gw
      real   :: maxkhpbl, tausurf_gw, fracover, cfc_aux
      real   :: aux1,aux2,aux3,hfs,hfl
      real(ESMF_KIND_R8)   :: autscx
      real   :: Nct, Wct, ksa1, Xscale
      integer :: K0, KCBLMIN
      integer  :: i_src_mode
      integer  :: i_dst_mode
      real :: USURF, RHEXCESS
      real, parameter :: r_air = 3.47d-3 !m3 Pa kg-1K-1
      integer,  parameter :: ncolmicro = 1

      type  (AerProps) :: AeroAux, AeroAux_b

      character(len=ESMF_MAXSTR)              :: IAm

    IAm = "MGB2_2M_Run"

    call ESMF_GridCompGet( GC, CONFIG=CF, __RC__ ) 
    

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__)
    

    call MAPL_TimerOn (MAPL,"--MGB2_2M",__RC__)

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_Get( MAPL, IM=IM, JM=JM, LM=LM,   &
         RUNALARM = ALARM,             &
         CF       = CF,                &
         LONS     = LONS,              &
         LATS     = LATS,              &
         INTERNAL_ESMF_STATE=INTERNAL, &
         __RC__ )
    

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
    allocate(ALPH_tmp(1,LM), __STAT__)
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

    allocate(KMIN_TROP(IM,JM), __STAT__)
    allocate(KLCL(IM,JM), __STAT__)
    allocate(KPBL(IM,JM), __STAT__)
    allocate(KCBL(IM,JM), __STAT__)
    allocate(NPRE_FRAC_2d(IM,JM), __STAT__)
    allocate(ZWS(IM,JM), __STAT__)
    allocate(ZPBL(IM,JM), __STAT__)
    allocate(KCT(IM,JM), __STAT__)

    allocate(QSS(IM,JM,LM ), __STAT__)
    allocate(DQSDT(IM,JM,LM ), __STAT__)
    allocate(DQST3(IM,JM,LM ), __STAT__)
    allocate(QST3(IM,JM,LM ), __STAT__)
    allocate(DZET(IM,JM,LM ), __STAT__)
    allocate(QDDF3(IM,JM,LM ), __STAT__)
    allocate(FQAI(IM,JM,LM ), __STAT__)
    allocate(FQAL(IM,JM,LM ), __STAT__)
    allocate(FQA(IM,JM,LM ), __STAT__)
    allocate(ZLO(IM,JM,LM ), __STAT__)
    allocate(GZLO(IM,JM,LM ), __STAT__)
    allocate(PLO(IM,JM,LM ), __STAT__)
    allocate(CNV_PLE(IM,JM,0:LM), __STAT__)
    allocate(TH1(IM,JM,LM ), __STAT__)
    allocate(TEMP(IM,JM,LM ), __STAT__)
    allocate(PK(IM,JM,LM ), __STAT__)
    allocate(Q1(IM,JM,LM ), __STAT__)
    allocate(QSNOW_AN(IM,JM,LM ), __STAT__)
    allocate(QRAIN_AN(IM,JM,LM ), __STAT__)
    allocate(DP(IM,JM,LM ), __STAT__)
    allocate(MASS(IM,JM,LM ), __STAT__)
    allocate(DQS(IM,JM,LM ), __STAT__)
    allocate(QCNTOT(IM,JM,LM), __STAT__)
    allocate(CFX(IM,JM,LM), __STAT__)
    allocate(U1(IM,JM,LM), __STAT__)
    allocate(V1(IM,JM,LM), __STAT__)

    allocate(QTOT(IM,JM,LM ), __STAT__)
    allocate(QL_TOT(IM,JM,LM ), __STAT__)
    allocate(QI_TOT(IM,JM,LM ), __STAT__)
    allocate(ACIL_AN_X(IM,JM,LM ), __STAT__)
    allocate(ACIL_LS_X(IM,JM,LM ), __STAT__)
    allocate(ACLL_AN_X(IM,JM,LM ), __STAT__)
    allocate(ACLL_LS_X(IM,JM,LM ), __STAT__)
    allocate(ZLE0(IM,JM,0:LM ), __STAT__)
    allocate(ZL0(IM,JM,LM ), __STAT__)
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
    allocate(ALPHT_X(IM,JM,LM ), __STAT__)
    allocate(VFALLRN_AN_X(IM,JM,LM ), __STAT__)
    allocate(VFALLSN_AN_X(IM,JM,LM ), __STAT__)
    allocate(TH(IM,JM,LM ), __STAT__)

    call ESMF_AlarmGet(ALARM, RingInterval=TINT, __RC__)
    call ESMF_TimeIntervalGet(TINT,   S_R8=DT_R8,__RC__)
    DT_MOIST = DT_R8

!#ifdef NODISABLE

    allocate( TMP3D(IM,JM,LM) )
    allocate( TMP2D(IM,JM)    )

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , __RC__)
    call MAPL_GetPointer(INTERNAL, QRAIN,    'QRAIN'   , __RC__)
    call MAPL_GetPointer(INTERNAL, QSNOW,    'QSNOW'   , __RC__)
    call MAPL_GetPointer(INTERNAL, QGRAUPEL, 'QGRAUPEL', __RC__)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , __RC__)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , __RC__)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , __RC__)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'    , __RC__)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL, 'NGRAUPEL'    , __RC__)

    call MAPL_GetPointer(IMPORT, SNOMAS, 'SNOMAS'     , __RC__)
    call MAPL_GetPointer(IMPORT, FRLANDICE, 'FRLANDICE'     , __RC__)
    call MAPL_GetPointer(IMPORT, FRLAND, 'FRLAND'     , __RC__)
    call MAPL_GetPointer(IMPORT, PREF, 'PREF'     , __RC__)
    call MAPL_GetPointer(IMPORT, KPBLIN, 'KPBL'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUGWX, 'TAUGWX'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUGWY, 'TAUGWY'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUX,   'TAUX'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUY,   'TAUY'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUOROX, 'TAUOROX'     , __RC__)
    call MAPL_GetPointer(IMPORT, TAUOROY, 'TAUOROY'     , __RC__)
    call MAPL_GetPointer(IMPORT, OMEGA,  'OMEGA'     , __RC__)
    call MAPL_GetPointer(IMPORT, ALH,    'ALH'     , __RC__)
    call MAPL_GetPointer(IMPORT, RADLW,  'RADLW'     , __RC__)
    call MAPL_GetPointer(IMPORT, RADSW,  'RADSW'     , __RC__)
    call MAPL_GetPointer(IMPORT, WSUB_NATURE,  'WSUB_NATURE'     , __RC__)
    call MAPL_GetPointer(IMPORT, PLE,  'PLE'     , __RC__)
    call MAPL_GetPointer(IMPORT, U,    'U'       , __RC__)
    call MAPL_GetPointer(IMPORT, V,    'V'       , __RC__)
    call MAPL_GetPointer(IMPORT, T,    'T'       , __RC__)
    call MAPL_GetPointer(IMPORT, TKE,    'TKE'       , __RC__)
    call MAPL_GetPointer(IMPORT, ZLE,  'ZLE'     , __RC__)
    call MAPL_GetPointer(IMPORT, KH,   'KH'      , __RC__)
    call MAPL_GetPointer(IMPORT, SH,  'SH'     , __RC__)
    call MAPL_GetPointer(IMPORT, EVAP,  'EVAP'     , __RC__)

    call MAPL_GetPointer(EXPORT, DTDT_micro, 'DTDT_micro' ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DQVDT_micro, 'DQVDT_micro' ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DQIDT_micro, 'DQIDT_micro' ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DQLDT_micro, 'DQLDT_micro' ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CFICE,   'CFICE'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CFLIQ,   'CFLIQ'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DT_RASP,   'DT_RASP'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_FICE,   'CNV_FICE'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CCNCOLUMN,   'CCNCOLUMN'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NDCOLUMN,   'NDCOLUMN'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, NCCOLUMN,   'NCCOLUMN'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RHLIQ,   'RHLIQ'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RHCmicro,   'RHCmicro'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QCVAR_EXP,   'QCVAR_EXP'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SC_NDROP,   'SC_NDROP'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SC_NICE,   'SC_NICE'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, MFD_SC,   'MFD_SC'   ,  ALLOC=.TRUE., __RC__)

    call MAPL_GetPointer(EXPORT, SC_ICE,      'SC_ICE'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFL,    'RL'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFI,    'RI'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFR,    'RR'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFS,    'RS'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFG,    'RG'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_CF,      'FCLD'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_QV,      'QV'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_QL,      'QL'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_QI,      'QI'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_QR,      'QR'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_QS,      'QS'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, RAD_QG,   'QG'   ,  ALLOC=.TRUE., __RC__)
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
    call MAPL_GetPointer(EXPORT, DNCCNV,      'DNCCNV'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_FRC,     'CNV_FRC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_UPDF,     'CNV_UPDF'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_CVW,     'CNV_CVW'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, DNHET_IMM,   'DNHET_IMM'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_MFD,   'CNV_MFD'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_DQCDT, 'CNV_DQCDT'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CNV_PRC3,   'CNV_PRC3'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SHLW_PRC3,   'SHLW_PRC3'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SHLW_SNO3,   'SHLW_SNO3'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QLDET_SC,   'QLDET_SC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, QIDET_SC,   'QIDET_SC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CUFRC_SC,   'CUFRC_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EVAPC,   'EVAPC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SUBLC,   'SUBLC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

    call init_Aer(AeroAux)
    call init_Aer(AeroAux_b)

   !============================================= Start Stratiform cloud processes==========================================

    DO J=1, JM
       DO I=1, IM
          if (KPBLIN(I,J) == 0) KPBLIN(I,J) = LM-1
       END DO
    END DO

      ! MAT The code below should use nint(KPBLIN)
      DO J=1, JM
         DO I=1, IM
            ZPBL(I, J) =ZLE(I, J,  NINT(KPBLIN(I, J )))
         END DO
      END DO

         ! find the minimun level for cloud micro calculations
         CNV_PLE  = PLE*.01
         PLO      = 0.5*(CNV_PLE(:,:,0:LM-1) +  CNV_PLE(:,:,1:LM  ) )
         call find_l(KMIN_TROP, PLO, pmin_trop, IM, JM, LM, 10, LM-2)
         ! Find Convective Cloud Top
         KCT = 20 !default upper limit. Less than 20 makes no difference
         call find_l(KCT, CNV_DQCDT, 1.0e-9, IM, JM, LM, 20, LM-2)

         U1       = U
         V1       = V
         Q1       = Q
         PK       = (100.0*PLO/MAPL_P00)**(MAPL_KAPPA)
         TH       = T/PK
         TH1      = TH
         DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
         MASS     = DP/MAPL_GRAV
         TEMP    = TH1*PK   

        DO K=0,LM
           ZLE0(:,:,K)= ZLE(:,:,K) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
        END DO
        ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
        DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
        TH       = T/PK
        DQST3    = GEOS_DQSAT(T, PLO, QSAT=QST3)
        WHERE ( ZL0 < 3000. )
           QDDF3 = -( ZL0-3000. ) * ZL0 * MASS
        ELSEWHERE
           QDDF3 = 0.
        END WHERE
        TMP2D = SUM(QDDF3, 3)
        DO K = 1,LM
           QDDF3(:,:,K) = QDDF3(:,:,K) / TMP2D
        END DO

        ZLO = ZL0
        GZLO = MAPL_GRAV*ZL0

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
         QSNOW_AN = 0.0
         QRAIN_AN = 0.0
         PFRZ= 0.0

      K0 = LM
      KCBLMIN  =       count(PREF < PMIN_CBL)

      ! Find estimated inversion strength (DONIF)

       DQS      = GEOS_DQSAT(TEMP, PLO, qsat=QSS)
       KLCL     = FINDLCL( TH1, Q1, PLO, PK, IM, JM, LM )
      !!    KPBL = FINDPBL( KH, IM, JM, LM )
      !! Set subcloud layer height to one level below PBL height level
      !!   make sure subcloud layer is at least 2 levels thick
      do j = 1,jm
         do i = 1,im
            if(nint(KPBLIN(i,j)).ne.0) then
               KPBL(i,j) = max(min(nint(KPBLIN(i,j))+1,LM-1), 1)
            else
               KPBL(i,j) = LM-1
            endif
         enddo
      enddo

      do J=1,JM
         do I=1,IM

            SELECT CASE( CBL_METHOD )

            CASE( 1 )
               KCBL(I,J)   =  K0 - KSTRAP
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )

            CASE( 2 )
               KCBL(I,J)   = KLCL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )

            CASE ( 3 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )

            CASE ( 4 )
               KCBL(I,J)   = KLCL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )

            CASE ( 5 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )

            CASE( 6 )
               KCBL(I,J)   = KPBL(I,J)
               KCBL(I,J)   =  MAX( KCBL(I,J), KCBLMIN )

            END SELECT

         end do
      end do

       call    FIND_EIS(TH1, QSS, TEMP, ZLE, PLO, KLCL, IM, JM, LM, LTS, EIS)

         ! Clean up any negative specific humidity before the microphysics scheme
         !-----------------------------------------
         call FILLQ2ZERO( Q1, MASS, TMP2D)  

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
             
                     pi_gw(1, 0:LM) = 100.0*CNV_PLE(I,J,0:LM)                     
                     theta_tr(1,1:LM) = TH1(I,J,1:LM)
                     rhoi_gw = 0.0  
                     ni_gw = 0.0 
                     ti_gw = 0.0                                         
                     ter8(1,1:LM) = TEMP(I,J,1:LM)   
                     pi_gw(1, 0:LM) = 100.0*CNV_PLE(I,J,0:LM) 
                     plevr8(1,1:LM) = 100.*PLO(I,J,:)
                     ndropr8(1,1:LM) = NCPL(I, J, 1:LM)
                     qir8(1,1:LM) =  QILS(I, J,1:LM)+QICN(I, J,1:LM)
                     qcr8(1,1:LM) =  QLLS(I, J,1:LM)+QLCN(I, J,1:LM)
                     npre8(1,1:LM) = NPRE_FRAC_2d(I,J)*NCPI(I,J,1:LM)
                     omegr8(1,1:LM) = OMEGA(I,J,1:LM)  
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
                    
               ! ==========================================================================================    
               ! ========================Activate the aerosols ============================================ 
           
               
               
                do K = KMIN_TROP(I, J), LM-1 !limit to troposphere and no activation at the surface
    
                ! find vertical velocity variance 
                       call vertical_vel_variance(omegr8(1, K), lc_turb(1, K), ter8(1, K), plevr8(1, K), rad_cooling(1,K),  uwind_gw(1,K), &
                                                         tausurf_gw, nm_gw(1, K), LCCIRRUS, Nct, Wct, &
                                                         ksa1, fcn(1, K), KH(I, J, K), FRLAND(I, J), ZPBL(I, J), ZLE(I, J, k), maxkhpbl, &
                                                            wparc_ls(1, K), wparc_gw(1, K), wparc_cgw(1, K), wparc_turb(1, K), EIS(I, J), TKE(I, J, K))
                                        
                   if (FRLAND(I, J) .lt. 0.1) then 
                    if (LTS(I, J) .gt. LTS_LOW) then                     
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
                        tauxr8 = ter8(1, K)
                             
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

    
        
         !=============================================End cloud particle nucleation=====================================
         !===============================================================================================================

         call MAPL_TimerOff(MAPL,"---ACTIV", __RC__)

    ! Exports to be filled 
    call MAPL_GetPointer(EXPORT, ER_PRCP,  'ER_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CN_PRCP,  'CN_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CN_SNR,   'CN_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, LS_PRCP,  'LS_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, LS_SNR,   'LS_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CN_ARF,   'CN_ARF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, LS_ARF,   'LS_ARF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RHX   ,   'RHX'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_AN,   'REV_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_AN,   'RSU_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_LS,   'REV_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_LS,   'RSU_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_AN,   'PFL_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_AN,   'PFI_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_LS,   'PFL_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_LS,   'PFI_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Unused Exports (foreced to 0.0)
    call MAPL_GetPointer(EXPORT, PTR2D,  'AN_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
    call MAPL_GetPointer(EXPORT, PTR2D,  'SC_PRCP'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS); PTR2D=0.0
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
    if (associated( DUDT_macro))  DUDT_macro=U1
    if (associated( DVDT_macro))  DVDT_macro=V1
    if (associated( DTDT_macro))  DTDT_macro=TEMP
    if (associated(DQVDT_macro)) DQVDT_macro=Q
    if (associated(DQLDT_macro)) DQLDT_macro=QLCN+QLLS
    if (associated(DQIDT_macro)) DQIDT_macro=QICN+QILS
    if (associated(DQADT_macro)) DQADT_macro=CLCN+CLLS
    if (associated(DQRDT_macro)) DQRDT_macro=QRAIN
    if (associated(DQSDT_macro)) DQSDT_macro=QSNOW
    if (associated(DQGDT_macro)) DQGDT_macro=QGRAUPEL

         !==========================================================================================================
         !===================================Cloud Macrophysics ====================================================
         !==========================================================================================================

         CFX=INC_NUC + NHET_IMM 
      
  call  macro_cloud (                    &
              IM*JM, LM         , &
              DT_MOIST          , &
              PLO               , &
              CNV_PLE           , &
              PK                , &
              FRLAND            , &   ! <- surf
              CNV_FRC           , &   ! <- convective fraction
              CNV_FRC           , &   ! <- convective fraction
              CNV_DQCDT         , &   ! <- dpcu              
              CNV_PRC3          , &   ! <- dpcu   
              CNV_UPDF          , &   ! <- dpcu
              QLDET_SC          , &   ! <- shcu   
              QIDET_SC          , &   ! <- shcu   
              SHLW_PRC3         , &   ! <- shcu   
              SHLW_SNO3         , &   ! <- shcu   
              CUFRC_SC          , &   ! <- shcu 
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
              CN_PRCP           , &            
              CN_ARF            , &
              CN_SNR            , &
              QST3              , &
              DZET              , &
              QDDF3             , &
                                ! Diagnostics
              RHX               , &
              REV_AN            , &
              RSU_AN            , &
              ACLL_AN_X,ACIL_AN_X   , &
              PFL_AN,PFI_AN     , &
              DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
              DCNVL_X, DCNVI_X,       &
              ALPHT_X, &
              VFALLSN_AN_X,  &
              VFALLRN_AN_X,  &
              EVAPC , SUBLC,  &
                                ! End diagnostics
            !!====2-Moment============
              SC_ICE,    &
              NCPL, &
              NCPI, &
              PFRZ, &
              DNDCNV, &
              DNCCNV, &
              DT_RASP , &
              QRAIN_AN, & !grid av
              QSNOW_AN, &
              KCBL)

         TEMP    = TH1*PK

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
          call FILLQ2ZERO( Q1, MASS, TMP2D) 
          call FILLQ2ZERO( QGRAUPEL, MASS, TMP2D) 
          call FILLQ2ZERO( QRAIN, MASS, TMP2D) 
          call FILLQ2ZERO( QSNOW, MASS, TMP2D) 
          call FILLQ2ZERO( QLLS, MASS, TMP2D)
          call FILLQ2ZERO( QLCN, MASS, TMP2D)  
          call FILLQ2ZERO( QILS, MASS, TMP2D)
          call FILLQ2ZERO( QICN, MASS, TMP2D)
         
         !=============================================End cloud macrophysics=====================================
         !======================================================================================================================
         !

         FQAI = 0.0
         FQAL = 0.0
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

       ! Update TH
       TH1 = TEMP/PK

    ! Update macrophysics tendencies
    if (associated( DUDT_macro))  DUDT_macro=( U1        - DUDT_macro)/DT_MOIST
    if (associated( DVDT_macro))  DVDT_macro=( V1        - DVDT_macro)/DT_MOIST
    if (associated( DTDT_macro))  DTDT_macro=( TEMP      - DTDT_macro)/DT_MOIST
    if (associated(DQVDT_macro)) DQVDT_macro=( Q1        -DQVDT_macro)/DT_MOIST
    if (associated(DQLDT_macro)) DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
    if (associated(DQIDT_macro)) DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
    if (associated(DQADT_macro)) DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
    if (associated(DQRDT_macro)) DQRDT_macro=( QRAIN     -DQRDT_macro)/DT_MOIST
    if (associated(DQSDT_macro)) DQSDT_macro=( QSNOW     -DQSDT_macro)/DT_MOIST
    if (associated(DQGDT_macro)) DQGDT_macro=( QGRAUPEL  -DQGDT_macro)/DT_MOIST
    call MAPL_TimerOff(MAPL,"---CLDMACRO")

    call MAPL_TimerOn (MAPL,"---CLDMICRO", __RC__)
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
    if (associated(DQVDT_micro)) DQVDT_micro = Q1
    if (associated(DQLDT_micro)) DQLDT_micro = QLLS + QLCN
    if (associated(DQIDT_micro)) DQIDT_micro = QILS + QICN
    if (associated(DQRDT_micro)) DQRDT_micro = QRAIN
    if (associated(DQSDT_micro)) DQSDT_micro = QSNOW
    if (associated(DQGDT_micro)) DQGDT_micro = QGRAUPEL
    if (associated(DQADT_micro)) DQADT_micro = CLLS + CLCN
    if (associated( DUDT_micro))  DUDT_micro = U1
    if (associated( DVDT_micro))  DVDT_micro = V1
    if (associated( DTDT_micro))  DTDT_micro = TEMP

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
     
         accre_enhanr8= ACC_ENH
         accre_enhan_icer8= ACC_ENH_ICE
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

                  cldfr8(1,1:LM)  = RAD_CF(I,J,1:LM) !Assume minimum overlap 
               liqcldfr8(1,1:LM)  =  CFLIQ(I,J,1:LM) 
               icecldfr8(1,1:LM)  =  CFICE(I,J,1:LM) 
                 
                  cldor8          = cldfr8  
               ter8(1,1:LM)       = TEMP(I,J,1:LM)
               qvr8(1,1:LM)       =   Q1(I,J,1:LM)

               qcr8(1,1:LM)        =     QL_TOT(I,J,1:LM)
               qir8(1,1:LM)        =     QI_TOT(I,J,1:LM)
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
                 ELSE
                    DO K=LM, 1, -1       
                         if (PLO(I,J,K) .lt. 500.0) exit  
                         HSMOIST_500 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + QST3(I, J, K)*MAPL_ALHL
                    END DO 
                    DO K=LM, 1, -1       
                         if (PLO(I,J,K) .lt. 950.0) exit  
                         HMOIST_950 = MAPL_CP*TEMP(I, J, K) + GZLO(I, J, K) + Q1(I, J, K)*MAPL_ALHL
                    END DO                                          
                    SINST = (HMOIST_950 -  HSMOIST_500)/45000.0                  
                  ENDIF
               
                  xscale = (36000.0/imsize)**(-0.666)
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
                  ALPH_tmp(1,1:LM)  = ALPHT_X(I, J, 1:LM)
                   
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
        
         call  micro_mg_tend_interface ( DT_MICRO, INT(CLDPARAMS%PDFSHAPE), ALPH_tmp, SCICE_tmp, FQA_tmp, &
                             ncolmicro,             LM,               dt_r8,       & 
                             CNV_FRC(I,J), CNV_FRC(I,J), &
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
                             prer_evap, &
                             frzimmr8,             frzcntr8,              frzdepr8,  & ! contact is not passed since it depends on the droplet size dist
                             nsootr8, rnsootr8,  & ! soot for contact IN
                             npccnor8, npsacwsor8,npraor8,nsubcor8, nprc1or8, &  ! Number tendencies for liquid
                             npraior8, nnucctor8, nnucccor8, nnuccdor8, nsubior8, nprcior8, nsacwior8,  &  ! Number tendencies for ice
                             ts_autice, ui_scale, autscx , disp_liu, nbincontactdust, urscale)

    end if 

        IF (MGVERSION > 1) then 

                   QRAIN   (I,J,1:LM)  = max(REAL(qroutr8(1,1:LM)), 0.0)
                   QSNOW   (I,J,1:LM)  = max(REAL(qsoutr8(1,1:LM)), 0.0)
                   QGRAUPEL(I,J,1:LM)  = max(REAL(qgoutr8(1,1:LM)), 0.0)
                   NRAIN   (I,J,1:LM)  = max(REAL(nroutr8(1,1:LM)), 0.0)
                   NSNOW   (I,J,1:LM)  = max(REAL(nsoutr8(1,1:LM)), 0.0)
                   NGRAUPEL(I,J,1:LM)  = max(REAL(ngoutr8(1,1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM)  = REAL(reff_rainr8(1,1:LM))
                   CLDREFFS(I,J,1:LM)  = REAL(reff_snowr8(1,1:LM))/scale_ri
                   CLDREFFG(I,J,1:LM)  = REAL(reff_graur8(1,1:LM))/scale_ri
            
        else
                    
                   QRAIN(I,J,1:LM)  = max(REAL(qrout2r8(1,1:LM)), 0.0) ! grid average 
                   QSNOW(I,J,1:LM)  = max(REAL(qsout2r8(1,1:LM)), 0.0)                      
                   NRAIN(I,J,1:LM)  = max(REAL(nrout2r8(1,1:LM)), 0.0)
                   NSNOW(I,J,1:LM)  = max(REAL(nsout2r8(1,1:LM)), 0.0)
                   CLDREFFR(I,J,1:LM) = REAL(drout2r8(1,1:LM))/2.0        
                   CLDREFFS(I,J,1:LM) = REAL(dsout2r8(1,1:LM))/2.0/scale_ri
                 
         end if          
         
               PFL_LS(I, J, 1:LM) = rflxr8(1,1:LM) !+ lflxr8(1,1:LM)
               PFI_LS(I, J, 1:LM) = sflxr8(1,1:LM) !+ gflxr8(1,1:LM) +  iflxr8(1,1:LM)
              
               !Update state after microphysisc
               LS_PRCP(I,J)     = max(1000.*REAL((prectr8(1)-precir8(1))), 0.0)
               LS_SNR(I,J)      = max(1000.*REAL(precir8(1)), 0.0)          
               QL_TOT(I,J,1:LM) = max(QL_TOT(I,J,1:LM)   + REAL(qctendr8(1,1:LM)) * DT_R8, 0.0)
               QI_TOT(I,J,1:LM) = max(QI_TOT(I,J,1:LM)   + REAL(qitendr8(1,1:LM)) * DT_R8, 0.0)    
               Q1(I,J,1:LM)   = MAX(Q1(I,J,1:LM)     + REAL(qvlatr8(1,1:LM)) * DT_R8, 0.0)
               TEMP(I,J,1:LM) = TEMP(I,J,1:LM)   + REAL(tlatr8(1,1:LM)) * DT_R8 / (MAPL_CP)  
               NCPL(I,J,1:LM) = MAX(NCPL(I,J,1:LM)   + REAL(nctendr8(1,1:LM)) * DT_R8, 0.0) 
               NCPI(I,J,1:LM) = MAX(NCPI(I,J,1:LM)   + REAL(nitendr8(1,1:LM)) * DT_R8, 0.0)  

               LS_ARF(I,J)     = maxval( REAL(cldfr8(1,1:LM)) )
                            
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

    DLPDF_X=  QLLS +QLCN
    DIPDF_X=  QILS +QICN
    do K= 1, LM
       do J=1,JM
            do I=1,IM
                  call update_cld( &
                         DT_MOIST                , &
                         ALPHT_X(I, J, K)        , &
                         INT(CLDPARAMS%PDFSHAPE) , &
                         CNV_FRC(I, J)           , &
                         CNV_FRC(I, J)           , &
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
                         RHCmicro(I, J, K), &
                         .TRUE.)
              
           end do 
       end do
    end do 
    DLPDF_X=((QLLS+QLCN) - DLPDF_X)/DT_MOIST
    DIPDF_X=((QILS+QICN) - DIPDF_X)/DT_MOIST

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
        !to m-3 s-1
         DNHET_CT    = DNHET_CT*CFX 
         DNHET_IMM   = DNHET_IMM*CFX 
         DNCNUC      = DNCNUC*CFX 
         DNCHMSPLIT  = DNCHMSPLIT*CFX
         DNCSUBL     = DNCSUBL*CFX 
         DNCACRIS    = DNCACRIS*CFX 
         DNCAUTICE   = DNCAUTICE*CFX 
         DNCCNV      = DNCCNV*CFX

         DNDCCN       = DNDCCN*CFX   
         DNDACRLS     = DNDACRLS*CFX
         DNDACRLR     = DNDACRLR*CFX    
         DNDEVAPC     = DNDEVAPC*CFX   
         DNDAUTLIQ    = DNDAUTLIQ*CFX 
         DNDCNV       = DNDCNV*CFX

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

       ! !Set rain water for radiation to 0 if preciprad flag is off (set to 0)
       ! if(CLDPARAMS%PRECIPRAD .eq. 0.) then
       !    RAD_QR = 0.
       !    RAD_QS = 0.
       !    RAD_QG = 0.      
       ! endif

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

      !=====Tune area cloud fraction of extended PBL clouds. Area cloud fraction may be different from Volume cloud fraction 
         FQA= 0.0
         where (RAD_CF .gt. 0.0)
              FQA =  CLCN/RAD_CF
         end where

#ifdef SKIP 
         DO I=1, IM
            DO J = 1, JM    

               if (FRLAND(I, J) .lt. 0.1) then 

                  
                  cfc_aux =(TEMP(I, J, LM) - TMAXCFCORR)/2.0
                  cfc_aux =  min(max(cfc_aux,-20.0), 20.0)
                  cfc_aux=   1.0/(1.0+exp(-cfc_aux))
                  
                  DO K=LM-1, 2, -1
                     if ((RAD_CF(I, J, K) .gt. 0.01) .and. (RAD_CF(I, J, K) .lt. 0.99)) then  

                           USURF=1.0                        
                           USURF= (LTS_UP-LTS(I, J))/(LTS_UP - LTS_LOW) 
                           USURF=min(max(USURF, MIN_EXP), MAX_EXP)                            
                           
                           fracover=min(max((TEMP(I, J, K) -TMAXLL)/2.0, -20.0), 20.0)
                           fracover = 1.0/(1.0+exp(-fracover))
                          
                           USURF = USURF*fracover + 1.0-fracover   !only near the surface                                  
                           USURF = USURF*cfc_aux + 1.0-cfc_aux !only for the subtropics                          
                           USURF =  usurf+ (1.0-USURF)*CNV_FRC(I, J) !only non-convective
                                            
                           RAD_CF(I, J, K)=RAD_CF(I, J, K)**USURF                                                                
                     END IF

                  END DO
               END IF

            end do
         end do
#endif
         
      CLCN =       FQA *RAD_CF
      CLLS =  (1.0-FQA)*RAD_CF   

      WHERE (QTOT .gt. 1.0e-12) 
           CFLIQ=RAD_CF*QL_TOT/QTOT
           CFICE=RAD_CF*QI_TOT/QTOT
      END WHERE
         
      ! Clean up Relative Humidity where RH > 110%
      !---------------------------------------------
      if (CLEANUP_RH == 1)  then

         RHEXCESS = 1.1

         !-----If using 2-moment microphysics, allow for supersaturation w.r.t ice ---

         QSS=GEOS_QsatICE (TH1*PK, PLO*100.0)
         where (CFICE .lt. 0.99 .and. QSS .gt. 1.0e-20)
            CFX = max(( Q1 - QSS*CFICE), 0.0)/(1.0-CFICE)
            SAT_RAT = min(CFX/QSS, 2.0) ! this is just diagnostic
         elsewhere
            SAT_RAT = 1.0
         end where

         if (associated(RHICE    )) RHICE  =  Q1/QSS
         where (TEMP>273.15)
            RHICE=0.0
         end where

         QSS  = GEOS_QsatLQU(         &
              TEMP   , &
              PLO*100.0 , DQ=DQSDT     ) !clean up only with respect to liquid water  DONIF
         if (associated(RHLIQ    ))   RHLIQ     =  Q1/QSS  !DONIF

         where ( Q1 > RHEXCESS*QSS )
            DQS = (Q1 - RHEXCESS*QSS)/( 1.0 + RHEXCESS*DQSDT*MAPL_ALHL/MAPL_CP )
         elsewhere
            DQS = 0.0
         endwhere

         Q1      =  Q1    - DQS
         TEMP    =  TEMP  + (MAPL_ALHL/MAPL_CP)*DQS
         TH1     =  TEMP  / PK

         ER_PRCP =  SUM( DQS * MASS, 3)/DT_MOIST
         LS_PRCP =  LS_PRCP + ER_PRCP
      ELSE
         ER_PRCP =  0.0
      ENDIF

      if (associated(CCNCOLUMN))   CCNCOLUMN  = SUM(    CCN1*MASS/(100.*PLO*r_air/TEMP) , 3)
      if (associated(NDCOLUMN ))    NDCOLUMN  = SUM(NCPL_VOL*MASS/(100.*PLO*r_air/TEMP) , 3)
      if (associated(NCCOLUMN ))    NCCOLUMN  = SUM(NCPI_VOL*MASS/(100.*PLO*r_air/TEMP) , 3)

      ! Update microphysics tendencies
      if (associated(DQVDT_micro)) DQVDT_micro = ( Q1         - DQVDT_micro) / DT_MOIST
      if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
      if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
      if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
      if (associated(DQRDT_micro)) DQRDT_micro = ( QRAIN      - DQRDT_micro) / DT_MOIST
      if (associated(DQSDT_micro)) DQSDT_micro = ( QSNOW      - DQSDT_micro) / DT_MOIST
      if (associated(DQGDT_micro)) DQGDT_micro = ( QGRAUPEL   - DQGDT_micro) / DT_MOIST
      if (associated( DUDT_micro))  DUDT_micro = ( U1         -  DUDT_micro) / DT_MOIST
      if (associated( DVDT_micro))  DVDT_micro = ( V1         -  DVDT_micro) / DT_MOIST
      if (associated( DTDT_micro))  DTDT_micro = ( TEMP       -  DTDT_micro) / DT_MOIST
      call MAPL_TimerOff (MAPL,"---CLDMICRO", __RC__)

      ! Exports

        call MAPL_GetPointer(EXPORT, PTR3D, 'DQRL', RC=STATUS); VERIFY_(STATUS)
        if(associated(PTR3D)) PTR3D = DQRDT_macro + DQRDT_micro

        ! Compute DBZ radar reflectivity
        call MAPL_GetPointer(EXPORT, PTR3D, 'DBZ'    , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, PTR2D, 'DBZ_MAX', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D) .OR. associated(PTR2D)) then
           call CALCDBZ(TMP3D,100*PLO,TEMP,Q1,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,0)
           if (associated(PTR3D)) PTR3D = TMP3D
           if (associated(PTR2D)) then
              PTR2D=-9999.0
              DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                 PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
              END DO ; END DO ; END DO
           endif
        endif

        ! Other Exports
        call MAPL_GetPointer(EXPORT, PTR2D, 'CWP', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS+QICN+QILS )*MASS , 3 )
        call MAPL_GetPointer(EXPORT, PTR2D, 'CLWP', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS ) *MASS , 3 )
        call MAPL_GetPointer(EXPORT, PTR2D, 'LWP', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS+QRAIN ) *MASS , 3 )
        call MAPL_GetPointer(EXPORT, PTR2D, 'IWP', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM( ( QICN+QILS+QSNOW+QGRAUPEL ) *MASS , 3 )
        call MAPL_GetPointer(EXPORT, PTR2D, 'TPW', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR2D)) PTR2D = SUM(   Q1*MASS , 3 )

   call MAPL_TimerOff(MAPL,"--MGB2_2M",__RC__)

end subroutine MGB2_2M_Run

end module GEOS_MGB2_2M_InterfaceMod
