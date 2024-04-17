! $Id$

#include "MAPL_Generic.h"
!#define PDFDIAG 1

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
  real    :: MINRHCRIT  
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
  logical :: USE_AV_V 


  real  :: DCS, WBFFACTOR, NC_CST, NI_CST, NG_CST, MUI_CST, &
           LCCIRRUS, UISCALE, LIU_MU, NPRE_FRAC, QCVAR_CST, &
           AUT_SCALE, TS_AUTO_ICE, CCN_PARAM, IN_PARAM, &
           FDROP_DUST, FDROP_SOOT, WSUB_OPTION, &
           DUST_INFAC, ORG_INFAC, BC_INFAC, SS_INFAC, RRTMG_IRRAD, RRTMG_SORAD,&
           MTIME,MINCDNC, Immersion_param, ACC_ENH, ACC_ENH_ICE, DT_MICRO, URSCALE, &
           CNV_GSC, CNV_BSC


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
    
    call ESMF_ConfigGetAttribute( CF, MGVERSION, Label="MGVERSION:",  default=3, __RC__)
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

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, CLLS, CLCN
    real, pointer, dimension(:,:,:)     :: NCPL, NCPI, NRAIN, NSNOW, NGRAUPEL
    
    logical  :: nccons, nicons, ngcons, do_graupel
    real(ESMF_KIND_R8)  Dcsr8, micro_mg_berg_eff_factor_in, ncnstr8, ninstr8, ngnstr8, mui_cnstr8

 

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: nn
    real                       :: tmprhL, tmprhO
  

    IAm = "MGB2_2M_Initialize"

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )
   
    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM.GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4
    
    
    
    
    call MAPL_GetResource( MAPL, LPHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default=.TRUE., RC=STATUS)
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
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NCPL,     'NCPL'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NCPI,     'NCPI'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NRAIN,    'NRAIN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NSNOW,    'NSNOW'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NGRAUPEL,  'NGRAUPEL' , RC=STATUS); VERIFY_(STATUS)

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
    call MAPL_GetResource( MAPL, MINRHCRIT, 'MINRHCRIT:', DEFAULT = 0.9, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, TURNRHCRIT, 'TURNRHCRIT:', DEFAULT = 884., RC=STATUS); VERIFY_(STATUS)

     
    !2M==tuning and options======
    
    
    call MAPL_GetResource(MAPL, UISCALE,        'UISCALE:',        DEFAULT= 1.0,    __RC__) !Scaling factor for sed vel of ice      
    call MAPL_GetResource(MAPL, LIU_MU,         'LIU_MU:',         DEFAULT= 2.0,    __RC__) !Liu autoconversion parameter
    call MAPL_GetResource(MAPL, NPRE_FRAC,      'NPRE_FRAC:',      DEFAULT= -1.0,   __RC__) !Fraction of preexisting ice affecting ice nucleationn            
    call MAPL_GetResource(MAPL, USE_AV_V,       'USE_AV_V:',       DEFAULT= .TRUE.,    __RC__) !Set to > 0 to use an average velocity for activation
    call MAPL_GetResource(MAPL, AUT_SCALE,      'AUT_SCALE:',      DEFAULT= 0.5,    __RC__) !scale factor for critical size for drizzle
    call MAPL_GetResource(MAPL, TS_AUTO_ICE,    'TS_AUTO_ICE:',    DEFAULT= 360.,    __RC__) !Ice autoconversion time scale
    call MAPL_GetResource(MAPL, CCN_PARAM,      'CCNPARAM:',       DEFAULT= 2.0,    __RC__) !CCN activation param
    call MAPL_GetResource(MAPL, IN_PARAM,       'INPARAM:',        DEFAULT= 6.0,    __RC__) !IN param
    call MAPL_GetResource(MAPL, Immersion_param,'ImmersionPARAM:', DEFAULT= 6.0,    __RC__) !Immersion param
    call MAPL_GetResource(MAPL, ACC_ENH,        'ACC_ENH:',        DEFAULT= 1.0,    __RC__) !accretion rain-liquid scaling for MG2
    call MAPL_GetResource(MAPL, ACC_ENH_ICE,    'ACC_ENH_ICE:',    DEFAULT= 1.0,    __RC__) !accretion snow-ice scaling for MG2
    call MAPL_GetResource(MAPL, FDROP_DUST,     'FDROP_DUST:',     DEFAULT= 0.5,    __RC__) !Fraction of dust within droplets for immersion freezing
    call MAPL_GetResource(MAPL, FDROP_SOOT,     'FDROP_SOOT:',     DEFAULT= 0.05,   __RC__) !Fraction of soot within droplets for immersion freezing        
    call MAPL_GetResource(MAPL, MINCDNC,        'MINCDNC:',        DEFAULT= 25.0,    __RC__) !min nucleated droplet conc. cm-3
    call MAPL_GetResource(MAPL, MTIME,          'MTIME:',          DEFAULT= -1.0,   __RC__) !Mixing time scale for aerosol within the cloud. Default is time step
    call MAPL_GetResource(MAPL, LCCIRRUS,       'LCCIRRUS:',       DEFAULT= 500.0,  __RC__) !Characteristic Length (m) of high freq gravity waves    
     call MAPL_GetResource(MAPL, QCVAR_CST,     'QCVAR_CST:',    DEFAULT= -1.,  __RC__) !Characteristic Length (m) of high freq gravity waves 
    !============

    call MAPL_GetResource(MAPL, DUST_INFAC,     'DUST_INFAC:',     DEFAULT= 1.0,    __RC__) !scalings for the INP concentrations
    call MAPL_GetResource(MAPL, BC_INFAC,       'BC_INFAC:',       DEFAULT= 0.1,    __RC__)
    call MAPL_GetResource(MAPL, ORG_INFAC,      'ORG_INFAC:',      DEFAULT= 1.0,    __RC__)
    call MAPL_GetResource(MAPL, SS_INFAC,       'SS_INFAC:',       DEFAULT= 1.0,    __RC__)    
    call MAPL_GetResource(MAPL, DT_MICRO,       'DT_MICRO:',       DEFAULT= 300.0,  __RC__)  ! time step of the microphysics substepping (s) (MG2) (5 min)
    call MAPL_GetResource(MAPL, URSCALE,       'URSCALE:',        DEFAULT= 1.0,    __RC__)  !Scaling factor for sed vel of rain    
    call MAPL_GetResource(MAPL, RRTMG_IRRAD ,  'USE_RRTMG_IRRAD:',DEFAULT=1.0,     __RC__)
    call MAPL_GetResource(MAPL, RRTMG_SORAD ,  'USE_RRTMG_SORAD:',DEFAULT=1.0,     __RC__)      
    call MAPL_GetResource(MAPL, CNV_GSC,   'CNV_GSC:', DEFAULT= 5.0e-5 ,RC=STATUS) !linear scaling for NCPL of conv detrainment 
    call MAPL_GetResource(MAPL, CNV_BSC,   'CNV_BSC:', DEFAULT= 0.3, RC=STATUS) !scaling for N=B*Nad for conv detrainment     
    call MAPL_GetResource(MAPL, DCS,      'DCS:'    , DEFAULT=200.0e-6, __RC__ ) !ice/snow separation diameter   
    Dcsr8 = DCS
    
    call MAPL_GetResource(MAPL, WBFFACTOR,   'WBFFACTOR:', DEFAULT= 0.1 ,__RC__) !scaling for the Bergeron-Findeinsen process rate    
    
    micro_mg_berg_eff_factor_in = WBFFACTOR
    call MAPL_GetResource(MAPL, NC_CST ,  'NC_CST:' , DEFAULT=  0.0 ,__RC__) !constant nd (set if greather than zero)      
    call MAPL_GetResource(MAPL, NI_CST ,  'NI_CST:' , DEFAULT=  0.0 ,__RC__) !constant nd (set if greather than zero)     
    call MAPL_GetResource(MAPL, NG_CST ,  'NG_CST:' , DEFAULT=  0.0 ,__RC__) !constant ng (set if greather than zero)     
    call MAPL_GetResource(MAPL, MUI_CST,  'MUI_CST:', DEFAULT= -1.0 ,__RC__) !constant ng (set if greather than zero) 
    
    call MAPL_GetResource(MAPL, WSUB_OPTION,  'WSUB_OPTION:',   DEFAULT= 1.0,    __RC__) !0- param 1- Use Wsub climatology 2-Wnet
    
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
                       nccons, nicons, ncnstr8, ninstr8, 2.0)
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
    real, pointer, dimension(:,:)   :: AREA, FRLAND, TS, DTSX, SH, EVAP, KPBL_SC
    real, pointer, dimension(:,:,:) :: SL2, SL3, QT2, QT3, W2, W3, SLQT, WQT, WQL, WSL
    real, pointer, dimension(:,:,:) :: WTHV2
    real, pointer, dimension(:,:,:) :: OMEGA
    
    real, pointer, dimension(:,:)   :: TAUOROX, TAUOROY
    real, pointer, dimension(:,:,:) :: ALH, RADLW, RADSW, WSUB_CLIM
    
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
    real, pointer, dimension(:,:,:) :: RHCRIT
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D
#ifdef PDFDIAG
    real, pointer, dimension(:,:,:) :: PDF_W1, PDF_W2, PDF_SIGW1, PDF_SIGW2,     &
                                       PDF_QT1, PDF_QT2, PDF_SIGQT1, PDF_SIGQT2, &
                                       PDF_TH1, PDF_TH2, PDF_SIGTH1, PDF_SIGTH2, &
                                       PDF_RQTTH, PDF_RWTH, PDF_RWQT
#endif    
    
    !2m
    real, pointer, dimension(:,:,:) :: SC_ICE, CDNC_NUC, INC_NUC, PFRZ, &
       CFICE, CFLIQ, DT_RASP, SMAX_LIQ, SMAX_ICE, WSUB, CCN01, CCN04, CCN1, &
       NHET_NUC, NLIM_NUC, SO4, ORG, BCARBON, DUST, SEASALT, NCPL_VOL, NCPI_VOL, &
       SAT_RAT, RHICE, RL_MASK, RI_MASK, &
       NHET_IMM, NHET_DEP, DUST_IMM, DUST_DEP, SIGW_GW, SIGW_CNV, SIGW_TURB, &
       SIGW_RC, BERG, BERGS, MELT, DNHET_CT, QCRES, QIRES, AUTICE, FRZPP_LS, &
       SNOWMELT_LS, DNCNUC, DNCSUBL, DNCHMSPLIT, DNCAUTICE, DNCACRIS, DNDCCN, &
       DNDACRLS, DNDACRLR, DNDEVAPC, DNDAUTLIQ, DNDCNV, DNICNV, &
       CNV_UPDF, CNV_CVW, DNHET_IMM, CNV_MFD, CNV_DQCDT, KAPPA, RHCmicro, RHLIQ, &
       CNV_NICE, CNV_NDROP, NWFA, CNV_FICE
       
     real, pointer, dimension(:,:)   :: EIS, LTS, QCVAR, &
       CCNCOLUMN, NDCOLUMN, NCCOLUMN
       
    
    
    real, allocatable, dimension(:,:,:) :: QCNTOT, CFX,  QTOT, &
       QL_TOT, QI_TOT, ACIL_LS_X, ACIL_AN_X, ACLL_LS_X, ACLL_AN_X, DLPDF_X, DIPDF_X, DLFIX_X, DIFIX_X, &
       AUT_X, SDM_X, FRZ_TT_X, FRZ_PP_X,  AIRDEN, TH1, FQA !check how much of these we are actually using
    
    integer, allocatable, dimension(:, :)   ::   KLCL
    real, allocatable, dimension(:, :)  :: CLDREFFI_TOP_X, CLDREFFL_TOP_X,  NCPL_TOP_X, NCPI_TOP_X, NCPL_CLDBASEX, uwind_gw

    ! Local variables
    real    :: ALPHA
    integer :: IM,JM,LM
    integer :: I, J, L, K
  

    integer :: num_steps_micro,  pcnst, n_modes, kbmin, kcldtop, kcldbot , &
                NAUX, kcldtopcvn, nbincontactdust, index, K0, KCBLMIN, i_src_mode, i_dst_mode
                  
    real, parameter :: pmin_trop = 10.0 !mbar minimum pressure to do cloud microphysics
    integer, parameter :: KMIN_TROP =  25

    REAL, allocatable, dimension(:,:) :: SCICE_tmp, FQA_tmp, cfaux 
                        
    real (ESMF_KIND_R8), dimension(3)       :: ccn_diag
    real(ESMF_KIND_R8), allocatable, dimension(:,:,:) :: rndstr8,naconr8  !Assume maximum 5 dust bins
    real(ESMF_KIND_R8), dimension(1)       :: prectr8, precir8            
    real (ESMF_KIND_R8)  :: disp_liu, ui_scale, dcrit, tfreez,  &
                            ts_autice, dcsr8,  scale_ri, mtimesc, ur_scale
      
    
    real(ESMF_KIND_R8), allocatable, dimension(:,:)  :: ttendr8, qtendr8, cwtendr8, &
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
      
      
      real(ESMF_KIND_R8)   :: autscx   
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
    call MAPL_GetPointer(IMPORT, KPBL_SC,  'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)
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
    call MAPL_GetPointer(EXPORT, QCVAR,       'QCVAR'   ,  ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SC_ICE,      'SC_ICE'      , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFR,    'RR'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFS,    'RS'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CLDREFFG,    'RG'          , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, CDNC_NUC,    'CDNC_NUC'    , ALLOC=.TRUE., __RC__) 
    call MAPL_GetPointer(EXPORT, INC_NUC,     'INC_NUC'     , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, PFRZ,        'PFRZ'        , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, LTS,         'LTS'         , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, EIS,         'EIS'         , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SMAX_LIQ,       'SMAX_LIQ'    , ALLOC=.TRUE., __RC__)
    call MAPL_GetPointer(EXPORT, SMAX_ICE,       'SMAX_ICE'    , ALLOC=.TRUE., __RC__)
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
    PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
    TH1       = T/PK
    AIRDEN = 100.*PLmb/T/MAPL_RGAS
    GZLO = MAPL_GRAV*ZL0
    
    ! Lowe tropospheric stability and estimated inversion strength
    call MAPL_GetPointer(EXPORT, LTS,   'LTS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EIS,   'EIS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    KLCL = FIND_KLCL( T, Q, PLmb, IM, JM, LM )    

    call FIND_EIS(TH1, QST3, T, ZL0, PLEmb, KLCL, IM, JM, LM, LTS, EIS)
    
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
    
        xscale = 8.7475*(real(imsize)**(-0.328)) ! scale for resolutions =! 50 km for WSUB_OPTION >= 1                   
        SIGW_RC  =  -OMEGA/AIRDEN/MAPL_GRAV + (RADLW + RADSW)*MAPL_CP/MAPL_GRAV
	    QL_TOT = QLCN+QLLS 
   		QI_TOT = QICN+QILS 
        QTOT   = QL_TOT+QI_TOT  

     !!=============== find vertical velocity variance


	   if (WSUB_OPTION .lt. 1.) then ! use parameterization from Barahona et al. GMD. 2014 (Appendix) 

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

         else  !WSUB climatology 

          WSUB  =  WSUB_CLIM    
          SIGW_TURB = WSUB                    
                 !call WRITE_PARALLEL ('Using Wclim***************') 

        end if      
                         
   ! ==========================================================================================    
   ! ========================Activate the aerosols ============================================ 


        do J=1,JM
          do I=1,IM
          
            kbmin= min(NINT(KPBL_SC(I, J)), LM-1)-2
            npre =   NPRE_FRAC   
            dpre=  1.0e-9      
            if (NPRE_FRAC < 0.0) npre =  CNV_FRC(I,J)*ABS(NPRE_FRAC) + (1-CNV_FRC(I,J))*0.05                       
                            
        	do K = KMIN_TROP, LM-1 !limit to troposphere and no activation at the surface
            
                 npre =  npre*NCPI(I,J,K)
                 if ((npre .gt. 0.0)  .and. (QI_TOT(I, J, K).gt. 0.)) dpre = ( QI_TOT(I, J, K)/(5400.0*npre*MAPL_PI))**(0.33) !Assume exponential distribution    

                 !!Subroutine aerosol_activate contains the CCN activation and ice nucleation parameterizations. Lives in aer_cloud.F90.

                 call   aerosol_activate(T(I, J, K), 100.*PLmb(I, J, K), WSUB(I, J, K), SIGW_RC(I, J, K), AeroProps(I, J, K), &
                      npre, dpre, ccn_diag, &
                      nact, SMAX_LIQ(I, J, K), INC_NUC (I, J, K), SMAX_ICE(I, J, K) ,  NHET_NUC(I, J, K), &
                      NHET_IMM(I, J, K), DNHET_IMM(I, J, K) , NHET_DEP(I, J, K) , SC_ICE(I, J, K) , &
                      DUST_IMM(I, J, K), DUST_DEP(I, J, K), NLIM_NUC(I, J, K), USE_AV_V, int(CCN_PARAM), int(IN_PARAM),  &
                      SO4(I, J, K), SEASALT(I, J, K), DUST(I, J, K), ORG(I, J, K), BCARBON(I, J, K), &                          
                                  FDROP_DUST, FDROP_SOOT, DUST_INFAC, BC_INFAC, ORG_INFAC, SS_INFAC, int(Immersion_PARAM))


                  CCN01(I, J, K) = max(ccn_diag(1), 0.0)
                  CCN04(I, J, K) = max(ccn_diag(2), 0.0)
                  CCN1 (I, J, K) = max(ccn_diag(3), 0.0)

                  if (K .ge. kbmin-4) nact = max(nact, (1.0-CNV_FRC(I, J))*MINCDNC*1.e6)

                  CDNC_NUC(I, J, K) =  nact


               end do 
            enddo
         enddo
      
        WSUB  =  SIGW_RC + 0.8*WSUB !diagnostic	

        where (T .gt. 238.0)
        SC_ICE  =  1.0
        end where 
        SC_ICE  =  MIN(MAX(SC_ICE, 1.0), 1.8)


             call MAPL_TimerOff(MAPL,"---ACTIV", __RC__)
        
 !=============================================End cloud particle nucleation=====================================
 !===============================================================================================================

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
		
         
        
      !====== Add convective detrainment of number concentration 
        
      call MAPL_GetPointer(EXPORT, CNV_NICE,  'CNV_NICE',  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      call MAPL_GetPointer(EXPORT, CNV_NDROP, 'CNV_NDROP', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      ! CNV_MFD includes Deep+Shallow mass flux

      call MAPL_GetPointer(EXPORT, CNV_MFD, 'CNV_MFD', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
      
      DO I=  1, IM
          DO J =  1, JM 
           kbmin =  max(min(NINT(KPBL_SC(I,J)), LM-1), NINT(0.8*LM))
           aux2= ZL0(I, J, kbmin ) !assume cldbase as PBLheight
           aux3  =  CDNC_NUC(I, J, kbmin)
           Do K  =  1, LM       
            call   make_cnv_ice_drop_number(CNV_NDROP(I, J, K), CNV_NICE(I, J, K), NHET_IMM(I, J, K),  \
                             aux3, ZL0(I, J, K), aux2, T(I, J, K), CNV_FICE(I, J, K), CNV_GSC, CNV_BSC)
           
           end do
          end do
         end do                            

      DNDCNV =  CNV_NDROP*CNV_MFD*iMASS
      DNICNV =  CNV_NICE*CNV_MFD*iMASS
            
      !update Number concentrations   
      NCPL =  NCPL + DNDCNV*DT_MOIST
      NCPI = NCPI  + DNICNV*DT_MOIST
         
    !==========================================================================================================
    !===================================Cloud Macrophysics ====================================================
    !==========================================================================================================

   
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
    
#ifdef PDFDIAG
   call MAPL_GetPointer(EXPORT,  PDF_W1,  'PDF_W1' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_W2,  'PDF_W2' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_SIGW1,  'PDF_SIGW1' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_SIGW2,  'PDF_SIGW2' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_QT1,  'PDF_QT1' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_QT2,  'PDF_QT2' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_SIGQT1,  'PDF_SIGQT1' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_SIGQT2,  'PDF_SIGQT2' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_TH1,  'PDF_TH1' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_TH2,  'PDF_TH2' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_SIGTH1,  'PDF_SIGTH1' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_SIGTH2,  'PDF_SIGTH2' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_RQTTH,  'PDF_RQTTH' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_RWTH,  'PDF_RWTH' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
   call MAPL_GetPointer(EXPORT,  PDF_RWQT,  'PDF_RWQT' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
#endif

      ! Include shallow precip condensates if present
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          QRAIN = QRAIN + PTR3D*DT_MOIST
        endif
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then 
          QSNOW = QSNOW + PTR3D*DT_MOIST
        endif
      
      
       !=========== evap/subl/pdf
         
        call MAPL_TimerOn(MAPL,"----hystpdf")
         
        call MAPL_GetPointer(EXPORT, RHCRIT,  'RHCRIT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        do I=1,IM
          do J=1,JM
           do L=1,LM
			
            DLPDF_X(I, J, L)=  QLLS(I, J, L) +QLCN(I, J, L)
            DIPDF_X(I, J, L)=  QILS(I, J, L) +QICN(I, J, L)
            
            call pdf_alpha(PLmb(I, J, L),PLmb(I, J, LM), ALPHA, FRLAND(I, J),  &
                              MINRHCRIT, TURNRHCRIT, EIS(I, J), 0) !0 uses old slingo formulation 
     		
            !include area scaling and limit RHcrit to > 70%
            ALPHA = min( 0.30, ALPHA*SQRT(SQRT(max(AREA(I,J), 0.0)/1.e10)) )
            RHCRIT(I, J, L) =  1.0 - ALPHA
            
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
                      NCPL(I,J,L)   , &
                      NCPI(I,J,L)   , &
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
#ifdef PDFDIAG
                      PDF_SIGW1(I,J,L),  &
                      PDF_SIGW2(I,J,L),  &
                      PDF_W1(I,J,L),     &
                      PDF_W2(I,J,L),     &
                      PDF_SIGTH1(I,J,L), &
                      PDF_SIGTH2(I,J,L), &
                      PDF_TH1(I,J,L),    &
                      PDF_TH2(I,J,L),    &
                      PDF_SIGQT1(I,J,L), &
                      PDF_SIGQT2(I,J,L), &
                      PDF_QT1(I,J,L),    &
                      PDF_QT2(I,J,L),    &
                      PDF_RQTTH(I,J,L),  &
                      PDF_RWTH(I,J,L),   &
                      PDF_RWQT(I,J,L),   &
#endif
                      WTHV2(I,J,L)   , &
                      WQL(I,J,L)     , &
                      .false.        , & 
                      .true., &
                      SC_ICE(I, J, L))
                      
         	DLPDF_X(I, J, L)=((QLLS(I, J, L)+QLCN(I, J, L)) - DLPDF_X(I, J, L))/DT_MOIST
         	DIPDF_X(I, J, L)=((QILS(I, J, L)+QICN(I, J, L)) - DIPDF_X(I, J, L))/DT_MOIST
         
           end do ! IM loop
         end do ! JM loop
       end do ! LM loop
       
     
      call MAPL_GetPointer(EXPORT, PTR3D,     'DIPDF'     , ALLOC=.TRUE., __RC__)
      PTR3D= DIPDF_X
      call MAPL_GetPointer(EXPORT, PTR3D,     'DLPDF'     , ALLOC=.TRUE., __RC__)
      PTR3D= DLPDF_X
    
      call MAPL_TimerOff(MAPL,"----hystpdf")
       
         do I=1,IM
          do J=1,JM
           do L=1,LM

         
       ! evaporation for CN/LS
             EVAPC(I,J,L) = Q(I,J,L)
             call EVAP3 (         &
                  DT_MOIST      , &
                  CCW_EVAP_EFF  , &
                  RHCRIT(I, J, L) , &
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

             SUBLC(I,J,L) =   Q(I,J,L)
             call SUBL3 (        &
                  DT_MOIST      , &
                  CCI_EVAP_EFF  , &
                  RHCRIT(I, J, L) , &
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
    PFL_LS =  0.0
    PFL_AN =  0.0
    PFI_LS =  0.0
    PFI_AN =  0.0
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
    naconr8   = 0.
    scale_ri =  1.3 ! scaling factor to account for the different definition of Ri in Chao and Suarez

    if ((RRTMG_SORAD .gt. 0.0) .or. (RRTMG_IRRAD .gt. 0.0)) then 
    scale_ri =  1.0
    end if 

    ! Update TH
    TH1 = T/PK
    
    !initialize MG variables
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
    autscx = 1.0
    disp_liu = LIU_MU
    ui_scale = UISCALE
    ur_scale  = URSCALE
    ts_autice = DT_MOIST
    mtimesc  = DT_MOIST     
    if (TS_AUTO_ICE .gt. 0.) ts_autice= TS_AUTO_ICE    
    if (MTIME .gt. 0.0) mtimesc=MTIME    
    xscale = (9000.0/real(imsize))**(-0.666)
	IF (QCVAR_CST .gt. 0.) then 
    	QCVAR =  QCVAR_CST
    else    
	    call estimate_qcvar(QCVAR, IM, JM, LM, PLmb, T, GZLO, Q, QST3, xscale)
    end if    
      
   
   
    do I=1,IM
	    do J=1,JM       
               kbmin =1                          
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
                frzcntr8 = frzimmr8 *0.0  
                frzdepr8(1,1:LM) = NHET_DEP(I, J, 1:LM)/DT_MOIST
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
        
         call  micro_mg_tend_interface ( DT_MICRO, INT(PDFSHAPE), 1.-RHCRIT(I, J, 1:LM), SCICE_tmp, FQA_tmp, &
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
                             ts_autice, ui_scale, autscx , disp_liu, nbincontactdust, ur_scale)


    end if 

        IF (MGVERSION > 1) then 

                   QRAIN(I,J,1:LM)  = max(QRAIN(I,J,1:LM) + REAL(qrtendr8(1, 1:LM)*DT_R8), 0.0) ! grid average 
                   QSNOW(I,J,1:LM)  = max(QSNOW(I,J,1:LM) + REAL(qstendr8(1, 1:LM)*DT_R8), 0.0) ! grid average                     
                   NRAIN(I,J,1:LM)  = max(NRAIN(I,J,1:LM) + REAL(nrtendr8(1, 1:LM)*DT_R8), 0.0)
                   NSNOW(I,J,1:LM)  = max(NSNOW(I,J,1:LM) + REAL(nstendr8(1, 1:LM)*DT_R8), 0.0)
                   CLDREFFR(I,J,1:LM)  = REAL(reff_rainr8(1,1:LM))
                   CLDREFFS(I,J,1:LM)  = REAL(reff_snowr8(1,1:LM))/scale_ri
                   CLDREFFG(I,J,1:LM)  = REAL(reff_graur8(1,1:LM))/scale_ri
                   QGRAUPEL(I,J,1:LM)  = max(QGRAUPEL(I,J,1:LM) + REAL(qgtendr8(1, 1:LM)*DT_R8), 0.0) ! grid average 
                   NGRAUPEL(I,J,1:LM)  = max(NGRAUPEL(I,J,1:LM) + REAL(ngtendr8(1, 1:LM)*DT_R8), 0.0)
                   
            
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
      
         
               PFL_LS(I, J, 1:LM) = rflxr8(1,2:LM+1) !+ lflxr8(1,1:LM)
               PFI_LS(I, J, 1:LM) = sflxr8(1,2:LM+1) + gflxr8(1,2:LM+1) !+  iflxr8(1,1:LM)
              
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

               
            IF (MGVERSION < 2) then   
               !normalize precip flux
               if (PFL_LS(I, J, LM) .gt. 1.0e-7) PFL_LS(I, J, 1:LM) = PFL_LS(I, J, 1:LM)*LS_PRCP(I,J)/PFL_LS(I, J, LM)                    
               if (PFI_LS(I, J, LM) .gt. 1.0e-7) PFI_LS(I, J, 1:LM) = PFI_LS(I, J, 1:LM)*LS_SNR(I,J)/PFI_LS(I, J, LM)                    
			end if 
               
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

          do I=1,IM
			  do J=1,JM
    			do K= 1, LM
        
                  
                  call update_cld( &
                         DT_MOIST                , &
                         1.- RHCRIT(I, J, K)        , &
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
            CLDREFFI =  MAPL_UNDEF
         end where

         where (QL_TOT .le. 0.0)
            CFLIQ =0.0
            NCPL  =0.0
            CLDREFFL =  MAPL_UNDEF
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
         TMP3D =  QLCN/(QLCN + QLLS+1.E-14)         
         PFL_AN(:,:,1:LM) = PFL_LS(:,:,1:LM) * TMP3D
         PFL_LS = PFL_LS - PFL_AN
         TMP3D =   QICN/(QICN + QILS + 1.E-14)
         PFI_AN(:,:,1:LM) = PFI_LS(:,:,1:LM) * TMP3D
         PFI_LS = PFI_LS - PFI_AN
      ! cleanup suspended precipitation condensates
         call FIX_NEGATIVE_PRECIP(RAD_QR, RAD_QS, RAD_QG)
         
         !=================================================================================
         !    Fill up diagnostics

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

         TH1 = T / PK

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
      if (associated(PTR2D)) PTR2D =  CLDREFFL_TOP_X
      
      call MAPL_GetPointer(EXPORT, PTR2D, 'NCPL_CLDBASE', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D=  NCPL_CLDBASEX
      
      call MAPL_GetPointer(EXPORT, PTR2D, 'NCPL_TOP', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D=  NCPL_TOP_X
      
      call MAPL_GetPointer(EXPORT, PTR2D, 'NCPI_TOP', RC=STATUS); VERIFY_(STATUS)
      if (associated(PTR2D)) PTR2D=  NCPI_TOP_X
      
         


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
           call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,1)
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
