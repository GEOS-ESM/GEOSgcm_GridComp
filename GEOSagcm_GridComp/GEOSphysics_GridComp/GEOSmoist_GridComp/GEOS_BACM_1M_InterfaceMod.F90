! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_BACM_1M_InterfaceMod -- A Module to interface with the
!   BACM_1M cloud microphysics

module GEOS_BACM_1M_InterfaceMod

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use Aer_Actv_Single_Moment,only: Aer_Actv_1M_interface, USE_AEROSOL_NN

  implicit none

      type CLDPARAM_TYPE
           real               :: CNV_BETA              ! 1
           real               :: ANV_BETA              ! 2
           real               :: LS_BETA               ! 3
           real               :: RH_CRIT               ! 4
           real               :: AUTOC_LS              ! 5
           real               :: QC_CRIT_LS            ! 6
           real               :: ACCRETION             ! 7
           real               :: RAIN_REVAP_FAC        ! 8
           real               :: VOL_TO_FRAC           ! 9
           real               :: SUPERSAT              ! 10
           real               :: SHEAR_EVAP_FAC        ! 11
           real               :: MIN_ALLOW_CCW         ! 12
           real               :: CCW_EVAP_EFF          ! 13
           real               :: CCI_EVAP_EFF          ! 13
           real               :: NSUB_AUTOCONV         ! 14
           real               :: LS_SUND_INTER         ! 15
           real               :: LS_SUND_COLD          ! 16
           real               :: LS_SUND_TEMP1         ! 17
           real               :: ANV_SUND_INTER        ! 18
           real               :: ANV_SUND_COLD         ! 19
           real               :: ANV_SUND_TEMP1        ! 20
           real               :: ANV_TO_LS_TIME        ! 21
           real               :: CCN_OCEAN             ! 22
           real               :: CCN_LAND              ! 23
           real               :: NCCN_ANVIL_NULL       ! 24
           real               :: NCCN_PBL_NULL         ! 25
           real               :: DISABLE_RAD           ! 26
           real               :: ICE_SETTLE            ! 27
           real               :: ANV_ICEFALL           ! 28
           real               :: LS_ICEFALL            ! 29
           real               :: REVAP_OFF_P           ! 30
           real               :: CNV_ENVF              ! 31
           real               :: ANV_ENVF              ! 31
           real               :: SC_ENVF               ! 31
           real               :: LS_ENVF               ! 31
           real               :: WRHODEP               ! 32
           real               :: ICE_RAMP              ! 33
           real               :: CNV_ICEPARAM          ! 34
           real               :: CNV_ICEFRPWR          ! 35
           real               :: CNV_DDRF              ! 36
           real               :: ANV_DDRF              ! 37
           real               :: LS_DDRF               ! 38
           real               :: AUTOC_ANV             ! 39
           real               :: QC_CRIT_ANV           ! 40
           real               :: TANHRHCRIT            ! 41
           real               :: MINRHCRIT             ! 42
           real               :: MAXRHCRIT             ! 43
           real               :: PRECIPRAD             ! 44
           real               :: TURNRHCRIT            ! 45
           real               :: MAXRHCRITLAND         ! 46
           real               :: FR_LS_WAT             ! 47
           real               :: FR_LS_ICE             ! 48
           real               :: FR_AN_WAT             ! 49
           real               :: FR_AN_ICE             ! 50
           real               :: MIN_RL                ! 51
           real               :: MIN_RI                ! 52
           real               :: MAX_RL                ! 53
           real               :: MAX_RI                ! 54
           real               :: FAC_RL                ! 55
           real               :: FAC_RI                ! 56
           real               :: SNOW_REVAP_FAC        ! 57
           real               :: PDFSHAPE              ! 58
           real               :: TURNRHCRIT_UP         ! 59
           real               :: SLOPERHCRIT           ! 60
           real               :: MIN_LTS               ! 61
           integer            :: CFPBL_EXP             ! 62
           real               :: DISP_FACTOR_LIQ       ! 63
           real               :: DISP_FACTOR_ICE       ! 63
           real               :: SCLM_SHALLOW          ! 63
           real               :: SCLM_DEEP             ! 63
      endtype CLDPARAM_TYPE
      type (CLDPARAM_TYPE) :: CLDPARAMS

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
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  character(len=ESMF_MAXSTR)        :: COMP_NAME

  public :: BACM_1M_Setup, BACM_1M_Initialize, BACM_1M_Run

contains

subroutine BACM_1M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code

    IAm = "GEOS_BACM_1M_InterfaceMod"
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)       
    Iam = trim(COMP_NAME) // Iam

    ! !INTERNAL STATE:

      FRIENDLIES%QV       = trim(COMP_NAME)//":DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS"
      FRIENDLIES%QW       = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%CLLS     = trim(COMP_NAME)//":DYNAMICS"
      FRIENDLIES%CLCN     = trim(COMP_NAME)//":DYNAMICS"
      FRIENDLIES%QLLS     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QLCN     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QILS     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"
      FRIENDLIES%QICN     = trim(COMP_NAME)//":DYNAMICS:TURBULENCE"

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
         SHORT_NAME = 'QW',                                        &
         LONG_NAME  = 'mass_fraction_of_wet_air',                  &
         UNITS      = 'kg kg-1',                                   &
         RESTART    = MAPL_RestartSkip,                            &
         FRIENDLYTO = trim(FRIENDLIES%QW),                         &
         DIMS       = MAPL_DimsHorzVert,                           &
         VLOCATION  = MAPL_VLocationCenter,             RC=STATUS  )  
    VERIFY_(STATUS)                                                      

    call MAPL_TimerAdd(GC, name="--BACM_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine BACM_1M_Setup

subroutine BACM_1M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State),         pointer  :: GIM(:)
    type (ESMF_State),         pointer  :: GEX(:)
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QW

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: imsize,nn
    integer                    :: LM
    real                       :: tmprh, TMP_ICEFALL

    call MAPL_Get ( MAPL, LM=LM, GIM=GIM, GEX=GEX, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QW,       'QW'      , RC=STATUS); VERIFY_(STATUS)
    QW = Q+QLLS+QLCN+QILS+QICN

    call MAPL_GetResource(MAPL,USE_AEROSOL_NN,'USE_AEROSOL_NN:',default=.TRUE., RC=STATUS )
    VERIFY_(STATUS)
    call aer_cloud_init()
    call WRITE_PARALLEL ("INITIALIZED aer_cloud_init for BACM_1M")

    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4

    call MAPL_GetResource( MAPL, CLDPARAMS%MIN_LTS,        'LTS_LOW:',        DEFAULT= 20.0    ) ! lower LTS for morphology correction
    call MAPL_GetResource( MAPL, CLDPARAMS%DISP_FACTOR_LIQ,'DISP_FACTOR_LIQ:',DEFAULT= 1.0     ) ! Scales the droplet/ice crystal number in convective detrainment 
    call MAPL_GetResource( MAPL, CLDPARAMS%DISP_FACTOR_ICE,'DISP_FACTOR_ICE:',DEFAULT= 1.0     ) ! Scales the droplet/ice crystal number in convective detrainment 
    tmprh = CEILING(100.0*(1.0-min(0.20, max(0.01, 0.1*SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))/100.0
    call MAPL_GetResource( MAPL, CLDPARAMS%MINRHCRIT,      'MINRHCRIT:'    ,  DEFAULT=tmprh    )
    call MAPL_GetResource( MAPL, CLDPARAMS%MAXRHCRIT,      'MAXRHCRIT:'    ,  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%MAXRHCRITLAND,  'MAXRHCRITLAND:',  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%TANHRHCRIT,     'TANHRHCRIT:'   ,  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%TURNRHCRIT,     'TURNRHCRIT:'   ,  DEFAULT= -999.0  )
    call MAPL_GetResource( MAPL, CLDPARAMS%TURNRHCRIT_UP,  'TURNRHCRIT_UP:',  DEFAULT= 300.0   )
    call MAPL_GetResource( MAPL, CLDPARAMS%SLOPERHCRIT,    'SLOPERHCRIT:',    DEFAULT= 20.0    )
    call MAPL_GetResource( MAPL, CLDPARAMS%CCW_EVAP_EFF,   'CCW_EVAP_EFF:',   DEFAULT= 4.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%CCI_EVAP_EFF,   'CCI_EVAP_EFF:',   DEFAULT= 1.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_ICEPARAM,   'CNV_ICEPARAM:',   DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SCLM_DEEP,      'SCLM_DEEP:',      DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SCLM_SHALLOW,   'SCLM_SHALLOW:',   DEFAULT= 2.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%PDFSHAPE,       'PDFSHAPE:',       DEFAULT= 2.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_BETA,       'CNV_BETA:',       DEFAULT= 10.0    )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_BETA,       'ANV_BETA:',       DEFAULT= 4.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_BETA,        'LS_BETA:',        DEFAULT= 4.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%RH_CRIT,        'RH_CRIT:',        DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%QC_CRIT_LS,     'QC_CRIT_LS:',     DEFAULT= 8.0e-4  )
    call MAPL_GetResource( MAPL, CLDPARAMS%ACCRETION,      'ACCRETION:',      DEFAULT= 2.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%RAIN_REVAP_FAC, 'RAIN_REVAP_FAC:', DEFAULT= 1.00    )
    call MAPL_GetResource( MAPL, CLDPARAMS%VOL_TO_FRAC,    'VOL_TO_FRAC:',    DEFAULT= -1.0    )
    call MAPL_GetResource( MAPL, CLDPARAMS%SUPERSAT,       'SUPERSAT:',       DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SHEAR_EVAP_FAC, 'SHEAR_EVAP_FAC:', DEFAULT= 1.3     )
    call MAPL_GetResource( MAPL, CLDPARAMS%MIN_ALLOW_CCW,  'MIN_ALLOW_CCW:',  DEFAULT= 1.0e-9  )
    call MAPL_GetResource( MAPL, CLDPARAMS%AUTOC_LS,       'AUTOC_LS:',       DEFAULT= 1.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%AUTOC_ANV,      'AUTOC_ANV:',      DEFAULT= 1.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%NSUB_AUTOCONV,  'NSUB_AUTOCONV:',  DEFAULT= 20.     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_SUND_INTER,  'LS_SUND_INTER:',  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_SUND_COLD,   'LS_SUND_COLD:',   DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_SUND_TEMP1,  'LS_SUND_TEMP1:',  DEFAULT= 230.    )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_SUND_INTER, 'ANV_SUND_INTER:', DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_SUND_COLD,  'ANV_SUND_COLD:',  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_SUND_TEMP1, 'ANV_SUND_TEMP1:', DEFAULT= 230.    )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_TO_LS_TIME, 'ANV_TO_LS_TIME:', DEFAULT= 14400.  )
    call MAPL_GetResource( MAPL, CLDPARAMS%CCN_OCEAN,      'NCCN_OCEAN:',     DEFAULT= 300.    )
    call MAPL_GetResource( MAPL, CLDPARAMS%CCN_LAND,       'NCCN_LAND:',      DEFAULT= 100.    )
    call MAPL_GetResource( MAPL, CLDPARAMS%DISABLE_RAD,    'DISABLE_RAD:',    DEFAULT= 0.      )
    call MAPL_GetResource( MAPL, CLDPARAMS%REVAP_OFF_P,    'REVAP_OFF_P:',    DEFAULT= 2000.   )
    call MAPL_GetResource( MAPL, CLDPARAMS%ICE_RAMP,       'ICE_RAMP:',       DEFAULT= -27.0   )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_ICEFRPWR,   'CNV_ICEFRPWR:',   DEFAULT= 4.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_DDRF,       'CNV_DDRF:',       DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_DDRF,       'ANV_DDRF:',       DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_DDRF,        'LS_DDRF:',        DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%QC_CRIT_ANV,    'QC_CRIT_ANV:',    DEFAULT= 8.0e-4  )
    call MAPL_GetResource( MAPL, CLDPARAMS%ICE_SETTLE,     'ICE_SETTLE:',     DEFAULT= 1.      )
    SELECT CASE ( LM )
       CASE ( 72 )
           TMP_ICEFALL = 0.5
       CASE ( 91 )
           TMP_ICEFALL = 0.25
       CASE ( 181 )
           TMP_ICEFALL = 0.125
       CASE DEFAULT
           TMP_ICEFALL = 1.0
    END SELECT
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_ICEFALL,    'ANV_ICEFALL:',    DEFAULT= TMP_ICEFALL )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_ICEFALL,     'LS_ICEFALL:',     DEFAULT= TMP_ICEFALL )
    call MAPL_GetResource( MAPL, CLDPARAMS%FAC_RI,         'FAC_RI:',         DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%MIN_RI,         'MIN_RI:',         DEFAULT=   5.e-6 )
    call MAPL_GetResource( MAPL, CLDPARAMS%MAX_RI,         'MAX_RI:',         DEFAULT= 140.e-6 )
    call MAPL_GetResource( MAPL, CLDPARAMS%FAC_RL,         'FAC_RL:',         DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%MIN_RL,         'MIN_RL:',         DEFAULT=  2.5e-6 )
    call MAPL_GetResource( MAPL, CLDPARAMS%MAX_RL,         'MAX_RL:',         DEFAULT= 60.0e-6 )
    call MAPL_GetResource( MAPL, CLDPARAMS%PRECIPRAD,      'PRECIPRAD:',      DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SNOW_REVAP_FAC, 'SNOW_REVAP_FAC:', DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_ENVF,       'CNV_ENVF:',       DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_ENVF,       'ANV_ENVF:',       DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SC_ENVF,        'SC_ENVF:',        DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_ENVF,        'LS_ENVF:',        DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_LS_WAT,      'FR_LS_WAT:',      DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_AN_WAT,      'FR_AN_WAT:',      DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_LS_ICE,      'FR_LS_ICE:',      DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_AN_ICE,      'FR_AN_ICE:',      DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%CFPBL_EXP,      'CFPBL_EXP:',      DEFAULT= 1       )

end subroutine BACM_1M_Initialize


subroutine BACM_1M_Run (GC, IMPORT, EXPORT, CLOCK, RC)
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
    integer                         :: I, J, L
    real, pointer, dimension(:,:)   :: LONS
    real, pointer, dimension(:,:)   :: LATS

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QW

    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, TH, U, V, KH
    real, pointer, dimension(:,:)   :: FRLAND, TROPP, CAPE, KPBLSC
    real, pointer, dimension(:,:,:) :: HL2,       &
                                       HL3,       &
                                       QT2,       &
                                       QT3,       &
                                       W2,        &
                                       W3,        &
                                       HLQT,      &
                                       WQT,       &
                                       WQL,       &
                                       WHL,       &
                                       EDMF_FRC
    ! Local
    real, allocatable, dimension(:,:,:) :: PLEmb, PKE, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  PK,  ZL0
    real, allocatable, dimension(:,:,:) :: DZET,  T,   MASS
    real, allocatable, dimension(:,:,:) :: QDDF3, DQST3, QST3
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: turnrhcrit2D
    real, allocatable, dimension(:,:)   :: minrhcrit2D
    real, allocatable, dimension(:,:)   :: maxrhcrit2D
    real, allocatable, dimension(:,:)   :: CNV_FRACTION
    real, allocatable, dimension(:,:)   :: TMP2D
    real                                :: CNV_FRACTION_MIN
    real                                :: CNV_FRACTION_MAX
    real                                :: CNV_FRACTION_EXP
    character(len=ESMF_MAXSTR)          :: GRIDNAME
    character(len=4)                    :: imchar
    character(len=2)                    :: dateline
    integer                             :: imsize,nn
    real                                :: tmprhO, tmprhL
    ! Exports
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro
    real, pointer, dimension(:,:,:) :: DUDT_micro, DVDT_micro, DTDT_micro
 
    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ) 
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--BACM_1M")

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

    ! Internal State
    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLCN,     'CLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, CLLS,     'CLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)

    ! Import State
    call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TH,      'TH'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, U,       'U'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, V,       'V'       , RC=STATUS); VERIFY_(STATUS)
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
    call MAPL_GetPointer(IMPORT, TROPP,   'TROPP'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, CAPE,    'CAPE'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBLSC,  'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,LM+1) )
    ALLOCATE ( PLEmb(IM,JM,LM+1) )
    ALLOCATE ( PKE  (IM,JM,LM+1) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PLmb (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( DZET (IM,JM,LM  ) )
    ALLOCATE ( T    (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( QDDF3(IM,JM,LM  ) )
    ALLOCATE ( DQST3(IM,JM,LM  ) )
    ALLOCATE (  QST3(IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( CNV_FRACTION (IM,JM) )
    ALLOCATE ( turnrhcrit2D (IM,JM) )
    ALLOCATE ( minrhcrit2D  (IM,JM) )
    ALLOCATE ( maxrhcrit2D  (IM,JM) )
    ALLOCATE ( TMP2D        (IM,JM) )

    ! Derived States
    PLEmb    =  PLE*.01
    PKE      = (PLE/MAPL_P00)**(MAPL_KAPPA)
    PLmb     = 0.5*(PLEmb(:,:,1:LM) + PLEmb(:,:,2:LM+1))
    PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
    DO L=1,LM+1
       ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
    END DO
    ZL0      = 0.5*(ZLE0(:,:,1:LM) + ZLE0(:,:,2:LM+1) ) ! Layer Height (m) above the surface
    DZET     =     (ZLE0(:,:,1:LM) - ZLE0(:,:,2:LM+1) ) ! Layer thickness (m)
    T        = TH*PK
    DQST3    = GEOS_DQSAT(T, PLmb, QSAT=QST3)
    MASS     = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )/MAPL_GRAV

    WHERE ( ZL0 < 3000. )
       QDDF3 = -( ZL0-3000. ) * ZL0 * MASS
    ELSEWHERE
       QDDF3 = 0.
    END WHERE
    TMP2D = SUM(QDDF3, 3)
    DO L = 1,LM
       QDDF3(:,:,L) = QDDF3(:,:,L) / TMP2D
    END DO

    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4
    tmprhL = min(0.99,(1.0-min(0.20, max(0.01, 0.035*SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))
    tmprhO = min(0.99,(1.0-min(0.20, max(0.01, 0.140*SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))
    do J=1,JM
       do I=1,IM
          if (CLDPARAMS%TURNRHCRIT .LT. 0) then
             turnrhcrit2D(I,J) = PLmb(I, J, NINT(KPBLSC(I,J)))-50.  ! 50mb above KHSFC top
             minrhcrit2D(I,J) = tmprhO             *(1.0-FRLAND(I,J)) + tmprhL                 *FRLAND(I,J)
             maxrhcrit2D(I,J) = CLDPARAMS%MAXRHCRIT*(1.0-FRLAND(I,J)) + CLDPARAMS%MAXRHCRITLAND*FRLAND(I,J)
          else
             turnrhcrit2D(I,J) = MIN( CLDPARAMS%TURNRHCRIT , CLDPARAMS%TURNRHCRIT-(1020-PLEmb(i,j,LM)) )
             minrhcrit2D(I,J) = CLDPARAMS%MINRHCRIT
             maxrhcrit2D(I,J) = CLDPARAMS%MAXRHCRIT*(1.0-FRLAND(I,J)) + CLDPARAMS%MAXRHCRITLAND*FRLAND(I,J)
          endif
       end do
    end do

    call MAPL_GetResource(MAPL,CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=  500.0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource(MAPL,CNV_FRACTION_EXP, 'CNV_FRACTION_EXP:', DEFAULT=    1.0, RC=STATUS)
    VERIFY_(STATUS)
    if( CNV_FRACTION_MAX > CNV_FRACTION_MIN ) then
       ! CAPE
       WHERE (CAPE .ne. MAPL_UNDEF)
          CNV_FRACTION =(MAX(1.e-6,MIN(1.0,(CAPE-CNV_FRACTION_MIN)/(CNV_FRACTION_MAX-CNV_FRACTION_MIN))))
       END WHERE
    endif
    if (CNV_FRACTION_EXP > 1.0) then
        CNV_FRACTION = CNV_FRACTION**CNV_FRACTION_EXP
    elseif (CNV_FRACTION_EXP > 0.0) then
        CNV_FRACTION = 1.0-(1.0-CNV_FRACTION)**(1.0/CNV_FRACTION_EXP)
    else
        CNV_FRACTION = 1
    endif

    ! Export Tendencies
    call MAPL_GetPointer(EXPORT, DQVDT_micro, 'DQVDT_micro'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDT_micro, 'DQIDT_micro'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDT_micro, 'DQLDT_micro'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQADT_micro, 'DQADT_micro'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DUDT_micro, 'DUDT_micro'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DVDT_micro, 'DVDT_micro'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTDT_micro, 'DTDT_micro'      , RC=STATUS); VERIFY_(STATUS)
    if (associated(DQVDT_micro)) DQVDT_micro = Q
    if (associated(DQLDT_micro)) DQLDT_micro = QLLS + QLCN
    if (associated(DQIDT_micro)) DQIDT_micro = QILS + QICN
    if (associated(DQADT_micro)) DQADT_micro = CLLS + CLCN
    if (associated(DUDT_micro) ) DUDT_micro  = U
    if (associated(DVDT_micro) ) DVDT_micro  = V
    if (associated(DTDT_micro) ) DTDT_micro  = T

    ! Export and/or sratch Variable
    call MAPL_GetPointer(EXPORT, RAD_CF,   'FCLD', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QV,   'QV'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QL,   'QL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RAD_QI,   'QI'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CLDREFFL, 'RL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CLDREFFI, 'RI'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

         SC_SNR = 0.0
         SC_PRC2 = 0.0
         SC_ARFX = 0.0
         REV_SC_X = 0.0
         RSU_SC_X = 0.0
         PFL_SC_X = 0.0
         PFI_SC_X = 0.0

         call  PROGNO_CLOUD (                    &
              IM*JM, LM         , &
              DT_MOIST          , &
              LATS              , &
              PLmb              , &
              ZL0               , &
              PLEmb             , &
              PK                , &
              FRLAND            , &   ! <- surf
              KH                , &   ! <- turb
              EDMF_FRC          , &   ! <- turb
              WQT               , &   ! <- turb
              WHL               , &   ! <- turb
              QT2               , &   ! <- turb
              HL2               , &   ! <- turb
              HLQT              , &   ! <- turb
              W2                , &   ! <- turb
              W3                , &   ! <- turb
              QT3               , &
              HL3               , &
              DTS               , &
              CNV_MFD           , &   ! <- ras
              CNV_DQLDT         , &   ! <- ras              
              CNV_PRC3          , &   ! <- ras   
              CNV_UPDF          , &   ! <- ras
              MFD_SC            , &   ! <- shlw   
              QLDET_SC          , &   ! <- shlw   
              QIDET_SC          , &   ! <- shlw   
              SHLW_PRC3         , &   ! <- shlw   
              SHLW_SNO3         , &   ! <- shlw   
              UFRC_SC           , &   ! <- shlw   
              U                 , &
              V                 , &
              TH                , &
              Q                 , &
              QLLS              , &
              QLCN              , &
              QILS              , &
              QICN              , &
              CLCN              , &
              CLLS              , &
              RAD_CF            , &
              RAD_QV            , &
              RAD_QL            , &
              RAD_QI            , &
              QRN               , &
              QSN               , &
              RAD_QG            , &
              QPLS              , &
              CLDREFFL          , &
              CLDREFFI          , &
              LS_PRC2           , &
              CN_PRC2           , &
              AN_PRC2           , &
              SC_PRC2           , &
              LS_ARFX           , &
              CN_ARFX           , &
              AN_ARFX           , &
              SC_ARFX           , &
              LS_SNR            , &
              CN_SNR            , &
              AN_SNR            , &
              SC_SNR            , &
              CLDPARAMS         , &
              minrhcrit2D       , &
              maxrhcrit2D       , &
              turnrhcrit2D      , &
              QST3              , &
              DZET              , &
              QDDF3             , &
              CNV_FRACTION      , &
              TROPP             , &
                                ! Diagnostics
              RHX_X             , &
              REV_LS_X          , &
              REV_AN_X          , &
              REV_CN_X          , &
              REV_SC_X          , &
              RSU_LS_X          , &
              RSU_AN_X          , &
              RSU_CN_X          , &
              RSU_SC_X          , &
              ACLL_CN_X,ACIL_CN_X   , &
              ACLL_AN_X,ACIL_AN_X   , &
              ACLL_LS_X,ACIL_LS_X   , &
              ACLL_SC_X,ACIL_SC_X   , &
              PFL_CN_X,PFI_CN_X     , &
              PFL_AN_X,PFI_AN_X     , &
              PFL_LS_X,PFI_LS_X     , &
              PFL_SC_X,PFI_SC_X     , &
              DLPDF_X,DIPDF_X,DLFIX_X,DIFIX_X,    &
              AUT_X, EVAPC_X , SDM_X , SUBLC_X ,  &
              FRZ_TT_X, DCNVL_X, DCNVi_X,       &
              ALPHT_X, ALPH1_X, ALPH2_X, CFPDF_X, &
              RHCLR_X,                      &
              DQRL_X, FRZ_PP_X,               &
              VFALLICE_AN_X,VFALLICE_LS_X,    &
              VFALLWAT_AN_X,VFALLWAT_LS_X,    &
              VFALLSN_AN_X,VFALLSN_LS_X,VFALLSN_CN_X,VFALLSN_SC_X,  &
              VFALLRN_AN_X,VFALLRN_LS_X,VFALLRN_CN_X,VFALLRN_SC_X,  &
              PDF_A, &
#ifdef PDFDIAG
              PDF_SIGW1, PDF_SIGW2, PDF_W1, PDF_W2, &
              PDF_SIGTH1, PDF_SIGTH2, PDF_TH1, PDF_TH2, &
              PDF_SIGQT1, PDF_SIGQT2, PDF_QT1, PDF_QT2, &
              PDF_RQTTH, PDF_RWTH, PDF_RWQT,            &
#endif
              WTHV2, WQL, &
              adjustl(SHALLOW_OPTION)=="UW",   &
              NACTL,    &
              NACTI,    &
              CONVPAR_OPTION )

         VERIFY_(STATUS)

         if (associated(PDF_AX)) PDF_AX = PDF_A

         if (USE_AEROSOL_NN) then
           CFX =100.*PLmb*r_air/T !density times conversion factor
           NCPL = NACTL/CFX ! kg-1
           NCPI = NACTI/CFX ! kg-1
         else
           NCPL = 0.
           NCPI = 0.
         endif

         !-----------------------------------------
         !! kluge in some skewness for PBL clouds 
         !! where DTS:=TS-TSFCAIR > 0. i.e., unstable
         RAD_QV   = max( Q1 , 0. )

         if (CLDPARAMS%CFPBL_EXP > 1.0) then
          CFPBL = RAD_CF
          do L = 1,LM
             where( (KH(:,:,L-1) > 5.) .AND. (DTS > 0.) )
                CFPBL(:,:,L) = RAD_CF(:,:,L)**CLDPARAMS%CFPBL_EXP
             endwhere
          enddo
          where( CFPBL > 0.0 )
             RAD_QL = RAD_QL*(RAD_CF/CFPBL)
             RAD_QI = RAD_QI*(RAD_CF/CFPBL)
             QRN    = QRN   *(RAD_CF/CFPBL)
             QSN    = QSN   *(RAD_CF/CFPBL)
          endwhere
          RAD_CF = CFPBL
         endif

         RAD_QL = MIN( RAD_QL , 0.001 )  ! Still a ridiculously large
         RAD_QI = MIN( RAD_QI , 0.001 )  ! value.
         QRN    = MIN( QRN   , 0.01 )  ! value.
         QSN    = MIN( QSN   , 0.01 )  ! value.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! Set rain water for radiation to 0 if preciprad flag is off (set to 0)
         if(CLDPARAMS%PRECIPRAD.eq.0.) then
            RAD_QR = 0.
            RAD_QS = 0.
            RAD_QG = 0.
         else
            RAD_QR = QRN
            RAD_QS = QSN
         endif

         if (associated(QRTOT)) QRTOT = QRN*RAD_CF
         if (associated(QSTOT)) QSTOT = QSN*RAD_CF
         if (associated(QGTOT)) QGTOT = 0.0

         if (associated(QPTOTLS)) QPTOTLS = QPLS

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

         if (associated(DQVDT_micro)) DQVDT_micro = (Q1 - DQVDT_micro         ) / DT_MOIST
         if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
         if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
         if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
         if (associated(DUDT_micro) ) DUDT_micro  = (U1 - DUDT_micro - U1) / DT_MOIST
         if (associated(DVDT_micro) ) DVDT_micro  = (V1 - DVDT_micro - V1) / DT_MOIST
         if (associated(DTDT_micro) ) DTDT_micro  = (TH1*PK - DTDT_micro) / DT_MOIST

    call MAPL_TimerOff (MAPL,"--BACM_1M")

end subroutine BACM_1M_Run

end module GEOS_BACM_1M_InterfaceMod
