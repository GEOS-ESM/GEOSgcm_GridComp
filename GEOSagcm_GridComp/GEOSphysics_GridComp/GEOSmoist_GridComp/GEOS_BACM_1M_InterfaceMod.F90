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
  use GEOSmoist_Process_Library
  use CLOUDNEW, only: CLDPARAMS, PROGNO_CLOUD

  implicit none

  private

  character(len=ESMF_MAXSTR)        :: COMP_NAME
  character(len=ESMF_MAXSTR)        :: IAm
  integer                           :: STATUS

  ! specify how to handle friendlies with DYN:TRB:CHM:ANA
  type FRIENDLIES_TYPE
         character(len=ESMF_MAXSTR) :: QV
         character(len=ESMF_MAXSTR) :: CLLS
         character(len=ESMF_MAXSTR) :: CLCN
         character(len=ESMF_MAXSTR) :: QLLS
         character(len=ESMF_MAXSTR) :: QLCN
         character(len=ESMF_MAXSTR) :: QILS
         character(len=ESMF_MAXSTR) :: QICN
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  ! Local resource variables
  integer :: imsize
  integer :: CFPBL_EXP
  real    :: MINRHCRITLND
  real    :: MINRHCRITOCN
  real    :: MAXRHCRITLND
  real    :: MAXRHCRITOCN
  real    :: TURNRHCRIT

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

      FRIENDLIES%QV       = "DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS"
      FRIENDLIES%CLLS     = "DYNAMICS"
      FRIENDLIES%CLCN     = "DYNAMICS"
      FRIENDLIES%QLLS     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QLCN     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QILS     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%QICN     = "DYNAMICS:TURBULENCE"

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

    call MAPL_TimerAdd(GC, name="--BACM_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine BACM_1M_Setup

subroutine BACM_1M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QW

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: nn
    integer                    :: LM
    real                       :: tmprhL, tmprhO, TMP_ICEFALL

    call MAPL_Get ( MAPL, LM=LM, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)

    call MAPL_GetPointer(INTERNAL, Q,        'Q'       , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLLS,     'QLLS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QLCN,     'QLCN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QW,       'QW'      , RC=STATUS); VERIFY_(STATUS)
    QW = Q+QLLS+QLCN+QILS+QICN

    call MAPL_GetResource( MAPL, CLDPARAMS%CCW_EVAP_EFF,   'CCW_EVAP_EFF:',   DEFAULT= 4.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%CCI_EVAP_EFF,   'CCI_EVAP_EFF:',   DEFAULT= 4.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%PDFSHAPE,       'PDFSHAPE:',       DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_BETA,       'CNV_BETA:',       DEFAULT= 10.0    )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_BETA,       'ANV_BETA:',       DEFAULT= 4.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_BETA,        'LS_BETA:',        DEFAULT= 4.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%RH_CRIT,        'RH_CRIT:',        DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%QC_CRIT_LS,     'QC_CRIT_LS:',     DEFAULT= 8.0e-4  )
    call MAPL_GetResource( MAPL, CLDPARAMS%ACCRETION,      'ACCRETION:',      DEFAULT= 2.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%RAIN_REVAP_FAC, 'RAIN_REVAP_FAC:', DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SNOW_REVAP_FAC, 'SNOW_REVAP_FAC:', DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%AUTOC_LS,       'AUTOC_LS:',       DEFAULT= 1.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%AUTOC_ANV,      'AUTOC_ANV:',      DEFAULT= 1.0e-3  )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_SUND_INTER,  'LS_SUND_INTER:',  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_SUND_COLD,   'LS_SUND_COLD:',   DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_SUND_TEMP1,  'LS_SUND_TEMP1:',  DEFAULT= 230.    )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_SUND_INTER, 'ANV_SUND_INTER:', DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_SUND_COLD,  'ANV_SUND_COLD:',  DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_SUND_TEMP1, 'ANV_SUND_TEMP1:', DEFAULT= 230.    )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_TO_LS_TIME, 'ANV_TO_LS_TIME:', DEFAULT= 14400.  )
    call MAPL_GetResource( MAPL, CLDPARAMS%DISABLE_RAD,    'DISABLE_RAD:',    DEFAULT= 0.      )
    call MAPL_GetResource( MAPL, CLDPARAMS%REVAP_OFF_P,    'REVAP_OFF_P:',    DEFAULT= 2000.   )
    call MAPL_GetResource( MAPL, CLDPARAMS%ICE_RAMP,       'ICE_RAMP:',       DEFAULT= -27.0   )
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
    call MAPL_GetResource( MAPL, CLDPARAMS%MIN_RI,         'MIN_RI:',         DEFAULT=  15.e-6 )
    call MAPL_GetResource( MAPL, CLDPARAMS%MAX_RI,         'MAX_RI:',         DEFAULT= 150.e-6 )
    call MAPL_GetResource( MAPL, CLDPARAMS%FAC_RL,         'FAC_RL:',         DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%MIN_RL,         'MIN_RL:',         DEFAULT=  5.e-6  )
    call MAPL_GetResource( MAPL, CLDPARAMS%MAX_RL,         'MAX_RL:',         DEFAULT= 21.e-6  )
    call MAPL_GetResource( MAPL, CLDPARAMS%CNV_ENVF,       'CNV_ENVF:',       DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%ANV_ENVF,       'ANV_ENVF:',       DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%SC_ENVF,        'SC_ENVF:',        DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%LS_ENVF,        'LS_ENVF:',        DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_LS_WAT,      'FR_LS_WAT:',      DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_AN_WAT,      'FR_AN_WAT:',      DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_LS_ICE,      'FR_LS_ICE:',      DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CLDPARAMS%FR_AN_ICE,      'FR_AN_ICE:',      DEFAULT= 0.0     )
    call MAPL_GetResource( MAPL, CFPBL_EXP,                'CFPBL_EXP:',      DEFAULT= 1       )

    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4

    tmprhL = CEILING(100.0*(1.0-min(0.20, max(0.01, 0.1*SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))/100.0 ! roundup by 0.01s
    tmprhL = min(0.99,tmprhL)
    tmprhO = min(0.99,tmprhL)
    call MAPL_GetResource( MAPL, MINRHCRITLND,             'MINRHCRITLND:',   DEFAULT=tmprhL   )
    call MAPL_GetResource( MAPL, MINRHCRITOCN,             'MINRHCRITOCN:',   DEFAULT=tmprhO   )
    call MAPL_GetResource( MAPL, MAXRHCRITLND,             'MAXRHCRITOCN:',   DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, MAXRHCRITOCN,             'MAXRHCRITLND:',   DEFAULT= 1.0     )
    call MAPL_GetResource( MAPL, TURNRHCRIT,               'TURNRHCRIT:',     DEFAULT= 750.0  )

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
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, KH
    real, pointer, dimension(:,:)   :: FRLAND, TS, DTSX, TROPP, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:) :: HL2, HL3, QT2, QT3, W2, W3, HLQT, WQT, WQL, WHL, EDMF_FRC
    real, pointer, dimension(:,:,:) :: OMEGA
    real, pointer, dimension(:,:,:) :: CNV_MFD, CNV_DQCDT, CNV_PRC3, CNV_UPDF
    real, pointer, dimension(:,:,:) :: MFD_SC, QLDET_SC, QIDET_SC, SHLW_PRC3, SHLW_SNO3, CUFRC_SC
    ! Local
    real, allocatable, dimension(:,:,:) :: PLEmb, PKE, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  PK,  ZL0
    real, allocatable, dimension(:,:,:) :: DZET,  TH,  MASS
    real, allocatable, dimension(:,:,:) :: QDDF3, DQST3, QST3
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: turnrhcrit2D
    real, allocatable, dimension(:,:)   :: minrhcrit2D
    real, allocatable, dimension(:,:)   :: maxrhcrit2D
    real, allocatable, dimension(:,:)   :: TMP2D
    type(ESMF_State)                    :: AERO
    ! Exports
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro
    real, pointer, dimension(:,:,:) ::  DUDT_micro,  DVDT_micro,  DTDT_micro, DTHDT_micro
    real, pointer, dimension(:,:,:) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, pointer, dimension(:,:,:) :: CLDREFFL, CLDREFFI 
    real, pointer, dimension(:,:,:) :: RHX
    real, pointer, dimension(:,:,:) :: REV_CN, REV_SC, REV_AN, REV_LS
    real, pointer, dimension(:,:,:) :: RSU_CN, RSU_SC, RSU_AN, RSU_LS
    real, pointer, dimension(:,:,:) :: ACLL_CN, ACLL_SC, ACLL_AN, ACLL_LS
    real, pointer, dimension(:,:,:) :: ACIL_CN, ACIL_SC, ACIL_AN, ACIL_LS
    real, pointer, dimension(:,:,:) :: PFL_CN, PFL_SC, PFL_AN, PFL_LS
    real, pointer, dimension(:,:,:) :: PFI_CN, PFI_SC, PFI_AN, PFI_LS
    real, pointer, dimension(:,:,:) :: DLPDF, DIPDF, DLFIX, DIFIX
    real, pointer, dimension(:,:,:) :: AUT, EVAPC, SUBLC, SDM
    real, pointer, dimension(:,:,:) :: VFALLICE_AN, VFALLICE_LS
    real, pointer, dimension(:,:,:) :: VFALLWAT_AN, VFALLWAT_LS
    real, pointer, dimension(:,:,:) :: VFALLRN_AN, VFALLRN_LS, VFALLRN_CN, VFALLRN_SC
    real, pointer, dimension(:,:,:) :: VFALLSN_AN, VFALLSN_LS, VFALLSN_CN, VFALLSN_SC
    real, pointer, dimension(:,:,:) :: FRZ_TT, FRZ_PP
    real, pointer, dimension(:,:,:) :: DCNVL, DCNVI
    real, pointer, dimension(:,:,:) :: ALPHT, ALPH1, ALPH2
    real, pointer, dimension(:,:,:) :: CFPDF, RHCLR
    real, pointer, dimension(:,:,:) :: DQRL, WTHV2
    real, pointer, dimension(:,:,:) :: PDF_A
    real, pointer, dimension(:,:  ) :: CNV_FRC
    real, pointer, dimension(:,:  ) :: LS_PRCP, CN_PRCP, AN_PRCP, SC_PRCP
    real, pointer, dimension(:,:  ) :: LS_ARF,  CN_ARF,  AN_ARF,  SC_ARF
    real, pointer, dimension(:,:  ) :: LS_SNR,  CN_SNR,  AN_SNR,  SC_SNR
    real, pointer, dimension(:,:,:) :: NCPL_VOL, NCPI_VOL
    real, pointer, dimension(:,:,:) :: CFICE, CFLIQ
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D

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
    call MAPL_GetPointer(INTERNAL, NACTL,    'NACTL'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, NACTI,    'NACTI'   , RC=STATUS); VERIFY_(STATUS)

    ! Import State
    call MAPL_GetPointer(IMPORT, ZLE,     'ZLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, PLE,     'PLE'     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, T,       'T'       , RC=STATUS); VERIFY_(STATUS)
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
    call MAPL_GetPointer(IMPORT, TS,      'TS'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, TROPP,   'TROPP'   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, KPBLSC,  'KPBL_SC' , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, SH,      'SH'      , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, EVAP,    'EVAP'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(IMPORT, OMEGA,   'OMEGA'   , RC=STATUS); VERIFY_(STATUS)
    call   ESMF_StateGet(IMPORT, 'AERO', AERO, __RC__)

    ! Allocatables
     ! Edge variables 
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PLEmb(IM,JM,0:LM) )
    ALLOCATE ( PKE  (IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( ZL0  (IM,JM,LM  ) )
    ALLOCATE ( PLmb (IM,JM,LM  ) )
    ALLOCATE ( PK   (IM,JM,LM  ) )
    ALLOCATE ( DZET (IM,JM,LM  ) )
    ALLOCATE ( TH   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( QDDF3(IM,JM,LM  ) )
    ALLOCATE ( DQST3(IM,JM,LM  ) )
    ALLOCATE (  QST3(IM,JM,LM  ) )
    ALLOCATE ( TMP3D(IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( turnrhcrit2D (IM,JM) )
    ALLOCATE ( minrhcrit2D  (IM,JM) )
    ALLOCATE ( maxrhcrit2D  (IM,JM) )
    ALLOCATE ( TMP2D        (IM,JM) )

    ! Derived States
    PLEmb    =  PLE*.01
    PKE      = (PLE/MAPL_P00)**(MAPL_KAPPA)
    PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
    PK       = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
    DO L=0,LM
       ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
    END DO
    ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
    DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
    TH       = T/PK
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

    do J=1,JM
       do I=1,IM
          if (TURNRHCRIT .LT. 0) then
             turnrhcrit2D(I,J) = PLmb(I, J, NINT(KPBLSC(I,J)))-50.  ! 50mb above KHSFC top
          else
             turnrhcrit2D(I,J) = TURNRHCRIT
          endif
          minrhcrit2D(I,J) = MINRHCRITOCN*(1.0-FRLAND(I,J)) + MINRHCRITLND*FRLAND(I,J)     
          maxrhcrit2D(I,J) = MAXRHCRITOCN*(1.0-FRLAND(I,J)) + MAXRHCRITLND*FRLAND(I,J)
       end do
    end do

    ! Export Tendencies
    call MAPL_GetPointer(EXPORT, DQVDT_micro, 'DQVDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQIDT_micro, 'DQIDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQLDT_micro, 'DQLDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQADT_micro, 'DQADT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DUDT_micro,  'DUDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DVDT_micro,  'DVDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT,  DTDT_micro,  'DTDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DTHDT_micro, 'DTHDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    if (associated(DQVDT_micro)) DQVDT_micro = Q
    if (associated(DQLDT_micro)) DQLDT_micro = QLLS + QLCN
    if (associated(DQIDT_micro)) DQIDT_micro = QILS + QICN
    if (associated(DQADT_micro)) DQADT_micro = CLLS + CLCN
    if (associated( DUDT_micro))  DUDT_micro  = U
    if (associated( DVDT_micro))  DVDT_micro  = V
    if (associated( DTDT_micro))  DTDT_micro  = T
    if (associated(DTHDT_micro)) DTHDT_micro  = TH

    ! Imports which are Exports from local siblings
    ! Required export MUST have been filled in GridComp
    call MAPL_GetPointer(EXPORT, CNV_FRC, 'CNV_FRC', RC=STATUS); VERIFY_(STATUS)
    ! DeepCu : default to 0.0 in MAPL if not running DeepCu
    call MAPL_GetPointer(EXPORT, CNV_MFD,    'CNV_MFD '  ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_DQCDT,  'CNV_DQCDT' ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_PRC3,   'CNV_PRC3'  ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CNV_UPDF,   'CNV_UPDF'  ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! ShallowCu : default to 0.0 in MAPL if not running ShallowCu
    call MAPL_GetPointer(EXPORT, MFD_SC,     'MFD_SC'    ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QLDET_SC,   'QLDET_SC'  ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, QIDET_SC,   'QIDET_SC'  ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SHLW_PRC3,  'SHLW_PRC3' ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SHLW_SNO3,  'SHLW_SNO3' ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CUFRC_SC,   'CUFRC_SC'  ,  ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
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
    ! Exports passed to be filled in progno_cloud
    call MAPL_GetPointer(EXPORT, LS_PRCP,  'LS_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CN_PRCP,  'CN_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AN_PRCP,  'AN_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_PRCP,  'SC_PRCP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, LS_ARF,   'LS_ARF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CN_ARF,   'CN_ARF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AN_ARF,   'AN_ARF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_ARF,   'SC_ARF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, LS_SNR,   'LS_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CN_SNR,   'CN_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AN_SNR,   'AN_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SC_SNR,   'SC_SNR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RHX   ,   'RHX'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_CN,   'REV_CN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_SC,   'REV_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_AN,   'REV_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, REV_LS,   'REV_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_CN,   'RSU_CN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_SC,   'RSU_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_AN,   'RSU_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RSU_LS,   'RSU_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACLL_CN,  'ACRLL_CN', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACLL_SC,  'ACRLL_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACLL_AN,  'ACRLL_AN', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACLL_LS,  'ACRLL_LS', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACIL_CN,  'ACRIL_CN', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACIL_SC,  'ACRIL_SC', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACIL_AN,  'ACRIL_AN', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ACIL_LS,  'ACRIL_LS', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_CN,   'PFL_CN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_SC,   'PFL_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_AN,   'PFL_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFL_LS,   'PFL_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_CN,   'PFI_CN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_SC,   'PFI_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_AN,   'PFI_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PFI_LS,   'PFI_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DLPDF,    'DLPDF'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DIPDF,    'DIPDF'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DLFIX,    'DLFIX'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DIFIX,    'DIFIX'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, AUT,      'AUT'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, EVAPC,    'EVAPC'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SDM,      'SDM'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLICE_AN, 'VFALLICE_AN' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLICE_LS, 'VFALLICE_LS' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLWAT_AN, 'VFALLWAT_AN' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLWAT_LS, 'VFALLWAT_LS' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLRN_AN,  'VFALLRN_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLRN_LS,  'VFALLRN_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLRN_CN,  'VFALLRN_CN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLRN_SC,  'VFALLRN_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLSN_AN,  'VFALLSN_AN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLSN_LS,  'VFALLSN_LS'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLSN_CN,  'VFALLSN_CN'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, VFALLSN_SC,  'VFALLSN_SC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SUBLC,     'SUBLC'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZ_TT,    'FRZ_TT' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, FRZ_PP,    'FRZ_PP' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DCNVL,     'DCNVL'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DCNVI,     'DCNVI'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ALPHT,     'ALPHT'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ALPH1,     'ALPH1'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, ALPH2,     'ALPH2'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, CFPDF,     'CFPDF'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, RHCLR,     'RHCLR'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, DQRL,      'DQRL'   , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, PDF_A,     'PDF_A'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WTHV2,     'WTHV2'  , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, WQL,       'WQL'    , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)

         call  PROGNO_CLOUD (     &
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
              QT3               , &   ! <- turb
              HL3               , &   ! <- turb
              CNV_MFD           , &   ! <- deep
              CNV_DQCDT         , &   ! <- deep              
              CNV_PRC3          , &   ! <- deep   
              CNV_UPDF          , &   ! <- deep
              MFD_SC            , &   ! <- shlw   
              QLDET_SC          , &   ! <- shlw   
              QIDET_SC          , &   ! <- shlw   
              SHLW_PRC3         , &   ! <- shlw   
              SHLW_SNO3         , &   ! <- shlw   
              CUFRC_SC*0.5      , &   ! <- shlw   
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
              RAD_QR            , &
              RAD_QS            , &
              RAD_QG            , &
              CLDREFFL          , &
              CLDREFFI          , &
              LS_PRCP           , &
              CN_PRCP           , &
              AN_PRCP           , &
              SC_PRCP           , &
              LS_ARF            , &
              CN_ARF            , &
              AN_ARF            , &
              SC_ARF            , &
              LS_SNR            , &
              CN_SNR            , &
              AN_SNR            , &
              SC_SNR            , &
              minrhcrit2D       , &
              maxrhcrit2D       , &
              turnrhcrit2D      , &
              QST3              , &
              DZET              , &
              QDDF3             , &
              CNV_FRC           , &
              TROPP             , &
                                ! Diagnostics
              RHX               , &
              REV_LS            , &
              REV_AN            , &
              REV_CN            , &
              REV_SC            , &
              RSU_LS            , &
              RSU_AN            , &
              RSU_CN            , &
              RSU_SC            , &
              ACLL_CN  ,ACIL_CN , &
              ACLL_AN  ,ACIL_AN , &
              ACLL_LS  ,ACIL_LS , &
              ACLL_SC  ,ACIL_SC , &
               PFL_CN  , PFI_CN , &
               PFL_AN  , PFI_AN , &
               PFL_LS  , PFI_LS , &
               PFL_SC  , PFI_SC , &
              DLPDF    ,DIPDF   , &
              DLFIX    ,DIFIX   , &
              AUT, EVAPC, SDM, SUBLC,  &
              FRZ_TT, DCNVL, DCNVI, &
              ALPHT, ALPH1, ALPH2,  &
              CFPDF, RHCLR, DQRL, FRZ_PP,  &
              VFALLICE_AN  ,VFALLICE_LS  , &
              VFALLWAT_AN  ,VFALLWAT_LS  , &
               VFALLSN_AN  , VFALLSN_LS  ,VFALLSN_CN  ,VFALLSN_SC  ,  &
               VFALLRN_AN  , VFALLRN_LS  ,VFALLRN_CN  ,VFALLRN_SC  ,  &
              PDF_A, &
              WTHV2, WQL, &
              .TRUE., &
              NACTL,    &
              NACTI,    &
              "GF" )

         call MAPL_GetPointer(EXPORT, DTSX, 'DTSX', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
         DTSX = TS - TH(:,:,LM)*PK(:,:,LM) + (MAPL_GRAV/MAPL_CP)*ZL0(:,:,LM)
         if (CFPBL_EXP > 1) then
           TMP3D = RAD_CF
           do L = 1,LM
              where( (KH(:,:,L-1) > 5.) .AND. (DTSX> 0.) )
                 TMP3D(:,:,L) = RAD_CF(:,:,L)**CFPBL_EXP
              endwhere
           enddo
           where( TMP3D > 0.0 )
              RAD_QL = RAD_QL*(RAD_CF/TMP3D)
              RAD_QI = RAD_QI*(RAD_CF/TMP3D)
              RAD_QR = RAD_QR*(RAD_CF/TMP3D)
              RAD_QS = RAD_QS*(RAD_CF/TMP3D)
              RAD_QG = RAD_QG*(RAD_CF/TMP3D)
           endwhere
           RAD_CF = TMP3D
         endif
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
         ! Cloud fraction exports
         call MAPL_GetPointer(EXPORT, CFICE, 'CFICE', RC=STATUS); VERIFY_(STATUS)
         if (associated(CFICE)) then
           CFICE=0.0
           WHERE (RAD_QI .gt. 1.0e-12)
              CFICE=RAD_CF*RAD_QI/(RAD_QL+RAD_QI)
           END WHERE
           CFICE=MAX(MIN(CFICE, 1.0), 0.0)
         endif
         call MAPL_GetPointer(EXPORT, CFLIQ, 'CFLIQ', RC=STATUS); VERIFY_(STATUS)
         if (associated(CFLIQ)) then
           CFLIQ=0.0
           WHERE (RAD_QL .gt. 1.0e-12)
              CFLIQ=RAD_CF*RAD_QL/(RAD_QL+RAD_QI)
           END WHERE
           CFLIQ=MAX(MIN(CFLIQ, 1.0), 0.0)
         endif
         ! Rain-out and Relative Humidity where RH > 110%
         DQST3 = GEOS_DQSAT(TH*PK, PLmb, QSAT=QST3)
         where ( Q > 1.1*QST3 )
            TMP3D = (Q - 1.1*QST3)/( 1.0 + 1.1*DQST3*MAPL_ALHL/MAPL_CP )
         elsewhere
            TMP3D = 0.0
         endwhere
         LS_PRCP = LS_PRCP + SUM(TMP3D*MASS,3)/DT_MOIST
         Q  =  Q - TMP3D
         TH = TH + (MAPL_ALHL/MAPL_CP)*TMP3D/PK
         T  = TH*PK

         ! cleanup any negative QV/QC/CF
         call FILLQ2ZERO(Q       , MASS, TMP2D)
         call MAPL_GetPointer(EXPORT, PTR2D, 'FILLNQV', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR2D)) PTR2D = TMP2D
         call FILLQ2ZERO(QLLS    , MASS, TMP2D)
         call FILLQ2ZERO(QLCN    , MASS, TMP2D)
         call FILLQ2ZERO(QILS    , MASS, TMP2D)
         call FILLQ2ZERO(QICN    , MASS, TMP2D)
         call FILLQ2ZERO(CLLS    , MASS, TMP2D)
         call FILLQ2ZERO(CLCN    , MASS, TMP2D)
         where (QILS+QICN .le. 0.0)
            CLDREFFI = 36.0e-6
         end where
         where (QLLS+QLCN .le. 0.0)
            CLDREFFL = 14.0e-6
         end where

         if (associated(DQVDT_micro)) DQVDT_micro = ( Q          - DQVDT_micro) / DT_MOIST
         if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
         if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
         if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
         if (associated( DUDT_micro))  DUDT_micro = ( U          -  DUDT_micro) / DT_MOIST
         if (associated( DVDT_micro))  DVDT_micro = ( V          -  DVDT_micro) / DT_MOIST
         if (associated( DTDT_micro))  DTDT_micro = (TH*PK       -  DTDT_micro) / DT_MOIST
         if (associated(DTHDT_micro)) DTHDT_micro = (TH          - DTHDT_micro) / DT_MOIST

         ! Compute DBZ radar reflectivity
         call MAPL_GetPointer(EXPORT, PTR3D, 'DBZ'    , RC=STATUS); VERIFY_(STATUS)
         call MAPL_GetPointer(EXPORT, PTR2D, 'DBZ_MAX', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D) .OR. associated(PTR2D)) then
            call CALCDBZ(TMP3D,100*PLmb,TH*PK,Q,RAD_QR,RAD_QS,RAD_QG,IM,JM,LM,1,0,0)
            if (associated(PTR3D)) PTR3D = TMP3D
            if (associated(PTR2D)) then
               PTR2D=-9999.0
               DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                  PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
               END DO ; END DO ; END DO
            endif
         endif

    call MAPL_TimerOff (MAPL,"--BACM_1M")

end subroutine BACM_1M_Run

end module GEOS_BACM_1M_InterfaceMod
