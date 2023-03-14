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
  use GEOSmoist_Process_Library
  use Aer_Actv_Single_Moment
  use gfdl2_cloud_microphys_mod

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
  integer :: imsize
  real    :: TURNRHCRIT
  real    :: MINRHCRITLND
  real    :: MINRHCRITOCN
  real    :: MAXRHCRITLND
  real    :: MAXRHCRITOCN
  real    :: CCW_EVAP_EFF
  real    :: CCI_EVAP_EFF
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

    call MAPL_TimerAdd(GC, name="--GFDL_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine GFDL_1M_Setup

subroutine GFDL_1M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, QW

    character(len=ESMF_MAXSTR) :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: nn
    real                       :: tmprhL, tmprhO

    type(ESMF_VM) :: VM
    integer :: comm

    call MAPL_GetResource( MAPL, LHYDROSTATIC, Label="HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
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
    call MAPL_GetPointer(INTERNAL, QW,       'QW'      , RC=STATUS); VERIFY_(STATUS)
    QW = Q+QLLS+QLCN+QILS+QICN+QRAIN+QSNOW+QGRAUPEL

    call ESMF_VMGetCurrent(VM, _RC)
    call ESMF_VMGet(VM, mpiCommunicator=comm, _RC)

    call gfdl_cloud_microphys_init(comm)
    call WRITE_PARALLEL ("INITIALIZED GFDL_1M microphysics in non-generic GC INIT")

    call MAPL_GetResource(MAPL, GRIDNAME, 'AGCM_GRIDNAME:', RC=STATUS)
    VERIFY_(STATUS)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4

    call MAPL_GetResource( MAPL, TURNRHCRIT      , 'TURNRHCRIT:'      , DEFAULT= -999.0 , RC=STATUS); VERIFY_(STATUS)
    tmprhL = min(0.99,(1.0-min(0.30, max(0.01, dw_land  * SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))
    tmprhO = min(0.99,(1.0-min(0.30, max(0.01, dw_ocean * SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))
    call MAPL_GetResource( MAPL, MINRHCRITLND    , 'MINRHCRITLND:'    , DEFAULT=tmprhL , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MINRHCRITOCN    , 'MINRHCRITOCN:'    , DEFAULT=tmprhO , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAXRHCRITLND    , 'MAXRHCRITOCN:'    , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAXRHCRITOCN    , 'MAXRHCRITLND:'    , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CCW_EVAP_EFF    , 'CCW_EVAP_EFF:'    , DEFAULT= 2.0e-3, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CCI_EVAP_EFF    , 'CCI_EVAP_EFF:'    , DEFAULT= 1.0e-3, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, PDFSHAPE        , 'PDFSHAPE:'        , DEFAULT= 2     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ANV_ICEFALL     , 'ANV_ICEFALL:'     , DEFAULT= 0.8   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LS_ICEFALL      , 'LS_ICEFALL:'      , DEFAULT= 0.8   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RI          , 'FAC_RI:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RI          , 'MIN_RI:'          , DEFAULT=  5.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RI          , 'MAX_RI:'          , DEFAULT=140.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RL          , 'FAC_RL:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RL          , 'MIN_RL:'          , DEFAULT= 2.5e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RL          , 'MAX_RL:'          , DEFAULT=60.0e-6, RC=STATUS); VERIFY_(STATUS)

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

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL, QW
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W, KH
    real, pointer, dimension(:,:)   :: AREA, FRLAND, TS, DTSX, TROPP, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:) :: HL2, HL3, QT2, QT3, W2, W3, HLQT, WQT, WQL, WHL, EDMF_FRC
    real, pointer, dimension(:,:,:) :: WTHV2
    real, pointer, dimension(:,:,:) :: OMEGA
    ! Local
    real, allocatable, dimension(:,:,:) :: PLEmb, PKE, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  PK,  ZL0
    real, allocatable, dimension(:,:,:) :: DZ, DZET, DP, MASS, iMASS
    real, allocatable, dimension(:,:,:) :: DQST3, QST3, TH
    real, allocatable, dimension(:,:,:) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                                           DQSDTmic, DQGDTmic, DQADTmic, &
                                            DUDTmic,  DVDTmic,  DTDTmic
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: turnrhcrit2D
    real, allocatable, dimension(:,:)   :: minrhcrit2D
    real, allocatable, dimension(:,:)   :: maxrhcrit2D
    real, allocatable, dimension(:,:)   :: TMP2D
    ! Exports
    real, pointer, dimension(:,:  ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ) :: LS_PRCP, LS_SNR, CNV_FRC
    real, pointer, dimension(:,:,:) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:) ::  DUDT_macro,  DVDT_macro,  DTDT_macro, DTHDT_macro
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro
    real, pointer, dimension(:,:,:) ::  DUDT_micro,  DVDT_micro,  DTDT_micro, DTHDT_micro
    real, pointer, dimension(:,:,:) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, pointer, dimension(:,:,:) :: CLDREFFL, CLDREFFI
    real, pointer, dimension(:,:,:) :: EVAPC, SUBLC
    real, pointer, dimension(:,:,:) :: RHX, REV_LS, RSU_LS
    real, pointer, dimension(:,:,:) :: PFL_LS, PFL_AN
    real, pointer, dimension(:,:,:) :: PFI_LS, PFI_AN
    real, pointer, dimension(:,:,:) :: PDF_A
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D

    ! Local variables
    real    :: turnrhcrit_up
    real    :: ALPHAl, ALPHAu, ALPHA, RHCRIT
    integer :: IM,JM,LM
    integer :: I, J, L

    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS )
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--GFDL_1M",RC=STATUS)

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
    ALLOCATE ( DZ   (IM,JM,LM  ) )
    ALLOCATE ( TH   (IM,JM,LM  ) )
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
    DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
    MASS     = DP/MAPL_GRAV
    iMASS    = 1.0/MASS

    do J=1,JM
       do I=1,IM
         if ( (TURNRHCRIT .LT. 0) .or. (FRLAND(I,J) .GT. 0.0) ) then
             turnrhcrit2D(I,J) = PLmb(I, J, NINT(KPBLSC(I,J)))-50.  ! 50mb above KHSFC top
          else
             turnrhcrit2D(I,J) = MIN( TURNRHCRIT , TURNRHCRIT-(1020-PLEmb(i,j,LM)) )
          endif
          minrhcrit2D(I,J) = MINRHCRITOCN*(1.0-FRLAND(I,J)) + MINRHCRITLND*FRLAND(I,J)
          maxrhcrit2D(I,J) = MAXRHCRITOCN*(1.0-FRLAND(I,J)) + MAXRHCRITLND*FRLAND(I,J)
       end do
    end do

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
    call MAPL_GetPointer(EXPORT, DTHDT_macro, 'DTHDT_macro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    if (associated( DUDT_macro))  DUDT_macro=U
    if (associated( DVDT_macro))  DVDT_macro=V
    if (associated( DTDT_macro))  DTDT_macro=T
    if (associated(DTHDT_macro)) DTHDT_macro=TH
    if (associated(DQVDT_macro)) DQVDT_macro=Q
    if (associated(DQLDT_macro)) DQLDT_macro=QLCN+QLLS
    if (associated(DQIDT_macro)) DQIDT_macro=QICN+QILS
    if (associated(DQADT_macro)) DQADT_macro=CLCN+CLLS
    if (associated(DQRDT_macro)) DQRDT_macro=QRAIN
    if (associated(DQSDT_macro)) DQSDT_macro=QSNOW
    if (associated(DQGDT_macro)) DQGDT_macro=QGRAUPEL
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
        do L=1,LM
          do J=1,JM
           do I=1,IM
             if (.not. do_qa) then ! if not doing the evap/subl/pdf inside of GFDL-MP
       ! Send the condensates through the pdf after convection
       !  Use Slingo-Ritter (1985) formulation for critical relative humidity
             ALPHA = maxrhcrit2D(I,J)
           ! lower turn from maxrhcrit
             if (PLmb(i,j,l) .le. turnrhcrit2D(I,J)) then
                ALPHAl = minrhcrit2D(I,J)
             else
                if (L.eq.LM) then
                   ALPHAl = maxrhcrit2D(I,J)
                else
                   ALPHAl = minrhcrit2D(I,J) + (maxrhcrit2D(I,J)-minrhcrit2D(I,J))/(19.) * &
                           ((atan( (2.*(PLmb(i,j,l)-turnrhcrit2D(I,J))/min(100., PLEmb(i,j,LM)-turnrhcrit2D(I,J))-1.) * &
                           tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                endif
             endif
           ! upper turn back to maxrhcrit
             turnrhcrit_up = TROPP(i,j)/100.0
             IF (turnrhcrit_up == MAPL_UNDEF) turnrhcrit_up = 100.
             if (PLmb(i,j,l) .le. turnrhcrit_up) then
                ALPHAu = maxrhcrit2D(I,J)
             else
                ALPHAu = maxrhcrit2D(I,J) - (maxrhcrit2D(I,J)-minrhcrit2D(I,J))/(19.) * &
                        ((atan( (2.*(PLmb(i,j,l)-turnrhcrit_up)/( turnrhcrit2D(I,J)-turnrhcrit_up)-1.) * &
                        tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
             endif
           ! combine and limit
             ALPHA = min( 0.30, 1.0 - min(max(ALPHAl,ALPHAu),1.) ) ! restrict RHcrit to > 70%
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
                  NACTL(I,J,L)  , &
                  NACTI(I,J,L)  , &
                   QST3(I,J,L)  )
             EVAPC(I,J,L) = ( Q(I,J,L) - EVAPC(I,J,L) ) / DT_MOIST
       ! sublimation for CN/LS
             RHCRIT = 1.0 - ALPHA
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
                  NACTL(I,J,L)  , &
                  NACTI(I,J,L)  , &
                   QST3(I,J,L)  )
             SUBLC(I,J,L) = ( Q(I,J,L) - SUBLC(I,J,L) ) / DT_MOIST
       ! Put condensates in touch with the PDF
             call hystpdf( &
                      DT_MOIST       , &
                      ALPHA          , &
                      PDFSHAPE       , &
                      CNV_FRC(I,J)   , &
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
                      WTHV2(I,J,L)   , &
                      WQL(I,J,L)     , &
                      USE_AEROSOL_NN)
       ! cleanup clouds
             call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
             RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
             endif
           end do ! IM loop
         end do ! JM loop
       end do ! LM loop

       ! Update TH
       TH = T/PK

    ! Update macrophysics tendencies
    if (associated( DUDT_macro))  DUDT_macro=( U         - DUDT_macro)/DT_MOIST
    if (associated( DVDT_macro))  DVDT_macro=( V         - DVDT_macro)/DT_MOIST
    if (associated( DTDT_macro))  DTDT_macro=( T         - DTDT_macro)/DT_MOIST
    if (associated(DTHDT_macro)) DTHDT_macro=( TH        -DTHDT_macro)/DT_MOIST
    if (associated(DQVDT_macro)) DQVDT_macro=( Q         -DQVDT_macro)/DT_MOIST
    if (associated(DQLDT_macro)) DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
    if (associated(DQIDT_macro)) DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
    if (associated(DQADT_macro)) DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
    if (associated(DQRDT_macro)) DQRDT_macro=( QRAIN     -DQRDT_macro)/DT_MOIST
    if (associated(DQSDT_macro)) DQSDT_macro=( QSNOW     -DQSDT_macro)/DT_MOIST
    if (associated(DQGDT_macro)) DQGDT_macro=( QGRAUPEL  -DQGDT_macro)/DT_MOIST
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
    call MAPL_GetPointer(EXPORT, DTHDT_micro, 'DTHDT_micro' , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    if (associated(DQVDT_micro)) DQVDT_micro = Q
    if (associated(DQLDT_micro)) DQLDT_micro = QLLS + QLCN
    if (associated(DQIDT_micro)) DQIDT_micro = QILS + QICN
    if (associated(DQRDT_micro)) DQRDT_micro = QRAIN
    if (associated(DQSDT_micro)) DQSDT_micro = QSNOW
    if (associated(DQGDT_micro)) DQGDT_micro = QGRAUPEL
    if (associated(DQADT_micro)) DQADT_micro = CLLS + CLCN
    if (associated( DUDT_micro))  DUDT_micro = U
    if (associated( DVDT_micro))  DVDT_micro = V
    if (associated( DTDT_micro))  DTDT_micro = T
    if (associated(DTHDT_micro)) DTHDT_micro = TH

        ! Delta-Z layer thickness (gfdl expects this to be negative)
         DZ = -1.0*DZET
        ! Zero-out local microphysics tendencies
         DQVDTmic = 0.
         DQLDTmic = 0.
         DQRDTmic = 0.
         DQIDTmic = 0.
         DQSDTmic = 0.
         DQGDTmic = 0.
         DQADTmic = 0.
          DUDTmic = 0.
          DVDTmic = 0.
          DTDTmic = 0.
       ! Zero-out 3D Precipitation Fluxes
        ! Ice
         PFI_LS = 0.
        ! Liquid
         PFL_LS = 0.
        ! Cloud
         RAD_CF = MIN(CLCN+CLLS,1.0)
        ! Liquid
         RAD_QL = QLCN+QLLS
        ! Ice
         RAD_QI = QICN+QILS
        ! VAPOR
         RAD_QV = Q
        ! RAIN
         RAD_QR = QRAIN
        ! SNOW
         RAD_QS = QSNOW
        ! GRAUPEL
         RAD_QG = QGRAUPEL
        ! Run the driver
         call gfdl_cloud_microphys_driver( &
                             ! Input water/cloud species and liquid+ice CCN [NACTL+NACTI]
                               RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI)/1.e6, &
                             ! Output tendencies
                               DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                               DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
                             ! Input fields
                               T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
                             ! constant inputs
                               AREA, DT_MOIST, FRLAND, CNV_FRC, &
                               ANV_ICEFALL, LS_ICEFALL, &
                             ! Output rain re-evaporation and sublimation
                               REV_LS, RSU_LS, &
                             ! Output precipitates
                               PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
                             ! Output mass flux during sedimentation (Pa kg/kg)
                               PFL_LS(:,:,1:LM), PFI_LS(:,:,1:LM), &
                             ! constant grid/time information
                               LHYDROSTATIC, LPHYS_HYDROSTATIC, &
                               1,IM, 1,JM, 1,LM, 1, LM)
     ! Apply tendencies
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
         RAD_CF = RAD_CF + DQADTmic * DT_MOIST
     ! Redistribute CN/LS CF/QL/QI
         call REDISTRIBUTE_CLOUDS(RAD_CF, RAD_QL, RAD_QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, T)
     ! Convert precip diagnostics from mm/day to kg m-2 s-1
         PRCP_RAIN    = MAX(PRCP_RAIN    / 86400.0, 0.0)
         PRCP_SNOW    = MAX(PRCP_SNOW    / 86400.0, 0.0)
         PRCP_ICE     = MAX(PRCP_ICE     / 86400.0, 0.0)
         PRCP_GRAUPEL = MAX(PRCP_GRAUPEL / 86400.0, 0.0)
     ! Fill GEOS precip diagnostics
         LS_PRCP = PRCP_RAIN
         LS_SNR  = PRCP_SNOW + PRCP_ICE + PRCP_GRAUPEL
     ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
         PFL_LS = PFL_LS/(MAPL_GRAV*DT_MOIST)
         PFI_LS = PFI_LS/(MAPL_GRAV*DT_MOIST)
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
              ! Get new TH
               TH(I,J,L) = T(I,J,L)/PK(I,J,L)
              ! get radiative properties
               call RADCOUPLE ( T(I,J,L), PLmb(I,J,L), CLLS(I,J,L), CLCN(I,J,L), &
                     Q(I,J,L), QLLS(I,J,L), QILS(I,J,L), QLCN(I,J,L), QICN(I,J,L), QRAIN(I,J,L), QSNOW(I,J,L), QGRAUPEL(I,J,L), NACTL(I,J,L), NACTI(I,J,L), &
                     RAD_QV(I,J,L), RAD_QL(I,J,L), RAD_QI(I,J,L), RAD_QR(I,J,L), RAD_QS(I,J,L), RAD_QG(I,J,L), RAD_CF(I,J,L), &
                     CLDREFFL(I,J,L), CLDREFFI(I,J,L), &
                     FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
               if (do_qa) RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
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
         ! Cloud fraction exports
         call MAPL_GetPointer(EXPORT, PTR3D, 'CFICE', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           PTR3D=0.0
           WHERE (RAD_QI .gt. 1.0e-12)
              PTR3D=RAD_CF*RAD_QI/(RAD_QL+RAD_QI)
           END WHERE
           PTR3D=MAX(MIN(PTR3D, 1.0), 0.0)
         endif
         call MAPL_GetPointer(EXPORT, PTR3D, 'CFLIQ', RC=STATUS); VERIFY_(STATUS)
         if (associated(PTR3D)) then
           PTR3D=0.0
           WHERE (RAD_QL .gt. 1.0e-12)
              PTR3D=RAD_CF*RAD_QL/(RAD_QL+RAD_QI)
           END WHERE
           PTR3D=MAX(MIN(PTR3D, 1.0), 0.0)
         endif
         ! Rain-out and Relative Humidity where RH > 110%
         DQST3 = GEOS_DQSAT(T, PLmb, QSAT=QST3)
         where ( Q > 1.1*QST3 )
            TMP3D = (Q - 1.1*QST3)/( 1.0 + 1.1*DQST3*MAPL_ALHL/MAPL_CP )
         elsewhere
            TMP3D = 0.0
         endwhere
         LS_PRCP = LS_PRCP + SUM(TMP3D*MASS,3)/DT_MOIST
         Q  =  Q - TMP3D
         T  =  T + (MAPL_ALHL/MAPL_CP)*TMP3D
         TH =  T / PK

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
         call FILLQ2ZERO(QRAIN   , MASS, TMP2D)
         call FILLQ2ZERO(QSNOW   , MASS, TMP2D)
         call FILLQ2ZERO(QGRAUPEL, MASS, TMP2D)

         ! Update microphysics tendencies
         if (associated(DQVDT_micro)) DQVDT_micro = ( Q          - DQVDT_micro) / DT_MOIST
         if (associated(DQLDT_micro)) DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
         if (associated(DQIDT_micro)) DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
         if (associated(DQADT_micro)) DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
         if (associated(DQRDT_micro)) DQRDT_micro = ( QRAIN      - DQRDT_micro) / DT_MOIST
         if (associated(DQSDT_micro)) DQSDT_micro = ( QSNOW      - DQSDT_micro) / DT_MOIST
         if (associated(DQGDT_micro)) DQGDT_micro = ( QGRAUPEL   - DQGDT_micro) / DT_MOIST
         if (associated( DUDT_micro))  DUDT_micro = ( U          -  DUDT_micro) / DT_MOIST
         if (associated( DVDT_micro))  DVDT_micro = ( V          -  DVDT_micro) / DT_MOIST
         if (associated( DTDT_micro))  DTDT_micro = ( T          -  DTDT_micro) / DT_MOIST
         if (associated(DTHDT_micro)) DTHDT_micro = ( TH         - DTHDT_micro) / DT_MOIST
        call MAPL_TimerOff(MAPL,"---CLDMICRO")

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
        if (associated(PTR2D)) PTR2D = SUM(   Q*MASS , 3 )

     call MAPL_TimerOff(MAPL,"--GFDL_1M",RC=STATUS)

end subroutine GFDL_1M_Run

   subroutine REDISTRIBUTE_CLOUDS(CF, QL, QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, QV, TE)
      real, dimension(:,:,:), intent(inout) :: CF, QL, QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, QV, TE

      WHERE (QL < 1.e-8)
        QV   = QV + QL
        TE   = TE - (MAPL_ALHL/MAPL_CP)*QL
        QL   = 0.0
        QLLS = 0.0
        QLCN = 0.0
      ELSE WHERE
        QLCN = MAX(0.0,MIN(QL,QLCN))
        QLLS = MAX(0.0,QL-QLCN)
      END WHERE

      WHERE (QI < 1.e-8)
        QV   = QV + QI
        TE   = TE - (MAPL_ALHS/MAPL_CP)*QI
        QI   = 0.0
        QILS = 0.0
        QICN = 0.0
      ELSE WHERE
        QICN = MAX(0.0,MIN(QI,QICN))
        QILS = MAX(0.0,QI-QICN)
      END WHERE

      WHERE ( (CF < 1.e-5) .or. (QL+QI < 1.e-8) )
        CF   = 0.0
        CLLS = 0.0
        CLCN = 0.0
        QV   = QV + QL + QI
        TE   = TE - (MAPL_ALHL/MAPL_CP)*QL - (MAPL_ALHS/MAPL_CP)*QI
        QL   = 0.0
        QLLS = 0.0
        QLCN = 0.0
        QI   = 0.0
        QILS = 0.0
        QICN = 0.0
      ELSE WHERE
        CLCN = MAX(0.0,MIN(CF,CLCN,1.0))
        CLLS = MAX(0.0,MIN(CF-CLCN,1.0))
      END WHERE

   end subroutine REDISTRIBUTE_CLOUDS

   subroutine FIX_NEGATIVE_PRECIP(QRAIN, QSNOW, QGRAUPEL)
      real, dimension(:,:,:), intent(inout) :: QRAIN, QSNOW, QGRAUPEL

      WHERE (QRAIN < 1.e-8)
        QRAIN = 0.0
      END WHERE

      WHERE (QSNOW < 1.e-8)
        QSNOW = 0.0
      END WHERE

      WHERE (QGRAUPEL < 1.e-8)
        QGRAUPEL = 0.0
      END WHERE

   end subroutine FIX_NEGATIVE_PRECIP

end module GEOS_GFDL_1M_InterfaceMod
