#define ACC_PREFIX !$acc
! $Id$

#include "MAPL_Generic.h"
!#define PDFDIAG 1

!=============================================================================
!BOP

! !MODULE: GEOS_GFDL_1M_InterfaceMod -- A Module to interface with the
!   GFDL_1M cloud microphysics

module GEOS_GFDL_1M_InterfaceMod

#ifdef SERIALIZE
USE m_serialize, ONLY: &
  fs_add_savepoint_metainfo, &
  fs_create_savepoint, &
  fs_read_field, &
  fs_write_field
USE utils_ppser, ONLY:  &
  ppser_get_mode, &
  ppser_savepoint, &
  ppser_serializer, &
  ppser_serializer_ref, &
  ppser_intlength, &
  ppser_reallength, &
  ppser_realtype, &
  ppser_zrperturb, &
  ppser_get_mode
USE savepoint_helpers
USE utils_ppser_ijkbuff
USE utils_ppser_kbuff
#endif


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
  logical :: LHYDROSTATIC
  logical :: LPHYS_HYDROSTATIC
  logical :: LMELTFRZ

  integer, save :: timestep = 0

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

    call MAPL_TimerAdd(GC, name="--GFDL_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine GFDL_1M_Setup

subroutine GFDL_1M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    type (ESMF_Alarm   )                :: ALARM
    type (ESMF_TimeInterval)            :: TINT
    real(ESMF_KIND_R8)                  :: DT_R8
    real                                :: DT_MOIST

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL

    type(ESMF_VM) :: VM
    integer :: comm

    call MAPL_GetResource( MAPL, LHYDROSTATIC, Label="HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LPHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LMELTFRZ, Label="MELTFRZ:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_Get ( MAPL, INTERNAL_ESMF_STATE=INTERNAL, RC=STATUS )
    VERIFY_(STATUS)
    call MAPL_Get( MAPL, &
         RUNALARM = ALARM,             &
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
    call MAPL_GetPointer(INTERNAL, QILS,     'QILS'    , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(INTERNAL, QICN,     'QICN'    , RC=STATUS); VERIFY_(STATUS)

    call ESMF_VMGetCurrent(VM, _RC)
    call ESMF_VMGet(VM, mpiCommunicator=comm, _RC)

    call gfdl_cloud_microphys_init(comm)
    call WRITE_PARALLEL ("INITIALIZED GFDL_1M microphysics in non-generic GC INIT")

    call MAPL_GetResource( MAPL, SH_MD_DP        , 'SH_MD_DP:'        , DEFAULT= .TRUE., RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, TURNRHCRIT_PARAM, 'TURNRHCRIT:'      , DEFAULT= -9999., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, PDFSHAPE        , 'PDFSHAPE:'        , DEFAULT= 1     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ANV_ICEFALL     , 'ANV_ICEFALL:'     , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LS_ICEFALL      , 'LS_ICEFALL:'      , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LIQ_RADII_PARAM , 'LIQ_RADII_PARAM:' , DEFAULT= 2     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ICE_RADII_PARAM , 'ICE_RADII_PARAM:' , DEFAULT= 1     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RI          , 'FAC_RI:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RI          , 'MIN_RI:'          , DEFAULT=  5.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RI          , 'MAX_RI:'          , DEFAULT=100.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RL          , 'FAC_RL:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RL          , 'MIN_RL:'          , DEFAULT= 2.5e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RL          , 'MAX_RL:'          , DEFAULT=60.0e-6, RC=STATUS); VERIFY_(STATUS)

                                 CCW_EVAP_EFF = 1.e-2
                    if (do_evap) CCW_EVAP_EFF = 0.0 ! Evap done inside GFDL-MP
    call MAPL_GetResource( MAPL, CCW_EVAP_EFF, 'CCW_EVAP_EFF:', DEFAULT= CCW_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

                                 CCI_EVAP_EFF = 1.e-2
                    if (do_subl) CCI_EVAP_EFF = 0.0 ! Subl done inside GFDL-MP
    call MAPL_GetResource( MAPL, CCI_EVAP_EFF, 'CCI_EVAP_EFF:', DEFAULT= CCI_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=  500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_EXP, 'CNV_FRACTION_EXP:', DEFAULT=    1.0, RC=STATUS); VERIFY_(STATUS)

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
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, W, KH
    real, pointer, dimension(:,:)   :: AREA, FRLAND, TS, DTSX, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:) :: SL2, SL3, QT2, QT3, W2, W3, SLQT, WQT, WQL, WSL, PDF_A
    real, pointer, dimension(:,:,:) :: WTHV2
    real, pointer, dimension(:,:,:) :: OMEGA
    ! Local
    real, allocatable, dimension(:,:,:) :: U0, V0
    real, allocatable, dimension(:,:,:) :: PLEmb, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  ZL0
    real, allocatable, dimension(:,:,:) :: DZ, DZET, DP, MASS, iMASS
    real, allocatable, dimension(:,:,:) :: DQST3, QST3
    real, allocatable, dimension(:,:,:) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                                           DQSDTmic, DQGDTmic, DQADTmic, &
                                            DUDTmic,  DVDTmic,  DTDTmic
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D
    integer, allocatable, dimension(:,:) :: KLCL
    ! Exports
    real, pointer, dimension(:,:  ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ) :: LS_PRCP, LS_SNR, ICE, FRZR, CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:) ::  DUDT_macro,  DVDT_macro,  DTDT_macro
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro
    real, pointer, dimension(:,:,:) ::  DUDT_micro,  DVDT_micro,  DTDT_micro
    real, pointer, dimension(:,:,:) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, pointer, dimension(:,:,:) :: CLDREFFL, CLDREFFI
    real, pointer, dimension(:,:,:) :: EVAPC, SUBLC
    real, pointer, dimension(:,:,:) :: RHX, REV_LS, RSU_LS
    real, pointer, dimension(:,:,:) :: PFL_LS, PFL_AN
    real, pointer, dimension(:,:,:) :: PFI_LS, PFI_AN
    real, pointer, dimension(:,:,:) :: PDFITERS
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

    ! Local variables
    real    :: facEIS
    real    :: minrhcrit, turnrhcrit, ALPHA, RHCRIT
    integer :: IM,JM,LM
    integer :: I, J, L

    integer :: comm, rank, mpierr
    
    type(ESMF_VM) :: vm
    call ESMF_VMGetCurrent(vm,rc=status)
    call ESMF_VMGet(VM, mpiCommunicator=comm, rc=status)
    call MPI_COMM_rank(comm,rank,mpierr)

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
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          QRAIN = QRAIN + PTR3D*DT_MOIST
        endif
        call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then 
          QSNOW = QSNOW + PTR3D*DT_MOIST
        endif
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #587
call fs_create_savepoint('Evap_subl_pdf-In', ppser_savepoint)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #588
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'EIS', EIS)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'EIS', EIS)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'EIS', EIS, ppser_zrperturb)
END SELECT
#endif
       ! evap/subl/pdf
        call MAPL_GetPointer(EXPORT, RHCRIT3D,  'RHCRIT', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
#ifdef SERIALIZE
call ser_set_nz(LM)
call ser_set_ny(JM)
call ser_set_nx(IM)
#endif
        do L=1,LM
          do J=1,JM
           do I=1,IM
#ifdef SERIALIZE
call ser_set_i(I)
call ser_set_j(J)
call ser_set_k(L)
#endif
           
           ! Send the condensates through the pdf after convection
             facEIS = MAX(0.0,MIN(1.0,EIS(I,J)/10.0))**2
           ! determine combined minrhcrit in stable/unstable regimes
             minrhcrit  = (1.0-dw_ocean)*(1.0-facEIS) + (1.0-dw_land)*facEIS
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
             if (.not. do_qa) then ! if not doing cloud pdf inside of GFDL-MP
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
                      USE_BERGERON)
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #681
call fs_create_savepoint('IJKBUFF_TEST-In', ppser_savepoint)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #682
    call fs_write_ijkbuff(ppser_serializer, ppser_savepoint, 'RHX', RHX(I,J,L), i=ser_i, i_size=ser_nx, j=ser_j, j_size=ser_ny, k=ser_k, k_size=ser_nz, mode=ppser_get_mode())
    call fs_write_ijkbuff(ppser_serializer, ppser_savepoint, 'Q', Q(I,J,L), i=ser_i, i_size=ser_nx, j=ser_j, j_size=ser_ny, k=ser_k, k_size=ser_nz, mode=ppser_get_mode())
    call fs_write_ijkbuff(ppser_serializer, ppser_savepoint, 'T', T(I,J,L), i=ser_i, i_size=ser_nx, j=ser_j, j_size=ser_ny, k=ser_k, k_size=ser_nz, mode=ppser_get_mode())
    call fs_write_ijkbuff(ppser_serializer, ppser_savepoint, 'PLmb', PLmb(I,J,L), i=ser_i, i_size=ser_nx, j=ser_j, j_size=ser_ny, k=ser_k, k_size=ser_nz, mode=ppser_get_mode())
#endif
             RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #684
call fs_create_savepoint('IJKBUFF_TEST-Out', ppser_savepoint)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #685
    call fs_write_ijkbuff(ppser_serializer, ppser_savepoint, 'RHX', RHX(I,J,L), i=ser_i, i_size=ser_nx, j=ser_j, j_size=ser_ny, k=ser_k, k_size=ser_nz, mode=ppser_get_mode())
#endif
            endif
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
             if (CCW_EVAP_EFF > 0.0) then ! else evap done inside GFDL
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
             endif
           ! sublimation for CN
             if (CCI_EVAP_EFF > 0.0) then ! else subl done inside GFDL
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
                  NACTL(I,J,L)  , &
                  NACTI(I,J,L)  , &
                   QST3(I,J,L)  )
             SUBLC(I,J,L) = ( Q(I,J,L) - SUBLC(I,J,L) ) / DT_MOIST
             endif
           ! cleanup clouds
             call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
           end do ! IM loop
         end do ! JM loop
       end do ! LM loop
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #745
call fs_create_savepoint('Evap_subl_pdf-Out', ppser_savepoint)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #746
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'SUBLC', SUBLC)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'SUBLC', SUBLC)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'SUBLC', SUBLC, ppser_zrperturb)
END SELECT
#endif

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
                             ! Input water/cloud species and liquid+ice CCN [NACTL+NACTI (#/m^3)]
                               RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI), &
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
         RAD_CF = MIN(1.0,MAX(0.0,RAD_CF + DQADTmic * DT_MOIST))
         ! Serialize data for RedistributeClouds
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #850
call fs_create_savepoint('RedistributeClouds-In', ppser_savepoint)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #851
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLCN', CLCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLLS', CLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLCN', QLCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLLS', QLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QICN', QICN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QILS', QILS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'T', T)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T, ppser_zrperturb)
END SELECT
#endif
     ! Redistribute CN/LS CF/QL/QI
         call REDISTRIBUTE_CLOUDS(RAD_CF, RAD_QL, RAD_QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, T)
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #854
call fs_create_savepoint('RedistributeClouds-Out', ppser_savepoint)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #855
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLCN', CLCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLLS', CLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLCN', QLCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLLS', QLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QICN', QICN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QILS', QILS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'T', T)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T, ppser_zrperturb)
END SELECT
#endif
     ! Convert precip diagnostics from mm/day to kg m-2 s-1
         PRCP_RAIN    = MAX(PRCP_RAIN    / 86400.0, 0.0)
         PRCP_SNOW    = MAX(PRCP_SNOW    / 86400.0, 0.0)
         PRCP_ICE     = MAX(PRCP_ICE     / 86400.0, 0.0)
         PRCP_GRAUPEL = MAX(PRCP_GRAUPEL / 86400.0, 0.0)
     ! Fill GEOS precip diagnostics
         LS_PRCP = PRCP_RAIN
         LS_SNR  = PRCP_SNOW
         ICE     = PRCP_ICE + PRCP_GRAUPEL
         FRZR    = 0.0
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

       timestep = timestep + 1
      !  print*, "TIMESTEP = ", timestep, ' rank = ', rank, ' do_qa = ', do_qa, "6 * rank + TIMESTEP = ", 6*rank+timestep
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #886
call fs_create_savepoint('RadCouple-In', ppser_savepoint)
call fs_add_savepoint_metainfo(ppser_savepoint, 'timestep', timestep)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #887
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'Q', Q)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'T', T)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLLS', QLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QILS', QILS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLLS', CLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLCN', QLCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QICN', QICN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLCN', CLCN)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'Q', Q)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'Q', Q, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN, ppser_zrperturb)
END SELECT
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #888
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'PLmb', PLmb)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QRAIN', QRAIN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QSNOW', QSNOW)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QGRAUPEL', QGRAUPEL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'NACTL', NACTL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'NACTI', NACTI)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PLmb', PLmb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QRAIN', QRAIN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QSNOW', QSNOW)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QGRAUPEL', QGRAUPEL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'NACTL', NACTL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'NACTI', NACTI)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'PLmb', PLmb, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QRAIN', QRAIN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QSNOW', QSNOW, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QGRAUPEL', QGRAUPEL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'NACTL', NACTL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'NACTI', NACTI, ppser_zrperturb)
END SELECT
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #889
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'FAC_RL', FAC_RL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'MIN_RL', MIN_RL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'MAX_RL', MAX_RL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'FAC_RI', FAC_RI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'MIN_RI', MIN_RI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'MAX_RI', MAX_RI)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'FAC_RL', FAC_RL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MIN_RL', MIN_RL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MAX_RL', MAX_RL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'FAC_RI', FAC_RI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MIN_RI', MIN_RI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MAX_RI', MAX_RI)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'FAC_RL', FAC_RL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MIN_RL', MIN_RL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MAX_RL', MAX_RL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'FAC_RI', FAC_RI, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MIN_RI', MIN_RI, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'MAX_RI', MAX_RI, ppser_zrperturb)
END SELECT
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #890
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QR', RAD_QR)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QS', RAD_QS)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QR', RAD_QR)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QS', RAD_QS)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QR', RAD_QR, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QS', RAD_QS, ppser_zrperturb)
END SELECT
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #891
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QG', RAD_QG)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLDREFFL', CLDREFFL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLDREFFI', CLDREFFI)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QG', RAD_QG)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFL', CLDREFFL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFI', CLDREFFI)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QG', RAD_QG, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFL', CLDREFFL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFI', CLDREFFI, ppser_zrperturb)
END SELECT
#endif

     ! Radiation Coupling
         do L = 1, LM
           do J = 1, JM
             do I = 1, IM
              ! cleanup clouds
               call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
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
#ifdef SERIALIZE
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #909
call fs_create_savepoint('RadCouple-Out', ppser_savepoint)
call fs_add_savepoint_metainfo(ppser_savepoint, 'timestep', timestep)
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #910
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'Q', Q)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'T', T)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLLS', QLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QILS', QILS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLLS', CLLS)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QLCN', QLCN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'QICN', QICN)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLCN', CLCN)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'Q', Q)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'Q', Q, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'T', T, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLLS', QLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QILS', QILS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLLS', CLLS, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QLCN', QLCN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'QICN', QICN, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLCN', CLCN, ppser_zrperturb)
END SELECT
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #911
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QR', RAD_QR)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QS', RAD_QS)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QR', RAD_QR)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QS', RAD_QS)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QV', RAD_QV, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QL', RAD_QL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QI', RAD_QI, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QR', RAD_QR, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QS', RAD_QS, ppser_zrperturb)
END SELECT
! file: /home/fgdeconi/work/git/fp/geos/src/Components/@GEOSgcm_GridComp/GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSmoist_GridComp/GEOS_GFDL_1M_InterfaceMod.F90.SER lineno: #912
SELECT CASE ( ppser_get_mode() )
  CASE(0)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_QG', RAD_QG)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLDREFFL', CLDREFFL)
    call fs_write_field(ppser_serializer, ppser_savepoint, 'CLDREFFI', CLDREFFI)
  CASE(1)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QG', RAD_QG)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFL', CLDREFFL)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFI', CLDREFFI)
  CASE(2)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_QG', RAD_QG, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'RAD_CF', RAD_CF, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFL', CLDREFFL, ppser_zrperturb)
    call fs_read_field(ppser_serializer_ref, ppser_savepoint, 'CLDREFFI', CLDREFFI, ppser_zrperturb)
END SELECT
#endif

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

            call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,1)
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

     call MAPL_TimerOff(MAPL,"--GFDL_1M",RC=STATUS)

end subroutine GFDL_1M_Run

end module GEOS_GFDL_1M_InterfaceMod

