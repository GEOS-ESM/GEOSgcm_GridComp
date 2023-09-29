! $Id$

#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: GEOS_THOM_1M_InterfaceMod -- A Module to interface with the
!   THOM_1M cloud microphysics

module GEOS_THOM_1M_InterfaceMod

  use ESMF
  use MAPL
  use GEOS_UtilsMod
  use GEOSmoist_Process_Library
  use Aer_Actv_Single_Moment
  use module_mp_thompson

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
         character(len=ESMF_MAXSTR) :: NCPL                     
         character(len=ESMF_MAXSTR) :: NCPI
         character(len=ESMF_MAXSTR) :: NRAIN
  end type FRIENDLIES_TYPE
  type (FRIENDLIES_TYPE) FRIENDLIES

  character(len=ESMF_MAXSTR)        :: COMP_NAME

  ! Local resource variables
  real    :: DT_THOM
  real    :: TURNRHCRIT_PARAM
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

  public :: THOM_1M_Setup, THOM_1M_Initialize, THOM_1M_Run

contains

subroutine THOM_1M_Setup (GC, CF, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(ESMF_Config),   intent(inout) :: CF
    integer, optional                  :: RC  ! return code
    character(len=ESMF_MAXSTR)         :: COMP_NAME

    IAm = "GEOS_THOM_1M_InterfaceMod"
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
      FRIENDLIES%NCPL     = "DYNAMICS:TURBULENCE"                  
      FRIENDLIES%NCPI     = "DYNAMICS:TURBULENCE"
      FRIENDLIES%NRAIN    = "DYNAMICS:TURBULENCE"

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

    call MAPL_TimerAdd(GC, name="--THOM_1M", RC=STATUS)
    VERIFY_(STATUS)

end subroutine THOM_1M_Setup

subroutine THOM_1M_Initialize (MAPL, RC)
    type (MAPL_MetaComp), intent(inout) :: MAPL
    integer, optional                   :: RC  ! return code

    type (ESMF_Grid )                   :: GRID
    type (ESMF_State)                   :: INTERNAL

    CHARACTER(len=ESMF_MAXSTR) :: errmsg

    real, pointer, dimension(:,:,:)     :: Q, QLLS, QLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL

    call MAPL_GetResource( MAPL, LHYDROSTATIC, Label="HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LPHYS_HYDROSTATIC, Label="PHYS_HYDROSTATIC:",  default=.TRUE., RC=STATUS)
    VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, DT_THOM, Label="DT_THOM:",  default=300.0, RC=STATUS)
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

    call thompson_init(.false., USE_AEROSOL_NN, MAPL_am_I_root() , 1, errmsg, STATUS)
    _ASSERT( STATUS==0, errmsg )
    call WRITE_PARALLEL ("INITIALIZED THOM_1M microphysics in non-generic GC INIT")

    call MAPL_GetResource( MAPL, TURNRHCRIT_PARAM, 'TURNRHCRIT:'      , DEFAULT= -9999., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, PDFSHAPE        , 'PDFSHAPE:'        , DEFAULT= 1     , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, ANV_ICEFALL     , 'ANV_ICEFALL:'     , DEFAULT= 0.8   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, LS_ICEFALL      , 'LS_ICEFALL:'      , DEFAULT= 0.8   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RI          , 'FAC_RI:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RI          , 'MIN_RI:'          , DEFAULT=  5.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RI          , 'MAX_RI:'          , DEFAULT=140.e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, FAC_RL          , 'FAC_RL:'          , DEFAULT= 1.0   , RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MIN_RL          , 'MIN_RL:'          , DEFAULT= 2.5e-6, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, MAX_RL          , 'MAX_RL:'          , DEFAULT=60.0e-6, RC=STATUS); VERIFY_(STATUS)

                                 CCW_EVAP_EFF = 4.e-3
    call MAPL_GetResource( MAPL, CCW_EVAP_EFF, 'CCW_EVAP_EFF:', DEFAULT= CCW_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

                                 CCI_EVAP_EFF = 4.e-3
    call MAPL_GetResource( MAPL, CCI_EVAP_EFF, 'CCI_EVAP_EFF:', DEFAULT= CCI_EVAP_EFF, RC=STATUS); VERIFY_(STATUS)

    call MAPL_GetResource( MAPL, CNV_FRACTION_MIN, 'CNV_FRACTION_MIN:', DEFAULT=    0.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_MAX, 'CNV_FRACTION_MAX:', DEFAULT= 1500.0, RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetResource( MAPL, CNV_FRACTION_EXP, 'CNV_FRACTION_EXP:', DEFAULT=    0.5, RC=STATUS); VERIFY_(STATUS)

end subroutine THOM_1M_Initialize

subroutine THOM_1M_Run (GC, IMPORT, EXPORT, CLOCK, RC)
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
    real, pointer, dimension(:,:,:) :: NCPL, NCPI, NRAIN, NACTL, NACTI
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
    real, allocatable, dimension(:,:,:) :: VVEL, DELZ, DP, MASS, iMASS
    real, allocatable, dimension(:,:,:) :: DQST3, QST3
    real, allocatable, dimension(:,:,:) :: AIRDEN
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: TMP2D
    integer, allocatable, dimension(:,:) :: KLCL
    integer, allocatable, dimension(:,:) :: iLand2D
    ! Exports
    real, pointer, dimension(:,:  ) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, SNOW_RATIO
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
    real, pointer, dimension(:,:,:) :: DBZ3D
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D
    ! Thompson Pointers for inputs
    real, dimension(:,:,:), allocatable, target  :: inputs
    real, dimension(:,:,:), pointer :: qv   => null()
    real, dimension(:,:,:), pointer :: qc   => null()
    real, dimension(:,:,:), pointer :: qr   => null()
    real, dimension(:,:,:), pointer :: qi   => null()
    real, dimension(:,:,:), pointer :: qs   => null()
    real, dimension(:,:,:), pointer :: qg   => null()
    real, dimension(:,:,:), pointer :: tt   => null()
    real, dimension(:,:,:), pointer :: pp   => null()
    real, dimension(:,:,:), pointer :: ww   => null()
    real, dimension(:,:,:), pointer :: dz   => null()
    real, dimension(:,:,:), pointer :: nc   => null()
    real, dimension(:,:,:), pointer :: ni   => null()
    real, dimension(:,:,:), pointer :: nr   => null()
    real, dimension(:,:,:), pointer :: nwfa => null()
    real, dimension(:,:,:), pointer :: nifa => null()
    ! Thompson Pointers for extended diags
    real, dimension(:,:,:), allocatable, target  :: diag3d
    real, dimension(:,:,:), pointer :: prw_vcdc   => null()
    real, dimension(:,:,:), pointer :: prw_vcde   => null()
    real, dimension(:,:,:), pointer :: tpri_inu   => null()
    real, dimension(:,:,:), pointer :: tpri_ide_d => null()
    real, dimension(:,:,:), pointer :: tpri_ide_s => null()
    real, dimension(:,:,:), pointer :: tprs_ide   => null()
    real, dimension(:,:,:), pointer :: tprs_sde_d => null()
    real, dimension(:,:,:), pointer :: tprs_sde_s => null()
    real, dimension(:,:,:), pointer :: tprg_gde_d => null()
    real, dimension(:,:,:), pointer :: tprg_gde_s => null()
    real, dimension(:,:,:), pointer :: tpri_iha   => null()
    real, dimension(:,:,:), pointer :: tpri_wfz   => null()
    real, dimension(:,:,:), pointer :: tpri_rfz   => null()
    real, dimension(:,:,:), pointer :: tprg_rfz   => null()
    real, dimension(:,:,:), pointer :: tprs_scw   => null()
    real, dimension(:,:,:), pointer :: tprg_scw   => null()
    real, dimension(:,:,:), pointer :: tprg_rcs   => null()
    real, dimension(:,:,:), pointer :: tprs_rcs   => null()
    real, dimension(:,:,:), pointer :: tprr_rci   => null()
    real, dimension(:,:,:), pointer :: tprg_rcg   => null()
    real, dimension(:,:,:), pointer :: tprw_vcd_c => null()
    real, dimension(:,:,:), pointer :: tprw_vcd_e => null()
    real, dimension(:,:,:), pointer :: tprr_sml   => null()
    real, dimension(:,:,:), pointer :: tprr_gml   => null()
    real, dimension(:,:,:), pointer :: tprr_rcg   => null()
    real, dimension(:,:,:), pointer :: tprr_rcs   => null()
    real, dimension(:,:,:), pointer :: tprv_rev   => null()
    real, dimension(:,:,:), pointer :: tten3      => null()
    real, dimension(:,:,:), pointer :: qvten3     => null()
    real, dimension(:,:,:), pointer :: qrten3     => null()
    real, dimension(:,:,:), pointer :: qsten3     => null()
    real, dimension(:,:,:), pointer :: qgten3     => null()
    real, dimension(:,:,:), pointer :: qiten3     => null()
    real, dimension(:,:,:), pointer :: niten3     => null()
    real, dimension(:,:,:), pointer :: nrten3     => null()
    real, dimension(:,:,:), pointer :: ncten3     => null()
    real, dimension(:,:,:), pointer :: qcten3     => null()

    ! Local variables
    real    :: facEIS
    real    :: minrhcrit, turnrhcrit, ALPHA, RHCRIT
    integer :: IM,JM,LM
    integer :: I, J, L
    CHARACTER(len=ESMF_MAXSTR) :: errmsg

    ! Effective cloud radii - turned off in CCPP (taken care off in radiation)
    logical, parameter :: do_effective_radii = .false.
    integer, parameter :: has_reqc = 0
    integer, parameter :: has_reqi = 0
    integer, parameter :: has_reqs = 0
    ! Stochastic physics off
    integer, parameter :: kme_stoch = 1
    integer            :: spp_mp_opt = 0 
    integer            :: spp_mp = 0
    integer            :: n_var_spp
    real, dimension(:,:), pointer :: spp_wts_mp => null()
    real, dimension(:)  , pointer :: spp_prt_list => null()
    real, dimension(:)  , pointer :: spp_stddev_cutoff => null()
    character(len=3),dimension(:)  , pointer :: spp_var_list => null()

    call ESMF_GridCompGet( GC, CONFIG=CF, RC=STATUS ) 
    VERIFY_(STATUS)

    ! Get my internal MAPL_Generic state
    !-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    VERIFY_(STATUS)

    call MAPL_TimerOn (MAPL,"--THOM_1M",RC=STATUS)

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
    ALLOCATE ( DELZ (IM,JM,LM  ) )
    ALLOCATE ( DP   (IM,JM,LM  ) )
    ALLOCATE ( MASS (IM,JM,LM  ) )
    ALLOCATE ( iMASS(IM,JM,LM  ) )
    ALLOCATE ( DQST3(IM,JM,LM  ) )
    ALLOCATE (  QST3(IM,JM,LM  ) )
    ALLOCATE (AIRDEN(IM,JM,LM  ) )
    ALLOCATE (  VVEL(IM,JM,LM  ) )
    ALLOCATE ( TMP3D(IM,JM,LM  ) )
     ! 2D Variables
    ALLOCATE ( iLand2D      (IM,JM) ) 
               iLand2D = NINT(FRLAND)
    ALLOCATE ( KLCL         (IM,JM) )
    ALLOCATE ( TMP2D        (IM,JM) )

    ! INOUT Arrays
    ALLOCATE ( inputs(IM*JM,LM,15) )
    inputs = 0.0
    qv => inputs(:,:,1:1)
    qc => inputs(:,:,2:2)
    qr => inputs(:,:,3:3)
    qi => inputs(:,:,4:4)
    qs => inputs(:,:,5:5)
    qg => inputs(:,:,6:6)
    tt => inputs(:,:,7:7)
    pp => inputs(:,:,8:8)
    ww => inputs(:,:,9:9)
    dz => inputs(:,:,10:10)
    nc => inputs(:,:,11:11)
    ni => inputs(:,:,12:12)
    nr => inputs(:,:,13:13)
    nwfa => inputs(:,:,14:14)
    nifa => inputs(:,:,15:15)

    ! Extended diagnostics
    ALLOCATE ( diag3d(IM*JM,LM,37) )
    diag3d = 0.0
    prw_vcdc   => diag3d(:,:,1:1)
    prw_vcde   => diag3d(:,:,2:2)
    tpri_inu   => diag3d(:,:,3:3)
    tpri_ide_d => diag3d(:,:,4:4)
    tpri_ide_s => diag3d(:,:,5:5)
    tprs_ide   => diag3d(:,:,6:6)
    tprs_sde_d => diag3d(:,:,7:7)
    tprs_sde_s => diag3d(:,:,8:8)
    tprg_gde_d => diag3d(:,:,9:9)
    tprg_gde_s => diag3d(:,:,10:10)
    tpri_iha   => diag3d(:,:,11:11)
    tpri_wfz   => diag3d(:,:,12:12)
    tpri_rfz   => diag3d(:,:,13:13)
    tprg_rfz   => diag3d(:,:,14:14)
    tprs_scw   => diag3d(:,:,15:15)
    tprg_scw   => diag3d(:,:,16:16)
    tprg_rcs   => diag3d(:,:,17:17)
    tprs_rcs   => diag3d(:,:,18:18)
    tprr_rci   => diag3d(:,:,19:19)
    tprg_rcg   => diag3d(:,:,20:20)
    tprw_vcd_c => diag3d(:,:,21:21)
    tprw_vcd_e => diag3d(:,:,22:22)
    tprr_sml   => diag3d(:,:,23:23)
    tprr_gml   => diag3d(:,:,24:24)
    tprr_rcg   => diag3d(:,:,25:25)
    tprr_rcs   => diag3d(:,:,26:26)
    tprv_rev   => diag3d(:,:,27:27)
    tten3      => diag3d(:,:,28:28)
    qvten3     => diag3d(:,:,29:29)
    qrten3     => diag3d(:,:,30:30)
    qsten3     => diag3d(:,:,31:31)
    qgten3     => diag3d(:,:,32:32)
    qiten3     => diag3d(:,:,33:33)
    niten3     => diag3d(:,:,34:34)
    nrten3     => diag3d(:,:,35:35)
    ncten3     => diag3d(:,:,36:36)
    qcten3     => diag3d(:,:,37:37)

    ! Derived States
    PLEmb    =  PLE*.01
    PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
    DO L=0,LM
       ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
    END DO
    ZL0      = 0.5*(ZLE0(:,:,1:LM) + ZLE0(:,:,0:LM-1) ) ! Layer Height (m) above the surface
    DELZ     =     (ZLE0(:,:,1:LM) - ZLE0(:,:,0:LM-1) ) ! Layer thickness (m)
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
    call MAPL_GetPointer(EXPORT, DBZ3D,    'DBZ'     , ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    call MAPL_GetPointer(EXPORT, SNOW_RATIO,'SNOW_RATIO', ALLOC=.TRUE., RC=STATUS); VERIFY_(STATUS)
    ! Unused Exports (foreced to 0.0)
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
             minrhcrit  = (1.0-0.1)*(1.0-facEIS) + (1.0-0.05)*facEIS
             if (TURNRHCRIT_PARAM <= 0.0) then
              ! determine the turn pressure using the LCL
                turnrhcrit  = PLmb(I, J, KLCL(I,J)) - 250.0 ! 250mb above the LCL
             else
                turnrhcrit  = TURNRHCRIT_PARAM
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
                      WTHV2(I,J,L)   , &
                      WQL(I,J,L)     , &
                      .false.        , & 
                      USE_BERGERON)
             RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
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
           ! evaporation for CN
             if (CCW_EVAP_EFF > 0.0) then ! else evap done inside MP
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
             if (CCI_EVAP_EFF > 0.0) then ! else subl done inside MP
             RHCRIT = 1.0 - ALPHA
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
         RAD_QV = Q/(1.0-Q)
        ! RAIN
         RAD_QR = QRAIN
        ! SNOW
         RAD_QS = QSNOW
        ! GRAUPEL
         RAD_QG = QGRAUPEL
        ! Air Density
         AIRDEN   = MAPL_EPSILON*100.*PLmb/(T*MAPL_RGAS*(Q+MAPL_EPSILON))
        ! Vertical velocity
         if (LHYDROSTATIC) then 
           VVEL = -OMEGA/(MAPL_GRAV*AIRDEN)
         else 
           VVEL = W                
         endif                  
        ! RESHAPE
         qv = RESHAPE(RAD_QV,(/IM*JM,LM,1/))
         qc = RESHAPE(RAD_QL,(/IM*JM,LM,1/))
         qr = RESHAPE(RAD_QR,(/IM*JM,LM,1/))
         qi = RESHAPE(RAD_QI,(/IM*JM,LM,1/))
         qs = RESHAPE(RAD_QS,(/IM*JM,LM,1/))
         qg = RESHAPE(RAD_QG,(/IM*JM,LM,1/))
         tt = RESHAPE(T,(/IM*JM,LM,1/))
         pp = RESHAPE(PLmb*100.0,(/IM*JM,LM,1/))
         ww = RESHAPE(VVEL,(/IM*JM,LM,1/))
         dz = RESHAPE(DELZ,(/IM*JM,LM,1/))
        ! aerosols
         nwfa = RESHAPE(NACTL/AIRDEN,(/IM*JM,LM,1/))
         nifa = RESHAPE(NACTI/AIRDEN,(/IM*JM,LM,1/))
        ! Ensure we have 1st guess ice number where mass non-zero but no number.
         if (all(NCPL .eq. 0.0)) NCPL  = make_DropletNumber (RAD_QL*AIRDEN, NACTL) / AIRDEN
         nc = RESHAPE(NCPL,(/IM*JM,LM,1/))
        ! Ensure we have 1st guess ice number where mass non-zero but no number.
         if (all(NCPI .eq. 0.0)) NCPI  = make_IceNumber (RAD_QI*AIRDEN, T) / AIRDEN
         ni = RESHAPE(NCPI,(/IM*JM,LM,1/))
        ! Ensure we have 1st guess rain number where mass non-zero but no number.
         if (all(NRAIN .eq. 0.0)) NRAIN = make_RainNumber(RAD_QR*AIRDEN, T) / AIRDEN
         nr = RESHAPE(NRAIN,(/IM*JM,LM,1/))
        ! Run the driver
         call mp_gt_driver(qv=qv, qc=qc, qr=qr, qi=qi, qs=qs, qg=qg, nc=nc, ni=ni, nr=nr, &
                              nwfa=nwfa, nifa=nifa, &
                            ! nwfa2d=nwfa2d, nifa2d=nifa2d,     &
                              tt=tt, p=pp, w=ww, dz=dz, dt_in=DT_MOIST, dt_inner=min(DT_THOM,DT_MOIST),  &
                              sedi_semi=.TRUE., decfl=10, lsm=iLand2D,                  &
                                 rainnc=PRCP_RAIN,      &
                                 snownc=PRCP_SNOW,      &
                                  icenc=PRCP_ICE,       &
                              graupelnc=PRCP_GRAUPEL, sr=SNOW_RATIO, &
                              ids=1, ide=IM*JM, jds=1, jde=1, kds=1, kde=LM,   &
                              ims=1, ime=IM*JM, jms=1, jme=1, kms=1, kme=LM,   &
                              its=1, ite=IM*JM, jts=1, jte=1, kts=1, kte=LM,   &
                              istep=1, nsteps=1, &
                              diagflag=.TRUE., &
                              do_radar_ref=.TRUE., fullradar_diag=.TRUE., first_time_step=.TRUE., &
                              refl_10cm=DBZ3D, &
                              has_reqc=has_reqc, has_reqi=has_reqi, has_reqs=has_reqs,       &
                              rand_perturb_on=spp_mp_opt, kme_stoch=kme_stoch, n_var_spp=n_var_spp,&
                          !   rand_pert=spp_wts_mp, spp_var_list=spp_var_list,               &
                          !   spp_prt_list=spp_prt_list,                                     &
                          !   spp_stddev_cutoff=spp_stddev_cutoff,                           &
                              ! Extended diagnostics
                              ext_diag=.TRUE., &
                              tten3=tten3, &
                              qvten3=qvten3, qrten3=qrten3, qsten3=qsten3, qgten3=qgten3, &
                              qiten3=qiten3, qcten3=qcten3, pfils=PFI_LS, pflls=PFL_LS, &
                              ! Unfilled
                              ! vts1=vts1, txri=txri, txrc=txrc,                             &
                              prw_vcdc=prw_vcdc,                                             &
                              prw_vcde=prw_vcde, tpri_inu=tpri_inu, tpri_ide_d=tpri_ide_d,   &
                              tpri_ide_s=tpri_ide_s, tprs_ide=tprs_ide,                      &
                              tprs_sde_d=tprs_sde_d,                                         &
                              tprs_sde_s=tprs_sde_s, tprg_gde_d=tprg_gde_d,                  &
                              tprg_gde_s=tprg_gde_s, tpri_iha=tpri_iha,                      &
                              tpri_wfz=tpri_wfz, tpri_rfz=tpri_rfz, tprg_rfz=tprg_rfz,       &
                              tprs_scw=tprs_scw, tprg_scw=tprg_scw, tprg_rcs=tprg_rcs,       &
                              tprs_rcs=tprs_rcs,                                             &
                              tprr_rci=tprr_rci, tprg_rcg=tprg_rcg, tprw_vcd_c=tprw_vcd_c,   &
                              tprw_vcd_e=tprw_vcd_e, tprr_sml=tprr_sml, tprr_gml=tprr_gml,   &
                              tprr_rcg=tprr_rcg, tprr_rcs=tprr_rcs,                          &
                              tprv_rev=tprv_rev,                                             &
                              niten3=niten3, nrten3=nrten3, ncten3=ncten3,    &
                              ! Error handling
                              errmsg=errmsg, errflg=STATUS)
          _ASSERT( STATUS==0, errmsg )
     ! RESHAPE
         RAD_QV = RESHAPE(qv/(1.0+qv),(/IM,JM,LM/))
         RAD_QL = RESHAPE(qc,(/IM,JM,LM/))
         RAD_QR = RESHAPE(qr,(/IM,JM,LM/))
         RAD_QI = RESHAPE(qi,(/IM,JM,LM/))
         RAD_QS = RESHAPE(qs,(/IM,JM,LM/))
         RAD_QG = RESHAPE(qg,(/IM,JM,LM/))
         T      = RESHAPE(tt,(/IM,JM,LM/))
         NCPL   = RESHAPE(nc,(/IM,JM,LM/))*AIRDEN
         NCPI   = RESHAPE(ni,(/IM,JM,LM/))*AIRDEN
         NRAIN  = RESHAPE(nr,(/IM,JM,LM/))*AIRDEN
     ! Redistribute CN/LS CF/QL/QI
         call REDISTRIBUTE_CLOUDS(RAD_CF, RAD_QL, RAD_QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, T)
     ! Convert precip diagnostics from mm to kg m-2 s-1
         PRCP_RAIN    = MAX(PRCP_RAIN    / DT_MOIST, 0.0)
         PRCP_SNOW    = MAX(PRCP_SNOW    / DT_MOIST, 0.0)
         PRCP_ICE     = MAX(PRCP_ICE     / DT_MOIST, 0.0)
         PRCP_GRAUPEL = MAX(PRCP_GRAUPEL / DT_MOIST, 0.0)
     ! Fill GEOS precip diagnostics
         LS_PRCP = PRCP_RAIN
         LS_SNR  = PRCP_SNOW
         ICE     = PRCP_ICE + PRCP_GRAUPEL
         FRZR    = 0.0
     ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
         PFL_LS = PFL_LS/DT_MOIST
         PFI_LS = PFI_LS/DT_MOIST
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
                     Q(I,J,L), QLLS(I,J,L), QILS(I,J,L), QLCN(I,J,L), QICN(I,J,L), QRAIN(I,J,L), QSNOW(I,J,L), QGRAUPEL(I,J,L), NACTL(I,J,L), NACTI(I,J,L), &
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
         where (RAD_QI .le. 0.0)
            CLDREFFI = MAPL_UNDEF
         end where
         where (RAD_QL .le. 0.0)
            CLDREFFL = MAPL_UNDEF
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
        call MAPL_GetPointer(EXPORT, DBZ_MAX , 'DBZ_MAX' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_1KM , 'DBZ_1KM' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_TOP , 'DBZ_TOP' , RC=STATUS); VERIFY_(STATUS)
        call MAPL_GetPointer(EXPORT, DBZ_M10C, 'DBZ_M10C', RC=STATUS); VERIFY_(STATUS)

        if (associated(DBZ3D) .OR. &
            associated(DBZ_MAX) .OR. associated(DBZ_1KM) .OR. associated(DBZ_TOP) .OR. associated(DBZ_M10C)) then

            call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,1)
            if (associated(DBZ3D)) DBZ3D = TMP3D

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

     call MAPL_TimerOff(MAPL,"--THOM_1M",RC=STATUS)

end subroutine THOM_1M_Run

end module GEOS_THOM_1M_InterfaceMod
