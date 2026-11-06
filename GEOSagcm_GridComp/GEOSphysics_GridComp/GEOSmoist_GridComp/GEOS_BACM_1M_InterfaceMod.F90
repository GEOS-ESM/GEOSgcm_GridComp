! $Id$

#include "MAPL.h"

!=============================================================================
!BOP

! !MODULE: GEOS_BACM_1M_InterfaceMod -- A Module to interface with the
!   BACM_1M cloud microphysics

module GEOS_BACM_1M_InterfaceMod

  use ESMF
  use MAPL
  use gFTL2_StringVector, only: StringVector
  use GEOS_UtilsMod
  use GEOSmoist_Process_Library
  use CLOUDNEW, only: CLDPARAMS, PROGNO_CLOUD
  use aer_cloud
  use Aer_Actv_Single_Moment

  implicit none

  private

  integer                           :: STATUS

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

subroutine BACM_1M_Setup (GC, ADVECTED, DIFFUSED, CHEM_UPDATED, ANALYZED, RC)
    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    type(StringVector),  intent(inout) :: ADVECTED     ! internals to subscribe to advection (TRADV)
    type(StringVector),  intent(inout) :: DIFFUSED     ! internals to subscribe to turbulent diffusion (TRDIFF)
    type(StringVector),  intent(inout) :: CHEM_UPDATED ! internals to subscribe to chemistry (TRCHEM)
    type(StringVector),  intent(inout) :: ANALYZED     ! internals to subscribe to analysis (TRANA)
    integer, optional                  :: RC  ! return code

    ! !INTERNAL STATE:

    ! MAPL2 FRIENDLYTO: Q "DYNAMICS:TURBULENCE:CHEMISTRY:ANALYSIS", CLLS and CLCN "DYNAMICS",
    ! QLLS, QLCN, QILS and QICN "DYNAMICS:TURBULENCE". BACM_1M has no QRAIN, QSNOW or QGRAUPEL.

      ! Updated by chemistry, i.e. PChem's water production/loss (MAPL2 FRIENDLYTO=CHEMISTRY)
      call CHEM_UPDATED%push_back('Q')

      ! Updated by the AGCM's analysis increments (MAPL2 FRIENDLYTO=ANALYSIS)
      call ANALYZED%push_back('Q')

      ! Diffused by turbulence (MAPL2 FRIENDLYTO=TURBULENCE)
      call DIFFUSED%push_back('Q')
      call DIFFUSED%push_back('QLLS')
      call DIFFUSED%push_back('QLCN')
      call DIFFUSED%push_back('QILS')
      call DIFFUSED%push_back('QICN')

      ! Advected by dynamics (MAPL2 FRIENDLYTO=DYNAMICS)
      call ADVECTED%push_back('Q')
      call ADVECTED%push_back('CLLS')
      call ADVECTED%push_back('CLCN')
      call ADVECTED%push_back('QLLS')
      call ADVECTED%push_back('QLCN')
      call ADVECTED%push_back('QILS')
      call ADVECTED%push_back('QICN')

    !BOS

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'Q',                                        &
         STANDARD_NAME    = 'specific_humidity',                        &
         UNITS            = 'kg kg-1',                                  &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         RESTART_MODE     = MAPL_RESTART_REQUIRED,                      &
         FILL_VALUE       = 1.0e-6,                                     &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'QLLS',                                     &
         STANDARD_NAME    = 'mass_fraction_of_large_scale_cloud_liquid_water', &
         UNITS            = 'kg kg-1',                                  &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'QLCN',                                     &
         STANDARD_NAME    = 'mass_fraction_of_convective_cloud_liquid_water', &
         UNITS            = 'kg kg-1',                                  &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'CLLS',                                     &
         STANDARD_NAME    = 'large_scale_cloud_area_fraction',          &
         UNITS            = '1',                                        &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'CLCN',                                     &
         STANDARD_NAME    = 'convective_cloud_area_fraction',           &
         UNITS            = '1',                                        &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'QILS',                                     &
         STANDARD_NAME    = 'mass_fraction_of_large_scale_cloud_ice_water', &
         UNITS            = 'kg kg-1',                                  &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'QICN',                                     &
         STANDARD_NAME    = 'mass_fraction_of_convective_cloud_ice_water', &
         UNITS            = 'kg kg-1',                                  &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'NACTL',                                    &
         STANDARD_NAME    = 'activ aero # conc liq phase for 1-mom',    &
         UNITS            = 'm-3',                                      &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         RESTART_MODE     = MAPL_RESTART_SKIP,                          &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    call MAPL_GridCompAddSpec(GC,                                       &
         STATE_INTENT     = ESMF_STATEINTENT_INTERNAL,                  &
         SHORT_NAME       = 'NACTI',                                    &
         STANDARD_NAME    = 'activ aero # conc ice phase for 1-mom',    &
         UNITS            = 'm-3',                                      &
         DIMS             = 'xyz',                                      &
         VERTICAL_STAGGER = MAPL_VERTICAL_STAGGER_CENTER,               &
         RESTART_MODE     = MAPL_RESTART_SKIP,                          &
         FILL_VALUE       = 0.0,                                        &
         _RC)

    _RETURN(_SUCCESS)
end subroutine BACM_1M_Setup

subroutine BACM_1M_Initialize (GC, CLOCK, IMPORT, EXPORT, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! gridded component
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    type(ESMF_State),    intent(inout) :: IMPORT
    type(ESMF_State),    intent(inout) :: EXPORT
    integer, optional                  :: RC     ! return code

    character(len=:), allocatable :: GRIDNAME
    character(len=4)           :: imchar
    character(len=2)           :: dateline
    integer                    :: nn
    integer                    :: LM
    real                       :: tmprhL, tmprhO, TMP_ICEFALL

    call MAPL_GridCompGet(GC, num_levels=LM, _RC)

    call MAPL_GridCompGetResource(GC, 'CCW_EVAP_EFF',   CLDPARAMS%CCW_EVAP_EFF,   default= 4.0e-3, _RC)
    call MAPL_GridCompGetResource(GC, 'CCI_EVAP_EFF',   CLDPARAMS%CCI_EVAP_EFF,   default= 4.0e-3, _RC)
    call MAPL_GridCompGetResource(GC, 'HYSTPDFOPT',     CLDPARAMS%HYSTPDFOPT,     default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'PDFSHAPE',       CLDPARAMS%PDFSHAPE,       default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'CNV_BETA',       CLDPARAMS%CNV_BETA,       default= 10.0  , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_BETA',       CLDPARAMS%ANV_BETA,       default= 4.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'LS_BETA',        CLDPARAMS%LS_BETA,        default= 4.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'RH_CRIT',        CLDPARAMS%RH_CRIT,        default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'QC_CRIT_LS',     CLDPARAMS%QC_CRIT_LS,     default= 8.0e-4, _RC)
    call MAPL_GridCompGetResource(GC, 'ACCRETION',      CLDPARAMS%ACCRETION,      default= 2.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'RAIN_REVAP_FAC', CLDPARAMS%RAIN_REVAP_FAC, default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'SNOW_REVAP_FAC', CLDPARAMS%SNOW_REVAP_FAC, default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'AUTOC_LS',       CLDPARAMS%AUTOC_LS,       default= 1.0e-3, _RC)
    call MAPL_GridCompGetResource(GC, 'AUTOC_ANV',      CLDPARAMS%AUTOC_ANV,      default= 1.0e-3, _RC)
    call MAPL_GridCompGetResource(GC, 'LS_SUND_INTER',  CLDPARAMS%LS_SUND_INTER,  default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'LS_SUND_COLD',   CLDPARAMS%LS_SUND_COLD,   default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'LS_SUND_TEMP1',  CLDPARAMS%LS_SUND_TEMP1,  default= 230.  , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_SUND_INTER', CLDPARAMS%ANV_SUND_INTER, default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_SUND_COLD',  CLDPARAMS%ANV_SUND_COLD,  default= 1.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_SUND_TEMP1', CLDPARAMS%ANV_SUND_TEMP1, default= 230.  , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_TO_LS_TIME', CLDPARAMS%ANV_TO_LS_TIME, default= 14400., _RC)
    call MAPL_GridCompGetResource(GC, 'DISABLE_RAD',    CLDPARAMS%DISABLE_RAD,    default= 0.    , _RC)
    call MAPL_GridCompGetResource(GC, 'REVAP_OFF_P',    CLDPARAMS%REVAP_OFF_P,    default= 2000. , _RC)
    call MAPL_GridCompGetResource(GC, 'ICE_RAMP',       CLDPARAMS%ICE_RAMP,       default= -27.0 , _RC)
    call MAPL_GridCompGetResource(GC, 'CNV_DDRF',       CLDPARAMS%CNV_DDRF,       default= 0.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_DDRF',       CLDPARAMS%ANV_DDRF,       default= 0.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'LS_DDRF',        CLDPARAMS%LS_DDRF,        default= 0.0   , _RC)
    call MAPL_GridCompGetResource(GC, 'QC_CRIT_ANV',    CLDPARAMS%QC_CRIT_ANV,    default= 8.0e-4, _RC)
    call MAPL_GridCompGetResource(GC, 'ICE_SETTLE',     CLDPARAMS%ICE_SETTLE,     default= 1.    , _RC)
    SELECT CASE ( LM )
    CASE ( 72 )
           TMP_ICEFALL = 1.0
    CASE ( 91 )
           TMP_ICEFALL = 0.5
    CASE ( 137 )
           TMP_ICEFALL = 0.33
    CASE ( 181 )
           TMP_ICEFALL = 0.25
    CASE DEFAULT
           TMP_ICEFALL = 1.0
    END SELECT
    if (CLDPARAMS%ICE_SETTLE==2) TMP_ICEFALL=TMP_ICEFALL*0.8
    call MAPL_GridCompGetResource(GC, 'ANV_ICEFALL',     CLDPARAMS%ANV_ICEFALL, default= TMP_ICEFALL, _RC)
    call MAPL_GridCompGetResource(GC, 'LS_ICEFALL',      CLDPARAMS%LS_ICEFALL,  default= TMP_ICEFALL, _RC)
    call MAPL_GridCompGetResource(GC, 'LIQ_RADII_PARAM', LIQ_RADII_PARAM,       default= 1          , _RC)
    call MAPL_GridCompGetResource(GC, 'ICE_RADII_PARAM', ICE_RADII_PARAM,       default= 1          , _RC)
    call MAPL_GridCompGetResource(GC, 'FAC_RI',          CLDPARAMS%FAC_RI,      default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'MIN_RI',          CLDPARAMS%MIN_RI,      default=  15.e-6    , _RC)
    call MAPL_GridCompGetResource(GC, 'MAX_RI',          CLDPARAMS%MAX_RI,      default= 150.e-6    , _RC)
    call MAPL_GridCompGetResource(GC, 'FAC_RL',          CLDPARAMS%FAC_RL,      default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'MIN_RL',          CLDPARAMS%MIN_RL,      default=  5.e-6     , _RC)
    call MAPL_GridCompGetResource(GC, 'MAX_RL',          CLDPARAMS%MAX_RL,      default= 21.e-6     , _RC)
    call MAPL_GridCompGetResource(GC, 'CNV_ENVF',        CLDPARAMS%CNV_ENVF,    default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'ANV_ENVF',        CLDPARAMS%ANV_ENVF,    default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'SC_ENVF',         CLDPARAMS%SC_ENVF,     default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'LS_ENVF',         CLDPARAMS%LS_ENVF,     default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'FR_LS_WAT',       CLDPARAMS%FR_LS_WAT,   default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'FR_AN_WAT',       CLDPARAMS%FR_AN_WAT,   default= 1.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'FR_LS_ICE',       CLDPARAMS%FR_LS_ICE,   default= 0.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'FR_AN_ICE',       CLDPARAMS%FR_AN_ICE,   default= 0.0        , _RC)
    call MAPL_GridCompGetResource(GC, 'CFPBL_EXP',       CFPBL_EXP,             default= 1          , _RC)

    call MAPL_GridCompGetResource(GC, 'DBZ_LIQUID_SKIN', DBZ_LIQUID_SKIN, default= 0, _RC)

    ! Grid name to determine IMSIZE, as GWD (GEOS_GwdGridComp.F90)
    call MAPL_GridCompGetResource(GC, 'AGCM.GRIDNAME', GRIDNAME, _RC)
    GRIDNAME =  AdjustL(GRIDNAME)
    nn = len_trim(GRIDNAME)
    dateline = GRIDNAME(nn-1:nn)
    imchar = GRIDNAME(3:index(GRIDNAME,'x')-1)
    read(imchar,*) imsize
    if(dateline.eq.'CF') imsize = imsize*4

    tmprhL = CEILING(100.0*(1.0-min(0.20, max(0.01, 0.1*SQRT(SQRT(((111000.0*360.0/FLOAT(imsize))**2)/1.e10))))))/100.0 ! roundup by 0.01s
    tmprhL = min(0.99,tmprhL)
    tmprhO = min(0.99,tmprhL)
    call MAPL_GridCompGetResource(GC, 'MINRHCRITLND', MINRHCRITLND, default=tmprhL, _RC)
    call MAPL_GridCompGetResource(GC, 'MINRHCRITOCN', MINRHCRITOCN, default=tmprhO, _RC)
    call MAPL_GridCompGetResource(GC, 'MAXRHCRITLND', MAXRHCRITLND, default= 1.0  , _RC)
    call MAPL_GridCompGetResource(GC, 'MAXRHCRITOCN', MAXRHCRITOCN, default= 1.0  , _RC)
    call MAPL_GridCompGetResource(GC, 'TURNRHCRIT',   TURNRHCRIT,   default= 750.0, _RC)

    call MAPL_GridCompGetResource(GC, 'CNV_FRACTION_MIN', CNV_FRACTION_MIN, default=  500.0, _RC)
    call MAPL_GridCompGetResource(GC, 'CNV_FRACTION_MAX', CNV_FRACTION_MAX, default= 1500.0, _RC)

    call MAPL_GridCompGetResource(GC, 'ICE_FRACTION_POLYNOMIAL', ICE_FRACTION_POLYNOMIAL, default=JASON_ICE_POLYNOMIAL, _RC)

    ! Self-defaults pass a copy: value is intent(inout) and must not alias default
    call MAPL_GridCompGetResource(GC, 'USE_AEROSOL_NN', USE_AEROSOL_NN, default=(USE_AEROSOL_NN), _RC)
    call MAPL_GridCompGetResource(GC, 'USE_BERGERON',   USE_BERGERON,   default=(USE_BERGERON)  , _RC)

    call MAPL_GridCompGetResource(GC, 'NN_MIN_ICE', NN_MIN_ICE, default= 100.0e6, _RC)
    call MAPL_GridCompGetResource(GC, 'NN_MAX_ICE', NN_MAX_ICE, default= 500.0e6, _RC)

    if (USE_AEROSOL_NN) then
      ! NOTE: For now we hard code in .false. for use_wnet as that is only an option with MG and will be handled there
      call aer_cloud_init(use_wnet = .false.)
      call WRITE_PARALLEL ("INITIALIZED aer_cloud_init")
    endif

    _RETURN(_SUCCESS)
end subroutine BACM_1M_Initialize


subroutine BACM_1M_Run (GC, IMPORT, EXPORT, CLOCK, RC)
    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code:

    ! Local derived type aliases

    type (ESMF_Grid    )            :: GRID
    type (ESMF_State   )            :: INTERNAL
    real                            :: DT_MOIST

    ! Local variables

    integer                         :: IM,JM,LM
    integer                         :: I, J, L
    real, allocatable, dimension(:,:) :: LONS
    real, allocatable, dimension(:,:) :: LATS
    real, pointer, contiguous, dimension(:,:,:) :: p3d ! scratch for the edge-field bounds remap

    ! Internals
    real, pointer, dimension(:,:,:) :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN
    real, pointer, dimension(:,:,:) :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:) :: ZLE, PLE, T, U, V, KH
    real, pointer, dimension(:,:)   :: FRLAND, TS, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:) :: SL2, SL3, QT2, QT3, W2, W3, SLQT, WQT, WSL, PDF_A
    real, pointer, dimension(:,:,:) :: OMEGA
    ! Internals re-exported (were ALLOC'd exports shared with the driver, GF and UW)
    real, pointer, dimension(:,:,:) :: CNV_MFD, CNV_DQCDT, CNV_PRC3, CNV_UPDF
    real, pointer, dimension(:,:,:) :: MFD_SC, QLDET_SC, QIDET_SC, SHLW_PRC3, SHLW_SNO3, CUFRC_SC
    real, pointer, dimension(:,:,:) :: UMST0, VMST0
    real, pointer, dimension(:,:,:) :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro
    real, pointer, dimension(:,:,:) :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro
    real, pointer, dimension(:,:,:) ::  DUDT_micro,  DVDT_micro,  DTDT_micro
    real, pointer, dimension(:,:,:) :: CLDREFFL, CLDREFFI
    real, pointer, dimension(:,:,:) :: REV_CN, RSU_CN
    real, pointer, dimension(:,:,:) :: PFL_CN, PFL_AN, PFL_LS
    real, pointer, dimension(:,:,:) :: PFI_CN, PFI_AN, PFI_LS
    real, pointer, dimension(:,:  ) :: CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:  ) :: LS_PRCP, CN_PRCP, AN_PRCP, SC_PRCP
    real, pointer, dimension(:,:  ) :: LS_SNR,  CN_SNR,  AN_SNR,  SC_SNR
    ! Local
    real, allocatable, dimension(:,:,:) :: U0, V0
    real, allocatable, dimension(:,:,:) :: PLEmb, ZLE0
    real, allocatable, dimension(:,:,:) :: PLmb,  PK,  ZL0
    real, allocatable, dimension(:,:,:) :: DZET,  TH,  MASS
    real, allocatable, dimension(:,:,:) :: QDDF3, DQST3, QST3
    real, allocatable, dimension(:,:,:) :: TMP3D
    real, allocatable, dimension(:,:)   :: turnrhcrit2D
    real, allocatable, dimension(:,:)   :: minrhcrit2D
    real, allocatable, dimension(:,:)   :: maxrhcrit2D
    real, allocatable, dimension(:,:)   :: IKEX, IKEX2
    real, allocatable, dimension(:,:)   :: TMP2D
    ! Exports computed in local arrays (were ALLOC'd), copied to the export if it is connected
    real, allocatable, dimension(:,:,:) :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, allocatable, dimension(:,:,:) :: RHX
    real, allocatable, dimension(:,:,:) :: REV_SC, REV_AN, REV_LS
    real, allocatable, dimension(:,:,:) :: RSU_SC, RSU_AN, RSU_LS
    real, allocatable, dimension(:,:,:) :: ACLL_CN, ACLL_SC, ACLL_AN, ACLL_LS
    real, allocatable, dimension(:,:,:) :: ACIL_CN, ACIL_SC, ACIL_AN, ACIL_LS
    real, allocatable, dimension(:,:,:) :: PFL_SC, PFI_SC
    real, allocatable, dimension(:,:,:) :: DLPDF, DIPDF, DLFIX, DIFIX
    real, allocatable, dimension(:,:,:) :: AUT, EVAPC, SUBLC, SDM
    real, allocatable, dimension(:,:,:) :: VFALLICE_AN, VFALLICE_LS
    real, allocatable, dimension(:,:,:) :: VFALLWAT_AN, VFALLWAT_LS
    real, allocatable, dimension(:,:,:) :: VFALLRN_AN, VFALLRN_LS, VFALLRN_CN, VFALLRN_SC
    real, allocatable, dimension(:,:,:) :: VFALLSN_AN, VFALLSN_LS, VFALLSN_CN, VFALLSN_SC
    real, allocatable, dimension(:,:,:) :: FRZ_TT, FRZ_PP
    real, allocatable, dimension(:,:,:) :: DCNVL, DCNVI
    real, allocatable, dimension(:,:,:) :: ALPHT, ALPH1, ALPH2
    real, allocatable, dimension(:,:,:) :: CFPDF, RHCLR
    real, allocatable, dimension(:,:,:) :: DQRL, WTHV2, WQL
    real, allocatable, dimension(:,:,:) :: PDFITERS
    real, allocatable, dimension(:,:  ) :: LS_ARF,  CN_ARF,  AN_ARF,  SC_ARF
    real, allocatable, dimension(:,:  ) :: DTSX
    real, pointer, dimension(:,:,:) :: PTR3D
    real, pointer, dimension(:,:  ) :: PTR2D

    call MAPL_GridCompTimerStart(GC, "--BACM_1M", _RC)

    ! Get parameters from generic state.
    !-----------------------------------

    call MAPL_GridCompGet(GC, num_levels=LM, grid=GRID, _RC)
    call MAPL_GridGet(GRID, IM=IM, JM=JM, _RC)
    call MAPL_GridGetCoordinates(GRID, longitudes=LONS, latitudes=LATS, _RC) ! radians
    call MAPL_GridCompGetInternalState(GC, INTERNAL, _RC)

    ! Moist's timestep (MAPL2: the RUNALARM ring interval)
    call MAPL_ClockGet(CLOCK, dt=DT_MOIST, _RC)

    ! Internal State
    call MAPL_StateGetPointer(INTERNAL, Q,        'Q'       , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLLS,     'QLLS'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLCN,     'QLCN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, CLCN,     'CLCN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, CLLS,     'CLLS'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QILS,     'QILS'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QICN,     'QICN'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, NACTL,    'NACTL'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, NACTI,    'NACTI'   , _RC)

    ! Import State
    call MAPL_StateGetPointer(IMPORT, ZLE,     'ZLE'     , _RC)
    call MAPL_StateGetPointer(IMPORT, PLE,     'PLE'     , _RC)
    call MAPL_StateGetPointer(IMPORT, T,       'T'       , _RC)
    call MAPL_StateGetPointer(IMPORT, U,       'U'       , _RC)
    call MAPL_StateGetPointer(IMPORT, V,       'V'       , _RC)
    call MAPL_StateGetPointer(IMPORT, FRLAND,  'FRLAND'  , _RC)
    call MAPL_StateGetPointer(IMPORT, KH,      'KH'      , _RC)
    call MAPL_StateGetPointer(IMPORT, PDF_A,   'PDF_A'   , _RC)
    call MAPL_StateGetPointer(IMPORT, W2,      'W2'      , _RC)
    call MAPL_StateGetPointer(IMPORT, W3,      'W3'      , _RC)
    call MAPL_StateGetPointer(IMPORT, WQT,     'WQT'     , _RC)
    call MAPL_StateGetPointer(IMPORT, WSL,     'WSL'     , _RC)
    call MAPL_StateGetPointer(IMPORT, SL2,     'SL2'     , _RC)
    call MAPL_StateGetPointer(IMPORT, SL3,     'SL3'     , _RC)
    call MAPL_StateGetPointer(IMPORT, QT2,     'QT2'     , _RC)
    call MAPL_StateGetPointer(IMPORT, QT3,     'QT3'     , _RC)
    call MAPL_StateGetPointer(IMPORT, SLQT,    'SLQT'    , _RC)
    call MAPL_StateGetPointer(IMPORT, TS,      'TS'      , _RC)
    call MAPL_StateGetPointer(IMPORT, KPBLSC,  'KPBL_SC' , _RC)
    call MAPL_StateGetPointer(IMPORT, SH,      'SH'      , _RC)
    call MAPL_StateGetPointer(IMPORT, EVAP,    'EVAP'    , _RC)
    call MAPL_StateGetPointer(IMPORT, OMEGA,   'OMEGA'   , _RC)
    ! MAPL2 allocated edge fields 0:LM; MAPL3 pointers are 1-based
    p3d => ZLE
    ZLE(1:IM, 1:JM, 0:LM) => p3d
    p3d => PLE
    PLE(1:IM, 1:JM, 0:LM) => p3d
    p3d => KH
    KH(1:IM, 1:JM, 0:LM) => p3d

    ! Allocatables
     ! Edge variables
    ALLOCATE ( ZLE0 (IM,JM,0:LM) )
    ALLOCATE ( PLEmb(IM,JM,0:LM) )
     ! Layer variables
    ALLOCATE ( U0   (IM,JM,LM  ) )
    ALLOCATE ( V0   (IM,JM,LM  ) )
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
    ALLOCATE ( IKEX         (IM,JM) )
    ALLOCATE ( IKEX2        (IM,JM) )
    ALLOCATE ( TMP2D        (IM,JM) )

    ! Derived States
    PLEmb    =  PLE*.01
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
    U0       = U
    V0       = V

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

    ! Tendencies (internals re-exported, mapl3-alloc-exports.md group A)
    call MAPL_StateGetPointer(INTERNAL, DQVDT_macro, 'DQVDT_macro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQIDT_macro, 'DQIDT_macro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQLDT_macro, 'DQLDT_macro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQADT_macro, 'DQADT_macro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQVDT_micro, 'DQVDT_micro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQIDT_micro, 'DQIDT_micro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQLDT_micro, 'DQLDT_micro' , _RC)
    call MAPL_StateGetPointer(INTERNAL, DQADT_micro, 'DQADT_micro' , _RC)
    call MAPL_StateGetPointer(INTERNAL,  DUDT_micro,  'DUDT_micro' , _RC)
    call MAPL_StateGetPointer(INTERNAL,  DVDT_micro,  'DVDT_micro' , _RC)
    call MAPL_StateGetPointer(INTERNAL,  DTDT_micro,  'DTDT_micro' , _RC)
    DQVDT_micro = Q
    DQLDT_micro = QLLS + QLCN
    DQIDT_micro = QILS + QICN
    DQADT_micro = CLLS + CLCN
     DUDT_micro = U
     DVDT_micro = V
     DTDT_micro = T

    ! Imports which are Exports from local siblings
    ! Required export MUST have been filled in GridComp
    call MAPL_StateGetPointer(INTERNAL, CNV_FRC,  'CNV_FRC' , _RC)
    call MAPL_StateGetPointer(INTERNAL, SRF_TYPE, 'SRF_TYPE', _RC)
    ! DeepCu : default to 0.0 in MAPL if not running DeepCu
    call MAPL_StateGetPointer(INTERNAL, CNV_MFD,    'CNV_MFD'   , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNV_DQCDT,  'CNV_DQCDT' , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNV_PRC3,   'CNV_PRC3'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, CNV_UPDF,   'CNV_UPDF'  , _RC)
    ! ShallowCu : default to 0.0 in MAPL if not running ShallowCu
    call MAPL_StateGetPointer(INTERNAL, MFD_SC,     'MFD_SC'    , _RC)
    call MAPL_StateGetPointer(INTERNAL, QLDET_SC,   'QLDET_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, QIDET_SC,   'QIDET_SC'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, SHLW_PRC3,  'SHLW_PRC3' , _RC)
    call MAPL_StateGetPointer(INTERNAL, SHLW_SNO3,  'SHLW_SNO3' , _RC)
    call MAPL_StateGetPointer(INTERNAL, CUFRC_SC,   'CUFRC_SC'  , _RC)
    ! Export and/or scratch Variable
    call MAPL_StateGetPointer(INTERNAL, CLDREFFL, 'RL'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, CLDREFFI, 'RI'  , _RC)
    ! Exports passed to be filled in progno_cloud
    call MAPL_StateGetPointer(INTERNAL, LS_PRCP,  'LS_PRCP' , _RC)
    call MAPL_StateGetPointer(INTERNAL, CN_PRCP,  'CN_PRCP' , _RC)
    call MAPL_StateGetPointer(INTERNAL, AN_PRCP,  'AN_PRCP' , _RC)
    call MAPL_StateGetPointer(INTERNAL, SC_PRCP,  'SC_PRCP' , _RC)
    call MAPL_StateGetPointer(INTERNAL, LS_SNR,   'LS_SNR'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, CN_SNR,   'CN_SNR'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, AN_SNR,   'AN_SNR'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, SC_SNR,   'SC_SNR'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, REV_CN,   'REV_CN'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, RSU_CN,   'RSU_CN'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFL_CN,   'PFL_CN'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFL_AN,   'PFL_AN'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFL_LS,   'PFL_LS'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFI_CN,   'PFI_CN'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFI_AN,   'PFI_AN'  , _RC)
    call MAPL_StateGetPointer(INTERNAL, PFI_LS,   'PFI_LS'  , _RC)
    ! Local to this routine (mapl3-alloc-exports.md group C), zero-filled as MAPL2's ALLOC did
    allocate(RAD_CF(IM,JM,LM), RAD_QV(IM,JM,LM), RAD_QL(IM,JM,LM), RAD_QI(IM,JM,LM), source=0.0)
    allocate(RAD_QR(IM,JM,LM), RAD_QS(IM,JM,LM), RAD_QG(IM,JM,LM), RHX(IM,JM,LM), source=0.0)
    allocate(REV_SC(IM,JM,LM), REV_AN(IM,JM,LM), REV_LS(IM,JM,LM), source=0.0)
    allocate(RSU_SC(IM,JM,LM), RSU_AN(IM,JM,LM), RSU_LS(IM,JM,LM), source=0.0)
    allocate(ACLL_CN(IM,JM,LM), ACLL_SC(IM,JM,LM), ACLL_AN(IM,JM,LM), ACLL_LS(IM,JM,LM), source=0.0)
    allocate(ACIL_CN(IM,JM,LM), ACIL_SC(IM,JM,LM), ACIL_AN(IM,JM,LM), ACIL_LS(IM,JM,LM), source=0.0)
    allocate(PFL_SC(IM,JM,0:LM), PFI_SC(IM,JM,0:LM), source=0.0) ! edge
    allocate(DLPDF(IM,JM,LM), DIPDF(IM,JM,LM), DLFIX(IM,JM,LM), DIFIX(IM,JM,LM), source=0.0)
    allocate(AUT(IM,JM,LM), EVAPC(IM,JM,LM), SDM(IM,JM,LM), SUBLC(IM,JM,LM), source=0.0)
    allocate(VFALLICE_AN(IM,JM,LM), VFALLICE_LS(IM,JM,LM), VFALLWAT_AN(IM,JM,LM), VFALLWAT_LS(IM,JM,LM), source=0.0)
    allocate(VFALLRN_AN(IM,JM,LM), VFALLRN_LS(IM,JM,LM), VFALLRN_CN(IM,JM,LM), VFALLRN_SC(IM,JM,LM), source=0.0)
    allocate(VFALLSN_AN(IM,JM,LM), VFALLSN_LS(IM,JM,LM), VFALLSN_CN(IM,JM,LM), VFALLSN_SC(IM,JM,LM), source=0.0)
    allocate(FRZ_TT(IM,JM,LM), FRZ_PP(IM,JM,LM), DCNVL(IM,JM,LM), DCNVI(IM,JM,LM), source=0.0)
    allocate(ALPHT(IM,JM,LM), ALPH1(IM,JM,LM), ALPH2(IM,JM,LM), CFPDF(IM,JM,LM), RHCLR(IM,JM,LM), source=0.0)
    allocate(DQRL(IM,JM,LM), PDFITERS(IM,JM,LM), WTHV2(IM,JM,LM), WQL(IM,JM,LM), source=0.0)
    allocate(LS_ARF(IM,JM), CN_ARF(IM,JM), AN_ARF(IM,JM), SC_ARF(IM,JM), DTSX(IM,JM), source=0.0)

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
              WQT               , &   ! <- turb
              WSL               , &   ! <- turb
              QT2               , &   ! <- turb
              SL2               , &   ! <- turb
              SLQT              , &   ! <- turb
              W2                , &   ! <- turb
              W3                , &   ! <- turb
              QT3               , &   ! <- turb
              SL3               , &   ! <- turb
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
              SRF_TYPE          , &
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
              PDF_A, PDFITERS, &
              DQVDT_macro, DQLDT_macro, DQIDT_macro, DQADT_macro, &
              WTHV2, WQL, &
              NACTL,    &
              NACTI,    &
              "GF" )
         ! update T
         T = TH*PK
         ! Exports
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
         call FILLQ2ZERO(RAD_QV, MASS, RC=STATUS); _VERIFY(STATUS)
         call FILLQ2ZERO(RAD_QL, MASS, RC=STATUS); _VERIFY(STATUS)
         call FILLQ2ZERO(RAD_QI, MASS, RC=STATUS); _VERIFY(STATUS)
         call FILLQ2ZERO(RAD_QR, MASS, RC=STATUS); _VERIFY(STATUS)
         call FILLQ2ZERO(RAD_QS, MASS, RC=STATUS); _VERIFY(STATUS)
         call FILLQ2ZERO(RAD_QG, MASS, RC=STATUS); _VERIFY(STATUS)
         call FILLQ2ZERO(RAD_CF, MASS, RC=STATUS); _VERIFY(STATUS)
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

         DQVDT_micro = ( Q          - DQVDT_micro) / DT_MOIST - DQVDT_macro
         DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST - DQLDT_macro
         DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST - DQIDT_macro
         DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST - DQADT_macro
          DUDT_micro = ( U          -  DUDT_micro) / DT_MOIST
          DVDT_micro = ( V          -  DVDT_micro) / DT_MOIST
          DTDT_micro = ( T          -  DTDT_micro) / DT_MOIST

#ifdef KEDISS_MICRO
        ! dissipative heating tendency from KE across the macro/micro physics
         call MAPL_StateGetPointer(EXPORT, PTR3D, 'DTDTFRIC', _RC)
         if(associated(PTR3D)) then
           call dissipative_ke_heating(IM,JM,LM, MASS,U0,V0, &
                                       DUDT_micro,DVDT_micro,PTR3D)
         endif
#else
        ! Cumulus Friction
        call MAPL_StateGetPointer(INTERNAL, UMST0, 'UMST0', _RC)
        call MAPL_StateGetPointer(INTERNAL, VMST0, 'VMST0', _RC)
        IKEX  = SUM( (0.5/DT_MOIST)*((V**2+U**2)  - (VMST0**2+UMST0**2))*MASS , 3 )
        IKEX2 = MAX(  SUM( 1.E-04 * MASS , 3 ) ,  1.0e-6 ) ! floor at 1e-6 W m-2
        !scaled 3D kinetic energy dissipation
        call MAPL_StateGetPointer(EXPORT, PTR3D, 'KEDISS', _RC)
        if (associated(PTR3D))  then
           do L=1,LM
              PTR3D(:,:,L) = (IKEX/IKEX2) * 1.E-04
           enddo
        end if
        call MAPL_StateGetPointer(EXPORT, PTR3D, 'DTDTFRIC', _RC)
        if(associated(PTR3D)) then
           do L=1,LM
              PTR3D(:,:,L) = -(1./MAPL_CP)*(IKEX/IKEX2) * 1.E-04 * (PLE(:,:,L)-PLE(:,:,L-1))
           end do
        end if
#endif

         ! Fill these with 0.0, diagnose in GridComp if desired
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'ICE',  _RC)
         PTR2D = 0.0
         call MAPL_StateGetPointer(INTERNAL, PTR2D, 'FRZR', _RC)
         PTR2D = 0.0

         ! Compute DBZ radar reflectivity
         call MAPL_StateGetPointer(EXPORT, PTR3D, 'DBZ'    , _RC)
         call MAPL_StateGetPointer(EXPORT, PTR2D, 'DBZ_MAX', _RC)
         if (associated(PTR3D) .OR. associated(PTR2D)) then
            call CALCDBZ(TMP3D,100*PLmb,T,Q,RAD_QR*RAD_CF,RAD_QS*RAD_CF,RAD_QG*RAD_CF,IM,JM,LM,1,0,DBZ_LIQUID_SKIN)
            if (associated(PTR3D)) PTR3D = TMP3D
            if (associated(PTR2D)) then
               PTR2D=-9999.0
               DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                  PTR2D(I,J) = MAX(PTR2D(I,J),TMP3D(I,J,L))
               END DO ; END DO ; END DO
            endif
         endif

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'QRTOT', _RC)
         if (associated(PTR3D)) PTR3D = RAD_QR*RAD_CF

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'QSTOT', _RC)
         if (associated(PTR3D)) PTR3D = RAD_QS*RAD_CF

         call MAPL_StateGetPointer(EXPORT, PTR3D, 'QGTOT', _RC)
         if (associated(PTR3D)) PTR3D = RAD_QG*RAD_CF

         call MAPL_StateGetPointer(EXPORT, PTR2D, 'LWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM( ( QLCN+QLLS ) *MASS , 3 )

         call MAPL_StateGetPointer(EXPORT, PTR2D, 'IWP', _RC)
         if (associated(PTR2D)) PTR2D = SUM( ( QICN+QILS ) *MASS , 3 )

         ! Copy the local arrays to their exports, if connected
         call copy_to_export_3d('FCLD', RAD_CF, _RC)
         call copy_to_export_3d('QV'  , RAD_QV, _RC)
         call copy_to_export_3d('QL'  , RAD_QL, _RC)
         call copy_to_export_3d('QI'  , RAD_QI, _RC)
         call copy_to_export_3d('QR'  , RAD_QR, _RC)
         call copy_to_export_3d('QS'  , RAD_QS, _RC)
         call copy_to_export_3d('QG'  , RAD_QG, _RC)
         call copy_to_export_3d('RHX'     , RHX    , _RC)
         call copy_to_export_3d('REV_SC'  , REV_SC , _RC)
         call copy_to_export_3d('REV_AN'  , REV_AN , _RC)
         call copy_to_export_3d('REV_LS'  , REV_LS , _RC)
         call copy_to_export_3d('RSU_SC'  , RSU_SC , _RC)
         call copy_to_export_3d('RSU_AN'  , RSU_AN , _RC)
         call copy_to_export_3d('RSU_LS'  , RSU_LS , _RC)
         call copy_to_export_3d('ACRLL_CN', ACLL_CN, _RC)
         call copy_to_export_3d('ACRLL_SC', ACLL_SC, _RC)
         call copy_to_export_3d('ACRLL_AN', ACLL_AN, _RC)
         call copy_to_export_3d('ACRLL_LS', ACLL_LS, _RC)
         call copy_to_export_3d('ACRIL_CN', ACIL_CN, _RC)
         call copy_to_export_3d('ACRIL_SC', ACIL_SC, _RC)
         call copy_to_export_3d('ACRIL_AN', ACIL_AN, _RC)
         call copy_to_export_3d('ACRIL_LS', ACIL_LS, _RC)
         call copy_to_export_3d('PFL_SC'  , PFL_SC , _RC)
         call copy_to_export_3d('PFI_SC'  , PFI_SC , _RC)
         call copy_to_export_3d('DLPDF'   , DLPDF  , _RC)
         call copy_to_export_3d('DIPDF'   , DIPDF  , _RC)
         call copy_to_export_3d('DLFIX'   , DLFIX  , _RC)
         call copy_to_export_3d('DIFIX'   , DIFIX  , _RC)
         call copy_to_export_3d('AUT'     , AUT    , _RC)
         call copy_to_export_3d('EVAPC'   , EVAPC  , _RC)
         call copy_to_export_3d('SDM'     , SDM    , _RC)
         call copy_to_export_3d('SUBLC'   , SUBLC  , _RC)
         call copy_to_export_3d('VFALLICE_AN', VFALLICE_AN, _RC)
         call copy_to_export_3d('VFALLICE_LS', VFALLICE_LS, _RC)
         call copy_to_export_3d('VFALLWAT_AN', VFALLWAT_AN, _RC)
         call copy_to_export_3d('VFALLWAT_LS', VFALLWAT_LS, _RC)
         call copy_to_export_3d('VFALLRN_AN' , VFALLRN_AN , _RC)
         call copy_to_export_3d('VFALLRN_LS' , VFALLRN_LS , _RC)
         call copy_to_export_3d('VFALLRN_CN' , VFALLRN_CN , _RC)
         call copy_to_export_3d('VFALLRN_SC' , VFALLRN_SC , _RC)
         call copy_to_export_3d('VFALLSN_AN' , VFALLSN_AN , _RC)
         call copy_to_export_3d('VFALLSN_LS' , VFALLSN_LS , _RC)
         call copy_to_export_3d('VFALLSN_CN' , VFALLSN_CN , _RC)
         call copy_to_export_3d('VFALLSN_SC' , VFALLSN_SC , _RC)
         call copy_to_export_3d('FRZ_TT'  , FRZ_TT , _RC)
         call copy_to_export_3d('FRZ_PP'  , FRZ_PP , _RC)
         call copy_to_export_3d('DCNVL'   , DCNVL  , _RC)
         call copy_to_export_3d('DCNVI'   , DCNVI  , _RC)
         call copy_to_export_3d('ALPHT'   , ALPHT  , _RC)
         call copy_to_export_3d('ALPH1'   , ALPH1  , _RC)
         call copy_to_export_3d('ALPH2'   , ALPH2  , _RC)
         call copy_to_export_3d('CFPDF'   , CFPDF  , _RC)
         call copy_to_export_3d('RHCLR'   , RHCLR  , _RC)
         call copy_to_export_3d('DQRL'    , DQRL   , _RC)
         call copy_to_export_3d('PDFITERS', PDFITERS, _RC)
         call copy_to_export_3d('WTHV2'   , WTHV2  , _RC)
         call copy_to_export_3d('WQL'     , WQL    , _RC)
         call copy_to_export_2d('LS_ARF'  , LS_ARF , _RC)
         call copy_to_export_2d('CN_ARF'  , CN_ARF , _RC)
         call copy_to_export_2d('AN_ARF'  , AN_ARF , _RC)
         call copy_to_export_2d('SC_ARF'  , SC_ARF , _RC)
         call copy_to_export_2d('DTSX'    , DTSX   , _RC)

    call MAPL_GridCompTimerStop(GC, "--BACM_1M", _RC)

    _RETURN(_SUCCESS)

contains

    subroutine copy_to_export_3d(name, local, rc)
       character(*), intent(in) :: name
       real, intent(in) :: local(:,:,:)
       integer, optional, intent(out) :: rc
       real, pointer :: ptr(:,:,:)
       integer :: status
       call MAPL_StateGetPointer(EXPORT, ptr, name, _RC)
       if (associated(ptr)) ptr = local
       _RETURN(_SUCCESS)
    end subroutine copy_to_export_3d

    subroutine copy_to_export_2d(name, local, rc)
       character(*), intent(in) :: name
       real, intent(in) :: local(:,:)
       integer, optional, intent(out) :: rc
       real, pointer :: ptr(:,:)
       integer :: status
       call MAPL_StateGetPointer(EXPORT, ptr, name, _RC)
       if (associated(ptr)) ptr = local
       _RETURN(_SUCCESS)
    end subroutine copy_to_export_2d

end subroutine BACM_1M_Run

end module GEOS_BACM_1M_InterfaceMod
