! $Id$
! $Name$

!#define PDFDIAG 1

module cloudnew

#ifndef _CUDA
   use GEOS_UtilsMod,     only: QSAT=>GEOS_Qsat, DQSAT=>GEOS_DQsat, &
         QSATLQ=>GEOS_QsatLQU, QSATIC=>GEOS_QsatICE
#else
   use cudafor
   ! NOTE: GPUs use the QSAT and DQSAT at the end of this module
#endif

   use MAPL_ConstantsMod, only: MAPL_TICE , MAPL_CP   , &
                                MAPL_GRAV , MAPL_ALHS , &
                                MAPL_ALHL , MAPL_ALHF , &
                                MAPL_RGAS , MAPL_H2OMW, &
                                MAPL_AIRMW, MAPL_RVAP , &
                                MAPL_PI   , MAPL_R8   , &
                                MAPL_R4   , MAPL_AVOGAD

   use MAPL_BaseMod,      only: MAPL_UNDEF
   use Aer_Actv_Single_Moment,only: USE_AEROSOL_NN
   use GEOSmoist_Process_Library

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
           real               :: SNOW_REVAP_FAC        ! 57
           real               :: CCW_EVAP_EFF          ! 13
           real               :: CCI_EVAP_EFF          ! 13
           real               :: LS_SUND_INTER         ! 15
           real               :: LS_SUND_COLD          ! 16
           real               :: LS_SUND_TEMP1         ! 17
           real               :: ANV_SUND_INTER        ! 18
           real               :: ANV_SUND_COLD         ! 19
           real               :: ANV_SUND_TEMP1        ! 20
           real               :: ANV_TO_LS_TIME        ! 21
           real               :: DISABLE_RAD           ! 26
           real               :: ICE_SETTLE            ! 27
           real               :: ANV_ICEFALL           ! 28
           real               :: LS_ICEFALL            ! 29
           real               :: REVAP_OFF_P           ! 30
           real               :: CNV_ENVF              ! 31
           real               :: ANV_ENVF              ! 31
           real               :: SC_ENVF               ! 31
           real               :: LS_ENVF               ! 31
           real               :: ICE_RAMP              ! 33
           real               :: CNV_DDRF              ! 36
           real               :: ANV_DDRF              ! 37
           real               :: LS_DDRF               ! 38
           real               :: AUTOC_ANV             ! 39
           real               :: QC_CRIT_ANV           ! 40
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
           real               :: PDFSHAPE              ! 58
      endtype CLDPARAM_TYPE
      type (CLDPARAM_TYPE) :: CLDPARAMS

#ifndef _CUDA
   private

   PUBLIC CLDPARAMS
   PUBLIC PROGNO_CLOUD
#endif

#ifdef _CUDA

   ! Inputs
   ! ------

   real, allocatable, dimension(:,:), device :: PP_dev
   real, allocatable, dimension(:,:), device :: EXNP_dev
   real, allocatable, dimension(:,:), device :: PPE_dev
   real, allocatable, dimension(:  ), device :: FRLAND_dev
   real, allocatable, dimension(:,:), device :: KH_dev
   real, allocatable, dimension(:,:), device :: RMFDTR_dev
   real, allocatable, dimension(:,:), device :: QLWDTR_dev
   real, allocatable, dimension(:,:), device :: U_dev
   real, allocatable, dimension(:,:), device :: V_dev
   real, allocatable, dimension(:,:), device :: QST3_dev
   real, allocatable, dimension(:,:), device :: DZET_dev
   real, allocatable, dimension(:,:), device :: QDDF3_dev
   real, allocatable, dimension(:  ), device :: CNVFRC_dev
   real, allocatable, dimension(:  ), device :: TROPP_dev

   ! Inoutputs
   ! ---------

   real, allocatable, dimension(:,:), device :: TH_dev
   real, allocatable, dimension(:,:), device :: Q_dev
   real, allocatable, dimension(:,:), device :: QRN_CU_dev
   real, allocatable, dimension(:,:), device :: CNV_UPDFRC_dev ! on edges, but dims=1:LM
   real, allocatable, dimension(:,:), device :: QLW_LS_dev  
   real, allocatable, dimension(:,:), device :: QLW_AN_dev
   real, allocatable, dimension(:,:), device :: QIW_LS_dev  
   real, allocatable, dimension(:,:), device :: QIW_AN_dev
   real, allocatable, dimension(:,:), device :: ANVFRC_dev
   real, allocatable, dimension(:,:), device :: CLDFRC_dev

   ! Outputs
   ! -------

   real, allocatable, dimension(:,:), device :: RAD_CLDFRC_dev
   real, allocatable, dimension(:,:), device :: RAD_QV_dev
   real, allocatable, dimension(:,:), device :: RAD_QL_dev
   real, allocatable, dimension(:,:), device :: RAD_QI_dev
   real, allocatable, dimension(:,:), device :: RAD_QR_dev
   real, allocatable, dimension(:,:), device :: RAD_QS_dev
   real, allocatable, dimension(:,:), device :: RAD_QG_dev
   real, allocatable, dimension(:,:), device :: CLDREFFL_dev
   real, allocatable, dimension(:,:), device :: CLDREFFI_dev
   real, allocatable, dimension(:  ), device :: PRELS_dev
   real, allocatable, dimension(:  ), device :: PRECU_dev
   real, allocatable, dimension(:  ), device :: PREAN_dev
   real, allocatable, dimension(:  ), device :: LSARF_dev
   real, allocatable, dimension(:  ), device :: CUARF_dev
   real, allocatable, dimension(:  ), device :: ANARF_dev
   real, allocatable, dimension(:  ), device :: SNRLS_dev
   real, allocatable, dimension(:  ), device :: SNRCU_dev
   real, allocatable, dimension(:  ), device :: SNRAN_dev

   real, allocatable, dimension(:,:), device :: PFL_CN_dev
   real, allocatable, dimension(:,:), device :: PFI_CN_dev
   real, allocatable, dimension(:,:), device :: PFL_AN_dev
   real, allocatable, dimension(:,:), device :: PFI_AN_dev
   real, allocatable, dimension(:,:), device :: PFL_LS_dev
   real, allocatable, dimension(:,:), device :: PFI_LS_dev

   real, allocatable, dimension(:,:), device :: RHX_dev
   real, allocatable, dimension(:,:), device :: REV_LS_dev
   real, allocatable, dimension(:,:), device :: REV_AN_dev
   real, allocatable, dimension(:,:), device :: REV_CN_dev
   real, allocatable, dimension(:,:), device :: RSU_LS_dev
   real, allocatable, dimension(:,:), device :: RSU_AN_dev
   real, allocatable, dimension(:,:), device :: RSU_CN_dev
   real, allocatable, dimension(:,:), device :: ACLL_CN_dev
   real, allocatable, dimension(:,:), device :: ACIL_CN_dev
   real, allocatable, dimension(:,:), device :: ACLL_AN_dev
   real, allocatable, dimension(:,:), device :: ACIL_AN_dev
   real, allocatable, dimension(:,:), device :: ACLL_LS_dev
   real, allocatable, dimension(:,:), device :: ACIL_LS_dev
   real, allocatable, dimension(:,:), device :: PDFL_dev
   real, allocatable, dimension(:,:), device :: PDFI_dev
   real, allocatable, dimension(:,:), device :: FIXL_dev
   real, allocatable, dimension(:,:), device :: FIXI_dev                          
   real, allocatable, dimension(:,:), device :: AUT_dev
   real, allocatable, dimension(:,:), device :: EVAPC_dev
   real, allocatable, dimension(:,:), device :: SDM_dev
   real, allocatable, dimension(:,:), device :: SUBLC_dev 
   real, allocatable, dimension(:,:), device :: FRZ_TT_dev
   real, allocatable, dimension(:,:), device :: FRZ_PP_dev
   real, allocatable, dimension(:,:), device :: DCNVL_dev
   real, allocatable, dimension(:,:), device :: DCNVi_dev
   real, allocatable, dimension(:,:), device :: ALPHT_dev
   real, allocatable, dimension(:,:), device :: ALPH1_dev
   real, allocatable, dimension(:,:), device :: ALPH2_dev
   real, allocatable, dimension(:,:), device :: CFPDF_dev
   real, allocatable, dimension(:,:), device :: RHCLR_dev
   real, allocatable, dimension(:,:), device :: DQRL_dev
   real, allocatable, dimension(:,:), device :: VFALLICE_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLICE_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLWAT_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLWAT_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLSN_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLSN_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLSN_CN_dev
   real, allocatable, dimension(:,:), device :: VFALLRN_AN_dev
   real, allocatable, dimension(:,:), device :: VFALLRN_LS_dev
   real, allocatable, dimension(:,:), device :: VFALLRN_CN_dev

   real, allocatable, dimension(:,:), device :: LIQANMOVE_dev
   real, allocatable, dimension(:,:), device :: ICEANMOVE_dev
   real, allocatable, dimension(:,:), device :: DANCLD_dev
   real, allocatable, dimension(:,:), device :: DLSCLD_dev
   real, allocatable, dimension(:,:), device :: CURAINMOVE_dev
   real, allocatable, dimension(:,:), device :: CUSNOWMOVE_dev

   ! Constants passed in from CLDPARAMS (from MAPL_GetResource)
   real,    constant :: CNV_BETA
   real,    constant :: ANV_BETA
   real,    constant :: LS_BETA
   real,    constant :: RH00
   real,    constant :: C_00
   real,    constant :: LWCRIT
   real,    constant :: C_ACC
   real,    constant :: C_EV_R
   real,    constant :: C_EV_S
   real,    constant :: CCW_EVP_EFF
   real,    constant :: CCI_EVP_EFF
   real,    constant :: LS_SDQV2
   real,    constant :: LS_SDQV3
   real,    constant :: LS_SDQVT1
   real,    constant :: ANV_SDQV2
   real,    constant :: ANV_SDQV3
   real,    constant :: ANV_SDQVT1
   real,    constant :: ANV_TO_LS
   integer, constant :: DISABLE_RAD
   integer, constant :: ICE_SETTLE
   real,    constant :: ANV_ICEFALL_C
   real,    constant :: LS_ICEFALL_C
   real,    constant :: REVAP_OFF_P
   real,    constant :: CNVENVFC
   real,    constant :: ANVENVFC
   real,    constant :: SCENVFC
   real,    constant :: LSENVFC
   real,    constant :: T_ICE_ALL
   real,    constant :: CNVDDRFC
   real,    constant :: ANVDDRFC
   real,    constant :: LSDDRFC
   integer, constant :: FR_LS_WAT
   integer, constant :: FR_LS_ICE
   integer, constant :: FR_AN_WAT
   integer, constant :: FR_AN_ICE
   real,    constant :: MIN_RL
   real,    constant :: MIN_RI
   real,    constant :: MAX_RL
   real,    constant :: MAX_RI
   real,    constant :: FAC_RL
   real,    constant :: FAC_RI
   integer, constant :: PDFFLAG

   ! Parameters for Internal DQSAT
   ! -----------------------------

   real, parameter :: MAX_MIXING_RATIO = 1.
   real, parameter :: ZEROC            = MAPL_TICE

   real, parameter :: TMINTBL   =  150.0
   real, parameter :: TMAXTBL   =  333.0
   real, parameter :: DEGSUBS   =  100
   real, parameter :: ERFAC     = (DEGSUBS/ESFAC)
   real, parameter :: DELTA_T   =  1.0 / DEGSUBS
   real, parameter :: TABLESIZE =  nint(TMAXTBL-TMINTBL)*DEGSUBS + 1
   real, parameter :: TMIX      = -20.

   real, parameter :: TMINSTR = -95.
   real, parameter :: TSTARR1 = -75.
   real, parameter :: TSTARR2 = -65.
   real, parameter :: TSTARR3 = -50.
   real, parameter :: TSTARR4 = -40.
   real, parameter :: TMAXSTR = +60.

   real(kind=MAPL_R8), parameter :: B6 = 6.136820929E-11*100.0
   real(kind=MAPL_R8), parameter :: B5 = 2.034080948E-8 *100.0
   real(kind=MAPL_R8), parameter :: B4 = 3.031240396E-6 *100.0
   real(kind=MAPL_R8), parameter :: B3 = 2.650648471E-4 *100.0
   real(kind=MAPL_R8), parameter :: B2 = 1.428945805E-2 *100.0
   real(kind=MAPL_R8), parameter :: B1 = 4.436518521E-1 *100.0
   real(kind=MAPL_R8), parameter :: B0 = 6.107799961E+0 *100.0
   real(kind=MAPL_R8), parameter :: BI6= 1.838826904E-10*100.0
   real(kind=MAPL_R8), parameter :: BI5= 4.838803174E-8 *100.0
   real(kind=MAPL_R8), parameter :: BI4= 5.824720280E-6 *100.0
   real(kind=MAPL_R8), parameter :: BI3= 4.176223716E-4 *100.0
   real(kind=MAPL_R8), parameter :: BI2= 1.886013408E-2 *100.0
   real(kind=MAPL_R8), parameter :: BI1= 5.034698970E-1 *100.0
   real(kind=MAPL_R8), parameter :: BI0= 6.109177956E+0 *100.0
   real(kind=MAPL_R8), parameter :: S16= 0.516000335E-11*100.0
   real(kind=MAPL_R8), parameter :: S15= 0.276961083E-8 *100.0
   real(kind=MAPL_R8), parameter :: S14= 0.623439266E-6 *100.0
   real(kind=MAPL_R8), parameter :: S13= 0.754129933E-4 *100.0
   real(kind=MAPL_R8), parameter :: S12= 0.517609116E-2 *100.0
   real(kind=MAPL_R8), parameter :: S11= 0.191372282E+0 *100.0
   real(kind=MAPL_R8), parameter :: S10= 0.298152339E+1 *100.0
   real(kind=MAPL_R8), parameter :: S26= 0.314296723E-10*100.0
   real(kind=MAPL_R8), parameter :: S25= 0.132243858E-7 *100.0
   real(kind=MAPL_R8), parameter :: S24= 0.236279781E-5 *100.0
   real(kind=MAPL_R8), parameter :: S23= 0.230325039E-3 *100.0
   real(kind=MAPL_R8), parameter :: S22= 0.129690326E-1 *100.0
   real(kind=MAPL_R8), parameter :: S21= 0.401390832E+0 *100.0
   real(kind=MAPL_R8), parameter :: S20= 0.535098336E+1 *100.0

   real(kind=MAPL_R8), parameter :: DI(0:3) = (/ 57518.5606E08, 2.01889049, 3.56654, 20.947031 /)
   real(kind=MAPL_R8), parameter :: CI(0:3) = (/ 9.550426, -5723.265, 3.53068, -.00728332 /)
   real(kind=MAPL_R8), parameter :: DL(1:6) = (/ -7.902980, 5.02808, -1.3816, 11.344, 8.1328, -3.49149 /)
   real(kind=MAPL_R8), parameter :: LOGPS   = 3.005714898  ! log10(1013.246)
   real(kind=MAPL_R8), parameter :: TS      = 373.16
   real(kind=MAPL_R8), parameter :: CL(0:9) = (/54.842763, -6763.22, -4.21000, .000367, &
                                       .0415, 218.8,  53.878000, -1331.22, -9.44523, .014025  /)

   real, parameter :: TMINLQU = MAPL_TICE - 40.0
   real, parameter :: TMINICE = MAPL_TICE + -95.

#else

!! Some parameters set by CLDPARAMS 

   real    :: CNV_BETA
   real    :: ANV_BETA
   real    :: LS_BETA
   real    :: RH00
   real    :: C_00
   real    :: LWCRIT
   real    :: C_ACC
   real    :: C_EV_R
   real    :: C_EV_S
   real    :: CCW_EVP_EFF
   real    :: CCI_EVP_EFF
   real    :: LS_SDQV2
   real    :: LS_SDQV3
   real    :: LS_SDQVT1
   real    :: ANV_SDQV2
   real    :: ANV_SDQV3
   real    :: ANV_SDQVT1
   real    :: ANV_TO_LS
   integer :: DISABLE_RAD
   integer :: ICE_SETTLE
   real    :: ANV_ICEFALL_C
   real    :: LS_ICEFALL_C
   real    :: REVAP_OFF_P
   real    :: CNVENVFC
   real    :: ANVENVFC
   real    :: SCENVFC
   real    :: LSENVFC
   real    :: T_ICE_ALL
   real    :: CNVDDRFC
   real    :: ANVDDRFC
   real    :: SCDDRFC
   real    :: LSDDRFC
   real    :: MIN_RI, MAX_RI, FAC_RI, MIN_RL, MAX_RL, FAC_RL
   integer :: FR_LS_WAT, FR_LS_ICE, FR_AN_WAT, FR_AN_ICE
   integer :: pdfflag
#endif

   real, parameter :: ESFAC        = MAPL_H2OMW/MAPL_AIRMW
   real, parameter :: T_ICE_MAX    = MAPL_TICE-10.0
   real, parameter :: RHO_I        =  916.8      ! Density of ice crystal in kg/m^3
   real, parameter :: RHO_W        = 1000.0      ! Density of liquid water in kg/m^3
   real, parameter :: MIN_CLD_FRAC = 1.0e-8
   real, parameter :: ZVIR = MAPL_RVAP/MAPL_RGAS - 1.
   real, parameter :: GORD = MAPL_GRAV/MAPL_RGAS
   real, parameter :: GFAC = 1.e5/MAPL_GRAV 
   real, parameter :: R_AIR = 3.47e-3 !m3 Pa kg-1K-1

  ! LDRADIUS4 constants
   real, parameter :: r13  = 1./3.
   real, parameter :: be   = r13 - 0.11
   real, parameter :: aewc = 0.13*(3./(4.*MAPL_PI*RHO_W*1.e3))**r13
   real, parameter :: aeic = 0.13*(3./(4.*MAPL_PI*RHO_I*1.e3))**r13

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains 

! GPU The GPU main routine call is smaller due to CUDA limit on
!     number of arguments permitted in call. Most inputs and outputs
!     are USE-associated in the GridComp

#ifdef _CUDA
   attributes(global) subroutine progno_cloud(IRUN,LM,DT)
#else
   subroutine progno_cloud( &
!!! first vars are (in) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IRUN, LM         , &
         DT               , &
         LATS_dev         , &
         PP_dev           , &
         ZZ_dev           , &
         PPE_dev          , &
         EXNP_dev         , &
         FRLAND_dev       , &
         KH_dev           , &
         mf_frc_dev       , &
!         wqtfac_dev       , &
!         whlfac_dev       , &
         wqt_dev          , &
         whl_dev          , &
         qt2_dev          , &
         hl2_dev          , &
         hlqt_dev         , &
         w2_dev           , &
         w3_dev           , &
         qt3_dev          , &
         hl3_dev          , &
         RMFDTR_dev       , &
         QLWDTR_dev       , &              
         QRN_CU_dev       , &
         CNV_UPDFRC_dev   , &
         SC_RMFDTR_dev    , &
         SC_QLWDTR_dev    , &
         SC_QIWDTR_dev    , &
         QRN_SC_dev       , &
         QSN_SC_dev       , &
         SC_UPDFRC_dev    , &
         U_dev            , &
         V_dev            , & 
         TH_dev           , &
         Q_dev            , &
         QLW_LS_dev       , &
         QLW_AN_dev       , &
         QIW_LS_dev       , &
         QIW_AN_dev       , &
         ANVFRC_dev       , &
         CLDFRC_dev       , &
         RAD_CLDFRC_dev   , &
         RAD_QV_dev       , &
         RAD_QL_dev       , &
         RAD_QI_dev       , &
         RAD_QR_dev       , &
         RAD_QS_dev       , &
         RAD_QG_dev       , &
         CLDREFFL_dev     , &
         CLDREFFI_dev     , &
         PRELS_dev        , &
         PRECU_dev        , &
         PREAN_dev        , &
         PRESC_dev        , &
         LSARF_dev        , &
         CUARF_dev        , &
         ANARF_dev        , &
         SCARF_dev        , &
         SNRLS_dev        , &
         SNRCU_dev        , &
         SNRAN_dev        , &
         SNRSC_dev        , &
         minrhcrit        , &
         maxrhcrit        , &
         turnrhcrit       , &
         QST3_dev         , &
         DZET_dev         , &
         QDDF3_dev        , &
         CNVFRC_dev       , &
         TROPP_dev        , &
         RHX_dev          , &
         REV_LS_dev       , &
         REV_AN_dev       , &
         REV_CN_dev       , &
         REV_SC_dev       , &
         RSU_LS_dev       , &
         RSU_AN_dev       , &
         RSU_CN_dev       , &
         RSU_SC_dev       , &
         ACLL_CN_dev,ACIL_CN_dev, &
         ACLL_AN_dev,ACIL_AN_dev, &
         ACLL_LS_dev,ACIL_LS_dev, &
         ACLL_SC_dev,ACIL_SC_dev, &
         PFL_CN_dev,PFI_CN_dev, &
         PFL_AN_dev,PFI_AN_dev, &
         PFL_LS_dev,PFI_LS_dev, &
         PFL_SC_dev,PFI_SC_dev, &
         PDFL_dev,PDFI_dev,FIXL_dev,FIXI_dev, &
         AUT_dev, EVAPC_dev, SDM_dev, SUBLC_dev,    &
         FRZ_TT_dev, DCNVL_dev, DCNVi_dev,      &
         ALPHT_dev, ALPH1_dev, ALPH2_dev,       &
         CFPDF_dev, RHCLR_dev,       &
         DQRL_dev,FRZ_PP_dev,               &
         VFALLICE_AN_dev,VFALLICE_LS_dev,   &
         VFALLWAT_AN_dev,VFALLWAT_LS_dev,   &
         VFALLSN_AN_dev,VFALLSN_LS_dev,VFALLSN_CN_dev,VFALLSN_SC_dev, &
         VFALLRN_AN_dev,VFALLRN_LS_dev,VFALLRN_CN_dev,VFALLRN_SC_dev,  &
         PDF_A_dev, & 
#ifdef PDFDIAG
         PDF_SIGW1_dev, PDF_SIGW2_dev, PDF_W1_dev, PDF_W2_dev, & 
         PDF_SIGTH1_dev, PDF_SIGTH2_dev, PDF_TH1_dev, PDF_TH2_dev, &
         PDF_SIGQT1_dev, PDF_SIGQT2_dev, PDF_QT1_dev, PDF_QT2_dev, &
         PDF_RQTTH_dev, PDF_RWTH_dev, PDF_RWQT_dev, &
#endif
         WTHV2_dev, wql_dev, &
         DOSHLW, &
         NACTL_dev,    &
         NACTI_dev,    &
         CONVPAR_OPTION )
#endif

      implicit none

#ifdef _CUDA
      integer, intent(in   ), value :: IRUN
      integer, intent(in   ), value :: LM
      real   , intent(in   ), value :: DT
#else
      real, intent(in   ), dimension(IRUN)      :: minrhcrit,maxrhcrit,turnrhcrit 

      integer, intent(in   )                    :: IRUN ! IM*JM
      integer, intent(in   )                    :: LM   ! LM
      real, intent(in   )                       :: DT   ! DT_MOIST
      real, intent(in   ), dimension(IRUN)      :: LATS_dev    ! LATS
      real, intent(in   ), dimension(IRUN,  LM) :: PP_dev      ! PLO
      real, intent(in   ), dimension(IRUN,  LM) :: ZZ_dev      ! ZLO
      real, intent(in   ), dimension(IRUN,0:LM) :: PPE_dev     ! CNV_PLE
      real, intent(in   ), dimension(IRUN,  LM) :: EXNP_dev    ! PK
      real, intent(in   ), dimension(IRUN     ) :: FRLAND_dev  ! FRLAND
      real, intent(in   ), dimension(IRUN,0:LM) :: KH_dev      ! KH
      real, intent(in   ), dimension(IRUN,  LM) :: mf_frc_dev  !
!      real, intent(in   ), dimension(IRUN,  LM) :: wqtfac_dev  !
!      real, intent(in   ), dimension(IRUN,  LM) :: whlfac_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: wqt_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: whl_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: qt2_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hl2_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hlqt_dev !
      real, intent(in   ), dimension(IRUN,  LM) :: w2_dev   !
      real, intent(in   ), dimension(IRUN,  LM) :: w3_dev   !
      real, intent(in   ), dimension(IRUN,  LM) :: qt3_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hl3_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: RMFDTR_dev  ! CNV_MFD
      real, intent(in   ), dimension(IRUN,  LM) :: QLWDTR_dev  ! CNV_DQLDT
      real, intent(inout), dimension(IRUN,  LM) :: QRN_CU_dev  ! CNV_PRC3 IS THIS INTENT IN?
      real, intent(inout), dimension(IRUN,  LM) :: CNV_UPDFRC_dev ! CNV_UPDF
      real, intent(in   ), dimension(IRUN,  LM) :: SC_RMFDTR_dev  ! MFDSHLW
      real, intent(in   ), dimension(IRUN,  LM) :: SC_QLWDTR_dev  ! DQLDTSHLW
      real, intent(in   ), dimension(IRUN,  LM) :: SC_QIWDTR_dev  ! DQIDTSHLW
      real, intent(inout), dimension(IRUN,  LM) :: QRN_SC_dev     ! SHLW_PRC3
      real, intent(inout), dimension(IRUN,  LM) :: QSN_SC_dev     ! SHLW_SNO3
      real, intent(in   ), dimension(IRUN,  LM) :: SC_UPDFRC_dev  ! UPDFSHLW
      real, intent(in   ), dimension(IRUN,  LM) :: U_dev  ! U1
      real, intent(in   ), dimension(IRUN,  LM) :: V_dev  ! V1
      real, intent(inout), dimension(IRUN,  LM) :: TH_dev ! TH1
      real, intent(inout), dimension(IRUN,  LM) :: Q_dev  ! Q1
      real, intent(inout), dimension(IRUN,  LM) :: QLW_LS_dev ! QLLS
      real, intent(inout), dimension(IRUN,  LM) :: QLW_AN_dev ! QLCN
      real, intent(inout), dimension(IRUN,  LM) :: QIW_LS_dev ! QILS
      real, intent(inout), dimension(IRUN,  LM) :: QIW_AN_dev ! QICN
      real, intent(inout), dimension(IRUN,  LM) :: ANVFRC_dev ! CLCN
      real, intent(inout), dimension(IRUN,  LM) :: CLDFRC_dev ! CLLS
      real, intent(  out), dimension(IRUN,  LM) :: RAD_CLDFRC_dev ! RAD_CF
      real, intent(  out), dimension(IRUN,  LM) :: RAD_QV_dev ! RAD_QV
      real, intent(  out), dimension(IRUN,  LM) :: RAD_QL_dev ! RAD_QL
      real, intent(  out), dimension(IRUN,  LM) :: RAD_QI_dev ! RAD_QI
      real, intent(  out), dimension(IRUN,  LM) :: RAD_QR_dev ! QRAIN
      real, intent(  out), dimension(IRUN,  LM) :: RAD_QS_dev ! QSNOW
      real, intent(  out), dimension(IRUN,  LM) :: RAD_QG_dev ! QGRAUPEL
      real, intent(  out), dimension(IRUN,  LM) :: CLDREFFL_dev ! CLDREFFL
      real, intent(  out), dimension(IRUN,  LM) :: CLDREFFI_dev ! CLDREFFI
      real, intent(  out), dimension(IRUN     ) :: PRELS_dev ! LS_PRC2
      real, intent(  out), dimension(IRUN     ) :: PRECU_dev ! CN_PRC2
      real, intent(  out), dimension(IRUN     ) :: PREAN_dev ! AN_PRC2
      real, intent(  out), dimension(IRUN     ) :: PRESC_dev ! SC_PRC2
      real, intent(  out), dimension(IRUN     ) :: LSARF_dev ! LS_ARFX
      real, intent(  out), dimension(IRUN     ) :: CUARF_dev ! CN_ARFX
      real, intent(  out), dimension(IRUN     ) :: ANARF_dev ! AN_ARFX
      real, intent(  out), dimension(IRUN     ) :: SCARF_dev ! SC_ARFX
      real, intent(  out), dimension(IRUN     ) :: SNRLS_dev ! LS_SNR
      real, intent(  out), dimension(IRUN     ) :: SNRCU_dev ! CN_SNR
      real, intent(  out), dimension(IRUN     ) :: SNRAN_dev ! AN_SNR
      real, intent(  out), dimension(IRUN     ) :: SNRSC_dev ! SC_SNR
      real, intent(in   ), dimension(IRUN,  LM) :: QST3_dev   ! QST3
      real, intent(in   ), dimension(IRUN,  LM) :: DZET_dev   ! DZET
      real, intent(in   ), dimension(IRUN,  LM) :: QDDF3_dev  ! QDDF3
      real, intent(in   ), dimension(IRUN)      :: CNVFRC_dev   ! CNV_FRACTION
      real, intent(in   ), dimension(IRUN)      :: TROPP_dev   ! TROPP

      real, intent(  out), dimension(IRUN,  LM) :: RHX_dev    ! RHX
      real, intent(  out), dimension(IRUN,  LM) :: REV_LS_dev ! REV_LS
      real, intent(  out), dimension(IRUN,  LM) :: REV_AN_dev ! REV_AN
      real, intent(  out), dimension(IRUN,  LM) :: REV_CN_dev ! REV_CN
      real, intent(  out), dimension(IRUN,  LM) :: REV_SC_dev ! REV_SC
      real, intent(  out), dimension(IRUN,  LM) :: RSU_LS_dev ! RSU_LS
      real, intent(  out), dimension(IRUN,  LM) :: RSU_AN_dev ! RSU_AN
      real, intent(  out), dimension(IRUN,  LM) :: RSU_CN_dev ! RSU_CN
      real, intent(  out), dimension(IRUN,  LM) :: RSU_SC_dev ! RSU_SC
      real, intent(  out), dimension(IRUN,  LM) :: ACLL_CN_dev ! ACLL_CN
      real, intent(  out), dimension(IRUN,  LM) :: ACIL_CN_dev ! ACIL_CN
      real, intent(  out), dimension(IRUN,  LM) :: ACLL_AN_dev ! ACLL_AN
      real, intent(  out), dimension(IRUN,  LM) :: ACIL_AN_dev ! ACIL_AN
      real, intent(  out), dimension(IRUN,  LM) :: ACLL_LS_dev ! ACLL_LS
      real, intent(  out), dimension(IRUN,  LM) :: ACIL_LS_dev ! ACIL_LS
      real, intent(  out), dimension(IRUN,  LM) :: ACLL_SC_dev ! ACLL_SC
      real, intent(  out), dimension(IRUN,  LM) :: ACIL_SC_dev ! ACIL_SC
      real, intent(  out), dimension(IRUN,0:LM) :: PFL_CN_dev ! PFL_CN
      real, intent(  out), dimension(IRUN,0:LM) :: PFI_CN_dev ! PFI_CN
      real, intent(  out), dimension(IRUN,0:LM) :: PFL_AN_dev ! PFL_AN
      real, intent(  out), dimension(IRUN,0:LM) :: PFI_AN_dev ! PFI_AN
      real, intent(  out), dimension(IRUN,0:LM) :: PFL_LS_dev ! PFL_LS
      real, intent(  out), dimension(IRUN,0:LM) :: PFI_LS_dev ! PFI_LS
      real, intent(  out), dimension(IRUN,0:LM) :: PFL_SC_dev ! PFL_SC
      real, intent(  out), dimension(IRUN,0:LM) :: PFI_SC_dev ! PFI_SC
      real, intent(  out), dimension(IRUN,  LM) :: PDFL_dev ! DlPDF
      real, intent(  out), dimension(IRUN,  LM) :: PDFI_dev ! DiPDF
      real, intent(  out), dimension(IRUN,  LM) :: FIXL_dev ! DlFIX
      real, intent(  out), dimension(IRUN,  LM) :: FIXI_dev ! DiFIX                  
      real, intent(  out), dimension(IRUN,  LM) :: AUT_dev   ! AUT
      real, intent(  out), dimension(IRUN,  LM) :: EVAPC_dev ! EVAPC
      real, intent(  out), dimension(IRUN,  LM) :: SDM_dev   ! SDM
      real, intent(  out), dimension(IRUN,  LM) :: SUBLC_dev ! SUBLC
      real, intent(  out), dimension(IRUN,  LM) :: FRZ_TT_dev ! FRZ_TT
      real, intent(  out), dimension(IRUN,  LM) :: FRZ_PP_dev ! FRZ_PP
      real, intent(  out), dimension(IRUN,  LM) :: DCNVL_dev ! DCNVL
      real, intent(  out), dimension(IRUN,  LM) :: DCNVi_dev ! DCNVi
      real, intent(  out), dimension(IRUN,  LM) :: ALPHT_dev ! ALPHT
      real, intent(  out), dimension(IRUN,  LM) :: ALPH1_dev ! ALPH1
      real, intent(  out), dimension(IRUN,  LM) :: ALPH2_dev ! ALPH2
      real, intent(  out), dimension(IRUN,  LM) :: CFPDF_dev ! CFPDF
      real, intent(  out), dimension(IRUN,  LM) :: RHCLR_dev ! RHCLR
      real, intent(  out), dimension(IRUN,  LM) :: DQRL_dev ! DQRL
      real, intent(  out), dimension(IRUN,  LM) :: VFALLICE_AN_dev ! VFALLICE_AN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLICE_LS_dev ! VFALLICE_LS
      real, intent(  out), dimension(IRUN,  LM) :: VFALLWAT_AN_dev ! VFALLWAT_AN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLWAT_LS_dev ! VFALLWAT_LS
      real, intent(  out), dimension(IRUN,  LM) :: VFALLSN_AN_dev ! VFALLSN_AN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLSN_LS_dev ! VFALLSN_LS
      real, intent(  out), dimension(IRUN,  LM) :: VFALLSN_CN_dev ! VFALLSN_CN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLSN_SC_dev ! VFALLSN_SC
      real, intent(  out), dimension(IRUN,  LM) :: VFALLRN_AN_dev ! VFALLRN_AN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLRN_LS_dev ! VFALLRN_LS
      real, intent(  out), dimension(IRUN,  LM) :: VFALLRN_CN_dev ! VFALLRN_CN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLRN_SC_dev ! VFALLRN_SC
      real, intent(inout), dimension(IRUN,  LM) :: PDF_A_dev
#ifdef PDFDIAG
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGW1_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGW2_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_W1_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_W2_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGTH1_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGTH2_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_TH1_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_TH2_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGQT1_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGQT2_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_QT1_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_QT2_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_RQTTH_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_RWTH_dev
      real, intent(  out), dimension(IRUN,  LM) :: PDF_RWQT_dev
#endif
      real, intent(  out), dimension(IRUN,  LM) :: WTHV2_dev
      real, intent(  out), dimension(IRUN,  LM) :: wql_dev
      real, intent(in   ), dimension(IRUN,  LM) :: NACTL_dev  ! NACTL
      real, intent(in   ), dimension(IRUN,  LM) :: NACTI_dev  ! NACTI


!!$      real, intent(  out), dimension(IRUN,  LM) :: LIQANMOVE_dev  ! LIQANMOVE
!!$      real, intent(  out), dimension(IRUN,  LM) :: ICEANMOVE_dev  ! ICEANMOVE
!!$      real, intent(  out), dimension(IRUN,  LM) :: DANCLD_dev     ! DANCLD
!!$      real, intent(  out), dimension(IRUN,  LM) :: DLSCLD_dev     ! DLSCLD
!!$      real, intent(  out), dimension(IRUN,  LM) :: CURAINMOVE_dev ! CURAINMOVE
!!$      real, intent(  out), dimension(IRUN,  LM) :: CUSNOWMOVE_dev ! CUSNOWMOVE

      LOGICAL, INTENT(IN) :: DOSHLW
      character(LEN=*), INTENT(IN)              :: CONVPAR_OPTION

#endif

! GPU The GPUs need to know how big local arrays are during compile-time
!     as the GPUs cannot allocate memory themselves. This command resets
!     this a priori size to LM for the CPU.
#ifndef GPU_MAXLEVS
#define GPU_MAXLEVS LM
#endif

      integer :: I , J , K , L, kd, ku

      integer :: FRACTION_REMOVAL

      real :: MASS, iMASS
      real :: TOTFRC
      real :: QRN_LS, QRN_AN, QRN_CU_1D, QRN_SC_1D
      real :: QSN_LS, QSN_AN, QSN_CU, QSN_SC_1D
      real :: QRN_ALL, QSN_ALL, QGR_ALL
      real :: QTMP1, QTMP2, QTMP3
      real :: TEMP
      real :: RHCRIT
      real :: AA3, BB3, ALPHA
      real :: VFALL, VFALLRN, VFALLSN

      real :: TOT_PREC_UPD
      real :: TOT_PREC_ANV
      real :: TOT_PREC_LS
      real :: TOT_PREC_SC
      real :: AREA_UPD_PRC
      real :: AREA_ANV_PRC
      real :: AREA_LS_PRC
      real :: AREA_SCUP_PRC

      real :: BETA
      real :: AREA_UPD_PRC_tolayer
      real :: AREA_ANV_PRC_tolayer
      real :: AREA_LS_PRC_tolayer
      real :: AREA_SCUP_PRC_tolayer

      real :: U_above,U_below
      real :: V_above,V_below
      real :: DZET_above,DZET_below

      real :: PRN_CU_above, PSN_CU_above
      real :: PRN_LS_above, PSN_LS_above
      real :: PRN_AN_above, PSN_AN_above
      real :: PRN_SC_above, PSN_SC_above

      real :: EVAP_DD_CU_above, SUBL_DD_CU_above
      real :: EVAP_DD_LS_above, SUBL_DD_LS_above
      real :: EVAP_DD_AN_above, SUBL_DD_AN_above
      real :: EVAP_DD_SC_above, SUBL_DD_SC_above

      logical :: use_autoconv_timescale

      real :: NI, NL, TROPICAL, EXTRATROPICAL

      real :: LSPDFLIQNEW, LSPDFICENEW, LSPDFFRACNEW

! These are in constant memory in CUDA and are set in the GridComp
#ifndef _CUDA
         CNV_BETA      = CLDPARAMS%CNV_BETA  ! Area factor for convective rain showers (non-dim)
         ANV_BETA      = CLDPARAMS%ANV_BETA  ! Area factor for anvil rain showers (non-dim)
         LS_BETA       = CLDPARAMS%LS_BETA   ! Area factor for Large Scale rain showers (non-dim)
         RH00          = CLDPARAMS%RH_CRIT   ! Critical relative humidity
         C_00          = CLDPARAMS%AUTOC_LS
         LWCRIT        = CLDPARAMS%QC_CRIT_LS
         C_ACC         = CLDPARAMS%ACCRETION
         C_EV_R        = CLDPARAMS%RAIN_REVAP_FAC
         C_EV_S        = CLDPARAMS%SNOW_REVAP_FAC
         CCW_EVP_EFF   = CLDPARAMS%CCW_EVAP_EFF
         CCI_EVP_EFF   = CLDPARAMS%CCI_EVAP_EFF
         LS_SDQV2      = CLDPARAMS%LS_SUND_INTER
         LS_SDQV3      = CLDPARAMS%LS_SUND_COLD
         LS_SDQVT1     = CLDPARAMS%LS_SUND_TEMP1
         ANV_SDQV2     = CLDPARAMS%ANV_SUND_INTER
         ANV_SDQV3     = CLDPARAMS%ANV_SUND_COLD
         ANV_SDQVT1    = CLDPARAMS%ANV_SUND_TEMP1
         ANV_TO_LS     = CLDPARAMS%ANV_TO_LS_TIME
         DISABLE_RAD   = INT( CLDPARAMS%DISABLE_RAD )
         ICE_SETTLE    = NINT( CLDPARAMS%ICE_SETTLE )
         ANV_ICEFALL_C = CLDPARAMS%ANV_ICEFALL
         LS_ICEFALL_C  = CLDPARAMS%LS_ICEFALL
         REVAP_OFF_P   = CLDPARAMS%REVAP_OFF_P
         CNVENVFC      = CLDPARAMS%CNV_ENVF
         ANVENVFC      = CLDPARAMS%ANV_ENVF
         SCENVFC       = CLDPARAMS%SC_ENVF
         LSENVFC       = CLDPARAMS%LS_ENVF
         T_ICE_ALL     = CLDPARAMS%ICE_RAMP + T_ICE_MAX
         CNVDDRFC      = CLDPARAMS%CNV_DDRF
         ANVDDRFC      = CLDPARAMS%ANV_DDRF
         LSDDRFC       = CLDPARAMS%LS_DDRF
         FR_LS_WAT     = INT( CLDPARAMS%FR_LS_WAT )
         FR_LS_ICE     = INT( CLDPARAMS%FR_LS_ICE )
         FR_AN_WAT     = INT( CLDPARAMS%FR_AN_WAT )
         FR_AN_ICE     = INT( CLDPARAMS%FR_AN_ICE )
         MIN_RL        = CLDPARAMS%MIN_RL
         MIN_RI        = CLDPARAMS%MIN_RI
         MAX_RL        = CLDPARAMS%MAX_RL
         MAX_RI        = CLDPARAMS%MAX_RI
         FAC_RL        = CLDPARAMS%FAC_RL
         FAC_RI        = CLDPARAMS%FAC_RI
         PDFFLAG       = INT(CLDPARAMS%PDFSHAPE)
#endif

      use_autoconv_timescale = .false.
      wthv2_dev = 0.0
      wql_dev = 0.0

#ifdef _CUDA
      i = (blockidx%x - 1) * blockdim%x + threadidx%x

      RUN_LOOP: IF ( I <= IRUN ) THEN
#else
      RUN_LOOP: DO I = 1, IRUN
#endif

         K_LOOP: DO K = 1, LM         
            if (K == 1) then
               TOT_PREC_UPD = 0.
               TOT_PREC_ANV = 0.
               TOT_PREC_LS  = 0.
               TOT_PREC_SC  = 0.

               AREA_UPD_PRC = 0.
               AREA_ANV_PRC = 0.
               AREA_LS_PRC  = 0.
               AREA_SCUP_PRC= 0.
            end if

            if (K == LM ) then
               !! ZERO DIAGNOSTIC OUTPUTS BEFORE SHOWERS !!
               
               PRELS_dev(I) = 0.
               PRECU_dev(I) = 0.
               PRESC_dev(I) = 0.
               PREAN_dev(I) = 0.

               SNRCU_dev(I) = 0. 
               SNRLS_dev(I) = 0. 
               SNRSC_dev(I) = 0. 
               SNRAN_dev(I) = 0. 

               LSARF_dev(I) = 0.
               CUARF_dev(I) = 0.
               ANARF_dev(I) = 0.
               SCARF_dev(I) = 0.
            end if

            !Zero out/initialize precips, except QRN_CU which comes from RAS 
            QRN_LS    = 0.
            QRN_AN    = 0.
            QRN_CU_1D = 0.
            QRN_SC_1D = 0.
            QSN_LS    = 0.
            QSN_AN    = 0.
            QSN_CU    = 0.
            QSN_SC_1D = 0.
            VFALL     = 0.

            RAD_QV_dev(I,K)     = 0.
            RAD_QL_dev(I,K)     = 0.
            RAD_QI_dev(I,K)     = 0.
            RAD_QR_dev(I,K)     = 0.
            RAD_QS_dev(I,K)     = 0.
            RAD_QG_dev(I,K)     = 0.
            RAD_CLDFRC_dev(I,K) = 0.
            CLDREFFL_dev(I,K)   = 0.
            CLDREFFI_dev(I,K)   = 0.

            PFL_CN_dev(I,K) = 0.
            PFI_CN_dev(I,K) = 0.
            PFL_SC_dev(I,K) = 0.
            PFI_SC_dev(I,K) = 0.
            PFL_AN_dev(I,K) = 0.
            PFI_AN_dev(I,K) = 0.
            PFL_LS_dev(I,K) = 0.
            PFI_LS_dev(I,K) = 0.

            IF (K == 1) THEN
               PFL_CN_dev(I,0) = 0.
               PFI_CN_dev(I,0) = 0.
               PFL_SC_dev(I,0) = 0.
               PFI_SC_dev(I,0) = 0.
               PFL_AN_dev(I,0) = 0.
               PFI_AN_dev(I,0) = 0.
               PFL_LS_dev(I,0) = 0.
               PFI_LS_dev(I,0) = 0.
            END IF

            ! Initialize other diagnostics 

            RHX_dev(I,K) = MAPL_UNDEF
            REV_LS_dev(I,K) = MAPL_UNDEF
            REV_AN_dev(I,K) = MAPL_UNDEF
            IF(CONVPAR_OPTION .ne. 'GF')REV_CN_dev(I,K) = MAPL_UNDEF
            REV_SC_dev(I,K) = MAPL_UNDEF
            RSU_LS_dev(I,K) = MAPL_UNDEF
            RSU_AN_dev(I,K) = MAPL_UNDEF
            IF(CONVPAR_OPTION .ne. 'GF')RSU_CN_dev(I,K) = MAPL_UNDEF
            RSU_SC_dev(I,K) = MAPL_UNDEF
            IF(CONVPAR_OPTION .ne. 'GF')ACLL_CN_dev(I,K) = MAPL_UNDEF
            IF(CONVPAR_OPTION .ne. 'GF')ACIL_CN_dev(I,K) = MAPL_UNDEF
            ACLL_SC_dev(I,K) = MAPL_UNDEF
            ACIL_SC_dev(I,K) = MAPL_UNDEF
            ACLL_AN_dev(I,K) = MAPL_UNDEF
            ACIL_AN_dev(I,K) = MAPL_UNDEF
            ACLL_LS_dev(I,K) = MAPL_UNDEF
            ACIL_LS_dev(I,K) = MAPL_UNDEF

            PDFL_dev(I,K) = MAPL_UNDEF
            PDFI_dev(I,K) = MAPL_UNDEF
            FIXL_dev(I,K) = MAPL_UNDEF
            FIXI_dev(I,K) = MAPL_UNDEF
            AUT_dev(I,K) = MAPL_UNDEF
            EVAPC_dev(I,K) = MAPL_UNDEF
            SDM_dev(I,K) = MAPL_UNDEF
            SUBLC_dev(I,K) = MAPL_UNDEF
            FRZ_TT_dev(I,K) = MAPL_UNDEF
            FRZ_PP_dev(I,K) = MAPL_UNDEF
            DCNVL_dev(I,K) = MAPL_UNDEF
            DCNVi_dev(I,K) = MAPL_UNDEF
            ALPHT_dev(I,K) = MAPL_UNDEF
            ALPH1_dev(I,K) = MAPL_UNDEF
            ALPH2_dev(I,K) = MAPL_UNDEF
            CFPDF_dev(I,K) = MAPL_UNDEF
            RHCLR_dev(I,K) = MAPL_UNDEF
            DQRL_dev(I,K) = MAPL_UNDEF
            VFALLICE_AN_dev(I,K) = MAPL_UNDEF
            VFALLICE_LS_dev(I,K) = MAPL_UNDEF
            VFALLWAT_AN_dev(I,K) = MAPL_UNDEF
            VFALLWAT_LS_dev(I,K) = MAPL_UNDEF
            VFALLSN_AN_dev(I,K) = MAPL_UNDEF
            VFALLSN_LS_dev(I,K) = MAPL_UNDEF
            VFALLSN_CN_dev(I,K) = MAPL_UNDEF
            VFALLSN_SC_dev(I,K) = MAPL_UNDEF
            VFALLRN_AN_dev(I,K) = MAPL_UNDEF
            VFALLRN_LS_dev(I,K) = MAPL_UNDEF
            VFALLRN_CN_dev(I,K) = MAPL_UNDEF
            VFALLRN_SC_dev(I,K) = MAPL_UNDEF


            ! Copy QRN_CU into a temp scalar
            !QRN_CU_1D = QRN_CU_dev(I,K)
            !- GF scheme handles its own conv precipitation.
            !- => QRN_CU_1D is equal to zero, in this case.
            !- Otherwise, the surface convective precip must be set to zero inside GF main routine.
            IF(CONVPAR_OPTION .ne. 'GF') QRN_CU_1D = QRN_CU_dev(I,K)

            MASS =  ( PPE_dev(I,K) - PPE_dev(I,K-1) )*100./MAPL_GRAV  ! layer-mass (kg/m**2)

            iMASS = 1.0 / MASS

            QRN_SC_1D = QRN_SC_dev(I,K)
            QSN_SC_1D = QSN_SC_dev(I,K)

            TEMP =  EXNP_dev(I,K) * TH_dev(I,K) 

            FRZ_PP_dev(I,K) = 0.00

            ! Buoyancy calculation for SHOC, if not using ADG PDF
            if (PDFFLAG .ne. 5) then
!              wthv2_dev(I,K) = -1.0*(TH_dev(I,K)*(1.+ESFAC*Q_dev(I,K))/MAPL_GRAV)*KH_dev(I,K)
              wthv2_dev(I,K) = -1.0*KH_dev(I,K)
              if (K.lt.LM .and. K.gt.1) then
                wthv2_dev(I,K) = wthv2_dev(I,K)*(TH_dev(I,K-1)*(1.+ESFAC*Q_dev(I,K-1))-TH_dev(I,K+1)*(1.+ESFAC*Q_dev(I,K+1))) &
                                               / (ZZ_dev(I,K-1)-ZZ_dev(I,K+1))
              else if (K.eq.LM) then
                wthv2_dev(I,K) = wthv2_dev(I,K)*(TH_dev(I,K-1)*(1.+ESFAC*Q_dev(I,K-1))-TH_dev(I,K)*(1.+ESFAC*Q_dev(I,K))) &
                                               / (ZZ_dev(I,K-1)-ZZ_dev(I,K))
              else if (K.eq.1) then
                wthv2_dev(I,K) = wthv2_dev(I,K)*(TH_dev(I,K)*(1.+ESFAC*Q_dev(I,K))-TH_dev(I,K+1)*(1.+ESFAC*Q_dev(I,K+1))) &
                                               / (ZZ_dev(I,K)-ZZ_dev(I,K+1))
              end if
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Total Condensate Source
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            FIXL_dev(I,K) = QLW_AN_dev(I,K) + QLW_LS_dev(I,K)
            FIXI_dev(I,K) = QIW_AN_dev(I,K) + QIW_LS_dev(I,K)

            CALL fix_up_clouds(    &
                  Q_dev(I,K)     , &
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QIW_LS_dev(I,K), & 
                  CLDFRC_dev(I,K), &   
                  QLW_AN_dev(I,K), & 
                  QIW_AN_dev(I,K), & 
                  ANVFRC_dev(I,K))

            FIXL_dev(I,K) = -( QLW_AN_dev(I,K) + QLW_LS_dev(I,K) - FIXL_dev(I,K) ) / DT 
            FIXI_dev(I,K) = -( QIW_AN_dev(I,K) + QIW_LS_dev(I,K) - FIXI_dev(I,K) ) / DT 
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            FRZ_TT_dev(I,K) = QIW_AN_dev(I,K) + QIW_LS_dev(I,K)

            IF(USE_AEROSOL_NN) THEN
            CALL meltfrz_inst (    &
                  DT             , &
                  CNVFRC_dev(I)  , &
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QLW_AN_dev(I,K), &
                  QIW_LS_dev(I,K), &
                  QIW_AN_dev(I,K))
            ELSE
            CALL meltfrz (         &
                  DT             , &
                  CNVFRC_dev(I)  , &
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QIW_LS_dev(I,K))
            CALL meltfrz (         &
                  DT             , &
                  CNVFRC_dev(I)  , &
                  TEMP           , &
                  QLW_AN_dev(I,K), & 
                  QIW_AN_dev(I,K))
            ENDIF

            FRZ_TT_dev(I,K) = ( QIW_AN_dev(I,K) + QIW_LS_dev(I,K) - FRZ_TT_dev(I,K) ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            DCNVi_dev(I,K) = QIW_AN_dev(I,K)
            DCNVL_dev(I,K) = QLW_AN_dev(I,K)
         !  cnvsrc is now handled inside convection codes
            DCNVi_dev(I,K) = ( QIW_AN_dev(I,K) - DCNVi_dev(I,K) ) / DT
            DCNVL_dev(I,K) = ( QLW_AN_dev(I,K) - DCNVL_dev(I,K) ) / DT
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            PDFL_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)
            PDFI_dev(I,K) = QIW_LS_dev(I,K)+QIW_AN_dev(I,K)

            if (k == 1 .or. k == lm) then
               U_above = 0.0
               U_below = 0.0
               V_above = 0.0
               V_below = 0.0
               DZET_above = 0.0
               DZET_below = 0.0
            else
               U_above = U_dev(i,k-1)
               U_below = U_dev(i,k+1)
               V_above = V_dev(i,k-1)
               V_below = V_dev(i,k+1)
               DZET_above = DZET_dev(i,k-1)
               DZET_below = DZET_dev(i,k+1)
            end if
  
            call pdf_spread ( &
                  PP_dev(I,K),PPE_dev(I,LM),TROPP_dev(I),&
                  minrhcrit(I), maxrhcrit(I), turnrhcrit(I),&
                  ALPHA, ALPHT_dev(I,K) ) 

            ! impose a minimum amount of variability
            ALPHA    = MAX(  ALPHA , 1.0 - RH00 )
 
            LSPDFLIQNEW = QLW_LS_dev(I,K)
            LSPDFICENEW = QIW_LS_dev(I,K)
            LSPDFFRACNEW= CLDFRC_dev(I,K)

            call hystpdf(      &
                  DT             , &
                  ALPHA          , &
                  PDFFLAG        , &
                  CNVFRC_dev(I)  , &                  
                  PP_dev(I,K)    , &
                  ZZ_dev(I,K)    , &
                  Q_dev(I,K)     , &
                  QLW_LS_dev(I,K), &
                  QLW_AN_dev(I,K), &
                  QIW_LS_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  TEMP           , &
                  CLDFRC_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  NACTL_dev(I,K),  &
                  NACTI_dev(I,K),  &
                  whl_dev(I,K),        &
                  wqt_dev(I,K),        &
!                  wqtfac_dev(I,K),     &
!                  whlfac_dev(I,K),     &
                  hl2_dev(I,K),        &
                  qt2_dev(I,K),        &
                  hlqt_dev(I,K),       & 
                  w3_dev(I,K),         &
                  w2_dev(I,K),         &
                  qt3_dev(I,K),        &
                  hl3_dev(I,K),        &
                  mf_frc_dev(I,K),     &
                  PDF_A_dev(I,K),      &  ! can remove these after development
#ifdef PDFDIAG
                  PDF_SIGW1_dev(I,K),  &
                  PDF_SIGW2_dev(I,K),  &
                  PDF_W1_dev(I,K),     &
                  PDF_W2_dev(I,K),     &
                  PDF_SIGTH1_dev(I,K), &
                  PDF_SIGTH2_dev(I,K), &
                  PDF_TH1_dev(I,K),    &
                  PDF_TH2_dev(I,K),    &
                  PDF_SIGQT1_dev(I,K), &
                  PDF_SIGQT2_dev(I,K), &
                  PDF_QT1_dev(I,K),    &
                  PDF_QT2_dev(I,K),    &
                  PDF_RQTTH_dev(I,K),  &
                  PDF_RWTH_dev(I,K),   &
                  PDF_RWQT_dev(I,K),   &
#endif
                  WTHV2_dev(I,K),      &
                  wql_dev(I,K),        &
                  USE_AEROSOL_NN)

            LSPDFLIQNEW = QLW_LS_dev(I,K) - LSPDFLIQNEW
            LSPDFICENEW = QIW_LS_dev(I,K) - LSPDFICENEW
            LSPDFFRACNEW= CLDFRC_dev(I,K) - LSPDFFRACNEW

            RHX_dev(I,K)   = Q_dev(I,K)/QSAT( TEMP, PP_dev(I,K) )
            CFPDF_dev(I,K) = CLDFRC_dev(I,K)
            PDFL_dev(I,K)  = ( QLW_LS_dev(I,K) + QLW_AN_dev(I,K) - PDFL_dev(I,K) ) / DT 
            PDFI_dev(I,K)  = ( QIW_LS_dev(I,K) + QIW_AN_dev(I,K) - PDFI_dev(I,K) ) / DT 

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            TOTFRC = CLDFRC_dev(I,K) + ANVFRC_dev(I,K)

            IF ( TOTFRC > 1.00 ) THEN
               CLDFRC_dev(I,k) = CLDFRC_dev(I,k)*(1.00 / TOTFRC )
               ANVFRC_dev(I,k) = ANVFRC_dev(I,k)*(1.00 / TOTFRC )
            END IF

            TOTFRC = CLDFRC_dev(I,K) + ANVFRC_dev(I,K)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! CONDENSATE/FRACTION SOURCES FINISHED. NOW LOSE CLOUD CONDENSATE !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Total Condensate Sink
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!       E  V  A  P  O  R  A  T  I  O  N
            !!                A  N  D 
            !!       S  U  B  L  I  M  A  T  I  O  N
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            EVAPC_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)
            SUBLC_dev(I,K) = QIW_LS_dev(I,K)+QIW_AN_dev(I,K)

            if (USE_AEROSOL_NN) then
               NL = NACTL_dev(I,K)
               NI = NACTI_dev(I,K)
            else
               NL = 50.e6
               NI =  5.e6
            endif

             ! 'Anvil' partition from RAS/Parameterized not done in hystpdf

            RHCRIT = 1.0
            call evap3(            &
                  DT             , &
                  CCW_EVP_EFF    , &
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  NL , &
                  NI , &
                  QST3_dev(I,K)  )  

            RHCRIT = 1.0 - ALPHA
            call subl3(            &
                  DT             , & 
                  CCI_EVP_EFF    , &
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  NL , &
                  NI , &
                  QST3_dev(I,K)  ) 

            EVAPC_dev(I,K) = ( EVAPC_dev(I,K) - (QLW_LS_dev(I,K)+QLW_AN_dev(I,K)) ) / DT
            SUBLC_dev(I,K) = ( SUBLC_dev(I,K) - (QIW_LS_dev(I,K)+QIW_AN_dev(I,K)) ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!       A U T O C O N V E R S I  O  N
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            !           FRACTION_REMOVAL = 0 -> none
            !                              1 -> constant in-cloud QC
            !                              2 -> trim high edge of PDF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            AUT_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)

            FRACTION_REMOVAL = fr_ls_wat

            call autocon3(         &
                  DT             , &
                  QLW_LS_dev(I,K), &
                  QRN_LS         , &
                  TEMP           , &
                  PP_dev(I,K)    , &
                  KH_dev(I,K-1)  , &
                  CLDFRC_dev(I,K), &
                  LS_SDQV2       , &
                  LS_SDQV3       , &   
                  LS_SDQVT1      , &
                  DZET_dev(I,K)  , &
                  VFALL          , &
                  FRACTION_REMOVAL )
         
            VFALLWAT_LS_dev(I,K) = VFALL

            FRACTION_REMOVAL = fr_an_wat

            call autocon3(         &
                  DT             , &
                  QLW_AN_dev(I,K), &
                  QRN_AN         , &
                  TEMP           , &
                  PP_dev(I,K)    , &
                  KH_dev(I,K)    , &
                  ANVFRC_dev(I,K), &
                  ANV_SDQV2      , &
                  ANV_SDQV3      , &
                  ANV_SDQVT1     , &
                  DZET_dev(I,K)  , &
                  VFALL          , &
                  FRACTION_REMOVAL )

            VFALLWAT_AN_dev(I,K) = VFALL
            AUT_dev(I,K) = ( AUT_dev(I,K) - ( QLW_AN_dev(I,K) + QLW_LS_dev(I,K) ) )/DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Ice cloud settling
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            SDM_dev(I,K) = QIW_AN_dev(I,K)+QIW_LS_dev(I,K)

            ! Parameterized (RAS) Ice Fall
            ! ----------------------------

            ! WMP: Adjustments to resolved scale ice fall speed options 
                                      FRACTION_REMOVAL = fr_ls_ice
            if (CNVFRC_dev(I) >= 0.5) FRACTION_REMOVAL = fr_an_ice

            SELECT CASE( ICE_SETTLE )
            CASE( 0 )
             ! MERRA-2 Formulation
              TROPICAL      = ANV_ICEFALL_C*1.0
              EXTRATROPICAL = ANV_ICEFALL_C*0.0
            CASE( 1 )
              TROPICAL      = ANV_ICEFALL_C
              EXTRATROPICAL =  LS_ICEFALL_C
            END SELECT

            CALL SETTLE_VEL(       &
                  QIW_AN_dev(I,K), &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  ANVFRC_dev(I,K), &
                  KH_dev(I,K-1)  , &
                  VFALL          , &
                  EXTRATROPICAL, TROPICAL, CNVFRC_dev(I) )
            VFALLICE_AN_dev(I,K) = VFALL

            CALL ICEFALL(          &
                  QIW_AN_dev(I,K), &
                  DZET_dev(I,K)  , &
                  QSN_AN         , &
                  VFALL          , &
                  ANVFRC_dev(I,K), &
                  DT             , &
                  FRACTION_REMOVAL )

            ! Resolved Scale Ice Fall
            ! -----------------------

            ! WMP: Adjustments to resolved scale ice fall speed options 
            SELECT CASE( ICE_SETTLE )
            CASE( 0 )
             ! MERRA-2 Formulation
              TROPICAL      =  LS_ICEFALL_C*0.0
              EXTRATROPICAL =  LS_ICEFALL_C*1.0
            CASE( 1 )
              TROPICAL      = ANV_ICEFALL_C
              EXTRATROPICAL =  LS_ICEFALL_C
            END SELECT

            CALL SETTLE_VEL(       &
                  QIW_LS_dev(I,K), &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  CLDFRC_dev(I,K), &
                  KH_dev(I,K-1)  , &
                  VFALL          , &
                  EXTRATROPICAL, TROPICAL, CNVFRC_dev(I) )
            VFALLICE_LS_dev(I,K) = VFALL

            CALL ICEFALL(          &
                  QIW_LS_dev(I,K), &
                  DZET_dev(I,K)  , &
                  QSN_LS         , &
                  VFALL          , &
                  CLDFRC_dev(I,K), &
                  DT             , &
                  FRACTION_REMOVAL )

            SDM_dev(I,K) = ( SDM_dev(I,K) - (QIW_LS_dev(I,K) + QIW_AN_dev(I,K)) )/DT
      
            DQRL_dev(I,K) = ( QRN_LS + QRN_AN + QSN_LS + QSN_AN ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Add in convective rain 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! CU-FREEZE 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Also "freeze" out any conv. precip that needs
            ! to be since this isnt done in RAS. This is
            ! precip w/ large particles, so freezing is 
            ! strict. Check up on this!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP2 = 0.

            if ( TEMP < MAPL_TICE ) then
               QTMP2     = QRN_CU_1D
               QSN_CU    = QRN_CU_1D
               QRN_CU_1D = 0.
               TEMP      = TEMP + QSN_CU*(MAPL_ALHS-MAPL_ALHL) / MAPL_CP
            end if
      
            FRZ_PP_dev(I,K) = FRZ_PP_dev(I,K) +  QTMP2/DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !----------------------------------------------------------------------------------------------
            ! Column will now be swept from top-down for precip accumulation/accretion/re-evaporation

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            AREA_LS_PRC_tolayer  = 0.0
            AREA_UPD_PRC_tolayer = 0.0
            AREA_SCUP_PRC_tolayer = 0.0
            AREA_ANV_PRC_tolayer = 0.0

            TOT_PREC_UPD  = TOT_PREC_UPD + ( ( QRN_CU_1D + QSN_CU ) * MASS )
            AREA_UPD_PRC  = AREA_UPD_PRC + ( CNV_UPDFRC_dev(I,K)* ( QRN_CU_1D + QSN_CU )* MASS )

            TOT_PREC_SC  = TOT_PREC_SC + ( ( QRN_SC_1D + QSN_SC_1D ) * MASS )
            AREA_SCUP_PRC  = AREA_SCUP_PRC + ( SC_UPDFRC_dev(I,K)* ( QRN_SC_1D + QSN_SC_1D )* MASS )

            TOT_PREC_ANV  = TOT_PREC_ANV + ( ( QRN_AN + QSN_AN ) * MASS )
            AREA_ANV_PRC  = AREA_ANV_PRC + ( ANVFRC_dev(I,K)* ( QRN_AN + QSN_AN )* MASS )

            TOT_PREC_LS   = TOT_PREC_LS  + ( ( QRN_LS + QSN_LS ) * MASS )
            AREA_LS_PRC   = AREA_LS_PRC  + ( CLDFRC_dev(I,K)* ( QRN_LS + QSN_LS )* MASS )

            if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC_tolayer = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
            if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC_tolayer = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
            if ( TOT_PREC_SC > 0.0 ) AREA_SCUP_PRC_tolayer = MAX( AREA_SCUP_PRC/TOT_PREC_SC, 1.E-6 )
            if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC_tolayer  = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

            AREA_LS_PRC_tolayer  = LS_BETA  * AREA_LS_PRC_tolayer
            AREA_UPD_PRC_tolayer = CNV_BETA * AREA_UPD_PRC_tolayer
            AREA_ANV_PRC_tolayer = ANV_BETA * AREA_ANV_PRC_tolayer

            IF (K == LM) THEN ! Weve accumulated over the whole column
               if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC  = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
               if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC  = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
               if ( TOT_PREC_SC  > 0.0 ) AREA_SCUP_PRC = MAX( AREA_SCUP_PRC/TOT_PREC_SC, 1.E-6 )
               if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC   = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

               AREA_LS_PRC  = LS_BETA  * AREA_LS_PRC
               AREA_UPD_PRC = CNV_BETA * AREA_UPD_PRC
               AREA_ANV_PRC = ANV_BETA * AREA_ANV_PRC

               !! "couple" to diagnostic areal fraction output 
               !! Intensity factor in PRECIP3 is floored at
               !! 1.0. So this is fair.

               LSARF_dev(I) = MIN( AREA_LS_PRC,   1.0 )
               CUARF_dev(I) = MIN( AREA_UPD_PRC,  1.0 )
               SCARF_dev(I) = MIN( AREA_SCUP_PRC, 1.0 )
               ANARF_dev(I) = MIN( AREA_ANV_PRC,  1.0 )
            END IF

            QRN_ALL = 0.
            QSN_ALL = 0.
            QGR_ALL = 0.

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! GET SOME MICROPHYSICAL QUANTITIES 

            CALL MICRO_AA_BB_3( TEMP,PP_dev(I,K),QST3_dev(I,K),AA3,BB3 )

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP1 = QLW_LS_dev(I,K) + QLW_AN_dev(I,K)
            QTMP2 = QIW_LS_dev(I,K) + QIW_AN_dev(I,K)

        !- Convection schemes must handle their own conv precipitation.
        !- Otherwise, the surface convective precip must be set to zero inside GF main routine.
        IF(CONVPAR_OPTION .ne. 'GF') then
            !  Convective
            ! ----------

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , &
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_CU_1D       , &
                  QSN_CU          , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , &
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_UPD_PRC_tolayer    , &
                  PRECU_dev(I)    , & 
                  SNRCU_dev(I)    , & 
                  PRN_CU_above    , &
                  PSN_CU_above    , &
                  EVAP_DD_CU_above, &
                  SUBL_DD_CU_above, &
                  REV_CN_dev(I,K) , &
                  RSU_CN_dev(I,K) , &
                  ACLL_CN_dev(I,K), &
                  ACIL_CN_dev(I,K), &
                  PFL_CN_dev(I,K) , &
                  PFI_CN_dev(I,K) , &
                  VFALLRN         , &
                  VFALLSN         , &
                  FRZ_PP_dev(I,K) , &
                  CNVENVFC, CNVDDRFC )

            VFALLSN_CN_dev(I,K) = VFALLSN
            VFALLRN_CN_dev(I,K) = VFALLRN

            if (.not.use_autoconv_timescale) then
               if (VFALLSN.NE.0.) then
                  QSN_ALL = QSN_ALL + PFI_CN_dev(I,K)/VFALLSN
               end if
               if (VFALLRN.NE.0.) then
                  QRN_ALL = QRN_ALL + PFL_CN_dev(I,K)/VFALLRN
               end if
            end if
        ENDIF

        if(DOSHLW) THEN 
            ! Shallow convective
            ! ------------------

            SCDDRFC = 0.

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , &
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_SC_1D       , &
                  QSN_SC_1D       , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , &
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_SCUP_PRC_tolayer    , &
                  PRESC_dev(I)    , & 
                  SNRSC_dev(I)    , & 
                  PRN_SC_above    , &
                  PSN_SC_above    , &
                  EVAP_DD_SC_above, &
                  SUBL_DD_SC_above, &
                  REV_SC_dev(I,K) , &
                  RSU_SC_dev(I,K) , &
                  ACLL_SC_dev(I,K), &
                  ACIL_SC_dev(I,K), &
                  PFL_SC_dev(I,K) , &
                  PFI_SC_dev(I,K) , &
                  VFALLRN         , &
                  VFALLSN         , &
                  FRZ_PP_dev(I,K) , &
                  SCENVFC, SCDDRFC )

            VFALLSN_SC_dev(I,K) = VFALLSN
            VFALLRN_SC_dev(I,K) = VFALLRN

            if (.not.use_autoconv_timescale) then
               if (VFALLSN.NE.0.) then
                  QSN_ALL = QSN_ALL + PFI_SC_dev(I,K)/VFALLSN
               end if
               if (VFALLRN.NE.0.) then
                  QRN_ALL = QRN_ALL + PFL_SC_dev(I,K)/VFALLRN
               end if
            end if

        ENDIF
            ! Anvil
            ! -----

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , & 
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_AN          , &
                  QSN_AN          , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , & 
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_ANV_PRC_tolayer    , &
                  PREAN_dev(I)    , & 
                  SNRAN_dev(I)    , &
                  PRN_AN_above    , &
                  PSN_AN_above    , &
                  EVAP_DD_AN_above, &
                  SUBL_DD_AN_above, &
                  REV_AN_dev(I,K) , &
                  RSU_AN_dev(I,K) , &
                  ACLL_AN_dev(I,K), &
                  ACIL_AN_dev(I,K), &
                  PFL_AN_dev(I,K) , &
                  PFI_AN_dev(I,K) , &
                  VFALLRN         , &
                  VFALLSN         , &
                  FRZ_PP_dev(I,K) , &
                  ANVENVFC, ANVDDRFC )

            VFALLSN_AN_dev(I,K) = VFALLSN
            VFALLRN_AN_dev(I,K) = VFALLRN

            if (.not.use_autoconv_timescale) then
               if (VFALLSN.NE.0.) then
                  QSN_ALL = QSN_ALL + PFI_AN_dev(I,K)/VFALLSN
               end if
               if (VFALLRN.NE.0.) then
                  QRN_ALL = QRN_ALL + PFL_AN_dev(I,K)/VFALLRN
               end if
            end if

            ! Largescale
            ! ----------

            call  PRECIP3(          &
                  K,LM            , &
                  DT              , &  
                  FRLAND_dev(I)   , &
                  RHCRIT          , &
                  QRN_LS          , &
                  QSN_LS          , &
                  QTMP1           , &
                  QTMP2           , &
                  TEMP            , &
                  Q_dev(I,K)      , &
                  mass            , &
                  imass           , &
                  PP_dev(I,K)     , &
                  DZET_dev(I,K)   , &
                  QDDF3_dev(I,K)  , &
                  AA3             , &
                  BB3             , &
                  AREA_LS_PRC_tolayer     , &
                  PRELS_dev(I)    , & 
                  SNRLS_dev(I)    , &
                  PRN_LS_above    , &
                  PSN_LS_above    , &
                  EVAP_DD_LS_above, &
                  SUBL_DD_LS_above, &
                  REV_LS_dev(I,K) , &
                  RSU_LS_dev(I,K) , &    
                  ACLL_LS_dev(I,K), &
                  ACIL_LS_dev(I,K), &
                  PFL_LS_dev(I,K) , &
                  PFI_LS_dev(I,K) , &
                  VFALLRN         , &
                  VFALLSN         , &
                  FRZ_PP_dev(I,K) , &
                  LSENVFC, LSDDRFC    )

            VFALLSN_LS_dev(I,K) = VFALLSN
            VFALLRN_LS_dev(I,K) = VFALLRN
      
            if (.not.use_autoconv_timescale) then
               if (VFALLSN.NE.0.) then
                  QSN_ALL = QSN_ALL + PFI_LS_dev(I,K)/VFALLSN
               end if
               if (VFALLRN.NE.0.) then
                  QRN_ALL = QRN_ALL + PFL_LS_dev(I,K)/VFALLRN
               end if
            end if

            IF ( (QLW_LS_dev(I,K)+QLW_AN_dev(I,K)) > tiny(0.00) ) THEN
               QTMP3 = 1./(QLW_LS_dev(I,K)+QLW_AN_dev(I,K))
            ELSE
               QTMP3 = 0.0
            END IF
            QLW_LS_dev(I,K) = QLW_LS_dev(I,K) * QTMP1 * QTMP3
            QLW_AN_dev(I,K) = QLW_AN_dev(I,K) * QTMP1 * QTMP3
     
            IF ( (QIW_LS_dev(I,K)+QIW_AN_dev(I,K)) > tiny(0.00) ) THEN
               QTMP3 = 1./(QIW_LS_dev(I,K)+QIW_AN_dev(I,K))
            ELSE
               QTMP3 = 0.0
            END IF
            QIW_LS_dev(I,K) = QIW_LS_dev(I,K) * QTMP2 * QTMP3
            QIW_AN_dev(I,K) = QIW_AN_dev(I,K) * QTMP2 * QTMP3

         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            TH_dev(I,K)  =  TEMP / EXNP_dev(I,K) 

            QRN_ALL = QRN_ALL / (100.*PP_dev(I,K) / (MAPL_RGAS*TEMP ))
            QSN_ALL = QSN_ALL / (100.*PP_dev(I,K) / (MAPL_RGAS*TEMP ))

            RHCLR_dev(I,K) = MIN( CLDFRC_dev(I,K) + ANVFRC_dev(I,K), 1.00 )
            IF ( RHCLR_dev(I,K) < 1.00 ) THEN
               RHCLR_dev(I,K) = ( MIN( Q_dev(I,K)/QSAT(TEMP,PP_dev(I,K)),1.0 ) - RHCLR_dev(I,K) ) / &
                                     (1. - RHCLR_dev(I,K))
               IF ( RHCLR_dev(I,K) < 0.00 ) THEN
                  RHCLR_dev(I,K) = MAPL_UNDEF
               END IF
            ELSE
               RHCLR_dev(I,K) = MAPL_UNDEF
            END IF

            IF (DISABLE_RAD==1) THEN
               RAD_QL_dev(I,K)     = 0.
               RAD_QI_dev(I,K)     = 0.
               RAD_QR_dev(I,K)     = 0.
               RAD_QS_dev(I,K)     = 0.
               RAD_QG_dev(I,K)     = 0.
               RAD_CLDFRC_dev(I,K) = 0.
               CLDREFFL_dev(I,K)   = 0.
               CLDREFFI_dev(I,K)   = 0.
            ELSE
               call RADCOUPLE ( TEMP, PP_dev(I,K), CLDFRC_dev(I,K), ANVFRC_dev(I,K), &
                     Q_dev(I,K), QLW_LS_dev(I,K), QIW_LS_dev(I,K), QLW_AN_dev(I,K), QIW_AN_dev(I,K), QRN_ALL, QSN_ALL, QGR_ALL, NACTL_dev(I,K), NACTI_dev(I,K), & 
                     RAD_QV_dev(I,K), RAD_QL_dev(I,K), RAD_QI_dev(I,K), RAD_QR_dev(I,K), RAD_QS_dev(I,K), RAD_QG_dev(I,K), RAD_CLDFRC_dev(I,K), & 
                     CLDREFFL_dev(I,K), CLDREFFI_dev(I,K), &
                     FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI )
            END IF

            QRN_CU_dev(I,K) = QRN_CU_1D

         end do K_LOOP

#ifndef _CUDA
      end do RUN_LOOP
#else
      end if RUN_LOOP
#endif

   END SUBROUTINE PROGNO_CLOUD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              P R O C E S S   S U B R O U T I N E S                 !!
!!                         * * * * *                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              P R O C E S S   S U B R O U T I N E S                 !!
!!                         * * * * *                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!                                                                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!              P R O C E S S   S U B R O U T I N E S                 !!
!!                         * * * * *                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine pdf_spread (PP,PPsfc,TROPP, &
         minrhcrit, maxrhcrit, turnrhcrit, &
         ALPHA, ALPHT_DIAG)

      real,    intent(in)  :: PP,PPsfc,TROPP
      real,    intent(in)  :: minrhcrit, maxrhcrit, turnrhcrit
      real,    intent(out) :: ALPHA
      real,    intent(out) :: ALPHT_DIAG

      real    :: a1, Al, Au, TURNRHCRIT_UP

#define OLD_RHCRIT
#ifdef OLD_RHCRIT
      ! alpha is the 1/2*width so RH_crit=1.0-alpha

         !  Use Slingo-Ritter (1985) formulation for critical relative humidity
         !  array a1 holds the critical rh, ranges from 0.8 to 1

         a1 = 1.0
         if (pp .le. turnrhcrit) then
            a1 = minrhcrit
         else
            a1 = minrhcrit + (maxrhcrit-minrhcrit)/(19.) * &
                  ((atan( (2.*(pp- turnrhcrit)/(1020.-turnrhcrit)-1.) * &
                  tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
         end if

         a1 = min(a1,1.)

         alpha = 1. - a1

         ALPHA = MIN( ALPHA , 0.25 )  ! restrict RHcrit to > 75% 
#else

       !  Use Slingo-Ritter (1985) formulation for critical relative humidity
           ! lower turn from maxrhcrit
             if (pp .le. TURNRHCRIT) then
                Al = minrhcrit
             else
                Al = minrhcrit + (maxrhcrit-minrhcrit)/(19.) * &
                         ((atan( (2.*(pp-TURNRHCRIT)/(PPsfc-TURNRHCRIT)-1.) * &
                         tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
             endif
           ! upper turn back to maxrhcrit
             TURNRHCRIT_UP = TROPP/100.0 ! Pa to mb
             IF (TURNRHCRIT_UP == MAPL_UNDEF) TURNRHCRIT_UP = 100.
             if (pp .le. TURNRHCRIT_UP) then
                Au = maxrhcrit
             else
                Au = maxrhcrit - (maxrhcrit-minrhcrit)/(19.) * &
                        ((atan( (2.*(pp-TURNRHCRIT_UP)/(TURNRHCRIT-TURNRHCRIT_UP)-1.) * &
                        tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
             endif
           ! combine and limit
             ALPHA = min( 0.25, 1.0 - min(max(Al,Au),1.) ) ! restrict RHcrit to > 75% 

#endif

      ALPHT_DIAG = ALPHA

   end subroutine pdf_spread

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine meltfrz( DT, CNVFRC, TE, QL, QI )

      real, intent(in)    :: DT, CNVFRC
      real, intent(inout) :: TE,QL,QI

      real  :: fQi,dQil

      real  ::  taufrz

      integer :: K

      ! freeze liquid
      if ( TE <= MAPL_TICE ) then
         fQi  = ice_fraction( TE, CNVFRC )
         taufrz = 1000.
         dQil = Ql *(1.0 - EXP( -Dt * fQi / taufrz ) )
         dQil = max(  0., dQil )
         Qi   = Qi + dQil
         Ql   = Ql - dQil
         TE   = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if

      ! melt ice instantly above 0^C
      if ( TE > MAPL_TICE ) then
         dQil = -Qi 
         dQil = min(  0., dQil )
         Qi   = Qi + dQil
         Ql   = Ql - dQil
         TE   = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if

   end subroutine meltfrz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine meltfrz_inst  (     &
         Dt       , &
         CNVFRC   , &
         TE       , &
         QCL      , &
         QAL      , &
         QCI      , &
         QAI      )

      real ,   intent(inout) :: TE,QCL,QCI,QAL,QAI
      real ,   intent(in   ) :: Dt, CNVFRC
      real                   :: fQi,dQil,DQmax, QLTOT, QITOT, FQA
      integer                :: n
      integer, parameter     :: MaxIterations=1
      logical                :: converged
      real                   :: taufrz=450 ! timescale for freezing (seconds)

    ! Make total condensate from LS and CN partitions and store fractions for later
      QITOT=QCI+QAI
      QLTOT=QCL+QAL
      FQA = 0.0
      if (QITOT+QLTOT .gt. 0.0) then
         FQA= (QAI+QAL)/(QITOT+QLTOT)
      end if

      convergence: do n=1,MaxIterations

      ! melt ice using ICE_FRACTION
      fQi = ice_fraction( TE, CNVFRC )
      if ( fQi < 1.0 ) then
         DQmax = (TE-MAPL_TICE)*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil  = QITOT*(1.0-fQi)
         dQil  = min(dQil, DQmax)
         dQil  = max(  0., dQil )
         QLTOT = max(QLTOT+dQil, 0.)
         QITOT = max(QITOT-dQil, 0.)
         TE = TE - (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if

      ! freeze liquid using ICE_FRACTION 
      fQi = ice_fraction( TE, CNVFRC )
      if ( fQi > 0.0 ) then
         DQmax = (MAPL_TICE-TE)*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil  = QLTOT *(1.0 - EXP( -Dt * fQi / taufrz ) )
         dQil  = min(dQil, DQmax)
         dQil  = max(  0., dQil )
         QLTOT = max(QLTOT-dQil, 0.)
         QITOT = max(QITOT+dQil, 0.)
         TE = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP
      end if
      
      if (qltot == 0.0) then
         converged=.true.
      else
         converged = (ABS(dQil/QLTOT) <= 0.01) ! Converge at <=1% of liquid changing phase 
      end if
      if (converged) exit convergence     
 
      end do convergence

    ! Parse into new LS and CN partitions based on previous fractions (FQA)
      QCI = QITOT*(1.0-FQA)
      QAI = QITOT*FQA
      QCL = QLTOT*(1.0-FQA)
      QAL = QLTOT*FQA

   end subroutine meltfrz_inst

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine autocon3( &
         DT       , &
         QC       , &
         QP       , &
         TE       , &
         PL       , &
         KH       , &
         F        , &
         SUNDQV2  , &
         SUNDQV3  , &
         SUNDQT1  , &
         DZET     , &
         VF       , &
         FRACTION_REMOVAL )

      integer, intent(in) :: FRACTION_REMOVAL

      real, intent(in   ) :: DT
      real, intent(in   ) :: TE
      real, intent(in   ) :: PL
      real, intent(in   ) :: KH
      real, intent(inout) :: QC
      real, intent(inout) :: QP
      real, intent(inout) :: F
      real, intent(in   ) :: DZET
      real, intent(inout) :: VF

      real, intent(in   ) :: SUNDQV2, SUNDQV3, SUNDQT1

      real    :: ACF0, ACF
      real    :: C00x, iQCcrx, F2, F3,RATE,dQP,QCm
      real    :: dqfac
      integer :: NS, K

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Precip. conversion from Smith (1990,    !
      !   QJRMS, 116, 435, Eq. 2.29). Similar    ! 
      !   to Del Genios Eq.(10).                 !
      !                                          !
      !   Coalesence term needs to be determined !
      !   through entire column and is done in   !
      !   subroutine "ACCRETE_EVAP_PRECIP"       !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      CALL SUNDQ3_ICE3(TE, SUNDQV2, SUNDQV3, SUNDQT1, F2, F3 )

      C00x  = C_00 * F2 * F3
      iQCcrx = F2 * F3 / LWCRIT

      if ( ( F > 0.) .and. ( QC > 0. ) )then
         QCm = QC/F
      else
         QCm = 0.
      end if


      RATE = C00x * ( 1.0 - EXP( - ( QCm * iQCcrx )**2 ) )

      !!! Make up a fall velocity for liquid precipitation - in analogy to falling ice,
      !!!    think of the fall velocity as the ratio of autoconverted to existing condensate 
      !!!    multiplied by delta z / delta t  
      !!!   (ie, autoconversion/sec is related to residence time in layer)

      VF  =  (DZET/DT) * ( 1.0 - EXP( -RATE * DT ) )

      !! temporary kluge until we can figure a better to make
      !! thicker low clouds ( reuse arrays F2 and F3 )
      F2 = 1.0
      F3 = 1.0 

#ifdef DONT_SKIP_KLUDGES
      ! Implement ramps for gradual change in autoconv
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thicken low high lat clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 775.  .AND. TE .LE.  275. ) then
!!!      F3 = max(-0.016 * PL + 13.4, 0.2)
         F3 = 0.2
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  282. ) then
!!!      F3 = max(0.11 * TE - 30.02, 0.2)
         F3 = 0.2
      end if
      if ( PL .GE. 775.  .AND. PL .LT. 825. .AND. TE .LE.  282. .AND. TE .GT. 275.) then
!!!      F3 = min(max(-0.016*PL + 0.11 * TE - 16.85, 0.2),1.)
         F3 = 0.2
      end if
      if ( PL .GE. 825.  .AND. TE .LE.  275. ) then
         F3 = 0.2
      end if
      if ( PL .LE. 775.  .OR. TE .GT.  282. ) then
         F3 = 1.
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Thin-out low tropical clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if ( PL .GE. 950.  .AND. TE .GE.  285. ) then
         F3 = min(0.2 * TE - 56, 2.)
      end if
      if ( PL .GE. 925.  .AND. TE .GE.  290. ) then
         F3 = min(0.04 * PL - 36., 2.)
      end if
      if ( PL .GE. 925.  .AND. PL .LT. 950. .AND. TE .GT.  285. .AND. TE .LT. 290.) then
         F3 = max(min(0.04*PL + 0.2 * TE - 94., 2.),1.)
      end if
      if ( PL .GE. 950.  .AND. TE .GE.  290. ) then
         F3 = 2.
      end if

      F3   = MAX( F3, 0.1 )

      RATE = F3 * RATE
#endif

      dQP  =  QC*( 1.0 - EXP( -RATE * DT ) )

      dQP  =  MAX( dQP , 0.0 )  ! Protects against floating point problems for tiny RATE

      QC   = QC - dQP  
      QP   = QP + dQP


      SELECT CASE( FRACTION_REMOVAL )

      CASE( 0 )
      ! do NOTHING

      CASE( 1 )
         if ( (QC + dQP) > 0. ) then
            F = QC * F / (QC + dQP )
         end if

      CASE( 2 )
         if ( (QC + dQP) > 0. .and. QC > 0. ) then
            F = F * SQRT( QC / (QC + dQP ) )
         end if

      CASE( 3 )
         if ( (QC + dQP) > 0. .and. QC > 0. ) then
            F = F * ( QC / (QC + dQP ) )**0.333
         end if


      END SELECT

   END SUBROUTINE AUTOCON3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device)   &
#endif
   subroutine PRECIP3(  &
         K,LM         , &
         DT           , & 
         FRLAND       , & 
         RHCR3        , &
         QPl          , &
         QPi          , &
         QCl          , &
         QCi          , &
         TE           , &
         QV           , &
         mass         , &
         imass        , &
         PL           , &
         dZE          , &
         QDDF3        , &
         AA           , &
         BB           , &
         AREA         , &
         RAIN         , & 
         SNOW         , &
         PFl_above    , &
         PFi_above    , &
         EVAP_DD_above, &
         SUBL_DD_above, &
         REVAP_DIAG   , &
         RSUBL_DIAG   , &
         ACRLL_DIAG        , &
         ACRIL_DIAG        , &
         PFL_DIAG     , &
         PFI_DIAG     , &
         VFALLRN      , &
         VFALLSN      , &
         FRZ_DIAG     , &
         ENVFC,DDRFC  )


      integer, intent(in) :: K,LM

      real, intent(in   ) :: DT

      real, intent(inout) :: QV,QPl,QPi,QCl,QCi,TE

      real, intent(in   ) :: mass,imass
      real, intent(in   ) :: PL
      real, intent(in   ) :: AA,BB
      real, intent(in   ) :: RHCR3
      real, intent(in   ) :: dZE
      real, intent(in   ) :: QDDF3
      real, intent(  out) :: RAIN,SNOW
      real, intent(in   ) :: AREA
      real, intent(in   ) :: FRLAND

      real, intent(inout) :: PFl_above, PFi_above
      real, intent(inout) :: EVAP_DD_above, SUBL_DD_above

      real, intent(  out) :: REVAP_DIAG
      real, intent(  out) :: RSUBL_DIAG
      real, intent(  out) :: ACRLL_DIAG,ACRIL_DIAG
      real, intent(  out) :: PFL_DIAG, PFI_DIAG
      real, intent(inout) :: FRZ_DIAG
      real, intent(  out) :: VFALLSN, VFALLRN

      real, intent(in   ) :: ENVFC,DDRFC


      real :: PFi,PFl,QS,dQS,ENVFRAC
      real :: TKo,QKo,QSTKo,DQSTKo,RH_BOX,T_ED,QPlKo,QPiKo
      real :: Ifactor,RAINRAT0,SNOWRAT0
      real :: FALLRN,FALLSN,VEsn,VErn,NRAIN,NSNOW,Efactor

      real :: TinLAYERrn,DIAMrn,DROPRAD
      real :: TinLAYERsn,DIAMsn,FLAKRAD

      real :: EVAP,SUBL,ACCR,MLTFRZ,EVAPx,SUBLx
      real :: EVAP_DD,SUBL_DD
      real :: LANDSEAF, TC, MAXMLT, iDT

      real, parameter :: TRMV_L  = 1.0         ! m/s
      real, parameter :: TAU_FRZ = 5000.0      ! sec       
      real, parameter :: FRZ_TAU = 1.0/TAU_FRZ ! sec^-1
      real, parameter :: MELT_T  = 5.0         ! degrees C
      real, parameter :: LFBYCP  = MAPL_ALHF/MAPL_CP
      real, parameter :: CPBYLF  = 1.0/LFBYCP
      real, parameter :: B_SUB   = 1.00
         
         
      integer :: NS, itr,L

      logical, parameter :: taneff = .false.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! fraction of precip falling through "environment" vs
      ! through cloud

      if(taneff) then
         envfrac = 1.00

         if (pl .le. 600.) then
            envfrac = 0.25
         else
            envfrac = 0.25 + (1.-0.25)/(19.) *                    &
                  ((atan( (2.*(pl-600.)/(900.-600.)-1.) *       &
                  tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
         end if

         envfrac = min(envfrac,1.)
      else
         ENVFRAC = ENVFC
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( AREA > 0. ) THEN
         Ifactor = 1./ ( AREA )
      ELSE
         Ifactor = 1.00
      END if

      Ifactor = MAX( Ifactor, 1.) !

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !   Start at top of precip column:
      !  
      !               a) Accrete                   
      !               b) Evaporate/Sublimate  
      !               c) Rain/Snow-out to next level down 
      !               d) return to (a)
      !
      !   ....................................................................
      !           
      !  Accretion formulated according to Smith (1990, Q.J.R.M.S., 116, 435
      !  Eq. 2.29)
      !  
      !  Evaporation (ibid. Eq. 2.32)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      iDT = 1.0/DT

      !!! INITIALIZE DIAGNOSTIC ARRAYS !!!!!!!!!!!!!!!!!!!!!

      PFL_DIAG   =  0.
      PFI_DIAG   =  0.
      ACRIL_DIAG =  0.
      ACRLL_DIAG =  0.
      REVAP_DIAG =  0.
      RSUBL_DIAG =  0.

      !!!!!!!!!!!!!! UPDATE SATURATED HUMIDITY  !!!!!!!!!!!!!

      dQS = DQSAT( TE, PL, QSAT = QS )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (K == 1) THEN
         PFl = QPl*MASS
         PFi = QPi*MASS

         EVAP_DD = 0.
         SUBL_DD = 0.

         VFALLRN = 0.0
         VFALLSN = 0.0
      ELSE 
         QPl = QPl + PFl_above * iMASS
         PFl = 0.00

         QPi = QPi + PFi_above * iMASS
         PFi = 0.00

         !! Accretion of liquid condensate by falling rain or ice/snow
         !! Liquid freezes when accreted by snow

         IF(QCL > 0.0) THEN
            IF(QPl > 0.0) THEN
               ACCR = min(C_ACC*(QPl*MASS)*QCl, QCl)
               QPl  = QPl + ACCR
               QCl  = QCl - ACCR

               ACRLL_DIAG = ACCR * iDT
            END IF

            IF(QPi > 0.0) THEN
               ACCR = min(C_ACC*(QPi*MASS)*QCl, QCl)
               QPi  = QPi + ACCR
               QCl  = QCl - ACCR
               TE   = TE + LFBYCP*ACCR

               ACRIL_DIAG = ACCR * iDT
            END IF
         END IF

         RAINRAT0 = Ifactor*QPl*MASS*iDT
         SNOWRAT0 = Ifactor*QPi*MASS*iDT

         call MARSHPALMQ2(RAINRAT0,PL,DIAMrn,NRAIN,FALLrn,VErn)
         call MARSHPALMQ2(SNOWRAT0,PL,DIAMsn,NSNOW,FALLsn,VEsn)

         VFALLRN = FALLrn
         VFALLSN = FALLsn

         TinLAYERrn = dZE / ( FALLrn+0.01 )
         TinLAYERsn = dZE / ( FALLsn+0.01 )


         !*****************************************
         !  Melting of Frozen precipitation      
         !*****************************************

         TC = TE - MAPL_TICE

         IF( QPi > 0.0 .AND. TC > 0.0) THEN

            MAXMLT = min(QPi, TC*CPBYLF)

            IF ( K < LM-3 .and. TC <= MELT_T) THEN
               MLTFRZ = min(TinLAYERsn*QPi*TC*FRZ_TAU, MAXMLT)
            else
               MLTFRZ = MAXMLT
            END IF

            TE       = TE  - LFBYCP*MLTFRZ
            QPl      = QPl + MLTFRZ
            QPi      = QPi - MLTFRZ
            FRZ_DIAG = FRZ_DIAG - MLTFRZ * iDT

         END IF

         !*****************************************
         !  Freezing of rain
         !*****************************************

         IF ( QPl > 0.0 .AND.  TC <= 0.0) THEN

            MLTFRZ = min(QPl,-TC*CPBYLF)
            TE     = TE  + LFBYCP*MLTFRZ
            QPi    = QPi + MLTFRZ
            QPl    = QPl - MLTFRZ
            FRZ_DIAG = FRZ_DIAG + MLTFRZ * iDT

         END IF


         ! ******************************************
         !   In the exp below, evaporation time 
         !   scale is determined "microphysically"
         !   from temp, press, and drop size. In this
         !   context C_EV becomes a dimensionless 
         !   fudge-fraction.
         !   Also remember that these microphysics 
         !   are still only for liquid.
         ! ******************************************

         QKo   = QV
         TKo   = TE
         QPlKo = QPl
         QPiKo = QPi

         do itr = 1,3

            DQSTKo = dQS
            QSTKo  = QS  + DQSTKo * ( TKo - TE )
            QSTKo  = MAX( QSTKo , 1.0e-7 )

            RH_BOX = QKo/QSTKo

            QKo   = QV
            TKo   = TE

            IF ( RH_BOX < RHCR3 ) THEN
               Efactor =  RHO_W * ( AA + BB )    / (RHCR3 - RH_BOX )
            else
               Efactor = 9.99e9
            end if

            if ( FRLAND < 0.1 ) then
               LANDSEAF = 0.5  ! Over Ocean
            else
               LANDSEAF = 0.5  ! Over Land
            end if

            LANDSEAF = 1.00

            !!!!! RAin falling !!!!!!!!!!!!!!!!!!!!!!!
            if ( ( RH_BOX < RHCR3 ) .AND. ( DIAMrn > 0.00 ) .AND. &
                  ( PL > 100. ) .AND. ( PL < REVAP_OFF_P ) ) then
               DROPRAD=0.5*DIAMrn
               T_ED =  Efactor * DROPRAD**2 
               T_ED =  T_ED * ( 1.0 + DQSTKo*MAPL_ALHL/MAPL_CP )
               EVAP =  QPl*(1.0 - EXP( -C_EV_R * VErn * LANDSEAF * ENVFRAC * TinLAYERrn / T_ED ) )
            ELSE
               EVAP = 0.0
            END if

            !!!!! Snow falling !!!!!!!!!!!!!!!!!!!!!!!
            if ( ( RH_BOX < RHCR3 ) .AND. ( DIAMsn > 0.00 ) .AND. &
                  ( PL > 100. ) .AND. ( PL < REVAP_OFF_P ) ) then
               FLAKRAD=0.5*DIAMsn
               T_ED =  Efactor * FLAKRAD**2   
               T_ED =  T_ED * ( 1.0 + DQSTKo*MAPL_ALHS/MAPL_CP )
               SUBL =  QPi*(1.0 - EXP( -C_EV_S * VEsn * LANDSEAF * ENVFRAC * TinLAYERsn / T_ED ) )
            ELSE
               SUBL = 0.0
            END IF

            if (itr == 1) then 
               EVAPx  = EVAP
               SUBLx  = SUBL
            else
               EVAP   = (EVAP+EVAPx) /2.0
               SUBL   = (SUBL+SUBLx) /2.0
            endif

            QKo=QV + EVAP + SUBL
            TKo=TE - EVAP * MAPL_ALHL / MAPL_CP - SUBL * MAPL_ALHS / MAPL_CP

         enddo
      
         QPi  = QPi - SUBL
         QPl  = QPl - EVAP

         !! Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
         EVAP_DD = EVAP_DD_above + DDRFC*EVAP*MASS 
         EVAP    = EVAP          - DDRFC*EVAP
         SUBL_DD = SUBL_DD_above + DDRFC*SUBL*MASS 
         SUBL    = SUBL          - DDRFC*SUBL
         ! -----

         QV   = QV  + EVAP + SUBL
         TE   = TE  - EVAP * MAPL_ALHL / MAPL_CP - SUBL * MAPL_ALHS / MAPL_CP

         REVAP_DIAG = EVAP / DT
         RSUBL_DIAG = SUBL / DT

         PFl  = QPl*MASS
         PFi  = QPi*MASS

         PFL_DIAG =  PFl/DT
         PFI_DIAG =  PFi/DT
      end if

      ! QDDF3 (<= QDDF3_dev) is calculated on the CPU in order to avoid
      ! the reverse loop on GPUs and thus save local memory use.
      EVAP = QDDF3*EVAP_DD/MASS
      SUBL = QDDF3*SUBL_DD/MASS
      QV   = QV  + EVAP + SUBL
      TE   = TE  - EVAP * MAPL_ALHL / MAPL_CP - SUBL * MAPL_ALHS / MAPL_CP
      REVAP_DIAG = REVAP_DIAG + EVAP / DT
      RSUBL_DIAG = RSUBL_DIAG + SUBL / DT

      IF (K == LM) THEN
         RAIN  = PFl/DT
         SNOW  = PFi/DT
      END IF

      QPi = 0.
      QPl = 0.

      PFl_above = PFl
      PFi_above = Pfi

      EVAP_DD_above = EVAP_DD
      SUBL_DD_above = SUBL_DD

   end subroutine precip3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine ICEFALL( QI, DZ, QP, VF, F, DT, FRACTION_REMOVAL )

      integer, intent(in) :: FRACTION_REMOVAL

      real, intent(inout) :: QI
      real, intent(in   ) :: DZ
      real, intent(inout) :: QP
      real, intent(in   ) :: VF
      real, intent(inout) :: F
      real, intent(in   ) :: DT

      real :: QIxP

      QIxP = QI * ( VF * DT / DZ ) 
      QIxP = MIN( QIxP , QI )

      QIxP = MAX( QIxP, 0.0 ) ! protects against precision problem

      QP = QP + QIxP
      QI = QI - QIxP

      SELECT CASE( FRACTION_REMOVAL )

      CASE( 0 )
         ! do NOTHING

      CASE( 1 )
         if ( (QI + QIxP) > 0. ) then
            F = QI * F / (QI + QIxP )
         end if

      CASE( 2 )
         if ( (QI + QIxP) > 0. .and. QI > 0. ) then
            F = F * SQRT( QI / (QI + QIxP ) )
         end if

      CASE( 3 )
         if ( (QI + QIxP) > 0. .and. QI > 0. ) then
            F = F * ( QI / (QI + QIxP ) )**0.333
         end if

      END SELECT

   end subroutine ICEFALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine SETTLE_VEL( QI, PL, TE, F, KH, VF, LARGESCALE, ANVIL, CNVFRC )

      real, intent(in   ) :: TE
      real, intent(in   ) :: QI, F, PL
      real, intent(in   ) :: KH
      real, intent(out  ) :: VF

      real, intent(in) :: ANVIL, LARGESCALE, CNVFRC
      
      real :: RHO, XIm, LXIm, VF_A, VF_L

      real :: PFAC

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Uses Eq. 1 Lawrence and Crutzen (1998, Tellus 50B, 263-289) 
      ! Except midlat form is taken to be for LS cloud, and tropical
      ! form is taken to be for ANVIL cloud
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      RHO = 1000.*100.*PL/(MAPL_RGAS*TE)  ! 1000 TAKES TO g m^-3 ; 100 takes mb TO Pa

      if ( ( F > 0.) .and. ( QI > 0. ) ) then
         XIm = (QI/F)*RHO
      else
         XIm = 0.
      end if

      if ( XIm > 0.) then
         LXIm = ALOG10(XIm)
      else
         LXIm = 0.0
      end if

    ! Pressure factor
    ! Reduce/increase fall speeds for high/low pressure (NOT in LC98!!! ) 
    ! Assume unmodified they represent situation at 100 mb
       PFAC = SIN( 0.5*MAPL_PI*MIN(1.0,100./PL))

    ! Convective anvil
       VF_A = (128.6 + 53.2*LXIm + 5.5*LXIm**2) * MIN(ANVIL,PFAC)

    ! Mid-latitude cirrus
       VF_L = (109.0*(XIm**0.16)) * MIN(LARGESCALE,PFAC)
 
    ! Combine the two and convert from cm/s to m/s
       VF = 0.01 * (CNVFRC*VF_A + (1.0-CNVFRC)*VF_L)

   end SUBROUTINE SETTLE_VEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine MARSHPALMQ2(RAIN,PR,DIAM3,NTOTAL,W,VE)

      real, intent(in ) :: RAIN,PR     ! in kg m**-2 s**-1, mbar
      real, intent(out) :: DIAM3,NTOTAL,W,VE

      real :: RAIN_DAY,LAMBDA,A,B,SLOPR,DIAM1

      real, parameter  :: N0   = 0.08  ! # cm**-3

      INTEGER :: IQD

      real :: RX(8) , D3X(8)

      ! Marshall-Palmer sizes at different rain-rates: avg(D^3)

      RX = (/ 0.   , 5.   , 20.  , 80.  , 320. , 1280., 4*1280., 16*1280. /)  ! rain per in mm/day
      D3X= (/ 0.019, 0.032, 0.043, 0.057, 0.076, 0.102, 0.137  ,  0.183   /)

      RAIN_DAY = RAIN * 3600. *24.

      IF ( RAIN_DAY <= 0.00 ) THEN
         DIAM1 = 0.00
         DIAM3 = 0.00
         NTOTAL= 0.00
         W     = 0.00
      END IF

      DO IQD = 1,7
         IF ( (RAIN_DAY <= RX(IQD+1)) .AND. (RAIN_DAY > RX(IQD) ) ) THEN
            SLOPR =( D3X(IQD+1)-D3X(IQD) ) / ( RX(IQD+1)-RX(IQD) )
            DIAM3 = D3X(IQD) + (RAIN_DAY-RX(IQD))*SLOPR
         END IF
      END DO

      IF ( RAIN_DAY >= RX(8) ) THEN
         DIAM3=D3X(8)
      END IF

      NTOTAL = 0.019*DIAM3

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DIAM3 = 0.664 * DIAM3  !!  DRYING/EVAP SHOULD PROBABLY GO AS          !!
      !!  D_1.5 == <<D^(3/2)>>^(2/3) NOT AS          !!
      !!  D_3   == <<D^3>>^(1/3)                     !!
      !!  RATIO D_1.5/D_3 =~ 0.66  (JTB 10/17/2002)  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      W      = (2483.8 * DIAM3 + 80.)*SQRT(1000./PR)
      !VE     = 1.0  + 28.0*DIAM3 
      VE     = MAX( 0.99*W/100. , 1.000 )

      DIAM1  = 3.0*DIAM3
      !  Change back to MKS units

      DIAM1  = DIAM1/100.
      DIAM3  = DIAM3/100.
      W      = W/100.
      NTOTAL = NTOTAL*1.0e6

   end subroutine MARSHPALMQ2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine MICRO_AA_BB_3(TEMP,PR,Q_SAT,AA,BB)

      real, intent(in ) :: TEMP,Q_SAT
      real, intent(in ) :: PR
      real, intent(out) :: AA,BB

      real :: E_SAT

      real, parameter  :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
      real, parameter  :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
      real, parameter  :: DIFFU   =  2.2e-5        ! m**2 s**-1

      E_SAT = 100.* PR * Q_SAT /( (EPSILON) + (1.0-(EPSILON))*Q_SAT )  ! (100 converts from mbar to Pa)
   
      AA  = ( GET_ALHX3(TEMP)**2 ) / ( K_COND*MAPL_RVAP*(TEMP**2) )
      ! AA  = ( MAPL_ALHL**2 ) / ( K_COND*MAPL_RVAP*(TEMP**2) )

      BB  = MAPL_RVAP*TEMP / ( DIFFU*(1000./PR)*E_SAT )

   end subroutine MICRO_AA_BB_3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   function GET_ALHX3(T) RESULT(ALHX3)

      real, intent(in) :: T
      real :: ALHX3

      real :: T_X

      T_X = T_ICE_MAX

      if ( T < T_ICE_ALL ) then
         ALHX3=MAPL_ALHS
      end if

      if ( T > T_X ) then
         ALHX3=MAPL_ALHL
      end if

      if ( (T <= T_X) .and. (T >= T_ICE_ALL) ) then
         ALHX3 = MAPL_ALHS + (MAPL_ALHL-MAPL_ALHS)*( T - T_ICE_ALL ) /( T_X - T_ICE_ALL )
      end if

   end function GET_ALHX3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   SUBROUTINE SUNDQ3_ICE3(TEMP,RATE2,RATE3,TE1, F2, F3)

      REAL, INTENT( IN) :: RATE2,RATE3,TE1

      REAL, INTENT( IN) :: TEMP
      REAL, INTENT(OUT) :: F2, F3

      REAL :: XX, YY,TE0,TE2,JUMP1  !,RATE2,RATE3,TE1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Ice - phase treatment totally invented
      !!  Sharp increase in autoconversion in range
      !!  ~~TE1 K ~< T < TE0 K .
      !!  (JTB, 3/25/2003)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      TE0=MAPL_TICE
      !TE1=263.
      TE2=200.
      !RATE2=  10.
      JUMP1=  (RATE2-1.0) / ( ( TE0-TE1 )**0.333 ) 
      !RATE3=  25.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Ice - phase treatment  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      IF ( TEMP .GE. TE0 ) THEN
         F2   = 1.0
         F3   = 1.0
      END IF
      IF ( ( TEMP .GE. TE1 ) .AND. ( TEMP .LT. TE0 ) ) THEN
         F2   = 1.0 + JUMP1 * (( TE0 - TEMP )**0.3333)
         F3   = 1.0
      END IF
      IF ( TEMP .LT. TE1 ) THEN
         F2   = RATE2 + (RATE3-RATE2)*(TE1-TEMP)/(TE1-TE2)
         F3   = 1.0
      END IF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      F2 = MIN(F2,27.0)


   end  subroutine sundq3_ice3

end module cloudnew
