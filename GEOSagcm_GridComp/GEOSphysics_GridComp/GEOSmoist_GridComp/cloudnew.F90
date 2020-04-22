! $Id$
! $Name$

module cloudnew

#ifndef _CUDA
   use GEOS_UtilsMod,     only: QSAT=>GEOS_Qsat, DQSAT=>GEOS_DQsat, &
         QSATLQ=>GEOS_QsatLQU, QSATIC=>GEOS_QsatICE
   use CLDPARAMS
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

   use partition_pdf

   implicit none

#ifndef _CUDA
   private

   PUBLIC PROGNO_CLOUD
   PUBLIC ICE_FRACTION
   PUBLIC T_CLOUD_CTL
   PUBLIC fix_up_clouds
#endif

   type T_CLOUD_CTL
      real  :: SCLMFDFR
      real  :: RSUB_RADIUS
   end type T_CLOUD_CTL

#ifdef _CUDA

   ! Inputs
   ! ------

   real, allocatable, dimension(:,:), device :: PP_dev
   real, allocatable, dimension(:,:), device :: EXNP_dev
   real, allocatable, dimension(:,:), device :: PPE_dev
   real, allocatable, dimension(:,:), device :: KH_dev
   real, allocatable, dimension(:  ), device :: DTS_dev
   real, allocatable, dimension(:  ), device :: SNOMAS_dev
   real, allocatable, dimension(:  ), device :: FRLANDICE_dev
   real, allocatable, dimension(:  ), device :: FRLAND_dev
   real, allocatable, dimension(:,:), device :: RMFDTR_dev
   real, allocatable, dimension(:,:), device :: QLWDTR_dev
   real, allocatable, dimension(:,:), device :: U_dev
   real, allocatable, dimension(:,:), device :: V_dev
   real, allocatable, dimension(:,:), device :: QST3_dev
   real, allocatable, dimension(:,:), device :: DZET_dev
   real, allocatable, dimension(:,:), device :: QDDF3_dev
   real, allocatable, dimension(:  ), device :: TEMPOR_dev
   real, allocatable, dimension(:  ), device :: CNV_FRACTION_dev
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
   real, allocatable, dimension(:,:), device :: QPLS_dev
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
   real,    constant :: CLDVOL2FRC
   real,    constant :: RHSUP_ICE
   real,    constant :: SHR_EVAP_FAC
   real,    constant :: MIN_CLD_WATER
   real,    constant :: CLD_EVP_EFF
   integer, constant :: NSMAX
   real,    constant :: LS_SDQV2
   real,    constant :: LS_SDQV3
   real,    constant :: LS_SDQVT1
   real,    constant :: ANV_SDQV2
   real,    constant :: ANV_SDQV3
   real,    constant :: ANV_SDQVT1
   real,    constant :: ANV_TO_LS
   real,    constant :: N_WARM
   real,    constant :: N_ICE
   real,    constant :: N_ANVIL
   real,    constant :: N_PBL
   integer, constant :: DISABLE_RAD
   integer, constant :: ICE_SETTLE
   real,    constant :: ANV_ICEFALL_C
   real,    constant :: LS_ICEFALL_C
   real,    constant :: REVAP_OFF_P
   real,    constant :: CNVENVFC
   real,    constant :: ANVENVFC
   real,    constant :: SCENVFC
   real,    constant :: LSENVFC
   real,    constant :: WRHODEP
   real,    constant :: T_ICE_ALL
   real,    constant :: CNVICEPARAM
   integer, constant :: ICEFRPWR
   real,    constant :: CNVDDRFC
   real,    constant :: ANVDDRFC
   real,    constant :: LSDDRFC
   integer, constant :: TANHRHCRIT
   real,             :: MINRHCRIT
   real,    constant :: MINRHCRIT_I
   real,    constant :: MAXRHCRIT
   real,    constant :: TURNRHCRIT
   real,    constant :: MAXRHCRITLAND
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
   real,    constant :: CFPBL_EXP
   integer, constant :: PDFFLAG

   ! Parameters for Internal DQSAT
   ! -----------------------------

   real, parameter :: ESFAC            = MAPL_H2OMW/MAPL_AIRMW
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
   real    :: CLDVOL2FRC
   real    :: RHSUP_ICE
   real    :: SHR_EVAP_FAC
   real    :: MIN_CLD_WATER
   real    :: CLD_EVP_EFF
   integer :: NSMAX
   real    :: LS_SDQV2
   real    :: LS_SDQV3
   real    :: LS_SDQVT1
   real    :: ANV_SDQV2
   real    :: ANV_SDQV3
   real    :: ANV_SDQVT1
   real    :: ANV_TO_LS
   real    :: N_WARM
   real    :: N_ICE
   real    :: N_ANVIL
   real    :: N_PBL
   integer :: DISABLE_RAD
   integer :: ICE_SETTLE
   real    :: ANV_ICEFALL_C
   real    :: LS_ICEFALL_C
   real    :: REVAP_OFF_P
   real    :: CNVENVFC
   real    :: ANVENVFC
   real    :: SCENVFC
   real    :: LSENVFC
   real    :: WRHODEP
   real    :: T_ICE_ALL
   real    :: CNVICEPARAM
   integer :: ICEFRPWR
   real    :: CNVDDRFC
   real    :: ANVDDRFC
   real    :: SCDDRFC
   real    :: LSDDRFC
   integer :: tanhrhcrit
   real    :: minrhcrit, minrhcrit_i
   real    :: maxrhcrit
   real    :: turnrhcrit
   real    :: MIN_RI, MAX_RI, FAC_RI, MIN_RL, MAX_RL, FAC_RL, CFPBL_EXP
   integer :: FR_LS_WAT, FR_LS_ICE, FR_AN_WAT, FR_AN_ICE
   real    :: maxrhcritland
   integer :: pdfflag
   real    :: thl2tune
   real    :: qw2tune
   real    :: qwthl2tune
   real    :: fac_cond, fac_fus
#endif

   real, parameter :: T_ICE_MAX    = MAPL_TICE-10.0
   real, parameter :: RHO_W        = 1.0e3      ! Density of liquid water in kg/m^3
   real, parameter :: MIN_CLD_FRAC = 1.0e-8
   real, parameter :: ZVIR = MAPL_RVAP/MAPL_RGAS - 1.
   real, parameter :: GORD = MAPL_GRAV/MAPL_RGAS
   real, parameter :: GFAC = 1.e5/MAPL_GRAV 
   real, parameter :: R_AIR = 3.47e-3 !m3 Pa kg-1K-1

  ! ICE_FRACTION constants
        ! In anvil/convective clouds
   real, parameter :: aT_ICE_ALL = 245.16
   real, parameter :: aT_ICE_MAX = 261.16
   real, parameter :: aICEFRPWR  = 2.0
        ! Over snow/ice
   real, parameter :: iT_ICE_ALL = MAPL_TICE-40.0
   real, parameter :: iT_ICE_MAX = MAPL_TICE
   real, parameter :: iICEFRPWR  = 4.0
        ! Over Land
   real, parameter :: lT_ICE_ALL = 239.16
   real, parameter :: lT_ICE_MAX = 261.16
   real, parameter :: lICEFRPWR  = 2.0
        ! Over Oceans
   real, parameter :: oT_ICE_ALL = 238.16
   real, parameter :: oT_ICE_MAX = 263.16
   real, parameter :: oICEFRPWR  = 4.0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains 

! GPU The GPU main routine call is smaller due to CUDA limit on
!     number of arguments permitted in call. Most inputs and outputs
!     are USE-associated in the GridComp

#ifdef _CUDA
   attributes(global) subroutine progno_cloud(IRUN,LM,DT,SCLMFDFR)
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
         SNOMAS_dev       , &
         FRLANDICE_dev    , &
         FRLAND_dev       , &
         KH_dev           , &
         ISOTROPY_dev     , &
         mf_frc_dev       , &
         au_dev           , &
         hle_dev          , &
         qte_dev          , &
         hl2u_dev         , &
         qt2u_dev         , &
         hlqtu_dev        , &    
         wqt_dev          , &
         whl_dev          , &
         qt2_dev          , &
         hl2_dev          , &
         hlqt_dev         , &
         w2_dev           , &
         w3_dev           , &
         qt3_dev          , &
         TKE_dev          , &
         DTS_dev          , &
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
         QPLS_dev         , &
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
         CLDPARAMS        , &
         SCLMFDFR         , &
         QST3_dev         , &
         DZET_dev         , &
         QDDF3_dev        , &
         CNV_FRACTION_dev , &
         TROPP_dev        , &
         A_cloud          , &
         B_cloud          , &
         qsat             , &
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
         PDF_A_dev, PDF_SIGW_dev, PDF_W1_dev, PDF_W2_dev, & 
         PDF_SIGTH1_dev, PDF_SIGTH2_dev, PDF_TH1_dev, PDF_TH2_dev, &
         PDF_SIGQT1_dev, PDF_SIGQT2_dev, PDF_QT1_dev, PDF_QT2_dev, &
         PDF_RQTTH_dev, WTHV2_dev, wql_dev, SKEW_QT_dev, &
         TEMPOR_dev, DOSHLW, &
         NACTL_dev,    &
         NACTI_dev,    &
         CONVPAR_OPTION)
#endif

      implicit none

#ifdef _CUDA
      integer, intent(in   ), value :: IRUN
      integer, intent(in   ), value :: LM
      real   , intent(in   ), value :: DT
      real   , intent(in   ), value :: SCLMFDFR   ! CLOUD_CTL%SCLMFDFR
#else
      type (CLDPARAM_TYPE), intent(in)          :: CLDPARAMS

      integer, intent(in   )                    :: IRUN ! IM*JM
      integer, intent(in   )                    :: LM   ! LM
      real, intent(in   )                       :: DT   ! DT_MOIST
      real, intent(in   ), dimension(IRUN)      :: LATS_dev    ! LATS
      real, intent(in   ), dimension(IRUN,  LM) :: PP_dev      ! PLO
      real, intent(in   ), dimension(IRUN,  LM) :: ZZ_dev      ! ZLO
      real, intent(in   ), dimension(IRUN,0:LM) :: PPE_dev     ! CNV_PLE
      real, intent(in   ), dimension(IRUN,  LM) :: EXNP_dev    ! PK
      real, intent(in   ), dimension(IRUN     ) :: SNOMAS_dev  ! SNOMAS
      real, intent(in   ), dimension(IRUN     ) :: FRLANDICE_dev  ! FRLANDICE
      real, intent(in   ), dimension(IRUN     ) :: FRLAND_dev  ! FRLAND
      real, intent(in   ), dimension(IRUN,0:LM) :: KH_dev      ! KH
      real, intent(in   ), dimension(IRUN,  LM) :: ISOTROPY_dev! ISOTROPY
      real, intent(in   ), dimension(IRUN,  LM) :: mf_frc_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: au_dev   !
      real, intent(in   ), dimension(IRUN,  LM) :: hle_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: qte_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hl2u_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: qt2u_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hlqtu_dev !
      real, intent(in   ), dimension(IRUN,  LM) :: wqt_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: whl_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: qt2_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hl2_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: hlqt_dev !
      real, intent(in   ), dimension(IRUN,  LM) :: w2_dev   !
      real, intent(in   ), dimension(IRUN,  LM) :: w3_dev   !
      real, intent(in   ), dimension(IRUN,  LM) :: qt3_dev  !
      real, intent(in   ), dimension(IRUN,  LM) :: tke_dev  !
      real, intent(in   ), dimension(IRUN     ) :: DTS_dev     ! DTS
      real, intent(in   ), dimension(IRUN,  LM) :: RMFDTR_dev  ! CNV_MFD
      real, intent(in   ), dimension(IRUN,  LM) :: QLWDTR_dev  ! CNV_DQLDT
      real, intent(inout), dimension(IRUN,  LM) :: QRN_CU_dev  ! CNV_PRC3 IS THIS INTENT IN?
      real, intent(inout), dimension(IRUN,  LM) :: CNV_UPDFRC_dev ! CNV_UPDF
      real, intent(in   ), dimension(IRUN,  LM) :: SC_RMFDTR_dev  ! MFDSHLW
      real, intent(in   ), dimension(IRUN,  LM) :: SC_QLWDTR_dev  ! DQLDTSHLW
      real, intent(in   ), dimension(IRUN,  LM) :: SC_QIWDTR_dev  ! DQIDTSHLW
      real, intent(inout), dimension(IRUN,  LM) :: QRN_SC_dev     ! SHLW_PRC3
      real, intent(inout), dimension(IRUN,  LM) :: QSN_SC_dev     ! SHLW_SNO3
      real, intent(inout), dimension(IRUN,  LM) :: SC_UPDFRC_dev  ! UPDFSHLW
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
      real, intent(  out), dimension(IRUN,  LM) :: QPLS_dev ! QPLS
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
      real, intent(in   )                       :: SCLMFDFR   ! CLOUD_CTL%SCLMFDFR
      real, intent(in   ), dimension(IRUN,  LM) :: QST3_dev   ! QST3
      real, intent(in   ), dimension(IRUN,  LM) :: DZET_dev   ! DZET
      real, intent(in   ), dimension(IRUN,  LM) :: QDDF3_dev  ! QDDF3
      real, intent(in   ), dimension(IRUN)      :: CNV_FRACTION_dev   ! CNV_FRACTION
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
      real, intent(  out), dimension(IRUN,  LM) :: PDF_SIGW_dev
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
      real, intent(  out), dimension(IRUN,  LM) :: WTHV2_dev
      real, intent(  out), dimension(IRUN,  LM) :: wql_dev
      real, intent(inout), dimension(IRUN,  LM) :: SKEW_QT_dev
      real, intent(in   ), dimension(IRUN,  LM) :: NACTL_dev  ! NACTL
      real, intent(in   ), dimension(IRUN,  LM) :: NACTI_dev  ! NACTI
      real, intent(  out), dimension(IRUN,  LM) :: A_cloud
      real, intent(  out), dimension(IRUN,  LM) :: B_cloud
      real, intent(  out), dimension(IRUN,  LM) :: qsat

!!$      real, intent(  out), dimension(IRUN,  LM) :: LIQANMOVE_dev  ! LIQANMOVE
!!$      real, intent(  out), dimension(IRUN,  LM) :: ICEANMOVE_dev  ! ICEANMOVE
!!$      real, intent(  out), dimension(IRUN,  LM) :: DANCLD_dev     ! DANCLD
!!$      real, intent(  out), dimension(IRUN,  LM) :: DLSCLD_dev     ! DLSCLD
!!$      real, intent(  out), dimension(IRUN,  LM) :: CURAINMOVE_dev ! CURAINMOVE
!!$      real, intent(  out), dimension(IRUN,  LM) :: CUSNOWMOVE_dev ! CUSNOWMOVE

      real, intent(in   ), dimension(IRUN     ) :: TEMPOR_dev  ! TEMPOR
      INTEGER, INTENT(IN) :: DOSHLW
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
      real :: QRN_ALL, QSN_ALL
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

      real :: TROPICAL, EXTRATROPICAL

      real :: LSPDFLIQNEW, LSPDFICENEW, LSPDFFRACNEW

      real, dimension(LM) :: hl_sec, qt_sec, hlqt_sec, wqt_sec, whl_sec, &
                             hl, total_water, w3var, w2var, &
                             hlsec, qtsec, hlqtsec, wqtsec, whlsec
      real :: wrk1, wrk2, wrk3, sm

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
         CLDVOL2FRC    = CLDPARAMS%VOL_TO_FRAC
         RHSUP_ICE     = CLDPARAMS%SUPERSAT
         SHR_EVAP_FAC  = CLDPARAMS%SHEAR_EVAP_FAC
         MIN_CLD_WATER = CLDPARAMS%MIN_ALLOW_CCW
         CLD_EVP_EFF   = CLDPARAMS%CCW_EVAP_EFF
         NSMAX         = INT( CLDPARAMS%NSUB_AUTOCONV  )
         LS_SDQV2      = CLDPARAMS%LS_SUND_INTER
         LS_SDQV3      = CLDPARAMS%LS_SUND_COLD
         LS_SDQVT1     = CLDPARAMS%LS_SUND_TEMP1
         ANV_SDQV2     = CLDPARAMS%ANV_SUND_INTER
         ANV_SDQV3     = CLDPARAMS%ANV_SUND_COLD
         ANV_SDQVT1    = CLDPARAMS%ANV_SUND_TEMP1
         ANV_TO_LS     = CLDPARAMS%ANV_TO_LS_TIME
         N_WARM        = CLDPARAMS%NCCN_WARM
         N_ICE         = CLDPARAMS%NCCN_ICE
         N_ANVIL       = CLDPARAMS%NCCN_ANVIL
         N_PBL         = CLDPARAMS%NCCN_PBL
         DISABLE_RAD   = INT( CLDPARAMS%DISABLE_RAD )
         ICE_SETTLE    = NINT( CLDPARAMS%ICE_SETTLE )
         ANV_ICEFALL_C = CLDPARAMS%ANV_ICEFALL
         LS_ICEFALL_C  = CLDPARAMS%LS_ICEFALL
         REVAP_OFF_P   = CLDPARAMS%REVAP_OFF_P
         CNVENVFC      = CLDPARAMS%CNV_ENVF
         ANVENVFC      = CLDPARAMS%ANV_ENVF
         SCENVFC       = CLDPARAMS%SC_ENVF
         LSENVFC       = CLDPARAMS%LS_ENVF
         WRHODEP       = CLDPARAMS%WRHODEP
         T_ICE_ALL     = CLDPARAMS%ICE_RAMP + T_ICE_MAX
         CNVICEPARAM   = CLDPARAMS%CNV_ICEPARAM
         ICEFRPWR      = INT( CLDPARAMS%CNV_ICEFRPWR )
         CNVDDRFC      = CLDPARAMS%CNV_DDRF
         ANVDDRFC      = CLDPARAMS%ANV_DDRF
         LSDDRFC       = CLDPARAMS%LS_DDRF
         TANHRHCRIT    = INT( CLDPARAMS%TANHRHCRIT )
         MINRHCRIT_I   = CLDPARAMS%MINRHCRIT
         MAXRHCRIT     = CLDPARAMS%MAXRHCRIT
         TURNRHCRIT    = CLDPARAMS%TURNRHCRIT
         MAXRHCRITLAND = CLDPARAMS%MAXRHCRITLAND
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
         CFPBL_EXP     = CLDPARAMS%CFPBL_EXP
         PDFFLAG       = INT(CLDPARAMS%PDFSHAPE)
         thl2tune      = 1.0
         qw2tune       = 1.0
         qwthl2tune    = 1.0
         fac_cond      = MAPL_ALHL/MAPL_CP
         fac_fus       = MAPL_ALHF/MAPL_CP
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

       ! Outside of coherent convective regions bring RHCRIT up to MAXRHCRIT 
       !   to remove some resolution sensitivity in low-level stratus clouds
!         MINRHCRIT  = MAXRHCRIT*(1.0-CNV_FRACTION_dev(I)) + MINRHCRIT_I*(CNV_FRACTION_dev(I))
         MINRHCRIT = MINRHCRIT_I


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
            QPLS_dev(I,K)       = 0.
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
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QLW_AN_dev(I,K), &
                  QIW_LS_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  CNV_FRACTION_dev(I), SNOMAS_dev(I), FRLANDICE_dev(I), FRLAND_dev(I))
            ELSE
            CALL meltfrz (         &
                  DT             , &
                  TEMP           , &
                  QLW_LS_dev(I,K), & 
                  QIW_LS_dev(I,K), &
                  CNV_FRACTION_dev(I), SNOMAS_dev(I), FRLANDICE_dev(I), FRLAND_dev(I))
            CALL meltfrz (         &
                  DT             , &
                  TEMP           , &
                  QLW_AN_dev(I,K), & 
                  QIW_AN_dev(I,K), &
                  CNV_FRACTION_dev(I), SNOMAS_dev(I), FRLANDICE_dev(I), FRLAND_dev(I))
            ENDIF

            FRZ_TT_dev(I,K) = ( QIW_AN_dev(I,K) + QIW_LS_dev(I,K) - FRZ_TT_dev(I,K) ) / DT

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            DCNVi_dev(I,K) = QIW_AN_dev(I,K)
            DCNVL_dev(I,K) = QLW_AN_dev(I,K)

            CALL cnvsrc (          &  
                  DT             , &
                  CNVICEPARAM    , &
                  SCLMFDFR       , &
                  MASS           , & 
                  iMASS          , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLWDTR_dev(I,K), &
                  RMFDTR_dev(I,K), &
                  SC_QLWDTR_dev(I,K), &
                  SC_QIWDTR_dev(I,K), &
                  SC_RMFDTR_dev(I,K), &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  CLDFRC_dev(I,K), & 
                  ANVFRC_dev(I,K), &
                  QST3_dev(I,K)  , &
                  CNV_FRACTION_dev(I), SNOMAS_dev(I), FRLANDICE_dev(I), FRLAND_dev(I), &
                  CONVPAR_OPTION )

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
  
            call pdf_spread (&
                  K,LM,&
                  U_dev(I,K),U_above,U_below,&
                  V_dev(I,K),V_above,V_below,&
                  KH_dev(I,K-1),DZET_above,DZET_below,&
                  CNV_UPDFRC_dev(I,K),PP_dev(I,K),ALPHA,&
                  ALPHT_dev(I,K),ALPH1_dev(I,K),ALPH2_dev(I,K), & 
                  FRLAND_dev(I),&
                  CONVPAR_OPTION   ) 

            ! impose a minimum amount of variability
            ALPHA    = MAX(  ALPHA , 1.0 - RH00 )
   
            RHCRIT = 1.0 - ALPHA

            LSPDFLIQNEW = QLW_LS_dev(I,K)
            LSPDFICENEW = QIW_LS_dev(I,K)
            LSPDFFRACNEW= CLDFRC_dev(I,K)

            IF(USE_AEROSOL_NN) THEN
            call hystpdf_new(      &
                  DT             , &
                  ALPHA          , &
                  PDFFLAG        , &
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
                  au_dev(I,K),     &
                  hle_dev(I,K),    &
                  qte_dev(I,K),    &
                  hl2u_dev(I,K),   &
                  qt2u_dev(I,K),   &
                  hlqtu_dev(I,K),  &
                  whl_dev(I,K),        &
                  wqt_dev(I,K),        &
                  hl2_dev(I,K),        &
                  qt2_dev(I,K),        &
                  hlqt_dev(I,K),       & 
                  w3_dev(I,K),         &
                  w2_dev(I,K),         &
                  qt3_dev(I,K),        &
                  mf_frc_dev(I,K),     &
                  PDF_A_dev(I,K),      &  ! can remove these after development
                  PDF_SIGW_dev(I,K),   &
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
                  WTHV2_dev(I,K),      &
                  wql_dev(I,K),        &
                  SKEW_QT_dev(I,K),    &
                  CNV_FRACTION_dev(I), SNOMAS_dev(I), FRLANDICE_dev(I), FRLAND_dev(I), &
                  A_cloud(I,K),        &
                  B_cloud(I,K),        &
                  qsat(I,K))
 
            else
            call hystpdf(          &
                  DT             , &
                  ALPHA          , &
                  PDFFLAG        , &
                  PP_dev(I,K)    , &
                  Q_dev(I,K)     , &
                  QLW_LS_dev(I,K), &
                  QLW_AN_dev(I,K), &
                  QIW_LS_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  TEMP           , &
                  CLDFRC_dev(I,K), & 
                  ANVFRC_dev(I,K), &
                  CNV_FRACTION_dev(I), SNOMAS_dev(I), FRLANDICE_dev(I), FRLAND_dev(I)  )
            endif

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

             ! 'Anvil' partition from RAS/Parameterized not done in hystpdf

            call evap3(            &
                  DT             , &
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  CLDFRC_dev(I,K), &
                  NACTL_dev(I,K) , &
                  NACTI_dev(I,K) , &
                  QST3_dev(I,K)  )  

            call subl3(            &
                  DT             , & 
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  CLDFRC_dev(I,K), &
                  NACTL_dev(I,K) , &
                  NACTI_dev(I,K) , &
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
            if (CNV_FRACTION_dev(I) >= 0.5) FRACTION_REMOVAL = fr_an_ice

            SELECT CASE( ICE_SETTLE )
            CASE( 0 )
             ! MERRA-2 Formulation
              TROPICAL      = ANV_ICEFALL_C*1.0
              EXTRATROPICAL = ANV_ICEFALL_C*0.0
            CASE( 1 )
              TROPICAL      = CNV_FRACTION_dev(I) 
              EXTRATROPICAL = 1.0-TROPICAL
              TROPICAL      = ANV_ICEFALL_C*TROPICAL
              EXTRATROPICAL =  LS_ICEFALL_C*EXTRATROPICAL
            END SELECT

            CALL SETTLE_VEL(       &
                  WRHODEP        , &
                  QIW_AN_dev(I,K), &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  ANVFRC_dev(I,K), &
                  KH_dev(I,K-1)  , &
                  VFALL          , &
                  EXTRATROPICAL, TROPICAL, TROPP_dev(I) )
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
              TROPICAL      = CNV_FRACTION_dev(I)
              EXTRATROPICAL = 1.0-TROPICAL
              TROPICAL      = ANV_ICEFALL_C*TROPICAL
              EXTRATROPICAL =  LS_ICEFALL_C*EXTRATROPICAL
            END SELECT

            CALL SETTLE_VEL(       &
                  WRHODEP        , &
                  QIW_LS_dev(I,K), &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  CLDFRC_dev(I,K), &
                  KH_dev(I,K-1)  , &
                  VFALL          , &
                  EXTRATROPICAL, TROPICAL, TROPP_dev(I) )
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
               if ( TOT_PREC_ANV > 0.0 ) AREA_ANV_PRC = MAX( AREA_ANV_PRC/TOT_PREC_ANV, 1.E-6 )
               if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )
               if ( TOT_PREC_SC > 0.0 ) AREA_SCUP_PRC = MAX( AREA_SCUP_PRC/TOT_PREC_SC, 1.E-6 )
               if ( TOT_PREC_LS  > 0.0 ) AREA_LS_PRC  = MAX( AREA_LS_PRC/TOT_PREC_LS,   1.E-6 )

               AREA_LS_PRC  = LS_BETA  * AREA_LS_PRC
               AREA_UPD_PRC = CNV_BETA * AREA_UPD_PRC
               AREA_ANV_PRC = ANV_BETA * AREA_ANV_PRC

               !! "couple" to diagnostic areal fraction output 
               !! Intensity factor in PRECIP3 is floored at
               !! 1.0. So this is fair.

               LSARF_dev(I) = MIN( AREA_LS_PRC,  1.0 )
               CUARF_dev(I) = MIN( AREA_UPD_PRC, 1.0 )
               SCARF_dev(I) = MIN( AREA_SCUP_PRC, 1.0 )
               ANARF_dev(I) = MIN( AREA_ANV_PRC, 1.0 )
            END IF

            QRN_ALL = 0.
            QSN_ALL = 0.

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! GET SOME MICROPHYSICAL QUANTITIES 

            CALL MICRO_AA_BB_3( TEMP,PP_dev(I,K),QST3_dev(I,K),AA3,BB3 )

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP1 = QLW_LS_dev(I,K) + QLW_AN_dev(I,K)
            QTMP2 = QIW_LS_dev(I,K) + QIW_AN_dev(I,K)

        !- GF scheme handles its own conv precipitation.
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

        if(DOSHLW==1) THEN 
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
               if (VFALLRN.NE.0. .AND. VFALLSN.NE.0.) then
                  QPLS_dev(I,K) = QPLS_dev(I,K) + PFL_LS_dev(I,K)/VFALLRN + PFI_LS_dev(I,K)/VFALLSN
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
               RAD_CLDFRC_dev(I,K) = 0.
               CLDREFFL_dev(I,K)   = 0.
               CLDREFFI_dev(I,K)   = 0.
            ELSE
               call RADCOUPLE ( TEMP, PP_dev(I,K), PPE_dev(I,K)-PPE_dev(I,K-1), KH_dev(I,K-1), DTS_dev(I), CLDFRC_dev(I,K), ANVFRC_dev(I,K), &
                     Q_dev(I,K), QLW_LS_dev(I,K), QIW_LS_dev(I,K), QLW_AN_dev(I,K), QIW_AN_dev(I,K), QRN_ALL, QSN_ALL, NACTL_dev(I,K), NACTI_dev(I,K), & 
                     RAD_QV_dev(I,K), RAD_QL_dev(I,K), RAD_QI_dev(I,K), RAD_QR_dev(I,K), RAD_QS_dev(I,K), RAD_CLDFRC_dev(I,K), & 
                     CLDREFFL_dev(I,K), CLDREFFI_dev(I,K), &
                     TEMPOR_dev(I), FRLAND_dev(I), CNV_FRACTION_dev(I), &
                     SCLMFDFR, FR_AN_WAT, CONVPAR_OPTION,RHX_DEV(I,K) )
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
   subroutine pdf_spread (K,LM,&
         U,U_above,U_below,&
         V,V_above,V_below,&
         KH,&
         DZ_above,DZ_below,&
         UPDF,PP,ALPHA,&
         ALPHT_DIAG, ALPH1_DIAG, ALPH2_DIAG,&
         FRLAND , & 
         CONVPAR_OPTION   ) 

      integer, intent(in)  :: k,lm
      real,    intent(in)  :: U,U_above,U_below
      real,    intent(in)  :: V,V_above,V_below
      real,    intent(in)  :: DZ_above,DZ_below
      real,    intent(in)  :: UPDF,PP
      real,    intent(in)  :: KH
      real,    intent(out) :: ALPHA
      real,    intent(out) :: ALPH1_DIAG, ALPH2_DIAG, ALPHT_DIAG
      real,    intent(in)  :: FRLAND
      character(LEN=*), INTENT(IN) :: CONVPAR_OPTION

      real    :: A1,A2,A3
      real    :: tempmaxrh

      ! alpha is the 1/2*width so RH_crit=1.0-alpha

      if (tanhrhcrit.eq.1) then

         !  Use Slingo-Ritter (1985) formulation for critical relative humidity
         !  array a1 holds the critical rh, ranges from 0.8 to 1

         tempmaxrh = maxrhcrit
         if (frland > 0.05) tempmaxrh = maxrhcritland
         a1 = 1.0
         if (pp .le. turnrhcrit) then
            a1 = minrhcrit
         else
            a1 = minrhcrit + (tempmaxrh-minrhcrit)/(19.) * &
                  ((atan( (2.*(pp- turnrhcrit)/(1020.-turnrhcrit)-1.) * &
                  tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
         end if

         a1 = min(a1,1.)

         alpha = 1. - a1

      else

         alpha = 0.001 ! 0.1% RH SLOP

         !! DIRECTIONAL SHEAR == ABS( e_normal dot [U_z,V_z] ) 
   
         A1 = 0.

         A3 = 1./SQRT( U**2 + V**2 + 0.01 )  ! inverse of wind mag 

         A2 = V * A3 ! x-component of unit normal to (U,V) 

         if (k > 1 .and. k < lm) then
            A1 = ( A2 * ( U_above - U_below ) &
                  / ( DZ_above+DZ_below ) )
         end if

         A2 = -U * A3 ! y-component of unit normal to (U,V) 

         if (k > 1 .and. k < lm) then
            A1 = ( A2 * ( V_above - V_below )  & 
                  / ( DZ_above+DZ_below ) )  + A1
         end if

         A1 = ABS( A1 )  ! A1 is now magnitude of veering shear at layers in (m/s) /m.  Thus, A1=.001  ==> 1 m/s/km

         ALPHA = ALPHA  +  10.*A1  

         ALPH1_DIAG = 10.*A1

         !! Total shear = SQRT( [U_z,V_z] dot [U_z,V_z] )

         A1  = 0.
         if (k > 1 .and. k < lm) then
            A1 = ( ( U_above - U_below )/ ( DZ_above+DZ_below ) )**2 & 
                  + ( ( V_above - V_below )/ ( DZ_above+DZ_below ) )**2  
         end if

         A1  = SQRT ( A1 )  ! A1 is now magnitude of TOTAL shear at layers in (m/s) /m.  Thus, A1=.001  ==> 1 m/s/km

         ALPHA = ALPHA  +  3.33*A1

         !! KH values ~100 m+2 s-1 typical of strong PBLs

         ALPHA = ALPHA  +  0.002*KH

         ALPH2_DIAG = 0.002*KH

      end if               ! end of slingo ritter if-sequence

      ALPHA = MIN( ALPHA , 0.25 )  ! restrict RHcrit to > 75% 
      ALPHT_DIAG = ALPHA

   end subroutine pdf_spread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine fix_up_clouds( &
         QV, &
         TE, &
         QLC,&
         QIC,&
         CF, &
         QLA,&
         QIA,&
         AF  )

      real, intent(inout) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA

      ! Fix if Anvil cloud fraction too small
      if (AF < 1.E-5) then
         QV  = QV + QLA + QIA
         TE  = TE - (MAPL_ALHL/MAPL_CP)*QLA - (MAPL_ALHS/MAPL_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if

      ! Fix if LS cloud fraction too small
      if ( CF < 1.E-5 ) then
         QV = QV + QLC + QIC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC - (MAPL_ALHS/MAPL_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if
      
      ! LS LIQUID too small
      if ( QLC  < 1.E-8 ) then
         QV = QV + QLC 
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC
         QLC = 0.
      end if
      ! LS ICE too small
      if ( QIC  < 1.E-8 ) then
         QV = QV + QIC 
         TE = TE - (MAPL_ALHS/MAPL_CP)*QIC
         QIC = 0.
      end if

      ! Anvil LIQUID too small
      if ( QLA  < 1.E-8 ) then
         QV = QV + QLA 
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLA
         QLA = 0.
      end if
      ! Anvil ICE too small
      if ( QIA  < 1.E-8 ) then
         QV = QV + QIA 
         TE = TE - (MAPL_ALHS/MAPL_CP)*QIA
         QIA = 0.
      end if

      ! Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
      if ( ( QLA + QIA ) < 1.E-8 ) then
         QV = QV + QLA + QIA
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLA - (MAPL_ALHS/MAPL_CP)*QIA
         AF  = 0.
         QLA = 0.
         QIA = 0.
      end if
      ! Ditto if LS cloud LIQUID+ICE too small
      if ( ( QLC + QIC ) < 1.E-8 ) then
         QV = QV + QLC + QIC
         TE = TE - (MAPL_ALHL/MAPL_CP)*QLC - (MAPL_ALHS/MAPL_CP)*QIC
         CF  = 0.
         QLC = 0.
         QIC = 0.
      end if

   end subroutine fix_up_clouds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine meltfrz( DT, TE, QL, QI, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND)

      real, intent(in)    :: DT, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND
      real, intent(inout) :: TE,QL,QI

      real  :: fQi,dQil

      real  ::  taufrz

      integer :: K

      ! freeze liquid
      if ( TE <= MAPL_TICE ) then
         fQi  = ice_fraction( TE, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
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
         TE       , &
         QCL      , &
         QAL      , &
         QCI      , &
         QAI      , &
         CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND)

      real ,   intent(inout) :: TE,QCL,QCI,QAL,QAI
      real ,   intent(in   ) :: Dt,CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND
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
      fQi = ice_fraction( TE, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
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
      fQi = ice_fraction( TE, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
      if ( fQi > 0.0 ) then
         DQmax = (MAPL_TICE-TE)*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)
         dQil  = QLTOT *(1.0 - EXP( -Dt * fQi / taufrz ) )
        !dQil  = QLTOT*fQi
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine hystpdf_new( &
         DT          , &
         ALPHA       , &
         PDFSHAPE    , &
         PL          , &
         ZL          , &
         QV          , &
         QCl         , &
         QAl         , &
         QCi         , &
         QAi         , &
         TE          , &
         CF          , &
         AF          , &
         NL          , &
         NI          , &
         AU          , &
         HLE         , &
         QTE         , &
         HL2U        , &
         QT2U        , &
         HLQTU       , &
         WHL         , &
         WQT         , &
         HL2         , &
         QT2         , &
         HLQT        , & 
         W3          , &
         W2          , &
         MFQT3       , &
         MF_FRC      , &
         PDF_A,      &  ! can remove these after development
         PDF_SIGW,   &
         PDF_W1,     &
         PDF_W2,     &
         PDF_SIGHL1, &
         PDF_SIGHL2, &
         PDF_HL1,    &
         PDF_HL2,    &
         PDF_SIGQT1, &
         PDF_SIGQT2, &
         PDF_QT1,    &
         PDF_QT2,    &
         PDF_RHLQT,  &
         WTHV2,      &
         WQL,        &
         SKEW_QT,    &
         CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND, &
         A_cloud,    &
         B_cloud,    &
         qsat)

      real, intent(in)    :: DT,ALPHA,PL,ZL
      integer, intent(in) :: pdfshape
      real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF,SKEW_QT,PDF_A
      real, intent(in)    :: NL,NI,CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND
!      real, intent(in)    :: HL,WHL,WQT,HL2,QT2,HLQT,W3,W2,MF_FRC,MFQT3
      real, intent(in)    :: WHL,WQT,HL2,QT2,HLQT,W3,W2,MF_FRC,MFQT3
      real, intent(in)    :: AU, HLE, QTE, HL2U, QT2U, HLQTU
      real, intent(out)   :: PDF_SIGW, PDF_W1, PDF_W2, &
                             PDF_SIGHL1, PDF_SIGHL2, PDF_HL1, PDF_HL2, &
                             PDF_SIGQT1, PDF_SIGQT2, PDF_QT1, PDF_QT2, &
                             PDF_RHLQT
      real, intent(out)   :: WTHV2, WQL
      real, intent(out)   :: A_cloud, B_cloud, qsat

      ! internal arrays
      real :: QCO, QVO, CFO, QAO, TAU,HL
      real :: QT, QMX, QMN, DQ, sigmaqt1, sigmaqt2

      real :: TEO,QSx,DQsx,QS,DQs

      real :: TEp, QSp, CFp, QVp, QCp
      real :: TEn, QSn, CFn, QVn, QCn

      real :: QCx, QVx, CFx, QAx, QC, QA, fQi
      real :: dQAi, dQAl, dQCi, dQCl, Nfac, NLv, NIv 

      real :: Tce, qle, ace, Tcu, qlu, acu, HLU, QTU

!      real :: fQip

      real :: tmpARR
      real :: ALHX, DQCALL
      ! internal scalars
      integer :: N, nmax

      pdfflag = PDFSHAPE

      QC = QCl + QCi
      QA = QAl + QAi
      QT  =  QC  + QA + QV  !Total water after microphysics
      tmpARR = 0.0
      nmax =  20
      QAx = 0.0

      if ( AF < 1.0 )  tmpARR = 1./(1.-AF)

      TEo = TE

      fQi = ice_fraction( TE, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
      DQSx  = DQSAT( TE, PL, QSAT=QSx )
      CFx = CF*tmpARR
      QCx = QC*tmpARR
      QVx = ( QV - QSx*AF )*tmpARR

      if ( AF >= 1.0 )    QVx = QSx*1.e-4 
      if ( AF > 0. )  QAx = QA/AF

      QT  = QCx + QVx

      TEp = TEo
      QSn = QSx
      TEn = TEo
      CFn = CFx
      QVn = QVx
      QCn = QCx
      DQS = DQSx

      do n=1,nmax

         QVp = QVn
         QCp = QCn
         CFp = CFn
         TEp = TEn
!         fQip= fQi

         if(pdfflag.lt.2) then

            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn

         elseif(pdfflag.eq.2) then
            ! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
            ! for triangular, skewed r : sigmaqt1 < sigmaqt2
            ! try: skewed right below 500 mb
!!!       if(pl.lt.500.) then
            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn
!!!       else
!!!       sigmaqt1  = 2*ALPHA*QSn*0.4
!!!       sigmaqt2  = 2*ALPHA*QSn*0.6
!!!       endif
         elseif(pdfflag .eq. 4) then !lognormal (sigma is dimmensionless)
            sigmaqt1 =  max(ALPHA/sqrt(3.0), 0.001)
         endif

         if (pdfflag.lt.5) then
           call pdffrac(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,CFn)
           call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)
         elseif (pdfflag.eq.5) then
!           if (ZL<400.) then
!             print *,'n=',n,'  ZL=',ZL,'  TEn=',TEn
!             print *,'HL=',HL,'  QT=',QT,'  QV=',QVn 
!             print *,'WHL=',WHL,'  WQT=',WQT,'  HL2=',HL2
!             print *,'QT2=',QT2,'  HLQT=',HLQT,'  W3=',W3
!           end if

           ! Update the liquid water static energy
           ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS
           HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn
!                fac_cond*QLW_LS_dev(I,:) - fac_fus*QIW_LS_dev(I,:) 
           QT = QVn+QCn

           call partition_dblgss(DT,           &
                                 TEn,          &
                                 QVn,          &
                                 QCn,          &
                                 0.0,          & ! assume OMEGA=0
                                 ZL,           &
                                 PL*100.,      &
!                                 qpl,         &
!                                 qpi,         &
                                 QT,           &
                                 HL,          &
                                 WHL,         &
                                 WQT,         &
                                 HL2,         &
                                 QT2,         &
                                 HLQT,        & 
                                 W3,           &
                                 W2,           &
                                 MFQT3,        &
                                 MF_FRC,       &
                                 SKEW_QT,      &
                                 PDF_A,        &
                                 PDF_SIGW,     &
                                 PDF_W1,       &
                                 PDF_W2,       &
                                 PDF_SIGHL1,   &
                                 PDF_SIGHL2,   &
                                 PDF_HL1,      &
                                 PDF_HL2,      &
                                 PDF_SIGQT1,   &
                                 PDF_SIGQT2,   &
                                 PDF_QT1,      &
                                 PDF_QT2,      &
                                 PDF_RHLQT,    &
                                 WTHV2,        &
                                 WQL,          &
                                 CFn)

           fQi = ice_fraction( TEn, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )


! TEMP:  overwrite ADG with uniform PDF
!           sigmaqt1  = ALPHA*QSn
!           sigmaqt2  = ALPHA*QSn
!           qt = qvp+qcp
!           call pdffrac(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,CFn)
!           call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)
! end TEMP
         elseif (pdfflag == 6) then ! Single gaussian
            ! Update the liquid water static energy
            ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS
            HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn
!                fac_cond*QLW_LS_dev(I,:) - fac_fus*QIW_LS_dev(I,:)
            QT = QVn + QCn
            
            call gaussian(ZL, 100.*PL, HL, QT, HL2, QT2, HLQT, &
                          TEn, QCn, CFn, &
                          A_cloud, B_cloud, qsat)
            
            fQi = ice_fraction( TEn, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
         elseif (pdfflag == 7) then ! Double Gaussian with consistent partitioning
            ! Update the liquid water static energy
            ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

            HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn
            QT = QVn + QCn

            !
            ! Updraft ensemble
            !
            if ( au > 0. ) then
               HLU = ( HL - ( 1. - au )*HLE )/au
               QTU = ( QT - ( 1. - au )*QTE )/au

               call gaussian(ZL, 100.*PL, HLU, QTU, HL2U, QT2U, HLQTU, &
                             Tcu, qlu, acu, &
                             A_cloud, B_cloud, qsat)

            end if

            !
            ! Environment
            !
            
            call gaussian(ZL, 100.*PL, HLE, QTE, HL2, QT2, HLQT, &
                          Tce, qle, ace, &
                          A_cloud, B_cloud, qsat)

            !
            ! Combine upddraft and environment
            !

            if ( au > 0. ) then
               TEn = au*Tcu + ( 1. - au )*Tce
               QCn = au*qlu + ( 1. - au )*qle
               CFn = au*acu + ( 1. - au )*ace
            else
               TEn = Tce
               QCn = qle
               CFn = ace
            end if

            fQi = ice_fraction( TEn, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )

            !
!            write(*,*) au, HL, HLU, HLE
!            write(*,*) au, QT, QTU, QTE
         endif

!         if (abs(QVn+QCn-QVp-QCp)>1e-6*(QVp+QCp) .and. QVp>0.0001) print *,'total water not conserved!'

!         if (pl>950.) print *,'hystpdf, af dblgss: wthv2=',wthv2

         IF(USE_AEROSOL_NN) THEN
           DQCALL = QCn - QCp
           CF  = CFn * ( 1.-AF)
           Nfac = 100.*PL*R_AIR/TEp !density times conversion factor
           NLv = NL/Nfac
           NIv = NI/Nfac
           call Bergeron_iter    (  &         !Microphysically-based partitions the new condensate
                 DT               , &
                 PL               , &
                 TEp              , &
                 QT               , &
                 QCi              , &
                 QAi              , &
                 QCl              , &
                 QAl              , &
                 CF               , &
                 AF               , &
                 NLv              , &
                 NIv              , &
                 DQCALL           , &
                 fQi                , & 
                 CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND, &
                 .false.)
         ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! These lines represent adjustments
         ! to anvil condensate due to the 
         ! assumption of a stationary TOTAL 
         ! water PDF subject to a varying 
         ! QSAT value during the iteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         if ( AF > 0. ) then
            QAo = QAx  ! + QSx - QS 
         else
            QAo = 0.
         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

         if(pdfflag.eq.1) then 
            QCn = QCp + ( QCn - QCp ) / ( 1. - (CFn * (ALPHA-1.) - (QCn/QSn))*DQS*ALHX/MAPL_CP)             
         elseif(pdfflag.eq.2 .or. pdfflag.eq.5) then
            ! This next line needs correcting - need proper d(del qc)/dT derivative for triangular
            ! for now, just use relaxation of 1/2.
            if (n.ne.nmax) QCn = QCp + ( QCn - QCp ) *0.5
         endif

         QVn = QVp - (QCn - QCp)
         TEn = TEp + (1.0-fQi)*(MAPL_ALHL/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
               +      fQi* (MAPL_ALHS/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

         if (abs(Ten - Tep) .lt. 0.00001) exit 
!         if (abs(Ten-Tep)>2.0) print *,'Ten-Tep large, Ten=',Ten,'  Tep=',Tep,'  PL=',PL,'  n=',n,'  QCn=',QCn,'  QCp=',QCp,'  AF=',AF

         DQS  = DQSAT( TEn, PL, QSAT=QSn )

      enddo ! qsat iteration

!      if (abs(TEo-TEn)>2.0) print *,'TEn=',TEn,' TEo=',TEo,' PL=',PL,' QCn=',QCn,' QCp=',QCp
!      if (abs(QVo-QVn)>0.1*QVo .and. Qvo>0.001) print *,'QVo-QVn large, QVo=',Qvo,'  QVn=',QVn

      CFo = CFn
      CF = CFn
      QCo = QCn
      QVo = QVn
      TEo = TEn

      ! Update prognostic variables.  Deal with special case of AF=1
      ! Temporary variables QCo, QAo become updated grid means.
      if ( AF < 1.0 ) then
         CF  = CFo * ( 1.-AF)
         QCo = QCo * ( 1.-AF)
         QAo = QAo *   AF  
      else

         ! Special case AF=1, i.e., box filled with anvil. 
         !   - Note: no guarantee QV_box > QS_box
         CF  = 0.          ! Remove any other cloud
         QAo = QA  + QC    ! Add any LS condensate to anvil type
         QCo = 0.          ! Remove same from LS   
         QT  = QAo + QV    ! Total water
         ! Now set anvil condensate to any excess of total water 
         ! over QSx (saturation value at top)
         QAo = MAX( QT - QSx, 0. )
      end if

      ! Now take {\em New} condensate and partition into ice and liquid
      ! taking care to keep both >=0 separately. New condensate can be
      ! less than old, so $\Delta$ can be < 0.

      dQCl = 0.0
      dQCi = 0.0
      dQAl = 0.0
      dQAi = 0.0

      !large scale   

      QCx   = QCo - QC
      if  (QCx .lt. 0.0) then  !net evaporation. Water evaporates first
         dQCl = max(QCx, -QCl)   
         dQCi = max(QCx - dQCl, -QCi)
      else
         dQCl  = (1.0-fQi)*QCx
         dQCi  =    fQi  * QCx
      end if

      !Anvil   
      QAx   = QAo - QA

      if  (QAx .lt. 0.0) then  !net evaporation. Water evaporates first
         dQAl = max(QAx, -QAl)   
         dQAi = max(QAx - dQAl, -QAi)
      else            
         dQAl  =  (1.0-fQi)*QAx
         dQAi  = QAx*fQi
      end if

      ! Clean-up cloud if fractions are too small
      if ( AF < 1.e-5 ) then
         dQAi = -QAi
         dQAl = -QAl
      end if
      if ( CF < 1.e-5 ) then
         dQCi = -QCi
         dQCl = -QCl
      end if

      QAi    = QAi + dQAi
      QAl    = QAl + dQAl
      QCi    = QCi + dQCi
      QCl    = QCl + dQCl
      QV     = QV  - ( dQAi+dQCi+dQAl+dQCl) 

      TE  = TE + (MAPL_ALHL*( dQAi+dQCi+dQAl+dQCl)+MAPL_ALHF*(dQAi+dQCi))/ MAPL_CP

      ! We need to take care of situations where QS moves past QA
      ! during QSAT iteration. This should be only when QA/AF is small
      ! to begin with. Effect is to make QAo negative. So, we 
      ! "evaporate" offending QA's
      !
      ! We get rid of anvil fraction also, although strictly
      ! speaking, PDF-wise, we should not do this.
      if ( QAo <= 0. ) then
         QV  = QV + QAi + QAl
         TE  = TE - (MAPL_ALHS/MAPL_CP)*QAi - (MAPL_ALHL/MAPL_CP)*QAl
         QAi = 0.
         QAl = 0.
         AF  = 0.  
      end if

   end subroutine hystpdf_new

   ! Single-gaussian cloud pdf
   subroutine gaussian(z, p, hl, qt, hl2, qt2, hlqt, &
                       T, ql, ac, &
                       A, B, qs)

     use MAPL_SatVaporMod,  only: MAPL_EQsat
     use MAPL_ConstantsMod, only: MAPL_CP, MAPL_ALHL, MAPL_GRAV, MAPL_RDRY, MAPL_RVAP, MAPL_PI

     real, intent(in)            :: z, p, hl, qt, hl2, qt2, hlqt
     real, intent(inout)         :: T, ql, ac
     real, intent(out), optional :: A, B, qs                                                                                                            
     real :: dqs, fac_cond, Tl, s, sigma_s, exner, Q

     exner    = (p*1.E-5)**(MAPL_RDRY/MAPL_CP) ! Exner function
     fac_cond = MAPL_ALHL/MAPL_CP              ! lv/cp 

     Tl = hl - (MAPL_GRAV/MAPL_CP)*z
!     Tl = T - fac_cond*ql        ! liquid water temperature                   
     qs = MAPL_EQsat(Tl, p, dqs) ! saturation specific humidity 
     s  = qt - qs                ! saturation excess/deficit 

     A = 1./( 1. + fac_cond*dqs )
     B = a*exner*dqs

     sigma_s = sqrt( A**2.*qt2 - 2*A*B*(hlqt/exner) + B**2.*(hl2/exner**2.) )

     ! Diagnose cloud properties 
     if (sigma_s > 0.) then
        Q = A*(qt - qs)/sigma_s

        ac = 0.5*( 1. + erf(Q/sqrt(2.)) )
        ql = sigma_s*( ac*Q + exp(-0.5*Q**2.)/sqrt(2.*MAPL_PI) )
     else
        ac = 0.
        ql = 0.
     end if

     T = hl + (MAPL_ALHL/MAPL_CP)*ql - (MAPL_GRAV/MAPL_CP)*z ! Update temperature 

   end subroutine gaussian

#ifdef _CUDA
   attributes(device) &
#endif
   subroutine hystpdf( &
         DT          , &
         ALPHA       , &
         PDFSHAPE    , &
         PL          , &
         QV          , &
         QCl         , &
         QAl         , &
         QCi         , &
         QAi         , &
         TE          , &
         CF          , &
         AF          , &
         CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND)

      real, intent(in)    :: DT,ALPHA,PL,CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND
      integer, intent(in) :: pdfshape
      real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF

      ! internal arrays
      real :: QCO, QVO, CFO, QAO, TAU
      real :: QT, QMX, QMN, DQ, QVtop, sigmaqt1, sigmaqt2

      real :: TEO,QSx,DQsx,QS,DQs

      real :: TEp, QSp, CFp, QVp, QCp
      real :: TEn, QSn, CFn, QVn, QCn

      real :: QCx, QVx, CFx, QAx, QC, QA, fQi, fQi_A
      real :: dQAi, dQAl, dQCi, dQCl 

      real :: tmpARR
      real :: ALHX
      ! internal scalars
      integer :: N

      QC = QCl + QCi
      QA = QAl + QAi

      if ( QA > 0.0 ) then
         fQi_A = QAi / QA 
      else
         fQi_A = 0.0
      end if

      TEo = TE

      DQSx  = DQSAT( TEo, PL, QSAT=QSx )

      if ( AF < 1.0 ) then
         tmpARR = 1./(1.-AF)
      else
         tmpARR = 0.0
      end if

      CFx = CF*tmpARR
      QCx = QC*tmpARR
      QVx = ( QV - QSx*AF )*tmpARR

      if ( AF >= 1.0 ) then
         QVx = QSx*1.e-4 
      end if


      if ( AF > 0. ) then
         QAx = QA/AF
      else
         QAx = 0.
      end if

      QT  = QCx + QVx

      TEp = TEo
      QSn = QSx
      TEn = TEo
      CFn = CFx
      QVn = QVx
      QCn = QCx

      DQS = DQSx

      do n=1,4

         QSp = QSn
         QVp = QVn
         QCp = QCn
         CFp = CFn

         DQS  = DQSAT( TEn, PL, QSAT=QSn )

         TEp = TEn
         fQi = ice_fraction( TEp, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )

         sigmaqt1  = ALPHA*QSn
         sigmaqt2  = ALPHA*QSn

!! WMP - 2018-05-17 - This section of code is repetitive (so skip it)....
#ifdef SKIP_WMP
         if(pdfflag.eq.2) then
! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
! for triangular, skewed r : sigmaqt1 < sigmaqt2
! try: skewed right below 500 mb
!!!       if(pl.lt.500.) then
          sigmaqt1  = ALPHA*QSn
          sigmaqt2  = ALPHA*QSn
!!!       else
!!!       sigmaqt1  = 2*ALPHA*QSn*0.4
!!!       sigmaqt2  = 2*ALPHA*QSn*0.6
!!!       endif
         endif
#endif

         call pdffrac(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,CFn)
         call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! These lines represent adjustments
         ! to anvil condensate due to the 
         ! assumption of a stationary TOTAL 
         ! water PDF subject to a varying 
         ! QSAT value during the iteration
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
         if ( AF > 0. ) then
            QAo = QAx  ! + QSx - QS 
         else
            QAo = 0.
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

         if(pdfflag.eq.1) then
          QCn = QCp + ( QCn - QCp ) / ( 1. - (CFn * (ALPHA-1.) - (QCn/QSn))*DQS*ALHX/MAPL_CP)
         elseif(pdfflag.eq.2) then
! This next line needs correcting - need proper d(del qc)/dT derivative for triangular
! for now, just use relaxation of 1/2.
          if (n.ne.4) QCn = QCp + ( QCn - QCp ) *0.5
         endif
         QVn = QVp - (QCn - QCp)

         TEn = TEp + (1.0-fQi)*(MAPL_ALHL/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
               +      fQi* (MAPL_ALHS/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

      enddo ! qsat iteration

         CFo = CFn
         CF = CFn
         QCo = QCn
         QVo = QVn
         TEo = TEn

      ! Update prognostic variables.  Deal with special case of AF=1
      ! Temporary variables QCo, QAo become updated grid means.
      if ( AF < 1.0 ) then
         CF  = CFo * ( 1.-AF)
         QCo = QCo * ( 1.-AF)
         QAo = QAo *   AF  
      else

         ! Special case AF=1, i.e., box filled with anvil. 
         !   - Note: no guarantee QV_box > QS_box

         CF  = 0.          ! Remove any other cloud
         QAo = QA  + QC    ! Add any LS condensate to anvil type
         QCo = 0.          ! Remove same from LS   

         QT  = QAo + QV    ! Total water

         ! Now set anvil condensate to any excess of total water 
         ! over QSx (saturation value at top)
         QAo = MAX( QT - QSx, 0. )
      end if

      ! Now take {\em New} condensate and partition into ice and liquid
      ! taking care to keep both >=0 separately. New condensate can be
      ! less than old, so $\Delta$ can be < 0.

      QCx   = QCo - QC
      dQCl  = (1.0-fQi)*QCx
      dQCi  =    fQi  * QCx

      if ((QCl+dQCl)<0.) then
         dQCi = dQCi + (QCl+dQCl)
         dQCl = -QCl !== dQCl - (QCl+dQCl)
      end if
      if ((QCi+dQCi)<0.) then
         dQCl = dQCl + (QCi+dQCi)
         dQCi = -QCi !== dQCi - (QCi+dQCi)
      end if

      QAx   = QAo - QA
      dQAl  = QAx ! (1.0-fQi)*QAx
      dQAi  = 0.  !  fQi  * QAx

      if ((QAl+dQAl)<0.) then
         dQAi = dQAi + (QAl+dQAl)
         dQAl = -QAl
      end if
      if ((QAi+dQAi)<0.) then
         dQAl = dQAl + (QAi+dQAi)
         dQAi = -QAi 
      end if

      ! Clean-up cloud if fractions are too small
      if ( AF < 1.e-5 ) then
         dQAi = -QAi
         dQAl = -QAl
      end if
      if ( CF < 1.e-5 ) then
         dQCi = -QCi
         dQCl = -QCl
      end if

      QAi    = QAi + dQAi
      QAl    = QAl + dQAl
      QCi    = QCi + dQCi
      QCl    = QCl + dQCl
      QV     = QV  - ( dQAi+dQCi+dQAl+dQCl) 


      !!TE  = TE + (MAPL_ALHS/MAPL_CP)*(dQAi+dQCi) + (MAPL_ALHL/MAPL_CP)*(dQAl+dQCl)
      TE  = TE + (MAPL_ALHL*( dQAi+dQCi+dQAl+dQCl)+MAPL_ALHF*(dQAi+dQCi))/ MAPL_CP

      ! We need to take care of situations where QS moves past QA
      ! during QSAT iteration. This should be only when QA/AF is small
      ! to begin with. Effect is to make QAo negative. So, we 
      ! "evaporate" offending QAs
      !
      ! We get rid of anvil fraction also, although strictly
      ! speaking, PDF-wise, we should not do this.
      if ( QAo <= 0. ) then
         QV  = QV + QAi + QAl
         TE  = TE - (MAPL_ALHS/MAPL_CP)*QAi - (MAPL_ALHL/MAPL_CP)*QAl
         QAi = 0.
         QAl = 0.
         AF  = 0.  
      end if

   end subroutine hystpdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
      subroutine pdffrac (flag,qtmean,sigmaqt1,sigmaqt2,qstar,clfrac)
      implicit none

      integer flag            ! flag to indicate shape of pdf
                              ! 1 for tophat, 2 for triangular, 3 for Gaussian
      real qtmean             ! Grid box value of q total
      real sigmaqt1           ! width of distribution (sigma)
      real sigmaqt2           ! width of distribution (sigma)
      real qstar              ! saturation q at grid box avg T
      real clfrac             ! cloud fraction (area under pdf from qs)

      real :: qtmode, qtmin, qtmax

      if(flag.eq.1) then
       if((qtmean+sigmaqt1).lt.qstar) then
        clfrac = 0.
       else
        if(sigmaqt1.gt.0.) then
        clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
        else
        clfrac = 1.
        endif
       endif
      elseif(flag.eq.2) then
       qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.
       qtmin = max(qtmode-sigmaqt1,0.)
       qtmax = qtmode + sigmaqt2
       if(qtmax.lt.qstar) then
        clfrac = 0.
       elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
        clfrac = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
       elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
        clfrac = 1. - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
       elseif ( qstar.le.qtmin ) then
        clfrac = 1.
       endif
      endif

      return
      end subroutine pdffrac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif

      subroutine pdfcondensate (flag,qtmean4,sigmaqt14,sigmaqt24,qstar4,condensate4)
      implicit none

      integer flag            ! flag to indicate shape of pdf
                              ! 1 for tophat, 2 for triangular
      real qtmean4            ! Grid box value of q total
      real sigmaqt14          ! width of distribution (to left)
      real sigmaqt24          ! width of distribution (to right)
      real qstar4             ! saturation q at grid box avg T
      real condensate4        ! condensate (area under (q*-qt)*pdf from qs)

      real *8 :: qtmode, qtmin, qtmax, constA, constB, cloudf
      real *8 :: term1, term2, term3
      real *8 :: qtmean, sigmaqt1, sigmaqt2, qstar, condensate

      qtmean = dble(qtmean4)
      sigmaqt1 = dble(sigmaqt14)
      sigmaqt2 = dble(sigmaqt24)
      qstar = dble(qstar4)

      if(flag.eq.1) then
       if(qtmean+sigmaqt1.lt.qstar) then
        condensate = 0.d0
       elseif(qstar.gt.qtmean-sigmaqt1)then
        if(sigmaqt1.gt.0.d0) then
         condensate = (min(qtmean + sigmaqt1 - qstar,2.d0*sigmaqt1)**2)/ (4.d0*sigmaqt1)
        else
         condensate = qtmean-qstar
        endif
       else
        condensate = qtmean-qstar
       endif
      elseif(flag.eq.2) then
       qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.d0
       qtmin = max(qtmode-sigmaqt1,0.d0)
       qtmax = qtmode + sigmaqt2
       if ( qtmax.lt.qstar ) then
        condensate = 0.d0
       elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
        constB = 2.d0 / ( (qtmax - qtmin)*(qtmax-qtmode) )
        cloudf = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
        term1 = (qstar*qstar*qstar)/3.d0
        term2 = (qtmax*qstar*qstar)/2.d0
        term3 = (qtmax*qtmax*qtmax)/6.d0
        condensate = constB * (term1-term2+term3) - qstar*cloudf
       elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
        constA = 2.d0 / ( (qtmax - qtmin)*(qtmode-qtmin) )
        cloudf = 1.d0 - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
        term1 = (qstar*qstar*qstar)/3.d0
        term2 = (qtmin*qstar*qstar)/2.d0
        term3 = (qtmin*qtmin*qtmin)/6.d0
        condensate = qtmean - ( constA * (term1-term2+term3) ) - qstar*cloudf
       elseif ( qstar.le.qtmin ) then
        condensate = qtmean-qstar
       endif
      endif
      condensate4 = real(condensate)

      return
      end subroutine pdfcondensate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) & 
#endif
   subroutine cnvsrc( & 
         DT      , &
         ICEPARAM, &
         SCLMFDFR, &
         MASS    , &
         iMASS   , &
         PL      , &
         TE      , &
         QV      , &
         DCF     , &
         DMF     , &
         DCLFshlw, &
         DCIFshlw, &
         DMFshlw , &
         QLA     , & 
         QIA     , & 
         CF      , &
         AF      , &
         QS      , &
         CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND, &
         CONVPAR_OPTION )

      !INPUTS:
      !
      !       ICEPARAM: 0-1  controls how strongly new conv condensate is partitioned in ice-liquid
      !                 1 means partitioning follows ice_fraction(TE,CNV_FRACTION,SNOMAS,FRLANDICE,FRLAND). 0 means all new condensate is
      !                 liquid 
      !
      !       SCLMFDFR: Scales detraining mass flux to a cloud fraction source - kludge. Thinly justified
      !                 by fuzziness of cloud boundaries and existence of PDF of condensates (for choices
      !                 0.-1.0) or by subgrid layering (for choices >1.0) 

      real, intent(in)    :: DT,ICEPARAM,SCLMFDFR
      real, intent(in)    :: MASS,iMASS,QS
      real, intent(in)    :: DMF,PL
      real, intent(in)    :: DCF,CF,DCIFshlw,DCLFshlw,DMFshlw
      real, intent(inout) :: TE
      real, intent(inout) :: AF,QV
      real, intent(inout) :: QLA, QIA
      real, intent(in)    :: CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND
      character(LEN=*), INTENT(IN) :: CONVPAR_OPTION

      real :: TEND,QVx,QCA,fQi

      integer  :: STRATEGY
      real     :: minrhx 

      STRATEGY = 1

      !Minimum allowed env RH
      minrhx    = 0.001  

      !Addition of condensate from RAS
      TEND = DCF*iMASS
      fQi  = 0.0 + ICEPARAM*ice_fraction( TE, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
      QLA  = QLA + (1.0-fQi)* TEND*DT
      QIA  = QIA +    fQi   * TEND*DT

      ! dont forget that conv cond has never frozen !!!!
      !TE   = TE +   (MAPL_ALHS-MAPL_ALHL) * fQi * TEND*DT/MAPL_CP
      IF(ADJUSTL(CONVPAR_OPTION) /= 'GF') TE   = TE +   (MAPL_ALHS-MAPL_ALHL) * fQi * TEND*DT/MAPL_CP

      ! add shallow convective ice/liquid source
      QLA = QLA + DCLFshlw*iMASS*DT
      QIA = QIA + DCIFshlw*iMASS*DT

      QCA  = QLA + QIA

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Tiedtke-style anvil fraction !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     TEND=(DMF+DMFshlw*1.2)*iMASS * SCLMFDFR
      TEND=(DMF+DMFshlw)*iMASS * SCLMFDFR
      AF = AF + TEND*DT
      AF = MIN( AF , 0.99 )

      ! where ( (AF+CF) > 1.00 )
      !    AF=1.00-CF
      ! endwhere
      !!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Check for funny (tiny, negative)
      ! external QV s, resulting from assumed
      ! QV=QSAT within anvil.
      !
      ! Two strategies to fix 
      !   1) Simply constrain AF assume condensate
      !      just gets "packed" in
      !   2) Evaporate QCA to bring up QVx leave AF alone
      !      Should include QSAT iteration, but ...        
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )

      if ( AF < 1.0 ) then
         QVx  = ( QV - QS * AF )/(1.-AF)
      else
         QVx  = QS
      end if

      if (STRATEGY==1) then 
         if ( (( QVx - minrhx*QS ) < 0.0 ) .and. (AF > 0.) ) then
            AF  = (QV  - minrhx*QS )/( QS*(1.0-minrhx) )
         end if
         if ( AF < 0. ) then  ! If still cant make suitable env RH then destroy anvil
            AF  = 0.
            QV  = QV + QLA + QIA
            TE  = TE - (MAPL_ALHL*QLA + MAPL_ALHS*QIA)/MAPL_CP        
            QLA = 0.
            QIA = 0.
         end if
      else if (STRATEGY==2) then
         if ( (( QVx - minrhx*QS ) < 0.0 ) .and. (AF > 0.) ) then
            QV  = QV  + (1.-AF)*( minrhx*QS - QVx )
            QCA = QCA - (1.-AF)*( minrhx*QS - QVx )
            TE  = TE  - (1.-AF)*( minrhx*QS - QVx )*MAPL_ALHL/MAPL_CP
         end if
      end if

   end subroutine cnvsrc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine evap3(&
         DT      , &
         RHCR    , &
         PL      , &
         TE      , &
         QV      , &
         QL      , &
         QI      , &
         F       , &
         XF      , & 
         NL      , &
         NI      , &
         QS        )

      real, intent(in   ) :: DT 
      real, intent(in   ) :: RHCR
      real, intent(in   ) :: PL
      real, intent(inout) :: TE
      real, intent(inout) :: QV
      real, intent(inout) :: QL,QI
      real, intent(inout) :: F
      real, intent(in   ) :: XF
      real, intent(in   ) :: NL,NI
      real, intent(in   ) :: QS

      real :: ES,NN,RADIUS,K1,K2,TEFF,QCm,EVAP,RHx,QC  !,QS

      real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
      real, parameter :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
      real, parameter :: DIFFU   =  2.2e-5        ! m**2 s**-1

      real :: A_eff

      A_EFF = CLD_EVP_EFF

      NN = 50.*1.0e6

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         EVAPORATION OF CLOUD WATER.             !!
      !!                                                 !!
      !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
      !!  formulation  (Eq.s 15-17)                      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )
      
      ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100's <-^ convert from mbar to Pa)

      RHx = MIN( QV/QS , 1.00 )


      K1 = (MAPL_ALHL**2) * RHO_W / ( K_COND*MAPL_RVAP*(TE**2))

      K2 = MAPL_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Here DIFFU is given for 1000 mb  !!
      !! so 1000./PR accounts for inc-    !!
      !! reased diffusivity at lower      !!
      !! pressure.                        !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      if ( ( F > 0.) .and. ( QL > 0. ) ) then
         QCm=QL/F
      else
         QCm=0.
      end if

      RADIUS = LDRADIUS4(PL,TE,QCm,NN,RHX,NL,NI,1)
      
      if ( (RHx < 1.0 ) .and.(RADIUS > 0.0) ) then
         TEFF   =   (1.0 - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
      else
         TEFF   = 0.0 ! -999.
      end if

      EVAP = a_eff*QL*DT*TEFF
      EVAP = MIN( EVAP , QL  )

      QC=QL+QI
      if (QC > 0.) then
         F    = F * ( QC - EVAP ) / QC
      end if

      QV   = QV   + EVAP
      QL   = QL   - EVAP
      TE   = TE   - (MAPL_ALHL/MAPL_CP)*EVAP

   end subroutine evap3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine subl3( &
         DT        , &
         RHCR      , &
         PL        , &
         TE        , &
         QV        , &
         QL        , &
         QI        , &
         F         , &
         XF        , & 
         NL        , &
         NI        , &
         QS        )

      real, intent(in   ) :: DT
      real, intent(in   ) :: RHCR
      real, intent(in   ) :: PL
      real, intent(inout) :: TE
      real, intent(inout) :: QV
      real, intent(inout) :: QL,QI
      real, intent(inout) :: F
      real, intent(in   ) :: XF
      real, intent(in   ) :: NL,NI
      real, intent(in   ) :: QS

      real :: ES,NN,RADIUS,K1,K2,TEFF,QCm,SUBL,RHx,QC !, QS

      real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
      real, parameter :: K_COND  =  2.4e-2        ! J m**-1 s**-1 K**-1
      real, parameter :: DIFFU   =  2.2e-5        ! m**2 s**-1

      real :: A_eff

      A_EFF = CLD_EVP_EFF

      NN = 5.*1.0e6

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!         SUBLORATION OF CLOUD WATER.             !!
      !!                                                 !!
      !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
      !!  formulation  (Eq.s 15-17)                      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !   QS  = QSAT(         &
      !               TE    , &
      !               PL      )

      ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100s <-^ convert from mbar to Pa)

      RHx = MIN( QV/QS , 1.00 )


      K1 = (MAPL_ALHL**2) * RHO_W / ( K_COND*MAPL_RVAP*(TE**2))

      K2 = MAPL_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Here DIFFU is given for 1000 mb  !!
      !! so 1000./PR accounts for inc-    !!
      !! reased diffusivity at lower      !!
      !! pressure.                        !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      if ( ( F > 0.) .and. ( QI > 0. ) ) then
         QCm=QI/F
      else
         QCm=0.
      end if

      RADIUS = LDRADIUS4(PL,TE,QCm,NN,RHX,NL,NI,2)
      
      if ( (RHx < RHCR ) .and.(RADIUS > 0.0) ) then
         TEFF   =   ( RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
      else
         TEFF   = 0.0 ! -999.
      end if

      SUBL = a_eff*QI*DT*TEFF
      SUBL = MIN( SUBL , QI  )

      QC=QL+QI
      if (QC > 0.) then
         F    = F * ( QC - SUBL ) / QC
      end if

      QV   = QV   + SUBL
      QI   = QI   - SUBL
      TE   = TE   - (MAPL_ALHS/MAPL_CP)*SUBL

   end subroutine subl3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

      integer :: NSMX

      real    :: ACF0, ACF, DTX
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


      NSMX = NSMAX
      DTX  = DT/NSMX

      CALL SUNDQ3_ICE3(TE, SUNDQV2, SUNDQV3, SUNDQT1, F2, F3 )

      C00x  = C_00 * F2 * F3
      !QCcrx = LWCRIT / ( F2 * F3 )
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
         
         
      integer :: NS, NSMX, itr,L

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
   subroutine SETTLE_VEL( WXR, QI, PL, TE, F, KH, VF, LARGESCALE, ANVIL, TROPP_Pa )

      real, intent(in   ) :: WXR 
      real, intent(in   ) :: TE
      real, intent(in   ) :: QI, F, PL
      real, intent(in   ) :: KH
      real, intent(out  ) :: VF

      real, intent(in) :: ANVIL, LARGESCALE
      real, intent(in) :: TROPP_Pa
      
      real :: RHO, XIm,LXIm, VF_A, VF_L, tpp_hPa, VF_PSC, wgt1, wgt2

! ------ For polar stratospheric clouds with single-moment microphysics ------

      REAL :: oneThird,ricecm,logsigicesq,ndensice,h2ocond,fluxcorr
      REAL :: rmedice,radius,rhoi,mdens,mfp,dynvis

      REAL, PARAMETER :: sigsq=1.3323e-19
      REAL, PARAMETER :: bet=1.458e-6
      REAL, PARAMETER :: s=110.4
      REAL, PARAMETER :: a=1.249
      REAL, PARAMETER :: b=0.42
      REAL, PARAMETER :: cc=0.87
      REAL, PARAMETER :: nice = 1.e-2
      REAL, PARAMETER :: sigice = 1.6 
      REAL, PARAMETER :: massh2o = 2.991e-23 
      REAL, PARAMETER :: densice = 1.

      REAL, PARAMETER :: BLEND_DEPTH_hPa = 50.
! ----------------------------------------------------------------------------

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

    ! Convective anvil
       VF_A = 128.6 + 53.2*LXIm + 5.5*LXIm**2

    ! Mid-latitude cirrus
       VF_L = 109.0*(XIm**0.16)
 
    ! Combine the two
       VF = ANVIL*VF_A + LARGESCALE*VF_L

    ! Convert from cm/s to m/s
       VF = 0.01 * VF

    ! Reduce/increase fall speeds for high/low pressure (NOT in LC98!!! ) 
    ! Assume unmodified they represent situation at 100 mb
      if (WXR > 0.) then
    !    VF = VF * ( 100./MAX(PL,10.) )**WXR
    !    VF = VF * SIN( 0.5*MAPL_PI*MIN(1.0,100./PL)**WXR )
         VF = VF * MIN(WXR,SIN( 0.5*MAPL_PI*MIN(1.0,100./PL)))
      endif

#ifdef DONT_SKIP_ICE_THICKEN
    ! More slowing at low levels in extratropics (Again NOT in LC98!!! ) 
      if ( KH > 2.0 ) then
        VF = VF * 0.01
      end if
#endif

#ifdef SKIP_PSC_FOR_NOW
! ------------------- Settling velocity for PSCs: -------------------

!    Change TROPP_Pa from Pa to hPa to match units of PL
      tpp_hPa = TROPP_Pa/100.  

!    If tropopause pressure is undefined set equal to 100. hPa
      IF (TROPP_Pa == MAPL_UNDEF) tpp_hPa = 100.

      IF( (PL < tpp_hPa) .AND.(QI > 0.) ) THEN

       mdens = (RHO/1000.)*MAPL_AVOGAD*1.00E-06/MAPL_AIRMW
       h2ocond = mdens*QI*MAPL_AIRMW/MAPL_H2OMW
       logsigicesq = LOG(sigice)*LOG(sigice)
       ndensice = densice/massh2o
       oneThird = 1./3.
       rmedice = (3.0*h2ocond/(ndensice*4.0*MAPL_PI*nice))**(oneThird)*EXP(-3.0/2.0*logsigicesq)

       IF(rmedice == 0.0 ) THEN

        VF_PSC = 0.0

       ELSE

! rmedice comes in as cm but need m so divide by 100
! densice is g/cm**3 but needs kg/m**3 so multiply by 1000
! -------------------------------------------------------------
        radius = 0.01*rmedice
        rhoi = densice*1000.
        mfp = .22508/(sigsq*mdens*1.e6)
        dynvis = bet*TE**1.5/(TE+s)
        fluxcorr = EXP(8.0*LOG(sigice)*LOG(sigice))

! took out divide by 100 so VF units m/s
! --------------------------------------
        VF_PSC = fluxcorr*(0.2222*rhoi*radius*radius*MAPL_GRAV/dynvis*(1.+ mfp/radius*(a+b*exp(-cc*radius/mfp))))

       END IF

      END IF

! --------- Modify VF with VF_PSC in the lower stratosphere ---------
     
      IF( (PL < tpp_hPa) .AND. (PL > (tpp_hPa-BLEND_DEPTH_hPa)) .AND. (QI > 0.) ) THEN
       wgt1 = ((tpp_hPa-BLEND_DEPTH_hPa)-PL)/(-BLEND_DEPTH_hPa)
       wgt2 = 1.-wgt1
      ELSE IF (PL <= (tpp_hPa-BLEND_DEPTH_hPa) .AND. (QI > 0.) ) THEN
       wgt1 = 0.0
       wgt2 = 1.0
      ELSE
      END IF
       

      IF( (PL < tpp_hPa) .AND. (QI > 0.) ) THEN
        VF = VF*wgt1+VF_PSC*wgt2
      ELSE
        VF = VF
      END IF

! -------------------------------------------------------------------
#endif

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
   function LDRADIUS3(PL,TE,QCL,NN) RESULT(RADIUS)

      real, intent(in) :: TE,PL,NN,QCL
      real :: RADIUS

      real :: MUU,RHO


      RHO = 100.*PL / (MAPL_RGAS*TE )
      MUU = QCL * RHO                       
      RADIUS = MUU/(NN*RHO_W*(4./3.)*MAPL_PI)
      RADIUS = RADIUS**(1./3.)    ! Equiv. Spherical Cloud Particle Radius in m


   end function LDRADIUS3
   function LDRADIUS4(PL,TE,QC,NN,RHX,NNL,NNI,ITYPE) RESULT(RADIUS)
   
       real, parameter :: betai     = -2.262e+3
       real, parameter :: gamai     =  5.113e+6
       real, parameter :: deltai    =  2.809e+3
       real, parameter :: densic    =  500.e0   !Ice crystal density in kgm-3

       real, parameter :: abeta = 0.07
       real, parameter :: bbeta = -0.14
       real, parameter :: bx = 100.* (3./(4.*MAPL_PI))**(1./3.)  
       real, parameter :: r13 = 1./3.  
       real, parameter :: r13bbeta = 1./3. - 0.14
       real, parameter :: premit= 750.e2            ! top pressure bound for mid level cloud
       real, parameter :: smnum = 1.e-1            !
       real, parameter :: rhminh= 0.80              ! minimum rh for high stable clouds      
       
       REAL   , INTENT(IN) :: TE,PL,NN,QC,RHX,NNL,NNI
       INTEGER, INTENT(IN) :: ITYPE
       REAL  :: RADIUS
       INTEGER, PARAMETER  :: CLOUD=1, ICE=2
       
       REAL :: NNX,RHO,IWL,BB, RIV,CLOUD_HIGH,RHDIF,WC
       INTEGER :: ICLD_HIGH
       
       !print*,"radius4=", PL,TE
       
       IF(ITYPE == CLOUD) THEN        
       !- liquid cloud effective radius ----- 
          !- [liu&daum, 2000 and 2005. liu et al 2008]
          !- air density (kg/m^3)
          RHO = 100.*PL / (MAPL_RGAS*TE )
   
          IF(USE_AEROSOL_NN) THEN
            !- cloud drop number concentration 
            !- from the aerosol model + ....
            NNX = NNL * 1.e-6  !#/cm3
          ELSE
            !- cloud drop number concentration :u[NNX]= #/cm^3
            NNX = NN * 1.e-6  !#/cm3
          ENDIF

          !- liquid water content : u[lwl] = g/m3
          WC = RHO*QC* 1.e+3  !g/m3
      
          !- radius in micrometers
          RADIUS= bx *  ( WC /NNX)**r13bbeta*abeta*6.92 !6.92=(1.e-6)**bbeta	      
   
          !- RADIUS is limited between 2.5 and 60 micrometers as 
          !- required by rrtm parameterization
          RADIUS = max(2.5, min( 60.0, RADIUS ) )
  
          !- convert to meter
          RADIUS = RADIUS*1.e-6
  
       ELSEIF(ITYPE == ICE) THEN

        IF (adjustl(CLDMICRO)=="GFDL") THEN

         RHO = 100.*PL / (MAPL_RGAS*TE )
         !- ice water content
         WC = RHO*QC  !kg/m3

         IF( (.not. USE_AEROSOL_NN) ) THEN 

            !------ice cloud effective radius ----- [klaus wyser, 1998]
             BB     =  -2. + log10(1000.*WC/50.)*(1.e-3*(273.15-TE)**1.5)
             BB     = MIN((MAX(BB,-6.)),-2.) 
             RADIUS =377.4 + 203.3 * bb+ 37.91 * bb **2 + 2.3696 * bb **3
             RADIUS =RADIUS * 1.e-6 !- convert to meter

         ELSE

            !--Mean volume and effective radius following Lohmann&Karcher (2002)
            !-- qice is the detrained ice water mixing ratio (kg/kg)
            !-- NNI  !#/m^3	 
            !-- RIV in micrometers

            RIV  = 1.E+6*((3.*WC)/(4.*MAPL_PI*densic*NNI))**0.33333  
            RIV  = MAX(RIV, 8.22)
            RADIUS= ((((RIV**3.-betai)**2.-gamai))/deltai)**0.33333
            
            !- convert to meter
            
            RADIUS = RADIUS*1.E-6

            !if((PL<300. .and. PL>200.) .and. (QC>1.e-6)) print*,"GOCRI=", RHO*QC,TE,NNI, RADIUS* 1.e+6
         ENDIF
  
        ELSE ! CLDMICRO =1MOMENT

            !------ice cloud effective radius ----- [klaus wyser, 1998]
            !- air density (kg/m^3)
            RHO = 100.*PL / (MAPL_RGAS*TE )
            !- ice water content
            iwl = RHO*QC* 1.e+3  !g/m3
            if(iwl<1.0e-6 .or. TE>273.0) then
             RADIUS = 5.0*1.e-6
            else
             BB     = -2. + log10(iwl/50.)*(1.e-3*(273.15-max(210.15,TE))**1.5)
             RADIUS =377.4 + 203.3 * bb+ 37.91 * bb **2 + 2.3696 * bb **3
             RADIUS =RADIUS * 1.e-6 !- convert to meter
             !print*,"bb=",temp,micro_g(ngrid)%rei(k,i,j),bb,iwl(k,i,j);call flush(6)
            endif

        ENDIF ! CLDMICRO

      ELSE
        STOP "WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2"
      ENDIF

   end function LDRADIUS4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef _CUDA
   attributes(device) &
#endif
   function ICE_FRACTION (TEMP,CNV_FRACTION,SNOMAS,FRLANDICE,FRLAND) RESULT(ICEFRCT)
      real, intent(in) :: TEMP,CNV_FRACTION,SNOMAS,FRLANDICE,FRLAND
      real             :: ICEFRCT
      real             :: ICEFRCT_0, ICEFRCT_1 

  ! Anvil-Convective sigmoidal function like figure 7(right) of Hu et al 2010 in tropical convective regimes
        ICEFRCT_0  = 0.00
        if ( TEMP <= aT_ICE_ALL ) then
           ICEFRCT_0 = 1.000
        else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
           ICEFRCT_0 = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
        end if
        ICEFRCT_0 = MIN(ICEFRCT_0,1.00)
        ICEFRCT_0 = MAX(ICEFRCT_0,0.00)
        ICEFRCT_0 = ICEFRCT_0**aICEFRPWR

  ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
      if ( (SNOMAS > 0.1) .OR. (FRLANDICE > 0.5) ) then
        ! Over snow/ice
        ICEFRCT_1  = 0.00
        if ( TEMP <= iT_ICE_ALL ) then
           ICEFRCT_1 = 1.000
        else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
           ICEFRCT_1 = 1.00 -  ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL )
        end if
        ICEFRCT_1 = MIN(ICEFRCT_1,1.00)
        ICEFRCT_1 = MAX(ICEFRCT_1,0.00)
        ICEFRCT_1 = ICEFRCT_1**iICEFRPWR
      else if (FRLAND > 0.1) then
        ! Over Land
        ICEFRCT_1  = 0.00
        if ( TEMP <= lT_ICE_ALL ) then
           ICEFRCT_1 = 1.000
        else if ( (TEMP > lT_ICE_ALL) .AND. (TEMP <= lT_ICE_MAX) ) then
           ICEFRCT_1 = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - lT_ICE_ALL ) / ( lT_ICE_MAX - lT_ICE_ALL ) ) )
        end if
        ICEFRCT_1 = MIN(ICEFRCT_1,1.00)
        ICEFRCT_1 = MAX(ICEFRCT_1,0.00)
        ICEFRCT_1 = ICEFRCT_1**lICEFRPWR
      else
        ! Over Oceans
        ICEFRCT_1  = 0.00
        if ( TEMP <= oT_ICE_ALL ) then
           ICEFRCT_1 = 1.000
        else if ( (TEMP > oT_ICE_ALL) .AND. (TEMP <= oT_ICE_MAX) ) then
           ICEFRCT_1 = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - oT_ICE_ALL ) / ( oT_ICE_MAX - oT_ICE_ALL ) ) )
        end if
        ICEFRCT_1 = MIN(ICEFRCT_1,1.00)
        ICEFRCT_1 = MAX(ICEFRCT_1,0.00)
        ICEFRCT_1 = ICEFRCT_1**oICEFRPWR
      endif

   ! Combine the Convective and Stratiform functions
      ICEFRCT  = ICEFRCT_1*(1.0-CNV_FRACTION) + ICEFRCT_0*(CNV_FRACTION)

   end function ICE_FRACTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
   real function ICEFRAC(T,T_TRANS,T_FREEZ)

      real, intent(in) :: T
      real, intent(in),optional :: T_TRANS
      real, intent(in),optional :: T_FREEZ

      real :: T_X,T_F

      if (present( T_TRANS )) then 
         T_X = T_TRANS
      else
         T_X = T_ICE_MAX
      endif
      if (present( T_FREEZ )) then 
         T_F = T_FREEZ
      else
         T_F = T_ICE_ALL
      endif


      if ( T < T_F ) ICEFRAC=1.000

      if ( T > T_X ) ICEFRAC=0.000

      if ( T <= T_X .and. T >= T_F ) then 
         ICEFRAC = 1.00 - ( T - T_F ) /( T_X - T_F )
      endif

   end function ICEFRAC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


      TE0=273.
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) &
#endif
   subroutine RADCOUPLE(  &
         TE,              & 
         PL,              & 
         DELP,            &
         KH,              &
         DTS,             &
         CF,              & 
         AF,              & 
         QV,              &
         QClLS,           & 
         QCiLS,           & 
         QClAN,           & 
         QCiAN,           & 
         QRN_ALL,         & 
         QSN_ALL,         & 
         NL,              &
         NI,              &
         RAD_QV,          &
         RAD_QL,          &  
         RAD_QI,          & 
         RAD_QR,          & 
         RAD_QS,          & 
         RAD_CF,          & 
         RAD_RL,          & 
         RAD_RI,          & 
         TEMPOR, FRLAND, CNV_FRACTION, SCLMFDFR, FR_AN_WAT, &
         CONVPAR_OPTION,  &
	 RHX)

      real, intent(in ) :: TE
      real, intent(in ) :: PL, DELP, KH, DTS
      real, intent(in ) :: AF,CF, QV, QClAN, QCiAN, QClLS, QCiLS
      real, intent(in ) :: QRN_ALL, QSN_ALL
      real, intent(in ) :: NL,NI
      real, intent(out) :: RAD_QV,RAD_QL,RAD_QI,RAD_QR,RAD_QS,RAD_CF,RAD_RL,RAD_RI

      real, intent(in )  :: tempor, FRLAND, CNV_FRACTION, SCLMFDFR,RHX
      integer, intent(in) :: FR_AN_WAT
      character(LEN=*), INTENT(IN) :: CONVPAR_OPTION
      real :: RElAN, REiAN, RElLS, REiLS, QCm, ss
      real :: QClANm, QCiANm, QClLSm, QCiLSm, QCtot, AFx
      real :: rampt, rampu, rampp

      real :: ALPH, CFPBL
      real :: CIP, TEM2D, TEM2, TEM3
      real :: NN, NN_LAND, NN_OCEAN

      ! Limits on Radii needed to ensure
      ! correct behavior of cloud optical
      ! properties currently calculated in 
      ! sorad and irrad (1e-6 m = micron)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adjust Anvil fractions for
      ! warm clouds
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (FR_AN_WAT == 0) then
        ALPH =  0.1
        SS   =  (280.-TE)/20.
        SS   =  MIN( 1.0 , SS )
        SS   =  MAX( 0.0 , SS )
        SS   =  ALPH + (SS**3) * ( 1.0 - ALPH )
        AFx  =  AF * SS * 0.5
      else
        AFx  =  AF
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Total cloud fraction
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RAD_CF = MIN( CF + AFx, 1.00 )

      ! Total In-cloud liquid
      if ( RAD_CF > 0. ) then
         RAD_QL = ( QClLS + QClAN ) / RAD_CF
      else
         RAD_QL = 0.0
      end if
      RAD_QL = MIN( RAD_QL, 0.01 )

      ! Total In-cloud ice
      if (  RAD_CF >0. ) then
         RAD_QI = ( QCiLS + QCiAN ) / RAD_CF
      else
         RAD_QI = 0.0
      end if
      RAD_QI = MIN( RAD_QI, 0.01 )


      ! Total In-cloud precipitation
      if (  RAD_CF >0. ) then
         RAD_QR = ( QRN_ALL ) / RAD_CF
         RAD_QS = ( QSN_ALL ) / RAD_CF
      else
         RAD_QR = 0.0
         RAD_QS = 0.0
      end if
      RAD_QL = MIN( RAD_QL, 0.01 )
      RAD_QI = MIN( RAD_QI, 0.01 )
      RAD_QR = MIN( RAD_QR, 0.01 )
      RAD_QS = MIN( RAD_QS, 0.01 )

     ! Number Concentration Assumptions
      NN_LAND  = 150.0e6
      NN_OCEAN =  30.0e6
     !          Over Land            Over Ocean
      NN = FRLAND*NN_LAND + (1.0-FRLAND)*NN_OCEAN

     ! LIQUID RADII
      !-BRAMS formulation     
      RAD_RL = LDRADIUS4(PL,TE,RAD_QL,NN,RHX,NL,NI,1)
     ! apply limits
      RAD_RL = RAD_RL*FAC_RL
      RAD_RL = MAX( MIN_RL, MIN(RAD_RL, MAX_RL) )

    ! ICE RADII
     !-BRAMS formulation  
      RAD_RI = LDRADIUS4(PL,TE,RAD_QI,NN,RHX,NL,NI,2)
    ! apply limits
      RAD_RI = RAD_RI*FAC_RI
      RAD_RI = MAX( MIN_RI, MIN(RAD_RI, MAX_RI) )

      if ( RAD_CF < 1.e-5 ) then
         RAD_QL = 0.
         RAD_QI = 0.
         RAD_CF = 0.
         RAD_QR = 0.
         RAD_QS = 0.
      end if

   end subroutine RADCOUPLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Parititions DQ into ice and liquid. Follows Barahona et al. GMD. 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine Bergeron_iter    (           &
         DTIME            , &
         PL               , &
         TE               , &
         QV               , &
         QILS             , &
         QICN             , &
         QLLS             , &
         QLCN             , &     
         CF               , &
         AF               , &
         NL               , &
         NI               , & 
         DQALL            , &
         FQI  ,    &
         CNV_FRACTION,SNOMAS,FRLANDICE,FRLAND, &
         needs_preexisting )

      real ,  intent(in   )    :: DTIME, PL, TE       !, RHCR
      real ,  intent(inout   )    ::  DQALL 
      real ,  intent(in)    :: QV, QLLS, QLCN, QICN, QILS
      real ,  intent(in)    :: CF, AF, NL, NI
      real, intent (out) :: FQI
      real, intent(in) :: CNV_FRACTION,SNOMAS,FRLANDICE,FRLAND
      logical, intent (in)  :: needs_preexisting
      
      real  :: DC, TEFF,QCm,DEP, &
            QC, QS, RHCR, DQSL, DQSI, QI, TC, &
            DIFF, DENAIR, DENICE, AUX, &
            DCF, QTOT, LHCORR,  QL, DQI, DQL, &
            QVINC, QSLIQ, CFALL,  new_QI, new_QL, &
            QSICE, fQI_0, QS_0, DQS_0, FQA, NIX

      DIFF = 0.0     
      DEP=0.0 
      QI = QILS + QICN !neccesary because NI is for convective and large scale 
      QL = QLLS +QLCN
      QTOT=QI+QL
      FQA = 0.0
      if (QTOT .gt. 0.0) FQA = (QICN+QILS)/QTOT
      NIX= (1.0-FQA)*NI

      DQALL=DQALL/DTIME                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      CFALL= min(CF+AF, 1.0)
      TC=TE-273.0
      fQI_0 = fQI

      !Completelely glaciated cloud:
      if (TE .ge. T_ICE_MAX) then   !liquid cloud
         FQI   = 0.0

      elseif(TE .le. T_ICE_ALL) then !ice cloud

         FQI   = 1.0

      else !mixed phase cloud

         FQI   = 0.0
         
          if (QILS .le. 0.0) then 
           
                    if (needs_preexisting) then
                   ! new 0518 this line ensures that only preexisting ice can grow by deposition.
                  ! Only works if explicit ice nucleation is available (2 moment muphysics and up)                        
                    else
                      fQi  =   ice_fraction( TE, CNV_FRACTION, SNOMAS, FRLANDICE, FRLAND )
                    end if                      
                  return 
         end if 
         
         
         QVINC=  QV 
         QSLIQ  = QSATLQ(         &
               TE   , &
               PL*100.0 , DQ=DQSL )

         QSICE  = QSATIC(         &
               TE   , &
               PL*100.0 , DQ=DQSI )

         QVINC =MIN(QVINC, QSLIQ) !limit to below water saturation 

         ! Calculate deposition onto preexisting ice 

         DIFF=(0.211*1013.25/(PL+0.1))*(((TE+0.1)/273.0)**1.94)*1e-4  !From Seinfeld and Pandis 2006
         DENAIR=PL*100.0/MAPL_RGAS/TE
         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
         LHcorr = ( 1.0 + DQSI*MAPL_ALHS/MAPL_CP) !must be ice deposition

         if  ((NIX .gt. 1.0) .and. (QILS .gt. 1.0e-10)) then 
            DC=max((QILS/(NIX*DENICE*MAPL_PI))**(0.333), 20.0e-6) !Assumme monodisperse size dsitribution 
         else
            DC = 20.0e-6
         end if

         TEFF= NIX*DENAIR*2.0*MAPL_PI*DIFF*DC/LHcorr ! 1/Dep time scale 

         DEP=0.0
         if ((TEFF .gt. 0.0) .and. (QILS .gt. 1.0e-14)) then 
            AUX =max(min(DTIME*TEFF, 20.0), 0.0)
            DEP=(QVINC-QSICE)*(1.0-EXP(-AUX))/DTIME
         end if

         DEP=MAX(DEP, -QILS/DTIME) !only existing ice can be sublimated

         !DEP=max(DEP, 0.0)

         DQI = 0.0
         DQL = 0.0
         FQI=0.0
         !QS_MIX=QSLIQ
         !DQS_MIX = DQSL
         !Partition DQALL accounting for Bergeron-Findensen process

         if  (DQALL .ge. 0.0) then !net condensation. Note: do not allow bergeron with QLCN

            if (DEP .gt. 0.0) then 
               DQI = min(DEP, DQALL + QLLS/DTIME)
               DQL = DQALL - DQI
            else
               DQL=DQALL ! could happen because the PDF allows condensation in subsaturated conditions
               DQI = 0.0 
            end if
         end if

         if  (DQALL .lt. 0.0) then  !net evaporation. Water evaporates first regaardless of DEP   
            DQL = max(DQALL, -QLLS/DTIME)   
            DQI = max(DQALL - DQL, -QILS/DTIME)        
         end if

         if (DQALL .ne. 0.0)  FQI=max(min(DQI/DQALL, 1.0), 0.0)

      end if !=====  

   end subroutine Bergeron_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   attributes(device) function QSAT(TL,PL,PASCALS)

     real,              intent(IN) :: TL, PL
     logical, optional, intent(IN) :: PASCALS
     real    :: QSAT

     real    :: URAMP, DD, QQ, TI, DQ, PP
     integer :: IT
 
     URAMP = TMIX

     if (present(PASCALS)) then 
        if (PASCALS) then
           PP = PL
        else
           PP = PL*100.
        end if
     else
        PP = PL*100.
     end if
 
     TI = TL - ZEROC

     if    (TI <= URAMP) then
        QSAT  =  QSATICE0(TL,PP,DQ)
     elseif(TI >= 0.0  ) then
        QSAT  =  QSATLQU0(TL,PP,DQ)
     else
        QSAT  =  QSATICE0(TL,PP,DQ)
        QQ    =  QSATLQU0(TL,PP,DQ)
        TI    =  TI/URAMP
        QSAT  =  TI*(QSAT - QQ) +  QQ
     end if

   end function QSAT

   attributes(device) function DQSAT(TL,PL,QSAT,PASCALS)

      real,              intent(IN) :: TL, PL
      real,              intent(OUT):: QSAT
      logical, optional, intent(IN ):: PASCALS
      real    :: DQSAT

      real    :: URAMP, TT, WW, DD, DQQ, QQ, TI, DQI, QI, PP, DQ
      integer :: IT

      URAMP = TMIX

      if (present(PASCALS)) then 
         if (PASCALS) then
            PP = PL
         else
            PP = PL*100.
         end if
      else
         PP = PL*100.
      end if

      TI = TL - ZEROC

      if    (TI <= URAMP) then
         QQ  = QSATICE0(TL,PP,DQ)
         QSAT  = QQ
         DQSAT = DQ
      elseif(TI >= 0.0  ) then
         QQ  = QSATLQU0(TL,PP,DQ)
         QSAT  = QQ
         DQSAT = DQ
      else
         QQ  = QSATLQU0(TL,PP,DQQ)
         QI  = QSATICE0(TL,PP,DQI)
         TI  = TI/URAMP
         DQSAT = TI*(DQI - DQQ) + DQQ
         QSAT  = TI*(QI - QQ) +  QQ
      end if

   end function DQSAT

   attributes(device) function QSATLQU0(TL,PL,DQ) result(QS)

      real, intent(IN) :: TL
      real, intent(IN) :: PL
      real, intent(OUT):: DQ
      real    :: QS

      real    :: TI,W
      real    :: DD
      real    :: TT
      real    :: DDQ
      integer :: IT

      integer, parameter :: TYPE = 1

#define TX TL
#define PX PL
#define EX QS
#define DX DQ


   if    (TX<TMINLQU) then
      TI = TMINLQU
   elseif(TX>TMAXTBL) then
      TI = TMAXTBL
   else
      TI = TX
   end if

#include "esatlqu.code"

   if    (TX<TMINLQU) then
      DDQ = 0.0
   elseif(TX>TMAXTBL) then
      DDQ = 0.0
   else
      if(PX>EX) then
         DD = EX
         TI = TX + DELTA_T
#include "esatlqu.code"
         DDQ = EX-DD
         EX  = DD
      endif
   end if

   if(PX > EX) then
      DD = ESFAC/(PX - (1.0-ESFAC)*EX)
      EX = EX*DD
      DX = DDQ*ERFAC*PX*DD*DD
   else
      EX = MAX_MIXING_RATIO
      DX = 0.0
   end if

#undef  DX
#undef  TX
#undef  EX
#undef  PX

      return
   end function QSATLQU0

   attributes(device) function QSATICE0(TL,PL,DQ) result(QS)

      real, intent(IN) :: TL
      real, intent(IN) :: PL
      real, intent(OUT):: DQ
      real    :: QS

      real    :: TI,W
      real    :: DD
      real    :: TT
      real    :: DDQ
      integer :: IT

      integer, parameter :: TYPE = 1

#define TX TL
#define PX PL
#define EX QS
#define DX DQ


   if    (TX<TMINICE) then
      TI = TMINICE
   elseif(TX>ZEROC  ) then
      TI = ZEROC
   else
      TI = TX
   end if

#include "esatice.code"

   if    (TX<TMINICE) then
      DDQ = 0.0
   elseif(TX>ZEROC  ) then
      DDQ = 0.0
   else
      if(PX>EX) then
         DD = EX
         TI = TX + DELTA_T
#include "esatice.code"
         DDQ = EX-DD
         EX  = DD
      endif
   end if

   if(PX > EX) then
      DD = ESFAC/(PX - (1.0-ESFAC)*EX)
      EX = EX*DD
      DX = DDQ*ERFAC*PX*DD*DD
   else
      EX = MAX_MIXING_RATIO
      DX = 0.0
   end if

#undef  DX
#undef  TX
#undef  EX
#undef  PX

         return
   end function QSATICE0

#endif

end module cloudnew
