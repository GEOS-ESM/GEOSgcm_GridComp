! $Id$
! $Name$

module cldmacro


   !This module handles large scale condesation and cloud fraction, convective precipitation, and makes 
   ! preliminary calculations for the two-moment cloud microphysics.  
   !=======================================================================

   use GEOS_UtilsMod,     only:QSAT=>GEOS_Qsat, DQSAT=>GEOS_DQsat, &
         QSATLQ=>GEOS_QsatLQU, QSATIC=>GEOS_QsatICE

   use MAPL_ConstantsMod, only: MAPL_TICE , MAPL_CP   , &
         MAPL_GRAV , MAPL_ALHS , &
         MAPL_ALHL , MAPL_ALHF , &
         MAPL_RGAS , MAPL_H2OMW, &
         MAPL_AIRMW, MAPL_RVAP , &
         MAPL_PI   , MAPL_R8   , &
         MAPL_R4

   use MAPL_BaseMod,      only: MAPL_UNDEF

   use GEOSmoist_Process_Library

   implicit none

   private

   public macro_cloud
   public update_cld
   public meltfrz_inst
   public CLDPARAMS, CLDPARAM_TYPE 
 
   type CLDPARAM_TYPE 
   real :: CNV_BETA
   real :: RH00
   real :: C_ACC
   real :: C_EV_R
   real :: C_EV_S
   real :: CCW_EVAP_EFF
   real :: CCI_EVAP_EFF
   real :: REVAP_OFF_P
   real :: CNVENVFC
   real :: T_ICE_ALL
   real :: CNVICEPARAM
   real :: CNVDDRFC
   integer :: PDFSHAPE
   real :: MINRHCRIT
   real :: TURNRHCRIT
   real :: TURNRHCRIT_UPPER
   real :: SLOPERHCRIT
   real :: DISP_FACTOR_ICE
   real :: DISP_FACTOR_LIQ
   real :: SCLM_DEEP, SCLM_SHALLOW
   endtype CLDPARAM_TYPE
   type (CLDPARAM_TYPE) :: CLDPARAMS

   real, parameter :: T_ICE_MAX    = MAPL_TICE  ! -7.0+MAPL_TICE
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
   real, parameter :: iT_ICE_ALL = 236.16
   real, parameter :: iT_ICE_MAX = 255.16
   real, parameter :: iICEFRPWR  = 6.0
        ! Over Land
   real, parameter :: lT_ICE_ALL = 239.16
   real, parameter :: lT_ICE_MAX = 261.16
   real, parameter :: lICEFRPWR  = 2.0
        ! Over Oceans
   real, parameter :: oT_ICE_ALL = 238.16
   real, parameter :: oT_ICE_MAX = 263.16
   real, parameter :: oICEFRPWR  = 4.0
   
   
   ! There are two PI's in this routine: PI_0 and MAPL_PI
   real, parameter :: PI_0 = 4.*atan(1.)
   logical, parameter ::  USE_AEROSOL_NN =  .TRUE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains 


   subroutine macro_cloud( &
!!! first vars are (in) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IRUN, LM         , &
         DT               , &
         PP_dev           , &
         PPE_dev          , &
         EXNP_dev         , &
         FRLAND_dev       , &
         CNVFRC_dev       , &
         SRFTYPE_dev      , &
         QLWDTR_dev       , &              
         QRN_CU_dev       , &
         CNV_UPDFRC_dev   , &
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
         PRECU_dev        , &
         CUARF_dev        , &
         SNRCU_dev        , &
         QST3_dev         , &
         DZET_dev         , &
         QDDF3_dev        , &
         RHX_dev          , &
         REV_AN_dev       , &
         RSU_AN_dev       , &
         ACLL_AN_dev,ACIL_AN_dev, &
         PFL_AN_dev,PFI_AN_dev, &
         PDFL_dev,PDFI_dev,FIXL_dev,FIXI_dev, &         
         DCNVL_dev, DCNVi_dev,      &
         ALPHT_dev,   &
         VFALLSN_AN_dev,  &
         VFALLRN_AN_dev,  &
         EVAPC_dev,  &
         SUBLC_dev,  &  
         SCICE_dev,  & 
         NCPL_dev,  &
         NCPI_dev,  &
         PFRZ_dev,  &
         DNDCNV_dev,  &
         DNCCNV_dev,  &
         RAS_DT_dev,  &
         QRAIN_CN, & 
         QSNOW_CN, &
         KCBL)


      integer, intent(in   )                    :: IRUN ! IM*JM
      integer, intent(in   )                    :: LM   ! LM
      real, intent(in   )                       :: DT   ! DT_MOIST
      real, intent(in   ), dimension(IRUN,  LM) :: PP_dev      ! PLO
      real, intent(in   ), dimension(IRUN,0:LM) :: PPE_dev     ! CNV_PLE
      real, intent(in   ), dimension(IRUN,  LM) :: EXNP_dev    ! PK
      real, intent(in   ), dimension(IRUN     ) :: FRLAND_dev  ! FRLAND
      real, intent(in   ), dimension(IRUN     ) :: CNVFRC_dev  ! CNVFRC
      real, intent(in   ), dimension(IRUN     ) :: SRFTYPE_dev
      real, intent(in   ), dimension(IRUN,  LM) :: QLWDTR_dev  ! CNV_DQLDT
      real, intent(inout), dimension(IRUN,  LM) :: QRN_CU_dev  ! CNV_PRC3 IS THIS INTENT IN?
      real, intent(inout), dimension(IRUN,  LM) :: CNV_UPDFRC_dev ! CNV_UPDF
      real, intent(in   ), dimension(IRUN,  LM) :: SC_QLWDTR_dev
      real, intent(in   ), dimension(IRUN,  LM) :: SC_QIWDTR_dev
      real, intent(inout), dimension(IRUN,  LM) :: QRN_SC_dev
      real, intent(inout), dimension(IRUN,  LM) :: QSN_SC_dev
      real, intent(inout), dimension(IRUN,  LM) :: SC_UPDFRC_dev
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
      real, intent(  out), dimension(IRUN     ) :: PRECU_dev ! CN_PRC2    
      real, intent(  out), dimension(IRUN     ) :: CUARF_dev ! CN_ARFX
      real, intent(  out), dimension(IRUN     ) :: SNRCU_dev ! CN_SNR
      real, intent(in   ), dimension(IRUN,  LM) :: QST3_dev   ! QST3
      real, intent(in   ), dimension(IRUN,  LM) :: DZET_dev   ! DZET
      real, intent(in   ), dimension(IRUN,  LM) :: QDDF3_dev  ! QDDF3
      real, intent(  out), dimension(IRUN,  LM) :: RHX_dev    ! RHX    
      real, intent(  out), dimension(IRUN,  LM) :: REV_AN_dev ! REV_CN
      real, intent(  out), dimension(IRUN,  LM) :: RSU_AN_dev ! RSU_CN
      real, intent(  out), dimension(IRUN,  LM) :: ACLL_AN_dev ! ACLL_CN
      real, intent(  out), dimension(IRUN,  LM) :: ACIL_AN_dev ! ACIL_CN     
      real, intent(  out), dimension(IRUN,0:LM) :: PFL_AN_dev ! PFL_CN
      real, intent(  out), dimension(IRUN,0:LM) :: PFI_AN_dev ! PFI_CN 
      real, intent(  out), dimension(IRUN,  LM) :: PDFL_dev ! DlPDF
      real, intent(  out), dimension(IRUN,  LM) :: PDFI_dev ! DiPDF
      real, intent(  out), dimension(IRUN,  LM) :: FIXL_dev ! DlFIX
      real, intent(  out), dimension(IRUN,  LM) :: FIXI_dev ! DiFIX                      
      real, intent(  out), dimension(IRUN,  LM) :: DCNVL_dev ! DCNVL
      real, intent(  out), dimension(IRUN,  LM) :: DCNVi_dev ! DCNVi
      real, intent(  out), dimension(IRUN,  LM) :: ALPHT_dev ! ALPHT
      real, intent(  out), dimension(IRUN,  LM) :: VFALLSN_AN_dev ! VFALLSN_CN
      real, intent(  out), dimension(IRUN,  LM) :: VFALLRN_AN_dev ! VFALLRN_CN
      real, intent(  out), dimension(IRUN,  LM) :: EVAPC_dev ! VFALLSN_CN
      real, intent(  out), dimension(IRUN,  LM) :: SUBLC_dev ! VFALLRN_CN
      
      !=====two_moment
      real, intent(inout), dimension(IRUN,  LM) :: SCICE_dev 
      real, intent(inout), dimension(IRUN,  LM) :: NCPL_dev
      real, intent(inout), dimension(IRUN,  LM) :: NCPI_dev
      real, intent(out), dimension(IRUN,  LM) :: PFRZ_dev
      real, intent(out), dimension(IRUN,  LM) :: DNDCNV_dev
      real, intent(out), dimension(IRUN,  LM) :: DNCCNV_dev
      real, intent(out), dimension(IRUN,  LM) :: RAS_DT_dev
      real, intent(out), dimension(IRUN,  LM) :: QRAIN_CN
      real, intent(out), dimension(IRUN,  LM) :: QSNOW_CN

      
      real, dimension(IRUN,  LM) :: FRZ_PP_dev ! FRZ_PP
      real :: TOT_UPDFRC
      integer, intent(in   ), dimension(IRUN     ) :: KCBL  ! RAS CLOUD BASE


      ! GPU The GPUs need to know how big local arrays are during compile-time
      !     as the GPUs cannot allocate memory themselves. This command resets
      !     this a priori size to LM for the CPU.


      integer :: I , K

      real :: MASS, iMASS
      real :: TOTFRC
      real :: QRN_CU_1D
      real :: QSN_CU
      real :: QRN_ALL, QSN_ALL
      real :: QTMP1, QTMP2, QTMP3, QTOT
      real :: TEMP
      real :: RHCRIT
      real :: AA3, BB3, ALPHA
      real :: VFALL, VFALLRN, VFALLSN
      real :: TOT_PREC_UPD   
      real :: AREA_UPD_PRC
      real :: AREA_UPD_PRC_tolayer
      real :: PRN_CU_above, PSN_CU_above
      real :: EVAP_DD_CU_above, SUBL_DD_CU_above

      real :: NIX, TOTAL_WATER, QRN_XS, QSN_XS

      logical :: use_autoconv_timescale

      use_autoconv_timescale = .false.
      QRN_XS = 0.0
      QSN_XS = 0.0

      RUN_LOOP: DO I = 1, IRUN

         K_LOOP: DO K = 1, LM         

            if (K == 1) then
               TOT_PREC_UPD = 0.
               AREA_UPD_PRC = 0.
            end if

            if (K == LM ) then
               !! ZERO DIAGNOSTIC OUTPUTS BEFORE SHOWERS !!
               PRECU_dev(I) = 0.
               SNRCU_dev(I) = 0. 
               CUARF_dev(I) = 0.
            end if

            !Zero out/initialize precips, except QRN_CU which comes from RAS        
            QRN_CU_1D = 0.
            QSN_CU    = 0.
            VFALL     = 0.

            PFL_AN_dev(I,K) = 0.
            PFI_AN_dev(I,K) = 0.
            IF (K == 1) THEN
               PFL_AN_dev(I,0) = 0.
               PFI_AN_dev(I,0) = 0.
            END IF

            ! Initialize other diagnostics 

            RHX_dev(I,K) = MAPL_UNDEF
            REV_AN_dev(I,K) = MAPL_UNDEF
            RSU_AN_dev(I,K) = MAPL_UNDEF
            ACLL_AN_dev(I,K) = MAPL_UNDEF
            ACIL_AN_dev(I,K) = MAPL_UNDEF
            PDFL_dev(I,K) = MAPL_UNDEF
            PDFI_dev(I,K) = MAPL_UNDEF
            FIXL_dev(I,K) = MAPL_UNDEF
            FIXI_dev(I,K) = MAPL_UNDEF
            DCNVL_dev(I,K) = MAPL_UNDEF
            DCNVi_dev(I,K) = MAPL_UNDEF
            ALPHT_dev(I,K) = MAPL_UNDEF
            VFALLSN_AN_dev(I,K) = MAPL_UNDEF
            VFALLRN_AN_dev(I,K) = MAPL_UNDEF
            
            EVAPC_dev(I,K) = 0.0
            SUBLC_dev(I,K) = 0.0
            
            !====two-moment

            DNDCNV_dev(I, K) = MAPL_UNDEF
            DNCCNV_dev(I, K) =  MAPL_UNDEF
            RAS_DT_dev(I, K) =  MAPL_UNDEF
            QRAIN_CN(I,K) = 0.0
            QSNOW_CN(I,K) = 0.0  
            NIX= 0.0

            ! Copy QRN_CU into a temp scalar
            QRN_CU_1D = QRN_CU_dev(I,K)

            MASS =  ( PPE_dev(I,K) - PPE_dev(I,K-1) )*100./MAPL_GRAV  ! layer-mass (kg/m**2)
            iMASS = 1.0 / MASS
            TEMP =  EXNP_dev(I,K) * TH_dev(I,K) 
            FRZ_PP_dev(I,K) = 0.00
            !PFRZ_dev(I, K) = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Condensate Source
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            TOTAL_WATER =  (QIW_AN_dev(I,K)+QLW_AN_dev(I,K) + QIW_LS_dev(I,K)+ QIW_LS_dev(I,K))*MASS +QLWDTR_dev(I,K)*DT+SC_QLWDTR_dev(I,K)*DT+SC_QIWDTR_dev(I,K)*DT

            DCNVi_dev(I,K) = QIW_AN_dev(I,K)
            DCNVL_dev(I,K) = QLW_AN_dev(I,K)
            DNDCNV_dev(I, K) = NCPL_dev(I, K)
            DNCCNV_dev(I, K) = NCPI_dev(I, K)

           ! cnvsrc now handled in convection routines

            DCNVi_dev(I,K) = ( QIW_AN_dev(I,K) - DCNVi_dev(I,K) ) / DT
            DCNVL_dev(I,K) = ( QLW_AN_dev(I,K) - DCNVL_dev(I,K) ) / DT
            DNDCNV_dev(I, K) = (NCPL_dev(I, K)-DNDCNV_dev(I, K))/DT
            DNCCNV_dev(I, K) = (NCPI_dev(I, K)-DNCCNV_dev(I, K))/DT




!!!!!!!!!!!!!!!!!!!!check consistency!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

            ! assume deep and shallow updraft fractions non-overlapping
            TOT_UPDFRC = CNV_UPDFRC_dev(I,K) + SC_UPDFRC_dev(I,K)
            TOT_UPDFRC = MAX( MIN( TOT_UPDFRC, 1.), 0.)

            call pdf_spread (&
                  PP_dev(I,K),ALPHA,&
                  ALPHT_dev(I,K), & 
                  FRLAND_dev(I))

            ! impose a minimum amount of variability
            ALPHA    = MAX(  ALPHA , 1.0 - CLDPARAMS%RH00 )

            RHCRIT = 1.0 - ALPHA


            !=================================new condensate ====================================
!!!!!!!!!Calculate probability of freezing to scale nucleated ice crystals !!
            !================================   


            call Pfreezing (  &       
                  CLDPARAMS%PDFSHAPE , &
                  ALPHA              , &
                  PP_dev(I,K)        , &
                  TEMP               , &
                  Q_dev(I,K)         , &
                  QLW_LS_dev(I,K)    , &
                  QLW_AN_dev(I,K)    , &
                  QIW_LS_dev(I,K)    , &
                  QIW_AN_dev(I,K)    , &    
                  SCICE_dev(I, K)    , &
                  CLDFRC_dev(I,K)    , & 
                  ANVFRC_dev(I,K)    , &
                  PFRZ_dev(I, K) )   



            !=============Collect convective precip==============


           PDFL_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)
           PDFI_dev(I,K) = QIW_LS_dev(I,K)+QIW_AN_dev(I,K)
  
            call hystpdf_new(      &
                  DT             , &
                  ALPHA          , &
                  CLDPARAMS%PDFSHAPE , &
                  CNVFRC_dev(I)  , &
                  SRFTYPE_dev(I)  , &
                  PP_dev(I,K)    , &
                  Q_dev(I,K)     , &
                  QLW_LS_dev(I,K), &
                  QLW_AN_dev(I,K), &
                  QIW_LS_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  TEMP           , &
                  CLDFRC_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  NCPL_dev(I,K)  , &
                  NCPI_dev(I,K)  , &
                  SCICE_dev(I, K))
            
            PDFL_dev(I,K)  = ( QLW_LS_dev(I,K) + QLW_AN_dev(I,K) - PDFL_dev(I,K) ) / DT 
            PDFI_dev(I,K)  = ( QIW_LS_dev(I,K) + QIW_AN_dev(I,K) - PDFI_dev(I,K) ) / DT 



            QTMP1 = QLW_LS_dev(I,K) + QLW_AN_dev(I,K)
            QTMP2 = QIW_LS_dev(I,K) + QIW_AN_dev(I,K)


            EVAPC_dev(I,K) = QLW_LS_dev(I,K)+QLW_AN_dev(I,K)
            SUBLC_dev(I,K) = QIW_LS_dev(I,K)+QIW_AN_dev(I,K)

             ! 'Anvil' partition from RAS/Parameterized not done in hystpdf

            call evap3(            &
                  DT             , &
                  CLDPARAMS%CCW_EVAP_EFF, &
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  NCPL_dev(I,K) , &
                  NCPI_dev(I,K) , &
                  QST3_dev(I,K)  )  

            call subl3(            &
                  DT             , & 
                  CLDPARAMS%CCI_EVAP_EFF, &
                  RHCRIT         , &
                  PP_dev(I,K)    , &
                  TEMP           , &
                  Q_dev(I,K)     , &
                  QLW_AN_dev(I,K), &
                  QIW_AN_dev(I,K), &
                  ANVFRC_dev(I,K), &
                  NCPL_dev(I,K) , &
                  NCPI_dev(I,K) , &
                  QST3_dev(I,K)  ) 

            EVAPC_dev(I,K) = ( EVAPC_dev(I,K) - (QLW_LS_dev(I,K)+QLW_AN_dev(I,K)) ) / DT
            SUBLC_dev(I,K) = ( SUBLC_dev(I,K) - (QIW_LS_dev(I,K)+QIW_AN_dev(I,K)) ) / DT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !  Add in convective rain 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! CU-FREEZE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Also "freeze" out any conv. precip that needs
            ! to be since this isn't done in RAS. This is
            ! precip w/ large particles, so freezing is 
            ! strict. Check up on this!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP1 = 0.
            QTMP2 = 0.
            QTMP3 = 0.
            QRN_ALL = 0.
            QSN_ALL = 0.


            if ( TEMP < MAPL_TICE ) then	    
               QTMP2     = QRN_CU_1D
               QSN_CU    = QRN_CU_1D
               QRN_CU_1D = 0.
               TEMP      = TEMP + QSN_CU*(MAPL_ALHS-MAPL_ALHL) / MAPL_CP
            end if

            QRN_CU_1D = QRN_CU_1D + QRN_SC_dev(I,K) + QRN_XS!! add any excess precip to convective 
            QSN_CU    = QSN_CU    + QSN_SC_dev(I,K) + QSN_XS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !----------------------------------------------------------------------------------------------
            ! Column will now be swept from top-down for precip accumulation/accretion/re-evaporation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            AREA_UPD_PRC_tolayer = 0.0


            TOT_PREC_UPD  = TOT_PREC_UPD + ( ( QRN_CU_1D + QSN_CU ) * MASS )
            AREA_UPD_PRC  = AREA_UPD_PRC + ( TOT_UPDFRC* ( QRN_CU_1D + QSN_CU )* MASS )

            if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC_tolayer = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )

            AREA_UPD_PRC_tolayer = CLDPARAMS%CNV_BETA * AREA_UPD_PRC_tolayer

            IF (K == LM) THEN ! We've accumulated over the whole column

               if ( TOT_PREC_UPD > 0.0 ) AREA_UPD_PRC = MAX( AREA_UPD_PRC/TOT_PREC_UPD, 1.E-6 )

               AREA_UPD_PRC = CLDPARAMS%CNV_BETA * AREA_UPD_PRC

               !! "couple" to diagnostic areal fraction output 
               !! Intensity factor in PRECIP3 is floored at
               !! 1.0. So this is fair.

               CUARF_dev(I) = MIN( AREA_UPD_PRC, 1.0 )

            END IF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! GET SOME MICROPHYSICAL QUANTITIES 

            CALL MICRO_AA_BB_3( TEMP,PP_dev(I,K),QST3_dev(I,K),AA3,BB3 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QTMP1 = QLW_LS_dev(I,K) + QLW_AN_dev(I,K)
            QTMP2 = QIW_LS_dev(I,K) + QIW_AN_dev(I,K)
            QTOT=QTMP1+QTMP2

            ! QTMP1 = 0.0
            ! QTMP2 = 0.0



            ! Convective
            ! ----------
            !RHCRIT=1.0

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
                  REV_AN_dev(I,K) , &
                  RSU_AN_dev(I,K) , &
                  ACLL_AN_dev(I,K), &
                  ACIL_AN_dev(I,K), &
                  PFL_AN_dev(I,K) , &
                  PFI_AN_dev(I,K) , &
                  VFALLRN         , &
                  VFALLSN         , &
                  FRZ_PP_dev(I,K) , &
                  CLDPARAMS%CNVENVFC, CLDPARAMS%CNVDDRFC, &
                  ANVFRC_dev(I,k), CLDFRC_dev(I,k),  &
                  PP_dev(I,KCBL(I)))

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

            if (.true.) then 

               IF ( (QLW_LS_dev(I,K)+QLW_AN_dev(I,K)) > 1.e-20 ) THEN
                  QTMP3 = 1./(QLW_LS_dev(I,K)+QLW_AN_dev(I,K))
               ELSE
                  QTMP3 = 0.0
               END IF
               QLW_LS_dev(I,K) = QLW_LS_dev(I,K) * QTMP1 * QTMP3
               QLW_AN_dev(I,K) = QLW_AN_dev(I,K) * QTMP1 * QTMP3
               NCPL_dev(I, K) = NCPL_dev(I, K)* QTMP1 * QTMP3

               IF ( (QIW_LS_dev(I,K)+QIW_AN_dev(I,K)) > 1.0e-20 ) THEN
                  QTMP3 = 1./(QIW_LS_dev(I,K)+QIW_AN_dev(I,K))
               ELSE
                  QTMP3 = 0.0
               END IF
               QIW_LS_dev(I,K) = QIW_LS_dev(I,K) * QTMP2 * QTMP3
               QIW_AN_dev(I,K) = QIW_AN_dev(I,K) * QTMP2 * QTMP3
               NCPI_dev(I, K) = NCPI_dev(I, K)* QTMP2 * QTMP3

               ! reduce cloud farction as well
               QTMP3 = QIW_LS_dev(I,K)+QIW_AN_dev(I,K) + QLW_LS_dev(I,K)+QLW_AN_dev(I,K)

               If (QTOT .gt. 0.0) then 
                  CLDFRC_dev(I,k) = CLDFRC_dev(I,k)*QTMP3/QTOT          
                  ANVFRC_dev(I,k) = ANVFRC_dev(I,k)*QTMP3/QTOT     
               end if

            end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            QRAIN_CN(I,K) = QRN_ALL / (100.*PP_dev(I,K) / (MAPL_RGAS*TEMP ))
            QSNOW_CN(I,K) = QSN_ALL / (100.*PP_dev(I,K) / (MAPL_RGAS*TEMP ))
            QRN_CU_dev(I,K) = QRN_CU_1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            TOTFRC = CLDFRC_dev(I,K) + ANVFRC_dev(I,K)

            IF ( TOTFRC > 1.00 ) THEN
               CLDFRC_dev(I,k) = CLDFRC_dev(I,k)*(1.00 / TOTFRC )
               ANVFRC_dev(I,k) = ANVFRC_dev(I,k)*(1.00 / TOTFRC )
            END IF

            TOTFRC = CLDFRC_dev(I,K) + ANVFRC_dev(I,K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



            TH_dev(I,K)  =  TEMP / EXNP_dev(I,K) 

         end do K_LOOP


      end do RUN_LOOP

   END SUBROUTINE MACRO_CLOUD


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

   subroutine pdf_spread ( &
         PP,ALPHA,&
         ALPHT_DIAG, &
         FRLAND )  

      real,    intent(in)  :: PP
      real,    intent(out) :: ALPHA
      real,    intent(out) :: ALPHT_DIAG
      real,    intent(in)  :: FRLAND

      real :: slope_up, aux1, aux2, maxalpha

      slope_up = 20.0
      maxalpha=1.0-CLDPARAMS%MINRHCRIT

      ! alpha is the 1/2*width so RH_crit=1.0-alpha

      !  Use Slingo-Ritter (1985) formulation for critical relative humidity
      !  array a1 holds the critical rh, ranges from 0.8 to 1
      !Reformulated by Donifan Barahona

      aux1 = min(max((pp- CLDPARAMS%TURNRHCRIT)/CLDPARAMS%SLOPERHCRIT, -20.0), 20.0) 
      aux2 = min(max((CLDPARAMS%TURNRHCRIT_UPPER - pp)/slope_up, -20.0), 20.0)


      if (frland > 0.05)  then           
         aux1=1.0
         !maxalpha=max(maxalpha-0.05, 0.001)
      else
         aux1 = 1.0/(1.0+exp(aux1)) !this function reproduces the old Sligo function. 
       !  aux2=min(max(2.0*(ltsx-min_lts), -20.0), 20.0)
         !aux2=0.5/(1.0+exp(aux2))
        ! aux1=max(aux2, aux1)
          
      end if

      !aux2= 1.0/(1.0+exp(aux2)) !this function reverses the profile at low P    
      aux2=1.0

      alpha  = maxalpha*aux1*aux2



      ALPHA = MIN( ALPHA , 0.4 )  ! restrict RHcrit to > 60% 
      ALPHT_DIAG = ALPHA

   end subroutine pdf_spread

   subroutine update_cld( &
         DT          , &
         ALPHA       , &
         PDFFLAG     , &
         CNVFRC      , &
         SRFTYPE     , &
         PL          , &
         QV          , &
         QCl         , &
         QAl         , &
         QCi         , &
         QAi         , &
         TE          , &
         CF          , &
         AF          , &
         SCICE       , &
         NI          , &
         NL          , &
         RHcmicro    , &
         DO_HYSTPDF)

      real, intent(in)    :: DT,ALPHA,PL,CNVFRC,SRFTYPE
      integer, intent(in) :: pdfflag
      real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF, NI, RHCmicro, NL,  SCICE

      ! internal arrays
      real :: CFO
      real :: QT

      real :: QSx,DQsx

      real :: QCx, QC, QA

      real :: QX, QSLIQ, QSICE, CFALL, DQx, FQA, DELQ

      real :: SHOM, maxalpha
      ! internal scalars
      logical :: DO_HYSTPDF

      maxalpha=1.0-CLDPARAMS%MINRHCRIT 

      QC = QCl + QCi
      QA = QAl + QAi
      QT  =  QC  + QA + QV  !Total water after microphysics
      CFALL  = AF+CF
      FQA = 0.0
      if (QA+QC .gt. tiny(1.0))  FQA=QA/(QA+QC)

      SHOM=2.349-(TE/259.0) !hom threeshold Si according to Ren & McKenzie, 2005

      !================================================
      ! First find the cloud fraction that would correspond to the current condensate 
      QSLIQ  = QSATLQ(         &
            TE   , &
            PL*100.0 , DQ=DQx )


      QSICE  = QSATIC(         &
            TE   , &
            PL*100.0 , DQ=DQx )

      if ((QC+QA) .gt. 1.0e-13) then 
         QSx=((QCl+QAl)*QSLIQ + QSICE*(QCi+QAi))/(QC+QA)               
      else
         DQSx  = DQSAT(         &
               TE   , &
               PL ,  35.0, QSAT=QSx     ) !use ramp to -40 
      end if


      
      if (TE .gt. CLDPARAMS%T_ICE_ALL)   SCICE = 1.0
      QCx=QC+QA   
      QX=QT-QSx*SCICE
      CFo=0.
      
      !====== recalculate QX if too low and SCICE<SHOM
        if ((QX .gt. QCx) .and. (QCx .gt. 0.0)) then     
           QX=QT-QSx*SHOM
       end if 
      
      !=======================

     DELQ=max(min(2.0*maxalpha*QSx, 0.5*QT), 1.0e-12)   
      
      if  ((QX .le. QCx)  .and. (QCx .gt. tiny(1.0)))  then          
         CFo =  (1.0+SQRT(1.0-(QX/QCx)))
         if (CFo .gt. 1.e-6) then
            CFo = min(1.0/CFo, 1.0)
            DELQ=  2.0*QCx/(CFo*CFo)
         else
            CFo = 0.0
         end if
      elseif (QCx .gt. tiny(1.0)) then  
         !   CFo = 1.0  !Outside of distribution but still with condensate         
        DELQ=max(min(2.0*maxalpha*QSx, 0.5*QT), 1.0e-12)          
        CFo = SQRT(2.0*QCx/DELQ)        
      else
        CFo = 0.0         
      end if

      if  (QSx .gt. tiny(1.0)) then 
         RHCmicro = SCICE - 0.5*DELQ/Qsx
      else
         RHCmicro = 0.0
      end if

      CFALL   = max(CFo, 0.0) 
      CFALL   = min(CFo, 1.0) 

      CF=CFALL*(1.0-FQA)
      AF=CFALL*FQA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.NOT. DO_HYSTPDF) return 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      !if ((TE .le. CLDPARAMS%T_ICE_ALL))  return !don't do anything else for cirrus


      !================================================
      !================================================
      ! Now find the equilibirum cloud fraction for mix-phase and liquid clouds. Create new condensate if neccesary

                  
                  
       call hystpdf_new(  &
            DT          , &
            ALPHA       , &
            PDFFLAG     , &
            CNVFRC      , &
            SRFTYPE     , &
            PL          , &
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
            SCICE  )


    CALL fix_up_clouds(    &
                  QV     , &
                  TE           , &
                  QCl, & 
                  QCi, & 
                  CF, &   
                  QAl, & 
                  QAi, & 
                  AF)



   end subroutine update_cld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hystpdf_new( &
         DT          , &
         ALPHA       , &
         PDFFLAG     , &
         CNVFRC      , &
         SRFTYPE     , &
         PL          , &
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
         SC_ICE)

      real, intent(in)    :: DT,ALPHA,PL,CNVFRC,SRFTYPE
      integer, intent(in) :: pdfflag
      real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF
      real, intent(in)    :: NL,NI
      real, intent(in)    :: SC_ICE
      
      
      ! internal arrays
      real :: QCO, QVO, CFO, QAO
      real :: QT, sigmaqt1, sigmaqt2

      real :: TEO,QSx,DQsx,DQs, qsnx

      real :: TEp, CFp, QVp, QCp
      real :: TEn, QSn, CFn, QVn, QCn

      real :: QCx, QVx, CFx, QAx, QC, QA, fQi
      real :: dQAi, dQAl, dQCi, dQCl, Nfac, NLv, NIv 

      real :: fQip

      real :: tmpARR
      real :: ALHX, DQCALL
      ! internal scalars
      integer :: N, nmax

      QC = QCl + QCi
      QA = QAl + QAi
      QT  =  QC  + QA + QV  !Total water after microphysics
      tmpARR = 0.0
      nmax =  20
      QAx = 0.0

      if ( AF < 1.0 )  tmpARR = 1./(1.-AF)

      TEo = TE

      fQi = ice_fraction( TE, CNVFRC,SRFTYPE )
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
         fQip= fQi

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


         qsnx= qsn*SC_ICE !
         if ((QCI .ge. 0.0) .and. (qsn .gt. qt))  qsnx=qsn !this way we do not evaporate preexisting ice but maintain supersat
         
         
         call pdffrac(PDFFLAG,qt,sigmaqt1,sigmaqt2,qsnx,CFn)
         call pdfcondensate(PDFFLAG,qt,sigmaqt1,sigmaqt2,qsnx,QCn)

       
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
                 CNVFRC,SRFTYPE   , &
                 DQCALL           , &
                 fQi              , & 
                 .true.)
       

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
            if (n.ne.nmax) QCn = QCp + ( QCn - QCp ) *0.5
         endif

         QVn = QVp - (QCn - QCp)
         TEn = TEp + (1.0-fQi)*(MAPL_ALHL/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
               +      fQi* (MAPL_ALHS/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

         if (abs(Ten - Tep) .lt. 0.00001) exit 

         DQS  = DQSAT( TEn, PL, QSAT=QSn )

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
         ENVFC,DDRFC, AF, CF, &
         PCBL  )


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

      real, intent(in   ) :: AF,CF, PCBL


      real :: PFi,PFl,ENVFRAC
      real :: TKo,QKo,QSTKo,DQSTKo,RH_BOX,T_ED,QPlKo,QPiKo
      real :: Ifactor,RAINRAT0,SNOWRAT0
      real :: FALLRN,FALLSN,VEsn,VErn,NRAIN,NSNOW,Efactor

      real :: TinLAYERrn,DIAMrn,DROPRAD
      real :: TinLAYERsn,DIAMsn,FLAKRAD

      real :: EVAP,SUBL,ACCR,MLTFRZ,EVAPx,SUBLx
      real :: EVAP_DD,SUBL_DD,DDFRACT
      real :: LANDSEAF, TC, MAXMLT, iDT

      real :: tmpARR, CFR, aux, RH_EVAP

    
      real :: QSICE, DQSI, Told, QKCLR

      integer :: itr

      logical, parameter :: taneff = .true.
      
      
      
       real, parameter :: TRMV_L  = 1.0         ! m/s
      real, parameter :: TAU_FRZ = 5000.0      ! sec       
      real, parameter :: FRZ_TAU = 1.0/TAU_FRZ ! sec^-1
      real, parameter :: MELT_T  = 5.0         ! degrees C
      real, parameter :: LFBYCP  = MAPL_ALHF/MAPL_CP
      real, parameter :: CPBYLF  = 1.0/LFBYCP
      real, parameter :: B_SUB   = 1.00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! fraction of precip falling through "environment" vs
      ! through cloud
     

      if(taneff) then
         !reproduces the atan profile but is less messy

         aux = min(max((pl- PCBL)/10.0, -20.0), 20.0) 
         aux = 1.0/(1.0+exp(-aux))
         envfrac = ENVFC + (1.0-ENVFC)*aux   !ENVFC is the minimum exposed area. Below cloud base envfrac becomes 1.    


         !if (pl .le. 600.) then
         !   envfrac = 0.25
         !else
         !   envfrac = 0.25 + (1.-0.25)/(19.) *                    &
         !         ((atan( (2.*(pl-600.)/(900.-600.)-1.) *       &
         !         tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
         ! end if


         envfrac = min(envfrac,1.)

      else
         ENVFRAC = ENVFC
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CFR= AF+CF
      if ( CFR < 0.99) then 
         tmpARR = 1./(1.-CFR)
      else
         tmpARR = 0.0
      end if
!!!!!!!!!!!!!!!!!!!

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
      PFL_DIAG =  0.
      PFI_DIAG =  0.
      ACRIL_DIAG    =  0.
      ACRLL_DIAG    =  0.
      REVAP_DIAG    =  0.
      RSUBL_DIAG    =  0.
      RH_EVAP= RHCR3
      
      !RH_EVAP= 1.0

      DDFRACT = DDRFC 

      IF (K == 1) THEN
         PFl=QPl*MASS
         PFi=QPi*MASS

         EVAP_DD = 0.
         SUBL_DD = 0.

         VFALLRN = 0.0
         VFALLSN = 0.0
      ELSE 
         QPl   = QPl + PFl_above * iMASS
         PFl = 0.00

         QPi   = QPi + PFi_above * iMASS
         PFi = 0.00


          IF(QCL > 0.0) THEN
            IF(QPi > 0.0) THEN
               ACCR = min(CLDPARAMS%C_ACC*(QPl*MASS)*QCl, QCl)
               QPl  = QPl + ACCR
               QCl  = QCl - ACCR

               ACRLL_DIAG = ACCR * iDT
            END IF

            IF(QPi > 0.0) THEN
               ACCR = min(CLDPARAMS%C_ACC*(QPi*MASS)*QCl, QCl)
               QPi  = QPi + ACCR
               QCl  = QCl - ACCR
               TE   = TE + LFBYCP*ACCR

               ACRIL_DIAG = ACCR * iDT
            END IF
         END IF
        

         RAINRAT0 = Ifactor*QPl*MASS/DT
         SNOWRAT0 = Ifactor*QPi*MASS/DT

         call MARSHPALMQ2(RAINRAT0,PL,DIAMrn,NRAIN,FALLrn,VErn)
         call MARSHPALMQ2(SNOWRAT0,PL,DIAMsn,NSNOW,FALLsn,VEsn)

         IF ( FRLAND < 0.1 ) THEN
            !!      DIAMsn = MAX(  DIAMsn, 1.0e-3 )   ! Over Ocean
         END IF

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

         EVAP = 0.0
         SUBL = 0.0

         !  if (TKo .gt. 240.0) then 
         do itr = 1,20 ! 

            DQSTKo = DQSAT  ( TKo , PL, QSAT=QSTko   ) !use for rain

            QSTKo  = MAX( QSTKo , 1.0e-7 )        

!!!!! RAin falling !!!!!!!!!!!!!!!!!!!!!!!
            if (tmpARR .gt. 0.0) then 
               QKCLR=(QKo -QSTKo*CFR)*tmpARR
               RH_BOX =QKCLR/QSTKo 
            else
               RH_BOX = QKo/QSTKo
            end if

            IF ( RH_BOX < RH_EVAP ) THEN
               Efactor =  RHO_W * ( AA + BB )    / (RH_EVAP - RH_BOX )
            else
               Efactor = 9.99e9
            end if


            LANDSEAF = 1.00


            if ( ( RH_BOX < RH_EVAP ) .AND. ( DIAMrn > 0.00 ) .AND. &
                  ( PL > 100. ) .AND. ( PL < CLDPARAMS%REVAP_OFF_P ) ) then
               DROPRAD=0.5*DIAMrn
               T_ED =  Efactor * DROPRAD**2 
               T_ED =  T_ED * ( 1.0 + DQSTKo*MAPL_ALHL/MAPL_CP )
               EVAP =  QPl*(1.0 - EXP( -CLDPARAMS%C_EV_R * VErn * LANDSEAF *ENVFRAC* TinLAYERrn / T_ED ) )
            ELSE
               EVAP = 0.0
            END if

!!!!! Snow falling !!!!!!!!!!!!!!!!!!!!!!! 

            !QSICE  = QSATIC(  min(TKo, T_ICE_MAX), PL*100.0 , DQ=DQSI ) ! use for snow

            DQSI = DQSAT( TKo , PL, 5.0, QSAT = QSICE ) !use for snow, small ramp to assure continuitiy at higher T 
            !DQSI = DQSAT( TKo , PL, QSAT = QSICE ) 
            QSICE  = MAX( QSICE , 1.0e-7 )
            if (tmpARR .gt. 0.0) then 
               QKCLR =(QKo -QSICE*CFR)*tmpARR !Snow only sublimates when QV<QICE  
               RH_BOX =QKCLR/QSICE 
            else
               RH_BOX =QKo/QSICE                                        
            end if



            IF ( RH_BOX < RH_EVAP ) THEN
               Efactor =  0.5*RHO_W * ( AA + BB )    / (RH_EVAP - RH_BOX )
            else
               Efactor = 9.99e9
            end if


            if ( ( RH_BOX < RH_EVAP ) .AND. ( DIAMsn > 0.00 ) .AND. &
                  ( PL > 100. ) .AND. ( PL < CLDPARAMS%REVAP_OFF_P ) ) then
               FLAKRAD=0.5*DIAMsn
               T_ED =  Efactor * FLAKRAD**2   
               T_ED =  T_ED * ( 1.0 + DQSI*MAPL_ALHS/MAPL_CP )
               SUBL =  QPi*(1.0 - EXP( -CLDPARAMS%C_EV_S * VEsn * LANDSEAF * ENVFRAC * TinLAYERsn / T_ED ) )
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

            EVAP= EVAP*(1.-CFR)! Can only evaporate in the clear parto of the cell DONIF
            SUBL = SUBL*(1.-CFR)

            Told = TKo
        	    !QKo=QKo + EVAP + SUBL
            TKo=TKo - EVAP * MAPL_ALHL / MAPL_CP - SUBL * MAPL_ALHS / MAPL_CP


            if (abs(Told-Tko) .le. 0.01) exit
         enddo
         ! end if       

         QPi  = QPi - SUBL
         QPl  = QPl - EVAP

         !! Put some re-evap/re-subl precip in to a \quote{downdraft} to be applied later
         EVAP_DD = EVAP_DD_above + DDFRACT*EVAP*MASS 
         EVAP    = EVAP          - DDFRACT*EVAP
         SUBL_DD = SUBL_DD_above + DDFRACT*SUBL*MASS 
         SUBL    = SUBL          - DDFRACT*SUBL
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
   subroutine MARSHPALMQ2(RAIN,PR,DIAM3,NTOTAL,W,VE)

      real, intent(in ) :: RAIN,PR     ! in kg m**-2 s**-1, mbar
      real, intent(out) :: DIAM3,NTOTAL,W,VE

      real :: RAIN_DAY,SLOPR,DIAM1

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
   !==========================================================

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

   function GET_ALHX3(T) RESULT(ALHX3)

      real, intent(in) :: T
      real :: ALHX3

      real :: T_X

      T_X = T_ICE_MAX

      if ( T < CLDPARAMS%T_ICE_ALL ) then
         ALHX3=MAPL_ALHS
      end if

      if ( T > T_X ) then
         ALHX3=MAPL_ALHL
      end if

      if ( (T <= T_X) .and. (T >= CLDPARAMS%T_ICE_ALL) ) then
         ALHX3 = MAPL_ALHS + (MAPL_ALHL-MAPL_ALHS)*( T - CLDPARAMS%T_ICE_ALL ) /( T_X - CLDPARAMS%T_ICE_ALL )
      end if

   end function GET_ALHX3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Partitions DQ into ice and liquid. Follows Barahona et al. GMD. 2014
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
         CNVFRC,SRFTYPE   , &  
         DQALL            , &
         FQI              , &
         needs_preexisting )

      real ,  intent(in   )    :: DTIME, PL, TE       !, RHCR
      real ,  intent(inout   )    ::  DQALL 
      real ,  intent(in)    :: QV, QLLS, QLCN, QICN, QILS
      real ,  intent(in)    :: CF, AF, NL, NI, CNVFRC,SRFTYPE
      real, intent (out) :: FQI
      logical, intent (in)  :: needs_preexisting
      
      real  :: DC, TEFF,DEP, &
            DQSL, DQSI, QI, TC, &
            DIFF, DENAIR, DENICE, AUX, &
            QTOT, LHCORR,  QL, DQI, DQL, &
            QVINC, QSLIQ, CFALL,  &
            QSICE, fQI_0, FQA, NIX

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

      elseif(TE .le. CLDPARAMS%T_ICE_ALL) then !ice cloud

         FQI   = 1.0

      else !mixed phase cloud

         FQI   = 0.0
         
          if (QILS .le. 0.0) then 
           
                    if (needs_preexisting) then
                   ! new 0518 this line ensures that only preexisting ice can grow by deposition.
                  ! Only works if explicit ice nucleation is available (2 moment muphysics and up)                        
                    else
                      fQi  =   ice_fraction( TE, CNVFRC,SRFTYPE )
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



   !=============================================================================
   ! Subroutine Pfreezing: calculates the probability of finding a supersaturated parcel in the grid cell
   !SC_ICE is the effective freezing point for ice (Barahona & Nenes. 2009)
   ! Modified 02/19/15. in situ nucleation only occurs in the non_convective part of the grid cell


   subroutine Pfreezing (     &             
         PDFFLAG  , &
         ALPHA    , &
         PL       , &
         TE       , &
         QV       , &              
         QCl      , &
         QAl      , & 
         QCi      , &
         QAi      , &    
         SC_ICE   , &
         CF       , &
         AF       , &          
         PF       )



      integer, intent(IN)  :: PDFFLAG
      real ,   intent(in)  :: PL,ALPHA, QV, SC_ICE, AF, TE, &
                              QCl, QCi, QAl, QAi, CF  
      real ,   intent(out) :: PF

      real   :: qt, QCx, QSn, tmpARR, CFALL, QVx, CFio, QA, QC, DQSx
      real   :: sigmaqt1, sigmaqt2, qsnx


      QA = QAl + QAi  
      QC=QCl+QCi 
      
      CFALL = AF

      if ( CFALL >= 1.0 ) then 
         PF = 0.0 
         return 
      end if

      QSn = QSATIC(         &
            TE   , &
            PL*100.0 , DQ=DQSx ) !only with respect to ice
      QSn  = MAX( QSn , 1.0e-9 )  

      tmpARR = 0.0    
      if  ( CFALL < 0.99 ) then 
         tmpARR = 1./(1.0-CFALL)
      end if

      QCx = QC*tmpARR
      QVx = ( QV - QSn*CFALL )*tmpARR

      qt =  QCx  +  QVx

      CFio=0.0

      if(pdfflag.lt.2) then
           sigmaqt1  = max(ALPHA, 0.01)*QSn
           sigmaqt2  = max(ALPHA, 0.01)*QSn 
      elseif(pdfflag.eq.2) then
         ! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
         ! for triangular, skewed r : sigmaqt1 < sigmaqt2
         sigmaqt1  = ALPHA*QSn
         sigmaqt2  = ALPHA*QSn
      elseif(pdfflag .eq. 4) then !lognormal (sigma is dimmensionless)
         sigmaqt1 =  max(ALPHA/sqrt(3.0), 0.001)     
      endif

      qsnx= Qsn*SC_ICE

     call pdffrac(pdfflag,qt,sigmaqt1,sigmaqt2,qsnx,CFio)

      PF  = CFio*(1.0-CFALL)        

      PF=min(max(PF, 0.0), 0.999)


   end subroutine Pfreezing




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Instantaneous freezing of condensate!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine meltfrz_inst  (     &
         IM,JM,LM , &
         TE       , &
         QCL       , &
         QAL       , &
         QCI        , &
         QAI        , &               
         NL       , &
         NI            )

      integer, intent(in)                             :: IM,JM,LM
      real ,   intent(inout), dimension(:,:,:)   :: TE,QCL,QCI, QAL, QAI, NI, NL

      real ,   dimension(im,jm,lm)              :: dQil, DQmax, QLTOT, QITOT, dNil, FQA

      QITOT= QCI+QAI
      QLTOT=QCL + QAL
      FQA = 0.0


      where (QITOT+QLTOT .gt. 0.0)
         FQA= (QAI+QAL)/(QITOT+QLTOT)
      end where


      dQil = 0.0
      dNil =0.0
      DQmax  = 0.0

      ! freeze liquid instantaneosly below -40 C
      where( TE <= CLDPARAMS%T_ICE_ALL )
         DQmax = (CLDPARAMS%T_ICE_ALL - TE)*MAPL_CP/(MAPL_ALHS-MAPL_ALHL)   
         dQil = min(QLTOT , DQmax)  
      end where

      where ((dQil .le. DQmax) .and. (dQil .gt. 0.0))
         dNil = NL
      end where

      where ((dQil .gt. DQmax) .and. (dQil .gt. 0.0)) 
         dNil  =  NL*DQmax/dQil 
      end where

      dQil = max(  0., dQil )
      QITOT = max(QITOT + dQil, 0.0)
      QLTOT= max(QLTOT -  dQil, 0.0)                  
      NL  = NL - dNil
      NI   = NI  + dNil
      TE   = TE + (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP

      dQil = 0.0
      dNil =0.0
      DQmax  = 0.0

      ! melt ice instantly above 0^C
      where( TE > T_ICE_MAX )
         DQmax =  (TE-T_ICE_MAX) *MAPL_CP/(MAPL_ALHS-MAPL_ALHL)   
         dQil = min(QITOT, DQmax)  
      endwhere

      where ((dQil .le. DQmax) .and. (dQil .gt. 0.0))
         dNil = NI
      end where
      where ((dQil .gt. DQmax) .and. (dQil .gt. 0.0)) 
         dNil  =  NI*DQmax/dQil 
      end where
      dQil = max(  0., dQil )
      QLTOT =  max(QLTOT+ dQil, 0.)
      QITOT = max(QITOT - dQil, 0.) 
      NL  = NL + dNil
      NI   = NI  - dNil

      TE   = TE - (MAPL_ALHS-MAPL_ALHL)*dQil/MAPL_CP

      QCI   = QITOT*(1.0-FQA)
      QAI = QITOT*FQA      
      QCL   = QLTOT*(1.0-FQA)
      QAL = QLTOT*FQA

   end subroutine meltfrz_inst



  
   !C=======================================================================
   !C
   !C *** REAL FUNCTION erf (overwrites previous versions)
   !C *** THIS SUBROUTINE CALCULATES THE ERROR FUNCTION USING A
   !C *** POLYNOMIAL APPROXIMATION
   !C
   !C=======================================================================
   !C
   REAL FUNCTION erf_app(x)
      REAL :: x
      REAL*8:: AA(4), axx, y
      DATA AA /0.278393d0,0.230389d0,0.000972d0,0.078108d0/

      y = dabs(dble(x))
      axx = 1.d0 + y*(AA(1)+y*(AA(2)+y*(AA(3)+y*AA(4))))
      axx = axx*axx
      axx = axx*axx
      axx = 1.d0 - (1.d0/axx)
      if(x.le.0.) then
         erf_app = sngl(-axx)
      else
         erf_app = sngl(axx)
      endif
      RETURN
   END FUNCTION erf_app

      
end module cldmacro
