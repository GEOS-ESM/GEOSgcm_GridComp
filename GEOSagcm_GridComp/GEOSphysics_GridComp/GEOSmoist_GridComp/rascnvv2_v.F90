! $Id$

      module module_ras

      use GEOS_Mod
      use GEOS_UtilsMod, only : DQSAT=>GEOS_DQsat

      implicit none
      SAVE  ! ????

      integer,        parameter :: kind_phys=8
      real(kind_phys), parameter :: grav = MAPL_GRAV
      real(kind_phys), parameter :: CP   = MAPL_CP
      real(kind_phys), parameter :: RGAS = MAPL_RGAS
      real(kind_phys), parameter :: ALHL = MAPL_ALHL
      real(kind_phys), parameter :: ALHF = MAPL_ALHF
      real(kind_phys), parameter :: RKAP = MAPL_KAPPA
      real(kind_phys), parameter :: NU   = .605

      integer,               parameter :: nrcmax=15
      real (kind=kind_phys), parameter :: delt_c=1800.0
      logical,               parameter :: fix_ncld_hr=.true.
!
      real(kind=kind_phys) ZERO, HALF, ONE, TWO
      real(kind=kind_phys) FOUR_P2,FOUR
      real(kind=kind_phys) ONE_M1,ONE_M2,ONE_M5,ONE_M6,ONE_M10
      PARAMETER (ZERO=0.0, HALF=0.5,  ONE=1.0, TWO=2.0)
      PARAMETER (FOUR_P2=4.E2,FOUR=4.,ONE_M10=1.E-10,ONE_M6=1.E-6       &
     &,          ONE_M5=1.E-5,ONE_M2=1.E-2,ONE_M1=1.E-1)
!
      real(kind=kind_phys), parameter :: cmb2pa = 100.0  ! Conversion from MB to PA
      real(kind=kind_phys) onebg, gravcon, gravfac, elocp, elfocp,      &
     &                     rkapi, rkpp1i,  zfac
!
      parameter (ONEBG   = ONE / GRAV,    GRAVCON = cmb2pa * ONEBG      &
     &,          GRAVFAC = GRAV / CMB2PA, ELOCP   = ALHL / CP           &
     &,          ELFOCP  = (ALHL+ALHF) / CP                             &
     &,          RKAPI   = ONE / RKAP,    RKPP1I  = ONE / (ONE+RKAP)    &
     &,          zfac    = 0.28888889E-4 * ONEBG)
!
      logical, parameter :: advcld=.true., advups=.true.
!
      real(kind=kind_phys), allocatable ::  RASAL(:)
      real(kind=kind_phys) RHMAX,  qudfac, QUAD_LAM, RHRAM, TESTMB,     &
     &                     TSTMBI, HCRIT,  DD_DP,   RKNOB, AFC, EKNOB

      PARAMETER (DD_DP=450.0,  RKNOB=1.5, EKNOB=1.0)
      PARAMETER (RHMAX=1.0   )  !  MAX RELATIVE HUMIDITY
      PARAMETER (QUAD_LAM=1.0)  !  MASK FOR QUADRATIC LAMBDA
      PARAMETER (RHRAM=0.05)    !  PBL RELATIVE HUMIDITY RAMP
      PARAMETER (HCRIT=4000.0)  !  Critical Moist Static Energy
      parameter (qudfac=quad_lam*half)
      parameter (testmb=0.1, tstmbi=one/testmb)
!
      real(kind=kind_phys) ALMIN1, ALMIN2, ALMAX
      real(kind=kind_phys) facdt
      PARAMETER (ALMIN1=0.00E-6, ALMIN2=2.50E-5, ALMAX=1.0E-2)
      PARAMETER (facdt=0.5)     !  relaxation time control?
      
      real(kind=kind_phys), parameter :: BLDMAX = 200.0
!
      INTEGER KBLMX
      real(kind=kind_phys) C0, C0I, QI0, QW0

      PARAMETER (QI0=1.0E-5, QW0=1.0E-5)
      PARAMETER (C0I=1.0E-3)
      parameter (c0=2.0e-3)

      real(kind=kind_phys) TF, TCR, TCRF, TCL
      parameter (TF=233.16, TCR=263.16, TCRF=1.0/(TCR-TF),TCL=2.0)

!     For Tilting Angle Specification

      real(kind=kind_phys) REFP(6), REFR(6), TLAC(8), PLAC(8), TLBPL(7) &
     &,                    drdp(5), VTP
!
      DATA PLAC/100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0/
      DATA TLAC/ 35.0,  25.0,  20.0,  17.5,  15.0,  12.5,  10.0,  7.5/
      DATA REFP/500.0, 300.0, 250.0, 200.0, 150.0, 100.0/
      DATA REFR/ 1.0,   2.0,  3.0,   4.0,   6.0,   8.0/
!
      real(kind=kind_phys) AC(16), AD(16)
!
      integer,parameter:: nqrp=500001
      real(kind=kind_phys) C1XQRP, C2XQRP, TBQRP(NQRP), TBQRA(NQRP),    &
     &                     TBQRB(NQRP)
!
      integer,parameter:: nvtp=10001
      real(kind=kind_phys) C1XVTP, C2XVTP, TBVTP(NVTP)
!
      contains

      subroutine set_ras_afc(dt)
      implicit none
      real(kind=kind_phys) DT
      AFC = -(1.01097E-4*DT)*(3600./DT)**0.57777778
      end subroutine set_ras_afc

      subroutine ras_init(levs,  me)
!
      Implicit none
!
      integer levs
!
      real(kind=kind_phys) actp,   facm, tem,  actop
      real(kind=kind_phys) rasalf, tem1, tem2
      integer              i, l, me

      PARAMETER (ACTP=1.7,   FACM=0.5)
!
      real(kind=kind_phys) PH(15),    A(15)
!
      DATA PH/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0    &
     &,       550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/
!
       DATA A/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677           &
     &,       0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664            &
     &,       0.0553, 0.0445, 0.0633/
!
      logical first
      data first/.true./
!
      if (first) then
!
        allocate (rasal(levs))
!                                   set critical workfunction arrays
        ACTOP = ACTP*FACM
        DO L=1,15
          A(L) = A(L)*FACM
        ENDDO
        DO L=2,15
          TEM   = 1.0 / (PH(L) - PH(L-1))
          AC(L) = (PH(L)*A(L-1) - PH(L-1)*A(L)) * TEM
          AD(L) = (A(L) - A(L-1)) * TEM
        ENDDO
        AC(1)  = ACTOP
        AC(16) = A(15)
        AD(1)  = 0.0
        AD(16) = 0.0

        CALL SETQRP
        CALL SETVTP

        RASALF  = 0.30

        DO L=1,LEVS
          RASAL(L) = RASALF
        ENDDO
!
        do i=1,7
          tlbpl(i) = (tlac(i)-tlac(i+1)) / (plac(i)-plac(i+1))
        enddo
        do i=1,5
          drdp(i)  = (REFR(i+1)-REFR(i)) / (REFP(i+1)-REFP(i))
        enddo
!
        VTP    = 36.34*SQRT(1.2)* (0.001)**0.1364
!
        if (me .eq. 0) print *,' NO DOWNDRAFT FOR CLOUD TYPES'          &
     &,        ' DETRAINING WITHIN THE BOTTOM ',DD_DP,' hPa LAYERS'
!
        first = .false.
      endif
!
      end subroutine ras_init

      end module module_ras


      SUBROUTINE CRTWRK(PL, ACR)
      use module_ras
      Implicit none
!
      real(kind=kind_phys) PL, ACR
      INTEGER IWK
!
      IWK = PL * 0.02 - 0.999999999
      IWK = MAX(1, MIN(IWK,16))
      ACR = AC(IWK) + PL * AD(IWK)
!
      RETURN
      END
      SUBROUTINE CLOUD(                                                 &
     &                  K, KD, M                                        &
     &,                 RASALF, FRACBL, MAX_NEG_BOUY                    &
     &,                 ALFINT, ALFINQ, RHFACL, RHFACS, garea           &
     &,                 alfind, rhc_ls                                  &

     &,                 TOI, QOI, ROI, PRS, PRSM, phil, phih            &
     &,                 QLI, QII, KPBL, DSFC                            &
     &,                 CD,lprnt, trcfac                                &

     &,                 TCU, QCU, RCU, PCU, FLX                         &
     &,                 CUP, REVAP, DT                                  &
     &,                 WFNC, WRKFUN, CALKBL, CRTFUN, TLA, DNDRFT, DPD)  

!
!***********************************************************************
!******************** Relaxed  Arakawa-Schubert ************************
!****************** Plug Compatible Scalar Version *********************
!************************ SUBROUTINE CLOUD  ****************************
!************************  October 2004     ****************************
!********************  VERSION 2.0  (modified) *************************
!************* Shrinivas.Moorthi@noaa.gov (301) 763 8000(X7233) ********
!***********************************************************************
!*Reference:
!-----------
!     NOAA Technical Report NWS/NCEP 99-01:
!     Documentation of Version 2 of Relaxed-Arakawa-Schubert
!     Cumulus Parameterization with Convective Downdrafts, June 1999.
!     by S. Moorthi and M. J. Suarez.
!
!***********************************************************************
!
!===>    UPDATES CLOUD TENDENCIES DUE TO A SINGLE CLOUD
!===>    DETRAINING AT LEVEL KD.
!
!***********************************************************************

!===>  K      INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
!===>  KD     INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
!===>  M      INPUT   NUMBER OF TRACERS. MAY BE ZERO.
!===>  RASALF INPUT   
!===>  FRACBL INPUT   MASS FLUX LIMIT
!===>  MAX_NEG_BOUY  INPUT  FRACTION OF NEG ALLOWED
!===>  ALFINT(K)  INPUT   INTERFACE UPWIND (=1) PARAMETER FOR T
!===>  ALFINQ(K)  INPUT   INTERFACE UPWIND (=1) PARAMETER FOR Q
!===>  RHFACL INPUT   CRITICAL RH AT PBL TOP (LAND)
!===>  RHFACS INPUT   CRITICAL RH AT PBL TOP (SEA) ???
!===>  GAREA  INPUT   AREA OF GRID BOX
!===>  ALFIND(K)  INPUT   INTERFACE UPWIND (=1) PARAMETER FOR DOWNDRAFT
!===>  RHC_LS(K)  INPUT   CRITICAL RH FOR RE-EVAP

!===>  TOI(K)     INOUT   TEMPERATURE             KELVIN
!===>  QOI(K)     INOUT   SPECIFIC HUMIDITY       NON-DIMENSIONAL
!===>  ROI(K,M)   INOUT   TRACER                  ARBITRARY
!===>  PRS(K+1)   INPUT   PRESSURE @ EDGES        MB
!===>  PRSM(K)    INPUT   PRESSURE @ LAYERS       MB
!===>  PHIL(K)    INPUT   GEOPOTENTIAL @ LAYERS IN MKS units
!===>  PHIH(K+1)  INPUT   GEOPOTENTIAL @ EDGES  IN MKS units
!===>  QLI(K)     INOUT   LIQUID WATER            NON-DIMENSIONAL
!===>  QII(K)     INOUT   ICE                     NON-DIMENSIONAL

!===>  KPBL       INOUT   PBL TOP
!===>  DSFC       INOUT   GUSTINESS (m/s)

!===>  CD         INOUT   DRAG COEF. (NON-DIMENSIONAL)
!===>  LPRNT      INOUT
!===>  TRCFAC     INOUT   tracer pg factor for u AND V 


!===>  TCU(K  )   UPDATE  TEMPERATURE TENDENCY       DEG
!===>  QCU(K  )   UPDATE  WATER VAPOR TENDENCY       (G/G)
!===>  RCU(K,M)   UPDATE  TRACER TENDENCIES          ND
!===>  PCU(K-1)   UPDATE  PRECIP @ BASE OF LAYER     KG/M^2
!===>  FLX(K  )   UPDATE  UPDRAFT MASS FLUX @ TOP OF LAYER   KG/M^2
!===>  CUP        UPDATE  PRECIPITATION AT THE SURFACE KG/M^2
!===>  REVAP      UPDATE  LOGICAL TO CONTROL REEVAP IN ENVIRONMENT
!===>  DT         UPDATE  TIME STEP USED FOR GUSTINESS
!===>  WFNC       UPDATE  WORK FUNCTION (in if CRTFUN = F, out if WRKFUN = T ) 
!===>  WRKFUN     UPDATE  LOGICAL TO DO WRKFUNC ONLY
!===>  CALKBP     UPDATE  LOGICAL TO COMPUTE KPBL INTERNALLY
!===>  CRTFUN     UPDATE  LOGICAL TO CONTROL CRITICAL WORK FUNCTION
!===>  TLA        UPDATE  TILTING ANGLE FOR DD (NEG USE TABLE)
!===>  DNDRFT     INPUT   LOGICAL TO CONTROL DOWNDRAFT
!===>  DPD        INPUT   Minumum Cloud Depth for DOWNDRFAT Computation hPa
!
      use module_ras
      IMPLICIT NONE
!
!  INPUT ARGUMENTS

      LOGICAL REVAP, DNDRFT, WRKFUN, CALKBL, CRTFUN, CALCUP
      logical lprnt
      INTEGER K, KD, M


      real(kind=kind_phys) TOI(K),    QOI(K ),  PRS(K+1), PRSM(K)       &
     &,                    QLI(K),    QII(K)                            &
     &,                    PHIH(K+1), ROI(K,M), PHIL(K)           



      real(kind=kind_phys) CD,        UFN,     DSFC
      INTEGER KPBL,   KBL,     KB1

      real(kind=kind_phys) RASALF, FRACBL, MAX_NEG_BOUY, ALFINT(K),     &
     &                     RHFACL, RHFACS, garea
      real(kind=kind_phys) ALFINQ(K), DPD, alfind(k), rhc_ls(k)
      real(kind=kind_phys) trcfac(M)
 
!  UPDATE ARGUMENTS

      real(kind=kind_phys) TCU(K),  QCU(K), RCU(K,M)                    &
     &,                    TCD(K),  QCD(K), PCU(K), FLX(K), CUP

!  TEMPORARY WORK SPACE

      real(kind=kind_phys) HOL(KD:K),  QOL(KD:K),   GAF(KD:K+1)         &
     &,                    HST(KD:K),  QST(KD:K),   TOL(KD:K)           &
     &,                    GMH(KD:K),  GMS(KD:K+1), GAM(KD:K+1)         &
     &,                    AKT(KD:K),  AKC(KD:K),   BKC(KD:K)           &
     &,                    LTL(KD:K),  RNN(KD:K),   FCO(KD:K)           &
     &,                                PRI(KD:K)                        &
     &,                    QIL(KD:K),  QLL(KD:K)                        &
     &,                    ZET(KD:K),  XI(KD:K),    RNS(KD:K)           &
     &,                    Q0U(KD:K),  Q0D(KD:K),   vtf(KD:K)           &
     &,                    DLB(KD:K+1),DLT(KD:K+1), ETA(KD:K+1)         &
     &,                    PRL(KD:K+1)                                  &
     &,                    CIL(KD:K),  CLL(KD:K),   ETAI(KD:K)

      real(kind=kind_phys) ALM,   DET,    HCC,  CLP                     &
     &,                    HSU,   HSD,    QTL,  QTV                     &
     &,                    AKM,   WFN,    HOS,  QOS                     &
     &,                    AMB,   TX1,    TX2,  TX3                     &
     &,                    TX4,   TX5,    QIS,  QLS                     &
     &,                    HBL,   QBL,    RBL(M)                        &
     &,                    QLB,   QIB,    PRIS                          &
     &,                    WFNC,  TX6,    ACR                           &
     &,                    TX7,   TX8,    TX9,  RHC                     &
     &,                    hstkd, qstkd,  ltlkd, q0ukd, q0dkd, dlbkd    &
     &,                    qtp, qw00, qi00, qrbkd                       &
     &,                    hstold, rel_fac


      LOGICAL UNSAT, ep_wfn

      LOGICAL LOWEST, SKPDD

      real(kind=kind_phys) TL, PL, QL, QS, DQS, ST1, SGN, TAU,          &
     &                     QTVP, HB, QB, TB, QQQ,                       &
     &                     HCCP, DS, DH, AMBMAX, X00, EPP, QTLP,        &
     &                     DPI, DPHIB, DPHIT, DEL_ETA, DETP,            &
     &                     TEM, TEM1, TEM2, TEM3, TEM4,                 &
     &                     ST2, ST3,                                    &
     &                     ERRMIN, ERRMI2, ERRH, ERRW, ERRE, TEM5,      &
     &                     TEM6, HBD, QBD, st1s
      parameter (ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
!     parameter (c0=1.0e-3, KBLMX=20, ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
      INTEGER I, L,  N,  KD1, II                                        &
     &,       KP1, IT, KM1, KTEM, KK, KK1, LM1, LL, LP1, kbls, kmxh

      real avt, avq, avr, avh
!
!     REEVAPORATION
      real(kind=kind_phys), parameter :: rainmin=1.0e-8
      real(kind=kind_phys), parameter :: oneopt9=1.0/0.09
      real(kind=kind_phys), parameter :: oneopt4=1.0/0.04

      real(kind=kind_phys) CLFRAC, DT, clf, clvfr

      real(kind=kind_phys) ACTEVAP,AREARAT,DELTAQ,MASS,MASSINV,POTEVAP  &
     &,                    TEQ,QSTEQ,DQDT,QEQ                        
!    &,                    ELOCP,GRAVCON, GRAVFAC,  AFC, RKNOB, ELFOCP
!
!     Temporary workspace and parameters needed for downdraft
!
      real(kind=kind_phys) TLA, GMF
!
      real(kind=kind_phys) BUY(KD:K+1), QRB(KD:K),   QRT(KD:K)          &
     &,                    ETD(KD:K+1), HOD(KD:K+1), QOD(KD:K+1)        &
     &,                    GHD(KD:K),   GSD(KD:K),   EVP(KD:K)          &
     &,                    ETZ(KD:K),   CLDFR(KD:K)                     &
     &,                    TRAIN, DOF, CLDFRD                           &
     &,                    FAC, RSUM1, RSUM2, RSUM3, dpneg
      INTEGER IDH
      LOGICAL DDFT, UPDRET
!
!***********************************************************************
!
!!CFPP$ EXPAND (QSATCN, CRTWRK)
!!CFPP$ NOCONCUR R
!
      do l=1,K
        tcd(L) = 0.0
        qcd(L) = 0.0
      enddo
!
      KP1     = K  + 1
      KM1     = K  - 1
      KD1     = KD + 1
      kblmx   = k / 2
!
!
      CLDFRD   = 0.0
      DOF      = 0.0
      PRL(KP1) = PRS(KP1)
!
      DO L=KD,K
        RNN(L) = 0.0
        ZET(L) = 0.0
        XI(L)  = 0.0
!
        TOL(L) = TOI(L)
        QOL(L) = QOI(L)
        PRL(L) = PRS(L)
        BUY(L) = 0.0
        CLL(L) = QLI(L)
        CIL(L) = QII(L)
      ENDDO
!
      DO L=KD, K
        DPI    = ONE / (PRL(L+1) - PRL(L))
        PRI(L) = GRAVFAC * DPI
!
        PL     = PRSM(L)
        TL     = TOL(L)


        AKT(L) = (PRL(L+1) - PL) * DPI

        CALL QSATCN(TL, PL, QS, DQS,lprnt)

        QST(L) = QS
        GAM(L) = DQS * ELOCP
        ST1    = ONE + GAM(L)
        GAF(L) = (ONE/ALHL) * (GAM(L)/(ONE + GAM(L)))
 
        QL     = MAX(MIN(QS*RHMAX,QOL(L)), ONE_M10)
        QOL(L) = QL
 
        TEM    = CP * TL
        LTL(L) = TEM * ST1 / (ONE+NU*(QST(L)+TL*DQS))
        vtf(L) = 1.0 + NU * QL
        ETA(L) = ONE / (LTL(L) * VTF(L))

        HOL(L) = TEM + QL * ALHL
        HST(L) = TEM + QS * ALHL
      ENDDO
!
      ETA(K+1) = ZERO
      GMS(K)   = ZERO
!
      AKT(KD)  = HALF
      GMS(KD)  = ZERO
!
      CLP      = ZERO
!
      GAM(K+1) = GAM(K)
      GAF(K+1) = GAF(K)
!
      DO L=K,KD1,-1
 
        DPHIB  = PHIL(L) - PHIH(L+1)
        DPHIT  = PHIH(L) - PHIL(L)
!
        DLB(L) = DPHIB * ETA(L)
        DLT(L) = DPHIT * ETA(L)
!
        QRB(L) = DPHIB
        QRT(L) = DPHIT
!
        ETA(L) = ETA(L+1) + DPHIB

!
        HOL(L) = HOL(L) + ETA(L)
        hstold = hst(l)
        HST(L) = HST(L) + ETA(L)
!

        ETA(L) = ETA(L) + DPHIT
      ENDDO
!
!     For the cloud top layer
!
      L = KD

      DPHIB  = PHIL(L) - PHIH(L+1)
!
      DLB(L) = DPHIB * ETA(L)
!
      QRB(L) = DPHIB
      QRT(L) = DPHIB
!
      ETA(L) = ETA(L+1) + DPHIB

      HOL(L) = HOL(L) + ETA(L)
      HST(L) = HST(L) + ETA(L)
!
!     To determine KBL internally -- If KBL is defined externally
!     the following two loop should be skipped

      IF (CALKBL) THEN
         KTEM = MAX(KD, K-KBLMX-2)
         kmxh = k


         DO L=kmxh,KTEM+1,-1
           kbls = l
           if (hst(l-1) .gt. hst(l)) exit
         ENDDO
         KBL   = Kmxh
         TX1   = ZERO
         UNSAT = .FALSE.
         DO L=kmxh-1,KTEM,-1
           TEM = HOL(K) - HOL(L)
           TX3 = (HOL(L) - HOL(L+1)) / (PRL(L+2) - PRL(L))

           IF (TX3 .LT. TX1 .AND. TEM .LT. HCRIT) THEN
             TX1   = TX3
             KBL   = L
             UNSAT = .TRUE.
           ELSEIF (UNSAT .AND.                                          &
     &           ( ((KBL .LT. K-1) .AND. TX3 .GT. 0.5*TX1)              &
     &              .OR. TEM .GT. HCRIT) ) THEN
             TX1 = -1.0E20
           ENDIF
         ENDDO

         ii = kbl
         do l=ktem,kmxh-1
           if (hol(kmxh) .gt. hst(l)) kbl = l+1    ! Commented on 09/20/04
         enddo

         if (prl(K+1) - prl(ii) .gt. 50.0 .and. ii .gt. kbl) kbl = ii
         if (kbl .ne. ii) then
           if (PRL(K+1)-PRL(KBL) .gt. bldmax) kbl = max(kbl,ii)
         endif
         KBL  = min(k, MAX(KBL,K-KBLMX))
         tem1 = min( (prl(kbl) - prl(kd))*0.05, 20.0d0   )
         tem1 = max( prl(k+1)-prl(k), tem1 )

         if (prl(k+1)-prl(kbl) .lt. tem1) then
           KTEM = MAX(KD+1, K-KBLMX)
           do l=k,KTEM,-1
             tem = prl(k+1) - prl(l)
             if (tem .gt. tem1) then
               kbl = min(kbl,l)
               exit
             endif
           enddo
         endif
!!!

         KPBL = KBL
      ELSE
         KBL  = KPBL
      ENDIF
      KBL      = MAX(KBL,KD)
      KB1      = KBL - 1
!!
      if(kb1 .le. kd)then
          !!!write(*,*) ' return con kb1=',kb1,' kd=',kd,' EXIT CLOUD'
        return
      endif
      if(PRL(K+1)-PRL(KBL) .gt. bldmax)then
         !!!write(*,*) " return con bldmax  exit"
        return
      endif
!
!
      PRIS     = ONE / (PRL(K+1)-PRL(KBL))
      TX1      = ETA(KBL)
!
      GMS(KBL) = 0.0
      XI(KBL)  = 0.0
      ZET(KBL) = 0.0
!
      DO L=K,KD,-1
        IF (L .GE. KBL) THEN
          ETA(L) = (PRL(K+1)-PRL(L)) * PRIS
        ELSE
          ZET(L) = (ETA(L) - TX1) * ONEBG
          XI(L)  =  ZET(L) * ZET(L) * QUDFAC
          ETA(L) =  ZET(L) - ZET(L+1)
          GMS(L) =  XI(L)  - XI(L+1)
        ENDIF
      ENDDO
!
      HBL = HOL(K) * ETA(K)
      QBL = QOL(K) * ETA(K)
      QLB = CLL(K) * ETA(K)
      QIB = CIL(K) * ETA(K)
      TX1 = QST(K) * ETA(K)
!
      DO L=KM1,KBL,-1
         TEM = ETA(L) - ETA(L+1)
         HBL = HBL + HOL(L) * TEM
         QBL = QBL + QOL(L) * TEM
         QLB = QLB + CLL(L) * TEM
         QIB = QIB + CIL(L) * TEM
         TX1 = TX1 + QST(L) * TEM
      ENDDO
!                                   Find Min value of HOL in TX2
      TX2 = HOL(KD)
      IDH = KD1
      DO L=KD1,KB1
        IF (HOL(L) .LT. TX2) THEN
           TX2 = HOL(L)
           IDH = L             ! Level of minimum moist static energy!
        ENDIF
      ENDDO
      IDH = 1
      IDH = MAX(KD1, IDH)
!
      TEM1 = HBL - HOL(KD)
      TEM  = HBL - HST(KD1)                                             &
     &             - LTL(KD1) *( NU *(QOL(KD1)-QST(KD1)))
      LOWEST = KD .EQ. KB1
!

      TX1   = RHFACS - QBL / TX1       !     Average RH
      UNSAT = (TEM .GT. ZERO .OR. (LOWEST .AND. TEM1 .GE. ZERO))        &
     &         .AND. (TX1 .LT. RHRAM)                                   &
     &         .AND. (KBL .GT. KD)


!===>  IF NO SOUNDING MEETS FIRST CONDITION, RETURN

      IF (.NOT. UNSAT) then 
            !!!write(*,*) " return con 1"
         RETURN
        endif
!     
      RHC    = MAX(ZERO, MIN(ONE, EXP(-20.0*TX1) ))
!
      DO N=1,M
        RBL(N) = ROI(K,N) * ETA(K)
      ENDDO
      DO N=1,M
        DO L=KM1,KBL,-1
          RBL(N) = RBL(N) + ROI(L,N)*(ETA(L)-ETA(L+1))
        ENDDO
      ENDDO
!
      TX4    = 0.0
      TX5    = 0.0
!
      TX3      = QST(KBL) - GAF(KBL) * HST(KBL)
      QIL(KBL) = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KBL))*TCRF))
!
      DO L=KB1,KD1,-1
        TEM      = QST(L) - GAF(L) * HST(L)
        TEM1     = (TX3 + TEM) * 0.5
        ST2      = (GAF(L)+GAF(L+1)) * 0.5
!
        FCO(L+1) =            TEM1 + ST2 * HBL


        RNN(L+1) = ZET(L+1) * TEM1 + ST2 * TX4
        GMH(L+1) = XI(L+1)  * TEM1 + ST2 * TX5
!
        TX3      = TEM
        TX4      = TX4 + ETA(L) * HOL(L)
        TX5      = TX5 + GMS(L) * HOL(L)
!
        QIL(L)   = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(L))*TCRF))
        QLL(L+1) = (0.5*ALHF) * ST2 * (QIL(L)+QIL(L+1)) + ONE
      ENDDO
!
!     FOR THE CLOUD TOP -- L=KD
!
      L = KD
!
      TEM      = QST(L) - GAF(L) * HST(L)
      TEM1     = (TX3 + TEM) * 0.5
      ST2      = (GAF(L)+GAF(L+1)) * 0.5
!
      FCO(L+1) =            TEM1 + ST2 * HBL
      RNN(L+1) = ZET(L+1) * TEM1 + ST2 * TX4
      GMH(L+1) = XI(L+1)  * TEM1 + ST2 * TX5
!
      FCO(L)   = TEM + GAF(L) * HBL
      RNN(L)   = TEM * ZET(L) + (TX4 + ETA(L)*HOL(L)) * GAF(L)
      GMH(L)   = TEM * XI(L)  + (TX5 + GMS(L)*HOL(L)) * GAF(L)
!
!   Replace FCO for the Bottom
!
      FCO(KBL) = QBL
      RNN(KBL) = 0.0
      GMH(KBL) = 0.0
!
      QIL(KD)  =  MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KD))*TCRF))
      QLL(KD1) = (0.5*ALHF) * ST2 * (QIL(KD) + QIL(KD1)) + ONE
      QLL(KD ) = ALHF * GAF(KD) * QIL(KD) + ONE
!
!
      st1  = qil(kd)
      st2  = c0i * st1
      tem  = c0  * (1.0-st1)
      tem2 = st2*qi0 + tem*qw0
!
      DO L=KD,KB1
         tx2    = akt(l) * eta(l)
         tx1    = tx2 * tem2
         q0u(l) = tx1
         FCO(L) = FCO(L+1) - FCO(L) + tx1
         RNN(L) = RNN(L+1) - RNN(L)                                     &
     &          + ETA(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*zet(l)
         GMH(L) = GMH(L+1) - GMH(L)                                     &
     &          + GMS(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*xi(l)
!
         tem1   = (1.0-akt(l)) * eta(l)


         AKT(L) = QLL(L)   + (st2 + tem) * tx2


         AKC(L) = 1.0 / AKT(L)
!
         st1    = 0.5 * (qil(l)+qil(l+1))
         st2    = c0i * st1
         tem    = c0  * (1.0-st1)
         tem2   = st2*qi0 + tem*qw0
!
         BKC(L) = QLL(L+1) - (st2 + tem) * tem1
!
         tx1    = tem1*tem2
         q0d(l) = tx1
         FCO(L) = FCO(L) + tx1
         RNN(L) = RNN(L) + tx1*zet(l+1)
         GMH(L) = GMH(L) + tx1*xi(l+1)
      ENDDO


      qw00 = qw0
      qi00 = qi0
      ii = 0
  777 continue
!
      ep_wfn = .false.
      RNN(KBL) = 0.0
      TX3      = bkc(kb1) * (QIB + QLB)
      TX4      = 0.0
      TX5      = 0.0
      DO L=KB1,KD1,-1
        TEM    = BKC(L-1)       * AKC(L)
        TX3    = (TX3 + FCO(L)) * TEM
        TX4    = (TX4 + RNN(L)) * TEM
        TX5    = (TX5 + GMH(L)) * TEM
      ENDDO
      IF (KD .LT. KB1) THEN
         HSD   = HST(KD1)                                               &
     &         + LTL(KD1) *  NU *(QOL(KD1)-QST(KD1))
      ELSE
         HSD   = HBL
      ENDIF
!
      TX3 = (TX3 + FCO(KD)) * AKC(KD)
      TX4 = (TX4 + RNN(KD)) * AKC(KD)
      TX5 = (TX5 + GMH(KD)) * AKC(KD)
      ALM = ALHF*QIL(KD) - LTL(KD) * VTF(KD)
!
      HSU = HST(KD) + LTL(KD) * NU * (QOL(KD)-QST(KD))

!
!===> VERTICAL INTEGRALS NEEDED TO COMPUTE THE ENTRAINMENT PARAMETER
!
      TX1 = ALM * TX4
      TX2 = ALM * TX5

      DO L=KD,KB1
        TAU = HOL(L) - HSU
        TX1 = TX1 + TAU * ETA(L)
        TX2 = TX2 + TAU * GMS(L)
      ENDDO
!
!     MODIFY HSU TO INCLUDE CLOUD LIQUID WATER AND ICE TERMS
!

      HSU   = HSU - ALM * TX3
!
      CLP   = ZERO
      ALM   = -100.0
      HOS   = HOL(KD)
      QOS   = QOL(KD)
      QIS   = CIL(KD)
      QLS   = CLL(KD)
      UNSAT = HBL .GT. HSU .and. abs(tx1) .gt. 1.0e-4



!***********************************************************************


       ST1  = HALF*(HSU + HSD)
       IF (UNSAT) THEN
!
!  STANDARD CASE:
!   CLOUD CAN BE NEUTRALLY BOUYANT AT MIDDLE OF LEVEL KD W/ +VE LAMBDA.
!   EPP < .25 IS REQUIRED TO HAVE REAL ROOTS.
!
       clp = 1.0
       st2 = hbl - hsu

       if (tx2 .eq. 0.0) then
         alm = - st2 / tx1
         if (alm .gt. almax) alm = -100.0
       else
         x00 = tx2 + tx2
         epp = tx1 * tx1 - (x00+x00)*st2
         if (epp .gt. 0.0) then
           x00  = 1.0 / x00
           tem  = sqrt(epp)
           tem1 = (-tx1-tem)*x00
           tem2 = (-tx1+tem)*x00
           if (tem1 .gt. almax) tem1 = -100.0
           if (tem2 .gt. almax) tem2 = -100.0
           alm  = max(tem1,tem2)


         endif
       endif

!
!  CLIP CASE:
!   NON-ENTRAINIG CLOUD DETRAINS IN LOWER HALF OF TOP LAYER.
!   NO CLOUDS ARE ALLOWED TO DETRAIN BELOW THE TOP LAYER.
!
       ELSEIF ( (HBL .LE. HSU) .AND.                                    &
     &          (HBL .GT. ST1   )     ) THEN
         ALM = ZERO
         CLP = (HBL-ST1) / (HSU-ST1)
       ENDIF
!
      UNSAT = .TRUE.
      IF (ALMIN1 .GT. 0.0) THEN
        IF (ALM .GE. ALMIN1) UNSAT = .FALSE.
      ELSE
        LOWEST   = KD .EQ. KB1
        IF ( (ALM .GT. ZERO) .OR.                                       &
     &      (.NOT. LOWEST .AND. ALM .EQ. ZERO) ) UNSAT = .FALSE.
      ENDIF
!
!
!===>  IF NO SOUNDING MEETS SECOND CONDITION, RETURN
!
      IF (UNSAT) THEN
         IF (ii .gt. 0 .or. (qw00 .eq. 0.0 .and. qi00 .eq. 0.0)) then 
                !!!write(*,*) " return con 2 "
           RETURN
          endif
         CLP = 1.0
         ep_wfn = .true.
         GO TO 888
      ENDIF
!

      IF(CLP.GT.ZERO .AND. CLP.LT.ONE) THEN
        ST1     = HALF*(ONE+CLP)
        ST2     = ONE - ST1
        st1s    = st1
        hstkd   = hst(kd)
        qstkd   = qst(kd)
        ltlkd   = ltl(kd)
        q0ukd   = q0u(kd)
        q0dkd   = q0d(kd)
        dlbkd   = dlb(kd)
        qrbkd   = qrb(kd)
!
        HST(KD) = HST(KD)*ST1 + HST(KD1)*ST2
        HOS     = HOL(KD)*ST1 + HOL(KD1)*ST2
        QST(KD) = QST(KD)*ST1 + QST(KD1)*ST2
        QOS     = QOL(KD)*ST1 + QOL(KD1)*ST2
        QLS     = CLL(KD)*ST1 + CLL(KD1)*ST2
        QIS     = CIL(KD)*ST1 + CIL(KD1)*ST2
        LTL(KD) = LTL(KD)*ST1 + LTL(KD1)*ST2
!
        DLB(KD) = DLB(KD)*CLP
        qrb(KD) = qrb(KD)*CLP
        ETA(KD) = ETA(KD)*CLP
        GMS(KD) = GMS(KD)*CLP
        Q0U(KD) = Q0U(KD)*CLP
        Q0D(KD) = Q0D(KD)*CLP
      ENDIF
!
!
!***********************************************************************
!
!    Critical workfunction is included in this version
!
      ACR = 0.0
      TEM = PRL(KD1) - (PRL(KD1)-PRL(KD)) * CLP * HALF
      tx1 = PRL(KBL) - TEM

!
      tx2 = min(800.0d0 ,max(tx1,100.0d0 ))
      tem1    = log(tx2*0.01) / log(8.0)
      rel_fac = (dt * facdt)  / (3600.0 * 1.5)
!
      rel_fac = min(one,rel_fac)
      
      IF (CRTFUN) THEN
        CALL CRTWRK(TEM, ST1)
!       ACR = (PRL(K) - TEM) * ST1
        ACR = TX1 * ST1
      ENDIF
!
!===>  NORMALIZED MASSFLUX
!
!  ETA IS THE THICKNESS COMING IN AND THE MASS FLUX GOING OUT.
!  GMS IS THE THICKNESS OF THE SQUARE; IT IS LATER REUSED FOR GAMMA_S
!
!     ETA(K) = ONE

      DO L=KB1,KD,-1
        ETA(L)  = ETA(L+1) + ALM * (ETA(L) + ALM * GMS(L))
      ENDDO
      DO L=KD,KBL
        ETAI(L) = 1.0 / ETA(L)
      ENDDO

!
!===>  CLOUD WORKFUNCTION
!
      WFN   = ZERO
      AKM   = ZERO
      DET   = ZERO
      HCC   = HBL
      UNSAT = .FALSE.
      QTL   = QST(KB1) - GAF(KB1)*HST(KB1)
      TX1   = HBL
!
      qtv   = qbl
      det   = qlb + qib
!
      tx2   = 0.0
      dpneg = 0.0
!
      DO L=KB1,KD1,-1
         DEL_ETA = ETA(L) - ETA(L+1)
         HCCP = HCC + DEL_ETA*HOL(L)
!
         QTLP = QST(L-1) - GAF(L-1)*HST(L-1)
         QTVP = 0.5 * ((QTLP+QTL)*ETA(L)                                &
     &              + (GAF(L)+GAF(L-1))*HCCP)
         ST1  = ETA(L)*Q0U(L) + ETA(L+1)*Q0D(L)
         DETP = (BKC(L)*DET - (QTVP-QTV)                                &
     &        + DEL_ETA*(QOL(L)+CLL(L)+CIL(L)) + ST1)  * AKC(L)

!
         TEM1   = AKT(L)   - QLL(L)
         TEM2   = QLL(L+1) - BKC(L)
         RNS(L) = TEM1*DETP  + TEM2*DET - ST1

         qtp    = 0.5 * (qil(L)+qil(L-1))
         tem2   = min(qtp*(detp-eta(l)*qw00),                           &
     &               (1.0-qtp)*(detp-eta(l)*qi00))
         st1    = min(tx2,tem2)
         tx2    = tem2
!
         IF (rns(l) .lt. zero .or. st1 .lt. zero) ep_wfn = .TRUE.
         IF (DETP .LE. ZERO) UNSAT = .TRUE.

         ST1  = HST(L) - LTL(L)*NU*(QST(L)-QOL(L))


         TEM2 = HCCP   + DETP   * QTP * ALHF
!

         ST2  = LTL(L) * VTF(L)
         TEM5 = CLL(L) + CIL(L)
         TEM3 = (TX1  - ETA(L+1)*ST1 - ST2*(DET-TEM5*eta(l+1))) * DLB(L)
         TEM4 = (TEM2 - ETA(L  )*ST1 - ST2*(DETP-TEM5*eta(l)))  * DLT(L)
!

         ST1  = TEM3 + TEM4


         WFN = WFN + ST1       
         AKM = AKM - min(ST1,ZERO)


         if (st1 .lt. zero .and. wfn .lt. zero) then
           dpneg = dpneg + prl(l+1) - prl(l)
         endif

         BUY(L) = 0.5 * (tem3/(eta(l+1)*qrb(l)) + tem4/(eta(l)*qrt(l)))
!
         HCC = HCCP
         DET = DETP
         QTL = QTLP
         QTV = QTVP
         TX1 = TEM2

      ENDDO

      DEL_ETA = ETA(KD) - ETA(KD1)
      HCCP    = HCC + DEL_ETA*HOS
!
      QTLP = QST(KD) - GAF(KD)*HST(KD)
      QTVP = QTLP*ETA(KD) + GAF(KD)*HCCP
      ST1  = ETA(KD)*Q0U(KD) + ETA(KD1)*Q0D(KD)
      DETP = (BKC(KD)*DET - (QTVP-QTV)                                  &
     &     + DEL_ETA*(QOS+QLS+QIS) + ST1) * AKC(KD)
!
      TEM1    = AKT(KD)  - QLL(KD)
      TEM2    = QLL(KD1) - BKC(KD)
      RNS(KD) = TEM1*DETP  + TEM2*DET - ST1
!
      IF (rns(kd) .lt. zero) ep_wfn = .TRUE.
      IF (DETP.LE.ZERO) UNSAT = .TRUE.
!
  888 continue


      if (ep_wfn) then
        IF ((qw00 .eq. 0.0 .and. qi00 .eq. 0.0)) then 
             !!!write(*,*) "return con 3"
           RETURN
         endif
        if (ii .eq. 0) then
          ii  = 1
          if (clp .gt. 0.0 .and. clp .lt. 1.0) then
            hst(kd) = hstkd
            qst(kd) = qstkd
            ltl(kd) = ltlkd
            q0u(kd) = q0ukd
            q0d(kd) = q0dkd
            dlb(kd) = dlbkd
            qrb(kd) = qrbkd
          endif
          do l=kd,kb1
            FCO(L) = FCO(L) - q0u(l) - q0d(l)
            RNN(L) = RNN(L) - q0u(l)*zet(l) - q0d(l)*zet(l+1)
            GMH(L) = GMH(L) - q0u(l)*xi(l)  - q0d(l)*zet(l+1)
            ETA(L) = ZET(L) - ZET(L+1)
            GMS(L) = XI(L)  - XI(L+1)
            Q0U(L) = 0.0
            Q0D(L) = 0.0
          ENDDO
          qw00 = 0.0
          qi00 = 0.0


          go to 777
        else
          unsat = .true.
        endif
      endif
!
!
!
      ST1 = HST(KD)  - LTL(KD)*NU*(QST(KD)-QOS)
      ST2 = LTL(KD)  * VTF(KD)
      TEM5 = (QLS + QIS) * eta(kd1)
      ST1  = HALF * (TX1-ETA(KD1)*ST1-ST2*(DET-TEM5))*DLB(KD)
!

      WFN = WFN + ST1
      AKM = AKM - min(ST1,ZERO)   ! Commented on 08/26/02 - does not include top
!

      BUY(KD) = ST1 / (ETA(KD1)*qrb(kd))

      DET = DETP
      HCC = HCCP
      AKM = AKM / WFN


!***********************************************************************
!
!     If only to calculate workfunction save it and return
!
      IF (WRKFUN) THEN
        IF (WFN .GE. 0.0) WFNC = WFN
             !!!write(*,*) " return con work fun"
        RETURN
      ELSEIF (.NOT. CRTFUN) THEN
        ACR = WFNC
      ENDIF
!
!===>  THIRD CHECK BASED ON CLOUD WORKFUNCTION
!
      CALCUP = .FALSE.

      TEM  =  MIN(CD*200.0, MAX_NEG_BOUY)
      IF (WFN .GT. ACR .AND.  (.NOT. UNSAT)                             &
     & .and. dpneg .lt. 150.0  .AND. AKM .LE. TEM) THEN

        CALCUP = .TRUE.
      ENDIF

!===>  IF NO SOUNDING MEETS THIRD CONDITION, RETURN

      IF (.NOT. CALCUP) then 
          !!!write(*,*) " return con calcup"
       RETURN
      endif
!
      CLP = CLP * RHC
      do l=kd,kb1
        rnn(l) = rns(l)
      enddo
      DO L=KBL,K 
        RNN(L) = 0.0 
      ENDDO
!
!     If downdraft is to be invoked, do preliminary check to see
!     if enough rain is available and then call DDRFT.
!
      DDFT = .FALSE.
      IF (DNDRFT) THEN
!
        TRAIN = 0.0
        IF (CLP .GT. 0.0) THEN
          DO L=KD,KB1
            TRAIN = TRAIN + RNN(L)
          ENDDO
        ENDIF

        PL = (PRL(KD1) + PRL(KD))*HALF
        TEM = PRL(K+1)*(1.0-DPD*0.001)
        IF (TRAIN .GT. 1.0E-4 .AND. PL .LE. TEM) DDFT  = .TRUE.
!
      ENDIF
!

      IF (DDFT) THEN
!
!     Call Downdraft scheme based on (Cheng and Arakawa, 1997)
!
        CALL DDRFT(                                                     &
     &              K, KD                                               &
     &,             TLA, ALFIND                                         &
     &,             TOL, QOL, HOL,   PRL, QST, HST, GAM, GAF, HBL, QBL  &
     &,             QRB, QRT, BUY,   KBL, IDH, ETA, RNN, ETAI           &
     &,             ALM, WFN, TRAIN, DDFT                               &
     &,             ETD, HOD, QOD,   EVP, DOF, CLDFR, ETZ               &
     &,             GMS, GSD, GHD,lprnt)               

      ENDIF
!
!  No Downdraft case (including case with no downdraft soln)
!  ---------------------------------------------------------
!
      IF (.NOT. DDFT) THEN
        DO L=KD,K+1
          ETD(L) = 0.0
          HOD(L) = 0.0
          QOD(L) = 0.0
        ENDDO
        DO L=KD,K
          EVP(L) = 0.0
          ETZ(L) = 0.0
        ENDDO

      ENDIF
!
!
!===> CALCULATE GAMMAS  i.e. TENDENCIES PER UNIT CLOUD BASE MASSFLUX
!           Includes downdraft terms!

      avh = 0.0

!
!     Fraction of detrained condensate evaporated
!
      tem1 = max(ZERO, min(HALF, (prl(kd)-FOUR_P2)*ONE_M2))
      tem1 = 0.0
!
      tem2    = 1.0 - tem1
      TEM = DET * QIL(KD)


      st1 = (HCC+ALHF*TEM-ETA(KD)*HST(KD)) / (1.0+gam(KD))
      DS  = ETA(KD1) * (HOS- HOL(KD)) - ALHL*(QOS - QOL(KD))
      DH  = ETA(KD1) * (HOS- HOL(KD))


      GMS(KD) = (DS + st1 - tem1*det*alhl-tem*alhf) * PRI(KD)
      GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOS + DH)


!
!      TENDENCY FOR SUSPENDED ENVIRONMENTAL ICE AND/OR LIQUID WATER
!

      QIL(KD) =     (tem2*TEM + ETA(KD1)*(QIS-CIL(KD))                  &
     &                        - ETA(KD)*QIS ) * PRI(KD)
      QLL(KD) = (tem2*(DET-TEM) + ETA(KD1)*(QLS-CLL(KD))                &
     &                          - ETA(KD)*QLS ) * PRI(KD)
!
      GHD(KD) = 0.0
      GSD(KD) = 0.0
!
      DO L=KD1,K
       ST1 = ONE - ALFINT(L)
       ST2 = ONE - ALFINQ(L)
       ST3 = ONE - ALFIND(L)

         HB       = ALFINT(L)*HOL(L-1) + ST1*HOL(L)
         QB       = ALFINT(L)*QOL(L-1) + ST1*QOL(L)

         TEM      = ALFINQ(L)*CIL(L-1) + ST2*CIL(L)
         TEM2     = ALFINQ(L)*CLL(L-1) + ST2*CLL(L)
 
         TEM1     = ETA(L) * (TEM - CIL(L))
         TEM3     = ETA(L) * (TEM2 - CLL(L))

         HBD      = ALFIND(L)*HOL(L-1) + ST3*HOL(L)
         QBD      = ALFIND(L)*QOL(L-1) + ST3*QOL(L)

         TEM5     = ETD(L) * (HOD(L) - HBD)
         TEM6     = ETD(L) * (QOD(L) - QBD)
!
         DH       = ETA(L) * (HB - HOL(L)) + TEM5
         DS       = DH - ALHL * (ETA(L) * (QB - QOL(L)) + TEM6)

         GMH(L)   = DH * PRI(L)
         GMS(L)   = DS * PRI(L)

!
         GHD(L)   = TEM5 * PRI(L)
         GSD(L)   = (TEM5 - ALHL * TEM6) * PRI(L)
!
         QIL(L)   = TEM1 * PRI(L)
         QLL(L)   = TEM3 * PRI(L)

         TEM1     = ETA(L) * (CIL(L-1) - TEM)
         TEM3     = ETA(L) * (CLL(L-1) - TEM2)

         DH       = ETA(L) * (HOL(L-1) - HB) - TEM5
         DS       = DH - ALHL * ETA(L) * (QOL(L-1) - QB)                &
     &                 + ALHL * (TEM6 - EVP(L-1))

         GMH(L-1) = GMH(L-1) + DH * PRI(L-1)
         GMS(L-1) = GMS(L-1) + DS * PRI(L-1)
!
!
         GHD(L-1) = GHD(L-1) - TEM5 * PRI(L-1)
         GSD(L-1) = GSD(L-1) - (TEM5-ALHL*(TEM6-EVP(L-1))) * PRI(L-1)

         QIL(L-1) = QIL(L-1) + TEM1 * PRI(L-1)
         QLL(L-1) = QLL(L-1) + TEM3 * PRI(L-1)

        avh = avh + gmh(l-1)*(prs(l)-prs(l-1))

      ENDDO
!

      HBD  = HOL(K)
      QBD  = QOL(K)
      TEM5 =  ETD(K+1) * (HOD(K+1) - HBD)
      TEM6 =  ETD(K+1) * (QOD(K+1) - QBD)
      DH   = - TEM5
      DS   = DH  + ALHL * TEM6
      TEM1 = DH * PRI(K)
      TEM2 = (DS - ALHL * EVP(K)) * PRI(K)
      GMH(K) = GMH(K) + TEM1
      GMS(K) = GMS(K) + TEM2
      GHD(K) = GHD(K) + TEM1
      GSD(K) = GSD(K) + TEM2

!
      avh = avh + gmh(K)*(prs(KP1)-prs(K))
!
      tem4   = - GRAVFAC * pris
      TX1    = DH * tem4
      TX2    = DS * tem4
!
      DO L=KBL,K
        GMH(L) = GMH(L) + TX1
        GMS(L) = GMS(L) + TX2
        GHD(L) = GHD(L) + TX1
        GSD(L) = GSD(L) + TX2
!
        avh = avh + tx1*(prs(l+1)-prs(l))
      ENDDO

!
!***********************************************************************
!***********************************************************************

!===>  KERNEL (AKM) CALCULATION BEGINS

!===>  MODIFY SOUNDING WITH UNIT MASS FLUX
!
!     TESTMB = 0.01

      DO L=KD,K

         TEM1   = GMH(L)
         TEM2   = GMS(L)
         HOL(L) = HOL(L) +  TEM1*TESTMB
         QOL(L) = QOL(L) + (TEM1-TEM2)  * (TESTMB/ALHL)
         HST(L) = HST(L) +  TEM2*(ONE+GAM(L))*TESTMB
         QST(L) = QST(L) +  TEM2*GAM(L)*(TESTMB/ALHL)
         CLL(L) = CLL(L) + QLL(L) * TESTMB
         CIL(L) = CIL(L) + QIL(L) * TESTMB
      ENDDO
!

      if (alm .gt. 0.0) then
      HOS = HOS + GMH(KD)  * TESTMB
      QOS = QOS + (GMH(KD)-GMS(KD)) * (TESTMB/ALHL)

        QLS     = QLS + QLL(KD) * TESTMB
        QIS     = QIS + QIL(KD) * TESTMB
      else
        st2 = 1.0 - st1s
        HOS = HOS + (st1s*GMH(KD)+st2*GMH(KD1))  * TESTMB
        QOS = QOS + (st1s * (GMH(KD)-GMS(KD))                           &
     &            +  st2  * (GMH(KD1)-GMS(KD1))) * (TESTMB/ALHL)
        HST(kd) = HST(kd) + (st1s*GMS(kd)*(ONE+GAM(kd))                 &
     &                    +  st2*gms(kd1)*(ONE+GAM(kd1))) * TESTMB
        QST(kd) = QST(kd) + (st1s*GMS(kd)*GAM(kd)                       &
     &                    +  st2*gms(kd1)*gam(kd1)) * (TESTMB/ALHL)

        QLS     = QLS + (st1s*QLL(KD)+st2*QLL(KD1)) * TESTMB
        QIS     = QIS + (st1s*QIL(KD)+st2*QIL(KD1)) * TESTMB
      endif

!
      TEM = PRL(K+1) - PRL(K)
      HBL = HOL(K) * TEM
      QBL = QOL(K) * TEM
      QLB = CLL(K) * TEM
      QIB = CIL(K) * TEM
      DO L=KM1,KBL,-1
        TEM = PRL(L+1) - PRL(L)
        HBL = HBL + HOL(L) * TEM
        QBL = QBL + QOL(L) * TEM
        QLB = QLB + CLL(L) * TEM
        QIB = QIB + CIL(L) * TEM
      ENDDO
      HBL = HBL * PRIS
      QBL = QBL * PRIS
      QLB = QLB * PRIS
      QIB = QIB * PRIS


!***********************************************************************

!===>  CLOUD WORKFUNCTION FOR MODIFIED SOUNDING, THEN KERNEL (AKM)
!
      AKM = ZERO
      TX1 = ZERO
      QTL = QST(KB1) - GAF(KB1)*HST(KB1)
      QTV = QBL
      HCC = HBL
      TX2 = HCC
      TX4 = (ALHF*0.5)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(KB1))*TCRF))
!
      qtv   = qbl
      tx1   = qib + qlb
!

      DO L=KB1,KD1,-1
         DEL_ETA = ETA(L) - ETA(L+1)
         HCCP = HCC + DEL_ETA*HOL(L)
!
         QTLP = QST(L-1) - GAF(L-1)*HST(L-1)
         QTVP = 0.5 * ((QTLP+QTL)*ETA(L)                                &
     &                +(GAF(L)+GAF(L-1))*HCCP)

         DETP = (BKC(L)*TX1 - (QTVP-QTV)                                &
     &        +  DEL_ETA*(QOL(L)+CLL(L)+CIL(L))                         &
     &        +  ETA(L)*Q0U(L) + ETA(L+1)*Q0D(L)) * AKC(L)
         IF (DETP .LE. ZERO) UNSAT = .TRUE.

         ST1  = HST(L) - LTL(L)*NU*(QST(L)-QOL(L))

         TEM2 = (ALHF*0.5)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(L-1))*TCRF))
         TEM1 = HCCP + DETP * (TEM2+TX4)

         ST2  = LTL(L) * VTF(L)
         TEM5 = CLL(L) + CIL(L)
         AKM  = AKM +                                                   &
     &     (  (TX2  -ETA(L+1)*ST1-ST2*(TX1-TEM5*eta(l+1))) * DLB(L)     &
     &      + (TEM1 -ETA(L  )*ST1-ST2*(DETP-TEM5*eta(l)))  * DLT(L) )
!
         HCC  = HCCP
         TX1  = DETP
         TX2  = TEM1
         QTL  = QTLP
         QTV  = QTVP
         TX4  = TEM2
      ENDDO
!
      if (unsat) then 
         !!!write(*,*) " return con AKM loop "
       return
      endif
!
!  Eventhough we ignore the change in lambda, we still assume
!  that the cLoud-top contribution is zero; as though we still
!  had non-bouyancy there.
!
!
      ST1 = HST(KD)  - LTL(KD)*NU*(QST(KD)-QOS)
      ST2 = LTL(KD)  * VTF(KD)
      TEM5 = (QLS + QIS) * eta(kd1)
      AKM  = AKM + HALF * (TX2-ETA(KD1)*ST1-ST2*(TX1-TEM5)) * DLB(KD)
!
      AKM = (AKM - WFN) * (ONE/TESTMB)


!***********************************************************************

!===>   MASS FLUX

      tem2 = rel_fac
!
      AMB = - (WFN-ACR) / AKM
!
      if(lprnt) print *,' wfn=',wfn,' acr=',acr,' akm=',akm             &
     &,' amb=',amb,' KD=',kd,' cldfrd=',cldfrd,' tem2=',tem2            &
     &,' rel_fac=',rel_fac,' prskd=',prs(kd)

!===>   RELAXATION AND CLIPPING FACTORS
!
      AMB = AMB * CLP * tem2

       
!===>   SUB-CLOUD LAYER DEPTH LIMIT ON MASS FLUX



        !!!write(*,*) 'PRL(KP1)-PRL(KBL)', PRL(KP1)-PRL(KBL)
        !!!write(*,*) 'FRACBL           ', FRACBL
        !!!write(*,*) 'GRAVCON          ', GRAVCON

      AMBMAX = (PRL(KP1)-PRL(KBL))*(FRACBL*GRAVCON)
      AMB    = MAX(MIN(AMB, AMBMAX),ZERO)


!***********************************************************************
!*************************RESULTS***************************************
!***********************************************************************

!===>  PRECIPITATION AND CLW DETRAINMENT
!
      avt = 0.0
      avq = 0.0
      avr = dof

!
      DSFC = DSFC + AMB * ETD(K) * (1.0/DT)
!
      DO L=K,KD,-1
          PCU(L) = PCU(L) + AMB*RNN(L)      !  (A40)
          avr = avr + rnn(l)
      ENDDO
!
!===> TEMPARATURE AND Q CHANGE AND CLOUD MASS FLUX DUE TO CLOUD TYPE KD
!
      TX1 = AMB * (ONE/CP)
      TX2 = AMB * (ONE/ALHL)
      DO L=KD,K
        ST1    = GMS(L)*TX1
        TOI(L) = TOI(L) + ST1
        TCU(L) = TCU(L) + ST1
        TCD(L) = TCD(L) + GSD(L) * TX1
!
        st1 = st1 - (alhl/cp) * (QIL(L) + QLL(L)) * AMB

        avt = avt + st1 * (prs(l+1)-prs(l))

        FLX(L) = FLX(L) + ETA(L)*AMB
!
        QII(L) = QII(L) + QIL(L) * AMB
        TEM    = 0.0

        QLI(L) = QLI(L) + QLL(L) * AMB + TEM

        ST1    = (GMH(L)-GMS(L)) * TX2

        QOI(L) = QOI(L) + ST1
        QCU(L) = QCU(L) + ST1
        QCD(L) = QCD(L) + (GHD(L)-GSD(L)) * TX2
!
        avq = avq + (st1+(QLL(L)+QIL(L))*amb) * (prs(l+1)-prs(l))

      ENDDO
      avr = avr * amb
!
!
      TX1 = 0.0
      TX2 = 0.0
!
!     REEVAPORATION OF FALLING CONVECTIVE RAIN
!
      IF (REVAP) THEN
!
       tem = 0.0
       do l=kd,kbl
!        tem = tem + pcu(l)
         IF (L .lt. IDH .or. (.not. DDFT)) THEN
           tem = tem + amb * rnn(l)
         endif
       enddo
       tem = tem + amb * dof
       tem = tem * (3600.0/dt)
       tem1 = max(1.0d0 , min(100.0d0  , (5.0E10/max(garea,one))))

       clfrac = max(ZERO, min(ONE, rknob*clf(tem)*tem1))

!
       DO L=KD,K
!                                                 for L=KD,K
         IF (L .GE. IDH .AND. DDFT) THEN
           TX2    = TX2 + AMB * RNN(L)
           CLDFRD = MIN(AMB*CLDFR(L), clfrac)
         ELSE
           TX1 = TX1 + AMB * RNN(L)
         ENDIF
         tx4 = zfac * phil(l)
         tx4 = (one - tx4 * (one - half*tx4)) * afc
!

         IF (TX1 .GT. 0. .OR. TX2 .GT. 0.0) THEN
          TEQ     = TOI(L)
          QEQ     = QOI(L)
          PL      = 0.5 * (PRL(L+1)+PRL(L))

          ST1     = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
          ST2     = ST1*ELFOCP + (1.0-ST1)*ELOCP

          CALL QSATCN ( TEQ,PL,QSTEQ,DQDT,.false.)
!
!
          DELTAQ = 0.5 * (QSTEQ*rhc_ls(l)-QEQ) / (1.+ST2*DQDT)
!
          QEQ    = QEQ + DELTAQ
          TEQ    = TEQ - DELTAQ*ST2
!
          TEM1   = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
          TEM2   = TEM1*ELFOCP + (1.0-TEM1)*ELOCP

          CALL QSATCN ( TEQ,PL,QSTEQ,DQDT,.false.)
!
!
          DELTAQ = (QSTEQ*rhc_ls(l)-QEQ) / (1.+TEM2*DQDT)
!
          QEQ    = QEQ + DELTAQ
          TEQ    = TEQ - DELTAQ*TEM2

          IF (QEQ .GT. QOI(L)) THEN
            POTEVAP = (QEQ-QOI(L))*(PRL(L+1)-PRL(L))*GRAVCON

            tem4    = 0.0
            if (tx1 .gt. 0.0)                                           &
     &      TEM4    = POTEVAP * (1. - EXP( tx4*TX1**0.57777778 ) )
            ACTEVAP = MIN(TX1, TEM4*CLFRAC)

            if (tx1 .lt. rainmin*dt) actevap = min(tx1, potevap)
!
            tem4    = 0.0
            if (tx2 .gt. 0.0)                                           &
     &      TEM4    = POTEVAP * (1. - EXP( tx4*TX2**0.57777778 ) )

            TEM4    = min(MIN(TX2, TEM4*CLDFRD), potevap-actevap)
            if (tx2 .lt. rainmin*dt) tem4 = min(tx2, potevap-actevap)
!
            TX1     = TX1 - ACTEVAP
            TX2     = TX2 - TEM4
            ST1     = (ACTEVAP+TEM4) * PRI(L)
            QOI(L)  = QOI(L) + ST1
            QCU(L)  = QCU(L) + ST1
!

            ST1     = ST1 * ELOCP
            TOI(L)  = TOI(L) - ST1 
            TCU(L)  = TCU(L) - ST1
          ENDIF
         ENDIF
       ENDDO
!
       CUP = CUP + TX1 + TX2 + DOF * AMB
      ELSE
       DO L=KD,K
         TX1 = TX1 + AMB * RNN(L)
       ENDDO
       CUP = CUP + TX1 + DOF * AMB
      ENDIF

!
!    MIXING OF PASSIVE TRACERS
!
      DO N=1,M

       DO L=KD,K
         HOL(L) = ROI(L,N)
       ENDDO
!
        HCC     = RBL(N)
        HOD(KD) = HOL(KD)
!      Compute downdraft properties for the tracer
       DO L=KD1,K
         ST1 = ONE - ALFIND(L)
         HB  = ALFIND(L)  * HOL(L-1) + ST1 * HOL(L)
         IF (ETZ(L-1) .NE. 0.0) THEN
           DEL_ETA = ETD(L) - ETD(L-1)
           TEM     = 1.0 / ETZ(L-1)
           IF (DEL_ETA .GT. 0.0) THEN
             HOD(L) = (ETD(L-1)*(HOD(L-1)-HOL(L-1))                     &
     &                +  ETD(L)  *(HOL(L-1)-HB)                         &
     &                +  ETZ(L-1)*HB) * TEM
           ELSE
             HOD(L) = (ETD(L-1)*(HOD(L-1)-HB) + ETZ(L-1)*HB) * TEM
           ENDIF
         ELSE
           HOD(L) = HB
         ENDIF
       ENDDO
             
       DO L=KB1,KD,-1
          HCC = HCC + (ETA(L)-ETA(L+1))*HOL(L)
       ENDDO
!
       GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOL(KD))
       DO L=KD1,K
        ST1 = ONE - ALFINT(L)
        ST2 = ONE - ALFIND(L)
           HB       = ALFINT(L) * HOL(L-1) + ST1 * HOL(L)
           HBD      = ALFIND(L) * HOL(L-1) + ST2 * HOL(L)
           TEM5     = ETD(L)    * (HOD(L) - HBD)
           DH       = ETA(L)    * (HB - HOL(L)) + TEM5
           GMH(L  ) = DH * PRI(L) * trcfac(n)
           DH       = ETA(L)    * (HOL(L-1) - HB) - TEM5
           GMH(L-1) = GMH(L-1)  + DH * PRI(L-1) * trcfac(n)
       ENDDO
!
       DO L=KD,K
         ST1      = GMH(L)*AMB
         ROI(L,N) = HOL(L)   + ST1
         RCU(L,N) = RCU(L,N) + ST1
       ENDDO
      ENDDO                             ! Tracer loop M

!***********************************************************************
!***********************************************************************
!***********************************************************************

             !!!write(*,*) " return con reached end of cloud "

      RETURN
      END

      SUBROUTINE DDRFT(                                                 &
     &                  K, KD                                           &
     &,                 TLA, ALFIND                                     &
     &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF, HBL, QBL&
     &,                 QRB, QRT, BUY, KBL, IDH, ETA, RNN, ETAI         &
     &,                 ALM, WFN, TRAIN, DDFT                           &
     &,                 ETD, HOD, QOD, EVP, DOF, CLDFRD, WCB            &
     &,                 GMS, GSD, GHD,lprnt)                   
!    &,                 GMS, GSD, GHD,lprnt)
!    &,                 TX1, TX2, TX3, TX4, TX5, TX6, TX7, TX8, TX9)

!
!***********************************************************************
!******************** Cumulus Downdraft Subroutine *********************
!****************** Based on Cheng and Arakawa (1997)  ****** **********
!************************ SUBROUTINE DDRFT  ****************************
!*************************  October 2004  ******************************
!***********************************************************************
!***********************************************************************
!************* Shrinivas.Moorthi@noaa.gov (301) 763 8000(X7233) ********
!***********************************************************************
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!===>  TOL(K)     INPUT   TEMPERATURE            KELVIN
!===>  QOL(K)     INPUT   SPECIFIC HUMIDITY      NON-DIMENSIONAL

!===>  PRL(K+1)   INPUT   PRESSURE @ EDGES       MB

!===>  K     INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
!===>  KD    INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
!     
      use module_ras
      IMPLICIT NONE
!
!  INPUT ARGUMENTS
!
      INTEGER K, KD
      real(kind=kind_phys) ALFIND(K)
      INTEGER KBL, KB1



      LOGICAL SKPDD, SKPUP

      real(kind=kind_phys) HOL(KD:K),   QOL(KD:K),   GAF(KD:K+1)        &
     &,                    HST(KD:K),   QST(KD:K),   TOL(KD:K)          &
     &,                    BUY(KD:K+1), QRB(KD:K),   QRT(KD:K)          &
     &,                    GAM(KD:K+1), RNN(KD:K),   RNS(KD:K)                       &
     &,                    ETA(KD:K+1), PRL(KD:K+1), ETAI(KD:K)
!
      real(kind=kind_phys)    HBL,     QBL,        PRIS                 &
     &,                       TRAIN,   WFN,        ALM
!
!     TEMPORARY WORK SPACE
!
      real(kind=kind_phys) GMS(KD:K+1)
      real(kind=kind_phys) TX1,    TX2,  TX3, TX4                       &
     &,                    TX5,    TX6,  TX7, TX8, TX9
      LOGICAL UNSAT

      real(kind=kind_phys) TL, PL, QL, QS, DQS, ST1,  HB, QB, TB        &
     &,                    QQQ, PICON, PIINV, DEL_ETA                   &
     &,                    TEM, TEM1, TEM2, TEM3, TEM4, ST2             &
     &,                    ERRMIN, ERRMI2, ERRH, ERRW, ERRE, TEM5       &
     &,                    TEM6, HBD, QBD
      INTEGER I, L,  N, IX, KD1, II                                     &
     &,       KP1, IT, KM1, KTEM, KK, KK1, LM1, LL, LP1                 &
     &,       IP1, JJ, ntla

!
      integer, parameter :: NUMTLA=2
!     integer, parameter :: NUMTLA=4
      parameter (ERRMIN=0.0001, ERRMI2=0.1*ERRMIN)
!     parameter (ERRMIN=0.00001, ERRMI2=0.1*ERRMIN)
!
      real(kind=kind_phys) TLA,    STLA,  CTL2, CTL3
      real(kind=kind_phys) GMF,    PI,    ONPG, CTLA, VTRM, VTPEXP      &
     &,    RPART,  QRMIN,  AA1,    BB1,   CC1,  DD1                     &
     &,    WC2MIN, WCMIN,  WCBASE, F2,    F3,   F5, GMF1, GMF5          &
     &,    QRAF,   QRBF,   CMPOR , del_tla               
!    &,    sialf
!
      parameter (ONPG=1.0+0.5, GMF=1.0/ONPG, RPART=0.0)
!     parameter (ONPG=1.0+0.5, GMF=1.0/ONPG, RPART=1.0)
!     parameter (ONPG=1.0+0.5, GMF=1.0/ONPG, RPART=0.5)
!     PARAMETER (AA1=1.0, BB1=1.5, CC1=1.1, DD1=0.85, F3=CC1, F5=2.5)
!     PARAMETER (AA1=2.0, BB1=1.5, CC1=1.1, DD1=0.85, F3=CC1, F5=2.5)
      PARAMETER (AA1=1.0, BB1=1.0, CC1=1.0, DD1=1.0, F3=CC1,  F5=1.0)
      parameter (QRMIN=1.0E-6, WC2MIN=0.01, GMF1=GMF/AA1, GMF5=GMF/F5)
!     parameter (QRMIN=1.0E-6, WC2MIN=1.00, GMF1=GMF/AA1, GMF5=GMF/F5)
!     parameter (sialf=0.5)
!
      PARAMETER (PI=3.1415926535897931, PIINV=1.0/PI)
      INTEGER ITR, ITRMU, ITRMD, KTPD, ITRMIN, ITRMND
!     PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=7)
      PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=12, ITRMND=12)
!     PARAMETER (ITRMU=25, ITRMD=25, ITRMIN=12)
!     PARAMETER (ITRMU=14, ITRMD=18, ITRMIN=7)
!     PARAMETER (ITRMU=10, ITRMD=10, ITRMIN=5)
      real(kind=kind_phys) QRP(KD:K+1), WVL(KD:K+1), AL2
      real(kind=kind_phys) WVLO(KD:K+1)
!
      real(kind=kind_phys) RNF(KD:K),   ETD(KD:K+1), WCB(KD:K)          &
     &,                    HOD(KD:K+1), QOD(KD:K+1), EVP(KD:K)          &
     &,                    ROR(KD:K+1), STLT(KD:K)                      &
     &,                    GHD(KD:K),   GSD(KD:K),   CLDFRD(KD:K)       &
     &,                    RNT,        RNB                              &
     &,                    ERRQ,       RNTP
      INTEGER IDW, IDH, IDN(K), idnm
      real(kind=kind_phys) ELM(K)
!     real(kind=kind_phys) EM(K*K), ELM(K)
      real(kind=kind_phys) EDZ, DDZ, CE, QHS, FAC, FACG, ASIN,          &
     &                     RSUM1, RSUM2, RSUM3, CEE
      LOGICAL DDFT, UPDRET, DDLGK
!
      real(kind=kind_phys) AA(KD:K,KD:K+1), QW(KD:K,KD:K)               &
     &,                    BUD(KD:K), VT(2), VRW(2), TRW(2)             &
     &,                    GQW(KD:K)                                    &
     &,                    QA(3),     WA(3),    DOF, DOFW               &
     &,                    QRPI(KD:K), QRPS(KD:K)
!    &,                    GQW(KD:K), WCB(KD:K)

!***********************************************************************

      real(kind=kind_phys) QRPF, VTPF
      logical lprnt
!!CFPP$ EXPAND (QRPF, QRABF, VTPF)
!!CFPP$ NOCONCUR R

!

!     if(lprnt) print *,' K=',K,' KD=',KD,' In Downdrft'

      KD1    = KD + 1
      KP1    = K  + 1
      KM1    = K  - 1
      KB1    = KBL - 1
!
      CMPOR  = CMB2PA / RGAS
!
!     VTP    = 36.34*SQRT(1.2)* (0.001)**0.1364
      VTPEXP = -0.3636
!     PIINV  = 1.0 / PI
      PICON  = PI * ONEBG * 0.5
!
!
!     Compute Rain Water Budget of the Updraft (Cheng and Arakawa, 1997)
!
      CLDFRD = 0.0
      RNTP   = 0.0
      DOF    = 0.0
      ERRQ   = 10.0
      RNB    = 0.0
      RNT    = 0.0
      TX2    = PRL(KBL)
!
      TX1      = (PRL(KD) + PRL(KD1)) * 0.5
      ROR(KD)  = CMPOR*TX1 / (TOL(KD)*(1.0+NU*QOL(KD)))
!     GMS(KD)  = VTP * ROR(KD) ** VTPEXP
      GMS(KD)  = VTP * VTPF(ROR(KD))
!
      QRP(KD)  = QRMIN
!
      TEM      = TOL(K) * (1.0 + NU * QOL(K))
      ROR(K+1) = 0.5 * CMPOR * (PRL(K+1)+PRL(K)) / TEM
      GMS(K+1) = VTP * VTPF(ROR(K+1))
      QRP(K+1) = QRMIN
!!    BUY(KD)  = MAX(BUY(KD),ONE_M1)
!     BUY(KD)  = MAX(BUY(KD), 0.1)
!     BUY(KD)  = MAX(BUY(KD), 0.0)
!
      kk = kbl
      DO L=KD1,K
        TEM = 0.5 * (TOL(L)+TOL(L-1))                                   &
     &      * (1.0 + (0.5*NU) * (QOL(L)+QOL(L-1)))
        ROR(L) = CMPOR * PRL(L) / TEM
!       GMS(L) = VTP * ROR(L) ** VTPEXP
        GMS(L) = VTP * VTPF(ROR(L))
        QRP(L) = QRMIN
!!      BUY(L) = MAX(BUY(L),ONE_M1)
!       BUY(L) = MAX(BUY(L), 0.1)
!       BUY(L) = MAX(BUY(L), 1.0E-5)
        if (buy(l) .le. 0.0 .and. kk .eq. KBL) then
          kk = l
!       if (buy(l) .le. 0.0) then
!         if (buy(l-1) .gt. 0.0 .and. buy(l+1) .gt. 0.0) then
!           buy(l) = 0.5 * (buy(l+1) + buy(l-1))
!         elseif (buy(l-1) .gt. 0.0) then
!           buy(l) = 0.5*buy(l-1)
!           buy(l) = 0.25 * buy(l-1)
!         else
!            BUY(L) = 1.0E-4
!            BUY(L) = 5.0E-4
!            BUY(L) = 1.0E-5
!         endif
        endif
!       BUY(L) = MAX(BUY(L), 1.0E-4)
!       BUY(L) = MAX(BUY(L), 1.0E-5)
!       BUY(L) = MAX(BUY(L), 5.0E-4)
      ENDDO
      if (kk .ne. kbl) then
        do l=kk,kbl
          buy(l) = 0.9 * buy(l-1)
        enddo
      endif
!
      do l=kd,k
        qrpi(l) = buy(l)
      enddo
      do l=kd1,kb1
        buy(l) = 0.25 * (qrpi(l-1)+qrpi(l)+qrpi(l)+qrpi(l+1))
!       tem = 0.5 * (eta(l)+eta(l+1))
!       buy(l) = buy(l) * tem * tem
      enddo
!     tem = 0.5 * (eta(KD)+eta(kd1))
!     buy(kd) = buy(kd) * tem * tem
      
!
!     CALL ANGRAD(TX1, ALM, STLA, CTL2, AL2, PI, TLA, TX2, WFN, TX3)
      tx1 = 1000.0 + tx1 - prl(k+1)
      CALL ANGRAD(TX1, ALM,  AL2, TLA, TX2, WFN, TX3)
!
!    Following Ucla approach for rain profile
!
      F2      = 2.0*BB1*ONEBG/(PI*0.2)
      WCMIN   = SQRT(WC2MIN)
      WCBASE  = WCMIN
!
!     del_tla = TLA * 0.2
!     del_tla = TLA * 0.25
      del_tla = TLA * 0.3
      TLA     = TLA - DEL_TLA
!
!     do ntla=1,numtla
!
!     if (errq .lt. 1.0 .or. tla .gt. 45.0) cycle
!
!     tla = tla + del_tla
!     STLA = SIN(TLA*PI/180.0)
!     CTL2 = 1.0 - STLA * STLA
!
!     if (lprnt) print *,' tla=',tla,' al2=',al2,' ptop='
!    &,0.5*(prl(kd)+prl(kd1)),' ntla=',ntla
!     if (lprnt) print *,' buy=',(buy(l),l=kd,kbl)
!
!     STLA = F2     * STLA * AL2
!     CTL2 = DD1    * CTL2
!     CTL3 = 0.1364 * CTL2
!
      DO L=KD,K
        RNF(L)   = 0.0
        RNS(L)   = 0.0
        WVL(L)   = 0.0
        STLT(L)  = 0.0
        GQW(L)   = 0.0
        QRP(L)   = QRMIN
        DO N=KD,K
          QW(N,L) = 0.0
        ENDDO
      ENDDO
!
!-----QW(N,L) = D(W(N)*W(N))/DQR(L)
!
      KK = KBL
!     WVL(KK)    = WCBASE
      QW(KD,KD)  = -QRB(KD)  * GMF1
      GHD(KD)    = ETA(KD)   * ETA(KD)
      GQW(KD)    = QW(KD,KD) * GHD(KD)
      GSD(KD)    = ETAI(KD)  * ETAI(KD)
!     GSD(KD)    = 1.0 / GHD(KD)
!
      GQW(KK)    = -  QRB(KK-1) * (GMF1+GMF1)
!
      WCB(KK)    = WCBASE * WCBASE
!     WVL(KK)    = WCBASE
!     STLT(KBL)  = 1.0 / WCBASE

      TX1        = WCB(KK)
      GSD(KK)    = 1.0
      GHD(KK)    = 1.0
!
      TEM        = GMF1 + GMF1
!!    TX1        = WCB(KK) + buy(kb1)*tem*qrb(kb1)
      DO L=KB1,KD1,-1
         GHD(L)  = ETA(L)  * ETA(L)
         GSD(L)  = ETAI(L) * ETAI(L)
!        GSD(L)  = 1.0 / GHD(L)
         GQW(L)  = - GHD(L) * (QRB(L-1)+QRT(L)) * TEM
         QW(L,L) = - QRT(L) * TEM
!
!        TX1     = TX1 + BUY(L) * TEM
!!       TX1     = TX1 + BUY(L) * TEM * (qrb(l-1)+qrt(l)) * ghd(l)
         st1     = 0.5 * (eta(l) + eta(l+1))
         TX1     = TX1 + BUY(L) * TEM * (qrb(l)+qrt(l)) * st1 * st1
         WCB(L)  = TX1 * GSD(L)
      ENDDO
!
      TEM1        = (QRB(KD) + QRT(KD1) + QRT(KD1)) * GMF1
      GQW(KD1)    = - GHD(KD1) * TEM1
!     QW(L,KD1)   = - QRT(KD1) * TEM
      QW(KD1,KD1) = - QRT(KD1) * TEM
!     WCB(KD)     = (TX1 + BUY(KD)*TEM) * GSD(KD)
      st1     = 0.5 * (eta(kd) + eta(kd1))
      WCB(KD)     = (TX1 + BUY(KD)*TEM*qrb(kd)*st1*st1) * GSD(KD)
!
      DO L=KD1,KBL
        DO N=KD,L-1
           QW(N,L) = GQW(L) * GSD(N)
        ENDDO
      ENDDO
      QW(KBL,KBL) = 0.0
!
!     WVL(KBL)    = WCBASE
!     STLT(KBL)   = 1.0 / WCBASE
!
!
      do ntla=1,numtla
!
!     if (errq .lt. 1.0 .or. tla .gt. 45.0) cycle
      if (errq .lt. 0.1 .or. tla .gt. 45.0) cycle
!
      tla = tla + del_tla
      STLA = SIN(TLA*PI/180.0)
      CTL2 = 1.0 - STLA * STLA
!
!     if (lprnt) print *,' tla=',tla,' al2=',al2,' ptop='
!    &,0.5*(prl(kd)+prl(kd1)),' ntla=',ntla,' f2=',f2,' stla=',stla
!     if (lprnt) print *,' buy=',(buy(l),l=kd,kbl)
!
      STLA = F2     * STLA * AL2
      CTL2 = DD1    * CTL2
      CTL3 = 0.1364 * CTL2
!
      DO L=KD,K
        RNF(L)   = 0.0
        WVL(L)   = 0.0
        STLT(L)  = 0.0
        QRP(L)   = QRMIN
      ENDDO
      WVL(KBL)    = WCBASE
      STLT(KBL)   = 1.0 / WCBASE
!
      DO L=KD,K+1
        DO N=KD,K
          AA(N,L) = 0.0
        ENDDO
      ENDDO
!
      SKPUP = .FALSE.
!
      DO ITR=1,ITRMU               ! Rain Profile Iteration starts!
        IF (.NOT. SKPUP) THEN
           wvlo = wvl
!
!-----CALCULATING THE VERTICAL VELOCITY
!
          TX1      = 0.0
          QRPI(KBL) = 1.0 / QRP(KBL)
          DO L=KB1,KD,-1
            TX1     = TX1    + QRP(L+1) * GQW(L+1)
            ST1     = WCB(L) + QW(L,L)  * QRP(L)                        &
     &                       + TX1      * GSD(L)
!           if (st1 .gt. 0.0) then
            if (st1 .gt. wc2min) then
!             WVL(L)  = SQRT(ST1)
              WVL(L)  = 0.5 * (SQRT(ST1) + WVL(L))
!             if (itr .eq. 1) wvl(l) = wvl(l) * 0.25
            else
!     if (lprnt)  print *,' l=',l,' st1=',st1,' wcb=',wcb(l),' qw='
!    &,qw(l,l),' qrp=',qrp(l),' tx1=',tx1,' gsd=',gsd(l),' ite=',itr
!             wvl(l) = 0.5*(wcmin+wvl(l))
              wvl(l) = 0.5 * (wvl(l) + wvl(l+1))
              qrp(l) = 0.5 * ((wvl(l)*wvl(l)-wcb(l)-tx1*gsd(l))/qw(l,l) &
     &                     + qrp(l))
!!            wvl(l) = 0.5 * (wvl(l) + wvl(l+1))
            endif
!           wvl(l)  = 0.5 * (wvl(l) + wvlo(l))
!           WVL(L)  = SQRT(MAX(ST1,WC2MIN))
            wvl(l)  = max(wvl(l), wcbase)
            STLT(L) = 1.0 / WVL(L)
            QRPI(L) = 1.0 / QRP(L)
          ENDDO
!         qrps = qrp
!         do l=kd1,kb1
!           qrp(l) = 0.25 * (qrps(l-1)+qrps(l)+qrps(l)+qrps(l+1))
!           qrpi(l) = 1.0 / qrp(l)
!         enddo
!         qrpi(kd) = 1.0 / qrp(kd)
!
!     if (lprnt) then
!     print *,' ITR=',ITR,' ITRMU=',ITRMU
!     print *,' WVL=',(WVL(L),L=KD,KBL)
!     print *,' qrp=',(qrp(L),L=KD,KBL)
!     print *,' qrpi=',(qrpi(L),L=KD,KBL)
!     print *,' rnf=',(rnf(L),L=KD,KBL)
!     endif
!
!-----CALCULATING TRW, VRW AND OF
!
!         VT(1)   = GMS(KD) * QRP(KD)**0.1364
          VT(1)   = GMS(KD) * QRPF(QRP(KD))
          TRW(1)  = ETA(KD) * QRP(KD) * STLT(KD)
          TX6     = TRW(1) * VT(1)
          VRW(1)  = F3*WVL(KD) - CTL2*VT(1)
          BUD(KD) = STLA * TX6 * QRB(KD) * 0.5
          RNF(KD) = BUD(KD)
          DOF     = 1.1364 * BUD(KD) * QRPI(KD)
          DOFW    = -BUD(KD) * STLT(KD)
!
          RNT     = TRW(1) * VRW(1)
          TX2     = 0.0
          TX4     = 0.0
          RNB     = RNT
          TX1     = 0.5
          TX8     = 0.0
!
          IF (RNT .GE. 0.0) THEN
            TX3 = (RNT-CTL3*TX6) * QRPI(KD)
            TX5 = CTL2 * TX6 * STLT(KD)
          ELSE
            TX3 = 0.0
            TX5 = 0.0
            RNT = 0.0
            RNB = 0.0
          ENDIF
!
          DO L=KD1,KB1
            KTEM    = MAX(L-2, KD)
            LL      = L - 1
!
!           VT(2)   = GMS(L) * QRP(L)**0.1364
            VT(2)   = GMS(L) * QRPF(QRP(L))
            TRW(2)  = ETA(L) * QRP(L) * STLT(L)
            VRW(2)  = F3*WVL(L) - CTL2*VT(2)
            QQQ     = STLA * TRW(2) * VT(2)
            ST1     = TX1  * QRB(LL)
            BUD(L)  = QQQ * (ST1 + QRT(L))
!
            QA(2)   = DOF
            WA(2)   = DOFW
            DOF     = 1.1364 * BUD(L) * QRPI(L)
            DOFW    = -BUD(L) * STLT(L)
!
            RNF(LL) = RNF(LL) + QQQ * ST1
            RNF(L)  =           QQQ * QRT(L)
!
            TEM3    = VRW(1) + VRW(2)
            TEM4    = TRW(1) + TRW(2)
!
            TX6     = .25 * TEM3 * TEM4
            TEM4    = TEM4 * CTL3
!
!-----BY QR ABOVE
!
!           TEM1    = .25*(TRW(1)*TEM3 - TEM4*VT(1))*TX7
            TEM1    = .25*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
            ST1     = .25*(TRW(1)*(CTL2*VT(1)-VRW(2))                   &
     &                  * STLT(LL) + F3*TRW(2))
!-----BY QR BELOW
            TEM2    = .25*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
            ST2     = .25*(TRW(2)*(CTL2*VT(2)-VRW(1))                   &
     &                 * STLT(L)  + F3*TRW(1))
!
!      From top to  the KBL-2 layer
!
            QA(1)   = TX2
            QA(2)   = QA(2) + TX3 - TEM1
            QA(3)   = -TEM2
!
            WA(1)   = TX4
            WA(2)   = WA(2) + TX5 - ST1
            WA(3)   = -ST2
!
            TX2     = TEM1
            TX3     = TEM2
            TX4     = ST1
            TX5     = ST2
!
            VT(1)   = VT(2)
            TRW(1)  = TRW(2)
            VRW(1)  = VRW(2)
!
            IF (WVL(KTEM) .EQ. WCMIN) WA(1) = 0.0
            IF (WVL(LL)   .EQ. WCMIN) WA(2) = 0.0
            IF (WVL(L)    .EQ. WCMIN) WA(3) = 0.0
            DO N=KTEM,KBL
              AA(LL,N) = (WA(1)*QW(KTEM,N) * STLT(KTEM)                 &
     &                 +  WA(2)*QW(LL,N)   * STLT(LL)                   &
     &                 +  WA(3)*QW(L,N)    * STLT(L) ) * 0.5
            ENDDO
            AA(LL,KTEM) = AA(LL,KTEM) + QA(1)
            AA(LL,LL)   = AA(LL,LL)   + QA(2)
            AA(LL,L)    = AA(LL,L)    + QA(3)
            BUD(LL)     = (TX8 + RNN(LL)) * 0.5                         &
     &                    - RNB + TX6 - BUD(LL)
            AA(LL,KBL+1) = BUD(LL)
            RNB = TX6
            TX1 = 1.0
            TX8 = RNN(LL)
          ENDDO
          L  = KBL
          LL = L - 1
!         VT(2)   = GMS(L) * QRP(L)**0.1364
          VT(2)   = GMS(L) * QRPF(QRP(L))
          TRW(2)  = ETA(L) * QRP(L) * STLT(L)
          VRW(2)  = F3*WVL(L) - CTL2*VT(2)
          ST1     = STLA * TRW(2) * VT(2) * QRB(LL)
          BUD(L)  = ST1

          QA(2)   = DOF
          WA(2)   = DOFW
          DOF     = 1.1364 * BUD(L) * QRPI(L)
          DOFW    = -BUD(L) * STLT(L)
!
          RNF(LL) = RNF(LL) + ST1
!
          TEM3    = VRW(1) + VRW(2)
          TEM4    = TRW(1) + TRW(2)
!
          TX6     = .25 * TEM3 * TEM4
          TEM4    = TEM4 * CTL3
!
!-----BY QR ABOVE
!
          TEM1    = .25*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
          ST1     = .25*(TRW(1)*(CTL2*VT(1)-VRW(2))                     &
     &                * STLT(LL) + F3*TRW(2))
!-----BY QR BELOW
          TEM2    = .25*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
          ST2     = .25*(TRW(2)*(CTL2*VT(2)-VRW(1))                     &
     &                 * STLT(L)  + F3*TRW(1))
!
!      For the layer next to the top of the boundary layer
!
          QA(1)   = TX2
          QA(2)   = QA(2) + TX3 - TEM1
          QA(3)   = -TEM2
!
          WA(1)   = TX4
          WA(2)   = WA(2) + TX5 - ST1
          WA(3)   = -ST2
!
          TX2     = TEM1
          TX3     = TEM2
          TX4     = ST1
          TX5     = ST2
!
          IDW     = MAX(L-2, KD)
!
          IF (WVL(IDW) .EQ. WCMIN) WA(1) = 0.0
          IF (WVL(LL)  .EQ. WCMIN) WA(2) = 0.0
          IF (WVL(L)   .EQ. WCMIN) WA(3) = 0.0
!
          KK = IDW
          DO N=KK,L
            AA(LL,N) = (WA(1)*QW(KK,N) * STLT(KK)                       &
     &               +  WA(2)*QW(LL,N) * STLT(LL)                       &
     &               +  WA(3)*QW(L,N)  * STLT(L) ) * 0.5

          ENDDO
!
          AA(LL,IDW) = AA(LL,IDW) + QA(1)
          AA(LL,LL)  = AA(LL,LL)  + QA(2)
          AA(LL,L)   = AA(LL,L)   + QA(3)
          BUD(LL)    = (TX8+RNN(LL)) * 0.5 - RNB + TX6 - BUD(LL)
!
          AA(LL,L+1) = BUD(LL)
!
          RNB        = TRW(2) * VRW(2)
!
!      For the top of the boundary layer
!
          IF (RNB .LT. 0.0) THEN
             KK    = KBL
             TEM   = VT(2) * TRW(2)
             QA(2) = (RNB - CTL3*TEM) * QRPI(KK)
             WA(2) = CTL2 * TEM * STLT(KK)
          ELSE
             RNB   = 0.0
             QA(2) = 0.0
             WA(2) = 0.0
          ENDIF
!
          QA(1) = TX2
          QA(2) = DOF + TX3 - QA(2)
          QA(3) = 0.0
!
          WA(1) = TX4
          WA(2) = DOFW + TX5 - WA(2)
          WA(3) = 0.0
!
          KK = KBL
          IF (WVL(KK-1) .EQ. WCMIN) WA(1) = 0.0
          IF (WVL(KK)   .EQ. WCMIN) WA(2) = 0.0
!
          DO II=1,2
             N = KK + II - 2
             AA(KK,N) = (WA(1)*QW(KK-1,N) * STLT(KK-1)                  &
     &                +  WA(2)*QW(KK,N)   * STLT(KK)) * 0.5
          ENDDO
          FAC = 0.5
          LL  = KBL
          L   = LL + 1
          LM1 = LL - 1
          AA(LL,LM1)  = AA(LL,LM1) + QA(1)
          AA(LL,LL)   = AA(LL,LL)  + QA(2)
          BUD(LL)     = 0.5*RNN(LM1) - TX6 + RNB - BUD(LL)
          AA(LL,LL+1) = BUD(LL)
!
!-----SOLVING THE BUDGET EQUATIONS FOR DQR
!
          DO L=KD1,KBL
            LM1  = L - 1
            UNSAT = ABS(AA(LM1,LM1)) .LT. ABS(AA(L,LM1))
            DO  N=LM1,KBL+1
               IF (UNSAT) THEN
                  TX1       = AA(LM1,N)
                  AA(LM1,N) = AA(L,N)
                  AA(L,N)   = TX1
               ENDIF
            ENDDO
            TX1 = AA(L,LM1) / AA(LM1,LM1)
            DO  N=L,KBL+1
               AA(L,N) = AA(L,N) - TX1 * AA(LM1,N)
            ENDDO
          ENDDO     
!
!-----BACK SUBSTITUTION AND CHECK IF THE SOLUTION CONVERGES
!
          KK = KBL
          KK1 = KK + 1
          AA(KK,KK1) = AA(KK,KK1) / AA(KK,KK)      !   Qr correction !
          TX2        = ABS(AA(KK,KK1)) * QRPI(KK)  !   Error Measure !
!     if (lprnt) print *,' tx2a=',tx2,' aa1=',aa(kk,kk1)
!    &,' qrpi=',qrpi(kk)
!
          KK = KBL + 1
          DO L=KB1,KD,-1
             LP1   = L + 1
             TX1  = 0.0
             DO N=LP1,KBL
               TX1  = TX1 + AA(L,N) * AA(N,KK)
             ENDDO
             AA(L,KK) = (AA(L,KK) - TX1) / AA(L,L)       ! Qr correction !
             TX2      = MAX(TX2, ABS(AA(L,KK))*QRPI(L))  ! Error Measure !
!     if (lprnt) print *,' tx2b=',tx2,' aa1=',aa(l,kk)
!    &,' qrpi=',qrpi(l),' L=',L
          ENDDO
!
!         tem = 0.5
          if (tx2 .gt. 1.0 .and. abs(errq-tx2) .gt. 0.1) then
            tem = 0.5
!!        elseif (tx2 .lt. 0.1) then
!!          tem = 1.2
          else
            tem = 1.0
          endif
!
          DO L=KD,KBL
!            QRP(L) = MAX(QRP(L)+AA(L,KBL+1), QRMIN)
             QRP(L) = MAX(QRP(L)+AA(L,KBL+1)*tem, QRMIN)
          ENDDO
!
!     if (lprnt) print *,' itr=',itr,' tx2=',tx2
          IF (ITR .LT. ITRMIN) THEN
             TEM = ABS(ERRQ-TX2) 
             IF (TEM .GE. ERRMI2 .AND. TX2 .GE. ERRMIN) THEN 
               ERRQ  = TX2                              ! Further iteration !
             ELSE 
               SKPUP = .TRUE.                           ! Converges      !
               ERRQ  = 0.0                              ! Rain profile exists!
!     print *,' here1',' tem=',tem,' tx2=',tx2,' errmi2=',
!    *errmi2,' errmin=',errmin
             ENDIF 
          ELSE
             TEM = ERRQ - TX2
!            IF (TEM .LT. ZERO .AND. ERRQ .GT. 0.1) THEN
             IF (TEM .LT. ZERO .AND. ERRQ .GT. 0.5) THEN
!            IF (TEM .LT. ZERO .and.                                    &
!    &          (ntla .lt. numtla .or. ERRQ .gt. 0.5)) THEN
!     if (lprnt) print *,' tx2=',tx2,' errq=',errq,' tem=',tem
               SKPUP = .TRUE.                           ! No convergence !
               ERRQ = 10.0                              ! No rain profile!
!!!!         ELSEIF (ABS(TEM).LT.ERRMI2 .OR. TX2.LT.ERRMIN) THEN
             ELSEIF (TX2.LT.ERRMIN) THEN
               SKPUP = .TRUE.                           ! Converges      !
               ERRQ = 0.0                               ! Rain profile exists!
!     print *,' here2'
             elseif (tem .lt. zero .and. errq .lt. 0.1) then
               skpup = .true.
!              if (ntla .eq. numtla .or. tem .gt. -0.003) then
                 errq  = 0.0
!              else
!                errq = 10.0
!              endif
             ELSE
               ERRQ = TX2                               ! Further iteration !
!     if (lprnt) print *,' itr=',itr,' errq=',errq
!              if (itr .eq. itrmu .and. ERRQ .GT. ERRMIN*10             &
!    &            .and. ntla .eq. 1) ERRQ = 10.0 
             ENDIF
          ENDIF
!
!         if (lprnt) print *,' ERRQ=',ERRQ

        ENDIF                                           ! SKPUP  ENDIF!
!
      ENDDO                                          ! End of the ITR Loop!!
!     enddo                                          ! End of ntla loop
!
!     if(lprnt) then
!       print *,' QRP=',(QRP(L),L=KD,KBL)
!       print *,'RNF=',(RNF(L),L=KD,KBL),' RNT=',RNT,' RNB=',RNB
!    &,' errq=',errq
!     endif
!
      IF (ERRQ .LT. 0.1) THEN
        DDFT = .TRUE.
        RNB  = - RNB
   !    do l=kd1,kb1-1
   !      if (wvl(l)-wcbase .lt. 1.0E-9) ddft = .false.
   !    enddo
      ELSE
        DDFT = .FALSE.
      ENDIF
!
!     Caution !! Below is an adjustment to rain flux to maintain
!                conservation of precip!
!
      IF (DDFT) THEN
        TX1 = 0.0
        DO L=KD,KB1
          TX1 = TX1 + RNF(L)
        ENDDO
!     if (lprnt) print *,' tx1+rnt+rnb=',tx1+rnt+rnb, ' train=',train
        TX1 = TRAIN / (TX1+RNT+RNB)
        IF (ABS(TX1-1.0) .LT. 0.2) THEN
           RNT = MAX(RNT*TX1,ZERO)
           RNB = RNB * TX1
        ELSE
           DDFT = .FALSE.
           ERRQ = 10.0
        ENDIF
      ENDIF
      enddo                                          ! End of ntla loop
!
      DOF = 0.0
      IF (.NOT. DDFT) RETURN     ! Rain profile did not converge!
!

      DO L=KD,KB1
         RNF(L) = RNF(L) * TX1

      ENDDO
!     if (lprnt) print *,' TRAIN=',TRAIN
!     if (lprnt) print *,' RNF=',RNF
!
!     Adjustment is over
!
!     Downdraft
!
      DO L=KD,K
        WCB(L) = 0.0
      ENDDO
!
      SKPDD = .NOT. DDFT
!
      ERRQ  = 10.0
      IF (.NOT. SKPDD) THEN
!
!     Calculate Downdraft Properties
!

        KK = MAX(KB1,KD1)
        DO L=KK,K
          STLT(L) = STLT(L-1)
        ENDDO
        TEM1 = 1.0 / BB1
!
        DO L=KD,K
          IF (L .LE. KBL) THEN
            TEM     = STLA * TEM1
            STLT(L) = ETA(L) * STLT(L) * TEM / ROR(L)
          ELSE
            STLT(L) = 0.0
          ENDIF
        ENDDO
!       if (lprnt) print *,' STLT=',stlt

        rsum1 = 0.0
        rsum2 = 0.0

!
        IDN      = 99
        DO L=KD,K+1
          ETD(L)  = 0.0
          WVL(L)  = 0.0
!         QRP(L)  = 0.0
        ENDDO
        DO L=KD,K
          EVP(L)   = 0.0
          BUY(L)   = 0.0
          QRP(L+1) = 0.0
        ENDDO
        HOD(KD)  = HOL(KD)
        QOD(KD)  = QOL(KD)
        TX1      = 0.0                               ! sigma at the top
!!!     TX1      = STLT(KD)*QRB(KD)*ONE              ! sigma at the top
!       TX1      = MIN(STLT(KD)*QRB(KD)*ONE, ONE)    ! sigma at the top
!       TX1      = MIN(STLT(KD)*QRB(KD)*0.5, ONE)    ! sigma at the top
        RNTP     = 0.0
        TX5      = TX1
        QA(1)    = 0.0
!     if(lprnt) print *,' stlt=',stlt(kd),' qrb=',qrb(kd)
!    *,' tx1=',tx1,' ror=',ror(kd),' gms=',gms(kd),' rpart=',rpart
!    *,' rnt=',rnt
!
!       Here we assume RPART of detrained rain RNT goes to Pd
!
        IF (RNT .GT. 0.0) THEN
          if (TX1 .gt. 0.0) THEN
            QRP(KD) = (RPART*RNT / (ROR(KD)*TX1*GMS(KD)))               &
     &                                          ** (1.0/1.1364)
           else
             tx1 = RPART*RNT / (ROR(KD)*GMS(KD)*QRP(KD)**1.1364)
           endif
            RNTP    = (1.0 - RPART) * RNT
            BUY(KD) = - ROR(KD) * TX1 * QRP(KD)
        ELSE
          QRP(KD) = 0.0
        ENDIF
!
!     L-loop for the downdraft iteration from KD1 to K+1 (bottom surface)
!
!     BUD(KD) = ROR(KD)
      idnm = 1
      DO L=KD1,K+1

          QA(1) = 0.0
          ddlgk = idn(idnm) .eq. 99
          if (.not. ddlgk) cycle
          IF (L .LE. K) THEN
            ST1   = 1.0 - ALFIND(L)
            WA(1) = ALFIND(L)*HOL(L-1) + ST1*HOL(L)
            WA(2) = ALFIND(L)*QOL(L-1) + ST1*QOL(L)
            WA(3) = ALFIND(L)*TOL(L-1) + ST1*TOL(L)
            QA(2) = ALFIND(L)*HST(L-1) + ST1*HST(L)
            QA(3) = ALFIND(L)*QST(L-1) + ST1*QST(L)
          ELSE
            WA(1) = HOL(K)
            WA(2) = QOL(K)
            WA(3) = TOL(K)
            QA(2) = HST(K)
            QA(3) = QST(K)
          ENDIF
!
          FAC = 2.0
          IF (L .EQ. KD1) FAC = 1.0

          FACG    = FAC * 0.5 * GMF5     !  12/17/97
!
!         DDLGK   =  IDN(idnm) .EQ. 99
          BUD(KD) = ROR(L)

!         IF (DDLGK) THEN
            TX1    = TX5
            WVL(L) = MAX(WVL(L-1),ONE_M1)

            QRP(L) = MAX(QRP(L-1),QRP(L))
!
!           VT(1)  = GMS(L-1) * QRP(L-1) ** 0.1364
            VT(1)  = GMS(L-1) * QRPF(QRP(L-1))
            RNT    = ROR(L-1) * (WVL(L-1)+VT(1))*QRP(L-1)
!     if(lprnt) print *,' l=',l,' qa=',qa(1), ' tx1RNT=',RNT*tx1,
!    *' wvl=',wvl(l-1)
!    *,' qrp=',qrp(l-1),' tx5=',tx5,' tx1=',tx1,' rnt=',rnt

!

!           TEM    = MAX(ALM, 2.5E-4) * MAX(ETA(L), 1.0)
            TEM    = MAX(ALM,ONE_M6) * MAX(ETA(L), ONE)
!           TEM    = MAX(ALM, 1.0E-5) * MAX(ETA(L), 1.0)
            TRW(1) = PICON*TEM*(QRB(L-1)+QRT(L-1))
            TRW(2) = 1.0 / TRW(1)
!
            VRW(1) = 0.5 * (GAM(L-1) + GAM(L))
            VRW(2) = 1.0 / (VRW(1) + VRW(1))
!
            TX4    =  (QRT(L-1)+QRB(L-1))*(ONEBG*FAC*500.00*EKNOB)
!
            DOFW   = 1.0 / (WA(3) * (1.0 + NU*WA(2)))      !  1.0 / TVbar!
!
            ETD(L) = ETD(L-1)
            HOD(L) = HOD(L-1)
            QOD(L) = QOD(L-1)
!
            ERRQ   = 10.0

!
            IF (L .LE. KBL) THEN
              TX3 = STLT(L-1) * QRT(L-1) * (0.5*FAC)
              TX8 = STLT(L)   * QRB(L-1) * (0.5*FAC)
              TX9 = TX8 + TX3
            ELSE
              TX3 = 0.0
              TX8 = 0.0
              TX9 = 0.0
            ENDIF
!
            TEM  = WVL(L-1) + VT(1)
            IF (TEM .GT. 0.0) THEN
              TEM1 = 1.0 / (TEM*ROR(L-1))
              TX3 = VT(1) * TEM1 * ROR(L-1) * TX3
              TX6 = TX1 * TEM1
            ELSE
              TX6 = 1.0
            ENDIF
!         ENDIF
!
          IF (L .EQ. KD1) THEN
            IF (RNT .GT. 0.0) THEN
              TEM    = MAX(QRP(L-1),QRP(L))
              WVL(L) = TX1 * TEM * QRB(L-1)*(FACG*5.0)
            ENDIF
            WVL(L) = MAX(ONE_M2, WVL(L))
            TRW(1) = TRW(1) * 0.5
            TRW(2) = TRW(2) + TRW(2)
          ELSE
            IF (DDLGK) EVP(L-1) = EVP(L-2)
          ENDIF
!
!       No downdraft above level IDH
!

          IF (L .LT. IDH) THEN

            ETD(L)   = 0.0
            HOD(L)   = WA(1)
            QOD(L)   = WA(2)
            EVP(L-1) = 0.0
            WVL(L)   = 0.0
            QRP(L)   = 0.0
            BUY(L)   = 0.0
            TX5      = TX9
            ERRQ     = 0.0
            RNTP     = RNTP + RNT * TX1
            RNT      = 0.0
            WCB(L-1) = 0.0
          ENDIF
!         BUD(KD) = ROR(L)
!
!       Iteration loop for a given level L begins
!
!         if (lprnt) print *,' tx8=',tx8,' tx9=',tx9,' tx5=',tx5
!    &,                      ' tx1=',tx1
          DO ITR=1,ITRMD
!
!           UNSAT =  DDLGK .AND. (ERRQ .GT. ERRMIN)
            UNSAT =  ERRQ .GT. ERRMIN
            IF (UNSAT) THEN
!
!             VT(1)  = GMS(L) * QRP(L) ** 0.1364
              VT(1)  = GMS(L) * QRPF(QRP(L))
              TEM    =  WVL(L) + VT(1)
!
              IF (TEM .GT. 0.0) THEN
                ST1    = ROR(L) * TEM * QRP(L) + RNT
                IF (ST1 .NE. 0.0) ST1 = 2.0 * EVP(L-1) / ST1
                TEM1   = 1.0 / (TEM*ROR(L))
                TEM2   = VT(1) * TEM1 * ROR(L) * TX8
              ELSE
                TEM1   = 0.0
                TEM2   = TX8
                ST1    = 0.0
              ENDIF
!     if (lprnt) print *,' st1=',st1,' tem=',tem,' ror=',ror(l)
!    &,' qrp=',qrp(l),' rnt=',rnt,' ror1=',ror(l-1),' wvl=',wvl(l)
!    &,' wvl1=',wvl(l-1),' tem2=',tem2,' vt=',vt(1),' tx3=',tx3
!
              st2 = tx5
              TEM = ROR(L)*WVL(L) - ROR(L-1)*WVL(L-1)
              if (tem .gt. 0.0) then
                TX5 = (TX1 - ST1 + TEM2 + TX3)/(1.0+tem*tem1)
              else
                TX5 = TX1 - tem*tx6 - ST1 + TEM2 + TX3
              endif
              TX5   = MAX(TX5,ZERO)
              tx5 = 0.5 * (tx5 + st2)
!
!             qqq = 1.0 + tem * tem1 * (1.0 - sialf)
!
!             if (qqq .gt. 0.0) then
!               TX5   = (TX1 - sialf*tem*tx6 - ST1 + TEM2 + TX3) / qqq
!             else
!               TX5   = (TX1 - tem*tx6 - ST1 + TEM2 + TX3)
!             endif
!
!     if(lprnt) print *,' tx51=',tx5,' tx1=',tx1,' st1=',st1,' tem2='
!     if(tx5 .le. 0.0 .and. l .gt. kd+2)
!    * print *,' tx51=',tx5,' tx1=',tx1,' st1=',st1,' tem2='
!    *,tem2,' tx3=',tx3,' tem=',tem,' tem1=',tem1,' wvl=',wvl(l-1),
!    &wvl(l),' l=',l,' itr=',itr,' evp=',evp(l-1),' vt=',vt(1)
!    *,' qrp=',qrp(l),' rnt=',rnt,' kd=',kd
!     if (lprnt) print *,' etd=',etd(l),' wvl=',wvl(l)
!    &,' trw=',trw(1),trw(2),' ror=',ror(l),' wa=',wa


!
              TEM1   = ETD(L)
              ETD(L) = ROR(L) * TX5 * MAX(WVL(L),ZERO)
!
              if (etd(l) .gt. 0.0) etd(l) = 0.5 * (etd(l) + tem1)
!

              DEL_ETA = ETD(L) - ETD(L-1)

!               TEM       = DEL_ETA * TRW(2)
!               TEM2      = MAX(MIN(TEM, 1.0), -1.0)
!               IF (ABS(TEM) .GT. 1.0 .AND. ETD(L) .GT. 0.0 ) THEN
!                 DEL_ETA = TEM2 * TRW(1)
!                 ETD(L)  = ETD(L-1) + DEL_ETA
!               ENDIF
!               IF (WVL(L) .GT. 0.0) TX5 = ETD(L) / (ROR(L)*WVL(L))
!
                ERRE  = ETD(L) - TEM1
!
                tem  = max(abs(del_eta), trw(1))
                tem2 = del_eta / tem
                TEM1 = SQRT(MAX((tem+DEL_ETA)*(tem-DEL_ETA),ZERO))
!               TEM1 = SQRT(MAX((TRW(1)+DEL_ETA)*(TRW(1)-DEL_ETA),0.0))

                EDZ  = (0.5 + ASIN(TEM2)*PIINV)*DEL_ETA + TEM1*PIINV

              DDZ   = EDZ - DEL_ETA
              WCB(L-1) = ETD(L) + DDZ
!
              TEM1  = HOD(L)
              IF (DEL_ETA .GT. 0.0) THEN
                QQQ    = 1.0 / (ETD(L) + DDZ)
                HOD(L) = (ETD(L-1)*HOD(L-1) + DEL_ETA*HOL(L-1)          &
     &                                            + DDZ*WA(1)) * QQQ
                QOD(L) = (ETD(L-1)*QOD(L-1) + DEL_ETA*QOL(L-1)          &
     &                                            + DDZ*WA(2)) * QQQ
              ELSEif((ETD(L-1) + EDZ) .gt. 0.0) then
                QQQ    = 1.0 / (ETD(L-1) + EDZ)
                HOD(L) = (ETD(L-1)*HOD(L-1) + EDZ*WA(1)) * QQQ
                QOD(L) = (ETD(L-1)*QOD(L-1) + EDZ*WA(2)) * QQQ
              ENDIF
              ERRH  = HOD(L) - TEM1
              ERRQ  = ABS(ERRH/HOD(L))  + ABS(ERRE/MAX(ETD(L),ONE_M5))
!     if (lprnt) print *,' ERRQP=',errq,' errh=',errh,' hod=',hod(l)
!    &,' erre=',erre,' etd=',etd(l),' del_eta=',del_eta
              DOF   = DDZ
              VT(2) = QQQ

!
              DDZ  = DOF
              TEM4 = QOD(L)
              TEM1 = VRW(1)
!
              QHS  = QA(3) + 0.5 * (GAF(L-1)+GAF(L))                    &
     &                           * (HOD(L)-QA(2))
!
!                                           First iteration       !
!
              ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
              TEM2 = ROR(L) * QRP(L)
              CALL QRABF(TEM2,QRAF,QRBF)
              TEM6 = TX5 * (1.6 + 124.9 * QRAF) * QRBF * TX4
!
              CE   = TEM6 * ST2 / ((5.4E5*ST2 + 2.55E6)*(ETD(L)+DDZ))
!
              TEM2   = - ((1.0+TEM1)*(QHS+CE) + TEM1*QOD(L))
              TEM3   = (1.0 + TEM1) * QHS * (QOD(L)+CE)
              TEM    = MAX(TEM2*TEM2 - 4.0*TEM1*TEM3,ZERO)
              QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
!

!
!                                            second iteration   !
!
              ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
              CE   = TEM6 * ST2 / ((5.4E5*ST2 + 2.55E6)*(ETD(L)+DDZ))
!             CEE  = CE * (ETD(L)+DDZ)
!


              TEM2   = - ((1.0+TEM1)*(QHS+CE) + TEM1*tem4)
              TEM3   = (1.0 + TEM1) * QHS * (tem4+CE)
              TEM    = MAX(TEM2*TEM2 - 4.0*TEM1*TEM3,ZERO)
              QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
!                                              Evaporation in Layer L-1
!

              EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)
!                                              Calculate Pd (L+1/2)
              QA(1)    = TX1*RNT + RNF(L-1) - EVP(L-1)
!
!     if(lprnt) print *,' etd=',etd(l),' tx5=',tx5,' rnt=',rnt
!    *,' rnf=',rnf(l-1),' evp=',evp(l-1),' itr=',itr,' L=',L

!
              if (qa(1) .gt. 0.0) then
              IF (ETD(L) .GT. 0.0) THEN
                TEM    = QA(1) / (ETD(L)+ROR(L)*TX5*VT(1))
                QRP(L) = MAX(TEM,ZERO)
              ELSEIF (TX5 .GT. 0.0) THEN
                QRP(L) = (MAX(ZERO,QA(1)/(ROR(L)*TX5*GMS(L))))           &
     &                                          ** (1.0/1.1364)
              ELSE
                QRP(L) = 0.0
              ENDIF
              else
                qrp(l) = 0.5 * qrp(l)
              endif
!                                              Compute Buoyancy
              TEM1   = WA(3)+(HOD(L)-WA(1)-ALHL*(QOD(L)-WA(2)))         &
     &                                                  * (1.0/CP)
!             if (lprnt) print *,' tem1=',tem1,' wa3=',wa(3),' hod='
!    &,hod(l),' wa1=',wa(1),' qod=',qod(l),' wa2=',wa(2),' alhl=',alhl
!    &,' cmpor=',cmpor,' dofw=',dofw,' prl=',prl(l),' qrp=',qrp(l)
              TEM1   = TEM1 * (1.0 + NU*QOD(L))
              ROR(L) = CMPOR * PRL(L) / TEM1
              TEM1   = TEM1 * DOFW
!!!           TEM1   = TEM1 * (1.0 + NU*QOD(L)) * DOFW

              BUY(L) = (TEM1 - 1.0 - QRP(L)) * ROR(L) * TX5
!                                              Compute W (L+1/2)

              TEM1   = WVL(L)
!             IF (ETD(L) .GT. 0.0) THEN
              WVL(L) = VT(2) * (ETD(L-1)*WVL(L-1) - FACG                &
     &                 * (BUY(L-1)*QRT(L-1)+BUY(L)*QRB(L-1)))
!
!             if (lprnt) print *,' wvl=',wvl(l),'vt2=',vt(2),' buy1='
!    &,buy(l-1),' buy=',buy(l),' qrt1=',qrt(l-1),' qrb1=',qrb(l-1)
!    &,' etd1=',etd(l-1),' wvl1=',wvl(l-1)
!             ENDIF
!
              if (wvl(l) .lt. 0.0) then
!               WVL(L) = max(wvl(l), 0.1*tem1)
!               WVL(L) = 0.5*tem1
!               WVL(L) = 0.1*tem1
!               WVL(L) = 0.0
                WVL(L) = 1.0e-10
              else
                WVL(L) = 0.5*(WVL(L)+TEM1)
              endif

!
!             WVL(L) = max(0.5*(WVL(L)+TEM1), 0.0)

              ERRW   = WVL(L) - TEM1
!
              ERRQ   = ERRQ + ABS(ERRW/MAX(WVL(L),ONE_M5))

!     if (lprnt) print *,' errw=',errw,' wvl=',wvl(l)
!     if(lprnt .or. tx5 .eq. 0.0) then
!     if(tx5 .eq. 0.0 .and. l .gt. kbl) then
!        print *,' errq=',errq,' itr=',itr,' l=',l,' wvl=',wvl(l)
!    &,' tx5=',tx5,' idnm=',idnm,' etd1=',etd(l-1),' etd=',etd(l)
!    &,' kbl=',kbl
!     endif
!
!     if(lprnt) print *,' itr=',itr,' itrmnd=',itrmnd,' itrmd=',itrmd
!             IF (ITR .GE. MIN(ITRMIN,ITRMD/2)) THEN
              IF (ITR .GE. MIN(ITRMND,ITRMD/2)) THEN
!     if(lprnt) print *,' itr=',itr,' etd1=',etd(l-1),' errq=',errq
                IF (ETD(L-1) .EQ. 0.0 .AND. ERRQ .GT. 0.2) THEN
!     if(lprnt) print *,' bud=',bud(kd),' wa=',wa(1),wa(2)
                  ROR(L)   = BUD(KD)
                  ETD(L)   = 0.0
                  WVL(L)   = 0.0
                  ERRQ     = 0.0
                  HOD(L)   = WA(1)
                  QOD(L)   = WA(2)
!                 TX5      = TX1 + TX9
                  if (L .le. KBL) then
                    TX5      = TX9
                  else
                    TX5 = (STLT(KB1) * QRT(KB1)                         &
     &                  +  STLT(KBL) * QRB(KB1)) * (0.5*FAC)
                  endif

!     if(lprnt) print *,' tx1=',tx1,' rnt=',rnt,' rnf=',rnf(l-1)
!    *,' evp=',evp(l-1),' l=',l
                  EVP(L-1) = 0.0
                  TEM      = MAX(TX1*RNT+RNF(L-1),ZERO)
                  QA(1)    = TEM - EVP(L-1)
!                 IF (QA(1) .GT. 0.0) THEN
!     if(lprnt) print *,' ror=',ror(l),' tx5=',tx5,' tx1=',tx1
!    *,' tx9=',tx9,' gms=',gms(l),' qa=',qa(1)
!     if(lprnt) call mpi_quit(13)
!     if (tx5 .eq. 0.0 .or. gms(l) .eq. 0.0)
!     if (lprnt) 
!    *  print *,' Atx5=',tx5,' gms=',gms(l),' ror=',ror(l)
!    *,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9
!    *,' kbl=',kbl,' etd1=',etd(l-1),' idnm=',idnm,' idn=',idn(idnm)
!    *,' errq=',errq
                  QRP(L)   = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
     &                                            ** (1.0/1.1364)
!                 endif
                  BUY(L)   = - ROR(L) * TX5 * QRP(L)
                  WCB(L-1) = 0.0
                ENDIF
!
                DEL_ETA = ETD(L) - ETD(L-1)
                IF(DEL_ETA .LT. 0.0 .AND. ERRQ .GT. 0.1) THEN
                  ROR(L)   = BUD(KD)
                  ETD(L)   = 0.0
                  WVL(L)   = 0.0
!!!!!             TX5      = TX1 + TX9
                  CLDFRD(L-1) = TX5
!
                  DEL_ETA  = - ETD(L-1)
                  EDZ      = 0.0
                  DDZ      = -DEL_ETA
                  WCB(L-1) = DDZ

!
                  HOD(L)   = HOD(L-1)
                  QOD(L)   = QOD(L-1)

!
                  TEM4     = QOD(L)
                  TEM1     = VRW(1)
!
                  QHS      = QA(3) + 0.5 * (GAF(L-1)+GAF(L))            &
     &                                   * (HOD(L)-QA(2))

!
!                                           First iteration       !
!
                  ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                  TEM2 = ROR(L) * QRP(L-1)
                  CALL QRABF(TEM2,QRAF,QRBF)
                  TEM6 = TX5 * (1.6 + 124.9 * QRAF) * QRBF * TX4
!
                  CE   = TEM6*ST2/((5.4E5*ST2 + 2.55E6)*(ETD(L)+DDZ))
!

                  TEM2   = - ((1.0+TEM1)*(QHS+CE) + TEM1*QOD(L))
                  TEM3   = (1.0 + TEM1) * QHS * (QOD(L)+CE)
                  TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                  QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
!
!                                            second iteration   !
!
                  ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                  CE   = TEM6*ST2/((5.4E5*ST2 + 2.55E6)*(ETD(L)+DDZ))
!                 CEE  = CE * (ETD(L)+DDZ)
!


                  TEM2   = - ((1.0+TEM1)*(QHS+CE) + TEM1*tem4)
                  TEM3   = (1.0 + TEM1) * QHS * (tem4+CE)
                  TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                  QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))

!                                              Evaporation in Layer L-1
!
                  EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)

!                                               Calculate Pd (L+1/2)
!                 RNN(L-1) = TX1*RNT + RNF(L-1) - EVP(L-1)
                  QA(1)    = TX1*RNT + RNF(L-1)
                  EVP(L-1) = min(EVP(L-1), QA(1))
                  QA(1)    = QA(1) - EVP(L-1)
                  qrp(l)   = 0.0
!
!     if (tx5 .eq. 0.0 .or. gms(l) .eq. 0.0)
!     if (lprnt)
!    *  print *,' Btx5=',tx5,' gms=',gms(l),' ror=',ror(l)
!    *,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9
!    *,' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
!    &,' evp=',evp(l-1)
!
!                 IF (QA(1) .GT. 0.0) THEN
!!                  RNS(L-1) = QA(1)
!!!                 tx5      = tx9
!                   QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
!    &                                         ** (1.0/1.1364)
!                 endif
!                 ERRQ   = 0.0
!                                              Compute Buoyancy
!                 TEM1   = WA(3)+(HOD(L)-WA(1)-ALHL*(QOD(L)-WA(2)))     &
!    &                                                  * (1.0/CP)
!                 TEM1   = TEM1 * (1.0 + NU*QOD(L)) * DOFW

!                 BUY(L) = (TEM1 - 1.0 - QRP(L)) * ROR(L) * TX5
!
!                 IF (QA(1) .GT. 0.0) RNS(L) = QA(1)
                  IF (L .LE. K) THEN
                     RNS(L) = QA(1)
                     QA(1)  = 0.0
                  ENDIF
                  tx5      = tx9
                  ERRQ     = 0.0
                  QRP(L)   = 0.0
                  BUY(L)   = 0.0
                                                                              
!
                ENDIF
              ENDIF
            ENDIF
!

          ENDDO                ! End of the iteration loop  for a given L!
!       if (kd .eq. 13 .and. .not. ddft) stop
          IF (L .LE. K) THEN
            IF (ETD(L-1) .EQ. 0.0                                       &
     &         .AND. ERRQ .GT. 0.1 .and. l .le. kbl) THEN
!!!  &         .AND. ERRQ .GT. ERRMIN*10.0 .and. l .le. kbl) THEN
!    &         .AND. ERRQ .GT. ERRMIN*10.0) THEN
               ROR(L)   = BUD(KD)
               HOD(L)   = WA(1)
               QOD(L)   = WA(2)
               TX5      =       TX9     ! Does not make too much difference!
!              TX5      = TX1 + TX9
               EVP(L-1) = 0.0
!              EVP(L-1) = CEE * (1.0 - qod(l)/qa(3))
               QA(1)    = TX1*RNT + RNF(L-1)
               EVP(L-1) = min(EVP(L-1), QA(1))
               QA(1)    = QA(1) - EVP(L-1)
!              QRP(L)   = 0.0
      if (tx5 .eq. 0.0 .or. gms(l) .eq. 0.0) then
        print *,' Ctx5=',tx5,' gms=',gms(l),' ror=',ror(l)              &
     &,' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9                     &
     &,' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
      endif
!              IF (QA(1) .GT. 0.0) THEN
                 QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))                 &
     &                                         ** (1.0/1.1364)
!              ENDIF
               ETD(L)   = 0.0
               WVL(L)   = 0.0
               ST1      = 1.0 - ALFIND(L)

               ERRQ     = 0.0
               BUY(L)   = - ROR(L) * TX5 * QRP(L)
               WCB(L-1) = 0.0
            ENDIF
          ENDIF

!
!
          LL = MIN(IDN(idnm), K+1)
          IF (ERRQ .LT. 1.0 .AND. L .LE. LL) THEN
            IF (ETD(L-1) .GT. 0.0 .AND. ETD(L) .EQ. 0.0) THEN
             IDN(idnm) = L
             wvl(l)    = 0.0
             if (L .lt. KBL .or. tx5 .gt. 0.0) idnm  = idnm + 1
             errq      = 0.0
            ENDIF
            if (etd(l) .eq. 0.0 .and. l .gt. kbl) then
              idn(idnm) = l
              if (tx5 .gt. 0.0) idnm  = idnm + 1
            endif
          ENDIF

!       if (lprnt) then
!       print *,' ERRQ=',ERRQ,' IDN=',IDN(idnm),' idnm=',idnm
!       print *,' L=',L,' QRP=',QRP(L),' ETD=',ETD(L),' QA=',QA(1)
!    *,' evp=',evp(l-1),' rnf=',rnf(l-1)
!       endif

! 
!     If downdraft properties are not obtainable, (i.e.solution does
!      not converge) , no downdraft is assumed
!
!          IF (ERRQ .GT. ERRMIN*100.0 .AND. IDN(idnm) .EQ. 99)          &
           IF (ERRQ .GT. 0.1 .AND. IDN(idnm) .EQ. 99)                   &
     &                          DDFT = .FALSE.
!
!
        DOF = 0.0
        IF (.NOT. DDFT) RETURN
!
!     if (ddlgk .or. l .le. idn(idnm)) then
!     rsum2 = rsum2 + evp(l-1)
!     print *,' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' qa=',qa(1)
!    *,' evp=',evp(l-1)
!     else
!     rsum1 = rsum1 + rnf(l-1)
!     print *,' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' rnf=',rnf(l-1)
!     endif

      ENDDO                      ! End of the L Loop of downdraft !

        TX1 = 0.0

        DOF = QA(1)
!
!     print *,' dof=',dof,' rntp=',rntp,' rnb=',rnb
!     print *,' total=',(rsum1+dof+rntp+rnb)

      ENDIF                       ! SKPDD endif
!

       RNN(KD) = RNTP
       TX1     = EVP(KD)
       TX2     = RNTP + RNB + DOF

!     if (lprnt) print *,' tx2=',tx2
       II = IDH
       IF (II .GE. KD1+1) THEN
          RNN(KD)   = RNN(KD) + RNF(KD)
          TX2       = TX2 + RNF(KD)
          RNN(II-1) = 0.0
          TX1       = EVP(II-1)
        ENDIF
!     if (lprnt) print *,' tx2=',tx2,' idnm=',idnm,' idn=',idn(idnm)
        DO L=KD,K
          II = IDH

          IF (L .GT. KD1 .AND. L .LT. II) THEN
            RNN(L-1) = RNF(L-1)
            TX2      = TX2 + RNN(L-1)

          ELSEIF (L .GE. II .AND. L .LT. IDN(idnm)) THEN
!!!       ELSEIF (L .GE. II .AND. L .LE. IDN(idnm)) THEN

!           do jj=2,idnm
!             if (l .ge. idn(jj-1) .and. l .lt. idn(jj)) then
!!!             RNN(L)   = 0.0
!!!             TX1      = TX1 + EVP(L)
!             endif
!           enddo
!
            rnn(l) = rns(l)
            tx2    = tx2 + rnn(l)
            TX1    = TX1 + EVP(L)

          ELSEIF (L .GE. IDN(idnm)) THEN
            ETD(L+1) = 0.0
            HOD(L+1) = 0.0
            QOD(L+1) = 0.0
            EVP(L)   = 0.0
            RNN(L)   = RNF(L) + RNS(L)
            TX2      = TX2    + RNN(L)
          ENDIF
!     if (lprnt) print *,' tx2=',tx2,' L=',L,' rnn=',rnn(l)
        ENDDO
!       IF (K+1 .GT. IDN(idnm)) THEN
!         ETD(K+1) = 0.0
!         HOD(K+1) = 0.0
!         QOD(K+1) = 0.0
!         EVP(K)   = 0.0
!         RNN(K)   = RNF(K)
!         TX2      = TX2 + RNN(K)
!       ENDIF
!
!      For Downdraft case the rain is that falls thru the bottom

        L = KBL

        RNN(L)    = RNN(L) + RNB
        CLDFRD(L) = TX5

!
!     Caution !! Below is an adjustment to rain flux to maintain
!                conservation of precip!

!
!     if (lprnt) print *,' train=',train,' tx2=',tx2,' tx1=',tx1

        IF (TX1 .GT. 0.0) THEN
          TX1 = (TRAIN - TX2) / TX1
        ELSE
          TX1 = 0.0
        ENDIF

!       TX5      = EVP(KBL)

!!      EVP(KBL) = EVP(KBL) * TX1

!       TX3      = RNN(KBL) + EVP(KBL) + DOF
!       TX2      = RNN(KBL)
!       TX4      = EVP(KBL)

!       DO L=KD,KB1
        DO L=KD,K

!         TX5    = TX5 + EVP(L)
          EVP(L) = EVP(L) * TX1
!         TX3    = TX3 + EVP(L) + RNN(L)
!         TX2    = TX2 + RNN(L)
!         TX4    = TX4 + EVP(L)
        ENDDO
!
!     if (lprnt .and. kd .eq. 52) stop
!***********************************************************************
!***********************************************************************

      RETURN
      END
      SUBROUTINE QSATCN(TT,P,Q,DQDT,lprnt)

      use module_ras
      real(kind=kind_phys) TT, P, Q, DQDT
      logical lprnt


      real TT4, P4, Q4, DQDT4

      TT4 =TT
      P4  =p
      

      DQDT4 = GEOS_DQsat(TT4,P4,qsat=q4)

      DQDT = DQDT4
      Q    = Q4
     

  !DQST3 = DQSAT  ( &
  !            TEMP   , &
  !            PP    , QSAT = QST3  )

      return
      end

      SUBROUTINE ANGRAD( PRES, ALM,  AL2, TLA, PRB, WFN, UFN)
!     SUBROUTINE ANGRAD( PRES, ALM, STLA, CTL2, AL2                     &
!    &,                  PI, TLA, PRB, WFN, UFN)
      use module_ras
      implicit none

!     real(kind=kind_phys) PRES, STLA, CTL2, pi,  pifac                 &
      real(kind=kind_phys) PRES                                         &
     &,                    ALM,  AL2,  TLA,  TEM, TEM1                  &
     &,                    PRB,  ACR,  WFN,  UFN
!
      integer i
!
!     pifac = pi / 180.0
!     print *,' pres=',pres
      IF (TLA .LT. 0.0) THEN
          IF (PRES .LE. PLAC(1)) THEN
            TLA = TLAC(1)
          ELSEIF (PRES .LE. PLAC(2)) THEN
            TLA = TLAC(2) + (PRES-PLAC(2))*tlbpl(1)
          ELSEIF (PRES .LE. PLAC(3)) THEN
            TLA = TLAC(3) + (PRES-PLAC(3))*tlbpl(2)
          ELSEIF (PRES .LE. PLAC(4)) THEN
            TLA = TLAC(4) + (PRES-PLAC(4))*tlbpl(3)
          ELSEIF (PRES .LE. PLAC(5)) THEN
            TLA = TLAC(5) + (PRES-PLAC(5))*tlbpl(4)
          ELSEIF (PRES .LE. PLAC(6)) THEN
            TLA = TLAC(6) + (PRES-PLAC(6))*tlbpl(5)
          ELSEIF (PRES .LE. PLAC(7)) THEN
            TLA = TLAC(7) + (PRES-PLAC(7))*tlbpl(6)
          ELSEIF (PRES .LE. PLAC(8)) THEN
            TLA = TLAC(8) + (PRES-PLAC(8))*tlbpl(7)
          ELSE
            TLA = TLAC(8)
          ENDIF
!         tla = tla * 1.5

!         STLA = SIN(TLA*PIFAC)
!         TEM1 = COS(TLA*PIFAC)
!         CTL2 = TEM1 * TEM1

      ELSE
!         STLA = SIN(TLA*PIFAC)
!         TEM1 = COS(TLA*PIFAC)
!         CTL2 = TEM1 * TEM1

      ENDIF
        IF (PRES .GE. REFP(1)) THEN
          TEM = REFR(1)
        ELSEIF (PRES .GE. REFP(2)) THEN
          TEM = REFR(1) + (PRES-REFP(1)) * drdp(1)
        ELSEIF (PRES .GE. REFP(3)) THEN
          TEM = REFR(2) + (PRES-REFP(2)) * drdp(2)
        ELSEIF (PRES .GE. REFP(4)) THEN
          TEM = REFR(3) + (PRES-REFP(3)) * drdp(3)
        ELSEIF (PRES .GE. REFP(5)) THEN
          TEM = REFR(4) + (PRES-REFP(4)) * drdp(4)
        ELSEIF (PRES .GE. REFP(6)) THEN
          TEM = REFR(5) + (PRES-REFP(5)) * drdp(5)
        ELSE
          TEM = REFR(6)
        ENDIF
!!      AL2 = min(ALMAX, MAX(ALM, 2.0E-4/TEM))
!       AL2 = min(2.0E-3, MAX(ALM, 2.0E-4/TEM))
!
        tem = 2.0E-4 / tem
        al2 = min(4.0*tem, max(alm, tem))
!
      RETURN
      END
      SUBROUTINE SETQRP
      use module_ras
      implicit none

      real(kind=kind_phys) tem2,tem1,x,xinc,xmax,xmin
      integer jx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!CFPP$ NOCONCUR R
!     XMIN=1.0E-6
      XMIN=0.0
      XMAX=5.0
      XINC=(XMAX-XMIN)/(NQRP-1)
      C1XQRP=1.-XMIN/XINC
      C2XQRP=1./XINC
      TEM1 = 0.001 ** 0.2046
      TEM2 = 0.001 ** 0.525
      DO JX=1,NQRP
        X         = XMIN + (JX-1)*XINC
        TBQRP(JX) =        X ** 0.1364
        TBQRA(JX) = TEM1 * X ** 0.2046
        TBQRB(JX) = TEM2 * X ** 0.525
      ENDDO    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      FUNCTION QRPF(QRP)
!
      use module_ras
      implicit none

      real(kind=kind_phys) QRP, QRPF, XJ, REAL_NQRP, ONE_2
      PARAMETER (ONE_2=1.)
      INTEGER JX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL_NQRP=REAL(NQRP)
      XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE_2),REAL_NQRP)
!     XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),FLOAT(NQRP))
      JX   = MIN(XJ,NQRP-ONE_2)
      QRPF = TBQRP(JX)  + (XJ-JX) * (TBQRP(JX+1)-TBQRP(JX))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE QRABF(QRP,QRAF,QRBF)
      use module_ras
      implicit none
!
      real(kind=kind_phys) QRP, QRAF, QRBF, XJ, REAL_NQRP, ONE_2
      PARAMETER (ONE_2=1.)
      INTEGER JX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL_NQRP=REAL(NQRP)
      XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE_2),REAL_NQRP)
      JX   = MIN(XJ,NQRP-ONE_2)
      XJ   = XJ - JX
      QRAF = TBQRA(JX)  + XJ * (TBQRA(JX+1)-TBQRA(JX))
      QRBF = TBQRB(JX)  + XJ * (TBQRB(JX+1)-TBQRB(JX))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      SUBROUTINE SETVTP
      use module_ras
      implicit none

      real(kind=kind_phys) vtpexp,xinc,x,xmax,xmin
      integer jx
      PARAMETER(VTPEXP=-0.3636)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!CFPP$ NOCONCUR R
      XMIN=0.05
      XMAX=1.5
      XINC=(XMAX-XMIN)/(NVTP-1)
      C1XVTP=1.-XMIN/XINC
      C2XVTP=1./XINC
      DO JX=1,NVTP
        X         = XMIN + (JX-1)*XINC
        TBVTP(JX) =        X ** VTPEXP
      ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      FUNCTION VTPF(ROR)
!
      use module_ras
      implicit none
      real(kind=kind_phys) ROR, VTPF, XJ, REAL_NVTP, ONE_2
      PARAMETER (ONE_2=1.)
      INTEGER JX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL_NVTP=REAL(NVTP)
      XJ   = MIN(MAX(C1XVTP+C2XVTP*ROR,ONE_2),REAL_NVTP)
      JX   = MIN(XJ,NVTP-ONE_2)
      VTPF = TBVTP(JX)  + (XJ-JX) * (TBVTP(JX+1)-TBVTP(JX))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
      FUNCTION CLF(PRATE)
!
      use MODULE_RAS
      implicit none
      real(kind=kind_phys) PRATE, CLF
!
      real (kind=kind_phys), parameter :: ccf1=0.30, ccf2=0.09          &
     &,                                   ccf3=0.04, ccf4=0.01          &
     &,                                   pr1=1.0,   pr2=5.0            &
     &,                                   pr3=20.0
!
      if (prate .lt. pr1) then
        clf = ccf1
      elseif (prate .lt. pr2) then
        clf = ccf2
      elseif (prate .lt. pr3) then
        clf = ccf3
      else
        clf = ccf4
      endif
!
      RETURN
      END
