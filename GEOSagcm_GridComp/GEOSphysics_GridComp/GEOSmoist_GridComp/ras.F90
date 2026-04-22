
! $Id$

MODULE RAS

   use ESMF
   use GEOS_Mod
   use GEOS_UtilsMod, only : DQSAT=>GEOS_DQsat
   use module_ras
   !use aer_cloud, only: AerProps, getINsubset
   !use aer_cloud
   use RASPARAMS


   IMPLICIT NONE

   PRIVATE
   PUBLIC RASE

   character(len=32), parameter, public :: RasCodes(-2:7) = (/ &
         "Default                         ", & ! -2
         "Invalid Code                    ", & ! -1
         "Cloud detraining here is active ", & ! 0
         "PBL h < layer's h*              ", & ! 1
         "No Valid lambda                 ", & ! 2
         "Lambda out of bounds            ", & ! 3
         "A < Acrit                       ", & ! 4
         "Negative Kernel                 ", & ! 5
         "Invalid Code                    ", & ! 6
         "RH Trigger not met              "  & ! 7
         /)

      integer, parameter :: nsmx_par = 20 !maximum number of modes allowed  
      type :: AerProps            
          sequence 
          real, dimension(nsmx_par)  :: num !Num conc m-3
          real, dimension(nsmx_par)  :: dpg !dry Geometric size, m
          real, dimension(nsmx_par)  :: sig  !logarithm (base e) of the dry geometric disp
	      real, dimension(nsmx_par)  :: den  !dry density , Kg m-3
  	      real, dimension(nsmx_par)  :: kap !Hygroscopicity parameter 
 	      real, dimension(nsmx_par)  :: fdust! mass fraction of dust 
	      real, dimension(nsmx_par)  :: fsoot ! mass fraction of soot
	      real, dimension(nsmx_par)  :: forg ! mass fraction of organics
	      integer   :: nmods  ! total number of modes (nmods<nmodmax)
      end type AerProps     
      

CONTAINS

   SUBROUTINE RASE(IDIM, IRUN, K0, ICMIN, DT ,                      &
         CPO,ALHLO,ALHL1,TICE,GRAVO,                      &
         SEEDRAS,IRAS,JRAS,SIGE,                          &
         KCBL,WGT0,WGT1,ZCBL,MXDIAM,TPERT,QPERT,          &
         THO, QHO, UHO, VHO,                              & 
         QSS, DQS, CNV_FRACTION, RASAL2_2d,               &
         CO_AUTO,                                         &
         pko, plo, phio, phie, qlo, qio,                  &
         PLE, PKE, CLW, FLX, FLXD, FLXC,                  &
         CNV_PRC3,                                        &
         CNV_UPDFRC,                                      &
         CNV_CVW,                                         &
         CNV_QC,                                          &
         ENTLAM,                                          &
         CLAN,                                            &
         HHO, HSO,PRECU,                                  &
         RASPARAMS,                                       &
         RAS_NO_NEG,                                      &
     !!  RAS Relaxation Diagnostics
         RAS_TIME, RAS_TRG, RAS_TOKI, RAS_PBL, RAS_WFN,   &
         RAS_TAU,                                         &
         
         
         !AEROPROPS,                                   & !!!!!AER_CLOUD
         CNV_FICE,                                        & 
         CNV_NICE,                                        & 
         CNV_NDROP,                                    &   !DONIF
         RAS_ALPHA,                                     &
         
         ITRCR, IRC, XHO,                                  &
         TRIEDLEV_DIAG,                                   &
         FSCAV , DISSKE                                   )


      !*********************************************************************
      !*********************************************************************
      !******************** Relaxed Arakawa-Schubert ***********************
      !************************ Parameterization ***************************
      !********************** SCALAR RAS-1 VERSION  ************************
      !************************* 31 DECEMBER 1999 **************************
      !*********************************************************************
      !************************** Developed By *****************************
      !*********************************************************************
      !************************ Shrinivas Moorthi **************************
      !******************************* and *********************************
      !************************** Max J. Suarez ****************************
      !*********************************************************************
      !******************** Laboratory for Atmospheres *********************
      !****************** NASA/GSFC, Greenbelt, MD 20771 *******************
      !*********************************************************************
      !*********************************************************************

      !  Input:
      !  ------
      ! 
      !     K0      : Number of vertical levels (increasing downwards)
      !
      !     DT      : Time step in seconds
      !
      !     RASAL   : Array of dimension K-1 containing relaxation parameters
      !               for cloud-types detraining at those levels
      !
      !     CPO     : Specific heat at constant pressure (J/kg/K)
      !
      !     ALHLO   : Latent Heat of condensation (J/kg)
      !
      !     ALHL1   : Latent Heat of condensation + fusion (J/kg)
      !
      !     GRAVO   : Acceleration due to gravity (m/s^2)
      !
      !     PLE     : 2D array of dimension (IDIM,K0+1) containing pressure
      !               in hPa at the interfaces of K-layers from top of the 
      !               atmosphere to the bottom  (mb)
      !
      !     PKE     : 2D array of dimension (IDIM,K0+1) containing (PRS/P00) **
      !               RKAP.  i.e. Exner function at layer edges.
      !
      !     PKL     : 2D array of dimension (IDIM,K0) ) containing the
      !               Exner function at the layers.
      !
      !     QSS     : 2D array of dimension (IDIM,K0  ) containing the
      !               saturation specific humidity at the layers. (kg/kg)
      !
      !     DQS     : 2D array of dimension (IDIM,K0  ) containing
      !               d(qss)/dt at the layers.  (1/K)
      !   
      !     CNV_FRACTION    : 1D array of dimension (IDIM) containing
      !               fraction of grid cell considered to be convective
      !   

      !  Update:
      !  -------
      !
      !     THO     : 2D array of dimension (IDIM,K0) containing potential
      !               temperature (K)
      !
      !     QHO     : 2D array of dimension (IDIM,K0) containing specific
      !               humidity (kg/kg)
      !
      !     UHO     : 2D array of dimension (IDIM,K0) containing u-wind (m/s)
      !
      !     VHO     : 2D array of dimension (IDIM,K0) containing v-wind (m/s)
      !
      !  Output:
      !  -------
      !!
      !     CLW     : 2D array of dimension (IDIM,K0) containing the
      !               detrained cloud liquid water.  (kg/m^2/s) 
      !
      !     FLX     : 2D array of dimension (IDIM,K0) containing the
      !               cloud-base mass flux for each cloud type ordered by
      !               detrainment level.   (kg/m^2/s) 
      !
      !     FLXD    : 2D array of dimension (IDIM,K0) containing the
      !               detrained  mass flux for each cloud type ordered by
      !               detrainment level.   (kg/m^2/s) 
      !
      !     FLXC    : 2D array of dimension (IDIM,K0+1) containing the
      !               total cloud mass flux for all cloud types through
      !               the top of each level. (e.g., FLXC(K)=SUM(FLX(ICMIN:K))
      !               and  FLXD(L) = FLXC(L+1)-FLSD(L) )
      !                (kg/m^2/s) 
      !
      !     PRECU   : 1D (IDIM) Locally-handled convective precip
      !               Zero if older version of RAS-1 is used. Nonzero if
      !               RAS-2 is used.

      !
      !   AEROPROPS, Structure containing aerosol propoerties (in)
      !   CNV_NICE, CNV_DROP.  Flux of ice crystals and droplet number at det level (1/m^2/s) (out)
      !   CNV_FICE: Ice fraction in the detrained condensate (out)

      !************************************************************************

      !  ARGUMENTS

      INTEGER,                     INTENT(IN   ) ::  IDIM, IRUN, K0, ICMIN
      REAL, DIMENSION (IDIM,K0  ), INTENT(INOUT) ::  THO, QHO, UHO, VHO, QLO, QIO, CLAN
      REAL, DIMENSION (IDIM,K0+1), INTENT(IN   ) ::  PLE, PKE,PHIE
      REAL, DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  QSS, DQS,PLO,PKO,PHIO
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  CNV_FRACTION, RASAL2_2d
      REAL, DIMENSION (     K0+1), INTENT(IN   ) ::  SIGE
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  FLX, CLW , FLXD
      REAL, DIMENSION (IDIM,K0+1), INTENT(  OUT) ::  FLXC
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CNV_PRC3
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  CNV_UPDFRC, CNV_QC, CNV_CVW
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  ENTLAM
      REAL, DIMENSION (IDIM,K0  ), INTENT(  OUT) ::  HHO, HSO
      REAL,                        INTENT(IN   ) ::  DT,  CPO, ALHLO, GRAVO
      REAL,                        INTENT(IN   ) ::  ALHL1, TICE
      INTEGER,                     INTENT(  IN)  ::  ITRCR
      INTEGER, DIMENSION (IDIM,2), INTENT(IN   ) ::  SEEDRAS
      INTEGER, DIMENSION (IDIM  ), INTENT(IN   ) ::  IRAS,JRAS,KCBL
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  ZCBL,TPERT,QPERT
      REAL, DIMENSION (IDIM     ), INTENT(IN   ) ::  CO_AUTO
      REAL, DIMENSION (IDIM     ), INTENT(  OUT) ::  MXDIAM
      REAL, DIMENSION (IDIM     ), INTENT(  OUT) ::  PRECU
      REAL, DIMENSION (IDIM,K0  ), INTENT(IN   ) ::  WGT0,WGT1

      !     REAL, DIMENSION(:),          INTENT(IN   ) ::  RASPARAMS
      type (RASPARAM_TYPE),        INTENT(IN   ) ::  RASPARAMS
      LOGICAL,                     INTENT(IN   ) ::  RAS_NO_NEG  ! Whether to guard against negatives

      REAL, DIMENSION (IDIM     ), INTENT(  OUT) :: RAS_TIME, RAS_TRG, RAS_TOKI, RAS_PBL, RAS_WFN 

      INTEGER, DIMENSION(IDIM,K0), INTENT(  OUT) ::  IRC
      REAL,    OPTIONAL,           INTENT(INOUT) ::  XHO(IDIM,K0,ITRCR)
      REAL,    OPTIONAL,           INTENT(  OUT) ::  TRIEDLEV_DIAG(IDIM,K0)
      REAL,    OPTIONAL,           INTENT(  OUT) ::  DISSKE(IDIM,K0)

      REAL,    OPTIONAL,           INTENT(IN)    :: FSCAV(ITRCR) ! Fraction scavenged per km
      ! = 0 (no scav), = 1 (full scav)

      !  LOCALS

      REAL, DIMENSION (IDIM,K0)  ::  DELP   ! MEM - needed for fill_z

      REAL,  DIMENSION(K0) :: POI_SV, QOI_SV, UOI_SV, VOI_SV
      REAL,  DIMENSION(K0) :: POI, QOI, UOI, VOI,  DQQ, BET, GAM, CLL
      REAL,  DIMENSION(K0) :: POI_c, QOI_c
      REAL,  DIMENSION(K0) :: PRH,  PRI,  GHT, DPT, DPB, PKI, DISSK0,DISSK1, CLANTND
      REAL,  DIMENSION(K0) :: TCU, QCU, UCU, VCU,  CLN, RNS, POL,DM
      REAL,  DIMENSION(K0) :: QST, SSL,  RMF, RNN, RN1, RMFC, RMFP
      REAL,  DIMENSION(K0) :: GMS, ETA, GMH, EHT,  GM1, HCC, RMFD
      REAL,  DIMENSION(K0) :: HOL, HST, QOL, ZOL, HCLD, CLL0,CLLX,CLLI,CLLB
      REAL,  DIMENSION(K0) :: WSP, LAMBDSV, BKE , CVW, UPDFRC
      REAL,  DIMENSION(K0) :: TAU, RASAL, MTKWI, UPDFRP,BK2,BK3,DLL0,DLLX

      REAL,  DIMENSION(ITRCR) :: XHT

      ! MEM - added XFF for unpacking in STRAP
      REAL,  DIMENSION(K0,ITRCR) :: XOI, XCU, XOI_SV, XFF

      REAL,  DIMENSION(K0+1) :: PRJ, PRS, QHT, SHT ,ZET, XYD, XYD0

      INTEGER,  DIMENSION(K0-1) :: RC

      INTEGER :: K,MY_PE

      REAL, DIMENSION(IDIM,K0) :: LAMBDSV2

      REAL TX2, TX3, UHT, VHT, AKM, ACR, ALM, TTH, QQH, SHTRG, WSPBL, DQX
      REAL WFN, TEM, TRG, TRGEXP, EVP, WLQ, QCC,MTKW_MAX !, BKE
      REAL SHTRG_FAC, SIGE_MINHOL, WFNOG

      INTEGER I, IC, L, ICL , ITR , ICL_C, N_DTL
      INTEGER NDTLEXPON

      INTEGER,  ALLOCATABLE, DIMENSION(:) :: ICL_V

      !  RASE GLOBAL CONSTANTS


      REAL GRAV, CP, ALHL, CPBG, ALHI, CPI, GRAVI, DDT, LBCP, OBG, AFC 

!!!!!!!!!
      REAL FRICFAC,DPTH_BL,WUPDRFT,PBLFRAC,AUTORAMPB,CO_ZDEP
      REAL RASAL1, RASAL2, RASAL2i, CO_T, RASNCL,FRICLAMBDA,SDQVT1,SDQV2 
      REAL LAMBDA_FAC,STRAPPING,ACRITFAC,HMINTRIGGER,LLDISAGGXP
      REAL LAMBMX_FAC, DIAMMN_MIN,RDTLEXPON, CLI_CRIT,SDQV3, MAXDALLOWED_D, MAXDALLOWED_S, MAXDALLOWED_E
      REAL RHMN, RHMX, CLDMICRO, FDROP_DUST, FDROP_SOOT, RASAL_SLOPE, FDROP_SEASALT
      INTEGER KSTRAP

      real cld_radius, areal_frac, spect_mflx, cvw_cbase
!!!!!!!!!

      REAL, PARAMETER :: ONEPKAP = 1.+ 2./7., DAYLEN = 86400.0
      !      REAL, PARAMETER :: PBLFRAC = 0.5
      REAL, PARAMETER :: RHMAX   = 0.9999

      !  LAMBDA LIMITS
      REAL            :: LAMBDA_MIN
      REAL            :: LAMBDA_MAX

      !  TRIGGER PARAMETERS
      real, parameter :: RHO_W  =  1.0e3  ! Density of liquid water in kg/m^3

      LOGICAL DYNA_STRAPPING, DO_TRACERS, SMOOTH_HST

      character(len=ESMF_MAXSTR)          :: CBL_STYLE

      real*8, dimension(k0) :: tcu8, qcu8, pcu, flx8
      real*8, dimension(k0,itrcr+2) :: rcu
      real*8    :: cup  !, dpd, tla
      logical :: revap, wrkfun, calkpb, crtfun, lprnt, dndrft

      real*8, dimension(k0) :: toi8, qoi8, prsm8, phil8, qli8, qii8, trcfac
      real*8, dimension(k0) :: ALFIND,ALFINT,ALFINQ,RHC_LS
      real*8, dimension(k0+1) :: prs8, phih8
      real*8, dimension(k0,itrcr+2) :: roi8
      real*8 :: FRACBL, dt8, rasalf
      integer :: KPBL



      real*8 :: MAX_NEG_BOUY = 1.0 ! no inhibition for =1.0
      !!real*8 :: ALFINT = 0.5
      !!real*8 :: ALFINQ = 0.5
      real*8 :: RHFACL = 0.0 ! not used
      real*8 :: RHFACS = 0.0 ! no inhibition
      real*8 :: GAREA  = 1.e10 ! 1 degree resolution
      !!real*8 :: ALFIND = 1.0
      !!real*8 :: RHC_LS = 0.80
      real*8 :: dsfc   = 0.001
      real*8 :: cd     = 1.e-3
      real*8 :: wfnc   = 0.0
      real*8 :: tla    = -1.0
      real*8 :: dpd    = 300.

      !  SCAVANGING RELATED PARAMETERS
      REAL                   :: DELZKM  ! layer thickness in km
      REAL                   :: FNOSCAV ! fraction of tracer *not* scavenged
      REAL                   :: FSCAV_(ITRCR) ! Fraction scavenged per km

      ! ************************AER_CLOUD *********************************************
      ! *****the AEROPROPS part of this code needs to be updated in this code is to be used again
      ! DONIF 3/10/25


    
      REAL, DIMENSION (IDIM,K0), INTENT( OUT) ::  CNV_NDROP, CNV_FICE, CNV_NICE
      REAL, DIMENSION (IDIM,K0), INTENT( OUT) ::  RAS_TAU,  RAS_ALPHA
      REAL,  DIMENSION(K0) :: CNVNDROP, CNVNICE, CNVFICE !DONIF

      TYPE(AerProps), DIMENSION(K0) ::  AERO   !AEROSOL VERTICAl ARRAY FOR ALL SPECIES    
      REAL ::  T_ICE_ALL, T_ICE_MAX, ASEASALT, f_seasalt 
      INTEGER, PARAMETER :: NDUSTMAX = 10

      INTEGER :: INDEX        
      
      TYPE(AerProps),    DIMENSION(IDIM, K0)    :: AEROPROPS !DONIF
      TYPE(AerProps) :: AERAUX, AER_BASE
      

      CNV_FICE  =0.0
      CNV_NDROP =0.0
      CNV_NICE  =0.0
      CNVFICE  =0.0
      CNVNDROP =0.0
      CNVNICE  =0.0
      T_ICE_ALL= 238.0
      T_ICE_MAX= MAPL_TICE
      CLDMICRO = 0.0
      FDROP_DUST = 0.1
      FDROP_SOOT = 0.01
      RAS_ALPHA = MAPL_UNDEF
      RAS_TAU = MAPL_UNDEF

      ! *********************************************************************

      IF ( PRESENT(FSCAV) ) then
         FSCAV_ = FSCAV
      ELSE
         FSCAV_ = 0.0 ! NO SCAVENGING BY DEFAULT
      ENDIF

      IF(IRUN <= 0) RETURN

      IF (PRESENT(TRIEDLEV_DIAG)) TRIEDLEV_DIAG=0.
      CNV_PRC3  =0.
      CNV_UPDFRC=0. 
      CNV_CVW   =0. 
      CNV_QC    =0.
      ENTLAM    =0.
      if ( PRESENT( DISSKE ))     DISSKE    =0.
      !!LAMBDSV2 = 0.
      IRC = -2



      !      SMOOTH_HST   = .TRUE. 
      SMOOTH_HST   = .FALSE.  

      FRICFAC      = RASPARAMS%CUFRICFAC
      SHTRG_FAC    = RASPARAMS%SHR_LAMBDA_FAC

      ! MAT CO_AUTO is now passed in from outside
      !     in order to allow this code to run over im*jm 
      !     columns
      !CO_AUTO      = RASPARAMS(3)     !  ---  3

      CLI_CRIT      = RASPARAMS%QC_CRIT_CN
      RASAL1        = RASPARAMS%RASAL1
      RASAL2        = RASPARAMS%RASAL2
      RASNCL        = RASPARAMS%RASNCL
      LAMBDA_FAC    = RASPARAMS%LAMBDA_FAC
      LAMBMX_FAC    = RASPARAMS%LAMBMX_FAC
      DIAMMN_MIN    = RASPARAMS%MIN_DIAMETER
      FRICLAMBDA    = RASPARAMS%CUFRICLAMBDA
      RDTLEXPON     = RASPARAMS%RDTLEXPON
      STRAPPING     = RASPARAMS%STRAPPING
      SDQV2         = RASPARAMS%SDQV2
      SDQV3         = RASPARAMS%SDQV3
      SDQVT1        = RASPARAMS%SDQVT1
      ACRITFAC      = RASPARAMS%ACRITFAC
      HMINTRIGGER   = RASPARAMS%HMINTRIGGER
      LLDISAGGXP    = RASPARAMS%LLDISAGGXP
      PBLFRAC       = RASPARAMS%PBLFRAC
      AUTORAMPB     = RASPARAMS%RASAUTORAMPB
      CO_ZDEP       = RASPARAMS%AUTOC_CN_ZDEP
      MAXDALLOWED_S = RASPARAMS%MAXDALLOWED_S
      MAXDALLOWED_D = RASPARAMS%MAXDALLOWED_D
      MAXDALLOWED_E = RASPARAMS%MAXDALLOWED_E
      RASAL_SLOPE     = RASPARAMS%RASAL_SLOPE
      RHMN          = RASPARAMS%RAS_RHMIN
      RHMX          = RASPARAMS%RAS_RHFULL
      CLDMICRO      = RASPARAMS%CLDMICRO
      FDROP_DUST    = RASPARAMS%FDROP_DUST
      FDROP_SOOT    = RASPARAMS%FDROP_SOOT
      FDROP_SEASALT    = RASPARAMS%FDROP_SEASALT

      IF ( STRAPPING <= 0.0 ) THEN
         DYNA_STRAPPING = .TRUE.
      ELSE
         DYNA_STRAPPING = .FALSE.
         KSTRAP      = INT( STRAPPING )
      ENDIF


      DO_TRACERS = present(XHO) .and. ITRCR>0

      WUPDRFT = 2.500

      GRAV  = GRAVO
      ALHL  = ALHLO
      CP    = CPO
      CPI   = 1.0/CP      
      ALHI  = 1.0/ALHL
      GRAVI = 1.0/GRAV
      CPBG  = CP*GRAVI
      DDT   = DAYLEN/DT
      AFC   = -1.04E-4*SQRT(DT*113.84)
      LBCP  = ALHL*CPI
      OBG   = 100.*GRAVI

      HHO = MAPL_UNDEF
      HSO = MAPL_UNDEF

#ifdef RAS2
      dt8    =  dt

      call set_ras_afc(dt8)
      call ras_init(k0,1)

      kpbl  = k0
      lprnt = .false.
      trcfac(:) = 1.0
      revap  = .true.
      wrkfun = .false.
      calkpb = .true.
      crtfun = .true.
      dndrft = .true.
      wfnc   = 0.0

#endif

      ! MEM
      DELP(:,ICMIN:K0) = PLE(:,ICMIN+1:K0+1) - PLE(:,ICMIN:K0)

      DO I=1,IRUN

         if (CLDMICRO > 0.0) then !===================AER_CLOUD

            AERO =  AEROPROPS(I, :)
            CNVFICE  =0.0
            CNVNDROP =0.0
            CNVNICE  =0.0

         end if

#ifndef RAS2

         !!CALL FINDBASE
         K = KCBL(I)

         rc(icmin) = 0

         CALL FINDDTLS



         IF ( K > 0 ) THEN 
            CALL STRAP( FINAL=0 )
            call HTEST

            HHO(I,:) = HOL
            HSO(I,:) = HST

            TAU = 0.0
            RASAL = 0.0
            DO ICL_C = 1,N_DTL
               ICL   = ICL_V( ICL_C )


               IF (DO_TRACERS) XCU(ICMIN:,:) = 0.

               IF (PRESENT(TRIEDLEV_DIAG)) TRIEDLEV_DIAG(I,ICL) = 1.

               UCU(ICMIN:) = 0.    ! This change makes cumulus friction
               VCU(ICMIN:) = 0.    ! correct.



               IF( ICL > ICMIN ) THEN
                  CALL CLOUDE(ICL)
               ENDIF

               ENTLAM(I,ICL) = ALM
               RAS_TAU(I, ICL) = TAU(ICL)
               RAS_ALPHA(I, ICL) =  RASAL(ICL) 
               CNV_NDROP   (I,ICL)  =    CNVNDROP(ICL) !DONIF
               CNV_NICE   (I,ICL)   =     CNVNICE(ICL) !DONIF
               CNV_FICE   (I,ICL)   =     CNVFICE(ICL) !DONIF

            ENDDO

            IF ( SUM( RMF(ICMIN:K) ) > 0.0 ) THEN

               CALL RNEVP

               CALL STRAP( FINAL=1 )

            ELSE

               CALL STRAP( FINAL=2 )

            ENDIF

         ELSE 

            CALL STRAP( FINAL=2 )

         ENDIF

         PRECU = 0.0  ! Zero out precip - TBD w/ in progno_cloud
#else


         K = K0
         CALL FINDDTLS

         !===>  K      INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
         !===>  KD     INPUT   DETRAINMENT LEVEL ( 1<= KD < K )
         !===>  M      INPUT   NUMBER OF TRACERS. MAY BE ZERO.
         !===>  RASALF INPUT
         !===>  FRACBL INPUT   MASS FLUX LIMIT
         !===>  MAX_NEG_BOUY  INPUT  FRACTION OF NEG ALLOWED
         !===>  ALFINT(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR T
         !===>  ALFINQ(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR Q
         !===>  RHFACL INPUT   CRITICAL RH AT PBL TOP (LAND)
         !===>  RHFACS INPUT   CRITICAL RH AT PBL TOP (SEA) ???
         !===>  GAREA  INPUT   AREA OF GRID BOX
         !===>  ALFIND(K) INPUT   INTERFACE UPWIND (=1) PARAMETER FOR DOWNDRAFT
         !===>  RHC_LS(K) INPUT   CRITICAL RH FOR RE-EVAP

         ALFINT = 0.5
         ALFINQ = 0.5
         ALFIND = 0.5
         RHC_LS = 0.8

         tcu8  = 0.0
         qcu8  = 0.0
         rcu  = 0.0
         pcu  = 0.0
         flx8 = 0.0
         cup  = 0.0

         cloudloop: DO ICL = ICMIN,K0-1
            IF (DO_TRACERS) XCU(ICMIN:,:) = 0.


            rasalf = rasal(ICL)
            FRACBL = pblfrac

            toi8 = pko(i,:)*tho(i,:) ! tbdone
            qoi8 = qho(i,:)
            roi8(:,1) = uho(i,:)
            roi8(:,2) = vho(i,:)
            if(present(xho)) roi8(:,3:size(xho,3)+2) = xho(i,:,:) 

            prs8  = ple(i,:)
            prsm8 = plo(i,:) ! tbdone
            phil8 = phio(i,:) ! tbdone
            phih8 = phie(i,:) ! tbdone
            qli8  = qlo(i,:) ! tbdone
            qii8  = qio(i,:) ! tbdone

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

            call CLOUD(                                                       &
                  &                  K0, ICL, size(xcu,2)                            &
                  &,                 RASALF, FRACBL, MAX_NEG_BOUY                    &
                  &,                 ALFINT, ALFINQ, RHFACL, RHFACS, garea           &
                  &,                 ALFIND, RHC_LS                                  &
                  
                  &,                 TOI8, QOI8, ROI8, PRS8, PRSM8, phil8, phih8     &
                  &,                 QLI8, QII8, KPBL, DSFC                            &
                  &,                 CD,lprnt, trcfac                                &
                  
                  &,                 TCU8, QCU8, RCU, PCU, FLX8                         &
                  &,                 CUP, REVAP, DT8                                 &
                  &,                 WFNC, WRKFUN, CALKPB, CRTFUN, TLA, DnDRFT, DPD)  

            tho(i,:)  = toi8/pko(i,:)
            qho(i,:)  = qoi8
            uho(i,:)  = roi8(:,1) 
            vho(i,:)  = roi8(:,2)
            qlo(i,:)  = qli8 ! tbdone
            qio(i,:)  = qii8  ! tbdone


            if(present(xho)) xho(i,:,:) =roi8(:,3:size(xho,3)+2)
         ENDDO cloudloop

         ! temporary zeroing of outputs needed by rest of moist
         !-----------------------------------------------------
         cnv_prc3(i,:)   = 0.0
         clw(i,:)        = 0.0
         HHO(i,:)        = 0.0
         HSO(i,:)        = 0.0
         cnv_updfrc(i,:) = 0.0001
         cnv_qc(i,:)     = 0.0
         flxc(i,1:k0)    = flx8 / DT
         flx(i,:)        = 0.0

         flxd(i,1:k0-1)  = MAX( flxc(i,2:k0)- flxc(i,1:k0-1) , 0. )
         flxd(i,k0)      = 0.0

         precu(i)        = cup / DT   ! 

         !
         ! Tiedtke-style prognostic anvil fraction. 
         ! Normally done in progno_cloud. Moved here
         ! since RAS-2 manages its own detraining 
         ! ice and liquid condensate.
         !---------------------------------------------
         dm           = 100.*(prs8(2:k0+1)-prs8(1:k0) )/grav
         CLANTND      = FLXD(i,:) / DM
         CLAN(i,:)    = CLAN(i,:) + CLANTND*DT
         CLAN(i,:)    = MIN( CLAN(i,:) , 0.99 )

         ! Now detraining mass flux has to be zeroed
         ! or clouds will be crated twice
         !--------------------------------
         flxd(i,:)    = 0.

#endif

      ENDDO

      IF(ALLOCATED(ICL_V)) DEALLOCATE(ICL_V)

      RETURN

   CONTAINS

      !*********************************************************************

      SUBROUTINE CLOUDE(IC)

         INTEGER, INTENT(IN ) :: IC
         REAL :: DEEP_FACT,CU_DIAM,WSCALE

         REAL :: CLI , TE_A, C00_X, CLI_CRIT_X, PETE, TOKI, GMHx, HSTx  !, dQx
         REAL :: DT_LYR, RATE, CVW_X, CLOSS, F2, F3, F4, F5
         INTEGER :: K700

         !=============================AER_CLOUD local variables ====================
         REAL :: WBASE, NDROP, NICE, FP_D, FF_A, FP_I, FICE, &
               NDROP_AMB, NSOOT_AMB, NSOOT, NIN, INSOOT, dCVW2, QICE, &
               dQICE, dQIG, FPICE, dNICE, dNDROP, DSOOT_AMB, DSOOT, QLIQ, dQLIQ, FPRECIP, AUX, QT, &
               MAXNICE, MAXNDROP, MINNICE, MINNDROP, NDROP_ACT, RIMM, FNDRIM, TminusTa, Tparcel , &
	      alph_e, beta_e, RH_AMB, ECRIT

         real :: lamb, minh, maxh, max_alpha, min_alpha

         REAL, DIMENSION (NDUSTMAX) :: NDUST, NDUST_AMB, INDUST, DDUST_AMB, DDUST 
         INTEGER :: INX, naux, INDEX


         T_ICE_ALL = 238.0 !A little higher so there is ice at the freezing level
         WBASE=1.0 
         FICE=0.0
         NICE=0.0
         NDROP=0.0
         NDROP_AMB=0.0
         NDROP_ACT=0.0
         NIN = 0.0
         NDUST=0.0
         NSOOT=0.0
         NDUST_AMB=0.0
         NSOOT_AMB = 0.0
         dCVW2=0.0
         QICE=0.0
         QLIQ=0.0
         FPICE = 0.0
         INDUST=0.0
         INSOOT= 0.0
         QT= 0.0
         FPRECIP=0.0   
         FNDRIM = 0.0
         RIMM = 0.0   
         TminusTa = 0.0

         f_seasalt = 0.0
         aseasalt = 0.0

         call init_Aer(AER_BASE)

         !AER_CLOUD=============================


         ALM   = 0.
         TRG   = AMIN1(1.,(QOI(K)/QST(K)-RHMN)/(RHMX-RHMN))

         F4  = MIN(   1.0,  MAX( 0.0 , (AUTORAMPB-SIGE(IC))/0.2 )  )  ! F4 should ramp from 0 at SIG=AUTORAMPB
         ! to 1 at SIG=AUTORAMPB-0.2

         if ( SIGE(IC) >= 0.5 ) then
            F5 = 1.0
         else
            F5 = 1.0 - 2.*CO_ZDEP *( 0.5 - SIGE(IC) )
            F5 = MAX( F5 , 0.0 )
         endif


         IF(TRG <= 1.0E-5) THEN    ! TRIGGER  =========>>
            RC(IC) = 7
            RETURN
         ENDIF

         !  RECOMPUTE SOUNDING UP TO DETRAINMENT LEVEL

         POI_c = POI
         QOI_c = QOI
         POI_c(K) =  POI_c(K) + TPERT(I)
         QOI_c(K) =  QOI_c(K) + QPERT(I)

         ZET(K+1) = 0.
         SHT(K+1) = CP*POI_c(K)*PRJ(K+1)
         DO L=K,IC,-1
            QOL(L)  = AMIN1(QST(L)*RHMAX,QOI_c(L))
            QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
            SSL(L)  = CP*PRJ(L+1)*POI_c(L) + GRAV*ZET(L+1)
            HOL(L)  = SSL(L) + QOL(L)*ALHL
            HST(L)  = SSL(L) + QST(L)*ALHL
            TEM     = POI_c(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
            ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI_c(L)*CPBG
         ENDDO

         DO L=IC+1,K
            TEM  = (PRJ(L)-PRH(L-1))/(PRH(L)-PRH(L-1))
            SHT(L)  = SSL(L-1) + TEM*(SSL(L)-SSL(L-1)) 
            QHT(L)  = .5*(QOL(L)+QOL(L-1))
         ENDDO

         ! SMOOTH HSTAR W/ 1-2-1 Filter
         if ( SMOOTH_HST ) then
            HSTx = HST(IC) ! save for later
            DO L=K-1,IC+1,-1
               HST(L) = 0.25*(HST(L+1)+HST(L-1))+0.5*HST(L)
            END DO
            DO L=IC,IC
               HST(L) = 0.5*HST(L+1)+0.5*HST(L)
            END DO
         endif



         !  CALCULATE LAMBDA, ETA, AND WORKFUNCTION


         LAMBDA_MIN = .2/MXDIAM(I)
         LAMBDA_MAX = .2/DIAMMN_MIN 

         !     LAMBDA_MIN = .2/(LAMBDA_FAC*DPTH_BL)
         !     LAMBDA_MAX = .2/( MAX( LAMBMX_FAC*DPTH_BL , DIAMMN_MIN ) )

         IF (HOL(K) <= HST(IC)) THEN   ! CANNOT REACH IC LEVEL  ======>>
            RC(IC) = 1
            RETURN
         ENDIF

         !  LAMBDA CALCULATION: MS-A18
         TEM  =       (HST(IC)-HOL(IC))*(ZOL(IC)-ZET(IC+1)) 
         DO L=IC+1,K-1
            TEM = TEM + (HST(IC)-HOL(L ))*(ZET(L )-ZET(L +1))
         ENDDO
         IF(TEM <= 0.0) THEN         ! NO VALID LAMBDA  ============>>
            RC(IC) = 2
            RETURN
         ENDIF

         ALM     = (HOL(K)-HST(IC)) / TEM
         IF(ALM > LAMBDA_MAX) THEN
            RC(IC) = 3
            RETURN
         ENDIF

         !   ALPHA CALCULATION 
         RASAL2i = RASAL2_2d(I)
         IF ( ZET(IC) <  2000. ) RASAL(IC) = RASAL1
         IF ( ZET(IC) >= 2000. ) THEN
            RASAL(IC) = RASAL1 + (RASAL2i-RASAL1)*MIN(1.0,(ZET(IC) - 2000.)/RASAL_SLOPE)
         ENDIF
         RASAL(IC) = DT / RASAL(IC)

         !   Tokioka Calculation
         TOKI=MIN(1.0,(ALM/LAMBDA_MIN)**2)

         !   RAS relaxation timescale
         IF (k0 > 96) THEN
          !! AMM kluge to run 132 levels for now -- multiply time scale by ratio of number of levels
          !! AMM between 900 and 30 mb in 72 and 132 level grids, 0.46
            TAU(IC) = TOKI*TRG*RASAL(IC) *0.46
         ELSE
            TAU(IC) = TOKI*TRG*RASAL(IC)
         ENDIF

        !IF (TAU(IC) < 1.0E-5) THEN
        !   RC(IC) = 6
        !   RETURN
        !ENDIF

         !LAMBDSV(IC) = ALM

         !  ETA CALCULATION: MS-A2

         DO L=IC+1,K
            ETA(L) = 1.0 + ALM * (ZET(L )-ZET(K))
         ENDDO
         ETA(IC) = 1.0 + ALM * (ZOL(IC)-ZET(K))

         !  WORKFUNCTION CALCULATION:  MS-A22

         WFN     = 0.0
         HCC(K)  = HOL(K)
         DO L=K-1,IC+1,-1
            HCC(L) = HCC(L+1) + (ETA(L) - ETA(L+1))*HOL(L)
            TEM    = HCC(L+1)*DPB(L) + HCC(L)*DPT(L)
            EHT(L) = ETA(L+1)*DPB(L) + ETA(L)*DPT(L)
            WFN    = WFN + (TEM - EHT(L)*HST(L))*GAM(L)
         ENDDO
         HCC(IC) = HST(IC)*ETA(IC)
         WFN     = WFN + (HCC(IC+1)-HST(IC)*ETA(IC+1))*GAM(IC)*DPB(IC)



         !  VERTICAL VELOCITY/KE CALCULATION (ADDED 12/2001 JTB)
         BK3(K)   = 0.0
         BK2(K)   = 0.0
         BKE(K)   = 0.0
         HCLD(K)  = HOL(K)
     
         DO L=K-1,IC,-1
            HCLD(L) = ( ETA(L+1)*HCLD(L+1)   +      & 
                  (ETA(L) - ETA(L+1))*HOL(L) ) / ETA(L)
            TEM     = (HCLD(L)-HST(L) )*(ZET(L)-ZET(L+1))/ (1.0+LBCP*DQQ(L))
       
       if (CLDMICRO .le. 0.0) then  
       
            BKE(L)  = BKE(L+1) + GRAV * TEM / ( CP*PRJ(L+1)*POI(L) )
            BK2(L)  = BK2(L+1) + GRAV * MAX(TEM,0.0) / ( CP*PRJ(L+1)*POI(L) )
            BK3(L)  = BK3(L+1) + GRAV * MIN(TEM,0.0) / ( CP*PRJ(L+1)*POI(L) )
            CVW(L) = SQRT(  2.0* MAX( BK2(L) , 0.0 )  )    
       end if 
         ENDDO


         !  HEAL/REPAIR/RATIONALIZE STUPID CALCULATION OF CVW
         !  IT IS STUPID BECAUSE IT IGNORES A NUMBER OF IMPORTANT EFFECTS
         !    1) CONVECTIVE PRESSURE PERT. COULD GO EITHER WAY + OR - 
         !    2) CONDENSATE LOADING/MECHANICAL DRAG
         !    3) ENVIRONMENTAL CONDENSATE COULD GO +
         !    4) SUBGRID SCALE CLOUD BASE ENHANCEMENT OF Q AND T
         !
         !  SO, HERE AFTER NAIVELY USING ITS VALUE TO TEST FOR CLOUD TOP
         !  W, WE FLOOR IT AT 1 M/S.  THIS ENSURES THAT IF A PLUME
         !  MAKES IT THROUGH ALL TEST CONDITIONS, I.E., GETS A MASS-FLUX 
         !  IT WILL ALSO HAVE SOME DIAGNOSED "UPWARD MOTION"

         CVW(IC:K) = MAX(  CVW(IC:K) , 1.00 )

         !  NOTE THIS "CENTRALIZES" A KLUGE PRESENT IN OTHER LOCATIONS.
         !  CLEAN UP SOME TIME.      -JTB 12/04/03



         !  TEST FOR CRITICAL WORK FUNCTION

         CALL ACRITN(POL(IC), PRS(K), ACR)




         IF(WFN <= ACR) THEN   ! SUB-CRITICAL WORK FUNCTION ======>>
            RC(IC) = 4
            RETURN
         ENDIF

         !  CLOUD TOP WATER AND MOMENTUM (TIMES ETA(IC)) MS-A16

         ! Tracer scavenging
         ! RAS loops over a series of plumes all having common cloud base level K 
         ! and different detrainment levels IC.  The plumes operate sequentially
         ! on the grid box mean quantities (wind, moisture, tracer) and so each
         ! subsequent plume is seeing the effects of previous plumes.  We parameterize
         ! scavenging following Liu et al. [JGR, 2001], their equation 1:
         !  AEROSOL FRACTION SCAVENGED = 1 - exp(-FSCAV*DZ)
         ! where FSCAV is a specified scavenging efficiency [km-1] and DZ is the
         ! distance [km] the tracer traverses in the plume from its entrainment
         ! level to its detrainment level.  We write the aerosol fraction surviving as:
         !  FNOSCAV = exp(- FSCAV_(ITR) * DZ)
         ! The total scavenging is proportional to the convective mass flux, which
         ! is not explicitly solved for at this point.

         IF (DO_TRACERS) THEN 
            DO ITR=1,ITRCR
               !           Scavenging of the below cloud tracer
               DELZKM = ( ZET(IC) - ZET(K) ) / 1000.
               FNOSCAV = max( min(exp(- FSCAV_(ITR) * DELZKM),1.), 0.)
               XHT(ITR) = XOI(K,ITR) * FNOSCAV
            END DO
         END IF

         WLQ     = QOL(K)
         UHT     = UOI(K)
         VHT     = VOI(K)
         RNN(K)  = 0.
         CLL0(K) = 0.

         !print *, '========================================='
         DO L=K-1,IC,-1
            TEM   = ETA(L) - ETA(L+1)
            WLQ   = WLQ + TEM * QOL(L)
            UHT   = UHT + TEM * UOI(L)
            VHT   = VHT + TEM * VOI(L)

            IF (DO_TRACERS)  THEN  
               DO ITR=1,ITRCR

                  !         Scavenging of the entrained tracer.  Updates transported tracer mass.
                  DELZKM = ( ZET(IC) - ZET(L+1) ) / 1000.
                  FNOSCAV = max( min(exp(- FSCAV_(ITR) * DELZKM),1.), 0.)
                  XHT(ITR) = XHT(ITR) + TEM * XOI(L,ITR) * FNOSCAV

               END DO
            END IF

!!!! How much condensate (CLI) is present here? 
            IF(L>IC) THEN
               TX2   = 0.5*(QST(L)+QST(L-1))*ETA(L)
               TX3   = 0.5*(HST(L)+HST(L-1))*ETA(L)
               QCC   = TX2 + GM1(L)*(HCC(L)-TX3)
               CLL0(L)  = (WLQ-QCC)
            ELSE
               CLL0(L)   = (WLQ-QST(IC)*ETA(IC))
            ENDIF

            CLL0(L)    = MAX( CLL0(L) , 0.00 )

            CLI  = CLL0(L) / ETA(L)  ! condensate (kg/kg)
            TE_A = POI(L)*PRH(L)     ! Temperature (K)


            !=====================================================================
            if (CLDMICRO .gt. 0.0) then  !AER_CLOUD MIcrophysics considering activation and nucleation 
               !recompute vertical velocity

               Tparcel = TE_A
               CVW(K) = 0.8  ! Assume a below cloud base  W of 0.8 m s-1            
               BK2(K)   = 0.0

     
               TEM     = (HCLD(L)-HST(L) )/ (1.0+LBCP*DQQ(L))  
               TminusTa = max(min(TEM/CP, 5.0), 0.0) !limit DT to 5 K. According to Wei, JAS, 1998   
	     TEM =0.33*TminusTa*CO_AUTO(I)/TE_A !Bouyancy term, effciency =0.5 mwr Roode et al    	     

               BK2(L)  = BK2(L+1) + GRAV * TEM*(ZET(L)-ZET(L+1)) 
               BK2(L) = BK2(L) - (ZET(L)-ZET(L+1))*(BK2(L+1)*ALM + CLI*GRAV)  !Account for drag from entrainment of stagnat air and condesate loading
               CVW(L) = max(SQRT(  2.0* MAX( BK2(L) , 0.0 )  ), 1.0) 


	    CVW_X = MIN(CVW(L), 50.0)
               DT_LYR  =  max(( ZET(L)-ZET(L+1) )/CVW_X, 1.0) !Sanity check 
               TEM   = ETA(L) - ETA(L+1)

               Tparcel  =  TE_A + TminusTa



!!!!!!!!!account for entrainment effects on activation !!!!!!!!!!!
               ! Barahona and Nenes, JGR, 2007
               alph_e = 2.8915e-8*Tparcel*Tparcel -2.1328e-5*Tparcel+4.2523e-3
               beta_e = MAPL_ALHL*TminusTa/MAPL_RVAP/Tparcel/Tparcel
               RH_AMB=QOI(L)/QST(L)
               ECRIT  = max(1.0-RH_AMB -beta_e, 1.0e-6) 
               ECRIT =  alph_e/ECRIT
               ! print *, L, Tparcel, RH_AMB, ECRIT, ALM
	           ECRIT =  ALM/ECRIT
               !Print *, ECRIT


               if (L .eq. K-1) then

                  FICE=0.0
                  NICE=0.0
                  NDROP=0.0
                  NIN =0.0
                  NDUST_AMB =0.0
                  NSOOT_AMB = 0.0
                  NSOOT=0.0
                  NDUST= 0.0
                  RATE=0.0
                  FPRECIP=0.0

                  AER_BASE%nmods = 0
                  AER_BASE%num   = 0.0
                  do INDEX = 1, AERO(L)%nmods
                      if (AERO(L)%num(INDEX) > 0.1) then
                          AER_BASE%nmods = AER_BASE%nmods + 1
                          naux = AER_BASE%nmods

                          AER_BASE%num(naux)   = AERO(L)%num(INDEX)
                          AER_BASE%dpg(naux)   = max(AERO(L)%dpg(INDEX), 1.0e-9)
                          AER_BASE%sig(naux)   = AERO(L)%sig(INDEX)
                          AER_BASE%den(naux)   = AERO(L)%den(INDEX)
                          AER_BASE%kap(naux)   = AERO(L)%kap(INDEX)
                          AER_BASE%fdust(naux) = AERO(L)%fdust(INDEX)
                          AER_BASE%fsoot(naux) = AERO(L)%fsoot(INDEX)
                          AER_BASE%forg(naux)  = AERO(L)%forg(INDEX)
                      end if
                  end do

                  !initial conditions     
                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L,  .true., DDUST_AMB, DSOOT_AMB, ECRIT) !cloud droplet number and INsource at cloud base 
                  NDUST=NDUST_AMB
                  NSOOT=NSOOT_AMB
                  DDUST=DDUST_AMB
                  DSOOT=DSOOT_AMB                                     

               else 
                  call ARGact(Tparcel, CVW_X, NDROP_ACT, NDROP_AMB, NDUST_AMB, NSOOT_AMB, L, .false., DDUST_AMB, DSOOT_AMB, ECRIT) !cloud droplet number above cloud base  

               end if

               QT = CLI
               RATE = 0.0
               FPRECIP = 0.0

               if (QT .gt. 0.0) then

                  ! if FICE is already >= 1.0 then the cloud is glaciated and there is no need to do anymore partitioning

                  if (FICE .ge. 1.0) then


                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
                           DT_LYR, RIMM, CO_AUTO(I)) 



                     dNICE = -NICE*FP_I 
                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice

                     MINNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/4.0e-8!assuming maximum vol radius 250 microns
                     MAXNICE = max(QICE*(1.0-RATE*DT_LYR), 0.0)/2.51e-12 !assuming minimum vol radius 10 microns

                     NICE=MIN(max(NICE, MINNICE), MAXNICE)

                     FICE = 1.0

                  else 

                     ! Cloud is not completely glaciated do the whole thing
                     ! ALL these subroutines return tendencies


                     CALL  INfreezing(QLIQ, NDROP, NIN, NDUST, NSOOT, INDUST, INSOOT, Tparcel, POL(L), CVW_X, DDUST, DSOOT)  !calculate the freezing fraction of the aerosol at this level

                     NIN = min(NIN, NDROP/DT_LYR)

                     call Qgrowth(Tparcel, POL(L), QICE, NICE, QT, NIN, dQIG, RIMM, FNDRIM)

                     CALL  Qremoval(RATE, FICE, FP_D, FP_I, Tparcel,  & 
                           POL(L), QT,  NICE, NDROP, CVW_X, FPICE, &
                           DT_LYR,  RIMM, CO_AUTO(I)) 



                     !ice number tendency: -precip + freezin
                     dNICE = -NICE*FP_I  + NIN     
                     NICE  =  (NICE +dNICE*DT_LYR)*ETA(L+1)/ETA(L) !ice
                     NICE =max(NICE, 0.0)


                     !ice mass tendency: growth - precip
                     dQICE = -QICE*FPICE + dQIG
                     QICE  =  min((QICE + dQICE*DT_LYR)*ETA(L+1)/ETA(L), QT) !ice
                     QICE=max(min(QT, QICE), 0.0)


                     ! Liquid Tendency: source/evap -  precip 
                     !dQLIQ = max((CLI-QICE), -QLIQ)/DT_LYR -QLIQ*max(RATE-FPICE, 0.0) 
                     ! dQLIQ = CLI*(1.0-RATE*DT_LYR)/DT_LYR -dQICE - QLIQ*max(RATE-FPICE, 0.0)           
                     !QLIQ  =  max((QLIQ + dQLIQ*DT_LYR)*ETA(L+1)/ETA(L), 0.0) !liquid. This is actually diagnostic
                     QLIQ=max((QT-QICE), 0.0)


                     !droplet number tendency: -precip - freezin + activation + activated entrained aerosol 


                     dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + max(NDROP_ACT-NDROP, 0.0)/DT_LYR          

                     !dNDROP =-NDROP*FP_D - NIN -  FNDRIM*NDROP/DT_LYR + NDROP_ACT/DT_LYR         

                     NDROP =  (NDROP + dNDROP*DT_LYR)*ETA(L+1)/ETA(L) + &
                           (ZET(L) - ZET(L+1))*ALM*MAX((NDROP_AMB-NDROP), 0.0)

                     !Aerosol tendency: Entrainment - freezing 

                     NDUST = (NDUST - INDUST*DT_LYR)*ETA(L+1)/ETA(L) + &
                           (ZET(L) - ZET(L+1))*ALM*MAX(NDUST_AMB-NDUST, 0.0) 

                     NSOOT =  (NSOOT - INSOOT*DT_LYR)*ETA(L+1)/ETA(L)  + &     
                           (ZET(L) - ZET(L+1))*ALM*MAX(NSOOT_AMB-NSOOT, 0.0)  


                           
                     !Update FICE and perform Sanity checks


                     MINNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/2.e-10    !assuming maximum vol radius 36 microns
                     MAXNDROP = (1.0-FICE)*QLIQ*max(1.0-RATE*DT_LYR, 0.0)/3.35e-14 !assuming minimum vol radius 2 microns
                     MINNICE = QICE/4.0e-8!assuming maximum vol radius 250 microns
                     MAXNICE = QICE/2.51e-12 !assuming minimum vol radius 10 microns

                     IF ((NICE .gt. MAXNICE) .or. (NICE .lt. MINNICE))   then    
                        !print *, 'nilim', NICE*1e-6, MINNICE*1e-6, MAXNICE*1e-6
                     END IF

                     IF ((NDROP .gt. MAXNDROP) .or. (NDROP .lt. MINNDROP))      then 
                        !print *, 'ndroplim', NDROP*1e-6, MINNDROP*1e-6, MAXNDROP*1e-6
                     end if


                     NSOOT=MAX(NSOOT, 0.0)
                     NDUST=MAX(NDUST, 0.0)              

                     NDROP=MIN(max(NDROP, MINNDROP), MAXNDROP)
                     NICE=MIN(max(NICE, MINNICE), MAXNICE)

                     FICE=max(min(QICE/QT, 1.0), 0.0)

                     IF (FICE .ge. 1.0) THEN !Complete glaciation 
                        NICE=NICE+NDROP 
                        NDROP = 0.0
                        QICE  = QT
                        QLIQ= 0.0
                     END IF

                     IF (Tparcel .LT. T_ICE_ALL) THEN !instantaneous freezing
                        NICE=NICE+NDROP 
                        NDROP = 0.0
                        FICE  = 1.0
                        QICE  = QT
                        QLIQ=0.0
                     END IF

                     IF (Tparcel .GT. T_ICE_MAX) THEN !instantaneous melting
                        NDROP=NICE+NDROP 
                        NICE = 0.0
                        FICE  = 0.0
                        QICE  = 0.0
                        QLIQ=QT
                     END IF

                  END IF

               else 

                  FICE =0.0 
                  QICE = 0.0
                  QLIQ = 0.0
                  NICE= 0.0 
                  NDROP = 0.0
                  RATE =0.0
               end if

               FPRECIP= RATE*DT_LYR

               !RATE=RATE*F4
               ! NDROP=NDROP*F4
               !NICE=NICE*(1.0-F4)

               !print *, TE_A, FICE, 'NICE', NICE*1e-6, 'NDROP', NDROP*1e-6, L 
               !print *, 'FPI', FP_I*DT_LYR, 'FPD', FP_D*DT_LYR, 'FPICE', FPICE, 'FPRE', FPRECIP, QT, QLIQ

            else !Bacmeister 2006 microphysics

               CALL SUNDQ3_ICE( TE_A,  SDQV2, SDQV3, SDQVT1, F2 , F3 )

               C00_X  =  CO_AUTO(I) * F2 * F3  * F4  ! * F5  ! F4 reduces AUTO for shallow clouds, F5 modifies auto for deep clouds
               CLI_CRIT_X =  CLI_CRIT / ( F2 * F3 )

               RATE = C00_X * ( 1.0 - EXP( -(CLI)**2 / CLI_CRIT_X**2 ) )
            endif
            !=============================================================================


            CVW_X     = MAX( CVW(L) , 1.00 )              ! Floor conv. vel. since we do not 
                                                          ! really trust it at low values

            DT_LYR  = ( ZET(L)-ZET(L+1) )/CVW_X           ! l.h.s. DT_LYR => time in layer (L,L+1)

            CLOSS   = CLL0(L) * RATE * DT_LYR
            CLOSS   = MIN( CLOSS , CLL0(L) )

            CLL0(L) = CLL0(L) - CLOSS
            DLL0(L) = CLOSS

            IF(CLOSS>0.) then
               WLQ = WLQ - CLOSS
               RNN(L) = CLOSS 
            else
               RNN(L) = 0.
            ENDIF




         ENDDO


         if (CLDMICRO .gt. 0.0) then !AER_CLOUD=======================================

            CNVNDROP(IC)=NDROP
            CNVNICE(IC)=NICE
            CNVFICE(IC)=FICE
         endif !=======================================



         WLQ = WLQ - QST(IC)*ETA(IC)

         !     CALCULATE GAMMAS AND KERNEL

         GMS(K) =          (SHT(K)-SSL(K))*PRI(K)          ! MS-A30 (W/O GRAV)
         GMH(K) = GMS(K) + (QHT(K)-QOL(K))*PRI(K)*ALHL     ! MS-A31 (W/O GRAV)
         AKM    = GMH(K)*GAM(K-1)*DPB(K-1)                 ! MS-A37 (W/O GRAV)

         TX2     = GMH(K)
         DO L=K-1,IC+1,-1
            GMS(L) = ( ETA(L  )*(SHT(L)-SSL(L  ))   & 
                  + ETA(L+1)*(SSL(L)-SHT(L+1)) )     *PRI(L)
            GMH(L) = GMS(L)                         &
                  + ( ETA(L  )*(QHT(L)-QOL(L  ))   &
                  + ETA(L+1)*(QOL(L)-QHT(L+1)) )*ALHL*PRI(L)
            TX2 = TX2 + (ETA(L) - ETA(L+1)) * GMH(L)
            AKM = AKM - GMS(L)*EHT(L)*PKI(L) + TX2*GHT(L)

         ENDDO

         GMS(IC) = ETA(IC+1)*(SSL(IC)-SHT(IC+1))*PRI(IC)
         AKM     = AKM - GMS(IC)*ETA(IC+1)*DPB(IC)*PKI(IC)

         GMH(IC) =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL &
               + ETA(IC  )*(HST(IC)-HOL(IC  ))     )*PRI(IC)

         if ( SMOOTH_HST ) then
            GMHx    =   GMS(IC) + ( ETA(IC+1)*(QOL(IC)-QHT(IC+1))*ALHL &
                  + ETA(IC  )*(HSTx   -HOL(IC  ))     )*PRI(IC)
         endif

         !    CLOUD BASE MASS FLUX

         IF (AKM >= 0.0 .OR. WLQ < 0.0)  THEN  !  =========>
            RC(IC) = 5
            RETURN
         ENDIF

         WFN = - (WFN-ACR)/AKM ! MS-A39 MASS-FLUX IN Pa/step
         RAS_WFN(I) = WFN ! WMP Store the full mass-flux to use in efficiency diagnostic
         WFN = MIN( TAU(IC)*WFN  ,   (PRS(K+1)-PRS(K) )*(100.*PBLFRAC))

     !! WMP RAS DIAGNOSTICS that make up TAU(IC)
       ! Fill the RAS timescale diagnostic
         RAS_TIME(I) = RASAL(IC)
       ! Fill the RAS RH trigger diagnostic
         RAS_TRG(I) = TRG
       ! Fill the RAS Tokioka diagnostic
         RAS_TOKI(I) = TOKI
       ! Fill the RAS PBL fraction diagnostic
         RAS_PBL(I) = (PRS(K+1)-PRS(K))*(100.*PBLFRAC)
       ! Fill the RAS efficiency diagnostic
       ! RAS_EFFICIENCY(I) = WFN/RAS_EFFICIENCY(I)
     !! WMP RAS DIAGNOSTICS

         !    CUMULATIVE PRECIP AND CLOUD-BASE MASS FLUX FOR OUTPUT

         WFNOG    = WFN*GRAVI
         TEM      = WFN*GRAVI
         CLL (IC) = CLL (IC) + WLQ*TEM           ! (kg/m^2/step)
         RMF (IC) = RMF (IC) +     TEM           ! (kg/m^2/step)
         RMFD(IC) = RMFD(IC) +     TEM * ETA(IC) ! (kg/m^2/step)

         DO L=IC+1,K
            RMFP(L) = TEM * ETA(L)                 ! (kg/m^2/step)
            RMFC(L) = RMFC(L)  +  RMFP(L)          ! (kg/m^2/step)

            DLLX(L) = DLLX(L)  +  TEM*DLL0(L)        


            IF ( CVW(L) > 0.0 ) THEN
               UPDFRP(L) = rmfp(L) * (DDT/DAYLEN) * 1000. / ( CVW(L) * PRS(L) )
            ELSE 
               UPDFRP(L) = 0.0       
            ENDIF

            CLLI(L) = CLL0(L)/ETA(L)  ! current cloud; incloud condensate        
            CLLB(L) = CLLB(L) +  UPDFRP(L) * CLLI(L) !  cumulative grid mean convective condensate        

            UPDFRC(L) =  UPDFRC(L) +  UPDFRP(L)      



         ENDDO



         !    THETA AND Q CHANGE DUE TO CLOUD TYPE IC

         DO L=IC,K
            RNS(L) = RNS(L) + RNN(L)*TEM ! (kg/m^2/step)
            GMH(L) = GMH(L) * WFN
            GMS(L) = GMS(L) * WFN
            QOI(L) = QOI(L) + (GMH(L) - GMS(L)) * ALHI
            POI(L) = POI(L) + GMS(L)*PKI(L)*CPI
            QST(L) = QST(L) + GMS(L)*BET(L)*CPI
         ENDDO

         if ( SMOOTH_HST ) then
            GMHx    = GMHx * WFN
            dQx     = (GMHx - GMH(IC)) * ALHI
            RNS(IC) = RNS(IC) + dQx / ( PRI(IC)*GRAV )
         endif

         IF (DO_TRACERS) THEN
            WFN     = WFN*0.5 *1.0           !*FRICFAC*0.5
            TEM     = WFN*PRI(K)

            DO ITR=1,ITRCR 
               XCU(K,ITR) =  XCU(K,ITR) + TEM * (XOI(K-1,ITR) - XOI(K,ITR))
            END DO

            DO ITR=1,ITRCR 
               DO L=K-1,IC+1,-1
                  TEM    = WFN*PRI(L)
                  ! MEM - note that TEM and ETA are always POSITIVE
                  XCU(L,ITR) = XCU(L,ITR) + TEM *                        &
                        ( (XOI(L-1,ITR) - XOI(L,ITR  )) * ETA(L)   &
                        + (XOI(L,ITR  ) - XOI(L+1,ITR)) * ETA(L+1) )
               ENDDO
            ENDDO

            TEM     = WFN*PRI(IC)

            DO ITR=1,ITRCR 
               XCU(IC,ITR) = XCU(IC,ITR) +  &
                     (2.*(XHT(ITR) - XOI(IC,ITR)*(ETA(IC)-ETA(IC+1))) & 
                     - (XOI(IC,ITR)+XOI(IC+1,ITR))*ETA(IC+1))*TEM
            ENDDO

            DO ITR=1,ITRCR 
               DO L=IC,K
                  XOI(L,ITR) = XOI(L,ITR) + XCU(L,ITR)
               ENDDO
            ENDDO

            ! MEM
            IF (RAS_NO_NEG) call fill_z(K-IC+1, ITRCR, XOI(IC:K,:), DELP(I,IC:K), IC, K )

         else 
            WFN     = WFN*0.5 *1.0           !*FRICFAC*0.5
         endif

         LAMBDSV(IC) = 1.000 


         !   CUMULUS FRICTION

         IF(FRICFAC <= 0.0) THEN
            RC(IC) = 0
            RETURN  !  NO CUMULUS FRICTION =========>>
         ENDIF

         WFN     = WFN*FRICFAC*EXP( -ALM / FRICLAMBDA )
         TEM     = WFN*PRI(K)

         UCU(K)  = UCU(K) + TEM * (UOI(K-1) - UOI(K))
         VCU(K)  = VCU(K) + TEM * (VOI(K-1) - VOI(K))

         DO L=K-1,IC+1,-1
            TEM    = WFN*PRI(L)
            UCU(L) = UCU(L) + TEM *                        &
                  ( (UOI(L-1) - UOI(L  )) * ETA(L  ) &
                  + (UOI(L  ) - UOI(L+1)) * ETA(L+1) )
            VCU(L) = VCU(L) + TEM *                        &
                  ( (VOI(L-1) - VOI(L  )) * ETA(L)   &
                  + (VOI(L  ) - VOI(L+1)) * ETA(L+1) )
         ENDDO

         TEM     = WFN*PRI(IC)

         UCU(IC) = UCU(IC) + (2.*(UHT - UOI(IC)*(ETA(IC)-ETA(IC+1))) & 
               - (UOI(IC)+UOI(IC+1))*ETA(IC+1))*TEM
         VCU(IC) = VCU(IC) + (2.*(VHT - VOI(IC)*(ETA(IC)-ETA(IC+1))) & 
               - (VOI(IC)+VOI(IC+1))*ETA(IC+1))*TEM


         DISSK0(IC) = ETA(IC)* GRAV * WFNOG * PRI(IC) * 0.5 & 
               *( (UHT/ETA(IC)  - UOI(IC) )**2    & 
               + (VHT/ETA(IC)  - VOI(IC) )**2 )  
         DO L=IC,K
            UOI(L) = UOI(L) + UCU(L)
            VOI(L) = VOI(L) + VCU(L)
         ENDDO

         RC(IC) = 0

         RETURN

      END SUBROUTINE CLOUDE



      SUBROUTINE ACRITN( PL, PLB, ACR)

         REAL, INTENT(IN ) :: PL, PLB
         REAL, INTENT(OUT) :: ACR

         INTEGER IWK

         !!REAL, PARAMETER :: FACM=0.5

         REAL, PARAMETER :: &
               PH(15)=(/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, &
               550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/)

         REAL, PARAMETER :: & 
               A(15)=(/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677, &
               0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664, &
               0.0553, 0.0445, 0.0633     /)   !!*FACM

         IWK = PL * 0.02 - 0.999999999

         IF (IWK .GT. 1 .AND. IWK .LE. 15) THEN
            ACR = A(IWK-1) + (PL-PH(IWK-1))*.02*(A(IWK)-A(IWK-1))
         ELSEIF(IWK > 15) THEN
            ACR = A(15)
         ELSE
            ACR = A(1)
         ENDIF

         ACR = ACRITFAC  * ACR * (PLB - PL)

         RETURN

      END SUBROUTINE ACRITN


      SUBROUTINE RNEVP


         ZET(K+1) = 0
         DO L=K,ICMIN,-1
            TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
         ENDDO

         DO L=ICMIN,K

            TEM    = PRI(L)*GRAV
            CNV_PRC3(I,L) = RNS(L)*TEM     
         ENDDO

         !! If hst is smoothed then adjusted precips may be negative
         if ( SMOOTH_HST ) then
            DO L=ICMIN,K
               if ( CNV_PRC3(I,L) < 0. ) then
                  QOI(L)        = QOI(L) +  CNV_PRC3(I,L)
                  POI(L)        = POI(L) -  CNV_PRC3(I,L) * (ALHL/CP) / PRJ(L+1)
                  CNV_PRC3(I,L) = 0.
               endif
            END DO
         endif

         RETURN

      END SUBROUTINE RNEVP



      SUBROUTINE HTEST
         REAL,  DIMENSION(K0) :: HOL1
         integer  :: lminhol
         real     :: minhol


         hol=0.            ! HOL initialized here in order not to confuse Valgrind debugger
         lminhol  = K+1
         MINHOL   = -999999.
         ZET(K+1) = 0
         SHT(K+1) = CP*POI(K)*PRJ(K+1)
         DO L=K,ICMIN,-1
            QOL(L)  = AMIN1(QST(L)*RHMAX,QOI(L))
            QOL(L)  = AMAX1( 0.000, QOL(L) )     ! GUARDIAN/NEG.PREC.
            SSL(L)  = CP*PRJ(L+1)*POI(L) + GRAV*ZET(L+1)
            HOL(L)  = SSL(L) + QOL(L)*ALHL
            HST(L)  = SSL(L) + QST(L)*ALHL
            TEM     = POI(L)*(PRJ(L+1)-PRJ(L))*CPBG
            ZET(L)  = ZET(L+1) + TEM
            ZOL(L)  = ZET(L+1) + (PRJ(L+1)-PRH(L))*POI(L)*CPBG
         ENDDO

         HOL1=HOL
         DO L=K-1,ICMIN+1,-1
            HOL1(L)  = (0.25*HOL(L+1)+0.50*HOL(L)+0.25*HOL(L-1))
            if ( ( MINHOL>=HOL1(L) ) .OR. (MINHOL<0.) ) THEN
               MINHOL   =  HOL1(L)
               LMINHOL  =  L
            end if
         ENDDO

         SIGE_MINHOL = SIGE(LMINHOL)




      end subroutine HTEST


      subroutine FINDDTLS
         real :: SIGDT0,sigmax,sigmin
         integer :: LL

         integer, allocatable :: THE_SEED(:)
         integer :: seed_len

         call random_seed(SIZE=seed_len)
         allocate(THE_SEED(seed_len))

         THE_SEED(1)=SEEDRAS(I,1)*IRAS(I) + SEEDRAS(I,2)*JRAS(I)
         THE_SEED(2)=SEEDRAS(I,1)*JRAS(I) + SEEDRAS(I,2)*IRAS(I)
         THE_SEED(1)=THE_SEED(1)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
         THE_SEED(2)=THE_SEED(2)*SEEDRAS(I,1)/( SEEDRAS(I,2) + 10)
         if(THE_SEED(1) == 0) THE_SEED(1) =  5
         if(THE_SEED(2) == 0) THE_SEED(2) = -5

         ! Gfortran uses longer seeds, so fill the rest with zero
         if (seed_len > 2) THE_SEED(3:) = 0

         call random_seed(PUT=THE_SEED)

         SIGMAX=SIGE(K)
         SIGMIN=SIGE(ICMIN)

         if ( RASNCL < 0.0 ) then 
            !! NO SHALLOW CONV   N_DTL = 56 - ICMIN 
            N_DTL = K - ICMIN 
         else 
            N_DTL = min( int( RASNCL ) , K-ICMIN )
         endif

         if(allocated(ICL_V)) deallocate(ICL_V)
         allocate(ICL_V(N_DTL))

         if ( ( RASNCL < 0.0 ) .and. ( RASNCL >=-100.) ) then 
            do L=1,N_DTL
               ICL_V(L) = ICMIN + L - 1
            enddo
         else if ( RASNCL < -100.0 ) then 
            do L=1,N_DTL
               ICL_V(L) = K - L
               !! NO SHALLOW CONV           ICL_V(L) = 56 - L
            enddo
         else
            do L=1,N_DTL
               call random_number ( SIGDT0 )
               SIGDT0 = 1.00 - ( SIGDT0**RDTLEXPON )
               SIGDT0 = SIGMIN+SIGDT0*(SIGMAX-SIGMIN)

               do LL=ICMIN,K
                  if ( (SIGE(LL+1)>=SIGDT0) .and. (SIGE(LL)<SIGDT0 ) ) ICL_V(L)=LL
               enddo
            end do
         endif

         deallocate(THE_SEED)

      end subroutine FINDDTLS




      SUBROUTINE STRAP(FINAL)

! MEM - these are inherited from the current scope:
!   K = KCBL  e.g. level of PBL
!   icmin = extreme level of conv detrainment ?  30 hPa
!   PRJ = PKE for the column

!   POI = THO for the column
!   QOI = QHO for the column
!   UOI = UHO for the column
!   VOI = VHO for the column

! MEM - added these:
!   XFF(L,ITR) = for range of layers KCBL to surface, record the weighted value / range total 
!   WW(L)  = level by level weighting, for unpacking

         INTEGER :: FINAL
         REAL , DIMENSION(K0)  :: WGHT, MASSF, WW

         REAL :: WGHT0, PRCBL

         integer, parameter :: nrands=1
         real ::    rndu(nrands)
         integer :: seedcbl(nrands)

         ! !DESCRIPTION: 
         !   {\tt STRAP} is called: FINAL=0, to compute cloud base layer CBL properties
         !   given a value K for the index of the upper {\em EDGE} of the CBL; FINAL=1
         !   to redistribute convective tendencies within CBL

         integer :: KK

         !  LOCAL VARIABLES FOR USE IN CLOUDE

         !!IF (.NOT. PRESENT(FINAL)) THEN
         IF (FINAL==0) THEN


            !!PRJ(ICMIN:K+1) = PKE(I,ICMIN:K+1)


            do kk=icmin,k+1
               PRJ(kk) = PKE(I,kk)
            enddo


            poi=0.        ! These initialized here in order not to confuse Valgrind debugger
            qoi=0.        ! Do not believe it actually makes any difference.
            uoi=0.
            voi=0.     

            PRS(ICMIN:K0+1) = PLE(I,ICMIN:K0+1)
            POI(ICMIN:K)   = THO(I,ICMIN:K)
            QOI(ICMIN:K)   = QHO(I,ICMIN:K)
            UOI(ICMIN:K)   = UHO(I,ICMIN:K)
            VOI(ICMIN:K)   = VHO(I,ICMIN:K)

            WSP(ICMIN:K)   = SQRT(   ( UOI(ICMIN:K) - UOI(K) )**2  & 
                  + ( VOI(ICMIN:K) - VOI(K) )**2    ) 

            QST(ICMIN:K) = QSS(I,ICMIN:K)
            DQQ(ICMIN:K) = DQS(I,ICMIN:K)

            IF (DO_TRACERS) THEN 
               DO ITR=1,ITRCR 
                  XOI(ICMIN:K,ITR) = XHO(I,ICMIN:K,ITR)  ! Init the column from 30 hPa down to KCBL
               END DO
            END IF


!!! Mass fraction of each layer below cloud base
!!! contributed to aggregate cloudbase layer (CBL) 
            MASSF(:) = WGT0(I,:)


!!! RESET PRESSURE at bottom edge of CBL 
            PRCBL = PRS(K)
            do l= K,K0
               PRCBL = PRCBL + MASSF(l)*( PRS(l+1)-PRS(l) )
            end do
            PRS(K+1) = PRCBL
            PRJ(K+1) = (PRS(K+1)/1000.)**(MAPL_RGAS/MAPL_CP)

            DO L=K,ICMIN,-1
               POL(L)  = 0.5*(PRS(L)+PRS(L+1))
               PRH(L)  = (PRS(L+1)*PRJ(L+1)-PRS(L)*PRJ(L)) &
                     / (ONEPKAP*(PRS(L+1)-PRS(L)))
               PKI(L)  = 1.0 / PRH(L)
               DPT(L)  = PRH(L  ) - PRJ(L)
               DPB(L)  = PRJ(L+1) - PRH(L)
               PRI(L)  = .01 / (PRS(L+1)-PRS(L))
            ENDDO


!!!!! RECALCULATE PROFILE QUAN. IN LOWEST STRAPPED LAYER
            if( K<=K0) then
               POI(K) = 0.
               QOI(K) = 0.
               UOI(K) = 0.
               VOI(K) = 0.

               !! SPECIFY WEIGHTS GIVEN TO EACH LAYER WITHIN SUBCLOUD "SUPERLAYER"
               WGHT = 0.
               DO L=K,K0
                  WGHT(L)   = MASSF(L) *                &
                        ( PLE(I,L+1) - PLE(I,L) ) &
                        /( PRS(K+1)   - PRS(K)  )
               END DO



               DO L=K,K0
                  POI(K) = POI(K) + WGHT(L)*THO(I,L)
                  QOI(K) = QOI(K) + WGHT(L)*QHO(I,L)
                  UOI(K) = UOI(K) + WGHT(L)*UHO(I,L)
                  VOI(K) = VOI(K) + WGHT(L)*VHO(I,L)
               ENDDO


               IF (DO_TRACERS) THEN 
                  IF (RAS_NO_NEG) THEN 
                     XOI(K,:)=0.                                           ! Init for accumulation
                     XFF(:,:)=0.
                     DO ITR=1,ITRCR
                        DO L=K,K0                                          ! From KCBL down to surface
                           XFF(L,ITR) =              WGHT(L)*XHO(I,L,ITR)  ! Record weighted tracer value
                           XOI(K,ITR) = XOI(K,ITR) + WGHT(L)*XHO(I,L,ITR)  ! Accumulate values in KCBL
                        END DO

                        IF ( XOI(K,ITR) .LT. 1.0e-25 ) THEN
!                         Cannot divide by the very small total
                          XFF(K:K0,ITR) = 1.0 / (K0-K+1)                     ! Divide equally among levels
                        ELSE
                          DO L=K,K0                                          ! From KCBL down to surface
                             XFF(L,ITR) = XFF(L,ITR) / XOI(K,ITR)            ! Divide weighted tracers by total
                          END DO
                        END IF
                     END DO
                     ! MEM - at this point, there are no negatives in XOI
                  ELSE
                     XOI(K,:)=0.
                     DO ITR=1,ITRCR
                        DO L=K,K0
                           XOI(K,ITR) = XOI(K,ITR) + WGHT(L)*XHO(I,L,ITR)
                        END DO
                     END DO
                  END IF
               END IF

               DQQ(K) = DQSAT( POI(K)*PRH(K) , POL(K), qsat=QST(K) )

            endif


            !!DPTH_BL = CPBG*POI(K)*( PRJ(K+1)-PRJ(K) )

            ! seedras(1,2) are both integers passed from
            ! from GEOS_Moist w/ values 0 - 1000000
            ! rndu(:) = 1.0*( seedras(1)+seedras(2) )/2000000.
            rndu(:) = max( seedras(I,1)/1000000., 1e-6 )

            !!call congvec( npoints , seedcbl , rndu )
            DPTH_BL   = ZCBL(I)

            MXDIAM(I) = CNV_FRACTION(I)*ABS(MAXDALLOWED_D) + (1-CNV_FRACTION(I))*ABS(MAXDALLOWED_S)
            if (MAXDALLOWED_E < 0) then
               ! Make MXDIAM stochastic
               MXDIAM(I) = MXDIAM(I)*( rndu(1)**(MAXDALLOWED_E) )
            endif

            DO L=K,ICMIN,-1
               BET(L)  = DQQ(L)*PKI(L)  !*
               GAM(L)  = PKI(L)/(1.0+LBCP*DQQ(L)) !*
               IF(L<K) THEN
                  GHT(L+1) = GAM(L)*DPB(L) + GAM(L+1)*DPT(L+1)
                  GM1(L+1) = 0.5*LBCP*(DQQ(L  )/(ALHL*(1.0+LBCP*DQQ(L  ))) + &
                        DQQ(L+1)/(ALHL*(1.0+LBCP*DQQ(L+1))) )
               ENDIF
            ENDDO


            TCU(ICMIN:K) = -POI(ICMIN:K)*PRH(ICMIN:K)
            QCU(ICMIN:K) = -QOI(ICMIN:K)

            RNS  = 0.
            CLL  = 0.
            RMF  = 0.
            RMFD = 0.
            RMFC = 0.
            RMFP = 0.
            CLL0 = 0.
            DLL0 = 0.
            CLLX = 0.
            DLLX = 0.
            CLLI = 0.
            CLLB = 0.

            POI_SV = POI
            QOI_SV = QOI
            UOI_SV = UOI
            VOI_SV = VOI
            IF (DO_TRACERS) XOI_SV = XOI

            LAMBDSV = 0.0
            CVW     = 0.0
            UPDFRC  = 0.0
            UPDFRP  = 0.0
            DISSK0  = 0.0

         END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!    IF (PRESENT(FINAL)) THEN


         IF (FINAL==1) THEN 


            THO(I,ICMIN:K-1) = POI(ICMIN:K-1)
            QHO(I,ICMIN:K-1) = QOI(ICMIN:K-1)
            UHO(I,ICMIN:K-1) = UOI(ICMIN:K-1)
            VHO(I,ICMIN:K-1) = VOI(ICMIN:K-1)
            CNV_UPDFRC(I,ICMIN:K-1)   =  UPDFRC(ICMIN:K-1)
            CNV_CVW   (I,ICMIN:K-1)   =     CVW(ICMIN:K-1)
            CNV_QC(I,ICMIN:K-1)       =  CLLB(ICMIN:K-1)

            if (CLDMICRO .gt. 0.0) then!======================AER_CLOUD=============
               CNV_NDROP   (I,ICMIN:K-1)  =    CNVNDROP(ICMIN:K-1) !DONIF
               CNV_NICE   (I,ICMIN:K-1)   =     CNVNICE(ICMIN:K-1) !DONIF
               CNV_FICE   (I,ICMIN:K-1)   =     CNVFICE(ICMIN:K-1) !DONIF
            endif
            !! De-strap tendencies from RAS
            !! specify weighting "SHAPE"
            WGHT   = WGT1(I,:)

            !! Scale properly by layer masses
            wght0 = 0.
            WW(:) = 0.
            DO L=K,K0 
               wght0 = wght0 +              WGHT(L)* ( PLE(I,L+1) - PLE(I,L) )
               WW(L) = (PRS(K+1) - PRS(K))/(WGHT(L)* ( PLE(I,L+1) - PLE(I,L) ))
            END DO

            wght0 = ( PRS(K+1)   - PRS(K)  )/wght0

            WGHT  = wght0 * WGHT


            DO L=K,K0 
               THO(I,L) =  THO(I,L) + WGHT(L)*(POI(K) - POI_SV(K))
               QHO(I,L) =  QHO(I,L) + WGHT(L)*(QOI(K) - QOI_SV(K))
               UHO(I,L) =  UHO(I,L) + WGHT(L)*(UOI(K) - UOI_SV(K))
               VHO(I,L) =  VHO(I,L) + WGHT(L)*(VOI(K) - VOI_SV(K))
            END DO

            IF (DO_TRACERS) THEN 
               XHO(I,ICMIN:K-1,:) = XOI(ICMIN:K-1,:) 
               IF ( RAS_NO_NEG ) THEN
!                 MEM - proportionally distributed:
                  DO ITR=1,ITRCR
                     DO L=K,K0                       !  For levels KCBL down to the surface
                        XHO(I,L,ITR) =  XFF(L,ITR) * XOI(K,ITR) * WW(L)
                     END DO
                  END DO
               ELSE
!                 Previous approach:
                  DO ITR=1,ITRCR
                     DO L=K,K0                       !  For levels KCBL down to the surface
                        XHO(I,L,ITR) =  XHO(I,L,ITR) + WGHT(L)*(XOI(K,ITR) - XOI_SV(K,ITR))
                     END DO
                  END DO
               END IF
            END IF




            FLX (I,ICMIN:K) = RMF (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD BASE)
            FLXD(I,ICMIN:K) = RMFD(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
            FLXC(I,ICMIN:K) = RMFC(ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s @ CLOUD TOP)
            CLW (I,ICMIN:K) = CLL (ICMIN:K) * DDT/DAYLEN  !  (KG/m^2/s )

            if ( PRESENT( DISSKE ))  DISSKE(I,ICMIN:K-1) = DISSK0(ICMIN:K-1) * DDT/DAYLEN

            FLX (I,1:ICMIN-1) = 0.
            FLXD(I,1:ICMIN-1) = 0.
            FLXC(I,1:ICMIN-1) = 0.
            CLW (I,1:ICMIN-1) = 0.


            IF ( K < K0 ) THEN 
               FLX (I,K:K0) = 0.
               FLXD(I,K:K0) = 0.
               FLXC(I,K:K0) = 0.
               CLW (I,K:K0) = 0.
            END IF

            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)

            IF(ALLOCATED(ICL_V)) DEALLOCATE( ICL_V )

         ENDIF

         IF (FINAL==2) THEN 


            FLX (I,:) = 0.
            FLXD(I,:) = 0.
            FLXC(I,:) = 0.
            CLW (I,:) = 0.

            IRC (I,ICMIN:K-1) = RC(ICMIN:K-1)

         ENDIF

         RETURN

      END SUBROUTINE STRAP


!!!!!!!!!!!!!!======================================     
!!!!!!!!!   Subroutine ARGact: finds the activated droplet number from Abdul_Razzak and Ghan 2000.
      ! Tailored for GOCART AEROSOL and works with AERO input 
      !Written by Donifan Barahona
      !!donifan.o.barahona@nasa.gov
!!!!!!!!!!!!!!====================================

      SUBROUTINE ARGact (TEMP, WX, NCPL_ACT, NCPL_AMB,  CDUST, CSOOT, LEV, ISBASE, DDUSTAMB, DSOOTAMB, ENT_PARAM)
         !
         Integer, intent(in)     ::  LEV    
         LOGICAL,  intent(in)     ::   ISBASE

         REAL, intent(inout)     ::   TEMP, WX, ENT_PARAM
         REAL, intent(out)     ::   NCPL_ACT, NCPL_AMB, CSOOT, DSOOTAMB
         REAL, DIMENSION(NDUSTMAX), INTENT(OUT) :: CDUST, DDUSTAMB
         integer                 :: INDEX, NMODES, naux       

         type(AerProps)  :: AER, auxaer      


         real     ::      kappa, alfa, beta, Akoh, G, T, smax, fi, gi, nui, &
                       citai, ui, aux1, PACT,  Ntot, auxx, aux, auxconc, W, alph, aseasalt_aux, f_seasalt1
         real, dimension (30) ::  SMI, TPI, SIGI 


         SMI=0.0      
         TPI = 0.0
         SIGI =2.0
         NCPL_ACT=0.0
         NCPL_AMB=0.0
         CDUST=0.0
         CSOOT=0.0
         DDUSTAMB =1.0e-9
         DSOOTAMB= 1.0e-9
         W=MIN(WX*(1.0-ENT_PARAM), 20.0)    

         PACT=0.0 !activation probability of entrained aerosol 

         AER%nmods = 0
         AER%num   = 0.0

         do INDEX = 1, AERO(LEV)%nmods
             if (AERO(LEV)%num(INDEX) > 0.1) then
                 AER%nmods = AER%nmods + 1
                 naux = AER%nmods

                 AER%num(naux)   = AERO(LEV)%num(INDEX)
                 AER%dpg(naux)   = max(AERO(LEV)%dpg(INDEX), 1.0e-9)
                 AER%sig(naux)   = AERO(LEV)%sig(INDEX)
                 AER%den(naux)   = AERO(LEV)%den(INDEX)
                 AER%kap(naux)   = AERO(LEV)%kap(INDEX)
                 AER%fdust(naux) = AERO(LEV)%fdust(INDEX)
                 AER%fsoot(naux) = AERO(LEV)%fsoot(INDEX)
                 AER%forg(naux)  = AERO(LEV)%forg(INDEX)
             end if
         end do

         !if (AER%nmods == 0) then
         !    call init_Aer(AER)
         !end if

!!!!!!!!!!activate aerosol transported from cloud base 
             NMODES =  AER_BASE%nmods
             TPI(1:nmodes) = AER_BASE%num(1:nmodes)
             SIGI(1:nmodes) = AER_BASE%sig(1:nmodes)                          

             
             Ntot= 0.0
              do index = 1, nmodes 
	              if (AER_BASE%kap(index) .gt. 0.1) Ntot =  Ntot + TPI(index)  
              end do
         

         if ((Ntot .lt. 1.0e4) .or. (TEMP .lt. 245.0) .or. (W .lt. 0.01)) then !no activation if aerosol < 1e-4 1/cm3            
            NCPL_ACT  = 0.0
         else

            ! Calculate constants. These fits were obtained from detailed correlations of physical properties. G is actually 1/G
            T = min(max(TEMP, 243.0), 323.0)     
            alfa=2.8915E-08*(T**2) - 2.1328E-05*T + 4.2523E-03
            beta=exp(3.49996E-04*T**2 - 2.27938E-01*T + 4.20901E+01)
            G=exp(-2.94362E-06*T**3 + 2.77941E-03*T**2 - 8.92889E-01*T + 1.18787E+02)
            Akoh= 0.66e-6/T  !from Seinfeld and Pandis (1998)
     
            !=======================================================
            !Activate droplets   
            !=======================================================
            !Calculate maximum supersaturation according to ARG2002

            auxx=0.0 
            
            
            DO INDEX = 1, NMODES            
                
               kappa=  max(AER_BASE%kap(INDEX), 0.001)
             
                  SMI (INDEX) = ((0.667*Akoh/AER_BASE%dpg(INDEX))**1.5)/SQRT(2.0*kappa)   ! Critical supersat for mode I   
                  SMI=MAX(SMI, 1.0e-5)   
                   
              if ((TPI(INDEX) .gt. 1e4) .and.  (kappa .gt. 0.1)) then                       
                  fi=0.5*exp(2.5*SIGI(INDEX)) !sigi is now log(sigi)
                  gi=1.0+0.25*SIGI(INDEX)
                  nui=((alfa*W*G)**1.5)/(2.0*MAPL_PI*980.0*beta*TPI(INDEX))
                  citai = 0.667*Akoh*SQRT(alfa*W*G)
                  aux1=fi*((citai/nui)**1.5) + gi*(SMI(INDEX)*SMI(INDEX)/(nui+(3.0*citai)))**0.75
                  aux1=aux1/(SMI(INDEX)*SMI(INDEX))      
                  auxx=auxx+aux1                  
                end if
            end do

  !Calculate number of activated droplets
            if (auxx .gt. 0.0) then
               smax = 1/sqrt(auxx)
               auxx=0.0

                   DO INDEX = 1, NMODES
                        if ((TPI(INDEX) .gt. 1e4) .and. (AER_BASE%kap(index) .gt. 0.1)) then
                           ui=sqrt(2.0)*log(SMI(INDEX)/smax)/3.0
                           aux1=0.5*TPI(INDEX)*(1.0-ERFAPP(ui))
                           auxx=auxx+aux1
                           AER_BASE%num(index) = max(TPI(INDEX) -aux1, 0.0) !remove already activated aerosol
                        end if                    
                   END DO
                  NCPL_ACT=auxx             
            else
                  NCPL_ACT = 0.0             
            end if

     
 
       end if

         !now filllup dust and soot number
         NMODES =  AER%nmods

         call getINsubset(1, AER,  auxaer)
         CDUST(1:auxaer%nmods)= auxaer%num(1:auxaer%nmods)
         DDUSTAMB(1:auxaer%nmods)= auxaer%dpg(1:auxaer%nmods)
         call getINsubset(2, AER,  auxaer)
         naux = max(auxaer%nmods, 1)
         CSOOT= sum(auxaer%num) 
         DSOOTAMB= sum(auxaer%dpg)/naux


         PACT=1.0 ! fraction of entrained aerosol that is activated
         auxconc =0.0
         aseasalt_aux  = 0.0

         do index = 1, nmodes 
	       if (AER%kap(index) .gt. 0.8)  auxconc = AER%num(index) + auxconc
           if (AER_BASE%kap(index) .gt. 0.8)    aseasalt_aux  = aseasalt_aux  + AER_BASE%num(index)*AER_BASE%dpg(index)*AER_BASE%dpg(index)*1.61*MAPL_PI !assumes a fixed sigma = 2.0

         end do
       aseasalt = max(aseasalt, aseasalt_aux)
     
	  NCPL_AMB=auxconc !Activate  entrained aerosol with kappa>0.8   



      END SUBROUTINE ARGact

!!!!!!!!!!!!!!!!!!!!!!!!!!
      !=====Subroutine INfreezing=========
      ! Freeze droplets in immersion and contact ice nucleation modes, according to Barahona et al. GMD (2014)
!!!!!!!!!!!!!!!!!!!!!

      subroutine INfreezing(QL, NL, NIN, NDUST, NSOOT, INDUST, INSOOT, TEMP, PRE, WSUB_, DDUST, DSOOT)  !NIN freezing tendency
         REAL, INTENT( IN) :: TEMP, NSOOT, WSUB_, PRE, QL, NL, DSOOT
         REAL, DIMENSION(NDUSTMAX), INTENT (IN) ::  NDUST, DDUST
         REAL, DIMENSION(NDUSTMAX), INTENT (OUT) ::  INDUST
         REAL, INTENT(OUT) ::  NIN, INSOOT 
         REAL, DIMENSION(NDUSTMAX) ::  INDUST_C

         real :: a, b, c , d, Tx, n05, ui, aux, nssoot, nsdust, ninbkg, SI, acorr, &
               dnsd, dnss, coolr, WSUB, ahet, INSOOT_C

         real :: nssoot_c, nsdust_c, mfp, nslip, ndfaer, viscosity, lam_r, taux, rho_a, &
                fdust_drop, fsoot_drop, min_ns_dust, min_ns_soot, nsss, INsea, dnsss, min_ns_seasalt 

         logical :: demott, Drop_mediated
         integer :: ix

         min_ns_dust= 3.75e6 !limits ice nucleation to -12 !new 02/10/14 
         min_ns_soot= 3.75e9 !limits ice nucleation to -18
         min_ns_seasalt = 4.0e2 !limits ice nucleation to -5
         
         demott=.false.
         INDUST=0.0
         INSOOT=0.0
         INDUST_C=0.0
         INSOOT_C=0.0   
         NIN=0.0
         Drop_mediated = .false.
       INsea = 0.0 ! sea salt only in immersion
   
   
   ! note for sea salt we just assume that it is equal to the current droplet concentration and take the area from the calculation at cloud base

         ! fraction of dust and soot within droplets
         fdust_drop= FDROP_DUST
         fsoot_drop = FDROP_SOOT


         WSUB=MAX(MIN(WSUB_, 10.0), 0.8)
         coolr=5.0e-3*WSUB  !Approximation to saturated cooling rate 
         n05=sum(NDUST)+NSOOT

         if (TEMP .gt. T_ICE_MAX) then 
            return
         end if

         if (TEMP .lt. T_ICE_ALL   ) then 
            return
         end if


         if ((QL .le. 1.0e-10) .or. (NL .le. 1.0)) then 
            return
         end if

         !Background IN
         ! SI at water saturation

         rho_a = PRE*100.0/MAPL_RGAS/TEMP 
         ninbkg=0.0
         SI = -1.2379e-2+3.3595 !Ice supersat at water sat. Derived from Murphy and Koop 2005
         !if (TEMP .lt. 260.0)  ninbkg=coolr*42.8*exp(3.88*si)*0.1 !tendency in IN from background IN. Derived from Phillips et al. 2007


         Tx = max(TEMP-273.16, -38.0 )

         lam_r=min((MAPL_PI*950.0*NL/rho_a/QL)**(1./3.), 1.0e8)
         viscosity=1.8e-5*(TEMP/298.0)**0.85    ! Viscosity (kg/m/s)
         mfp=2.0*viscosity/(PRE  &                   ! Mean free path (m)
               *sqrt(8.0*28.96e-3/(MAPL_PI*8.314409*TEMP)))        

         if ((n05 .gt.1.0) .and. (TEMP .lt. 272.0)) then

            nsdust=  max(exp(-0.517*Tx + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
            nssoot= max(1.0e4*exp(-0.0101*Tx*Tx - 0.8525*Tx + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)
            dnsd  = 0.517*nsdust
            dnss  = max(-(-2.0*0.0101*Tx -0.8525)*nssoot, 0.0)

            !ns in  contact. It is assumed that in contact is T-3 immersion
            taux=max(Tx-3.0, -35.0)
            nsdust_c= max(exp(-0.517*taux + 8.934)-min_ns_dust, 0.0) !From Niemand 2012
            nssoot_c= max(1.0e4*exp(-0.0101*taux*taux - 0.8525*taux + 0.7667)-min_ns_soot, 0.0) !Murray (review_ 2012)

            aux=0.0        
            acorr=2.7e7 !m2/m3 correction to the area due to non sphericity and aggregation. Assumes 10 m2/g (Murray 2011)


            DO ix=1, NDUSTMAX
               !Immersion
               ahet=0.52*DDUST(ix)*DDUST(ix)*DDUST(ix)*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))    !this needs to be improved

               INDUST(ix) = NDUST(ix)*exp(-nsdust*ahet)* &
                     dnsd*coolr*ahet*fdust_drop
               !Contact   
               nslip =1.0+(2.0*mfp/DDUST(ix))*(1.257+(0.4*exp(-(1.1*DDUST(ix)*0.5/mfp))))! Slip correction factor                    
               ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nsdust_c*ahet)) /(12.*MAPL_PI*viscosity*DDUST(ix))             
               INDUST_C(ix) = 2.0*MAPL_PI*ndfaer*NDUST(ix)*NL/lam_r

            END DO


            acorr=8.0e7 !m2/m3 correction to the area due to non sphericity and aggregation  Assumes 50 m2/g (Popovicheva 2003)       
            ahet =0.52*DSOOT*DSOOT*DSOOT*acorr*exp(4.5*log(2.0)*log(2.0)*log(2.0))             
            INSOOT=fsoot_drop*NSOOT*exp(-nssoot*ahet)*dnss*ahet*coolr !

            nslip =1.0+(2.0*mfp/DSOOT)*(1.257+(0.4*exp(-(1.1*DSOOT*0.5/mfp))))! Slip correction factor                    
            ndfaer =1.381e-23*TEMP*nslip*(1.0-exp(-nssoot_c*ahet)) /(12.*MAPL_PI*viscosity*DSOOT)             
            INSOOT_c= 2.0*MAPL_PI*ndfaer*NSOOT*NL/lam_r

         ! sea salt
         nsss =  -0.459*TEMP +128.6235 ! from Demott et al. PNAS, 2015  
         nsss=  max(exp(nsss)-min_ns_seasalt, 0.0)           
    	 dnsss=  max(0.459*nsss, 0.0)
         INsea= aseasalt*dnsss*coolr*FDROP_SEASALT

 
         end if

	    NIN =ninbkg+ INSOOT + SUM(INDUST) + INSOOT_C + SUM(INDUST_C) + INsea!
         INSOOT=INSOOT +INSOOT_C
         INDUST =INDUST + INDUST_C

	    
      end subroutine INfreezing




!!!!!!!!!!!!!!!!!!!!!!!!!!
      !===================Subroutine Qgrowth======================
      !Partitions water and ice according to Korolev and Mazin 2005. Assume spheres by now.
!!!!!!!!!!!!!!!!!!!!
      subroutine Qgrowth(TEMP, PRE, QICE0, NIPRE, QT, NINUC, DQIG, RIM, FNDRIM) !freezing of IN according to Demott et al 2010 (everything SI)
         REAL, INTENT(IN) :: TEMP, QICE0, NIPRE, NINUC,  PRE, QT
         REAL, INTENT(INOUT) :: DQIG, RIM, FNDRIM

         !real :: A, denI, aux, Dco, SI, denA, DQold, DQnew 
         real :: DIFF, DENICE, DENAIR, K1, K2, K3, SI, AUX, DC, TEFF, PLo, TEo, TC, &
               DQnew, DQold, rho_a, LWC, IWC, qmin_rim


         if (TEMP .gt. 272.15) then
            DQIG =0.0
            return
         end if

         TC=TEMP-273.0 
         PLo = max(PRE, 10.0) !limits  of the correlations 
         TEo = max(190.0, TEMP)

         qmin_rim = 1.0e-12


         DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
         DENAIR= PLo*100.0/MAPL_RGAS/TEMP
         DIFF=(0.211*1013.25/(PLo+0.1))*(((TEo+0.1)/273.0)**1.94)*1e-4  !From Seinfeld and Pandis 2006

         K1 = EXP(7.1170e-4*TEo*TEo-0.43563*TEo+78.744) 
         K2 = EXP(-9.933e-3*TEo+25.26)
         K3 = EXP(7.1772e-4*TEo*TEo-0.44055*TEo+73.996)


         AUX= 210368.0 + 131.438*TEMP - (3.32373E6/TEMP)- (41729.1*LOG(TEMP)) !From Murphy and Koop 2005
         SI=exp(-aux/8.314/TEMP)-1.0 !ratio of pw/pi-1
         rho_a = PRE*100.0/MAPL_RGAS/TEMP 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if  ((NIPRE .gt. 1.0) .and. (QICE0 .gt. 1.0e-10)) then 
            DC=max((QICE0/(NIPRE*500.0*MAPL_PI))**(0.333), 40.0e-6) !Assumme monodisperse size distribution about size distribution is made. 
         else   
            DC = 40.0e-6
         end if

         AUX=  NIPRE*DENICE*MAPL_PI*DC*DC
         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))

         if  (AUX .gt. 1.0e-12) then 
            TEFF=min(TEFF/AUX, 1.0e20)   
            DQold= SI/TEFF  
         else
            DQold=0.0
         end if

         ! Calculate rimming fraction
         IWC =  QICE0*rho_a
         LWC =  max(QT-QICE0, 0.0)*rho_a
         aux = DQold ! only due to deposition


         !Account for rimming

         if ((LWC .gt. qmin_rim)  .and. (IWC .gt. qmin_rim))  then 
            RIM = 6.0e-5/(LWC*(IWC**0.17)) !Fom Lin and Colle, NRW, 2011
            RIM = 1.0/(1.0+RIM)
            RIM  = min (0.95, RIM)
            DQold =  DQold*(1 + RIM/(1.0-RIM))
            FNDrim =  max(min(rho_a*(DQold -aux)/LWC, 1.0), 0.0) !Fraction of liquid condensate removed due to riming   
         END if




!!!!!!!!!!!!!!!!!!!!!recently nucleated!!!!!!!!!!!!!!!!!!!!!!!!


         AUX=  NINUC*DENICE*MAPL_PI*20.0e-6*20.0e-6
         TEFF = DENAIR*2.0*((K1*DIFF+K2)*DC+(K3/0.1))

         if  (AUX .gt. 1.0e-12)  then 
            TEFF=min(TEFF/AUX, 1.0e10)
            DQnew= SI/TEFF
         else
            DQnew = 0.0
         end if

         DQIG = DQold+DQnew     

      end subroutine Qgrowth



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !*************************************************************
      ! Function PDG07 (simplified background ice nucleation 
      !                     spectra according to Phillips et. al. 2007).  
      ! si is supersaturation wrt ice and T is in K 
      !************************************************************  

      subroutine PDG07_ice(si, Tx, N)     

         real, intent(in) :: si, Tx
         real, intent(out)  :: N 
         N=0.0

         !if (Tx .le. 243.0)then


         N=1000.0*exp(-0.388)*(exp(3.88*si)-1.0)/0.76
         !elseif (Tx .le. 260.0) then
         !  N=60.0*exp(-0.639)*(exp(12.96*si)-1.0)/0.76   
         !end if      
      end subroutine PDG07_ice


      !*********************************************************************



   END SUBROUTINE RASE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE SUNDQ3_ICE( TEMP,RATE2,RATE3,TE1, F2, F3)

      REAL, INTENT( IN) :: TEMP,RATE2,RATE3,TE1
      REAL, INTENT(OUT) :: F2, F3

      REAL :: XX, YY,TE0,TE2,JUMP1  !,RATE2,RATE3,TE1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  Ice - phase treatment totally invented
      !!  Sharp increase in autoconversion in range
      !!  ~~TE1 K ~< T < TE0 K .
      !!  (JTB, 3/25/2003)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      TE0=273.
      TE2=200.
      JUMP1=  (RATE2-1.0) / ( ( TE0-TE1 )**0.333 ) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Ice - phase treatment  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      IF   ( TEMP .GE. TE0 ) THEN 
         F2   = 1.0
         F3   = 1.0
      ENDIF
      IF ( ( TEMP .GE. TE1 ) .AND. ( TEMP .LT. TE0 ) )THEN 
         F2   = 1.0 + JUMP1 * (( TE0 - TEMP )**0.3333)
         F3   = 1.0
      ENDIF
      IF   ( TEMP .LT. TE1 ) THEN 
         F2   = RATE2 + (RATE3-RATE2)*(TE1-TEMP)/(TE1-TE2)
         F3   = 1.0
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF ( F2 .GT. 27.0 ) F2=27.0


   end  subroutine sundq3_ice

   subroutine congvec( npoints , seed, ran )

      INTEGER, intent(in   ) :: npoints
      INTEGER, intent(inout) :: seed(npoints)
      REAL,    intent(inout) :: ran(npoints)

      INTEGER irand,i2_16,overflow_32  ! variables for RNG
      INTEGER, PARAMETER :: huge32=2147483647
      i2_16=65536
      !i2_16=111

      ! Num. Recipes
      do irand = 1, npoints
         ! Marsaglia CONG algorithm
         seed(irand)=1664525*seed(irand)+1013904223
         ! mod 32 bit overflow
         seed(irand)=mod(seed(irand),2**32)
         ran(irand)=seed(irand)*0.931322574615479E-09
      enddo

#if 0
      do irand = 1, npoints
         ! Marsaglia CONG algorithm
         seed(irand)=69069*seed(irand)+1234567
         ! mod 32 bit overflow
         seed(irand)=mod(seed(irand),2**30)
         ran(irand)=seed(irand)*0.931322574615479E-09
      enddo

      ! convert to range 0-1 (32 bit only)
      overflow_32=i2_16*i2_16
      if ( overflow_32 .le. huge32 ) then
         do irand = 1, npoints
            ran(irand)=ran(irand)+1
            ran(irand)=(ran(irand))-int(ran(irand))
         enddo
      endif
#endif


   end subroutine congvec

   !==Subroutine Qremoval ============
   !Physically based parameterization for removal of condensate (Barahona et al. GMD 2014)
   !Written by Donifan Barahona according to DelGenio et al. 2005. J Climate with some modifications
   !donifan.barahona@nasa.gov 
   !================================================================

   subroutine Qremoval(RATE_Q, FICE_, F_NL, F_NI, TE, PRE, QC, NI, &
         NL, VCNx, FPICE_, DTL, RIM, COAUTO) 

      !This calculations are needed since we do not have explicit microphysics in RAS. All Fs are tendencies

      REAL, INTENT( IN) ::  FICE_,  TE, PRE, QC, VCNx, NI, NL, DTL, COAUTO
      REAL, INTENT (INOUT) :: RIM
      REAL, INTENT(OUT) :: RATE_Q, F_NL, F_NI, FPICE_

      REAL :: Dc_r, Dc_g, Dc_i, lamr, lamg, lami, fqi, &
            fqr, fqg, fpi, fpr, fpg, No, rho_r, rho_g, rho_i, &
            qr, qg, qi, fqig, rho_a, psc, QPT, aux, VCN, NG, fni, fnr, fng, NI_, &
            beta6, xs, prc, SI, Ai, Ah, Bh, Tc, L, N

      real  :: mui


      VCN=1.0
      VCN=max(min(VCNx, 10.0), 0.8) 

      Ai=exp(7.48416E-04*TE*TE - 4.38424E-01*TE + 8.62639E+01) !From detailed calculations
      Ai=1.0/Ai   
      aux= 210368.0 + 131.438*TE - (3.32373E6/TE)- (41729.1*LOG(TE)) !From Murphy and Koop 2005
      SI=exp(-aux/8.314/TE)-1.0 !Assumes a saturated environment 

      fqr=0.0
      fpr=0.0
      fqi=0.0
      fpi=0.0
      fqg=0.0
      fpg=0.0
      RATE_Q = 0.0
      FPICE_=0.0
      xs =  0.8
      !Some constants
      rho_r = 1000.0 !Kg/m3 rain
      rho_g = 500.0 !Kg/m3 graupel
      rho_i = 900.0 !Kg/m3 ice/snow
      rho_a = PRE*100.0/MAPL_RGAS/TE 


      psc=(PRE*273.0/1000.0/TE)**0.54

      Tc=max(TE-273.15, -80.0)

      !fqi  = min(0.25*(1.0-exp((TE-273.15)/10.0)), 1.0) !From delGenio 2005
      !fqg  = 1.0-fqi

      fqg = RIM
      fqi = 1.0-fqg
      NG=NI*fqg
      NI_=NI-NG

      qr=QC*(1.0-FICE_)
      qg=QC*fqg*FICE_
      qi=QC*fqi*FICE_


      Dc_r =  min(max((VCN/842.0/psc)**1.25, 1e-6), 1.e-4) !Morrison2008, limit to 100 microns
      Dc_g =  min(max(1.0e-3*(VCN/19.3/psc)**2.70, 1.e-6), 1.e-2)  !Locatelli and Hobbs 1974

      Ah= 2.0*psc  
      Bh= 0.244-0.0049*min(Tc, 0.0)
      Bh=1/Bh!
      Dc_i =  0.01*((VCN/Ah)**Bh)  
      Dc_i =  min(max(Dc_i, 4.0e-6), 5.e-4)


      !rain===============Follow Liu 2006====================

      if ((qr .gt. 1.e-9) .and. (NL .gt. 1.0)) then

         lamr=min((MAPL_PI*rho_r*NL/rho_a/qr)**(1./3.), 1.0e8)
         aux=max(min(Dc_r*lamr, 10.0), 0.1)

         L =  1.0e-3*qr*rho_a !g cm-3
         N =  NL*1.e-6*rho_a !cm-3

         ! relative dispersion from Yangang Liu et al 2008 Environ. Res. Lett. 3 045021 doi:10.1088/1748-9326/3/4/045021
         !Corrected

         xs=0.070*((L/N)**(-0.140))   
         xs=max(min(xs, 1.7), 1.0001)      !maintains 0.0<e<1.0
         xs=xs*xs*xs   
         xs = (xs + sqrt(xs+8.0)*sqrt(xs) - 4.0)/8.0 ! from Eq. 2 (gives e^2)

         ! Autoconvesion rate from Liu, Yangang, et al. Journal of the atmospheric sciences 63.3 (2006): 1103-1109.      

         beta6=(1.0+3.0*xs)*(1.0+4.0*xs)*(1.0+5.0*xs) / (1.0+ xs)/(1.0+2.0*xs)

         !xs =aux*aux*aux
         !xs = max(min(1/xs, 20.0), 1.e-6)

         xs =  1.03e16*(L*L)/(N*SQRT(N))    !ratio of mean mass to critical mass


         xs=    max(min(20.0, xs), 1.0e-6)



         !Using miu =2.0  
         prc =  1.1e10*beta6*L*L*L*(1.0 -exp(-min(xs*xs, 20.0)))/N                         
         prc =  prc*1.0e3/rho_a !return to mixing ratio


         fpr=prc/qr  !COAUTO is a tuning factor that accounts for the increase in rain by accretion 


         ! NUmber autoconversion rate is approximated according to Liu, Yangang, et al. Geophysical Research Letters 34.16 (2007).
         ! Assumes e =  0.33. Replaces previous formulatin by Barahona et al. 2014 GMD

         !xs= 9.7e-17*N*sqrt(N)/L/L !Eq. 8c        
         !fnr=fpr*6.0/(aux*aux*aux+3.0*aux*aux+6.0*aux + 6.)


         fnr=fpr/(1.0+ (1./xs)) ! from Eqs. 13a and 13b 

         !print *, '====f===',xs, fpr, fnr 

      else
         fpr=0.0
         fnr=0.0 
      end if

      !graupel
      if (.false.) then 
         if ((qg .gt. 1.e-9) .and. (NG .gt. 1.0)) then

            lamg=min((MAPL_PI*rho_g*NG/rho_a/qg)**(1./3.), 1.0e6)
            aux=Dc_g*lamg
            fpg=exp(-aux)*(aux*aux*aux+3.0*aux*aux+6.0*aux + 6.)/6.0  !precipitated fraction
            fpg=min(fpg, 1.0)/DTL
            fng=exp(-aux)/DTL  
         else
            fpg=0.0
            fng=0.0
         end if


      else !use diffusion equation instead

         if ((qg .gt. 1.e-9) .and. (NG .gt. 1.0)) then
            lamg=min((MAPL_PI*rho_g*NG/rho_a/qg)**(1./3.), 1.0e6)
            aux=Dc_g*lamg
            !fpi=exp(-aux)*(aux*aux*aux+3.0*aux*aux+6.0*aux + 6.)/6.0  !precipitated fraction
            !fpi=min(fpi, 1.0)
            !fni=exp(-aux)
            fng=Ai*SI/Dc_g/Dc_g !Assume that ice crystals only grow by diffusion
            fng=fng*(1.0-exp(-aux)) !Only sizes with D<Dc can grow to Dc
            fpg=fng*(aux*aux*aux+3.0*aux*aux+6.0*aux + 6.)/6.0
         else
            fpg=0.0
            fng=0.0
         end if

      end if


      !ice
      if ((qi .gt. 1.e-9) .and. (NI_ .gt. 1.0)) then
         lami=(MAPL_PI*rho_i*NI_/rho_a/qi)**(1./3.)
         aux=Dc_i*lami
         !fpi=exp(-aux)*(aux*aux*aux+3.0*aux*aux+6.0*aux + 6.)/6.0  !precipitated fraction
         !fpi=min(fpi, 1.0)/DTL
         !fni=exp(-aux)/DTL
         fni=Ai*SI/Dc_i/Dc_i !Assume that ice crystals only grow by diffusion 
         fni=fni*(1.0-exp(-aux)) !Only sizes with D<Dc can grow to Dc 
         fpi=fni*(aux*aux*aux+3.0*aux*aux+6.0*aux + 6.)/6.0


      else
         fpi=0.0
         fni=0.0
      end if

      !COAUTO tuning constant since we do not consider accretion by precip
 fng=min(fng, 1.0/DTL)
 fpg=min(fpg, 1.0/DTL)
 fni=min(fni, 1.0/DTL)
 fpi=min(fpi, 1.0/DTL)
 fnr=min(fnr, 1.0/DTL)
 fpr=min(fpr, 1.0/DTL)


      QPT  =  fpr*qr*(1.0-FICE_)+fpg*qg*FICE_+fpi*qi*FICE_

      F_NL =  fnr
      F_NI =  (fqg*fng+fqi*fni)
      FPICE_ = fqi*fpi+fqg*fpg



      if (QC .gt. 1e-10) then
         RATE_Q=QPT/QC
      end if


      RATE_Q =  min(RATE_Q, 1.0/DTL)
      FPICE_=min(FPICE_, 1.0/DTL)

   end subroutine Qremoval


   !=================================================================

   function ICE_FRACTION (TEMP) RESULT(ICEFRCT)
      real,  intent(in) :: TEMP
      real           :: ICEFRCT, T_ICE_ALL, T_ICE_MAX

      T_ICE_ALL=238.0
      T_ICE_MAX=273.0

      ICEFRCT  = 0.00
      if ( TEMP <= T_ICE_ALL ) then 
         ICEFRCT = 1.000
      end if

      if ( (TEMP > T_ICE_ALL) .AND. (TEMP <= T_ICE_MAX) ) then 
         ICEFRCT = 1.00 -  ( TEMP - T_ICE_ALL ) / ( T_ICE_MAX - T_ICE_ALL ) 
      end if

      ICEFRCT = MIN(ICEFRCT,1.00)
      ICEFRCT = MAX(ICEFRCT,0.00)


   end function ICE_FRACTION



   !*************************************************************     
   !Approximation to the error function
   !*************************************************************
   real*8 function ERFAPP(x)

      real,  intent(in) :: x 
      real :: a  
      a=x*x*(1.27324d0+(0.147d0*x*x))/(1d0+(0.147d0*x*x))
      ERFAPP=SQRT(1d0-exp(-a))

      if (x .lt. 0.0) then  
         ERFAPP=-ERFAPP
      end if

   end function ERFAPP

! Manyin - Adapted from: fv_fill.F90
!    Modified it to work bottom-up
!>@brief The subroutine 'fill_z' is for mass-conservative filling of nonphysical negative values in the tracers. 
!>@details This routine takes mass from adjacent cells in the same column to fill negatives, if possible.
 subroutine fill_z(km, nq, q, dp, iic, iik)
   integer,  intent(in   ):: km                !< No. of levels
   integer,  intent(in   ):: nq                !< Total number of tracers
   real ,    intent(inout)::  q(km,nq)         !< tracer mixing ratio
   real ,    intent(in   ):: dp(km)            !< pressure thickness
   integer,  intent(in   ):: iic, iik          !< Top level, bottom level [in the range of 1 to 72 (or 132)]
! LOCAL VARIABLES:
   logical :: zfix
   real    :: dm(km)
   integer :: i, k, ic
   real    :: dq, sum0, sum1, fac

   do ic=1,nq

      zfix = .false.

! Top layer
      if( q(1,ic) < 0. ) then
          q(2,ic) = q(2,ic) + q(1,ic)*dp(1)/dp(2)
          q(1,ic) = 0.
      endif

! Bottom layer
      k = km
!     if( q(k,ic)<0. .and. q(k-1,ic)>0.) then
!         zfix = .true.
! Borrow from above
!         dq = min ( q(k-1,ic)*dp(k-1), -q(k,ic)*dp(k) )
!         q(k-1,ic) = q(k-1,ic) - dq/dp(k-1)
!         q(k  ,ic) = q(k  ,ic) + dq/dp(k  )
!     endif
      if( q(k  ,ic) < 0. ) then
          q(k-1,ic) = q(k-1,ic) + q(k,ic)*dp(k)/dp(k-1)
          q(k  ,ic) = 0.
      endif

! Interior

#if 0
    IF ( SUM(q(2:2+(km/4),ic)*dp(2:2+(km/4))) > SUM(q(km-(1+km/4):km-1,ic)*dp(km-(1+km/4):km-1)) ) THEN

! Top-down
      do k=2,km-1

         if ( q(k,ic)<0.0) zfix = .true.

         if ( q(k,ic)<0.0 .and. q(k-1,ic)>0. ) then
! Borrow from above
            dq =  min ( q(k-1,ic)*dp(k-1), -q(k,ic)*dp(k) )
            q(k-1,ic) = q(k-1,ic) - dq/dp(k-1)
            q(k  ,ic) = q(k  ,ic) + dq/dp(k  )
         endif

         if ( q(k,ic)<0.0 .and. q(k+1,ic)>0. ) then
! Borrow from below:
            dq =  min ( q(k+1,ic)*dp(k+1), -q(k,ic)*dp(k) )
            q(k+1,ic) = q(k+1,ic) - dq/dp(k+1)
            q(k  ,ic) = q(k  ,ic) + dq/dp(k  )
         endif

      enddo

    ELSE
#endif

! Bottom-up
      do k=km-1,2,-1

         if ( q(k,ic)<0.0) zfix = .true.

         if ( q(k,ic)<0.0 .and. q(k+1,ic)>0. ) then
! Borrow from below:
            dq =  min ( q(k+1,ic)*dp(k+1), -q(k,ic)*dp(k) )
            q(k+1,ic) = q(k+1,ic) - dq/dp(k+1)
            q(k  ,ic) = q(k  ,ic) + dq/dp(k  )
         endif

         if ( q(k,ic)<0.0 .and. q(k-1,ic)>0. ) then
! Borrow from above
            dq =  min ( q(k-1,ic)*dp(k-1), -q(k,ic)*dp(k) )
            q(k-1,ic) = q(k-1,ic) - dq/dp(k-1)
            q(k  ,ic) = q(k  ,ic) + dq/dp(k  )
         endif

      enddo

#if 0
    END IF
#endif


! Perform final check and non-local fix if needed
         if ( zfix ) then

           sum0 = 0.
           do k=km,1,-1
              dm(k) = q(k,ic)*dp(k)
              sum0 = sum0 + dm(k)
           enddo

           if ( sum0 > 0. ) then
! PRINT*,'RAS NEGATIVES SPREAD IN Z', iic, iik
             sum1 = 0.
             do k=km,1,-1
                sum1 = sum1 + max(0., dm(k))
             enddo
             fac = sum0 / sum1
             do k=km,1,-1
                q(k,ic) = max(0., fac*dm(k)/dp(k))
             enddo
           endif

         endif

   enddo
 end subroutine fill_z

subroutine init_Aer(aerout)

    type (AerProps), intent(inout) :: aerout
   
           aerout%num = 0.0
	   aerout%dpg =  1.0e-9
	   aerout%sig =  2.0
	   aerout%kap =  0.2
	   aerout%den = 2200.0
	   aerout%fdust  =  0.0
           aerout%fsoot  =  0.0
	   aerout%forg   =  0.0
	   aerout%nmods = 1
	   
   end subroutine init_Aer

END MODULE RAS
