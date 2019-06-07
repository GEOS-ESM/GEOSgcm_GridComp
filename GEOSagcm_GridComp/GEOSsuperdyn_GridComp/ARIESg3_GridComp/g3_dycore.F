      subroutine dycore ( IM,JM,LM,SIG,PTOP,KM,DT,
     .                    OMEGA, CP, RGAS, AE,
     .                    ROTATION, TILT, DLAM, DPHI,
     .                    PHS,  PKH,
     .                    PIB,  UOB,  VOB,  POB,  QOB,
     .                    PIM,  UOM,  VOM,  POM,  QOM,
     .                    PII,  UOI,  VOI,  POI,  QOI,
     .                    OMG,  VOR,  PSIGDOT, ASSLN, lattice )

C **********************************************************************
C  PURPOSE
C     Update Time-Tendencies of Prognostic Fields 
C     due to Hydrodynamical Processes
C
C  INPUT ARGUMENT DESCRIPTION
C
C     IM .... Number of Grid Intervals in Zonal      Direction
C     JM .... Number of Grid Intervals in Meridional Direction
C     LM .... Number of Vertical Levels
C     SIG ... (LM+1): Sigma at Interfaces. SIG(1)=0; SIG(LM+1)=1
C     PTOP .. Model Top Pressure
C     KM .... Number of Scalars, Including H2O, but not Theta.
C     DT .... Time-Step fron n-1 to n+1 (Seconds)
C
C     OMEGA.. Rotation rate (rad/sec)
C     CP..... Specific heat at constant pressure (J/(kg K))
C     RGAS... Gas contant (J/(kg K))
C     AE..... 'Earth' radius (meters)
C
C     ROTATION ... Longitude of the Computational Pole (Degrees) 
C     TILT ....... Latitude  of the Computational Pole (Degrees) 
C     DLAM ....... Longitude Grid Increments
C     DPHI ....... Latitude  Grid Increments
C
C     PHS ... (IM,JM): Surface Geopotential (m * m/sec**2)
C     PKH ... (IM,JM,LM+1): (P/P00)**KAPPA
C
C     PIB ... (IM,JM):    Mass (Psurf-Ptop) mb  at Current  Time-Level
C     UOB ... (IM,JM,LM): Zonal      Wind   m/s at Current  Time-Level
C     VOB ... (IM,JM,LM): Meridional Wind   m/s at Current  Time-Level
C     POB ... (IM,JM,LM): Potential Temp.   K   at Current  Time-Level
C     QOB ... (IM,JM,LM,KM): Scalar Fields      at Current  Time-Level
C
C     PIM ... (IM,JM):    Mass (Psurf-Ptop)     at Previous Time-Level
C     UOM ... (IM,JM,LM): Zonal      Wind       at Previous Time-Level
C     VOM ... (IM,JM,LM): Meridional Wind       at Previous Time-Level
C     POM ... (IM,JM,LM): Potential Temperature at Previous Time-Level
C     QOM ... (IM,JM,LM,KM): Scalar Fields      at Previous Time-Level
C
C  OUTPUT ARGUMENT DESCRIPTION-- Tendencies are in per second.
c
C     PII ... (IM,JM):       Updated Surface Pressure   Time-Tendency
C     UOI ... (IM,JM,LM):    Updated Zonal Wind         Time-Tendency
C     VOI ... (IM,JM,LM):    Updated Meridional Wind    Time-Tendency
C     POI ... (IM,JM,LM):    Updated PI-Weighted Theta  Time-Tendency
C     QOI ... (IM,JM,LM,KM): Updated PI-Weighted Scalar Time-Tendency
C
C     OMG ... (IM,JM,LM): Omega       Diagnostic (mb/sec)
C     VOR ... (IM,JM,LM): Vorticity   Diagnostic (1/sec)
C     PSIGDOT (IM,JM,LM): Pi*Sigdot   Diagnostic (mb/sec)
C     ASSLN.. Asselin Filter Coefficient
C
C  NOTES:
C     (1) JDIM is no longer used
C         Over-indexing vectorization technique no longer used.
C     (2) The Vertical Layers are numbered from TOP(1) to BOTTOM(LM).
C     (3) All Time-Tendencies are INCREMENTED (bumped).
C         The Momentum Time-Tendencies ARE NOT mass-weighted.
C         The Potential Temperature and Scalar Time-Tendencies ARE 
C         mass-weighted (by PI).
C     (4) JM is 180 degrees divided by the meridional grid size.
C     (5) UXX(I,J) are located half a grid interval EAST of PXX(I,J).
C         VXX(I,J) are located half a grid interval SOUTH of PXX(I,J).
C     (6) If PTOP>0, the PKH MUST be defined. 
C     (7) The previous time level fields (PIM,UOM,etc) are used for the
C         economical explicit calculation done in conjunction with
C         leap-frog steps. If you are not doing leap-frog or do not
C         wish to have economical explicit tendencies, pass the current
C         time-level fields twice (i.e., in PIB,UOB,etc and again in
C         PIM,UOM,etc.). 
C
C  SPACE REQUIREMENTS:
C     (1) Takes IM*JM*19+4*JM+2*IM words from the heap for STATIC storage;
C         these are kept throughout the run.
C     (2) Takes IM*JM*(LM+25) + 3*LM + 1 words from the heap for DYNAMIC
C         storage when PTOP=0; for PTOP!=0, add IM*JM*LM words.
C         All of this storage is freed before returning.
C
C **********************************************************************

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      REAL(kind=8) ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX
      REAL(kind=8) HALF, FOURTH, THIRD
      REAL(kind=8) pi

      PARAMETER ( ZERO=0.0)
      PARAMETER (  ONE=1.0)
      PARAMETER (  TWO=2.0)
      PARAMETER (THREE=3.0)
      PARAMETER ( FOUR=4.0)
      PARAMETER ( FIVE=5.0)
      PARAMETER ( SIX =6.0)
      PARAMETER (  PI =3.1415926535898)

      PARAMETER (FOURTH=ONE/FOUR )
      PARAMETER (  HALF=ONE/TWO  )
      PARAMETER ( THIRD=ONE/THREE)

      INTEGER    TMPVRS
      PARAMETER (TMPVRS=12)

      REAL(kind=8) ALPHA, ASSLN, thbar

c ARGUMENTS
c ---------
      INTEGER IM
      INTEGER JM
      INTEGER LM
      INTEGER KM

      REAL(kind=8) OMEGA
      REAL(kind=8) AE
      REAL(kind=8) CP
      REAL(kind=8) RGAS

      REAL(kind=8) PTOP
      REAL(kind=8) DT

      REAL(kind=8) SIG(LM+1)

      REAL(kind=8) UOB(IM,JM,LM)
      REAL(kind=8) VOB(IM,JM,LM)
      REAL(kind=8) POB(IM,JM,LM)
      REAL(kind=8) QOB(IM,JM,LM,1)
      REAL(kind=8) PIB(IM,JM)

      REAL(kind=8) UOM(IM,JM,LM)
      REAL(kind=8) VOM(IM,JM,LM)
      REAL(kind=8) POM(IM,JM,LM)
      REAL(kind=8) QOM(IM,JM,LM,1)
      REAL(kind=8) PIM(IM,JM)

      REAL(kind=8) UOI(IM,JM,LM)
      REAL(kind=8) VOI(IM,JM,LM)
      REAL(kind=8) POI(IM,JM,LM)
      REAL(kind=8) QOI(IM,JM,LM,1)
      REAL(kind=8) PII(IM,JM)

      REAL(kind=8) TILT, PHI_NP
      REAL(kind=8) ROTATION, LAMBDA_NP
      REAL(kind=8) DLAM(lattice%imglobal), DPHI(lattice%jmglobal)

      REAL(kind=8) PHS(IM,JM)
      REAL(kind=8) PKH(IM,JM,LM+1)

      REAL(kind=8)     OMG(IM,JM,LM)
      REAL(kind=8)     VOR(IM,JM,LM)
      REAL(kind=8) psigdot(IM,JM,LM)

c THESE ARE DYNAMIC (RESET EACH TIME) LOCALS
c ------------------------------------------
      REAL(kind=8) DSG(LM)

c LOCAL VECTOR TEMPORARIES
c ------------------------
      REAL(kind=8) ALF (0:im+1,0:jm+1)
      REAL(kind=8) BET (0:im+1,0:jm+1)
      REAL(kind=8) GAM (0:im+1,0:jm+1)
      REAL(kind=8) DEL (0:im+1,0:jm+1)
      REAL(kind=8) EPS (0:im+1,0:jm+1)
      REAL(kind=8) LAM (0:im+1,0:jm+1)
      REAL(kind=8) MUU (0:im+1,0:jm+1)
      REAL(kind=8) PSD (0:im+1,0:jm+1,1:lm)

      REAL(kind=8) VT1 (0:im+1,0:jm+1)
      REAL(kind=8) VT2 (0:im+1,0:jm+1)
      REAL(kind=8) VT3 (0:im+1,0:jm+1)
      REAL(kind=8) VT4 (0:im+1,0:jm+1)

      REAL(kind=8) USB (0:im+1, 0:jm+1)
      REAL(kind=8) VSB (0:im+1, 0:jm+1)
      REAL(kind=8) ZOB (0:im+1,-1:jm+2)

c GHOSTED LOCALS
c --------------
      real(kind=8) pgx(0:im+1,jm)
      real(kind=8)  gx(0:im+1,jm)
      real(kind=8) pgy(im,0:jm+1)

c VECTOR TEMPORARIES
c ------------------
      REAL(kind=8) PKL  (0:im+1,0:jm+1,LM)
      REAL(kind=8) PHI  (0:im+1,0:jm+1,LM)
      REAL(kind=8) theta(0:im+1,0:jm+1,LM)
      REAL(kind=8) pkht (0:im+1,0:jm+1,LM+1)

      REAL(kind=8) PBI (IM,JM)
      REAL(kind=8) PBJ (IM,JM)
      REAL(kind=8) BIP (IM,JM)
      REAL(kind=8) BJP (IM,JM)
      REAL(kind=8) ACH (IM,JM)
      REAL(kind=8) DDX (IM,JM)
      REAL(kind=8) DDY (IM,JM)
      REAL(kind=8) DPX (IM,JM)
      REAL(kind=8) DPY (IM,JM)
      REAL(kind=8) PIV (IM,JM)

c SCALAR TEMPORARIES
c ------------------
      REAL(kind=8) PSUMS
      REAL(kind=8) PSUMN

      REAL(kind=8) SUMSO
      REAL(kind=8) SUMNO

      REAL(kind=8) ST1
      REAL(kind=8) ST2

      LOGICAL LEAP

c STATIC LOCALS
c -------------
      real(kind=8),    allocatable, save :: DXUIJ(:,:)
      real(kind=8),    allocatable, save :: DYUIJ(:,:)
      real(kind=8),    allocatable, save :: DXVIJ(:,:)
      real(kind=8),    allocatable, save :: DYVIJ(:,:)
      real(kind=8),    allocatable, save :: D2PIJ(:,:)
      real(kind=8),    allocatable, save :: D2UIJ(:,:)
      real(kind=8),    allocatable, save :: D2VIJ(:,:)
      real(kind=8),    allocatable, save :: D2ZIJ(:,:)
      real(kind=8),    allocatable, save :: DYVIN(:,:)
      real(kind=8),    allocatable, save :: DXUIN(:,:)
      real(kind=8),    allocatable, save :: D2PIN(:,:)
      real(kind=8),    allocatable, save :: D2ZIN(:,:)
      real(kind=8),    allocatable, save :: FFFIJ(:,:)

c Uniform-Grid Standard Filter Arrays
c -----------------------------------
      real(kind=8), allocatable, save :: su (:,:)
      real(kind=8), allocatable, save :: sv (:,:)

c Stretched-Grid Convolution Filter Arrays
c ----------------------------------------
c     real, allocatable, save :: wp(:,:,:)   ! Convolution Weights for P-Point
c     real, allocatable, save :: wu(:,:,:)   ! Convolution Weights for U-Point
c     real, allocatable, save :: wv(:,:,:)   ! Convolution Weights for V-Point

c     logical stretched
c     save    stretched


C  THESE DETERMINE THE REINITIALIZATION OF THE HORIZONTAL GRID
c ------------------------------------------------------------
      INTEGER IM0, JM0
      DATA IM0/0/, JM0/0/
      REAL(kind=8) AE0
      DATA AE0/0./

      INTEGER I,j, L, K, LL
      integer p1,p2,p3,jmm1,jmm2,jmm3


c Initialization
c --------------
      theta(1:im,1:jm,1:lm)   = pob(1:im,1:jm,1:lm)
      pkht (1:im,1:jm,1:lm+1) = pkh(1:im,1:jm,1:lm+1)

      do L=1,lm
      call ghostx ( theta(1:im  ,1:jm,L)   ,theta(0:im+1,1:jm  ,L)   ,im  ,jm,1,    1,lattice,'east'  )
      call ghosty ( theta(0:im+1,1:jm,L)   ,theta(0:im+1,0:jm+1,L)   ,im+2,jm,1,0,1,1,lattice,'north' )
      call ghostx (  pkht(1:im  ,1:jm,L)   , pkht(0:im+1,1:jm  ,L)   ,im  ,jm,1,    1,lattice,'east'  )
      call ghosty (  pkht(0:im+1,1:jm,L)   , pkht(0:im+1,0:jm+1,L)   ,im+2,jm,1,0,1,1,lattice,'north' )
      enddo
      call ghostx (  pkht(1:im  ,1:jm,lm+1), pkht(0:im+1,1:jm  ,lm+1),im  ,jm,1,    1,lattice,'east'  )
      call ghosty (  pkht(0:im+1,1:jm,lm+1), pkht(0:im+1,0:jm+1,lm+1),im+2,jm,1,0,1,1,lattice,'north' )

            omg(:,:,:) = 0.0
        psigdot(:,:,:) = 0.0

c Set J-Index Range
c -----------------
      if ( lattice%pej.eq.0 ) then
           p1 = 1
           p2 = 2
           p3 = 3
      else
           p1 = 1
           p2 = 1
           p3 = 1
      endif
      if ( lattice%pej.eq.lattice%ny-1 ) then
           jmm1 = jm-1
           jmm2 = jm-2
           jmm3 = jm-3
      else
           jmm1 = jm
           jmm2 = jm
           jmm3 = jm
      endif

c Compute Optimum Brown-Campana Coefficent
c ----------------------------------------
      ALPHA = 0.965*( (ASSLN**2+4.0)*(ASSLN+2.0)/32.0 )
      
      IF( IM.NE.IM0 .OR. JM.NE.JM0 .OR. AE.NE.AE0 ) THEN

      IF(IM0.NE.0) THEN
        deallocate (DXUIJ)
        deallocate (DYUIJ)
        deallocate (DXVIJ)
        deallocate (DYVIJ)
        deallocate (D2PIJ)
        deallocate (D2UIJ)
        deallocate (D2VIJ)
        deallocate (D2ZIJ)
        deallocate (DYVIN)
        deallocate (DXUIN)
        deallocate (D2PIN)
        deallocate (D2ZIN)
        deallocate (FFFIJ)
        deallocate (SU   )
        deallocate (SV   )
c       deallocate (wp)
c       deallocate (wu)
c       deallocate (wv)
      ENDIF

        allocate ( DXUIJ(IM,JM) )
        allocate ( DYUIJ(IM,JM) )
        allocate ( DXVIJ(IM,JM) )
        allocate ( DYVIJ(IM,JM) )
        allocate ( D2PIJ(IM,JM) )
        allocate ( D2UIJ(IM,JM) )
        allocate ( D2VIJ(IM,JM) )
        allocate ( D2ZIJ(IM,JM) )
        allocate ( DYVIN(IM,JM) )
        allocate ( DXUIN(IM,JM) )
        allocate ( D2PIN(IM,JM) )
        allocate ( D2ZIN(IM,JM) )
        allocate ( FFFIJ(IM,JM) )
        allocate ( SU(lattice%imglobal+2 ,jm) )
        allocate ( SV(lattice%imglobal+2 ,jm) )
c       allocate ( wp(lattice%imglobal,lattice%imglobal,jm) )
c       allocate ( wu(lattice%imglobal,lattice%imglobal,jm) )
c       allocate ( wv(lattice%imglobal,lattice%imglobal,jm) )

C *********************************************************************
C *****                 Compute Grid Metric Factors                ****
C *********************************************************************

        phi_np    =     tilt*pi/180.
        lambda_np = rotation*pi/180.

        CALL GRIDH( IM,JM, AE, OMEGA
     *,             SIN(PHI_NP), LAMBDA_NP, DLAM, DPHI
     *,             DXUIJ, DYUIJ, DXVIJ, DYVIJ
     *,             D2PIJ, D2UIJ, D2VIJ, D2ZIJ
     *,             DYVIN, DXUIN, D2PIN, D2ZIN, FFFIJ, SU, SV
     *,                                  lattice )
c    *,             stretched, wp,wu,wv, lattice )

         IM0 = IM
         JM0 = JM
         AE0 = AE
      ENDIF

      CALL GRIDV ( LM,SIG,DSG )

C *********************************************************************
C *****               L-Independent Quantities                     ****
C *********************************************************************

      LEAP = DT .GT. ZERO

      call ghostx ( pib,pgx,im,jm,1,1,lattice,'east' )

c Mass-Points
c -----------
      DO j = p2,jmm1
      DO I = 1,IM
       PBI(I,j) = (PIB(I,j) + pgx(i+1,j)) * HALF
       BIP(I,j) = ONE / PBI(I,j)
       PIV(I,j) = ONE / PIB(I,j)
       VT2(I,j) = D2PIJ(I,j) * PIB(I,j)
      ENDDO
      ENDDO

      call ghostx ( vt2(1:im,1:jm),vt2(0:im+1,1:jm),im,jm,1,1,lattice,'east' )

      DO j = p2,jmm1
      DO I = 1,IM
       VT1(I,j) = (VT2(I,j) + VT2(i+1,j)) * HALF   ! Mass-Point
      ENDDO
      ENDDO

      call ghosty ( vt1(1:im,1:jm),vt1(1:im,0:jm+1),im,jm,1,0,1,1,lattice,'north' )
      call ghosty ( pib,pgy,                        im,jm,1,0,1,1,lattice,'north' )

c Vort-Points
c -----------
      DO j = p2,jmm2
      DO I = 1,IM
       PBJ(I,j) = (pgy(I,j+1) + PIB(I,j)) * HALF
       BJP(I,j) = ONE / PBJ(I,j )
       ACH(I,j) = D2ZIJ(I,j) / (HALF* ( VT1(I,j+1) + VT1(I,j)))
      ENDDO
      ENDDO

c Vort-Point at North Pole
c ------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
           call par_dot ( pib(1,jm-1),d2pij(1,jm-1),psumn,im,1,lattice )
           psumn = psumn * (HALF/FLOAT(lattice%imglobal))
           DO I =1,IM
             ACH(I,JM-1) = D2ZIJ(I,JM-1) * (ONE/PSUMN)
           ENDDO
      endif

c Vort-Point at South Pole
c ------------------------
      if ( lattice%pej.eq.0 ) then
           call par_dot ( pib(1,2),d2pij(1,2),psums,im,1,lattice )
           psums = psums * (HALF/FLOAT(lattice%imglobal))
           DO I =1,IM
             ACH(I,1) = D2ZIJ(I,1) * (ONE/PSUMS)
           ENDDO
      endif

C *********************************************************************
C *****                   Begin First L-Loop                       ****
C *********************************************************************

      DO 1000 L=1,LM

c P TO THE KAPPA
c --------------
       DO j = p2,jmm1
       DO I = 1,IM
        ST1 = PIB(I,j)*SIG(L  )+PTOP
        ST2 = PIB(I,j)*SIG(L+1)+PTOP
        PKL(I,j,L) = (ST2*PKH(I,j,L+1)-ST1*PKH(I,j,L)) / ( (ST2-ST1)*((RGAS/CP)+ONE) )
       ENDDO
       ENDDO

c MASS FLUXES
c -----------
      DO j = p2,jmm1
      DO I = 1,IM
       USB(I,j) = DYUIJ(I,j) * PBI(I,j) * UOB(I,j,L)   ! UWND-Point
      ENDDO
      ENDDO
      DO j = p2,jmm2
      DO I = 1,IM
       VSB(I,j) = DXVIJ(I,j) * PBJ(I,j) * VOB(I,j,L)   ! VWND-Point
      ENDDO
      ENDDO


c V-Mass Flux at North Pole
c -------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
           DO I = 1,IM
           VSB(I,JM-1) = ZERO
           ENDDO
      endif

c V-Mass Flux at South Pole
c -------------------------
      if ( lattice%pej.eq.0 ) then
           DO I = 1,IM
           VSB(I,1) = ZERO
           ENDDO
      endif

      call ghostx ( usb(1:im,  1:jm),usb(0:im+1,1:jm  ),im,  jm,1,     1,lattice,'both' )
      call ghosty ( usb(0:im+1,1:jm),usb(0:im+1,0:jm+1),im+2,jm,1,0,-1,1,lattice,'both' )

      call ghostx ( vsb(1:im,  1:jm),vsb(0:im+1,1:jm  ),im,  jm,1,     1,lattice,'both' )
      call ghosty ( vsb(0:im+1,1:jm),vsb(0:im+1,0:jm+1),im+2,jm,1,1,-1,1,lattice,'both' )

c HORIZONTAL ADVECTION OF TRACERS
c -------------------------------
      IF(KM.GT.0) THEN
        DO K = 1,KM
          CALL HADVECT( USB(0:im+1,0:jm+1),VSB(0:im+1,0:jm+1),QOB(1,1,L,K),IM,JM,
     .                  VT3(0:im+1,0:jm+1),VT1(0:im+1,0:jm+1),VT2(0:im+1,0:jm+1),lattice )
          DO j = p2,jmm1
          DO I = 1,IM
           QOI(I,j,L,K) =  QOI(I,j,L,K) - VT3(I,j)*D2PIN(I,j)
          ENDDO
          ENDDO
        ENDDO
      ENDIF

c HORIZONTAL ADVECTION OF POTENTIAL TEMPERATURE
c ---------------------------------------------
      CALL HADVECT( USB(0:im+1,0:jm+1),VSB(0:im+1,0:jm+1),POB(1,1,L),IM,JM,
     .              VT3(0:im+1,0:jm+1),VT1(0:im+1,0:jm+1),VT2(0:im+1,0:jm+1),lattice )
      DO j = p2,jmm1
      DO I = 1,IM
       POI(I,j,L) =  POI(I,j,L) - VT3(I,j)*D2PIN(I,j)
      ENDDO
      ENDDO

c COMPUTE CONVERGENCE
c -------------------
      DO j = p2,jmm1
      DO I = 1,IM
       PSD(I,j,L) = ( (USB(i-1,j)-USB(I,j))
     *              + (VSB(I,j-1)-VSB(I,j)) ) * (D2PIN(I,j)*DSG(L))
      ENDDO
      ENDDO

c COMPUTE VORTICITY
c -----------------
      DO j = p2,jmm1
      DO I = 1,IM
       VT2(I,j) = DXUIJ(I,j) * UOB(I,j,L)   ! UWND-Point
      ENDDO
      ENDDO

      call ghostx ( vob(1,1,L),    gx,              im,jm,1,     1,lattice,'east'  )
      call ghosty ( vt2(1:im,1:jm),vt2(1:im,0:jm+1),im,jm,1,0,-1,1,lattice,'north' )
      DO j = p2,jmm2
      DO I = 1,IM
       ZOB(I,j) = ( (gx(i+1,j) - VOB(I,j,L))   * DYVIJ(I,j)
     *          +   ( VT2(I,j) - VT2(I,j+1)) ) * D2ZIN(I,j)   ! Vort-Point
      ENDDO
      ENDDO

c Vort-Point at North Pole
c ------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
           call par_dot ( uob(1,jm-1,L),dxuij(1,jm-1),sumno,im,1,lattice )
           sumno = sumno * ( ONE/FLOAT(lattice%imglobal))
           DO I =1,IM
             ZOB(I,JM-1) = SUMNO * D2ZIN(1,JM-1)   ! Vort-Point at Pole
           ENDDO
      endif

c Vort-Point at South Pole
c ------------------------
      if ( lattice%pej.eq.0 ) then
           call par_dot ( uob(1,2,L),dxuij(1,2),sumso,im,1,lattice )
           sumso = sumso * (-ONE/FLOAT(lattice%imglobal))
           DO I =1,IM
             ZOB(I,1) = SUMSO * D2ZIN(1,1)   ! Vort-Point at Pole
           ENDDO
      endif

       DO j = p1,jmm1
       DO I = 1,IM
            VOR(I,j,L) = ZOB(I,j)
       ENDDO
       ENDDO

      DO j = p1,jmm1
      DO I = 1,IM
       ZOB(I,j) = (ZOB(I,j)+FFFIJ(I,j)) * ACH(I,j)   ! Vort-Point
      ENDDO
      ENDDO

      call ghostx ( zob(1:im,  1:jm),zob(0:im+1, 1:jm  ),im,  jm,1,    1,lattice,'both' )
      call ghosty ( zob(0:im+1,1:jm),zob(0:im+1,-1:jm+2),im+2,jm,1,1,1,2,lattice,'both' )

c COMPUTE VORTICITY COEFFICIENTS
c ------------------------------
      DO j = p2,jmm1
      DO I = 1,IM
       EPS(I,j) = ( ZOB(i  ,j-1) + ZOB(i  ,j  ) )
       LAM(I,j) = ( ZOB(i-1,j-1) + ZOB(i+1,j  ) )
       MUU(I,j) = ( ZOB(i-1,j  ) + ZOB(i+1,j-1) )
      ENDDO
      ENDDO
       
      DO j = p2,jmm2
      DO I = 1,IM
       ALF(I,j) = ( EPS(I,j) + ZOB(i+1,j  ) )
       BET(I,j) = ( EPS(I,j) + ZOB(i-1,j  ) )
       ALF(I,j) = (THREE*HALF)*ALF(I,j) - HALF*( MUU(I,j)+ZOB(I,j+1) )
       BET(I,j) = (THREE*HALF)*BET(I,j) - HALF*( LAM(I,j)+ZOB(I,j+1) )
      ENDDO
      ENDDO
      DO j = p3,jmm1
      DO I = 1,IM
       GAM(I,j) = ( EPS(I,j) + ZOB(i-1,j-1) )
       DEL(I,j) = ( EPS(I,j) + ZOB(i+1,j-1) )
       GAM(I,j) = (THREE*HALF)*GAM(I,j) - HALF*( MUU(I,j)+ZOB(I,j-2) )
       DEL(I,j) = (THREE*HALF)*DEL(I,j) - HALF*( LAM(I,j)+ZOB(I,j-2) )
      ENDDO
      ENDDO


      DO j = p2,jmm2
      DO I = 1,IM
       MUU(I,j) = HALF*(ZOB(I,j-1) - ZOB(I,j+1))  ! Vort-Point
      ENDDO
      ENDDO
      DO j = p3,jmm3
      DO I = 1,IM
       LAM(I,j) = HALF*(ZOB(i+1,j) - ZOB(i-1,j))  ! Vort-Point
      ENDDO
      ENDDO


c Lambda Averaging near North Pole
c --------------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
           DO I =1,IM
            LAM(I,JM-2) = (THREE/FIVE)*LAM(I,JM-3)
           ENDDO
      endif

c Lambda Averaging near South Pole
c --------------------------------
      if ( lattice%pej.eq.0 ) then
           DO I =1,IM
            LAM(I,2) = (THREE/FIVE)*LAM(I,3)
           ENDDO
      endif

c Ghost ALF,BET,GAM,DEL,LAM,MUU
c -----------------------------
      call ghostx ( alf(1:im,1:jm),alf(0:im+1,1:jm),im,jm,1,1,lattice,'west' )
      call ghostx ( del(1:im,1:jm),del(0:im+1,1:jm),im,jm,1,1,lattice,'west' )
      call ghostx ( muu(1:im,1:jm),muu(0:im+1,1:jm),im,jm,1,1,lattice,'west' )

      call ghosty ( gam(1:im  ,1:jm),gam(1:im  ,0:jm+1),im  ,jm,1,0,1,1,lattice,'north' )
      call ghosty ( del(0:im+1,1:jm),del(0:im+1,0:jm+1),im+2,jm,1,0,1,1,lattice,'north' )
      call ghosty ( lam(1:im  ,1:jm),lam(1:im  ,0:jm+1),im  ,jm,1,1,1,1,lattice,'south' )

c U INCREMENT
c -----------
      DO j = p2,jmm2
      DO I = 1,IM
       EPS(I,j) = BET(I,j) * VSB(i  ,j  )
     *          + ALF(I,j) * VSB(i+1,j  )
     *          - LAM(I,j) * USB(i  ,j+1)
      ENDDO
      ENDDO

c Zero EPS near North Pole
c ------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
           DO I = 1,IM
           EPS(I,JM-1) = ZERO
           ENDDO
      endif

      DO j = p3,jmm1
      DO I = 1,IM
       EPS(I,j) = EPS(I,j  )
     *          + GAM(I,j  ) * VSB(i  ,j-1)
     *          + DEL(I,j  ) * VSB(i+1,j-1)
     *          + LAM(I,j-1) * USB(i  ,j-1)
      ENDDO
      ENDDO

      DO j = p2,jmm1
      DO I = 1,IM
       UOI(I,j,L) = UOI(I,j,L) + (THIRD*FOURTH)*EPS(I,j)*DXUIN(I,j)
      ENDDO
      ENDDO

c V INCREMENT
c -----------
      DO j = p2,jmm2
      DO I = 1,IM
       ST1 = BET(I  ,j  ) * USB(i  ,j  )
     *     + GAM(I  ,j+1) * USB(i  ,j+1)
     *     + DEL(i-1,j+1) * USB(i-1,j+1)
     *     + ALF(i-1,j  ) * USB(i-1,j  )
     *     + MUU(I  ,j  ) * VSB(i+1,j  )
     *     - MUU(i-1,j  ) * VSB(i-1,j  )
       VOI(I,j,L) = VOI(I,j,L) - (THIRD*FOURTH)*ST1*DYVIN(I,j)
      ENDDO
      ENDDO

1000  CONTINUE

C *********************************************************************
C *****               Begin L-Critical Section                     ****
C *********************************************************************

c VERTICAL INTERGRAL OF CONTINUITY EQUATION
c -----------------------------------------
      DO L = 2,LM
      DO j = p2,jmm1
      DO I = 1,IM
      PSD(I,j,L) = PSD(I,j,L-1) + PSD(I,j,L)
      ENDDO
      ENDDO
      ENDDO

c SURFACE PRESSURE TENDENCY
c -------------------------
      DO j = p2,jmm1
      DO I = 1,IM
      PII(I,j) = PII(I,j) + PSD(I,j,LM)
      ENDDO
      ENDDO

c DDX AND DDY ARE USED IN OMEGA CALCULATION
c -----------------------------------------
      DO j = p2,jmm1
      DO I = 1,IM
      DDX(I,j) = (pgx(i+1,j)-PIB(I,j))*DYUIJ(I,j)*PBI(I,j)
      ENDDO
      ENDDO

      DO j = p2,jmm2
      DO I = 1,IM
      DDY(I,j) = (pgy(I,j+1)-PIB(I,j))*DXVIJ(I,j)*PBJ(I,j)
      ENDDO
      ENDDO

C *********************************************************************
C *****                  Begin Second L-Loop                       ****
C *********************************************************************

      DO 2000 LL = 1,2
      DO 2000  L = LL,LM-1,2

C VERTICAL ADVECTION
C ------------------
      DO j = p2,jmm1
      DO I = 1,IM
          PSD(I,j,L) = PSD(I,j,L) - PSD(I,j,LM)*SIG(L+1)
      psigdot(i,j,L) = psd(i,j,L)   ! PSIGDOT Diagnostic
      ENDDO
      ENDDO

      call ghostx ( psd(1:im,  1:jm,L),psd(0:im+1,1:jm  ,L),im,  jm,1,    1,lattice,'east' )
      call ghosty ( psd(0:im+1,1:jm,L),psd(0:im+1,0:jm+1,L),im+2,jm,1,0,1,1,lattice,'both' )

       DO j = p2,jmm1
       DO I = 1,IM
        ALF(I,j) = (PKH(I,j,L+1)-PKL(I,j,L)) / (PKL(I,j,L+1)-PKL(I,j,L))
       ENDDO
       ENDDO

      DO j = p2,jmm1
      DO I = 1,IM
       ST1 = ( (POB(I,j,L)-POB(I,j,L+1))*ALF(I,j) + POB(I,j,L+1) )*PSD(I,j,L)
       POI(I,j,L  ) = POI(I,j,L  ) - ST1 * (ONE/DSG(L  ))
       POI(I,j,L+1) = POI(I,j,L+1) + ST1 * (ONE/DSG(L+1))
      ENDDO
      ENDDO

      IF(KM.GT.0) THEN
       DO K =1,KM
        DO j = p2,jmm1
        DO I = 1,IM
         ST1 = ( (QOB(I,j,L,K)+QOB(I,j,L+1,K))*HALF )*PSD(I,j,L)
         QOI(I,j,L  ,K) = QOI(I,j,L  ,K) - ST1 * (ONE/DSG(L  ))
         QOI(I,j,L+1,K) = QOI(I,j,L+1,K) + ST1 * (ONE/DSG(L+1))
        ENDDO
        ENDDO
       ENDDO
      ENDIF

      DO j = p2,jmm1
      DO I = 1,IM
       ST1 = ( UOB(I,j,L+1)-UOB(I,j,L) )*( PSD(i+1,j,L)+PSD(I,j,L) )*BIP(I,j)
       UOI(I,j,L  ) = UOI(I,j,L  ) - ST1 * (FOURTH/DSG(L  ))
       UOI(I,j,L+1) = UOI(I,j,L+1) - ST1 * (FOURTH/DSG(L+1))
      ENDDO
      ENDDO

      DO j = p2,jmm2
      DO I = 1,IM
       ST1 = ( VOB(I,j,L+1)-VOB(I,j,L) )*( PSD(I,j+1,L)+PSD(I,j,L) )*BJP(I,j)
       VOI(I,j,L  ) = VOI(I,j,L  ) - ST1 * (FOURTH/DSG(L  ))
       VOI(I,j,L+1) = VOI(I,j,L+1) - ST1 * (FOURTH/DSG(L+1))
      ENDDO
      ENDDO

c Omega Diagnostic
c ----------------
       DO j = p2,jmm1
       DO I = 1,IM
        VT1(I,j) = UOB(I,j,L  )*DDX(I,j)   ! UWND-Point
        VT3(I,j) = UOB(I,j,L+1)*DDX(I,j)   ! UWND-Point
       ENDDO
       ENDDO

       DO j = p2,jmm2
       DO I = 1,IM
        VT2(I,j) = VOB(I,j,L  )*DDY(I,j)   ! VWND-Point
        VT4(I,j) = VOB(I,j,L+1)*DDY(I,j)   ! VWND-Point
       ENDDO
       ENDDO

c Zero VT2 and VT4 at North Pole
c ------------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
       DO I = 1,IM
        VT2(I,JM-1) = ZERO
        VT4(I,JM-1) = ZERO
       ENDDO
      endif

c Zero VT2 and VT4 at South Pole
c ------------------------------
      if ( lattice%pej.eq.0 ) then
       DO I = 1,IM
        VT2(I,1) = ZERO
        VT4(I,1) = ZERO
       ENDDO
      endif

      call ghostx ( vt1(1:im,1:jm),vt1(0:im+1,1:jm),im,jm,1,     1,lattice,'west'  )
      call ghostx ( vt3(1:im,1:jm),vt3(0:im+1,1:jm),im,jm,1,     1,lattice,'west'  )
      call ghosty ( vt2(1:im,1:jm),vt2(1:im,0:jm+1),im,jm,1,1,-1,1,lattice,'south' )
      call ghosty ( vt4(1:im,1:jm),vt4(1:im,0:jm+1),im,jm,1,1,-1,1,lattice,'south' )

       DO j = p2,jmm1
       DO I = 1,IM
        ALF(I,j) = (PKH(I,j,L+1)-PKL(I,j,L  )) / (PKH(I,j,L+1)-PKH(I,j,L  ))
        BET(I,j) = (PKL(I,j,L+1)-PKH(I,j,L+1)) / (PKH(I,j,L+2)-PKH(I,j,L+1))
       ENDDO
       ENDDO

       DO j = p2,jmm1
       DO I = 1,IM
                 ST1 = SIG(L+1)*PSD(I,j,LM) + PSD(I,j,L)
                 ST2 = HALF*(VT1(i-1,j)+VT1(I,j)+VT2(I,j)+VT2(I,j-1)) * SIG(L+1)*D2PIN(I,j)*PIV(I,j)
        OMG(I,j,L)   = OMG(I,j,L  ) + ALF(I,j) * (ST1 + ST2)
                 ST2 = HALF*(VT3(i-1,j)+VT3(I,j)+VT4(I,j)+VT4(I,j-1)) * SIG(L+1)*D2PIN(I,j)*PIV(I,j)
        OMG(I,j,L+1) = OMG(I,j,L+1) + BET(I,j) * (ST1 + ST2)
       ENDDO
       ENDDO

       IF(L.EQ.LM-1) THEN
          DO j = p2,jmm1
          DO I = 1,IM
           BET(I,j) = (PKH(I,j,L+2)-PKL(I,j,L+1)) / (PKH(I,j,L+2)-PKH(I,j,L+1))
          ENDDO
          ENDDO
        DO j = p2,jmm1
        DO I = 1,IM
                  ST1 = PSD(I,j,LM)
                  ST2 = HALF*(VT3(i-1,j)+VT3(I,j)+VT4(I,j)+VT4(I,j-1)) * D2PIN(I,j)*PIV(I,j)
         OMG(I,j,L+1) = OMG(I,j,L+1) + BET(I,j) * (ST1 + ST2)
        ENDDO
        ENDDO
       ENDIF

2000  CONTINUE

C *********************************************************************
C *****               Begin L-Critical Section                     ****
C *********************************************************************

C Filtering of Scalars on Stretched-Grid
c --------------------------------------
c     if( stretched ) then
c         call avrxg ( poi,wp,im,jm,lm,lattice )           ! Potential Temperature Tendency
c         call avrxg ( pii,wp,im,jm,1 ,lattice )           ! Surface   Pressure    Tendency
c         do k=1,km
c         call avrxg ( qoi(1,1,1,k),wp,im,jm,lm,lattice )  ! Tracer Tendencies
c         enddo
c         call avrxg ( omg    ,wp,im,jm,lm,lattice )       ! Vertical Velocity (Pressure)
c         call avrxg ( psigdot,wp,im,jm,lm,lattice )       ! Vertical Velocity (Sigma)
c     else
C Filtering of Scalars on Uniform-Grid
c ------------------------------------
          call avrx ( poi,im,jm,lm,su,lattice )
          call avrx ( pii,im,jm,1 ,su,lattice )
          do k = 1,km
          call avrx ( qoi(1,1,1,K),im,jm,lm,su,lattice )
          enddo
          call avrx ( omg    ,im,jm,lm,su,lattice )
          call avrx ( psigdot,im,jm,lm,su,lattice )
c     endif

c AT THIS POINT THETA, Q, AND PS INCREMENTS ARE COMPLETE
c AND CAN BE USED FOR COMPUTING THE PRESSURE GRADIENT
c IN THE ECONOMICAL EXPLICIT SCHEME
c ---------------------------------
      IF(LEAP) THEN
       DO j = p2,jmm1
       DO I = 1,IM
        ST1 = PIM(I,j) + DT*PII(I,j)
        VT3(I,j) = ONE / ST1
        VT1(I,j) = ALPHA*(PIM(I,j)+ST1) + (ONE-TWO*ALPHA)*PIB(I,j)
       ENDDO
       ENDDO
      ELSE
       DO j = p2,jmm1
       DO I = 1,IM
        VT1(I,j) = PIB(I,j)
       ENDDO
       ENDDO
      ENDIF

      call ghostx ( vt1(1:im,  1:jm),vt1(0:im+1,1:jm  ),im,  jm,1,    1,lattice,'east'  )
      call ghosty ( vt1(0:im+1,1:jm),vt1(0:im+1,0:jm+1),im+2,jm,1,0,1,1,lattice,'north' )

c GRADIENT OF SURFACE PRESSURE
c ----------------------------
      DO j = p2,jmm1
      DO I = 1,IM
       DPX(I,j) = (VT1(i+1,j) - VT1(I,j))
      ENDDO
      ENDDO

      DO j = p2,jmm2
      DO I = 1,IM
       DPY(I,j) = (VT1(I,j+1) - VT1(I,j))
      ENDDO
      ENDDO

c HYDROSTATIC EQUATION
c --------------------
      DO j = p2,jmm1
      DO I = 1,IM
       PHI(I,j,LM) = PHS(I,j)
      ENDDO
      ENDDO

      DO L = LM,2,-1
         IF(LEAP) THEN
          DO j = p2,jmm1
          DO I = 1,IM
           VT1(I,j) = ( POM(I,j,L)*PIM(I,j) + DT*POI(I,j,L) )*VT3(I,j)
           VT1(I,j) = ALPHA*(POM(I,j,L) + VT1(I,j)) + (ONE-TWO*ALPHA)*POB(I,j,L)
          ENDDO
          ENDDO
         ELSE
          DO j = p2,jmm1
          DO I = 1,IM
           VT1(I,j) = POB(I,j,L)
          ENDDO
          ENDDO
         ENDIF
         DO j = p2,jmm1
         DO I = 1,IM
                 ST1  = VT1(I,j)
         PHI(I,j,L-1) = PHI(I,j,L) + ST1*(CP*(PKH(I,j,L+1)-PKH(I,j,L)))
         PHI(I,j,L  ) = PHI(I,j,L) + ST1*(CP*(PKH(I,j,L+1)-PKL(I,j,L)))
         ENDDO
         ENDDO
      ENDDO

      IF(LEAP) THEN
        DO j = p2,jmm1
        DO I = 1,IM
         VT1(I,j) = ( POM(I,j,1)*PIM(I,j) + DT*POI(I,j,1) )*VT3(I,j)
         VT1(I,j) = ALPHA*(POM(I,j,1) + VT1(I,j)) + (ONE-TWO*ALPHA)*POB(I,j,1)
        ENDDO
        ENDDO
      ELSE
        DO j = p2,jmm1
        DO I = 1,IM
         VT1(I,j) = POB(I,j,1)
        ENDDO
        ENDDO
      ENDIF

      DO j = p2,jmm1
      DO I = 1,IM
               ST1  = VT1(I,j)
         PHI(I,j,1) = PHI(I,j,1) + ST1*(CP*(PKH(I,j,2)-PKL(I,j,1)))
      ENDDO
      ENDDO

C *********************************************************************
C *****                   Begin Third L-Loop                       ****
C *********************************************************************

      DO 3000 L = 1,LM

      call ghostx ( vob(1:im,1:jm,L),vsb(0:im+1,1:jm),im,jm,1,     1,lattice,'both' )
      call ghosty ( uob(1:im,1:jm,L),usb(1:im,0:jm+1),im,jm,1,0,-1,1,lattice,'both' )

c KINETIC ENERGY  (SADOURNEY SICKPROOF), THERE IS A 1/2 IN DXDY[UV]
c -----------------------------------------------------------------
      ST2 = FIVE/SIX

      DO j = p3,jmm2
      DO I = 1,IM
       ST1      =  HALF*( USB(i,j-1)+USB(i,j+1) )
       VT1(I,j) = ( (ST2)*UOB(i,j,L)*UOB(i,j,L) + (ONE-ST2)*ST1*ST1 ) *D2UIJ(I,j)
      ENDDO
      ENDDO

      DO j = p2,jmm2
      DO I = 1,IM
       ST1      =  HALF*( VSB(i+1,j)+VSB(i-1,j) )
       VT2(I,j) = ( (ST2)*VOB(i,j,L)*VOB(i,j,L) + (ONE-ST2)*ST1*ST1 ) *D2VIJ(I,j)
      ENDDO
      ENDDO

c Sickproof at North Pole
c -----------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       VT1(I,JM-1) = UOB(I,JM-1,L)*UOB(I,JM-1,L)*D2UIJ(I,JM-1)
       VT2(I,JM-1) = ZERO
      ENDDO
      endif

c Sickproof at South Pole
c -----------------------
      if ( lattice%pej.eq.0 ) then
      DO I = 1,IM
       VT1(I,2) = UOB(I,2,L)*UOB(I,2,L)*D2UIJ(I,2)
       VT2(I,1) = ZERO
      ENDDO
      endif

      call ghostx ( vt1(1:im,1:jm),vt1(0:im+1,1:jm),im,jm,1,    1,lattice,'west'  )
      call ghosty ( vt2(1:im,1:jm),vt2(1:im,0:jm+1),im,jm,1,1,1,1,lattice,'south' )

      DO  j = p2,jmm1
      DO  I = 1,IM
       PHI(I,j,L) = PHI(I,j,L)
     *            + HALF*( HALF*(VT1(i-1,j) + VT1(I,j  ))
     *            +        HALF*(VT2(I  ,j) + VT2(I,j-1)) )*D2PIN(I,j)
      ENDDO
      ENDDO

      call ghostx ( phi(1:im  ,1:jm,L),  phi(0:im+1,1:jm  ,L),im  ,jm,1,    1,lattice,'east'  )
      call ghosty ( phi(0:im+1,1:jm,L),  phi(0:im+1,0:jm+1,L),im+2,jm,1,0,1,1,lattice,'north' )
      call ghostx ( pkl(1:im  ,1:jm,L),  pkl(0:im+1,1:jm  ,L),im  ,jm,1,    1,lattice,'east'  )
      call ghosty ( pkl(0:im+1,1:jm,L),  pkl(0:im+1,0:jm+1,L),im+2,jm,1,0,1,1,lattice,'north' )

c PRESSURE GRADIENT FORCE, INCLUDING GRAD OF KE,  GAM IS (CP THETA DP/DPI)
c ------------------------------------------------------------------------
#if junk
c conservation of irrotationality method
c --------------------------------------
      DO j = p2,jmm1
      DO I = 1,IM
        GAM(i,j) = ( SIG(L  )*(PKL(i,j,L  )-PKH(i,j,L))
     *             + SIG(L+1)*(PKH(i,j,L+1)-PKL(i,j,L)) )*(CP/DSG(L))*PIV(i,j)*POB(i,j,L)
      ENDDO
      ENDDO
      call ghostx ( gam(1:im  ,1:jm  ),gam(0:im+1,1:jm  )  ,im  ,jm,1,    1,lattice,'east'  )
      call ghosty ( gam(0:im+1,1:jm  ),gam(0:im+1,0:jm+1)  ,im+2,jm,1,0,1,1,lattice,'north' )

      DO j = p2,jmm1
      DO I = 1,IM
       ST1        = ( phi(i+1,j,L)-phi(i,j,L) ) + DPX(i,j) * HALF*( gam(i+1,j)+gam(i,j) )
       UOI(I,j,L) = UOI(I,j,L) - ST1*DXUIN(i,j)
      ENDDO
      ENDDO

      DO j = p2,jmm2
      DO I = 1,IM
       ST1        = ( phi(i,j+1,L)-phi(i,j,L) ) + DPY(i,j) * HALF*( gam(i,j+1)+gam(i,j) )
       VOI(I,j,L) = VOI(I,j,L) - ST1*DYVIN(i,j)
      ENDDO
      ENDDO
#endif

c conservation of energy method
c -----------------------------
      do j=p2,jmm1
      do i=1,im

c finite-volume theta_bar
c -----------------------
c      thbar      = ( theta(i+1,j,L)*( pkht(i+1,j,L+1)-pkht(i+1,j,L) )
c    .              + theta(i  ,j,L)*( pkht(i  ,j,L+1)-pkht(i  ,j,L) ) )
c    .                              /( pkht(i+1,j,L+1)-pkht(i+1,j,L) 
c    .                               + pkht(i  ,j,L+1)-pkht(i  ,j,L) )
 
c conservation theta_bar
c ----------------------
       thbar      = ( theta(i+1,j,L)+theta(i,j,L) )*0.5

       ST1        = ( phi(i+1,j,L)-phi(i,j,L) ) + cp*thbar*( pkl(i+1,j,L)-pkl(i,j,L) )
       UOI(I,j,L) = UOI(I,j,L) - ST1*DXUIN(i,j)
      enddo
      enddo

      do j=p2,jmm2
      do i=1,im

c finite-volume theta_bar
c -----------------------
c      thbar      = ( theta(i,j+1,L)*( pkht(i,j+1,L+1)-pkht(i,j+1,L) )
c    .              + theta(i  ,j,L)*( pkht(i  ,j,L+1)-pkht(i  ,j,L) ) )
c    .                              /( pkht(i,j+1,L+1)-pkht(i,j+1,L) 
c    .                               + pkht(i  ,j,L+1)-pkht(i  ,j,L) )
 
c conservation theta_bar
c ----------------------
       thbar      = ( theta(i,j+1,L)+theta(i,j,L) )*0.5

       ST1        = ( phi(i,j+1,L)-phi(i,j,L) ) + cp*thbar*( pkl(i,j+1,L)-pkl(i,j,L) )
       VOI(I,j,L) = VOI(I,j,L) - ST1*DYVIN(i,j)
      enddo
      enddo

3000  CONTINUE

C************************************************************
C***************** END THIRD L LOOP *************************
C************************************************************

C Filtering of Momentum Fields on Stretched-Grid
c ----------------------------------------------
c     if( stretched ) then
c         call avrxg ( uoi,wu,im,jm,lm,lattice )  ! U-Wind Tendency
c         call avrxg ( voi,wv,im,jm,lm,lattice )  ! V-Wind Tendency
c     else
C Filtering of Momentum Fields on Uniform-Grid
c --------------------------------------------
          call avrx  ( uoi,im,jm,lm,su,lattice )  ! U-Wind Tendency
          call avrx  ( voi,im,jm,lm,sv,lattice )  ! V-Wind Tendency
c     endif

      return
      end

      subroutine gridh ( IM,JM, AE, OMEGA
     *,                  SINPNP, LAMNP,DLAM,DPHI
     *,                  DXUIJ, DYUIJ, DXVIJ, DYVIJ
     *,                  D2PIJ, D2UIJ, D2VIJ, D2ZIJ
     *,                  DYVIN, DXUIN, D2PIN, D2ZIN, FFFIJ, SU, SV
     *,                                       lattice )
c    *,                  stretched, wp,wu,wv, lattice )

      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      REAL(kind=8) ONE
      REAL(kind=8) TWO
      REAL(kind=8) HALF
      REAL(kind=8) FOURTH
      REAL(kind=8) ZERO
      REAL(kind=8) PI
      REAL(kind=8) UNDEF,GETCON

      PARAMETER (   ONE=1.0 )
      PARAMETER (   TWO=2.0 )
      PARAMETER (  ZERO=0.0 )
      PARAMETER (  HALF=0.5 )
      PARAMETER (FOURTH=0.25)
      PARAMETER (    PI=3.1415926535898)

      REAL(kind=8) SINPNP
      REAL(kind=8) LAMNP

      INTEGER IM, JM

      REAL(kind=8) AE
      REAL(kind=8) OMEGA


      REAL(kind=8) DXUIJ(IM,JM)
      REAL(kind=8) DYUIJ(IM,JM)
      REAL(kind=8) DYVIJ(IM,JM)
      REAL(kind=8) DXVIJ(IM,JM)

      REAL(kind=8) DYVIN(IM,JM)
      REAL(kind=8) DXUIN(IM,JM)
      REAL(kind=8) D2ZIN(IM,JM)
      REAL(kind=8) D2PIN(IM,JM)

      REAL(kind=8) D2PIJ(IM,JM)
      REAL(kind=8) D2UIJ(IM,JM)
      REAL(kind=8) D2VIJ(IM,JM)
      REAL(kind=8) D2ZIJ(IM,JM)

      REAL(kind=8) FFFIJ(IM,JM)

      REAL(kind=8) SU(lattice%imglobal+2,JM)
      REAL(kind=8) SV(lattice%imglobal+2,JM)

      REAL(kind=8) DLAM (lattice%imglobal)
      REAL(kind=8) DPHI (lattice%jmglobal)

      real(kind=8) gostx (0:im+1,1:jm  )
      real(kind=8) gosty (1:im  ,0:jm+1)
      real(kind=8) gostxy(0:im+1,0:jm+1)

c Stretched-Grid Variables
c ------------------------
c     real(kind=8) wp(lattice%imglobal,lattice%imglobal,jm)
c     real(kind=8) wu(lattice%imglobal,lattice%imglobal,jm)
c     real(kind=8) wv(lattice%imglobal,lattice%imglobal,jm)
c     logical stretched


C   LOCALS

      REAL(kind=8) COSPNP, FS, FC, SUM,DOT

      REAL(kind=8) LAMBDA(lattice%imglobal)
      REAL(kind=8) PHI   (lattice%jmglobal)

      REAL(kind=8) DXPIJ(IM,JM)
      REAL(kind=8) DYPIJ(IM,JM)
      REAL(kind=8)   SL (lattice%imglobal+1)
      REAL(kind=8)   CL (lattice%imglobal+1)
      REAL(kind=8)   CC (lattice%jmglobal+1)
      REAL(kind=8)   SC (lattice%jmglobal+1)
      REAL(kind=8)   CV (lattice%jmglobal+1)

      INTEGER I,J
      REAL(kind=8) SUMIN, SVMIN
      integer p1,p2,p3,jmm1,jmm2,jmm3

      UNDEF  = getcon('UNDEF')
      COSPNP = SQRT( ONE-SINPNP**2 )


c Set J-Index Range
c -----------------
      if ( lattice%pej.eq.0 ) then
           p1 = 1
           p2 = 2
           p3 = 3
      else
           p1 = 1
           p2 = 1
           p3 = 1
      endif
      if ( lattice%pej.eq.lattice%ny-1 ) then
           jmm1 = jm-1
           jmm2 = jm-2
           jmm3 = jm-3
      else
           jmm1 = jm
           jmm2 = jm
           jmm3 = jm
      endif

c Check for Stretched-Grid
c ------------------------
c     stretched = .false.
c     i = 2
c     do while( .not.stretched .and. i.le.lattice%imglobal )
c     stretched = dlam(i).ne.dlam(1)
c     i = i+1
c     enddo
c     j = 2
c     do while( .not.stretched .and. j.le.lattice%jmglobal )
c     stretched = dphi(j).ne.dphi(1)
c     j = j+1
c     enddo

C Global Longitude and Latitude
c -----------------------------
      LAMBDA(1) = -PI
      DO I = 2,lattice%imglobal
      LAMBDA(I) = LAMBDA(I-1) + DLAM(I-1)
      ENDDO

      PHI(1)  = -PI*HALF
      DO J = 2,lattice%jmglobal-1
      PHI(J)  = PHI(J-1) + DPHI(J-1)
      ENDDO
      PHI(lattice%jmglobal) = PI*HALF

C HORIZONTAL GRID
c ---------------
      DO I = 1,lattice%imglobal
       SL(I) = SIN( LAMBDA(i)-LAMNP )
       CL(I) = COS( LAMBDA(i)-LAMNP )
      ENDDO
       SL(lattice%imglobal+1) = SL(1)
       CL(lattice%imglobal+1) = CL(1)

      DO J = 1,lattice%jmglobal
        CC(J) = COS( PHI(j) )
        SC(J) = SIN( PHI(j) )
      ENDDO
       CC(1) = ZERO
       SC(1) = -ONE
       CC(lattice%jmglobal)   = ZERO
       SC(lattice%jmglobal)   = ONE
       CC(lattice%jmglobal+1) = undef
       SC(lattice%jmglobal+1) = undef

      DO J = 2,lattice%jmglobal-2
       CV(J)  = HALF*(CC(J+1)+CC(J))
      ENDDO
       CV(1) = ZERO
       CV(lattice%jmglobal-1) = ZERO
       CV(lattice%jmglobal  ) = undef
       CV(lattice%jmglobal+1) = undef


c DXU & DYV ARE DEFINED FROM THE GRID:  ALL OTHER FACTORS ARE DEFINED IN TERMS OF THEM
c ------------------------------------------------------------------------------------
      DO J = 1,JM
      DO I = 1,IM
       DXUIJ(I,J)   = AE*DLAM( lattice%iglobal(i) )*CC( lattice%jglobal(j) )
      ENDDO
      ENDDO

      DO J = 1,JM
      DO I = 1,IM
       DYVIJ(I,J)   = AE*DPHI( lattice%jglobal(j) )
      ENDDO
      ENDDO


c DXP AND DXV
c -----------
      call ghostx ( dxuij,gostx,im,jm,1,1,lattice,'west' )
      DO j = 1,JM
      DO I = 1,IM
       DXPIJ(I,j) = HALF*( DXUIJ(I,j)+gostx(i-1,j) )
      ENDDO
      ENDDO

      call ghosty ( dxpij,gosty,im,jm,1,0,1,1,lattice,'north' )
      DO j = p2,jmm2
      DO I = 1,IM
       DXVIJ(I,j) = HALF*( gosty(I,j+1)+DXPIJ(I,j) )
      ENDDO
      ENDDO

      if( lattice%pej.eq.0 ) then
      DO I = 1,IM
       DXVIJ(I,1) = ZERO
      ENDDO
      endif
      if( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       DXVIJ(I,JM-1) = ZERO
      ENDDO
      endif


c DYU & DYP
c ---------
      call ghosty ( dyvij,gosty,im,jm,1,1,1,1,lattice,'south' )
      DO j = p3,jmm2
      DO I = 1,IM
       DYPIJ(I,j) = HALF*( gosty(I,j-1)+DYVIJ(I,j) )
      ENDDO
      ENDDO

      if( lattice%pej.eq.0 ) then
      DO I = 1,IM
       DYPIJ(I,1) = ZERO
       DYPIJ(I,2) = AE*( DPHI(1) + HALF*DPHI(2) )
      ENDDO
      endif
      if( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       DYPIJ(I,JM-1) = AE*( DPHI(lattice%jmglobal-1) + HALF*DPHI(lattice%jmglobal-2) )
       DYPIJ(I,JM  ) = ZERO
      ENDDO
      endif

      DO j = 1,JM
      DO I = 1,IM
       DYUIJ(I,j) = DYPIJ(I,j)
      ENDDO
      ENDDO


c AREA FACTORS
c ------------
      call ghosty ( dxvij,gosty,im,jm,1,1,1,1,lattice,'south' )
      DO j = p3,jmm2
      DO I = 1,IM
       D2PIJ(I,j) = HALF*( gosty(I,j-1)+DXVIJ(I,j) ) * DYPIJ(I,j)
       D2UIJ(I,j) = DXUIJ(I,j) * DYUIJ(I,j)
      ENDDO
      ENDDO

      if( lattice%pej.eq.0 ) then
          DO I = 1,IM
          D2PIJ(I,2) = DXPIJ(I,2) * DYVIJ(I,1)
          D2UIJ(I,2) = D2PIJ(I,2)
          ENDDO
      endif
      if( lattice%pej.eq.lattice%ny-1 ) then
          DO I = 1,IM
          D2PIJ(I,JM-1) = DXPIJ(I,JM-1) * DYVIJ(I,JM-1)
          D2PIJ(I,JM  ) = ZERO
          D2UIJ(I,JM-1) = D2PIJ(I,JM-1)
          D2UIJ(I,JM  ) = ZERO
          ENDDO
      endif

      call ghostx (  d2pij(1:im,  1:jm),gostxy(0:im+1,1:jm  ),im,  jm,1,    1,lattice,'east'  )
      call ghosty ( gostxy(0:im+1,1:jm),gostxy(0:im+1,0:jm+1),im+2,jm,1,0,1,1,lattice,'north' )
 
      DO j = p2,jmm2
      DO I = 1,IM
       D2VIJ(I,j) = DXVIJ(I,j) * DYVIJ(I,j)
       D2ZIJ(I,j) = FOURTH*(  gostxy(i  ,j)+gostxy(i  ,j+1)
     *                      + gostxy(i+1,j)+gostxy(i+1,j+1) )
      ENDDO
      ENDDO

      if( lattice%pej.eq.0 ) then
          call par_sum ( d2pij(1,2),sum,im,1,lattice )
          DO I = 1,IM
          D2VIJ(I,1) = ZERO
          D2ZIJ(I,1) = (HALF/FLOAT( lattice%imglobal ))*sum
          ENDDO
      endif
      if( lattice%pej.eq.lattice%ny-1 ) then
          call par_sum ( d2pij(1,jm-1),sum,im,1,lattice )
          DO I = 1,IM
          D2VIJ(I,JM-1) = ZERO
          D2ZIJ(I,JM-1) = (HALF/FLOAT( lattice%imglobal ))*sum
          ENDDO
      endif


c CORIOLIS PARAMETER
c ------------------
      call ghosty ( dxuij,gosty,im,jm,1,0,1,1,lattice,'north' )
      DO J = p2,jmm2
      DO I = 1,IM
       FS = (CC( lattice%jglobal(j)+1 )*gosty(i,j+1)-CC( lattice%jglobal(j) )*DXUIJ(I,J))
       FC = ( DYVIJ(I,J)*(SL(lattice%iglobal(i)+1)-SL(lattice%iglobal(i))) 
     .            + HALF*(CL(lattice%iglobal(i)+1)+CL(lattice%iglobal(i)))
     .            * (gosty(i,j+1)*SC( lattice%jglobal(j)+1 )-DXUIJ(I,J)*SC( lattice%jglobal(j) )) )
       FFFIJ(I,J) = -(AE*OMEGA/D2ZIJ(I,J)) * (SINPNP*FS - COSPNP*FC)
      ENDDO
      ENDDO

      if( lattice%pej.eq.0 ) then
          call par_sum ( dxuij(1,2),sum,im,1,lattice )
          FS = CC(2)*(ONE/FLOAT( lattice%imglobal ))*sum
          call par_dot ( dxuij(1,2),cl(lattice%iglobal(1)),dot,im,1,lattice )
          FC = SC(2)*(ONE/FLOAT( lattice%imglobal ))*dot
          DO I = 1,IM
          FFFIJ(I,1) = -(AE*OMEGA/D2ZIJ(I,1))*(SINPNP*FS - COSPNP*FC)
          ENDDO
      endif

      if( lattice%pej.eq.lattice%ny-1 ) then
          call par_sum ( dxuij(1,jm-1),sum,im,1,lattice )
          FS = CC(lattice%jmglobal-1)*(ONE/FLOAT( lattice%imglobal ))*sum
          call par_dot ( dxuij(1,JM-1),cl(lattice%iglobal(1)),dot,im,1,lattice )
          FC = SC(lattice%jmglobal-1)*(ONE/FLOAT( lattice%imglobal ))*dot
          DO I = 1,IM
          FFFIJ(I,JM-1) = (AE*OMEGA/D2ZIJ(I,JM-1))*(SINPNP*FS - COSPNP*FC)
          ENDDO
      endif


c PRE-COMPUTE INVERSES
c --------------------
      DO J = p2,jmm1
      DO I = 1,IM
       DXUIN(I,J) = ONE / DXUIJ(I,J)
       D2PIN(I,J) = ONE / D2PIJ(I,J)
      ENDDO
      ENDDO

      DO J = p2,jmm2
      DO I = 1,IM
       DYVIN(I,J) = ONE / DYVIJ(I,J)
      ENDDO
      ENDDO

      DO J = p1,jmm1
      DO I = 1,IM
       D2ZIN(I,J) = ONE / D2ZIJ(I,J)
      ENDDO
      ENDDO


C AVRX Filter Arrays for the Stretched-Grid
c -----------------------------------------
c     if( stretched ) then
c         if( lattice%myid.eq.0 ) then
c         print *
c         print *, 'Computing Convolution Weights for Stretched-Grid'
c         print *, '------------------------------------------------'
c         print *
c         endif
c         call weights_p ( wp,dlam,dphi,lattice%imglobal,lattice%jmglobal,lattice )
c         call weights_u ( wu,dlam,dphi,lattice%imglobal,lattice%jmglobal,lattice )
c         call weights_v ( wv,dlam,dphi,lattice%imglobal,lattice%jmglobal,lattice )
c     else

C AVRX Filter Arrays for the Uniform-Grid
c ---------------------------------------
c         if( lattice%myid.eq.0 ) then
c         print *
c         print *, 'Computing Spectral Filters for Uniform-Grid'
c         print *, '-------------------------------------------'
c         print *
c         endif
          DO j = 1,JM
          DO I = 1,lattice%imglobal+2
           SU(I,j) = ONE
           SV(I,j) = ONE
          ENDDO
          ENDDO
          FC = COS(45.*PI/180.)
          DO J = p2,jmm1
          DO I = 3,lattice%imglobal+2
            SU(I,J) = MIN ( (CC(lattice%jglobal(j))/FC)/SIN(INT((I-1)/2)*PI/FLOAT(lattice%imglobal)),ONE  ) ** 2
          ENDDO
          ENDDO
          DO J = p2,jmm2
          DO I = 3,lattice%imglobal+2
            SV(I,J) = MIN ( (CV(lattice%jglobal(j))/FC)/SIN(INT((I-1)/2)*PI/FLOAT(lattice%imglobal)),ONE  ) ** 2
          ENDDO
          ENDDO
c     endif

      return
      end

      subroutine gridv ( LM,SIG,DSG )
      IMPLICIT NONE
      INTEGER  LM
      REAL(kind=8) SIG(LM+1)
      REAL(kind=8) DSG(LM)
      INTEGER L
      DO L=1,LM
       DSG(L) =  SIG(L+1) - SIG(L)
      ENDDO
      RETURN
      END

      subroutine avrx ( U,IM,JM,LM,S,lattice )
      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      include 'mpif.h'
      integer  status(mpi_status_size)
      integer   stats(mpi_status_size,0:lattice%nx-1)
      integer   statr(mpi_status_size,0:lattice%nx-1)

      integer  sendquest(0:lattice%nx-1)
      integer  recvquest(0:lattice%nx-1)
      integer  ierror,peid

      logical first,  flag
      data    first /.true./

      INTEGER    FORWARD,    BACKWARD
      PARAMETER (FORWARD=-1, BACKWARD=1)

      REAL(kind=8) ONE
      REAL(kind=8) ZERO
      PARAMETER ( ONE=1.0)
      PARAMETER (ZERO=0.0)

c ARGUMENTS
c ---------
      INTEGER IM,JM,LM
      REAL(kind=8)  U(IM,JM,LM)
      REAL(kind=8)  S(lattice%imglobal+2,JM)

c SCRATCH SPACE
c -------------
      real(kind=8),    allocatable :: sendbuf(:,:)
      real(kind=8),    allocatable :: recvbuf(:,:)
      real(kind=8),    allocatable ::       z(:,:)
      real(kind=8),    allocatable ::       b(:,:)
      integer, allocatable ::      j2(:)

      INTEGER I,J,L, J1(JM*LM), L1(JM*LM)
      INTEGER NFFTS
      INTEGER IM0
      DATA    IM0/0/

      integer n,num,rem,isum,lsum,len,len0
      integer p1,p2,p3,jmm1,jmm2,jmm3
      integer,           save :: ix(19)
      real(kind=8), allocatable, save :: tr(:)

      if (lattice%imglobal.ne.im0) then
        if(im0.ne.0) deallocate ( tr )
                       allocate ( tr(lattice%imglobal*2) )
        call fftfax (lattice%imglobal,ix,tr)
        im0=lattice%imglobal
      endif

c Set J-Index Range
c -----------------
      if ( lattice%pej.eq.0 ) then
           p1 = 1
           p2 = 2
           p3 = 3
      else
           p1 = 1
           p2 = 1
           p3 = 1
      endif
      if ( lattice%pej.eq.lattice%ny-1 ) then
           jmm1 = jm-1
           jmm2 = jm-2
           jmm3 = jm-3
      else
           jmm1 = jm
           jmm2 = jm
           jmm3 = jm
      endif

c Compute Number of FFTs to Perform and Load Data into Buffer
c -----------------------------------------------------------
      allocate ( sendbuf(im,jm*lm) )

      nffts = 0
      do j=p1,jmm1
        if( s(lattice%imglobal+2,j).lt.0.9999 ) then
          do L=1,lm
             nffts  = nffts + 1
          j1(nffts) = j
          L1(nffts) = L
             do i=1,im
             sendbuf(i,nffts) = u(i,j,L)
             enddo
          enddo
        endif
      enddo

      num = nffts/lattice%nx
      rem = nffts-lattice%nx*num

                                 len0 = num    ! Define number of FFTs on myid
      if( lattice%pei.le.rem-1 ) len0 = num+1  ! Define number of FFTs on myid

#if 0
      if( first ) then
          if( lattice%myid.eq.0 ) then
          print *, 'FFT Load Balance for Upper-Air Field:'
          print *, '-------------------------------------'
          endif
            do n=0,lattice%nx*lattice%ny-1
            if( n.eq.lattice%myid ) then
            write(6,1000) lattice%myid,lattice%pei,lattice%pej,nffts,num,rem,len0
            if( mod(n+1,lattice%nx).eq.0 ) print *
            endif
            call my_barrier (lattice%comm)
            enddo
            if( lattice%myid.eq.lattice%nx*lattice%ny-1 ) print *
      first = .false.
      endif
 1000 format(1x,'absolute PE id: ',i3,'  relative (pei,pej): ',i2,',',i2,'  nffts: ',i6,2x,
     .          'num: ',i6,2x,'rem: ',i6,2x,'len0: ',i6)
#endif

c Distribute Data Across PEs in X-direction
c -----------------------------------------
      if( nffts.ne.0 ) then

      allocate ( j2(len0) )
      allocate ( z(lattice%imglobal+2,len0  ) )
      allocate ( b(lattice%imglobal+2,len0*2) )

      isum = 0
      lsum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
                       len = num
      if( n.le.rem-1 ) len = num+1
              sendquest(n) = mpi_request_null
      if( len.ne.0   ) then
          if( peid.ne.lattice%myid ) then

              call mpi_isend ( sendbuf(1,1+lsum),im*len,mpi_double_precision,peid,peid,lattice%comm,sendquest(n),ierror )

          else
              do L=1,len0
                 j2(L) = j1(L+lsum)
                 do i=1,im
                 z(i+isum,L) = sendbuf(i,L+lsum)
                 enddo
              enddo
          endif
      endif
      isum = isum + lattice%im(n)
      lsum = lsum + len
      enddo

c Receive Data and Perform FFT
c ----------------------------
      if( len0.ne.0 ) then

      allocate ( recvbuf( lattice%imglobal*len0,0:lattice%nx-1 ) )
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
      if( peid.ne.lattice%myid ) then

          call mpi_irecv ( recvbuf(1,n),lattice%im(n)*len0,mpi_double_precision,peid,lattice%myid,lattice%comm,
     .                     recvquest(n),ierror )

      else
          recvquest(n) = mpi_request_null
      endif
      enddo

      call mpi_waitall ( lattice%nx,sendquest(0:lattice%nx-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx,recvquest(0:lattice%nx-1),statr(1,0),ierror )

      isum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
      if( peid.ne.lattice%myid ) then
           lsum = 0
           do L=1,len0
           do i=1,lattice%im(n)
           lsum = lsum+1
            z(i+isum,L) = recvbuf(lsum,n)
           enddo
           enddo
      endif
      isum = isum + lattice%im(n)
      enddo

      do L=1,len0
      z(lattice%imglobal+1,L) = 0.0
      z(lattice%imglobal+2,L) = 0.0
      enddo

c Perform FFT
c -----------
      call rfftmlt( Z,B,TR,IX,1,lattice%imglobal+2,lattice%imglobal,len0,FORWARD )
      do L=1,len0
      do i=1,lattice%imglobal+2
       z(i,L) = z(i,L)*s(i,j2(L))
      enddo
      enddo
      call rfftmlt( Z,B,TR,IX,1,lattice%imglobal+2,lattice%imglobal,len0,BACKWARD )

c Distribute Filtered Data Back to source PEs
c -------------------------------------------
      deallocate ( sendbuf )
      deallocate ( recvbuf )
        allocate ( sendbuf( lattice%imglobal*len0,0:lattice%nx-1 ) )
      isum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
      if( peid.ne.lattice%myid ) then

          lsum = 0
          do L=1,len0
          do i=1,lattice%im(n)
          lsum = lsum+1
          sendbuf(lsum,n) = z(i+isum,L)
          enddo
          enddo
          call mpi_isend ( sendbuf(1,n),lsum,mpi_double_precision,peid,peid,lattice%comm,sendquest(n),ierror )

      else
          sendquest(n) = mpi_request_null
      endif
      isum = isum + lattice%im(n)
      enddo

      endif  ! End len0.ne.0 check

c Receive Filtered Data
c ---------------------
      allocate ( recvbuf(im,nffts) )
      isum = 0
      lsum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
                       len = num
      if( n.le.rem-1 ) len = num+1
              recvquest(n) = mpi_request_null
      if( len.ne.0   ) then
          if( peid.ne.lattice%myid ) then

              call mpi_irecv ( recvbuf(1,1+lsum),im*len,mpi_double_precision,peid,lattice%myid,lattice%comm,recvquest(n),ierror )

          else
              do L=1,len0
              do i=1,im
                recvbuf(i,L+lsum) = z(i+isum,L)
              enddo
              enddo
          endif
      endif
      isum = isum + lattice%im(n)
      lsum = lsum + len
      enddo

      call mpi_waitall ( lattice%nx,sendquest(0:lattice%nx-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx,recvquest(0:lattice%nx-1),statr(1,0),ierror )

c Reconstruct Filtered Field
c --------------------------
      do n=1,nffts
      do i=1,im
      u(i,j1(n),L1(n)) = recvbuf(i,n)
      enddo
      enddo

      deallocate ( z,b,j2  )
      deallocate ( sendbuf )
      deallocate ( recvbuf )
      else
      deallocate ( sendbuf )
      endif

      return
      end

      subroutine hadvect (U,V,P,IM,JM,DIV,VT1,VT2,lattice)
      use g3_dynamics_lattice_module
      implicit none
      type ( dynamics_lattice_type ) lattice

      INTEGER  IM,JM
      REAL(kind=8)   P(IM,JM)

      REAL(kind=8)   U(0:im+1,0:jm+1)
      REAL(kind=8)   V(0:im+1,0:jm+1)
      REAL(kind=8)   Q(0:im+1,0:jm+1)
      REAL(kind=8) DIV(0:im+1,0:jm+1)
      REAL(kind=8) VT1(0:im+1,0:jm+1)
      REAL(kind=8) VT2(0:im+1,0:jm+1)

      REAL(kind=8) ST1
      REAL(kind=8) ST2

      integer i,j,p1,p2,p3,jmm1,jmm2,jmm3

c Set J-Index Range
c -----------------
      if ( lattice%pej.eq.0 ) then
           p1 = 1
           p2 = 2
           p3 = 3
      else
           p1 = 1
           p2 = 1
           p3 = 1
      endif
      if ( lattice%pej.eq.lattice%ny-1 ) then
           jmm1 = jm-1
           jmm2 = jm-2
           jmm3 = jm-3
      else
           jmm1 = jm
           jmm2 = jm
           jmm3 = jm
      endif

      call ghostx ( p(1:im,  1:jm),q(0:im+1,1:jm  ),im,  jm,1,    1,lattice,'both' )
      call ghosty ( q(0:im+1,1:jm),q(0:im+1,0:jm+1),im+2,jm,1,0,1,1,lattice,'both' )

c CORNER FLUXES
c -------------
      DO j = p2,jmm1
      DO I = 1,IM
       DIV(I,j) = (U(I,j)-U(i-1,j)) - (V(I,j)-V(I,j-1))   ! Mass-Point
      ENDDO
      ENDDO

      call ghostx ( div(1:im,  1:jm),div(0:im+1,1:jm  ),im,  jm,1,    1,lattice,'east' )
      call ghosty ( div(0:im+1,1:jm),div(0:im+1,0:jm+1),im+2,jm,1,0,1,1,lattice,'both' )

      DO j = p2,jmm2
      DO I = 1,IM
       VT1(I,j) = 0.5*( Q(i,j  )+Q(i+1,j+1) ) * ( DIV(i  ,j+1)-DIV(i+1,j) )   ! VWND-Point
       VT2(I,j) = 0.5*( Q(i,j+1)+Q(i+1,j  ) ) * ( DIV(i+1,j+1)-DIV(i  ,j) )   ! VWND-Point
      ENDDO
      ENDDO

c Corner Flux at North Pole
c -------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       VT1(I,JM-1) = 0.   ! VWND-Point at Pole
       VT2(I,JM-1) = 0.   ! VWND-Point at Pole
      ENDDO
      endif

c Corner Flux at South Pole
c -------------------------
      if ( lattice%pej.eq.0 ) then
      DO I = 1,IM
       VT1(I,1) = 0.   ! VWND-Point at Pole
       VT2(I,1) = 0.   ! VWND-Point at Pole
      ENDDO
      endif

      call ghostx ( vt1(1:im,  1:jm),vt1(0:im+1,1:jm  ),im,  jm,1,    1,lattice,'west'  )
      call ghosty ( vt1(0:im+1,1:jm),vt1(0:im+1,0:jm+1),im+2,jm,1,1,1,1,lattice,'south' )
      call ghostx ( vt2(1:im,  1:jm),vt2(0:im+1,1:jm  ),im,  jm,1,    1,lattice,'west'  )
      call ghosty ( vt2(0:im+1,1:jm),vt2(0:im+1,0:jm+1),im+2,jm,1,1,1,1,lattice,'south' )

      DO j = p2,jmm1
      DO I = 1,IM
       DIV(I,j) = (1./48.)*( (VT1(i  ,j) - VT1(i-1,j-1))
     *                     + (VT2(i-1,j) - VT2(i  ,j-1)) )   ! Mass-Point
      ENDDO
      ENDDO
 
c NEAR FLUXES
c -----------
      DO j = p3,jmm2
      DO I = 1,IM
       ST1 = 28.*U(I,j) + 3.*(U(i+1,j)+U(i-1,j)) - (U(I,j-1)+U(I,j+1))
       VT1(I,j) = 0.5*(Q(i+1,j) + Q(i,j))*ST1   ! UWND-Point
      ENDDO
      ENDDO
      DO j = p3,jmm3
      DO I = 1,IM
       ST2 = 28.*V(I,j) + 3.*(V(I,j+1)+V(I,j-1)) - (V(i+1,j)+V(i-1,j))
       VT2(I,j) = 0.5*(Q(I,j+1) + Q(I,j))*ST2   ! VWND-Point
      ENDDO
      ENDDO
 
 
c Near Flux Near North Pole
c -------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       ST1 = 0.5*(U(i+1,JM-1)+U(i-1,JM-1)-V(I,JM-2)+V(i+1,JM-2))
       ST1 = 28.*U(I,JM-1) + 3.*(U(i+1,JM-1)+U(i-1,JM-1)) - (ST1 + U(I,JM-2))
       VT1(I,JM-1) = 0.5*(Q(i+1,JM-1) + Q(i,JM-1))*ST1                         ! UWND-Point Near Pole

       ST2 = 26.*V(I,JM-2) + 3.*(V(I,JM-3)) - (V(i+1,JM-2)+V(i-1,JM-2))
       VT2(I,JM-2) = 0.5*(Q(I,JM-1) + Q(I,JM-2))*ST2                           ! VWND-Point Near Pole
       VT2(I,JM-1) = 0.                                                        ! VWND-Point at   Pole
      ENDDO
      endif

c Near Flux Near South Pole
c -------------------------
      if ( lattice%pej.eq.0 ) then
      DO I = 1,IM
       ST1 = 0.5*(U(i+1,2)+U(i-1,2)+V(I,2)-V(i+1,2))
       ST1 = 28.*U(I,2) + 3.*(U(i+1,2)+U(i-1,2)) - (ST1 + U(I,3))
       VT1(I,2) = 0.5*(Q(i+1,2) + Q(i,2))*ST1                                  ! UWND-Point Near Pole

       ST2 = 26.*V(I,2) + 3.*(V(I,3)) - (V(i+1,2)+V(i-1,2))
       VT2(I,2) = 0.5*(Q(I,3) + Q(I,2))*ST2                                    ! VWND-Point Near Pole
       VT2(I,1) = 0.                                                           ! VWND-Point at Pole
      ENDDO
      endif

      call ghostx ( vt1(1:im,1:jm),vt1(0:im+1,1:jm  ),im,jm,1,     1,lattice,'west'  )
      call ghosty ( vt2(1:im,1:jm),vt2(1:im  ,0:jm+1),im,jm,1,1,-1,1,lattice,'south' )

      DO j = p2,jmm1
      DO I = 1,IM
       DIV(I,j) = DIV(I,j) + (1./24.)*( (VT1(I,j) - VT1(i-1,j))
     *                                + (VT2(I,j) - VT2(i,j-1)) )   ! Mass-Point
      ENDDO
      ENDDO

c FAR FLUXES
c ----------
      DO j = p2,jmm1
      DO I = 1,IM
       ST1 = U(I,j) + U(i-1,j)
       VT1(I,j) = 0.5*(Q(i+1,j) + Q(i-1,j))*ST1   ! Mass-Point
      ENDDO
      ENDDO
      DO j = p3,jmm2
      DO I = 1,IM
       ST2 = V(I,j) + V(i,j-1)
       VT2(I,j) = 0.5*(Q(i,j+1) + Q(i,j-1))*ST2   ! Mass-Point
      ENDDO
      ENDDO
 
c Far Flux Near North Pole
c ------------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       VT2(I,JM-1) = 0.   ! Mass-Point Near Pole
      ENDDO
      endif

c Far Flux Near South Pole
c ------------------------
      if ( lattice%pej.eq.0 ) then
      DO I = 1,IM
       VT2(I,2) = 0.      ! Mass-Point Near Pole
      ENDDO
      endif

      call ghostx ( vt1(1:im,1:jm),vt1(0:im+1,1:jm  ),im,jm,1,     1,lattice,'both' )
      call ghosty ( vt2(1:im,1:jm),vt2(1:im  ,0:jm+1),im,jm,1,0,-1,1,lattice,'both' )
 
      DO j = p3,jmm2
      DO I = 1,IM
       DIV(I,j) = DIV(I,j) - (1./12.)*( (VT1(i+1,j) - VT1(i-1,j))
     *                                + (VT2(i,j+1) - VT2(i,j-1)) )   ! Mass-Point
      ENDDO
      ENDDO
 
c DIV Near North Pole
c -------------------
      if ( lattice%pej.eq.lattice%ny-1 ) then
      DO I = 1,IM
       DIV(I,JM-1) = DIV(I,JM-1) - (1./12.)*( VT1(i+1,JM-1)-VT1(i-1,JM-1)-VT2(i,JM-2) )   ! Mass-Point Near Pole
      ENDDO
      endif

c DIV Near South Pole
c -------------------
      if ( lattice%pej.eq.0 ) then
      DO I = 1,IM
       DIV(I,2) = DIV(I,2) - (1./12.)*( VT1(i+1,2)-VT1(i-1,2)+VT2(i,3) )   ! Mass-Point Near Pole
      ENDDO
      endif

      RETURN
      END
