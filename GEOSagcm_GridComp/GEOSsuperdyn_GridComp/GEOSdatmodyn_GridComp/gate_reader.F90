           sUBROUTINE GATE(NT, LM,                         &
                     PREF_IN,                            &
                     time,                               &
                     yy,                                 &
                     mo,                                 &
                     dd,                                 &
                     hh,                                 &
                     mm,                                 &
                     PCP_OBS,                            &
                     TS_AIR,                             &
                     TG_SOIL,                            &
                     TSKIN,                              &
                     QSFCAIR,                            &
                     QSKIN,                              &
                     PSFC,                               &
                     LHF,                                &
                     SHF,                                &
                     tt,                                 &
                     qq,                                 &
                     uu,                                 &
                     vv,                                 &
                     T_H_adv,                            &
                     T_V_adv,                            &
                     Q_H_adv,                            &
                     Q_V_adv,                            &
                     Q1,                                 &
                     Q2,                                 & 
                     P_MODEL_E                            )

        use GEOS_Mod
 

       IMPLICIT NONE

!
!  lm is number of sigma layers in SCM
! ARGUMENTS
      INTEGER UNIT, LMM1   ! MOdel
!     INTEGER, parameter :: LMP1 = LM + 1
      INTEGER, parameter ::  NB=16, NTM1= 159  !
!     INTEGER, parameter ::  NT=160
! INPUT

      INTEGER,                       INTENT(IN   ) :: NT, LM
      REAL,  DIMENSION(0:LM),     INTENT(IN   ) :: PREF_IN
      
      
! OUTPUT
! Singli-layer:
      real, DIMENSION(nt       ), INTENT(  OUT) :: time        !Calenday day
      real, DIMENSION(nt       ), INTENT(  OUT) ::  yy        !Year
      real, DIMENSION(nt       ), INTENT(  OUT) ::  mo        !Month
      real, DIMENSION(nt       ), INTENT(  OUT) ::  dd        !Day
      real, DIMENSION(nt       ), INTENT(  OUT) ::  hh        !Hour
      real, DIMENSION(nt       ), INTENT(  OUT) ::  mm        !Minutes
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: PCP_OBS
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TS_AIR 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TG_SOIL 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: PSFC, LHF, SHF
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: TSKIN
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: QSKIN 
      REAL, DIMENSION(NT       ), INTENT(  OUT) :: QSFCAIR 

! Multiple-layer

      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: tt
      REAL, DIMENSION(NT, LM ), INTENT(  OUT) :: qq
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: uu
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: vv
!     REAL, DIMENSION(NT, LM), INTENT(  OUT) :: ww
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: T_H_adv
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: T_V_adv
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: Q_H_adv
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: Q_V_adv
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: Q1
      REAL, DIMENSION(NT, LM), INTENT(  OUT) :: Q2

      REAL, DIMENSION(NT, 0:LM), INTENT(  OUT) :: P_MODEL_E


! LOCAL


      REAL,    PARAMETER :: CP=1003.5, GRAV=9.81, RKAP=0.286
      REAL,    PARAMETER :: PI=3.1415926535898, AE=6.371E6
      REAL,    PARAMETER :: HEATL=2.52E6
      REAL     GRAV2
     

!  gate data: 38 pressure levels, 160 time periods, 16 variables
!  starting time, 0Z30aug1974, interval 3hr. gate region centered 9N, 23.5W

      INTEGER, PARAMETER :: nlayr = 38       ! gate pressure levels
      INTEGER, PARAMETER :: NP=41  !! Layers og GCM 
      INTEGER, PARAMETER :: NVAR3=16         !, NVAR2=16
       real, DIMENSION(LM) :: p_model 
      integer, parameter :: ipp = 3
      integer ilv, itim

      REAL, DIMENSION(nlayr, NT,NVAR3) :: FLD    ! (38,160,16) 
      INTEGER, DIMENSION(ipp    , NT   ) :: IME    ! (3,160)
      REAL, DIMENSION(nlayr           ) :: PRL    ! (38), 
      REAL, DIMENSION(nlayr, NB       ) :: avg    ! (38,16)
      REAL, DIMENSION(nlayr           ) :: DQTRMS ! (38)     
      REAL, DIMENSION(nlayr           ) :: DTTRMH ! (38) 
      REAL, DIMENSION(nlayr           ) :: DTTRM2 ! (38)
!
      REAL, DIMENSION(NT, LM)           :: TEMP   ! (NT,LM), 
      REAL, DIMENSION(NT, LM)           :: Q      ! (NT,LM), 
      REAL, DIMENSION(NT, LM)           :: u      ! (nt,lm), 
      REAL, DIMENSION(NT, LM)           :: v      ! (nt,lm)
      REAL, DIMENSION(NT)               :: PS     ! (NT), 
!     REAL, DIMENSION(LMP1)             :: SGE    ! (LMP1), 
!     REAL, DIMENSION(LM)               :: SGO    ! (LM), 
!     REAL, DIMENSION(LM)               :: DSG    ! (LM)
!     REAL, DIMENSION(NT)               :: SGEOP  ! (NT)
      REAL, DIMENSION(NT, LM)           :: VDQ    ! (NT,LM), 
      REAL, DIMENSION(NT, LM)           :: VDT    ! (NT,LM), 
      REAL, DIMENSION(NT, LM)           :: VDQH   !
      REAL, DIMENSION(NT, LM)           :: div
      REAL, DIMENSION(NT, LM)           :: omega
      REAL, DIMENSION(NT, LM)           :: VDQV
      REAL, DIMENSION(NT, LM)           :: VDTH
      REAL, DIMENSION(NT, LM)           :: VDTV
!
      REAL, DIMENSION(NT)               :: PRE1   ! (NT),   
      REAL, DIMENSION(NT)               :: PRE2   ! (NT),   
      REAL, DIMENSION(NT, LM)           :: QR     ! (NT,LM), 
      REAL, DIMENSION(NT, LM)           :: VDT2   ! (NT,LM)
      REAL, DIMENSION(NT    )           :: DQDTB  ! (NT),  
      REAL, DIMENSION(NT    )           :: PRE3   ! (NT)
      REAL, DIMENSION(NT    )           :: EVAP   ! (NT),   
      REAL, DIMENSION(NT    )           :: SST
!
      REAL, DIMENSION(NT, nlayr)        :: QRX    ! (NT,38), 
      REAL, DIMENSION(NT, nlayr)        :: QR1    ! (NT,38), 
      REAL, DIMENSION(NT, nlayr)        :: QR2    ! (NT,38)
      REAL, DIMENSION(NT, nlayr)        :: HBR    ! (NT,38), 
      REAL, DIMENSION(NT, nlayr)        :: QBR    ! (NT,38)
      REAL, DIMENSION(NT, nlayr)        :: UQX    ! (NT,38), 
      REAL, DIMENSION(NT, nlayr)        :: VQY    ! (NT,38), 
      REAL, DIMENSION(NT, nlayr)        :: WQP    ! (NT,38)
      REAL, DIMENSION(NT    )           :: PRQ1   ! (NT),   
      REAL, DIMENSION(NT    )           :: SENH   ! (NT), 
      REAL, DIMENSION(NT    )           :: EVAP2  ! (NT), 
      REAL, DIMENSION(NT    )           :: SENH2  ! (NT)

      REAL, DIMENSION(NLAYR )           :: DTF    ! (38),  
      REAL, DIMENSION(NLAYR )           :: DQF    ! (38), ! TX1(50), TX2(50)
      REAL, DIMENSION(NLAYR )           :: DTH    ! 
      REAL, DIMENSION(NLAYR )           :: DTV    !
      REAL, DIMENSION(NLAYR )           :: DQV    !
      REAL, DIMENSION(NLAYR )           :: DQH    !
      REAL, DIMENSION(NLAYR )           :: fake
  !   REAL, DIMENSION(NLAYR )           :: VV
      REAL, DIMENSION(NLAYR )           :: DQTRMV
      REAL, DIMENSION(NLAYR )           :: DQTRMH
      REAL, DIMENSION(NLAYR )           :: DTTRMV
      REAL, DIMENSION(NLAYR )           :: DTTRMS
      REAL, DIMENSION(NLAYR )           :: PLEVS_GATE

      
      
      
      real UTTRM, VTTRM, WETRM, WQTRM, VQTRM, UQTRM, DPMX, PIFAC,PTOP
      real ELGRAV, CPOR, RHOAIR, TEM, TEMU, TEMD, PRESK, PDWNK
      real PDWN, PUPP, PRGAT, dummy, RGAS 
      integer IC, n, IGTLEV, IGU, IGD, KL, KR, LR  
      real DPMY, DELP 
      real QZERO, PRES, WSTRM, FAC1, FAC2, SUM1, PUPPK  

      real TDEGK, WNDMAG, TEM2, TEM1, QDIFF, SUM2

      real CSW(10,20), CLW(10,20), PRSR(10)

      REAL * 8 PQ1, PQ2, PQA, DTEM, Q1MQ2(NT), QRBR(NT), Q1MQ2M, QRBM

      character*16 desc
      data desc /'GATE Data       '/
      real, parameter ::  xlat = 9. 
      real, parameter ::  xlon = -23.5
      integer L, ND, NL, JJ, I, K, IFIELD
!
!  SPE!IFY PHOENIX SIGMA EDGE VALUES
!  ---------------------------------
!     DATA SGE /0.000000,0.025000,0.070000,0.124000,0.184000,0.248000, &
!              0.316000,0.387500,0.462000,0.539000,0.617779,0.696448,  &
!              0.772512,0.843153,0.905120,0.954730,0.987871,1.000000/

!  the following is LW & SW forcing terms. i don't know the source
!  looks like daily values (20 days), on the 10 pressure levels

      DATA PRSR/150.,250.,350.,450.,550.,650.,750.,850.,950.,1006./
      DATA CSW/3.1,6.8,11.9,11.6,10.8,8.8,5.9,4.2,2.6,0.1,   &
               2.8,5.8,11.,11.7,12.0,11.,9.2,7.2,5.1,0.2,   &
               2.8,5.8,11.2,11.6,11.,9.6,8.3,7.2,5.8,0.3,  &
               6.4,9.9,9.6,7.3,6.2,3.5,2.1,1.4,0.9,0.0,     &
               3.1,6.3,11.7,12.8,13.,8.3,4.9,3.,1.8,0.0,   &
               5.1,9.0,9.7,8.3,7.5,6.3,5.5,3.9,2.5,0.1,    &
               5.9,9.4,8.9,7.8,7.3,6.3,4.6,3.2,2.1,0.1,    &
               3.3,7.1,11.3,11.1,10.8,9.1,6.5,4.6,2.5,0.1, &
               2.9,6.2,11.4,11.9,12.2,11.4,7.8,5.2,3.2,0.1,  &
               3.4,7.4,11.0,10.8,10.7,8.8,6.6,4.4,3.2,0.1,  &
               3.4,8.7,11.8,9.9,8.9,7.2,4.3,2.3,1.2,0.,     &
               3.1,6.8,11.4,11.1,10.7,9.1,6.8,5.2,3.8,0.2,  &
               4.1,6.3,10.3,10.5,10.,8.9,7.7,6.1,4.9,0.3,   &
               5.6,8.,9.7,8.7,7.9,6.5,5.,4.,3.,0.2,         &
               3.4,7.6,10.6,10.3,9.9,8.2,6.5,5.3,3.9,0.2,   &
               3.6,7.7,11.1,10.6,10.2,8.1,5.5,3.7,2.7,0.1,    &
               2.9,5.9,11.,11.6,11.8,10.9,8.,6.1,5.2,0.3,    &
               4.4,8.5,10.8,9.4,8.4,6.4,4.6,3.4,2.,0.1,     &
               4.6,8.6,11.,9.1,7.3,5.8,4.3,3.0,1.9,0.1,      &
               2.9,5.9,10.9,11.4,11.2,10.1,8.7,7.9,5.9,0.3/

      DATA CLW/-9.0,-18.6,-28.3,-31.5,-28.7,-22.0,-17.9,-14.1,-8.7,-4.0,     &
            -5.9,-15.6,-24.1,-28.1,-32.4,-30.3,-30.1,-24.,-13.6,-4.0,     &
            -5.9,-15.6,-24.2,-26.,-29.2,-29.3,-32.4,-28.1,-15.6,-4.0,                &
            -12.7,-21.6,-26.5,-26.8,-23.3,-16.1,-14.3,-12.1,-8.3,-4.0,        &
            -6.8,-16.5,-25.3,-38.3,-35.6,-23.4,-19.9,-15.2,-9.9,-4.0,            &
            -11.8,-19.6,-24.7,-25.6,-24.,-21.8,-22.,-16.3,-9.4,-4.0,           &
            -11.5,-20.5,-24.9,-25.4,-25.9,-22.0,-19.,-13.8,-8.5,-4.0,         &
            -7.1,-17.2,-25.7,-29.0,-33.2,-29.1,-24.7,-17.3,-8.8,-3.9,         &
            -6.4,-16.5,-25.4,-29.,-37.5,-31.2,-23.7,-16.6,-9.2,-3.9,        &
            -7.6,-18.8,-26.1,-28.1,-31.9,-26.5,-21.9,-15.3,-10.,-4.0,        &
            -7.0,-20.2,-30.6,-30.4,-29.1,-22.9,-16.9,-12.3,-8.6,-4.0,          &
            -6.3,-17.,-26.9,-28.9,-30.4,-27.5,-25.3,-20.5,-12.4,-4.0,         &
            -9.2,-17.1,-24.3,-25.5,-26.2,-26.9,-27.,-21.2,-12.9,-4.1,      &
            -10.4,-18.5,-25.4,-25.2,-24.4,-23.7,-23.2,-18.6,-11.2,-4.,         &
            -8.3,-18.8,-24.9,-26.,-27.7,-24.7,-23.8,-20.3,-12.3,-4.1,           &
            -9.,-19.,-25.7,-28.1,-30.,-24.4,-20.2,-15.5,-9.7,-4.0,           &
            -6.5,-16.5,-24.6,-26.9,-36.2,-29.6,-24.4,-19.6,-13.4,-4.1,      &
            -10.6,-19.8,-25.9,-26.7,-27.1,-22.9,-19.2,-13.6,-8.4,-4.0,                 &
            -12.0,-19.0,-27.3,-26.9,-23.6,-21.4,-19.4,-13.6,-8.6,-4.0,          &
            -6.5,-16.2,-24.9,-26.9,-29.9,-29.9,-29.4,-24.7,-12.4,-3.9/

       RGAS=CP*RKAP
       GRAV2=GRAV/100.0
   !   LMP1=LM+1  
   !   LMM1=LM-1   ! MOdel

! 
! Average edge Pressures to layer-mids and convert to hPa

      !!p_model(1:LM) = ( p_model_e(0:LM-1) + p_model_e(1:LM) )/2.0/100.
      p_model(1:LM) = ( pref_in(0:LM-1) + pref_in(1:LM) )/2.0/100.

!  SET THE SIGMA LEVELS OF THE MODEL
!  ---------------------------------
!        DO L=1,LM
!        DSG(L)   = (SGE(L+1) - SGE(L))
!        SGO(L)   = (SGE(L+1) + SGE(L)) * 0.5
!     END DO

!     WRITE(6,1001) SGE
!     WRITE(6,1002) SGO
!     WRITE(6,1003) DSG

 1001 FORMAT(' SGE:',10E12.4)
 1002 FORMAT(' SGO:',10E12.4)
 1003 FORMAT(' DSG:',10E12.4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! no surface flux forcing data
            LHF = -9999.
            SHF = -9999.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      FAC1 = (GRAV/CP) * 0.0001
      FAC2 = GRAV/(CP*1200.0)
      DO  ND=1,20
      SUM1 = 0.0
      SUM2 = 0.0

      DO 6005 NL=1,10
      SUM1 = SUM1 + CSW(NL,ND)
      SUM2 = SUM2 + CLW(NL,ND)
 6005 CONTINUE
!
   !  WRITE(6,6006) ND, SUM1, SUM2
 6006 FORMAT(' VERT. MEAN AT DAY=',I2,' AND SW=',F6.1,' AND LW=',F7.1)
!
      DO 6007 NL=1,9
      CSW(NL,ND) = (CSW(NL,ND) + CLW(NL,ND)) * FAC1
 6007 CONTINUE
      CSW(10,ND) = (CSW(10,ND) + CLW(10,ND)) * FAC2

      ENDDO
!
       RHOAIR = 1.275
       CPOR   = CP/RGAS
       ELGRAV = HEATL*GRAV
       PTOP   = 10.
!
!   DEGREES LONGITUDE PER METER (  1./ (5.*PI/180*RADEARTH*COS(DEG.LAT))
!
!      DPMX   = .18E-5
       PIFAC = PI / 180.0
       DPMX = AE * PIFAC * COS(xlat*PIFAC)
       DPMX = 1.0 / DPMX
!
!   DEGREES LATUTUDE PER METER (  1./ (5.*PI/180*RADEARTH))
!
!      DPMY   = .22E-5
       DPMY = AE * PIFAC
       DPMY = 1.0 / DPMY
!
   !   WRITE(6,*)' DPMX=',DPMX,' DPMY=',DPMY
!
!  READ DATA
!  ---------

     UNIT = GETFILE('gateorig.data.ascii.80', form="formatted")

!     open(888,file='gateorig.data.ascii.80',form='formatted', status='old')  

!  note: i seem to recall that there might be something suspicious in the 
!        first 3 time periods of the forcing terms. check before using.

       DO JJ = 1,160
         SST (JJ) = 27.2 + 273.16   !  Frank add
         READ(UNIT,'(1X,I2,I3,I3)') (IME(I,JJ), I = 1,3)
         time(jj) = 999. 
         yy(jj) = 1974.
         mo(JJ) = IME(1,JJ)
         dd(jj) = IME(2,jj)
         hh(jj) = IME(3,jj)
         mm(jj) = 0.
   !     print *, 'time',yy(jj),dd(jj),mo(jj),hh(jj),mm(jj)
         DO IFIELD = 1,16
           READ(UNIT,'(5(1X,E15.8))')  (FLD(I,JJ,IFIELD), I = 1,38)
         end do
       end do

   call FREE_FILE(UNIT)

!  field  1 is temperature (C)
!  field  2 is specific humidity (g/kg)
!  field  3 is u wind (m/sec)
!  field  4 is v wind (m/sec)
!  field  5 is geopotential height (m)
!  field  6 is omega (1e-3 mb/sec)
!  field  7 is dT/dx advection term (K/deg lat)
!  field  8 is dT/dy advection term (K/deg lon)
!  field  9 is dT/dt tendency (K/3 hrs)
!  field 10 is dq/dx advection term (g/kg/deg lat)
!  field 11 is dq/dy advection term (g/kg/deg lon)
!  field 12 is dq/dt tendency (g/kg/3 hrs)
!  field 13 is dq/dp vertical advection term (g/kg/mb)
!  field 14 is dS/dp vertical advection term (1e3 J/kg/mb)
!  field 15 is Q1 (1e-2 J/kg/sec)
!  field 16 is Q2 (1e-2 J/kg/sec)

          DO IGTLEV=2,38
           PRGAT=1000.0 - (IGTLEV-2)*25.0
           PLEVS_GATE(IGTLEV) = PRGAT
          ENDDO

       DO L=1,160


         DO K = 1,38
           FLD(K,L,1) = FLD(K,L,1) + 273.16    ! T     CONVERT TO DEG K
           FLD(K,L,6) = FLD(K,L,6)*1.E-3       ! OMEGA CONVERT TO MB/SEC
           FLD(K,L,14) = FLD(K,L,14) * 1000.0  ! CONVERT DS/DP
           FLD(K,L,15) = FLD(K,L,15)*1.E-2     ! CONVERT Q1
           FLD(K,L,16) = FLD(K,L,16)*1.E-2     ! CONVERT Q2
         end do

!   SURFACE GEOPOTENTIAL FROM SURFACE HEIGHT (MULTIPLY BY GRAV).
!   OBTAIN SURFACE PRESSURE FROM TEMPERATURES, USING EQ. OF STATE,
!          HYDROSTATIC RELATION AND NEUTRAL LAPSE RATE
!                PS = 1000. * ( (T1000/TS)**-CP/R )
!   --------------------------------------------------------------
     !   SGEOP(L) = FLD(1,L,5)*GRAV

         PS(L) = 1000.*((FLD(2,L,1)/FLD(1,L,1))**(-CPOR))

! Some "surface air" things         

         TS_AIR(L)  = FLD(1,L,1)
         QSFCAIR(L) = FLD(1,L,2)
         PSFC(L)    = PS(L)

! "Skin" quantities
         TSKIN(L)   = 273.15 + 27.5 
         QSKIN(L)   = 1000.*GEOS_QSAT (TSKIN(L) , PSFC(L) )

         TG_SOIL(L) = -999.9  !No "soil"

! 
! Average edge Pressures to layer-mids and convert to hPa
         p_model_e(L,0:LM) = pref_in(0:LM)*(PS(L)*100./pref_in(LM))
         p_model(1:LM) = ( p_model_e(L,0:LM-1) + p_model_e(L,1:LM) )/2.0/100.
         PLEVS_GATE(1) = PS(L)
!
!  COMPUTE DQ AND DT DUE TO HORIZONTAL AND VERTICAL ADVECTION
!   ----------------------------------------------------------

         DO K = 1,38
           UQTRM  = FLD(K,L,3)*FLD(K,L,10)*DPMX
           VQTRM  = FLD(K,L,4)*FLD(K,L,11)*DPMY
           WQTRM  = FLD(K,L,6)*FLD(K,L,13)

           DQH(K) = UQTRM + VQTRM    ! 
           DQV(K) = WQTRM

           UTTRM  = FLD(K,L,3)*FLD(K,L,7)*DPMX
           VTTRM  = FLD(K,L,4)*FLD(K,L,8)*DPMY
           WSTRM  = FLD(K,L,6)*FLD(K,L,14) / CP

           DTH(K) = UTTRM + VTTRM
           DTV(K) = WSTRM
         end do

!  DELP IS TH PS-PTOP IN KG/SEC**2-M( MB*1.E2 = KG/SEC**2-M )
         DELP=(PS(L)-PTOP)*1.E2
!
!   OBTAIN THE GLAS PRESSURE AT DIFFERENT LEVELS FROM PS AND SIGMA
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

         DO KL=1,LM    ! LM: after pressure convert

!         PRES = PS(L)*SGO(KL)
!         PRES = 1000.0*SGO(KL)

          PRES = p_model(KL)

          IF (PRES .GT. 1000.0) THEN
           IGD  = 1
           IGU  = 2
           PUPP = 1000.
           PDWN = PS(L)
           GO TO 205
          ENDIF

          IF (PRES .LE. PLEVS_GATE(38)) THEN
           IGD  = 38
           IGU  = 38
           PUPP = 0.
           PDWN = PLEVS_GATE(38)
           GO TO 205
          ENDIF
!
!   INTERPOLATE T,Q,DQ AT 38 LEVELS TO GLA 9 SIGMA LEVELS
!   +++++++++++++++++++++++++++++++++++++++++++++++++++
            PUPP = 10.
            PDWN = 1000.0
            IGU  = 38
            IGD  = 2
!
          DO IGTLEV=1,38
           PRGAT=PLEVS_GATE(IGTLEV)   
           IF (PRGAT .LE. PRES) THEN
            PUPP = PRGAT
            IGU  = IGTLEV
            EXIT
           ENDIF
           IF (PRGAT .GT. PRES) THEN
            PDWN = PRGAT
            IGD  = IGTLEV
           ENDIF
          ENDDO
!  900     CONTINUE
!
  205   CONTINUE
!
          PUPPK = PUPP**0.286
          PDWNK = PDWN**0.286
          PRESK = PRES**0.286
!
          TEM   = PUPPK - PDWNK
          TEMU  = (PUPPK-PRESK) / TEM
          TEMD  = 1.0 - TEMU

          TEMP(L,KL) = FLD(IGU,L,1) * TEMD + FLD(IGD,L,1) * TEMU
          Q(L,KL)    = FLD(IGU,L,2) * TEMD + FLD(IGD,L,2) * TEMU
          U(L,KL)    = FLD(IGU,L,3) * TEMD + FLD(IGD,L,3) * TEMU
          V(L,KL)    = FLD(IGU,L,4) * TEMD + FLD(IGD,L,4) * TEMU
! output
          TT(L,KL) = FLD(IGU,L,1) * TEMD + FLD(IGD,L,1) * TEMU
          QQ(L,KL)    = FLD(IGU,L,2) * TEMD + FLD(IGD,L,2) * TEMU
          UU(L,KL)    = FLD(IGU,L,3) * TEMD + FLD(IGD,L,3) * TEMU
          VV(L,KL)    = FLD(IGU,L,4) * TEMD + FLD(IGD,L,4) * TEMU

! I think the following block is pointless - JTB
!          VDQ(L,KL)  = DQTRMS(IGU)  * TEMD + DQTRMS(IGD)  * TEMU
!          VDQH(L,KL) = DQTRMH(IGU)  * TEMD + DQTRMH(IGD)  * TEMU
!          VDQV(L,KL) = DQTRMV(IGU)  * TEMD + DQTRMV(IGD)  * TEMU
!          VDT(L,KL)  = DTTRMS(IGU)  * TEMD + DTTRMS(IGD)  * TEMU
!          VDT2(L,KL) = DTTRM2(IGU)  * TEMD + DTTRM2(IGD)  * TEMU
!          VDTH(L,KL) = DTTRMH(IGU)  * TEMD + DTTRMH(IGD)  * TEMU
!          VDTV(L,KL) = DTTRMV(IGU)  * TEMD + DTTRMV(IGD)  * TEMU

         omega(L,KL)=  FLD(IGU,L,6) * TEMD + FLD(IGD,L,6) * TEMU
            div(L,KL)=   fake(IGU)  * TEMD +   fake(IGD)  * TEMU

! Frank_add_5-10-04

          Q_H_adv(L,KL) = DQH(IGU)  * TEMD + DQH(IGD)  * TEMU
          Q_V_adv(L,KL) = DQV(IGU)  * TEMD + DQV(IGD)  * TEMU
          T_H_adv(L,KL) = DTH(IGU)  * TEMD + DTH(IGD)  * TEMU
          T_V_adv(L,KL) = DTV(IGU)  * TEMD + DTV(IGD)  * TEMU

! Frank_add_end
!
          Q1(L,KL)   = FLD(IGU,L,15)* TEMD + FLD(IGD,L,15)* TEMU
          Q2(L,KL)   = FLD(IGU,L,16)* TEMD + FLD(IGD,L,16)* TEMU
!
        ENDDO    ! big time loop

!------ End of vertical interplation ---!




!  ESTIMATE SURFACE EVAPORATION (BULK METHOD)
!  ------------------------------------------
       TDEGK  = FLD(1,L,1)
       QZERO  = GEOS_QSAT(TDEGK,PS(L))
       QDIFF  = QZERO-(FLD(1,L,2)*1.E-3)
       TEM1   = ((FLD(1,L,3)**2)+(FLD(1,L,4)**2))**0.5
       TEM2   = ((FLD(2,L,3)**2)+(FLD(2,L,4)**2))**0.5
       WNDMAG =  0.5 * (TEM1+TEM2)
       EVAP(L)=     1.4E-3*RHOAIR*QDIFF*WNDMAG
       SENH(L)=     1.5E-3*RHOAIR*WNDMAG*(FLD(1,L,1)-FLD(2,L,1)) * CP

! EVAP IS IN MM/SEC - TO ADD TO VDQ WE WANT TO DIVIDE BY THE MASS IN
! LOWEST LAYER -- MASS = (SURF P)*(DELTA SIGMA,IE,0.11)/(GRAV*10**-2)
! FACTOR 1000 TO CONVERT FROM KG/KG / SEC TO G/KG / SEC

       EVAP2(L) = EVAP(L)*1000.0 / ( (p_model(LM))/GRAV2 )
       SENH2(L) = SENH(L) / (CP * (p_model(LM))/GRAV2 )



!  PRECIP ESTIMATE FROM Q2 & Q1
!  ----------------------------
      TEM = (PS(L) - 1000.0) * 100.0 / 2.0                            &
             + (FLD(2,L,15) - FLD(2,L,16))*(TEM+1250.0)
      PRE1(L) = FLD(1,L,16)*TEM + FLD(2,L,16)*(TEM+1250.0)
      DQDTB(L) = FLD(1,L,12)*TEM + FLD(2,L,12)*(TEM+1250.0)
      
      !PRE2(L) = -(FLD(1,L,12)/10800.0 + DQTRMS(1)) * TEM              &
      !         -(FLD(2,L,12)/10800.0 + DQTRMS(2)) * (TEM+1250.0)

      DO 4005 K=3,38
        PRE1(L) = PRE1(L)  + FLD(K,L,16) * 2500.0
        DQDTB(L) = DQDTB(L)  + FLD(K,L,12) * 2500.0
        !! PRE2(L) = PRE2(L)  - (FLD(K,L,12)/10800.0 + DQTRMS(K)) * 2500.0
 4005 CONTINUE
      PRE1(L) = (PRE1(L)/(HEATL*GRAV) + EVAP(L)) * 86400.0
      DQDTB(L) = DQDTB(L) / (GRAV * 125.0)
      PRE3(L) = PRE1(L) - DQDTB(L)
      !! PRE2(L) = (PRE2(L)/GRAV + EVAP(L)*1000.0) * 86.40

      if (PRE1(L) .ge. 0.) then
      PCP_OBS(L) = PRE1(L)   ! 
      else
      PCP_OBS(L) = 0.
      endif


!  INTERPOLATE RADIATION TO MODEL LEVELS
!  -------------------------------------
      LR = (L-1) / 8 + 1
!
      DO 6050 KL=1,LM
!     PRES = PS(L) * SGO(KL)
!     PRES = 1000.0 * SGO(KL)
          PRES = p_model(KL)
      IF (PRES .GT. PRSR(1)) THEN
         DO 6020 KR=1,10
         IF (PRES .LE. PRSR(KR)) THEN
            TEM1 = PRSR(KR) - PRSR(KR-1)
            TEM1 = (PRSR(KR) - PRES) / TEM1
            QR(L,KL) = CSW(KR-1,LR)*TEM1 + CSW(KR,LR)*(1.0-TEM1)
            GO TO 6021
         ENDIF
 6020    CONTINUE
         QR(L,KL) = CSW(10,LR)
 6021    CONTINUE
      ELSE
         QR(L,KL) = CSW(1,LR)
      ENDIF
 6050 CONTINUE
!
      DO 6060 KL=1,38
      IF (KL .NE. 38) THEN
         PRES = 100.0 + (KL-1)*25.0
      ELSE
         PRES = PS(L)
      ENDIF
      IF (PRES .GT. PRSR(1)) THEN
         DO 6055 KR=1,10
         IF (PRES .LE. PRSR(KR)) THEN
            TEM1 = PRSR(KR) - PRSR(KR-1)
            TEM1 = (PRSR(KR) - PRES) / TEM1
            QRX(L,KL) = CSW(KR-1,LR)*TEM1 + CSW(KR,LR)*(1.0-TEM1)
            GO TO 6056
         ENDIF
 6055    CONTINUE
         QRX(L,KL) = CSW(10,LR)
 6056    CONTINUE
      ELSE
         QRX(L,KL) = CSW(1,LR)
      ENDIF
!
      QR2(L,KL) = FLD(39-KL,L,16)
      QR1(L,KL) = FLD(39-KL,L,15) - QRX(L,KL)*CP
!
 6060 CONTINUE

      ENDDO   ! NT loop


      END SUBROUTINE GATE
