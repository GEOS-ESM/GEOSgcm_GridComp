
C
C      do 100 j=1,12
C      print *
C      do 100 i=1,24
C
C      hour = i
C
C      call solar(0.,0.,j*30,i*3600,zth,slr)
C
C      print *, j, i, zth, slr
C
C100   continue
C
C      stop
C      end
C
C
C
      SUBROUTINE SOLAR(CLON,SLON,CLAT,SLAT,DAY0,SEC0,ZTH,SLR)

      INTEGER  DAYSPRYR, EQNX


      PARAMETER
     *         (DAYLEN=86400.,  DAYSPRYR=365,  AE=6.371E6
     *,         OB=23.45,       PER=102.0,     ECC=.0167
     *,         EQNX=80,        SO=1365.0)

      PARAMETER
     *          (ZERO=0.0, ONE=1.0,         TWO=2.0,     THREE=3.0
     *,          FOUR=4.0, FIVE=5.0,        SIX=6.0,     ONE80=180.0
     *,          HALF=0.5, THIRD=ONE/THREE, FOURTH=0.25
     *,          PI=3.1415926535898,        UNDEF=-999.)



C*********************************************************************
C*********************** ARIES  MODEL ********************************
C********************* SUBROUTINE SOLAR ******************************
C********************* 3 SEPTEMBER 1992 ******************************
C*********************************************************************

C     ARGUMENTS

      REAL      CLAT,         CLON
      REAL      SLAT,         SLON
      INTEGER   DAY0,         SEC0
      REAL      ZTH,          SLR

C     LOCALS

      DIMENSION TH(DAYSPRYR)

      LOGICAL   FIRST
      DATA      FIRST/ .TRUE. /

      SAVE TH

C     reichle, 10 July 2002: why is astro only called in the first
C     call to solar? shouldn't "first" be saved to achieve that?

      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL ASTRO(TH)
      ENDIF


      C = CLAT
      S = SLAT

      CH = CLON
      SH = SLON

c --------------------------------
c
c next line modified by reichle, 10 July 2002
c (prevent iday from reaching 366 in a leap year, which
c  does not work because there are only 365 days per year
c  hard-wired into this routine...)
c
c      IDAY  = DAY0
      iday = min( day0, dayspryr)
c     
c ---------------------------------      


      IDYP1 = MOD(IDAY,DAYSPRYR) + 1
      FAC   = FLOAT(SEC0)/DAYLEN

        CFAC  = COS(FAC*TWO*PI)
        SFAC  = SIN(FAC*TWO*PI)

        HC =  CFAC*CH - SFAC*SH

        THNOW = TH(IDAY)*(ONE-FAC) + TH(IDYP1)*FAC
        ZS = SIN(THNOW) * SIN(OB*PI/ONE80)
        ZC = SQRT(ONE-ZS*ZS)
        PP = (ONE-ECC*COS(THNOW - PER*PI/ONE80 )) / (1.-ECC*ECC)
        PP = PP * PP

        SZS = S*ZS
        CZC = C*ZC
        ZTH = AMAX1(SZS + CZC*HC,ZERO)
        SLR = ZTH*(SO*PP)

      RETURN
      END



      SUBROUTINE ASTRO(TH)

C   SCCS VERSION %W% %G%

C   SCCS VERSION @(#)mlsize.h	1.1 8/31/92


      INTEGER  DAYSPRYR, EQNX

      PARAMETER
     *         (DAYLEN=86400.,  DAYSPRYR=365,  AE=6.371E6
     *,         OB=23.45,       PER=102.0,     ECC=.0167
     *,         EQNX=80,        SO=1365.0)

      PARAMETER
     *          (ZERO=0.0, ONE=1.0,         TWO=2.0,     THREE=3.0
     *,          FOUR=4.0, FIVE=5.0,        SIX=6.0,     ONE80=180.0
     *,          HALF=0.5, THIRD=ONE/THREE, FOURTH=0.25
     *,          PI=3.1415926535898,        UNDEF=-999.)

      DIMENSION  TH(DAYSPRYR)

c     reichle, 10 July 2002: are these function definitions???

      PYTH(X) = SQRT(ONE-X*X)
      FUN(Y)  = (TWO*PI/(PYTH(ECC)**3))*(ONE/DAYSPRYR)
     *          *(ONE - ECC*COS(Y-PER*PI/ONE80) ) ** 2


      KM     = EQNX
      TH(KM) = ZERO
      DO 200 K=2,DAYSPRYR

      T1 = FUN(TH(KM)       )
      T2 = FUN(TH(KM)+T1/TWO)
      T3 = FUN(TH(KM)+T2/TWO)
      T4 = FUN(TH(KM)+T3    )

      KP     = MOD(EQNX-2+K,DAYSPRYR) + 1
      TH(KP) = TH(KM) + (T1 + TWO*(T2 + T3) + T4) / SIX
      KM     = KP

200   CONTINUE


      RETURN
      END

