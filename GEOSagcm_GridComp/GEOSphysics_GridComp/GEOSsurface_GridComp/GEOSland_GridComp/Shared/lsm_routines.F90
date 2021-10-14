MODULE lsm_routines

! The module contains subroutines that are shared by Catchment and Catchment-CN models.
! Sarith, 10 Nov 2015  - The first version
!                      - moved   RZDRAIN, INTERC, BASE, PARTITION, RZEQUIL, gndtp0
!                        SIBALB, catch_calc_soil_moist, catch_calc_subtile2tile
!                        gndtmp, catch_calc_tp,  catch_calc_ght, catch_calc_FT,
!                        catch_calc_wtotl, dampen_tc_oscillations and catch_echo_constants from
!                        land models
!                      - moved DZTC, FWETL, FWETC, DZGT, PHIGT, ALHMGT, FSN, CATCH_FT_WEIGHT_TP1,
!                        CATCH_FT_THRESHOLD_TEFF, CATCH_FT_THRESHOLD_ASNOW, ZERO, and ONE from
!                        from land models.
!                      - removed ITYP from argument list in subroutines INTERC, PARTITION, RZEQUIL
!                      - removed CSOIL from arguments to INTERC
!                      - removed CDCR2 from arguments to BASE
!                      - change catch_echo_constants to lsmroutines_echo_constants
! Justin, 16 Apr 2018  - replaced LAND_UPD ifdef with LAND_FIX from SurfParams, CSOIL_2 now called
!                        from SurfParams, as well as others
! Sarith, 14 Aug 2018  - Added irrigation routines, considered experimental
! Sarith, 22 Apr 2020  - moved SUBROUTINE SRUNOFF here and modified to account for separate convective and 
!                        large-scale throughfalls. FWETC and FWETL are now passed through the resource file.

  USE MAPL_BaseMod,      ONLY:                &
       NTYPS             => MAPL_NumVegTypes, &
       MAPL_Land,                             &
       MAPL_UNDEF

  USE MAPL_ConstantsMod, ONLY:           &
       PIE               => MAPL_PI,     &  ! -
       ALHE              => MAPL_ALHL,   &  ! J/kg  @15C
       ALHM              => MAPL_ALHF,   &  ! J/kg
       ALHS              => MAPL_ALHS,   &  ! J/kg
       TF                => MAPL_TICE,   &  ! K
       RGAS              => MAPL_RGAS,   &  ! J/(kg K)
       SHW               => MAPL_CAPWTR, &  ! J/kg/K  spec heat of wat
       SHI               => MAPL_CAPICE, &  ! J/kg/K  spec heat of ice
       EPSILON           => MAPL_EPSILON,&
       MAPL_H2OMW,                       &
       MAPL_AIRMW

  USE CATCH_CONSTANTS,   ONLY:                   &
       N_SNOW            => CATCH_N_SNOW,        &
       N_GT              => CATCH_N_GT,          &
       RHOFS             => CATCH_SNWALB_RHOFS,  &
       SNWALB_VISMAX     => CATCH_SNWALB_VISMAX, &
       SNWALB_NIRMAX     => CATCH_SNWALB_NIRMAX, &
       SLOPE             => CATCH_SNWALB_SLOPE,  &
       MAXSNDEPTH        => CATCH_MAXSNDEPTH,    &
       DZ1MAX            => CATCH_DZ1MAX,        &
       SHR, N_SM, SCONST, CSOIL_1,               &
       C_CANOP, SATCAPFR

  USE SURFPARAMS,        ONLY:                   &
       LAND_FIX, CSOIL_2, WEMIN, AICEV, AICEN,   &
       FLWALPHA, ASTRFR, STEXP, RSWILT

  USE SIBALB_COEFF,  ONLY: coeffsib

  USE STIEGLITZSNOW, ONLY: &
       snowrt, StieglitzSnow_calc_asnow, StieglitzSnow_calc_tpsnow, get_tf0d

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: INTERC, BASE, PARTITION, RZEQUIL, gndtp0
  PUBLIC :: SIBALB, catch_calc_soil_moist, catch_calc_subtile2tile
  PUBLIC :: gndtmp, catch_calc_tp,  catch_calc_ght, catch_calc_FT, catch_calc_wtotl
  PUBLIC :: dampen_tc_oscillations, lsmroutines_echo_constants, irrigation_rate, SRUNOFF

  ! layer depth associated with snow-free land temperatures
  !
  ! Note: DZTC = .05 is a hardwired setting of the depth of the bottom of
  ! the surface soil layer.  It should be made a parameter that is tied to
  ! the heat capacity CSOIL, which had been set to either CSOIL_1 or
  ! CSOIL_2 based on vegetation type.  For now we leave
  ! it set to 0.05 despite an apparent inconsistency with CSOIL as
  ! currently used.  We do this (again, for now) because the effects of the
  ! inconsistency are drowned out by our arbitrary assumption, in computing
  ! the thermal conductivities, that the unsaturated soil has a degree of
  ! saturation of 50%.  For the flux calculation, setting the depth to .05m
  ! here provides approximately the same fluxes as setting the depth to much
  ! closer to 0 (as the value of CSOIL_2 suggests) and assuming a degree of
  ! saturation of about 25%, which is no less realistic an assumption.  There
  ! are other impacts in wet climates regarding the effect of
  ! the depth of the water table on the thermal conductivity; these impacts
  ! are presumably very small.

  REAL,    PARAMETER, PUBLIC :: DZTC     = 0.05   ! m  layer depth for tc1, tc2, tc4

  ! ---------------------------------------------------------------------------
  !
 
   REAL,    PARAMETER :: TIMFRL = 1.0
   REAL,    PARAMETER :: TIMFRC = 0.333
  
  ! ---------------------------------------------------------------------------
  !
  ! constants for ground temperature routine (gndtp0() and gndtmp())

  REAL,    PARAMETER, DIMENSION(N_gt), PUBLIC :: DZGT = &  ! m  layer depths
       (/ 0.0988, 0.1952, 0.3859, 0.7626, 1.5071, 10.0 /)

  ! PHIGT and ALHMGT are needed for backward compatibility with
  !  off-line (land-only) MERRA replay:
  !
  ! PHIGT = porosity used in gndtp0() and gndtmp()
  !         if neg,  POROS(n) from soil moisture submodel will be used
  !
  !               |   PHIGT      ALHMGT
  ! ------------------------------------------------
  !  MERRA        |      0.45    3.34e+5
  !  Fortuna-2_3  |  -9999.       ALHM

  REAL,    PARAMETER, PUBLIC :: PHIGT   = -9999.
  REAL,    PARAMETER, PUBLIC :: ALHMGT  = ALHM

  REAL,    PARAMETER, PUBLIC :: FSN     = 1.e3*ALHMGT ! unit change J/kg/K -> J/m/K

  ! ---------------------------------------------------------------------------
  !
  ! constants for "landscape" freeze/thaw (FT) state (see subroutine catch_calc_FT())

  REAL, PARAMETER  :: CATCH_FT_WEIGHT_TP1      = 0.5   !
  REAL, PARAMETER  :: CATCH_FT_THRESHOLD_TEFF  = TF    ! [Kelvin]
  REAL, PARAMETER  :: CATCH_FT_THRESHOLD_ASNOW = 0.2   !

  REAL,    PARAMETER :: ZERO     = 0.
  REAL,    PARAMETER :: ONE      = 1.
  
  CONTAINS

!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN INTERC ]
!****

      SUBROUTINE INTERC (                                                      &
                         NCH, DTSTEP, FWETC, FWETL, TRAINL, TRAINC,SMELT,      &
                         SATCAP,BUG,                                           &
                         CAPAC,                                                &
                         THRUL, THRUC                                          &
                        )
!****
!**** THIS ROUTINE USES THE PRECIPITATION FORCING TO DETERMINE
!**** CHANGES IN INTERCEPTION AND SOIL MOISTURE STORAGE.
!****
      IMPLICIT NONE

!****
      INTEGER, INTENT(IN) ::  NCH
      REAL, INTENT(IN) :: DTSTEP, FWETC, FWETL
      REAL, INTENT(IN), DIMENSION(NCH) :: TRAINL, TRAINC, SMELT, SATCAP
      LOGICAL, INTENT(IN) :: BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CAPAC

      REAL, INTENT(OUT), DIMENSION(NCH) :: THRUC, THRUL

      INTEGER CHNO
      REAL WETINT, WATADD, CAVAIL, THRU1, THRU2, XTCORR,SMPERS

!****
!**** ------------------------------------------------------------------
!**** LOOP OVER CHIPS:
      DO 100 CHNO = 1, NCH

!**** =======================================================
!****
!**** LOAD INTERCEPTION RESERVOIR.  STEP 1: LARGE SCALE CONDENSATION.
!****
!**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
!**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.
!**** (TIME SCALE TIMFRL FOR LARGE SCALE STORMS SET TO ONE FOR FWETL=1
!**** TO REFLECT THE EFFECTIVE LOSS OF "POSITION MEMORY" WHEN STORM
!**** COVERS ENTIRE GRID SQUARE.)

      XTCORR= (1.-TIMFRL) *                                                    &
            AMIN1( 1.,(CAPAC(CHNO)/SATCAP(CHNO))/FWETL )

!****
!**** FILL INTERCEPTION RESERVOIR WITH PRECIPITATION.
!**** THRU1 IS FIRST CALCULATED AS THE AMOUNT FALLING THROUGH THE
!****    CANOPY UNDER THE ASSUMPTION THAT ALL RAIN FALLS RANDOMLY.
!****    ONLY A FRACTION 1-XTCORR FALLS RANDOMLY, THOUGH, SO THE RESULT
!****    IS MULTIPLIED BY 1-XTCORR.
!****
      WATADD = TRAINL(CHNO)*DTSTEP + SMELT(CHNO)*DTSTEP
      CAVAIL = ( SATCAP(CHNO) - CAPAC(CHNO) ) * FWETL
      WETINT = CAPAC(CHNO)/SATCAP(CHNO)
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

!**** THRU2 IS THE AMOUNT THAT FALLS IMMEDIATELY THROUGH THE CANOPY DUE
!**** TO 'POSITION MEMORY'.

      THRU2=XTCORR*WATADD

      THRUL(CHNO)=THRU1+THRU2

      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2

!****
!**** ---------------------------------------------------
!****
!**** STEP 2: MOIST CONVECTIVE PRECIPITATION.
!****
!**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
!**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.

      XTCORR= (1.-TIMFRC) *                                                    &
           AMIN1( 1.,(CAPAC(CHNO)/SATCAP(CHNO))/FWETC )

!****
!**** FILL INTERCEPTION RESERVOIR WITH PRECIPITATION.
!**** THRU1 IS FIRST CALCULATED AS THE AMOUNT FALLING THROUGH THE
!****    CANOPY UNDER THE ASSUMPTION THAT ALL RAIN FALLS RANDOMLY.
!****    ONLY A FRACTION 1-XTCORR FALLS RANDOMLY, THOUGH, SO THE RESULT
!****    IS MULTIPLIED BY 1-XTCORR.
!****
      WATADD = TRAINC(CHNO)*DTSTEP
      CAVAIL = ( SATCAP(CHNO) - CAPAC(CHNO) ) * FWETC
      WETINT = CAPAC(CHNO)/SATCAP(CHNO)
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

!**** THRU2 IS THE AMOUNT THAT FALLS IMMEDIATELY THROUGH THE CANOPY DUE
!**** TO 'POSITION MEMORY'.

      THRU2=XTCORR*WATADD

      THRUC(CHNO)=THRU1+THRU2
      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2
!****
      IF (THRUL(CHNO)+THRUC(CHNO) .LT. -1.e-8) WRITE(*,*) 'THRU= ',                        &
          THRUL(CHNO), THRUC(CHNO), TRAINC(CHNO), TRAINL(CHNO), SMELT(CHNO)
      THRUL(CHNO)=AMAX1(0., THRUL(CHNO))
      THRUC(CHNO)=AMAX1(0., THRUC(CHNO))
      
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE INTERC

!****
!**** [ END INTERC ]
!****
!****
!**** ===================================================
!**** ///////////////////////////////////////////////////
!**** ===================================================

      SUBROUTINE SRUNOFF (                                                  &
           NCH,DTSTEP,UFW4RO, FWETC, FWETL, AR1,ar2,ar4, THRUL,THRUC,       &
           frice,tp1,srfmx, BUG,                                            &
           SRFEXC,RUNSRF,                                                   &
           QINFIL                                                           &
           )

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN)    :: DTSTEP, FWETC, FWETL
      LOGICAL, INTENT (IN):: UFW4RO 
      REAL, INTENT(IN), DIMENSION(NCH) :: AR1, ar2, ar4, frice, tp1,     &
             srfmx, THRUL, THRUC
      LOGICAL, INTENT(IN) :: BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) ::  SRFEXC ,RUNSRF

      REAL, INTENT(OUT), DIMENSION(NCH) :: QINFIL

      INTEGER N
      REAL deficit,srun0,frun,qin, qinfil_l, qinfil_c, qcapac, excess_infil, srunc, srunl, ptotal

!**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO N=1,NCH

         if(.not.UFW4RO) then
            
            PTOTAL=THRUL(N) + THRUC(N)
            frun=AR1(N)
            srun0=PTOTAL*frun
            
            !**** Comment out this line in order to allow moisture
            !**** to infiltrate soil:
            !       if(tp1(n) .lt. 0.) srun0=ptotal
            
            if(ptotal-srun0 .gt. srfmx(n)-srfexc(n))                               &
                 srun0=ptotal-(srfmx(n)-srfexc(n)) 
            
            if (srun0 .gt. ptotal) srun0=ptotal
            
            RUNSRF(N)=RUNSRF(N)+srun0
            QIN=PTOTAL-srun0
            
         endif

         if(UFW4RO) then

            !**** Compute runoff from large-scale and convective storms separately:
 
            deficit=srfmx(n)-srfexc(n)
            srunl=AR1(n)*THRUL(n)
            qinfil_l=(1.-ar1(n))*THRUL(n)
            qcapac=deficit*FWETL
            
            if(qinfil_l .gt. qcapac) then
               excess_infil=qinfil_l-qcapac
               srunl=srunl+excess_infil
               qinfil_l=qinfil_l-excess_infil
            endif
         
            srunc=AR1(n)*THRUC(n)
            qinfil_c=(1.-ar1(n))*THRUC(n)
            qcapac=deficit*FWETC
            
            if(qinfil_c .gt. qcapac) then
               excess_infil=qinfil_c-qcapac
               srunc=srunc+excess_infil
               qinfil_c=qinfil_c-excess_infil
            endif
         
            !**** Comment out this line in order to allow moisture
            !**** to infiltrate soil:
            !       if(tp1(n) .lt. 0.) srun0=ptotal
         
            if (srunl .gt. THRUL(n)) srunl=THRUL(n)
         
            if (srunc .gt. THRUC(n)) srunc=THRUC(n)
            
            RUNSRF(N)=RUNSRF(N)+srunl+srunc
            QIN=THRUL(n)+THRUC(n)-(srunl+srunc)

         endif

         SRFEXC(N)=SRFEXC(N)+QIN
         RUNSRF(N)=RUNSRF(N)/DTSTEP
         QINFIL(N)=QIN/DTSTEP
     
      END DO

      RETURN

      END SUBROUTINE SRUNOFF

!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
      SUBROUTINE BASE (                                                        &
                       NCH,DTSTEP,BF1,BF2,BF3,CDCR1,FRICE,COND,GNU,            &
                       CATDEF,                                                 &
                       BFLOW                                                   &
                      )

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: BF1, BF2, BF3, CDCR1, FRICE, COND,   &
          GNU

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CATDEF

      REAL, INTENT(OUT), DIMENSION(NCH) :: BFLOW


      INTEGER N
      REAL ZBAR, ashift

      data ashift/0./


      DO N=1,NCH
         ! note intentionally opposite sign w.r.t. zbar defined above, - reichle, 16 Nov 2015
        ZBAR=SQRT(1.e-20+catdef(n)/bf1(n))-bf2(n)
        BFLOW(N)=(1.-FRICE(N))*1000.*                                          &
              cond(n)*exp(-(bf3(n)-ashift)-gnu(n)*zbar)/gnu(n)
! *1000 is to convert from m/s to mm/s
        IF (CATDEF(N) .GE. CDCR1(N)) BFLOW(N)=0.
!#ifdef LAND_UPD
	IF (LAND_FIX) THEN
	      bflow(n)=amin1(1000.*cond(n),bflow(n))
	ELSE
!#else
	      bflow(n)=amin1(cond(n),bflow(n))
	END IF
!#endif
        CATDEF(N)=CATDEF(N)+BFLOW(N)*dtstep
        ENDDO


      RETURN
      END SUBROUTINE BASE


!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****

      SUBROUTINE PARTITION (                                                   &
                            NCH,DTSTEP,DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,     &
                            PSIS,BEE,poros,WPWET,                              &
                            ars1,ars2,ars3,ara1,ara2,ara3,ara4,                &
                            arw1,arw2,arw3,arw4,BUG,                           &
                            srfexc,catdef,runsrf,                              &
                            AR1, AR2, AR4, srfmx, srfmn,                       &
                            SWSRF1,SWSRF2,SWSRF4,RZI                           &
                           )

      IMPLICIT NONE

! -------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NCH

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,  &
                                          PSIS,BEE,poros,WPWET,                &
                                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,  &
                                          arw1,arw2,arw3,arw4

      LOGICAL, INTENT(IN) :: BUG
! -------------------------------------------------------------------
      REAL, INTENT(INOUT), DIMENSION(NCH) :: srfexc,catdef,runsrf
! -------------------------------------------------------------------
      REAL, INTENT(OUT), DIMENSION(NCH) :: AR1, AR2, AR4, srfmx, srfmn,        &
                                           SWSRF1, SWSRF2, SWSRF4, RZI
! -------------------------------------------------------------------
      INTEGER :: N

      REAL :: cor, A150, W150, WMIN, AX, WMNEW, WRZ, TERM1, TERM2, TERM3,      &
              AREA0, AREA1, AREA2, AREA3, AREA4, ASCALE, WILT, D1, D2, CDI,    &
              DELTA1, DELTA2, DELTA4, MULTAR, CATDEFX, RZEQX, RZEQW, FACTOR,   &
              X0, RZEQY, CATDEFW, AR1W, ASUM, RZEQYI, RZEQWI, RZEQXI, AR20,    &
              ARG1, EXPARG1, ARG2, EXPARG2, ARG3, EXPARG3  !, surflay

      LOGICAL :: LSTRESS


      DATA LSTRESS/.FALSE./    !,surflay/20./

!****
!**** --------------------------------------------------

!rr   next line for debugging, sep 23, 2003, reichle
!rr
!rr   write (*,*) 'entering partition()'

      DO N=1,NCH

        WILT=WPWET(N)
        WRZ=RZEXC(N)/VGWMAX(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

        AR1(N)= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFX)                         &
                 /(1.+ars2(n)*CATDEFX+ars3(n)*CATDEFX*CATDEFX)))

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*                           &
                 (1.+arw1(n)*CATDEFX)                                          &
                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

!**** CRITICAL VALUE 1: AVERAGE MOISTURE IN ROOT ZONE AT WMIN
!**** ASSOCIATED WITH CATDEF.
        ARG1=AMAX1(-40., AMIN1(40., -AX*(1.-WMIN)))
        EXPARG1=EXP(ARG1)
        RZEQX=(WMIN-1.-(2./AX))*EXPARG1 + WMIN + (2./AX)
        RZEQXI=AX*EXPARG1 *                                                    &
          ( -1. -(2./AX) - (2./(AX*AX)) + WMIN + (WMIN/AX) )                   &
          + WMIN + 2./AX
        AR20=1.+(-AX-1.+AX*WMIN)*EXPARG1
        RZEQXI=RZEQXI/(AR20+1.E-20)

!**** CRITICAL VALUE 2: AVERAGE MOISTURE IN ROOT ZONE WHEN WMIN
!**** IS EXACTLY AT WILTING POINT.
        ARG2=AMAX1(-40., AMIN1(40., -AX*(1.-WILT)))
        EXPARG2=EXP(ARG2)
        RZEQW=(WILT-1.-(2./AX))*EXPARG2 + WILT + (2./AX)
        RZEQWI=AX*EXPARG2 *                                                    &
         ( -1. -(2./AX) - (2./(AX*AX)) + WILT + (WILT/AX) )                    &
         + WILT + 2./AX
        AR20=1.+(-AX-1.+AX*WILT)*EXPARG2
        RZEQWI=RZEQWI/(AR20+1.E-20)

!**** SITUATION 1: CATDEF LE CDCR1
        IF(CATDEF(N) .LE. CDCR1(N)) THEN
          RZEQY=RZEQX+WRZ
          RZEQYI=RZEQXI+WRZ
          WMNEW=WMIN+WRZ
          ARG3=AMAX1(-40., AMIN1(40., -AX*(1.-WMNEW)))
          EXPARG3=EXP(ARG3)
          AREA1=(1.+AX-AX*WMIN)*EXPARG1
          AREA2=(1.+AX-AX*WMNEW)*EXPARG3
          IF(WMNEW .GE. WILT) THEN
            AR1(N)=AR1(N)+AREA2-AREA1
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF
          IF(WMNEW .LT. WILT) THEN
            AREA3=(1.+AX-AX*WILT)*EXPARG2
            AR1(N)=AR1(N)+AREA3-AREA1
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQX+WRZ-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF
          ENDIF

!**** SITUATION 2: CATDEF GT CDCR1
        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQY=WILT+(RZEQX-WILT)*FACTOR+WRZ
          RZEQYI=WILT+(RZEQXI-WILT)*FACTOR+WRZ

          IF(RZEQY .LT. WILT) THEN
            IF(RZEQY .LT. WILT-.001) THEN
!rr                WRITE(*,*) 'RZEXC WAY TOO LOW!  N=',N,' RZEQY=',RZEQY
!rr                WRITE(*,*) 'SRFEXC=',SRFEXC(N),'RZEXC=',RZEXC(N),
!rr     &                     'CATDEF=',CATDEF(N)
!             ELSE
!               WRITE(*,*) 'RZEXC TOO LOW  N=',N
              ENDIF
            RZEQY=WILT
            RZEQYI=WILT
            ENDIF

          IF(RZEQY .GE. RZEQX) THEN  ! RZEXC BRINGS MOISTURE ABOVE CDCR1 POINT
            WMNEW=WMIN+(RZEQY-RZEQX)
            ARG3=AMAX1(-40., AMIN1(40., -AX*(1.-WMNEW)))
            EXPARG3=EXP(ARG3)
            AREA1=(1.+AX-AX*WMIN)*EXPARG1
            AREA2=(1.+AX-AX*WMNEW)*EXPARG3
            AR1(N)=AR1(N)+(AREA2-AREA1)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQX .AND. RZEQY .GE. RZEQW) THEN
            CATDEFW=CDCR2(N)+((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            FACTOR=(RZEQY-RZEQW)/(RZEQX-RZEQW)
            AR1(N)=AR1W+FACTOR*(AR1(N)-AR1W)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQW) THEN
            CATDEFW=CDCR2(N)+((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            AR1(N)=AR1W
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQY-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF

          ENDIF

        RZI(N)=RZEQYI

        SWSRF1(N)=1.
!mjs: changed .001 temporarily because of large bee.
        SWSRF2(N)=AMIN1(1., AMAX1(0.01, RZEQYI))
        SWSRF4(N)=AMIN1(1., AMAX1(0.01, WILT))

!**** EXTRAPOLATION OF THE SURFACE WETNESSES

! 1st step: surface wetness in the unstressed fraction without considering
!           the surface excess; we just assume an equilibrium profile from
!           the middle of the root zone to the surface.

        SWSRF2(N)=((SWSRF2(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))
        SWSRF4(N)=((SWSRF4(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))

! srfmx is the maximum amount of water that can be added to the surface layer
! The choice of defining SWSRF4 like SWSRF2 needs to be better examined.
        srfmx(n)=ar2(n)*(1.-swsrf2(n))*(dzsf(n)*poros(n))
        srfmx(n)=srfmx(n)+ar4(n)*(1.-swsrf4(n))*(dzsf(n)*poros(n))
!**** For calculation of srfmn, assume surface moisture associated with
!**** AR1 is constantly replenished by water table.
        srfmn(n)=-(ar2(n)*swsrf2(n)+ar4(n)*swsrf4(n))*(dzsf(n)*poros(n))

        if(srfexc(n).gt.srfmx(n)) then
            cor=srfexc(n)-srfmx(n)     !  The correction is here
            srfexc(n)=srfmx(n)
            catdef(n)=catdef(n)-cor
            if(catdef(n).lt.0.) then
              runsrf(n)=runsrf(n)-catdef(n)/dtstep
              catdef(n)=0.
              endif
          else if(srfexc(n).lt.srfmn(n)) then
            cor=srfexc(n)-srfmn(n)
            catdef(n)=catdef(n)-cor
            srfexc(n)=srfmn(n)
          else
            cor=0.
          endif

        SWSRF2(N)=SWSRF2(N)+SRFEXC(N)/(dzsf(n)*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF2(N)=AMIN1(1., AMAX1(1.E-5, SWSRF2(N)))
        swsrf4(n)=swsrf4(n)+srfexc(n)/(dzsf(n)*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF4(N)=AMIN1(1., AMAX1(1.E-5, SWSRF4(N)))

        IF (AR1(N) .ge. 1.-1.E-5) then
          AR1(N)=1.
          AR2(N)=0.
          AR4(N)=0.
          SWSRF2(N)=1.
          SWSRF4(N)=wilt
          ENDIF

        IF (AR1(N) .LT. 0.) then
!rr          IF(AR1(N) .LT. -1.E-3) WRITE(*,*) 'AR1 TOO LOW: AR1=',AR1(N)
          AR1(N)=0.
          ENDIF
        ar1(n)=amax1(0., amin1(1., ar1(n)))
        ar2(n)=amax1(0., amin1(1., ar2(n)))
        ar4(n)=amax1(0., amin1(1., ar4(n)))
        asum=ar1(n)+ar2(n)+ar4(n)
        if(asum .lt. .9999 .or. asum .gt. 1.0001) then
          write(*,*) 'Areas do not add to 1: sum=',asum,'N=',n
       endif


        ENDDO


      RETURN
      END SUBROUTINE PARTITION

!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------

      SUBROUTINE RZEQUIL (                                                     &
                          NCH,CATDEF,VGWMAX,CDCR1,CDCR2,WPWET,                 &
                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,                  &
                          arw1,arw2,arw3,arw4,                                 &
                          RZEQ                                                 &
                         )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: CATDEF, VGWMAX, CDCR1, CDCR2,        &
                   WPWET, ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1,      &
                   arw2, arw3, arw4

      REAL, INTENT(OUT), DIMENSION(NCH) :: RZEQ

      INTEGER N
      REAL AX,WMIN,ASCALE,cdi,wilt,catdefx,factor,ARG1,EXPARG1

! ----------------------------------------------------------------------

      DO N=1,NCH

        WILT=WPWET(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*(1.+arw1(n)*CATDEFX)       &
                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

        ARG1=AMAX1(-40., AMIN1(40., -AX*(1.-WMIN)))
        EXPARG1=EXP(ARG1)
        RZEQ(N)=(WMIN-1.-(2./AX))*EXPARG1 + WMIN + (2./AX)

        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQ(N)=WILT+(RZEQ(N)-WILT)*FACTOR
          ENDIF

! scaling:
        RZEQ(N)=AMIN1(1.,AMAX1(0.,RZEQ(N)))
        RZEQ(N)=RZEQ(N)*VGWMAX(N)

      ENDDO

      RETURN
      END SUBROUTINE RZEQUIL

!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------


      subroutine gndtp0(t1,phi,zbar,thetaf,ht,fh21w,fh21i,fh21d,                   &
                      dfh21w,dfh21i,dfh21d,tp)

! using a diffusion equation this code generates ground temperatures
! with depth given t1
!            *****************************************
!        input
!        dts     timestep in seconds
!        t1      terrestrial (layer 1) temperature in deg C
!        phi     porosity
!        zbar    mean depth to the water table.
!        thetaf  mean vadose zone soil moisture factor (0-1)
!        output,
!        ht      heat content in layers 2-7
!        tp      ground temperatures in layers 2-7
!        tdeep   the temperature of the "deep"
!        f21     heat flux between layer 2 and the terrestrial layer (1)
!        df21    derivative of f21 with respect to temperature
!        xfice   a total soil column ice factor (0-1)
!             ***********************************

      REAL, INTENT(IN) :: phi, ZBAR, THETAF
      REAL, INTENT(IN), DIMENSION(*) :: HT
      REAL, INTENT(IN), DIMENSION(N_SM) :: T1

      REAL, INTENT(OUT) :: FH21W, FH21I, FH21D, DFH21W, DFH21I, DFH21D
      REAL, INTENT(OUT), DIMENSION(*) :: TP

      INTEGER L, K
      REAL, DIMENSION(N_GT) :: FICE, SHC, ZC, XKLH
      REAL, DIMENSION(N_GT+1) :: FH, ZB
      REAL SHW0, SHI0, SHR0, WS, XW, A1, TK1, A2, TK2, TK3, TKSAT,         &
           XWI, XD1, XD2, DENOM, XKLHW, TKDRY

      !data dz/0.0988,0.1952,0.3859,0.7626,1.5071,10.0/
      !DATA PHI/0.45/, FSN/3.34e+8/, SHR/2.4E6/

!     initialize parameters
      shw0=SHW*1000. ! PER M RATHER THAN PER KG
      shi0=SHI*1000. ! PER M RATHER THAN PER KG
      shr0=SHR*1000. ! PER M RATHER THAN PER KG [kg of water equivalent density]

! calculate the boundaries, based on the layer thicknesses(DZGT)

      zb(1)=-DZTC
      zb(2)=zb(1)-DZGT(1)
      shc(1)=shr0*(1.-phi)*DZGT(1)
      zc(1)=0.5*(zb(1)+zb(2))

! evaluates the temperatures in the soil layers based on the heat values.
!             ***********************************
!             input:
!             xw - water in soil layers, m
!             ht - heat in soil layers
!             fsn - heat of fusion of water   J/m
!             shc - specific heat capacity of soil
!             shi - specific heat capacity of ice
!             shw - specific heat capcity of water
!             snowd - snow depth, equivalent water m
!             output:
!             tp - temperature of layers, c
!             fice - fraction of ice of layers
!             pre - extra precipitation, i.e. snowmelt, m s-1
!             snowd - snow depth after melting, equivalent water m.
!             ***********************************
! determine fraction of ice and temp of soil layers based on layer
! heat and water content

      ws=phi*DZGT(1)  ! PORE SPACE IN LAYER 2
      xw=0.5*ws     ! ASSUME FOR THESE CALCULATIONS THAT THE PORE SPACE
                    ! IS ALWAYS HALF FILLED WITH WATER.  XW IS THE
                    ! AMOUNT OF WATER IN THE LAYER.

      tp(1)=0.
      FICE(1) = AMAX1( 0., AMIN1( 1., -ht(1)/(fsn*xw) ) )

      IF(FICE(1) .EQ. 1.) THEN
          tp(1)=(ht(1)+xw*fsn)/(shc(1)+xw*shi0)
        ELSEIF(FICE(1) .EQ. 0.) THEN
          tp(1)=ht(1)/(shc(1)+xw*shw0)
        ELSE
          TP(1)=0.
        ENDIF

! evaluates:  layer thermal conductivities
! *****************************************
!             from farouki(cold regions sci and tech, 5, 1981,
!             67-75) the values for tk1,tk2,tk3 are as follows:
!             tk2=2.2**(phi(l,ibv)-xw), tk3=.57**xw, and
!             tk1=3**(1-phi(l,ibv) for %sand<.5, and
!             for computation purposes i have fit these eqs.
!             to 2nd order polynomials.
!             ***********************************
!             input:
!             sklh - soil heat conductivities of layers
!             zb - soil layer boundaries, m
!             zc - soil layer centers, m
!             DZGT - layer thickness, m
!             w - soil water content, m
!             phi - soil porosity, dimensionless
!             q - % sand, silt, clay, peat
!             fice - fraction of ice in layers
!             output:
!             xklh - thermal conductivity, w m-2 k-1
!             ***********************************
! lets get the thermal conductivity for the layers

      a1=1-phi
      tk1=1.01692+a1*(0.89865+1.06211*a1)
      xw=phi*(1.-fice(1))
      a2=phi-xw
      tk2=1.00543+a2*(0.723371+.464342*a2)
      tk3=0.998899+xw*(-0.548043+0.120291*xw)
      tksat=tk1*tk2*tk3

      xwi=1.0
      if (zbar .le. zb(2))then
            xwi=thetaf
         elseif (zbar .ge. zb(2) .and. zbar .le. zb(1))then
            xd1=zb(1)-zbar
            xd2=zbar-zb(2)
            xwi=((xd1*thetaf)+xd2)/(xd1+xd2)
         endif

      xwi=min(xwi,1.)
      tkdry=0.226 ! = .039*0.45^(-2.2), from Farouki, p. 71
      xklh(1)=(tksat-tkdry)*xwi + tkdry
      xklhw=tksat

      denom=-(DZTC*0.5)-zc(1)
      fh21w=-xklhw  *(t1(1)-TF-tp(1))/denom
      fh21i=-xklh(1)*(t1(2)-TF-tp(1))/denom
      fh21d=-xklh(1)*(t1(3)-TF-tp(1))/denom
      dfh21w=-xklhw/denom
      dfh21i=-xklh(1)/denom
      dfh21d=dfh21i


      return
      end subroutine gndtp0

!
! -----------------------------------------------------------------------
!
      SUBROUTINE SIBALB (NCH, ITYP, VLAI, VGRN, ZTH,                 &
                 SCALVDR, SCALVDF, SCALIDR, SCALIDF,                 &
                 AVISDR, ANIRDR, AVISDF, ANIRDF, MODIS_SCALE         &
                 )

      USE sibalb_coeff

      IMPLICIT NONE

! OUTPUTS:
! AVISDR:   visible, direct albedo.
! ANIRDR:   near infra-red, direct albedo.
! AVISDF:   visible, diffuse albedo.
! ANIRDF:   near infra-red, diffuse albedo.

! INPUTS:
! SCALVDR:  MODIS scale factor for visible, direct.
! SCALVDF:  MODIS scale factor for visible, diffuse.
! SCALIDR:  MODIS scale factor for NIR, direct.
! SCALIDF:  MODIS scale factor for NIR, diffuse.
! VLAI:     the leaf area index.
! VGRN:     the greenness index.
      INTEGER :: M,J,K

! *********************************************************************

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) ::  ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: SCALVDR,SCALVDF,SCALIDR,SCALIDF,     &
                                             VLAI, VGRN,  ZTH

      REAL, INTENT(OUT), DIMENSION(NCH) :: AVISDR, ANIRDR, AVISDF,             &
                        ANIRDF
      LOGICAL, INTENT(IN), OPTIONAL :: MODIS_SCALE


      REAL, PARAMETER :: ALVDRS = 0.100
      REAL, PARAMETER :: ALIDRS = 0.200
      REAL, PARAMETER :: ALVDRD = 0.300
      REAL, PARAMETER :: ALIDRD = 0.350
      REAL, PARAMETER :: ALVDRI = 0.700
      REAL, PARAMETER :: ALIDRI = 0.700


!      REAL, PARAMETER :: WEMIN  = 13.0   ! [KG/M2]

! ALVDRS:  Albedo of soil for visible   direct  solar radiation.
! ALIDRS:  Albedo of soil for infra-red direct  solar radiation.
! ALVDFS:  Albedo of soil for visible   diffuse solar radiation.
! ALIDFS:  Albedo of soil for infra-red diffuse solar radiation.

      INTEGER, PARAMETER :: NLAI = 14

      REAL, PARAMETER :: EPSLN = 1.E-6
      REAL, PARAMETER :: BLAI = 0.5
      REAL, PARAMETER :: DLAI = 0.5

      REAL, PARAMETER :: ALATRM = BLAI + (NLAI - 1) * DLAI - EPSLN

      INTEGER, PARAMETER :: NTYPS_SIB=9



! ITYP: Vegetation type as follows:
!                  1:  BROADLEAF EVERGREEN TREES
!                  2:  BROADLEAF DECIDUOUS TREES
!                  3:  NEEDLELEAF TREES
!                  4:  GROUND COVER
!                  5:  BROADLEAF SHRUBS
!                  6:  DWARF TREES (TUNDRA)
!                  7:  BARE SOIL
!                  8:  DESERT
!                  9:  ICE
!  NCH: Chip index
!

	INTEGER I, LAI
	REAL FAC, GAMMA, BETA, ALPHA, DX, DY, ALA, FVEG
        REAL, DIMENSION(2) :: GRN
        REAL, DIMENSION(4,NTYPS_SIB) :: SNWALB (4, NTYPS_SIB)
        REAL, DIMENSION(NTYPS_SIB) :: SNWMSK
! by Teppei J. Yasunari
        REAL ,DIMENSION(N_snow,NCH) :: DENEL

      DATA GRN /0.33, 0.67/
! moved to catch_constants - SM 10/02/2012
!      REAL, PARAMETER :: SNWALB_VISMAX = 0.7
!      REAL, PARAMETER :: SNWALB_VISMIN = 0.5
!      REAL, PARAMETER :: SNWALB_NIRMAX = 0.5
!      REAL, PARAMETER :: SNWALB_NIRMIN = 0.3
      REAL, DIMENSION(NTYPS_SIB) :: SNWMID


! [ Definition of Functions: ]
!
!	REAL COEFFSIB

! --------------------------------------------------

      LOGICAL :: MODIS_SCALE_


!   Constants used in albedo calculations:

      REAL ALVDR (NLAI, 2, NTYPS_SIB)
      REAL BTVDR (NLAI, 2, NTYPS_SIB)
      REAL GMVDR (NLAI, 2, NTYPS_SIB)
      REAL ALIDR (NLAI, 2, NTYPS_SIB)
      REAL BTIDR (NLAI, 2, NTYPS_SIB)
      REAL GMIDR (NLAI, 2, NTYPS_SIB)

!  (Data statements for ALVDR described in full; data statements for
!   other constants follow same framework.)


!    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 1), I = 1, 14)                                      &
      	  /0.0808, 0.0796, 0.0792, 0.0790, 10*0.0789/

!    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 1), I = 1, 14)                                      &
      	  /0.0788, 0.0775, 0.0771, 0.0769, 10*0.0768/

!    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 2), I = 1, 14)                                      &
      	  /0.0803, 0.0790, 0.0785, 0.0784, 3*0.0783, 7*0.0782/

!    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 2), I = 1, 14)                                      &
      	  /0.0782, 0.0770, 0.0765, 0.0763, 10*0.0762/

!    NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
	DATA (ALVDR (I, 1, 3), I = 1, 14)                                      &
      	  /0.0758, 0.0746, 0.0742, 0.0740, 10*0.0739/

!    NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 3), I = 1, 14)                                      &
      	  /0.0683, 0.0672, 0.0667, 2*0.0665, 9*0.0664/

!    GROUNDCOVER (ITYP=4); GREEN=0.33; LAI=.5-7
	DATA (ALVDR (I, 1, 4), I = 1, 14)                                      &
      	  /0.2436, 0.2470, 0.2486, 0.2494, 0.2498, 0.2500, 2*0.2501,           &
      		6*0.2502 /

!    GROUNDCOVER (ITYP=4); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 4), I = 1, 14) /14*0.1637/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 5), I = 1, 14)                                      &
          /0.0807, 0.0798, 0.0794, 0.0792, 0.0792, 9*0.0791/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 5), I = 1, 14)                                      &
          /0.0787, 0.0777, 0.0772, 0.0771, 10*0.0770/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 6), I = 1, 14)                                      &
          /0.0802, 0.0791, 0.0787, 0.0786, 10*0.0785/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 6), I = 1, 14)                                      &
          /0.0781, 0.0771, 0.0767, 0.0765, 0.0765, 9*0.0764/


!    BARE SOIL
	DATA (ALVDR (I, 1, 7), I = 1, 14) /14*ALVDRS/
	DATA (ALVDR (I, 2, 7), I = 1, 14) /14*ALVDRS/

!    DESERT
	DATA (ALVDR (I, 1, 8), I = 1, 14) /14*ALVDRD/
	DATA (ALVDR (I, 2, 8), I = 1, 14) /14*ALVDRD/

!    ICE
	DATA (ALVDR (I, 1, 9), I = 1, 14) /14*ALVDRI/
	DATA (ALVDR (I, 2, 9), I = 1, 14) /14*ALVDRI/
!****
!**** -------------------------------------------------
	DATA (BTVDR (I, 1, 1), I = 1, 14)                                      &
        /0.0153, 0.0372, 0.0506, 0.0587, 0.0630, 0.0652, 0.0663,               &
      	0.0668, 0.0671, 0.0672, 4*0.0673 /
	DATA (BTVDR (I, 2, 1), I = 1, 14)                                      &
     	  /0.0135, 0.0354, 0.0487, 0.0568, 0.0611, 0.0633, 0.0644,             &
     	0.0650, 0.0652, 0.0654, 0.0654, 3*0.0655 /
	DATA (BTVDR (I, 1, 2), I = 1, 14)                                      &
      	  /0.0148, 0.0357, 0.0462, 0.0524, 0.0554, 0.0569, 0.0576,             &
      	0.0579, 0.0580, 0.0581, 0.0581, 3*0.0582 /
	DATA (BTVDR (I, 2, 2), I = 1, 14)                                      &
      	  /0.0131, 0.0342, 0.0446, 0.0508, 0.0539, 0.0554, 0.0560,             &
      	0.0564, 0.0565, 5*0.0566 /
	DATA (BTVDR (I, 1, 3), I = 1, 14)                                      &
      	  /0.0108, 0.0334, 0.0478, 0.0571, 0.0624, 0.0652, 0.0666,             &
      	0.0673, 0.0677, 0.0679, 4*0.0680 /
	DATA (BTVDR (I, 2, 3), I = 1, 14)                                      &
      	  /0.0034, 0.0272, 0.0408, 0.0501, 0.0554, 0.0582, 0.0597,             &
      		0.0604, 0.0608, 0.0610, 4*0.0611 /
	DATA (BTVDR (I, 1, 4), I = 1, 14)                                      &
      	  /0.2050, 0.2524, 0.2799, 0.2947, 0.3022, 0.3059, 0.3076,             &
      		0.3085, 0.3088, 0.3090, 4*0.3091 /
	DATA (BTVDR (I, 2, 4), I = 1, 14)                                      &
      	  /0.1084, 0.1404, 0.1617, 0.1754, 0.1837, 0.1887, 0.1915,             &
      		0.1931, 0.1940, 0.1946, 0.1948, 0.1950, 2*0.1951  /
        DATA (BTVDR (I, 1, 5), I = 1, 14)                                      &
          /0.0203, 0.0406, 0.0548, 0.0632, 0.0679, 0.0703, 0.0716,             &
           0.0722, 0.0726, 0.0727, 0.0728, 0.0728, 0.0728, 0.0729 /
        DATA (BTVDR (I, 2, 5), I = 1, 14)                                      &
          /0.0184, 0.0385, 0.0526, 0.0611,  0.0658, 0.0683, 0.0696,            &
           0.0702, 0.0705, 0.0707, 4*0.0708 /
        DATA (BTVDR (I, 1, 6), I = 1, 14)                                      &
          /0.0199, 0.0388, 0.0494,  0.0554, 0.0584, 0.0599, 0.0606,            &
           0.0609, 0.0611, 5*0.0612  /
        DATA (BTVDR (I, 2, 6), I = 1, 14)                                      &
          /0.0181, 0.0371, 0.0476, 0.0537,  0.0568, 0.0583, 0.0590,            &
           0.0593, 0.0595, 0.0595, 4*0.0596 /
	DATA (BTVDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 9), I = 1, 14) /14*0./

!****
!**** -----------------------------------------------------------
	DATA (GMVDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.0814, 0.1361, 0.2078, 0.2650, 0.2986, 0.3169,  0.3265,            &
        	   0.3313, 0.3337, 0.3348, 0.3354, 0.3357, 2*0.3358 /
	DATA (GMVDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.0760, 0.1336, 0.2034, 0.2622, 0.2969, 0.3159,  0.3259,            &
       	   0.3309, 0.3333, 0.3346, 0.3352, 0.3354, 2*0.3356 /
	DATA (GMVDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.0834, 0.1252, 0.1558, 0.1927, 0.2131,   0.2237, 0.2290,           &
       	   0.2315, 0.2327, 0.2332, 0.2335, 2*0.2336, 0.2337 /
	DATA (GMVDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.0789, 0.1235, 0.1531, 0.1912, 0.2122, 0.2232,  0.2286,            &
      	   0.2312, 0.2324, 0.2330, 0.2333, 0.2334, 2*0.2335 /
	DATA (GMVDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.0647, 0.1342, 0.2215, 0.2968, 0.3432, 0.3696, 0.3838,             &
       	   0.3912, 0.3950, 0.3968, 0.3978, 0.3982, 0.3984, 0.3985 /
	DATA (GMVDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.0258, 0.1227, 0.1999, 0.2825, 0.3339, 0.3634, 0.3794,             &
       	   0.3877, 0.3919, 0.3940, 0.3950, 0.3956, 0.3958, 0.3959 /
	DATA (GMVDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.3371, 0.5762, 0.7159, 0.7927, 0.8324, 0.8526,  0.8624,            &
       	   0.8671, 0.8693, 0.8704, 0.8709, 0.8710, 2*0.8712 /
	DATA (GMVDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.2634, 0.4375, 0.5532, 0.6291, 0.6763, 0.7048, 0.7213,             &
       	   0.7310, 0.7363, 0.7395, 0.7411, 0.7420, 0.7426, 0.7428 /
        DATA (GMVDR (I, 1, 5), I = 1, 14)                                      &
           /0.0971, 0.1544, 0.2511, 0.3157, 0.3548, 0.3768, 0.3886,            &
            0.3948, 0.3978, 0.3994, 0.4001, 0.4006, 0.4007, 0.4008 /
        DATA (GMVDR (I, 2, 5), I = 1, 14)                                      &
           /0.0924, 0.1470, 0.2458, 0.3123, 0.3527, 0.3756, 0.3877,            &
            0.3942, 0.3974, 0.3990, 0.3998, 0.4002, 0.4004, 0.4005 /
        DATA (GMVDR (I, 1, 6), I = 1, 14)                                      &
           /0.0970, 0.1355, 0.1841, 0.2230, 0.2447,  0.2561, 0.2617,           &
            0.2645, 0.2658, 0.2664, 0.2667, 3*0.2669 /
        DATA (GMVDR (I, 2, 6), I = 1, 14)                                      &
           /0.0934, 0.1337, 0.1812, 0.2213, 0.2437, 0.2554, 0.2613,            &
            0.2642, 0.2656, 0.2662, 0.2665, 0.2667, 0.2667, 0.2668 /
	DATA (GMVDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 9), I = 1, 14) /14*1./

!****
!****  -----------------------------------------------------------

	DATA (ALIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.2867,  0.2840, 0.2828, 0.2822, 0.2819, 0.2818, 2*0.2817,          &
       	   6*0.2816 /
	DATA (ALIDR (I, 2, 1), I = 1, 14)                                      &
        	  /0.3564, 0.3573, 0.3577, 0.3580, 2*0.3581, 8*0.3582 /
	DATA (ALIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.2848, 0.2819, 0.2804, 0.2798, 0.2795, 2*0.2793, 7*0.2792 /
	DATA (ALIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.3544, 0.3550, 0.3553, 2*0.3555, 9*0.3556 /
	DATA (ALIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.2350, 0.2311, 0.2293, 0.2285, 0.2281, 0.2280, 8*0.2279 /
	DATA (ALIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.2474, 0.2436, 0.2418, 0.2410, 0.2406, 0.2405, 3*0.2404,           &
       	   5*0.2403 /
	DATA (ALIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.5816, 0.6157, 0.6391, 0.6556, 0.6673, 0.6758, 0.6820,             &
       	   0.6866, 0.6899, 0.6924, 0.6943, 0.6956, 0.6966, 0.6974 /
	DATA (ALIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.5489, 0.5770, 0.5955, 0.6079, 0.6163, 0.6221, 0.6261,             &
       	   0.6288, 0.6308, 0.6321, 0.6330, 0.6337, 0.6341, 0.6344 /
        DATA (ALIDR (I, 1, 5), I = 1, 14)                                      &
           /0.2845, 0.2837, 0.2832, 0.2831, 0.2830, 9*0.2829 /
        DATA (ALIDR (I, 2, 5), I = 1, 14)                                      &
           /0.3532, 0.3562, 0.3578,  0.3586, 0.3590, 0.3592, 0.3594,           &
            0.3594, 0.3594, 5*0.3595 /
        DATA (ALIDR (I, 1, 6), I = 1, 14)                                      &
           /0.2825, 0.2812, 0.2806, 0.2803, 0.2802, 9*0.2801 /
        DATA (ALIDR (I, 2, 6), I = 1, 14)                                      &
           /0.3512, 0.3538,  0.3552, 0.3559, 0.3562, 0.3564, 0.3565,           &
            0.3565, 6*0.3566 /
	DATA (ALIDR (I, 1, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 2, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 1, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 2, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 1, 9), I = 1, 14) /14*ALIDRI/
	DATA (ALIDR (I, 2, 9), I = 1, 14) /14*ALIDRI/

!****
!**** -----------------------------------------------------------
	DATA (BTIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.1291, 0.1707, 0.1969, 0.2125, 0.2216,   0.2267, 0.2295,           &
       	   0.2311, 0.2319, 0.2323, 0.2326, 2*0.2327, 0.2328 /
	DATA (BTIDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.1939, 0.2357, 0.2598, 0.2735, 0.2810,  0.2851, 0.2874,            &
       	   0.2885, 0.2892, 0.2895, 0.2897, 3*0.2898 /
	DATA (BTIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.1217, 0.1522, 0.1713, 0.1820,   0.1879,  0.1910, 0.1926,          &
      	   0.1935, 0.1939, 0.1942, 2*0.1943, 2*0.1944 /
	DATA (BTIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.1781, 0.2067, 0.2221, 0.2301,   0.2342,  0.2363, 0.2374,          &
       	   0.2379, 0.2382, 0.2383, 2*0.2384, 2*0.2385 /
	DATA (BTIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.0846, 0.1299, 0.1614, 0.1814, 0.1935,   0.2004, 0.2043,           &
           0.2064, 0.2076, 0.2082, 0.2085, 2*0.2087, 0.2088 /
	DATA (BTIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.0950, 0.1410, 0.1722, 0.1921, 0.2042, 0.2111,  0.2151,            &
       	   0.2172, 0.2184, 0.2191, 0.2194, 0.2196, 2*0.2197 /
	DATA (BTIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.5256, 0.7444, 0.9908, 1.2700, 1.5680, 1.8505, 2.0767,             &
       	   2.2211, 2.2808, 2.2774, 2.2362, 2.1779, 2.1160, 2.0564 /
	DATA (BTIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.4843, 0.6714, 0.8577, 1.0335, 1.1812, 1.2858, 1.3458,             &
       	   1.3688, 1.3685, 1.3546, 1.3360, 1.3168, 1.2989, 1.2838 /
	DATA (BTIDR (I, 1, 5), I = 1, 14)                                      &
           /0.1498, 0.1930, 0.2201, 0.2364, 0.2460, 0.2514, 0.2544,            &
            0.2560, 0.2569, 0.2574, 0.2577, 0.2578, 0.2579, 0.2579 /
        DATA (BTIDR (I, 2, 5), I = 1, 14)                                      &
           /0.2184, 0.2656, 0.2927, 0.3078, 0.3159,  0.3202, 0.3224,           &
            0.3235, 0.3241, 0.3244, 0.3245, 3*0.3246 /
        DATA (BTIDR (I, 1, 6), I = 1, 14)                                      &
           /0.1369, 0.1681, 0.1860, 0.1958, 0.2010,  0.2038, 0.2053,           &
            0.2060, 0.2064, 0.2066, 0.2067, 3*0.2068 /
        DATA (BTIDR (I, 2, 6), I = 1, 14)                                      &
           /0.1969, 0.2268, 0.2416,  0.2488, 0.2521, 0.2537, 0.2544,           &
            0.2547, 0.2548, 5*0.2549 /
	DATA (BTIDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 9), I = 1, 14) /14*0./

!****
!**** --------------------------------------------------------------
	DATA (GMIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.1582, 0.2581, 0.3227, 0.3635, 0.3882, 0.4026, 0.4108,             &
       	   0.4154, 0.4179, 0.4193, 0.4200, 0.4204, 0.4206, 0.4207 /
	DATA (GMIDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.1934, 0.3141, 0.3818, 0.4200, 0.4415, 0.4533, 0.4598,             &
       	   0.4633, 0.4651, 0.4662, 0.4667, 0.4671, 2*0.4672 /
	DATA (GMIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.1347, 0.1871, 0.2277, 0.2515, 0.2651, 0.2727, 0.2768,             &
       	   0.2790, 0.2801, 0.2808, 0.2811, 0.2812, 0.2813, 0.2814 /
	DATA (GMIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.1440, 0.2217, 0.2629, 0.2839, 0.2947, 0.3003, 0.3031,             &
       	   0.3046, 0.3054, 0.3058, 0.3060, 2*0.3061, 0.3062 /
	DATA (GMIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.1372, 0.2368, 0.3235, 0.3839, 0.4229, 0.4465, 0.4602,             &
       	   0.4679, 0.4722, 0.4745, 0.4758, 0.4764, 0.4768, 0.4770 /
	DATA (GMIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.1435, 0.2524, 0.3370, 0.3955, 0.4332, 0.4563, 0.4697,             &
       	   0.4773, 0.4815, 0.4839, 0.4851, 0.4858, 0.4861, 0.4863 /
	DATA (GMIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.4298, 0.9651, 1.6189, 2.4084, 3.2992, 4.1928, 4.9611,             &
       	   5.5095, 5.8085, 5.9069, 5.8726, 5.7674, 5.6346, 5.4944 /
	DATA (GMIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.4167, 0.8974, 1.4160, 1.9414, 2.4147, 2.7803, 3.0202,             &
      	   3.1468, 3.1954, 3.1932, 3.1676, 3.1328, 3.0958, 3.0625 /
        DATA (GMIDR (I, 1, 5), I = 1, 14)                                      &
           /0.1959, 0.3203, 0.3985, 0.4472, 0.4766, 0.4937, 0.5034,            &
            0.5088, 0.5117, 0.5134, 0.5143, 0.5147, 0.5150, 0.5152 /
        DATA (GMIDR (I, 2, 5), I = 1, 14)                                      &
           /0.2328, 0.3859, 0.4734, 0.5227, 0.5498, 0.5644, 0.5720,            &
            0.5761, 0.5781, 0.5792, 0.5797, 0.5800, 0.5802, 0.5802 /
        DATA (GMIDR (I, 1, 6), I = 1, 14)                                      &
           /0.1447, 0.2244, 0.2698, 0.2953, 0.3094, 0.3170, 0.3211,            &
            0.3233, 0.3244, 0.3250, 0.3253, 0.3255, 0.3256, 0.3256 /
        DATA (GMIDR (I, 2, 6), I = 1, 14)                                      &
           /0.1643, 0.2624, 0.3110, 0.3347, 0.3461, 0.3517, 0.3543,            &
            0.3556, 0.3562, 0.3564, 0.3565, 0.3566, 0.3566, 0.3566 /
	DATA (GMIDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 9), I = 1, 14) /14*1./

!**** -----------------------------------------------------------

      IF (present(MODIS_SCALE)) THEN
         MODIS_SCALE_ = MODIS_SCALE
      ELSE
         MODIS_SCALE_ = .FALSE.
      END IF

!FPP$ EXPAND (COEFFSIB)

      DO I=1,NCH

        ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
        LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
        DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
        DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))

        ALPHA = COEFFSIB (ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        BETA  = COEFFSIB (BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        GAMMA = COEFFSIB (GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

        GAMMA = MAX(GAMMA,0.01)

        AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
        AVISDF(I) = ALPHA-BETA                                                 &
                 + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

        ALPHA = COEFFSIB (ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        BETA  = COEFFSIB (BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        GAMMA = COEFFSIB (GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

        GAMMA = MAX(GAMMA,0.01)

        ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
        ANIRDF(I) = ALPHA-BETA                                                 &
                 + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))


     IF (MODIS_SCALE_) THEN

        ! SCALE TO MODIS VALUES (SNOW-FREE)

        fac = min(VLAI(I),1.)

        AVISDR(I) = AVISDR(I)*fac + SCALVDR(I)*(1.-fac)
             ANIRDR(I) = ANIRDR(I)*fac + SCALIDR(I)*(1.-fac)
        AVISDF(I) = AVISDF(I)*fac + SCALVDF(I)*(1.-fac)
             ANIRDF(I) = ANIRDF(I)*fac + SCALIDF(I)*(1.-fac)

     ELSE

        AVISDR(I) = AVISDR(I) * SCALVDR(I)
                  ANIRDR(I) = ANIRDR(I) * SCALIDR(I)
        AVISDF(I) = AVISDF(I) * SCALVDF(I)
                  ANIRDF(I) = ANIRDF(I) * SCALIDF(I)

     END IF

! PROTECT AGAINST BAD SCALING

	  AVISDR(I) = AMIN1( 1., AMAX1( 0., AVISDR(I) ) )
          ANIRDR(I) = AMIN1( 1., AMAX1( 0., ANIRDR(I) ) )
	  AVISDF(I) = AMIN1( 1., AMAX1( 0., AVISDF(I) ) )
          ANIRDF(I) = AMIN1( 1., AMAX1( 0., ANIRDF(I) ) )



        ENDDO

      RETURN
      END SUBROUTINE SIBALB


  ! ================================================================================
  !
  !                      Catchment diagnostic routines
  !
  ! moved to here from file "catch_diagn_routines.F90" of "lana" directory (LDAS)
  ! (except calc_soil_moist), renamed for consistency, and revised for efficiency
  ! - reichle, 3 Apr 2012
  !
  ! ================================================================================

  subroutine catch_calc_soil_moist( &
       NTILES,vegcls,dzsf,vgwmax,cdcr1,cdcr2,psis,bee,poros,wpwet, &
       ars1,ars2,ars3,ara1,ara2, &
       ara3,ara4,arw1,arw2,arw3,arw4, &
       srfexc,rzexc,catdef, &
       ar1, ar2, ar4, &
       sfmc, rzmc, prmc,  &
       werror, sfmcun, rzmcun, prmcun, &
       swsrf1out, swsrf2out, swsrf4out )

    ! Calculate diagnostic soil moisture content from prognostic
    ! excess/deficit variables.
    !
    ! On input, also check validity of prognostic excess/deficit variables
    ! and modify if necessary.  Perturbed or updated excess/deficit variables
    ! in data assimilation integrations may be unphysical.
    ! Optional output "werror" contains excess or missing water related
    ! to inconsistency.
    !
    ! Optional outputs "smfcun", "rzmcun", "prmcun" are surface,
    ! root zone, and profile moisture content for unsaturated areas only,
    ! ie. excluding the saturated area of the catchment.
    !
    ! NOTE: When calling with optional output arguments, use keywords
    !       unless arguments are in proper order!
    !
    !       Example:
    !       (don't want "werror" as output, but want "*mcun" output)
    !
    !       call catch_calc_soil_moist(         &
    !            NTILES, ...                &
    !            sfmc, rzmc, prmc,        &
    !            sfmcun=sfmc_unsat,   &
    !            rzmcun=rzmc_unsat,   &
    !            prmcun=prmc_unsat )
    !
    ! replaces moisture_sep_22_2003.f (and older moisture.f)
    !
    ! koster+reichle, Feb 5, 2004
    !
    ! revised - koster+reichle, Mar 19, 2004
    !
    ! added optional *un output - koster+reichle, Apr 6, 2004
    !
    ! removed parameter "KSNGL" ("kind=single")
    ! added output of "ar1", "ar2", "ar4"
    ! changed output arguments "sfmc", "rzmc", "prmc" to optional
    ! - reichle, 2 Apr 2012
    !
    ! ----------------------------------------------------------------


    implicit none

    integer,                    intent(in) :: NTILES
    integer, dimension(NTILES), intent(in) :: vegcls

    real,    dimension(NTILES), intent(in) :: dzsf,vgwmax,cdcr1,cdcr2
    real,    dimension(NTILES), intent(in) :: wpwet,poros,psis
    real,    dimension(NTILES), intent(in) :: bee,ars1
    real,    dimension(NTILES), intent(in) :: ars2,ars3,ara1,ara2,ara3
    real,    dimension(NTILES), intent(in) :: ara4,arw1,arw2,arw3,arw4

    real,    dimension(NTILES), intent(inout) :: srfexc, rzexc, catdef

    real,    dimension(NTILES), intent(out) :: ar1, ar2, ar4

    real,    dimension(NTILES), intent(out), optional :: sfmc, rzmc, prmc

    real,    dimension(NTILES), intent(out), optional :: werror

    real,    dimension(NTILES), intent(out), optional :: sfmcun
    real,    dimension(NTILES), intent(out), optional :: rzmcun
    real,    dimension(NTILES), intent(out), optional :: prmcun
    real,    dimension(NTILES), intent(out), optional :: swsrf1out, swsrf2out, swsrf4out    

    ! ----------------------------
    !
    ! local variables

    integer :: n

    real, parameter :: dtstep_dummy = -9999.

    real, dimension(NTILES) :: rzeq, runsrf_dummy, catdef_dummy
    real, dimension(NTILES) :: prmc_orig
    real, dimension(NTILES) :: srfmn, srfmx, swsrf1, swsrf2, swsrf4, rzi


    ! --------------------------------------------------------------------
    !
    ! compute soil water storage upon input [mm]

    do n=1,NTILES
       prmc_orig(n) =                                                 &
            (cdcr2(n)/(1.-wpwet(n))-catdef(n)+rzexc(n)+srfexc(n))
    enddo

    ! -----------------------------------
    !
    ! check limits of catchment deficit
    !
    ! increased minimum catchment deficit from 0.01 to 1. to make the
    ! check work with perturbed parameters and initial condition
    ! reichle, 16 May 01
    !
    ! IT REALLY SHOULD WORK WITH catdef > 0 (rather than >1.) ????
    ! reichle, 5 Feb 2004

    do n=1,NTILES
       catdef(n)=max(1.,min(cdcr2(n),catdef(n)))
    end do

    ! ------------------------------------------------------------------
    !
    ! check limits of root zone excess
    !
    ! calculate root zone equilibrium moisture for given catchment deficit

    call rzequil( &
         NTILES, catdef, vgwmax,    &
         cdcr1, cdcr2, wpwet, &
         ars1, ars2, ars3, ara1, ara2, ara3, ara4, &
         arw1, arw2, arw3, arw4, &
         rzeq)

    ! assume srfexc=0 and constrain rzexc appropriately
    ! (iteration would be needed to constrain srfexc and rzexc simultaneously)

    do n=1,NTILES
       rzexc(n)=max(wpwet(n)*vgwmax(n)-rzeq(n),min(vgwmax(n)-rzeq(n),rzexc(n)))
    end do

    ! this translates into:
    !
    ! wilting level < rzmc < porosity
    !
    ! or more precisely:  wpwet*vgwmax < rzeq+rzexc < vgwmax
    !
    ! NOTE: root zone moisture is not allowed to drop below wilting level

    ! -----------------------------------------------------------------
    !
    ! Call partition() for computation of surface moisture content.
    !
    ! Call to partition() also checks limits of surface excess.
    !
    ! Call partition with dtstep_dummy:
    !  In partition, dtstep is only used for a correction that
    !  puts water into runsrf (for which runsrf_dummy is used here).
    !  Also use catdef_dummy because partition() updates catdef
    !  whenever srfexc exceeds physical bounds, but this is not desired here.

    runsrf_dummy = 0.
    catdef_dummy = catdef

    call partition( &
         NTILES,dtstep_dummy,dzsf,rzexc, &
         rzeq,vgwmax,cdcr1,cdcr2, &
         psis,bee,poros,wpwet, &
         ars1,ars2,ars3, &
         ara1,ara2,ara3,ara4, &
         arw1,arw2,arw3,arw4,.false., &
         srfexc,catdef_dummy,runsrf_dummy, &
         ar1, ar2, ar4,srfmx,srfmn, &
         swsrf1,swsrf2,swsrf4,rzi &
         )

     if(present(swsrf1out)) swsrf1out = swsrf1
     if(present(swsrf2out)) swsrf2out = swsrf2
     if(present(swsrf4out)) swsrf4out = swsrf4

    ! compute surface, root zone, and profile soil moisture

    if (present(sfmc) .and. present(rzmc) .and. present(prmc)) then

       do n=1,NTILES

       sfmc(n) = poros(n) *                                           &
            (swsrf1(n)*ar1(n) + swsrf2(n)*ar2(n) + swsrf4(n)*ar4(n))

       rzmc(n) = (rzeq(n)+rzexc(n)+srfexc(n))*poros(n)/vgwmax(n)

       ! compute revised soil water storage [mm]

       prmc(n) =                                                               &
            (cdcr2(n)/(1.-wpwet(n))-catdef(n)+rzexc(n)+srfexc(n))

       ! compute error in soil water storage [mm] (if argument is present)

       if (present(werror))  werror(n)=(prmc(n)-prmc_orig(n))

       ! convert to volumetric soil moisture
       ! note: dzpr = (cdcr2/(1-wpwet)) / poros

       prmc(n) = prmc(n)*poros(n) / (cdcr2(n)/(1.-wpwet(n)))


       ! check for negative soil moisture

       if ( (sfmc(n)<.0) .or. (rzmc(n)<.0) .or. (prmc(n)<.0) ) then

          write (*,*) 'FOUND NEGATIVE SOIL MOISTURE CONTENT.... stopping'
          write (*,*) n, sfmc(n), rzmc(n), prmc(n)
          stop
       end if

       ! compute moisture content in unsaturated areas [m3/m3] (if arg present)

       if (ar1(n)<1.) then

          if (present(prmcun))  prmcun(n)=(prmc(n)-poros(n)*ar1(n))/(1.-ar1(n))
          if (present(rzmcun))  rzmcun(n)=(rzmc(n)-poros(n)*ar1(n))/(1.-ar1(n))
          if (present(sfmcun))  sfmcun(n)=(sfmc(n)-poros(n)*ar1(n))/(1.-ar1(n))

       else

          if (present(prmcun))  prmcun(n)=poros(n)
          if (present(rzmcun))  rzmcun(n)=poros(n)
          if (present(sfmcun))  sfmcun(n)=poros(n)

       end if

    enddo

    end if

  return

  end subroutine catch_calc_soil_moist

  ! *******************************************************************

  subroutine catch_calc_subtile2tile( NTILES, ar1, ar2, asnow, subtile_data, tile_data )

    ! average from subtile space to tile-average
    !
    ! subtile areas correspond to subtile_data as follows:
    !
    !   ar1    <--->  subtile_data(:,1)    [saturated]
    !   ar2    <--->  subtile_data(:,2)    [transpiring]
    !   ar4    <--->  subtile_data(:,3)    [wilting]
    !   asnow  <--->  subtile_data(:,4)    [snow-covered]
    !
    ! reichle, Feb  2, 2011

    implicit none

    integer,                      intent(in)  :: NTILES

    real,    dimension(NTILES),   intent(in)  :: ar1, ar2, asnow

    real,    dimension(NTILES,4), intent(in)  :: subtile_data

    real,    dimension(NTILES),   intent(out) :: tile_data

    ! ----------------------------

    ! compute non-snow average

    tile_data = ar1*subtile_data(:,1) + ar2*subtile_data(:,2) &
         + (1. - ar1 - ar2)*subtile_data(:,3)

    ! mix in snow-covered area

    tile_data = asnow*subtile_data(:,4) + (1. - asnow)*tile_data

  end subroutine catch_calc_subtile2tile

  ! *******************************************************************

      subroutine gndtmp(dts,phi,zbar,thetaf,fh21,ht,xfice,tp, FICE)
! using a diffusion equation this code generates ground temperatures
! with depth given t1
!            *****************************************
!        input
!        dts     timestep in seconds
!        phi     porosity
!        t1      terrestrial (layer 1) temperature in deg C
!        zbar    mean depth to the water table.
!        thetaf  mean vadose zone soil moisture factor (0-1)
!        output,
!        ht      heat content in layers 2-7
!        tp      ground temperatures in layers 2-7
!        tdeep   the temperature of the "deep"
!        f21     heat flux between layer 2 and the terrestrial layer (1)
!        df21    derivative of f21 with respect to temperature
!        xfice   a total soil column ice factor (0-1)
!             ***********************************

      REAL, INTENT(IN) :: DTS, phi, ZBAR, THETAF, FH21

      REAL, INTENT(INOUT), DIMENSION(*) :: HT

      REAL, INTENT(OUT) :: XFICE
      REAL, INTENT(OUT), DIMENSION(*) :: TP, FICE

      INTEGER L, LSTART, K
      REAL, DIMENSION(N_GT) :: ZC, SHC, XKLH
      REAL, DIMENSION(N_GT+1) :: FH, ZB
      REAL SHW0, SHI0, SHR0, WS, XW, A1, TK1, A2, TK2, TK3, TKSAT,      &
           XWI, XD1, XD2, XKTH, TKDRY

      !data dz/0.0988,0.1952,0.3859,0.7626,1.5071,10.0/
      !DATA PHI/0.45/, FSN/3.34e+8/, SHR/2.4E6/

! initialize parameters

      shw0=SHW*1000. ! PER M RATHER THAN PER KG
      shi0=SHI*1000. ! PER M RATHER THAN PER KG
      shr0=SHR*1000. ! PER M RATHER THAN PER KG [kg of water equivalent density]

!----------------------------------
! initialize fice in ALL components
! reichle, July 8, 2002

      do l=1,N_GT
        fice(l)=0.0
        enddo
!----------------------------------

! calculate the boundaries, based on the layer thicknesses(DZGT)

      zb(1)=-DZTC    ! Bottom of surface layer, which is handled outside
                    ! this routine.
      do l=1,N_GT
        zb(l+1)=zb(l)-DZGT(l)
        shc(l)=shr0*(1-phi)*DZGT(l)
        enddo
      do l=1,N_GT
        zc(l)=0.5*(zb(l)+zb(l+1))
        enddo


! evaluates the temperatures in the soil layers based on the heat values.
!             ***********************************
!             input:
!             xw - water in soil layers, m
!             ht - heat in soil layers
!             fsn - heat of fusion of water   J/m
!             shc - specific heat capacity of soil
!             shi - specific heat capacity of ice
!             shw - specific heat capcity of water
!             snowd - snow depth, equivalent water m
!             output:
!             tp - temperature of layers, c
!             fice - fraction of ice of layers
!             pre - extra precipitation, i.e. snowmelt, m s-1
!             snowd - snow depth after melting, equivalent water m.
!             ***********************************
! determine fraction of ice and temp of soil layers based on layer
! heat and water content

      do 10 k=1,N_GT
        ws=phi*DZGT(k)  ! PORE SPACE IN LAYER
        xw=0.5*ws   ! ASSUME FOR THESE CALCULATIONS THAT THE PORE SPACE
                    ! IS ALWAYS HALF FILLED WITH WATER.  XW IS THE
                    ! AMOUNT OF WATER IN THE LAYER.
        tp(k)=0.
        fice(k)= AMAX1( 0., AMIN1( 1., -ht(k)/(fsn*xw) ) )

        IF(FICE(K) .EQ. 1.) THEN
            tp(k)=(ht(k)+xw*fsn)/(shc(k)+xw*shi0)
          ELSEIF(FICE(K) .EQ. 0.) THEN
            tp(k)=ht(k)/(shc(k)+xw*shw0)
          ELSE
            TP(K)=0.
          ENDIF

    10  continue

! evaluates:  layer thermal conductivities
! *****************************************
!             from farouki(cold regions sci and tech, 5, 1981,
!             67-75) the values for tk1,tk2,tk3 are as follows:
!             tk2=2.2**(phi(l,ibv)-xw), tk3=.57**xw, and
!             tk1=3**(1-phi(l,ibv) for %sand<.5, and
!             for computation purposes i have fit these eqs.
!             to 2nd order polynomials.
!             ***********************************
!             input:
!             sklh - soil heat conductivities of layers
!             zb - soil layer boundaries, m
!             zc - soil layer centers, m
!             DZGT - layer thickness, m
!             w - soil water content, m
!             phi - soil porosity, dimensionless
!             q - % sand, silt, clay, peat
!             fice - fraction of ice in layers
!             output:
!             xklh - thermal conductivity, w m-2 k-1
!             ***********************************



! lets get the thermal conductivity for the layers

      do k=1,N_GT

         a1=1-phi              ! ROCK FRACTION
         tk1=1.01692+a1*(0.89865+1.06211*a1)
         xw=phi*(1.-fice(k))   ! FOR SATURATED SOIL, XW HERE IS
                               !   THE LIQUID WATER FRACTION
         a2=phi-xw             ! FOR SATURATED SOIL, THE ICE FRACTION

         tk2=1.00543+a2*(0.723371+.464342*a2)
         tk3=0.998899+xw*(-0.548043+0.120291*xw)
         tksat=tk1*tk2*tk3

         xwi=1.0
         if (zbar .le. zb(k+1))then
             xwi=thetaf
           elseif (zbar .ge. zb(k+1) .and. zbar .le. zb(k))then
             xd1=zb(k)-zbar
             xd2=zbar-zb(k+1)
             xwi=((xd1*thetaf)+xd2)/(xd1+xd2)
           endif

         xwi=min(xwi,1.)
         tkdry=0.226 ! = .039*0.45^(-2.2), from Farouki, p. 71
         xklh(k)=(tksat-tkdry)*xwi + tkdry

         enddo

! evaluates heat flux between layers due to heat diffussion
!             ***********************************
!             input:
!             zb - soil layer boundaries, m
!             zc - soil layer centers, m
!             DZGT - layer thickness, m
!             fice - fraction of ice in layers
!             tp - temperature of layers, c
!             shw - specific heat of water
!             shi - specific heat of ice
!             fsn - heat of fusion    J/m
!             output:
!             fh - heat flux between layers
!             ***********************************


! total heat flux is via diffusion along the temperature gradient
      fh(N_GT+1)=0.
      fh(1)=fh21
      do k=2,N_GT
! THIS xkth is NEW (ie., Agnes corrected) - it should be fixed in all
! codes I'm using
         xkth=((zb(k)-zc(k-1))*xklh(k-1)+(zc(k)-zb(k))*xklh(k))                &
               /(zc(k)-zc(k-1))
         fh(k)=-xkth*(tp(k-1)-tp(k))/(zc(k-1)-zc(k))
         enddo

! update the heat contents in the model layers; ht(l)
! IF THERE'S SNOW THIS WILL HAVE TO BE MODIFIED L=1,N

      do k=1,N_GT
         ht(k)=ht(k)+(fh(k+1)-fh(k))*dts
         enddo



! evaluates the temperatures in the soil  layers based on the heat
! values.
!             ***********************************
!             input:
!             xw - water in soil layers, m
!             ht - heat in soil layers
!             fsn - heat of fusion of water    J/m
!             shc - specific heat capacity of soil
!             shi - specific heat capacity of ice
!             shw - specific heat capcity of water
!             snowd - snow depth, equivalent water m
!             output:
!             tp - temperature of layers, c
!             fice - fraction of ice of layers
!             pre - extra precipitation, i.e. snowmelt, m s-1
!             snowd - snow depth after melting, equivalent water m.
!             ***********************************
! determine fraction of ice and temp of soil layers based on layer
! heat and water content
      do 1000 k=1,N_GT
         ws=phi*DZGT(k)          ! saturated water content
!            xl=l
!            xw=(1/(7-xl))*ws
         xw=0.5*ws             ! For calculations here, assume soil
                               ! is half full of water.
         fice(k)=AMAX1( 0., AMIN1( 1., -ht(k)/(fsn*xw) ) )

         IF(FICE(K) .EQ. 1.) THEN
               tp(k)=(ht(k)+xw*fsn)/(shc(k)+xw*shi0)
            ELSEIF(FICE(K) .EQ. 0.) THEN
               tp(k)=ht(k)/(shc(k)+xw*shw0)
            ELSE
               TP(K)=0.
            ENDIF

 1000    continue

! determine the value of xfice
      xfice=0.0

      lstart=N_GT
      DO L=N_GT,1,-1
         IF(ZBAR .GE. ZB(L+1))THEN
            LSTART=L
            ENDIF
         ENDDO

      do l=lstart,N_GT
         xfice=xfice+fice(l)
         enddo
      xfice=xfice/((N_GT+1)-lstart)

      Return

      end subroutine gndtmp

 ! *******************************************************************

  subroutine catch_calc_tp( NTILES, poros, ghtcnt, tp, fice )

    ! Diagnose soil temperatures tp (all_layers, all tiles) from
    ! prognostic ground heat contents
    !
    ! return temperature in degree CELSIUS!!! [for consistency w/ catchment.F90]
    !
    ! (could also return fractional ice content fice)
    !
    ! GDL, 25Mar11, based on code by Randy Koster
    !
    ! reichle, 13 May 2011
    ! reichle,  2 Apr 2012 - revised for use without catch_types structures
    !
    !-----------------------------------------------------------------

    implicit none

    integer,                                      intent(in)            :: NTILES

    real,                 dimension(     NTILES), intent(in)            :: poros

    real,                 dimension(N_gt,NTILES), intent(in)            :: ghtcnt

    real,                 dimension(N_gt,NTILES), intent(out)           :: tp

    real,                 dimension(N_gt,NTILES), intent(out), optional :: fice

    ! Local variables

    real, parameter              :: shw0 = SHW   *1000. ! unit change J/kg/K -> J/m/K
    real, parameter              :: shi0 = SHI   *1000. ! unit change J/kg/K -> J/m/K
    real, parameter              :: shr0 = SHR   *1000. ! unit change J/kg/K -> J/m/K [where "per kg" is something like "per kg of water equivalent density"]

    integer                      :: n, k

    real                         :: phi, WS, XW
    real, dimension(N_gt)        :: SHC
    real, dimension(N_gt,NTILES) :: FICE_TMP

    ! ---------------------------------------------------------------------------

    ! initialize

    tp = 0.

    do n=1,NTILES

       if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
          phi=poros(n)
       else
          phi=PHIGT
       end if

       do k=1,N_gt

          shc(k) = shr0*(1-phi)*DZGT(k)

          ws = phi*DZGT(k) ! PORE SPACE IN LAYER

          xw = 0.5*ws      ! ASSUME FOR THESE CALCULATIONS THAT THE PORE SPACE
                           ! IS ALWAYS HALF FILLED WITH WATER.  XW IS THE
                           ! AMOUNT OF WATER IN THE LAYER.

          FICE_TMP(k,n)= AMAX1( 0., AMIN1( 1., -ghtcnt(k,n)/(fsn*xw) ) )

          IF     (FICE_TMP(K,n) .EQ. 1.) THEN

             tp(k,n) = (ghtcnt(k,n)+xw*fsn)/(shc(k)+xw*shi0)   ! Celsius

          ELSEIF (FICE_TMP(K,n) .EQ. 0.) THEN

             tp(k,n) = (ghtcnt(k,n)       )/(shc(k)+xw*shw0)   ! Celsius

          ELSE

             tp(k,n) = 0.                                              ! Celsius

          ENDIF

       end do

    end do

    if (present(fice))  fice = FICE_TMP

  end subroutine catch_calc_tp

  ! *******************************************************************

  subroutine catch_calc_wtotl( NTILES,                         &
       cdcr2, wpwet, srfexc, rzexc, catdef, capac, wesnn,      &
       wtotl )

    ! compute total water stored in land tiles
    !
    ! reichle,  4 Jan 2012
    ! reichle,  2 Apr 2012 - revised for use without catch_types structures
    !
    ! ----------------------------------------------------------------

    implicit none

    integer,                          intent(in)  :: NTILES

    real,    dimension(       NTILES), intent(in)  :: cdcr2, wpwet
    real,    dimension(       NTILES), intent(in)  :: srfexc, rzexc, catdef, capac
    real,    dimension(N_snow,NTILES), intent(in)  :: wesnn

    real,    dimension(       NTILES), intent(out) :: wtotl

    ! ----------------------------
    !
    ! local variables

    integer :: n

    ! ----------------------------------------------------------------

    do n=1,NTILES

       ! total water = soil water holding capacity +/- soil water excess/deficit
       !               + canopy interception + snow water

       wtotl(n) =                                  &
            ( cdcr2(n)/(1.-wpwet(n))               &
            - catdef(n)                            &
            + rzexc(n)                             &
            + srfexc(n)                            &
            + capac(n)                             &
            + sum( wesnn(1:N_snow,n) ) )

    end do

  end subroutine catch_calc_wtotl

  ! *******************************************************************

  ! *******************************************************************

  subroutine catch_calc_ght( dzgt, poros, tp, fice, ghtcnt )

    ! Invert (model diagnostic) soil temperature and ice fraction
    ! into (model prognostic) ground heat content
    !
    ! Input soil temperature is in deg CELSIUS!!!
    !
    ! subroutine is for a single layer only and single tile only!!!
    !
    ! reichle, 13 Oct 2014
    !
    !------------------------------------------------------------------

    implicit none

    real,    intent(in)    :: dzgt      ! soil temperature layer depth [m?]
    real,    intent(in)    :: poros     ! soil porosity
    real,    intent(in)    :: tp        ! soil temperature [deg CELSIUS]
    real,    intent(in)    :: fice      ! soil ice fraction

    real,    intent(out)   :: ghtcnt    ! ground heat content [J?]

    ! Local variables

    real    :: phi, ws, xw, shc, shw0, shi0, shr0

    ! initialize parameters

    shw0=SHW*1000. ! PER M RATHER THAN PER KG
    shi0=SHI*1000. ! PER M RATHER THAN PER KG
    shr0=SHR*1000. ! PER M RATHER THAN PER KG [kg of water equivalent density]

    ! ---------------------------------------------------------------------------

    if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
       phi=poros
    else
       phi=PHIGT
    end if

    shc = shr0*(1-phi)*DZGT

    ws  = phi*DZGT           ! PORE SPACE IN LAYER

    xw  = 0.5*ws             ! ASSUME FOR THESE CALCULATIONS THAT THE PORE SPACE
                             ! IS ALWAYS HALF FILLED WITH WATER.  XW IS THE
                             ! AMOUNT OF WATER IN THE LAYER.

    IF     (tp<0.0) THEN     ! water in soil is fully frozen

       ghtcnt = tp*(shc + xw*shi0) - FSN*xw

    ELSEIF (tp>0.0) THEN     ! water in soil is fully thawed

       ghtcnt = tp*(shc + xw*shw0)

    ELSE                     ! water in soil is partially frozen

       ghtcnt = -fice*(FSN*xw)

    END IF

  end subroutine catch_calc_ght

  ! ********************************************************************

  subroutine catch_calc_FT( NTILES, asnow, tp1, tsurf_excl_snow, FT, Teff )

    ! Diagnose "landscape" freeze/thaw (FT) state of the Catchment/StieglitzSnow model.
    !
    ! The landscape FT state is determined via the snow cover fraction and
    ! an "effective" temperature (Teff), computed as the weighted average
    ! of the top-layer soil temperature and the surface temperature
    ! in the absence of snow (see Farhadi et al., 2014, JHM, section 3).
    !
    ! Input soil temperature is in CELSIUS!!!
    !
    ! Constants:
    !
    !  CATCH_FT_WEIGHT_TP1      : weight for tp1 vs. Tsurf_excl_snow in Teff,
    !                             determines effective depth associated with FT state
    !                             = 1. - alpha, with alpha as in Farhadi et al. 2014,
    !
    !  CATCH_FT_THRESHOLD_TEFF  : FT threshold for Teff              [Kelvin]
    !  CATCH_FT_THRESHOLD_ASNOW : FT threshold for asnow
    !
    ! Inputs:
    !
    !  asnow              : snow cover fraction
    !  tp1                : top-layer soil temperature               [Celsius]
    !  tsurf_excl_snow    : surface temperature in absence of snow   [Kelvin]
    !  Teff               : effective temperature (optional output)  [Kelvin]
    !
    ! Outputs:
    !
    !  FT                 : FT state (0.=thawed, 1.=frozen)
    !
    ! reichle, 14 Oct 2014
    !
    !------------------------------------------------------------------

    implicit none

    integer,                 intent(in)  :: NTILES

    real, dimension(NTILES), intent(in)  :: asnow           ! snow cover fraction
    real, dimension(NTILES), intent(in)  :: tsurf_excl_snow ! ar1*tc1+ar2*tc2+ar4*tc4
    real, dimension(NTILES), intent(in)  :: tp1             ! top-layer soil temp [CELSIUS]

    real, dimension(NTILES), intent(out) :: FT

    real, dimension(NTILES), intent(out), optional :: Teff

    ! Local variables

    real                    :: w

    real, dimension(NTILES) :: Teff_tmp                     ! [Kelvin]

    ! ------------------------------------------------------------------

    w = CATCH_FT_WEIGHT_TP1

    ! compute snow-free "effective" temperature  [Kelvin]

    Teff_tmp = w*(tp1 + TF) + (1.-w)*Tsurf_excl_snow

    ! Note: For the thawed state use Teff_threshold with ">" and
    !        not ">=" as in Farhadi et al, 2014 because at exactly
    !        0 deg Celsius the soil may still be partially frozen.

    where (                                                   &
         (Teff_tmp  > CATCH_FT_THRESHOLD_TEFF ) .and.         &
         (asnow     < CATCH_FT_THRESHOLD_ASNOW)           )

       FT = 0.   ! thawed

    elsewhere

       FT = 1.   ! frozen

    end where

    if (present(Teff))  Teff = Teff_tmp

  end subroutine catch_calc_FT

  ! ********************************************************************

  subroutine dampen_tc_oscillations( dtstep, tair, tc_old, tc_new_in, &
       tc_new_out, dtc_new )

    implicit none

    ! apply corretion to TC to dampen surface energy balance oscillations

    ! reichle, 4 Apr 2014

    real, intent(in)  :: dtstep     ! model time step [s]
    real, intent(in)  :: tair       ! air temperature
    real, intent(in)  :: tc_old     ! Tc at beginning of time step
    real, intent(in)  :: tc_new_in  ! proposed  Tc at end of time step
    real, intent(out) :: tc_new_out ! corrected Tc at end of time step
    real, intent(out) :: dtc_new    ! change in tc_new from this subroutine

    ! local variables

    real, parameter :: xover_frac             = 0.1   ! [dimensionless]
    real, parameter :: xover_max_dTc_per_hour = 1.    ! [Kelvin/hour]

    real :: xover_max_dTc, dTc_tmp1, dTc_tmp2

    ! --------------------------------------------------------------------

    ! maximum Tc change in Kelvin that is allowed if tc1 changes across tm
    !   [this might better be done only once, at the beginning of subroutine
    !    catchment() and outside of the loop through tiles]

    xover_max_dTc = xover_max_dTc_per_hour*dtstep/3600.

    ! establish if tc changed across tair

    if ( ((tc_old < tair) .and. (tair < tc_new_in)) .or.                 &
         ((tc_old > tair) .and. (tair > tc_new_in))         ) then

       ! limit amount of tc change beyond tair (dTc_tmp) to the smaller of
       ! (i) a fraction of tcNew_minus_tm and (ii) a max value in Kelvin

       dTc_tmp1    = tc_new_in - tair

       dTc_tmp2    = min( abs(xover_frac*dTc_tmp1), xover_max_dTc )  ! never negative

       tc_new_out  = tair + sign( dTc_tmp2, dTc_tmp1 )

       dtc_new     = tc_new_out - tc_new_in

    else

       ! tc_new remains unchanged

       tc_new_out  = tc_new_in

       dtc_new     = 0.

    end if

  end subroutine dampen_tc_oscillations

  ! ********************************************************************

  subroutine lsmroutines_echo_constants(logunit)

    ! moved to here from catch_constants.F90, reichle, 14 Aug 2014

    ! reichle, 10 Oct 2008

    implicit none

    integer, intent(in) :: logunit

    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)
    write (logunit,*) 'lsmroutines_echo_constants():'
    write (logunit,*)
    write (logunit,*) 'PIE      = ', PIE
    write (logunit,*) 'ALHE     = ', ALHE
    write (logunit,*) 'ALHM     = ', ALHM
    write (logunit,*) 'ALHS     = ', ALHS
    write (logunit,*) 'TF       = ', TF
    write (logunit,*) 'RGAS     = ', RGAS
    write (logunit,*) 'SHW      = ', SHW
    write (logunit,*) 'SHI      = ', SHI
    write (logunit,*)
    write (logunit,*) 'N_snow        = ', N_snow
    write (logunit,*) 'N_gt          = ', N_gt
    write (logunit,*)
    write (logunit,*) 'RHOFS         = ', RHOFS
    write (logunit,*) 'SNWALB_VISMAX = ', SNWALB_VISMAX
    write (logunit,*) 'SNWALB_NIRMAX = ', SNWALB_NIRMAX
    write (logunit,*) 'SLOPE         = ', SLOPE
    write (logunit,*) 'MAXSNDEPTH    = ', MAXSNDEPTH
    write (logunit,*) 'DZ1MAX        = ', DZ1MAX
    write (logunit,*)
    write (logunit,*) 'SHR           = ', SHR
    write (logunit,*) 'EPSILON       = ', EPSILON
    write (logunit,*)
    write (logunit,*) 'N_sm          = ', N_sm
    write (logunit,*)
    write (logunit,*) 'SCONST        = ', SCONST
    write (logunit,*) 'CSOIL_2       = ', CSOIL_2
    write (logunit,*) 'LAND_FIX      = ', LAND_FIX
    write (logunit,*) 'WEMIN         = ', WEMIN
    write (logunit,*) 'AICEV         = ', AICEV
    write (logunit,*) 'AICEN         = ', AICEN
    write (logunit,*) 'FLWALPHA      = ', FLWALPHA
    write (logunit,*) 'ASTRFR        = ', ASTRFR
    write (logunit,*) 'STEXP         = ', STEXP
    write (logunit,*) 'RSWILT        = ', RSWILT

    write (logunit,*) 'C_CANOP (catchCN only)  = ', C_CANOP
    write (logunit,*)
    write (logunit,*) 'DZTC          = ', DZTC
    write (logunit,*)
    write (logunit,*) 'SATCAPFR      = ', SATCAPFR
    write (logunit,*)
    write (logunit,*) 'DZGT          = ', DZGT
    write (logunit,*) 'PHIGT         = ', PHIGT
    write (logunit,*) 'ALHMGT        = ', ALHMGT
    write (logunit,*)
    write (logunit,*) 'end lsmroutines_echo_constants()'
    write (logunit,*)
    write (logunit,*) '-----------------------------------------------------------'
    write (logunit,*)

  end subroutine lsmroutines_echo_constants

  ! ********************************************************************

  SUBROUTINE irrigation_rate (IRRIG_METHOD,                                    &
                NTILES, AGCM_HH, AGCM_MI, AGCM_S, lons, IRRIGFRAC, PADDYFRAC,  &
                CLMPT,CLMST, CLMPF, CLMSF, LAIMAX, LAIMIN, LAI,                &
                POROS, WPWET, VGWMAX, RZMC, IRRIGRATE)

    ! !DESCRIPTION:
    !
    ! NOTE: This is an experimental feature under development.
    !
    ! Calculate water requirement and apply the amount to precipitation.
    !
    ! Irrigate when available root zone soil moisture falls below tunable
    !        irrigation threshold parameter.
    ! Below GRIPC irrigated data provide fractions of croplands and paddy croplands.
    ! The irrigation model is applied on a tile if:
    ! (1) the irrigated fraction of the tile is greater than 0. AND
    ! (2) primary or secondary type in the tile is CLM4 type 16 (cropland) AND
    ! (3) LAI exceeds the LAI theshhold  (60% of LAI range)
    !
    ! GRIPC croplands and paddy croplands fractions determine whether to apply
    !    either sprinkler or flood OR both irrigation methods. Each method has
    !    its own local start times, durations and irrigation threshold parameters.
    !
    ! We assume plants need available soil moisture stay above 1/3 of soil moisture range
    !    [ wilting - saturation]
    ! Irrigation amount is scaled to grid total crop fraction when intensity
    ! is less than the fraction.  Irrigation is expanded to non-crop, non-forest,
    ! non-baresoil/urban tiles if intensity exceeds grid total crop fraction.
    ! In latter case, scaled irrigation is applied to grassland first,
    ! then further applied over the rest of tiles equally if the intensity
    ! exceeds grassland fraction as well.
    !
    ! Optionally efficiency correction is applied to account for field loss.
    !
    ! REVISION HISTORY:
    !
    ! Aug 2018: Sarith Mahanama ; Version 1 adapted from LIS subroutine clsmf25_getirrigationstates.F90

    implicit none

    ! INPUTS
    ! ------
    integer, intent (in)                       :: IRRIG_METHOD, NTILES, AGCM_HH, AGCM_MI, AGCM_S
    real   , intent (in), dimension (ntiles)   :: lons, IRRIGFRAC, PADDYFRAC, LAIMAX,  &
                  LAIMIN, LAI, CLMPT,CLMST, CLMPF, CLMSF, POROS, WPWET, VGWMAX, RZMC
    ! IRRIG_METHOD : 0 sprinkler and flood combined; 1 sprinkler irrigation ; 2 flood irrigation
    ! AGCM_HH / AGCM_MI / AGCM_S/ lons : Current hour, minute, second (UTC) and longitude

    ! Irrigation hotspots : Using the Global Rain-Fed, Irrigated, and Paddy Croplands (GRIPC) Dataset  (Salmon et al., 2015)
    !       Salmon JM, Friedl MA, Frolking S, Wisser D and Douglas EM: Global rain-fed, irrigated,
    !             and paddy croplands: A new high resolution map derived from remote sensing, crop
    !             inventories and climate data, Int. J. Appl. Earth Obs. Geoinf, 38, 321–334,
    !             doi:10.1016/j.jag.2015.01.014, 2015.

    ! IRRIGFRAC                        : Fraction of irrigated croplands [-] = total number of 500m irrigated croplands pixels in the tile /
    !                                                                          total number of 500m pixels in the tile
    ! PADDYFRAC                        : Fraction of paddy croplands [-]     = total number of 500m paddy croplands pixels in the tile /
    !                                                                          total number of 500m pixels in the tile
    ! LAIMAX / LAIMIN / LAI            : Maximum, minimum and current Leaf Area Index
    ! CLMPT / CLMST                    : CLM4 primary and secondary types (Note type 16 is cropland)
    ! CLMPF / CLMSF                    : CLM4 fractions of primary and secondary types
    ! POROS / WPWET / VGWMAX / RZMC    : porosity [m3/m3], wilting point wetness [-], maximum and current root zone soil moisture content [m3/m3]

    ! ONLY output
    ! -----------
    real   , intent (out), dimension (ntiles)  :: IRRIGRATE

    real, parameter      :: efcor = 76.0           ! Efficiency Correction (%)

    ! Sprinkler parameters
    ! --------------------
    real, parameter      :: otimess          = 6.0 ! local trigger check start time [hour]
    real, parameter      :: irrhrs           = 4.0 ! duration of irrigation hours
    real, parameter      :: sprinkler_thersh = 0.5 ! soil moisture threshhold to trigger sprinkler irrigation

    ! Drip parameters (not currently implemented)
    ! -------------------------------------------
    real, parameter      :: otimeds = 6.0 ! local trigger check start time [hour]
    real, parameter      :: irrhrd = 12.0 ! duration of irrigation hours

    ! Flood parameters
    ! ----------------
    real, parameter      :: otimefs      = 6.0  ! local trigger check start time [hour]
    real, parameter      :: irrhrf       = 1.0  ! duration of irrigation hours
    real, parameter      :: flood_thersh = 0.25 ! soil moisture threshhold to trigger flood irrigation

    ! local vars
    ! ----------
    real    :: smcwlt, smcref, smcmax, asmc, laithresh, laifac, RZDEP, vfrac, ma,  &
         otimee, irrig_thresh, IrrigScale, s_irate, f_irate, local_long, local_hour
    integer :: n, t, vtyp

    IRRIGRATE (:) = 0.

    TILE_LOOP : DO N = 1, NTILES

       local_long = 180. * lons(n) / PIE                                             ! local logitude [degrees]
       local_hour = AGCM_HH + AGCM_MI / 60. + AGCM_S / 3600. + 12.* local_long /180. ! local time [hours]
       if (local_hour >= 24.) local_hour =  local_hour - 24.
       if (local_hour <   0.) local_hour =  local_hour + 24.

       laithresh  = laimin (n) + 0.60 * (laimax (n) - laimin (n))
       if(laimax (n) /= laimin (n)) then
          laifac = (lai(n) - laimin (n)) / (laimax(n) - laimin(n))
       else
          laifac = 0.
       endif

       RZDEP      = laifac * VGWMAX (n) / poros (n)                                  ! root zone depth [mm]
       smcwlt     = RZDEP * wpwet (n) * poros (n)                                    ! RZ soil moisture content at wilting point [mm]
       smcref     = RZDEP * (wpwet (n) + 0.333 * (1. - wpwet (n))) * poros(n)        ! RZ reference soil moisture content [mm]
       smcmax     = RZDEP * poros (n)                                                ! RZ soil moisture at saturatopm [mm]
       asmc       = RZDEP * rzmc (n)                                                 ! actual RZ soil moisture content [mm]

       CHECK_IRRIG_INTENSITY : IF ((IRRIGFRAC(N) + PADDYFRAC(N)) > 0.) THEN

          s_irate = 0.
          f_irate = 0.

          TWO_CLMTYPS : DO t = 1, 2

             if (t == 1) then
                ! Primary CLM fraction
                vtyp  = NINT (CLMPT (n))
                vfrac = CLMPF (n)
             endif

             if (t == 2) then
                ! Secondary CLM fraction
                vtyp  = NINT (CLMST (n))
                vfrac = CLMSF (n)
             endif

             CHECK_CROP_LAITHRESH : IF ((vtyp == 16).and.(vfrac > 0.).and.(lai(n) >= laithresh).and.(laifac > 0.)) THEN

                !-----------------------------------------------------------------------------
                !     Compute irrigation scale parameter :
                !     Scale the irrigation intensity to the crop % when intensity < crop%.
                !     Expand irrigation for non-crop, non-forest when intensity > crop %
                !     in preference order of grassland first then rest.
                !-----------------------------------------------------------------------------

                IF ((IRRIGFRAC(N) + PADDYFRAC(N)) <  vfrac) THEN
                   IrrigScale = vfrac / (IRRIGFRAC(N) + PADDYFRAC(N))
                ELSE
                   IrrigScale = 1.
                ENDIF

                !-----------------------------------------------------------------------------
                !     Get the root zone moisture availability to the plant
                !-----------------------------------------------------------------------------

                if(smcref.ge.smcwlt) then
                   ma = (asmc - smcwlt) /(smcref - smcwlt)
                else
                   ma = -1
                endif

                SELECT CASE (IRRIG_METHOD)

                !--------------------------------------------------------------------------------------------------------------------------
                !     IRRIGRATE : irrigation rate required to fill up water deficit before END OF IRRIGATION PERIOD  (otimee - local_hour)
                !--------------------------------------------------------------------------------------------------------------------------

                CASE (0)
                   ! SPRINKLER AND FLOOD IRRIGATION COMBINED
                   ! ---------------------------------------
                   C_SPRINKLER : IF((IRRIGFRAC (N) > 0.).and.(local_hour >= otimess).and. (local_hour < otimess + irrhrs)) THEN
                      otimee = otimess + irrhrs ;  irrig_thresh = sprinkler_thersh
                      IF ((ma  <= irrig_thresh).and.(ma.ge.0)) THEN
                         s_irate = crop_water_deficit (IRRIGFRAC (N) * irrigScale, asmc, smcref, efcor) / (otimee - local_hour) /3600.0
                      ENDIF
                   ENDIF C_SPRINKLER

                   C_FLOOD : IF((PADDYFRAC (N) > 0.).and.(local_hour >= otimefs).and. (local_hour <= otimefs + irrhrf)) THEN
                      otimee = otimefs + irrhrf ; irrig_thresh = flood_thersh
                      IF ((ma  <= irrig_thresh).and.(ma.ge.0)) THEN
                         f_irate = crop_water_deficit (PADDYFRAC (N) * irrigScale, asmc, smcref, efcor) / (otimee - local_hour) /3600.0
                      ENDIF
                   ENDIF C_FLOOD

                   IRRIGRATE (N) = (s_irate * IRRIGFRAC (N) + f_irate * PADDYFRAC (N)) / (IRRIGFRAC(N) + PADDYFRAC(N)) ! weighted averaged sprinkler + flood

                CASE (1)
                   ! SPRINKLER IRRIGATION ONLY
                   ! -------------------------
                   SPRINKLER : IF(((IRRIGFRAC (N) + PADDYFRAC (N)) > 0.).and.(local_hour >= otimess).and. (local_hour < otimess + irrhrs)) THEN
                      otimee = otimess + irrhrs ;  irrig_thresh = sprinkler_thersh
                      IF ((ma  <= irrig_thresh).and.(ma.ge.0)) THEN
                         IRRIGRATE (N) = crop_water_deficit ((IRRIGFRAC (N) + PADDYFRAC (N)) * irrigScale, asmc, smcref, efcor) / &
                              (otimee - local_hour) /3600.0
                      ENDIF
                   ENDIF SPRINKLER

                CASE (2)
                   ! FLOOD IRRIGATION ONLY
                   ! ---------------------
                   FLOOD : IF(((IRRIGFRAC (N) + PADDYFRAC (N)) > 0.).and.(local_hour >= otimefs).and. (local_hour <= otimefs + irrhrf)) THEN
                      otimee = otimefs + irrhrf ; irrig_thresh = flood_thersh
                      IF ((ma  <= irrig_thresh).and.(ma.ge.0)) THEN
                         IRRIGRATE (N) = crop_water_deficit ((IRRIGFRAC (N) +PADDYFRAC (N)) * irrigScale, asmc, smcref, efcor) / &
                              (otimee - local_hour) /3600.0
                      ENDIF
                   ENDIF FLOOD

                CASE DEFAULT
                   print *, 'IN IRRIGATION_RATE : IRRIGATION_METHOD can  be 0, 1, or 2'
                   call exit(1)
                END SELECT
             END IF CHECK_CROP_LAITHRESH
          END DO TWO_CLMTYPS
         END IF CHECK_IRRIG_INTENSITY
    END DO TILE_LOOP

  END SUBROUTINE irrigation_rate

  ! ********************************************************************

  REAL FUNCTION crop_water_deficit (IrrigScale, asmc, smcref, efcor)

    implicit none

    real, intent (in)  :: IrrigScale, asmc, smcref, efcor
    real               :: twater

    twater             = smcref - asmc
    twater             = twater * IrrigScale            ! Scale the irrigation intensity
    crop_water_deficit = twater*(100.0/(100.0-efcor))   ! Apply efficiency correction

  END FUNCTION crop_water_deficit

  ! ********************************************************************

END MODULE lsm_routines
