MODULE lsm_routines

  ! The module contains subroutines that are shared by Catchment and CatchmentCN.
  ! Sarith,  10 Nov 2015 - The first version
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
  ! Justin,  16 Apr 2018 - replaced LAND_UPD ifdef with LAND_FIX from SurfParams, CSOIL_2 now called
  !                        from SurfParams, as well as others
  ! Sarith,  14 Aug 2018 - Added irrigation routines, considered experimental
  ! Sarith,  22 Apr 2020 - moved SUBROUTINE SRUNOFF here and modified to account for separate convective and 
  !                        large-scale throughfalls. FWETC and FWETL are now passed through the resource file.
  ! reichle, 27 Jan 2022 - moved "public" constants & subroutine echo_catch_constants() to catch_constants.f90
  
  use MAPL, ONLY:                                &
       MAPL_UNDEF

  USE MAPL_ConstantsMod, ONLY:                   &
       PIE               => MAPL_PI,             &  ! -
       TF                => MAPL_TICE,           &  ! K
       SHW               => MAPL_CAPWTR,         &  ! J/kg/K  spec heat of water
       SHI               => MAPL_CAPICE             ! J/kg/K  spec heat of ice

  USE CATCH_CONSTANTS,   ONLY:                   &
       N_SNOW            => CATCH_N_SNOW,        &
       N_GT              => CATCH_N_GT,          &
       DZTSURF           => CATCH_DZTSURF,       &
       DZGT              => CATCH_DZGT,          &
       PHIGT             => CATCH_PHIGT,         &
       FSN               => CATCH_FSN,           &       
       SHR               => CATCH_SHR,           &
       N_SM              => CATCH_N_ZONES,       &
       PEATCLSM_POROS_THRESHOLD,                 &
       PEATCLSM_ZBARMAX_4_SYSOIL
  
  USE SURFPARAMS,        ONLY:                   &
       LAND_FIX, FLWALPHA
  
  USE SIBALB_COEFF,      ONLY:                   &
       coeffsib

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: INTERC, SRUNOFF, RZDRAIN, BASE, PARTITION, RZEQUIL, gndtp0
  PUBLIC :: SIBALB, catch_calc_soil_moist, catch_calc_zbar, catch_calc_peatclsm_waterlevel
  PUBLIC :: catch_calc_subtile2tile
  PUBLIC :: gndtmp, catch_calc_tp, catch_calc_wtotl,  catch_calc_ght, catch_calc_FT
  PUBLIC :: dampen_tc_oscillations

  INTERFACE catch_calc_zbar
     MODULE PROCEDURE catch_calc_zbar_scalar
     MODULE PROCEDURE catch_calc_zbar_vector
  END INTERFACE catch_calc_zbar

  ! ---------------------------------------------------------------------------
  !
  !    ***** Do not define *public* Catchment model constants here. *****
  !
  ! Public Catchment model constants should be defined in ./catch_constants.f90
  !  (or in ../../GEOSsurface_GridComp/Shared/SurfParams.F90 if necessary).
  !
  !
  ! local constants:
  
  REAL,    PARAMETER :: CATCH_TIMFRL             = 1.0
  REAL,    PARAMETER :: CATCH_TIMFRC             = 0.333
  
  ! constants for "landscape" freeze/thaw (FT) state (see subroutine catch_calc_FT())
  
  REAL,    PARAMETER :: CATCH_FT_WEIGHT_TP1      = 0.5   !
  REAL,    PARAMETER :: CATCH_FT_THRESHOLD_TEFF  = TF    ! [Kelvin]
  REAL,    PARAMETER :: CATCH_FT_THRESHOLD_ASNOW = 0.2   !

  REAL,    PARAMETER :: ZERO                     = 0.
  REAL,    PARAMETER :: ONE                      = 1.
  
CONTAINS

!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN INTERC ]
!****

      SUBROUTINE INTERC (                             &
                         NCH, DTSTEP, FWETC, FWETL,   &
                         TRAINL, TRAINC, SMELT,       &  ! [kg m-2 s-1]
                         SATCAP,BUG,                  &
                         CAPAC,                       &
                         THRUL_VOL, THRUC_VOL         &  ! [kg m-2] !!!!
                        )
!****
!**** THIS ROUTINE USES THE PRECIPITATION FORCING TO DETERMINE
!**** CHANGES IN INTERCEPTION AND SOIL MOISTURE STORAGE.
!****
!****
!**** NOTE: Input precip & snow melt are in flux units [kg m-2 s-1]
!****       Output throughfall is in volume units      [kg m-2]
!****       (added comments for clarity but did not change units to preserve 0-diff, reichle, 6 Feb 2022)
        
      IMPLICIT NONE

!****
      INTEGER, INTENT(IN)                    :: NCH
      REAL,    INTENT(IN)                    :: DTSTEP, FWETC, FWETL
      REAL,    INTENT(IN),    DIMENSION(NCH) :: TRAINL, TRAINC, SMELT  ! [kg m-2 s-1]  (flux units)
      REAL,    INTENT(IN),    DIMENSION(NCH) :: SATCAP
      LOGICAL, INTENT(IN)                    :: BUG

      REAL,    INTENT(INOUT), DIMENSION(NCH) :: CAPAC

      REAL,    INTENT(OUT),   DIMENSION(NCH) :: THRUL_VOL, THRUC_VOL   ! [kg m-2]      ("volume" units) !!!
      
      ! --------------------------
      
      INTEGER :: CHNO
      REAL    :: WETINT, WATADD, CAVAIL, THRU1, THRU2, XTCORR, SMPERS, THRUL, THRUC

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

      XTCORR= (1.-CATCH_TIMFRL) *                                                    &
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

      THRUL=THRU1+THRU2

      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2

!****
!**** ---------------------------------------------------
!****
!**** STEP 2: MOIST CONVECTIVE PRECIPITATION.
!****
!**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
!**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.

      XTCORR= (1.-CATCH_TIMFRC) *                                                    &
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

      THRUC=THRU1+THRU2
      
      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2
!****
      IF (THRUL+THRUC .LT. -1.e-8) WRITE(*,*) 'THRU= ',                        &
          THRUL, THRUC, TRAINC(CHNO), TRAINL(CHNO), SMELT(CHNO)
      THRUL_VOL(CHNO)=AMAX1(0., THRUL)
      THRUC_VOL(CHNO)=AMAX1(0., THRUC)
      
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

      SUBROUTINE SRUNOFF (                                    &
           NCH, DTSTEP, UFW4RO, FWETC, FWETL, AR1, AR2, AR4,  &
           THRUL_VOL, THRUC_VOL,                              &  ! [kg m-2]      ("volume" units) !!!
           FRICE, TP1, SRFMX, BUG,                            &
           VGWMAX, RZEQ, POROS,                               &
           SRFEXC, RZEXC,                                     &
           RUNSRF,                                            &  ! [kg m-2 s-1]  (flux units)
           QINFIL                                             &  ! [kg m-2 s-1]  (flux units)
           )

!**** NOTE: Input throughfall is in volume units, as are calcs throughout this subroutine  [kg m-2]
!****       Input-output surface runoff and output infiltration are in flux units          [kg m-2 s-1]
!****         (added comments and clarified variable names, left throughfall units 
!****          unchanged to preserve 0-diff, reichle, 6 Feb 2022)
        
      IMPLICIT NONE

      INTEGER, INTENT(IN)                    :: NCH
      REAL,    INTENT(IN)                    :: DTSTEP, FWETC, FWETL
      LOGICAL, INTENT(IN)                    :: UFW4RO 
      REAL,    INTENT(IN),    DIMENSION(NCH) :: AR1, AR2, AR4, FRICE, TP1, SRFMX
      REAL,    INTENT(IN),    DIMENSION(NCH) :: VGWMAX, RZEQ, POROS
      REAL,    INTENT(IN),    DIMENSION(NCH) :: THRUL_VOL, THRUC_VOL              ! [kg m-2] 
      LOGICAL, INTENT(IN)                    :: BUG
      REAL,    INTENT(INOUT), DIMENSION(NCH) :: SRFEXC, RZEXC
      REAL,    INTENT(INOUT), DIMENSION(NCH) :: RUNSRF                            ! [kg m-2 s-1]
      REAL,    INTENT(OUT),   DIMENSION(NCH) :: QINFIL                            ! [kg m-2 s-1]

      ! ---------------------------
      
      INTEGER              :: N
      REAL                 :: deficit,srun0,frun,qin, qinfil_l, qinfil_c, qcapac, excess_infil
      REAL                 :: srunc, srunl, ptotal, excess, totcapac, watadd
      REAL, DIMENSION(NCH) :: THRUL, THRUC

      ! constants for PEATCLSM piecewise linear relationship between surface runoff and AR1
      
      REAL, PARAMETER      :: SRUN_AR1_MIN   =  0.5
      REAL, PARAMETER      :: SRUN_AR1_SLOPE = 10.
      
!**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      ! calculations throughout srunoff() are in "volume" units [kg m-2]
      
      THRUL  = THRUL_VOL        ! introduced *_VOL variables for clarity 
      THRUC  = THRUC_VOL   
      
      RUNSRF = RUNSRF * DTSTEP  ! convert input surface runoff from flux to "volume" units
      
      DO N=1,NCH

         if(.not.UFW4RO) then

            PTOTAL=THRUL(N) + THRUC(N)

            IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) THEN
               ! Non-peatland
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
               SRFEXC(N)=SRFEXC(N)+QIN
            ELSE
               ! Peatland
               ! MB: no Hortonian surface runoff
               !Note (rdk, from discussion w/MB; email 01/04/2021): 
               ! forcing all the rain to fall onto 
               ! the soil (1-AR1 fraction) rather than
               ! onto both the soil and the free water surface is a simple shortcut;
               ! the key function of this code is to retain all rainwater in the system
               ! and *only* to produce surface runoff when the 
               ! ground is already ridiculously wet.
               ! This prevents problems (e.g., numerical instabilities) found in 
               ! discharge calculations elsewhere in the code.

               srun0 = 0.
               ! handling numerical instability due to exceptional snow melt events at some pixels
               ! avoid AR1 to increase much higher than > 0.5 by enabling runoff
               !Added ramping to avoid potential oscillations (rdk, 09/18/20)
               IF (AR1(N)>SRUN_AR1_MIN) srun0=PTOTAL*amin1(1.,(ar1(n)-SRUN_AR1_MIN)*SRUN_AR1_SLOPE)

               ! MB: even no surface runoff when srfmx is exceeded (activating macro-pore flow)
               ! Rewrote code to determine excess over capacity all at once (rdk, 09/18/20)

               totcapac=(srfmx(n)-srfexc(n))+(vgwmax(n)-(rzeq(n)+rzexc(n)))
               watadd=ptotal-srun0
               if (watadd .gt. totcapac) then
                  excess=watadd-totcapac
                  srun0=srun0+excess
                  srfexc(n)=srfmx(n)
                  rzexc(n)=vgwmax(n)-rzeq(n)
               elseif(watadd .gt. srfmx(n)-srfexc(n)) then
                  excess=watadd-(srfmx(n)-srfexc(n))
                  srfexc(n)=srfmx(n)
                  rzexc(n)=rzexc(n)+excess
               else
                  srfexc(n)=srfexc(n)+watadd
               endif
               ! MB: check if VGWMAX is exceeded
               !IF(RZEQ(N) + RZEXC(N) .GT. (VGWMAX(N))) THEN
               !  srun0 = srun0 + RZEQ(N)+RZEXC(N)-VGWMAX(N)
               !  RZEXC(N)=VGWMAX(N)-RZEQ(N)
               !  ENDIF
               !(Commented out following lines to retain water balance -- rdk, 9/18/20)
               !if (srun0 .gt. ptotal) then
               !   srun0=ptotal
               !   endif
               RUNSRF(N)=RUNSRF(N)+srun0
               QIN=PTOTAL-srun0
               !SRFEXC(N)=amin1(SRFEXC(N)+QIN,srfmx(n))               
            ENDIF

         endif

         if(UFW4RO) then

            !**** Compute runoff from large-scale and convective storms separately:
            IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) THEN
                !non-peatland
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
               SRFEXC(N)=SRFEXC(N)+QIN

            else
               ! peatland
               ! MB: no Hortonian surface runoff
               !Note (rdk, from discussion w/MB; email 01/04/2021): 
               ! forcing all the rain to fall onto 
               ! the soil (1-AR1 fraction) rather than
               ! onto both the soil and the free water surface is a simple shortcut;
               ! the key function of this code is to retain all rainwater in the system
               ! and *only* to produce surface runoff when the 
               ! ground is already ridiculously wet.
               ! This prevents problems (e.g., numerical instabilities) found in 
               ! discharge calculations elsewhere in the code.

               srunl = 0.
               srunc = 0.
               ! handling numerical instability due to exceptional snow melt events at some pixels
               ! avoid AR1 to increase much higher than > 0.5 by enabling runoff
               IF (AR1(N)>SRUN_AR1_MIN) THEN
                  !Added ramping to avoid potential oscillations (rdk, 09/18/20)
                  srunl = THRUL(n)*amin1(1.,(ar1(n)-SRUN_AR1_MIN)*SRUN_AR1_SLOPE)
                  srunc = THRUC(n)*amin1(1.,(ar1(n)-SRUN_AR1_MIN)*SRUN_AR1_SLOPE)
               ENDIF
               PTOTAL = THRUL(N) + THRUC(N)
               SRUN0  = srunl + srunc
               ! MB: even no surface runoff when srfmx is exceeded (activating macro-pore flow)
               ! Rewrote code to determine excess over capacity all at once (rdk, 09/18/20)
               totcapac=(srfmx(n)-srfexc(n))+(vgwmax(n)-(rzeq(n)+rzexc(n)))
               watadd=ptotal-srun0
               if (watadd .gt. totcapac) then
                  excess=watadd-totcapac
                  srun0=srun0+excess
                  srfexc(n)=srfmx(n)
                  rzexc(n)=vgwmax(n)-rzeq(n)
               elseif(watadd .gt. srfmx(n)-srfexc(n)) then
                  excess=watadd-(srfmx(n)-srfexc(n))
                  srfexc(n)=srfmx(n)
                  rzexc(n)=rzexc(n)+excess
               else
                  srfexc(n)=srfexc(n)+watadd
               endif
               !if (ptotal-srun0 .gt. srfmx(n)-srfexc(n)) then
               !  excess=(ptotal-srun0)-(srfmx(n)-srfexc(n))
               !  rzexc(n)=rzexc(n) + excess
               !  ptotal=ptotal-excess
               !  endif                   
               ! MB: check if VGWMAX is exceeded
               !IF(RZEQ(N) + RZEXC(N) .GT. (VGWMAX(N))) THEN
               !  srun0 = srun0 + RZEQ(N)+RZEXC(N)-VGWMAX(N)
               !  RZEXC(N)=VGWMAX(N)-RZEQ(N)
               !  ENDIF
               RUNSRF(N)=RUNSRF(N)+srun0
               QIN=PTOTAL-srun0
               ! SRFEXC(N)=amin1(SRFEXC(N)+QIN,srfmx(n))  
            endif

         endif

         ! convert outputs to flux units [kg m-2 s-1]
         RUNSRF(N)=RUNSRF(N)/DTSTEP
         QINFIL(N)=QIN/DTSTEP

      END DO

      RETURN

      END SUBROUTINE SRUNOFF

!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****


      SUBROUTINE RZDRAIN (                                                   &
                          NCH, DTSTEP, VGWMAX, SATCAP, RZEQ, AR1, WPWET,     &
                          TSA1, TSA2, TSB1, TSB2, ATAU, BTAU, CDCR2, POROS,  &
                          BF1, BF2, ARS1, ARS2, ARS3, BUG,                   &
                          CAPAC, RZEXC, SRFEXC, CATDEF,                      &
                          RUNSRF                                             &  ! [kg m-2 s-1]
                          )

!-----------------------------------------------------------------
!        defines drainage timescales:
!             - tsc0, between srfex and rzex
!             - tsc2, between rzex and catdef
!        then defines correponding drainages
!        and updates the water contents
!-----------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN)                    :: NCH
      REAL,    INTENT(IN)                    :: DTSTEP
      REAL,    INTENT(IN),    DIMENSION(NCH) :: VGWMAX, SATCAP, RZEQ, AR1, wpwet
      REAL,    INTENT(IN),    DIMENSION(NCH) :: TSA1, TSA2, TSB1, TSB2, ATAU, BTAU, CDCR2, POROS
      REAL,    INTENT(IN),    DIMENSION(NCH) :: BF1,  BF2, ARS1, ARS2, ARS3
      LOGICAL, INTENT(IN)                    :: BUG

      REAL,    INTENT(INOUT), DIMENSION(NCH) :: RZEXC, SRFEXC, CATDEF, CAPAC   
      REAL,    INTENT(INOUT), DIMENSION(NCH) :: RUNSRF                          ! [kg m-2 s-1]

      ! --------------------

      INTEGER N
      REAL srflw,rzflw,FLOW,EXCESS,TSC0,tsc2,rzave,rz0,wanom,rztot,          &
           rzx,btaux,ax,bx,rzdif, rzavemin,ZBAR1,SYSOIL,RZFLW_CATDEF,        &
           EXCESS_CATDEF, CATDEF_PEAT_THRESHOLD, RZFLW_AR1, AR1eq



!**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO 100 N=1,NCH

!****   Compute equivalent of root zone excess in non-saturated area:
        rztot=rzeq(n)+rzexc(n)
        if(ar1(n).ne.1.) then
        !!! rzave=(rztot-ar1(n)*vgwmax(n))/(1.-ar1(n))
        !!! rzave=rzave*poros(n)/vgwmax(n)
            rzave=rztot*poros(n)/vgwmax(n)
          else
            rzave=poros(n)
          endif
  
! updated warning statement, reichle+koster, 12 Aug 2014
!
! Impose minimum of 1.e-4, rather than leaving positive values <1.e-4 unchanged.
! -reichle, 15 Jan 2016

        if (LAND_FIX) then
	   rzavemin = 1.e-4
	else
	   rzavemin = 0.
	end if

        if (rzave .le. rzavemin) then  ! JP: could put rzavemin in catch_constants
          rzave=1.e-4
          print*,'problem: rzave <= 1.e-4 in catchment',n
          end if

        btaux=btau(n)
        if (srfexc(n) .lt. 0.) btaux=btau(n)*(poros(n)/rzave)
        rz0=amax1(0.001,rzave-srfexc(n)/(1000.*(-btaux)))
        tsc0=atau(n)/(rz0**3.)

        tsc0=tsc0*3600.
        if(tsc0.lt.dtstep) tsc0=dtstep

! ---------------------------------------------------------------------

        SRFLW=SRFEXC(N)*DTSTEP/TSC0

        IF(SRFLW < 0.    ) SRFLW = FLWALPHA * SRFLW ! C05 change

!rr   following inserted by koster Sep 22, 2003
        rzdif=rzave/poros(n)-wpwet(n)
!**** No moisture transport up if rz at wilting; employ ramping.
        if(rzdif.le.0. .and. srflw.lt.0.)  srflw=0.
        if(rzdif.gt.0. .and. rzdif.lt.0.01                                     &
                   .and. srflw.lt.0.) srflw=srflw*(rzdif/0.01)
        RZEXC(N)=RZEXC(N)+SRFLW
        SRFEXC(N)=SRFEXC(N)-SRFLW

!**** Topography-dependent tsc2, between rzex and catdef

        rzx=rzexc(n)/vgwmax(n)

        if(rzx .gt. .01) then
            ax=tsa1(n)
            bx=tsb1(n)
          elseif(rzx .lt. -.01) then
            ax=tsa2(n)
            bx=tsb2(n)
          else
            ax=tsa2(n)+(rzx+.01)*(tsa1(n)-tsa2(n))/.02
            bx=tsb2(n)+(rzx+.01)*(tsb1(n)-tsb2(n))/.02
          endif

        tsc2=exp(ax+bx*catdef(n))
        rzflw=rzexc(n)*tsc2*dtstep/3600.

        IF (CATDEF(N)-RZFLW .GT. CDCR2(N)) then
          RZFLW=CATDEF(N)-CDCR2(N)
        end if

       IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) then
          ! mineral soil
          CATDEF(N)=CATDEF(N)-RZFLW
          RZEXC(N)=RZEXC(N)-RZFLW
       else
          !MB2021: use AR1eq, equilibrium assumption between water level in soil hummocks and surface water level in hollows
          AR1eq = (1.+ars1(n)*(catdef(n)))/(1.+ars2(n)*(catdef(n))+ars3(n)*(catdef(n))**2)
          ! PEAT
          ! MB: accounting for water ponding on AR1
          ! RZFLOW is partitioned into two flux components: (1) going in/out ponding water volume and (1) going in/out unsaturated soil storage
          ! Specific yield of ponding water surface fraction is 1.0
          ! calculate SYSOIL (see Dettmann and Bechtold, VZJ, 2016, for detailed theory)
          ! SYSOIL in CLSM can be derived from first derivative of 
          ! f_catdef(zbar) = ((zbar + bf2)^2 +1.0E-20)*bf1
          ! division by 1000 to convert from m to mm gives (Note: catdef in PEATCLSM remains
          ! the soil profile deficit, i.e. does not include the ponding water storage).
          ! SYSOIL = (2*bf1*zbar + 2*bf1*bf2)/1000
          ! Note: zbar defined here positive below ground.
          ! For the SYSOIL estimation zbar must be constrained to 0.0 to 0.45 m,
          ! to avoid extrapolation errors due to the non-optimal
          ! (linear) approximation with the bf1-bf2-CLSM function,
          ! theoretical SYSOIL curve levels off approximately at 0 m and 0.45 m.
          ZBAR1 = catch_calc_zbar( BF1(N), BF2(N), CATDEF(N) )  
          SYSOIL = (2.*bf1(n)*amin1(amax1(zbar1,0.),PEATCLSM_ZBARMAX_4_SYSOIL) + 2.*bf1(n)*bf2(n))/1000.
          SYSOIL = amin1(SYSOIL,poros(n))
          ! Calculate fraction of RZFLW removed/added to catdef
          RZFLW_CATDEF = (1.-AR1eq)*SYSOIL*RZFLW/(1.*AR1eq+SYSOIL*(1.-AR1eq))
          CATDEF(N)=CATDEF(N)-RZFLW_CATDEF
          ! MB: remove all RZFLW from RZEXC because the other part 
          ! flows into the surface water storage (microtopgraphy)
          RZEXC(N)=RZEXC(N)-RZFLW

       ENDIF
 

!****   REMOVE ANY EXCESS FROM MOISTURE RESERVOIRS:

        IF(CAPAC(N) .GT. SATCAP(N)) THEN
          RZEXC(N)=RZEXC(N)+CAPAC(N)-SATCAP(N)
          CAPAC(N)=SATCAP(N)
        ENDIF

        IF(RZEQ(N) + RZEXC(N) .GT. VGWMAX(N)) THEN
          EXCESS=RZEQ(N)+RZEXC(N)-VGWMAX(N)
          RZEXC(N)=VGWMAX(N)-RZEQ(N)

          IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) THEN
             CATDEF(N)=CATDEF(N)-EXCESS
          ELSE
             ! PEAT
             ! MB: like for RZFLW --> EXCESS_CATDEF is the fraction in/out of catdef
             EXCESS_CATDEF=(1.-AR1eq)*SYSOIL*EXCESS/(1.*AR1eq+SYSOIL*(1.-AR1eq))
             CATDEF(N)=CATDEF(N)-EXCESS_CATDEF
          ENDIF
       ENDIF

       IF (POROS(N) >= PEATCLSM_POROS_THRESHOLD) THEN
          ! MB: CATDEF Threshold at zbar=0
          ! water table not allowed to rise higher (numerically instable) 
          ! zbar<0 only occurred due to extreme infiltration rates
          ! (noticed this only snow melt events, very few locations and times)
          ! (--> NOTE: PEATCLSM has no Hortonian runoff for zbar > 0)            
          CATDEF_PEAT_THRESHOLD = ((BF2(N))**2.0)*BF1(N)
          IF(CATDEF(N) .LT. CATDEF_PEAT_THRESHOLD) THEN
             ! RUNSRF(N)=RUNSRF(N) + (CATDEF_PEAT_THRESHOLD - CATDEF(N))
             ! runoff from AR1 for zbar>0
             ! RZFLW_AR1 = RZFLW - RZFLW_CATDEF + (CATDEF_PEAT_THRESHOLD - CATDEF(N))
             ! AR1=0.5 at zbar=0
             ! SYsurface=0.5 at zbar=0
             ! RUNSRF(N) = RUNSRF(N) + amax1(0.0, RZFLW_AR1 - 0.5*1000.*ZBAR1)
             ! 
             ! revised (rdk, 1/04/2021): take excess water from both 
             ! soil and free standing water, the latter assumed to cover area AR1=0.5
             RUNSRF(N) = RUNSRF(N) + (CATDEF_PEAT_THRESHOLD-CATDEF(N) + 0.5*1000.*(-ZBAR1))/DTSTEP
             CATDEF(N)=CATDEF_PEAT_THRESHOLD
          ENDIF
       ENDIF

       IF(CATDEF(N) .LT. 0.) THEN
          ! bug fix: RUNSRF in flux units [kg m-2 s-1] for consistency with partition()
          ! reichle, 6 Feb 2022
          RUNSRF(N)=RUNSRF(N)-CATDEF(N)/DTSTEP
          CATDEF(N)=0.
       ENDIF

  100 ENDDO

      RETURN
      END SUBROUTINE RZDRAIN


!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****

      SUBROUTINE BASE (                                                        &
                       NCH,DTSTEP,BF1,BF2,BF3,CDCR1,FRICE,COND,GNU,            &
                       AR1, POROS, ARS1, ARS2, ARS3,                           &           
                       CATDEF,                                                 &
                       BFLOW                                                   &
                      )

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: BF1, BF2, BF3, CDCR1, FRICE, COND,   &
          GNU, AR1, POROS, ars1, ars2, ars3

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CATDEF

      REAL, INTENT(OUT), DIMENSION(NCH) :: BFLOW


      INTEGER N
      REAL ZBAR, ashift, CFRICE,Ksz_zero,m_Ivanov,v_slope,Ta,dztmp,SYSOIL,BFLOW_CATDEF,ICERAMP,AR1eq

      data ashift/0./

      DO N=1,NCH
         ZBAR = catch_calc_zbar( bf1(n), bf2(n), catdef(n) )  
         IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) THEN
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
         ELSE
            ! PEAT
            ! MB:  
            IF (FRICE(N) .GE. 0.9999) THEN
               CFRICE = 1.
            ELSE
               CFRICE = FRICE(N)
            ENDIF
            ! BFLOW in mm/s
            ! based on Ivanov 1981
            ! Ksz0  in  m/s
            ! m_Ivanov [-]  value depends on unit of Ksz0 and z
            ! v_slope in m^(-1)
            Ksz_zero=10.
            m_Ivanov=3.0
            v_slope = 1.5e-05
            ! Ta in m2/s, BFLOW in mm/s
            Ta = (Ksz_zero*(1.+100.*amax1(0.,ZBAR))**(1.-m_Ivanov))/(100.*(m_Ivanov-1.))
            BFLOW(N)=v_slope*Ta*1000.
            ! handling numerical instability due to extreme snow melt events on partly frozen ground
            ! --> allow BFLOW/DISCHARGE for zbar .LE. 0.05
            ICERAMP= AMAX1(0., AMIN1(1., ZBAR/0.05))
            ICERAMP= 1.-ICERAMP*CFRICE
            BFLOW(N)=ICERAMP*BFLOW(N)

            ! MB: Remove water from CATDEF and surface water storage
            IF (BFLOW(N) .NE. 0.0) THEN
               ! PEAT
               ! MB: accounting for water ponding on AR1
               ! same approach as for RZFLW (see subroutine RZDRAIN for
               ! comments)
               SYSOIL = (2.*bf1(N)*amin1(amax1(zbar,0.),PEATCLSM_ZBARMAX_4_SYSOIL) + 2.*bf1(N)*bf2(N))/1000.
               SYSOIL = amin1(SYSOIL,poros(n))
               !MB2021: use AR1eq, equilibrium assumption between water level in soil hummocks and surface water level in hollows
               AR1eq = (1.+ars1(n)*(catdef(n)))/(1.+ars2(n)*(catdef(n))+ars3(n)*(catdef(n))**2)
               BFLOW_CATDEF = (1.-AR1eq)*SYSOIL*BFLOW(N)/(1.*AR1eq+SYSOIL*(1.-AR1eq))
               CATDEF(N)=CATDEF(N)+BFLOW_CATDEF*dtstep
            ENDIF

         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE BASE


!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****

      SUBROUTINE PARTITION (                                                &
                            NCH,DTSTEP,DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,  &
                            PSIS,BEE,poros,WPWET,                           &
                            BF1, BF2,                                       &
                            ars1,ars2,ars3,ara1,ara2,ara3,ara4,             &
                            arw1,arw2,arw3,arw4,BUG,                        &
                            srfexc,catdef,                                  &
                            runsrf,                                         & ! [kg m-2 s-1]
                            AR1, AR2, AR4, srfmx, srfmn,                    &
                            SWSRF1,SWSRF2,SWSRF4,RZI                        &
                           )

      IMPLICIT NONE

! -------------------------------------------------------------------
      INTEGER, INTENT(IN)                    :: NCH

      REAL,    INTENT(IN)                    :: DTSTEP
      REAL,    INTENT(IN),    DIMENSION(NCH) :: DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,  &
                                                PSIS,BEE,poros,WPWET,BF1,BF2,        &
                                                ars1,ars2,ars3,ara1,ara2,ara3,ara4,  &
                                                arw1,arw2,arw3,arw4

      LOGICAL, INTENT(IN)                    :: BUG
! -------------------------------------------------------------------
      REAL,    INTENT(INOUT), DIMENSION(NCH) :: srfexc,catdef
      REAL,    INTENT(INOUT), DIMENSION(NCH) :: runsrf                                 ! [kg m-2 s-1]
! -------------------------------------------------------------------
      REAL,    INTENT(OUT),   DIMENSION(NCH) :: AR1, AR2, AR4, srfmx, srfmn,         &
                                                SWSRF1, SWSRF2, SWSRF4, RZI
! -------------------------------------------------------------------
      INTEGER :: N

      REAL :: cor, A150, W150, WMIN, AX, WMNEW, WRZ, TERM1, TERM2, TERM3,      &
              AREA0, AREA1, AREA2, AREA3, AREA4, ASCALE, WILT, D1, D2, CDI,    &
              DELTA1, DELTA2, DELTA4, MULTAR, CATDEFX, RZEQX, RZEQW, FACTOR,   &
              X0, RZEQY, CATDEFW, AR1W, ASUM, RZEQYI, RZEQWI, RZEQXI, AR20,    &
              ARG1, EXPARG1, ARG2, EXPARG2, ARG3, EXPARG3  !, surflay

      LOGICAL :: LSTRESS
      REAL    :: ZBAR, ARREST


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

        IF (POROS(N) >= PEATCLSM_POROS_THRESHOLD) THEN
           ! peat
           ! MB: AR4 (wilting fraction) for peatland depending on water table depth
           !ZBAR defined here positive below ground and in meter
           ZBAR = catch_calc_zbar( BF1(N), BF2(N), CATDEF(N) )  
           AR4(N)=amax1(0.,amin1(1.0,(ZBAR-0.30)/(1.0)))
           ARREST = 1.0 - AR1(N)
           AR4(N)=amin1(ARREST,AR4(N))
           AR2(N)=1.0-AR4(n)-AR1(N)
           ENDIF

        RZI(N)=RZEQYI

        SWSRF1(N)=1.
!mjs: changed .001 temporarily because of large bee.
        IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) THEN
        SWSRF2(N)=AMIN1(1., AMAX1(0.01, RZEQYI))
        SWSRF4(N)=AMIN1(1., AMAX1(0.01, WILT))

!**** EXTRAPOLATION OF THE SURFACE WETNESSES

! 1st step: surface wetness in the unstressed fraction without considering
!           the surface excess; we just assume an equilibrium profile from
!           the middle of the root zone to the surface.

        SWSRF2(N)=((SWSRF2(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))
        SWSRF4(N)=((SWSRF4(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))
        ELSE

             ! PEAT
             ! MB: for peatlands integrate across surface soil moisture distribution
             ! coefficients fitted for equilibrium conditions
             ! SWSRF2 and SWSRF4 as wetness (not moisture)
             ! MB: bug April 2018, AMIN1 function due to problems when spin up from
             ! scratch (i.e. dry soil at time=0)
             SWSRF2(N)=0.79437 - 0.99996*AMIN1(1.5,ZBAR) + 0.68801*(AMIN1(1.5,ZBAR))**2 + &
                     0.04186*(AMIN1(1.5,ZBAR))**3 - 0.15042*(AMIN1(1.5,ZBAR))**4
             SWSRF4(N)=SWSRF2(N)
        ENDIF

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
                          NCH,CATDEF,VGWMAX,CDCR1,CDCR2,WPWET,POROS,           &
                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,                  &
                          arw1,arw2,arw3,arw4,                                 &
                          RZEQ                                                 &
                         )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: CATDEF, VGWMAX, CDCR1, CDCR2,        &
                   WPWET, ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1,      &
                   arw2, arw3, arw4, POROS

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
! with depth given t1 (soil surface temperature; see below)
! gndtp0(): calculate heat fluxes between TSURF components and TP1; calculate derivatives
! gndtmp(): calculate soil heat content and temperature profiles
!             (layers 2-7: TP1, ..., TP6; GHTCNT1, ... GHTCNT6)
!
!            *****************************************
!      input:
!        t1      terrestrial (layer 1) surface temperature in deg C
!                (TC1, TC2, TC4 in Catchment and TG1, TG2, TG4 in CatchmentCN)
!        phi     porosity
!        zbar    mean depth to the water table.
!        thetaf  mean vadose zone soil moisture factor (0-1)
!      output:
!        ht      heat content in layers 2-7
!        tp      ground temperatures in layers 2-7
!        f21     heat flux between layer 2 and the terrestrial layer (1)
!                   "w"=wet         =saturated   <=> AR1
!                   "i"=intermediate=transpiring <=> AR2
!                   "d"=dry         =wilting     <=> AR4
!        df21    derivative of f21 with respect to temperature
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

      zb(1)=-DZTSURF
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

      a1=1.-phi
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

      denom=-(DZTSURF*0.5)-zc(1)
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
       NTILES,dzsf,vgwmax,cdcr1,cdcr2,psis,bee,poros,wpwet, &
       ars1,ars2,ars3,ara1,ara2, &
       ara3,ara4,arw1,arw2,arw3,arw4,bf1, bf2, &
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

    real,    dimension(NTILES), intent(in) :: dzsf,vgwmax,cdcr1,cdcr2
    real,    dimension(NTILES), intent(in) :: wpwet,poros,psis
    real,    dimension(NTILES), intent(in) :: bee,ars1,bf1, bf2
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
         cdcr1, cdcr2, wpwet, poros, &
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
         psis,bee,poros,wpwet,bf1, bf2, &
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
  
  ! Calculate zbar for Catchment[CN] model.
  !
  ! For mineral tiles, zbar is a fitted function that approximates the  
  !   water table depth EXCEPT for wet conditions.  For low values of 
  !   catdef, negative values on the order of -1 meter can be encountered,
  !   which would imply crazy water tables well above the surface.
  !   For this reason, zbar is not suitable for general estimation 
  !   of water table depth.
  !
  ! For PEATCLSM, zbar is a well-fitted function that describes the 
  !   water table depth for all wetness conditions.  At zbar=0, half of 
  !   the microtopography is flooeded.  Slightly negative values of
  !   up to -14 cm (-bf2) are theoretically possible but are not realized.
  !   Lesser negative values would represent slightly elevated water levels
  !   that imply more than half of the microtopography is flooded. 
  !
  ! Convention: zbar positive below ground (downward).
  !             
  ! This convention applies to water calculations, incl. subroutines RZDRAIN(),
  !   WUPDAT(), BASE(), PEATCLSM
  !
  ! WARNING:
  !   Opposite convention applies when zbar is used in ground heat
  !   diffusion model, incl. subroutines GNDTP0(), GNDTMP(), GNDTMP_CN().
  !
  ! - reichle, 29 Jan 2022
  ! - reichle,  3 Jun 2022 (updated documentation above)

  function catch_calc_zbar_scalar( bf1, bf2, catdef ) result(zbar)
    
    implicit none
    
    real,                       intent(in) :: bf1, bf2, catdef
    real                                   :: zbar

    zbar = SQRT(1.e-20 + catdef/bf1) - bf2
    
  end function catch_calc_zbar_scalar
  
  ! --------------------------
  
  function catch_calc_zbar_vector( bf1, bf2, catdef ) result(zbar)
    
    ! vector version of catch_calc_zbar 

    implicit none
    
    real, dimension(:),         intent(in) :: bf1, bf2, catdef
    real, dimension(size(bf1))             :: zbar
    
    zbar = SQRT(1.e-20 + catdef/bf1) - bf2
    
  end function catch_calc_zbar_vector
  
  ! *******************************************************************

  function catch_calc_peatclsm_waterlevel( bf1, bf2, cdcr2, poros, wpwet, catdef ) result(waterlevel)
    
    ! calculate water level (a.k.a. water table depth) for PEATCLSM only [m]
    !
    ! Convention: water leve positive above ground (opposite of zbar convention!)

    implicit none
    
    real, dimension(:),         intent(in) :: bf1, bf2, cdcr2, poros, wpwet, catdef
    real, dimension(size(bf1))             :: waterlevel
    
    WHERE (POROS >= PEATCLSM_POROS_THRESHOLD)
       
       ! note change of sign from zbar

       waterlevel = -1.*MIN( catch_calc_zbar(BF1,BF2,CATDEF), CDCR2/(1.-WPWET)/POROS/1000. )
       
    ELSEWHERE
       
       waterlevel = MAPL_UNDEF

    ENDWHERE       
       
  end function catch_calc_peatclsm_waterlevel
  
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

  subroutine gndtmp(phi,zbar,ht,xfice,tp,FICE,dts,thetaf,fh21)

    ! Added functionality for "short" and "full" versions of gndtmp():
    !   The "short" version of gndtmp() was originally introduced by Greg Walker 
    !   for CatchmentCN by duplicating code in catchmentCN.F90.
    !   It is now part of gndtmp() when called without the optional arguments.
    !   - reichle+jkolassa, 10 Feb 2022
    
! using a diffusion equation this code generates ground temperatures
! with depth given t1 (soil surface temperature; see below)
! gndtp0(): calculate heat fluxes between TSURF components and TP1; calculate derivatives
! gndtmp(): calculate soil heat content and temperature profiles
!             (layers 2-7: TP1, ..., TP6; GHTCNT1, ... GHTCNT6)
!
! t1 = terrestrial (layer 1) surface temperature in deg C
!                (TC1, TC2, TC4 in Catchment and TG1, TG2, TG4 in CatchmentCN)
!
!            *****************************************
!      input:
!        dts     timestep in seconds
!        phi     porosity
!        zbar    mean depth to the water table.
!        thetaf  mean vadose zone soil moisture factor (0-1)
!        fh21    ground heat flux (between layer 1 and layer 2)
!      input/output:
!        ht      heat content in layers 2-7
!      output:
!        tp      ground temperatures in layers 2-7
!        xfice   a total soil column ice factor (0-1)  [in saturated portion of layers 2-7 only???]
!        FICE    soil ice fraction in layers 2-7
!             ***********************************

      REAL, INTENT(IN)                            :: phi, ZBAR
      REAL, INTENT(IN),                  OPTIONAL :: DTS, THETAF, FH21

      REAL, INTENT(INOUT), DIMENSION(*)           :: HT                        ! dimension(N_GT)

      REAL, INTENT(OUT) :: XFICE
      REAL, INTENT(OUT),   DIMENSION(*)           :: TP, FICE                  ! dimension(N_GT)

      ! ----------------------------------
      
      INTEGER                    :: L, LSTART, K
      REAL, DIMENSION(N_GT) :: ZC, SHC, XKLH
      REAL, DIMENSION(N_GT+1) :: FH, ZB
      REAL                       :: SHW0, SHI0, SHR0, WS, XW, A1, TK1, A2, TK2, TK3, TKSAT
      REAL                       :: XWI, XD1, XD2, XKTH, TKDRY
      LOGICAL                    :: full_version

      !data dz/0.0988,0.1952,0.3859,0.7626,1.5071,10.0/
      !DATA PHI/0.45/, FSN/3.34e+8/, SHR/2.4E6/

! initialize parameters

      shw0=SHW*1000. ! PER M RATHER THAN PER KG
      shi0=SHI*1000. ! PER M RATHER THAN PER KG
      shr0=SHR*1000. ! PER M RATHER THAN PER KG [kg of water equivalent density]

! ---------------------------------

      ! process optional arguments (all or none must be present)
      
      if     (       present(DTS)  .and.        present(THETAF)  .and.        present(FH21) ) then
         full_version = .true.
      elseif ((.not. present(DTS)) .and. (.not. present(THETAF)) .and. (.not. present(FH21))) then
         full_version = .false.
      else
         ! should be replaced with better error handling...
         write (*,*) 'ERROR with optional input arguments for gndtmp()... stopping'
         stop
      endif

!----------------------------------
! initialize fice in ALL components
! reichle, July 8, 2002

      do l=1,N_GT
        fice(l)=0.0
        enddo
!----------------------------------

! calculate the boundaries, based on the layer thicknesses(DZGT)

      zb(1)=-DZTSURF  ! Bottom of surface layer, which is handled outside
                      ! this routine.

      do l=1,N_GT
        zb(l+1)=zb(l)-DZGT(l)
        shc(l)=shr0*(1.-phi)*DZGT(l)
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

      if (full_version) then

         ! full version of gndtmp: additional updates of tp and fice
         
         do l=1,N_GT
           zc(l)=0.5*(zb(l)+zb(l+1))
           enddo
    
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

         a1=1.-phi             ! ROCK FRACTION
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

      end if ! full version of gndtmp()

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

      IF (phi < PEATCLSM_POROS_THRESHOLD) THEN
         xfice=xfice/((N_GT+1)-lstart)
      ELSE
         !PEAT
         !MB: only first layer for total runoff reduction
        xfice=AMIN1(1.0,fice(1))
      ENDIF

      Return

      end subroutine gndtmp

 ! *******************************************************************

  

  





    

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

          shc(k) = shr0*(1.-phi)*DZGT(k)

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

    shc = shr0*(1.-phi)*DZGT

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

  ! *******************************************************************
  


END MODULE lsm_routines
