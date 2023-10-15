!**** "Cleaned-up" version: July 2002 (rk)
!rr   ---------------------------
!rr
!rr   this version:!rr   catchment.f as from /land/koster/src on May 29, 2003
!rr   snow_may2003.f as from Stephen Dery by email on May 29, 2003
!rr   
!rr   several snow code bug fixes
!rr   snow water equivalent in [mm]=[kg/m2] *throughout*
!rr   order of arguments changed in snowrt
!rr   
!rr   additional modifications:
!rr   removed "mxnch" AND "mxchp"
!rr   commented out a bunch of write statements
!rr   - reichle, 29 May 03
!rr   ---------------------------

! koster+reichle, 13 Aug 2008: added optional output of tcX and qaX as 
!                               computed before AGCM adjustments
! reichle, 28 Oct 2010 - moved DZ, SHR, PHI, FSN, FWETL, FWETC to module catch_constants
!                      - renamed N_gndtmp -> N_gt
! reichle, 28 Oct 2010 - moved SURFLAY to GEOS_CatchGridComp, pass into catchment()
! reichle, 23 Nov 2010 - replaced PHIGT with POROS(N), ALHMGT with ALHM
! reichle, 30 Nov 2010 - zero-diff revisions and clean-up for off-line (land-only) MERRA
!                        replay capability
!                         - restored PHIGT, ALHMGT 
!                         - moved MIN_SNOW_MASS->MINSWE, DZ1MAX, and SATCAPFR to 
!                            catch_constants
!                         - moved "small" back from catch_constants() into snowrt()
! reichle,  2 Apr 2012 - moved Catchment diagnostics routines to here from 
!                        file "catch_diagn_routines.F90" of "lana" directory (LDAS), 
!                        renamed for consistency, and revised for efficiency
! reichle, 12 Aug 2014 - added engineering fix for surface energy balance oscillations 
!                        (so far only used for off-line land modeling)
!                      - moved "catch_calc_tpsnow()" and "catch_calc_asnow()" to StieglitzSnow.F90
!                      - use "get_tf0d()" from StieglitzSnow.F90 [later to be unified with
!                         "StieglitzSnow_calc_tpsnow()"]
!                      - moved constants that are only needed in subroutine catchment() 
!                         from catch_constants.f90 to catchment.F90 where they are now private 
!                      - renamed remaining public constants by prefacing with "catch_*"
!                      - added "catch_echo_constants()" (moved and renamed from catch_constants.f90)
!                      - minor clean-up
! reichle, 20 Oct 2014 - added diagnostic subroutines catch_calc_ght(), catch_calc_FT(), 
!                         and catch_calc_tsurf_excl_snow()
! reichle, 24 Nov 2015 - changed CSOIL_2 back to pre-MERRA value (70,000 J/K)
!                      - use engineering fix for all veg types (dampen oscillations in off-line mode)
!                      - "zbar" bug fix
! Sarith, 10 Nov 2015  - moved   RZDRAIN, INTERC, BASE, PARTITION, RZEQUIL, gndtp0
!                        SIBALB, catch_calc_soil_moist, catch_calc_subtile2tile
!                        gndtmp, catch_calc_tp,  catch_calc_ght, catch_calc_FT, 
!                        catch_calc_wtotl, dampen_tc_oscillations and catch_echo_constants to
!                        ../Shared/lsm_routines.F90 during adding GEOScatchCN_GridComp and 
!                        reorganizing GEOSland_gridComp. 
!                      - moved DZTC, FWETL, FWETC, DZGT, PHIGT, ALHMGT, FSN, CATCH_FT_WEIGHT_TP1,
!                        CATCH_FT_THRESHOLD_TEFF, CATCH_FT_THRESHOLD_ASNOW, ZERO, and ONE to 
!                        ../Shared/lsm_routines.F90 
!                      - moved SHR, EPSILON, SCONST, CSOIL_1,  CSOIL_2, N_sm, and SATCAPFR to
!                        ../Shared/catch_constants.f90
! Justin, 19 Apr 2018  - removed LAND_UPD ifdefs, use SurfParams
! Justin, 11 Dec 2018  - put in ASNOW fix affecting AGCM only
! Sarith, 20 Apr 2020  - introducing USE_FWET_FOR_RUNOFF and passing FWETL and FWETC via GEOS_SurfaceGridComp.rc
! Reichle, 14 Jan 2022 - removed redundant qa constraint; removed commented-out #ifdef LAND_UPD directives

      MODULE CATCHMENT_MODEL

      USE MAPL_BaseMod,      ONLY:               &
           NTYPS             => MAPL_NumVegTypes, &
           MAPL_Land,                             &
           MAPL_UNDEF
      
      USE MAPL_ConstantsMod, ONLY:          &
           PIE               => MAPL_PI,     &  ! -                       
           ALHE              => MAPL_ALHL,   &  ! J/kg  @15C              
           ALHM              => MAPL_ALHF,   &  ! J/kg                    
           ALHS              => MAPL_ALHS,   &  ! J/kg                    
           TF                => MAPL_TICE,   &  ! K                       
           RGAS              => MAPL_RGAS,   &  ! J/(kg K)                
           SHW               => MAPL_CAPWTR, &  ! J/kg/K  spec heat of wat
           SHI               => MAPL_CAPICE, &  ! J/kg/K  spec heat of ice
           EPSILON           => MAPL_EPSILON
      
      USE CATCH_CONSTANTS,   ONLY:                   &
           N_SNOW            => CATCH_N_SNOW,        &
           N_GT              => CATCH_N_GT,          &
           CATCH_SNOW_RHOFS,                         &
           CATCH_SNOW_MAXDEPTH,                      &
           CATCH_SNOW_DZPARAM,                       &
           CSOIL_1           => CATCH_CSOIL_1,       &
           N_sm              => CATCH_N_ZONES,       &
           SATCAPFR          => CATCH_SATCAPFR,      &
           PHIGT             => CATCH_PHIGT,         &
           DZTSURF           => CATCH_DZTSURF,       &
           PEATCLSM_POROS_THRESHOLD,                 &
           PEATCLSM_ZBARMAX_4_SYSOIL

      USE SURFPARAMS,       ONLY:                    &
	   LAND_FIX, ASTRFR, STEXP, RSWILT,          &
	   FLWALPHA, CSOIL_2 

      USE lsm_routines, only :                       &
           INTERC, RZDRAIN, BASE, PARTITION, RZEQUIL,&
           gndtp0, gndtmp,                           &
           catch_calc_soil_moist, catch_calc_zbar,   &
           catch_calc_wtotl, dampen_tc_oscillations, &
           SRUNOFF 
      
      USE SIBALB_COEFF,  ONLY: coeffsib
      
      USE STIEGLITZSNOW, ONLY:                                                               &
           StieglitzSnow_snowrt,                     & 
           StieglitzSnow_calc_asnow,                 &
           StieglitzSnow_calc_tpsnow,                &
           N_constit
      
      IMPLICIT NONE

      private

      public :: catchment
      public :: catch_calc_tsurf
      public :: catch_calc_tsurf_excl_snow
      public :: catch_calc_etotl

      ! -----------------------------------------------------------------------------
      
      ! moved all "private" constants to here from "catch_constants.f90"
      ! - reichle, 14 Aug 2014
            
      REAL,    PARAMETER :: ZERO     = 0.
      REAL,    PARAMETER :: ONE      = 1.
       
      CONTAINS

      SUBROUTINE CATCHMENT (                                                   &
                     NCH, LONS, LATS, DTSTEP, UFW4RO, FWETC, FWETL,            &
                     cat_id,ITYP,DZSF,TRAINC,TRAINL, TSNOW, TICE, TFRZR, UM,   &
                     ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,         &
                     ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,         &
                     ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,         &
                     ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,         &
                     TM, QM, ra1, ra2, ra4, raS, SUNANG, PARDIR, PARDIF,       &
                     SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,  Z2,        &  ! HLWDWN = *absorbed* longwave only (excl reflected)
                     SQSCAT, RSOIL1, RSOIL2,   RDC,                            &
                     QSAT1, DQS1, ALW1, BLW1,  QSAT2, DQS2, ALW2, BLW2,        &
                     QSAT4, DQS4, ALW4, BLW4,  QSATS, DQSS, ALWS, BLWS,        &
                     BF1, BF2, BF3,VGWMAX,                                     &
                     CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,          &
                     ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,   &
                     tsa1,tsa2,tsb1,tsb2,atau,btau,BUG,                        &
                     TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,                      &
                     CATDEF, RZEXC, srfexc, GHTCNT, TSURF,                     &
                     WESNN, HTSNNN, SNDZN,                                     &
                     EVAP, SHFLUX, RUNOFF,                                     &
                     EINT, ESOI, EVEG, ESNO,  BFLOW,RUNSRF,SMELT,              &
                     HLWUP,SWLAND,HLATN,QINFIL,AR1, AR2, RZEQ,                 &  ! HLWUP = *emitted* longwave only (excl reflected)
                     GHFLUX, GHFLUXSNO, GHTSKIN, TPSN1, ASNOW0,                &
                     TP1, TP2, TP3, TP4, TP5, TP6,                             &
                     sfmc, rzmc, prmc, entot, wtot, WCHANGE, ECHANGE, HSNACC,  &
                     EVACC, SHACC,                                             &
                     SH_SNOW, AVET_SNOW, WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL,   &
                     LH_SNOW, LWUP_SNOW, LWDOWN_SNOW, NETSW_SNOW,              &
                     TCSORIG, TPSN1IN, TPSN1OUT, FSW_CHANGE, FICESOUT,         &
                     lonbeg,lonend,latbeg,latend,                              &
                     TC1_0, TC2_0, TC4_0, QA1_0, QA2_0, QA4_0, EACC_0,         &  ! OPTIONAL
                     RCONSTIT, RMELT, TOTDEPOS,  LHACC )                          ! OPTIONAL

      IMPLICIT NONE

! -----------------------------------------------------------
!     INPUTS

      INTEGER, INTENT(IN)                 :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP, cat_id

      REAL,    INTENT(IN)                 :: DTSTEP, FWETC, FWETL
      LOGICAL, INTENT(IN)                 :: UFW4RO

      REAL,    INTENT(IN), DIMENSION(NCH) :: DZSF, TRAINC, TRAINL,             &
                     TSNOW, TICE, TFRZR,  UM,                                  &
                     ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,         &
                     ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,         &
                     ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,         &
                     ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,         &
                     TM, QM, ra1, ra2, ra4, raS, SUNANG, PARDIR, PARDIF,       &
                     SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,  Z2,        &
                     SQSCAT, RSOIL1, RSOIL2,   RDC,                            &
                     QSAT1, DQS1, ALW1, BLW1,  QSAT2, DQS2, ALW2, BLW2,        &
                     QSAT4, DQS4, ALW4, BLW4,  QSATS, DQSS, ALWS, BLWS,        &
                     BF1, BF2, BF3,VGWMAX,                                     &
                     CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,          &
                     ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,   &
                     tsa1,tsa2,tsb1,tsb2,atau,btau
      
      REAL,    INTENT(IN), DIMENSION(NCH)                      :: LONS, LATS

      REAL,    INTENT(IN), DIMENSION(NCH, N_Constit), OPTIONAL :: TOTDEPOS

      LOGICAL, INTENT(IN) :: BUG

      REAL,    INTENT(IN) :: lonbeg, lonend, latbeg, latend

! -----------------------------------------------------------
!     PROGNOSTIC VARIABLES

      REAL, INTENT(INOUT), DIMENSION(NCH) ::                                   &
                     TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,                      &
                     CATDEF, RZEXC, SRFEXC
 
      REAL, INTENT(INOUT), DIMENSION(N_GT, NCH) ::  GHTCNT
 
      REAL, INTENT(INOUT), DIMENSION(N_SNOW, NCH) :: WESNN, HTSNNN, SNDZN

      REAL, INTENT(INOUT), DIMENSION(NCH, N_SNOW, N_Constit),                  &
               OPTIONAL :: RCONSTIT

! -----------------------------------------------------------
!     DIAGNOSTIC OUTPUT VARIABLES

      REAL, INTENT(OUT), DIMENSION(NCH) ::      EVAP, SHFLUX, RUNOFF,          &
                     EINT, ESOI, EVEG, ESNO, BFLOW,RUNSRF,SMELT,               &
                     HLWUP,SWLAND,HLATN,QINFIL,AR1, AR2, RZEQ,                 &
                     GHFLUX, TPSN1, ASNOW0, TP1, TP2, TP3, TP4, TP5, TP6,      &
                     sfmc, rzmc, prmc, entot, wtot, tsurf, WCHANGE, ECHANGE,   &
                     HSNACC, EVACC, SHACC
      REAL, INTENT(OUT), DIMENSION(NCH) :: GHFLUXSNO, GHTSKIN

      REAL, INTENT(OUT), DIMENSION(NCH) :: SH_SNOW, AVET_SNOW,         &
                     WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL
      REAL, INTENT(OUT), DIMENSION(NCH) :: LH_SNOW, LWUP_SNOW,         &
                     LWDOWN_SNOW, NETSW_SNOW
      REAL, INTENT(OUT), DIMENSION(NCH) :: TCSORIG, TPSN1IN, TPSN1OUT, &
                     FSW_CHANGE

      REAL, INTENT(OUT), DIMENSION(N_SNOW, NCH)   :: FICESOUT
      
      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: LHACC

      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: TC1_0,TC2_0,TC4_0
      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: QA1_0,QA2_0,QA4_0 	
      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: EACC_0 

      REAL, INTENT(OUT), DIMENSION(NCH, N_Constit), OPTIONAL :: RMELT	

! -----------------------------------------------------------
!     LOCAL VARIABLES

      INTEGER I,K,N,LAYER

      REAL, DIMENSION(NCH) ::  CSOIL, ASNOW, traincx, trainlx, EMAXRT,         &
            RC, SATCAP, SNWFRC, POTFRC,  ESNFRC, EVSNOW, SHFLUXS, HLWUPS,      &
            HFTDS1, HFTDS2, HFTDS4, DHFT1, DHFT2, DHFT4, TPSNB,                &
            QSATTC, DQSDTC, SWSRF1, SWSRF2, SWSRF4, AR4, RX11, RX21, RX12,     &
            RX14, RX24, RX22, EIRFRC, FCAN, THRUL_VOL, THRUC_VOL,              &
            RZEQOL, frice, srfmx,                                              &
            srfmn, RCST, EVAPFR, RCUN, PAR, PDIR, RDCX, EVAP1, EVAP2,          &
            EVAP4, SHFLUX1, SHFLUX2, SHFLUX4, HLWUP1, HLWUP2, HLWUP4,          &
            GHFLUX1, GHFLUX2, GHFLUX4, RZI, TC1SF, TC2SF, TC4SF, ar1old,       &
            ar2old, ar4old, GHFLUXS, DEDQA1X, DEDTC1X,                         &
            DHSDQA1X, DHSDTC1X, DEDQA2X, DEDTC2X, DHSDQA2X, DHSDTC2X,          &
            DEDQA4X, DEDTC4X, DHSDQA4X, DHSDTC4X, werror, sfmcun, rzmcun,      &
            prmcun,WTOT_ORIG,ENTOT_ORIG,HSNACC1,HSNACC2,HSNACC4,               &
            TC1_00, TC2_00, TC4_00, EACC_00,                                   &
              qa1_orig,qa2_orig,qa4_orig,tc1_orig,tc2_orig,tc4_orig,           &
              tcs_orig


      REAL, DIMENSION(N_GT) :: HT, TP, soilice

      REAL, DIMENSION(N_SNOW) :: TPSN, WESN, HTSNN, SNDZ, fices,               &
             wesnperc,wesndens,wesnrepar,excs,drho0,tksno, tmpvec_Nsnow

      REAL, DIMENSION(N_SNOW, N_Constit) :: RCONSTIT1

      REAL, DIMENSION(N_Constit)         :: RMELT1, TOTDEP1

      REAL, DIMENSION(N_SM) :: T1, AREA, tkgnd, fhgnd

      REAL :: EPFRC1, EPFRC2, EPFRC4, SUMEP, SUME, TC1SN, TC2SN, TC4SN,        &
              DTC1SN,DTC2SN,DTC4SN, RTBS1, RTBS2, RTBS4, ZBAR, THETAF,         &
              XFICE, FH21, FH21W, FH21I, FH21D, DFH21W, DFH21I, DFH21D,        &
              EVSN, SHFLS, HUPS, HCORR, SWNET0, HLWDWN0, TMPSNW, HLWTC,        &
              DHLWTC, HSTURB, DHSDEA, DHSDTC, ESATTC, ETURB, DEDEA, DEDTC,     &
              SNOWF, TS, fh31w, fh31i, fh31d, pr, ea, desdtc, areasc,          &
              pre, dummy1, dummy2, dummy3, areasc0, EDIF, EINTX,               &
              SCLAI, tsn1, tsn2, tsn3, hold, hnew, emaxrz, dedtc0,             &
              dhsdtc0, alhfsn, ADJ, raddn, zc1, tsnowsrf, dum, tsoil,          &
              QA1X, QA2X, QA4X, TC1X, TC2X, TC4X, TCSX,                        &
              EVAPX1,EVAPX2,EVAPX4,SHFLUXX1,SHFLUXX2,SHFLUXX4,EVEGFRC,         &
              EVAPXS,SHFLUXXS,phi,rho_fs,sumdepth,                             &   
              sndzsc, wesnprec, sndzprec,  sndz1perc,                          &   
              mltwtr, wesnbot, dtss



      LOGICAL :: ldum

      REAL    :: dtc1, dtc2, dtc4

      integer  numout
      integer  n_out
      integer  n_outs(20)

      numout =  0

! choose output point by lon and lat Input lons and lats are in radians
! EXAMPLE:
! 0.120643381534E+03  0.650233927779E+02       in degrees
! 2.1056              1.13487                  in radians (= deg * pi/180)
!
! User beware: if the range of lats and lons spans processors the output will be a mess!
!
      if( lonbeg.ne.MAPL_UNDEF .and. &
          lonend.ne.MAPL_UNDEF .and. &
          latbeg.ne.MAPL_UNDEF .and. &
          latend.ne.MAPL_UNDEF       ) then
      do i = 1,nch
       if ((lons(i).ge.lonbeg*PIE/180.).and.(lons(i).le.lonend*PIE/180.).and.(lats(i).ge.latbeg*PIE/180.).and.(lats(i).le.latend*PIE/180.)) then
        numout = numout + 1
        n_outs(numout) = i
        write (*,*) 'found point ',n_out,' ',lons(i)*180./pie,' ',lats(i)*180./pie
       endif
      enddo
      endif

      if(numout.ne.0) then
       do i = 1,numout
         n_out = n_outs(i)

         write (*,*) 'INPUT catchment arguments for n_out: ',n_out

         write (*,*) 'qsats = ', qsats(n_out)  
         write (*,*)  'TC1 = ', tc1(n_out)    
         write (*,*)  'TC2 = ', tc2(n_out)  
         write (*,*)  'TC4 = ', tc4(n_out)
         write (*,*)  'TPSN1 = ', tpsn1(n_out)
         write (*,*)  'TSURF = ',    TSURF(n_out)  
         write (*,*)  'WESNN(1), = ',    WESNN(1,n_out)    
         write (*,*)  'HTSNNN(1), = ',    HTSNNN(1,n_out)    
         write (*,*)  'SNDZN(1), = ',  SNDZN(1,n_out)  
         
         write (*,*) NCH  
         write (*,*) DTSTEP  
         write (*,*) UFW4RO
         write (*,*) ITYP(n_out)  
         write (*,*) TRAINC(n_out)    
         write (*,*) TRAINL(n_out)    
         write (*,*) TSNOW(n_out)    
         write (*,*) TICE(n_out)    
         write (*,*) TFRZR(n_out)    
         write (*,*) UM(n_out)  
         write (*,*) ETURB1(n_out)    
         write (*,*) DEDQA1(n_out)    
         write (*,*) DEDTC1(n_out)    
         write (*,*)  HSTURB1(n_out)    
         write (*,*) DHSDQA1(n_out)    
         write (*,*)  DHSDTC1(n_out)  
         write (*,*)  ETURB2(n_out)    
         write (*,*)  DEDQA2(n_out)    
         write (*,*)  DEDTC2(n_out)    
         write (*,*)  HSTURB2(n_out)    
         write (*,*)  DHSDQA2(n_out)    
         write (*,*)  DHSDTC2(n_out)  
         write (*,*)  ETURB4(n_out)    
         write (*,*)  DEDQA4(n_out)    
         write (*,*)  DEDTC4(n_out)    
         write (*,*)  HSTURB4(n_out)    
         write (*,*)  DHSDQA4(n_out)    
         write (*,*)  DHSDTC4(n_out)  
         write (*,*)  ETURBS(n_out)    
         write (*,*)  DEDQAS(n_out)    
         write (*,*)  DEDTCS(n_out)    
         write (*,*)  HSTURBS(n_out)    
         write (*,*)  DHSDQAS(n_out)    
         write (*,*)  DHSDTCS(n_out)  
         write (*,*)  TM(n_out)    
         write (*,*)  QM(n_out)    
         write (*,*)  ra1(n_out)    
         write (*,*)  ra2(n_out)    
         write (*,*)  ra4(n_out)    
         write (*,*)  raS(n_out)    
         write (*,*)  SUNANG(n_out)    
         write (*,*)  PARDIR(n_out)    
         write (*,*)  PARDIF(n_out)  
         write (*,*)  SWNETF(n_out)    
         write (*,*)  SWNETS(n_out)    
         write (*,*)   HLWDWN(n_out)    
         write (*,*)  PSUR(n_out)    
         write (*,*)   ZLAI(n_out)    
         write (*,*)    GREEN(n_out)    
         write (*,*)   Z2(n_out)  
         write (*,*)  SQSCAT(n_out)    
         write (*,*)  RSOIL1(n_out)    
         write (*,*)  RSOIL2(n_out)    
         write (*,*)    RDC(n_out)    
!         write (*,*)     U2FAC(n_out)  
         write (*,*)  QSAT1(n_out)    
         write (*,*)  DQS1(n_out)    
         write (*,*)  ALW1(n_out)    
         write (*,*)  BLW1(n_out)  
         write (*,*)  QSAT2(n_out)    
         write (*,*)  DQS2(n_out)    
         write (*,*)  ALW2(n_out)    
         write (*,*)  BLW2(n_out)  
         write (*,*)  QSAT4(n_out)    
         write (*,*)  DQS4(n_out)    
         write (*,*)  ALW4(n_out)    
         write (*,*)  BLW4(n_out)  
         write (*,*)  QSATS(n_out)    
         write (*,*)  DQSS(n_out)    
         write (*,*)  ALWS(n_out)    
         write (*,*)  BLWS(n_out)  
         write (*,*)  BF1(n_out)    
         write (*,*)  BF2(n_out)    
         write (*,*)  BF3(n_out)    
         write (*,*)  VGWMAX(n_out)  
         write (*,*)  CDCR1(n_out)    
         write (*,*)  CDCR2(n_out)    
         write (*,*)  psis(n_out)    
         write (*,*)  bee(n_out)    
         write (*,*)  poros(n_out)    
         write (*,*)  wpwet(n_out)    
         write (*,*)  cond(n_out)    
         write (*,*)  'gnu=',gnu(n_out)  
         write (*,*)  'ars1=',ARS1(n_out)    
         write (*,*)  'ars2=',ARS2(n_out)    
         write (*,*)  'ars3=',ARS3(n_out)    
         write (*,*)  'ara1=',ARA1(n_out)    
         write (*,*)  'ara2=',ARA2(n_out)    
         write (*,*)  'ara3=',ARA3(n_out)    
         write (*,*)  'ara4=',ARA4(n_out)  
         write (*,*)  'arw1=',ARW1(n_out)    
         write (*,*)  'arw2=',ARW2(n_out)    
         write (*,*)  'arw3=',ARW3(n_out)    
         write (*,*)  'arw4=',ARW4(n_out)  
         write (*,*)  'tsa1=',tsa1(n_out)    
         write (*,*)  'tsa2=',tsa2(n_out)    
         write (*,*)  'tsb1=',tsb1(n_out)    
         write (*,*)  'tsb2=',tsb2(n_out)    
         write (*,*)  'atau=',atau(n_out)    
         write (*,*)  'btau=',btau(n_out)    
         write (*,*)  'BUG=',BUG
         write (*,*)  TC1(n_out)    
         write (*,*)  TC2(n_out)    
         write (*,*)  TC4(n_out)    
         write (*,*)  QA1(n_out)    
         write (*,*)  QA2(n_out)    
         write (*,*)  QA4(n_out)    
         write (*,*)  CAPAC(n_out)  
         write (*,*)  CATDEF(n_out)    
         write (*,*)  RZEXC(n_out)    
         write (*,*)  srfexc(n_out)    
         write (*,*)  GHTCNT(:,n_out)    
         write (*,*)  TSURF(n_out)  
         write (*,*)  WESNN(:,n_out)    
         write (*,*)  HTSNNN(:,n_out)    
         write (*,*)  SNDZN(:,n_out)  
         write (*,*)  EVAP(n_out)    
         write (*,*)  SHFLUX(n_out)    
         write (*,*)  RUNOFF(n_out)  
         write (*,*)  EINT(n_out)    
         write (*,*)    ESOI(n_out)    
         write (*,*)    EVEG(n_out)    
         write (*,*)  ESNO(n_out)  
         write (*,*)  BFLOW(n_out)    
         write (*,*)  RUNSRF(n_out)    
         write (*,*)  SMELT(n_out)  
         write (*,*)  HLWUP(n_out)    
         write (*,*)  HLATN(n_out)    
         write (*,*)  QINFIL(n_out)    
         write (*,*)  AR1(n_out)    
         write (*,*)  AR2(n_out)    
         write (*,*)  RZEQ(n_out)  
         write (*,*)  GHFLUX(n_out)    
         write (*,*)  TPSN1(n_out)    
         write (*,*)  ASNOW0(n_out)    
         write (*,*)  TP1(n_out)    
         write (*,*)  TP2(n_out)  
         write (*,*)  TP3(n_out)   
         write (*,*)  TP4(n_out)    
         write (*,*)  TP5(n_out)    
         write (*,*)  TP6(n_out)   
       enddo 
      endif 
       
!rr ------------------------------------------------------------------      

      do n=1,nch
         rcst(n) = 1.E10
      end do
      
      do n=1,nch
       LH_SNOW(N)=0.
       SH_SNOW(N)=0.
       LWUP_SNOW(N)=0.
       LWDOWN_SNOW(N)=0.
       NETSW_SNOW(N)=0.
      end do

!**** ---------------------------------------------------
!**** PRE-PROCESS DATA AS NECESSARY:
!****
 
!rr      if (nch .gt. mxnch) then
!rr        write(*,*) 'catchment.f mxnch exceed: must be greater than',nch
!rr        stop
!rr        end if

      DO N=1,NCH
!       SATCAP(N) = 0.1 * ZLAI(N)
! change for convergence towards MOSAIC
!        SATCAP(N) = 0.2 * ZLAI(N) + 1.e-5
        SATCAP(N) = SATCAPFR * ZLAI(N) + 1.e-5
!!AMM        SATCAP(N) = 1.0 * ZLAI(N) + 1.e-5
        CSOIL(N)  = CSOIL_1

        if ( ityp(n) .ne. 1) CSOIL(N)  = CSOIL_2
        FCAN(N) = AMIN1( 1., AMAX1(0.,CAPAC(N)/SATCAP(N)) )
        SCLAI=amin1( 1., zlai(n)/2. )
        POTFRC(N)=FCAN(N)*SCLAI
        if(fcan(n) .lt. .1) POTFRC(N)=POTFRC(N)*(10.*fcan(n))

! Correction to RDC formulation -Randy Koster, 4/1/2011
!        RDCX(N)    = RDC(N)*SCLAI
        RDCX(N)    = RDC(N)

        DEDQA1X(N)  = AMAX1( DEDQA1(N), 500./ALHE )
        DEDTC1X(N)  = AMAX1( DEDTC1(N),   0. )
        DHSDQA1X(N) = AMAX1( DHSDQA1(N),   0. )
        DHSDTC1X(N) = AMAX1( DHSDTC1(N), -10. )

        DEDQA2X(N)  = AMAX1( DEDQA2(N), 500./ALHE )
        DEDTC2X(N)  = AMAX1( DEDTC2(N),   0. )
        DHSDQA2X(N) = AMAX1( DHSDQA2(N),   0. )
        DHSDTC2X(N) = AMAX1( DHSDTC2(N), -10. )

        DEDQA4X(N)  = AMAX1( DEDQA4(N), 500./ALHE )
        DEDTC4X(N)  = AMAX1( DEDTC4(N),   0. )
        DHSDQA4X(N) = AMAX1( DHSDQA4(N),   0. )
        DHSDTC4X(N) = AMAX1( DHSDTC4(N), -10. )


        qa1_orig(n)=qa1(n)
        qa2_orig(n)=qa2(n)
        qa4_orig(n)=qa4(n)
        tc1_orig(n)=tc1(n)
        tc2_orig(n)=tc2(n)
        tc4_orig(n)=tc4(n)
        
        ! Removed an #ifdef LAND_UPD [if (LAND_FIX)] block from here that constrained qa
        ! between qm and qsat.
        ! For energy conservation reasons, the same code was added in GEOS_CatchGridComp.F90 
        ! by Andrea Molod ca. 2016.
        ! Because the legacy LDASsa did not use GEOS_CatchGridComp.F90, the constraint was
        ! still needed here until the offline land software migrated to GEOSldas in 2020.
        ! - reichle, 14 Jan 2022

        if(ityp(n) .ge. 7) potfrc(n)=0.
!$$$        RA1(N)     = ONE / ( CD1(N) * max(UM(N),1.) )
!$$$        RA2(N)     = ONE / ( CD2(N) * max(UM(N),1.) )
!$$$        RA4(N)     = ONE / ( CD4(N) * max(UM(N),1.) )
!$$$        RAS(N)     = ONE / ( CDS(N) * max(UM(N),1.) ) 


!     HSNACC is an energy accounting term designed to account (among other,
!     lesser things) for the fact that snow is deposited at the snowpack
!     surface temperature while the atmosphere does not account for variations
!     in the heat content of deposited snow.
 
        HSNACC(N)=0.
        EVACC(N)=0.
        SHACC(N)=0.
        RUNSRF(N)=0.


!! !****   RESET LAND ICE VARIABLES, MAINTAINING TEMPS. AT EACH LAYER
!!         IF(ITYP(N) .EQ. 9) THEN
!! 
!!           ! This block of the code should no longer be used.
!!           ! If it is, Randy wants to know about it.
!!           ! reichle+koster, 12 Aug 2014
!!           write (*,*) 'catchment() encountered ityp==9. STOPPING.'
!!           stop 
!! 
!!           if(sum(htsnnn(:,n)+wesnn(:,n))==0.) then
!!               TSN1=tc1(n)-TF
!!               TSN2=tc1(n)-TF
!!               TSN3=tc1(n)-TF
!!             else
!!               TSN1=(HTSNNN(1,N)+WESNN(1,N)*ALHM)/(SCONST*WESNN(1,N)+1.e-5)
!!               TSN2=(HTSNNN(2,N)+WESNN(2,N)*ALHM)/(SCONST*WESNN(2,N)+1.e-5)
!!               TSN3=(HTSNNN(3,N)+WESNN(3,N)*ALHM)/(SCONST*WESNN(3,N)+1.e-5)
!!             endif
!!           WESNN(1,N)=.1
!!           WESNN(2,N)=.2
!!           WESNN(3,N)=.1
!!           HTSNNN(1,N)=-ALHM*WESNN(1,N)+TSN1*SCONST*WESNN(1,N)
!!           HTSNNN(2,N)=-ALHM*WESNN(2,N)+TSN1*SCONST*WESNN(2,N)
!!           HTSNNN(3,N)=-ALHM*WESNN(3,N)+TSN1*SCONST*WESNN(3,N)
!!           SNDZN(1,N)=WESNN(1,N)/.9
!!           SNDZN(2,N)=WESNN(2,N)/.9
!!           SNDZN(3,N)=WESNN(3,N)/.9
!!           POTFRC(N)=1.
!! 
!!           ENDIF

!****   RESET LAKE VARIABLES
        IF(ITYP(N) .EQ. 10) THEN
          CATDEF(N)=0.
          RZEXC(N)=0.
          SRFEXC(N)=0.
          SATCAP(N)=1000.
          CAPAC(N)=SATCAP(N)
          POTFRC(N)=1.
          ENDIF

        ENDDO

!**** ---------------------------------------------------
!**** DETERMINE INITIAL VALUE OF RZEQ:

      CALL RZEQUIL (                                                           &
                    NCH, CATDEF, VGWMAX,CDCR1,CDCR2,WPWET,POROS,               &
                    ars1,ars2,ars3,ara1,ara2,ara3,ara4,                        &
                    arw1,arw2,arw3,arw4,                                       &
                    RZEQOL                                                     &
                   )

      IF (BUG) THEN
        WRITE(*,*) 'RZEQUIL OK'
        ENDIF


!rr   switched order of call to partition and do-loop for emaxrz & emaxrt
!rr   because srfmn was not initialized in the computation of emaxrt
!rr   reichle, Oct 22, 2003

!**** PARTITION CATCHMENT INTO EVAPORATION SUBREGIONS:
!****
      CALL PARTITION (                                                         &
                      NCH,DTSTEP,DZSF,RZEXC,  RZEQOL,VGWMAX,CDCR1,CDCR2,       &
                      PSIS,BEE,poros,WPWET,                                    &
                      bf1, bf2,                                                &
                      ars1,ars2,ars3,ara1,ara2,ara3,ara4,                      &
                      arw1,arw2,arw3,arw4,BUG,                                 &
                      SRFEXC,CATDEF,RUNSRF,                                    &
                      AR1, AR2, AR4,srfmx,srfmn,  SWSRF1,SWSRF2,SWSRF4,RZI     &
                     )


      DO N=1,NCH
         TSOIL=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)

         ENTOT_ORIG(N) =                                                       &
              sum(HTSNNN(1:N_snow,N)) + TSOIL*CSOIL(N) + sum(GHTCNT(1:N_gt,N))

      ENDDO

      ! reichle, 1 May 2013, fix TWLAND<0 bug, use correct calculation in 
      
      CALL CATCH_CALC_WTOTL( NCH,                                              &
           CDCR2, WPWET, SRFEXC, RZEXC, CATDEF, CAPAC, WESNN,                  &
           WTOT_ORIG )
      
      DO N=1,NCH
         emaxrz=amax1(0.,RZEQOL(N)+RZEXC(N)-WPWET(N)*VGWMAX(N))
         EMAXRT(N)=(CAPAC(N)+emaxrz+(SRFEXC(N)-SRFMN(N)))/DTSTEP
         ENDDO
      
      do n=1,nch
        ar1old(n)=ar1(n)
        ar2old(n)=ar2(n)
        ar4old(n)=ar4(n)
        enddo

      IF (BUG) THEN
        WRITE(*,*) 'PARTITION OK'
        ENDIF

!**** ========================================================
!**** ENERGY BALANCES.

!**** COMPUTE "INITIAL ESTIMATE" OF HEAT FLUX TO DEEP SOIL (HFTDS)
!**** AND ITS DERIVATIVE WITH RESPECT TO TEMPERATURE (DHFTDS):

      DO N=1,NCH
        T1(1)=TC1(N)
        T1(2)=TC2(N)
        T1(3)=TC4(N)
        if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
           phi=POROS(N)
        else
           phi=PHIGT
        end if

	if (LAND_FIX) then
           ! zbar bug fix, - reichle, 16 Nov 2015
           ! zbar function - reichle, 29 Jan 2022 (minus sign applied in call to GNDTP0)
           ZBAR = catch_calc_zbar( bf1(n), bf2(n), catdef(n) )  
        else
           ! zbar minus sign applied in call to GNDTP0
	   ZBAR = SQRT(1.e-20+catdef(n)/bf1(n))+bf2(n)  ! old bug is wrong sign for bf2 here
	end if

        THETAF=.5
        DO LAYER=1,6
          HT(LAYER)=GHTCNT(LAYER,N)
          ENDDO

        CALL GNDTP0(                                               &
                    T1,phi,-1.*ZBAR,THETAF,                        &   ! note minus sign for zbar
                    HT,                                            &
                    fh21w,fH21i,fh21d,dfh21w,dfh21i,dfh21D,tp      &
                   )

        HFTDS1(N)=-FH21W
        HFTDS2(N)=-FH21I
        HFTDS4(N)=-FH21D
        DHFT1(N)=-DFH21W
        DHFT2(N)=-DFH21I
        DHFT4(N)=-DFH21D

        ENDDO 

      IF (BUG) THEN
        WRITE(*,*) 'HEAT FLUX INITIAL ESTIMATE OK'
        ENDIF

!**** -------------------------------------------------------------
!**** A. SNOW-FREE FRACTION.
!**** DETERMINE EVAPORATION, SENSIBLE HEAT FLUXES; UPDATE TEMPS:

      DO N=1,NCH
        PAR(N)    = PARDIR(N) + PARDIF(N) + 1.E-20
        PDIR(N)   = PARDIR(N) / PAR(N)
        TC1SF(N)  = TC1(N)
        TC2SF(N)  = TC2(N)
        TC4SF(N)  = TC4(N)
        ENDDO

      CALL RCUNST (                                                            &
                   NCH, ITYP, SUNANG, SQSCAT, PDIR, PAR, ZLAI, GREEN, BUG,     &
                   RCUN                                                        &
                  )

      IF (BUG) THEN
         WRITE(*,*) 'RCUNST OK'
      ENDIF

!**** 1. SATURATED FRACTION

      CALL ENERGY1 (                                                           &
                   NCH, DTSTEP, ITYP, UM, RCUN,                                &
                   ETURB1, DEDQA1X, DEDTC1X, HSTURB1, DHSDQA1X, DHSDTC1X,      &
                   QM,     RA1,   SWNETF,  HLWDWN, PSUR,                       &
                   RDCX,    HFTDS1, DHFT1,  QSAT1, DQS1, ALW1, BLW1,           &
                   EMAXRT,CSOIL,SWSRF1,POTFRC,.false.,                         &
                   TC1SF, QA1,                                                 &
                   EVAP1, SHFLUX1, HLWUP1, RX11, RX21, GHFLUX1, HSNACC1        &
                  )

      IF (BUG) THEN
        WRITE(*,*) 'ENERGY1 OK'
        ENDIF

!**** 2. SUBSATURATED BUT UNSTRESSED FRACTION


!CC    print*,'energy2'
      CALL ENERGY2 (                                                           &
                   NCH, DTSTEP, ITYP, UM, RCUN,                                &
                   ETURB2, DEDQA2X, DEDTC2X, HSTURB2, DHSDQA2X, DHSDTC2X,      &
                   QM,     RA2,   SWNETF,  HLWDWN, PSUR,                       &
                   RDCX,    HFTDS2, DHFT2, QSAT2, DQS2, ALW2, BLW2,            &
                   EMAXRT,CSOIL,SWSRF2,POTFRC,.false., RZI, WPWET,             &
                   TC2SF, QA2,                                                 &
                   EVAP2, SHFLUX2, HLWUP2, RX12, RX22, GHFLUX2, HSNACC2        &
                  )

      IF (BUG) THEN
        WRITE(*,*) 'ENERGY2 OK'
        ENDIF

!**** 3. WILTING FRACTION
!CC    print*,'energy4'
      CALL ENERGY4 (                                                           &
                   NCH, DTSTEP, ITYP,POROS, UM, RCST,                          &
                   ETURB4, DEDQA4X, DEDTC4X, HSTURB4, DHSDQA4X, DHSDTC4X,      &
                   QM,     RA4,   SWNETF,  HLWDWN, PSUR,                       &
                   RDCX,   HFTDS4, DHFT4, QSAT4, DQS4, ALW4, BLW4,             &
                   EMAXRT,CSOIL,SWSRF4,POTFRC,.false., WPWET,                  &
                   TC4SF, QA4,                                                 &
                   EVAP4, SHFLUX4, HLWUP4, RX14, RX24, GHFLUX4, HSNACC4        &
                  )

      IF (BUG) THEN
        WRITE(*,*) 'ENERGY4 OK'
        ENDIF

!**** COMPUTE EIRFRC
      DO N=1,NCH
         
        RTBS1=RX11(N)*RX21(N)/(RX11(N)+RX21(N)+1.E-20)
        EPFRC1=POTFRC(N) * ( RA1(N) + RTBS1 ) / ( RA1(N) + POTFRC(N)*RTBS1 )
         
        RTBS2=RX12(N)*RX22(N)/(RX12(N)+RX22(N)+1.E-20)
        EPFRC2=POTFRC(N) * ( RA2(N) + RTBS2 ) / ( RA2(N) + POTFRC(N)*RTBS2 )
         
        RTBS4=RX14(N)*RX24(N)/(RX14(N)+RX24(N)+1.E-20)
        EPFRC4=POTFRC(N) * ( RA4(N) + RTBS4 ) / ( RA4(N) + POTFRC(N)*RTBS4 )
         
        SUMEP=EPFRC1*EVAP1(N)*AR1(N)+EPFRC2*EVAP2(N)*AR2(N)+                   &
              EPFRC4*EVAP4(N)*AR4(N)
        SUME=EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N)
        
        !   "quick fix" gone wrong in global AMSR-E assimilation
        !   trying to correct while staying as close as possible to past fix
        !   30 July 2007, reichle: 
        !   

!           EIRFRC(N)=SUMEP/(SUME+1.E-20)

        if (SUME/=-1.e-20) then
            EIRFRC(N)=SUMEP/(SUME+1.E-20)
          else
            EIRFRC(N)=SUMEP/(SUME+2.E-20)
        end if

        ENDDO

      IF (BUG) THEN
        WRITE(*,*) 'EIRFRC OK'
        ENDIF

!**** --------------------------------------------------------
!**** B. SNOW-COVERED FRACTION.

      DO N=1,NCH

        TS     = TM(N) 
        T1(1)  = TC1(N)-TF 
        T1(2)  = TC2(N)-TF 
        T1(3)  = TC4(N)-TF 
        ! MB: avoid division by zero (AR1=0) in PEATCLSM equations
        IF(POROS(N) >= PEATCLSM_POROS_THRESHOLD) THEN
           AREA(1)= amax1(AR1(N),2.E-20)
        ELSE
           AREA(1) = AR1(N)
        END IF
        AREA(2)= AR2(N) 
        AREA(3)= AR4(N) 
        pr     = trainc(n)+trainl(n)+tsnow(n)+tice(n)+tfrzr(n)
        snowf  = tsnow(n)+tice(n)+tfrzr(n)
        dedea  = dedqas(n)*epsilon/psur(n) 
        dhsdea = dhsdqas(n)*epsilon/psur(n) 
        ea     = qm(n)*psur(n)/epsilon 
        esattc = qsats(n)*psur(n)/epsilon 
        desdtc = dqss(n)*psur(n)/epsilon 
        dedtc0  = dedtcs(n) + dedea*desdtc 
        dhsdtc0 = dhsdtcs(n) + dhsdea*desdtc 
        hsturb=hsturbs(n) 
        tkgnd(1)=1.8      !STEPH  
        tkgnd(2)=1.8 
        tkgnd(3)=1.8 
        raddn=hlwdwn(n)+swnets(n) 
        zc1=-(DZTSURF*0.5)
        hups=0.0 
 
!**** 1. RUN SNOW MODEL: 
 
        do i=1,N_SNOW
          wesn(i)=wesnn(i,n)
          htsnn(i)=htsnnn(i,n) 
          sndz(i)=sndzn(i,n) 
          tpsn(i)=0.0
          if (present(rconstit)) then
                DO K=1,N_Constit
                   RCONSTIT1(I,K)=RCONSTIT(N,I,K)
                   TOTDEP1 (K)   = totdepos(N,K)
                   ENDDO
             else
                DO K=1,N_Constit
                   RCONSTIT1(I,K)=0.
		   RMELT1(K)   = 0.
                   TOTDEP1 (K) = 0.
                   ENDDO
             endif
          enddo 

!     TPSN1 is used as input here, contradicts "declaration" as output only.
!     EnKF has been using tpsn1 as part of state vector, which should fix this.
!     reichle, 18 Nov 02

!     Removed tpsn1 from state vector, now compute tpsn1 from prognostics
!     in process
!     reichle, 29 May 03

        call StieglitzSnow_calc_tpsnow(htsnn(1),wesn(1),tsnowsrf,dum,ldum,ldum,.true.)
        tcs_orig(n)=tsnowsrf+tf
        if(wesn(1)+wesn(2)+wesn(3) .eq. 0.) tcs_orig(n)=                       &
                  amin1( tf, tc1_orig(n)*ar1(n)+tc2_orig(n)*ar2(n)+            &
                  tc4_orig(n)*(1.-ar1(n)-ar2(n)) )

        hlwtc=ALWS(N) + BLWS(N)*(TSNOWSRF+TF) 
        dhlwtc=BLWS(N)
        hcorr=0.

!! the field tpsn1(n) contains the value of TC(snow tile) that
!! came in from the catch grid comp, and catch internal state
!! spit it out here as a diagnostic
!! also spit out here the value of tsnowsrf+tf which is used
!! as the "original" value of TC for purposes of LW and turb fluxes

        tcsorig(N) = tcs_orig(n)
        tpsn1in(n) = tpsn1(n)    ! tpsn1 is "intent(out)", should NOT be used here, use catch_calc_tpsnow instead?  shouldn't this be the same as tcs_orig?  - reichle, 8/8/2014

        sumdepth=sum(sndz)

        CALL StieglitzSnow_snowrt(                                             &
                   N_sm, N_snow, MAPL_Land,                                    &  ! in   
                   CATCH_SNOW_MAXDEPTH, CATCH_SNOW_RHOFS, CATCH_SNOW_DZPARAM,  &  ! in   
                   t1, area, tkgnd, pr, snowf, ts, DTSTEP,                     &  ! in   
                   eturbs(n), dedtc0, hsturb, dhsdtc0, hlwtc, dhlwtc,          &  ! in   
                   raddn, zc1, totdep1,                                        &  ! in   
                   wesn, htsnn, sndz, RCONSTIT1,                               &  ! inout
                   hups, fices, tpsn, RMELT1,                                  &  ! out  
                   areasc, areasc0, pre, fhgnd,                                &  ! out  
                   EVSN, SHFLS, alhfsn, hcorr, ghfluxsno(n),                   &  ! out  
                   sndzsc, wesnprec, sndzprec, sndz1perc,                      &  ! out     
                   wesnperc, wesndens, wesnrepar, mltwtr,                      &  ! out     
                   excs, drho0, wesnbot, tksno, dtss                   )          ! out     


        FICESOUT(:,N)  = fices

        LH_SNOW(N)     = areasc*EVSN*ALHS
        SH_SNOW(N)     = areasc*SHFLS
        LWUP_SNOW(N)   = areasc*HUPS
        LWDOWN_SNOW(N) = areasc*HLWDWN(N)
        NETSW_SNOW(N)  = areasc*SWNETS(N)
 
        TPSN1(N) = TPSN(1)+TF 

        tpsn1out(n) = tpsn1(n)     ! why is tpsn1out needed?  same as tpsn1, reichle, 8/8/2014

        ! removed TPSN2 (never needed);
        ! renamed TPSN3 to TPSNB (bottom-layer snow temperature) 
        !  (for consistency with use of "N_snow" layers)
        ! reichle+koster, 12 Aug 2014 

        TPSNB(N) = TPSN(N_snow)+TF 
        SMELT(N) = PRE+sum(EXCS)            
        fh31w=fhgnd(1) 
        fh31i=fhgnd(2) 
        fh31d=fhgnd(3) 
        asnow(n) = areasc 
        asnow0(n)= areasc0 
        HSNACC(N) = HSNACC(N) + (1.-ASNOW(N))*                                 &
             (HSNACC1(N)*AR1(N)+HSNACC2(N)*AR2(N)+HSNACC4(N)*AR4(N))           &
             + hcorr
 

        do i=1,N_SNOW 
          wesnn(i,n)=wesn(i) 
          htsnnn(i,n)=htsnn(i) 
          sndzn(i,n)=sndz(i)
          if(present(rconstit)) then
             DO K=1,N_Constit
                RCONSTIT(N,I,K)=RCONSTIT1(I,K)
		RMELT (N,K) = RMELT1(K)
                ENDDO
             endif
          enddo 
 
        traincx(n)= trainc(n)*(1.-areasc) 
        trainlx(n)= trainl(n)*(1.-areasc)

!**** 2. UPDATE SURFACE TEMPERATURE

        DTC1SN=((-(FH31W/(area(1)+1.e-20))-HFTDS1(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT1(N)*DTSTEP/CSOIL(N))
        DTC2SN=((-(FH31I/(area(2)+1.e-20))-HFTDS2(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT2(N)*DTSTEP/CSOIL(N))
        DTC4SN=((-(FH31D/(area(3)+1.e-20))-HFTDS4(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT4(N)*DTSTEP/CSOIL(N))

        TC1SN=TC1(N)+DTC1SN
        IF((TC1SN-TPSNB(N))*(TC1(N)-TPSNB(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(1)*                                 &
                 (TC1SN-TPSNB(N))*CSOIL(N)/DTSTEP
          TC1SN=TPSNB(N)
          ENDIF

        TC2SN=TC2(N)+DTC2SN
        IF((TC2SN-TPSNB(N))*(TC2(N)-TPSNB(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(2)*                                 &
                 (TC2SN-TPSNB(N))*CSOIL(N)/DTSTEP
          TC2SN=TPSNB(N)
          ENDIF

        TC4SN=TC4(N)+DTC4SN
        IF((TC4SN-TPSNB(N))*(TC4(N)-TPSNB(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(3)*                                 &
                 (TC4SN-TPSNB(N))*CSOIL(N)/DTSTEP
          TC4SN=TPSNB(N)
          ENDIF



        TC1(N)=TC1SF(N)*(1.-AREASC)+TC1SN*AREASC
        TC2(N)=TC2SF(N)*(1.-AREASC)+TC2SN*AREASC
        TC4(N)=TC4SF(N)*(1.-AREASC)+TC4SN*AREASC
        
        EVSNOW(N)=EVSN
        esno(n)=evsnow(n)*asnow(n)*DTSTEP ! to have esno in mm/20min (03-17-99)
        SHFLUXS(N)=SHFLS 
        HLWUPS(N) =HUPS 
        GHFLUXS(N)=AREA(1)*(HFTDS1(N)+DHFT1(N)*DTC1SN) +                       &
                   AREA(2)*(HFTDS2(N)+DHFT2(N)*DTC2SN) +                       &
                   AREA(3)*(HFTDS4(N)+DHFT4(N)*DTC4SN)
        ENDDO 
 


      IF (BUG) THEN 
        WRITE(*,*) 'SNOW FRACTION OK' 
        ENDIF 
 
      DO N=1,NCH 
        HLATN(N)=(1.-ASNOW(N))*                                                &
              (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))*ALHE           &
              +ASNOW(N)*EVSNOW(N)*ALHS 
        EVAP(N)=(1.-ASNOW(N))*                                                 &
              (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))                &
              +ASNOW(N)*EVSNOW(N) 
        EVAPFR(N)=(1.-ASNOW(N))*                                               &
              (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))
        SHFLUX(N)=(1.-ASNOW(N))*                                               &
              (SHFLUX1(N)*AR1(N)+SHFLUX2(N)*AR2(N)+SHFLUX4(N)*AR4(N))          &
              +ASNOW(N)*SHFLUXS(N) 
        HLWUP(N)=(1.-ASNOW(N))*                                                &
              (HLWUP1(N)*AR1(N)+HLWUP2(N)*AR2(N)+HLWUP4(N)*AR4(N))             &
              +ASNOW(N)*HLWUPS(N) 
        SWLAND(N)=(1.-ASNOW(N))*SWNETF(N) + ASNOW(N)*SWNETS(N) 
        GHFLUX(N)=(1.-ASNOW(N))*                                               &
              (GHFLUX1(N)*AR1(N)+GHFLUX2(N)*AR2(N)+GHFLUX4(N)*AR4(N))          &
              +ASNOW(N)*GHFLUXS(N) 
        GHTSKIN(N)=(1.-ASNOW(N))*                                              &
              (GHFLUX1(N)*AR1(N)+GHFLUX2(N)*AR2(N)+GHFLUX4(N)*AR4(N))          &
              -ASNOW(N)*ghfluxsno(N)
        ENDDO 


      IF (BUG) THEN
        WRITE(*,*) 'ENERGY FLUXES OK'
        ENDIF


!****
!**** NOW ALLOW DEEPER SOIL TEMPERATURES TO BE UPDATED:

      DO N=1,NCH
        if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
           phi=POROS(N)
        else
           phi=PHIGT
        end if
        if (LAND_FIX) then
           ! zbar bug fix, - reichle, 16 Nov 2015
           ! zbar function - reichle, 29 Jan 2022 (minus sign applied in call to GNDTMP)
           ZBAR = catch_calc_zbar( bf1(n), bf2(n), catdef(n) )  
        else
           ! zbar minus sign applied in call to GNDTMP
	   ZBAR = SQRT(1.e-20+catdef(n)/bf1(n))+bf2(n)  ! old bug is wrong sign for bf2 here
        end if	
        THETAF=.5
        DO LAYER=1,6
          HT(LAYER)=GHTCNT(LAYER,N)
          ENDDO
        FH21=-GHFLUX(N)

        CALL GNDTMP(                                   &
              phi, -1.*zbar,                           &   ! note minus sign for zbar
              ht,                                      &
              xfice, tp, soilice,                      &
              DTS=dtstep, THETAF=thetaf, FH21=fh21)

        DO LAYER=1,6
          GHTCNT(LAYER,N)=HT(LAYER)
          ENDDO
        tp1(n)=tp(1)
        tp2(n)=tp(2)
        tp3(n)=tp(3)
        tp4(n)=tp(4)
        tp5(n)=tp(5)
        tp6(n)=tp(6)
        frice(n)=xfice
         
        ENDDO


      IF (BUG) THEN
        WRITE(*,*) 'DEEPER SOIL TEMPERATURES UPDATE OK'
        ENDIF

!**** ========================================================

!**** REMOVE EVAPORATED WATER FROM SURFACE RESERVOIRS:
!****
!**** (FIRST CORRECT FOR EXCESSIVE INTERCEPTION LOSS)

      DO N=1,NCH
        EINTX=EIRFRC(N)*EVAPFR(N)*DTSTEP
        IF(EINTX .GT. CAPAC(N)) THEN
          EDIF=(EINTX-CAPAC(N))/DTSTEP
!          EVACC(N)=EVACC(N)-EDIF
          EVAPFR(N)=EVAPFR(N)-EDIF
          EVAP(N)=EVAP(N)-EDIF
          HLATN(N)=HLATN(N)-EDIF*ALHE
          SHFLUX(N)=SHFLUX(N)+EDIF*ALHE
!          SHACC(N)=SHACC(N)+EDIF*ALHE
!          HSNACC(N)=HSNACC(N)+EDIF*ALHE
          EIRFRC(N)=CAPAC(N)/((EVAPFR(N)+1.E-20)*DTSTEP)
          ENDIF
        ENDDO


      CALL WUPDAT (                                                            &
                     NCH, DTSTEP, EVAPFR, SATCAP, TC1, RA1, RC,                &
                     RX11,RX21,RX12,RX22,RX14,RX24,                            &
                     AR1,AR2,AR4,CDCR1,EIRFRC,RZEQOL,srfmn,WPWET,VGWMAX,POROS, &
                     BF1, BF2, ARS1, ARS2, ARS3,                               &
                     CAPAC, RZEXC, CATDEF, SRFEXC,                             &
                     EINT, ESOI, EVEG                                          &
                    )

! ---------------------------------------------------------------------

      IF (BUG) THEN
        WRITE(*,*) 'WUPDAT OK'
        ENDIF

!**** REDISTRIBUTE MOISTURE BETWEEN RESERVOIRS:

      CALL RZDRAIN (                                                           &
                    NCH,DTSTEP,VGWMAX,SATCAP,RZEQOL,AR1,WPWET,                 &
                    tsa1,tsa2,tsb1,tsb2,atau,btau,CDCR2,poros,                 &
                    BF1, BF2, ARS1, ARS2, ARS3, BUG,                           &
                    CAPAC,RZEXC,SRFEXC,CATDEF,RUNSRF                           &
                    )

! ---------------------------------------------------------------------

      IF (BUG) THEN
        WRITE(*,*) 'RZDRAIN OK'
        ENDIF

!**** COMPUTE BASEFLOW FROM TOPMODEL EQUATIONS

      CALL BASE (                                                              &
                 NCH, DTSTEP,BF1, BF2, BF3, CDCR1, FRICE, COND, GNU,AR1, POROS,&
                 ARS1, ARS2, ARS3,                                             &
                 CATDEF, BFLOW                                                 &
                )

! ---------------------------------------------------------------------

      IF (BUG) THEN
        WRITE(*,*) 'BASE OK'
        ENDIF

!**** UPDATE CANOPY INTERCEPTION; DETERMINE THROUGHFALL RATES.

        CALL INTERC (                                                    &
             NCH, DTSTEP, FWETC, FWETL, TRAINLX, TRAINCX, SMELT,         &
             SATCAP, BUG,                                                &
             CAPAC,                                                      &
             THRUL_VOL, THRUC_VOL                                        &
             )

      IF (BUG) THEN
        WRITE(*,*) 'INTERC OK'
        ENDIF

!**** DETERMINE SURFACE RUNOFF AND INFILTRATION RATES:

        CALL SRUNOFF ( NCH, DTSTEP, UFW4RO, FWETC, FWETL,               &
             AR1, AR2, AR4, THRUL_VOL, THRUC_VOL,                       &
             FRICE, TP1, SRFMX, BUG,                                    & 
             VGWMAX, RZEQOL, POROS,                                     &
             SRFEXC, RZEXC, RUNSRF,                                     &
             QINFIL                                                     &
             )

      IF (BUG) THEN
        WRITE(*,*) 'SRUNOFF'
        ENDIF

!**** (ADD CHECK TO ENSURE RZEXC KEPT WITHIN BOUNDS IN SRUNOFF)
      
!**** RECOMPUTE RZEXC:

      CALL RZEQUIL (                                                           &
                    NCH, CATDEF, VGWMAX,CDCR1,CDCR2,WPWET,POROS,               &
                    ars1,ars2,ars3,ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4,    &
                    RZEQ                                                       &
                   )

      IF (BUG) THEN
        WRITE(*,*) 'RZEQUIL'
        ENDIF

      DO N=1,NCH
        ADJ=0.5*(RZEQOL(N)-RZEQ(N))
        RZEXC(N)=RZEXC(N)+ADJ
        CATDEF(N)=CATDEF(N)+ADJ
        ! make sure catdef does not become negative
        ! reichle, Aug 16, 2002
        IF(CATDEF(N) .LT. 0.) THEN
           RUNSRF(N)=RUNSRF(N)-CATDEF(N)/DTSTEP
           CATDEF(N)=0.
           ENDIF
         ENDDO

!**** Correct energy imbalance due to changing areas:
  
      ! note revised interface - reichle, 3 Apr 2012

      CALL CATCH_CALC_SOIL_MOIST (                                             &
          nch,dzsf,vgwmax,cdcr1,cdcr2,psis,bee,poros,wpwet,                    &
          ars1,ars2,ars3,ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4,bf1,bf2,      &
          srfexc,rzexc,catdef,                                                 &
          AR1, AR2, AR4,                                                       &
          sfmc, rzmc, prmc,                                                    &
          werror, sfmcun, rzmcun, prmcun  )

! Add differences due to adjustments to land moisture prognostics
      do n=1,nch
         if(werror(n) .le. 0.) runsrf(n)=runsrf(n)-werror(n)/dtstep
         if(werror(n) .gt. 0.) then
           edif=werror(n)/dtstep
           EVAP(N)=EVAP(N)-EDIF
           HLATN(N)=HLATN(N)-EDIF*ALHE
           EVEGFRC=EVEG(N)/(EVEG(N)+ESOI(N)+1.E-20)
           EVEG(N)=EVEG(N)-EDIF*EVEGFRC*DTSTEP
           ESOI(N)=ESOI(N)-EDIF*(1.-EVEGFRC)*DTSTEP
           SHFLUX(N)=SHFLUX(N)+EDIF*ALHE
!           EVACC(N)=EVACC(n)-EDIF
!           SHACC(N)=SHACC(N)+EDIF*ALHE
           endif
         enddo


    ! after revisions of calc_soil_moist() the call to partition is now obsolete 
    ! - reichle, 3 Apr 2012     
    !
    !  CALL PARTITION (                                                         &
    !                  NCH,DTSTEP,ITYP,DZSF,RZEXC,  RZEQOL,VGWMAX,CDCR1,CDCR2,  &
    !                  PSIS,BEE,poros,WPWET,                                    &
    !                  ars1,ars2,ars3,ara1,ara2,ara3,ara4,                      &
    !                  arw1,arw2,arw3,arw4,BUG,                                 &
    !                  SRFEXC,CATDEF,RUNSRF,                                    &
    !                  AR1, AR2, AR4,srfmx,srfmn,  SWSRF1,SWSRF2,SWSRF4,RZI     &
    !                 )

      do n=1,nch
        hold=csoil(n)*(ar1old(n)*tc1(n)+ar2old(n)*tc2(n)+ar4old(n)*tc4(n))
        hnew=csoil(n)*(ar1(n)*tc1(n)+ar2(n)*tc2(n)+ar4(n)*tc4(n))
        shflux(n)=shflux(n)-(hnew-hold)/dtstep
!        SHACC(N)=SHACC(N)-(hnew-hold)/dtstep
        enddo

!**** ---------------------------------------------------
!**** PROCESS DATA AS NECESSARY PRIOR TO RETURN:
!****
!**** ---------------------------------------------------


      DO N=1,NCH

        RUNOFF(N) = RUNSRF(N)+BFLOW(N)
        IF(CAPAC(N).LT.1.E-10) THEN
           RUNOFF(N) = RUNOFF(N)+CAPAC(N)/DTSTEP
           CAPAC(N) = 0.0
           endif

        EINT(N) = EINT(N) * ALHE / DTSTEP
        ESOI(N) = ESOI(N) * ALHE / DTSTEP
        EVEG(N) = EVEG(N) * ALHE / DTSTEP
        ESNO(N) = ESNO(N) * ALHS / DTSTEP
         
        TSOIL=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)
        TSURF(N)=(1.-ASNOW0(N))*TSOIL+ASNOW0(N)*TPSN1(N)

        if(asnow0(n) .eq. 0) then
          tpsn1(n)=amin1( TSURF(N),tf )
          tpsnb(n)=amin1( TSURF(N),tf )
        endif
        
        ! reichle, 1 May 2013, fix TWLAND<0 bug, use correct calculation in 

        CALL CATCH_CALC_WTOTL( 1, CDCR2(N), WPWET(N),                          &
             SRFEXC(N), RZEXC(N), CATDEF(N), CAPAC(N), WESNN(:,N),             &
             WTOT(N) )

        ENTOT(N) =                                                             &
             sum(HTSNNN(1:N_snow,N)) + TSOIL*CSOIL(N) + sum(GHTCNT(1:N_gt,N))
                
        WCHANGE(N) = (WTOT(N)-WTOT_ORIG(N))/DTSTEP
        ECHANGE(N) = (ENTOT(N)-ENTOT_ORIG(N))/DTSTEP

        !FSW_CHANGE IS THE CHANGE IN THE FREE-STANDING WATER, RELEVANT FOR PEATLAND ONLY
        FSW_CHANGE(N) = 0.
        IF(POROS(N) >= PEATCLSM_POROS_THRESHOLD) THEN
           pr = trainc(n)+trainl(n)+tsnow(n)+tice(n)+tfrzr(n)
           FSW_CHANGE(N) = PR - EVAP(N) - RUNOFF(N) - WCHANGE(N)
        ENDIF

! Perform check on sum of AR1 and AR2, to avoid calculation of negative 
! wilting fraction due to roundoff, outside of catchment:
        IF(AR1(N)+AR2(N) .GT. 1.) THEN
          IF(AR1(N) .GE. AR2(N)) AR1(N)=1.-AR2(N)
          IF(AR1(N) .LT. AR2(N)) AR2(N)=1.-AR1(N)
          ENDIF


       ! Engineering fix to dampen surface energy balance oscillations.
       ! Dampened temperatures are in tc[X]_00.  So far only used for 
       ! off-line land modeling.
       ! reichle, 12 Aug 2014

       ! Apply engineering fix for all vegetation classes.
       ! - reichle, 24 Nov 2015

       IF (LAND_FIX .OR. (ityp(N) .ne. 1)) THEN           
          call dampen_tc_oscillations(dtstep,tm(N),tc1_orig(N),tc1(N),     &
               tc1_00(N),dtc1)
          call dampen_tc_oscillations(dtstep,tm(N),tc2_orig(N),tc2(N),     &
               tc2_00(N),dtc2)
          call dampen_tc_oscillations(dtstep,tm(N),tc4_orig(N),tc4(N),     &
               tc4_00(N),dtc4)

          ! energy correction term resulting from resetting tc[X] to tc[X]_00          

          EACC_00(N) =                                                     &
               (dtc1*AR1(N) + dtc2*AR2(N) + dtc4*AR4(N))*CSOIL(N)/DTSTEP

       ELSE
          
          ! do not apply engineering fix for ityp=1 
          ! (because CSOIL for is much larger)

          TC1_00(N)=TC1(N)
          TC2_00(N)=TC2(N)
          TC4_00(N)=TC4(N)
          
          EACC_00(N) = 0.
          
       ENDIF

       ! To try dampening for AGCM output, uncomment the following
       ! *and* make sure EACC_00 is properly tracked. 
       !
       !TC1(N) = TC1_00(N)
       !TC2(N) = TC2_00(N)
       !TC4(N) = TC4_00(N)


       ! Optional output of TC and QA before they are changed for the
       ! benefit of the GCM.  Use for off-line land modeling, includes
       ! engineering fix for dampening of oscillations (except for 
       ! broadleaf evergreen tiles because for these tiles the fix
       ! would result in very large EACC_0 garbage terms)
       ! - reichle, 12 Aug 2014
                 

       IF (PRESENT(TC1_0 ))  TC1_0( N)=TC1_00( N)
       IF (PRESENT(TC2_0 ))  TC2_0( N)=TC2_00( N)
       IF (PRESENT(TC4_0 ))  TC4_0( N)=TC4_00( N)
       IF (PRESENT(QA1_0 ))  QA1_0( N)=QA1(    N)
       IF (PRESENT(QA2_0 ))  QA2_0( N)=QA2(    N)
       IF (PRESENT(QA4_0 ))  QA4_0( N)=QA4(    N)    
       IF (PRESENT(EACC_0))  EACC_0(N)=EACC_00(N)    

       
!       Revise values of qa and tc prior to return to AGCM, to ensure that 
!       AGCM uses evaporation and sensible heat flux rates consistent
!       with those of land model (in situations of imposed water limits, 
!       etc.).

        QA1X=QA1_ORIG(N)
        QA2X=QA2_ORIG(N)
        QA4X=QA4_ORIG(N)
        TC1X=TC1_ORIG(N)
        TC2X=TC2_ORIG(N)
        TC4X=TC4_ORIG(N)
        TCSX=TCS_ORIG(N)

        IF(DEDQA1(N) .NE. 0.) QA1X = QA1_ORIG(N)+(EVAP1(N)-ETURB1(N))/DEDQA1(N)
        IF(DEDQA2(N) .NE. 0.) QA2X = QA2_ORIG(N)+(EVAP2(N)-ETURB2(N))/DEDQA2(N)
        IF(DEDQA4(N) .NE. 0.) QA4X = QA4_ORIG(N)+(EVAP4(N)-ETURB4(N))/DEDQA4(N)

        IF(DHSDTC1(N) .NE. 0.) TC1X = TC1_ORIG(N)+(SHFLUX1(N)-HSTURB1(N))/     &
                           DHSDTC1(N)
        IF(DHSDTC2(N) .NE. 0.) TC2X = TC2_ORIG(N)+(SHFLUX2(N)-HSTURB2(N))/     &
                           DHSDTC2(N)
        IF(DHSDTC4(N) .NE. 0.) TC4X = TC4_ORIG(N)+(SHFLUX4(N)-HSTURB4(N))/     &
                           DHSDTC4(N)

!       Ensure that modifications made to QA and TC are not too large:

! JP added patch to exlude correction under snow:
        IF ( (.NOT. LAND_FIX) .OR. (ASNOW0(N) .EQ. 0. ) ) THEN          
           IF(ABS(QA1X-QA1(N)) .LE. 0.5*QA1(N)) THEN
              QA1(N)=QA1X
           ELSE
              QA1(N)=QA1(N)+SIGN(0.5*QA1(N),QA1X-QA1(N))
           ENDIF

           IF(ABS(QA2X-QA2(N)) .LE. 0.5*QA2(N)) THEN
              QA2(N)=QA2X
           ELSE
              QA2(N)=QA2(N)+SIGN(0.5*QA2(N),QA2X-QA2(N))
           ENDIF

           IF(ABS(QA4X-QA4(N)) .LE. 0.5*QA4(N)) THEN
              QA4(N)=QA4X
           ELSE
              QA4(N)=QA4(N)+SIGN(0.5*QA4(N),QA4X-QA4(N))
           ENDIF

           IF(ABS(TC1X-TC1(N)) .LE. 10.) THEN
              TC1(N)=TC1X
           ELSE
              TC1(N)=TC1(N)+SIGN(10.,TC1X-TC1(N))
           ENDIF

           IF(ABS(TC2X-TC2(N)) .LE. 10.) THEN
              TC2(N)=TC2X
           ELSE
              TC2(N)=TC2(N)+SIGN(10.,TC2X-TC2(N))
           ENDIF

           IF(ABS(TC4X-TC4(N)) .LE. 10.) THEN
              TC4(N)=TC4X
           ELSE
              TC4(N)=TC4(N)+SIGN(10.,TC4X-TC4(N))
           ENDIF
	ENDIF ! LAND_FIX and ASNOW=0

! EVACC and SHACC are the differences ("errors") between what the land surface
! model computes for the evaporation and sensible heat flux and what the AGCM
! will compute, since the latter is forced to compute these fluxes
! based on changes in near-surface humidity and temperature only -- and 
! because the fractions (AR1, AR2, AR4, and ASNOW0) provided back to the
! AGCM are not the same as those that went into computing the land model's
! fluxes, since the land model has to update those areas (based on the fluxes)
! as a matter of course.

        EVAPX1=ETURB1(N)+DEDQA1(N)*(QA1(N)-QA1_ORIG(N))
        EVAPX2=ETURB2(N)+DEDQA2(N)*(QA2(N)-QA2_ORIG(N))
        EVAPX4=ETURB4(N)+DEDQA4(N)*(QA4(N)-QA4_ORIG(N))
        EVAPXS=ETURBS(N)+DEDQAS(N)*DQSS(N)*(TPSN1(N)-TCS_ORIG(N))
        EVACC(N)=        (1.-ASNOW0(N))*                                       &
                        ( AR1(N)*EVAPX1+                                       &
                          AR2(N)*EVAPX2+                                       &
                          AR4(N)*EVAPX4 )                                      &
                      +  ASNOW0(N)*EVAPXS
        EVACC(N)=EVAP(N)-EVACC(N)


        ! added term for latent heat flux correction, reichle+qliu, 9 Oct 2008

        if(present(lhacc)) then
           LHACC(N)=  ALHE*(1.-ASNOW0(N))*                           &
                ( AR1(N)*EVAPX1+                                     &
                  AR2(N)*EVAPX2+                                     &
                  AR4(N)*EVAPX4 )                                    &
                + ALHS*ASNOW0(N)*EVAPXS
           LHACC(N)=HLATN(N)-LHACC(N)
           end if

        SHFLUXX1=HSTURB1(N)+DHSDTC1(N)*(TC1(N)-TC1_ORIG(N))
        SHFLUXX2=HSTURB2(N)+DHSDTC2(N)*(TC2(N)-TC2_ORIG(N))
        SHFLUXX4=HSTURB4(N)+DHSDTC4(N)*(TC4(N)-TC4_ORIG(N))
        SHFLUXXS=HSTURBS(N)+DHSDTCS(N)*(TPSN1(N)-TCS_ORIG(N))
        SHACC(N)=         (1.-ASNOW0(N))*                                      &
                        ( AR1(N)*SHFLUXX1+                                     &
                          AR2(N)*SHFLUXX2+                                     &
                          AR4(N)*SHFLUXX4 )                                    &
                      + ASNOW0(N)*SHFLUXXS
        SHACC(N)=SHFLUX(N)-SHACC(N)


! **** SPECIAL DIAGNOSTICS FOR AR5 DECADAL RUNS

        ! the following assumes that fices is unchanged, otherwise may need to update FICESOUT
        ! - reichle,  4 Oct 2023

        CALL STIEGLITZSNOW_CALC_TPSNOW(N_SNOW, HTSNNN(:,N), WESNN(:,N), TPSN, FICES)

        !AVET_SNOW(N)=(TPSN(1)+TF)*WESNN(1,N) + (TPSN(2)+TF)*WESNN(2,N) +       &
        !     (TPSN(3)+TF)*WESNN(3,N)

        tmpvec_Nsnow = (tpsn(1:N_snow)+tf)*wesnn(1:N_snow,N)

        AVET_SNOW(N) = sum(tmpvec_Nsnow(1:N_snow))
           
        WAT_10CM(N)=0.1*(RZEQ(N)+RZEXC(N))+SRFEXC(N)

        TOTWAT_SOIL(N)=(CDCR2(N)/(1.-WPWET(N))-CATDEF(N)+RZEXC(N)+SRFEXC(N))

        TOTICE_SOIL(N)=TOTWAT_SOIL(N)*FRICE(N)


        ENDDO  ! N=1,NCH (PROCESS DATA AS NECESSARY PRIOR TO RETURN)   

      if(numout.ne.0) then
       do i = 1,numout
         n_out = n_outs(i)
         write (*,*) 'OUTPUT catchment arguments for n_out: ',n_out
         write (*,*) 'qsats = ', qsats(n_out)  
         write (*,*)  'TC1 = ', tc1(n_out)    
         write (*,*)  'TC2 = ', tc2(n_out)  
         write (*,*)  'TC4 = ', tc4(n_out)  
         write (*,*)  'TPSN1 = ', tpsn1(n_out)
         write (*,*)  'TSURF = ',    TSURF(n_out)  
         write (*,*)  'WESNN(1), = ',    WESNN(1,n_out)    
         write (*,*)  'HTSNNN(1), = ',    HTSNNN(1,n_out)    
         write (*,*)  'SNDZN(1), = ',  SNDZN(1,n_out)  
         write (*,*)  'TP1 = ', TP(1)         
         
         write (*,*) NCH  
         write (*,*) DTSTEP  
         write (*,*) UFW4RO
         write (*,*) ITYP(n_out)  
         write (*,*) TRAINC(n_out)    
         write (*,*) TRAINL(n_out)    
         write (*,*) TSNOW(n_out)    
         write (*,*) TICE(n_out)    
         write (*,*) TFRZR(n_out)    
         write (*,*) UM(n_out)  
         write (*,*) ETURB1(n_out)    
         write (*,*) DEDQA1(n_out)    
         write (*,*) DEDTC1(n_out)    
         write (*,*)  HSTURB1(n_out)    
         write (*,*) DHSDQA1(n_out)    
         write (*,*)  DHSDTC1(n_out)  
         write (*,*)  ETURB2(n_out)    
         write (*,*)  DEDQA2(n_out)    
         write (*,*)  DEDTC2(n_out)    
         write (*,*)  HSTURB2(n_out)    
         write (*,*)  DHSDQA2(n_out)    
         write (*,*)  DHSDTC2(n_out)  
         write (*,*)  ETURB4(n_out)    
         write (*,*)  DEDQA4(n_out)    
         write (*,*)  DEDTC4(n_out)    
         write (*,*)  HSTURB4(n_out)    
         write (*,*)  DHSDQA4(n_out)    
         write (*,*)  DHSDTC4(n_out)  
         write (*,*)  ETURBS(n_out)    
         write (*,*)  DEDQAS(n_out)    
         write (*,*)  DEDTCS(n_out)    
         write (*,*)  HSTURBS(n_out)    
         write (*,*)  DHSDQAS(n_out)    
         write (*,*)  DHSDTCS(n_out)  
         write (*,*)  TM(n_out)    
         write (*,*)  QM(n_out)    
         write (*,*)  ra1(n_out)    
         write (*,*)  ra2(n_out)    
         write (*,*)  ra4(n_out)    
         write (*,*)  raS(n_out)    
         write (*,*)  SUNANG(n_out)    
         write (*,*)  PARDIR(n_out)    
         write (*,*)  PARDIF(n_out)  
         write (*,*)  SWNETF(n_out)    
         write (*,*)  SWNETS(n_out)    
         write (*,*)   HLWDWN(n_out)    
         write (*,*)  PSUR(n_out)    
         write (*,*)   ZLAI(n_out)    
         write (*,*)    GREEN(n_out)    
         write (*,*)   Z2(n_out)  
         write (*,*)  SQSCAT(n_out)    
         write (*,*)  RSOIL1(n_out)    
         write (*,*)  RSOIL2(n_out)    
         write (*,*)    RDC(n_out)    
!         write (*,*)     U2FAC(n_out)  
         write (*,*)  QSAT1(n_out)    
         write (*,*)  DQS1(n_out)    
         write (*,*)  ALW1(n_out)    
         write (*,*)  BLW1(n_out)  
         write (*,*)  QSAT2(n_out)    
         write (*,*)  DQS2(n_out)    
         write (*,*)  ALW2(n_out)    
         write (*,*)  BLW2(n_out)  
         write (*,*)  QSAT4(n_out)    
         write (*,*)  DQS4(n_out)    
         write (*,*)  ALW4(n_out)    
         write (*,*)  BLW4(n_out)  
         write (*,*)  QSATS(n_out)    
         write (*,*)  DQSS(n_out)    
         write (*,*)  ALWS(n_out)    
         write (*,*)  BLWS(n_out)  
         write (*,*)  BF1(n_out)    
         write (*,*)  BF2(n_out)    
         write (*,*)  BF3(n_out)    
         write (*,*)  VGWMAX(n_out)  
         write (*,*)  CDCR1(n_out)    
         write (*,*)  CDCR2(n_out)    
         write (*,*)  psis(n_out)    
         write (*,*)  bee(n_out)    
         write (*,*)  poros(n_out)    
         write (*,*)  wpwet(n_out)    
         write (*,*)  cond(n_out)    
         write (*,*)  gnu(n_out)  
         write (*,*)  ARS1(n_out)    
         write (*,*)  ARS2(n_out)    
         write (*,*)  ARS3(n_out)    
         write (*,*)  ARA1(n_out)    
         write (*,*)  ARA2(n_out)    
         write (*,*)  ARA3(n_out)    
         write (*,*)  ARA4(n_out)  
         write (*,*)  ARW1(n_out)    
         write (*,*)  ARW2(n_out)    
         write (*,*)  ARW3(n_out)    
         write (*,*)  ARW4(n_out)  
         write (*,*)  tsa1(n_out)    
         write (*,*)  tsa2(n_out)    
         write (*,*)  tsb1(n_out)    
         write (*,*)  tsb2(n_out)    
         write (*,*)  atau(n_out)    
         write (*,*)  btau(n_out)    
         write (*,*)  BUG
         write (*,*)  TC1(n_out)    
         write (*,*)  TC2(n_out)    
         write (*,*)  TC4(n_out)    
         write (*,*)  QA1(n_out)    
         write (*,*)  QA2(n_out)    
         write (*,*)  QA4(n_out)    
         write (*,*)  CAPAC(n_out)  
         write (*,*)  CATDEF(n_out)    
         write (*,*)  RZEXC(n_out)    
         write (*,*)  srfexc(n_out)    
         write (*,*)  GHTCNT(:,n_out)    
         write (*,*)  TSURF(n_out)  
         write (*,*)  WESNN(:,n_out)    
         write (*,*)  HTSNNN(:,n_out)    
         write (*,*)  SNDZN(:,n_out)  
         write (*,*)  EVAP(n_out)    
         write (*,*)  SHFLUX(n_out)    
         write (*,*)  RUNOFF(n_out)  
         write (*,*)  EINT(n_out)    
         write (*,*)    ESOI(n_out)    
         write (*,*)    EVEG(n_out)    
         write (*,*)  ESNO(n_out)  
         write (*,*)  BFLOW(n_out)    
         write (*,*)  RUNSRF(n_out)    
         write (*,*)  SMELT(n_out)  
         write (*,*)  HLWUP(n_out)    
         write (*,*)  HLATN(n_out)    
         write (*,*)  QINFIL(n_out)    
         write (*,*)  AR1(n_out)    
         write (*,*)  AR2(n_out)    
         write (*,*)  RZEQ(n_out)  
         write (*,*)  GHFLUX(n_out)    
         write (*,*)  TPSN1(n_out)    
         write (*,*)  ASNOW0(n_out)    
         write (*,*)  TP1(n_out)    
         write (*,*)  TP2(n_out)  
         write (*,*)  TP3(n_out)   
         write (*,*)  TP4(n_out)    
         write (*,*)  TP5(n_out)    
         write (*,*)  TP6(n_out)   
       enddo
      endif
    

      RETURN
      END SUBROUTINE CATCHMENT

!**** ===================================================
!**** ///////////////////////////////////////////////////
!**** ===================================================
!****
!**** [ BEGIN RCUNST ]
!****
      SUBROUTINE RCUNST (                                                      &
                         NCH, ITYP, SUNANG, SQSCAT, PDIR,                      &
                         PAR, ZLAI, GREEN,BUG,                                 &
                         RCUN                                                  &
                        )
!****
!****     This subroutine calculates the unstressed canopy resistance.
!**** (p. 1353, Sellers 1985.)  Extinction coefficients are computed first.
!****
      IMPLICIT NONE
!****
      LOGICAL, INTENT(IN) :: BUG
      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN), DIMENSION(NCH) :: SUNANG, PDIR, PAR, ZLAI,             &
            SQSCAT, GREEN

      REAL, INTENT(OUT), DIMENSION(NCH) :: RCUN

      REAL, DIMENSION(NTYPS) :: VGCHIL, VGZMEW, VGRST1, VGRST2, VGRST3


      INTEGER CHNO
      REAL  RHO4, EXTK1, EXTK2, RCINV, GAMMA, EKAT, DUM1, DUM2, DUM3,          &
              AA, BB, ZK, CC


      DATA VGCHIL /        0.1,        0.25,        0.01,        -0.3,         &
                          0.01,        0.20/

      DATA VGZMEW/      0.9809,      0.9638,      0.9980,      1.0773,         &
                        0.9980,      0.9676/

      DATA VGRST1 /     2335.9,      9802.2,      2869.7,      2582.0,         &
                       93989.4,      9802.2/

      DATA VGRST2 /        0.0,        10.6,         3.7,         1.1,         &
                          0.01,        10.6/

      DATA VGRST3 /      153.5,       180.0,       233.0,       110.0,         &
                         855.0,       180.0/



      DO 100 ChNo = 1, NCH

!**** First compute optical parameters.
!**** (Note: CHIL is constrained to be >= 0.01, as in SiB calcs.)

      AA = 0.5 - (0.633 + 0.330*VGCHIL(ITYP(ChNo)))*VGCHIL(ITYP(ChNo))
      BB = 0.877 * ( ONE - 2.*AA )
      CC =  ( AA + BB*SUNANG(ChNo) ) / SUNANG(ChNo)

      EXTK1 =  CC * SQSCAT(ChNo)
      EXTK2 = (ONE / VGZMEW(ITYP(ChNo))) * SQSCAT(ChNo)

      DUM1 =      PDIR(ChNo)  *   CC
      DUM2 = (ONE-PDIR(ChNo)) * ( BB*(ONE/3.+PIE/4.) + AA*1.5 )

!**** Bound extinction coefficient by 50./ZLAI:

      ZK =     PDIR(ChNo) *AMIN1( EXTK1, 50./ZLAI(ChNo) ) +                    &
          (ONE-PDIR(ChNo))*AMIN1( EXTK2, 50./ZLAI(ChNo) )

!**** Now compute unstressed canopy resistance:

      GAMMA = VGRST1(ITYP(ChNo)) / VGRST3(ITYP(ChNo)) + VGRST2(ITYP(ChNo))

      EKAT = EXP( ZK*ZLAI(ChNo) )
      RHO4 = GAMMA / (PAR(ChNo) * (DUM1 + DUM2))

      DUM1 = (VGRST2(ITYP(ChNo)) - GAMMA) / (GAMMA + 1.E-20)
      DUM2 = (RHO4 * EKAT + ONE) / (RHO4 + ONE)
      DUM3 = ZK * VGRST3(ITYP(ChNo))

      RCINV = ( DUM1*ALOG(DUM2) + ZK*ZLAI(ChNo) ) / DUM3         
      rcinv = amax1(rcinv,0.)

      RCUN(ChNo) = ONE / (RCINV * GREEN(ChNo) + 1.E-10)

 100  CONTINUE


      RETURN
      END SUBROUTINE RCUNST
!****
!**** [ END RCUNST ]


!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****

      SUBROUTINE energy1 (                                                     &
                       NCH, DTSTEP, ITYP, UM, RCIN,                            &
                       ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,         &
                       QM,     RA,   SWNET,  HLWDWN, PSUR,                     &
                       RDC,    HFTDS, DHFTDS,                                  &
                       QSATTC, DQSDTC, ALWRAD, BLWRAD,                         &
                       EMAXRT,CSOIL,SWSRF,POTFRC,BUG,                          &
                       TC, QA,                                                 &
                       EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX, HSNACC           &
                       )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RCIN, ETURB, HSTURB, QM, RA,     &
                    SWNET, HLWDWN, PSUR, RDC, HFTDS, DHFTDS, QSATTC, DQSDTC,   &
                    ALWRAD, BLWRAD, EMAXRT, CSOIL, SWSRF, POTFRC, DEDQA,       &
                    DEDTC, DHSDQA, DHSDTC
      LOGICAL, INTENT(IN) :: BUG


      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, RX1, RX2,      &
                    GHFLUX, HSNACC


      INTEGER ChNo, N
      REAL, DIMENSION(NCH) :: VPDSTR, ESATTX, VPDSTX, FTEMP, RC, EAX, TX,      &
                    RCX, DRCDTC, DUMMY,  FTEMPX, DRCDEA, DEDEA, DHSDEA, EM,    &
                    ESATTC, DESDTC, EA
      REAL  DELTC, DELEA
 

!**** - - - - - - - -
!****
      DATA DELTC /0.01/, DELEA /0.001/
!****
!**** --------------------------------------------------------------------- 
!****
!**** Expand data as specified by ITYP into arrays of size NCH.
!****

!****
!**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

!      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
!      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
!      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
!      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)


 100  CONTINUE

!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
!****          ( VPDFAC  computes the vapor pressure deficit stress term,
!****           TMPFAC  computes the temperature stress term,
!****           RSURFP  computes rc given a parallel resist. from the
!****                   surface,
!****           RCANOP  computes rc corrected for snow and interception.)
!****

      CALL VPDFAC (                                                            &
                   NCH,  ITYP,  ESATTC, EA,                                    &
                   VPDSTR                                                      &
                  )


      CALL TMPFAC (                                                            &
                   NCH,  ITYP, TC,                                             &
                   FTEMP                                                       &
                  )


      DO N=1,NCH
        RC(N)=RCIN(N)/(VPDSTR(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (                                                           &
                   NCH, UM, RDC, SWSRF,ESATTC, EA,                             &
                   RC,                                                         &
                   RX1, RX2                                                    &
                  )


      CALL RCANOP (                                                            &
                   NCH, RA, ETURB, POTFRC,                                     &
                   RC                                                          &
                  )

!****
!**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
!****
!**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
!****

      DO ChNo = 1, NCH
        TX(ChNo) = TC(ChNo) + DELTC
        ESATTX(ChNo) = ESATTC(ChNo) + DESDTC(CHNO) * DELTC
        EAX(ChNo) = EA(ChNo) + DELEA
        ENDDO

!****
!**** temperature:
      CALL VPDFAC (NCH, ITYP, ESATTX, EA, VPDSTX)
      CALL TMPFAC (NCH, ITYP, TX, FTEMPX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMPX(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (NCH, UM, RDC, SWSRF, ESATTX, EA,                           &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO  ChNo = 1, NCH
        DRCDTC(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELTC
        ENDDO
!****

!**** vapor pressure:
      CALL VPDFAC (NCH, ITYP, ESATTC, EAX, VPDSTX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (NCH, UM, RDC, SWSRF, ESATTC, EAX,                          &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO ChNo = 1, NCH
         DRCDEA(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELEA
         ENDDO
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 2: Solve the energy balance at the surface.
!****
      CALL FLUXES (                                                            &
                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                     &
                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,             &
                       RC, DRCDEA, DRCDTC, SWNET, HLWDWN, ALWRAD, BLWRAD,      &
                       EM,  CSOIL,   PSUR, EMAXRT,  HFTDS, DHFTDS,             &
                       TC,     EA,                                             &
                     EVAP, SHFLUX,  HLWUP, GHFLUX,  HSNACC                     &
                  )

!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****
!****
!**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE

!****
      RETURN
      END SUBROUTINE energy1

!****
!**** [ END CHIP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
      SUBROUTINE energy2 (                                                     &
                       NCH, DTSTEP, ITYP, UM, RCIN,                            &
                       ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,         &
                       QM,     RA,   SWNET,  HLWDWN, PSUR,                     &
                       RDC,    HFTDS, DHFTDS,                                  &
                       QSATTC, DQSDTC, ALWRAD, BLWRAD,                         &
                       EMAXRT,CSOIL,SWSRF,POTFRC,BUG,RZI, WPWET,               &
                       TC, QA,                                                 &
                       EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX, HSNACC           &
                       )


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) ::  ITYP

      REAL, INTENT(IN) ::  DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) ::  UM, RCIN, ETURB, HSTURB, QM, RA,    &
                    SWNET, HLWDWN, PSUR, RDC, HFTDS, DHFTDS, QSATTC, DQSDTC,   &
                    ALWRAD, BLWRAD, EMAXRT, CSOIL, SWSRF, POTFRC, RZI, WPWET,  &
                    DEDQA, DEDTC, DHSDQA, DHSDTC

      LOGICAL, INTENT(IN) ::   BUG


      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, RX1, RX2,      &
                    GHFLUX, HSNACC



      INTEGER ChNo, N
      REAL, DIMENSION(NCH) :: VPDSTR, ESATTX, VPDSTX, FTEMP, RC, EAX, TX,      &
                    RCX, DRCDTC, DUMMY, FTEMPX, DRCDEA, DEDEA, DHSDEA, EM,     &
                    ESATTC, DESDTC, EA, RSTFAC
!      REAL DELTC, DELEA, STEXP, ATRANS, ASTRFR
      REAL DELTC, DELEA, ATRANS

      DATA DELTC /0.01/, DELEA /0.001/


!****
!**** --------------------------------------------------------------------- 
!****
!**** Expand data as specified by ITYP into arrays of size NCH.
!****

!****
!**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

!      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
!      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
!      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
!      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
!****          ( VPDFAC  computes the vapor pressure deficit stress term,
!****           TMPFAC  computes the temperature stress term,
!****           RSURFP  computes rc given a parallel resist. from the
!****                   surface,
!****           RCANOP  computes rc corrected for snow and interception.)
!****


!**** Compute water stress effect (RDK 03/30/06)

      CALL VPDFAC (                                                            &
                   NCH,  ITYP,  ESATTC, EA,                                    &
                   VPDSTR                                                      &
                  )


      CALL TMPFAC (                                                            &
                   NCH,  ITYP, TC,                                             &
                   FTEMP                                                       &
                  )

      DO N=1,NCH
        RC(N)=RCIN(N)/(VPDSTR(N)*FTEMP(N)+1.E-20)
        ENDDO


      DO CHNO = 1, NCH
        ATRANS = WPWET(CHNO)+ASTRFR*(1.-WPWET(CHNO))
        RSTFAC(CHNO)=AMAX1( 1.E-3, (RZI(CHNO)-WPWET(CHNO))/                    &
                                         (ATRANS-WPWET(CHNO)) )
        RSTFAC(CHNO)=AMIN1( 1., RSTFAC(CHNO))
        RC(ChNo) = RC(CHNO) / RSTFAC(CHNO)**STEXP
        RC(CHNO) = AMIN1 (RC(CHNO) , 1.E10)
        ENDDO


      CALL RSURFP2 (                                                           &
                   NCH, UM, RDC, SWSRF, ESATTC, EA, WPWET,                     &
                   RC,                                                         &
                   RX1, RX2                                                    &
                  )


      CALL RCANOP (                                                            &
                   NCH, RA, ETURB, POTFRC,                                     &
                   RC                                                          &
                  )

!****
!**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
!****
!**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
!****

      DO ChNo = 1, NCH
        TX(ChNo) = TC(ChNo) + DELTC
        ESATTX(ChNo) = ESATTC(ChNo) + DESDTC(CHNO) * DELTC
        EAX(ChNo) = EA(ChNo) + DELEA
        ENDDO

!****
!**** temperature:
      CALL VPDFAC (NCH, ITYP, ESATTX, EA, VPDSTX)
      CALL TMPFAC (NCH, ITYP, TX, FTEMPX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMPX(N)+1.E-20)
        RCX(N) = RCX(N) / RSTFAC(N)**STEXP
        RCX(N) = AMIN1 (RCX(N) , 1.E10)
        ENDDO

      CALL RSURFP2 (NCH, UM, RDC, SWSRF,ESATTX, EA, WPWET,                     &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO  ChNo = 1, NCH
        DRCDTC(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELTC
        ENDDO
!****

!**** vapor pressure:
      CALL VPDFAC (NCH, ITYP, ESATTC, EAX, VPDSTX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMP(N)+1.E-20)
        RCX(N) = RCX(N) / RSTFAC(N)**STEXP
        RCX(N) = AMIN1 (RCX(N) , 1.E10)
        ENDDO

      CALL RSURFP2 (NCH, UM, RDC, SWSRF, ESATTC, EAX, WPWET,                   &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO ChNo = 1, NCH
         DRCDEA(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELEA
         ENDDO
!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 2: Solve the energy balance at the surface.
!****
      CALL FLUXES (                                                            &
                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                     &
                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,             &
                       RC, DRCDEA, DRCDTC,SWNET, HLWDWN, ALWRAD, BLWRAD,       &
                       EM,  CSOIL,   PSUR, EMAXRT, HFTDS, DHFTDS,              &
                       TC,     EA,                                             &
                     EVAP, SHFLUX,  HLWUP, GHFLUX, HSNACC                      &
                  )

!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****
!****
!**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE
!****

      RETURN
      END SUBROUTINE energy2

!****
!**** [ END CHIP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
      SUBROUTINE energy4 (                                                     &
                       NCH, DTSTEP, ITYP, POROS,UM, RCIN,                      &
                       ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,         &
                       QM,     RA,   SWNET,  HLWDWN, PSUR,                     &
                       RDC,    HFTDS, DHFTDS,                                  &
                       QSATTC, DQSDTC, ALWRAD, BLWRAD,                         &
                       EMAXRT,CSOIL,SWSRF,POTFRC,BUG,WPWET,                    &
                       TC, QA,                                                 &
                       EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX, HSNACC           &
                         )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RCIN, ETURB, HSTURB, QM, RA,     &
                SWNET, HLWDWN, PSUR, RDC, HFTDS, DHFTDS, QSATTC, DQSDTC,       &
                ALWRAD, BLWRAD, EMAXRT, CSOIL, SWSRF, POTFRC, WPWET, DEDQA,    &
                DEDTC, DHSDQA, DHSDTC, POROS
      LOGICAL, INTENT(IN) ::  BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, RX1, RX2,      &
                GHFLUX, HSNACC


      INTEGER ChNo, N
      REAL, DIMENSION(NCH) :: DEDEA, DHSDEA, EM, ESATTC, DESDTC, EA, RC,       &
                DRCDTC, DRCDEA, SWSRF4
      REAL  DELTC, DELEA

!****
      DATA DELTC /0.01/, DELEA /0.001/
!****
!**** --------------------------------------------------------------------- 
!****
!**** Expand data as specified by ITYP into arrays of size NCH.
!****

!****
!**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

!      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
!      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
!      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
!      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

      IF (POROS(CHNO) < PEATCLSM_POROS_THRESHOLD) THEN
         ! mineral soil
         SWSRF4(CHNO) = SWSRF(CHNO)
      ELSE
         ! PEAT
         ! MB: For ET calculation, AR4 surface wetness is set to WPWET
         SWSRF4(CHNO) = WPWET(CHNO)
      ENDIF

 100  CONTINUE

!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
!****          ( VPDFAC  computes the vapor pressure deficit stress term,
!****           TMPFAC  computes the temperature stress term,
!****           RSURFP  computes rc given a parallel resist. from the
!****                   surface,
!****           RCANOP  computes rc corrected for snow and interception.)
!****

      DO N=1,NCH
        RC(N)=RCIN(N)
        ENDDO

      CALL RSURFP2 (                                                           &
                   NCH, UM, RDC, SWSRF4, ESATTC, EA, WPWET,                     &
                   RC,                                                         &
                   RX1, RX2                                                    &
                  )

      CALL RCANOP (                                                            &
                   NCH, RA, ETURB, POTFRC,                                     &
                   RC                                                          &
                  )

!****
!**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
!****
!**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
!****

      DO ChNo = 1, NCH
        DRCDTC(CHNO)=0.
        DRCDEA(CHNO)=0.
        ENDDO

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 2: Solve the energy balance at the surface.
!****
      CALL FLUXES (                                                            &
                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                     &
                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,             &
                       RC, DRCDEA, DRCDTC,SWNET, HLWDWN, ALWRAD, BLWRAD,       &
                       EM,  CSOIL,   PSUR, EMAXRT, HFTDS, DHFTDS,              &
                       TC,     EA,                                             &
                     EVAP, SHFLUX,  HLWUP, GHFLUX, HSNACC                      &
                   )

!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****
!****
!**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE
!****

      RETURN
      END SUBROUTINE energy4

!****
!**** [ END CHIP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!**** [ BEGIN FLUXES ]
!****
      SUBROUTINE FLUXES (                                                      &
                           NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                &
                          ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,       &
                          RC, DRCDEA, DRCDTC,SWNET, HLWDWN, ALWRAD, BLWRAD,    &
                             EM,  CSOIL,   PSUR, EMAXRT,                       &
                          HFTDS, DHFTDS,                                       &
                             TC,     EA,                                       &
                           EVAP, SHFLUX,  HLWUP, GHFLUX, HSNACC                &
                         )
!****
!**** This subroutine computes the fluxes of latent and sensible heat
!**** from the surface through an energy balance calculation.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: ESATTC, DESDTC, ETURB,  DEDEA,       &
                 DEDTC, HSTURB, DHSDEA, DHSDTC, RC, DRCDEA, DRCDTC, SWNET,     &
                 HLWDWN, ALWRAD, BLWRAD, EM,  CSOIL, PSUR, EMAXRT, HFTDS,      &
                 DHFTDS

      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, EA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, GHFLUX,        &
           HSNACC


      INTEGER ChNo
      REAL HLWTC,  CDEEPS, Q0,  RHOAIR,   CONST,  DHLWTC, EPLANT, A11, A12,    &
               A21, A22, F0, DEA, DTC, EANEW,  ESATNW,  EHARMN, DETERM, DENOM, &
               EDIF
      LOGICAL DEBUG, CHOKE


      DATA DEBUG /.FALSE./

!****
!**** -------------------------------------------------------------------

      HSNACC = 0.0  ! Initialize INTENT OUT variable before use -- LLT 

      DO 200 ChNo = 1, NCH
!****
      HLWTC = ALWRAD(CHNO) + BLWRAD(CHNO) * TC(CHNO)
      RHOAIR = PSUR(ChNo) * 100. / (RGAS * TC(ChNo))
      CONST = RHOAIR * EPSILON / PSUR(ChNo)
      DHLWTC = BLWRAD(CHNO)
!****
!**** Compute matrix elements A11, A22, AND Q0 (energy balance equation).
!****
      A11 = CSOIL(ChNo)/DTSTEP +                                               &
              DHLWTC +                                                         &
              DHSDTC(ChNo) +                                                   &
              ALHE*DEDTC(ChNo) +                                               &
              DHFTDS(CHNO)
      A12 = DHSDEA(ChNo) + ALHE * DEDEA(ChNo)
      Q0 =  SWNET(ChNo) +                                                      &
              HLWDWN(ChNo) -                                                   &
              HLWTC -                                                          &
              HSTURB(ChNo) -                                                   &
              ALHE * ETURB(ChNo) -                                             &
              HFTDS(CHNO)
!****
!**** Compute matrix elements A21, A22, and F0 (canopy water budget  
!**** equation) and solve for fluxes.  Three cases are considered:
!****
!**** 1. Standard case: RC>0.
!**** 2. RC = 0.  Can only occur if CIR is full or ETURB is negative.
!****
      CHOKE = .TRUE.

      IF( RC(CHNO) .GT. 0.) THEN 
          EPLANT = CONST * (ESATTC(ChNo) - EA(ChNo)) / RC(ChNo)
          IF(EPLANT*ETURB(ChNo).GT.0.) THEN
              EHARMN = 2.*EPLANT*ETURB(CHNO) / (EPLANT + ETURB(ChNo))
            ELSE
              EHARMN=0.
            ENDIF
!****
!****            Some limitations to A21 and A22 are applied:
!****            we assume that the increase in plant evaporation
!****            due to an increase in either TC or EA balances 
!****            or outweighs any decrease due to RC changes.
!****

          A21 =  -DEDTC(ChNo)*RC(ChNo) +                                       &
            amax1(0., CONST*DESDTC(ChNo) - EHARMN*DRCDTC(ChNo) )
          A22 = -( RC(ChNo)*DEDEA(ChNo) +                                      &
                     amax1( 0., CONST + EHARMN*DRCDEA(ChNo) )   )

          F0 = RC(ChNo) * (ETURB(ChNo) - EPLANT)
          DETERM = AMIN1( A12*A21/(A11*A22) - 1., -0.1 )
          DEA = ( Q0*A21 - A11*F0 ) / ( DETERM * A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
          EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
          SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
          DENOM = DETERM * A11*A22
        ELSE
          CHOKE = .FALSE.
          A21 = -DESDTC(ChNo)
          A22 = 1.
          F0 = ESATTC(ChNo) - EA(ChNo)
          DEA = ( Q0*A21 - A11*F0 ) / ( A12*A21 - A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
          EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
          SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
          DENOM = A12 * A21 - A11*A22
        ENDIF

!**** - - - - - - - - - - - - - - - - - - - - - - -
!**** Adjustments

!**** 1. Adjust deltas and fluxes if all available water evaporates
!****    during time step:
!**** NOTE: SOME CALCS BELOW ASSUME CROSS DERIVATIVE TERMS ARE ZERO
!****
      IF( EVAP(CHNO) .GT. EMAXRT(CHNO) ) THEN
        CHOKE = .FALSE.
        DEA=SIGN(ETURB(CHNO),EA(CHNO))
        IF(DEDEA(CHNO) .NE. 0.) DEA = (EMAXRT(CHNO)-ETURB(CHNO))/DEDEA(CHNO) 
        DEA = EM(CHNO) - EA(CHNO)
        DTC =  (Q0 + ALHE*(ETURB(ChNo)-EMAXRT(CHNO)) - DHSDEA(CHNO)*DEA)       &
                /  ( A11 - ALHE*DEDTC(ChNo) )
        EVAP(CHNO) = EMAXRT(CHNO)
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
        ENDIF


!**** 2. Pathological cases. 

!**** CHECK THAT EVAP AND ETURB ARE STILL CONSISTENT
!**** NOTE: SOME CALCS BELOW ASSUME CROSS DERIVATIVE TERMS ARE ZERO

      IF( EVAP(CHNO)*(ETURB(CHNO)+DEDEA(CHNO)*DEA) .LT. 0. )  THEN 
        CHOKE = .FALSE.
        DEA=SIGN(ETURB(CHNO),EA(CHNO))
        IF(DEDEA(CHNO) .NE. 0.) DEA = -ETURB(CHNO)/DEDEA(CHNO) 
        DTC = ( Q0 + ALHE*ETURB(ChNo) - DHSDEA(CHNO)*DEA ) /                   &
                  ( A11 - ALHE*DEDTC(ChNo) )
        EVAP(CHNO) = 0.
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
        ENDIF



!**** 3. Excessive dea change: apply "choke".
! -- (03.09.98) : changed to correct the conservation.

      IF( CHOKE .AND. ABS(DEA) .GT. 0.5*EA(CHNO) ) THEN
        DEA = SIGN(.5*EA(CHNO),DEA)
        DTC = ( Q0 - A12*DEA ) / A11
        EVAP(ChNo)   = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC

        IF(EVAP(CHNO) .GT. EMAXRT(CHNO)) THEN
          EDIF=EVAP(CHNO)-EMAXRT(CHNO)
          EVAP(CHNO) = EMAXRT(CHNO)
          SHFLUX(ChNo) = SHFLUX(CHNO)+EDIF*ALHE
!          HSNACC(ChNo) = HSNACC(CHNO)+EDIF*ALHE
          ENDIF

        ENDIF

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      TC(ChNo) = TC(ChNo) + DTC
      EA(ChNo) = EA(ChNo) + DEA
      HLWUP(CHNO) = HLWTC + DHLWTC*DTC

      HLWTC = ALWRAD(CHNO) + BLWRAD(CHNO) * TC(CHNO)
! warning: this ghflux is the real ground heat flux, and does not include
! the temperature variation
      GHFLUX(CHNO)=HFTDS(CHNO)+DHFTDS(CHNO)*DTC

!**** Make sure EA remains positive

      EA(CHNO) = AMAX1(EA(CHNO), 0.0)

  200 CONTINUE

      RETURN
      END SUBROUTINE FLUXES
!****
!**** [ END FLUXES ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN VPDFAC ]
!****
      SUBROUTINE VPDFAC (                                                      &
                         NCH, ITYP, ESATTC, EA,                                &
                         VPDSTR                                                &
                        )
!****
!**** This subroutine computes the vapor pressure deficit stress.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN):: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: ESATTC, EA
      REAL, INTENT(OUT), DIMENSION(NCH) :: VPDSTR

      INTEGER :: ChNo
      REAL, DIMENSION(NTYPS) :: VGDFAC
!****
      DATA VGDFAC /   .0273,    .0357,    .0310,    .0238,                     &
                      .0275,    .0275/
!****
!**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
!      VPDSTR(ChNo) = 1. - (ESATTC(ChNo)-EA(ChNo)) * VGDFAC(ITYP(ChNo))
!      VPDSTR (ChNo) = AMIN1( 1., AMAX1( VPDSTR(ChNo), 1.E-10 ) )
      VPDSTR(CHNO) = 1.
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE VPDFAC
!****
!**** [ END VPDFAC ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN TMPFAC ]
!****
      SUBROUTINE TMPFAC (                                                      &
                         NCH,  ITYP, TC,                                       &
                         FTEMP                                                 &
                        )
!****
!**** Compute temperature stress factor.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: TC
      REAL, INTENT(OUT), DIMENSION(NCH) :: FTEMP

      INTEGER ChNo, TypPtr
      INTEGER, PARAMETER :: MEMFAC = 5
      REAL, DIMENSION(MEMFAC*NTYPS) :: VGTLL, VGTU, VGTCF1, VGTCF2, VGTCF3

      DATA VGTLL /MemFac*273., MemFac*273., MemFac*268., MemFac*283.,          &
                  MemFac*283., MemFac*273./
      DATA VGTU /MemFac*318., MemFac*318., MemFac*313., MemFac*328.,           &
                 MemFac*323., MemFac*323./
      DATA VGTCF1 / MemFac*-1.43549E-06,  MemFac*-6.83584E-07,                 &
                    MemFac* 1.67699E-07,  MemFac*-1.43465E-06,                 &
                    MemFac*-2.76097E-06,  MemFac*-1.58094E-07/
      DATA VGTCF2 / MemFac* 7.95859E-04,  MemFac* 3.72064E-04,                 &
                    MemFac*-7.65944E-05,  MemFac* 8.24060E-04,                 &
                    MemFac* 1.57617E-03,  MemFac* 8.44847E-05/
      DATA VGTCF3 / MemFac*-1.11575E-01,  MemFac*-5.21533E-02,                 &
                    MemFac* 6.14960E-03,  MemFac*-1.19602E-01,                 &
                    MemFac*-2.26109E-01,  MemFac*-1.27272E-02/
!****
!**** ----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
      TypPtr = MOD(ChNo,MemFac) + (ITYP(ChNo)-1)*MemFac + 1
      FTEMP(ChNo) = (TC(ChNo) - VGTLL(TypPtr)) * (TC(ChNo) - VGTU(TypPtr)) *   &
                          ( VGTCF1(TypPtr)*TC(ChNo)*TC(ChNo) +                 &
                            VGTCF2(TypPtr)*TC(ChNo) + VGTCF3(TypPtr) )
      IF ( TC(ChNo) .LE. VGTLL(TypPtr) .OR. TC(ChNo) .GE. VGTU(TypPtr) )       &
            FTEMP (ChNo) = 1.E-10
      FTEMP(CHNO) = AMIN1( 1., AMAX1( FTEMP(ChNo), 1.E-10 ) )
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE TMPFAC
!****
!**** [ END TMPFAC ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN WUPDAT ]
!****
      SUBROUTINE WUPDAT (                                                      &
                           NCH, DTSTEP, EVAP, SATCAP, TC, RA, RC,              &
                           RX11,RX21,RX12,RX22,RX14,RX24, AR1,AR2,AR4,CDCR1,   &
                           EIRFRC,RZEQ,srfmn,WPWET,VGWMAX, POROS,              &
                           BF1, BF2, ARS1, ARS2, ARS3,                         &
                           CAPAC, RZEXC, CATDEF, SRFEXC,                       &
                           EINT, ESOI, EVEG                                    &
                          )
!****
!**** THIS SUBROUTINE ALLOWS EVAPOTRANSPIRATION TO ADJUST THE WATER
!**** CONTENTS OF THE INTERCEPTION RESERVOIR AND THE SOIL LAYERS.
!****
      IMPLICIT NONE
!****
      INTEGER, INTENT(IN) :: NCH

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: EVAP, SATCAP, TC, RA, RC, RX11,      &
             RX21, RX12, RX22, RX14, RX24, AR1, AR2, AR4, CDCR1, EIRFRC,       &
             RZEQ, srfmn, WPWET, VGWMAX, POROS, BF1, BF2, ARS1, ARS2, ARS3

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CAPAC, CATDEF, RZEXC, SRFEXC

      REAL, INTENT(OUT), DIMENSION(NCH) :: EINT, ESOI, EVEG


      INTEGER CHNO
      REAL EGRO, CNDSAT, CNDUNS, ESATFR, cndv, cnds, WILT, egromx,rzemax
      REAL :: ZBAR1,SYSOIL,ET_CATDEF,AR1eq

!****
!**** -----------------------------------------------------------------
      DO 100 CHNO = 1, NCH

!**** COMPUTE EFFECTIVE SURFACE CONDUCTANCES IN SATURATED AND UNSATURATED
!**** AREAS:

      CNDSAT=(AR1(CHNO)/RX11(CHNO)) + (AR1(CHNO)/RX21(CHNO))
      CNDUNS=(AR2(CHNO)/RX12(CHNO)) + (AR2(CHNO)/RX22(CHNO)) +                 &
             (AR4(CHNO)/RX14(CHNO)) + (AR4(CHNO)/RX24(CHNO))
     
      ESATFR=CNDSAT/(CNDSAT+CNDUNS)

!****
!**** PARTITION EVAP BETWEEN INTERCEPTION AND GROUND RESERVOIRS.
!****

      EINT(CHNO)=EIRFRC(CHNO)*EVAP(CHNO)*DTSTEP
      EGRO = EVAP(CHNO)*DTSTEP - EINT(CHNO)

!**** ENSURE THAT INDIVIDUAL CAPACITIES ARE NOT EXCEEDED.

      IF(EINT(CHNO) .GT. CAPAC(CHNO)) THEN
        EGRO=EGRO+EINT(CHNO)-CAPAC(CHNO)
        EINT(CHNO)=CAPAC(CHNO)
        ENDIF

! RK 09/16/03
      WILT=WPWET(CHNO)*VGWMAX(CHNO)
      rzemax=amax1(0.,RZEXC(CHNO)+RZEQ(CHNO)-WILT )
      egromx= rzemax + (srfexc(chno)-srfmn(chno))
      IF(EGRO .GT. egromx) THEN
! 06.02.98: the minimum is designed to prevent truncation errors 
        EINT(CHNO)=AMIN1(CAPAC(CHNO),EINT(CHNO)+EGRO-egromx)
        EGRO=egromx
        ENDIF
! RK 09/16/03
! RK: the above test ensures in particular that rzexc+rzeq never < 0.



      ESOI(CHNO)=AMIN1(SRFEXC(CHNO)-SRFMN(CHNO),                               &
              EGRO*(AR1(CHNO)*RX11(CHNO)/(RX11(CHNO)+RX21(CHNO)+1.E-20)        &
              +AR2(CHNO)*RX12(CHNO)/(RX12(CHNO)+RX22(CHNO)+1.E-20)             &
              +AR4(CHNO)*RX14(CHNO)/(RX14(CHNO)+RX24(CHNO)+1.E-20)))
      EVEG(CHNO)=EGRO-ESOI(CHNO)

!rdk   following inserted by koster Oct. 16, 2007
!**** special case if soil is at wilting point and rx14=rx24
!      if(rzemax .eq. 0.) then
!        esoi(chno)=AMIN1(SRFEXC(CHNO)-SRFMN(CHNO),EGRO)
!        EVEG(CHNO)=EGRO-ESOI(CHNO)
!        endif

      if(esoi(chno) .gt. SRFEXC(CHNO)-SRFMN(CHNO)) then
        esoi(chno)=SRFEXC(CHNO)-SRFMN(CHNO)
        EVEG(CHNO)=EGRO-ESOI(CHNO)
        endif

      if(eveg(chno) .gt. rzemax) then
        EVEG(CHNO)=rzemax
        esoi(chno)=egro-eveg(chno)
        endif


!****
!**** SPECIAL CASE FOR CONDENSATION:
      IF(EVAP(CHNO) .LT. 0.) THEN
        EINT(CHNO)=EVAP(CHNO)*DTSTEP
! 05.20.98: to prevent negative throughfall due to truncation errors 
        EINT(CHNO)=AMIN1(0.,EINT(CHNO))
        ESOI(CHNO)=0.
        EVEG(CHNO)=0.
        ENDIF

!****
!**** REMOVE MOISTURE FROM RESERVOIRS:
!****

        IF (CATDEF(CHNO) .LT. CDCR1(CHNO)) THEN
          CAPAC(CHNO) = AMAX1(0., CAPAC(CHNO) - EINT(CHNO))
          RZEXC(CHNO) = RZEXC(CHNO) - EVEG(CHNO)*(1.-ESATFR)
          SRFEXC(CHNO) = SRFEXC(CHNO) - ESOI(CHNO)*(1.-ESATFR)

          IF (POROS(CHNO) < PEATCLSM_POROS_THRESHOLD) THEN
             CATDEF(CHNO) = CATDEF(CHNO) + (ESOI(CHNO) + EVEG(CHNO))*ESATFR
          ELSE
             ! PEAT
             ! MB: accounting for water ponding on AR1
             ! same approach as for RZFLW (see subroutine RZDRAIN for
             ! comments)
             ZBAR1  = catch_calc_zbar( BF1(CHNO), BF2(CHNO), CATDEF(CHNO) )  
             SYSOIL = (2.*bf1(CHNO)*amin1(amax1(zbar1,0.),PEATCLSM_ZBARMAX_4_SYSOIL) + 2.*bf1(CHNO)*bf2(CHNO))/1000.
             SYSOIL = amin1(SYSOIL,poros(CHNO))
             ET_CATDEF = SYSOIL*(ESOI(CHNO) + EVEG(CHNO))*ESATFR/(1.*AR1(CHNO)+SYSOIL*(1.-AR1(CHNO)))
             AR1eq = (1.+ars1(chno)*(catdef(chno)))/(1.+ars2(chno)*(catdef(chno))+ars3(chno)*(catdef(chno))**2)
             CATDEF(CHNO) = CATDEF(CHNO) + (1.-AR1eq)*ET_CATDEF
          ENDIF
! 05.12.98: first attempt to include bedrock
        ELSE
          CAPAC(CHNO) = AMAX1(0., CAPAC(CHNO) - EINT(CHNO))
          RZEXC(CHNO) = RZEXC(CHNO) -  EVEG(CHNO)
          SRFEXC(CHNO) = SRFEXC(CHNO) - ESOI(CHNO)
        ENDIF

!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE WUPDAT
!****
!**** [ END WUPDAT ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RSURFP ]
!****
      SUBROUTINE RSURFP1 (                                                     &
                         NCH, UM, RDC, WET, ESATTC, EA,                        &
                         RC,                                                   &
                         RX1, RX2                                              &
                         )
!****
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RDC, WET, ESATTC, EA
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RC
      REAL, INTENT(OUT), DIMENSION(NCH) :: RX1, RX2

      INTEGER ChNo
      REAL  U2, RSURF, HESAT
!****
!**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
      U2 = UM(ChNo)
!      RSURF = RDC(ChNo) / U2 + 26. + 6. / (1.E-10 + WET(ChNo))**2
      RSURF = RDC(ChNo) / U2

      RX1(CHNO)=RC(CHNO)
      RX2(CHNO)=RSURF

      RC(ChNo) = RC(CHNO) * RSURF / ( RC(ChNo) + RSURF )
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE RSURFP1
!****
!**** [ END RSURFP ]
!****
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RSURFP ]
!****
      SUBROUTINE RSURFP2 (                                                     &
                         NCH, UM, RDC, WET, ESATTC, EA, WPWET,                 &
                         RC,                                                   &
                         RX1, RX2                                              &
                        )
!****
      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RDC, WET, ESATTC, EA, WPWET
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RC
      REAL, INTENT(OUT), DIMENSION(NCH) :: RX1, RX2


      INTEGER ChNo
!      REAL U2, RSURF, HESAT, RSWILT, RSSAT, ATERM, BTERM
      REAL U2, RSURF, HESAT, RSSAT, ATERM, BTERM

! RDK 04/04/06
!  VALUES OF BARE SOIL SURFACE RESISTANCE AT WILTING POINT, SATURATION

     PARAMETER (RSSAT=25.)

!****
!**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
      U2 = UM(ChNo)
      BTERM=(RSWILT-RSSAT) / (1./WPWET(CHNO)**2 -1.)
      ATERM=RSSAT-BTERM
      RSURF = RDC(ChNo) / U2 + ATERM + BTERM / (1.E-10 + WET(ChNo))**2

!**** Account for subsaturated humidity at soil surface:
!****
      HESAT = ESATTC(CHNO) * MIN( 1., WET(CHNO)*2. )
      IF( EA(CHNO) .LT. HESAT ) THEN
          RSURF=RSURF*( 1. + (ESATTC(CHNO)-HESAT)/(HESAT-EA(CHNO)) )
        ELSE
          RSURF=1.E10
        ENDIF


      RX1(CHNO)=RC(CHNO)
      RX2(CHNO)=RSURF

      RC(ChNo) = RC(CHNO) * RSURF / ( RC(ChNo) + RSURF )
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE RSURFP2
!****
!**** [ END RSURFP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RCANOP ]
!****
      SUBROUTINE RCANOP (                                                      &
                         NCH, RA, ETURB, POTFRC,                               &
                         RC                                                    &
                        )
!****
!**** The effective latent heat resistance RC depends on the quantity 
!**** of interception reservoir water and the snow cover.  POTFRC
!**** is the fraction of the tile from which potential evaporation
!**** occurs.
!****
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: RA, ETURB, POTFRC
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RC

      INTEGER N
!      REAL ETCRIT,RAMPFC

!**** (Note: ETCRIT arbitrarily set to ~-5 W/m2, or -2.e-6 mm/sec.)
!      DATA ETCRIT/ -2.E-6 /
!****
!**** -----------------------------------------------------------------

      DO N = 1, NCH

        RC(N)=RC(N)*(1.-POTFRC(N)) / ( 1.+POTFRC(N)*RC(N)/RA(N) )

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****   Assume RC=0 for condensation (dew).
!****   RAMPFC is used to ensure continuity in RC.

!_RDK   Remove zeroing of resistance for dew, due to stability problems
!        RAMPFC=ETURB(N)/ETCRIT
!        IF ( RAMPFC .GE. 0. ) RC(N) = RC(N)*(1.-RAMPFC)
!        IF ( RAMPFC .GT. 1. ) RC(N) = 0.
!****
        ENDDO

      RETURN
      END SUBROUTINE RCANOP
!****
!**** [ END RCANOP ]
!****

 ! *******************************************************************

  subroutine catch_calc_tsurf( NTILES, tc1, tc2, tc4, wesnn, htsnn,    &
       ar1, ar2, ar4, tsurf )
        
    ! Calculate diagnostic surface temperature "tsurf" from prognostics
    !
    ! reichle, Aug 31, 2004
    ! reichle, Jan  4, 2012 - optionally "ignore_snow"
    ! reichle, Apr  2, 2012 - revised for use without catch_types structures and
    !                          to avoid duplicate calls to rzequil() and partition()
    ! reichle, Oct 20, 2014 - removed option to "ignore_snow"; 
    !                          use subroutine catch_calc_tsurf_excl_snow() instead
    !
    ! ----------------------------------------------------------------
    
    implicit none
    
    integer,                           intent(in)           :: NTILES
    
    real,    dimension(       NTILES), intent(in)           :: tc1, tc2, tc4
    real,    dimension(N_snow,NTILES), intent(in)           :: wesnn, htsnn

    real,    dimension(       NTILES), intent(in)           :: ar1, ar2, ar4
    
    real,    dimension(       NTILES), intent(out)          :: tsurf
    
    ! ----------------------------
    !    
    ! local variables
    
    integer                    :: n
    
    real,    dimension(NTILES) :: asnow
    
    real                       :: tpsn1, real_dummy

    logical                    :: ice1, tzero

    logical, parameter         :: use_threshold_fac = .false.

    ! ------------------------------------------------------------------

    ! Compute tsurf excluding snow

    call catch_calc_tsurf_excl_snow( NTILES, tc1, tc2, tc4, ar1, ar2, ar4, tsurf )

    ! Compute snow covered area
    
    call StieglitzSnow_calc_asnow( N_snow, NTILES, wesnn, asnow )
    
    ! Add contribution of snow temperature 
    
    do n=1,NTILES
       
       if (asnow(n)>0.) then
                    
          ! StieglitzSnow_calc_tpsnow() returns snow temperature in deg Celsius
          
          call StieglitzSnow_calc_tpsnow( htsnn(1,n), wesnn(1,n), tpsn1, real_dummy,  &
               ice1, tzero, use_threshold_fac ) 
          
          tsurf(n) = (1. - asnow(n))*tsurf(n) + asnow(n)*(tpsn1 + TF)
          
       end if
       
    end do
        
  end subroutine catch_calc_tsurf

  ! *******************************************************************
   
  subroutine catch_calc_tsurf_excl_snow( NTILES, tc1, tc2, tc4, ar1, ar2, ar4, &
       tsurf_excl_snow )
    
    ! Calculate diagnostic surface temperature "tsurf" ignoring snow
    !
    ! reichle, 20 Oct 2014
    !
    ! ----------------------------------------------------------------
    
    implicit none
    
    integer,                           intent(in)           :: NTILES
    real,    dimension(       NTILES), intent(in)           :: tc1, tc2, tc4
    real,    dimension(       NTILES), intent(in)           :: ar1, ar2, ar4    
    real,    dimension(       NTILES), intent(out)          :: tsurf_excl_snow
        
    ! ------------------------------------------------------------------

    tsurf_excl_snow = ar1*tc1 + ar2*tc2 + ar4*tc4
    
  end subroutine catch_calc_tsurf_excl_snow
  
  ! *******************************************************************
 
  subroutine catch_calc_etotl( NTILES, vegcls, dzsf, vgwmax, cdcr1, cdcr2, &
       psis, bee, poros, wpwet, bf1, bf2,                                  &
       ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1, arw2, arw3, arw4,   &
       srfexc, rzexc, catdef, tc1, tc2, tc4, wesnn, htsnn, ghtcnt,         &
       etotl )
    
    ! compute total energy stored in land tiles
    !
    ! reichle,  4 Jan 2012
    ! reichle,  2 Apr 2012 - revised for use without catch_types structures
    !
    ! ----------------------------------------------------------------
    
    implicit none
    
    integer,                           intent(in)  :: NTILES
    
    integer, dimension(       NTILES), intent(in)  :: vegcls
    real,    dimension(       NTILES), intent(in)  :: dzsf
    real,    dimension(       NTILES), intent(in)  :: vgwmax
    real,    dimension(       NTILES), intent(in)  :: cdcr1, cdcr2, bf1, bf2
    real,    dimension(       NTILES), intent(in)  :: psis, bee, poros, wpwet    
    real,    dimension(       NTILES), intent(in)  :: ars1, ars2, ars3
    real,    dimension(       NTILES), intent(in)  :: ara1, ara2, ara3, ara4
    real,    dimension(       NTILES), intent(in)  :: arw1, arw2, arw3, arw4
    real,    dimension(       NTILES), intent(in)  :: srfexc, rzexc, catdef
    real,    dimension(       NTILES), intent(in)  :: tc1, tc2, tc4
    real,    dimension(N_snow,NTILES), intent(in)  :: wesnn, htsnn
    real,    dimension(N_gt  ,NTILES), intent(in)  :: ghtcnt
    
    real,    dimension(       NTILES), intent(out) :: etotl
    
    ! ----------------------------
    !    
    ! local variables
    
    integer :: n

    real    :: tot_htsn, tot_ght, csoil
    
    real, dimension(NTILES) :: srfexc_tmp, rzexc_tmp, catdef_tmp
        
    real, dimension(NTILES) :: ar1, ar2, ar4, avg_tc

    ! ----------------------------------------------------------------
    !
    ! diagnose ar1, ar2, ar4 prior to catch_calc_tsurf()
    
    srfexc_tmp = srfexc   ! srfexc is "inout" in catch_calc_soil_moist()
    rzexc_tmp  = rzexc    ! rzexc  is "inout" in catch_calc_soil_moist()
    catdef_tmp = catdef   ! catdef is "inout" in catch_calc_soil_moist()
    
    call catch_calc_soil_moist(                                                       &
         NTILES, dzsf, vgwmax, cdcr1, cdcr2, psis, bee, poros, wpwet,                 &
         ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1, arw2, arw3, arw4, bf1, bf2,  &
         srfexc_tmp, rzexc_tmp, catdef_tmp, ar1, ar2, ar4 )
    
    ! compute snow-free tsurf
   
    call catch_calc_tsurf_excl_snow(                                           &
         NTILES, tc1, tc2, tc4, ar1, ar2, ar4, avg_tc ) 
    
    do n=1,NTILES
       
       ! total snow heat content
       
       tot_htsn = sum( htsnn(1:N_snow,n) )
       
       ! total ground heat content
       
       tot_ght  = sum( ghtcnt(1:N_gt,n))
       
       ! total energy
       
       if (vegcls(n)==1) then
          csoil = CSOIL_1
       else
          csoil = CSOIL_2
       end if
       
       etotl(n) = csoil*avg_tc(n) + tot_htsn + tot_ght 
       
    end do
    
  end subroutine catch_calc_etotl

  ! ********************************************************************

END MODULE CATCHMENT_MODEL
