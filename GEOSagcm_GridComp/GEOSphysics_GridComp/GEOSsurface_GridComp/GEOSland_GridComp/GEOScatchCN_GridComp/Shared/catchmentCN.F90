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
! koster, 2010-2013    - unified model, latest revision rc3f with matrix_calc fix
! walker, 2010-2015    - updated to unified land model with carbon physics (koster version xxxxx)
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
! Justin, 20 Jun 2018  - moved CSOIL_2 to SURFPARAMS
! Sarith , 20 Apr 2020  - introducing USE_FWET_FOR_RUNOFF and passing FWETL and FWETC via GEOS_SurfaceGridComp.rc

MODULE CATCHMENT_CN_MODEL
  
  USE MAPL_BaseMod,      ONLY: MAPL_Land
  
  USE MAPL_ConstantsMod, ONLY:           &
       PIE               => MAPL_PI,     &  ! -                       
       ALHE              => MAPL_ALHL,   &  ! J/kg  @15C              
       ALHM              => MAPL_ALHF,   &  ! J/kg                    
       ALHS              => MAPL_ALHS,   &  ! J/kg                    
       TF                => MAPL_TICE,   &  ! K                       
       RGAS              => MAPL_RGAS,   &  ! J/(kg K)                
       CPAIR             => MAPL_CP,     &  ! J/(kg K)
       EPSILON           => MAPL_EPSILON 
  
  USE CATCH_CONSTANTS,   ONLY:                   &
       N_SNOW            => CATCH_N_SNOW,        &
       N_GT              => CATCH_N_GT,          &
       CATCH_SNOW_RHOFS,                         &
       CATCH_SNOW_MAXDEPTH,                      &
       CATCH_SNOW_DZPARAM,                       &
       C_CANOP           => CATCH_C_CANOP,       &
       N_sm              => CATCH_N_ZONES,       &
       SATCAPFR          => CATCH_SATCAPFR,      &
       PHIGT             => CATCH_PHIGT,         &
       DZTSURF           => CATCH_DZTSURF,       &
       DZGT              => CATCH_DZGT,          &
       FSN               => CATCH_FSN,           &
       PEATCLSM_POROS_THRESHOLD,                 &
       PEATCLSM_ZBARMAX_4_SYSOIL

  USE SURFPARAMS,       ONLY: CSOIL_2, RSWILT,   &
       LAND_FIX, FLWALPHA
  
  USE lsm_routines, only :                       &
       INTERC, RZDRAIN, BASE, PARTITION, RZEQUIL,&
       gndtp0, gndtmp,                           &   
       catch_calc_soil_moist, catch_calc_zbar,   &
       catch_calc_wtotl, dampen_tc_oscillations, &
       SRUNOFF
  
  USE SIBALB_COEFF,  ONLY: coeffsib
  
  USE STIEGLITZSNOW, ONLY: &
       StieglitzSnow_snowrt,                     &
       StieglitzSnow_calc_asnow,                 &
       StieglitzSnow_calc_tpsnow,                &
       N_constit

  
  IMPLICIT NONE
  
  private
  
  public :: catchcn

  ! -----------------------------------------------------------------------------
  
  ! moved all "private" constants to here from "catch_constants.f90"
  ! - reichle, 14 Aug 2014
  
  REAL,    PARAMETER :: ZERO     = 0.
  REAL,    PARAMETER :: ONE      = 1.
  
  REAL,    PARAMETER :: DTCCRIT  = 8.
  REAL,    PARAMETER :: DTGCRIT  = 3.
  REAL,    PARAMETER :: DEACRIT  = 10.  
  
CONTAINS
  
  !****
  !
  ! This code represents a merge of the catchment model (as of April 28, 2009)
  ! and the most recent version of the "unified" model (from Sept. 20, 2006).
  
  SUBROUTINE CATCHCN (                                           &
       NCH, LONS, LATS, DTSTEP, UFW4RO, FWETC, FWETL, cat_id,    &
       ITYP1,ITYP2,FVEG1,FVEG2,                                  &
       DZSF, TRAINC,TRAINL, TSNOW, TICE, TFRZR, UM,              &
       ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,         &
       ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,         &
       ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,         &
       ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,         &
       TM, QM, ra1, ra2, ra4, raS, SUNANG,                       &
       SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,             &  ! HLWDWN = *absorbed* longwave only (excl reflected)
       SQSCAT, RSOIL1, RSOIL2,   RDC,                            &
       QSAT1, DQS1, ALW1, BLW1,  QSAT2, DQS2, ALW2, BLW2,        &
       QSAT4, DQS4, ALW4, BLW4,  QSATS, DQSS, ALWS, BLWS,        &
       RCSAT,DRCSDT,DRCSDQ,  RCUNS,DRCUDT,DRCUDQ,                &
       BF1, BF2, BF3,VGWMAX,                                     &
       CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,          &
       ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,   &
       tsa1,tsa2,tsb1,tsb2,atau,btau,BUG,                        &
       TG1, TG2, TG4,                                            &
       TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,                      &
       CATDEF, RZEXC, srfexc, GHTCNT,                            &
       WESNN, HTSNNN, SNDZN,     EVAP, SHFLUX, RUNOFF,           &
       EINT, ESOI, EVEG, ESNO,  BFLOW,RUNSRF,SMELT,              &
       HLWUP,SWLAND,HLATN,QINFIL,AR1, AR2, RZEQ,                 &  ! HLWUP = *emitted* longwave only (excl reflected)
       GHFLUX, GHFLUXSNO, GHTSKIN, TPSN1, ASNOW0,                &
       TP1, TP2, TP3, TP4, TP5, TP6,                             &
       sfmc, rzmc, prmc, entot, wtot, WCHANGE, ECHANGE, HSNACC,  &
       EVACC, SHACC, TSURF,                                      &
       SH_SNOW, AVET_SNOW, WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL,   &
       LH_SNOW, LWUP_SNOW, LWDOWN_SNOW, NETSW_SNOW,              &
       TCSORIG, TPSN1IN, TPSN1OUT, FSW_CHANGE, FICESOUT,         &
       TC1_0, TC2_0, TC4_0, QA1_0, QA2_0, QA4_0, EACC_0,         &  ! OPTIONAL
       RCONSTIT, RMELT, TOTDEPOS, LHACC                          &  ! OPTIONAL
       )
    
    IMPLICIT NONE
    
    ! -----------------------------------------------------------
    !     INPUTS
    
    INTEGER, INTENT(IN) :: NCH
    INTEGER, INTENT(IN), DIMENSION(:) :: ITYP1, ITYP2, cat_id
    
    REAL, INTENT(IN)    :: DTSTEP, FWETC, FWETL
    LOGICAL, INTENT(IN) :: UFW4RO
    REAL, INTENT(IN), DIMENSION(:) :: DZSF, TRAINC, TRAINL, TSNOW, &
         TICE, TFRZR, UM, FVEG1,   FVEG2,                          &
         ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,         &
         ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,         &
         ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,         &
         ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,         &
         TM, QM, ra1, ra2, ra4, raS, SUNANG,                       &
         SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,             &
         SQSCAT, RSOIL1, RSOIL2,   RDC,                            &
         QSAT1, DQS1, ALW1, BLW1,  QSAT2, DQS2, ALW2, BLW2,        &
         QSAT4, DQS4, ALW4, BLW4,  QSATS, DQSS, ALWS, BLWS,        &
         BF1, BF2, BF3,VGWMAX,                                     &
         CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,          &
         ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,   &
         tsa1,tsa2,tsb1,tsb2,atau,btau,                            &
         RCSAT,DRCSDT,DRCSDQ, RCUNS,DRCUDT,DRCUDQ
    REAL, INTENT(IN), DIMENSION(:) :: LONS, LATS

    REAL, INTENT(IN), DIMENSION(:, :), OPTIONAL :: TOTDEPOS
    
    LOGICAL, INTENT(IN) :: BUG
    
    ! -----------------------------------------------------------
    !     PROGNOSTIC VARIABLES
    
    REAL, INTENT(INOUT), DIMENSION(:) ::                         &
         TG1, TG2, TG4, TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,       &
         CATDEF, RZEXC, SRFEXC
 
    REAL, INTENT(INOUT), DIMENSION(:,:) ::  GHTCNT
    
    REAL, INTENT(INOUT), DIMENSION(:, :) :: WESNN, HTSNNN, SNDZN
    
    REAL, INTENT(INOUT), DIMENSION(:, :, :),                  &
               OPTIONAL :: RCONSTIT
    
    ! -----------------------------------------------------------
    !     DIAGNOSTIC OUTPUT VARIABLES
    
    REAL, INTENT(OUT), DIMENSION(:) ::      EVAP, SHFLUX, RUNOFF,          &
         EINT, ESOI, EVEG, ESNO, BFLOW,RUNSRF,SMELT,               &
         HLWUP,SWLAND,HLATN,QINFIL,AR1, AR2, RZEQ,                 &
         GHFLUX, TPSN1, ASNOW0, TP1, TP2, TP3, TP4, TP5, TP6,      &
         sfmc, rzmc, prmc, entot, wtot, tsurf, WCHANGE, ECHANGE,   &
         HSNACC, EVACC, SHACC
    REAL, INTENT(OUT), DIMENSION(:) :: GHFLUXSNO, GHTSKIN
    
    REAL, INTENT(OUT), DIMENSION(:) :: SH_SNOW, AVET_SNOW,         &
         WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL
    REAL, INTENT(OUT), DIMENSION(:) :: LH_SNOW, LWUP_SNOW,         &
         LWDOWN_SNOW, NETSW_SNOW
    REAL, INTENT(OUT), DIMENSION(:) :: TCSORIG, TPSN1IN, TPSN1OUT, &
         FSW_CHANGE
    
    REAL, INTENT(OUT), DIMENSION(:, :) :: FICESOUT
    
    REAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: LHACC
    
    REAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: TC1_0,TC2_0,TC4_0
    REAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: QA1_0,QA2_0,QA4_0 	
    REAL, INTENT(OUT), DIMENSION(:), OPTIONAL :: EACC_0
    REAL, INTENT(OUT), DIMENSION(:, :), OPTIONAL :: RMELT	
    
    ! -----------------------------------------------------------
    !     LOCAL VARIABLES
    
    INTEGER I,K,N,LAYER
    
    REAL, DIMENSION(NCH) :: CSOIL, CCANOP, ASNOW, traincx, trainlx,       &
         RC, SATCAP, SNWFRC, POTFRC,  ESNFRC, EVSNOW, SHFLUXS, HLWUPS,    &
         HFTDS1, HFTDS2, HFTDS4, DHFT1, DHFT2, DHFT4, TPSNB,              &
         QSATTC, DQSDTC, SWSRF1, SWSRF2, SWSRF4, AR4,                     &
         FCAN, THRUL_VOL, THRUC_VOL, RZEQOL, frice, srfmx,                &
         srfmn, RCST1, RCST2, EVAPFR, RDCX, EVAP1, EVAP2,                 &
         EVAP4, SHFLUX1, SHFLUX2, SHFLUX4, HLWUP1, HLWUP2, HLWUP4,        &
         GHFLUX1, GHFLUX2, GHFLUX4, RZI, TC1SF, TC2SF, TC4SF, ar1old,     &
         ar2old, ar4old, GHFLUXS, DEDQA1X, DEDTC1X,                       &
         DHSDQA1X, DHSDTC1X, DEDQA2X, DEDTC2X, DHSDQA2X, DHSDTC2X,        &
         DEDQA4X, DEDTC4X, DHSDQA4X, DHSDTC4X, werror, sfmcun, rzmcun,    &
         prmcun,WTOT_ORIG,ENTOT_ORIG,                                     &
         TC1_00, TC2_00, TC4_00, EACC_00,                                 &
         qa1_orig,qa2_orig,qa4_orig,tc1_orig,tc2_orig,tc4_orig,           &
         tgs_orig,TG1SF,TG2SF,TG4SF,RCUN1,RCUN2,                          &
         tg1_orig,tg2_orig,tg4_orig,                                      &
         EVROOT1, EVROOT2, EVROOT4, EVSURF1, EVSURF2, EVSURF4,            &
         EVINT1, EVINT2, EVINT4, ESATFR, ECORR, DRCST1DT, DRCST1DQ,       &
         DRCST2DT, DRCST2DQ, FVEG, RD, RCST, DRCSTDT, DRCSTDQ, RSURF
    
    
    REAL, DIMENSION(N_gt) :: HT, TP, soilice
    
    REAL, DIMENSION(N_SNOW) :: TPSN, WESN, HTSNN, SNDZ, fices,            &
         wesnperc,wesndens,wesnrepar,excs,drho0,tksno, tmpvec_Nsnow
    
    REAL, DIMENSION(N_SNOW, N_Constit) :: RCONSTIT1
    REAL, DIMENSION(N_Constit)         :: RMELT1, TOTDEP1
    
    REAL, DIMENSION(N_SM) :: T1, AREA, tkgnd, fhgnd
    
    REAL :: TG1SN, TG2SN, TG4SN, DTG1SN,DTG2SN,DTG4SN, ZBAR, THETAF,      &
         XFICE, FH21, FH21W, FH21I, FH21D, DFH21W, DFH21I, DFH21D,        &
         EVSN, SHFLS, HUPS, HCORR, SWNET0, HLWDWN0, TMPSNW, HLWTC,        &
         DHLWTC, HSTURB, DHSDEA, DHSDTC, ESATTC, ETURB, DEDEA, DEDTC,     &
         SNOWF, TS, fh31w, fh31i, fh31d, pr, ea, desdtc, areasc,          &
         pre, dummy1, dummy2, dummy3, areasc0, EDIF, EINTX,               &
         SCLAI, tsn1, tsn2, tsn3, hold, hnew, dedtc0,                     &
         dhsdtc0, alhfsn, ADJ, raddn, zc1, tsnowsrf, dum, tsoil,          &
         QA1X, QA2X, QA4X, TC1X, TC2X, TC4X, TCSX,                        &
         EVAPX1,EVAPX2,EVAPX4,SHFLUXX1,SHFLUXX2,SHFLUXX4,EVEGFRC,         &
         EVAPXS,SHFLUXXS,DTC1SN,DTC2SN,DTC4SN,TCANOP,                     &
         ZLAI0, phi,rho_fs,WSS,sumdepth,                                  &   
         sndzsc, wesnprec, sndzprec,  sndz1perc,                          &   
         mltwtr, wesnbot, dtss
    
    
    
    LOGICAL :: ldum

    REAL    :: dtc1, dtc2, dtc4
    
    real, parameter :: zlai_min=.001  !  Assume there's some interception 
    !  (if only ponding or on 
    !  dead vegetation) at all times.
    REAL, PARAMETER :: RSSAT=25.       !  Resistance to bare soil evaporation
    !  under saturated conditions
    
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
    !     do i = 1,nch
    !      if ((lons(i).gt.2.105).and.(lons(i).lt.2.106).and.(lats(i).gt.1.134).and.(lats(i).lt.1.135)) then
    !       numout = numout + 1
    !       n_outs(numout) = i
    !       write (*,*) 'found point ',n_out,' ',lons(i),' ',lats(i)
    !      endif
    !     enddo
    
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
          write (*,*) ITYP1(n_out)  
          write (*,*) ITYP2(n_out)  
          write (*,*) FVEG1(n_out)  
          write (*,*) FVEG2(n_out)  
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
          write (*,*)  SWNETF(n_out)    
          write (*,*)  SWNETS(n_out)    
          write (*,*)   HLWDWN(n_out)    
          write (*,*)  PSUR(n_out)    
          write (*,*)   ZLAI(n_out)    
          write (*,*)    GREEN(n_out)    
          write (*,*)  SQSCAT(n_out)    
          write (*,*)  RSOIL1(n_out)    
          write (*,*)  RSOIL2(n_out)    
          write (*,*)    RDC(n_out)    
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
          write (*,*)  RCSAT(n_out)
          write (*,*)  DRCSDT(n_out)
          write (*,*)  DRCSDQ(n_out)
          write (*,*)  RCUNS(n_out)
          write (*,*)  DRCUDT(n_out)
          write (*,*)  DRCUDQ(n_out)
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
         write (*,*)  TG1(n_out)    
         write (*,*)  TG2(n_out)    
         write (*,*)  TG4(n_out)    
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
       LH_SNOW(N)=0.
       SH_SNOW(N)=0.
       LWUP_SNOW(N)=0.
       LWDOWN_SNOW(N)=0.
       NETSW_SNOW(N)=0.
      end do

!**** ---------------------------------------------------
!**** PRE-PROCESS DATA AS NECESSARY:
!****

      DO N=1,NCH
!        ZLAI0 = AMAX1(ZLAI(N), ZLAI_MIN)
        ZLAI0 = ZLAI(N)
        SATCAP(N) = SATCAPFR * ZLAI0 + 1.e-5
        CSOIL(N)  = CSOIL_2
        CCANOP(N) = C_CANOP
!        if ( ityp1(n) .ne. 1) CSOIL(N)  = CSOIL_2
        FCAN(N) = AMIN1( 1., AMAX1(0.,CAPAC(N)/SATCAP(N)) )
!gkw    SCLAI=amin1( 1., ZLAI0/2. )
        SCLAI=amax1(0.001,  amin1( 1., ZLAI0/2. ))
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
        tg1_orig(n)=tg1(n)
        tg2_orig(n)=tg2(n)
        tg4_orig(n)=tg4(n)
        tc1_orig(n)=tc1(n)
        tc2_orig(n)=tc2(n)
        tc4_orig(n)=tc4(n)

        IF (LAND_FIX) THEN
	   ! Andrea Molod (Oct 21, 2016):
           qa1(n) = min(max(qm(N),qsat1(N)),qa1(N))
           qa1(n) = max(min(qm(N),qsat1(N)),qa1(N))
           qa2(n) = min(max(qm(N),qsat2(N)),qa2(N))
           qa2(n) = max(min(qm(N),qsat2(N)),qa2(N))
           qa4(n) = min(max(qm(N),qsat4(N)),qa4(N))
           qa4(n) = max(min(qm(N),qsat4(N)),qa4(N))
        END IF
!        if(ityp1(n) .ge. 7) potfrc(n)=0.

!     HSNACC is an energy accounting term designed to account (among other,
!     lesser things) for the fact that snow is deposited at the snowpack
!     surface temperature while the atmosphere does not account for variations
!     in the heat content of deposited snow.
 
        HSNACC(N)=0.
        EVACC(N)=0.
        SHACC(N)=0.
        RUNSRF(N)=0.


!! !****   RESET LAND ICE VARIABLES, MAINTAINING TEMPS. AT EACH LAYER
!!        IF(ITYP1(N) .EQ. 9) THEN
!!
!!          ! This block of the code should no longer be used.
!!          ! If it is, Randy wants to know about it.
!!          ! reichle+koster, 12 Aug 2014
!!          write (*,*) 'catchment() encountered ityp==9. STOPPING.'
!!          stop 
!!
!!          if(sum(htsnnn(:,n)+wesnn(:,n))==0.) then
!!              TSN1=tc1(n)-TF
!!              TSN2=tc1(n)-TF
!!              TSN3=tc1(n)-TF
!!            else
!!              TSN1=(HTSNNN(1,N)+WESNN(1,N)*ALHM)/(SCONST*WESNN(1,N)+1.e-5)
!!              TSN2=(HTSNNN(2,N)+WESNN(2,N)*ALHM)/(SCONST*WESNN(2,N)+1.e-5)
!!              TSN3=(HTSNNN(3,N)+WESNN(3,N)*ALHM)/(SCONST*WESNN(3,N)+1.e-5)
!!            endif
!!          WESNN(1,N)=.1
!!          WESNN(2,N)=.2
!!          WESNN(3,N)=.1
!!          HTSNNN(1,N)=-ALHM*WESNN(1,N)+TSN1*SCONST*WESNN(1,N)
!!          HTSNNN(2,N)=-ALHM*WESNN(2,N)+TSN1*SCONST*WESNN(2,N)
!!          HTSNNN(3,N)=-ALHM*WESNN(3,N)+TSN1*SCONST*WESNN(3,N)
!!          SNDZN(1,N)=WESNN(1,N)/.9
!!          SNDZN(2,N)=WESNN(2,N)/.9
!!          SNDZN(3,N)=WESNN(3,N)/.9
!!          POTFRC(N)=1.
!!
!!          ENDIF

!****   RESET LAKE VARIABLES
        IF(ITYP1(N) .EQ. 10) THEN
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
         TCANOP=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)
         TSOIL=AR1(N)*TG1(N)+AR2(N)*TG2(N)+AR4(N)*TG4(N)

         ENTOT_ORIG(N) =                                                       &
            sum(HTSNNN(1:N_snow,N)) + TSOIL*CSOIL(N) + TCANOP*CCANOP(N) + sum(GHTCNT(1:N_gt,N))

      ENDDO

      ! reichle, 1 May 2013, fix TWLAND<0 bug, use correct calculation in 
      
      CALL CATCH_CALC_WTOTL( NCH,                                              &
           CDCR2, WPWET, SRFEXC, RZEXC, CATDEF, CAPAC, WESNN,                  &
           WTOT_ORIG )

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
        T1(1)=TG1(N)
        T1(2)=TG2(N)
        T1(3)=TG4(N)
        if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
           phi=POROS(N)
        else
           phi=PHIGT
        end if
        ! zbar function - reichle, 29 Jan 2022 (minus sign applied in call to GNDTP0)
        ZBAR = catch_calc_zbar( bf1(n), bf2(n), catdef(n) )  
        THETAF=.5
        DO LAYER=1,6
          HT(LAYER)=GHTCNT(LAYER,N)
          ENDDO

        CALL GNDTP0(                                                &
                    T1,phi,-1.*ZBAR,THETAF,                         &   ! note minus sign for zbar
                    HT,                                             &
                    fh21w,fH21i,fh21d,dfh21w,dfh21i,dfh21D,tp       &
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
        TG1SF(N)  = TG1(N)
        TG2SF(N)  = TG2(N)
        TG4SF(N)  = TG4(N)
        TC1SF(N)  = TC1(N)
        TC2SF(N)  = TC2(N)
        TC4SF(N)  = TC4(N)
        FVEG(N)   = FVEG1(N) + FVEG2(N)
        RD(N) = RDCX(N) / UM(N)  !  UM sent down as U2
        rcst(n) = 1.E10
        DRCSTDT(N) = 0.
        DRCSTDQ(N) = 0.

        ENDDO

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!**** 1. SATURATED FRACTION

      DO N=1,NCH
        RSURF(N)=RSSAT
        ENDDO

      CALL FLUXES (                                                            &
                     NCH,   FVEG, DTSTEP, QSAT1, DQS1,                         &
                     ETURB1,  DEDQA1X,  DEDTC1X, HSTURB1, DHSDQA1X, DHSDTC1X,  &
                     RCSAT, DRCSDQ, DRCSDT,                                    &
                     SWNETF, HLWDWN, ALW1, BLW1,                               &
                     QM,  CSOIL, CCANOP,   PSUR,                               &
                     HFTDS1, DHFT1, RD, RSURF, POTFRC,                         &
                     TG1SF, TC1SF, QA1,                                        &
                     EVAP1, EVROOT1, EVSURF1, EVINT1,                          &
                     SHFLUX1,  HLWUP1, GHFLUX1                                 &
                   )

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!**** 2. SUBSATURATED BUT UNSTRESSED FRACTION

      CALL RSURFP2 ( NCH, RSSAT, SWSRF2, QSAT2, QA2, PSUR, WPWET, RSURF )

      CALL FLUXES (                                                            &
                     NCH,   FVEG, DTSTEP, QSAT2, DQS2,                         &
                     ETURB2,  DEDQA2X,  DEDTC2X, HSTURB2, DHSDQA2X, DHSDTC2X,  &
                     RCUNS, DRCUDQ, DRCUDT,                                    &
                     SWNETF, HLWDWN, ALW2, BLW2,                               &
                     QM,  CSOIL, CCANOP,   PSUR,                               &
                     HFTDS2, DHFT2, RD, RSURF, POTFRC,                         &
                     TG2SF,  TC2SF,     QA2,                                   &
                     EVAP2, EVROOT2, EVSURF2, EVINT2,                          &
                     SHFLUX2,  HLWUP2, GHFLUX2                                 &
                  )

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!**** 3. WILTING FRACTION

      ! MB: For ET calculation, AR4 surface wetness is set to WPWET 
      WHERE (POROS > PEATCLSM_POROS_THRESHOLD)  SWSRF4 = WPWET  ! PEATCLSM
      
      CALL RSURFP2 ( NCH, RSSAT, SWSRF4, QSAT4, QA4, PSUR, WPWET, RSURF  )

      CALL FLUXES (                                                            &
                     NCH,   FVEG, DTSTEP, QSAT4, DQS4,                         &
                     ETURB4,  DEDQA4X,  DEDTC4X, HSTURB4, DHSDQA4X, DHSDTC4X,  &
                     RCST, DRCSTDQ, DRCSTDT,                                   &
                     SWNETF, HLWDWN, ALW4, BLW4,                               &
                     QM,  CSOIL, CCANOP,   PSUR,                               &
                     HFTDS4, DHFT4, RD, RSURF, POTFRC,                         &
                     TG4SF, TC4SF, QA4,                                        &
                     EVAP4, EVROOT4, EVSURF4, EVINT4,                          &
                     SHFLUX4,  HLWUP4, GHFLUX4                                 &
                   )

!**** --------------------------------------------------------
!**** B. SNOW-COVERED FRACTION.

      DO N=1,NCH

        WSS    = UM(N)
        TS     = TM(N) 
        T1(1)  = TG1(N)-TF
        T1(2)  = TG2(N)-TF
        T1(3)  = TG4(N)-TF
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
 		   RMELT1(K) = 0.
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
        tgs_orig(n)=tsnowsrf+tf
        if(wesn(1)+wesn(2)+wesn(3) .eq. 0.) tgs_orig(n)=                       &
                  amin1( tf, tg1_orig(n)*ar1(n)+tg2_orig(n)*ar2(n)+            &
                  tg4_orig(n)*(1.-ar1(n)-ar2(n)) )

        hlwtc=ALWS(N) + BLWS(N)*(TSNOWSRF+TF) 
        dhlwtc=BLWS(N)
        hcorr=0.

!! the field tpsn1(n) contains the value of TC(snow tile) that
!! came in from the catch grid comp, and catch internal state
!! spit it out here as a diagnostic
!! also spit out here the value of tsnowsrf+tf which is used
!! as the "original" value of TC for purposes of LW and turb fluxes

        tcsorig(N) = tgs_orig(n)
        tpsn1in(n) = tpsn1(n)    ! tpsn1 is "intent(out)", should NOT be used here, use catch_calc_tpsnow instead?  shouldn't this be the same as tgs_orig?  - reichle, 8/8/2014

        sumdepth=sum(sndz)

        CALL StieglitzSnow_snowrt(                                             &
                   N_sm, N_snow, MAPL_Land,                                    &
                   t1,area,tkgnd,pr,snowf,ts,DTSTEP,                           &
                   eturbs(n),dedtc0,hsturb,dhsdtc0,hlwtc,dhlwtc,               &
                   desdtc,hups,raddn,zc1, totdep1, wss,                        &
                   wesn,htsnn,sndz,fices,tpsn,RCONSTIT1, RMELT1,               &
                   areasc,areasc0,pre,fhgnd,                                   &
                   EVSN,SHFLS,alhfsn,hcorr, ghfluxsno(n),                      &
                   sndzsc, wesnprec, sndzprec,  sndz1perc,                     &   
                   wesnperc, wesndens, wesnrepar, mltwtr,                      &
                   excs, drho0, wesnbot, tksno, dtss,                          &
                   CATCH_SNOW_MAXDEPTH, CATCH_SNOW_RHOFS, CATCH_SNOW_DZPARAM )

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
        HSNACC(N) = HSNACC(N) + hcorr
 
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

        DTG1SN=((-(FH31W/(area(1)+1.e-20))-HFTDS1(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT1(N)*DTSTEP/CSOIL(N))
        DTG2SN=((-(FH31I/(area(2)+1.e-20))-HFTDS2(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT2(N)*DTSTEP/CSOIL(N))
        DTG4SN=((-(FH31D/(area(3)+1.e-20))-HFTDS4(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT4(N)*DTSTEP/CSOIL(N))

        TG1SN=TG1(N)+DTG1SN
        IF((TG1SN-TPSNB(N))*(TG1(N)-TPSNB(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(1)*                                 &
                 (TG1SN-TPSNB(N))*CSOIL(N)/DTSTEP
          TG1SN=TPSNB(N)
          ENDIF

        TG2SN=TG2(N)+DTG2SN
        IF((TG2SN-TPSNB(N))*(TG2(N)-TPSNB(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(2)*                                 &
                 (TG2SN-TPSNB(N))*CSOIL(N)/DTSTEP
          TG2SN=TPSNB(N)
          ENDIF

        TG4SN=TG4(N)+DTG4SN
        IF((TG4SN-TPSNB(N))*(TG4(N)-TPSNB(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(3)*                                 &
                 (TG4SN-TPSNB(N))*CSOIL(N)/DTSTEP
          TG4SN=TPSNB(N)
          ENDIF

        TG1(N)=TG1SF(N)*(1.-AREASC)+TG1SN*AREASC
        TG2(N)=TG2SF(N)*(1.-AREASC)+TG2SN*AREASC
        TG4(N)=TG4SF(N)*(1.-AREASC)+TG4SN*AREASC


        
!****   SET CANOPY TEMPERATURE ABOVE SNOWPACK TO SURFACE SNOW
!****   TEMPERATURE (TPSN1)
! Problem: as it's now coded, TC1, TC2, and TC4 cannot change if areasc=1
        
        DTC1SN=TPSN1(N)-TC1(N)
        DTC2SN=TPSN1(N)-TC2(N)
        DTC4SN=TPSN1(N)-TC4(N)

        TC1(N)=TC1SF(N)*(1.-AREASC)+TPSN1(N)*AREASC
        TC2(N)=TC2SF(N)*(1.-AREASC)+TPSN1(N)*AREASC
        TC4(N)=TC4SF(N)*(1.-AREASC)+TPSN1(N)*AREASC

!        TC1(N)=TC1SF(N)*(1-AREASC)+TC1_ORIG(N)*AREASC
!        TC2(N)=TC2SF(N)*(1-AREASC)+TC2_ORIG(N)*AREASC
!        TC4(N)=TC4SF(N)*(1-AREASC)+TC4_ORIG(N)*AREASC

        HSNACC(N)=HSNACC(N)-AREASC*CCANOP(N)*                                  &
              (DTC1SN*AR1(N)+DTC2SN*AR2(N)+DTC4SN*AR4(N))/DTSTEP



        EVSNOW(N)=EVSN
        esno(n)=evsnow(n)*asnow(n) ! to have esno in mm/20min (03-17-99)
        SHFLUXS(N)=SHFLS 
        HLWUPS(N) =HUPS 
        GHFLUXS(N)=AREA(1)*(HFTDS1(N)+DHFT1(N)*DTG1SN) +                       &
                   AREA(2)*(HFTDS2(N)+DHFT2(N)*DTG2SN) +                       &
                   AREA(3)*(HFTDS4(N)+DHFT4(N)*DTG4SN)
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
        GHTSKIN(N)=(1.-ASNOW(N))*                                               &
              (GHFLUX1(N)*AR1(N)+GHFLUX2(N)*AR2(N)+GHFLUX4(N)*AR4(N))          &
              -ASNOW(N)*ghfluxsno(N)
        ESOI(N)=(1.-ASNOW(N))*                                                 &
             (EVSURF1(N)*AR1(N)+EVSURF2(N)*AR2(N)+EVSURF4(N)*AR4(N)) 
        EVEG(N)=(1.-ASNOW(N))*                                                 &
             (EVROOT1(N)*AR1(N)+EVROOT2(N)*AR2(N)+EVROOT4(N)*AR4(N)) 
        EINT(N)=(1.-ASNOW(N))*                                                 &
             (EVINT1(N)*AR1(N)+EVINT2(N)*AR2(N)+EVINT4(N)*AR4(N)) 
        ESATFR(N)=(EVROOT1(N)+EVSURF1(N))*AR1(N)
        ESATFR(N)=ESATFR(N)/(ESATFR(N) +                                       &
                  AR2(N)*(EVROOT2(N) +  EVSURF2(N)) +                          &
                  AR4(N)*(EVROOT4(N) + EVSURF4(N)) + 1.E-20)
        ESATFR(N)=AMIN1(1., AMAX1(0., ESATFR(N)))
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
        ! zbar function - reichle, 29 Jan 2022 (minus sign applied in call to GNDTMP)
        ZBAR = catch_calc_zbar( bf1(n), bf2(n), catdef(n) )  
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

      CALL WUPDAT (                                                            &
                     NCH, DTSTEP,                                              &
                     EVAPFR, SATCAP, TG1, RA1, RC,                             &
                     AR1,AR2,AR4,CDCR1, ESATFR,                                &
                     RZEQOL,SRFMN,WPWET,VGWMAX, POROS,                         &
                     BF1, BF2, ARS1, ARS2, ARS3,                               &
                     CAPAC, RZEXC, CATDEF, SRFEXC,                             &
                     ESOI, EVEG, EINT,                                         &
                     ECORR                                                     &
                    )

!**** UPDATE SENSIBLE HEAT IF WATER LIMITATIONS WERE IMPOSED:
      DO N=1,NCH
        IF(ECORR(N) .GT. 0.) THEN
          SHFLUX(N)=SHFLUX(N)+ECORR(N)*ALHE
          EVAP(N)=EVAP(N)-ECORR(N)
          HLATN(N)=ESOI(N)*ALHE + EVEG(N)*ALHE + EINT(N)*ALHE + ESNO(N)*ALHS
          ENDIF
        ENDDO


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

! Add differences due to adjustments to land mosture prognostics
      do n=1,nch
         if(werror(n) .le. 0.) runsrf(n)=runsrf(n)-werror(n)/dtstep
         if(werror(n) .gt. 0.) then
           edif=werror(n)/dtstep
           EVAP(N)=EVAP(N)-EDIF
           HLATN(N)=HLATN(N)-EDIF*ALHE
           EVEGFRC=EVEG(N)/(EVEG(N)+ESOI(N)+1.E-20)
           EVEG(N)=EVEG(N)-EDIF*EVEGFRC
           ESOI(N)=ESOI(N)-EDIF*(1.-EVEGFRC)
           SHFLUX(N)=SHFLUX(N)+EDIF*ALHE
           endif
         enddo


    ! after revisions of calc_soil_moist() the call to partition is now obsolete 
    ! - reichle, 3 Apr 2012     
    !
    !  CALL PARTITION (                                                         &
    !                  NCH,DTSTEP,DZSF,RZEXC,  RZEQOL,VGWMAX,CDCR1,CDCR2,       &
    !                  PSIS,BEE,poros,WPWET,                                    &
    !                  ars1,ars2,ars3,ara1,ara2,ara3,ara4,                      &
    !                  arw1,arw2,arw3,arw4,BUG,                                 &
    !                  SRFEXC,CATDEF,RUNSRF,                                    &
    !                  AR1, AR2, AR4,srfmx,srfmn,  SWSRF1,SWSRF2,SWSRF4,RZI     &
    !                 )

      do n=1,nch
        hold=csoil(n)*(ar1old(n)*tg1(n)+ar2old(n)*tg2(n)+ar4old(n)*tg4(n)) +   &
           ccanop(n)*(ar1old(n)*tc1(n)+ar2old(n)*tc2(n)+ar4old(n)*tc4(n))
        hnew=csoil(n)*(ar1(n)*tg1(n)+ar2(n)*tg2(n)+ar4(n)*tg4(n)) +            &
           ccanop(n)*(ar1(n)*tc1(n)+ar2(n)*tc2(n)+ar4(n)*tc4(n))
        shflux(n)=shflux(n)-(hnew-hold)/dtstep
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

        EINT(N) = EINT(N) * ALHE
        ESOI(N) = ESOI(N) * ALHE
        EVEG(N) = EVEG(N) * ALHE
        ESNO(N) = ESNO(N) * ALHS
         
        TCANOP=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)
        TSOIL=AR1(N)*TG1(N)+AR2(N)*TG2(N)+AR4(N)*TG4(N)
        TSURF(N)=(1.-ASNOW0(N))*TCANOP+ASNOW0(N)*TPSN1(N) ! gkw: assumes fveg1+fveg2=1

        if(asnow(n) .eq. 0) then
          tpsn1(n)=amin1( TSURF(N),tf )
          tpsnb(n)=amin1( TSURF(N),tf )
        endif

        ! reichle, 1 May 2013, fix TWLAND<0 bug, use correct calculation in 

        CALL CATCH_CALC_WTOTL( 1, CDCR2(N:N), WPWET(N:N),                          &
             SRFEXC(N:N), RZEXC(N:N), CATDEF(N:N), CAPAC(N:N), WESNN(:,N:N),             &
             WTOT(N:N) )
        
        ENTOT(N) =                                                             &
            sum(HTSNNN(1:N_snow,N)) + TSOIL*CSOIL(N) + TCANOP*CCANOP(N) + sum(GHTCNT(1:N_gt,N))

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

!!!    IF (ityp1(N) .ge. 1) THEN ! gkw: would do TC correction on all types in unified 
!!!       IF (ityp1(N) .eq.-1) THEN ! gkw: using Helfand with viscous sublayer; don't need TC correction; using -1 skips this block
     IF (LAND_FIX .OR. (ityp1(N) .ne. 1)) THEN   ! jkolassa Oct 2020: changed to be equivalent to catchment.F90; was previously set to always be false
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
          ! (because CSOIL for is much larger) gkw: does not apply for unified code; would apply corrction for type 1

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
        TCSX=TGS_ORIG(N)

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
        ENDIF ! IF ( (.NOT. LAND_FIX) .OR. (ASNOW0(N) .EQ. 0. ) ) 

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
        EVAPXS=ETURBS(N)+DEDQAS(N)*DQSS(N)*(TPSN1(N)-TGS_ORIG(N))
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
        SHFLUXXS=HSTURBS(N)+DHSDTCS(N)*(TPSN1(N)-TGS_ORIG(N))
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
         write (*,*) ITYP1(n_out)  
         write (*,*) ITYP2(n_out)  
         write (*,*) FVEG1(n_out)  
         write (*,*) FVEG2(n_out)  
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
         write (*,*)  SWNETF(n_out)    
         write (*,*)  SWNETS(n_out)    
         write (*,*)   HLWDWN(n_out)    
         write (*,*)  PSUR(n_out)    
         write (*,*)   ZLAI(n_out)    
         write (*,*)    GREEN(n_out)    
         write (*,*)  SQSCAT(n_out)    
         write (*,*)  RSOIL1(n_out)    
         write (*,*)  RSOIL2(n_out)    
         write (*,*)    RDC(n_out)    
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
    END SUBROUTINE CATCHCN

!**** ===================================================
!**** ///////////////////////////////////////////////////
!**** ===================================================
!****
!**** [ BEGIN FLUXES ]
!****
      SUBROUTINE FLUXES (                                                      &
                            NCH,   FVEG, DTSTEP, QSATTC, DQSDTC,               &
                          ETURB,  DEDQA,  DEDTC, HSTURB, DHSDQA, DHSDTC,       &
                             RC, DRCDQA, DRCDTC,                               &
                          SWNET, HLWDWN, ALWRAD, BLWRAD,                       &
                             QM,  CSOIL, CCANOP,   PSUR,                       &
                          HFTDS, DHFTDS, RD, RSURF, POTFRC,                    &
                             TG,  TC,     QA,                                  &
                           EVAP, EVROOT, EVSURF, EVINT,                        &
                           SHFLUX,  HLWUP, GHFLUX                              &
                        )



!****
!**** THIS SUBROUTINE COMPUTES THE FLUXES OF LATENT AND SENSIBLE HEAT
!**** FROM THE SURFACE THROUGH AN ENERGY BALANCE CALCULATION.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: FVEG, QSATTC, DQSDTC, ETURB,      &
            DEDQA, DEDTC, HSTURB, DHSDQA, DHSDTC, DRCDQA, DRCDTC,       &
            SWNET, HLWDWN, ALWRAD, BLWRAD, QM, CSOIL, CCANOP, PSUR,         &
            HFTDS, DHFTDS, RD, RSURF, POTFRC, RC

      REAL, INTENT(INOUT), DIMENSION(NCH) :: TG, TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, EVROOT, EVSURF, EVINT,     &
            SHFLUX, HLWUP, GHFLUX

      INTEGER  N
      REAL, DIMENSION(NCH) :: CDEEPS, CONST, EPLANT,EANEW, ESATNW, EHARMN,  &
            DETERM, DENOM, EDIF, CCAPAC, EM, EA, ESATTC, DESDTC, DEDEA,     &
            DHSDEA, DRCDEA
      REAL EVAPW, EVAPD, SHFLW, SHFLD, GHFLW, GHFLD, HLWUPW, HLWUPD, ADUM1, &
           ADUM2, ADUM3, ADUM4, ADUM5, ADUM6, ADUM7, ADUM8, ADUM9, ADUM10,  &
           ADUM11,ADUM12, DEADRY, DTCDRY, DTGDRY, DEAWET, DTCWET, DTGWET,   &
           DEA, DTC, DTG,  A11, A12, A13, A21, A22, A23, A31, A32, A33,     &
           B11, B12, B21, B22, Q1, Q2, Q3, R1, R2, DHLWTC, HLWTG, DEDTG,    &
           DHDTG, RHOTERM, ESATTG, RHOAIR, HLWTC, ENBAL, ENBAL1, ENBAL2,    &
           WETFRC, EVSUM, TCOLD, TGOLD, ESATTC_NEW, ESATTG_NEW, EAOLD,      &
           TERM1, TERM2, DRCDT, DRCDE, RNEGTEST
      LOGICAL DEBUG, CHOKED, CHOKEW

      REAL, PARAMETER :: SMALL = 1.E-20   ! (ARBITRARY VERY SMALL NUMBER)

      DATA DEBUG /.FALSE./
!****
!**** -------------------------------------------------------------------

      DO 200 N = 1, NCH

      EM(N)     = QM(N) * PSUR(N) / EPSILON
      EA(N)     = QA(N) * PSUR(N) / EPSILON
      ESATTC(N) = QSATTC(N) * PSUR(N) / EPSILON
      DESDTC(N) = DQSDTC(N) * PSUR(N) / EPSILON
      DEDEA(N)  = DEDQA(N) * EPSILON / PSUR(N)
      DHSDEA(N) = DHSDQA(N) * EPSILON / PSUR(N)
      DRCDEA(N) = DRCDQA(N) * EPSILON / PSUR(N)


!**** REMOVE THIS LINE:
!      RC(N)=AMIN1(AMAX1(RC(N),1.), 1.E10)
!****
      HLWTC = ALWRAD(N) + BLWRAD(N) * TC(N)
      HLWTG = ALWRAD(N) + BLWRAD(N) * TG(N)
      RHOAIR = PSUR(N) * 100. / (RGAS * TC(N))
      CONST = RHOAIR * EPSILON / PSUR(N)
      DHLWTC = BLWRAD(N)
      DHDTG = RHOAIR*CPAIR/(RD(N)+small)
      DEDTG = (RHOAIR*EPSILON/PSUR(N))*DESDTC(N)/(RD(N)+RSURF(N))
      RHOTERM = RHOAIR*EPSILON/PSUR(N)
      ESATTG = ESATTC(N)+DESDTC(N)*(TG(N)-TC(N))
      DRCDT = DRCDTC(N)
      DRCDE = DRCDEA(N)


      CALL MATRIX_CALC (DTSTEP, CSOIL(N), DHLWTC, DHDTG, DEDTG,             &
                        HSTURB(N), DHFTDS(N), FVEG(N), BLWRAD(N),           &
                        RHOTERM, RD(N), RSURF(N),                           &
                        CCANOP(N), DHSDTC(N), DESDTC(N), RC(N), ESATTC(N),  &
                        ESATTG, DRCDT, EA(N), DHSDEA(N), DRCDE,             &
                        SWNET(N), HLWDWN(N), TG(N), TC(N), HLWTG,           &
                        HFTDS(N), DEDTC(N), HLWTC, ETURB(N), DEDEA(N),      &
                        DTCDRY, DTGDRY, DEADRY )

! Ensure that resistances haven't gone negative; if they have, set
! derivatives of resistances to zero and recompute

      RNEGTEST = RC(N) + DRCDTC(N)*DTCDRY + DRCDEA(N)*DEADRY
      IF(RNEGTEST .LE. 0.) THEN
         DRCDT=0.
         DRCDE=0.
         CALL MATRIX_CALC (DTSTEP, CSOIL(N), DHLWTC, DHDTG, DEDTG,          &
                        HSTURB(N), DHFTDS(N), FVEG(N), BLWRAD(N),           &
                        RHOTERM, RD(N), RSURF(N),                           &
                        CCANOP(N), DHSDTC(N), DESDTC(N), RC(N), ESATTC(N),  &
                        ESATTG, DRCDT, EA(N), DHSDEA(N), DRCDE,             &
                        SWNET(N), HLWDWN(N), TG(N), TC(N), HLWTG,           &
                        HFTDS(N), DEDTC(N), HLWTC, ETURB(N), DEDEA(N),      &
                        DTCDRY, DTGDRY, DEADRY )
         ENDIF



      CHOKED=.FALSE.
      IF(DTCDRY > DTCCRIT) THEN 
        DTCDRY=DTCCRIT
        CHOKED=.TRUE.
        ENDIF
      IF(DTCDRY .LT. -DTCCRIT) THEN
        DTCDRY=-DTCCRIT
        CHOKED=.TRUE.
        ENDIF
      IF(DTGDRY .GT. DTGCRIT) THEN
        DTGDRY=DTGCRIT
        CHOKED=.TRUE.
        ENDIF
      IF(DTGDRY .LT. -DTGCRIT) THEN
        DTGDRY=-DTGCRIT
        CHOKED=.TRUE.
        ENDIF
      IF(DEADRY .GT. EA(N)/2.) THEN
        DEADRY=EA(N)/2.
        CHOKED=.TRUE.
        ENDIF
      IF(DEADRY .LT. -EA(N)/2.) THEN
        DEADRY=-EA(N)/2.
        CHOKED=.TRUE.
        ENDIF


!**** NOW SOLVE MATRIX GIVEN PRESENCE OF INTERCEPTION WATER.  ASSUME
!**** ALL EVAPORATION IS FROM (OR TO) INTERCEPTION RESERVOIR.
!****
      A11 = CSOIL(N)/DTSTEP +                                               &
              DHLWTC +                                                      &
              DHDTG +                                                       &
              DHFTDS(N)
      A12 = -FVEG(N)*BLWRAD(N)-DHDTG

      A21 = -FVEG(N)*BLWRAD(N)-DHDTG
      A22 = CCANOP(N)/DTSTEP +                                              &
             2.*FVEG(N)*BLWRAD(N)+DHDTG+DHSDTC(N)+                          &
             DHSDEA(N)*DESDTC(N) +                                          &
             ALHE*DEDTC(N) + ALHE*DEDEA(N)*DESDTC(N)


      Q1 =  SWNET(N)*(1.-FVEG(N)) +                                         &
              FVEG(N)*HLWTC +                                               &
              (1.-FVEG(N))*HLWDWN(N) -                                      &
              DHDTG*(TG(N)-TC(N)) -                                         &
              HLWTG -                                                       &
              HFTDS(N)
      Q2 =  SWNET(N)*FVEG(N) +                                              &
              FVEG(N)*HLWDWN(N) -                                           &
              2.*FVEG(N)*HLWTC +                                            &
              FVEG(N)*HLWTG +                                               &
              DHDTG*(TG(N)-TC(N)) -                                         &
              HSTURB(N) -                                                   &
              ALHE*ETURB(N) -                                               &
              DHSDEA(N)*(ESATTC(N)-EA(N)) -                                 &
              ALHE*DEDEA(N)*(ESATTC(N)-EA(N))

        ADUM1=SWNET(N)*FVEG(N)
        ADUM2=FVEG(N)*HLWDWN(N)
        ADUM3=2.*FVEG(N)*HLWTC
        ADUM4=FVEG(N)*HLWTG
        ADUM5=DHDTG*(TG(N)-TC(N))
        ADUM6=HSTURB(N)
        ADUM7=ALHE*ETURB(N)
        ADUM8=DHSDEA(N)*(ESATTC(N)-EA(N))
        ADUM9=ALHE*DEDEA(N)*(ESATTC(N)-EA(N))
        ADUM10=ADUM1+ADUM2-ADUM3+ADUM4+ADUM5-ADUM6-ADUM7-ADUM8-ADUM9



      DTCWET = ( Q1*A21 - A11*Q2 ) / ( A12*A21 - A11*A22 )
      DTGWET = ( Q1 - A12*DTCWET ) / A11
      DEAWET = ESATTC(N)+DESDTC(N)*DTCWET-EA(N)



!**** CHECK:
      ADUM11=A11*DTGWET+A12*DTCWET
      ADUM12=A21*DTGWET+A22*DTCWET


      CHOKEW=.FALSE.
      IF(DTCWET .GT. DTCCRIT) THEN 
        DTCWET=3.
        DEAWET = ESATTC(N)+DESDTC(N)*DTCWET-EA(N)
        CHOKEW=.TRUE.
        ENDIF
      IF(DTCWET .LT. -DTCCRIT) THEN
        DTCWET=-3.
        DEAWET = ESATTC(N)+DESDTC(N)*DTCWET-EA(N)
        CHOKEW=.TRUE.
        ENDIF
      IF(DTGWET .GT. DTGCRIT) THEN
        DTGWET=3.
        CHOKEW=.TRUE.
        ENDIF
      IF(DTGWET .LT. -DTGCRIT) THEN
        DTGWET=-3.
        CHOKEW=.TRUE.
        ENDIF

!    PUT IN LINES TO CONSERVE ENERGY, IF CHOKE = .TRUE.

      TCOLD=TC(N)
      TGOLD=TG(N)
      EAOLD=EA(N)

      WETFRC=POTFRC(N)
      IF(ETURB(N)+DEDTC(N)*DTCDRY+DEDEA(N)*DEADRY .LT. 0. .AND.             &
         ETURB(N)+DEDTC(N)*DTCWET+DEDEA(N)*DEAWET .LT. 0.) WETFRC=1.

      EVAPW=ETURB(N)+DEDTC(N)*DTCWET+DEDEA(N)*DEAWET
      EVAPD=ETURB(N)+DEDTC(N)*DTCDRY+DEDEA(N)*DEADRY
      EVAP(N)=EVAPW*WETFRC+EVAPD*(1.-WETFRC)

      HLWUPW=FVEG(N)*(HLWTC + DHLWTC*DTCWET) +                              &
           (1.-FVEG(N))*(HLWTG + DHLWTC*DTGWET)
      HLWUPD=FVEG(N)*(HLWTC + DHLWTC*DTCDRY) +                              &
           (1.-FVEG(N))*(HLWTG + DHLWTC*DTGDRY)
      HLWUP(N) = HLWUPW*WETFRC + HLWUPD*(1.-WETFRC)

      GHFLW=HFTDS(N)+DHFTDS(N)*DTGWET
      GHFLD=HFTDS(N)+DHFTDS(N)*DTGDRY
      GHFLUX(N) = GHFLW*WETFRC + GHFLD*(1.-WETFRC)

      SHFLW=HSTURB(N)+DHSDTC(N)*DTCWET+DHSDEA(N)*DEAWET
      IF(CHOKEW) THEN
        SHFLW = SWNET(N) + HLWDWN(N) - HLWUPW - ALHE*EVAPW -                &
              GHFLW - CCANOP(N)*DTCWET/DTSTEP - CSOIL(N)*DTGWET/DTSTEP
        ENDIF
      SHFLD=HSTURB(N)+DHSDTC(N)*DTCDRY+DHSDEA(N)*DEADRY
      IF(CHOKED) THEN
        SHFLD = SWNET(N) + HLWDWN(N) - HLWUPD - ALHE*EVAPD -                &
              GHFLD - CCANOP(N)*DTCDRY/DTSTEP - CSOIL(N)*DTGDRY/DTSTEP
        ENDIF
      SHFLUX(N)=SHFLW*WETFRC + SHFLD*(1.-WETFRC)



      TG(N) = TG(N) + DTGDRY*(1.-WETFRC) + DTGWET*WETFRC
      TC(N) = TC(N) + DTCDRY*(1.-WETFRC) + DTCWET*WETFRC
      EA(N) = EA(N) + DEADRY*(1.-WETFRC) + DEAWET*WETFRC


      ENBAL=SWNET(N) + HLWDWN(N) - HLWUP(N) - ALHE*EVAP(N) - SHFLUX(N) -    &
              GHFLUX(N) - CCANOP(N)*(TC(N)-TCOLD)/DTSTEP -                  &
              CSOIL(N)*(TG(N)-TGOLD)/DTSTEP


      DTC=DTCDRY+WETFRC*(DTCWET-DTCDRY)
      DTG=DTGDRY+WETFRC*(DTGWET-DTGDRY)
      DEA=DEADRY+WETFRC*(DEAWET-DEADRY)
      ENBAL1=(1.-FVEG(N))*SWNET(N) +  FVEG(N)*(HLWTC + DHLWTC*DTC) +        &
          (1.-FVEG(N))*HLWDWN(N) - (HLWTG + DHLWTC*DTG) -                   &
          DHDTG*(TG(N)-TC(N)) -                                             &
          RHOTERM*ALHE*(ESATTC(N)+DESDTC(N)*(TG(N)-TCOLD)                   &
          -EA(N))/(RD(N)+RSURF(N)) - GHFLUX(N) - CSOIL(N)*DTG/DTSTEP
      ENBAL2=FVEG(N)*SWNET(N) + FVEG(N)*HLWDWN(N) -                         &
         2.*FVEG(N)*(HLWTC + DHLWTC*DTC) + FVEG(N)*(HLWTG + DHLWTC*DTG) +   &
         DHDTG*(TG(N)-TC(N)) - SHFLUX(N) -                                  &
         RHOTERM*ALHE*(ESATTC(N)+DESDTC(N)*(TC(N)-TCOLD)                    &
         - EA(N))/RC(N) - CCANOP(N)*DTC/DTSTEP


! WARNING: THIS GHFLUX IS THE REAL GROUND HEAT FLUX, AND DOES NOT INCLUDE
! THE TEMPERATURE VARIATION

!**** MAKE SURE EA REMAINS POSITIVE

      EA(N) = AMAX1(EA(N), 0.0)

!**** - - - - - - - - - - - - - - - - - - - - - - -
!**** ADJUSTMENTS


!**** ADJUST DEA AND DTC IF SOLUTIONS WERE PATHOLOGICAL:
!****
!      ESATNW = ESATTC(N)+DESDTC(N)*DTC
!      EANEW = EA(N) + DEA



!**** PATHOLOGICAL CASES. 

!**** CASE 1:EVAP AND ETURB HAVE OPPOSITE SIGN, IMPLYING THAT 
!**** CALCULATED RC VALUE IS INAPPROPRIATE.
!**** ASSUME ZERO EVAPORATION FOR THE TIME STEP; INCLUDE EXCESS
!**** OR DEFICIT OF HEAT WITH SENSIBLE HEAT FLUX..

!      IF( EVAP(N)*ETURB(N) .LT. 0. )  THEN
!        SHFLUX(N)=SHFLUX(N)+ALHE*EVAP(N)
!        EVAP(N) = 0.
!        EA(N)=EM(N)
!        ENDIF

!      IF( EVAP(N) .LT. 0. .AND. EM(N).LT.ESATNW )  THEN 
!        CHOKE = .FALSE.
!        DEA = EM(N) - EA(N)
!        DTC = ( Q1 + ALHE*ETURB(N) - DHSDEA(N)*DEA ) / 
!     &            ( A11 - ALHE*DEDTC(N) )
!        EVAP(N) = 0.
!        SHFLUX(N) = HSTURB(N) + DHSDEA(N)*DEA +
!     &                                            DHSDTC(N)*DTC
!        ENDIF



!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ESATTC_NEW=ESATTC(N)+DESDTC(N)*DTCDRY
      ESATTG_NEW=ESATTG+DESDTC(N)*DTGDRY
      EVSURF(N)=(1.-WETFRC)*RHOTERM*(ESATTG_NEW-(EAOLD+DEADRY))/(RD(N)+RSURF(N))

      TERM1=RC(N) + DRCDTC(N)*DTCDRY + DRCDEA(N)*DEADRY
      IF(TERM1 .LE. .01) TERM1=.01
      EVROOT(N)=(1.-WETFRC)*RHOTERM*(ESATTC_NEW-(EAOLD+DEADRY))/TERM1

      EVINT(N)=WETFRC*(ETURB(N)+DEDTC(N)*DTCWET+DEDEA(N)*DEAWET)

!**** CORRECT INDIVIDUAL EVAP FLUXES IN CASES OF "CHOKE" OR
!**** OTHER CORRECTION:

      IF(EVAP(N) .LE. 0.) THEN
        EVINT(N)=EVAP(N)
        EVSURF(N)=0.
        EVROOT(N)=0.
        ENDIF

      IF(EVAP(N) .GT. 0.) THEN

        IF(EVINT(N).LT.SMALL) EVINT(N)=SMALL
        IF(EVSURF(N).LT.SMALL) EVSURF(N)=SMALL
        IF(EVROOT(N).LT.SMALL) EVROOT(N)=SMALL

        EVSUM=EVSURF(N)+EVROOT(N)+EVINT(N)

        EVSURF(N)=EVSURF(N)*EVAP(N)/EVSUM
        EVROOT(N)=EVROOT(N)*EVAP(N)/EVSUM
        EVINT(N)=EVINT(N)*EVAP(N)/EVSUM

        ENDIF




      QA(N) = EA(N) * EPSILON / PSUR(N)

 200  CONTINUE

      RETURN
      END SUBROUTINE FLUXES
!****
!**** [ END FLUXES ]


!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN MATRIX_CALC ]
!****
      SUBROUTINE MATRIX_CALC (DTSTEP, CSOIL, DHLWTC, DHDTG, DEDTG, HSTURB,  &
                              DHFTDS, FVEG, BLWRAD, RHOTERM, RD, RSURF,     &
                              CCANOP, DHSDTC, DESDTC, RC, ESATTC,           &
                              ESATTG, DRCDTC, EA, DHSDEA, DRCDEA, SWNET,    &
                              HLWDWN, TG, TC, HLWTG, HFTDS, DEDTC, HLWTC,   &
                              ETURB, DEDEA,                                 &
                              DTC, DTG, DEA )

      REAL, INTENT(IN) ::     DTSTEP, CSOIL, DHLWTC, DHDTG, DEDTG,          &
                              DHFTDS, FVEG, BLWRAD, RHOTERM, RD, RSURF,     &
                              CCANOP, DHSDTC, DESDTC, RC, ESATTC,           &
                              DRCDTC, EA, DHSDEA, DRCDEA, SWNET, HLWDWN,    &
                              TG, TC, HLWTG, HFTDS, ESATTG, HSTURB,         &
                              DEDTC, HLWTC, ETURB, DEDEA

      REAL, INTENT(OUT) ::    DTC, DTG, DEA

      REAL*8 :: A11, A12, A13, A21, A22, A23, A31, A32, A33, Q1, Q2, Q3,      &
              B11, B12, R1, B21, B22, R2, DETERM


! ---------------------------------------------------
!
!     Solve 3x3 matrix equation


      A11 = CSOIL/DTSTEP +                                                  &
              DHLWTC +                                                      &
              DHDTG +                                                       &
              ALHE*DEDTG +                                                  &
              DHFTDS
      A12 = -FVEG*BLWRAD-DHDTG
      A13 = -RHOTERM*ALHE/(RD+RSURF)

      A21 = -FVEG*BLWRAD-DHDTG
      A22 = CCANOP/DTSTEP +                                                 &
             2.*FVEG*BLWRAD+DHDTG+DHSDTC+                                   &
             RHOTERM*ALHE*DESDTC/RC-                                        &
             RHOTERM*ALHE*ESATTC*DRCDTC/                                    &
                            (RC*RC) +                                       &
             RHOTERM*ALHE*EA*DRCDTC/(RC*RC)
      A23 =  DHSDEA-RHOTERM*ALHE/RC -                                       &
             RHOTERM*ALHE*ESATTC*DRCDEA/                                    &
                            (RC*RC) +                                       &
             RHOTERM*ALHE*EA*DRCDEA/(RC*RC)

      A31 = RHOTERM*DESDTC/(RD+RSURF)
      A32 = -DEDTC+RHOTERM*DESDTC/RC -                                      &
             RHOTERM*ESATTC*DRCDTC/(RC*RC) +                                &
             RHOTERM*EA*DRCDTC/(RC*RC)
      A33 = -DEDEA - RHOTERM/(RD+RSURF) -                                   &
             RHOTERM/RC -                                                   &
             RHOTERM*ESATTC*DRCDEA/(RC*RC) +                                &
             RHOTERM*EA*DRCDEA/(RC*RC)

      DETERM = A11*A22*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31 -      &
               A12*A21*A33 - A11*A32*A23

      IF(DETERM .GT. -1.E-8) THEN
         IF(ABS(A22*A33-A32*A23) .GT. 1.E-8) THEN
                      A11=(-1.E-8 - A12*A23*A31 - A13*A21*A32               &
                                  + A13*A22*A31 + A12*A21*A33) /            &
                                    (A22*A33 - A32*A23)
         ELSEIF(ABS(A11*A33-A13*A31) .GT. 1.E-8) THEN
                      A22=(-1.E-8 - A12*A23*A31 - A13*A21*A32               &
                                  + A12*A21*A33 + A11*A32*A23) /            &
                                    (A11*A33 - A13*A31)
         ELSE
                      A33=(-1.E-8 - A12*A23*A31 - A13*A21*A32               &
                                  + A13*A22*A31 + A11*A32*A23) /            &
                                    (A11*A22 - A12*A21)
         ENDIF
         DETERM = A11*A22*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31 -   &
                  A12*A21*A33 - A11*A32*A23
      ENDIF

      Q1 =  SWNET*(1.-FVEG) +                                               &
              FVEG*HLWTC +                                                  &
              (1.-FVEG)*HLWDWN -                                            &
              DHDTG*(TG-TC) -                                               &
              RHOTERM*ALHE*(ESATTG-EA)/(RD+RSURF) -                         &
              HLWTG -                                                       &
              HFTDS
      Q2 =  SWNET*FVEG +                                                    &
              FVEG*HLWDWN -                                                 &
              2.*FVEG*HLWTC +                                               &
              FVEG*HLWTG +                                                  &
              DHDTG*(TG-TC) -                                               &
              HSTURB -                                                      &
              RHOTERM*ALHE*(ESATTC-EA)/RC
      Q3 = ETURB -                                                          &
              RHOTERM*(ESATTG-EA)/(RD+RSURF) -                              &
              RHOTERM*(ESATTC-EA)/RC

      B11=A12*A21-A11*A22
      B12=A13*A21-A11*A23
      R1=Q1*A21-Q2*A11

      B21=A12*A31-A11*A32
      B22=A13*A31-A11*A33
      R2=Q1*A31-Q3*A11

      DETERM=B11*B22-B12*B21
      IF(ABS(DETERM) .LT. 1.E-8) B11=(SIGN(1.d-8,DETERM)+B12*B21)/B22

      if(DETERM .ne. 0.) then
        DEA = ( R1*B21 - B11*R2 ) / ( B12*B21 - B11*B22 )
       else
        DEA = 0.
      endif
      if(B11 .ne. 0.) then
        DTC = ( R1 - B12*DEA ) / B11
       else
	DTC = 0.
      endif
      if(A11 .ne. 0.) then
        DTG = (Q1 - A12*DTC - A13*DEA)/A11
       else
        DTG = 0.
      endif

      RETURN
      END SUBROUTINE MATRIX_CALC


!**** [ END MATRIX_CALC ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN WUPDAT ]
!****
      SUBROUTINE WUPDAT (                                                   &
                           NCH, DTSTEP,                                     &
                           EVAP, SATCAP, TC, RA, RC,                        &
                           AR1,AR2,AR4,CDCR1, ESATFR,                       &
                           RZEQ,SRFMN,WPWET,VGWMAX, POROS,                  &
                           BF1, BF2, ARS1, ARS2, ARS3,                      &
                           CAPAC, RZEXC, CATDEF, SRFEXC,                    &
                           EVROOT, EVSURF, EVINT,                           &
                           ECORR                                            &
                          )
!****
!**** THIS SUBROUTINE ALLOWS EVAPOTRANSPIRATION TO ADJUST THE WATER
!**** CONTENTS OF THE INTERCEPTION RESERVOIR AND THE SOIL LAYERS.
!****
      IMPLICIT NONE
!****
      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: EVAP,  SATCAP, TC, RA, RC, AR1,   &
            AR2, AR4, CDCR1, ESATFR, RZEQ, SRFMN, WPWET, VGWMAX,            &
            POROS, BF1, BF2, ARS1, ARS2, ARS3

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CAPAC, RZEXC, CATDEF,          &
            SRFEXC, EVROOT, EVSURF, EVINT

      REAL, INTENT(OUT), DIMENSION(NCH) :: ECORR

      INTEGER  N
      REAL  EGRO,CNDSAT,CNDUNS,CNDV,CNDS,WILT,EGROMX,RZEMAX
      REAL :: ZBAR1,SYSOIL,ET_CATDEF,AR1eq
!****
!**** -----------------------------------------------------------------

      DO 100 N = 1, NCH
!****
!**** PARTITION EVAP BETWEEN INTERCEPTION AND GROUND RESERVOIRS.
!****

!**** ENSURE THAT INDIVIDUAL CAPACITIES ARE NOT EXCEEDED.
      ECORR(N)=0.
      IF(EVINT(N) .GT. CAPAC(N)/DTSTEP) THEN
        ECORR(N)=EVINT(N)-CAPAC(N)/DTSTEP
        EVINT(N)=CAPAC(N)/DTSTEP
        ENDIF

      EVSURF(N)=EVSURF(N) + ECORR(N)
      ECORR(N)=0.
      IF(EVSURF(N) .GT. (SRFEXC(N)-SRFMN(N))/DTSTEP) THEN
        ECORR(N)=EVSURF(N)-(SRFEXC(N)-SRFMN(N))/DTSTEP
        EVSURF(N)=(SRFEXC(N)-SRFMN(N))/DTSTEP
        ENDIF

      EVROOT(N)=EVROOT(N) + ECORR(N)
      ECORR(N)=0.
      WILT=WPWET(N)*VGWMAX(N)
      RZEMAX=AMAX1(0.,RZEXC(N)+RZEQ(N)-WILT ) 
      IF(EVROOT(N) .GT. RZEMAX/DTSTEP) THEN
        ECORR(N)=EVROOT(N)-RZEMAX/DTSTEP
        EVROOT(N)=RZEMAX/DTSTEP
        ENDIF


!****
!**** SPECIAL CASE FOR CONDENSATION:
      IF(EVAP(N) .LT. 0.) THEN
        EVINT(N)=EVAP(N)
        EVSURF(N)=0.
        EVROOT(N)=0.
        ENDIF

!****
!**** REMOVE MOISTURE FROM RESERVOIRS:
!****

       IF (CATDEF(N) .LT. CDCR1(N)) THEN
          CAPAC(N) = AMAX1(0., CAPAC(N) - EVINT(N)*DTSTEP)
          RZEXC(N) = RZEXC(N) - EVROOT(N)*(1.-ESATFR(N))*DTSTEP
          SRFEXC(N) = SRFEXC(N) - EVSURF(N)*(1.-ESATFR(N))*DTSTEP
          IF (POROS(N) < PEATCLSM_POROS_THRESHOLD) THEN
             CATDEF(N) = CATDEF(N) + (EVSURF(N) + EVROOT(N))*ESATFR(N)*DTSTEP
             ! 05.12.98: FIRST ATTEMPT TO INCLUDE BEDROCK
          ELSE
             ! PEAT
             ! MB: accounting for water ponding on AR1
             ! same approach as for RZFLW (see subroutine RZDRAIN for
             ! comments)
             ZBAR1  = catch_calc_zbar( BF1(N), BF2(N), CATDEF(N) )  
             SYSOIL = (2.*bf1(N)*amin1(amax1(zbar1,0.),PEATCLSM_ZBARMAX_4_SYSOIL) + 2.*bf1(N)*bf2(N))/1000.
             SYSOIL = amin1(SYSOIL,poros(N))
             ET_CATDEF = SYSOIL*(EVSURF(N) + EVROOT(N))*ESATFR(N)/(1.*AR1(N)+SYSOIL*(1.-AR1(N)))
             AR1eq = (1.+ars1(N)*(catdef(N)))/(1.+ars2(N)*(catdef(N))+ars3(N)*(catdef(N))**2)
             CATDEF(N) = CATDEF(N) + (1.-AR1eq)*ET_CATDEF
          ENDIF
       ELSE
          CAPAC(N) = AMAX1(0., CAPAC(N) - EVINT(N)*DTSTEP)
          RZEXC(N) = RZEXC(N) -  EVROOT(N)*DTSTEP
          SRFEXC(N) = SRFEXC(N) - EVSURF(N)*DTSTEP
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
!**** [ BEGIN RSURFP2 ]
!****
      SUBROUTINE RSURFP2 (                                                     &
                         NCH, RSSAT, WET, QSATTC, QA, PSUR, WPWET,             &
                         RSURF                                                 &
                        )
!****
      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: RSSAT
      REAL, INTENT(IN), DIMENSION(NCH) :: WET, QSATTC, QA, WPWET,              &
                                            PSUR
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RSURF


      INTEGER N
      REAL, DIMENSION(NCH) :: EA, ESATTC
      REAL HESAT, ATERM, BTERM

! RDK 04/04/06
!  VALUES OF BARE SOIL SURFACE RESISTANCE AT WILTING POINT, SATURATION
!!!   PARAMETER (RSWILT=500.) ! gkw: 61u
!!!   PARAMETER (RSWILT=1.e3) ! gkw: 62u
!!!   PARAMETER (RSWILT=5.e3) ! gkw: 63u
!!!   PARAMETER (RSWILT=1.e9) ! gkw: 64u
!!!   PARAMETER (RSWILT=2.e3) ! gkw: 65u
!!!   PARAMETER (RSWILT=2000.)! gkw: 66u & 67u
!****
!**** -----------------------------------------------------------------

      DO 100 N = 1, NCH
!****
      EA(N)     = QA(N) * PSUR(N) / EPSILON
      ESATTC(N) = QSATTC(N) * PSUR(N) / EPSILON
      BTERM=(RSWILT-RSSAT) / (1./WPWET(N)**2 -1.)
      ATERM=RSSAT-BTERM
      RSURF(N) = ATERM + BTERM / (1.E-10 + WET(N))**2

!**** Account for subsaturated humidity at soil surface:
!****
      HESAT = ESATTC(N) * MIN( 1., WET(N)*2. )
      IF( EA(N) .LT. HESAT ) THEN
          RSURF(N)=RSURF(N)*( 1. + (ESATTC(N)-HESAT)/(HESAT-EA(N)) )
        ELSE
          RSURF(N)=1.E10
        ENDIF

 100  CONTINUE

      RETURN
      END SUBROUTINE RSURFP2

END MODULE CATCHMENT_CN_MODEL

! ================================== EOF ============================================
