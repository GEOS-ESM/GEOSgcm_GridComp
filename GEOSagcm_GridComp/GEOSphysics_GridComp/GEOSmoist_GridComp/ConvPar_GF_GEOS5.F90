MODULE ConvPar_GF_GEOS5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- This convective parameterization is build to attempt                 !
!  a smooth transition to cloud resolving scales as proposed            !
!  by Arakawa et al (2011, ACP). The scheme is  described               !
!  in the paper Grell and Freitas (ACP, 2014).                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- Implemented in GEOS5 AGCM by Saulo Freitas (July 2016)               !
!- Use the following reference for this implementation:                 !
!- Freitas et al (2018, JAMES/AGU, https://doi.org/10.1029/2017MS001251)!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
USE module_gate
USE MAPL
!
USE Henrys_law_ConstantsMod, ONLY: get_HenrysLawCts
!.. USE GTMP_2_GFCONVPAR, only : GTMP_2_GFCONVPAR_interface

 IMPLICIT NONE
 PRIVATE
 PUBLIC  GF_GEOS5_INTERFACE, maxiens, icumulus_gf, closure_choice, deep, shal, mid &
        ,DEBUG_GF,USE_SCALE_DEP,DICYCLE,TAU_DEEP,TAU_MID,Hcts &
        ,USE_TRACER_TRANSP, USE_TRACER_SCAVEN&
	,USE_FLUX_FORM, USE_FCT, USE_TRACER_EVAP,ALP1

 !- for internal debugging
 PUBLIC GF_GEOS5_DRV       !- for debugging purposes (set private for normal runs)
 LOGICAL :: wrtgrads = .false.
 INTEGER :: nrec = 0, ntimes = 0
 INTEGER, PARAMETER :: DEBUG_GF=1       != 0/1  !- debug flag
 !-
 !- plume spectral size
 INTEGER ,PARAMETER  :: maxiens = 3, deep=1 ,shal=2 , mid = 3
 CHARACTER(LEN=10),PARAMETER,DIMENSION(maxiens)  :: cumulus_type = (/ &
                                                   'deep      ' &
                                                  ,'shallow   ' &
                                                  ,'mid       ' &
                                                                    /)
 !------------------- namelist variables
 !-- plume to be activated (1 true, 0 false): deep, shallow, congestus
 INTEGER, DIMENSION(maxiens) :: icumulus_gf

 !-- choice for the closures:
 !--  deep   : 0 ensemble (all)          , 1 GR, 4 ll omega, 7 moist conv, 10 PB
 !--  shallow: 0 ensemble (Wstar/BLQE)   , 1 Wstar, 4 heat-engine or 7 BLQE
 !--  mid    : 0 ensemble (Wstar/BLQE/PB), 1 Wstar, 2 BLQE, 3 PB, 4 PB_BL
 INTEGER, DIMENSION(maxiens) :: closure_choice ! deep, shallow, congestus

 INTEGER :: USE_TRACER_TRANSP  != 0/1     - default 1

 INTEGER :: USE_TRACER_SCAVEN  != 0/1/2/3 - default 1

 INTEGER :: USE_FLUX_FORM      != 1/2/3   - default 1

 INTEGER :: USE_FCT            != 0/1     - default 1 (only for USE_FLUX_FORM     = 2)

 INTEGER :: USE_TRACER_EVAP    != 0/1     - default 1 (only for USE_TRACER_SCAVEN > 0)


 REAL    :: ALP1               !- 0/0.5/1: apply subsidence transport of LS/anvil cloud fraction using 
                               !-          time implicit discretization

 INTEGER :: USE_SCALE_DEP  != 0/1  !- scale dependence flag
 INTEGER :: DICYCLE        != 0/1  !- diurnal cycle flag

 REAL :: tau_deep = 5400. ! deep      convective timescale
 REAL :: tau_mid  = 3600. ! congestus convective timescale
 
 !------------------- internal variables
 !-- turn ON/OFF deep/shallow/mid plumes
 INTEGER, PARAMETER :: ON=1, OFF=0
 !-- gross entraiment rate: deep, shallow, congestus
 REAL   , PARAMETER,DIMENSION(maxiens) :: cum_entr_rate=(/&
                                          1.00e-4  & !deep
                                         ,1.40e-3  & !P9_P6 & !shallow
 !                                       ,2.00e-3  & !shallow
                                         ,9.00e-4  & !mid
                                                        /)
 !-- max cloud fraction allowed: deep, shallow, congestus
 REAL   , PARAMETER,DIMENSION(maxiens) :: mxclfr=(/&
                                          0.6  & !deep
                                         ,0.2  & !shallow
                                         ,0.4  & !mid
                                                         /)
!-- General control for the diverse options in GF
 LOGICAL, PARAMETER :: ENTRNEW        = .TRUE.  !- new entr formulation
 LOGICAL, PARAMETER :: COUPL_MPHYSICS = .TRUE.  !- coupling with cloud microphysics
                                                ! (do not change  to false)
 LOGICAL, PARAMETER :: MELT_GLAC      = .TRUE.  !- turn ON/OFF ice phase/melting
 LOGICAL, PARAMETER :: DOWNDRAFT      = .TRUE.  !- turn ON/OFF downdrafts
 LOGICAL, PARAMETER :: MOMENTUM       = .TRUE.  !- turn ON/OFF conv transp of momentum

 LOGICAL, PARAMETER :: FEED_3DMODEL   = .TRUE.  !- set "false" to not feedback the AGCM with the
                                                !- heating/drying/transport conv tendencies
 LOGICAL, PARAMETER :: USE_C1D        = .TRUE.  !- turn ON/OFF the 'c1d' detrainment approach

 INTEGER, PARAMETER :: VERT_DISCR = 1

 !-autonversion formulation: (1) original , (2) MP_GT
 INTEGER, PARAMETER :: CLOUDMP = 1

 !-autonversion formulation: (1) Kessler, (2) Berry, (3) NOAA, (4) Sundvisqt
 INTEGER, PARAMETER :: autoconv = 1 
 !-rainfall evaporation(1) orig (2) mix orig+new (3) new
 INTEGER, PARAMETER :: aeroevap = 1

 INTEGER, PARAMETER ::               &
                       maxens  = 1,  & ! 1  ensemble one on cap_max
                       maxens2 = 1,  & ! 1  ensemble two on precip efficiency
                       maxens3 = 16, & !16 ensemble three done in cup_forcing_ens16 for G3d
                       ensdim  = maxens*maxens2*maxens3,&
                       ens4    = 1
 
 !

 !- physical constants
 REAL, PARAMETER ::  &
  rgas    = 287.,    & ! J K-1 kg-1
  cp      = 1004.,   & ! J K-1 kg-1
  rv      = 461.,    & ! J K-1 kg-1
  p00     = 1.e5,    & ! hPa
  tcrit   = 258.,    & ! K
  g       = MAPL_GRAV,& ! m s-2
  cpor    = cp/rgas, &
  xlv     = 2.5e6,   & ! J kg-1
  akmin   = 1.0,     & ! #
  tkmin   = 1.e-5,   & ! m+2 s-2
  ccnclean= 250.,    & ! # cm-3
  c1      = 0.001,   & !
  T_0     = 273.16,  & ! K
  T_ice   = 250.16,  & ! K
  xlf     = 0.333e6, & ! latent heat of freezing (J K-1 kg-1)
  c0_mid  = 0.002,   & ! conversion rate (cloud to rain, m-1)
  c0_deep = 0.002,   &
  qrc_crit= 2.e-4      ! 0.1e-3 kg/kg

 !- proportionality constant to estimate pressure
 !- gradient of updraft (Zhang and Wu, 2003, JAS)
 REAL, PARAMETER ::    pgcd=1., pgcon=0.

 !- numerical constraints
 REAL, PARAMETER ::       &
  xmbmaxshal  =  0.05,    &  ! kg/m2/s
  mintracer   =  tiny(1.),&  ! kg/kg - tiny(x)
  smallerQV   =  1.e-16      ! kg/kg  

 INTEGER, PARAMETER :: MAX_NSPEC=200
 INTEGER, DIMENSION(MAX_NSPEC) :: ind_chem
 CHARACTER(len=100),DIMENSION(MAX_NSPEC)    ::  CHEM_NAME
 INTEGER           ,DIMENSION(MAX_NSPEC)    ::  CHEM_NAME_MASK,CHEM_NAME_MASK_EVAP
 INTEGER :: ispc_CO
 TYPE Hcts_vars
   REAL :: hstar,dhr,ak0,dak
 END TYPE Hcts_vars
 TYPE (Hcts_vars), ALLOCATABLE :: Hcts(:)
 
 CHARACTER(LEN=10),PARAMETER  :: host_model = 'NEW_GEOS5'
!
!--for GATE soundings only
!CHARACTER(LEN=10),PARAMETER  :: host_model = 'OLD_GEOS5'

 INTEGER :: whoami_all

!---------------------------------------------------------------------------------------------------

 interface 
    module SUBROUTINE GF_GEOS5_INTERFACE(mxp,myp,mzp,mtp,ITRCR,LONS,LATS,dt_moist          &
                               ,T, PLE, PLO, ZLE, ZLO, PK,  U, V, OMEGA           &
                               ,TH1, Q1, U1, V1 ,QLCN ,QICN,QLLS,QILS, CNPCPRATE  &
                               ,CNV_MF0, CNV_PRC3, CNV_MFD, CNV_DQLDT ,ENTLAM     &
                               ,CNV_MFC, CNV_UPDF, CNV_CVW, CNV_QC, CLCN          &
                               ,QV_DYN_IN,PLE_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                               ,RADSW   ,RADLW ,DQDT_BL  ,DTDT_BL                 &
                               ,FRLAND  ,AREA  ,USTAR ,TSTAR ,QSTAR ,T2M ,Q2M     &
                               ,TA      ,QA    ,SH    ,EVAP  ,PHIS                &
                               ,KPBLIN         &
                               ,MAPL_GRAV      &
                               ,STOCHASTIC_SIG, SIGMA_DEEP, SIGMA_MID             &
                               ,DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD                 &
                               ,MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD                  &
                               ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC      &
                               ,DTDTDYN,DQVDTDYN                                  &
                               ,NCPL, NCPI, CNV_NICE, CNV_NDROP, CNV_FICE, CLDMICRO &
                               ,qcrit, C0_auto &
                               ,TRACER,FSCAV,CNAMES,QNAMES,DTRDT_GF               &
                               ,RSU_CN,REV_CN, PFI_CN, PFL_CN                     )

    IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN) :: CLDMICRO !set two-moment microphysics

    INTEGER ,INTENT(IN) :: mxp,myp,mzp,mtp,ITRCR

    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(IN)   :: PLE,ZLE,PLE_DYN_IN
    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(IN)   :: T ,U ,V ,ZLO ,PLO ,PK ,OMEGA         &
                                                     ,RADSW  ,RADLW  ,DQDT_BL  ,DTDT_BL    &
                                                     ,QV_DYN_IN,U_DYN_IN,V_DYN_IN,T_DYN_IN   &
                                                     ,DTDTDYN,DQVDTDYN

    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(INOUT):: TH1,Q1,U1,V1,QLCN ,QICN, NCPL, NCPI
    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(OUT)  :: CNV_MFC
    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(OUT)  :: CNV_MF0 ,CNV_PRC3,CNV_MFD,CNV_DQLDT  &
                                                     ,CNV_UPDF, CNV_CVW, CNV_QC,CLCN,ENTLAM&
                                                     ,QLLS,QILS,  CNV_NICE, CNV_NDROP, CNV_FICE
    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(INOUT):: RSU_CN,REV_CN
    REAL   ,DIMENSION(mxp,myp,0:mzp) ,INTENT(INOUT):: PFI_CN, PFL_CN

    !-for debug purposes
    REAL   ,DIMENSION(mxp,myp,mzp)   ,INTENT(INOUT)  :: DQDT_GF,DTDT_GF,MUPDP,MUPSH,MUPMD ,DTRDT_GF
    REAL   ,DIMENSION(mxp,myp)       ,INTENT(INOUT)  :: MFDP,MFSH,MFMD,ERRDP,ERRSH,ERRMD
    REAL   ,DIMENSION(mxp,myp)       ,INTENT(INOUT)  :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC

    REAL   ,DIMENSION(MXP,MYP)       ,INTENT(IN)   :: FRLAND ,AREA ,USTAR ,TSTAR ,QSTAR &
                                                     ,T2M ,Q2M ,TA ,QA ,SH ,EVAP ,PHIS  &
                                                     ,KPBLIN,LONS,LATS                  &
                                                     ,STOCHASTIC_SIG
    REAL   ,DIMENSION(MXP,MYP)       ,INTENT(OUT)  :: SIGMA_DEEP, SIGMA_MID
    REAL   ,DIMENSION(MXP,MYP)       ,INTENT(OUT)  :: CNPCPRATE ! kg m-2 s-1

    REAL                             ,INTENT(IN)   :: DT_moist ,MAPL_GRAV, qcrit, c0_auto

    REAL   ,DIMENSION(mxp,myp,mzp,itrcr) ,INTENT(INOUT)   :: TRACER  !=XHO in grid_moist_comp.f90
    REAL   ,DIMENSION(itrcr)             ,INTENT(IN   )   :: FSCAV
    CHARACTER(len=*)  ,DIMENSION(mtp)    ,INTENT(IN   )   :: CNAMES,QNAMES

END SUBROUTINE GF_GEOS5_INTERFACE
end interface 
!---------------------------------------------------------------------------------------------------

interface
   module SUBROUTINE GF_GEOS5_DRV(mxp,myp,mzp,mtp             &
              ,ims,ime, jms,jme, kms,kme              &
              ,its,ite, jts,jte, kts,kte              &
	      ,flip                                   &
              ,FSCAV                                  &
              ,mynum                 &
              ,dt                    &
              ,dx2d                  &
              ,stochastic_sig        &
              ,zm                    &
              ,zt                    &
	      ,dm                    &

              ,lons                  &
              ,lats                  &
              ,aot500                &
              ,temp2m                &
              ,sflux_r               &
              ,sflux_t               &
              ,topt                  &
              ,xland                 &
              ,sfc_press             &
              ,kpbl                  &

              ,u                     &
              ,v                     &
              ,w                     &
              ,temp                  &
              ,press                 &
              ,rvap                  &
              ,curr_rvap             &
              ,TRACER                &!-note: uses GEOS-5 data structure

              ,rthften               &!gsf_t
              ,rqvften               &!gsf_q
              ,rthblten              &!sgsf_t
              ,rqvblten              &!sgsf_q
              !---- output ----
              ,CONPRR                &
              ,rthcuten              &
              ,rqvcuten              &
              ,rqccuten              &
              ,rucuten               &
              ,rvcuten               &
              ,rchemcuten            &
              ,revsu_gf              & 
              ,prfil_gf              & 
              !
              ,do_this_column        &
              ,ierr4d                &
              ,jmin4d                &
              ,klcl4d                &
              ,k224d                 &
              ,kbcon4d               &
              ,ktop4d                &
              ,kstabi4d              &
              ,kstabm4d              &
              ,cprr4d                &
              ,xmb4d                 &
              ,edt4d                 &
              ,pwav4d                &
              ,sigma4d               &
              ,pcup5d                &
              ,up_massentr5d         &
              ,up_massdetr5d         &
              ,dd_massentr5d         &
              ,dd_massdetr5d         &
              ,zup5d                 &
              ,zdn5d                 &
              ,prup5d                &
              ,prdn5d                &
              ,clwup5d               &
              ,tup5d                 &
              ,conv_cld_fr5d         &
              !-- for debug/diagnostic
              ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC &
              )

   IMPLICIT NONE
   INTEGER, PARAMETER,DIMENSION(maxiens) :: use_excess=(/1,1,1/)
  !------------------------------------------------------------------------
   INTEGER, INTENT(IN) :: ims,ime, jms,jme, kms,kme,    &
                          its,ite, jts,jte, kts,kte,    &
                          mynum,mzp,mxp,myp,mtp

   REAL,    INTENT(IN) :: DT

   INTEGER, INTENT(IN), dimension(mzp) :: flip
   
   REAL:: FSCAV(mtp)

   REAL,    DIMENSION(kts:kte,its:ite,jts:jte), INTENT(IN)  :: &
                                                          zm,   &
                                                          zt,   &
                                                          u,    &
                                                          v,    &
                                                          w,    &
                                                          rvap, &
                                                          temp, &
                                                          press,&
							  dm,   &
							  curr_rvap

   INTEGER, DIMENSION(its:ite,jts:jte), INTENT(IN) :: kpbl
   REAL,    DIMENSION(its:ite,jts:jte), INTENT(IN) :: topt ,aot500 ,temp2m ,sfc_press &
                                                     ,sflux_r ,sflux_t ,xland,lons,lats,dx2d &
                                                     ,stochastic_sig
   REAL,    DIMENSION(kts:kte,its:ite,jts:jte), INTENT(IN) ::        &
                                                         rthften  &
                                                        ,rqvften  &
                                                        ,rthblten &
                                                        ,rqvblten

   REAL,    DIMENSION(its:ite,jts:jte),         INTENT(OUT) ::   CONPRR
   REAL,    DIMENSION(kts:kte,its:ite,jts:jte), INTENT(OUT) ::    &
                                                    rthcuten   &
                                                   ,rqvcuten   &
                                                   ,rqccuten   &
                                                   ,rucuten    &
                                                   ,rvcuten    &
						   ,revsu_gf   &
						   ,prfil_gf   
						   
   !-***** TRACER has different data structure   (i,j,k,ispc) *********
   REAL,    DIMENSION(its:ite,jts:jte,kts:kte,mtp), INTENT(IN )  :: TRACER
   !-***** rchemcuten uses the GF data structure (ispc,k,i,j) *********
   REAL,    DIMENSION(mtp,kts:kte,its:ite,jts:jte), INTENT(OUT)  :: rchemcuten
 
   INTEGER, DIMENSION(its:ite,jts:jte), INTENT(INOUT) :: do_this_column

!- for convective transport and cloud/radiation (OUT)
   INTEGER,DIMENSION(mxp,myp,maxiens)::  &
!  integer, dimension(its:ite,jts:jte,maxiens) , INTENT(OUT) ::    &
               ierr4d                    &
              ,jmin4d                    &
              ,klcl4d                    &
              ,k224d                     &
              ,kbcon4d                   &
              ,ktop4d                    &
              ,kstabi4d                  &
              ,kstabm4d

   REAL,DIMENSION(mxp,myp,maxiens)::     &
!   real,dimension(its:ite,jts:jte,maxiens)    , INTENT(OUT) ::    &
               cprr4d                    &
              ,xmb4d                     &
              ,edt4d                     &
              ,pwav4d                    &
              ,sigma4d
   REAL,DIMENSION(mxp,myp,mzp,maxiens):: &
!   real,dimension(its:ite,jts:jte,kts:kte,maxiens), INTENT(OUT) ::    &
               pcup5d                    &
              ,up_massentr5d             &
              ,up_massdetr5d             &
              ,dd_massentr5d             &
              ,dd_massdetr5d             &
              ,zup5d                     &
              ,zdn5d                     &
              ,prup5d                    &
              ,prdn5d                    &
              ,clwup5d                   &
              ,tup5d                     &
              ,conv_cld_fr5d
!--for debug
   REAL   ,DIMENSION(mxp,myp)  ,INTENT(INOUT)  :: AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_EC

  END SUBROUTINE GF_GEOS5_DRV
end interface 
END MODULE ConvPar_GF_GEOS5
